% MoKsm test harness for t-distribution parameter estimation

%% Example 1: univariate, single-component, stationary

% Eliminate the regularizing ridge for this example
Df = 3;
DTmu = 60e3;
CovRidge = 0;

% Generate samples with mu = 0, Sigma = 1, and all in the same time bin
N = 1e4;
Y = trnd(Df, 1, N);
t = (DTmu / 2) * ones(1,N);

% Plotting parameters
xlimits = [-5 5];
legend_fmt = '%s: \\mu = %0.2f, \\Sigma = %0.2f, \\nu = %0.2f';

% Plot the histogram and the generating PDF
figure
% histogram
bins = linspace(xlimits(1), xlimits(2), 101)';
binwidth = mean(diff(bins));
bins = [bins(1) - binwidth; bins; bins(end) + binwidth];
counts = hist(Y', bins);
h_hist = bar(bins, counts/(N * binwidth), 1);
str_hist = 'Synthesized data';
set(h_hist,'FaceColor',[.7 .7 .7],'EdgeColor','none');
xlim(xlimits);
hold on
title('1-component, 1-dimensional, stationary, CovRidge=0');
% true PDF
x = linspace(xlimits(1), xlimits(2), 150)';
pdf_true = tpdf(x, Df);
h_true = plot(x, pdf_true, 'k-');
str_true = sprintf(legend_fmt, 'Truth', 0, 1, Df);

% Fit using stats toolbox
pd = fitdist(Y', 'tLocationScale');
pdf_fitdist = pd.pdf( x );
h_fitdist = plot(x, pdf_fitdist, 'g-');
str_fitdist = sprintf(legend_fmt, 'fitdist', pd.mu, pd.sigma, pd.nu);

% Fit using MoKsm
model = MoKsm('Df',Df, 'DTmu',DTmu, 'CovRidge',CovRidge);
model = model.fit(Y, t);
% Initial fit only does 1 EM step before split-and-merge; use refit() to
% run until convergerce
fprintf('Calling refit() to run EM until convergence');
model = model.refit();
fprintf('done\n');
% make sure we only got the one component
assert(size(model.mu,3)==1, 'Multiple components fitted');
% compute and plot the PDF
pdf_moksm = MoKsm.mvt( x' - model.mu, model.C, model.df);
h_moksm = plot(x, pdf_moksm, 'r-');
str_moksm = sprintf(legend_fmt, 'MoKsm', model.mu, model.C, model.df);

% Create legend
h_legend = legend([h_true h_hist h_fitdist h_moksm], ...
    str_true, str_hist, str_fitdist, str_moksm);
% Make sure the legend doesn't overlap the data
heights = [0 1 0 0 0 -1 0 0; 0 0 0 0 0 0 0 1] ...
    * [get(h_legend, 'Position'), get(gca,'Position')]';
ylim(ylim() * heights(2)/heights(1));


%% Example 2: something more realistic

% Ground truth
Df = 3;
n_t = 20;
DTmu = 60e3;
Tmu = (0:n_t) * DTmu; % length = n_t+1 for the last fencepost
mu_1 = [linspace(100, 200, n_t) ; linspace(300, 400, n_t)];
mu_2 = repmat([400; 100], [1 n_t]);
C_1 = [30 -15; -15 30]^2;
C_2 = [30 5; 5 30]^2;
N_1 = 2000;
N_2 = 1000;
ndims = 2;

% Generate points
mu = cat(3, mu_1, mu_2);
C = cat(3, C_1, C_2);
N = cat(2, N_1, N_2);
n_cl = size(mu,3);
t_arr = cell(n_cl,1);
Y_arr = cell(n_cl,1);
for k = 1:n_cl
    % Time
    t = linspace(1, Tmu(end)-1, N(k));
    % Spikes
    % mvtrnd is a bit too paternalistic; it doesn't let you specify a scale
    % matrix (probably because it's worried about the user confusing the
    % scale matrix for the expected covariance).
    R = chol(C(:,:,k));
    u = chi2rnd(Df/2, [N(k) 1]);
    X = bsxfun(@rdivide, randn(N(k), ndims)*R, sqrt(u)); % row vectors
    Y = X' + mu(:,ceil(t/DTmu),k);
    % Collect and go on to the next component
    t_arr{k} = t;
    Y_arr{k} = Y;
end
% Combine the two sources
t = cat(2, t_arr{:});
Y = cat(2, Y_arr{:});
[t, sortperm] = sort(t);
Y = Y(:,sortperm);

% Figure parameters
n_skip = 5;
xlimits = [0 500];
ylimits = [0 500];

% Plot the raw data and the training points
figure
subplot(1,2,1)
% Synthetic data
h_data = plot(Y(1,:), Y(2,:), 'k.', 'MarkerSize',2);
str_data = 'Data points';
xlim(xlimits);
ylim(ylimits);
title(sprintf('%d non-stationary clusters, df=%.1f', n_cl, Df));
hold on
% Ground truth (ellipses)
theta = linspace(0, 2*pi, 100);
unit_circle = [sin(theta); cos(theta)];
for k = 1:n_cl
    % For each component
    R = chol(C(:,:,k));
    for j = ceil(n_skip/2):n_skip:n_t
        % For each selected time step
        ellipse = bsxfun(@plus, R' * unit_circle, mu(:,j,k));
        h_truth = plot(ellipse(1,:), ellipse(2,:), 'k-');
        % Go on to the next time step
    end
    % Go on to the next component
end
str_truth = 'Ground truth';

% Fit using MoKsm
model = MoKsm('Df',Df, 'DTmu',DTmu, 'ClusterCost',.05);
model = model.fit(Y, t);
% Plot
n_cl_moksm = size(model.mu, 3);
colors = lines(n_cl_moksm);
h_moksm2 = zeros(1,n_cl_moksm);
str_moksm2 = cell(1,n_cl_moksm);
for k = 1:size(model.mu,3)
    % For each component
    R = chol(model.C(:,:,k));
    for j = ceil(n_skip/2):n_skip:size(model.mu,2)
        % For each selected time step
        ellipse = bsxfun(@plus, R' * unit_circle, model.mu(:,j,k));
        h_moksm2(k) = plot(ellipse(1,:), ellipse(2,:), '-', ...
            'Color',colors(k,:));
        str_moksm2{k} = sprintf('MoKsm cl %d',k);
        % Go on to the next time step
    end
    % Go on to the next component
end

% Legend
h_legend = legend([h_truth, h_data, h_moksm2], ...
    str_truth, str_data, str_moksm2{:});

% Use the MoKsm model to generate data
subplot(1,2,2);
% Generate points
N_moksm = round(sum(N) * model.priors);
for k = 1:n_cl_moksm
    % Time
    t = linspace(1, model.mu_t(end)-1, N_moksm(k));
    % Spikes
    R = chol(model.C(:,:,k));
    u = chi2rnd(model.df/2, [N_moksm(k) 1]);
    X = bsxfun(@rdivide, randn(N_moksm(k), ndims)*R, sqrt(u));
    Y = X' + model.mu(:,ceil(t/DTmu),k);
    % Plot
    plot(Y(1,:), Y(2,:), '.', 'MarkerSize',2, 'Color',colors(k,:));
    hold on
end
% Make it pretty
xlim(xlimits);
ylim(ylimits);
title('Points generated from Moksm-fitted model');

