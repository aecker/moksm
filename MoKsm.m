function model = MoKsm(data, varargin)
% MoK to 12d data (adapted from Calabrese & Paninski 2011)
% AE 2012-02-03

params.MaxTrainSpikes = 5000;
params.MaxTestSpikes = 20000;
params.TrainFrac = 0.8;
params.verbose = false;
params = parseVarArgs(params,varargin{:});

if nargin < 3, verbose = false; end
randn('state', 1)
rand('state', 1)

% warning off MATLAB:nearlySingularMatrix

tol = 0.0002;

t = data.SpikeTimes.data;
Y = data.Features.data;

if size(Y,1) == length(t)
    Y = Y';
elseif size(Y,2) ~= length(t)
    error('Time dimension doesn''t match dataset');
end

t = t(:)';

T = size(Y, 2);

assert(size(Y,1) <= 50, 'Dimensionality way too high');

% split into training & test data
rnd = randperm(T);
nTrain = fix(params.TrainFrac * T);
train = sort(rnd(1 : min(nTrain,params.MaxTrainSpikes)));
test = sort(rnd(nTrain + 1 : min(end, nTrain + params.MaxTestSpikes)));
Ytrain = Y(:, train);
Ytest = Y(:, test);
ttrain = t(train);
ttest = t(test);

nTrain = length(train);
T = length(train) + length(test);

% Initialize model using 1 component
fprintf('Running initial Kalman filter model with one cluster ')
model.mu = repmat(mean(Ytrain, 2), [1 nTrain]); % cluster means
model.C = cov(Ytrain');                                  % observation covariance
model.Cmu = diag(median(abs(Ytrain), 2) / 0.6745 / 1000);  % cluster mean drift [TODO: make it dependent on average firing rate]
model.priors = 1;                       % cluster weights
model.post = ones(1, nTrain);
model.pk = mvn(Ytrain - model.mu, model.C + model.Cmu);
model.logLike = sum(mylog(model.pk));
model = fullEM(Ytrain, model, tol, params.verbose);
if params.verbose, plotData(model,Ytrain,ttrain), end
fprintf(' done (likelihood: %.5g)\n', model.logLike(end))

% Run split & merge 
% We alternate between trying to split and trying to merge. Both are done
% until no candidate leads to success. If both splitting and merging didn't
% lead to success we terminate
op = {@trySplit, @tryMerge};
i = 1;
success = true(1, 2);
while any(success)
    [model, success(i)] = op{i}(Ytrain, Ytest, ttest, ttrain, model, tol, params.verbose);
    if ~success(i)
        i = 3 - i;
    end
end
 
fprintf('Done with split & merge\n')

fprintf('Performing fit on entire dataset ...');
[logLike p] = evalTestSet(Y, t, ttrain, model)
model.postAll = p;
fprintf(' Done\n');

fprintf('--\n')
fprintf('Number of clusters: %d\n', size(model.mu, 3))
fprintf('Log-likelihoods\n')
fprintf('  training set: %.8g\n', model.logLike(end))
fprintf('      test set: %.8g\n', evalTestSet(Ytest, ttest, ttrain, model))
fprintf('\n\n')

model.Ytrain = Ytrain;
model.Ytest = Ytest;
model.ttrain = ttrain;
model.ttest = ttest;


   
    
function [model, success] = tryMerge(Ytrain, Ytest, ttest, ttrain, model, tol, verbose)

success = false;
cands = getMergeCandidates(model.post);
logLikeTest = evalTestSet(Ytest, ttest, ttrain, model);
for ij = cands'
    fprintf('Trying to merge clusters %d and %d ', ij(1), ij(2))
    newModel = mergeClusters(model, ij(1), ij(2));
    newModel = fullE(Ytrain, newModel);
    newModel = fullEM(Ytrain, newModel, tol, verbose);
    newLogLikeTest = evalTestSet(Ytest, ttest, ttrain, newModel);
    if newLogLikeTest > logLikeTest
        fprintf(' success (likelihood improved by %.5g)\n', newLogLikeTest - logLikeTest)
        if verbose, plotData(newModel,Ytrain,ttrain), end
        model = newModel;
        success = true;
        break
    else
        fprintf(' aborted\n')
    end
end


function [model, success] = trySplit(Ytrain, Ytest, ttest, ttrain, model, tol, verbose)

success = false;
splitCands = getSplitCandidates(model.post, model.pk, model.priors);
logLikeTest = evalTestSet(Ytest, ttest, ttrain, model);
for i = splitCands'
    try
        fprintf('Trying to split cluster %d ', i)
        newModel = splitCluster(model, i);
        newModel = fullE(Ytrain, newModel);
        newModel = fullEM(Ytrain, newModel, tol, verbose);
        newLogLikeTest = evalTestSet(Ytest, ttest, ttrain, newModel);
        if newLogLikeTest > logLikeTest
            fprintf(' success (likelihood improved by %.5g)\n', newLogLikeTest - logLikeTest)
            if verbose, plotData(newModel,Ytrain,ttrain), end
            model = newModel;
            success = true;
            break
        else
            fprintf(' aborted\n')
        end
    catch err
        if strcmp(err.identifier, 'MoKsm:starvation')
            fprintf(' aborted due to component starvation\n')
        else
            retrhow(err)
        end
    end    
end



function model = fullEM(Y, model, tol, verbose)
% EM recursion

[mu, C, Cmu, priors, post, pk, logLike] = expand(model);
[D, T, K] = size(mu);
Cf = zeros(D, D, T);                            % state covariances
iCfCmu = zeros(D, D, T);
iter = 0;
logLikeBase = logLike(end);
while iter < 2 || (logLike(end) - logLike(end - 1)) / (logLike(end - 1) - logLikeBase) > tol

    if ~mod(iter, 10)
        fprintf('.')
    end
    iter = iter + 1;
    
    for k = 1 : K
        
        muk = mu(:, :, k);
        Ck = C(:, :, k);
    
        % Forward step for updating the means (Eq. 9)
        Cf(:, :, 1) = Ck;
        iCk = inv(Ck);
        for t = 2 : T
            piCk = post(k, t - 1) * iCk; %#ok
            iCfCmu(:, :, t - 1) = inv(Cf(:, :, t - 1) + Cmu);
            Cf(:, :, t) = inv(iCfCmu(:, :, t - 1) + piCk);
            muk(:, t) = Cf(:, :, t) * (iCfCmu(:, :, t - 1) * muk(:, t - 1) + piCk * Y(:, t - 1));
        end

        % Backward step for updating the means (Eq. 10)
        for t = T - 1 : -1 : 1
            muk(:, t) = muk(:, t) + Cf(:, :, t) * (iCfCmu(:, :, t) * (muk(:, t + 1) - muk(:, t)));
        end

        % Update observation covariance (Eq. 11)
        Ymu = Y - muk;
        Ck = (bsxfun(@times, post(k, :), Ymu) * Ymu') / sum(post(k, :));

        % Estimate (unnormalized) probabilities
        pk(k, :) = mvn(Ymu, Ck + Cmu);
        post(k, :) = priors(k) * pk(k, :);
        
        mu(:, :, k) = muk;
        C(:, :, k) = Ck;
    end
    
    % calculate log-likelihood
    p = sum(post, 1);
    logLike(end + 1) = sum(mylog(p)); %#ok
    if verbose
        figure(1)
        plot(logLike, '.-k')
        ylim(prctile(logLike, [10 100]))
        drawnow
    end
    
    % normalize probabilities
    post = bsxfun(@rdivide, post, p);
    post(:, p == 0) = 0;
    
    % update class priors
    priors = sum(post, 2) / T;
    
    % check for starvation
    if any(priors * T < 2 * D)
        error('MoKsm:starvation', 'Component starvation: cluster %d', find(priors * T < 2 * D, 1))
    end
end

model = collect(mu, C, Cmu, priors, post, pk, logLike);



function model = fullE(Y, model)
% Do one full E-step

[~, T, K] = size(model.mu);
for k = 1 : K
    model.pk(k, :) = mvn(Y - model.mu(:, :, k), model.C(:, :, k) + model.Cmu);
    model.post(k, :) = model.priors(k) * model.pk(k, :);
end
p = sum(model.post, 1);
model.logLike(end + 1) = sum(mylog(p));
model.post = bsxfun(@rdivide, model.post, p);
model.post(:, p == 0) = 0;
model.priors = sum(model.post, 2) / T;



function [logLike p] = evalTestSet(Ytest, ttest, ttrain, model)
% Evaluate log-likelihood on test set by interpolating cluster means from
% training set.

ttrain(1) = min(ttest(1), ttrain(1));
ttrain(end) = max(ttest(end), ttrain(end));
K = size(model.mu, 3);
T = numel(ttest);
p = zeros(K, T);
for k = 1 : K
    muk = interp1(ttrain, model.mu(:, :, k)', ttest)';
    Ymu = Ytest - muk;
    p(k, :) = model.priors(k) * mvn(Ymu, model.C(:, :, k) + model.Cmu);
end
logLike = sum(mylog(sum(p, 1)));


function model = splitCluster(model, k)
% Split cluster k

[mu, C, Cmu, priors, post, pk, logLike] = expand(model);
[D, ~, K] = size(mu);
deltaMu  = chol(C(:, :, k))' * randn(D, 1);
mu(:, :, k) = bsxfun(@plus, mu(:, :, k), deltaMu);
mu(:, :, K + 1) = bsxfun(@minus, mu(:, :, k), deltaMu);
C(:, :, k) = det(C(:, :, k))^(1 / D) * eye(D);
C(:, :, K + 1) = C(:, :, k);
priors(k) = priors(k) / 2;
priors(K + 1) = priors(k);
model = collect(mu, C, Cmu, priors, post, pk, logLike);


function model = mergeClusters(model, i, j)
% Merge clusters i and j

[mu, C, Cmu, priors, post, pk, logLike] = expand(model);
mu(:, :, i) = (priors(i) * mu(:, :, i) + priors(j) * mu(:, :, j)) / (priors(i) + priors(j));
C(:, :, i) = (priors(i) * C(:, :, i) + priors(j) * C(:, :, j)) / (priors(i) + priors(j));
priors(i) = priors(i) + priors(j);
mu(:, :, j) = [];
C(:, :, j) = [];
priors(j) = [];
post(j, :) = [];
pk(j, :) = [];
model = collect(mu, C, Cmu, priors, post, pk, logLike);


function cand = getSplitCandidates(post, pk, priors)

fk = bsxfun(@rdivide, post, sum(post, 2));
Jsplit = sum(fk .* (mylog(fk) - mylog(pk)), 2);
[~, cand] = sort(Jsplit, 'descend');
[D, T] = size(post);
cand = cand(priors(cand) * T > 4 * D); % don't split small clusters


function cand = getMergeCandidates(post)

K = size(post, 1);
maxCandidates = ceil(K * sqrt(K) / 2);
np = sqrt(sum(post .* post, 2));
Jmerge = zeros(K * (K - 1) / 2, 1);
cand = zeros(K * (K - 1) / 2, 2);
k = 0;
for i = 1:K
    for j = i+1:K
        k = k + 1;
        Jmerge(k) = post(i, :) * post(j, :)' / (np(i) * np(j));
        cand(k, :) = [i j];
    end
end
[~, order] = sort(Jmerge, 'descend');
cand = cand(order(1:min(end, maxCandidates)), :);


function [mu, C, Cmu, priors, post, pk, logLike] = expand(model)

mu = model.mu;
C = model.C;
Cmu = model.Cmu;
priors = model.priors;
post = model.post;
pk = model.pk;
logLike = model.logLike;


function model = collect(mu, C, Cmu, priors, post, pk, logLike)

model.mu = mu;
model.C = C;
model.Cmu = Cmu;
model.priors = priors;
model.post = post;
model.pk = pk;
model.logLike = logLike;


function y = mylog(x)
% Natural logarithm excluding zeros

y = reallog(x);
y(x == 0) = 0;



function p = mvn(X, C)
% Zero-mean multivariate normal probability density
%   p = clus_mvn(X, C) calculates the density of the multivariate normal
%   distribution with zero mean and covariance matrix C at X. X is assumed
%   to be a column vector or a matrix of multiple column vectors, in which
%   case the result, p, is a row vector.
%
% AE 2012-02-05

D = size(C, 1);
const = (2*pi)^(-D/2);
[Ch, ~] = chol(C);
p = const / prod(diag(Ch)) * exp(-1/2 * sum((Ch' \ X).^2, 1));




function plotData(model,Y,t)

if nargin == 3, d1 = 2; d2 = 1; end
[~, j] = max(model.post, [], 1);
K = size(model.post, 1);
c = lines;
figure(2), clf, hold all
hdl = zeros(1, K);

subplot(211)
cla
hold on
for i = 1:K
    plot(Y(d1, j == i), Y(d2, j == i), '.', 'markersize', 1, 'color', c(i, :))
    hdl(i) = plot(model.mu(d1, :, i), model.mu(d2, :, i), '*-', 'color', c(i, :));
end
legend(hdl, arrayfun(@(x) sprintf('Cluster %d', x), 1:K, 'UniformOutput', false))

subplot(212)
cla
hold on
for i = 1:K
    plot(t(j==i),Y(d1, j == i), '.', 'markersize', 1, 'color', c(i, :))
    hdl(i) = plot(t,model.mu(d1, :, i), '*-', 'color', c(i, :));
end
legend(hdl, arrayfun(@(x) sprintf('Cluster %d', x), 1:K, 'UniformOutput', false))

