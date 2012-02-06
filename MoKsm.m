function model = MoKsm(Y, t, K)
% MoK to 12d data (adapted from Calabrese & Paninski 2011)
% AE 2012-02-03

randn('state', 1)
rand('state', 1)

tol = 0.0002;
trainFraction = 0.8;
T = size(Y, 2);

% split into training & test data
rnd = randperm(T);
nTrain = fix(trainFraction * T);
train = sort(rnd(1 : nTrain));
test = sort(rnd(nTrain + 1 : end));
Ytrain = Y(:, train);
Ytest = Y(:, test);
ttrain = t(train);
ttest = t(test);

% Model initialization using MoG with 5 components
fprintf('Initializing using mixture of Gaussians with %d clusters\n', K)
mog = gmdistribution.fit(Ytrain', K, 'options', statset('MaxIter', 500));
model.mu = repmat(permute(mog.mu, [2 3 1]), [1 nTrain 1]); % cluster means
model.C = mog.Sigma;                                  % observation covariance
model.Cmu = diag(median(abs(Ytrain), 2) / 0.6745 / 1000);  % cluster mean drift [TODO: make it dependent on average firing rate]
model.priors = mog.PComponents;                       % cluster weights
model.post = mog.posterior(Ytrain')';
model.pk = zeros(size(model.post));
model.logLike = -mog.NlogL;

% run full EM
fprintf('Running initial mixture of Kalman filters with %d clusters ', K)
model = fullEM(Ytrain, model, tol);
plotData(Ytrain, model.mu, model.post, 1, 7)
fprintf(' done (likelihood: %.5g)\n', model.logLike(end))

% try split & merge
done = false;
while ~done
    
    merged = false;
    mergeCands = getMergeCandidates(model.post);
    logLikeTest = evalTestSet(Ytest, ttest, ttrain, model);
    for ij = mergeCands'
        fprintf('Trying to merge clusters %d and %d ', ij(1), ij(2))
        newModel = mergeClusters(model, ij(1), ij(2));
        newModel = fullE(Ytrain, newModel);
        newModel = fullEM(Ytrain, newModel, tol);
        newLogLikeTest = evalTestSet(Ytest, ttest, ttrain, newModel);
        if newLogLikeTest > logLikeTest
            fprintf(' success (likelihood improved by %.5g)\n', newLogLikeTest - logLikeTest)
            plotData(Ytrain, newModel.mu, newModel.post, 1, 7)
            model = newModel;
            merged = true;
            break
        else
            fprintf(' aborted\n')
        end
    end
    
    split = false;
    splitCands = getSplitCandidates(model.post, model.pk, model.priors);
    logLikeTest = evalTestSet(Ytest, ttest, ttrain, model);
    for i = splitCands'
        fprintf('Trying to split cluster %d ', i)
        newModel = splitCluster(model, i);
        newModel = fullE(Ytrain, newModel);
        newModel = fullEM(Ytrain, newModel, tol);
        newLogLikeTest = evalTestSet(Ytest, ttest, ttrain, newModel);
        if newLogLikeTest > logLikeTest
            fprintf(' success (likelihood improved by %.5g)\n', newLogLikeTest - logLikeTest)
            plotData(Ytrain, newModel.mu, newModel.post, 1, 7)
            model = newModel;
            split = true;
            break
        else
            fprintf(' aborted\n')
        end
    end
    
    done = ~merged && ~split;
end

fprintf('Done with split & merge\n')
fprintf('--\n')
fprintf('Number of clusters: %d\n', size(model.mu, 3))
fprintf('Log-likelihoods\n')
fprintf('  training set: %.8g\n', model.logLike(end))
fprintf('      test set: %.8g\n', logLikeTest)
fprintf('\n\n')

model.Ytrain = Ytrain;
model.Ytest = Ytest;
model.ttrain = ttrain;
model.ttest = ttest;





function model = fullEM(Y, model, tol)
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
    logLike(end + 1) = sum(reallog(p)); %#ok
    figure(1), plot(logLike, '.-k'), drawnow
    
    % normalize probabilities
    post = bsxfun(@rdivide, post, p);
    
    % update class priors
    priors = sum(post, 2) / T;
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
model.logLike(end + 1) = sum(reallog(p));
model.post = bsxfun(@rdivide, model.post, p);
model.priors = sum(model.post, 2) / T;



function logLike = evalTestSet(Ytest, ttest, ttrain, model)
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
logLike = sum(reallog(sum(p, 1)));


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
Jsplit = sum(fk .* (reallog(fk) - reallog(pk)), 2);
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




function plotData(Y, mu, p, d1, d2)

if nargin == 3, d1 = 1; d2 = 7; end
[~, j] = max(p);
K = size(p, 1);
c = lines;
figure(2), clf, hold all
hdl = zeros(1, K);
for i = 1:K
    plot(Y(d1, j == i), Y(d2, j == i), '.', 'markersize', 1, 'color', c(i, :))
    hdl(i) = plot(mu(d1, :, i), mu(d2, :, i), '*-', 'color', c(i, :));
end
legend(hdl, arrayfun(@(x) sprintf('Cluster %d', x), 1:K, 'UniformOutput', false))








