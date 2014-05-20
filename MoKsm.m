classdef MoKsm
    % Mixture of Kalman filters model adapted from Calabrese & Paninski
    % 2011.
    %
    % This class fits a Mixture of Kalman filters model using a split &
    % merge approach. We optimize a penalized average likelihood, where the
    % penalizer is a constant cost per cluster.
    %
    % To improve efficiency, the model is slightly modified so that it
    % tracks the cluster means in time blocks, each of which can contain
    % multiple spikes. Within each block the means are assumed to be
    % constant.
    %
    % In addition, we use a mixture of t distributions (with fixed degrees
    % of freedom) instead of Gaussians. Thanks to Kevin Shan from Caltech
    % for working out and implementing the correct update equations for the
    % case of t distributions!
    %
    % Furthermore, the covariance matrices are regularized by adding a
    % smnall ridge (diagonal matrix) during the M step.
    %
    % Alexander S. Ecker, R. James Cotton and Kevin Shan
    % 2014-05-20
    
    properties
        params      % parameters for fitting
        mu_t        % times at which cluster means are updated
        mu          % cluster means
        C           % cluster covariances
        Cmu         % covariance of cluster mean drift
        priors      % cluster priors (mixing proportions)
        df          % degress of freedom for t distribution
        Y           % data
        t           % times
        train       % indices of training set: Y(:, train), t(train)
        test        % indices of test set
        runtime     % run time of the fitting (sec)
    end
    
    properties (Access = private)
        like        % cache for model likelihoods
        post        % cache for posteriors
        logLike     % log-likelihood curve during fitting
        blockId     % mapping from data point to time blocks
        spikeId     % mapping from time blocks to data points (training data)
    end
    
    methods

        function self = MoKsm(varargin)
            % MoKsm constructor
            %   mok = MoKsm('param1', value1, 'param2', value2, ...)
            %   constructs a MoKsm object with the following optional
            %   parameters:
            %
            %   MaxTrainSpikes  max. number of spikes for training data
            %   MaxTestSpikes   max. number of spikes for test data
            %   TrainFrac       fraction of data points used for training
            %   Seed            initial seed for random number generator
            %   Df              degrees of freedom for t distribution
            %   DriftRate       drift rate per unit of time
            %   DTmu            block size for means in units of time
            %   Tolerance       tolerance for determining convergence
            %   Verbose         verbose output
            %   CovRidge        independent variance added to cluster covariances
            %   ClusterCost     penalizer for adding additional clusters
            
            % create from struct
            if nargin && isstruct(varargin{1}) && isfield(varargin{1}, 'params')
                s = varargin{1};
                self.params = s.params;
                self = self.initialize(s.Y, s.t, s.train, s.test, s.mu_t, s.mu, s.C, s.Cmu, s.priors, s.df);
                return
            end
            
            % parse optional parameters
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('MaxTrainSpikes', 50000);
            p.addOptional('MaxTestSpikes', 50000);
            p.addOptional('TrainFrac', 0.8);
            p.addOptional('Seed', 1);
            p.addOptional('Df', 2);
            p.addOptional('DriftRate', 400 / 3600 / 1000);
            p.addOptional('DTmu', 60 * 1000);
            p.addOptional('Tolerance', 0.0005);
            p.addOptional('Verbose', false);
            p.addOptional('CovRidge', 1.5);
            p.addOptional('ClusterCost', 0.0025);
            p.parse(varargin{:});
            self.params = p.Results;
        end
        
        
        function self = initialize(self, Y, t, train, test, mu_t, mu, C, Cmu, priors, df)
            % Initialize model
            %   self = initialize(self, Y, t, train, test, mu_t, mu, C, Cmu, priors, df)
            %
            %   Inputs (D: dimensions, N: samples, T: blocks, K: clusters)
            %     Y         data (D x N)
            %     t         times of data points (1 x N)
            %     train     indices of training data
            %     test      indices of test data
            %     mu_t      time blocks (bin edges, 1 x T+1)
            %     mu        cluster means (D x T x K)
            %     C         cluster covariances (D x D x K)
            %     Cmu       covariance of mean drift (D x D)
            %     priors    cluster priors (1 x K)
            %     df        degrees of freedom for t distribution (scalar)
            
            % model parameters
            self.mu_t = mu_t;
            self.mu = mu;
            self.C = C;
            self.Cmu = Cmu;
            self.priors = priors;
            self.df = df;
            
            % training & test data
            self.Y = Y;
            self.t = t;
            self.train = train;
            self.test = test;
            if ~isempty(Y)
                [~, self.blockId] = histc(t, self.mu_t);
                self.spikeId = arrayfun(@(x) find(self.blockId(self.train) == x), 1 : numel(mu_t) - 1, 'UniformOutput', false);
                self = self.updateCache();
            end
        end
        
        
        function self = fit(self, Y, t)
            % Fit the model
            %   self = fit(self, Y, t) fits the model to data (Y, t). Y is
            %   a matrix of size (#dimensions x #samples) and t is a vector
            %   of length #samples.
            %
            %   See MoKsm for optional parameters to use for fitting.
            
            self.runtime = now();
            
            % make sure dimensions of input are correct
            if size(Y, 1) == length(t)
                Y = Y';
            elseif size(Y, 2) ~= length(t)
                error('Time dimension doesn''t match dataset');
            end
            assert(size(Y, 1) <= 50, 'Dimensionality way too high');
            t = reshape(t, 1, []);
            
            % sort by time
            [t, order] = sort(t);
            Y = Y(:, order);

            % ensure deterministic behavior
            rng('default')
            rng(self.params.Seed);
            
            % split into training & test data
            T = size(Y,2);
            rnd = randperm(T);
            nTrain = fix(self.params.TrainFrac * T);
            train = sort(rnd(1 : min(nTrain, self.params.MaxTrainSpikes))); %#ok<*PROP>
            test = sort(rnd(nTrain + 1 : min(end, nTrain + self.params.MaxTestSpikes)));
            Ytrain = Y(:, train);
            
            % Assign spikes to time blocks
            mu_t = t(1) : self.params.DTmu : t(end) + self.params.DTmu;
            nTime = numel(mu_t) - 1;
            fprintf('Estimating cluster means at %d time points\n', nTime);

            % fit initial model with one component
            fprintf('Running initial Kalman filter model with one cluster ')
            mu = repmat(mean(Ytrain, 2), [1 nTime]);
            C = cov(Ytrain');
            Cmu = eye(size(mu, 1)) * self.params.DTmu * self.params.DriftRate;
            df = self.params.Df;
            self = self.initialize(Y, t, train, test, mu_t, mu, C, Cmu, 1, df);
            if isinf(df), maxIter = 1; else maxIter = Inf; end
            self = self.EM(maxIter);
            fprintf(' done (likelihood: %.5g)\n', self.logLikelihood())
            
            % Run split & merge
            % We alternate between trying to split and trying to merge. Both are done
            % until no candidate leads to success. If both splitting and merging didn't
            % lead to success we terminate
            op = {@trySplit, @tryMerge};
            i = 1;
            success = true(1, 2);
            while any(success)
                [self, success(i)] = op{i}(self);
                if success(i)
                    success(:) = true;
                else
                    i = 3 - i;
                end
            end
            
            % calculate likelihoods and posteriors
            fprintf('Calculating likelihoods and posteriors for entire dataset...')
            self = self.updateCache();
            fprintf(' done\n\n')
            
            fprintf('--\n')
            fprintf('Number of clusters: %d\n', size(self.mu, 3))
            fprintf('Log-likelihoods\n')
            fprintf('  training set: %.8g\n', self.logLikelihood(self.train))
            fprintf('      test set: %.8g\n', self.logLikelihood(self.test))
            
            self.runtime = (now() - self.runtime) * 24 * 60 * 60; % convert to sec
            fprintf('Total run time: %.1f sec.\n\n\n', self.runtime)
        end
        
        
        function self = refit(self)
            % Refit model.
            
            self = self.EM();
            self = self.updateCache();
        end
        
        
        function self = splitCluster(self, k)
            % Split cluster k.
            %   Randomly perturb cluster means followed by partial EM
            %   update.

            partial = self.getPartial(k);
            [mu, C, Cmu, priors, df] = partial.expand();
            D = size(mu, 1);
            deltaMu  = chol(C)' * randn(D, 1) * 0.2;
            mu(:, :, 1) = bsxfun(@plus, mu(:, :, 1), deltaMu);
            mu(:, :, 2) = bsxfun(@minus, mu(:, :, 1), deltaMu);
            C(:, :, 1) = det(C(:, :, 1))^(1 / D) * eye(D);
            C(:, :, 2) = C(:, :, 1);
            priors(1 : 2) = priors(1) / 2;
            partial = partial.collect(mu, C, Cmu, priors, df);
            partial = partial.EM(20);
            self.mu(:, :, [k end+1]) = partial.mu;
            self.C(:, :, [k end+1]) = partial.C;
            self.priors([k end+1]) = partial.priors * self.priors(k);
            self = self.updateCache();
        end
        
        
        function self = mergeClusters(self, ids)
            % Merge clusters by partial EM update.
            
            ids = sort(ids);
            partial = self.getPartial(ids);
            p = permute(partial.priors, [3 2 1]);
            partial.mu(:, :, 1) = sum(bsxfun(@times, partial.mu, p), 3) / sum(p);
            partial.C(:, :, 1) = sum(bsxfun(@times, partial.C, p), 3) / sum(p);
            partial.mu(:, :, 2 : end) = [];
            partial.C(:, :, 2 : end) = [];
            partial.priors = 1;
            partial = partial.EM(20);
            self.mu(:, :, ids(1)) = partial.mu;
            self.C(:, :, ids(1)) = partial.C;
            self.priors(ids(1)) = sum(self.priors(ids));
            self.mu(:, :, ids(2 : end)) = [];
            self.C(:, :, ids(2 : end)) = [];
            self.priors(ids(2 : end)) = [];
            self = self.updateCache();
        end
        
        
        function self = deleteCluster(self, k)
            % Delete cluster k.
            
            self.mu(:, :, k) = [];
            self.C(:, :, k) = [];
            self.priors(k) = [];
            self.priors = self.priors / sum(self.priors);
            self = self.updateCache();
        end
        
        
        function like = likelihood(self, index)
            % Likelihood of the data for each cluster
            %   p = likelihood(self) returns the likelihoods of the entire
            %   dataset.
            %
            %   p = likelihood(self, index) returns the likelihood for the
            %   given data points.
            
            if nargin == 1 && ~isempty(self.like)
                like = self.like;
            else
                if nargin == 1
                    index = ':';
                end
                [mu, C, ~, priors, df] = self.expand();
                K = numel(priors);
                like = zeros(K, size(self.Y(:, index), 2));
                for k = 1 : K
                    muk = mu(:, self.blockId(index), k);
                    like(k, :) = priors(k) * MoKsm.mixtureDistribution(self.Y(:, index) - muk, C(:, :, k), df);
                end
            end
        end
        
        
        function post = posterior(self, varargin)
            % Posterior of class membership for each spike
            %   post = posterior(self) returns the posteriors on the entire
            %   dataset.
            %   
            %   post = posterior(self, index) returns the posteriors for
            %   the given data points.
            
            if nargin == 1 && ~isempty(self.post)
                post = self.post;
            else
                like = self.likelihood(varargin{:});
                p = sum(like, 1);
                post = bsxfun(@rdivide, like, p);
                post(:, p == 0) = 0;
            end
        end
        
        
        function logl = logLikelihood(self, varargin)
            % Penalized average log-likelihood
            %   logl = logLikelihood(self) returns the penalized average
            %   log-likelihood of the entire dataset.
            %   
            %   logl = logLikelihood(self, index) returns the penalized
            %   average log-likelihood for the given subset.
            
            [D, ~, K] = size(self.mu);
            logl = mean(MoKsm.mylog(sum(self.likelihood(varargin{:}), 1)));
            logl = logl - K * D * self.params.ClusterCost;
        end
        
        
        function ids = cluster(self, varargin)
            % Return cluster ids for all spikes
            %   ids = cluster(self) returns the cluster ids for the entire
            %   dataset.
            %
            %   ids = cluster(self, index) returns the cluster ids for the
            %   given subset.
            
            [~, ids] = max(self.posterior(varargin{:}), [], 1);
        end
        
        
        function [pairwise, n] = overlap(self, varargin)
            % Return matrix of cluster overlaps.
            %   [pairwise, n] = overlap(self)
            %   pairwise(i, j) contains the number of spikes assigned to
            %   cluster i that were generated by cluster j. n(i) is the
            %   number of spikes assigned to cluster i.
            
            post = self.posterior(varargin{:});
            [~, assignment] = max(post);
            K = size(post, 1);
            pairwise = zeros(K);
            n = hist(assignment, 1 : K);
            for i = 1 : K
                for j = 1 : K
                    pairwise(i, j) = sum(post(i, assignment == j));
                end
            end
        end
        
        
        function fp = falsePositives(self, varargin)
            % Return false positive rate for all clusters.
            %   fp = falsePositives(self)
            
            [pairwise, n] = self.overlap(varargin{:});
            fp = (sum(pairwise, 2) - diag(pairwise))' ./ n;
        end
        
        
        function fn = falseNegatives(self, varargin)
            % Return false negative rate for all clusters.
            %   fn = falseNegatives(self)
            
            [pairwise, n] = self.overlap(varargin{:});
            fn = 1 - diag(pairwise)' ./ n;
        end
        
        
        function plot(self, d)
            % Plot mixture model
            
            Ytrain = self.Y(:, self.train);
            ttrain = self.t(self.train);
            if nargin < 2
                if size(Ytrain, 1) > 3
                    d = [7 1 4 10]; % tetrodes
                else
                    d = 1 : 3;      % single electrodes
                end
            end
            d(d > size(Ytrain, 1)) = [];
            
            j = self.cluster(self.train);
            K = max(j);
            figure(2), clf, hold all
            c = lines;
            hdl = zeros(1, K);
            
            dd = combnk(d, 2);
            
            if size(Ytrain, 1) > 1
                N = size(dd, 1);
                for k = 1 : N
                    subplot(2, N, k)
                    cla
                    hold on
                    for i = 1:K
                        plot(Ytrain(dd(k, 1), j == i), Ytrain(dd(k, 2), j == i), '.', 'markersize', 1, 'color', c(i, :))
                        hdl(i) = plot(self.mu(dd(k, 1), :, i), self.mu(dd(k, 2), :, i), '-', 'color', c(i, :),'LineWidth',3);
                    end
                    xlim(quantile(Ytrain(dd(k, 1), :), [0.001 0.999]));
                    ylim(quantile(Ytrain(dd(k, 2), :), [0.001 0.999]));
                end
                subplot(2, N, N + (1 : N))
            else
                subplot(111);
            end
            
            
            cla
            hold on
            for i = 1:K
                plot(ttrain(j == i), Ytrain(d(1), j == i), '.', 'markersize', 1, 'color', c(i, :))
                mu_t = self.mu_t(1 : end - 1) + min(diff(self.mu_t)) / 2;
                hdl(i) = plot(mu_t, self.mu(d(1), :, i), '-', 'color', c(i, :), 'LineWidth', 3);
            end
            legend(hdl, arrayfun(@(x) sprintf('Cluster %d', x), 1 : K, 'UniformOutput', false))
            ylim(quantile(Ytrain(d(1),:), [0.001 0.999]));
        end
        
        
        function self = compress(self)
            % Compress MoKsm object
            
            self.Y = [];
            self.t = [];
            self.like = [];
            self.post = [];
            self.blockId = [];
            self.spikeId = {};
        end
        
        
        function self = uncompress(self, Y, t)
            % Uncompress MoKsm object.
            
            self = self.initialize(Y, t, self.train, self.test, ...
                self.mu_t, self.mu, self.C, self.Cmu, self.priors, self.df);
        end
        
    end
    
    
    methods (Access = protected)
        
        function [mu, C, Cmu, priors, df] = expand(self)
            mu = self.mu;
            C = self.C;
            Cmu = self.Cmu;
            priors = self.priors;
            df = self.df;
        end
        
        
        function self = collect(self, mu, C, Cmu, priors, df)
            self.mu = mu;
            self.C = C;
            self.Cmu = Cmu;
            self.priors = priors;
            self.df = df;
        end
        
    end
    
    
    methods (Access = private)
        
        function self = EM(self, maxIter)
            % Run EM.
            %   self = EM(self) runs the EM iteration until convergence.
            %
            %   self = EM(self, maxIter) runs the EM iteration until
            %   convergence but at most maxIter iterations.
            
            if nargin < 2, maxIter = Inf; end
            
            [mu, C, Cmu, priors, df] = self.expand(); %#ok<*PROP>
            Ytrain = self.Y(:, self.train);
            N = size(Ytrain, 2);
            [D, T, K] = size(mu);
            Cf = zeros(D, D, T);
                
            % Initial E step
            like = zeros(K, N);
            for k = 1 : K
                muk = mu(:, self.blockId(self.train), k);
                like(k, :) = priors(k) * MoKsm.mixtureDistribution(Ytrain - muk, C(:, :, k), df);
            end
            p = sum(like, 1);
            post = bsxfun(@rdivide, like, p);
            post(:, p == 0) = 0;
            
            % Perform EM iterations until convergence or maxIter
            iter = 0;
            logLikeBase = NaN;
            while iter < maxIter && (iter < 2 || (self.logLike(end) - self.logLike(end - 1)) / (self.logLike(end - 1) - logLikeBase) > self.params.Tolerance)
                
                if ~mod(iter, 5), fprintf('.'), end
                iter = iter + 1;

                % Perform M step
                for k = 1 : K
                    
                    postk = post(k, :);
                    muk = mu(:, :, k);
                    Ck = C(:, :, k);
                    iCk = inv(Ck);
                    
                    % Additional latent variable for mixture of t-distributions
                    if ~isinf(df)
                        % u = (df + D) / (df + (Y-mu)'*Ck^-1*(Y-mu))
                        Ymu = Ytrain - muk(:, self.blockId(self.train));
                        [R, ~] = chol(Ck);
                        mahal_sq_dist = sum((R' \ Ymu).^2, 1);
                        uk = (df + D) ./ (df + mahal_sq_dist);
                    else
                        uk = ones(1, size(Ytrain,2));
                    end
                    
                    % Initialize Kalman update 
                    Cf_0 = Ck; % Initial uncertainty on mean
                    idx = self.spikeId{1};
                    if isempty(idx)
                        Cf(:, :, 1) = Cf_0;
                    else
                        pred_infomat = inv(Cf_0);
                        meas_infomat = sum(postk(idx) .* uk(idx)) * iCk; %#ok<MINV>
                        meas_infovec = (iCk * Ytrain(:,idx)) * (postk(idx) .* uk(idx))'; %#ok<MINV>
                        Cf(:, :, 1) = inv(pred_infomat + meas_infomat);
                        muk(:, 1) = Cf(:,:,1) * (pred_infomat * muk(:,1) + meas_infovec); %#ok<MINV>
                    end
                    
                    % Forward iteration for updating the means (Eq. 9)
                    iCfCmu = zeros(D, D, T);
                    for tt = 2 : T
                        idx = self.spikeId{tt};
                        iCfCmu(:, :, tt - 1) = inv(Cf(:, :, tt - 1) + Cmu);
                        if isempty(idx)
                            Cf(:, :, tt) = Cf(:, :, tt - 1) + Cmu;
                            muk(:, tt) = Cf(:, :, tt) * (iCfCmu(:, :, tt - 1) * muk(:, tt - 1));
                        else
                            piCk = sum(postk(idx) .* uk(idx)) * iCk; %#ok
                            Cf(:, :, tt) = inv(iCfCmu(:, :, tt - 1) + piCk);
                            muk(:, tt) = Cf(:, :, tt) * (iCfCmu(:, :, tt - 1) * muk(:, tt - 1) + ...
                                (iCk * Ytrain(:, idx)) * (postk(idx) .* uk(idx))'); %#ok
                        end
                    end
                    
                    % Backward iteration for updating the means (Eq. 10)
                    for tt = T-1 : -1 : 1
                        muk(:, tt) = muk(:, tt) + Cf(:, :, tt) * (iCfCmu(:, :, tt) * (muk(:, tt + 1) - muk(:, tt)));
                    end
                    assert(~any(isnan(muk(:))), 'Got nan');
                    
                    % Update cluster covariances (Eq. 11)
                    Ymu = Ytrain - muk(:, self.blockId(self.train));
                    Ck = (bsxfun(@times, postk .* uk, Ymu) * Ymu') / sum(postk);
                    Ck = Ck + eye(D) * self.params.CovRidge; % add ridge to regularize
                    
                    mu(:, :, k) = muk;
                    C(:, :, k) = Ck;
                end
                
                % update class priors
                priors = sum(post, 2) / N;
                
                % Perform E step
                like = zeros(K, N);
                for k = 1 : K
                    muk = mu(:, self.blockId(self.train), k);
                    like(k, :) = priors(k) * MoKsm.mixtureDistribution(Ytrain - muk, C(:, :, k), df);
                end
                p = sum(like, 1);
                post = bsxfun(@rdivide, like, p);
                post(:, p == 0) = 0;
                
                % check for starvation
                [~, assignments] = max(post, [], 1);
                if any(priors * N < 2 * D) || any(~ismember(1 : K, assignments))
                    error('MoKsm:starvation', 'Component starvation: cluster %d', ...
                        find((priors' * N < 2 * D) | ~ismember(1 : K, assignments), 1))
                end
                
                % calculate log-likelihood
                ll_spikes = MoKsm.mylog(sum(like, 1)); % [1 x N]
                if T > 1
                    drift = diff(mu,1,2);
                    ll_drift = MoKsm.mylog(MoKsm.mvn(drift(:,:), Cmu)); % [1 x (T-1)*K]
                    ll_drift = sum(reshape(ll_drift, T-1, K), 1); % [1 x K]
                else
                    ll_drift = zeros(1, K);
                end
                self.logLike(end + 1) = (sum(ll_spikes) + sum(ll_drift))/N;
                if self.params.Verbose && numel(self.logLike) > 1
                    figure(1)
                    plot(self.logLike, '.-k')
                    ylim(prctile(self.logLike, [10 100]))
                    drawnow
                end
            
                if iter == 1
                    logLikeBase = self.logLike(end);
                end
            end
            
            self = self.collect(mu, C, Cmu, priors, df);
       end
        
        
        function [self, success] = tryMerge(self)
            % Try merging clusters
            %   Merge is accepted if penalized average likelihood improved.
            
            verbose = self.params.Verbose;
            success = false;
            cands = self.getMergeCandidates();
            logLikeTest = self.logLikelihood(self.test);
            for ij = cands'
                try
                    fprintf('Trying to merge clusters %d and %d ', ij(1), ij(2))
                    newSelf = self.mergeClusters(ij);
                    newSelf = newSelf.EM();
                    newLogLikeTest = newSelf.logLikelihood(self.test);
                    if newLogLikeTest > logLikeTest
                        fprintf(' success (likelihood improved by %.5g)\n', newLogLikeTest - logLikeTest)
                        self = newSelf;
                        self.logLike(end + 1) = NaN;
                        success = true;
                        if verbose, plot(self), end
                        break
                    else
                        fprintf(' aborted\n')
                    end
                catch err
                    if strcmp(err.identifier, 'MoKsm:starvation')
                        fprintf(' aborted due to component starvation\n')
                    else
                        rethrow(err)
                    end
                end
            end
        end
        
        
        function [self, success] = trySplit(self)
            % Try splitting clusters.
            %   Split is accepted if penalized average likelihood improved.
            
            verbose = self.params.Verbose;
            success = false;
            splitCands = self.getSplitCandidates();
            logLikeTest = self.logLikelihood(self.test);
            for i = splitCands'
                try
                    fprintf('Trying to split cluster %d ', i)
                    newSelf = self.splitCluster(i);
                    newSelf = newSelf.EM();
                    newLogLikeTest = newSelf.logLikelihood(self.test);
                    if newLogLikeTest > logLikeTest
                        fprintf(' success (likelihood improved by %.5g)\n', newLogLikeTest - logLikeTest)
                        self = newSelf;
                        self.logLike(end + 1) = NaN;
                        success = true;
                        if verbose, plot(self), end
                        break
                    else
                        fprintf(' aborted\n')
                    end
                 catch err
                     if strcmp(err.identifier, 'MoKsm:starvation')
                         fprintf(' aborted due to component starvation\n')
                     else
                         rethrow(err)
                     end
                 end
            end
        end
        
        
        function partial = getPartial(self, ids)
            % Return partial model.
            %   The partial model contains only the given clusters and the
            %   spikes assigned to them.
            
            [~, assignments] = max(self.posterior(self.train), [], 1);
            ndx = ismember(assignments, ids);
            Yp = self.Y(:, self.train(ndx));
            tp = self.t(self.train(ndx));
            partial = MoKsm(self.params, 'Verbose', false);
            partial = partial.initialize(Yp, tp, 1 : numel(tp), [], self.mu_t, ...
                self.mu(:, :, ids), self.C(:, :, ids), self.Cmu, self.priors(ids), self.df);
        end%
        
      
        function cand = getSplitCandidates(self)
            p = self.likelihood(self.train);
            pk = bsxfun(@rdivide, p, self.priors);
            post = self.posterior(self.train);
            fk = bsxfun(@rdivide, post, sum(post, 2));
            Jsplit = sum(fk .* (MoKsm.mylog(fk) - MoKsm.mylog(pk)), 2);
            [~, cand] = sort(Jsplit, 'descend');
            [D, T] = size(post);
            [~,assignments] = max(post,[],1);
            cand = cand(ismember(cand, assignments) & self.priors(cand) * T > 4 * D);  % don't split small clusters
        end
        
        
        function cand = getMergeCandidates(self)
            post = self.posterior(self.train);
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
        end
        
        
        function self = updateCache(self)
            % Update cached likelihoods and posteriors
            
            self.like = self.likelihood(':');
            self.post = self.posterior(':');
        end
        
    end

    
    methods(Static)
        
        function y = mylog(x)
            % Natural logarithm excluding zeros
            
            y = reallog(x);
            y(x == 0) = 0;
        end
        
        
        function p = mixtureDistribution(X, C, df)
            % Probability distribution of mixture components (normal or t)
            
            if isinf(df)
                p = MoKsm.mvn(X, C);
            else
                p = MoKsm.mvt(X, C, df);
            end
        end
        
        
        function p = mvn(X, C)
            % Zero-mean multivariate normal probability density
            %   p = mvn(X, C) calculates the density of the multivariate normal
            %   distribution with zero mean and covariance matrix C at X. X is assumed
            %   to be a column vector or a matrix of multiple column vectors, in which
            %   case the result, p, is a row vector.
            %
            % AE 2012-02-05
            
            D = size(C, 1);
            const = (2*pi)^(-D/2);
            [Ch, ~] = chol(C);
            p = const / prod(diag(Ch)) * exp(-1/2 * sum((Ch' \ X).^2, 1));
        end
        
        
        function p = mvt(X, C, df)
            % Zero-mean multivariate Student's t probability density
            %   p = mvt(X, C, df) calculates the density of the multivariate t
            %   distribution with scale parameter C and df degrees of freedom
            %   at X. X is assumed to be a column vector or a matrix of multiple
            %   column vectors, in which case the result, p, is a row vector.
            D = size(C, 1);
            [Ch, ~] = chol(C);
            delta = sum((Ch' \ X).^2, 1);
            p = exp(gammaln((df + D) / 2) - gammaln(df / 2) ...
                - ((df + D) / 2) .* log(1 + delta / df) ...
                - sum(log(diag(Ch))) - (D / 2) * log(df * pi));
        end
        
    end
    
end
