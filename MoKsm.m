classdef MoKsm
    % Mixture of Kalman filters model adapted from Calabrese & Paninski
    % 2011.
    % 
    % To improve efficiency the algorithm is slightly modified so that it
    % doesn't update the cluster means at each spike but rather uses longer
    % blocks of time, whithin which the means are assumed to be constant.
    %
    % We use a split & merge strategy to optimize a penalized average
    % log-likelihood.
    %
    % Alexander S. Ecker & R. James Cotton
    % 2012-07-04
    
    properties
        params      % parameters for fitting
        model       % mixture model
        post        % posteriors for class membership
        logLike     % log-likelihood curve during fitting
        Y           % data
        t           % times
        train       % indices of training set: Y(:, train), t(train)
        test        % indices of test set
        blockId     % mapping from data point to time blocks
        spikeId     % spike ids of the training data for each time block
    end
    
    methods

        function self = MoKsm(varargin)
            
            % parse optional parameters
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('MaxTrainSpikes', 20000); % max. number of spikes for training data
            p.addOptional('MaxTestSpikes', 50000);  % max. number of spikes for test data
            p.addOptional('TrainFrac', 0.8);        % fraction of data points used for training
            p.addOptional('Tolerance', 0.0002);     % tolerance for determining convergence
            p.addOptional('Verbose', false);        % verbose output
            p.addOptional('Seed', 1);               % seed for random number generator
            p.addOptional('Df', 2);                 % degrees of freedom for t distribution
            p.addOptional('CovRidge', 1.5);         % independent variance added to cluster covariances
            p.addOptional('DriftRate', 10 / 3600 / 1000);  % drift rate per ms
            p.addOptional('DTmu', 60 * 1000);              % block size for means in ms
            p.addOptional('ClusterCost', 0.0025);   % penalizer for adding additional clusters
            p.addOptional('IterPerBlock', 10);      % log-likelihood is evaluated every n iterations
            p.parse(varargin{:});
            self.params = p.Results;
        end
        
        
        function self = fit(self, Y, t)
            % Fit the model
            
            % make sure dimensions of input are correct
            if size(Y, 1) == length(t)
                Y = Y';
            elseif size(Y, 2) ~= length(t)
                error('Time dimension doesn''t match dataset');
            end
            assert(size(Y, 1) <= 50, 'Dimensionality way too high');
            t = reshape(t, 1, []);
            
            % sort by time
            [self.t, order] = sort(t);
            self.Y = Y(:, order);

            % ensure deterministic behavior
            rng(self.params.Seed);
            
            % split into training & test data
            T = size(self.Y,2);
            rnd = randperm(T);
            nTrain = fix(self.params.TrainFrac * T);
            self.train = sort(rnd(1 : min(nTrain,self.params.MaxTrainSpikes)));
            self.test = sort(rnd(nTrain + 1 : min(end, nTrain + self.params.MaxTestSpikes)));
            
            % Assign spikes to time blocks
            self.model.mu_t = self.t(1) : self.params.DTmu : self.t(end) + self.params.DTmu;
            nTime = length(self.model.mu_t) - 1;
            [~, blockIdTrain] = histc(self.t(self.train), self.model.mu_t);

            self.spikeId = arrayfun(@(x) find(blockIdTrain == x), 1 : nTime, 'UniformOutput', false);
            unsupported = find(cellfun(@isempty, self.spikeId)) + 1;
            self.model.mu_t(unsupported) = [];
            nTime = length(self.model.mu_t) - 1;
            self.spikeId(unsupported) = [];
            
            % block ids for all data points
            [~, self.blockId] = histc(self.t, self.model.mu_t);
            
            fprintf('Estimating cluster means at %d time points\n', nTime);

            % fit initial model with one component
            fprintf('Running initial Kalman filter model with one cluster ')
            Ytrain = self.Y(:, self.train);
            self.model.mu = repmat(mean(Ytrain, 2), [1 nTime]); % cluster means
            self.model.C = cov(Ytrain');                                  % observation covariance
            self.model.Cmu = eye(size(self.model.mu, 1)) * self.params.DTmu * self.params.DriftRate;
            self.model.df = self.params.Df;
            self.model.priors = 1;                       % cluster weights
            self = self.EM();
            fprintf(' done (likelihood: %.5g)\n', self.logLikelihood(self.test))
            
            % Run split & merge
            % We alternate between trying to split and trying to merge. Both are done
            % until no candidate leads to success. If both splitting and merging didn't
            % lead to success we terminate
            op = {@trySplit, @tryMerge};
            i = 1;
            success = true(1, 2);
            while any(success)
                [self, success(i)] = op{i}(self);
                if ~success(i)
                    i = 3 - i;
                end
            end
            
            fprintf('Done with split & merge\n')
            
            fprintf('--\n')
            fprintf('Number of clusters: %d\n', size(self.model.mu, 3))
            fprintf('Log-likelihoods\n')
            fprintf('  training set: %.8g\n', self.logLikelihood(self.train))
            fprintf('      test set: %.8g\n', self.logLikelihood(self.test))
            fprintf('\n\n')
        end

        
        function [self, success] = tryMerge(self)
            verbose = self.params.Verbose;
            
            success = false;
            cands = self.getMergeCandidates();
            logLikeTest = self.logLikelihood(self.test);
            for ij = cands'
                try
                    fprintf('Trying to merge clusters %d and %d ', ij(1), ij(2))
                    newSelf = self.merge(ij);
                    newSelf = newSelf.EM(2);
                    newLogLikeTest = newSelf.logLikelihood(self.test);
                    if newLogLikeTest > logLikeTest
                        fprintf(' success (likelihood improved by %.5g)\n', newLogLikeTest - logLikeTest)
                        self = newSelf;
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
            verbose = self.params.Verbose;
            
            success = false;
            splitCands = self.getSplitCandidates();
            logLikeTest = self.logLikelihood(self.test);
            for i = splitCands'
                try
                    fprintf('Trying to split cluster %d ', i)
                    newSelf = self.split(i);
                    newSelf = newSelf.EM(2);
                    newLogLikeTest = newSelf.logLikelihood(self.test);
                    if newLogLikeTest > logLikeTest
                        fprintf(' success (likelihood improved by %.5g)\n', newLogLikeTest - logLikeTest)
                        self = newSelf;
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
        
        
        function self = EM(self, maxIter)
            % Run EM until convergence.
            
            if nargin < 2, maxIter = Inf; end
            iter = 0;
            logLikeBase = NaN;
            while iter < maxIter && (iter < 2 || (self.logLike(end) - self.logLike(end - 1)) / (self.logLike(end - 1) - logLikeBase) > self.params.Tolerance)
                
                fprintf('.')
                iter = iter + 1;
                
                % run a couple of EM iterations
                for i = 1 : self.params.IterPerBlock
                    self = self.EStep();
                    self = self.MStep();
                end
                
                % calculate log-likelihood
                self.logLike(end + 1) = self.logLikelihood(self.train);
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
        end

        
        function self = EStep(self)
            % Perform E step
            
            self.post = self.posterior(self.train);
        end
        
        
        function self = MStep(self)
            % Perform M step
            
            Ytrain = self.trainingData();
            N = size(Ytrain, 2);
            
            % EM recursion
            [mu, C, Cmu, ~, df] = self.expand();
            [D, T, K] = size(mu);
                        
            for k = 1 : K
                
                postk = self.post(k, :);
                muk = mu(:, :, k);
                Ck = C(:, :, k);
                
                % Forward step for updating the means (Eq. 9)
                Cf = zeros(D, D, T);                % state covariances
                iCfCmu = zeros(D, D, T);
                Cf(:, :, 1) = Ck;
                iCk = inv(Ck);
                for tt = 2 : T
                    idx = self.spikeId{tt-1};
                    piCk = sum(postk(idx)) * iCk; %#ok
                    
                    % hacky, for now just hopping along the time axis
                    iCfCmu(:, :, tt - 1) = inv(Cf(:, :, tt - 1) + Cmu);
                    Cf(:, :, tt) = inv(iCfCmu(:, :, tt - 1) + piCk);
                    muk(:, tt) = Cf(:, :, tt) * (iCfCmu(:, :, tt - 1) * muk(:, tt - 1) + ...
                        (iCk * Ytrain(:, idx)) * postk(idx)'); %#ok
                end
                
                % Backward step for updating the means (Eq. 10)
                for tt = T-1 : -1 : 1
                    muk(:, tt) = muk(:, tt) + Cf(:, :, tt) * (iCfCmu(:, :, tt) * (muk(:, tt + 1) - muk(:, tt)));
                end
                assert(~any(isnan(muk(:))), 'Got nan');
                
                % Update observation covariance (Eq. 11)
                Ymu = Ytrain - muk(:, self.blockId(self.train));
                Ck = (bsxfun(@times, postk, Ymu) * Ymu') / sum(postk);
                Ck = Ck + eye(D) * self.params.CovRidge; % add ridge to regularize
                                
                mu(:, :, k) = muk;
                C(:, :, k) = Ck;
            end
            
            % update class priors
            priors = sum(self.post, 2) / N;
                
            % check for starvation
            if any(priors * N < 2 * D)
                error('MoKsm:starvation', 'Component starvation: cluster %d', find(priors * N < 2 * D, 1))
            end
            
            self = self.collect(mu, C, Cmu, priors, df);
        end
        
        
        function self = split(self, k)
            % Split cluster k
            
            [mu, C, Cmu, priors, df] = self.expand();
            [D, ~, K] = size(mu);
            deltaMu  = chol(C(:, :, k))' * randn(D, 1);
            mu(:, :, k) = bsxfun(@plus, mu(:, :, k), deltaMu);
            mu(:, :, K + 1) = bsxfun(@minus, mu(:, :, k), deltaMu);
            C(:, :, k) = det(C(:, :, k))^(1 / D) * eye(D);
            C(:, :, K + 1) = C(:, :, k);
            priors(k) = priors(k) / 2;
            priors(K + 1) = priors(k);
            self = self.collect(mu, C, Cmu, priors, df);
        end
        
        
        function self = merge(self, ids)
            % Merge clusters by comuting prior-weighted parameter averages.
            
            [mu, C, Cmu, priors, df] = self.expand();
            p = permute(priors(ids), [3 2 1]);
            mu(:, :, ids(1)) = sum(bsxfun(@times, mu(:, :, ids), p), 3) / sum(p);
            C(:, :, ids(1)) = sum(bsxfun(@times, C(:, :, ids), p), 3) / sum(p);
            priors(ids(1)) = sum(p);
            mu(:, :, ids(2 : end)) = [];
            C(:, :, ids(2 : end)) = [];
            priors(ids(2 : end)) = [];
            self = self.collect(mu, C, Cmu, priors, df);
        end
        
        
        function [Ytrain, ttrain] = trainingData(self)
            % Return training data
            
            Ytrain = self.Y(:, self.train);
            ttrain = self.t(self.train);
        end
        
        
        function [Ytest, ttest] = testData(self)
            % Return test data
            
            Ytest = self.Y(:, self.test);
            ttest = self.t(self.test);
        end

        
        function p = likelihood(self, index)
            % likelihood of the data for each cluster
            
            if nargin < 2, index = ':'; end
            [mu, C, Cmu, priors, df] = self.expand();
            K = numel(priors);
            p = zeros(K, size(self.Y(:, index), 2));
            for k = 1 : K
                muk = mu(:, self.blockId(index), k);
                p(k, :) = priors(k) * MoKsm.mixtureDistribution(self.Y(:, index) - muk, C(:, :, k) + Cmu, df);
            end
        end
        
        
        function post = posterior(self, varargin)
            % posterior of class membership for each spike
            
            likelihood = self.likelihood(varargin{:});
            p = sum(likelihood, 1);
            post = bsxfun(@rdivide, likelihood, p);
            post(:, p == 0) = 0;            
        end
        
        
        function logl = logLikelihood(self, varargin)
            % Penalized average log-likelihood
            
            [D, ~, K] = size(self.model.mu);
            logl = mean(MoKsm.mylog(sum(self.likelihood(varargin{:}), 1)));
            logl = logl - K * D * self.params.ClusterCost;
        end
        
        
        function ids = cluster(self, varargin)
            % return cluster ids for all spikes
            
            [~, ids] = max(self.posterior(varargin{:}), [], 1);
        end
        
        
        function [pairwise, n] = overlap(self)
            % Return matrix of cluster overlaps.
            %   [pairwise, n] = overlap(self)
            %   pairwise(i, j) contains the number of spikes assigned to
            %   cluster i that were generated by cluster j. n(i) is the
            %   number of spikes assigned to cluster i.
            
            [~, assignment] = max(self.post);
            K = size(self.post, 1);
            pairwise = zeros(K);
            n = hist(assignment, 1 : K);
            for i = 1 : K
                for j = 1 : K
                    pairwise(i, j) = sum(self.post(i, assignment == j));
                end
            end
        end
        
        
        function fp = falsePositives(self)
            % Return false positive rate for all clusters.
            %   fp = falsePositives(self)
            
            [pairwise, n] = overlap(self);
            fp = (sum(pairwise, 2) - diag(pairwise))' ./ n;
        end
        
        
        function fn = falseNegatives(self)
            % Return false negative rate for all clusters.
            %   fn = falseNegatives(self)
            
            [pairwise, n] = overlap(self);
            fn = 1 - diag(pairwise)' ./ n;
        end

        
        function plot(self, d)
            % Plot mixture model
            
            [Ytrain, ttrain] = self.trainingData();

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
                        hdl(i) = plot(self.model.mu(dd(k, 1), :, i), self.model.mu(dd(k, 2), :, i), '-', 'color', c(i, :),'LineWidth',3);
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
                mu_t = self.model.mu_t(1 : end - 1) + min(diff(self.model.mu_t)) / 2;
                hdl(i) = plot(mu_t, self. model.mu(d(1), :, i), '-', 'color', c(i, :), 'LineWidth', 3);
            end
            legend(hdl, arrayfun(@(x) sprintf('Cluster %d', x), 1 : K, 'UniformOutput', false))
            ylim(quantile(Ytrain(d(1),:), [0.001 0.999]));
        end
        
    end
    
    
    methods (Access = private)
        
        function cand = getSplitCandidates(self)
            pk = bsxfun(@rdivide, self.likelihood(self.train), self.model.priors);
            fk = bsxfun(@rdivide, self.post, sum(self.post, 2));
            Jsplit = sum(fk .* (MoKsm.mylog(fk) - MoKsm.mylog(pk)), 2);
            [~, cand] = sort(Jsplit, 'descend');
            [D, T] = size(self.post);
            cand = cand(self.model.priors(cand) * T > 4 * D); % don't split small clusters
        end
        
        
        function cand = getMergeCandidates(self)
            K = size(self.post, 1);
            maxCandidates = ceil(K * sqrt(K) / 2);
            np = sqrt(sum(self.post .* self.post, 2));
            Jmerge = zeros(K * (K - 1) / 2, 1);
            cand = zeros(K * (K - 1) / 2, 2);
            k = 0;
            for i = 1:K
                for j = i+1:K
                    k = k + 1;
                    Jmerge(k) = self.post(i, :) * self.post(j, :)' / (np(i) * np(j));
                    cand(k, :) = [i j];
                end
            end
            [~, order] = sort(Jmerge, 'descend');
            cand = cand(order(1:min(end, maxCandidates)), :);
        end
        
        
        function [mu, C, Cmu, priors, df] = expand(self)
            mu = self.model.mu;
            C = self.model.C;
            Cmu = self.model.Cmu;
            priors = self.model.priors;
            df = self.model.df;
        end
        
        
        function self = collect(self, mu, C, Cmu, priors, df)
            self.model.mu = mu;
            self.model.C = C;
            self.model.Cmu = Cmu;
            self.model.priors = priors;
            self.model.df = df;
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
        end
        
        
        function p = mvt(X, C, df)
            % Multivariate Student's t probability density.
            %   p = mvt(x, mu, Sigma) calculates the density of the multivariate t
            %   distribution mean mu, covariance matrix Sigma, and df degrees of
            %   freedom at x. x is assumed to be a row vector or a matrix of multiple
            %   row vectors, in which case the result, p, is a column vector.
            D = size(C, 1);
            [Ch, ~] = chol(C);
            delta = sum((Ch' \ X).^2, 1);
            p = exp(gammaln((df + D) / 2) - gammaln(df / 2) ...
                - ((df + D) / 2) .* log(1 + delta / df) ...
                - sum(log(diag(Ch))) - (D / 2) * log(df * pi));
        end
        
    end
end
