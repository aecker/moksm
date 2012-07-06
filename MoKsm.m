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
        like        % likelihoods (training data)
        post        % posteriors for class membership (training data)
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
                    newSelf = self.mergeClusters(ij);
                    newSelf = newSelf.EM(20);
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
                    newSelf = self.splitCluster(i);
                    newSelf = newSelf.EM(20);
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
        
        
            
            % update class priors
            priors = sum(self.post, 2) / N;
                
            % check for starvation
            if any(priors * N < 2 * D)
                error('MoKsm:starvation', 'Component starvation: cluster %d', find(priors * N < 2 * D, 1))
            end
            
            self = self.collect(mu, C, Cmu, priors, df);
        end
        
        
        function self = splitCluster(self, k)
            % Split cluster k
            
            [mu, C, Cmu, priors, df] = self.expand();
            [D, ~, K] = size(mu);
            deltaMu  = chol(C(:, :, k))' * randn(D, 1) * 0.2;
            mu(:, :, k) = bsxfun(@plus, mu(:, :, k), deltaMu);
            mu(:, :, K + 1) = bsxfun(@minus, mu(:, :, k), deltaMu);
            C(:, :, k) = det(C(:, :, k))^(1 / D) * eye(D);
            C(:, :, K + 1) = C(:, :, k);
            priors(k) = priors(k) / 2;
            priors(K + 1) = priors(k);
            self = self.collect(mu, C, Cmu, priors, df);
        end
        
        
        function self = mergeClusters(self, ids)
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
        
    end
        
end
