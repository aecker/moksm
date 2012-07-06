classdef MoKsm < MoK
    % Mixture of Kalman filters model adapted from Calabrese & Paninski
    % 2011.
    %
    % This class fits a Mixture of Kalman filters model using a split &
    % merge approach. We optimize a penalized average likelihood, where the
    % penalizer is a constant cost per cluster.
    %
    % Alexander S. Ecker & R. James Cotton
    % 2012-07-04
    
    properties
        logLike     % log-likelihood curve during fitting
        Y           % data
        t           % times
        train       % indices of training set: Y(:, train), t(train)
        test        % indices of test set
        blockId     % mapping from data point to time blocks
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
            
            self = self@MoK(varargin{:});
            
            % parse optional parameters
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('MaxTrainSpikes', 20000);
            p.addOptional('MaxTestSpikes', 50000);
            p.addOptional('TrainFrac', 0.8);
            p.addOptional('Seed', 1);
            p.addOptional('Df', 2);
            p.addOptional('DriftRate', 10 / 3600 / 1000);
            p.addOptional('DTmu', 60 * 1000);

            % MoK params (overwrites defaults there)
            p.addOptional('Tolerance', 0.0002);
            p.addOptional('Verbose', false);
            p.addOptional('CovRidge', 1.5);
            p.addOptional('ClusterCost', 0.0025);
            
            p.parse(varargin{:}, self.params);
            self.params = p.Results;
        end
        
        
        function self = fit(self, Y, t)
            % Fit the model
            %   self = fit(self, Y, t) fits the model to data (Y, t). Y is
            %   a matrix of size (#dimensions x #samples) and t is a vector
            %   of length #samples.
            %
            %   See MoKsm for optional parameters to use for fitting.
            
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
            Ytrain = self.Y(:, self.train);
            ttrain = self.t(self.train);
            
            % Assign spikes to time blocks
            self.mu_t = self.t(1) : self.params.DTmu : self.t(end) + self.params.DTmu;
            nTime = length(self.mu_t) - 1;
            [~, self.blockId] = histc(self.t, self.mu_t);
            fprintf('Estimating cluster means at %d time points\n', nTime);

            % fit initial model with one component
            fprintf('Running initial Kalman filter model with one cluster ')
            mu = repmat(mean(Ytrain, 2), [1 nTime]);
            C = cov(Ytrain');
            Cmu = eye(size(mu, 1)) * self.params.DTmu * self.params.DriftRate;
            df = self.params.Df;
            self = self.initialize(Ytrain, ttrain, self.mu_t, mu, C, Cmu, 1, df);
            self = self.EM(10);
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
                if ~success(i)
                    i = 3 - i;
                end
            end
            
            fprintf('Final model fit ')
            self = self.EM();
            fprintf('\nDone with split & merge\n')
            
            fprintf('--\n')
            fprintf('Number of clusters: %d\n', size(self.mu, 3))
            fprintf('Log-likelihoods\n')
            fprintf('  training set: %.8g\n', self.logLikelihood())
            fprintf('      test set: %.8g\n', self.logLikelihood(self.Y(:, self.test), self.blockId(self.test)))
            fprintf('\n\n')
        end

        
        function [self, success] = tryMerge(self)
            % Try merging clusters
            %   Merge is accepted if penalized average likelihood improved.
            
            verbose = self.params.Verbose;
            success = false;
            cands = self.getMergeCandidates();
            Ytest = self.Y(:, self.test);
            blockTest = self.blockId(self.test);
            logLikeTest = self.logLikelihood(Ytest, blockTest);
            for ij = cands'
                try
                    fprintf('Trying to merge clusters %d and %d ', ij(1), ij(2))
                    newSelf = self.mergeClusters(ij);
                    newSelf = newSelf.EM(20);
                    newLogLikeTest = newSelf.logLikelihood(Ytest, blockTest);
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
                    if strcmp(err.identifier, 'MoK:starvation')
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
            Ytest = self.Y(:, self.test);
            blockTest = self.blockId(self.test);
            logLikeTest = self.logLikelihood(Ytest, blockTest);
            for i = splitCands'
                try
                    fprintf('Trying to split cluster %d ', i)
                    newSelf = self.splitCluster(i);
                    newSelf = newSelf.EM(20);
                    newLogLikeTest = newSelf.logLikelihood(Ytest, blockTest);
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
                     if strcmp(err.identifier, 'MoK:starvation')
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
            
            [~, assignments] = max(self.posterior(), [], 1);
            ndx = ismember(assignments, ids);
            Yp = self.Ytrain(:, ndx);
            tp = self.ttrain(ndx);
            partial = MoK(self.params);
            partial = partial.initialize(Yp, tp, self.mu_t, self.mu(:, :, ids), ...
                self.C(:, :, ids), self.Cmu, self.priors(ids), self.df);
        end%
        
      
        function self = splitCluster(self, k)
            % Split cluster k.
            %   Randomly perturb cluster means followed by partial EM
            %   update.

            partial = self.getPartial(k);
            partial = partial.splitCluster(1);
            partial = partial.EM(20);
            self.mu(:, :, [k end+1]) = partial.mu;
            self.C(:, :, [k end+1]) = partial.C;
            self.priors([k end+1]) = partial.priors;
        end
        
        
        function self = mergeClusters(self, ids)
            % Merge clusters.
            %   Merge by comuting prior-weighted parameter averages
            %   followed by partial EM update.
            
            partial = self.getPartial(ids);
            partial = partial.mergeClusters(1 : numel(ids));
            partial = partial.EM(20);
            self.mu(:, :, ids(1)) = partial.mu;
            self.C(:, :, ids(1)) = partial.C;
            self.priors(ids(1)) = partial.priors;
            self.mu(:, :, ids(2 : end)) = [];
            self.C(:, :, ids(2 : end)) = [];
            self.priors(ids(2 : end)) = [];
        end
        
    end
    
    
    methods (Access = private)
        
        function cand = getSplitCandidates(self)
            p = self.likelihood();
            pk = bsxfun(@rdivide, p, self.priors);
            post = self.posterior();
            fk = bsxfun(@rdivide, post, sum(post, 2));
            Jsplit = sum(fk .* (MoKsm.mylog(fk) - MoKsm.mylog(pk)), 2);
            [~, cand] = sort(Jsplit, 'descend');
            [D, T] = size(post);
            cand = cand(self.priors(cand) * T > 4 * D); % don't split small clusters
        end
        
        
        function cand = getMergeCandidates(self)
            post = self.posterior();
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
        
    end
        
end
