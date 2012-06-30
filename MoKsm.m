classdef MoKsm
    % MoK to 12d data (adapted from Calabrese & Paninski 2011)
    % AE 2012-02-03
    % JC 2012-02-15
    
    properties
        params;
        model;
        Y;
        t;
        Ytrain;
        Ytest
        ttrain;
        ttest;
        blockId = []; % maps the training data to mu's
    end
    
    methods

        function self = MoKsm(Y, t, varargin)
            
            % parse optional parameters
            p = inputParser;
            p.addOptional('MaxTrainSpikes', 20000); % max. number of spikes for training data
            p.addOptional('MaxTestSpikes', 50000);  % max. number of spikes for test data
            p.addOptional('TrainFrac', 0.8);        % max. number of spikes for test data
            p.addOptional('Tolerance', 0.0002);     % tolerance for determining convergence
            p.addOptional('Verbose', false);        % verbose output
            p.addOptional('Seed', 1);               % seed for random number generator
            p.addOptional('Df', 2);                 % degrees of freedom for t distribution
            p.addOptional('CovRidge', 1.5);         % independent variance in muV
            p.addOptional('DriftRate', 10 / 3600 / 1000);  % drift rate in muV/h
            p.addOptional('DTmu', 60 * 1000);              % block size for means (sec)
            p.addOptional('ClusterCost', 0.03);     % penalizer for adding additional clusters
            p.parse(varargin{:});
            self.params = p.Results;
            
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
        end
        
        
        % Fit the model
        function self = fitModel(self)
            
            rng(self.params.Seed);
            
            % split into training & test data
            T = size(self.Y,2);
            rnd = randperm(T);
            nTrain = fix(self.params.TrainFrac * T);
            train = sort(rnd(1 : min(nTrain,self.params.MaxTrainSpikes)));
            test = sort(rnd(nTrain + 1 : min(end, nTrain + self.params.MaxTestSpikes)));
            self.Ytrain = self.Y(:, train);
            self.Ytest = self.Y(:, test);
            self.ttrain = self.t(train);
            self.ttest = self.t(test);
            
            % Initialize model using 1 component
            self.model.mu_t = self.t(1):self.params.DTmu:self.t(end);
            nTime = length(self.model.mu_t);

            % Assign spikes to blocks for the M step
            blockId = zeros(size(self.ttrain));
            % Compute boundaries between clusters
            mp = [0 (self.model.mu_t(1:end-1) + self.model.mu_t(2:end)) / 2];
            % So spikes between mp(i) and mp(i+1) should be assigned to
            % mu(i)
            lastIdx = 1;
            for i = 1:length(mp)-1
                startIdx = lastIdx - 1 + ...
                    find(self.ttrain(lastIdx:end) >= mp(i) & ...
                    self.ttrain(lastIdx:end) < mp(i+1), 1, 'first');
                endIdx = lastIdx - 1 + find(self.ttrain(lastIdx:end) >= mp(i+1), 1, 'first');
                if isempty(startIdx) % Found no spikes in this block
                    continue;
                end
                if isempty(endIdx)
                    blockId(startIdx:end) = i;
                    lastIdx = length(blockId);
                    break;
                end
                blockId(startIdx:endIdx-1) = i;
                lastIdx = endIdx;
            end
            blockId(lastIdx+1:end) = i+1;

            self.blockId = arrayfun(@(x) find(blockId == x), 1:nTime, 'UniformOutput', false);
            unsupported_blocks = find(cellfun(@isempty, self.blockId));
            self.model.mu_t(unsupported_blocks) = [];
            nTime = length(self.model.mu_t);
            self.blockId(unsupported_blocks) = [];
            
            fprintf('Estimating cluster means at %d time points\n', nTime);
            fprintf('Running initial Kalman filter model with one cluster ')

            self.model.mu = repmat(mean(self.Ytrain, 2), [1 nTime]); % cluster means
            self.model.C = cov(self.Ytrain');                                  % observation covariance
            self.model.Cmu = eye(size(self.model.mu, 1)) * self.params.DTmu * self.params.DriftRate;
            self.model.df = self.params.Df;
            self.model.priors = 1;                       % cluster weights
            self.model.post = ones(1, size(self.Ytrain,2));
            self.model.pk = MoKsm.mixtureDistribution(bsxfun(@minus,self.Ytrain,mean(self.Ytrain,2)), self.model.C + self.model.Cmu, self.model.df);
            self.model.logLike = sum(MoKsm.mylog(self.model.pk));
            
            
            self = fullEM(self);
            if self.params.Verbose, plot(self), end
            fprintf(' done (likelihood: %.5g)\n', self.model.logLike(end))
            
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
            
            fprintf('Performing fit on entire dataset ...');
            [self.model.logLike p] = self.evalTestSet();
            self.model.postAll = p;
            fprintf(' Done\n');
            
            fprintf('--\n')
            fprintf('Number of clusters: %d\n', size(self.model.mu, 3))
            fprintf('Log-likelihoods\n')
            fprintf('  training set: %.8g\n', self.model.logLike(end))
            fprintf('      test set: %.8g\n', self.evalTestSet())
            fprintf('\n\n')
        end

        function self = fullEM(self)
            % warning off MATLAB:nearlySingularMatrix
            
            Y = self.Ytrain;
            blockId = self.blockId;
            
            % EM recursion
            [mu, C, Cmu, priors, post, pk, logLike, mu_t, df] = MoKsm.expand(self.model);
            [D, T, K] = size(mu);
            
            iter = 0;
            logLikeBase = logLike(end);
            tic
            while iter < 2 || (logLike(end) - logLike(end - 1)) / (logLike(end - 1) - logLikeBase) > self.params.Tolerance
                
                if ~mod(iter, 10)
                    fprintf('.')
                end
                iter = iter + 1;
                
                parfor k = 1 : K
                    
                    this_post = post(k, :);
                    muk = mu(:, :, k);
                    Ck = C(:, :, k);
                    
                    % Forward step for updating the means (Eq. 9)
                    Cf = zeros(D, D, T);                % state covariances
                    iCfCmu = zeros(D, D, T);
                    Cf(:, :, 1) = Ck;
                    iCk = inv(Ck);
                    for t = 2:T
                        idx = blockId{t-1};
                        piCk = sum(this_post(idx)) * iCk; %#ok
                        %piCk = iCk * post(k,t-1);
                        
                        % hacky, for now just hopping along the time axis
                        iCfCmu(:, :, t - 1) = inv(Cf(:, :, t - 1) + Cmu);
                        Cf(:, :, t) = inv(iCfCmu(:, :, t - 1) + piCk);
                        muk(:, t) = Cf(:, :, t) * (iCfCmu(:, :, t - 1) * muk(:, t - 1) + ...
                            (iCk * Y(:, idx)) * this_post(idx)');
                        %muk(:, t) = Cf(:, :, t) * (iCfCmu(:, :, t - 1) * muk(:, t - 1) + ...
                        %    piCk * Y(:, t-1));
                    end
                    
                    assert(~any(isnan(muk(:))), 'Got nan');
                    % Backward step for updating the means (Eq. 10)
                    for t = T-1 : -1 : 1
                        muk(:, t) = muk(:, t) + Cf(:, :, t) * (iCfCmu(:, :, t) * (muk(:, t + 1) - muk(:, t)));
                    end
                    
                    % interpolate muk back out
                    muk_interp = interp1(self.model.mu_t,muk',self.ttrain,'linear','extrap')';

                    % Update observation covariance (Eq. 11)
                    Ymu = Y - muk_interp;
                    Ck = (bsxfun(@times, this_post, Ymu) * Ymu') / sum(this_post);
                    Ck = Ck + eye(D) * self.params.CovRidge; % add ridge to regularize
                    
                    % Estimate (unnormalized) probabilities
                    pk(k, :) = MoKsm.mixtureDistribution(Ymu, Ck + Cmu, df);
                    post(k, :) = priors(k) * pk(k, :);
                    
                    mu(:, :, k) = muk;
                    C(:, :, k) = Ck;
                end
                
                % calculate log-likelihood
                p = sum(post, 1);
                %logLike(end + 1) = sum(MoKsm.mylog(p)); %#ok
                model = MoKsm.collect(mu, C, Cmu, priors, post, pk, logLike, mu_t, df);
                logLike(end+1) = self.evalTestSet();
                if self.params.Verbose
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
            time = toc;
            %fprintf('Time per iter per cluster : %g\n', time / iter / K);
            
            self.model = MoKsm.collect(mu, C, Cmu, priors, post, pk, logLike, mu_t, df);
        end

        function [self, success] = tryMerge(self)
            tol = self.params.Tolerance;
            verbose = self.params.Verbose;
            
            success = false;
            cands = MoKsm.getMergeCandidates(self.model.post);
            logLikeTest = self.evalTestSet();
            for ij = cands'
                try
                    fprintf('Trying to merge clusters %d and %d ', ij(1), ij(2))
                    newSelf = self;
                    newSelf.model = MoKsm.mergeClusters(newSelf.model, ij(1), ij(2));
                    newSelf = fullE(newSelf);
                    newSelf = fullEM(newSelf);
                    newLogLikeTest = newSelf.evalTestSet();
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
            Ytrain = self.Ytrain;
            Ytest = self.Ytest;
            ttest = self.ttest;
            ttrain = self.ttrain;
            model = self.model;
            tol = self.params.Tolerance;
            verbose = self.params.Verbose;
            
            success = false;
            splitCands = MoKsm.getSplitCandidates(self.model.post, model.pk, model.priors);
            logLikeTest = self.evalTestSet();
            for i = splitCands'
                try
                    fprintf('Trying to split cluster %d ', i)
                    newSelf = self;
                    newSelf.model = MoKsm.splitCluster(self.model, i);
                    newSelf = fullE(newSelf);
                    newSelf = fullEM(newSelf);
                    newLogLikeTest = newSelf.evalTestSet();
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
        
        function self = fullE(self)
            % Do one full E-step
            
            [~, T, K] = size(self.model.mu);
            for k = 1 : K
                % interpolate muk back out
                muk_interp = interp1(self.model.mu_t,self.model.mu(:, :, k)',self.ttrain,'linear','extrap')';

                self.model.pk(k, :) = MoKsm.mixtureDistribution(self.Ytrain - muk_interp, ...
                    self.model.C(:, :, k) + self.model.Cmu, self.model.df);
                self.model.post(k, :) = self.model.priors(k) * self.model.pk(k, :);
            end
            p = sum(self.model.post, 1);
            self.model.logLike(end + 1) =  self.evalTestSet();
            self.model.post = bsxfun(@rdivide, self.model.post, p);
            self.model.post(:, p == 0) = 0;
            self.model.priors = sum(self.model.post, 2) / T;
        end
                
        function plot(self, varargin)
            MoKsm.plotModel(self.model, self.Ytrain, self.ttrain, varargin{:})
        end
        
    end
    
    
    methods (Access = private)
        
        function [logLike, p] = evalTestSet(self)
            % Evaluate log-likelihood on test set by interpolating cluster means from
            % training set.
            
            K = size(self.model.mu, 3);
            T = numel(self.ttest);
            p = zeros(K, T);
            for k = 1 : K
                muk = interp1(self.model.mu_t, self.model.mu(:, :, k)', self.ttest, 'linear', 'extrap')';
                Ymu = self.Ytest - muk;
                p(k, :) = self.model.priors(k) * MoKsm.mixtureDistribution(Ymu, self.model.C(:, :, k) + self.model.Cmu, self.model.df);
            end
            logLike = mean(MoKsm.mylog(sum(p, 1)));
            logLike = logLike - K * self.params.ClusterCost;
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
        
        function model = splitCluster(model, k)
            % Split cluster k
            
            [mu, C, Cmu, priors, post, pk, logLike, mu_t, df] = MoKsm.expand(model);
            [D, ~, K] = size(mu);
            deltaMu  = chol(C(:, :, k))' * randn(D, 1);
            mu(:, :, k) = bsxfun(@plus, mu(:, :, k), deltaMu);
            mu(:, :, K + 1) = bsxfun(@minus, mu(:, :, k), deltaMu);
            C(:, :, k) = det(C(:, :, k))^(1 / D) * eye(D);
            C(:, :, K + 1) = C(:, :, k);
            priors(k) = priors(k) / 2;
            priors(K + 1) = priors(k);
            model = MoKsm.collect(mu, C, Cmu, priors, post, pk, logLike, mu_t, df);
        end
        
        function model = mergeClusters(model, i, j)
            % Merge clusters i and j
            
            [mu, C, Cmu, priors, post, pk, logLike, mu_t, df] = MoKsm.expand(model);
            mu(:, :, i) = (priors(i) * mu(:, :, i) + priors(j) * mu(:, :, j)) / (priors(i) + priors(j));
            C(:, :, i) = (priors(i) * C(:, :, i) + priors(j) * C(:, :, j)) / (priors(i) + priors(j));
            priors(i) = priors(i) + priors(j);
            mu(:, :, j) = [];
            C(:, :, j) = [];
            priors(j) = [];
            post(j, :) = [];
            pk(j, :) = [];
            model = MoKsm.collect(mu, C, Cmu, priors, post, pk, logLike, mu_t, df);
        end
        
        function cand = getSplitCandidates(post, pk, priors)
            fk = bsxfun(@rdivide, post, sum(post, 2));
            Jsplit = sum(fk .* (MoKsm.mylog(fk) - MoKsm.mylog(pk)), 2);
            [~, cand] = sort(Jsplit, 'descend');
            [D, T] = size(post);
            cand = cand(priors(cand) * T > 4 * D); % don't split small clusters
        end
        
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
        end
        
        function [mu, C, Cmu, priors, post, pk, logLike, mu_t, df] = expand(model)
            mu = model.mu;
            C = model.C;
            Cmu = model.Cmu;
            priors = model.priors;
            post = model.post;
            pk = model.pk;
            logLike = model.logLike;
            mu_t = model.mu_t;
            df = model.df;
        end
        
        function model = collect(mu, C, Cmu, priors, post, pk, logLike, mu_t, df)
            model.mu = mu;
            model.C = C;
            model.Cmu = Cmu;
            model.priors = priors;
            model.post = post;
            model.pk = pk;
            model.logLike = logLike;
            model.mu_t = mu_t;
            model.df = df;
        end
        
        function plotModel(model, Y, t, d)
            if nargin < 4
                if size(Y, 1) > 3
                    d = [7 1 4 10]; % tetrodes
                else
                    d = 1 : 3;      % single electrodes
                end
            end
            [~, j] = max(model.post, [], 1);
            K = size(model.post, 1);
            c = lines;
            figure(2), clf, hold all
            hdl = zeros(1, K);
            
            dd = combnk(d, 2);
            N = size(dd, 1);
            for k = 1 : N
                subplot(2, N, k)
                cla
                hold on
                for i = 1:K
                    plot(Y(dd(k, 1), j == i), Y(dd(k, 2), j == i), '.', 'markersize', 1, 'color', c(i, :))
                    hdl(i) = plot(model.mu(dd(k, 1), :, i), model.mu(dd(k, 2), :, i), '-', 'color', c(i, :),'LineWidth',3);
                end
                xlim(quantile(Y(dd(k, 1),:), [0.001 0.999]));
                ylim(quantile(Y(dd(k, 2),:), [0.001 0.999]));
            end
            
            subplot(2, N, N + (1 : N))
            cla
            hold on
            for i = 1:K
                plot(t(j==i),Y(d(1), j == i), '.', 'markersize', 1, 'color', c(i, :))
                hdl(i) = plot(model.mu_t,model.mu(d(1), :, i), '-', 'color', c(i, :),'LineWidth',3);
            end
            legend(hdl, arrayfun(@(x) sprintf('Cluster %d', x), 1:K, 'UniformOutput', false))
            ylim(quantile(Y(d(1),:), [0.001 0.999]));

        end        
        
    end
end
