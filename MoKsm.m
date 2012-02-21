classdef MoKsm < SpikeSortingHelper
    % MoK to 12d data (adapted from Calabrese & Paninski 2011)
    % AE 2012-02-03
    % JC 2012-02-15
    
    properties
        params = struct('MaxTrainSpikes', 5000, ...20000, ...
            'MaxTestSpikes', 5000, ... 20000, ...
            'TrainFrac', 0.8, ...'
            'verbose', false, ...
            'tol', 0.0002, ...
            'Feature', 'Points', ...
            'dTmu',5000, ...
            'driftRate', 1 / 5000000 ...
            );
        model = struct;
        Ytrain = [];
        Ytest = [];
        ttrain = [];
        ttest = [];
        blockId = []; % maps the training data to mu's
    end
    
    methods
        % Smart initialization to detect either explict data to use
        % or the standard data format
        function self = MoKsm(data,varargin)
            
            %             if ismatrix(data) && numel(varargin) > 0 && ...
            %                     any(size(data) == numel(varargin{1}))
            %                 warning('Passing in data this way will be deprecated');
            %                 if numel(varargin) > 1
            %                     self = parseParams(self, varargin{2:end});
            %                 end
            %                 self = initializeRaw(self, data, varargin{1});
            %             else
            self = self@SpikeSortingHelper(data, varargin{:});
            % Hacky.  For creating with just features and time
            % need to skip this
            if length(varargin) > 1
                self = parseParams(self, varargin{:});
            end
            self = getFeatures(self, 'Points');
        end
        
        % Process the input arguments into the structure
        function self = parseParams(self, varargin)
            self.params = parseVarArgs(self.params, varargin{:});
        end
        
        % Fit the model
        function self = fitModel(self)
            if nargin < 3, verbose = false; end
            randn('state', 1)
            rand('state', 1)
            
            % warning off MATLAB:nearlySingularMatrix
            
            t = self.data.SpikeTimes.data;
            Y = self.data.Features.data;
            
            if size(Y,1) == length(t)
                Y = Y';
            elseif size(Y,2) ~= length(t)
                error('Time dimension doesn''t match dataset');
            end
            
            t = t(:)'; % force time to row vector
            
            assert(size(Y,1) <= 50, 'Dimensionality way too high');
            
            % split into training & test data
            T = size(Y,2);
            rnd = randperm(T);
            nTrain = fix(self.params.TrainFrac * T);
            train = sort(rnd(1 : min(nTrain,self.params.MaxTrainSpikes)));
            test = sort(rnd(nTrain + 1 : min(end, nTrain + self.params.MaxTestSpikes)));
            self.Ytrain = Y(:, train);
            self.Ytest = Y(:, test);
            self.ttrain = t(train);
            self.ttest = t(test);
            
            % Initialize model using 1 component
            self.model.mu_t = t(1):self.params.dTmu:t(end);
            nTime = length(self.model.mu_t);

            % Assign spikes to blocks for the M step
            blockId = zeros(size(self.ttrain));
            lastIdx = 0;
            for i = 1:length(self.model.mu_t)-1
                nextIdx = lastIdx + find(abs(self.ttrain(lastIdx+1:end)-self.model.mu_t(i)) > ...
                    abs(self.ttrain(lastIdx+1:end)-self.model.mu_t(i+1)),1,'first');
                blockId(lastIdx+1:nextIdx) = i;
                lastIdx = nextIdx;
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
            self.model.Cmu = diag(median(abs(bsxfun(@minus,self.Ytrain,mean(self.Ytrain,2))),2).^2 / 0.6745^2) * self.params.dTmu * self.params.driftRate;
            self.model.priors = 1;                       % cluster weights
            self.model.post = ones(1, size(self.Ytrain,2));
            self.model.pk = MoKsm.mvn(bsxfun(@minus,self.Ytrain,mean(self.Ytrain,2)), self.model.C + self.model.Cmu);
            self.model.logLike = sum(MoKsm.mylog(self.model.pk));
            
            
            self = fullEM(self);
            if self.params.verbose, plot(self), end
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
            [self.model.logLike p] = MoKsm.evalTestSet(Y, t, self.model);
            self.model.postAll = p;
            fprintf(' Done\n');
            
            fprintf('--\n')
            fprintf('Number of clusters: %d\n', size(self.model.mu, 3))
            fprintf('Log-likelihoods\n')
            fprintf('  training set: %.8g\n', self.model.logLike(end))
            fprintf('      test set: %.8g\n', MoKsm.evalTestSet(self.Ytest, self.ttest, self.model))
            fprintf('\n\n')
            
            model.Ytrain = self.Ytrain;
            model.Ytest = self.Ytest;
            model.ttrain = self.ttrain;
            model.ttest = self.ttest;
        end

        function self = fullEM(self)
            % warning off MATLAB:nearlySingularMatrix
            
            Y = self.Ytrain;
            blockId = self.blockId;
            
            % EM recursion
            [mu, C, Cmu, priors, post, pk, logLike, mu_t] = MoKsm.expand(self.model);
            [D, T, K] = size(mu);
            Cf = zeros(D, D, T);                            % state covariances
            iCfCmu = zeros(D, D, T);
            iter = 0;
            logLikeBase = logLike(end);
            tic
            while iter < 2 || (logLike(end) - logLike(end - 1)) / (logLike(end - 1) - logLikeBase) > self.params.tol
                
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
                    for t = 2:T
                        idx = blockId{t-1};
                        piCk = sum(post(k, idx)) * iCk; %#ok
                        %piCk = iCk * post(k,t-1);
                        
                        % hacky, for now just hopping along the time axis
                        iCfCmu(:, :, t - 1) = inv(Cf(:, :, t - 1) + Cmu);
                        Cf(:, :, t) = inv(iCfCmu(:, :, t - 1) + piCk);
                        muk(:, t) = Cf(:, :, t) * (iCfCmu(:, :, t - 1) * muk(:, t - 1) + ...
                            (iCk * Y(:, idx)) * post(k,idx)');
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
                    Ck = (bsxfun(@times, post(k, :), Ymu) * Ymu') / sum(post(k, :));
                    
                    % Estimate (unnormalized) probabilities
                    pk(k, :) = MoKsm.mvn(Ymu, Ck + Cmu);
                    post(k, :) = priors(k) * pk(k, :);
                    
                    mu(:, :, k) = muk;
                    C(:, :, k) = Ck;
                end
                
                % calculate log-likelihood
                p = sum(post, 1);
                %logLike(end + 1) = sum(MoKsm.mylog(p)); %#ok
                model = MoKsm.collect(mu, C, Cmu, priors, post, pk, logLike, mu_t);
                logLike(end+1) = MoKsm.evalTestSet(self.Ytrain, self.ttrain, model);
                if self.params.verbose
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
            
            self.model = MoKsm.collect(mu, C, Cmu, priors, post, pk, logLike, mu_t);
        end

        function [self, success] = tryMerge(self)
            tol = self.params.tol;
            verbose = self.params.verbose;
            
            success = false;
            cands = MoKsm.getMergeCandidates(self.model.post);
            logLikeTest = MoKsm.evalTestSet(self.Ytest, self.ttest, self.model);
            for ij = cands'
                try
                    fprintf('Trying to merge clusters %d and %d ', ij(1), ij(2))
                    newSelf = self;
                    newSelf.model = MoKsm.mergeClusters(newSelf.model, ij(1), ij(2));
                    newSelf = fullE(newSelf);
                    newSelf = fullEM(newSelf);
                    newLogLikeTest = MoKsm.evalTestSet(newSelf.Ytest, newSelf.ttest, newSelf.model);
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
            tol = self.params.tol;
            verbose = self.params.verbose;
            
            success = false;
            splitCands = MoKsm.getSplitCandidates(self.model.post, model.pk, model.priors);
            logLikeTest = MoKsm.evalTestSet(Ytest, ttest, model);
            for i = splitCands'
                try
                    fprintf('Trying to split cluster %d ', i)
                    newSelf = self;
                    newSelf.model = MoKsm.splitCluster(self.model, i);
                    newSelf = fullE(newSelf);
                    newSelf = fullEM(newSelf);
                    newLogLikeTest = MoKsm.evalTestSet(Ytest, ttest, newSelf.model);
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
        
        function self = fullE(self);
            % Do one full E-step
            
            [~, T, K] = size(self.model.mu);
            for k = 1 : K
                % interpolate muk back out
                muk_interp = interp1(self.model.mu_t,self.model.mu(:, :, k)',self.ttrain,'linear','extrap')';

                self.model.pk(k, :) = MoKsm.mvn(self.Ytrain - muk_interp, ...
                    self.model.C(:, :, k) + self.model.Cmu);
                self.model.post(k, :) = self.model.priors(k) * self.model.pk(k, :);
            end
            p = sum(self.model.post, 1);
            self.model.logLike(end + 1) =  MoKsm.evalTestSet(self.Ytrain, self.ttrain, self.model);
            self.model.post = bsxfun(@rdivide, self.model.post, p);
            self.model.post(:, p == 0) = 0;
            self.model.priors = sum(self.model.post, 2) / T;
        end
                
        function plot(self)
            MoKsm.plotModel(self.model,self.Ytrain,self.ttrain)
        end
        
    end
    
    methods(Static)
        function y = mylog(x)
            % Natural logarithm excluding zeros
            
            y = reallog(x);
            y(x == 0) = 0;
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
        
        function [logLike p] = evalTestSet(Ytest, ttest, model)
            % Evaluate log-likelihood on test set by interpolating cluster means from
            % training set.
            
            clusterCost = 0.001;
            
            K = size(model.mu, 3);
            T = numel(ttest);
            p = zeros(K, T);
            for k = 1 : K
                muk = interp1(model.mu_t, model.mu(:, :, k)', ttest,'linear','extrap')';
                Ymu = Ytest - muk;
                p(k, :) = model.priors(k) * MoKsm.mvn(Ymu, model.C(:, :, k) + model.Cmu);
            end
            logLike = mean(MoKsm.mylog(sum(p, 1)));
            %logLike = logLike - K * clusterCost;
        end
        
        function model = splitCluster(model, k)
            % Split cluster k
            
            [mu, C, Cmu, priors, post, pk, logLike, mu_t] = MoKsm.expand(model);
            [D, ~, K] = size(mu);
            deltaMu  = chol(C(:, :, k))' * randn(D, 1);
            mu(:, :, k) = bsxfun(@plus, mu(:, :, k), deltaMu);
            mu(:, :, K + 1) = bsxfun(@minus, mu(:, :, k), deltaMu);
            C(:, :, k) = det(C(:, :, k))^(1 / D) * eye(D);
            C(:, :, K + 1) = C(:, :, k);
            priors(k) = priors(k) / 2;
            priors(K + 1) = priors(k);
            model = MoKsm.collect(mu, C, Cmu, priors, post, pk, logLike, mu_t);
        end
        
        function model = mergeClusters(model, i, j)
            % Merge clusters i and j
            
            [mu, C, Cmu, priors, post, pk, logLike, mu_t] = MoKsm.expand(model);
            mu(:, :, i) = (priors(i) * mu(:, :, i) + priors(j) * mu(:, :, j)) / (priors(i) + priors(j));
            C(:, :, i) = (priors(i) * C(:, :, i) + priors(j) * C(:, :, j)) / (priors(i) + priors(j));
            priors(i) = priors(i) + priors(j);
            mu(:, :, j) = [];
            C(:, :, j) = [];
            priors(j) = [];
            post(j, :) = [];
            pk(j, :) = [];
            model = MoKsm.collect(mu, C, Cmu, priors, post, pk, logLike, mu_t);
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
        
        function [mu, C, Cmu, priors, post, pk, logLike, mu_t] = expand(model)
            mu = model.mu;
            C = model.C;
            Cmu = model.Cmu;
            priors = model.priors;
            post = model.post;
            pk = model.pk;
            logLike = model.logLike;
            mu_t = model.mu_t;
        end
        
        function model = collect(mu, C, Cmu, priors, post, pk, logLike, mu_t)
            model.mu = mu;
            model.C = C;
            model.Cmu = Cmu;
            model.priors = priors;
            model.post = post;
            model.pk = pk;
            model.logLike = logLike;
            model.mu_t = mu_t;
        end
        
        function plotModel(model,Y,t)
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
                hdl(i) = plot(model.mu(d1, :, i), model.mu(d2, :, i), '-', 'color', c(i, :),'LineWidth',3);
            end
            legend(hdl, arrayfun(@(x) sprintf('Cluster %d', x), 1:K, 'UniformOutput', false))
            xlim(quantile(Y(d1,:),[0.001 0.999]));
            ylim(quantile(Y(d2,:),[0.001 0.999]));
            
            subplot(212)
            cla
            hold on
            for i = 1:K
                plot(t(j==i),Y(d1, j == i), '.', 'markersize', 1, 'color', c(i, :))
                hdl(i) = plot(model.mu_t,model.mu(d1, :, i), '-', 'color', c(i, :),'LineWidth',3);
            end
            legend(hdl, arrayfun(@(x) sprintf('Cluster %d', x), 1:K, 'UniformOutput', false))
            ylim(quantile(Y(d1,:),[0.001 0.999]));

        end        
        
    end
end
