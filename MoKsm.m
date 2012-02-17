classdef MoKsm < SpikeSortingHelper
    % MoK to 12d data (adapted from Calabrese & Paninski 2011)
    % AE 2012-02-03
    % JC 2012-02-15
    
    properties
        params = struct('MaxTrainSpikes', 5000, ...
            'MaxTestSpikes', 20000, ...
            'TrainFrac', 0.8, ...'
            'verbose', false, ...
            'tol', 0.0002, ...
            'Feature', 'Points');
        model = struct;
        Ytrain = [];
        Ytest = [];
        ttrain = [];
        ttest = [];
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
            
            t = t(:)';
            
            T = size(Y, 2);
            
            assert(size(Y,1) <= 50, 'Dimensionality way too high');
            
            % split into training & test data
            rnd = randperm(T);
            nTrain = fix(self.params.TrainFrac * T);
            train = sort(rnd(1 : min(nTrain,self.params.MaxTrainSpikes)));
            test = sort(rnd(nTrain + 1 : min(end, nTrain + self.params.MaxTestSpikes)));
            self.Ytrain = Y(:, train);
            self.Ytest = Y(:, test);
            self.ttrain = t(train);
            self.ttest = t(test);
            
            nTrain = length(train);
            T = length(train) + length(test);
            
            % Initialize model using 1 component
            fprintf('Running initial Kalman filter model with one cluster ')
            self.model.mu = repmat(mean(self.Ytrain, 2), [1 nTrain]); % cluster means
            self.model.C = cov(self.Ytrain');                                  % observation covariance
            self.model.Cmu = diag(median(abs(self.Ytrain), 2) / 0.6745 / 1000);  % cluster mean drift [TODO: make it dependent on average firing rate]
            self.model.priors = 1;                       % cluster weights
            self.model.post = ones(1, nTrain);
            self.model.pk = MoKsm.mvn(self.Ytrain - self.model.mu, self.model.C + self.model.Cmu);
            self.model.logLike = sum(MoKsm.mylog(self.model.pk));
            self.model = MoKsm.fullEM(self.Ytrain, self.model, self.params.tol, self.params.verbose);
            if self.params.verbose, plotData(self), end
            fprintf(' done (likelihood: %.5g)\n', self.model.logLike(end))
            
            % Run split & merge
            % We alternate between trying to split and trying to merge. Both are done
            % until no candidate leads to success. If both splitting and merging didn't
            % lead to success we terminate
            op = {@MoKsm.trySplit, @MoKsm.tryMerge};
            i = 1;
            success = true(1, 2);
            while any(success)
                [self.model, success(i)] = op{i}(self.Ytrain, self.Ytest, self.ttest, self.ttrain, self.model, self.params.tol, self.params.verbose);
                if ~success(i)
                    i = 3 - i;
                end
            end
            
            fprintf('Done with split & merge\n')
            
            fprintf('Performing fit on entire dataset ...');
            [self.model.logLike p] = MoKsm.evalTestSet(Y, t, self.ttrain, self.model)
            self.model.postAll = p;
            fprintf(' Done\n');
            
            fprintf('--\n')
            fprintf('Number of clusters: %d\n', size(self.model.mu, 3))
            fprintf('Log-likelihoods\n')
            fprintf('  training set: %.8g\n', self.model.logLike(end))
            fprintf('      test set: %.8g\n', MoKsm.evalTestSet(self.Ytest, self.ttest, self.ttrain, self.model))
            fprintf('\n\n')
            
            model.Ytrain = self.Ytrain;
            model.Ytest = self.Ytest;
            model.ttrain = self.ttrain;
            model.ttest = self.ttest;
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
        function [model, success] = tryMerge(Ytrain, Ytest, ttest, ttrain, model, tol, verbose)
            success = false;
            cands = MoKsm.getMergeCandidates(model.post);
            logLikeTest = MoKsm.evalTestSet(Ytest, ttest, ttrain, model);
            for ij = cands'
                fprintf('Trying to merge clusters %d and %d ', ij(1), ij(2))
                newModel = MoKsm.mergeClusters(model, ij(1), ij(2));
                newModel = MoKsm.fullE(Ytrain, newModel);
                newModel = MoKsm.fullEM(Ytrain, newModel, tol, verbose);
                newLogLikeTest = MoKsm.evalTestSet(Ytest, ttest, ttrain, newModel);
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
        end
        
        
        function [model, success] = trySplit(Ytrain, Ytest, ttest, ttrain, model, tol, verbose)
            
            success = false;
            splitCands = MoKsm.getSplitCandidates(model.post, model.pk, model.priors);
            logLikeTest = MoKsm.evalTestSet(Ytest, ttest, ttrain, model);
            for i = splitCands'
                try
                    fprintf('Trying to split cluster %d ', i)
                    newModel = MoKsm.splitCluster(model, i);
                    newModel = MoKsm.fullE(Ytrain, newModel);
                    newModel = MoKsm.fullEM(Ytrain, newModel, tol, verbose);
                    newLogLikeTest = MoKsm.evalTestSet(Ytest, ttest, ttrain, newModel);
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
                        rethrow(err)
                    end
                end
            end
        end
        
        function model = fullEM(Y, model, tol, verbose)
            % EM recursion
            
            [mu, C, Cmu, priors, post, pk, logLike] = MoKsm.expand(model);
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
                    pk(k, :) = MoKsm.mvn(Ymu, Ck + Cmu);
                    post(k, :) = priors(k) * pk(k, :);
                    
                    mu(:, :, k) = muk;
                    C(:, :, k) = Ck;
                end
                
                % calculate log-likelihood
                p = sum(post, 1);
                logLike(end + 1) = sum(MoKsm.mylog(p)); %#ok
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
            
            model = MoKsm.collect(mu, C, Cmu, priors, post, pk, logLike);
        end
        
        function model = fullE(Y, model)
            % Do one full E-step
            
            [~, T, K] = size(model.mu);
            for k = 1 : K
                model.pk(k, :) = MoKsm.mvn(Y - model.mu(:, :, k), model.C(:, :, k) + model.Cmu);
                model.post(k, :) = model.priors(k) * model.pk(k, :);
            end
            p = sum(model.post, 1);
            model.logLike(end + 1) = sum(MoKsm.mylog(p));
            model.post = bsxfun(@rdivide, model.post, p);
            model.post(:, p == 0) = 0;
            model.priors = sum(model.post, 2) / T;
        end
        
        
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
                p(k, :) = model.priors(k) * MoKsm.mvn(Ymu, model.C(:, :, k) + model.Cmu);
            end
            logLike = sum(MoKsm.mylog(sum(p, 1)));
        end
        
        function model = splitCluster(model, k)
            % Split cluster k
            
            [mu, C, Cmu, priors, post, pk, logLike] = MoKsm.expand(model);
            [D, ~, K] = size(mu);
            deltaMu  = chol(C(:, :, k))' * randn(D, 1);
            mu(:, :, k) = bsxfun(@plus, mu(:, :, k), deltaMu);
            mu(:, :, K + 1) = bsxfun(@minus, mu(:, :, k), deltaMu);
            C(:, :, k) = det(C(:, :, k))^(1 / D) * eye(D);
            C(:, :, K + 1) = C(:, :, k);
            priors(k) = priors(k) / 2;
            priors(K + 1) = priors(k);
            model = MoKsm.collect(mu, C, Cmu, priors, post, pk, logLike);
        end
        
        function model = mergeClusters(model, i, j)
            % Merge clusters i and j
            
            [mu, C, Cmu, priors, post, pk, logLike] = MoKsm.expand(model);
            mu(:, :, i) = (priors(i) * mu(:, :, i) + priors(j) * mu(:, :, j)) / (priors(i) + priors(j));
            C(:, :, i) = (priors(i) * C(:, :, i) + priors(j) * C(:, :, j)) / (priors(i) + priors(j));
            priors(i) = priors(i) + priors(j);
            mu(:, :, j) = [];
            C(:, :, j) = [];
            priors(j) = [];
            post(j, :) = [];
            pk(j, :) = [];
            model = MoKsm.collect(mu, C, Cmu, priors, post, pk, logLike);
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
        
        function [mu, C, Cmu, priors, post, pk, logLike] = expand(model)
            mu = model.mu;
            C = model.C;
            Cmu = model.Cmu;
            priors = model.priors;
            post = model.post;
            pk = model.pk;
            logLike = model.logLike;
        end
        
        function model = collect(mu, C, Cmu, priors, post, pk, logLike)
            model.mu = mu;
            model.C = C;
            model.Cmu = Cmu;
            model.priors = priors;
            model.post = post;
            model.pk = pk;
            model.logLike = logLike;
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
        end        
        
    end
end
