classdef MoK
    % Mixture of Kalman filters model adapted from Calabrese & Paninski
    % 2011.
    %
    % Alexander S. Ecker & R. James Cotton
    % 2012-07-04
    
    properties 
        params      % parameters for fitting
        mu_t        % times at which cluster means are updated
        mu          % cluster means
        C           % cluster covariances
        Cmu         % covariance of cluster mean drift
        priors      % cluster priors (mixing proportions)
        df          % degress of freedom for t distribution
    end
    
    properties (Access = private)
        logLike     % log-likelihood curve during fitting
        Y           % data
        t           % times
        blockId     % mapping from data point to time blocks
        spikeId     % spike ids of the training data for each time block
    end
    
    methods

        function self = MoK(Y, t, mu_t, mu, C, Cmu, priors, df, varargin)
            
            % model parameters
            self.mu_t = mu_t;
            self.mu = mu;
            self.C = C;
            self.Cmu = Cmu;
            self.priors = priors;
            self.df = df;
            
            % training data
            self.Y = Y;
            self.t = t;
            [~, self.blockId] = histc(t, self.mu_t);
            self.spikeId = arrayfun(@(x) find(self.blockId == x), 1 : numel(mu_t), 'UniformOutput', false);
            
            % parse optional parameters
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('Tolerance', 0.0002);     % tolerance for determining convergence
            p.addOptional('Verbose', false);        % verbose output
            p.addOptional('CovRidge', 1.5);         % independent variance added to cluster covariances
            p.addOptional('ClusterCost', 0.0025);   % penalizer for adding additional clusters
            p.parse(varargin{:});
            self.params = p.Results;
        end
        
        
        function self = EM(self, maxIter)
            % Run EM until convergence.
            
            if nargin < 2, maxIter = Inf; end
            iter = 0;
            logLikeBase = NaN;
            while iter < maxIter && (iter < 2 || (self.logLike(end) - self.logLike(end - 1)) / (self.logLike(end - 1) - logLikeBase) > self.params.Tolerance)
                
                if ~mod(iter, 5), fprintf('.'), end
                iter = iter + 1;

                % Perform E step
                like = self.likelihood();
                p = sum(like, 1);
                post = bsxfun(@rdivide, like, p);
                post(:, p == 0) = 0;

                % Perform M step
                N = size(self.Y, 2);
                [mu, C, Cmu, ~, df] = self.expand(); %#ok<*PROP>
                [D, T, K] = size(mu);
                Cf = zeros(D, D, T);
                for k = 1 : K
                    
                    postk = post(k, :);
                    muk = mu(:, :, k);
                    Ck = C(:, :, k);
                    
                    % Forward iteration for updating the means (Eq. 9)
                    iCfCmu = zeros(D, D, T);
                    Cf(:, :, 1) = Ck;
                    iCk = inv(Ck);
                    for tt = 2 : T
                        idx = self.spikeId{tt-1};
                        iCfCmu(:, :, tt - 1) = inv(Cf(:, :, tt - 1) + Cmu);
                        if isempty(idx)
                            Cf(:, :, tt) = Cf(:, :, tt - 1) + Cmu;
                            muk(:, tt) = Cf(:, :, tt) * (iCfCmu(:, :, tt - 1) * muk(:, tt - 1));
                        else
                            piCk = sum(postk(idx)) * iCk; %#ok
                            Cf(:, :, tt) = inv(iCfCmu(:, :, tt - 1) + piCk);
                            muk(:, tt) = Cf(:, :, tt) * (iCfCmu(:, :, tt - 1) * muk(:, tt - 1) + ...
                                (iCk * self.Y(:, idx)) * postk(idx)'); %#ok
                        end
                    end
                    
                    % Backward iteration for updating the means (Eq. 10)
                    for tt = T-1 : -1 : 1
                        muk(:, tt) = muk(:, tt) + Cf(:, :, tt) * (iCfCmu(:, :, tt) * (muk(:, tt + 1) - muk(:, tt)));
                    end
                    assert(~any(isnan(muk(:))), 'Got nan');
                    
                    % Update cluster covariances (Eq. 11)
                    Ymu = self.Y - muk(:, self.blockId);
                    Ck = (bsxfun(@times, postk, Ymu) * Ymu') / sum(postk);
                    Ck = Ck + eye(D) * self.params.CovRidge; % add ridge to regularize
                    
                    mu(:, :, k) = muk;
                    C(:, :, k) = Ck;
                end
                
                % update class priors
                priors = sum(post, 2) / N;
                
                % check for starvation
                if any(priors * N < 2 * D)
                    error('MoK:starvation', 'Component starvation: cluster %d', find(priors * N < 2 * D, 1))
                end
                
                self = self.collect(mu, C, Cmu, priors, df);
                
                % calculate log-likelihood
                self.logLike(end + 1) = mean(MoK.mylog(sum(like, 1)));
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
        
        
        function self = deleteCluster(self, k)
            % Delete cluster k.
            
            [mu, C, Cmu, priors, df] = self.expand();
            mu(:, :, k) = [];
            C(:, :, k) = [];
            priors(k) = [];
            priors = priors / sum(priors);
            self = self.collect(mu, C, Cmu, priors, df);            
        end

        
        function p = likelihood(self, Y, block)
            % likelihood of the data for each cluster
            %   p = likelihood(self) for the likelihood of the training set
            %
            %   p = likelihood(self, Y, tBlockId) for a test set. tBlockId
            %   is the index of the time block for each column in Y.
            
            if nargin < 2
                Y = self.Y;
                block = self.blockId;
            end
            [mu, C, Cmu, priors, df] = self.expand();
            K = numel(priors);
            p = zeros(K, size(Y, 2));
            for k = 1 : K
                muk = mu(:, block, k);
                p(k, :) = priors(k) * MoK.mixtureDistribution(Y - muk, C(:, :, k) + Cmu, df);
            end
        end
        
        
        function post = posterior(self, varargin)
            % posterior of class membership for each spike
            %   post = posterior(self) for the posteriors on the training
            %   set 
            %   
            %   post = posterior(self, Y, t) for a test set. tBlockId is
            %   the index of the time block for each column in Y.
            
            likelihood = self.likelihood(varargin{:});
            p = sum(likelihood, 1);
            post = bsxfun(@rdivide, likelihood, p);
            post(:, p == 0) = 0;
        end
        
        
        function logl = logLikelihood(self, varargin)
            % Penalized average log-likelihood
            %   logl = logLikelihood(self) for the penalized log-likelihood
            %   on the training set.
            %   
            %   logl = logLikelihood(self, Y, t) for a test set. tBlockId
            %   is the index of the time block for each column in Y.
            
            [D, ~, K] = size(self.mu);
            logl = mean(MoK.mylog(sum(self.likelihood(varargin{:}), 1)));
            logl = logl - K * D * self.params.ClusterCost;
        end
        
        
        function ids = cluster(self, varargin)
            % return cluster ids for all spikes
            %   ids = cluster(self) to cluster the training set.
            %
            %   ids = cluster(self, Y, t) to cluster a test set. tBlockId
            %   is the index of the time block for each column in Y.
            
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
            
            if nargin < 2
                if size(self.Y, 1) > 3
                    d = [7 1 4 10]; % tetrodes
                else
                    d = 1 : 3;      % single electrodes
                end
            end
            d(d > size(self.Y, 1)) = [];
            
            j = self.cluster();
            K = max(j);
            figure(2), clf, hold all
            c = lines;
            hdl = zeros(1, K);
            
            dd = combnk(d, 2);
            
            if size(self.Y, 1) > 1
                N = size(dd, 1);
                for k = 1 : N
                    subplot(2, N, k)
                    cla
                    hold on
                    for i = 1:K
                        plot(self.Y(dd(k, 1), j == i), self.Y(dd(k, 2), j == i), '.', 'markersize', 1, 'color', c(i, :))
                        hdl(i) = plot(self.mu(dd(k, 1), :, i), self.mu(dd(k, 2), :, i), '-', 'color', c(i, :),'LineWidth',3);
                    end
                    xlim(quantile(self.Y(dd(k, 1), :), [0.001 0.999]));
                    ylim(quantile(self.Y(dd(k, 2), :), [0.001 0.999]));
                end
                subplot(2, N, N + (1 : N))
            else
                subplot(111);
            end
            
            
            cla
            hold on
            for i = 1:K
                plot(self.t(j == i), self.Y(d(1), j == i), '.', 'markersize', 1, 'color', c(i, :))
                mu_t = self.mu_t(1 : end - 1) + min(diff(self.mu_t)) / 2;
                hdl(i) = plot(mu_t, self.mu(d(1), :, i), '-', 'color', c(i, :), 'LineWidth', 3);
            end
            legend(hdl, arrayfun(@(x) sprintf('Cluster %d', x), 1 : K, 'UniformOutput', false))
            ylim(quantile(self.Y(d(1),:), [0.001 0.999]));
        end
        
    end
    
    
    methods (Access = private)
        
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
    
    
    methods(Static)
        
        function y = mylog(x)
            % Natural logarithm excluding zeros
            
            y = reallog(x);
            y(x == 0) = 0;
        end
        
        
        function p = mixtureDistribution(X, C, df)
            % Probability distribution of mixture components (normal or t)
            
            if isinf(df)
                p = MoK.mvn(X, C);
            else
                p = MoK.mvt(X, C, df);
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
