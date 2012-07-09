classdef ClusteringHelper

	properties
        ClusterAssignment = [];
        GroupingAssignment = [];
        ClusterTags = [];
        ContaminationMatrix = [];
	end

    methods (Abstract)
        % This method updates the ClusterAssignment and ContaminationMatrix
        self = updateInformation(self);
        
        % fit the model
        self = fit(self);

        % Delete a cluster by ID.  Should call updateInformation
        % afterwards.
        self = delete(self, id);

        % Splits a cluster by ID.  Should call updateInformation
        % afterwards.
        self = split(self, id);

        % Splits a cluster by ID.  Should call updateInformation
        % afterwards.
        self = merge(self, ids);

        % Group the selected clusters
        self = group(self, ids);

        % Refit the complete data set again
        self = refit(self);
        
        % Remove any information that can be recomputed and doesn't
        % need to be stored on disk
        self = compress(self);
        
        % Recreate any information that compress strips out
        self = uncompress(self);
    end

	methods
		function self = ClusteringHelper()
        end
        
        function s = saveStructure(self)
            warning off MATLAB:structOnObject
            s = struct(self);
            warning on MATLAB:structOnObject
        end
        
        function self = singleUnit(self, clusterIds)
            % Toggles the SingleUnit flag for the selected ids
            for i = 1:length(clusterIds)
                tags = self.ClusterTags.data{clusterIds(i)};
                su_flags = strcmp(tags,'SingleUnit');
                if ~any(su_flags) % If flag not found, add it
                    self.ClusterTags.data{clusterIds(i)}{end+1} = 'SingleUnit';
                else              % Delete any found flags
                    self.ClusterTags.data{clusterIds(i)}(su_flags) = [];
                end
            end
        end
        
        function cm = getContamination(self, ids)
            % Returns the full or a subset of the contamination matrix
            
            groups = self.GroupingAssignment.data;
            if nargin < 2
                ids = 1 : numel(groups);
            end
            groups = groups(ids);
            
            pairwise = self.ContaminationMatrix.data.pairwise;
            n = self.ContaminationMatrix.data.n;
            K = numel(ids);
            cm = zeros(K);
            for i = 1 : K
                for j = 1 : K
                    pij = pairwise(groups{i}, groups{j});
                    if i == j
                        cm(i, j) = 1 - sum(pij(:)) / sum(n(groups{i}));
                    else
                        cm(i, j) = sum(pij(:)) / sum(n(groups{i}));
                    end
                end
            end
        end
        
        function [ids group] = getClusterIds(self)            
            ids = 1:length(self.GroupingAssignment.data);
            group = cellfun(@(x) length(x) > 1, self.GroupingAssignment.data);
        end
                        
        function spikeIds = getSpikesByClusIds(self,clusterIds)
            modelIds = cat(2,self.GroupingAssignment.data{clusterIds});
            spikeIds = cat(2,self.ClusterAssignment.data{modelIds});
            
            % Replace any groups with their spike ids
            while any(spikeIds < 0)
                groups = find(spikeIds < 0);
                groupIds = unique(spikeIds(groups));
                
                % Remove the negative numbers and replace with their sets
                % of spike ids
                spikeIds(groups) = [];                   %#ok<*AGROW>
                spikeIds = [spikeIds, cat(2,self.ClusterAssignment.data{groupIds})];
            end
            
            spikeIds = sort(unique(spikeIds));
        end
        
        % TODO: Probably there should be a private method for actually
        % adding and deleting clusters from updateInformation so the
        % tags follow the clusters around properly
    
        function [corrs time] = getCrossCorrs(self, varargin)
            params.clusIds = getClusterIds(self);
            params.binSize = 0.5;
            params.nBins = 80;
            params = parseVarArgs(params,varargin{:});
            
            for i = 1:length(params.clusIds)
                ids1 = getSpikesByClusIds(self,params.clusIds(i));
                t1 = getSpikeTimes(self, ids1);
                
                for j = i:length(params.clusIds)
                    ids2 = getSpikesByClusIds(self,params.clusIds(j));
                    t2 = getSpikeTimes(self,ids2);
                    % units for time in CrossCorr are 1/10000 sec!!
                    % indices are such that CCGs are relative to the spikes
                    % of the cell in the column
                    [corrs{j,i} time] = CrossCorr(10 * t1, 10 * t2, params.binSize, params.nBins);
                    if i == j
                        corrs{j,i}(params.nBins / 2 + 1) = 0;
                    end
                    corrs{i,j} = flipud(corrs{j,i});
                end
            end
        end
        
        function bool = hasTag(self,tag,varargin)
            % Return if a cluster has a tag
            %
            % hasTag(data, tag, varargin)
            %     clusIds [getActiveClusters(data)]
            %
            % JC 2009-09-27
            
            params.clusIds = getClusterIds(self);
            params = parseVarArgs(params,varargin{:});
            
            if isprop(self,'ClusterTags') && isfield(self.ClusterTags,'data')
                bool = cellfun(@(x) any(strcmp(x,tag)), self.ClusterTags.data);
            else
                bool = zeros(1,length(params.clusIds));
            end
            bool = reshape(bool,1,[]);
        end
        
        function varargout = plotProjections(self,varargin)
            % Plot the prjoections
            %
            % [axHdl, plotHdl] = plotProjections(data, varargin)
            %     clusIds [getActiveClusters(data)]
            %     maxPoints [20000]
            %     scatter [0] - use scatter plot so colors more mixed, slower
            %
            % JC 2009-09-27
            
            params.maxPoints = 5000;
            params.clusIds = getClusterIds(self);
            params.pointSize = 1;
            params.scatter = 0;
            params.position = [];
            params = parseVarArgs(params,varargin{:});
            
            if isempty(params.position)
                params.position = [0 0 1 1];
            end
            
            X = self.Features.data;
            color = getClusColor(self, params.clusIds);
            
            % select features to plot
            if size(X, 2) > 3
                feat = 1 : 3 : 10;      % tetrodes
            else
                feat = 1 : size(X, 2);  % single channel data
            end
            c = combnk(feat, 2);
            N = size(c, 1);
            n = fix(sqrt(N));
            m = ceil(N / n);
            
            % select subset of points to show
            K = numel(params.clusIds);
            show = cell(1, K);
            for k = 1 : K
                ids = getSpikesByClusIds(self, params.clusIds(k));
                r = randperm(numel(ids));
                prior = sum(self.priors(self.GroupingAssignment.data{params.clusIds(k)}));
                show{k} = ids(r(1 : min(end, round(params.maxPoints * prior))));
            end

            % plot 2D scatter plots
            axHdl = zeros(1, N);
            plotHdl = zeros(K, N);
            ij = 1;
            for i = 1 : n
                for j = 1 : m
                    pos = params.position;
                    pos(1) = pos(1) + (i - 1) * pos(3) / n;
                    pos(2) = pos(2) + (m - j) * pos(4) / m;
                    pos(3) = 0.98 * pos(3) / n;
                    pos(4) = 0.98 * pos(4) / m;
                    axHdl(ij) = axes('Position', pos, 'Color', 0.5 * ones(1, 3)); %#ok
                    hold on
                    ax = Inf * [1, -1, 1, -1];
                    for k = 1 : K
                        X1 = X(show{k}, c(ij, 1));
                        X2 = X(show{k}, c(ij, 2));
                        plotHdl(k, ij) = plot(X1, X2, '.', 'MarkerSize', params.pointSize, 'Color', color(k, :));
                        ax = [ClusteringHelper.updateLimits(ax(1 : 2), X1), ...
                              ClusteringHelper.updateLimits(ax(3 : 4), X2)];
                    end
                    axis(ax)
                    set(axHdl(ij), 'box', 'on', 'xtick', [], 'ytick', [])
                    ij = ij + 1;
                end
            end
            if nargout > 0, varargout{1} = axHdl; end
            if nargout > 1, varargout{2} = plotHdl; end
        end
        
        function varargout = plotContaminations(self,varargin)
            % Plot the contamination matrix
            %
            % [axesHdl, imgHdl] = plotContaminations(data, varargin)
            %     clusIds [getActiveClusters(data)]
            %     maxPoints [20]
            %
            % JC 2009-09-27
            
            params.clusIds = getClusterIds(self);
            params = parseVarArgs(params,varargin{:});
            
            cla
            hdl = imagesc(getContamination(self, params.clusIds));
            xlabel('Source'); ylabel('Classified');
            
            set(gca,'XTick',1:length(params.clusIds));
            set(gca,'XTickLabel',params.clusIds);
            set(gca,'YTick',1:length(params.clusIds));
            set(gca,'YTickLabel',params.clusIds);
            
            set(gca,'CLim',[0 0.2]);
            
            if nargout > 0, varargout{1} = gca; end
            if nargout > 1, varargout{2} = hdl; end
            
        end
        
        function varargout = plotWaveforms(self, varargin)
            % Plot the waveforms
            %
            % plotWaveforms(data, varargin)
            %     clusIds [getActiveClusters(data)]
            %     maxPoints [20]
            %
            % JC 2009-09-27
            
            params.maxPoints = 50;
            params.clusIds = getClusterIds(self);
            params = parseVarArgs(params,varargin{:});
            
            k = numel(params.clusIds);
            chans = numel(self.Waveforms.data);
            hdl = zeros(k, chans);
            yl = [-100 100];
            
            for i = 1 : k
                color = getClusColor(self, params.clusIds(i));
                ids = getSpikesByClusIds(self, params.clusIds(i));
                if isfield(self.Waveforms.meta, 'subset')
                    [~, ids] = ismember(ids, self.Waveforms.meta.subset);
                    ids = ids(ids > 0);
                end
                r = randperm(length(ids));
                ids = sort(ids(r(1 : min(end, params.maxPoints))));
                for j = 1 : chans
                    hdl(i, j) = axes('Position', [(i - 1) / k, (j - 1) / chans, 1 / k, 1 / chans]);
                    if isempty(ids), axis off, continue; end
                    plot(self.Waveforms.data{j}(:, ids), 'Color', color)
                    hold on
                    axis tight off
                    m = median(self.Waveforms.data{j}(:, ids), 2);
                    plot(m, 'k', 'linewidth', 2)
                    yl = [min(yl(1), 1.5 * min(m)), max(yl(2), 1.5 * max(m))];
                end
            end
            
            linkaxes(hdl(:));
            set(hdl(:), 'YLim', yl);
            
            if nargout
                varargout{1} = hdl;
            end
        end
        
        function plotTimeFeatures(self, varargin)
            % Plot the features against time (for stability)
            %
            % plotTimeFeatures(data, varargin)
            %     clusIds [getActiveClusters(data)]
            %
            % JC 2009-09-29
            
            params.clusIds = getClusterIds(self);
            params = parseVarArgs(params,varargin{:});
            
            assert(length(params.clusIds) >= 1, 'Can only do for multiple clusters');
            
            clf
            X = self.Features.data;
            
            if size(X, 2) > 3
                feat = 1 : 3 : 10;
            else
                feat = 1 : size(X, 2);
            end
            n = numel(feat);
            
            for i = 1 : n
                subplot(n, 1, i);
                hold on
            end
            
            yl = [Inf(n, 1), -Inf(n, 1)];
            for i = 1 : length(params.clusIds)
                color = getClusColor(self, params.clusIds(i));
                ids = getSpikesByClusIds(self, params.clusIds(i));
                for j = 1 : n
                    subplot(n, 1, j);
                    Xij = X(ids, feat(j));
                    plot(self.SpikeTimes.data(ids), Xij, '.', 'color', color, 'markersize', 1);
                    axis tight
                    yl(j, :) = ClusteringHelper.updateLimits(yl(j, :), Xij);
                    ylim(yl(j, :))
                end
            end
        end
        
        function varargout = plotTimeAxes(self, handle)
            % Plot time axes (largest variance feature against time)
            %   plotTimeAxes(self) plots in gca.
            %   plotTimeAxes(self, handle) plots in the given handle.
            %       handle:   axes handle to plot in
            %
            % AE 2012-07-09
            
            if nargin < 2, handle = gca; end
            
            [~, ndx] = max(var(self.Features.data)); % largest variance feature
            clusIds = getClusterIds(self);
            n = numel(clusIds);
            plotAxes = zeros(1, n);
            xl = [Inf -Inf];
            set(handle, 'NextPlot', 'add', 'YGrid', 'on', 'Box', 'on')
            for i = 1 : n
                color = getClusColor(self, clusIds(i));
                ids = getSpikesByClusIds(self, clusIds(i));
                x = self.Features.data(ids, ndx);
                t = self.SpikeTimes.data(ids) / 1000 / 60; % convert to minutes
                plotAxes(i) = plot(handle, x, t, '.', 'Color', color, 'MarkerSize', 1);
                xl = ClusteringHelper.updateLimits(xl, x);
            end
            axis(handle, 'tight')
            xlim(xl)
            
            if nargout
                varargout{1} = plotAxes;
            end
        end
        
        function varargout = plotCrossCorrs(self, varargin)
            % Plot cross-correlograms
            %
            % hdl = plotCrossCorrs(self, varargin)
            %     clusIds [getActiveClusters(data)]
            %
            % AE 2012-07-02
            
            params.clusIds = getClusterIds(self);
            params = parseVarArgs(params, varargin{:});
            
            [corrs time] = getCrossCorrs(self, 'clusIds', params.clusIds);
            N = numel(params.clusIds);
            hdl = zeros(N);
            for i = 1 : N
                color = getClusColor(self, params.clusIds(i));
                for j = 1 : N
                    hdl(i, j) = axes('Position', [i-1 N-j 1 1] / N);
                    bb = bar(time, corrs{i, j}, 1, 'FaceColor', 'k', 'LineStyle', 'none');
                    if i == j
                        set(bb, 'FaceColor', color)
                    else
                        hold on
                        plot(0, 0, '*', 'Color', color)
                    end
                    axis tight off
                    axis([xlim * 1.05, ylim * 1.1])
                end
            end
            
            if nargout
                varargout{1} = hdl;
            end
        end
        
        function [fp fn snr frac] = getStats(self,varargin)
            % Get statistics for a given cluster or all
            %
            % [fp fn snr frac] = getStats(data,varargin)
            %      clusIds [getActiveClusters(data)]
            %
            % JC 2009-09-27
            
            params.clusIds = getClusterIds(self);
            params = parseVarArgs(params,varargin{:});
            
            cm = getContamination(self);
            fp = (sum(cm, 2) - diag(cm))';
            fn = diag(cm)';

            fp = fp(params.clusIds);
            fn = fn(params.clusIds);
            
            for i = 1:length(params.clusIds)
                index = getSpikesByClusIds(self, params.clusIds(i));
                ids{i} = index;
                if isempty(ids{i}), snr(i) = NaN; continue; end
                if isfield(self.Waveforms.meta, 'subset')
                    [~, index] = ismember(index, self.Waveforms.meta.subset);
                    index = index(index > 0);
                end
                wf = self.Waveforms.data{1}(:, index);
                snr(i) = range(mean(wf, 2)) / mean(std(wf, [], 2)) / 2;
            end
            
            frac = cellfun(@length, ids) / sum(cellfun(@length, ids));
        end
        
        function plotLDAs(self, varargin)
            % Plot LDA
            %
            % plotLDAs(data, varargin)
            %     clusIds [getActiveClusters(data)]
            %
            % AE 2012-07-09
            % JC 2009-09-29
            
            params.clusIds = getClusterIds(self);
            params = parseVarArgs(params,varargin{:});
            
            assert(length(params.clusIds) > 1, 'Can only do for multiple clusters');
            
            clf
            for i = 1:length(params.clusIds)
                for j = i+1:length(params.clusIds)
                    subplot(length(params.clusIds)-1,length(params.clusIds)-1,i + (length(params.clusIds)-1)*(j-2));
                    
                    % get spike for that cluster
                    ids1 = getSpikesByClusIds(self, abs(params.clusIds(i)));
                    ids2 = getSpikesByClusIds(self, abs(params.clusIds(j)));
                    
                    % project using linear discriminant analysis
                    dat1 = self.Features.data(ids1,:);
                    dat2 = self.Features.data(ids2,:);
                    w1 = (cov(dat1)+cov(dat2))^-1 * (mean(dat1,1)'-mean(dat2,1)');
                    
                    p1 = dat1*w1;
                    p2 = dat2*w1;
                    
                    bins = linspace(min([p1;p2]),max([p1;p2]),50);
                    n1 = hist(p1,bins);
                    n2 = hist(p2,bins);
                    
                    bar(bins,[n1;n2]',1,'stacked','linestyle','none');
                end
            end
        end
            
        
        function color = getClusColor(self,clusId)
            % Get the color for a cluster
            color = jet(length(self.ClusterAssignment.data));
            color = color(clusId,:);
        end
    end
    
    methods (Static)
        
        function lim = updateLimits(lim, X)
            range = prctile(X, 100 * normcdf([-1 1])); % 1 sigma percentiles
            range = range + [-2.5 2.5] * diff(range) / 2;  % +/- 3 sigma
            lim = [min(lim(1), range(1)), max(lim(2), range(2))];
        end
        
    end
end
