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
        % Creates a ClusteringHelper class.  Doesn't need to initialize
        % anything yet
%            self = updateInformation(self);
            %self.GroupingAssignment = struct('data',{},'meta',struct);
            %self.ClusterAssignment = struct('data',{},'meta',struct);
        end
        
        function s = saveStructure(self)
            f = properties(self);
            for i = 1:length(f)
                s.(f{i}) = self.(f{i});
            end
        end
        
        function self = singleUnit(self, clusterIds)
            % Toggles the SingleUnit flag for the selected ids
            for i = 1:length(clusterIds)
                tags = self.ClusterTags.data(clusterIds(i));
                su_flags = strcmp(tags,'SingleUnit');
                if ~any(su_flags) % If flag not found, add it
                    self.ClusterTags.data{clusterIds(i)}{end+1} = 'SingleUnit';
                else              % Delete any found flags
                    self.ClusterTags.data{clusterIds(i)}(su_flags) = [];
                end
            end
        end
        
        function cm = getContamination(self,ids)
            % Returns the full or a subset of the contamination matrix
            if nargin < 2
                cm = self.ContaminationMatrix;
            else
                cm = self.ContaminationMatrix(ids,ids);
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
            m = fix(sqrt(N));
            n = ceil(N / m);
            
            % select subset of points to show
            K = numel(params.clusIds);
            show = cell(1, K);
            for k = 1 : K
                ids = getSpikesByClusIds(self, params.clusIds(k));
                r = randperm(numel(ids));
                perCluster = round(params.maxPoints * self.model.priors(params.clusIds(k)));
                show{k} = ids(r(1 : min(end, perCluster)));
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
                    axHdl(ij) = axes('Position', pos, 'Color', 0.35 * ones(1, 3));
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
            hdl = imagesc(self.ContaminationMatrix.data(params.clusIds,params.clusIds));
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
                ids = getSpikesByClusIds(self, abs(params.clusIds(i)));
                if isempty(ids), continue; end
                r = randperm(length(ids));
                ids = sort(ids(r(1 : min(end, params.maxPoints))));
                for j = 1 : chans
                    hdl(i, j) = axes('Position', [(i - 1) / k, (j - 1) / chans, 1 / k, 1 / chans]);
                    plot(self.Waveforms.data{j}(:, ids), 'Color', color)
                    axis tight off
                    m = median(self.Waveforms.data{j}(:, ids), 2);
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
                    ylj = prctile(Xij, 100 * normcdf([-1 1])); % 1 sigma percentiles
                    ylj = ylj + [-3 3] * diff(ylj) / 2;        % plot +/- 3 sigma
                    yl(j, :) = [min(yl(j, 1), ylj(1)), max(yl(j, 2), ylj(2))];
                    ylim(yl(j, :))
                end
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
                    axis([xlim * 1.1, ylim * 1.2])
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
            
            %[cm fp fn] = getContamination(self.Model);
            cm = zeros(length(params.clusIds));
            fp = zeros(1,length(params.clusIds));
            fn = zeros(1,length(params.clusIds));
            
            fp = fp(params.clusIds);
            fn = fn(params.clusIds);
            cm = cm(params.clusIds,params.clusIds);
            % fp = ones(size(params.clusIds));
            % fn = ones(size(params.clusIds));
            % cm = ones(length(params.clusIds));
            
            for i = 1:length(params.clusIds)
                ids{i} = getSpikesByClusIds(self, params.clusIds(i));
                
                if isempty(ids{i}), snr(i) = NaN; continue; end
                
                wf = self.Waveforms.data{1}(:,ids{i});
                snr(i) = range(mean(wf,2)) / mean(std(wf,[],2)) / 2;
            end
            
            frac = cellfun(@length, ids) / sum(cellfun(@length, ids));
        end
        
        function plotLDAs(self, varargin)
            % Plot the waveforms
            %
            % plotLDAs(data, varargin)
            %     clusIds [getActiveClusters(data)]
            %
            % JC 2009-09-29
            
            params.clusIds = getClusterIds(self);
            params = parseVarArgs(params,varargin{:});
            
            assert(length(params.clusIds) > 1, 'Can only do for multiple clusters');
            
            clf
            wf = cat(1,self.Waveforms.data{:});
            
            for i = 1:length(params.clusIds)
                for j = i+1:length(params.clusIds)
                    subplot(length(params.clusIds)-1,length(params.clusIds)-1,i + (length(params.clusIds)-1)*(j-2));
                    
                    % get spike for that cluster
                    ids1 = getSpikesByClusIds(self, abs(params.clusIds(i)));
                    ids2 = getSpikesByClusIds(self, abs(params.clusIds(j)));
                    
                    % project using linear discriminant analysis
                    dat1 = wf(:,ids1)';
                    dat2 = wf(:,ids2)';
                    w1 = (cov(dat1)+cov(dat2))^-1 * (mean(dat1,1)'-mean(dat2,1)');
                    
                    p1 = dat1*w1;
                    p2 = dat2*w1;
                    
                    bins = linspace(min([p1;p2]),max([p1;p2]),50);
                    n1 = hist(p1,bins);
                    n2 = hist(p2,bins);
                    
                    bar(bins,[n1;n2]',1,'stacked');
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
