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
            params.nBins = 40;
            params = parseVarArgs(params,varargin{:});
            
            for i = 1:length(params.clusIds)
                ids1 = getSpikesByClusIds(self,params.clusIds(i));
                t1 = getSpikeTimes(self, ids1);
                
                for j = i:length(params.clusIds)
                    ids2 = getSpikesByClusIds(self,params.clusIds(j));
                    t2 = getSpikeTimes(self,ids2);
                    [corrs{i,j} time] = CrossCorr(t1,t2,params.binSize,params.nBins);
                    if i == j
                        corrs{i,j}(params.nBins / 2 + 1) = 0;
                    end
                    corrs{j,i} = flipud(corrs{i,j});
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
        
        function plotProjections(self,varargin)
            % Plot the prjoections
            %
            % plotProjections(data, varargin)
            %     clusIds [getActiveClusters(data)]
            %     maxPoints [20000]
            %     scatter [0] - use scatter plot so colors more mixed, slower
            %
            % JC 2009-09-27
            
            params.maxPoints = 5000;
            params.clusIds = getClusterIds(self);
            params.pointSize = 1;
            params.scatter = 0;
            params = parseVarArgs(params,varargin{:});
            
            c = getClusColor(self, params.clusIds);
            
            cla
            hold on
            set(gca,'Color',[.1 .1 .1]);
            
            X = self.Features.data;
            
            if params.scatter
                
                for i = 1:length(params.clusIds)
                    ids{i} = getSpikesByClusIds(self, params.clusIds(i));
                    colors{i} = repmat(c(i,:),length(ids{i}),1);
                end
                
                ids = cat(2,ids{:});
                colors = cat(1,colors{:});
                
                r = randperm(length(ids));
                r = r(1:min(end,params.maxPoints));
                
                ids = ids(r);
                colors = colors(r,:);
                
                scatter(X(ids,1),X(ids,2),ones(size(ids)) * params.pointSize,colors,'filled');
            else
                show = getSpikesByClusIds(self, params.clusIds);
                r = randperm(length(show));
                show = show(r(1:min(end,params.maxPoints)));
                
                if size(X,2) <= 3
                    for i = 1:length(params.clusIds)
                        ids = getSpikesByClusIds(self, params.clusIds(i));
                        ids = ids(ismember(ids,show));
                        plot(X(ids,1),X(ids,2),'.','color',c(i,:),'MarkerSize',params.pointSize);
                    end
                else
                    spacing = 400;
                    for i = 1:length(params.clusIds)
                        ids = getSpikesByClusIds(self, params.clusIds(i));
                        ids = ids(ismember(ids,show));
                        plot(X(ids,1),X(ids,4),'.','color',c(i,:),'MarkerSize',params.pointSize);
                        plot(X(ids,1),X(ids,7)+spacing,'.','color',c(i,:),'MarkerSize',params.pointSize);
                        plot(X(ids,1),X(ids,10)+spacing*2,'.','color',c(i,:),'MarkerSize',params.pointSize);
                        plot(X(ids,4)+spacing,X(ids,7),'.','color',c(i,:),'MarkerSize',params.pointSize);
                        plot(X(ids,4)+spacing,X(ids,10)+spacing,'.','color',c(i,:),'MarkerSize',params.pointSize);
                        plot(X(ids,7)+spacing,X(ids,10)+spacing*2,'.','color',c(i,:),'MarkerSize',params.pointSize);
                    end
                end
            end
        end
        
        function plotContaminations(self,varargin)
            % Plot the contamination matrix
            %
            % plotWaveforms(data, varargin)
            %     clusIds [getActiveClusters(data)]
            %     maxPoints [20]
            %
            % JC 2009-09-27
            
            params.clusIds = getClusterIds(self);
            params = parseVarArgs(params,varargin{:});
            
            cla
            imagesc(self.ContaminationMatrix.data(params.clusIds,params.clusIds));
            xlabel('Source'); ylabel('Classified');
            
            set(gca,'XTick',1:length(params.clusIds));
            set(gca,'XTickLabel',params.clusIds);
            set(gca,'YTick',1:length(params.clusIds));
            set(gca,'YTickLabel',params.clusIds);
            
            set(gca,'CLim',[0 1]);
        end
        
        function varargout = plotWaveforms(self, varargin)
            % Plot the waveforms
            %
            % plotWaveforms(data, varargin)
            %     clusIds [getActiveClusters(data)]
            %     maxPoints [20]
            %
            % JC 2009-09-27
            
            params.maxPoints = 20;
            params.clusIds = getClusterIds(self);
            params.figure = [];
            params = parseVarArgs(params,varargin{:});
            
            if isempty(params.figure)
                figure
            else
                figure(params.figure)
            end
            
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
                    yli = ylim;
                    yl = [min(yli(1), yl(1)), max(yli(2), yl(2))];
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
            
            for i = 1:size(X,2)
                subplot(size(X,2),1,i);
                hold on
            end
            
            for i = 1:length(params.clusIds)
                color = getClusColor(self, params.clusIds(i));
                ids = getSpikesByClusIds(self,params.clusIds(i));
                for j = 1:size(X,2)
                    subplot(size(X,2),1,j);
                    plot(self.SpikeTimes.data(ids),X(ids,j),'.','color',color);
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
            params.figure = [];
            params = parseVarArgs(params, varargin{:});
            
            if isempty(params.figure)
                figure
            else
                figure(params.figure)
            end
            
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
end
