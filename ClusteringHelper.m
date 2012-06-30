classdef ClusteringHelper

	properties
        ClusterAssignment = [];
        ClusterTags = [];
        ContaminationMatrix = [];
	end

	methods
		function self = ClusteringHelper()
        % Creates a ClusteringHelper class.  Doesn't need to initialize
        % anything yet
            self = updateInformation(self);
        end
        
        function self = singleUnit(self, clusterIds)
            % Toggles the SingleUnit flag for the selected ids
            for i = 1:length(clusterIds)
                tags = data.ClusterTags.data{clusterIds(i)};
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
        
        function ids = getClusterIds(self)
            % Returns the valid cluster ids
            ids = 1:length(self.ClusterAssignment);
            containGroups = find(cellfun(@(x) any(x < 0), self.ClusterAssignment));
            
            toRemove = [];
            for i = 1:length(containGroups)
                % Find the negative numbers (groups) in a cluster
                % assignment
                groups = self.ClusterAssignment{containGroups(i)};
                groups = -groups(groups < 0);
                toRemove = [toRemove groups];
            end
            
            % Remove any cluster ids contained by other cluster ids
            ids(unique(toRemove)) = [];
        end
                        
        function spikeIds = getSpikesByClusIds(self,clusterIds)
            spikeIds = cat(1,self.ClusterAssignment{clusterIds});
            
            % Replace any groups with their spike ids
            while any(spikeIds < 0)
                groups = find(spikeIds < 0);
                groupIds = unique(spikeIds(groups));
                
                % Remove the negative numbers and replace with their sets
                % of spike ids
                spikeIds(groups) = [];                   %#ok<*AGROW>
                spikeIds = [spikeIds, cat(1,self.ClusterAssignment{groupIds})];
            end
            
            spikeIds = sort(unique(spikeIds));
        end
        
        % TODO: Probably there should be a private method for actually
        % adding and deleting clusters from updateInformation so the
        % tags follow the clusters around properly
    
        function [corrs time] = getCrossCorrs(self, varargin)
            params.clusIds = getClusterIds(self);
            params.binSize = 1;
            params.nBins = 30;
            params = parseVarArgs(params,varargin{:});
            
            for i = 1:length(params.clusIds)
                ids1 = getSpikesByClusIds(self,params.clusIds(i));
                t1 = getSpikeTimes(self, ids1);
                
                for j = i:length(params.clusIds)
                    ids2 = getSpikesByClusIds(self,params.clusIds(j));
                    t2 = getSpikeTimes(self,ids2);
                    [corrs{i,j} time] = CrossCorr(t1,t2,params.binSize,params.nBins);
                end
            end
        end
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

        % Refit the complete data set again
        self = refit(self);
        
        % Remove any information that can be recomputed and doesn't
        % need to be stored on disk
        self = compress(self);
        
        % Recreate any information that compress strips out
        self = uncompress(self);
    end

end
