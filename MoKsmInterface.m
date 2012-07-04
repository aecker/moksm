classdef MoKsmInterface < SpikeSortingHelper & ClusteringHelper & MoKsm

	properties
	end

	methods
		function self = MoKsmInterface(varargin)
        % Creates a ClusteringHelper class.  Doesn't need to initialize
        % anything yet
        
            self = self@SpikeSortingHelper(varargin{:});
            self = self@ClusteringHelper();
            first = find(cellfun(@ischar, varargin), 1, 'first');
            self = self@MoKsm(varargin{first : end});
        end
        
        function self = updateInformation(self)
            % This method updates the ClusterAssignment and ContaminationMatrix
            ids = cluster(self);
            N = length(unique(ids));
            self.ClusterAssignment.data = cell(N,1);
            for i = 1:length(unique(ids))
                self.ClusterAssignment.data(i) = {find(ids == i)};
            end
            
            [pairwise, n] = overlap(self);
            self.ContaminationMatrix.data.pairwise = pairwise;
            self.ContaminationMatrix.data.n = n;
            
            if isempty(self.GroupingAssignment)
                self.GroupingAssignment(1).data = num2cell(1:length(self.ClusterAssignment.data));
                self.ClusterTags.data = cell(1,length(self.ClusterAssignment.data));
            end
        end
        
        function self = fit(self)
            % Fits the model
            self = fit@MoKsm(self, self.Features.data, self.SpikeTimes.data);
        end

        function self = delete(self, id)
            % Delete a cluster by ID.  Should call updateInformation afterwards.
            assert(all(cellfun(@length, self.GroupingAssignment.data(id)) == 1), ...
                'Can only merge ungrouped clusters');
        
            remove_id = cat(1,self.GroupingAssignment.data{id});
                           
            % remove unused components
            self.model.mu(:, :, remove_id) = [];
            self.model.C(:, :, remove_id) = [];
            self.model.priors(remove_id) = [];
            
            self = EStep(self);
            
            % Remove the pointer to the deleted cluster and decrement all
            % others that are greater than it.  Need to go from back to
            % front in case something has multiple clusters before it
            remove_id = sort(remove_id,'descend');
            for i = 1:length(remove_id)
                
                assert(~any(cellfun(@(x) any(x == remove_id(i)) && length(x) > 1, self.GroupingAssignment.data)), ...
                    'The deleted cluster was in a group.  Should not have happened');

                % Must delete the original pointer to the removed cluster
                % before shifting others or we delete two things
                idx = find(cellfun(@(x) all(x == remove_id(i)), self.GroupingAssignment.data));
                self.GroupingAssignment.data(idx) = [];
                self.ClusterTags.data(idx) = [];

                % Decrement all points higher than the removed element in
                % the list
                for j = 1:length(self.GroupingAssignment.data)
                    idx = find(self.GroupingAssignment.data{j} > remove_id(i));
                    self.GroupingAssignment.data{j}(idx) = self.GroupingAssignment.data{j}(idx) - 1;
                end
            end
            
            % Reclassify all the points into the remaining clusters
            self = updateInformation(self);
        end
        
        function self = split(self, id)
            % Splits a cluster by ID.  Should call updateInformation
            % afterwards.
        
            assert(length(id) == 1, 'Only split one cluster at a time');
            group = length(self.GroupingAssignment.data{id}) > 1;
            
            switch(group)
                case true % For groups simply split apart into raw clusters
                    clusterIds = self.GroupingAssignment.data{id};
                    self.GroupingAssignment.data(id) = [];
                    self.ClusterTags.data(id) = [];
                    
                    self.GroupingAssignment.data(end+1:end+length(clusterIds)) = num2cell(clusterIds);
                    self.ClusterTags.data(end+1:end+length(clusterIds)) = cell(1,length(clusterIds));
                    % No need to rerun updateInformation as cluster
                    % assignments unchanged
                case false
                    error('Not implemented yet to split a single cluster');
                    self = updateInformation(self);
            end

        end

        function self = merge(self, ids)
        % Splits a cluster by ID.  Should call updateInformation
        % afterwards.
        
            assert(all(cellfun(@length, self.GroupingAssignment.data(ids)) == 1), ...
                'Can only merge ungrouped clusters');
        
            ids = cat(1,self.GroupingAssignment.data{ids});
            dest_id = min(ids);
            remove_id = setdiff(ids,dest_id);

            % remove unused components
            self.model.mu(:, :, remove_id) = [];
            self.model.C(:, :, remove_id) = [];
            self.model.priors(remove_id) = [];

            % update model
            self = EStep(self);
            self = MStep(self);
            
            % Remove the pointer to the deleted cluster and decrement all
            % others that are greater than it.  Need to go from back to
            % front in case something has multiple clusters before it
            remove_id = sort(remove_id,'descend');
            for i = 1:length(remove_id)
                
                assert(~any(cellfun(@(x) any(x == remove_id(i)) && length(x) > 1, self.GroupingAssignment.data)), ...
                    'The deleted cluster was in a group.  Should not have happened');

                % Must delete the original pointer to the removed cluster
                % before shifting others or we delete two things
                idx = find(cellfun(@(x) all(x == remove_id(i)), self.GroupingAssignment.data));
                self.GroupingAssignment.data(idx) = [];
                self.ClusterTags.data(idx) = [];

                % Decrement all points higher than the removed element in
                % the list
                for j = 1:length(self.GroupingAssignment.data)
                    idx = find(self.GroupingAssignment.data{j} > remove_id(i));
                    self.GroupingAssignment.data{j}(idx) = self.GroupingAssignment.data{j}(idx) - 1;
                end
            end
            
            % Reclassify all the points into the remaining clusters
            self = updateInformation(self);
        end
        
        function self = group(self, ids)
            % TODO: Support grouping when a group is included
            
            finalGroup = cat(2, self.GroupingAssignment.data{ids});
            
            % Delete previous groups and tags
            self.GroupingAssignment.data(ids) = [];
            self.ClusterTags.data(ids) = [];
            
            % Create new one
            self.GroupingAssignment.data(end+1) = {finalGroup};
            self.ClusterTags.data(end+1) = {[]};
        end
        
        function self = refit(self)
        % Refit the complete data set again
            self = EM(self);
        end
        
        function self = compress(self)
        % Remove any information that can be recomputed and doesn't
        % need to be stored on disk
        
            % Remove the Waveforms, SpikeTimes, Features, tt
            self = compress@SpikeSortingHelper(self);
            
            % Remove the intermediate data used for clustering
            self.Y = [];
            self.t = [];
        end
        
        function self = uncompress(self)
        % Recreate any information that compress strips out
            self = uncompress@SpikeSortingHelper(self);
            
            % Restore the MoKsm stuff
            self.Y = self.Features.data';
            self.t = self.SpikeTimes.data;
        end
        
    end

end
