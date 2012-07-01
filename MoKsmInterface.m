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
            
            self.ContaminationMatrix.data = eye(N);
        end
        
        function self = fit(self)
            % Fits the model
            self = fit@MoKsm(self, self.Features.data, self.SpikeTimes.data);
        end

        function self = delete(self, id)
            % Delete a cluster by ID.  Should call updateInformation afterwards.
        end
        
        function self = split(self, id)
            % Splits a cluster by ID.  Should call updateInformation
            % afterwards.
        
        end

        function self = merge(self, ids)
        % Splits a cluster by ID.  Should call updateInformation
        % afterwards.
        
        end
        
        function self = refit(self)
        % Refit the complete data set again
        
        end
        
        function self = compress(self)
        % Remove any information that can be recomputed and doesn't
        % need to be stored on disk
        
        end
        
        function self = uncompress(self)
        % Recreate any information that compress strips out
        
        end
        
    end

end
