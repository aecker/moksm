classdef SpikeSortingHelper

	properties
        Waveforms = [];
        SpikeTimes = [];
        Features = [];
        tt = [];
        dataSource = [];
	end

	methods
		function [self args] = SpikeSortingHelper(electrode,varargin)
            if isstruct(electrode) && isfield(electrode, 't')
                self.dataSource = struct('type','tt');
                self.tt = electrode;
                args = varargin;
            elseif isstruct(electrode) && count(detect.Electrodes(electrode)) > 0
                self.dataSource = struct('type','DataJoint', 'key', electrode);
                self = loadTT(self);
                args = varargin;
            elseif ismatrix(electrode) && nargin > 1 && any(size(electrode) == length(varargin{1}))
                warning('Construct fake spike structure.  Only use for debugging.');
                self.dataSource = struct('type','Raw');
                self.SpikeTimes = struct('data', varargin{1}, 'meta', struct);
                self.Features = struct('data', electrode, 'meta', struct);
                args = varargin(2:end);
                return; % Don't try and get waveforms and spike times
            else
                error('Could not construct data for the SpikeSortingHelper');
            end
            self = getWaveforms(self);
            self = getTimes(self);
        end        

        function self = loadTT(self)
            % Load the TT file
            assert(count(detect.Electrodes(self.electrode)) == 1, 'Only call this for one VC');

            de = fetch(detect.Electrodes(self.electrode), 'detect_electrode_file');
            fn = getLocalPath(de.detect_electrode_file);

            self.tt = ah_readTetData(fn);
        end

        % Get the waveforms (and scale them) from the tt file
        function self = getWaveforms(self)
            % Get and scale the waveforms

            if max(mean(self.tt.w{1},2)) > 1  % data originally was in raw values
                wf = cellfun(@(x) x / 2^23*317000, self.tt.w, 'UniformOutput',false);
            else % new data is in volts, convert to 
                wf = cellfun(@(x) x * 1e6, self.tt.w, 'UniformOutput', false);
            end

            self.Waveforms = struct('data',{wf},'meta',struct('units', 'muV'));
        end

        % Get the times from the tt file
        function self = getTimes(self)
            self.SpikeTimes = struct('data', self.tt.t, 'meta', struct);
        end

        % Extract features from spike waveforms
        function self = getFeatures(self,feature)
            if strcmp(feature, 'Points') == 1
                dat = cat(1,self.data.Waveforms.data{:});
                X = dat([25 15 10],:)';
            elseif strcmp(feature, 'PCA') == 1
                X = [];
                for i = 1:length(self.data.Waveforms.data)
                    [~,P] = princomp(self.data.Waveforms.data{i}');
                    X = [X P(:,1:3)];
                end
            else
                error('Unsupported feature');
            end
            self.Features = struct('data',X,'meta',struct('Feature',feature));
        end
    end
end
