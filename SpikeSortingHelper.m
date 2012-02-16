classdef SpikeSortingHelper

	properties
        electrode = [];
        tt = [];
        data = [];
	end

	methods
		function self = SpikeSortingHelper(electrode)
            self.electrode = electrode;
            self = loadTT(self);
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

            self.data.Waveforms = struct('data',{wf},'meta',struct('units', 'muV'));
        end

        % Get the times from the tt file
        function self = getTimes(self)
            self.data.SpikeTimes = struct('data', self.tt.t, 'meta', struct);
        end

        % Extract features from spike waveforms
        function self = getFeatures(self,feature)
            dat = cat(1,self.data.Waveforms.data{:});
            if strcmp(feature, 'Points') == 1
                X = dat([25 15 10],:)';
            else
                error('Unsupported feature');
            end
            self.data.Features = struct('data',X,'meta',struct('Feature',feature));
        end
    end
end
