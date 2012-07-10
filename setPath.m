function setPath
% Set Matlab path.

base = fileparts(mfilename('fullpath'));
addpath(base)
addpath(fullfile(base, 'lib'))
