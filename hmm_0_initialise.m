function config = hmm_0_initialise
% hmm_0_initialise - Set paths and run sanity check for the HMM burst tutorial
%
% The following toolboxes will be added to the path based on diretories defined in the top part of this script
% - HMM-MAR
% - SPM
% - distributionPlot
%
% Checks will be run to make sure that
% - matlab version is greater than 2018a (scripts may run on older versions but haven't been tested')
% - the Signal Processing and Wavelet toolboxes are present
% - HMM-MAR, spm and distributionPlot can all be found on the path


%% User specified paths
% Please set the following paths based on your system set-up

config = [];
% Set the path to this download on the following line
config.basepath = '';
% Set the path to the HMM-MAR toolbox on the following line
config.hmmmarpath = '';
% Set the path to distributionPlot here
config.distPlotpath = '';

config.figpath = fullfile(config.basepath,'figures');

%% Sanity check paths are specified and exist

flds = fields(config);
for ii = 1:numel(flds)
    if isempty( config.(flds{ii}) )
        warning( [ flds{ii} ' path not specified']);
    end
    if exist( config.(flds{ii}),'dir' ) ~= 7 % not sure why this doesn't return bool?
        warning( [ flds{ii} ' path not found']);
    end
end


%% Add toolboxes to path

% Add HMM-MAR to path
addpath(genpath(config.hmmmarpath))
addpath(genpath(config.distPlotpath))


%% Sanity check all toolboxes are present
allclear = true;

% Is matlab recent enough?
if verLessThan('matlab','9.4')
    warning('Matlab is older than R2018a - the tutorial might not run...');
    allclear = false;
end

% Do we have the signal processing toolbox and EMD functions?
if ~license( 'test', 'Signal_Toolbox' )
    warning('Signal Processing Toolbox is missing');
    warning('This is available as a MatLab toolbox or add-on');
    allclear = false;
end

if isempty( which('cwt') )
    warning('Wavelet toolbox is missing!');
    warning('This is available as a MatLab toolbox or add-on');
    allclear = false;
end

if isempty( which('hmmmar') )
    warning('hmmmar toolbox is missing!');
    warning('This can be downloaded from: https://github.com/OHBA-analysis/HMM-MAR');
    allclear = false;
end

if isempty( which('ROInets.remove_source_leakage') )
    warning('ROInets is missing!');
    warning('This can be downloaded from: ');
    allclear = false;
end


% Final message
if allclear == true
    disp('Burst HMM tutorial checks passed');
else
    disp('Sanity check failed - please check warnings above');
end

end
