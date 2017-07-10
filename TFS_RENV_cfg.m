% LPF_TFS_cfg 
%
% Configuration file for experiment to test VCV identification
% with various forms of TFS and RENV processed speech.  With and without noise (contin/interrupted).
%
% afc_main('lpf_tfs','subj',param1,param2,param3,param4,param5,param6,param7,param8,param9,param10,param11,param12,param13,param14, param15,comment)
%
% param1 = 'booth1','booth2','booth3','booth4','test'
% param2 = attenuator setting (incorporate HB gain)
% param3 = processing type
% 'INTACT','EEQ1','EEQ4'
% param4 = test ear 'Left','Right' 'Both'
% param5 = test mode 'Train', 'Test', 'Train(3)/Test'
% param6 = number of 64-trial runs
% param7 = source level in dB SPL
% param8 = noise type '30c', 'SNRc', 'SNRi', 'SNRi-SAM', 'VC-1', 'VC-2', 'VC-4'
% param9 = noise SNR in dB
% param10 = name of hearing loss file
% param11 = 0 for no NAL-RP and 1 for NAL-RP
% param12 = 'no' for no hearing loss simulation and 'yes' for simulation.

close all;

% Extract parameter 9 to set the number of sets to run
for k = 1:12,
	evalstr = sprintf('param%i = varargin{%i};',k,k);
	eval(evalstr);
end

def=struct(...
'expname','TFS_RENV',		...		% name of experiment   
'headphone','HD580',		...		% Headphone type (ER2 or HD580)   
'intervalnum',16,			...		% number of intervals (one per consonant
'measurementProcedure','constantStimuli',	...	% 'constantStimuli' or 'transformedUpDown' (default)
'expvar',[0],				...		% Everything presented at one level
'expvarunit','N/A',			...		% units of the tracking variable
'expvarnum', 64,		...		% number of presentations of variable, same index as in expvar 
'practice',0,				...		% enables practice presentations
'expvarord',0,				...		% order of presentation 0 = random, others not implemented yet
'exppar1',[1:param6],			...		% vector containing experimental parameters for which the exp is performed
'exppar1unit','N/A',			...		% units of experimental parameter
'exppar2',[(randperm(4)+1)],          ...
                            ...     % In this experiment, exppar2 is:
                            ...     % 1 = 30
                            ...     % 2 = SNRc
                            ...     % 3 = SNRi-SQW5
                            ...     % 4 = SNRi-SQW10
                            ...     % 5 = SNRi-SQW20
                            ...     % 6 = randomize SNRc + x
'exppar2unit','N/A',			...		% units of experimental parameter
'exppar3', [-12 -18], ...
'exppar3unit', 'N/A', ...
'repeatnum',1,				...		% number of repeatitions of the experiment
'parrand',[0 0 1],				...		% toggles random presentation of the elements in "exppar" on (1), off(0)
'mouse',1,				...		% enables mouse control (1), or disables mouse control (0)  
'feedback',0,				...		% visual feedback after response: 0 = no feedback, 1 = correct/false/measurement phase
'samplerate',44100,			...		% sampling rate in Hz
'intervallen',1,			...		% length of each signal-presentation interval in samples (might be overloaded in 'expname_set')
'pauselen',1,			...		% length of pauses between signal-presentation intervals in samples (might be overloaded in 'expname_set')
'presiglen',44100,			...		% length of signal leading the first presentation interval in samples (might be overloaded in 'expname_set')
'postsiglen',1,			...		% length of signal following the last presentation interval in samples (might be overloaded in 'expname_set')
'result_path','',			...		% where to save results
'control_path','',			...		% where to save control files
'markinterval',0,			...		% toggles visual interval marking on (1), off(0)
'messages','autoSelect',			...		% message configuration file
'savefcn','TFS_RENV',			...		% function which writes results to disk
'interleaved',0,			...		% toggles block interleaving on (1), off (0)
'interleavenum',3,		...		% number of interleaved runs
'debug',0,					...		% set 1 for debugging (displays all changible variables during measurement)
'dither',0,					...		% 1 = enable +- 0.5 LSB uniformly distributed dither, 0 = disable dither
'bits',32,				...		% output bit depth: 8 or 16
'sethandle',@set,            ...    % JGD -- Store a handle to the set function
'msghandle',0 ,              ...    % JGD -- Store a handle to message text (assigned in *_set.m)
'checkOutputClip', 1, ...
'backgroundsig',1			...		% allows a backgroundsignal during output: 0 = no bgs, 1 = bgs is added to the other signals, 2 = bgs and the other signals are multiplied
);

def.feedback = 0;
def.markcorrect = 0;

switch lower(param8)
    case '30c'
        def.exppar2 = 1;
    case 'snrc'
        def.exppar2 = 2;
    case 'snri-sqw5'
        def.exppar2 = 3;
    case 'snri-sqw10'
        def.exppar2 = 4;
    case 'snri-sqw20'
        def.exppar2 = 5;
end


%switch param9
 %   case -6
  %      def.exppar3 = 1;
   % case -12
    %    def.exppar3 = 2;
%end

def.consonants = {'B','D','F','G','J','K','L','M','N','P','R','S','SH','T','V','Z'};

% Make sure the hearing loss functions are in the path
currentdir = pwd;
cd ..
hl_path = lower([pwd,'\HL_functions']);
cd(currentdir);
p = lower(path);
if isempty(strfind(p,hl_path)), path(hl_path,path); end


def.control_filename = ['control_',def.expname,'_',vpname,'_',varargin{end},'.dat'];
if exist(def.control_filename,'file')
    def.new_expt = 0;
else
    def.new_expt = 1;
end

% eof
