% exampleConstantStimuli_user - stimulus generation function of experiment 'exampleConstantStimuli' -
%
% This function is called by afc_main when starting
% the experiment 'exampleConstantStimuli'. It generates the stimuli which
% are presented during the experiment.
% The stimuli must be elements of the structure 'work' as follows:
%
% work.signal = def.intervallen by 2 times def.intervalnum matrix.
%               The first two columns must contain the test signal
%               (column 1 = left, column 2 = right) ...
% 
% work.presig = def.presiglen by 2 matrix.
%               The pre-signal. 
%               (column 1 = left, column 2 = right).
%
% work.postsig = def.postsiglen by 2 matrix.
%                The post-signal. 
%               ( column 1 = left, column 2 = right).
%
% work.pausesig = def.pausesiglen by 2 matrix.
%                 The pause-signal. 
%                 (column 1 = left, column 2 = right).
% 
% To design an own experiment, e.g., 'myexperiment',
% make changes in this file and save it as 'myexperiment_user.m'.
% Ensure that the referenced elements of structure 'work' are existing.
%
% See also help exampleConstantStimuli_cfg, exampleConstantStimuli_set, afc_main

function TFS_RENV_user

global def
global work
global set
global speaker_choices;


% Assign the correct interval for the current stimulus
def.ranpos = 17-set.stim_consonant(set.stim_indx);

% Load the current stimulus
[s,Fs] = audioread([work.stim_directory, set.stim_list{set.stim_indx}]);
s = s(:,1);
s = resample(s, def.samplerate, Fs);
s = fftfilt(set.lpf,s);

% Keep track of the stimulus that was presented to the user
work.test_stimuli{set.stim_indx} = set.stim_list{set.stim_indx};

if strcmpi(work.userpar12,'yes')
    s = [zeros(round(def.samplerate/2),1); s; zeros(round(def.samplerate/2),1)];
end

% Scale to desired level
base_level = 10*log10(mean(s.^2));
s = s*(10^((set.desiredLevel-base_level)/20));

% extract noise clip
% Extract a block of noise same length as stimulus
start_ind = 5000 + round((length(set.n_base)-length(s)-10000)*rand(1));
n = set.n_base(start_ind + [1:length(s)]);

baseline_start = max(1,round((length(set.contra_masker)-length(s)-10)*rand(1)));
baseline_end = baseline_start + length(s) - 1;
noise = set.baseline_noise(baseline_start:baseline_end);

chosenNoise = lower(deblank(work.exppar3)); 
if chosenNoise == 1
    % 30 dB continuous noise.
    n = n * (10^((30)/20));	% Scale noise (originally at 0 dB SPL) to desired level
elseif chosenNoise == 2
    % Continuous noise.
    n = n * (10^((set.noiseLevel)/20));	% Scale noise (originally at 0 dB SPL) to desired level
elseif chosenNoise == 3
    % Interrupted noise.
    % Create 4 Hz square-wave window that alternates between 0 and 1.
    win = 0.5 + 0.5*square(2*pi*4*[0:2*size(n,1)-1]/def.samplerate);
    win = win(max(1,round(size(n,1)*rand(1)))+[0:size(n,1)-1]);
    % Level in gaps is always 30 dB.  Determine level at peaks
    % to yield desired overall level
    peak_level = 10*log10(2*(10^(set.noiseLevel/10)) - 10^(30/10));
    win = 30 + (peak_level-30)*win(:);
    % Convert to linear and scale the noise
    n = n .* (10.^(win/20));
elseif chosenNoise == 4
    % Interrupted noise.
    % Create 10 Hz square-wave window that alternates between 0 and 1.
    win = 0.5 + 0.5*square(2*pi*10*[0:2*size(n,1)-1]/def.samplerate);
    win = win(max(1,round(size(n,1)*rand(1)))+[0:size(n,1)-1]);
    % Level in gaps is always 30 dB.  Determine level at peaks
    % to yield desired overall level
    peak_level = 10*log10(2*(10^(set.noiseLevel/10)) - 10^(30/10));
    win = 30 + (peak_level-30)*win(:);
    % Convert to linear and scale the noise
    n = n .* (10.^(win/20));
elseif chosenNoise == 5
    % Interrupted noise.
    % Create 32 Hz square-wave window that alternates between 0 and 1.
    win = 0.5 + 0.5*square(2*pi*32*[0:2*size(n,1)-1]/def.samplerate);
    win = win(max(1,round(size(n,1)*rand(1)))+[0:size(n,1)-1]);
    % Level in gaps is always 30 dB.  Determine level at peaks
    % to yield desired overall level
    peak_level = 10*log10(2*(10^(set.noiseLevel/10)) - 10^(30/10));
    win = 30 + (peak_level-30)*win(:);
    % Convert to linear and scale the noise
    n = n .* (10.^(win/20));
elseif chosenNoise == 6
    % Interrupted noise.
    % Create 256 Hz square-wave window that alternates between 0 and 1.
    win = 0.5 + 0.5*square(2*pi*256*[0:2*size(n,1)-1]/def.samplerate);
    win = win(max(1,round(size(n,1)*rand(1)))+[0:size(n,1)-1]);
    % Level in gaps is always 30 dB.  Determine level at peaks
    % to yield desired overall level
    peak_level = 10*log10(2*(10^(set.noiseLevel/10)) - 10^(30/10));
    win = 30 + (peak_level-30)*win(:);
    % Convert to linear and scale the noise
    n = n .* (10.^(win/20));
end

% Add noise and signal
s = s+n;

% Increment stim index
set.stim_indx = set.stim_indx + 1;


% Process the simulus
s1 = s;
s = resample(s1,24000,def.samplerate);
switch lower(deblank(work.userpar3(1:3)))
	case {'eeq'},
        % Called when running an experiment
        i = find((double(work.userpar3)>=48)&(double(work.userpar3)<=57),1,'first');
        eeq_bands = eval(work.userpar3(i:end));
        Tfast = 5;
        Tslow = 200;
        gain_bounds_dB = [0,20];
        [s] = LD_EnergyEqualize_NB(s,24000,eeq_bands,Tfast,Tslow,true,gain_bounds_dB);
end
s = resample(s,def.samplerate,24000);

% Convert to stereo signal using left/right presentation flags
s = s(:) * [set.left set.right];

% Apply NAL gain if desired
if (work.userpar11 == 1)
    s = [s;0*set.NALGainFilter];
    for k = 1:2,
        s(:,k) = fftfilt(set.NALGainFilter(:,k),s(:,k));
    end
end

% Simulate hearing loss
work.bgsig = 0*[s; set.headphone_comp*[0 0]];
set.ear_scale = [set.left set.right];
for k = 1:2,
    if set.MBE_active(k)==1,
    % If the MBE is required, filter signal and do the MBE processing
        s_filt = mg_HL_filtering(s(:,k), set.filter_data{k});

        % Use newer MG processing
        s(:,k) = mg_HL_processing_v4(s_filt,set.mg_data{k},set.filter_data{k});
    end
    % Generate a block of masking noise of appropriate length plus length of
    % headphone comp filter. Apply compensation filter and truncate to desire
    % length.
    if strcmpi(work.userpar12,'yes')
		% If the loss is being simulated, create non-zero noise
		HL_noise = masking_noise(set.fc,set.HL_TN(:,k),(size(s,1)+length(set.headphone_comp))/def.samplerate,def.samplerate); 
	else
		% Else, the masking noise is for a 0 dB loss
		HL_noise = 0*masking_noise(set.fc,0*set.HL_TN(:,k),(size(s,1)+length(set.headphone_comp))/def.samplerate,def.samplerate); 
	end
	HL_noise = fftfilt(set.headphone_comp,[HL_noise(:)]); % Now, apply the filter
    HL_noise = HL_noise([1:size(s,1)+length(set.headphone_comp)])';
	%HL_noise = [0.5-0.5*cos(2*pi*[0:999]'*pi/1000); ones(length(HL_noise)-2000,1); 0.5+0.5*cos(2*pi*[0:999]'*pi/1000)].*HL_noise(:);
	HL_noise = hann(HL_noise(:)',10,def.samplerate)';
	work.bgsig(:,k) = scale(HL_noise,-set.max_level);
end
work.bgsig = [work.bgsig;zeros(def.intervalnum*def.intervallen + (def.intervalnum-1)*def.pauselen + def.postsiglen,2)];


% scale for headphone presentation
s = [s; set.headphone_comp*[0 0]];
for k = 1:2,
	s(:,k) = scale(s(:,k),-set.max_level);
    s(:,k) = fftfilt(set.headphone_comp,s(:,k));
    s(:,k) = hann(s(:,k).', 10, def.samplerate).';
end

% Formulate outptus
def.presiglen = size(s,1);
work.signal = zeros(1,32);
work.presig = min(0.999,max(-0.999,s));
work.pausesig = [0 0];
work.postsig = [0 0];


% filename = ['stim',int2str(set.stim_num),'.wav'];
% wavwrite(s, def.samplerate, filename);
% set.stim_num = set.stim_num + 1;

	

% eof