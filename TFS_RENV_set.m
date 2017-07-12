% exampleConstantStimuli_set - setup function of experiment 'exampleConstantStimuli' -
%
% This function is called by afc_main when starting
% the experiment 'exampleConstantStimuli'. It defines elements
% of the structure 'set'. The elements of 'set' might be used 
% by the function 'exampleConstantStimuli_user.m'.
% 
% If an experiments can be started with different experimental 
% conditions, e.g, presentation at different sound preasure levels,
% one might switch between different condition dependent elements 
% of structure 'set' here.
%
% For most experiments it is also suitable to pregenerate some 
% stimuli here.
% 
% To design an own experiment, e.g., 'myexperiment',
% make changes in this file and save it as 'myexperiment_set.m'.
% Ensure that this function does exist, even if absolutely nothing 
% is done here.
%
% See also help exampleConstantStimuli_cfg, exampleConstantStimuli_user, afc_main

function TFS_RENV_set

global def
global work
global set

% Reconfigure window
try
	H = get(gcf,'Children');
    Ht = findobj(H,'Tag','afc_message');
	feval(def.sethandle,Ht,'Position',[0.1 0.55 0.8 0.4]);
	for k1 = 1:4,
		for k2 = 1:4,
			indx = (k1-1)*4 + k2;
            Ht = findobj(H,'Tag',sprintf('afc_button%i',17-indx));
			feval(def.sethandle,Ht,'Position',[0.1+(k2-1)*0.21 0.05+(4-k1)*0.12 0.17 0.09],'String',def.consonants{indx},'FontWeight','bold','FontSize',0.55);
		end
	end
	pause(0.1);
	drawnow;
	pause(0.1);
end

% If we are a new experiment of INTACT speech, then reconfigure the control
% file so that there is just one repetition of quiet (instead of 8)
switch lower(work.userpar5)
	case 'train'
		use_test = 0;
	case 'test'
		use_test = 1;
	case 'train(3)test'
		if (work.exppar1 <= 3)
			use_test = 0;
		else
			use_test = 1;
		end
end
% work.exppar1 = number of 64 trial runs
% if strcmpi(work.userpar3,'intact')&&(length(def.exppar2)==3)&&(work.exppar1==1)&&(work.exppar2==1)
% 	use_test = 1;
% end

% Commented since this experiment does not use test stimuli
% if (def.new_expt)&&strcmpi(work.userpar3,'intact')&&(length(def.exppar2)==3)
%     fid = fopen(def.control_filename,'rt');
%     s = textscan(fid,'%s','Delimiter','\n');
%     fclose(fid);
%     fid = fopen(def.control_filename,'wt');
%     fprintf(fid,'%s\n',s{1}{1:4});
%     fprintf(fid,'%s\n',s{1}{12:end});
%     fclose(fid);
%     def.new_expt = 0;
% 	use_test = 1;
% end


% Get the max infromation for he booth in use.
set.max_level = findmaxlevel(work.userpar1,work.userpar2,def.headphone);

% Load headphone compensation
if strcmp(def.headphone,'HD580'),
    load ..\HeadphoneCorr\HD580_correction;
    if def.samplerate ~= Fs_g,
        g = resample(g,def.samplerate,Fs_g);
    end
else
    g = [zeros(63,1);1;zeros(63,1)];
end
set.headphone_comp = g;
% Incorporate a bandlimiting lowpass filter into the headphone correction
% cutoff is upper end of 8kHz 3rd oct filter
fcut = 1000*(2^(9.5/3));
set.lpf = firpm(256,[0 fcut-300 fcut+300 def.samplerate/2]/(def.samplerate/2),[1 1 0 0]);
set.lpf = set.lpf(:);
%set.headphone_comp = conv(set.headphone_comp(:),lpf(:));
[G_headphone,f_headphone] = freqz(set.headphone_comp,1,8192,def.samplerate);

% Get hearing loss info
work.userpar10 = ['..\Losses\' work.userpar10];
hl_filename = work.userpar10;
load(hl_filename);
fc_HL = fc;
[m,n] = size(HL);
if n>2,
    HL = HL';
    [m,n] = size(HL);
end
if n==1, HL = HL*[1 1]; end
if n>2,
    error('Hearing loss is specified for a maximum of two ears!');
end
if m~=length(HL(:,1)),
    error('The vector fc must contain the frequencies at which the hearing loss is specified!');
end
if m<3, 
    error('Hearing loss must be specified at a minimum of three frequencies');
end
set.fc = fc_HL;
set.HL = HL;

% Now, calculate the HL masking noise parameters and a dummy set of masking
% noise.  First, determine how long noise must be to cover entire
% presentation
set.HL_noise_duration_samp = def.presiglen + def.intervalnum*def.intervallen + (def.intervalnum-1)*def.pauselen + def.postsiglen;
set.HL_noise_duration = set.HL_noise_duration_samp/def.samplerate;
set.HL_noise_counter = 0;
for k = 1:2,
	if 0
		% Use TN to simulate up to 30 dB of threshold shift.
		set.HL_TN(:,k) = min(30,set.HL(:,k));
	else
		% We will use to simulate thresholds of up to 40 dB SPL and MBE
		% otherwise.
		load thresh_normal
		set.thresh_interp = interpolate_vals(log(fc),thresh,log(set.fc));
		set.HL_TN(:,k) = max(0,min(40,set.HL(:,k) + set.thresh_interp) - set.thresh_interp);
	end
	% Remaining threshold will be simulated with MBE
	MBE_shift = set.HL(:,k) - set.HL_TN(:,k);
	if any(MBE_shift > 0) && strcmpi(work.userpar12,'yes')
		% We need to use MBE
		set.MBE_active(k) = 1;
		set.filter_data{k} = generate_HL_filters(def.samplerate);
		disp('MBE/MN in use!!!');

		set.HL_noise_delta(k) = 0;

		slope_MBE = (100-set.HL_TN(:,k))./(100-HL(:,k));
		HL_MBE = 100 - 100./slope_MBE;
		set.mg_data{k} = parameter_generate(set.filter_data{k}.f_HL,def.samplerate,fc_HL,HL_MBE,2);

	else
		set.MBE_active(k) = 0;
		set.HL_noise_delta(k) = 0;
	end
end

% Calculate the NAL gain filters
set.NALGainFilter = [nal_rp_filter(set.HL(:,1),set.fc,def.samplerate) nal_rp_filter(set.HL(:,2),set.fc,def.samplerate)];

% Create a list of stimuli from either training or testing stimuli
% Create a list of stimuli from either training or testing stimuli
set.stim_list = [];
set.stim_consonant = [];
if (use_test == 0),
    % If we are in runs 1-3, use training stimuli
    work.stim_directory = '..\TFS_Stimuli\consonants\train_stim\';
else
    % else, use testing stimuli
    work.stim_directory = '..\TFS_Stimuli\consonants\test_stim\';
end
d = dir([work.stim_directory,'*.wav']);
for k=1:length(d)
    % ldaquila change - don't add SPACED.wav
    if ~strcmp('SPACED.wav', d(k).name)
        set.stim_list{end+1,1} = d(k).name;
    end
    for k2 = 1:length(def.consonants),
        str = ['A',def.consonants{k2},'A'];
        if ~isempty(findstr(str,d(k).name))
            set.stim_consonant(end+1,1) = k2;
        end
    end
end
p = randperm(length(set.stim_list));
set.stim_list = set.stim_list(p);
% Prints the correct answers in order to a textfile "curranswers.txt"
fileID = fopen('curranswers.txt', 'w');
for i=1:length(set.stim_list)
    fprintf(fileID, '%s\n', set.stim_list{i});
end
fprintf(fileID, '\nEND OF RUN\N\N');
fclose(fileID);
set.stim_consonant = set.stim_consonant(p);
set.stim_indx = 1;

% Set up the left and right parameter
switch lower(deblank(work.userpar4)),
	case 'left',
        set.left = 1;
        set.right = 0;
    case 'right',
        set.left = 0;
        set.right = 1;
	case 'both',
        set.left = 1;
        set.right = 1;
end

% Desired level
set.desiredLevel = work.userpar7;

% Design bandlimiting filter for 

% load base noise signal and scale to 0 dB
[n_base,Fs] = audioread('..\TFS_Stimuli\consonants\SSN.wav');
n_base = n_base(:,1);
n_base = resample(n_base, def.samplerate, Fs);
n_base = fftfilt(set.lpf,n_base);
set.n_base = n_base*(10^(0/20))/sqrt(mean(n_base.^2));
clear n_base

% noise type and level
set.noiseType = work.userpar8;
set.noiseLevel = set.desiredLevel - work.exppar2;

% Make contra masker from SSN
[noise, Fs] = audioread('Sentences\SSN_male.wav');
noise = resample(noise, def.samplerate, Fs);
set.contra_level = -1000;
set.contra_masker = scale(noise(def.samplerate+[1:10*def.samplerate]),set.contra_level);
set.baseline_noise = scale(noise,30);
set.noise_level_with_baseline = 10*log10(10^(set.noiseLevel/10)-10^(30/10));
set.SSN = scale(noise,set.noise_level_with_baseline);

set.stim_num = 1;
% eof
