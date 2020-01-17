%% AM_FST
% Script to run FST for AM Ported from my previous isss_multi script. 
% Author - Matt Heard

% MM/DD/YY -- CHANGELOG
% 08/07/17 -- Started changelog. MH
% 08/07/17 -- Found error in "prepare timing keys" that overwrote eventStartKey
%   and stimStartKey every time code completed a run. Fixed! MH
% 08/09/17 -- Preparing for testing, making sure code looks pretty. MH
% 07/09/19 -- Cloned for new project. Lots of updates to make. MH
% 10/04/19 -- Updates begin in earnest with a new stimuli set. 
% 10/10/19 -- Debug mode added

%% Startup
sca; DisableKeysForKbCheck([]); KbQueueStop; 

try 
    PsychPortAudio('Close'); 
catch
    disp('PsychPortAudio is already closed.')
end

InitializePsychSound
clearvars; clc; 
DEBUG = 0; 
codeStart = GetSecs(); 

if DEBUG
    AudioDevice = PsychPortAudio('GetDevices'); 
    Screen('Preference','VisualDebugLevel', 0); 
end

%% Parameters
if DEBUG
    warning('USING DEBUG DEFAULTS!')
    dlg_ans = {'xxxx', '1', '9', '0'}; 
else
    prompt = {...
        'Subject number (####YL):', ...
        'First run (1)', ... 
        'Last run (9)', ... 
        'RTBox connected (0/1):', ...
        }; 
    dlg_ans = inputdlg(prompt); 
end


subj.Num  = dlg_ans{1};
subj.firstRun = str2double(dlg_ans{2}); 
subj.lastRun  = str2double(dlg_ans{3}); 
ConnectedToRTBox = str2double(dlg_ans{4}); 

% Scan paradigms
% Will likely change after meeting with Xiangrui
p.TR     = 1.000; 
p.epiNum = 8; 

% Number of stimuli. Likely to change. 
% NumSpStim  = 192; 
% NumStim    = 200;
% p.runsMax = 6; 

% Timing
p.runsMax = 9; % Now enough for 9?
p.events = 24; % Events per block, may change
% 4 hard OR
% 4 hard SR
% 4 easy OR
% 4 easy SR
% 4 noise
% 4 silence
p.sentences  = 16; % how many sentence stimuli per block?
p.structures = 144; % how many sentence structures?

p.presTime = 4.000; % 4 seconds
p.jitter   = 0.500; % Maximum stimuli duration is ~2.75 seconds
% Jitter should be half of TR to prevent interference?

p.rxnWindow = 3.000;  % 3 seconds

p.epiTime   = p.TR * p.epiNum;  % 8 seconds because I said so
p.eventTime = p.presTime + p.epiTime;
p.runDuration = p.epiTime + ...   % After first pulse
    p.eventTime * p.events + ...  % Each event
    p.eventTime;                  % After last acquisition

% p.vocode = '15ch'; % how many channels of vocoding?
p.snr = 3; % SNR between babble and speech (hard condition, easy condition tests clear speech)
%  1
% -1
% -3

%% Paths
cd ..
dir_exp = pwd; 

dir_stim = fullfile(dir_exp, 'stimuli');
dir_stim_clear  = fullfile(dir_stim, 'clear'); 
dir_stim_base   = fullfile(dir_stim, 'noi_sil'); 

dir_scripts = fullfile(dir_exp, 'scripts');
dir_results = fullfile(dir_exp, 'results');
dir_funcs   = fullfile(dir_scripts, 'functions');

Instructions = 'instructions_lang.txt';

%% Preallocating timing variables
real_eventStart = NaN(p.events, p.runsMax);
real_stimStart  = NaN(p.events, p.runsMax); 
real_stimEnd    = NaN(p.events, p.runsMax); 
real_eventEnd   = NaN(p.events, p.runsMax); 

real_respTime = cell(p.events, p.runsMax); 
real_respKey  = cell(p.events, p.runsMax); 

abs_eventStart = NaN(p.events, p.runsMax); 
abs_stimStart  = NaN(p.events, p.runsMax); 
abs_stimEnd    = NaN(p.events, p.runsMax); 
abs_rxnEnd     = NaN(p.events, p.runsMax); 
abs_eventEnd   = NaN(p.events, p.runsMax); 

firstPulse = NaN(1, p.runsMax); 
runEnd     = NaN(1, p.runsMax); 

%% File names
results_xlsx = [subj.Num '_lang.xlsx']; 
results_mat  = [subj.Num '_lang.mat']; 

%% Load stimuli
ad_count = 1; 
files_clear = dir(fullfile(dir_stim_clear, '*.wav')); 
fs_clear = zeros(1, length(files_clear)); 
ad_all = cell(1, 2*length(files_clear) + 4*p.runsMax); % clear, vocode, noise, no silence yet

disp('loading clear stimuli...')
for ii = 1:length(files_clear)
    thisfile = fullfile(dir_stim_clear, files_clear(ii).name); 
    [tempAudio, fs_clear(ii)] = audioread(thisfile); 
    ad_all{ad_count} = [tempAudio'; tempAudio']; 
    ad_count = ad_count + 1; 
end

disp('done!')

if ~all(fs_clear(1) == fs_clear)
    error('fs clear is not equal across all stimuli!')
end

fs = fs_clear(1); 

% Add babble to clear stimuli
cfg.noisefile = fullfile(dir_stim, 'babble_track_330m_mono_44100.wav'); 
cfg.prestim  = 0.120;
cfg.poststim = 0.120;
cfg.snrs = p.snr; 
cfg.fs = fs; 
cd(dir_stim)
disp('adding babble...')
ad_babble = jp_addnoise_hwk_mh(dir_stim_clear, cfg); % DOES NOT SET RMS
for ii = 1:length(ad_babble)
    ad_all{ad_count} = [ad_babble{ii}, ad_babble{ii}]'; 
    ad_count = ad_count + 1; 
end

disp('done!')

% Noise 
disp('making noise stimuli...')
noise_stim = randi(length(files_clear), [1, 4*p.runsMax]); 
for ii = noise_stim % number of noise trials
    tempAudio = jp_vocode_mh(ad_all{ii}(1, :), 1, fs); 
    ad_all{ad_count} = [tempAudio; tempAudio]; 
    ad_count = ad_count + 1;     
end

disp('done!')

% Equalize RMS
ad_all_rms = jp_equalizerms_mh(ad_all, 'verbose'); 

% Silence
disp('loading silence...')
tempAudio = audioread('silence01.wav'); 
ad_all_rms{ad_count} = [tempAudio, tempAudio]'; 
disp('done!')

% Clean up
dur_all = cellfun((@(x) length(x)/fs), ad_all_rms); 
if any(dur_all > 3)
    error('stim are too long!')
end

clear ad_all ad_babble tempAudio

%% Create keys
% Each block consists of 16 sentences, 8 high quality and 8 degraded. There
% will also be 4 noise and 4 silence events. 
% key_stim: what stim to play during what event
% key_jitter: how much jitter?
% key_eventStart: when does each event start?
% key_stimStart: when does each stimuli start?
% key_stimEnd
% key_rxnEnd
% key_eventEnd

% Events
events_clear   = 1:length(files_clear); 
events_babble  = length(files_clear)+1:2*length(files_clear);
events_noise   = 2*length(files_clear)+1:2*length(files_clear)+4*p.runsMax; 
matrix_noise   = reshape(events_noise, [4, p.runsMax]); 
events_silence = 2*length(files_clear)+4*p.runsMax+1; 

if events_silence ~= length(ad_all_rms)
    error('stimuli are not loaded correctly?')
end

key_sentence = reshape(Shuffle(1:p.structures), p.sentences, p.runsMax); 
key_sentence = (4*(key_sentence-1))+1; 

key_stim   = nan(p.events, p.runsMax); 
key_answer = nan(p.events, p.runsMax); 
for rr = 1:p.runsMax
    thesesent = key_sentence(:, rr); 
    orsr_mf = repelem([0 1 2 3]', 4); 
    % 0 is OF, actually 1, mod 1, answer 1
    % 1 is OM, actually 2, mod 0, answer 2
    % 2 is SF, actually 3, mod 1, answer 1
    % 3 is SM, actually 4, mod 0, answer 2
    
    clearbabble = repmat([0 length(files_clear)]', [8 1]); 
    noi_sil = [matrix_noise(:, rr); repmat(events_silence, [4 1])]; 
    
    order_all = Shuffle(orsr_mf + clearbabble) + thesesent; 
    order_all = Shuffle([order_all; noi_sil]); 
    
    key_stim(:, rr) = order_all; 
    
    noise   = ismember(order_all, events_noise); 
    silence = ismember(order_all, events_silence);  
    sent_mf = find(~ismember(order_all, [events_noise events_silence])); 
    temp = mod(order_all(sent_mf), 2); 
    male   = temp == 0;
    female = temp == 1;
    
    key_answer(noise, rr) = 3; 
    key_answer(silence, rr) = 0; 
    key_answer(sent_mf(male), rr) = 2; 
    key_answer(sent_mf(female), rr) = 1; 
end

% Now is a good time to run check_cb
% check_cb

% Timing
key_jitter = p.jitter * rand(p.events, p.runsMax); 

temp = p.epiTime + [0:p.eventTime:((p.events-1)*p.eventTime)]'; %#ok<NBRAK>
key_eventStart = repmat(temp, [1, p.runsMax]); 
key_stimStart  = key_eventStart + key_jitter; 
key_stimDur    = dur_all(key_stim); 
key_stimEnd    = key_stimStart  + key_stimDur;
key_rxnEnd     = key_stimEnd    + p.rxnWindow; 
key_eventEnd   = key_eventStart + p.eventTime;

%% PTB
if DEBUG
    [wPtr, rect] = Screen('OpenWindow', 1, 185, [0 0 1280 720]);
else
    [wPtr, rect] = Screen('OpenWindow', 1, 185);
end

DrawFormattedText(wPtr, 'Please wait, preparing experiment...');
Screen('Flip', wPtr);
centerX = rect(3)/2;
centerY = rect(4)/2;
crossCoords = [-30, 30, 0, 0; 0, 0, -30, 30]; 
HideCursor(); 

pahandle = PsychPortAudio('Open', 5, [], [], fs); % check at scanner
RTBox('fake', ~ConnectedToRTBox); 
RTBox('UntilTimeout', 1);
Screen('TextSize', wPtr, 42); 

%% Prepare test
try
    for rr = subj.firstRun:subj.lastRun
        % Wait for first pulse
        DrawFormattedText(wPtr, ['Waiting for first pulse. Run ' num2str(rr)], 'center', 'center'); 
        Screen('Flip', wPtr); 
        
        RTBox('Clear'); 
        firstPulse(rr) = RTBox('WaitTR'); 

        % Draw onto screen after recieving first pulse
        Screen('DrawLines', wPtr, crossCoords, 2, 0, [centerX, centerY]);
        Screen('Flip', wPtr); 

        % Generate absolute time keys
        abs_eventStart(:, rr) = firstPulse(rr) + key_eventStart(:, rr); 
        abs_stimStart(:, rr)  = firstPulse(rr) + key_stimStart(:, rr); 
        abs_stimEnd(:, rr)    = firstPulse(rr) + key_stimEnd(:, rr); 
        abs_rxnEnd(:, rr)     = firstPulse(rr) + key_rxnEnd(:, rr); 
        abs_eventEnd(:, rr)   = firstPulse(rr) + key_eventEnd(:, rr); 
        
        WaitTill(firstPulse(rr) + p.epiTime); 
        
        %% Present audio stimuli
        for ev = 1:p.events
            real_eventStart(ev, rr) = GetSecs(); 

            PsychPortAudio('FillBuffer', pahandle, ad_all_rms{key_stim(ev, rr)});
            WaitTill(abs_stimStart(ev, rr)-0.1); 
            
            real_stimStart(ev, rr) = PsychPortAudio('Start', pahandle, ... 
                1, abs_stimStart(ev, rr), 1);
            WaitTill(abs_stimEnd(ev, rr)); 
            real_stimEnd(ev, rr) = GetSecs(); 
            RTBox('Clear'); 

            [real_respTime{ev, rr}, real_respKey{ev, rr}] = RTBox(abs_rxnEnd(ev, rr)); 

            WaitTill(abs_eventEnd(ev, rr));    
            real_eventEnd(ev, rr) = GetSecs(); 
        end

        WaitSecs(p.eventTime); 
        runEnd(rr) = GetSecs(); 

        if rr ~= subj.lastRun
            endstring = 'End of run. Press any button when ready to continue.'; 
            DrawFormattedText(wPtr, endstring, 'center', 'center'); 
            Screen('Flip', wPtr); 
            RTBox('Clear'); 
            RTBox(inf); 
        end 
        
    end
    
catch err
    sca; 
    runEnd(rr) = GetSecs();  %#ok<NASGU>
    cd(dir_scripts)
    OutputData
    PsychPortAudio('Close'); 
    rethrow(err)
end

%% Closing down
Screen('CloseAll');
PsychPortAudio('Close'); 
DisableKeysForKbCheck([]); 

%% Save data
cd(dir_scripts)
disp('Please wait, saving data...')
OutputData
disp('All done!')
