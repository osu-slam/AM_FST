%% AM_FST_snr_pilot
% Script to run FST for AM Ported from my previous isss_multi script. 
% Cloned on 11/14/19 to pilot SNR in seniors. 
% Author - Matt Heard

% MM/DD/YY -- CHANGELOG
% 08/07/17 -- Started changelog. MH
% 08/07/17 -- Found error in "prepare timing keys" that overwrote eventStartKey
%   and stimStartKey every time code completed a run. Fixed! MH
% 08/09/17 -- Preparing for testing, making sure code looks pretty. MH
% 07/09/19 -- Cloned for new project. Lots of updates to make. MH
% 10/04/19 -- Updates begin in earnest with a new stimuli set. 
% 10/10/19 -- Debug mode added
% 11/14/19 -- Cloned from AM_FST_MRI to pilot the correct SNR for seniors.
%   MH
% 11/20/19 -- Version 2 tackles the memory issues that happen when we load
%   1000 audio files. MH

% TODO:

%% Startup
sca; DisableKeysForKbCheck([]); KbQueueStop; 

try 
    PsychPortAudio('Close'); 
catch
    disp('PsychPortAudio is already closed.')
end

InitializePsychSound
clearvars; clc; 
DEBUG = 1; 
codeStart = GetSecs(); 

if DEBUG
    AudioDevice = PsychPortAudio('GetDevices'); 
    Screen('Preference','VisualDebugLevel', 0); 
end

%% Parameters
% We have about 45 minutes to test this task as we are not collecting any
% other metrics. Also, it's a behavioral task so the timing is much easier.
% 

if DEBUG
    warning('USING DEBUG DEFAULTS!')
    dlg_ans = {'xxxx', '1', '4'}; 
else
    prompt = {...
        'Subject number (####YL):', ...
        'First run (1)', ... 
        'Last run (4)', ... 
        }; 
    dlg_ans = inputdlg(prompt); 
end


subj.Num  = dlg_ans{1};
subj.firstRun = str2double(dlg_ans{2}); 
subj.lastRun  = str2double(dlg_ans{3}); 
ConnectedToRTBox = 0; 

p.snr = [-3, -1, 1, 3]; % SNR between babble and speech 
% Piloting to determine proper SNR for artificial speech and naturalistic
% babble. Let's check many values, starting with [-3, -1, 1, 3]--too hard,
% now on -2, -1, 0, 1
% Four blocks means we could test up to 36 sentences in each block. 

% Testing 144 sentences, maaaybe 6 second per trial... that's a total of 15
% minutes in the worst case scenario!
p.structures = 144; % how many sentence structures?

p.runsMax = length(p.snr); 
p.events = floor(p.structures/length(p.snr)); % Events per block, may change

p.jitter    = 2; % Maximum stimuli duration is ~2.75 seconds
p.rxnWindow = 10.000; % 10 seconds, starting after the sentence starts to play
p.fs = 44100; 

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
% real_stimEnd    = NaN(p.events, p.runsMax); 
real_eventEnd   = NaN(p.events, p.runsMax); 

real_respTime = cell(p.events, p.runsMax); 
real_respKey  = cell(p.events, p.runsMax); 

% abs_eventStart = NaN(p.events, p.runsMax); 
% abs_stimStart  = NaN(p.events, p.runsMax); 
% abs_stimEnd    = NaN(p.events, p.runsMax); 
% abs_rxnEnd     = NaN(p.events, p.runsMax); 
% abs_eventEnd   = NaN(p.events, p.runsMax); 

firstPulse = NaN(1, p.runsMax); 
runEnd     = NaN(1, p.runsMax); 

%% File names
results_xlsx = [subj.Num '_lang_pilot.xlsx']; 
results_mat  = [subj.Num '_lang_pilot.mat']; 

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
% Dramatically simplified since we are just assessing babble!
key_sentence = reshape(Shuffle(1:p.structures), p.events, p.runsMax); 
key_sentence_struct = (4*(key_sentence-1))+1; 
key_sentence_block  = mod(key_sentence-1, p.events)+1; 

key_stim   = nan(p.events, p.runsMax); 
key_answer = cell(p.events, p.runsMax); 
for rr = 1:p.runsMax
    thesesent = key_sentence_struct(:, rr); 
    orsr_mf = Shuffle(repelem([0 1 2 3]', p.events/4)); 
    % 0 is OF, actually 1, mod 1, answer 1
    % 1 is OM, actually 2, mod 0, answer 2
    % 2 is SF, actually 3, mod 1, answer 1
    % 3 is SM, actually 4, mod 0, answer 2
    order_all = orsr_mf + thesesent;
    
    key_stim(:, rr) = order_all; 
    
    temp = mod(order_all, 2); 
    male   = temp == 0;
    female = temp == 1;
    
    for ev = 1:p.events
        if male(ev) == 1
            key_answer{ev, rr} = 'right'; 
        elseif female(ev) == 1
            key_answer{ev, rr} = 'left'; 
        end
        
    end
    
end

% Timing
% So much easier now!
key_jitter = 2 + p.jitter * rand(p.events, p.runsMax); 
whichBlock = Shuffle(1:4); 

%% Load stimuli
these_sent = unique(key_stim); % find which stim to load
files_clear = dir(fullfile(dir_stim_clear, '*.wav')); 
these_files = files_clear(these_sent); 
these_files = cellfun((@(x, y) [x filesep y]), {these_files(:).folder}, {these_files(:).name}, 'UniformOutput', false)'; 

ad_all = cell(length(p.snr), p.events); % babble only
disp('done!')

% Add babble to clear stimuli
% When creating many SNR, things are different...
% cfg.noisefile = fullfile(dir_stim, 'babble_track_330m_mono_44100.wav'); 
disp('adding babble...')
cfg.noisefile = fullfile(dir_stim, 'babble_track_330m_mono_44100.wav'); 
cfg.fs = p.fs; 
cfg.stereo = 1; 
cd(dir_stim)

for ii = 1:p.runsMax
    this_run = p.events*(ii-1)+1:p.events*ii; 
    cfg.snrs = p.snr(ii); 

    ad_babble = jp_addnoise_hwk_mh(these_files(this_run), cfg); % DOES NOT SET RMS
    ad_all(ii, :) = ad_babble; 
    clear ad_babble
end

disp('done!')

% Equalize RMS
if DEBUG
    ad_all_rms = jp_equalizerms_mh_v2(ad_all, 'verbose'); 
else
    ad_all_rms = jp_equalizerms_mh_v2(ad_all); 
end

% Clean up
dur_all = cellfun((@(x) length(x)/p.fs), ad_all_rms); 

clear ad_all files_clear

%% PTB
if DEBUG
    [wPtr, rect] = Screen('OpenWindow', 1, 185, [0 0 1280 720]);
else
    [wPtr, rect] = Screen('OpenWindow', 1, 185);
    HideCursor(); 
end

DrawFormattedText(wPtr, 'Please wait, preparing experiment...');
Screen('Flip', wPtr);
centerX = rect(3)/2;
centerY = rect(4)/2;
% crossCoords = [-30, 30, 0, 0; 0, 0, -30, 30]; 
maleX   = rect(3)/3; 
femaleX = 2*rect(3)/3; 

pahandle = PsychPortAudio('Open', 5, [], [], p.fs); % check at scanner
RTBox('fake', ~ConnectedToRTBox); 
RTBox('UntilTimeout', 1);
Screen('TextSize', wPtr, 42); 

RTBox('ButtonNames',{'left' 'right' 'space' '4'}); 

%% Prepare test
try
    for rr = subj.firstRun:subj.lastRun
        thisBlock = ad_all_rms(whichBlock(rr), :); 
        
        % Wait for first pulse
        DrawFormattedText(wPtr, ['Press space to begin block ' num2str(rr)], 'center', 'center'); 
        Screen('Flip', wPtr); 
        
        RTBox('Clear'); 
        firstPulse(rr) = RTBox(inf); 

        % Draw onto screen after recieving first pulse
%         Screen('DrawLines', wPtr, crossCoords, 2, 0, [centerX, centerY]);
        DrawFormattedText(wPtr, 'Female', 'center', 'center', [], [], [], [], [], [], [0 0 centerX rect(4)]);
        DrawFormattedText(wPtr, 'Male', 'center', 'center', [], [], [], [], [], [], [centerX 0 rect(3) rect(4)]);
        Screen('Flip', wPtr); 

        % Generate absolute time keys
        WaitTill(firstPulse(rr) + 5); 
        
        %% Present audio stimuli
        for ev = 1:p.events
            real_eventStart(ev, rr) = GetSecs(); 
            
            PsychPortAudio('FillBuffer', pahandle, thisBlock{key_sentence_block(ev, rr)});
            WaitTill(GetSecs() + key_jitter(ev, rr)); 
            
            real_stimStart(ev, rr) = PsychPortAudio('Start', pahandle, ... 
                1, [], 1);
            RTBox('Clear'); 

            [real_respTime{ev, rr}, real_respKey{ev, rr}] = RTBox(GetSecs() + p.rxnWindow); 
            real_eventEnd(ev, rr) = GetSecs(); 
        end

        runEnd(rr) = GetSecs(); 

        if rr ~= subj.lastRun
            endstring = 'End of run. Press space to continue.'; 
            DrawFormattedText(wPtr, endstring, 'center', 'center'); 
            Screen('Flip', wPtr); 
            RTBox('Clear'); 
            RTBox(inf); 
            Screen('Flip', wPtr);
            WaitTill(GetSecs()+2); 
        end 
        
    end
    
catch err
    sca; 
    runEnd(rr) = GetSecs();  %#ok<NASGU>
    cd(dir_scripts)
    OutputData_snr_pilot
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
OutputData_snr_pilot
disp('All done!')
