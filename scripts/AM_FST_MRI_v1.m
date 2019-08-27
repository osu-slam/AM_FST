%% AM_FST
% Script to run FST for AM Ported from my previous isss_multi script. 
% Author - Matt Heard

% MM/DD/YY -- CHANGELOG
% 08/07/17 -- Started changelog. MH
% 08/07/17 -- Found error in "prepare timing keys" that overwrote eventStartKey
%   and stimStartKey every time code completed a run. Fixed! MH
% 08/09/17 -- Preparing for testing, making sure code looks pretty. MH
% 07/09/19 -- Cloned for new project. Lots of updates to make. MH

%% Startup
sca; DisableKeysForKbCheck([]); KbQueueStop; 
Screen('Preference','VisualDebugLevel', 0); 

try 
    PsychPortAudio('Close'); 
catch
    disp('PsychPortAudio is already closed.')
end

InitializePsychSound
clearvars; clc; 
codeStart = GetSecs(); 
% AudioDevice = PsychPortAudio('GetDevices'); % For debugging

%% Parameters
prompt = {...
    'Subject number (####YL):', ...
    'First run (TBD)', ... 
    'Last run (TBD)', ... 
    'RTBox connected (0/1):', ...
    }; 
dlg_ans = inputdlg(prompt); 

subj.Num  = dlg_ans{1};
subj.firstRun = str2double(dlg_ans{2}); 
subj.lastRun  = str2double(dlg_ans{3}); 
ConnectedToRTBox = str2double(dlg_ans{4}); 

% Scan paradigms
% Will likely change after meeting with Xiangrui/Tatiana
scan.TR     = 1.000; 
scan.epiNum = 8; 

% Number of stimuli. Likely to change. 
% NumSpStim  = 192; 
% NumStim    = 200;
% p.runsMax = 6; 

% Timing
p.runsMax = 4; % ONLY ENOUGH SENTENCES FOR 4 BLOCKS RIGHT NOW
p.events = 24; % Events per block, may change
% 4 high OR
% 4 high SR
% 4 low OR
% 4 low SR
% 4 noise
% 4 silence
p.sentences = 16; % how many sentence stimuli per block?

p.presTime = 4.000;  % 4 seconds
p.jitter   = 1.000; % Maximum stimuli duration is ~2.75 seconds

p.rxnWindow = 3.000;  % 3 seconds

p.epiTime   = 8.000;  % 8 seconds because I said so
p.eventTime = p.presTime + p.epiTime;
p.runDuration = p.epiTime + ...   % After first pulse
    p.eventTime * p.events + ...  % Each event
    p.eventTime;                  % After last acquisition

p.vocode = '15ch'; % how many channels of vocoding?
    
%% Paths
cd ..
dir_exp = pwd; 

dir_stim = fullfile(dir_exp, 'stimuli');
dir_stim_clear  = fullfile(dir_stim, 'clear_resample'); 
dir_stim_vocode = fullfile(dir_stim, p.vocode); 
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
% Clear
files_clear = dir(fullfile(dir_stim_clear, '*.wav')); 
ad_clear = cell(1, length(files_clear)); 
fs_clear = zeros(1, length(files_clear)); 

for ii = 1:length(files_clear)
    thisfile = fullfile(dir_stim_clear, files_clear(ii).name); 
    [ad_clear{ii}, fs_clear(ii)] = audioread(thisfile); 
    ad_clear{ii} = [ad_clear{ii}'; ad_clear{ii}']; 
end

if ~all(fs_clear(1) == fs_clear)
    error('fs clear is not equal across all stimuli!')
end

fs_clear = fs_clear(1); 

% Vocode
files_vocode = dir(fullfile(dir_stim_vocode, '*.wav')); 
ad_vocode = cell(1, length(files_vocode)); 
fs_vocode = zeros(1, length(files_vocode)); 

for ii = 1:length(files_vocode)
    thisfile = fullfile(dir_stim_vocode, files_vocode(ii).name); 
    [ad_vocode{ii}, fs_vocode(ii)] = audioread(thisfile); 
    ad_vocode{ii} = [ad_vocode{ii}'; ad_vocode{ii}']; 
end

if ~all(fs_vocode(1) == fs_vocode)
    error('fs vocode is not equal across all stimuli!')
end

fs_vocode = fs_vocode(1); 

% Noise and silence
files_base = dir(fullfile(dir_stim_base, '*.wav')); 
ad_base = cell(1, length(files_base)); 
fs_base = zeros(1, length(files_base)); 

for ii = 1:length(files_base)
    thisfile = fullfile(dir_stim_base, files_base(ii).name); 
    [ad_base{ii}, fs_base(ii)] = audioread(thisfile); 
    ad_base{ii} = [ad_base{ii}'; ad_base{ii}']; 
end

if ~all(fs_base(1) == fs_base)
    error('fs vocode is not equal across all stimuli!')
end

fs_base = fs_base(1); 

% Clean up
if any(fs_clear ~= [fs_vocode fs_base])
    error('fs clear vocode and base are not equal!')
end

fs = fs_clear; 
dur_clear  = cellfun((@(x) length(x)/fs), ad_clear); 
dur_vocode = cellfun((@(x) length(x)/fs), ad_vocode); 
dur_base   = cellfun((@(x) length(x)/fs), ad_base); 

if any(dur_clear > 3) || any(dur_vocode > 3) || any(dur_base > 3)
    error('stim are too long!')
end

stim_all = [files_clear; files_vocode; files_base]; stim_all = {stim_all.name};
num_files_clear = length(files_clear); % how many sentences
clear files_clear files_vocode files_base
ad_all   = [ad_clear, ad_vocode, ad_base]; 
clear ad_clear ad_vocode ad_base
dur_all  = [dur_clear, dur_vocode, dur_base]; 
clear dur_clear dur_vocode dur_base

%% Create keys
% Each block consists of 8 sentences, 4 high quality and 4 degraded. There
% will also be 4 noise and 4 silence events. 
% key_stim: what stim to play during what event
% key_jitter: how much jitter?
% key_eventStart: when does each event start?
% key_stimStart: when does each stimuli start?
% key_stimEnd
% key_rxnEnd
% key_eventEnd

% Events
events_noise   = cellfun((@(x) contains(x, '01ch')), stim_all);
events_silence = cellfun((@(x) contains(x, 'silence')), stim_all);

events_sent    = cellfun((@(x) x(1:2)), stim_all, 'UniformOutput', false); 
temp = cell(1, length(find(events_noise | events_silence))); 
for ii = 1:length(temp); temp{ii} = 'nan'; end
events_sent(events_noise | events_silence) = temp; 
events_sent = cellfun(@str2double, events_sent); 

order_sentence = reshape(1:max(events_sent), p.sentences, p.runsMax); 
order_sentence = (4*(order_sentence-1))+1; 

key_stim   = nan(p.events, p.runsMax); 
key_answer = nan(p.events, p.runsMax); 
for rr = 1:p.runsMax
    thesesent = Shuffle(order_sentence(:, rr)); 
    orsr_mf = repelem([0 1 2 3]', 4); 
    % 0 is OF, actually 1, mod 1, answer 1
    % 1 is OM, actually 2, mod 0, answer 2
    % 2 is SF, actually 3, mod 1, answer 1
    % 3 is SM, actually 4, mod 0, answer 2
    
    clearvocode = repmat([0 num_files_clear]', [8 1]); 
    noi_sil = find(events_noise | events_silence)'; 
    
    order_all = Shuffle(orsr_mf + clearvocode) + thesesent; 
    order_all = Shuffle([order_all; noi_sil]); 
    
    key_stim(:, rr) = order_all; 
    
    noise   = ismember(order_all, find(events_noise)); 
    silence = ismember(order_all, find(events_silence));  
    sent_mf = find(~ismember(order_all, [find(events_noise) find(events_silence)])); 
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
key_jitter = rand(p.events, p.runsMax); 

temp = p.epiTime + [0:p.eventTime:((p.events-1)*p.eventTime)]'; %#ok<NBRAK>
key_eventStart = repmat(temp, [1, p.runsMax]); 

key_stimStart = key_eventStart + key_jitter; 
key_stimDur   = dur_all(key_stim); 
key_stimEnd   = key_stimStart  + key_stimDur;
key_rxnEnd    = key_stimEnd    + p.rxnWindow; 
key_eventEnd  = key_eventStart + p.eventTime;

%% PTB
[wPtr, rect] = Screen('OpenWindow', 1, 185);
DrawFormattedText(wPtr, 'Please wait, preparing experiment...');
Screen('Flip', wPtr);
centerX = rect(3)/2;
centerY = rect(4)/2;
crossCoords = [-30, 30, 0, 0; 0, 0, -30, 30]; 
HideCursor(); 

pahandle = PsychPortAudio('Open', 5, [], [], fs); % check at scanner
RTBox('fake', ~ConnectedToRTBox); 
RTBox('UntilTimeout', 1);

%% Prepare test
try
    for rr = subj.firstRun:subj.lastRun
        % Wait for first pulse
        DrawFormattedText(wPtr, ['Waiting for first pulse. Run ' num2str(rr)]); 
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

            PsychPortAudio('FillBuffer', pahandle, ad_all{key_stim(ev, rr)});
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
            DrawFormattedText(wPtr, endstring); 
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
