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
AudioDevice = PsychPortAudio('GetDevices', 3); 

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
p.events = 24; % May change
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
eventStart    = NaN(p.events, p.runsMax);
stimStart     = NaN(p.events, p.runsMax); 
stimEnd       = NaN(p.events, p.runsMax); 
eventEnd      = NaN(p.events, p.runsMax); 
key_stimStart  = NaN(p.events, p.runsMax); 
key_stimEnd    = NaN(p.events, p.runsMax); 
stimDuration  = NaN(p.events, p.runsMax); 
key_eventStart = NaN(p.events, p.runsMax); 

AbsEvStart   = NaN(p.events, p.runsMax); 
AbsStimStart = NaN(p.events, p.runsMax); 
AbsStimEnd   = NaN(p.events, p.runsMax); 
AbsRxnEnd    = NaN(p.events, p.runsMax); 
AbsEvEnd     = NaN(p.events, p.runsMax); 

respTime = cell(p.events, p.runsMax); 
respKey  = cell(p.events, p.runsMax); 

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
ad_all   = [ad_clear, ad_vocode, ad_base]; 
dur_all  = [dur_clear, dur_vocode, dur_base]; 

%% Create keys
% Each block consists of 8 sentences, 4 high quality and 4 degraded. There
% will also be 4 noise and 4 silence events. 
% Events
key_sent = cellfun((@(x) x(1:2)), stim_all, 'UniformOutput', false); 
temp = cell(1, length(find(key_noise | key_silence))); 
for ii = 1:length(temp); temp{ii} = 'nan'; end
key_sent(key_noise | key_silence) = temp; 
key_sent = cellfun(@str2double, key_sent); 

order_sentence = reshape(1:max(key_sent), p.sentences, p.runsMax); 
order_sentence = (4*(order_sentence-1))+1; 

key_stim = nan(p.events, p.runsMax); 
for rr = 1:p.runsMax
    thesesent = Shuffle(order_sentence(:, rr)); 
    orsr_mf = repelem([0 1 2 3]', 4); 
    % 0 is OF
    % 1 is OM
    % 2 is SF
    % 3 is SM
    
    clearvocode = repmat([0 length(files_clear)]', [8 1]); 
    noi_sil = find(key_noise | key_silence)'; 
    
    order_all = Shuffle(orsr_mf + clearvocode) + thesesent; 
    order_all = Shuffle([order_all; noi_sil]); 
    
    key_stim(:, rr) = order_all; 
end

% Now is a good time to run check_cb
% check_cb

% Timing
key_jitter = rand(p.events, p.runsMax); 

temp = p.epiTime + [0:p.eventTime:((p.events-1)*p.eventTime)]'; %#ok<NBRAK>
key_eventStart = repmat(temp, [1, p.runsMax]); 

key_stimStart = key_eventStart + key_jitter; 
key_stimEnd   = key_stimStart  + rawStimDur(eventKey(:,rr))';
rxnEndKey   = key_stimEnd + p.rxnWindow; 
eventEndKey = key_eventStart + p.eventTime;

%% PTB
[wPtr, rect] = Screen('OpenWindow', 1, 185);
DrawFormattedText(wPtr, 'Please wait, preparing experiment...');
Screen('Flip', wPtr);
centerX = rect(3)/2;
centerY = rect(4)/2;
crossCoords = [-30, 30, 0, 0; 0, 0, -30, 30]; 
HideCursor(); 

pahandle = PsychPortAudio('Open', [], [], [], fs);

if ~Training
    cd(dir_funcs)
    stimulicheck(NumSpStim, eventKey); 
end
cd(dir_exp)

for ii = subj.firstRun:subj.lastRun
    stimDuration(:, ii) = rawStimDur(eventKey(:,ii))'; 
end

% Check if using RTBox or Keyboard
if ~ConnectedToRTBox
    cd(RTBoxLoc)
    RTBox('fake', 1); 
    cd(dir_exp)
end

%% Prepare test
try
    for rr = subj.firstRun:subj.lastRun
        % Wait for first pulse
        DrawFormattedText(wPtr, ['Waiting for first pulse. ', ... 
            protocol_order{rr}, num2str(rr)]); 
        Screen('Flip', wPtr); 
        
        cd(RTBoxLoc)
        RTBox('Clear'); 
        RTBox('UntilTimeout', 1);
        firstPulse(rr) = RTBox('WaitTR'); 

        % Draw onto screen after recieving first pulse
        Screen('DrawLines', wPtr, crossCoords, 2, 0, [centerX, centerY]);
        Screen('Flip', wPtr); 

        % Generate absolute time keys
        AbsEvStart(:, rr)   = firstPulse(rr) + key_eventStart(:,rr); 
        AbsStimStart(:, rr) = firstPulse(rr) + key_stimStart(:,rr); 
        AbsStimEnd(:, rr)   = firstPulse(rr) + key_stimEnd(:,rr); 
        AbsRxnEnd(:, rr)    = firstPulse(rr) + rxnEndKey(:,rr); 
        AbsEvEnd(:, rr)     = firstPulse(rr) + eventEndKey(:,rr); 

        WaitTill(firstPulse(rr) + p.epiTime); 

        %% Present audio stimuli
        for event = 1:p.events
            eventStart(event, rr) = GetSecs(); 

            PsychPortAudio('FillBuffer', pahandle, audio{eventKey(event, rr)});
            WaitTill(AbsStimStart(event, rr)); 

            stimStart(event, rr) = GetSecs; 
            PsychPortAudio('Start', pahandle, 1);
            WaitTill(AbsStimEnd(event, rr)); 
            stimEnd(event, rr) = GetSecs; 
            RTBox('Clear'); 

            [respTime{event, rr}, respKey{event, rr}] = RTBox(AbsRxnEnd(event, rr)); 

            WaitTill(AbsEvEnd(event, rr));    
            eventEnd(event, rr) = GetSecs(); 
        end

        WaitSecs(p.eventTime); 
        runEnd(rr) = GetSecs(); 

        if rr ~= subj.lastRun
            endstring = ['End of run. Press any button when ready to continue. Experimenters, prepare for '...
                protocol_order{rr+1}, num2str(rr+1)]; 
            DrawFormattedText(wPtr, endstring); 
            Screen('Flip', wPtr); 
            RTBox('Clear'); 
            RTBox(inf); 
        end 
            
        
    end
    
catch err
    sca; 
    runEnd(rr) = GetSecs();  %#ok<NASGU>
    cd(dir_funcs)
    OutputData
    cd(dir_scripts)
    PsychPortAudio('Close'); 
    rethrow(err)
end
%% Closing down
Screen('CloseAll');
PsychPortAudio('Close'); 
DisableKeysForKbCheck([]); 

%% Save data
cd(dir_funcs)
disp('Please wait, saving data...')
OutputData
disp('All done!')
cd(dir_scripts)