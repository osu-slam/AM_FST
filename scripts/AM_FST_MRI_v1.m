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
NumSpStim = 192; 
NumStim   = 200; 

% Timing
p.runs   = length(subj.firstRun:subj.lastRun); 
p.events = 16; % May change

p.presTime   = 4.000;  % 4 seconds
p.epiTime    = 10.000; % 10 seconds
p.eventTime  = p.presTime + p.epiTime;

p.runDuration = p.epiTime + ...   % After first pulse
    p.eventTime * p.events + ...  % Each event
    p.eventTime;                  % After last acquisition

p.rxnWindow = 3.000;  % 3 seconds
p.jitWindow = 0.900;  % 1 second, see notes below
    % For this experiment, the first second of the silent window will not
    % have stimuli presented. To code for this, I add an additional 1 s
    % to the jitterKey. So, the jitter window ranges from 1 s to 2 s.

% Training override
if Training
    NumSpStim = 4; 
    NumStim   = 5; 
    p.events  = 5; 
end
    
%% Paths
cd ..
expDir = pwd; 

StimuliLoc = [expDir, '\stimuli_lang'];
ScriptsLoc = [expDir, '\scripts'];
ResultsLoc = [expDir, '\results'];
FuncsLoc   = [ScriptsLoc, '\functions'];

cd ..
RTBoxLoc = [pwd, '\RTBox']; 

Instructions = 'instructions_lang.txt';

%% Preallocating timing variables
if Training
    maxNumRuns = 1;
else
    maxNumRuns = 6; 
end

eventStart    = NaN(p.events, maxNumRuns);
stimStart     = NaN(p.events, maxNumRuns); 
stimEnd       = NaN(p.events, maxNumRuns); 
eventEnd      = NaN(p.events, maxNumRuns); 
stimStartKey  = NaN(p.events, maxNumRuns); 
stimEndKey    = NaN(p.events, maxNumRuns); 
stimDuration  = NaN(p.events, maxNumRuns); 
eventStartKey = NaN(p.events, maxNumRuns); 

AbsEvStart   = NaN(p.events, maxNumRuns); 
AbsStimStart = NaN(p.events, maxNumRuns); 
AbsStimEnd   = NaN(p.events, maxNumRuns); 
AbsRxnEnd    = NaN(p.events, maxNumRuns); 
AbsEvEnd     = NaN(p.events, maxNumRuns); 

respTime = cell(p.events, maxNumRuns); 
respKey  = cell(p.events, maxNumRuns); 

firstPulse = NaN(1, maxNumRuns); 
runEnd     = NaN(1, maxNumRuns); 

%% File names
if Training
    filetag = [subj.Num '_' subj.Init '_practice_']; 
else
    filetag = [subj.Num '_' subj.Init '_']; 
end

ResultsTxt = [filetag 'lang_results.txt']; 
ResultsXls = [filetag 'lang_results.xlsx']; 
Variables  = [filetag 'lang_variables.mat']; 
    
%% Determine protocol
if Training
    protocol_order = {'multi'}; 
else
    cd(ScriptsLoc)
    load('master_scan_order.mat')
    protocol_order = cell(1, 6); 
    for s = 1:10
        if master{1, s} == subj.Num
            for k = 1:6
                protocol_order{k} = master{k+1, s}; 
            end
        end
    end
    cd(expDir)
end

%% Load PTB and stimuli
% PTB
[wPtr, rect] = Screen('OpenWindow', 1, 185);
DrawFormattedText(wPtr, 'Please wait, preparing experiment...');
Screen('Flip', wPtr);
centerX = rect(3)/2;
centerY = rect(4)/2;
crossCoords = [-30, 30, 0, 0; 0, 0, -30, 30]; 
HideCursor(); 

% Stimuli, check counterbalance
cd(FuncsLoc) 
[audio, fs, rawStimDur, jitterKey, eventKey, answerKey, speechKey] = ...
    LoadStimAndKeys(StimuliLoc, p.events, subj.firstRun, subj.lastRun, NumSpStim, maxNumRuns, Training);
fs = fs{1}; % Above func checks that all fs are the same.  

pahandle = PsychPortAudio('Open', [], [], [], fs);

if ~Training
    cd(FuncsLoc)
    stimulicheck(NumSpStim, eventKey); 
end
cd(expDir)

for i = subj.firstRun:subj.lastRun
    stimDuration(:, i) = rawStimDur(eventKey(:,i))'; 
end

% Check if using RTBox or Keyboard
if ~ConnectedToRTBox
    cd(RTBoxLoc)
    RTBox('fake', 1); 
    cd(expDir)
end

%% Prepare test
try
    for run = subj.firstRun:subj.lastRun

        DrawFormattedText(wPtr, 'Please wait, preparing run...');
        Screen('Flip', wPtr); 

        % Prepare timing keys
        eventStartKey(:, run) = p.epiTime + [0:p.eventTime:((p.events-1)*p.eventTime)]'; %#ok<NBRAK>
        stimStartKey(:, run)  = eventStartKey(:, run) + jitterKey(:, run); 

        if Training
            stimEndKey = stimStartKey + rawStimDur(eventKey)';
        else
            stimEndKey(:, run) = stimStartKey(:, run) + rawStimDur(eventKey(:,run))';
        end

        rxnEndKey   = stimEndKey + p.rxnWindow; 
        eventEndKey = eventStartKey + p.eventTime;

        % Display instructions
        if Training
            cd(FuncsLoc)
            DisplayInstructions_bkfw_rtbox(Instructions, wPtr, RTBoxLoc); 
            cd(expDir)
        end


        % Wait for first pulse
        DrawFormattedText(wPtr, ['Waiting for first pulse. ', ... 
            protocol_order{run}, num2str(run)]); 
        Screen('Flip', wPtr); 
        
        cd(RTBoxLoc)
        RTBox('Clear'); 
        RTBox('UntilTimeout', 1);
        firstPulse(run) = RTBox('WaitTR'); 

        % Draw onto screen after recieving first pulse
        Screen('DrawLines', wPtr, crossCoords, 2, 0, [centerX, centerY]);
        Screen('Flip', wPtr); 

        % Generate absolute time keys
        AbsEvStart(:, run)   = firstPulse(run) + eventStartKey(:,run); 
        AbsStimStart(:, run) = firstPulse(run) + stimStartKey(:,run); 
        AbsStimEnd(:, run)   = firstPulse(run) + stimEndKey(:,run); 
        AbsRxnEnd(:, run)    = firstPulse(run) + rxnEndKey(:,run); 
        AbsEvEnd(:, run)     = firstPulse(run) + eventEndKey(:,run); 

        WaitTill(firstPulse(run) + p.epiTime); 

        %% Present audio stimuli
        for event = 1:p.events
            eventStart(event, run) = GetSecs(); 

            PsychPortAudio('FillBuffer', pahandle, audio{eventKey(event, run)});
            WaitTill(AbsStimStart(event, run)); 

            stimStart(event, run) = GetSecs; 
            PsychPortAudio('Start', pahandle, 1);
            WaitTill(AbsStimEnd(event, run)); 
            stimEnd(event, run) = GetSecs; 
            RTBox('Clear'); 

            [respTime{event, run}, respKey{event, run}] = RTBox(AbsRxnEnd(event, run)); 

            WaitTill(AbsEvEnd(event, run));    
            eventEnd(event, run) = GetSecs(); 
        end

        WaitSecs(p.eventTime); 
        runEnd(run) = GetSecs(); 

        if run ~= subj.lastRun
            endstring = ['End of run. Press any button when ready to continue. Experimenters, prepare for '...
                protocol_order{run+1}, num2str(run+1)]; 
            DrawFormattedText(wPtr, endstring); 
            Screen('Flip', wPtr); 
            RTBox('Clear'); 
            RTBox(inf); 
        end 
            
        
    end
    
catch err
    sca; 
    runEnd(run) = GetSecs();  %#ok<NASGU>
    cd(FuncsLoc)
    OutputData
    cd(ScriptsLoc)
    PsychPortAudio('Close'); 
    rethrow(err)
end
%% Closing down
Screen('CloseAll');
PsychPortAudio('Close'); 
DisableKeysForKbCheck([]); 

%% Save data
cd(FuncsLoc)
disp('Please wait, saving data...')
OutputData
disp('All done!')
cd(ScriptsLoc)