%% AM_FST
% Script to run FST for AM Ported from my previous isss_multi script. 
% Author - Matt Heard

% MM/DD/YY -- CHANGELOG
% 08/07/17 -- Started changelog. MH
% 08/07/17 -- Found error in "prepare timing keys" that overwrote eventStartKey
%   and stimStartKey every time code completed a run. Fixed! MH
% 08/09/17 -- Preparing for testing, making sure code looks pretty. MH
% 07/09/19 -- Cloned for new project. Lots of updates to make. MH
%   Split into MRI and Behavioral versions. 

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
    'Subject number (#):', ...
    'First run (TBD)', ... 
    'Last run (TBD)', ... 
    'Instructions? (0/1)', ...
    }; 
dlg_ans = inputdlg(prompt); 

subj.Num  = dlg_ans{1}; 
subj.firstRun = str2double(dlg_ans{2}); 
subj.lastRun  = str2double(dlg_ans{3}); 
Instructions = str2double(dlg_ans{4}); 

% Number of stimuli. Likely to change. 
NumSpStim = 192; 
NumStim   = 200; 

% Timing
p.runs   = length(subj.firstRun:subj.lastRun); 
p.eventsPerRun = 16; % May change
    
%% Paths
dir_scripts = pwd; 
cd ..
dir_exp = pwd; 

dir_stim = fullfile(dir_exp, 'stimuli');
dir_results = fullfile(dir_exp, 'results');
dir_funcs   = fullfile(dir_scripts, 'functions');

%% Preallocating timing variables
respTime = cell(p.events, maxNumRuns); 
respKey  = cell(p.events, maxNumRuns); 

%% File names
if subj.Num < 10
    filetag = ['000' num2str(subj.Num)]; 
elseif subj.Num < 100
    filetag = ['00' num2str(subj.Num)]; 
elseif subj.Num < 1000
    filetag = ['0' num2str(subj.Num)]; 
elseif subj.Num < 10000
    filetag = num2str(subj.Num); 
else
    error('Subject number too high!')
end

if ~isempty(mfilename())
    results_xlsx = fullfile(dir_results, [filetag, '_', mfilename(), '.xlsx']); 
    results_mat  = fullfile(dir_results, [filetag, '_', mfilename(), '.mat']); 
else % if testing... 
    results_xlsx = fullfile(dir_results, [filetag, '_TEST.xlsx']); 
    results_mat  = fullfile(dir_results, [filetag, '_TEST.mat']); 
end

%% Load PTB and stimuli
% PTB
[wPtr, rect] = Screen('OpenWindow', 1, 185);
DrawFormattedText(wPtr, 'Please wait, preparing experiment...');
Screen('Flip', wPtr);
centerX = rect(3)/2;
centerY = rect(4)/2;
% crossCoords = [-30, 30, 0, 0; 0, 0, -30, 30]; 
HideCursor(); 
RTBox('fake', 1); 
RTBox('UntilTimeout', 1);

% Stimuli, check counterbalance
% Needs updating
cd(dir_funcs) 
% fs = fs{1}; % Above func checks that all fs are the same.  

pahandle = PsychPortAudio('Open', [], [], [], fs);

%% Instructions
if Instructions
    
    instructions = fullfile(dir_scripts, 'instructions_BEHAV.txt'); 
    % need to be modified based on whether sentence is vocoded or
    % speech-in-noise
    fid = fopen(instructions);
    ii = 1;
    while 1
        line = fgetl(fid);
        if line == -1
            break
        end

        inst_lines{ii} = line; %#ok<SAGROW>
        ii = ii + 1;
    end
    
    noClear = [0 0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 0];
    for ii = 1:18
        if any(ii == [4, 6, 9, 11, 14, 16])
            DrawFormattedText(wPtr, inst_lines{ii}, 'center', centerY + 200, 255);
        else
            DrawFormattedText(wPtr, inst_lines{ii}, 'center', 'center', 255);
        end

        Screen('Flip', wPtr, [], noClear(ii));
        WaitSecs(0.5);

        %%% NEEDS UPDATING FOR STIMULI
        if ii == 4
            PsychPortAudio('FillBuffer', pahandle, audio_pract{7});
            PsychPortAudio('Start', pahandle);
            % Brothers that tutor sisters are helpful. 15ch
        elseif ii == 6
            PsychPortAudio('FillBuffer', pahandle, audio_pract{8});
            PsychPortAudio('Start', pahandle);
            % Brothers that tutor sisters are helpful. 24ch
        elseif ii == 8
            PsychPortAudio('FillBuffer', pahandle, audio_pract{7});
            PsychPortAudio('Start', pahandle);
            % Brothers that tutor sisters are helpful. 15ch
        elseif ii == 10
            PsychPortAudio('FillBuffer', pahandle, audio_pract{6});
            PsychPortAudio('Start', pahandle);
            % Sisters that tutor brothers are helpful. 24ch
        elseif ii == 13
            PsychPortAudio('FillBuffer', pahandle, audio_pract{25});
            PsychPortAudio('Start', pahandle);
            % Nephews that nieces welcome are gracious. 15ch
        elseif ii == 15
            PsychPortAudio('FillBuffer', pahandle, audio_pract{28});
            PsychPortAudio('Start', pahandle);
            % Nieces that nephews welcome are gracious. 24ch
        end

        RTBox('Clear');
        RTBox(inf);
    end

    Screen('Flip', wPtr);
    
%% Practice
% STILL NEED TO EXTRACT ANSWER KEY FOR PRACTICE!
    while 1
        correct = 0;

        DrawFormattedText(wPtr, 'Start!', 'center', 'center', 255);
        Screen('Flip', wPtr);
        WaitTill(GetSecs() + 0.5);

        for evt = 1:p.stimPerBlock
            WaitTill(GetSecs() + 0.5);
            Screen('DrawLines', wPtr, crossCoords, 2, 255, [centerX, centerY]);
            Screen('Flip', wPtr); 

            stimEnd = GetSecs() + dur_pract(key_pract(evt));
            PsychPortAudio('FillBuffer', pahandle, audio_pract{key_pract(evt)});
            PsychPortAudio('Start', pahandle);

            DrawFormattedText(wPtr, 'female', centerX - 500, 'center', 255);
            DrawFormattedText(wPtr, 'male', centerX + 500, 'center', 255);
            WaitTill(stimEnd); 
            
            Screen('Flip', wPtr);

            RTBox('Clear'); 
            windowStart = GetSecs();
            [~, answer] = RTBox(windowStart + 5); 

            if strcmp('', answer) % If subject timed out
                DrawFormattedText(wPtr, 'Too slow! Be sure to respond quicker.', 'center', 'center', 255);
            elseif strcmp(key_pract_direction{evt}, answer) % If correct
                correct = correct + 1;
                DrawFormattedText(wPtr, 'You are correct! Good job!', 'center', 'center', 255);
            else % If wrong
                DrawFormattedText(wPtr, 'Oops, wrong answer!', 'center', 'center', 255);
            end

            Screen('Flip', wPtr);
            WaitTill(GetSecs() + 1);
            Screen('Flip', wPtr);
        end

        % Feedback
        correct_trials = sprintf(inst_lines{19}, num2str(correct));
        DrawFormattedText(wPtr, correct_trials, 'center', 'center', 255);
        Screen('Flip', wPtr);

        WaitSecs(0.5);
        RTBox('Clear');
        [~, cont] = RTBox(inf);
        if strcmp(cont, 'space')
            DrawFormattedText(wPtr, inst_lines{20}, 'center', 'center', 255);
            Screen('Flip', wPtr);
            WaitTill(GetSecs + 0.5);
            RTBox('Clear');
            RTBox(inf);
            % To prevent the practice condition from interfering with the
            % rest of the experiment, I've inserted a 20 second break where
            % the experiment is "loading". 
            DrawFormattedText(wPtr, 'Please wait...', 'center', 'center', 255);
            Screen('Flip', wPtr);
            WaitTill(GetSecs() + 20);
            break
        else
            Screen('Flip', wPtr);
        end

    end
    
end

%% Main trials
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
            cd(dir_funcs)
            DisplayInstructions_bkfw_rtbox(instructions, wPtr, RTBoxLoc); 
            cd(dir_exp)
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