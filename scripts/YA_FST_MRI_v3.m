%% YA_FST
% Script to run FST for YA, ported from my previous isss_multi script. 
% Author - Matt Heard (heardmatthew49@gmail.com)

% MM/DD/YY -- CHANGELOG
% 08/07/17 -- Started changelog. MH
% 08/07/17 -- Found error in "prepare timing keys" that overwrote eventStartKey
%   and stimStartKey every time code completed a run. Fixed! MH
% 08/09/17 -- Preparing for testing, making sure code looks pretty. MH
% 07/09/19 -- Cloned for new project. Lots of updates to make. MH
% 10/04/19 -- Updates begin in earnest with a new stimuli set. 
% 10/10/19 -- Debug mode added
% 01/13/20 -- Forked into YA version:
%   + Switched back to traditional ISSS
%   + Implementing fix for timing error
%   + Adjusted output column for easy analysis later
%   + Removed empty column (Real Stim Duration) because timing broke it
%   ~ Still need to counterbalance stimuli presentation!  
%   ~ New jitter which ensures it is no longer than (max allowed by
%     duration) or 1s long
% 01/15/20 -- Pulled to CCBBI scan computer
%   + New scheme for timing based on discussion w/ Xiangrui
%   + Had to add jp_scripts
% 01/24/20 -- New stimuli and counterbalancing. 
% 01/27/20 -- New stimuli loading done! Changed to female voice stim. 
% 01/28/20 -- Removed jitter. 
% 01/31/20 -- Forked for CCBBI. 
% 01/31/20 -- Brought back to SLAM. New stimuli with SNRs 2 and -2.
%   Modified timing at end. 

%% Startup
clearvars; clc; 
sca; DisableKeysForKbCheck([]); KbQueueStop; 

try 
    PsychPortAudio('Close'); 
catch
    disp('PsychPortAudio is already closed.')
end

InitializePsychSound
DEBUG = 1; 
codeStart = GetSecs(); 
Screen('Preference','SkipSyncTests', 1); 

if DEBUG
    AudioDevice = PsychPortAudio('GetDevices'); 
end

%% Parameters
if DEBUG
    warning('USING DEBUG DEFAULTS!')
    dlg_ans = {'TEST', '1', '1', '0'}; 
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

%% Scan paradigm
% Will likely change after meeting with Xiangrui
% Reverted back to traditional ISSS
p.TR     = 2.000; 
p.epiNum = 4; 
% 40 minutes of scan time!

% Timing
p.runsMax = 9; % Now enough for 9?
p.events = 20; % Events per block, needs updating. 
% Each sentence (O/S) is presented nine times. 3 speeds x 3 SNRs
p.sentences  = 9; % how many sentence stimuli per block?
p.structures = 9; % how many sentence structures?

p.presTime = 4.000; % 4 seconds
p.jitter   = 0; % No jitter because FIR (thanks Xiangrui

p.rxnWindow = 3.000;  % 3 seconds, should we expand this?

p.epiTime   = p.TR * p.epiNum;  % 8 seconds because I said so
p.eventTime = p.presTime + p.epiTime;
p.runDuration = p.epiTime + ...   % After first pulse
    p.eventTime * p.events + ...  % Each event
    p.eventTime;                  % After last acquisition

p.snr = [2 -2]; % SNR between babble and speech
                % (clear, easy, hard)
p.rate = [0.75 1 1.25]; % make sure to put in ascending order

%% Paths
cd ..
dir_exp = pwd; 

dir_stim = fullfile(dir_exp, 'stimuli');
dir_stim_all = fullfile(dir_stim, 'YA_FST_v3_select_babble_v7_norm'); 

dir_scripts = fullfile(dir_exp, 'scripts');
dir_results = fullfile(dir_exp, 'results');

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

abs_eventStart_delta = NaN(p.events, p.runsMax); 
abs_stimStart_delta  = NaN(p.events, p.runsMax); 
abs_stimEnd_delta    = NaN(p.events, p.runsMax); 
abs_rxnEnd_delta     = NaN(p.events, p.runsMax); 
abs_eventEnd_delta   = NaN(p.events, p.runsMax); 

t1 = NaN(p.events + 1, p.runsMax); 

firstPulse = NaN(1, p.runsMax); 
runEnd     = NaN(1, p.runsMax); 

%% File names
results_xlsx = ['YA_' subj.Num '_lang.xlsx']; 
results_mat  = ['YA_' subj.Num '_lang.mat']; 

%% Create keys
% Each block consists of 18 sentences. There will be 9 OR and 9 SR
% sentences (pre-selected) based on 9 sentence structures. There will be 6
% clear, 6 easy babble (high SNR), and 6 hard babble (low SNR) sentences.
% There will be 6 slow (, 6 normal, and 6 fast sentences. 

% KEYS TO THE CODE
% key_events: what to play, and at what point
% key_jitter: how much jitter?
% key_eventStart: when does each event start?
% key_stimStart: when does each stimuli start?
% key_stimEnd
% key_rxnEnd
% key_eventEnd

% Thanks Hyun!!!
runs = 1:p.runsMax;
syn = ["O","S"];
cla = ["", ['SNR' num2str(max(p.snr))], ['SNR' num2str(min(p.snr))]];
spr = string(p.rate); 
spr(strcmp(spr , "1")) = "1.0"; 

ntrial = 9*2*3*3;

cond = strings(ntrial,6);
in = 0;
for rr=1:9
    for ss=1:2
        for cc=1:3
            for pp=1:3
                in = in+1;
                cond(in,1:4) = [runs(rr), syn(ss), cla(cc), spr(pp)];
            end
        end
    end
end

for ss=1:2
    str = CreateLatinSquare(9);
    str = str(randperm(9),:);
    in = 0;
    for cc=1:3
        for pp=1:3
            in=in+1;
            idx = ( cond(:,2)==syn(ss) & cond(:,3)==cla(cc) & cond(:,4)==spr(pp) );
            cond(idx,5) = str(in,:);
        end
    end
end

gen = "*";
for tt=1:ntrial
    cond(tt,6) = strcat("00",cond(tt,5),"_f",cond(tt,2),gen,"_",cond(tt,4));
end

% Add noise/silence and shuffle within each block. 
key_events = []; 
for rr = runs
    thisblock = cond(strcmp(cond(:, 1), num2str(rr)), :); 
    % Randomly grab trial from next block to serve as noise
    
    noise = [num2str(rr), "N", "N", "N", "N", "noise"];
    silence = [num2str(rr), "S", "S", "S", "S", "silence01"];
    
    thisblock = [thisblock; noise; silence]; 
    idx = Shuffle(1:length(thisblock)); 
    
    key_events = [key_events; thisblock(idx, :)]; 
end

snr_easy = ['SNR' num2str(max(p.snr))];
snr_hard = ['SNR' num2str(min(p.snr))];

for ii = 1:length(key_events)
    if strcmp(key_events(ii, 3), "N") % if noise
        key_events(ii, 6) = strcat(key_events(ii, 6), '.wav');
    elseif strcmp(key_events(ii, 3), "S") % if silence
        key_events(ii, 6) = strcat(key_events(ii, 6), '.wav');
    elseif strcmp(key_events(ii, 3), snr_hard) % if hard babble
        key_events(ii, 6) = strcat(key_events(ii, 6), '_', snr_hard, '.wav');
    elseif strcmp(key_events(ii, 3), snr_easy) % if easy babble
        key_events(ii, 6) = strcat(key_events(ii, 6), '_', snr_easy, '.wav');
    elseif strcmp(key_events(ii, 3), "") % if clear
        key_events(ii, 6) = strcat(key_events(ii, 6), '.wav');
    else
        error('exception when preparing key_events!')
    end
    
end

clear cond thisblock

% Timing
% while 1 % use this for loop to ensure jitter is longer than 0.001
%     key_jitter = p.jitter * rand(p.events, p.runsMax); 
%     if all(all(key_jitter > 0.001))
%         break
%     end
% end
key_jitter = 0.1*ones(p.events, p.runsMax); % add slight delay
temp = p.epiTime + [0:p.eventTime:((p.events-1)*p.eventTime)]'; %#ok<NBRAK>
key_eventStart = repmat(temp, [1, p.runsMax]); 
key_stimStart  = key_eventStart + key_jitter; 
key_eventEnd   = key_eventStart + p.eventTime;

%% Load stimuli
% files_all = dir(fullfile(dir_stim_all, '*.wav')); 

ad_events = cell(p.events, p.runsMax); 
fs_events = zeros(p.events, p.runsMax); 
% We are loading audio data in the order of stimuli presentation!

disp('loading all stimuli...')
thisevent = 0; 
for ii = 1:size(key_events, 1)
    thisblock = double(key_events(ii, 1)); 
    if mod(ii, p.events) == 1
        thisevent = 1; 
    else
        thisevent = thisevent + 1; 
    end
    
    thisfile = dir(fullfile(dir_stim_all, char(key_events(ii, 6)))); 
    key_events(ii, 6) = thisfile.name; 
    thisfile = fullfile(thisfile.folder, thisfile.name);
    [tempAudio, fs_events(thisevent, thisblock)] = audioread(thisfile); 
    ad_events{thisevent, thisblock} = [tempAudio'; tempAudio']; 
end

disp('done!')

if ~all(all(fs_events(1) == fs_events))
    error('fs clear is not equal across all stimuli!')
end

fs = fs_events(1); 

% Clean up
dur_all = cellfun((@(x) length(x)/fs), ad_events); 
if any(dur_all > p.presTime)
    error('stim are too long!')
end

clear tempAudio

% Few more keys that depend on stim data!
key_stimEnd    = key_stimStart  + dur_all;
key_rxnEnd     = key_stimEnd    + p.rxnWindow; 

%% PTB
if DEBUG
    [wPtr, rect] = Screen('OpenWindow', 1, 185, [0 0 1280 720]);
else
    [wPtr, rect] = Screen('OpenWindow', 0, 185);
end

DrawFormattedText(wPtr, 'Please wait, preparing experiment...');
Screen('Flip', wPtr);
centerX = rect(3)/2;
centerY = rect(4)/2;
crossCoords = [-30, 30, 0, 0; 0, 0, -30, 30]; 
HideCursor(); 

pahandle = PsychPortAudio('Open', 5, [], [], fs); % 5 in lab, 18 at CCBBI
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
        
        WaitTill(firstPulse(rr) + p.epiTime - 1.5*p.TR); % edits Xiangrui
        t1(1, rr) = RTBox('WaitTR');
        
        %% Present audio stimuli
        for ev = 1:p.events
            PsychPortAudio('FillBuffer', pahandle, ... 
                ad_events{ev, rr});
            
            abs_stimStart_delta(ev, rr) = t1(ev, rr) + p.TR + key_jitter(ev, rr); 
            WaitTill(abs_stimStart_delta(ev, rr) - 0.1); 
            
            real_stimStart(ev, rr) = PsychPortAudio('Start', pahandle, ... 
                1, abs_stimStart_delta(ev, rr), 1);
            
            abs_rxnEnd_delta(ev, rr) = abs_stimStart_delta(ev, rr) + ... 
                dur_all(ev, rr) + p.rxnWindow; 
            
            RTBox('Clear'); 
            [real_respTime{ev, rr}, real_respKey{ev, rr}] = ... 
                RTBox(abs_rxnEnd_delta(ev, rr)); 
            
            % added delay because of scanner error
            WaitTill(t1(ev, rr) + p.eventTime - 0.5*p.TR); % can wait a while!
            if ev ~= 20
                t1(ev + 1, rr) = RTBox('WaitTR'); % does not need 'Clear'--TRs don't buffer
            end
        end

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
    OutputData_v3
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
OutputData_v3
disp('All done!')

%% how to counterbalance
function latinSquare=CreateLatinSquare(n,interleave)
% create a latinsquare with order: n
% by Niki ---2012/9/18
% the interleave style has a better order balance especially if n is
% large, to get a latinsquare with interleve style: CreatLatinSquare(n,1).
% by Niki---2013/3/18
% rewrite alghrithm to imporve speed
% by Niki 2014/7/1
% test:
if nargin==0
    CreateLatinSquare(10,1)
    CreateLatinSquare(9,0)
    CreateLatinSquare(10,0)
    CreateLatinSquare(9,1)
    return
end
if nargin<2
    interleave=1;
end
latinSquare=nan(n);
latinSquare(:,1)=1:n;
if interleave==1
    shiftsize=(.5-mod(1:n-1,2))/.5.*ceil((1:n-1)/2);
else
    shiftsize=n-1:-1:1;
end
for col=2:n
    latinSquare(:,col)=circshift((1:n)',shiftsize(col-1));
end
end