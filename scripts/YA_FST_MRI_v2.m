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
% 01/24/20 -- New stimuli and counterbalancing
%   TODO: new stimuli loading!

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
    dlg_ans = {'TEST', '1', '1', '1'}; 
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
p.TR     = 1.000; 
p.epiNum = 8; 

% 40 minutes of scan time!

% Number of stimuli. Likely to change. 
% NumSpStim  = 192; 55
% NumStim    = 200;
% p.runsMax = 6; 

% Timing
p.runsMax = 9; % Now enough for 9?
p.events = 20; % Events per block, needs updating. 
% Each sentence (O/S) is presented nine times. 3 speeds x 3 SNRs
p.sentences  = 9; % how many sentence stimuli per block?
p.structures = 9; % how many sentence structures?

p.presTime = 4.000; % 4 seconds
p.jitter   = 1.000; % Jitter should be half of TR to prevent interference?

p.rxnWindow = 3.000;  % 3 seconds, should we expand this?

p.epiTime   = p.TR * p.epiNum;  % 8 seconds because I said so
p.eventTime = p.presTime + p.epiTime;
p.runDuration = p.epiTime + ...   % After first pulse
    p.eventTime * p.events + ...  % Each event
    p.eventTime;                  % After last acquisition

p.snr = [3 0]; % SNR between babble and speech
               % (clear, easy, hard)
%  1 % -1 % -3
warning('We have not discussed SNR yet!')

%% Paths
cd ..
dir_exp = pwd; 

dir_stim = fullfile(dir_exp, 'stimuli');
cfg.noisefile = fullfile(dir_stim, 'babble_4speakers.wav'); 
dir_stim_clear  = fullfile(dir_stim, 'YA_FST_v2_norm'); 

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

% %% Load stimuli
% files_clear = dir(fullfile(dir_stim_clear, '*.wav')); 
% fs_clear = zeros(1, length(files_clear)); 
% ad_clear = cell(length(files_clear), 1); % just clear speech
% 
% disp('loading clear stimuli...')
% for ii = 1:length(files_clear)
%     thisfile = fullfile(dir_stim_clear, files_clear(ii).name); 
%     [tempAudio, fs_clear(ii)] = audioread(thisfile); 
%     ad_clear{ii} = [tempAudio'; tempAudio']; 
% end
% 
% disp('done!')
% 
% if ~all(fs_clear(1) == fs_clear)
%     error('fs clear is not equal across all stimuli!')
% end
% 
% fs = fs_clear(1); 
% 
% % Add babble to clear stimuli
% disp('adding babble...')
% cd(dir_stim)
% 
% ad_babble_easy = cell(length(files_clear), 1); % high SNR = easy
% ad_babble_hard = cell(length(files_clear), 1); % low SNR = hard
% 
% cfg.prestim  = 0.210; % 0.25 - 0.04, leading silence in each stim
% cfg.poststim = 0; % each stim has ~0.3 lagging silence!
% cfg.fs = fs;
% 
% cfg.snrs = max(p.snr); 
% temp = jp_addnoise_hwk_mh_edits_CCBBI(dir_stim_clear, cfg); % DOES NOT SET RMS
% for ii = 1:length(ad_babble_easy)
%     ad_babble_easy{ii} = [temp{ii}', temp{ii}']';
% end
% 
% cfg.snrs = min(p.snr); 
% temp = jp_addnoise_hwk_mh_edits_CCBBI(dir_stim_clear, cfg); % DOES NOT SET RMS
% for ii = 1:length(ad_babble_hard)
%     ad_babble_hard{ii} = [temp{ii}', temp{ii}']';
% end
% 
% disp('done!')
% 
% % Noise 
% disp('making noise stimuli...')
% ad_noise = cell(p.runsMax, 1); % randomly chooses normal-speed stimuli
% noise_stim = Shuffle(2:3:length(ad_clear)); 
% noise_stim = noise_stim(1:p.runsMax); 
% 
% for ii = 1:length(noise_stim) % number of noise trials
%     tempAudio = jp_vocode_mh(ad_clear{noise_stim(ii)}(1, :), 1, fs); 
%     ad_noise{ii} = [tempAudio; tempAudio]; 
% end
% 
% disp('done!')
% 
% % Combine into one cell for RMS normalization
% ad_all = [ad_clear; ad_babble_easy; ad_babble_hard; ad_noise];
% 
% % Equalize RMS
% ad_all_rms = jp_equalizerms_mh(ad_all, 'verbose'); 
% 
% % Silence
% disp('loading silence...')
% tempAudio = audioread('silence01.wav'); 
% ad_all_rms{end+1} = [tempAudio, tempAudio]'; 
% disp('done!')
% 
% % Clean up
% dur_all = cellfun((@(x) length(x)/fs), ad_all_rms); 
% if any(dur_all > p.presTime)
%     error('stim are too long!')
% end
% 
% clear ad_all ad_clear ad_babble_easy ad_babble_hard tempAudio

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
% events_clear       =                         1:length(files_clear); 
% events_babble_easy =   length(files_clear) + 1:2*length(files_clear);
% events_babble_hard = 2*length(files_clear) + 1:3*length(files_clear);
% events_noise       = 3*length(files_clear) + 1:3*length(files_clear) + p.runsMax; 
% events_silence     = 3*length(files_clear) + p.runsMax + 1; 
% 
% events_slow = 1:3:3*length(files_clear);
% events_med  = 2:3:3*length(files_clear);
% events_fast = 3:3:3*length(files_clear);
% 
% events_sentence = repmat(repelem(1:9, 6), [1 3]);
% 
% if events_silence ~= length(ad_all_rms)
%     error('stimuli are not loaded correctly?')
% end
% 
% key_stim   = nan(p.events, p.runsMax); 
% key_answer = nan(p.events, p.runsMax); 
% babbleIdx = length(files_clear);

% cb_square = CreateLatinSquare(p.sentences); 
% cb_square = cb_square(randperm(p.sentences), :);
% key_sentence = (cb_square-1)*3 + 1; 
% key_speed = repmat(repelem([0 1 2]', 3), [2 p.runsMax]); 
% key_babble = repmat(repmat([0 babbleIdx 2*babbleIdx]', [6, 1]), [1, p.runsMax]); 
% 
% foo = key_sentence + key_speed + key_babble; 

% Thanks Hyun!!!
run = 1:p.runsMax;
syn = ["O","S"];
cla = ["","SNR3","SNR0"];
spr = ["0.75","1.0","1.25"];
ntrial = 9*2*3*3;

cond = strings(ntrial,6);
in = 0;
for rr=1:9
    for ss=1:2
        for cc=1:3
            for pp=1:3
                in = in+1;
                cond(in,1:4) = [run(rr), syn(ss), cla(cc), spr(pp)];
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

for rr = 1:p.runsMax
    thesesent = key_sentence(:, rr); 
    
    % we need to assign each babble and speech rate across runs
    % we need to make sure each assignment only happens once!
    % M/F and OR/SR counterbalancing is already addressed. 
    
    % 18 sentences/block, therefore:
    % 6 slow   -- add 0
    % 6 medium -- add 1
    % 6 fast   -- add 2
    
    % 18 sentences/block, therefore:
    % 6 clear    -- add 0
    % 6 easy SNR -- add babbleIdx
    % 6 hard SNR -- add 2*babbleIdx
    
    clearbabble = repmat([0 length(files_clear)]', [8 1]); 
    noi_sil = [events_noise(:, rr); repmat(events_silence, [4 1])]; 
    
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
while 1 % use this for loop to ensure jitter is longer than 0.001
    key_jitter = p.jitter * rand(p.events, p.runsMax); 
    if all(all(key_jitter > 0.001))
        break
    end
end

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
    [wPtr, rect] = Screen('OpenWindow', 0, 185);
end

DrawFormattedText(wPtr, 'Please wait, preparing experiment...');
Screen('Flip', wPtr);
centerX = rect(3)/2;
centerY = rect(4)/2;
crossCoords = [-30, 30, 0, 0; 0, 0, -30, 30]; 
HideCursor(); 

pahandle = PsychPortAudio('Open', 5, [], [], fs); % 18 at scanner?
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
                ad_all_rms{key_stim(ev, rr)});
            
            abs_stimStart_delta(ev, rr) = t1(ev, rr) + p.TR + key_jitter(ev, rr); 
            WaitTill(abs_stimStart_delta(ev, rr) - 0.1); 
            
            real_stimStart(ev, rr) = PsychPortAudio('Start', pahandle, ... 
                1, abs_stimStart_delta(ev, rr), 1);
            
            abs_rxnEnd_delta(ev, rr) = abs_stimStart_delta(ev, rr) + ... 
                key_stimDur(ev, rr) + p.rxnWindow; 
            
            RTBox('Clear'); 
            [real_respTime{ev, rr}, real_respKey{ev, rr}] = ... 
                RTBox(abs_rxnEnd_delta(ev, rr)); 
            
            % added delay because of scanner error
            WaitTill(t1(ev, rr) + p.eventTime - 0.5*p.TR); % can wait a while!
            t1(ev + 1, rr) = RTBox('WaitTR'); % does not need 'Clear'--TRs don't buffer
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
    OutputData_v2_edits_XL
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
OutputData_v2_edits_XL
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