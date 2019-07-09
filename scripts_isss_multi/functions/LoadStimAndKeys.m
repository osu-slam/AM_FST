% LoadStimuli.m
function [audiodata, samplingrate, duration, jitkey, eventkey, anskey, speechkey] = ... 
    LoadStimAndKeys(fileloc, numEvents, firstRun, lastRun, numSpeech, maxRuns, train) 
% Loads audio data from stimuli.wav into two cell arrays; one for the
% waveform (audiodata) and one for the sampling rate (samplingrate). Now 
% generates the answer, jitter, and event keys at the same time. 
% Keep all stimuli in fileloc or training.Author - Matt H

% CHANGELOG
% 07/12/17 Started keeping record of changelog on file. -- MH
% 07/19/17 Changed to v6, allowing for multiple runs and training. -- MH
% 08/09/17 Made it pretty for experiments. -- MH

%% Set variables while debugging
%  fileloc   = StimuliLoc; 
%  numEvents = p.events; 
%  firstRun  = subj.firstRun;
%  lastRun   = subj.lastRun;
%  numSpeech = NumSpStim; 
%  maxRuns   = maxNumRuns; 
%  train     = Training; 

%% Preparing to load stimuli
% Load file names
if train
    cd([fileloc filesep 'training'])
else
    cd(fileloc)
end
files = dir('*.wav'); 

% Preallocate variables
anskey = NaN(numEvents, maxRuns); 

audiodata    = cell(1, length(files));
samplingrate = cell(1, length(files));

%% Load ALL audio stimuli
for i = 1:length(files)
    [tempAudio, tempFs] = audioread(files(i).name);
    audiodata{i}    = [tempAudio'; tempAudio']; % Convert mono to stereo
    samplingrate{i} = tempFs;
end

% Check samplingrate is same across files (assumption made in main code)
for i = 2:length(samplingrate)
    if samplingrate{i} ~= samplingrate{i-1}
        error('Your sampling rates are not all the same. Stimuli will not play correctly.')
    end
end

% Load duration of files
clear info
info(length(audiodata)) = audioinfo(files(end).name); % Preallocate struct
for i = 1:length(files)
    info(i) = audioinfo(files(i).name); 
end

duration = NaN(1, length(info)); 
for i = 1:length(duration)
    duration(i) = info(i).Duration; 
end

%% Make keys
% jitkey -- How much is the silent period jittered by?
if train
    jitkey = 0.9 + rand(numEvents, 1);
else
    jitkey = NaN(numEvents, 6);
    for i = firstRun:lastRun
        jitkey(:, i) = 0.9 + rand(numEvents, 1); % Add 1 because stimuli are short-ish
    end
end

% speechkey -- Which speech stimuli should we use this run?
% Carefully choose which stimuli to present to ensure stimuli are
% counterbalanced. 
% eventkey -- In what order will stimuli be presented?
if train
    speechkey = [1; 2; 3; 4]; 
    eventkey  = Shuffle([speechkey; 5]); 
else
    randomstim   = NaN(numEvents/2, maxRuns); 
    sentence     = NaN(numEvents/2, maxRuns); 
    noisesilence = NaN(numEvents/2, maxRuns); 
    
    for i = firstRun:lastRun
        randomstim(:, i) = Shuffle(vertcat( ... 
            0 * ones(numEvents/8, 1), ... % Half of events are speech stim
            1 * ones(numEvents/8, 1), ... % Four conditions of speech means...
            2 * ones(numEvents/8, 1), ... % events/(2*4) = events/8
            3 * ones(numEvents/8, 1) ... 
            )); 
        sentence(:, i) = (((i-1)*32)+1:4:i*numEvents*2)'; 
        noisesilence(:, i) = (numSpeech+1:length(audiodata))'; 
    end
    
    speechkey = sentence + randomstim;
    eventkey  = Shuffle(vertcat(speechkey, noisesilence)); 
end

% anskey -- What should have subjects responded with?
for i = 1:numEvents
    for j = firstRun:lastRun
        if     eventkey(i, j) > numSpeech + 4 % Silence
            anskey(i, j) = 0; 
        elseif eventkey(i, j) > numSpeech     % Noise
            anskey(i, j) = 3; 
        elseif mod(eventkey(i, j), 2) == 0    % Male
            anskey(i, j) = 2; 
        elseif mod(eventkey(i, j), 2) == 1    % Female
            anskey(i, j) = 1; 
        end
    end
end

end