% LoadStimuli.m
function [audiodata, samplingrate, duration, jitkey, eventkey, anskey, rhythmkey] = ... 
    LoadStimAndKeys_rhythm(fileloc, numEvents, firstRun, lastRun, maxRuns) 
% Loads audio data from \stimuli into two cell arrays; one for the
% waveform (audiodata) and one for the sampling rate (samplingrate). Now 
% generates the answer, jitter, and event keys at the same time. 
% Author - Matt H

% CHANGELOG
% 07/17/17  Started keeping changelog on file. -- MH
% 07/17/17  Changed to v2, allowing for multiple runs. -- MH
% 08/09/17  Made sure it was pretty. -- MH

% Set variables while debugging
%  fileloc   = StimuliLoc; 
%  numEvents = p.events; 
%  firstRun  = subj.firstRun; 
%  lastRun   = subj.lastRun; 
%  maxRuns   = 2; 

%% Preparing to load stimuli
% Load file names
cd(fileloc)
files = dir('*.wav'); 

% Preallocate variables
anskey = NaN(numEvents, maxRuns); 

audiodata    = cell(1, length(files));
samplingrate = cell(1, length(files));

%% Load ALL audio stimuli
for i = 1:length(files)
    [tempAudio, tempFs] = audioread(files(i).name);
    audiodata{i} = [tempAudio'; tempAudio'];
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
% eventkey -- In what order will stimuli be presented?
% Carefully choose which stimuli to present to ensure stimuli are
% counterbalanced. 
eventkey = NaN(numEvents, maxRuns); 
for i = firstRun:lastRun
    eventkey(:, i) = Shuffle(vertcat( ... 
        1 * ones(4, 1), ... % Simple long
        2 * ones(4, 1), ... % Complex long
        3 * ones(4, 1), ... % Simple short
        4 * ones(4, 1), ... % Complex short
        5, 6, 7, 8 ...      % Oddball, oddball, silence, silence
        )); 
end

% jitkey -- How much was the jitter before stimulus presentation? This
% depends on which stimuli is being presented, as durations vary across one
% condition. 
shortstimduration = 1.700;
longstimduration  = 3.000;
jitkey = NaN(numEvents, 2); 
for i = firstRun:lastRun
    for j = 1:numEvents
        if ~isempty(find(eventkey(j, i)     == [1 2 5], 1)) % Long conditions
            jitkey(j, i) = (4-longstimduration)*rand(1); 
        elseif ~isempty(find(eventkey(j, i) == [3 4 6], 1)) % Short conditions
            jitkey(j, i) = (4-shortstimduration)*rand(1); 
        elseif ~isempty(find(eventkey(j, i) == [7 8], 1))   % Silent conditions
            jitkey(j, i) = 2*rand(1); 
        end
    end
end

% rhythmkey -- Which rhythms were presented to the subject? 
rhythmkey = sort(eventkey);

% anskey -- What should have subjects responded with?
for i = firstRun:lastRun
    for j = 1:numEvents
        if ~isempty(find(eventkey(j, i) == [5 6], 1)) % events that are odd
            anskey(j, i) = 1; 
        else % events that are not odd
            anskey(j, i) = 0; 
        end
    end
end

end