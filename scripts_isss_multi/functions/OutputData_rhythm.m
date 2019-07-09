%% OutputData.m
% Saves data into a txt and xlsx file. Run as part of the main experiment 
% script. Author -- Matt H

% CHANGELOG
% 08/08/17  Starter keeping changelog. --MH
% 08/08/17  Adjusted to allow multiple runs. --MH

%% Preallocate all variables
% Timing and subject response
stimDuration   = NaN(p.events, maxNumRuns); 
eventKeyPrint  = NaN(p.events, maxNumRuns); 
answerKeyPrint = NaN(p.events, maxNumRuns); 
reaction = NaN(p.events, maxNumRuns); 
response = NaN(p.events, maxNumRuns); 

% Xls file
headers = {'Jitter key', 'Actual jitter', ... 
        'Stim duration key', 'Actual stim duration', ...
        'Event duration', 'Event key', 'Answer key', ... 
        'Subj response', 'RT'}; 
data = cell(p.events + 1, 9); 
data(1,:)  = headers; 
sheetnames = {'run1', 'run2'}; 

% Txt file
dstring = '';
fstring = '';
for i = 1:p.events
    dstring = strcat(dstring, ' %d '); 
    fstring = strcat(fstring, ' %f ');
end

%% Convert relevant variables
% Convert to relative time, instead of system
runDuration  = runEnd - firstPulse; 
actualJitter = stimStart - eventStart; 
actualStimDuration  = stimEnd - stimStart; 
actualEventDuration = eventEnd - eventStart; 

% Convert variables to proper size. 
stimDuration(:, subj.firstRun:subj.lastRun) = ... 
    rawStimDur(eventKey(:, subj.firstRun:subj.lastRun)); 
%eventKeyPrint(:, subj.firstRun:subj.lastRun)     = eventKey; 
%answerKeyPrint(:, subj.firstRun:subj.lastRun)    = answerKey; 

% Convert keys from cells to vectors
for i = 1:p.events
    for j = 1:maxNumRuns
        if isempty(respKey{i, j})
            respKey{i, j}  = '0'; 
            respTime{i, j} = NaN; 
        end
        reaction(i, j) = respTime{i, j} - stimEnd(i, j); 
        response(i, j) = str2double(respKey{i, j});
    end
end

%% Path
cd(ResultsLoc)
if ~exist(subj.Num, 'file')
    mkdir(subj.Num); 
end
cd(subj.Num)

% Checks if files and variables already exists
while exist(ResultsXls, 'file') == 2
	ResultsXls = [ResultsXls(1:end-5), '_new', ResultsXls(end-4:end)]; 
end

while exist(ResultsTxt, 'file') == 2
	ResultsTxt = [ResultsTxt(1:end-4), '_new', ResultsTxt(end-3:end)]; 
end

while exist(Variables, 'file') == 2
	Variables = [Variables(1:end-4), '_new', Variables(end-3:end)]; 
end

%% Begin printing to txt file
fid = fopen(ResultsTxt, 'w');    
for run = subj.firstRun:subj.lastRun
    fprintf(fid, 'DATA FOR RUN ---------- \n');

    %% Timing data
    fprintf(fid, '# Timing data \n');
    
    fprintf(fid, 'Run started %6.2f after code started \n', ...
        firstPulse(run) - codeStart); 
    fprintf(fid, 'Run duration: %6.2f \n', runDuration(run));
    fprintf(fid, 'Expected run duration: %6.2f \n', p.runDuration); 

    jitterstring = ['Jitter key (msec): ', fstring, '\n']; 
    fprintf(fid, jitterstring, (jitterKey(:, run) * 1000)); 

    actualjitterstring = ['Actual jitter (msec): ', fstring, '\n']; 
    fprintf(fid, actualjitterstring, (actualJitter(:, run) * 1000)); 

    stimulistring = ['Stimuli duration key: ', fstring, '\n'];
    fprintf(fid, stimulistring, stimDuration(:, run));

    actualstimulistring = ['Actual stimuli duration: ', fstring, '\n']; 
    fprintf(fid, actualstimulistring, actualStimDuration(:, run)); 

    eventdurationstring = ['Event duration: ', fstring, '\n']; 
    fprintf(fid, eventdurationstring, actualEventDuration(:, run)); 

    fprintf(fid, 'Expected total duration: %f \n', p.eventTime); 

    %% Stimuli data
    fprintf(fid, '# Stimuli data \n'); 
    
    keystring = ['Event key: ', dstring, '\n'];
    fprintf(fid, keystring, eventKey(:, run));

    ansstring = ['Answer key: ', dstring, '\n'];
    fprintf(fid, ansstring, answerKey(:, run)); 

    %% Subject data
    fprintf(fid, '# Response data \n'); 

    respstring = ['Subject responses: ', dstring, '\n']; 
    fprintf(fid, respstring, response(:, run)); 

    reactionstring = ['Reaction time (msec): ', fstring, '\n'];
    fprintf(fid, reactionstring, (reaction(:, run) * 1000)); 

    %% Done printing!
    fprintf(fid, '\n'); 

    %% Prepare to print to xlsx file
    M    = horzcat(jitterKey(:, run), actualJitter(:, run), ...
        stimDuration(:, run), actualStimDuration(:, run), ... 
        actualEventDuration(:, run), ...
        eventKey(:, run), answerKey(:, run), ...
        response(:, run), reaction(:, run)); 
    
    for i = 1:p.events
        for j = 1:9
            data{i+1, j} = M(i, j); 
        end
    end

    %% Print to xlsx file
    warning off
    xlswrite(ResultsXls, data, sheetnames{run})
    warning on

end

%% All done!
fclose(fid); 
save(Variables); 