%% OutputData.m
% Saves data into a txt and xlsx file. Run as part of the main experiment 
% script. Author -- Matt H

% CHANGELOG
% 08/08/17  Started keeping changelog. --MH
% 08/08/17  Now allows for mulitple runs. --MH
% 08/10/17  Ready for subject 3. --MH
% 10/10/19  Updating for AM_FST. --MH

%% Preallocate all results_mat
% Timing and subject response
reaction = NaN(p.events, p.runsMax); 
response = NaN(p.events, p.runsMax); 

% Xls file
headers = {'BLOCK', 'Jitter key', 'Actual jitter', ... 
        'Stim duration key', 'Actual stim duration', ...
        'Actual Event duration', 'Event key', 'Answer key', ... 
        'Subj response', 'RT'}; 

%% Saving relevant timing information
% Convert to relative time, instead of system
runDur = runEnd - firstPulse; 

real_jitter   = real_stimStart - real_eventStart; 
real_stimDur  = real_stimEnd - real_stimStart; 
real_eventDur = real_eventEnd - real_eventStart; 

% Convert keys from cells to vectors
for ii = 1:p.events
    for jj = 1:size(real_respKey, 2)
        if isempty(real_respKey{ii, jj})
            real_respTime{ii, jj} = NaN; 
        end
        reaction(ii, jj) = real_respTime{ii, jj} - real_stimEnd(ii, jj); 
        response(ii, jj) = str2double(real_respKey{ii, jj});
    end
end

%% Checks if files already exists to prevent overwrite
cd(dir_results)

while exist(results_xlsx, 'file') == 2
	results_xlsx = [results_xlsx(1:end-5), '_new', results_xlsx(end-4:end)]; 
end

while exist(results_mat, 'file') == 2
	results_mat = [results_mat(1:end-4), '_new', results_mat(end-3:end)]; 
end

%% Begin printing to xlsx file
% Rather than do many sheets, let's just do one long table. 
data = cell(p.events*p.runsMax + 1, length(headers)); 
data(1,:) = headers; 

idx = 1; 
for rr = subj.firstRun:subj.lastRun
    M = horzcat(repmat(rr, [p.events 1]), ...
        key_jitter(:,rr), real_jitter(:,rr), ...
        key_stimDur(:,rr), real_stimDur(:,rr), ... 
        real_eventDur(:,rr), ...
        key_stim(:,rr), key_answer(:,rr), ...
        response(:,rr), reaction(:,rr)); 

    for ii = 1:p.events
        for jj = 1:length(headers)
            data{idx+1, jj} = M(ii, jj); 
        end
        
        idx = idx + 1; 
    end
    
end

%% Save data
xlswrite(results_xlsx, data)

clear data ad_all_rms
save(results_mat)
