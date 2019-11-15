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
correct = false(p.events, p.runsMax); 
orsr = cell(p.events, p.runsMax); 

% Xls file
headers = {'SNR', 'Event key', 'Answer key', 'Subj response', 'RT', 'Correct', 'ORSR'}; 

%% Get or/sr events
% 1 -- OF -- mod4 gives 1
% 2 -- OM -- mod4 gives 2
% 3 -- SF -- mod4 gives 3
% 4 -- SM -- mod4 gives 0

key_stim_orsr = mod(key_stim, 4); 
stim_or = key_stim_orsr == 1 | key_stim_orsr == 2;
stim_sr = key_stim_orsr == 3 | key_stim_orsr == 0;

%% Saving relevant timing information
% Convert to relative time, instead of system
runDur = runEnd - firstPulse;

% real_jitter   = real_stimStart - real_eventStart; 
% real_stimDur  = real_stimEnd - real_stimStart; 
% real_eventDur = real_eventEnd - real_eventStart; 

% Convert keys from cells to vectors
for ii = 1:p.events
    for jj = 1:size(real_respKey, 2)
        if isempty(real_respKey{ii, jj})
            real_respTime{ii, jj} = NaN; 
        end
        
        reaction(ii, jj) = real_respTime{ii, jj} - real_stimStart(ii, jj); 
        if strcmp(real_respKey{ii, jj}, key_answer{ii, jj})
            correct(ii, jj) = true; 
        end
        
        if stim_or(ii, jj)
            orsr{ii, jj} = 'OR'; 
        elseif stim_sr(ii, jj)
            orsr{ii, jj} = 'SR'; 
        end
        
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

idx = 2; 
for rr = subj.firstRun:subj.lastRun
    for ev = 1:p.events
        data{idx, 1} = p.snr(whichBlock(rr)); 
        data{idx, 2} = key_stim(ev, rr); 
        data{idx, 3} = key_answer{ev, rr}; 
        data{idx, 4} = real_respKey{ev, rr}; 
        data{idx, 5} = reaction(ev, rr); 
        data{idx, 6} = correct(ev, rr); 
        data{idx, 7} = orsr{ev, rr}; 
        idx = idx + 1; 
    end
        
end

%% Save data
xlswrite(results_xlsx, data)

clear data ad_all_rms files_clear thisBlock
save(results_mat)
