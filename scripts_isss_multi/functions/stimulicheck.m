%% Stimuli Check
% Makes sure stimuli are counterbalanced. Checks how many stimuli have a
% male actor, a female actor, a subject-, and a object-relative clause.
% Then checks how many are MO, MS, FO, and FS. Author -- Matt

% Set variables for debugging
% ns = NumSpeechStimuli; 
% eventkey = eventKey; 

function stimulicheck(ns, eventkey)
% Set up keys and preallocate when possible
speech  = eventkey(eventkey <= ns); 
s = length(speech); 

obj     = []; 
subj    = []; 

objKey  = sort(horzcat(1:4:ns, 2:4:ns)); 
subjKey = sort(horzcat(3:4:ns, 4:4:ns)); 

% Pull out which stimuli are male/fem or obj/subj
male = speech(logical(mod(speech + 1, 2))); 
fem  = speech(logical(mod(speech, 2))); 

for i = 1:length(speech)
    if     find(speech(i) == objKey)
        obj  = horzcat(obj, speech(i)); 
    elseif find(speech(i) == subjKey) 
        subj = horzcat(subj, speech(i)); 
    end
end

% Prepare to pull out OM, SM, OF, SF
objMale  = [];
subjMale = [];
objFem   = [];
subjFem  = [];

% Find OM, SM, OF, SF stimuli
for i = 1:length(male)
    if find(male(i) == objKey)
        objMale = horzcat(objMale, male(i)); 
    elseif find(male(i) == subjKey)
        subjMale = horzcat(subjMale, male(i)); 
    end
end

for i = 1:length(fem)
    if find(fem(i) == objKey)
        objFem = horzcat(objFem, fem(i)); 
    elseif find(fem(i) == subjKey)
        subjFem = horzcat(subjFem, fem(i)); 
    end
end

% Double check that we have correct number of stimuli
numMale = length(male); % Should be 1/2
numFem  = length(fem) ; % Should be 1/2
numObj  = length(obj) ; % Should be 1/2
numSubj = length(subj); % Should be 1/2

numObjMale  = length(objMale) ; % Should be 1/4
numObjFem   = length(objFem)  ; % Should be 1/4
numSubjMale = length(subjMale); % Should be 1/4
numSubjFem  = length(subjFem) ; % Should be 1/4

if (numMale     ~= s/2 || numFem     ~= s/2 || numObj      ~= s/2 || numSubj    ~= s/2 || ... 
    numObjMale  ~= s/4 || numObjFem  ~= s/4 || numSubjMale ~= s/4 || numSubjFem ~= s/4)
    error('Check your stimuli list, something is wrong')
end

end