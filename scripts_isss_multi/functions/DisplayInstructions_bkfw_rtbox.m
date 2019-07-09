% DisplayInstructions_bkfw_rtbox.m
% Subfunction called to display instructions on screen from a .txt file.
% Now runs with RTBox commands as of 6/30/17. 
% Keep the .txt file in the same folder as this script. Author - Matt H
function DisplayInstructions_bkfw_rtbox(filename, windowPointer, rtboxloc)

% Set properties of display
font = 'Cambria'; % The best font. 
size = 40;

% Read instructions into cell array
fid = fopen(filename, 'r'); 
i = 1; 
instructions{i} = fgetl(fid); % Creates a cell array with each line of 
                              % instructions as an element. 
while ischar(instructions{i})
    i = i + 1;
    instructions{i} = fgetl(fid); 
end
fclose(fid); 

% Display instructions
i = 1; 
while i < length(instructions) % Last element is empty
    WaitSecs(0.5);
    if i == 0 % Prevent script from quitting on wrong button press
        DrawFormattedText(windowPointer, 'Wrong button! Try the other one.', 'center', 'center');
        Screen('Flip',windowPointer);
    elseif i == 1 % Welcome screen
        Screen('TextFont', windowPointer, font); 
        Screen('TextSize', windowPointer, size);
        Screen('TextStyle', windowPointer, 2);
        DrawFormattedText(windowPointer, instructions{i}, 'center', 'center');
        Screen('Flip',windowPointer);
    elseif i == length(instructions)-1 % Final prompt
        Screen('TextFont', windowPointer, font); 
        Screen('TextSize', windowPointer, size);
        Screen('TextStyle', windowPointer, 1);
        DrawFormattedText(windowPointer, instructions{i}, 'center', 'center');
        Screen('Flip',windowPointer);
    else % Instructions
        Screen('TextFont', windowPointer, font); 
        Screen('TextSize', windowPointer, size);
        Screen('TextStyle', windowPointer, 0);
        DrawFormattedText(windowPointer, instructions{i}, 'center', 'center');
        Screen('Flip',windowPointer);
    end

% Collect subject response
    cd(rtboxloc)    
while 1
    [isDown, key] = RTBox(inf); 
    if isDown
        if strcmp(key, '2')
            i = i + 1; 
            break
        elseif strcmp(key, '1')
            i = i - 1; 
            break
        end
    end
end

end
    