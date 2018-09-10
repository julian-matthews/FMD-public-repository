%% flash_screen(window, flash_repeats, flash_time, flash_string)
% 2015-08-27 - Julian
% Flashes screen with text at a set rate a few times.
% Annoying enough to discourage participants from making errors.
%
% Requires input of:
% window - window for Psychtoolbox screen
% flash_repeats - number of times to flash from one to other (ie. 5)
% flash_time - time that screen remains one colour (ie 80ms)
% flash_string - message in centre of screen
%
% UPDATE: 2017-02-08 - simplifications for FMD dual-task

% Uses GetMouse.m to register mouse click

function flash_screen(window, flash_repeats, flash_time, flash_string)

% Define colours
black = [0 0 0];
white = [255 255 255];
red = [255 0 0];
yellow = [255 255 0];

Screen('TextSize', window, 40);

% Add some pointless variability to the flashing colour
coin_toss = randi(2);

if coin_toss == 1
    flash_colour = red;
else
    flash_colour = yellow;
end

% Cycle for n=flash_repeats
for repeats = 1:flash_repeats
    
    flash_remaining = flash_time;
    flashSecs = GetSecs;
    
    while flash_remaining > 0
        
        Screen('FillRect', window, flash_colour);
        DrawFormattedText(window,flash_string,'center','center',black);
        Screen('Flip', window, [], []);
        
        time_elapsed = GetSecs - flashSecs;
        flash_remaining = flash_time - time_elapsed;
        
    end
    
    flash_remaining = flash_time;
    flashSecs = GetSecs;
    
    while flash_remaining > 0
        
        Screen('FillRect', window, black);
        [~,y_cent,~] = DrawFormattedText(window,flash_string,'center','center',white);
        Screen('Flip', window, [], []);
        
        time_elapsed = GetSecs - flashSecs;
        flash_remaining = flash_time - time_elapsed;
        
    end
    
end

Screen('FillRect', window, black);
Screen('TextSize', window, 30);
DrawFormattedText(window,flash_string,'center','center',white);
DrawFormattedText(window,'<< click to continue >>','center',[y_cent + 50],white);
Screen('Flip', window, [], []);

WaitSecs(.5);

while (1)
    [~,~,buttons] = GetMouse(window);
    if buttons(1) || KbCheck
        break;
    end
end

end