function [Cfg] = screen_parameters
%% PARAMETERS FOR FMD-DUAL-TASK
% [ Cfg ] = screen_parameters
% Set up all experimental parameters for screen and stimuli

Cfg.aux_buffer = 1;

%% INITIALISE SCREEN

% Window size (blank is full screen)
Cfg.WinSize = [];
% Cfg.WinSize = [10 10 850 750];

% Select screen
Cfg.screens = Screen('Screens');
Cfg.screenNumber = min(Cfg.screens); % Main screen
% Cfg.screenNumber = max(Cfg.screens); % Laptop

% Define colours
Cfg.white = WhiteIndex(Cfg.screenNumber);
Cfg.black = BlackIndex(Cfg.screenNumber);
Cfg.gray = round((Cfg.white + Cfg.black)/2);

if round(Cfg.gray) == Cfg.white
    Cfg.gray = Cfg.black;
end

Cfg.WinColor = Cfg.gray;

% 2017-02-23 - JULIAN
% Experienced syncerrors on my laptop, included synctest skip input to
% debug code but have commented out for piloting & experiment.

% skip_sync = input('Skip sync test? (y/n) \n','s');
%
% if strcmp(skip_sync,'y')
%     Screen('Preference', 'SkipSyncTests', 1);
% end

%('OpenWin', WinPtr, WinColour, WinRect, PixelSize, AuxBuffers, Stereo)
[Cfg.windowPtr, Cfg.windowRect] = Screen('OpenWindow', ...
    Cfg.screenNumber, Cfg.WinColor, Cfg.WinSize, [], 2, 0);

% Find window size
[Cfg.width, Cfg.height] = Screen('WindowSize', Cfg.windowPtr);

% Define center X & Y
[Cfg.xCentre , Cfg.yCentre] = RectCenter(Cfg.windowRect);

% Font
Screen('TextFont', Cfg.windowPtr, 'Arial');

% Text size
Screen('TextSize', Cfg.windowPtr, 20);

Cfg.computer = Screen('Computer');
Cfg.version = Screen('Version');

% This should be '60', will be a possible source of crashes if not
Cfg.FrameRate = Screen('NominalFrameRate', Cfg.windowPtr); 

[x,y] = Screen('DisplaySize',Cfg.windowPtr);
Cfg.xDimCm = x/10;
Cfg.yDimCm = y/10;

Cfg.distanceCm = 60;    % Defined by Sherman, Seth, Barrett, Kanai (2015)

%DEG VISUAL ANGLE FOR SCREEN
% 2018-03-27 Confirmed that this operation is correct, has unintuitive order of operations that simultanously convert 
% to radians & degrees. Consider 'atand' function for future versions.
Cfg.visualAngleDegX = atan(Cfg.xDimCm/(2*Cfg.distanceCm))/pi*180*2;
Cfg.visualAngleDegY = atan(Cfg.yDimCm/(2*Cfg.distanceCm))/pi*180*2;

%DEG VISUAL ANGLE PER PIXEL
Cfg.visualAnglePixelPerDegX = Cfg.width/Cfg.visualAngleDegX;
Cfg.visualAnglePixelPerDegY = Cfg.height/Cfg.visualAngleDegY;

Cfg.pixelsPerDegree = mean([Cfg.visualAnglePixelPerDegX Cfg.visualAnglePixelPerDegY]); % Usually the mean is reported in papers

%% SET UP PARAMETERS FOR FIXATION CROSS

%Set colour, width, length etc.
Cfg.crossColour = Cfg.black;  %255 = white
Cfg.crossLength = Cfg.pixelsPerDegree * 0.38; % Occupied 0.38 degrees in Sherman
Cfg.crossWidth = 5;

%Set start and end points of lines
crossLines = [-Cfg.crossLength, 0; Cfg.crossLength, 0; 0, -Cfg.crossLength; 0, Cfg.crossLength];
Cfg.crossLines = crossLines';

%% LOAD INSTRUCTION IMAGES & MAKE TEXTURES

image_location = '../addons/instructions/';

train_peri = imread([image_location '01-peri-instruct.png']);
train_cent = imread([image_location '02-cent-instruct.png']);
train_dual = imread([image_location '03-dual-instruct.png']);

part_2 = imread([image_location '04-part-2.png']);

divert_25 = imread([image_location 'divert-25.png']);
divert_50 = imread([image_location 'divert-50.png']);
divert_75 = imread([image_location 'divert-75.png']);
full_25 = imread([image_location 'full-25.png']);
full_50 = imread([image_location 'full-50.png']);
full_75 = imread([image_location 'full-75.png']);

intermish = imread([image_location 'intermission.png']);

Cfg.instruct.train_peri = Screen('MakeTexture',Cfg.windowPtr,train_peri);
Cfg.instruct.train_cent = Screen('MakeTexture',Cfg.windowPtr,train_cent);
Cfg.instruct.train_dual = Screen('MakeTexture',Cfg.windowPtr,train_dual);

Cfg.instruct.part_2 = Screen('MakeTexture',Cfg.windowPtr,part_2);

Cfg.instruct.divert_25 = Screen('MakeTexture',Cfg.windowPtr,divert_25);
Cfg.instruct.divert_50 = Screen('MakeTexture',Cfg.windowPtr,divert_50);
Cfg.instruct.divert_75 = Screen('MakeTexture',Cfg.windowPtr,divert_75);
Cfg.instruct.full_25 = Screen('MakeTexture',Cfg.windowPtr,full_25);
Cfg.instruct.full_50 = Screen('MakeTexture',Cfg.windowPtr,full_50);
Cfg.instruct.full_75 = Screen('MakeTexture',Cfg.windowPtr,full_75);

Cfg.instruct.intermish = Screen('MakeTexture',Cfg.windowPtr,intermish);

end
