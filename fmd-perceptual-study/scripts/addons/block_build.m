function [ built_block ] = block_build( Cfg,Subj,expectation_condition )
%% BLOCK_BUILD block_build(Cfg,Subj,ExpStr)
% Constructs basic setup for a block of 'Subj.nTrial' trials
%
% Cfg = Screen parameters struct
% Subj = Subject parameters struct
%
% ExpStr = string for expectation (i.e 'lo_full' or 'med' or 'hi_divert')
% Examines for presence of 'lo', 'med', and 'hi' to determine probability
% of peripheral target presence (25, 50 or 75%)

%% SHUFFLE LETTER PRESENCE
% Define letter presence style (per trial or counterbalanced)
% We will counterbalance to ensure 50% proportion of each letter condition

letter_counter = 1;

if letter_counter == 1 % counterbalanced over block
    letter_presence = zeros(1,Subj.nTrials);
    letter_presence(1:Subj.nTrials/2) = 0;
    letter_presence(Subj.nTrials/2+1:end) = 1;
    letter_presence = Shuffle(letter_presence);
elseif letter_counter == 0 % per trial 50% chance
    letter_presence = randi(2,[1,Subj.nTrials]);
end

%% DEFINE FRAMERATE CONVERSION
% Cfg.FrameRate is defined in screen_parameters.m
% Usually 60Hz which makes frame conversion (milliseconds * 0.06)
frame_conversion = Cfg.FrameRate/1000;

%% FIXATION CROSS DURATION
% Random value from 500ms to 1500ms
fixation_temp = randi(1000,[1,Subj.nTrials]);
fixation_temp = fixation_temp + 500;

%% LETTER ROTATION ANGLES
% Random rotation angle (0 to 359 degrees) for each of the 4 letters
ang_temp = randi(360,[4,Subj.nTrials]);
ang_temp = ang_temp - 1;

%% SPECIFY EXPECTATION LEVEL
% Sherman et al. = 25% 50% 75% probability of target presence
% Search for 'lo','med', & 'hi' strings
if any(strfind(expectation_condition,'lo'))
    expectation = .25;
elseif any(strfind(expectation_condition,'med'))
    expectation = .5;
elseif any(strfind(expectation_condition,'hi'))
    expectation = .75;
else
    disp('Expectation unspecified, defaulting to 50%')
    expectation = .5;
end

%% SPECIFY PROPORTION OF PERIPHERAL STIMULUS BASED ON EXPECTATION
% Build variable containing specified proportions
% Assumes two stimulus conditions: presence vs. absence | male vs. female | Red-Green vs. Green-Red | etc

% Specified expectation (target presence in Sherman et al. version)
presence = ones(1,round(Subj.nTrials * expectation));

% Unspecified expectation (target absence in Sherman et al. version)
absence = zeros(1,round(Subj.nTrials * (1-expectation)));

% Concatenate them together and shuffle the order
peripheral_presence = Shuffle(horzcat(presence,absence));

%% BUILD TRIALS
% Clear TR ready for build
TR = [];

for tr = 1 : Subj.nTrials
    
    % The condition of the letter task
    % In our study, 4 conditions: all Ts, all Ls, Ts+L, Ls+T
    % Sherman task is different, either all Ls or Ls+T
    
    TR(tr).targetTrialType = letter_presence(tr); %#ok<*AGROW>
    
    % Sanity check
    if TR(tr).targetTrialType == 0
        TR(tr).letter_condition = 'Ls';
    elseif TR(tr).targetTrialType == 1
        TR(tr).letter_condition = 'Ls+T';
    end
    
    %% MAIN TRIAL SETUP
    
    % Initial values for SOAs (in frames)
    % These will be assigned by QUEST during training
    % In Sherman cSOA mean=254ms SD=75ms
    % pSOA fixed @ 388ms "with gradual onset and offset"
    TR(tr).cSOA = [];
    TR(tr).pSOA = round(388*frame_conversion);
    TR(tr).true_cSOA = [];
    TR(tr).true_pSOA = 388*frame_conversion;
    
    % Peripheral contrast QUEST
    % Separate for full and diverted attention conditions
    % TR(tr).pCON_full = [];
    % TR(tr).pCON_divert = [];
    
    % In frames: @ 60Hz 500ms is 30 frames
    % TR(tr).screenInterval = 30;
    
    % In frames
    TR(tr).fixation_duration = round(fixation_temp(tr)*frame_conversion);
    
    % TR(tr).intertrialInt = 50; % In frames! 60hz... therefore, 60 = 1sec
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%  CENTRAL TASK   %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% CENTRAL TASK
    
    % Central letter displacement (1.43 degrees x 1.43 degrees: visual angle)
    % 35 pixels in our original study, seems larger in Sherman
    TR(tr).letterHeight = round(Cfg.pixelsPerDegree * 1.43);
    
    TR(tr).letter_angles = ang_temp(:,tr)'; % 4 angles @ 0:359 degrees
    
    % Find centre of screen for central letter presentation
    centreLeft = Cfg.xCentre - TR(tr).letterHeight;
    centreTop = Cfg.yCentre - TR(tr).letterHeight;
    centreRight = Cfg.xCentre + TR(tr).letterHeight;
    centreBottom = Cfg.yCentre + TR(tr).letterHeight;
    
    % Centre of letter texture (x coordinates)
    TR(tr).centreLetterRect(:,1) = [centreLeft, Cfg.xCentre, centreRight, Cfg.xCentre];
    
    % Centre of letter texture (y coordinates)
    TR(tr).centreLetterRect(:,2) = [Cfg.yCentre, centreTop, Cfg.yCentre, centreBottom];
    
    % If this is a "different" letter trial, randomly select which location
    % will be different
    % LURD = Left Up Right Down
    if strcmp(TR(tr).letter_condition,'Ls+T')
        TR(tr).replace_LURD = randi(4,1);
    else
        TR(tr).replace_LURD = NaN;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%  PERIPHERAL TASK  %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% PERIPHERAL TASK
    
    % Sherman et al have a fairly complex display procedure
    % Define peripheral rectangle (8x10 degree of visual angle) on which
    % the stimulus will be displayed
    %
    % In the Sherman copy this will contain a gabor patch
    % We might add a face-gender discrimination but the location will be
    % the same
    
    % Include jitter from Sherman study (0.66 x 1.24 degrees)
    % xJit = 10/2 + (randi(124*2)-124)/100;
    % yJit = 8/2 + (randi(66*2)-66)/100;
    
    % Baseline position of (25.2 x 21.08 degrees)
    xJit = 25.2/3 + (randi(124*2)-124)/100;
    yJit = 21.08/3 + (randi(66*2)-66)/100;
    
    xptsRectPeriph = [-xJit xJit -xJit xJit] .* Cfg.pixelsPerDegree + Cfg.xCentre;
    yptsRectPeriph = [-yJit -yJit yJit yJit] .* Cfg.pixelsPerDegree + Cfg.yCentre;
    
    centrePeriRect = vertcat(xptsRectPeriph,yptsRectPeriph);
    
    % Randomly select one point, allocate this as the position for the trial
    usePoint = randi(4,1);
    TR(tr).centrePeriRect = round(centrePeriRect(:,usePoint)');
    
    % Set size of the item to be displayed in the periphery
    % Face in Matthews et al., gabor patch for Sherman et al.
    imageHeightPeriph = 4 * Cfg.pixelsPerDegree; % Approximately 2 degrees visible
    TR(tr).imageHeightPeriph = round(imageHeightPeriph);
    
    %% PERIPHERAL STIMULUS
    % Define whether face-gender (Matthews et al.) or gabor task (ala Sherman et al.)
    
    % Select peripheral condition (2 or 1: Present|Absent or Male|Female)
    % Expectation is baked in above
    TR(tr).peripheral_presence = peripheral_presence(tr);
    
    switch Subj.peripheral_task
        case 'face'
            % Set gender for main peripheral task
            if TR(tr).peripheral_presence == 0
                TR(tr).peripheral_type = 'f';
            elseif TR(tr).peripheral_presence == 1
                TR(tr).peripheral_type = 'm';
            end
            
            %Set face image number (here from 1 - 65) randomly
            picMainGen = randperm(65);
            TR(tr).peripheral_tex = picMainGen(33);
            
            %Set gender for MASK to male/female 50/50 ratio (as above) using
            %scrambled faces
            mask_type = rand(1);
            
            if mask_type <= 0.5
                TR(tr).peripheral_mask = 'f_sc';
            else
                TR(tr).peripheral_mask = 'm_sc';
            end
            
            %Set face image number (here from 1 - 20) randomly
            picMaskGen = randperm(20);
            TR(tr).mask_tex = picMaskGen(8);
            
        case 'gabor'
            
            % Set presence of gabor
            if TR(tr).peripheral_presence == 0
                TR(tr).peripheral_type = 'absent';
            elseif TR(tr).peripheral_presence == 1
                TR(tr).peripheral_type = 'present';
            end
            
            % Mebbe define some other stuff here
            
        otherwise
            disp('Peripheral task undefined, defaulting to gabor task')
            
            % Redefine peripheral_task to 'gabor' and restart
            Subj.peripheral_task = 'gabor';
            return
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%  OUTPUT/SAVE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    TR(tr).mouseResponsesCent = [];
    TR(tr).central_confid = [];
    TR(tr).central_RT = []; % Time taken to respond
    TR(tr).mouseResponsesPer = [];
    TR(tr).peripheral_confid = [];
    TR(tr).peripheral_RT = []; % Time taken to respond
    
end
%% RANDOMISE AND INPUT INTO STRUCT

% Randomize trials
random_allocation = randperm(Subj.nTrials);
TR = TR(random_allocation);

% Input into this block build
built_block = TR;

end

