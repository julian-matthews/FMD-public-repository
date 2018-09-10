%% FMD-DUAL-TASK
% Created 2016-11-21: Julian
%
% Dual-task + expectation experiment adapted from Sherman, Seth, Barrett &
% Kanai (2015). Central letter and peripheral "grill" tasks categorising
% targets present or absent + two-level confidence judgement.
%
% Expectation of "grill" presence manipulated per block.
%
% Updates from Sherman include QUEST thresholding, per trial expectation
% reminders, 4-level confidence + 2AFC w/ single click, continuous 
% performance staircasing, and flashed-screen "deterrent" for poor 
% central task performance.
%
% Participants should sit 60cm from screen for optimal presentation and
% (ideally) in darkened room for true replication of Sherman. 

function runExp
%% FIXATION CROSS
% 0.38° x 0.38° visual angle
% Random presentation interval 500-1500ms

%% LETTER TASK: 'T' PRESENT OR ABSENT
% 4 letters presented around fixation cross (0°, 90°, 180°, 270°)
% 1.43° x 1.43° visual angle
% Random rotation angle for each letter (0° to 359°)
% Masked by letter F at same rotation angle

%% PERIPHERAL GABOR TASK: GABOR PRESENT OR ABSENT
% Spatial frequency 2c/°, Gaussian SD = 2°

% Commenting out for experiment
% dbstop if error

% Add supporting functions to path
addpath('../addons');

% Collect participant details
Subj.subjID = input('Subject''s initials:\n','s');  % Subject initials
Subj.subjNo = input('Subject number (1 to 99):\n','s');  % Enter a subject number

% Load Psychtoolbox parameters and OpenWindow
Cfg = screen_parameters;

% Raw data save location
data_save = '../../data/raw/';

if ~exist([data_save Subj.subjNo '_' Subj.subjID '/' Subj.subjNo '_' Subj.subjID '_settings.mat'],'file')
    
    % Start clock
    tic
    
    % Creates trials (if not created) and assigns appropriate order
    % Runs training and contrast thresholding using run_training
    Subj = assign_order(Subj,Cfg);
    
    % Take a break
    Screen('FillRect', Cfg.windowPtr, Cfg.black);
    
    Screen('DrawTexture',Cfg.windowPtr,Cfg.instruct.intermish);
    
    Screen('Flip', Cfg.windowPtr, [], []);
    
    WaitSecs(4);
    
    while (1)
        [~,~,buttons] = GetMouse(Cfg.windowPtr);
        if buttons(1) || KbCheck
            break;
        end
    end
    
else
    
    % Assume experiment has crashed, continue from after training if
    % thresholding has been completed successfully
    load([data_save Subj.subjNo '_' Subj.subjID '/' Subj.subjNo '_' Subj.subjID '_settings'])
    
    if isempty(Subj.training_time)
        % Training attempted but not complete
        
        % Restart clock
        tic
        
        % Creates trials (if not created) and assigns appropriate order
        % Runs training and contrast thresholding using run_training
        Subj = assign_order(Subj,Cfg);
        
        % Take a break
        Screen('FillRect', Cfg.windowPtr, Cfg.black);
        
        Screen('DrawTexture',Cfg.windowPtr,Cfg.instruct.intermish);
        
        Screen('Flip', Cfg.windowPtr, [], []);
        
        WaitSecs(4);
        
        while (1)
            [~,~,buttons] = GetMouse(Cfg.windowPtr);
            if buttons(1) || KbCheck
                break;
            end
        end
    end
end

% Load QUEST
load([data_save Subj.subjNo '_' Subj.subjID '/' Subj.subjNo '_' Subj.subjID '_settings'])

%% EXPERIMENT CONTINUES:
% 36 blocks of 12 trials, split over 6 "run" files to save memory

% Restart clock
tic

% Experiment instructions
Screen('FillRect', Cfg.windowPtr, Cfg.black);

Screen('DrawTexture',Cfg.windowPtr,Cfg.instruct.part_2);

Screen('Flip', Cfg.windowPtr, [], []);

WaitSecs(8);

while (1)
    [~,~,buttons] = GetMouse(Cfg.windowPtr);
    if buttons(1) || KbCheck
        break;
    end
end

%% LOAD TRIALS & GET STARTED
% Location of subjects trials:
trial_location = ['../../trials/' Subj.subjNo '_' Subj.subjID '/'];

% Count "run" files in location (should be 6)
trial_files = dir([trial_location Subj.subjNo '_' Subj.subjID '_r*.mat']);

for this_run = 1:length(trial_files(~([trial_files.isdir])))
    
    % Load file: i.e. 99_JM_r1.mat
    load([trial_location Subj.subjNo '_' Subj.subjID '_r' num2str(this_run) '.mat'])
    
    for block = 1:length(run_file)
        
        % Check the condition
        this_condition = run_file(block).condition;
        
        % Select the appropriate instructions
        switch this_condition
            case 'lo_full'
                the_instr = Cfg.instruct.full_25;
            case 'med_full'
                the_instr = Cfg.instruct.full_50;
            case 'hi_full'
                the_instr = Cfg.instruct.full_75;
            case 'lo_divert'
                the_instr = Cfg.instruct.divert_25;
            case 'med_divert'
                the_instr = Cfg.instruct.divert_50;
            case 'hi_divert'
                the_instr = Cfg.instruct.divert_75;
        end
        
        %% PRESENT TRIALS & SAVE
        
        % Colour screen black
        Screen('FillRect', Cfg.windowPtr, Cfg.black);
        
        % Present appropriate instruction screen
        Screen('DrawTexture',Cfg.windowPtr,the_instr);
        
        Screen('Flip', Cfg.windowPtr, [], []);
        
        WaitSecs(2);
        
        while (1)
            [~,~,buttons] = GetMouse(Cfg.windowPtr);
            if buttons(1) || KbCheck
                break;
            end
        end
        
        % Load trials for this staircase (#1 = GABOR CONTRAST)
        TR = run_file(block).TR;
        
        for tr = 1:length(TR)
            %% RETRIEVE QUEST ESTIMATES
            TR(tr).cSOA = round(QuestMean(centQUEST));
            TR(tr).true_cSOA = QuestMean(centQUEST);
            
            if ~isempty(strfind(this_condition,'divert'))
                % Present dual-task contrast
                TR(tr).pCON = QuestMean(periQUEST.divert);
            else
                % Present single-task contrast
                TR(tr).pCON = QuestMean(periQUEST.full);
            end
            
            % Should be controlled in training but accounts for negative
            % contrast
            if TR(tr).pCON < 0
                TR(tr).pCON = 0;
            end
            
            %% PRESENT LETTER STIMULI & GABOR
            
            % Letter mask is presented for 300ms as per Sherman
            mask_frames = round(300*Cfg.FrameRate/1000);
            
            % Total number of frames to be presented
            nFrames = TR(tr).cSOA + mask_frames;
            
            % Letter display starts immediately and ends after cSOA frames
            % letter_start = 1;
            letter_dur = TR(tr).cSOA + 1;
            
            % Peripheral fades in immediately and ends after pSOA frames
            peri_dur = TR(tr).pSOA + 1;
            
            %% DEFINE GABOR LOCATION & MAKE IT
            if strcmp(TR(tr).peripheral_type,'present')
                
                % Creates gabor tex using details from GaborDemo
                % gabor_height is always the same, just grabbing a default from
                % our first trial
                gabor_height = TR(tr).imageHeightPeriph;
                
                % 50% in either phase 45 or 225
                flipcoin = randi(2);
                
                if flipcoin == 1
                    GABOR_PHASE = 45;
                else
                    GABOR_PHASE = 225;
                end
                
                % Angle of gabor
                GABOR_ROTATE = 45;
                
                % Using a default contrast for the moment (0.8)
                [gabortex,gaborprops] = make_gabor(Cfg,gabor_height,0.8,GABOR_PHASE);
                
                % Gabor located here on this trial
                GABORPOSITION = CenterRectOnPoint([0,0,gabor_height,gabor_height],...
                    TR(tr).centrePeriRect(1),TR(tr).centrePeriRect(2));
                
                % Proportion of frames of peripheral presentation that fade
                % begins and ends
                fade_factor = round(peri_dur*.25);
                
                % Build vector TR(tr).pSOA long that fades in and out
                faders = ones(1,TR(tr).pSOA);
                faders(:) = TR(tr).pCON;
                
                if TR(tr).pCON > 0
                    fadein = 0:TR(tr).pCON/fade_factor:TR(tr).pCON;
                
                    faders(1:(fade_factor+1)) = fadein;
                    faders(end-fade_factor:end) = fliplr(fadein);
                end
                
            end
            
            %% BUILD 'T', 'L', 'F' TEXTURES
            % Psychtoolbox does not provide a simple text rotation function,
            % this work-around is provided in RotatingTextDemo by Peter Scarfe
            % Essentially, textures are made for each character
            
            % Draw TLF text in the middle of the screen to define bounds
            Screen('TextSize', Cfg.windowPtr, TR(tr).letterHeight);
            [~, ~, T_bounds] = DrawFormattedText(Cfg.windowPtr, 'T', 'center', 'center', Cfg.white);
            [~, ~, L_bounds] = DrawFormattedText(Cfg.windowPtr, 'L', 'center', 'center', Cfg.white);
            [~, ~, F_bounds] = DrawFormattedText(Cfg.windowPtr, 'F', 'center', 'center', Cfg.white);
            
            % Over-write the screen in WinColor so that it is back to its
            % original state
            Screen('FillRect', Cfg.windowPtr, Cfg.WinColor);
            
            % Make a rectangular texture to hold our text. This has the same
            % background color to that of the screen. Note also, that we
            % increase the size of the text bounds slightly and round upwards
            % to the nearest pixel. This is to make sure the text fits in the
            % texture and because texture dimensions can only be to interger
            % pixels.
            
            % We'll do this for each character in case the font choice has
            % kerning differences
            T_Rect = ones(ceil((T_bounds(4) - T_bounds(2)) * 1.05),...
                ceil((T_bounds(3) - T_bounds(1)) * 1.05)) .* Cfg.WinColor;
            T_tex = Screen('MakeTexture', Cfg.windowPtr, T_Rect);
            
            L_Rect = ones(ceil((L_bounds(4) - L_bounds(2)) * 1.05),...
                ceil((L_bounds(3) - L_bounds(1)) * 1.05)) .* Cfg.WinColor;
            L_tex = Screen('MakeTexture', Cfg.windowPtr, L_Rect);
            
            F_Rect = ones(ceil((F_bounds(4) - F_bounds(2)) * 1.05),...
                ceil((F_bounds(3) - F_bounds(1)) * 1.05)) .* Cfg.WinColor;
            F_tex = Screen('MakeTexture', Cfg.windowPtr, F_Rect);
            
            % Set the text size for this texture and then draw our text to the
            % texture, just as if we were drawing it to the screen
            Screen('TextSize', T_tex, TR(tr).letterHeight);
            Screen('TextSize', L_tex, TR(tr).letterHeight);
            Screen('TextSize', F_tex, TR(tr).letterHeight);
            
            % Now draw our text, but here we draw it to a texture "pretending"
            % that it is the screen
            DrawFormattedText(T_tex, 'T', 'center', 'center', Cfg.white);
            DrawFormattedText(L_tex, 'L', 'center', 'center', Cfg.white);
            DrawFormattedText(F_tex, 'F', 'center', 'center', Cfg.white);
            
            % Preallocate variables for coordinates
            L_coords = nan(length(TR(tr).letter_angles));
            F_coords = nan(length(TR(tr).letter_angles));
            
            % Define coordinates of Ls & Fs
            for letter = 1:length(TR(tr).letter_angles)
                % Define 4 coordinates of L_Rect
                L_dims = [0,0,fliplr(size(L_Rect))];
                
                % Define 4 coordinates of F_Rect
                F_dims = [0,0,fliplr(size(F_Rect))];
                
                L_coords(:,letter) = CenterRectOnPoint(L_dims,...
                    TR(tr).centreLetterRect(letter,1),...
                    TR(tr).centreLetterRect(letter,2));
                
                F_coords(:,letter) = CenterRectOnPoint(F_dims,...
                    TR(tr).centreLetterRect(letter,1),...
                    TR(tr).centreLetterRect(letter,2));
            end
            
            % Define coordinates of T (if applicable)
            if strcmp(TR(tr).letter_condition,'Ls+T')
                % Define 4 coordinates of T_Rect
                T_dims = [0,0,fliplr(size(T_Rect))];
                
                T_coords = CenterRectOnPoint(T_dims,...
                    TR(tr).centreLetterRect(TR(tr).replace_LURD,1),...
                    TR(tr).centreLetterRect(TR(tr).replace_LURD,2));
            end
            
            %% PRESENT TRIAL
            
            HideCursor;
            
            % Clear screen
            Screen('FillRect', Cfg.windowPtr, Cfg.WinColor);
            
            % Draw the lines for 'fixation_duration' many frames
            Screen('DrawLines', Cfg.windowPtr, ...
                Cfg.crossLines, Cfg.crossWidth, Cfg.crossColour, [Cfg.xCentre, Cfg.yCentre]);
            
            for m = 1 : (TR(tr).fixation_duration)
                Screen('Flip', Cfg.windowPtr, [], Cfg.aux_buffer);
            end
            
            %% PRESENT TRIAL
            for f = 1:nFrames
                
                % Clear screen
                Screen('FillRect', Cfg.windowPtr, Cfg.WinColor);
                
                if f < letter_dur
                    %% PRESENT THE LETTER STIMULI
                    % Draw 4 letter textures with rotation/etc defined from TR
                    % Left-Up-Right-Down, all Ls by default
                    Screen('DrawTextures', Cfg.windowPtr, L_tex,[],...
                        L_coords, TR(tr).letter_angles);
                    
                    if strcmp(TR(tr).letter_condition,'Ls+T')
                        %% OVERWRITE WITH 'T' TEXTURE
                        Screen('DrawTexture', Cfg.windowPtr, T_tex,[],...
                            T_coords,TR(tr).letter_angles(TR(tr).replace_LURD));
                    end
                    
                else
                    %% PRESENT THE LETTER MASKS
                    Screen('DrawTextures', Cfg.windowPtr, F_tex,[],...
                        F_coords, TR(tr).letter_angles);
                    
                end
                
                if f < peri_dur && strcmp(TR(tr).peripheral_type,'present')
                    %% PRESENT PERIPHERAL STIMULI
                    % Fade in and out
                    
                    gaborprops(4) = faders(f); % Contrast
                    
                    Screen('DrawTextures', Cfg.windowPtr, gabortex, [], ...
                        GABORPOSITION, GABOR_ROTATE, [], [], [], [],...
                        kPsychDontDoRotation, gaborprops');
                    
                end
                
                % Draw fixation cross to screen
                Screen('DrawLines', Cfg.windowPtr, ...
                    Cfg.crossLines, Cfg.crossWidth, Cfg.crossColour, [Cfg.xCentre, Cfg.yCentre]);
                
                % Present stimuli on screen for this frame
                Screen('Flip', Cfg.windowPtr);
                
            end
            
            %% RESPONSE SCREEN
            question_string = sprintf('Was the Grill (P)resent or (A)bsent?');
            
            ShowCursor;
            
            Cfg = response_screen(Cfg,question_string,'P','A');
            
            WaitSecs(.3);
            
            %% CLOSE TEXTURES FOR THIS TRIAL
            % Preserves memory for other trials
            Screen('Close',F_tex);
            Screen('Close',T_tex);
            Screen('Close',L_tex);
            
            if strcmp(TR(tr).peripheral_type,'present')
                Screen('Close',gabortex);
            end
            
            %% COLLECT A RESPONSE
            
            clicks = 0;
            
            peri_start_time = GetSecs;
            
            % Wait until subject has given a response
            while clicks == 0
                
                [x,y] = getMouseResponse();
                
                % Check whether the click went inside a box area
                for m = 1 : size(Cfg.polyL, 1)
                    idxs_left(m) = inpolygon(x,y,squeeze(Cfg.polyL(m,1,:)),squeeze(Cfg.polyL(m,2,:))); %#ok<*AGROW>
                    
                    idxs_right(m) = inpolygon(x,y,squeeze(Cfg.polyR(m,1,:)),squeeze(Cfg.polyR(m,2,:)));
                end
                
                idx_pos_left = find(idxs_left == 1);
                idx_pos_right = find(idxs_right == 1);
                
                % Left boxes click
                if length(idx_pos_left) == 1 %~isempty(idx_pos_left)
                    keyid = -1;
                    keyid2 = idx_pos_left;
                    
                    clicks = 1;
                    
                    % Paint selected box blue
                    Screen('FillPoly', Cfg.windowPtr, [0 0 255], squeeze(Cfg.polyL(idx_pos_left,:,:))',1);
                    for wait = 1:10
                        Screen('Flip', Cfg.windowPtr,  [], Cfg.aux_buffer);
                    end
                    
                end
                
                if length(idx_pos_right) == 1 %~isempty(idx_pos_right)
                    keyid = 1;
                    keyid2 = idx_pos_right;
                    
                    clicks= 1;
                    
                    % Paint selected box blue
                    Screen('FillPoly', Cfg.windowPtr, [0 0 255], squeeze(Cfg.polyR(idx_pos_right,:,:))',1);
                    for wait = 1:10
                        Screen('Flip', Cfg.windowPtr,  [], Cfg.aux_buffer);
                    end
                    
                end
            end
            
            peri_end_time = GetSecs;
            TR(tr).peripheral_RT = peri_end_time - peri_start_time; % In seconds
            
            % Define response to TR
            if keyid == -1
                % Response on left: 'present'
                TR(tr).mouseResponsesPer = 'present';
                
            elseif keyid == 1
                % Response on right: 'absent'
                TR(tr).mouseResponsesPer = 'absent';
            end
            
            %% EXAMINE RESPONSE AND UPDATE QUEST
            if strcmp(TR(tr).mouseResponsesPer,TR(tr).peripheral_type)
                TR(tr).peripheral_accuracy = 1;
            else
                TR(tr).peripheral_accuracy = 0;
            end
            
            TR(tr).peripheral_confid = keyid2;
            TR(tr).peripheral_mouse_pos = [x y];
            
            if ~isempty(strfind(this_condition,'divert'))
                %% DUAL TASK CONDITION
                % RESPONSE SCREEN #2: CENTRAL TASK
                question_string = sprintf('Did the letter (T) appear or (Ls)?');
                
                ShowCursor;
                
                Cfg = response_screen(Cfg,question_string,'T','Ls');
                
                WaitSecs(.3);
                
                %% COLLECT RESPONSE #2: 'T' PRESENCE
                
                clicks = 0;
                
                cent_start_time = GetSecs;
                
                % Wait until subject has given a response
                while clicks == 0
                    
                    [x,y] = getMouseResponse();
                    
                    % Check whether the click went inside a box area
                    for m = 1 : size(Cfg.polyL, 1)
                        idxs_left(m) = inpolygon(x,y,squeeze(Cfg.polyL(m,1,:)),squeeze(Cfg.polyL(m,2,:)));
                        
                        idxs_right(m) = inpolygon(x,y,squeeze(Cfg.polyR(m,1,:)),squeeze(Cfg.polyR(m,2,:)));
                    end
                    
                    idx_pos_left = find(idxs_left == 1);
                    idx_pos_right = find(idxs_right == 1);
                    
                    % Left boxes click
                    if length(idx_pos_left) == 1 %~isempty(idx_pos_left)
                        keyid = -1;
                        keyid2 = idx_pos_left;
                        
                        clicks = 1;
                        
                        % Paint selected box blue
                        Screen('FillPoly', Cfg.windowPtr, [0 0 255], squeeze(Cfg.polyL(idx_pos_left,:,:))',1);
                        for wait = 1:10
                            Screen('Flip', Cfg.windowPtr,  [], Cfg.aux_buffer);
                        end
                        
                    end
                    
                    if length(idx_pos_right) == 1 %~isempty(idx_pos_right)
                        keyid = 1;
                        keyid2 = idx_pos_right;
                        
                        clicks= 1;
                        
                        % Paint selected box blue
                        Screen('FillPoly', Cfg.windowPtr, [0 0 255], squeeze(Cfg.polyR(idx_pos_right,:,:))',1);
                        for wait = 1:10
                            Screen('Flip', Cfg.windowPtr,  [], Cfg.aux_buffer);
                        end
                        
                    end
                end
                
                cent_end_time = GetSecs;
                TR(tr).central_RT = cent_end_time - cent_start_time; % In seconds
                
                % Define response to TR
                if keyid == -1
                    % Response on left: 'present'
                    TR(tr).mouseResponsesCent = 'Ls+T';
                    
                elseif keyid == 1
                    % Response on right: 'absent'
                    TR(tr).mouseResponsesCent = 'Ls';
                end
                
                %% EXAMINE CENTRAL RESPONSE
                if strcmp(TR(tr).mouseResponsesCent,TR(tr).letter_condition)
                    TR(tr).central_accuracy = 1;
                else
                    TR(tr).central_accuracy = 0;
                end
                
                TR(tr).central_confid = keyid2;
                TR(tr).central_mouse_pos = [x y];
                
            end
            
            % QUEST contrasts are updated during the experiment now but exclusively
            % in the 50% expectation blocks (like in training)
            
            % This is intended to prevent contrast see-sawing between the high & low
            % conditions
            
            if ~isempty(strfind(this_condition,'med'))
                if ~isempty(strfind(this_condition,'full'))
                     periQUEST.full = QuestUpdate(periQUEST.full,...
                        TR(tr).pCON,TR(tr).peripheral_accuracy);
                else
                     periQUEST.divert = QuestUpdate(periQUEST.divert,...
                        TR(tr).pCON,TR(tr).peripheral_accuracy);
                end
            end
            
            %% CLICK TO CONTINUE...
            if tr < length(TR)
                
                % Expectation reminder
                if ~isempty(strfind(this_condition,'lo'))
                    % 25% appearance
                    expect_string = 'The Grill appears in 25% of trials this block';
                elseif ~isempty(strfind(this_condition,'med'))
                    % 50% appearance
                    expect_string = 'The Grill appears in 50% of trials this block';
                elseif ~isempty(strfind(this_condition,'hi'))
                    % 75% appearance
                    expect_string = 'The Grill appears in 75% of trials this block';
                    
                end
                
                % Clear screen
                Screen('FillRect', Cfg.windowPtr, Cfg.WinColor);
                
                Screen('TextSize',Cfg.windowPtr,24);
                
                [~,y_cent,~] = DrawFormattedText(Cfg.windowPtr, ...
                    '<< click to continue >>', 'center', 'center',...
                    Cfg.white,50,[],[],1.5);
                
                DrawFormattedText(Cfg.windowPtr, expect_string, 'center', ...
                    (y_cent+50),Cfg.white,50,[],[],1.5);
                
                Screen('Flip', Cfg.windowPtr, [], []);
                
                WaitSecs(.3);
                
                while (1)
                    [~,~,buttons] = GetMouse(Cfg.windowPtr);
                    if buttons(1) || KbCheck
                        break;
                    end
                end
                
            elseif tr == length(TR) && ~isempty(strfind(this_condition,'divert'))
                %% CENTRAL PERFORMANCE CHECK
                % Examines accuracy of block and presents on-screen reminder if
                % below 60%
                
                accuracy_this_block = mean([TR(:).central_accuracy]);
                
                if accuracy_this_block < .6
                    
                    flash_screen(Cfg.windowPtr,4,.4,'Remember to prioritise the letters in the next dual task block!')
                    
                end
            end
        end
        %% SAVE TRIAL BLOCK
    
        run_file(block).TR = TR;
        
        filename = [Subj.subjNo '_' Subj.subjID '_block' num2str(this_run)];
        save([data_save Subj.subjNo '_' Subj.subjID '/' filename '.mat'],'run_file');
        
    end
    
    if this_run < length(trial_files(~([trial_files.isdir])))
        
        % Clear screen
        Screen('FillRect', Cfg.windowPtr, Cfg.black);
        
        Screen('TextSize',Cfg.windowPtr,40);
        
        end_block = sprintf(['Run complete\n\nTake a break\n\n'...
            '<< click to continue >>']);
        
        DrawFormattedText(Cfg.windowPtr, end_block, 'center', 'center',...
            Cfg.white,50,[],[],1.5);
        
        Screen('Flip', Cfg.windowPtr, [], []);
        
        % Append contrasts and SOAs to Subj struct and save in 'data' location
        QUESTfile_update = [Subj.subjNo '_' Subj.subjID '_settings'];
        QUEST_location = [data_save Subj.subjNo '_' Subj.subjID '/'];

        save([QUEST_location QUESTfile_update '.mat'],'periQUEST','centQUEST');
        
        WaitSecs(.8);
        
        while (1)
            [~,~,buttons] = GetMouse(Cfg.windowPtr);
            if buttons(1) || KbCheck
                break;
            end
        end
        
    end
end

%% FINAL SCREEN

% Clear screen
Screen('FillRect', Cfg.windowPtr, Cfg.black);

Screen('TextSize',Cfg.windowPtr,60);

DrawFormattedText(Cfg.windowPtr, 'Thanks for participating!', 'center', 'center',...
    Cfg.white,50,[],[],1.5);

Screen('Flip', Cfg.windowPtr, [], []);

%% SAVE & CLOSE

% Append contrasts and SOAs to Subj struct and save in 'data' location
filename = [Subj.subjNo '_' Subj.subjID '_settings'];
save_location = [data_save Subj.subjNo '_' Subj.subjID '/'];

% How long did the experiment take:
Subj.experiment_time = toc;
Subj.end_time = datestr(now);
Subj.exp_duration = datestr(datenum(datevec(Subj.end_time)) - ...
    datenum(datevec(Subj.start_time)),'HH:MM:SS');
Subj.task_time = datestr((Subj.training_time + Subj.experiment_time)/(24*60*60),'HH:MM:SS');

save([save_location filename '.mat'],'Subj','Cfg','periQUEST','centQUEST');

WaitSecs(2);

ShowCursor;
sca;

clear all

end
