function run_training( Cfg,Subj,periQUEST,centQUEST )
%% RUN_TRAINING
% Conducts training as per Sherman et al (2015), three staircase procedures
% 1. GABOR CONTRAST under full attention
% 2. LETTER SOA under full attention
% 3. GABOR CONTRAST under diverted attention
%
% Two interleaved identical staircases which terminate after 8 reversals
% Performance set @ 79.4% and no confidence judgment acquired
% Gabor expectation fixed @ 50%

% If saving screencaptures for presentation/etc
screengrab = 0; % Set to 0 or 1

% Staircase cycles: 8 sets of 12 trials in Sherman et al. (2015)
% This is equivalent to Matthews et al: 2x48 trials per condition
reversals = 8;

% Same location as regular data with '_training' suffix appended
training_save = '../../data/raw/';

% 'data' location
settingsfile = [Subj.subjNo '_' Subj.subjID '_settings'];
save_location = ['../../data/raw/' Subj.subjNo '_' Subj.subjID '/'];

if ~exist([training_save Subj.subjNo '_' Subj.subjID '/' Subj.subjNo '_' Subj.subjID '_training.mat'],'file')
    % Define structure of training trials if not created already
    create_training(Cfg,Subj,reversals,training_save)
end

% Load training struct, ready for implementation
filename = [Subj.subjNo '_' Subj.subjID '_training'];
load([training_save Subj.subjNo '_' Subj.subjID '/' filename '.mat'] ,'stair')

if isempty(stair(1).block) %#ok<*NODEF>
    %% STAIRCASE #1: GABOR CONTRAST (full attention)
    % Letter task SOA fixed @ 300ms, contrast titrated from 5%
    % "Was gabor present or absent?"
    
    % Colour screen black
    Screen('FillRect', Cfg.windowPtr, Cfg.black);
    
    Screen('DrawTexture',Cfg.windowPtr,Cfg.instruct.train_peri);
    
    Screen('Flip', Cfg.windowPtr, [], []);
    
    % Screenshots for presentation
    if screengrab
        capture_count = 1; %#ok<*UNRCH>
        imageArray = Screen('GetImage', Cfg.windowPtr);
        imwrite(imageArray, ['screenshot' mat2str(capture_count) '.jpg'])
    end
    
    WaitSecs(2);
    
    while (1)
        [~,~,buttons] = GetMouse(Cfg.windowPtr);
        if buttons(1) || KbCheck
            break;
        end
    end
    
    % Start looping training set 'reversals' number of times (default=8)
    for iteration = 1:reversals
        
        % Load trials for this staircase (#1 = GABOR CONTRAST)
        TR = stair(1).trial_data(iteration).TR;
        
        for tr = 1:length(TR)
            %% RETRIEVE QUEST ESTIMATES
            TR(tr).cSOA = round(QuestMean(centQUEST));
            TR(tr).true_cSOA = QuestMean(centQUEST);
            TR(tr).pCON_full = QuestMean(periQUEST.full);
            
            if TR(tr).pCON_full < 0
                TR(tr).pCON_full = 0;
            end
            
            %% PRESENT LETTER STIMULI & GABOR
            % Presentation parameters
            
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
                faders(:) = TR(tr).pCON_full;
                
                if TR(tr).pCON_full > 0
                    fadein = 0:TR(tr).pCON_full/fade_factor:TR(tr).pCON_full;
                
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
            
            %% PRESENT FIXATION CROSS
            HideCursor;
            
            % Clear screen
            Screen('FillRect', Cfg.windowPtr, Cfg.WinColor);
            
            % Draw the lines for 'fixation_duration' many frames
            Screen('DrawLines', Cfg.windowPtr, ...
                Cfg.crossLines, Cfg.crossWidth, Cfg.crossColour, [Cfg.xCentre, Cfg.yCentre]);
            
            for m = 1 : (TR(tr).fixation_duration)
                Screen('Flip', Cfg.windowPtr, [], Cfg.aux_buffer);
            end
            
            % Screenshots for presentation
            if screengrab
                capture_count = capture_count + 1;
                imageArray = Screen('GetImage', Cfg.windowPtr);
                imwrite(imageArray, ['screenshot' mat2str(capture_count) '.jpg'])
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
                
                if f == round(TR(tr).pSOA/2) || f == nFrames-1
                    % Screenshot with gabor (if present) and with letter mask
                    
                    % Screenshots for presentation
                    if screengrab
                        capture_count = capture_count + 1;
                        imageArray = Screen('GetImage', Cfg.windowPtr);
                        imwrite(imageArray, ['screenshot' mat2str(capture_count) '.jpg'])
                    end
                    
                end
                
            end
            
            %% RESPONSE SCREEN
            question_string = sprintf('Was the Grill (P)resent or (A)bsent?');
            
            ShowCursor;
            
            Cfg = response_screen(Cfg,question_string,'P','A');
            
            % Screenshots for presentation
            if screengrab
                capture_count = capture_count + 1;
                imageArray = Screen('GetImage', Cfg.windowPtr);
                imwrite(imageArray, ['screenshot' mat2str(capture_count) '.jpg'])
            end
            
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
            
            peri_end_time = GetSecs;
            TR(tr).peripheral_RT = peri_end_time - peri_start_time; % In seconds
            
            % Screenshots for presentation
            if screengrab
                capture_count = capture_count + 1;
                imageArray = Screen('GetImage', Cfg.windowPtr);
                imwrite(imageArray, ['screenshot' mat2str(capture_count) '.jpg'])
            end
            
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
            
            periQUEST.full = QuestUpdate(periQUEST.full,TR(tr).pCON_full,TR(tr).peripheral_accuracy);
            
            %% CLICK TO CONTINUE...
            
            if tr < length(TR)
                % Clear screen
                Screen('FillRect', Cfg.windowPtr, Cfg.WinColor);
                
                DrawFormattedText(Cfg.windowPtr, '<< click to continue >>', 'center', 'center',...
                    Cfg.white,50,[],[],1.5);
                
                Screen('Flip', Cfg.windowPtr, [], []);
                
                % Screenshots for presentation
                if screengrab
                    capture_count = capture_count + 1;
                    imageArray = Screen('GetImage', Cfg.windowPtr);
                    imwrite(imageArray, ['screenshot' mat2str(capture_count) '.jpg'])
                end
                
                WaitSecs(.3);
                
                while (1)
                    [~,~,buttons] = GetMouse(Cfg.windowPtr);
                    if buttons(1) || KbCheck
                        break;
                    end
                end
                
            end
            
        end
        
        if iteration < reversals
            
            % Clear screen
            Screen('FillRect', Cfg.windowPtr, Cfg.black);
            
            Screen('TextSize',Cfg.windowPtr,30);
            
            end_block = sprintf(['Time for a break\n\n\n'...
                '<< click to continue >>'],iteration,reversals);
            
            DrawFormattedText(Cfg.windowPtr, end_block, 'center', 'center',...
                Cfg.white,50,[],[],1.5);
            
            Screen('Flip', Cfg.windowPtr, [], []);
            WaitSecs(.8);
            
            while (1)
                [~,~,buttons] = GetMouse(Cfg.windowPtr);
                if buttons(1) || KbCheck
                    break;
                end
            end
            
        end
        
        %% SAVE TRIAL BLOCK TO 'stair'
        
        stair(1).trial_data(iteration).TR = TR;
        
        filename = [Subj.subjNo '_' Subj.subjID '_training'];
        save([training_save Subj.subjNo '_' Subj.subjID '/' filename '.mat'],'stair');
        
    end
    
    %% STAIR(1) COMPLETE + SAVE
    stair(1).block = 1;
    
    save([save_location settingsfile '.mat'],'Subj','Cfg','periQUEST','centQUEST');
    
    save([training_save Subj.subjNo '_' Subj.subjID '/' filename '.mat'],'stair');
    
else
    % Load the appropriate QUEST
    load([save_location settingsfile '.mat'],'periQUEST','centQUEST');
    
end

if isempty(stair(2).block)
    
    %% STAIRCASE #2: LETTER SOA
    % Gabor contrast fixed from Staircase #1, letter SOA titrated from 300ms
    % "Was letter target 'T' present or absent?"
    
    % Colour screen black
    Screen('FillRect', Cfg.windowPtr, Cfg.black);
    
    % Present instruction screen and click to start
    
    Screen('DrawTexture',Cfg.windowPtr,Cfg.instruct.train_cent);
    
    Screen('Flip', Cfg.windowPtr, [], []);
    
    WaitSecs(2);
    
    while (1)
        [~,~,buttons] = GetMouse(Cfg.windowPtr);
        if buttons(1) || KbCheck
            break;
        end
    end
    
    % Start looping training set 'reversals' number of times (default=8)
    for iteration = 1:reversals
        
        % Load trials for this staircase (#2 = LETTER SOA)
        TR = stair(2).trial_data(iteration).TR;
        
        for tr = 1:length(TR)
            %% RETRIEVE QUEST ESTIMATES
            TR(tr).cSOA = round(QuestMean(centQUEST));
            TR(tr).true_cSOA = QuestMean(centQUEST);
            TR(tr).pCON_full = QuestMean(periQUEST.full);
            
            if TR(tr).pCON_full < 0
                TR(tr).pCON_full = 0;
            end
            
            %% PRESENT LETTER STIMULI & GABOR
            % Presentation parameters
            
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
                faders(:) = TR(tr).pCON_full;
                
                if TR(tr).pCON_full > 0
                    fadein = 0:TR(tr).pCON_full/fade_factor:TR(tr).pCON_full;
                
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
            
            %% PRESENT FIXATION CROSS
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
            question_string = sprintf('Did the letter (T) appear or (Ls)?');
            
            ShowCursor;
            
            Cfg = response_screen(Cfg,question_string,'T','Ls');
            
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
            
            %% EXAMINE RESPONSE AND UPDATE QUEST
            if strcmp(TR(tr).mouseResponsesCent,TR(tr).letter_condition)
                TR(tr).central_accuracy = 1;
            else
                TR(tr).central_accuracy = 0;
            end
            
            TR(tr).central_confid = keyid2;
            TR(tr).central_mouse_pos = [x y];
            
            centQUEST = QuestUpdate(centQUEST,TR(tr).cSOA,TR(tr).central_accuracy);
            
            %% CLICK TO CONTINUE...
            
            if tr < length(TR)
                % Clear screen
                Screen('FillRect', Cfg.windowPtr, Cfg.WinColor);
                
                DrawFormattedText(Cfg.windowPtr, '<< click to continue >>', 'center', 'center',...
                    Cfg.white,50,[],[],1.5);
                
                Screen('Flip', Cfg.windowPtr, [], []);
                
                WaitSecs(.3);
                
                while (1)
                    [~,~,buttons] = GetMouse(Cfg.windowPtr);
                    if buttons(1) || KbCheck
                        break;
                    end
                end
                
            end
            
        end
        
        if iteration < reversals
            
            % Clear screen
            Screen('FillRect', Cfg.windowPtr, Cfg.black);
            
            Screen('TextSize',Cfg.windowPtr,30);
            
            end_block = sprintf(['Time for a break\n\n\n'...
                '<< click to continue >>'],iteration,reversals);
            
            DrawFormattedText(Cfg.windowPtr, end_block, 'center', 'center',...
                Cfg.white,50,[],[],1.5);
            
            Screen('Flip', Cfg.windowPtr, [], []);
            WaitSecs(.8);
            
            while (1)
                [~,~,buttons] = GetMouse(Cfg.windowPtr);
                if buttons(1) || KbCheck
                    break;
                end
            end
            
        end
        
        %% SAVE TRIAL BLOCK TO 'stair'
        
        stair(2).trial_data(iteration).TR = TR;
        
        filename = [Subj.subjNo '_' Subj.subjID '_training'];
        save([training_save Subj.subjNo '_' Subj.subjID '/' filename '.mat'],'stair');
        
    end
    
    %% STAIR(2) COMPLETE + SAVE
    stair(2).block = 1;
    
    save([save_location settingsfile '.mat'],'Subj','Cfg','periQUEST','centQUEST');
    
    save([training_save Subj.subjNo '_' Subj.subjID '/' filename '.mat'],'stair');
    
else
    % Load the appropriate QUESTs
    load([save_location settingsfile '.mat'],'periQUEST','centQUEST')
    
end

if isempty(stair(3).block)
    %% STAIRCASE #3: GABOR CONTRAST (diverted attention)
    % Letter SOA fixed from Staircase #2, gabor contrast titrated from 1.05 Staircase #1
    % "Was gabor present or absent?" "Was letter target 'T' present or absent?"
    %
    % If letter task below 60% for block:
    % "Maintain concentration on the visual search task"
    
    % Colour screen black
    Screen('FillRect', Cfg.windowPtr, Cfg.black);
    
    % Present instruction screen
    Screen('DrawTexture',Cfg.windowPtr,Cfg.instruct.train_dual);
    
    Screen('Flip', Cfg.windowPtr, [], []);
    
    %% CREATE QUEST FOR DIVERTED ATTENTION GABOR
    
    % Load previous gabor contrast
    pCON_full = QuestMean(periQUEST.full);
    
    % Initial contrast starts from contrast in full attention condition
    p_contrast_Guess = 1.05 * pCON_full; % From Sherman et al.
    
    % Uses defaults from begin_QUEST.m
    periQUEST.divert = QuestCreate(p_contrast_Guess,.05,.794,3,.05,.5,.01,20);
    periQUEST.divert.normalizePdf = 1;
    
    %% CLICK AND CONTINUE
    
    WaitSecs(2);
    
    while (1)
        [~,~,buttons] = GetMouse(Cfg.windowPtr);
        if buttons(1) || KbCheck
            break;
        end
    end
    
    % Start looping training set 'reversals' number of times (default=8)
    for iteration = 1:reversals
        
        % Load trials for this staircase (#2 = LETTER SOA)
        TR = stair(3).trial_data(iteration).TR;
        
        for tr = 1:length(TR)
            %% RETRIEVE QUEST ESTIMATES
            TR(tr).cSOA = round(QuestMean(centQUEST));
            TR(tr).true_cSOA = QuestMean(centQUEST);
            TR(tr).pCON_divert = QuestMean(periQUEST.divert);
            
            if TR(tr).pCON_divert < 0
                TR(tr).pCON_divert = 0;
            end
            
            %% PRESENT LETTER STIMULI & GABOR
            % Presentation parameters
            
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
                faders(:) = TR(tr).pCON_divert;
                
                if TR(tr).pCON_divert > 0
                    fadein = 0:TR(tr).pCON_divert/fade_factor:TR(tr).pCON_divert;
                
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
            
            %% PRESENT FIXATION CROSS
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
            
            %% RESPONSE SCREEN #1: PERIPHERAL TASK
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
            
            %% COLLECT RESPONSE #1: GABOR PRESENCE
            
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
            
            %% EXAMINE PERIPHERAL RESPONSE
            if strcmp(TR(tr).mouseResponsesPer,TR(tr).peripheral_type)
                TR(tr).peripheral_accuracy = 1;
            else
                TR(tr).peripheral_accuracy = 0;
            end
            
            TR(tr).peripheral_confid = keyid2;
            TR(tr).peripheral_mouse_pos = [x y];
            
            %% RESPONSE SCREEN #2: CENTRAL TASK
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
            
            %% UPDATE QUEST (only peripheral contrast)
            
            periQUEST.divert = QuestUpdate(periQUEST.divert,TR(tr).pCON_divert,TR(tr).peripheral_accuracy);
            
            %% CLICK TO CONTINUE...
            
            if tr < length(TR)
                % Clear screen
                Screen('FillRect', Cfg.windowPtr, Cfg.WinColor);
                
                DrawFormattedText(Cfg.windowPtr, '<< click to continue >>', 'center', 'center',...
                    Cfg.white,50,[],[],1.5);
                
                Screen('Flip', Cfg.windowPtr, [], []);
                
                WaitSecs(.3);
                
                while (1)
                    [~,~,buttons] = GetMouse(Cfg.windowPtr);
                    if buttons(1) || KbCheck
                        break;
                    end
                end
                
            else
                %% CENTRAL PERFORMANCE CHECK
                % Examines accuracy of block and presents on-screen reminder if
                % below 60%
                
                accuracy_this_block = mean([TR(:).central_accuracy]);
                
                if accuracy_this_block < .6
                    
                    flash_screen(Cfg.windowPtr,4,.4,'Remember to prioritise the letter task!')
                    
                end
            end
            
        end
        
        if iteration < reversals
            
            % Clear screen
            Screen('FillRect', Cfg.windowPtr, Cfg.black);
            
            Screen('TextSize',Cfg.windowPtr,30);
            
            end_block = sprintf(['Time for a break\n\n\n'...
                '<< click to continue >>'],iteration,reversals);
            
            DrawFormattedText(Cfg.windowPtr, end_block, 'center', 'center',...
                Cfg.white,50,[],[],1.5);
            
            Screen('Flip', Cfg.windowPtr, [], []);
            WaitSecs(.8);
            
            while (1)
                [~,~,buttons] = GetMouse(Cfg.windowPtr);
                if buttons(1) || KbCheck
                    break;
                end
            end
            
        end
        
        %% SAVE TRIAL BLOCK TO 'stair'
        
        stair(3).trial_data(iteration).TR = TR;
        
        filename = [Subj.subjNo '_' Subj.subjID '_training'];
        save([training_save Subj.subjNo '_' Subj.subjID '/' filename '.mat'],'stair');
        
    end
    
    %% SAVE
    % How long did training take:
    Subj.training_time = toc;
    
    save([save_location settingsfile '.mat'],'Subj','Cfg','periQUEST','centQUEST');
    
    %% STAIR(3) COMPLETE + SAVE
    stair(3).block = 1; %#ok<NASGU>
    
    save([training_save Subj.subjNo '_' Subj.subjID '/' filename '.mat'],'stair');
    
else
    % Load the appropriate QUEST
    load([save_location settingsfile '.mat'],'Subj','Cfg','periQUEST','centQUEST');
    
end

end

function create_training(Cfg,Subj,reversals,training_save)
%% FMD CREATE TRAINING
% Builds TR structure for 3x training/staircase procedures

% Three staircases of 8 blocks, 288 trials total
total_staircases = 3;
total_blocks = reversals;

% Create staircase structure for storing trial information
task_order = {'ST:gabor contrast' 'letter SOA' 'DT:gabor contrast'};
stair = struct('task',task_order,'block',[],'trial_data',[]);

for staircase = 1:total_staircases
    
    for block = 1:total_blocks
        
        % Create a block of trials using fixed .5 gabor probability
        all_trials = block_build( Cfg, Subj, 'med');
        
        % Save trials to staircase struct
        stair(staircase).trial_data(block).TR = all_trials;
        
    end
    
end

%% SAVE TO data/training/
filename = [Subj.subjNo '_' Subj.subjID '_training'];
save([training_save Subj.subjNo '_' Subj.subjID '/' filename '.mat'],'stair');

end
