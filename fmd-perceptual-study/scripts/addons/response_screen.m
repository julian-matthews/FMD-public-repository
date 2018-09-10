function [Cfg] = response_screen(Cfg, Question_string,present_char,absent_char)
%% FMD RESPONSE SCREEN
% Draws/flips 2AFC + 4-lvl confidence square
%
% Requires 'Cfg' data for configuration, 'Question_string' for presentation
% at the top of the screen, and characters for present and absent judgements
% arranged on the left and right respectively.
% 
% Custom version of confidence wheel employed by Matthews et al for RSVP
% and dual-task studies

% Centre mouse on screen
SetMouse(Cfg.xCentre,Cfg.yCentre);

% If parameters for drawing screen do not exist, create now
if ~isfield(Cfg,'rectFrame')
    Cfg = configure_response_screen(Cfg);
end

% Setting (P)resent and (A)bsent as default choices
if nargin == 4
    left_text = present_char;
    right_text = absent_char;
else
    left_text = 'P';
    right_text = 'A';
end

%Draw the response boxes
Screen('FillRect',  Cfg.windowPtr, Cfg.WinColor);

Screen('FillRect', Cfg.windowPtr, Cfg.white, Cfg.bigrect_{1} );
Screen('FillRect', Cfg.windowPtr, Cfg.gray, Cfg.smallrect_{1} );
Screen('FillRect', Cfg.windowPtr, Cfg.gray, Cfg.cleavage_{1});

Screen('DrawLine', Cfg.windowPtr,Cfg.black,Cfg.smallrect_{1}(1),Cfg.smallrect_{1}(4),Cfg.bigrect_{1}(1),Cfg.bigrect_{1}(4),2);
Screen('DrawLine', Cfg.windowPtr,Cfg.black,Cfg.smallrect_{1}(1),Cfg.smallrect_{1}(2),Cfg.bigrect_{1}(1),Cfg.bigrect_{1}(2),2);
Screen('DrawLine', Cfg.windowPtr,Cfg.black,Cfg.smallrect_{1}(1),Cfg.smallrect_{1}(2),Cfg.bigrect_{1}(1),Cfg.bigrect_{1}(2),2);
Screen('DrawLine', Cfg.windowPtr,Cfg.black,Cfg.smallrect_{1}(3),Cfg.smallrect_{1}(2),Cfg.bigrect_{1}(3),Cfg.bigrect_{1}(2),2);
Screen('DrawLine', Cfg.windowPtr,Cfg.black,Cfg.smallrect_{1}(3),Cfg.smallrect_{1}(4),Cfg.bigrect_{1}(3),Cfg.bigrect_{1}(4),2);

Screen('DrawLine', Cfg.windowPtr,Cfg.black,Cfg.xCentre,Cfg.smallrect_{1}(4),Cfg.xCentre,Cfg.bigrect_{1}(4),2);
Screen('DrawLine', Cfg.windowPtr,Cfg.black,Cfg.xCentre,Cfg.smallrect_{1}(2),Cfg.xCentre,Cfg.bigrect_{1}(2),2);

Screen('DrawLine', Cfg.windowPtr,Cfg.black,Cfg.smallrect_{1}(1),Cfg.yCentre,Cfg.bigrect_{1}(1),Cfg.yCentre,1.3);
Screen('DrawLine', Cfg.windowPtr,Cfg.black,Cfg.smallrect_{1}(3),Cfg.yCentre,Cfg.bigrect_{1}(3),Cfg.yCentre,1.3);

Cfg.edge_offset = 1; % Polygon edge offset (prevents overlap)
polyL(1,:,:)=[Cfg.smallrect_{1}(1)+Cfg.edge_offset Cfg.cleavage_{1}(1) Cfg.cleavage_{1}(1) Cfg.bigrect_{1}(1)+Cfg.edge_offset;Cfg.smallrect_{1}(4) Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4) Cfg.bigrect_{1}(4)]; %Left conf 1
polyL(2,:,:)=[Cfg.smallrect_{1}(1) Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1) Cfg.bigrect_{1}(1);Cfg.smallrect_{1}(4)-Cfg.edge_offset Cfg.yCentre+Cfg.edge_offset Cfg.yCentre+Cfg.edge_offset Cfg.bigrect_{1}(4)-Cfg.edge_offset]; %2
polyL(3,:,:)=[Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1) Cfg.bigrect_{1}(1) Cfg.smallrect_{1}(1);Cfg.yCentre-Cfg.edge_offset Cfg.yCentre-Cfg.edge_offset Cfg.bigrect_{1}(2)+Cfg.edge_offset Cfg.smallrect_{1}(2)+Cfg.edge_offset]; %3
polyL(4,:,:)=[Cfg.smallrect_{1}(1)+Cfg.edge_offset Cfg.bigrect_{1}(1)+Cfg.edge_offset Cfg.cleavage_{1}(1) Cfg.cleavage_{1}(1);Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2) Cfg.bigrect_{1}(2) Cfg.smallrect_{1}(2)]; %4

polyR(1,:,:)=[Cfg.cleavage_{1}(3) Cfg.smallrect_{1}(3)-Cfg.edge_offset Cfg.bigrect_{1}(3)-Cfg.edge_offset Cfg.cleavage_{1}(3);Cfg.smallrect_{1}(4) Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4) Cfg.bigrect_{1}(4)];
polyR(2,:,:)=[Cfg.smallrect_{1}(3) Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3) Cfg.bigrect_{1}(3);Cfg.smallrect_{1}(4)-Cfg.edge_offset Cfg.yCentre+Cfg.edge_offset Cfg.yCentre+Cfg.edge_offset Cfg.bigrect_{1}(4)-Cfg.edge_offset];
polyR(3,:,:)=[Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3) Cfg.bigrect_{1}(3) Cfg.smallrect_{1}(3);Cfg.yCentre-Cfg.edge_offset Cfg.yCentre-Cfg.edge_offset Cfg.bigrect_{1}(2)+Cfg.edge_offset Cfg.smallrect_{1}(2)+Cfg.edge_offset];
polyR(4,:,:)=[Cfg.smallrect_{1}(3)-Cfg.edge_offset Cfg.cleavage_{1}(3) Cfg.cleavage_{1}(3) Cfg.bigrect_{1}(3)-Cfg.edge_offset;Cfg.smallrect_{1}(2) Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2) Cfg.bigrect_{1}(2)];

% Save the poly matrices for registering responses
Cfg.polyL = polyL;
Cfg.polyR = polyR;

% Draw 'sureness'
Screen('TextSize',Cfg.windowPtr,28);
DrawFormattedText(Cfg.windowPtr, 'absolutely certain', 'center', Cfg.bigrect_{1}(2) - 50 , [255 255 255]);
DrawFormattedText(Cfg.windowPtr, 'complete guess', 'center',  Cfg.bigrect_{1}(4) + 10 , [255 255 255]);

% Draw 'Question'
Screen('TextSize',Cfg.windowPtr,37);
DrawFormattedText(Cfg.windowPtr, Question_string, 'center', Cfg.bigrect_{1}(2) - 100 , [255 255 255]);

% Boxes for text
text_box=[0 0 1 1];
offset=(Cfg.xCentre-mean([Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1)]))/2.;

global ptb_drawformattedtext_disableClipping; % needed for DrawFormattedText with winRect ... otherwise text is not drawn
ptb_drawformattedtext_disableClipping=1;

rect_1L=text_box+[mean([Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1)])+offset mean([Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4)]) mean([Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1)])+offset mean([Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4)])];
rect_2L=text_box+[mean([Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1)]) mean([Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4)])-offset mean([Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1)]) mean([Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4)])-offset];
rect_3L=text_box+[mean([Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1)]) mean([Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2)])+offset mean([Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1)]) mean([Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2)])+offset];
rect_4L=text_box+[mean([Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1)])+offset mean([Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2)]) mean([Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1)])+offset mean([Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2)])];

rect_1R=text_box+[mean([Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3)])-offset mean([Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4)]) mean([Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3)])-offset mean([Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4)])];
rect_2R=text_box+[mean([Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3)]) mean([Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4)])-offset mean([Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3)]) mean([Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4)])-offset];
rect_3R=text_box+[mean([Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3)]) mean([Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2)])+offset mean([Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3)]) mean([Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2)])+offset];
rect_4R=text_box+[mean([Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3)])-offset mean([Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2)]) mean([Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3)])-offset mean([Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2)])];

Screen('TextSize',Cfg.windowPtr,30);

DrawFormattedText(Cfg.windowPtr, '1', 'center', 'center',  [0 0 0],[],[],[],[],[],rect_1L); 
DrawFormattedText(Cfg.windowPtr, '2', 'center', 'center',  [0 0 0],[],[],[],[],[],rect_2L);
DrawFormattedText(Cfg.windowPtr, '3', 'center', 'center',  [0 0 0],[],[],[],[],[],rect_3L);
DrawFormattedText(Cfg.windowPtr, '4', 'center', 'center',  [0 0 0],[],[],[],[],[],rect_4L);

DrawFormattedText(Cfg.windowPtr, '1', 'center', 'center',  [0 0 0],[],[],[],[],[],rect_1R);
DrawFormattedText(Cfg.windowPtr, '2', 'center', 'center',  [0 0 0],[],[],[],[],[],rect_2R);
DrawFormattedText(Cfg.windowPtr, '3', 'center', 'center',  [0 0 0],[],[],[],[],[],rect_3R);
DrawFormattedText(Cfg.windowPtr, '4', 'center', 'center',  [0 0 0],[],[],[],[],[],rect_4R);

Screen('TextSize',Cfg.windowPtr,20);

%% DISPLAY 2AFC OPTIONS

% Draw text for the left boxes
move_x = 10;
move_y= 24;

Screen('TextSize',Cfg.windowPtr,40);

Screen('DrawText', Cfg.windowPtr, left_text, Cfg.xCentre-6.8*move_x, Cfg.yCentre-move_y, Cfg.white);
Screen('DrawText', Cfg.windowPtr, right_text,  Cfg.xCentre+4*move_x, Cfg.yCentre-move_y, Cfg.white);

% Present everything
Screen('Flip',Cfg.windowPtr, [], Cfg.aux_buffer);

end

function [Cfg] = configure_response_screen(Cfg)
%% 2AFC CONFIDENCE SQUARE SETUP
% Parameters for defining the look of the confidence wheel, ensures
% response function is mobile

% Define rectangles for response collection
Cfg.frameThickness = 20;

Cfg.rs = 400; % Formerly 250

Cfg.rect = [0 0 Cfg.rs Cfg.rs];
Cfg.smallrect = [0 0 Cfg.rs/1.5 Cfg.rs/4];
Cfg.bigrect = [0 0 2.25*Cfg.rs 2*Cfg.rs];
Cfg.cleavage = [0 0 Cfg.rs/4 2*Cfg.rs];

Cfg.yoff = 0;

Cfg.screensize_r = Cfg.windowRect(4);
Cfg.screensize_c = Cfg.windowRect(3);

Cfg.rectFrame = [0 0 Cfg.width Cfg.height]+[Cfg.screensize_c/2-Cfg.width/2 Cfg.screensize_r/2-Cfg.height/2 Cfg.screensize_c/2-Cfg.width/2 Cfg.screensize_r/2-Cfg.height/2];
Cfg.rectFrame = Cfg.rectFrame + [0 Cfg.yoff 0 Cfg.yoff];
Cfg.rect_{1} = Cfg.rect + [Cfg.screensize_c/2-Cfg.rs/2 Cfg.screensize_r/2-Cfg.height/2+Cfg.frameThickness Cfg.screensize_c/2-Cfg.rs/2 Cfg.screensize_r/2-Cfg.height/2+Cfg.frameThickness]; % top quadrant
Cfg.bigrect_{1}=Cfg.bigrect + [Cfg.screensize_c/2-1.125*Cfg.rs Cfg.screensize_r/2-Cfg.rs Cfg.screensize_c/2-1.125*Cfg.rs Cfg.screensize_r/2-Cfg.rs]; % top quadrant
Cfg.smallrect_{1}=Cfg.smallrect + [Cfg.screensize_c/2-Cfg.rs/3 Cfg.screensize_r/2-Cfg.rs/8 Cfg.screensize_c/2-Cfg.rs/3 Cfg.screensize_r/2-Cfg.rs/8]; % top quadrant

Cfg.cleavage_{1}=Cfg.cleavage + [Cfg.screensize_c/2-Cfg.rs/8 Cfg.screensize_r/2-Cfg.rs Cfg.screensize_c/2-Cfg.rs/8 Cfg.screensize_r/2-Cfg.rs];

end
