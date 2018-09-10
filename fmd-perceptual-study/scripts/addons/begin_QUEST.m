function [ periQUEST, centQUEST ] = begin_QUEST( Cfg, Subj )
%% BEGIN_QUEST
% Establish and assign the initial QUEST parameters for this experiment

%% DEFINE FRAMERATE CONVERSION
% Cfg.FrameRate is defined in screen_parameters.m
% Usually 60Hz which makes frame conversion (milliseconds * 0.06)
frame_conversion = Cfg.FrameRate/1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  QUEST  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% QUEST PARAMETERS
% See QuestCreate for details
% In Sherman final cSOA mean=254ms, SD=75ms
% pSOA fixed @ 388ms "with gradual onset and offset"

cSOA_time = 300; % Sherman started central SOA @ 300ms

% This is in frames, 18 for a 60Hz monitor
cSOAGuess = round(cSOA_time * frame_conversion);
cSOAGuessSd = 20; % Also in frames
pThreshold = 0.794; % Sherman et al used 79.4% performance thresholding, oddly specific

% Parameters of a Weibull psychometric function
beta = 3; % This was '2' in our original study, '3' is typical for contrast thresholding
delta = 0.05; % Fraction of trials observers press blindly, conservative
gamma = 0.5; % 0.5 for a 2AFC design

centQUEST = QuestCreate(cSOAGuess,cSOAGuessSd,pThreshold,beta,delta,gamma,1,50);
centQUEST.normalizePdf=1; % Prevents underflow errors

switch Subj.peripheral_task
    case 'face'
    % 2016-12-13: This will need work if we want to add a face task
    % SOA for peripheral face, values from Matthews et al.
    p_SOA_Guess = 10;
    p_SOA_GuessSd = 8;
    
    % Initial QUEST for single-task gabor task
    periQUEST = QuestCreate(p_SOA_Guess,p_SOA_GuessSd,pThreshold,beta,delta,gamma,.01,20);
    periQUEST.normalizePdf = 1;
    
    case 'gabor'
    
    % Gabor presence task, separate contrast QUESTs for SINGLE+DUAL
    
    % Contrast means were roughly ~0.05 for Sherman et al (SEM ~0.01)
    p_contrast_Guess = 0.2; % This is the starting contrast
    p_contrast_GuessSd = 0.075; % Conservative selection recommended
    
    % Initial QUEST for single-task gabor task
    periQUEST.full = QuestCreate(p_contrast_Guess,p_contrast_GuessSd,pThreshold,beta,delta,gamma,.01,20);
    periQUEST.full.normalizePdf = 1;
    
    % We're doing this in the training script now....
    % Separate QUEST for dual-task gabor task: 
    % periQUEST.divert = QuestCreate(p_contrast_Guess,p_contrast_GuessSd,pThreshold,beta,delta,gamma,.01,20);
    % periQUEST.divert.normalizePdf = 1;
    
end


end

