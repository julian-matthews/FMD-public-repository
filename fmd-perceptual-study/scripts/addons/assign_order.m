function [ Subj] = assign_order( Subj , Cfg )
%% FMD ASSIGN_ORDER
% Assigns block order
% Creates trials using the appropriate create_trials.m file
% Establishes QUEST using begin_QUEST.m
% Accounts for crashes by restarting from latest block

%% PERIPHERAL DISCRIMINATION
% Will likely remain unchanged but can specify here.
% Matthews et al., Face-Gender Task = 'face'
% Sherman et al., Gabor-Presence Task = 'gabor'
% As of 2016-12-14 - additional coding required to add Face-Gender task
Subj.peripheral_task = 'gabor';

%% BUILD subj STRUCT

% Builds directory for data ie. fmd-dual-task/data/raw/99_JM/
% If folder does not exist we assume this participant needs their task
% order assigned and trials created
if ~exist(['../../data/raw/' Subj.subjNo '_' Subj.subjID],'dir')
    mkdir('../../data/raw/', [Subj.subjNo '_' Subj.subjID]);
end

% Checks whether trial folder exists for this subject and creates
if ~exist(['../../trials/' Subj.subjNo '_' Subj.subjID],'dir')
    mkdir('../../trials/', [Subj.subjNo '_' Subj.subjID]);
end

trial_save = ['../../trials/' Subj.subjNo '_' Subj.subjID '/'];

% Reset random number generator by the clock time & save
t = clock;
Subj.RNG_seed = t(3)*t(4)*t(5);
rng(Subj.RNG_seed,'twister')

Subj.start_time = datestr(now);
Subj.end_time = [];
Subj.exp_duration = [];
Subj.task_time = [];
Subj.training_time = [];
Subj.experiment_time = [];

Subj.nTrials = 12;

all_conditions = {'lo_full','med_full','hi_full','lo_divert','med_divert','hi_divert'};

% Unique shuffled order for this participant
Subj.condition_order = Shuffle(all_conditions);

%% SUBJECT PARAMETERS
% Build and save trials for this participant (6 runs of 6 blocks)
create_trials(Cfg, Subj, trial_save);

% Establish QUEST parameters
[periQUEST,centQUEST] = begin_QUEST(Cfg, Subj);

%% TRAINING
% 3x staircase procedures establishing:
% 1. Gabor contrast (full attention)
% 2. Letter SOA (full attention)
% 2. Gabor contrast (diverted attention)
run_training(Cfg,Subj,periQUEST,centQUEST)

end