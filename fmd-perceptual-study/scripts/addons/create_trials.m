function create_trials( Cfg, Subj, save_location )
%% FMD CREATE_TRIALS
% Builds TR structure for the FMD-dual-task
%
% Two letter conditions (all Ls or Ls+T) + two peripheral conditions

% Counterbalance conditions using latin square (each row is a 6th of the exp)
% This counterbalances order such that every condition has immediately
% followed each alternative, checks for carryover effects
order = [1,2,6,3,5,4;2,3,1,4,6,5;3,4,2,5,1,6;4,5,3,6,2,1;5,6,4,1,3,2;6,1,5,2,4,3];

% Create 6 runs of 6 blocks, each with 12 trials: 432 total
total_runs = 6;
total_blocks = 6;

for this_run = 1:total_runs
    
    % Define order of block conditions from latin square
    this_block_order = order(this_run,:);
    this_block_string = Subj.condition_order(this_block_order);
    
    % Create Struct for this run, 6 blocks in the above order
    % Includes a few sanity checks to confirm order & condition
    run_file = struct('run',this_run,'condition',this_block_string,'block',[],'TR',[]);
    
    for this_block = 1:total_blocks
        
        run_file(this_block).block = this_block;
        
        % Create trials using the specified expectation condition
        all_trials = block_build(Cfg, Subj, run_file(this_block).condition);
        
        % Input into this block
        run_file(this_block).TR = all_trials;
        
    end
    
    %% SAVE THIS RUN
    % Creates total of 6 runs, each containing 6 blocks of 12 trials: 432 total
    
    filename = [Subj.subjNo '_' Subj.subjID '_r' num2str(this_run)];
    save([save_location filename '.mat'],'run_file')
    
end

end