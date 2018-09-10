function [ data, demographics ] = fmd_preprocessing (raw_location,preprocessing_location)
%FMD-PREPROCESSING Runs preprocessing on raw data
%   Examines raw_data_location for FMD folders and builds 'data' struct
%   Outputs the data struct and saves to preprocessing folder
%
%   Typical locations:
%   raw_location= '/Users/julian/Desktop/fmd-dual-task/data/raw'
%   preprocessing_location= '/Users/julian/Desktop/fmd-dual-task/data/preprocessed'

demographics = readtable([preprocessing_location '/demographics.csv']);

% Go into folder for controls, fmds, and organics
groupids = {'control' 'fmd' 'organic'};

% Open each block of trials and build data structure with following:
data = struct('subjID',[],'group',[],'subjnum',[],'block',[],'task',[],'condition',[],...
    'central_accuracy',[],'central_signal',[],'central_decision',[],'central_confid',[],...
    'peri_accuracy',[],'peri_signal',[],'peri_decision',[],'peri_confid',[],...
    'cSOA',[],'pSOA',[],'pCON',[]);

count = 0;

for group = 1:length(groupids)
    groupid = groupids{group};
    
    % Reads all containing files from specified location
    raw_data = [raw_location,'/',groupid,'/'];   
    
    % Determine subject names
    temp = dir([raw_data '*FMD*']);
    files = {temp.name};
    
    for subj = 1:length(files)
        
        % Input training if it exists
        if exist([raw_data files{subj} '/' files{subj} '_training.mat'],'file')
            load([raw_data files{subj} '/' files{subj} '_training.mat'],'stair')
            
            for training = 1:length(stair)
                for block = 1:length(stair(training).trial_data)
                    
                    TR = stair(training).trial_data(block).TR;
                    
                    for trial = 1:length(TR)
                        
                        count = count + 1;
                        
                        data(count).subjID = files{subj};
                        data(count).group = groupid;
                        data(count).subjnum = subj;
                        data(count).block = (training-1)*length(stair(training).trial_data)+block-25;
                        data(count).task = stair(training).task;
                        
                        % 'lo_full' 'med_full' 'hi_full' 'lo_divert' 'med_divert' 'hi_divert'
                        switch stair(training).task
                            case 'ST:gabor contrast'
                                % ST:gabor training is 'med_full'
                                data(count).condition = 2;
                            case 'DT:gabor contrast'
                                % DT:gabor training is 'med_divert'
                                data(count).condition = 5;
                            case 'letter SOA'
                                % letter SOA is condition #0
                                data(count).condition = 0;
                        end
                        
                        if data(count).condition == 0
                            % Central task only
                            data(count).pSOA = TR(trial).pSOA; % Sanity checks
                            data(count).pCON = TR(trial).pCON_full;
                            
                            data(count).central_accuracy = TR(trial).central_accuracy;
                            
                            switch TR(trial).letter_condition
                                case 'Ls+T'
                                    data(count).central_signal = -1;
                                case 'Ls'
                                    data(count).central_signal = 1;
                            end
                            
                            switch TR(trial).mouseResponsesCent
                                case 'Ls+T'
                                    data(count).central_decision = -1;
                                case 'Ls'
                                    data(count).central_decision = 1;
                            end
                            
                            data(count).central_confid = ...
                                TR(trial).central_confid * data(count).central_decision;
                            data(count).cSOA = TR(trial).cSOA;
                            
                        elseif data(count).condition <= 3
                            % Full attention, no central task
                            data(count).cSOA = TR(trial).cSOA; % Sanity check
                            
                            data(count).peri_accuracy = TR(trial).peripheral_accuracy;
                            
                            switch TR(trial).peripheral_type
                                case 'present'
                                    data(count).peri_signal = -1;
                                case 'absent'
                                    data(count).peri_signal = 1;
                            end
                            
                            switch TR(trial).mouseResponsesPer
                                case 'present'
                                    data(count).peri_decision = -1;
                                case 'absent'
                                    data(count).peri_decision = 1;
                            end
                            
                            data(count).peri_confid = ...
                                TR(trial).peripheral_confid * data(count).peri_decision;
                            
                            data(count).pSOA = TR(trial).pSOA;
                            data(count).pCON = TR(trial).pCON_full;
                            data(count).pRT = TR(trial).peripheral_RT;
                            
                        elseif data(count).condition >= 4
                            % Diverted attention, both tasks
                            % Central
                            data(count).central_accuracy = TR(trial).central_accuracy;
                            
                            switch TR(trial).letter_condition
                                case 'Ls+T'
                                    data(count).central_signal = -1;
                                case 'Ls'
                                    data(count).central_signal = 1;
                            end
                            
                            switch TR(trial).mouseResponsesCent
                                case 'Ls+T'
                                    data(count).central_decision = -1;
                                case 'Ls'
                                    data(count).central_decision = 1;
                            end
                            
                            data(count).central_confid = ...
                                TR(trial).central_confid * data(count).central_decision;
                            data(count).cSOA = TR(trial).cSOA;
                            
                            % Peripheral
                            data(count).peri_accuracy = TR(trial).peripheral_accuracy;
                            
                            switch TR(trial).peripheral_type
                                case 'present'
                                    data(count).peri_signal = -1;
                                case 'absent'
                                    data(count).peri_signal = 1;
                            end
                            
                            switch TR(trial).mouseResponsesPer
                                case 'present'
                                    data(count).peri_decision = -1;
                                case 'absent'
                                    data(count).peri_decision = 1;
                            end
                            
                            data(count).peri_confid = ...
                                TR(trial).peripheral_confid * data(count).peri_decision;
                            
                            data(count).pSOA = TR(trial).pSOA;
                            data(count).pCON = TR(trial).pCON_divert;
                            data(count).pRT = TR(trial).peripheral_RT;
                            
                        end
                    end
                end
            end
        end
        
        %% TRIAL DATA
        for runner = 1:6
            
            runnum = num2str(runner);
            
            if exist([raw_data files{subj} '/' files{subj} '_block' runnum '.mat'],'file')
                load([raw_data files{subj} '/' files{subj} '_block' runnum '.mat'],'run_file')
                
                for block = 1:length(run_file)
                    
                    TR = run_file(block).TR;
                    
                    for trial = 1:length(TR)
                        count = count + 1;
                        
                        data(count).subjID = files{subj};
                        data(count).group = groupid;
                        data(count).subjnum = subj;
                        data(count).block = (run_file(block).run-1)*6 + ...
                            run_file(block).block; % 1 through 36
                        data(count).task = run_file(block).condition;
                        
                        % 'lo_full' 'med_full' 'hi_full' 'lo_divert' 'med_divert' 'hi_divert'
                        switch run_file(block).condition
                            case 'lo_full'
                                data(count).condition = 1;
                            case 'med_full'
                                data(count).condition = 2;
                            case 'hi_full'
                                data(count).condition = 3;
                            case 'lo_divert'
                                data(count).condition = 4;
                            case 'med_divert'
                                data(count).condition = 5;
                            case 'hi_divert'
                                data(count).condition = 6;
                        end
                        
                        if data(count).condition <= 3
                            % Full attention, no central task
                            data(count).cSOA = TR(trial).cSOA; % Sanity check
                            
                            data(count).peri_accuracy = TR(trial).peripheral_accuracy;
                            
                            switch TR(trial).peripheral_type
                                case 'present'
                                    data(count).peri_signal = -1;
                                case 'absent'
                                    data(count).peri_signal = 1;
                            end
                            
                            switch TR(trial).mouseResponsesPer
                                case 'present'
                                    data(count).peri_decision = -1;
                                case 'absent'
                                    data(count).peri_decision = 1;
                            end
                            
                            data(count).peri_confid = ...
                                TR(trial).peripheral_confid * data(count).peri_decision;
                            
                            data(count).pSOA = TR(trial).pSOA;
                            data(count).pCON = TR(trial).pCON;
                            data(count).pRT = TR(trial).peripheral_RT;
                            
                        elseif data(count).condition >= 4
                            % Diverted attention, both tasks
                            % Central
                            data(count).central_accuracy = TR(trial).central_accuracy;
                            
                            switch TR(trial).letter_condition
                                case 'Ls+T'
                                    data(count).central_signal = -1;
                                case 'Ls'
                                    data(count).central_signal = 1;
                            end
                            
                            switch TR(trial).mouseResponsesCent
                                case 'Ls+T'
                                    data(count).central_decision = -1;
                                case 'Ls'
                                    data(count).central_decision = 1;
                            end
                            
                            data(count).central_confid = ...
                                TR(trial).central_confid * data(count).central_decision;
                            data(count).cSOA = TR(trial).cSOA;
                            
                            % Peripheral
                            data(count).peri_accuracy = TR(trial).peripheral_accuracy;
                            
                            switch TR(trial).peripheral_type
                                case 'present'
                                    data(count).peri_signal = -1;
                                case 'absent'
                                    data(count).peri_signal = 1;
                            end
                            
                            switch TR(trial).mouseResponsesPer
                                case 'present'
                                    data(count).peri_decision = -1;
                                case 'absent'
                                    data(count).peri_decision = 1;
                            end
                            
                            data(count).peri_confid = ...
                                TR(trial).peripheral_confid * data(count).peri_decision;
                            
                            data(count).pSOA = TR(trial).pSOA;
                            data(count).pCON = TR(trial).pCON;
                            data(count).pRT = TR(trial).peripheral_RT;
                            
                        end
                    end
                end
            end
        end
    end
    
end

%% SAVE TO PREPROCESSED DATA LOCATION
if ~exist(preprocessing_location,'dir')
    mkdir(preprocessing_location);
end

save([preprocessing_location '/preprocessed_data.mat'],'data');

end
