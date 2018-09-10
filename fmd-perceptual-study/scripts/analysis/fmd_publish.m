% FMD-DUAL-TASK ANALYSIS
% Concatenate data into workable format

% Preprocess data from specified raw path with bespoke script
[data, dems] = fmd_preprocessing('../../data/raw','../../data/preprocessed');

alldat = struct2cell(data);
groupings = unique({alldat{2,:,:}});

% Perform separately for each group
for gr = 1:length(groupings)
    
    % Build structured matrix of measures
    cdat = alldat(:,:,strcmp(alldat(2,1,:),groupings(gr)));
    subjs = unique({cdat{1,:,:}}); %#ok<*CCAT1>
    conditions = unique({cdat{5,:,:}},'stable');
    clear COND;
    for condition = 1:length(conditions)
        COND(condition).group = groupings{gr};
        condset = cdat(:,:,strcmp(cdat(5,1,:),conditions(condition)));
        COND(condition).task = conditions{condition};  %#ok<*SAGROW>
        COND(condition).subjs = subjs;
        for subject = 1:length(subjs)
            subset = condset(:,:,strcmp(condset(1,1,:),subjs(subject)));
            blocks = unique([subset{4,:,:}]);
            for block = 1:length(blocks)
                count = 0; % Reset count
                for trial = 1:size(subset,3)
                    if subset{4,1,trial}==blocks(block)
                        count = count + 1;
                        COND(condition).pCON(count,block,subject) = ...
                            subset{17,1,trial};
                        COND(condition).pSOA(count,block,subject) = ...
                            subset{16,1,trial};
                        COND(condition).cSOA(count,block,subject) = ...
                            subset{15,1,trial};
                        if ~isempty(subset{14,1,trial}) % Peripheral data
                            COND(condition).peri_confid(count,block,subject) = ...
                                subset{14,1,trial};
                            COND(condition).peri_decision(count,block,subject) = ...
                                subset{13,1,trial};
                            COND(condition).peri_signal(count,block,subject) = ...
                                subset{12,1,trial};
                            COND(condition).peri_accuracy(count,block,subject) = ...
                                subset{11,1,trial};
                            COND(condition).pRT(count,block,subject) = ...
                                subset{18,1,trial};
                        end
                        if ~isempty(subset{10,1,trial}) % Central data
                            COND(condition).cent_confid(count,block,subject) = ...
                                subset{10,1,trial};
                            COND(condition).cent_decision(count,block,subject) = ...
                                subset{9,1,trial};
                            COND(condition).cent_signal(count,block,subject) = ...
                                subset{8,1,trial};
                            COND(condition).cent_accuracy(count,block,subject) = ...
                                subset{7,1,trial};
                        end
                    end
                end
            end
        end
    end
    
    GROUP(gr).COND = COND;
    
end

%% DEMOGRAPHICS
% Input HADS (and/or other demographics) into GROUP struct

for gr = 1:3
    for subj = 1:length(GROUP(gr).COND(1).subjs)
        for demog = 1:height(dems)
            subject = dems.Subject{demog};
            
            % Find subject and record HADS
            if contains(GROUP(gr).COND(1).subjs{subj},subject(4:end))
                GROUP(gr).DEMS.HADS(subj) = dems.HADS(demog);
                GROUP(gr).DEMS.BDI(subj) = dems.BDI(demog);
                GROUP(gr).DEMS.MEDS(subj) = dems.Medicated(demog);
                GROUP(gr).DEMS.HADSd(subj) = dems.HADS_depression(demog);
                GROUP(gr).DEMS.HADSa(subj) = dems.HADS_anxiety(demog);
            end
        end
    end
end

%% CONTRAST THRESHOLDS
% Compare thresholds for gabor contrast under full and diverted attention.
% Also examine Reaction Times for peripheral task

for gr = 1:length(groupings)
    % Subject x block dimensions
    COND = GROUP(gr).COND;
    subjs = COND(1).subjs;
    
    clear pCON;
    clear pRT;
    
    subnum = length(subjs);
    for condition = 1:length(conditions)
        switch COND(condition).task
            case 'ST:gabor contrast'
                errors = size(COND(condition).pCON,1);
                pCON.full_tr = (reshape(mean(COND(condition).pCON),[],subnum))';
                pCON.full_SEM_tr = (reshape(std(COND(condition).pCON)/sqrt(errors),[],subnum))';
                
                pRT.full_tr = (reshape(mean(COND(condition).pRT),[],subnum))';
                pRT.full_SEM_tr = (reshape(std(COND(condition).pRT)/sqrt(errors),[],subnum))';
                
            case 'DT:gabor contrast'
                errors = size(COND(condition).pCON,1);
                pCON.divert_tr = (reshape(mean(COND(condition).pCON),[],subnum))';
                pCON.divert_SEM_tr = (reshape(std(COND(condition).pCON)/sqrt(errors),[],subnum))';
                
                pRT.divert_tr = (reshape(mean(COND(condition).pRT),[],subnum))';
                pRT.divert_SEM_tr = (reshape(std(COND(condition).pRT)/sqrt(errors),[],subnum))';
            case 'med_full'
                errors = size(COND(condition).pCON,1);
                pCON.full = (reshape(mean(COND(condition).pCON),[],subnum))';
                pCON.full_SEM = (reshape(std(COND(condition).pCON)/sqrt(errors),[],subnum))';
                
                pRT.full = (reshape(mean(COND(condition).pRT),[],subnum))';
                pRT.full_SEM = (reshape(std(COND(condition).pRT)/sqrt(errors),[],subnum))';
            case 'med_divert'
                errors = size(COND(condition).pCON,1);
                pCON.divert = (reshape(mean(COND(condition).pCON),[],subnum))';
                pCON.divert_SEM = (reshape(std(COND(condition).pCON)/sqrt(errors),[],subnum))';
                
                pRT.divert = (reshape(mean(COND(condition).pRT),[],subnum))';
                pRT.divert_SEM = (reshape(std(COND(condition).pRT)/sqrt(errors),[],subnum))';
        end
    end
    GROUP(gr).pCON = pCON;
    GROUP(gr).pRT = pRT;
end

%% Start plottin'

figure

% Controls
subplot(2,3,1)
y = GROUP(1).pCON.full';
e = GROUP(1).pCON.full_SEM';
subjs = GROUP(1).COND(1).subjs; subnum = length(subjs);

c = bone(subnum);
er = errorbar(y,e,'o-');
for subject = 1:subnum
    set(er(subject),'color',c(subject,:))
end
box off
set(gca,'TickDir','out','XTick',1:6)
xlim([0.5,6.5])
ylim([-.05,.75])
ylabel('Full Attention');
% xlabel('Experimental Block');
words = sprintf('CONTROLS (n=%s)',mat2str(subnum)); title(words);

subplot(2,3,4)
y = GROUP(1).pCON.divert';
e = GROUP(1).pCON.divert_SEM';

c = bone(subnum);
er = errorbar(y,e,'o-');
for subject = 1:subnum
    set(er(subject),'color',c(subject,:))
end
box off
set(gca,'TickDir','out','XTick',1:6)
xlim([0.5,6.5])
ylim([-.05,.75])
ylabel('Diverted Attention');
% xlabel('Experimental Block');
% title('Diverted Attention');

% FMD
subplot(2,3,2)
y = GROUP(2).pCON.full';
e = GROUP(2).pCON.full_SEM';
subjs = GROUP(2).COND(1).subjs; subnum = length(subjs);

c = cool(subnum);
er = errorbar(y,e,'o-');
for subject = 1:subnum
    set(er(subject),'color',c(subject,:))
end
box off
set(gca,'TickDir','out','XTick',1:6)
xlim([0.5,6.5])
ylim([-.05,.75])
% ylabel('Gabor Contrast');
% xlabel('Experimental Block');
words = sprintf('FMD (n=%s)',mat2str(subnum)); title(words);

subplot(2,3,5)
y = GROUP(2).pCON.divert';
e = GROUP(2).pCON.divert_SEM';

c = cool(subnum);
er = errorbar(y,e,'o-');
for subject = 1:subnum
    set(er(subject),'color',c(subject,:))
end
box off
set(gca,'TickDir','out','XTick',1:6)
xlim([0.5,6.5])
ylim([-.05,.75])
% ylabel('Gabor Contrast');
xlabel('Experimental Block');
% title('Diverted Attention');

% Organic
subplot(2,3,3)
y = GROUP(3).pCON.full';
e = GROUP(3).pCON.full_SEM';
subjs = GROUP(3).COND(1).subjs; subnum = length(subjs);

c = copper(subnum);
er = errorbar(y,e,'o-');
for subject = 1:subnum
    set(er(subject),'color',c(subject,:))
end
box off
set(gca,'TickDir','out','XTick',1:6)
xlim([0.5,6.5])
ylim([-.05,.75])
% ylabel('Gabor Contrast');
% xlabel('Experimental Block');
words = sprintf('ORGANIC (n=%s)',mat2str(subnum)); title(words);

subplot(2,3,6)
y = GROUP(3).pCON.divert';
e = GROUP(3).pCON.divert_SEM';

c = copper(subnum);
er = errorbar(y,e,'o-');
for subject = 1:subnum
    set(er(subject),'color',c(subject,:))
end
box off
set(gca,'TickDir','out','XTick',1:6)
xlim([0.5,6.5])
ylim([-.05,.75])
% ylabel('Gabor Contrast');
% xlabel('Experimental Block');
% title('Diverted Attention');

%% pCON statistics

for gr = 1:length(groupings)
    [h(gr),p(gr),ci(:,gr),stat{gr}]=ttest(mean(GROUP(gr).pCON.full,2),mean(GROUP(gr).pCON.divert,2));
end

disp(stat)

%% BETWEEN-GROUPS stats

% CNTRLS = [mean(GROUP(1).pCON.full,2);mean(GROUP(1).pCON.divert,2)];
% FMDs = [mean(GROUP(2).pCON.full,2);mean(GROUP(2).pCON.divert,2)];
% OMDs = [mean(GROUP(3).pCON.full,2);mean(GROUP(3).pCON.divert,2)];

CNTRLS = [mean(GROUP(1).pRT.full,2);mean(GROUP(1).pRT.divert,2)];
FMDs = [mean(GROUP(2).pRT.full,2);mean(GROUP(2).pRT.divert,2)];
OMDs = [mean(GROUP(3).pRT.full,2);mean(GROUP(3).pRT.divert,2)];

% Controls vs. FMD
[p,tbl,stats] = anova2([CNTRLS,FMDs,OMDs],20,'off');
disp(p)
disp(tbl)
c = multcompare(stats)

% And confirming with effect sizes:

treatment = reshape(repmat(1:3,40,1),[],1); % Groups 
attentions = reshape(repmat(1:2,20,3),[],1); % Attention
fectsiz = reshape([CNTRLS,FMDs,OMDs],[],1);

mes2way(fectsiz,[treatment attentions],'partialeta2'); % Significance testing

% Significant difference in contrast between-groups:
% F(2)=7.30, p=.001, pu2=0.11

% Tukey-Kramer adjusted post-hoc reveals sig difference between:
% Controls & FMD: p=0.011
% Controls & OMD: p=0.001
% But not, FMD & OMD: p=.791

%% pCON as BAR PLOTS

figure

for gr = 1:length(groupings)
    
    subnum = length(GROUP(gr).COND(1).subjs);
    
    % Measures [full and divert]
    measure = mean(GROUP(gr).pCON.full,2);
    FULL_mean(gr) = mean(measure); FULL_SEM(gr) = std(measure)/sqrt(subnum);
    measure2 = mean(GROUP(gr).pCON.divert,2);
    DIV_mean(gr) = mean(measure2); DIV_SEM(gr) = std(measure2)/sqrt(subnum);
    
    tempWSE(gr,:)=within_subject_error([measure measure2]);
    
end

% err = [FULL_SEM;DIV_SEM]';
err = tempWSE*1.9;

colormap(gray)
h = bar([FULL_mean;DIV_mean]',1);
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,1),'r','LineWidth',2)

ylim([0,.3])
xticklabels({'Controls' 'FMD' 'Organics'})
box off
set(gca,'Tickdir','out')

title('Mean Grill Contrast: Group x Attention')

legend({'full attention' 'diverted attention'}, 'Location','best')
legend('boxoff')

%% DOES SENSITIVITY CHANGE OVER BLOCKS?

for GR = 1:3 
    % Examine training sensitivity
    CD = 1;
    
    for subj = 1:size(GROUP(GR).COND(CD).subjs,2)
        for block = 1:size(GROUP(GR).COND(CD).peri_signal,2)
            [~,~,~,training(subj,block,GR)]=type1auc(GROUP(GR).COND(CD).peri_signal(:,block,subj),GROUP(GR).COND(CD).peri_decision(:,block,subj));
        end
    end
    
    % Examine experimental sensitivity
    CDm = 8; % med_full
    CDh = 4; % hi_full
    CDl = 6; % lo_full
    
    for subj = 1:size(GROUP(GR).COND(CD).subjs,2)
        for block = 1:6
            [~,~,~,exp25(subj,block,GR)]=type1auc(GROUP(GR).COND(CDl).peri_signal(:,block,subj),GROUP(GR).COND(CDl).peri_decision(:,block,subj));
            [~,~,~,exp50(subj,block,GR)]=type1auc(GROUP(GR).COND(CDm).peri_signal(:,block,subj),GROUP(GR).COND(CDm).peri_decision(:,block,subj));
            [~,~,~,exp75(subj,block,GR)]=type1auc(GROUP(GR).COND(CDh).peri_signal(:,block,subj),GROUP(GR).COND(CDh).peri_decision(:,block,subj));
        end
    end
    
end

%% IS THERE A CORRELATION BETWEEN pCON & HADS?

for gr = 1:3
    pCONs(:,gr) = mean(GROUP(gr).pCON.divert,2);
    HADS(:,gr) = GROUP(gr).DEMS.HADSd;
end

plot(pCONs(:,1),HADS(:,1),'ko')
hold
plot(pCONs(:,2),HADS(:,2),'rs')
plot(pCONs(:,3),HADS(:,3),'b^')

[rho,pval]=corr(pCONs,HADS,'Tail','both');

% Full attention: 
%       Control: rho=-.10, p>.25; FMD: rho=-.14, p>.25; OMD: rho=-.12, p>.25
% Diverted attention:
%       Control: rho=.24, p>.25; FMD: rho=-.14, p>.25; OMD: rho=-.10, p>.25

% DEMOGRAPHICS

[~,~,stats]=anova1(HADS);
mes1way(HADS,'partialeta2')

% Significant difference in HADS (depression subscale) between groups:
% F(2,57)=6.66, p=.003, p-eta2=.19

% Control & FMD sig (p=.002)
% Control & OMD sig (p=.040)
% FMD & OMD ns (p>.25)

% However, correlation is not significant

% All subjects
corr(reshape(pCONs,[],1),reshape(HADS,[],1))

% Full attention: rho=.02, p>.25
% Diverted attention: rho=.07, p>.25

%% DOES MEDICATION AFFECT CONTRAST PERCEPTION?

pCONs = [mean(GROUP(1).pCON.divert,2);mean(GROUP(2).pCON.divert,2);mean(GROUP(3).pCON.divert,2)];
MEDs = [GROUP(1).DEMS.MEDS,GROUP(2).DEMS.MEDS,GROUP(3).DEMS.MEDS];
HADS = [GROUP(1).DEMS.HADS,GROUP(2).DEMS.HADS,GROUP(3).DEMS.HADS];
GROUPS = reshape(repmat(1:3,20,1),[],1);

plot(pCONs(MEDs==0),HADS(MEDs==0),'ko')
hold on
plot(pCONs(MEDs==1),HADS(MEDs==1),'k.')

h1 = histogram(pCONs(MEDs==1));
hold on 
h2 = histogram(pCONs(MEDs==0));

% h1.Normalization = 'probability';
h1.BinWidth = 0.01;
h2.BinWidth = 0.01;

%% MEASURES

met_thresh = 4;

for gr = 1:length(groupings)
    clear EXP
    for expec = 1:3
        
        switch expec
            case 1
                EXP(expec).level = 75;
                full_expectation = 'hi_full';
                div_expectation = 'hi_divert';
            case 2
                EXP(expec).level = 50;
                full_expectation = 'med_full';
                div_expectation = 'med_divert';
            case 3
                EXP(expec).level = 25;
                full_expectation = 'lo_full';
                div_expectation = 'lo_divert';
        end
        
        for taco = 1:length(GROUP(gr).COND)
            if strcmp(GROUP(gr).COND(taco).task,full_expectation)
                single_tasks = GROUP(gr).COND(taco);
            end
            if strcmp(GROUP(gr).COND(taco).task,div_expectation)
                dual_tasks = GROUP(gr).COND(taco);
            end
        end

        % Calculate measures
        subnum = length(single_tasks.subjs);
        for subj = 1:subnum
            
            if expec == 2
                % CENTRAL
                GROUP(gr).CST.Accuracy(subj)=mean(mean(GROUP(gr).COND(2).cent_accuracy(:,:,subj)));
                
                signal = reshape(GROUP(gr).COND(2).cent_signal(:,:,subj),[],1);
                decision = reshape(GROUP(gr).COND(2).cent_confid(:,:,subj),[],1);
                [~,~,GROUP(gr).CST.OBJPER(subj)...
                    ,GROUP(gr).CST.dPRIME(subj),GROUP(gr).CST.CRIT(subj)...
                    ] = type1auc(signal,decision);
                
                correct = reshape(GROUP(gr).COND(2).cent_accuracy(:,:,subj),[],1);
                conf = reshape(GROUP(gr).COND(2).cent_confid(:,:,subj),[],1);
                [GROUP(gr).CST.METACOG(subj),GROUP(gr).CST.METACrit(subj)]...
                    = type2dprime(correct,conf,met_thresh); % type2roc(correct,abs(conf),4);
                [GROUP(gr).CST.METPrs(subj),GROUP(gr).CST.METCritPrs(subj)]... 
                    = type2dprime(correct(conf(:,:)<0),conf(conf(:,:)<0),met_thresh); % type2roc(correct(conf(:,:)<0),abs(conf(conf(:,:)<0)),4);
                [GROUP(gr).CST.METAbs(subj),GROUP(gr).CST.METCritAbs(subj)]...
                    = type2dprime(correct(conf(:,:)>0),conf(conf(:,:)>0),met_thresh); % type2roc(correct(conf(:,:)>0),abs(conf(conf(:,:)>0)),4);
                
                GROUP(gr).CST.CONFID(subj) = mean(conf);
                GROUP(gr).CST.CONFID_correct(subj) = mean(conf(correct==1)); GROUP(gr).CST.CONFID_cn(subj) = size(conf(correct==1),1);
                GROUP(gr).CST.CONFID_incorrect(subj) = mean(conf(correct==0)); GROUP(gr).CST.CONFID_in(subj) = size(conf(correct==0),1);
            end
            
            %% SINGLE
            
            EXP(expec).PST.Accuracy(subj)=mean(mean(single_tasks.peri_accuracy(:,:,subj)));
            
            signal = reshape(single_tasks.peri_signal(:,:,subj),[],1);
            decision = reshape(single_tasks.peri_confid(:,:,subj),[],1);
            [~,~,EXP(expec).PST.OBJPER(subj)...
                    ,EXP(expec).PST.dPRIME(subj),EXP(expec).PST.CRIT(subj)...
                    ] = type1auc(signal,decision);
            
            correct = reshape(single_tasks.peri_accuracy(:,:,subj),[],1);
            conf = reshape(single_tasks.peri_confid(:,:,subj),[],1);
            
            if ~isempty(abs(conf(conf(:,:)<0)))
                EXP(expec).PST.METAUCPrs(subj) = type2roc(correct(conf(:,:)<0),abs(conf(conf(:,:)<0)),4);
                EXP(expec).PST.METAUC2Prs(subj) = type2auc(correct(conf(:,:)<0),abs(conf(conf(:,:)<0)));
            else
                EXP(expec).PST.METAUCPrs(subj) = NaN;
                EXP(expec).PST.METAUC2Prs(subj) = NaN;
            end
            
            
            
            if ~isempty(abs(conf(conf(:,:)>0)))
                EXP(expec).PST.METAUCAbs(subj) = type2roc(correct(conf(:,:)>0),abs(conf(conf(:,:)>0)),4);
                EXP(expec).PST.METAUC2Abs(subj) = type2auc(correct(conf(:,:)>0),abs(conf(conf(:,:)>0)));
            else
                EXP(expec).PST.METAUCAbs(subj) = NaN;
                EXP(expec).PST.METAUC2Abs(subj) = NaN;
            end
            
            
            
            [EXP(expec).PST.METACOG(subj),EXP(expec).PST.METACrit(subj)]...
                = type2dprime(correct,conf,met_thresh); % type2roc(correct,abs(conf),4);
            [EXP(expec).PST.METPrs(subj),EXP(expec).PST.METCritPrs(subj)]...
                = type2dprime(correct(conf(:,:)<0),conf(conf(:,:)<0),met_thresh); % type2roc(correct(conf(:,:)<0),abs(conf(conf(:,:)<0)),4);
            [EXP(expec).PST.METAbs(subj),EXP(expec).PST.METCritAbs(subj)]...
                = type2dprime(correct(conf(:,:)>0),conf(conf(:,:)>0),met_thresh); % type2roc(correct(conf(:,:)>0),abs(conf(conf(:,:)>0)),4);
            
            conf = abs(conf);
            EXP(expec).PST.CONFID(subj) = mean(conf);
            EXP(expec).PST.CONFID_correct(subj) = mean(conf(correct==1)); EXP(expec).PST.CONFID_cn(subj) = size(conf(correct==1),1);
            EXP(expec).PST.CONFID_incorrect(subj) = mean(conf(correct==0)); EXP(expec).PST.CONFID_in(subj) = size(conf(correct==0),1);
            
            pres=decision(signal<0); abse=decision(signal>0);
            
            % Present "it's present!"
            pp=pres(pres<0); EXP(expec).PST.CONFsd_pp(subj) = abs(mean(pp));
            % Present "it's absent!"
            pa=pres(pres>0); EXP(expec).PST.CONFsd_pa(subj) = abs(mean(pa));
            % Absent "it's present!"
            ap=abse(abse<0); EXP(expec).PST.CONFsd_ap(subj) = abs(mean(ap));
            % absent "it's absent!"
            aa=abse(abse>0); EXP(expec).PST.CONFsd_aa(subj) = abs(mean(aa));
            %% DUAL
            
            EXP(expec).PDT.Accuracy(subj)=mean(mean(dual_tasks.peri_accuracy(:,:,subj)));
            
            signal = reshape(dual_tasks.peri_signal(:,:,subj),[],1);
            decision = reshape(dual_tasks.peri_confid(:,:,subj),[],1);
            [~,~,EXP(expec).PDT.OBJPER(subj)...
                    ,EXP(expec).PDT.dPRIME(subj),EXP(expec).PDT.CRIT(subj)...
                    ] = type1auc(signal,decision);
            
            correct = reshape(dual_tasks.peri_accuracy(:,:,subj),[],1);
            conf = reshape(dual_tasks.peri_confid(:,:,subj),[],1);
            
            if ~isempty(abs(conf(conf(:,:)<0)))
                EXP(expec).PDT.METAUCPrs(subj) = type2roc(correct(conf(:,:)<0),abs(conf(conf(:,:)<0)),4);
                EXP(expec).PDT.METAUC2Prs(subj) = type2auc(correct(conf(:,:)<0),abs(conf(conf(:,:)<0)));
            else
                EXP(expec).PDT.METAUCPrs(subj) = NaN;
                EXP(expec).PDT.METAUC2Prs(subj) = NaN;
            end
            
            if ~isempty(abs(conf(conf(:,:)>0)))
                EXP(expec).PDT.METAUCAbs(subj) = type2roc(correct(conf(:,:)>0),abs(conf(conf(:,:)>0)),4);
                EXP(expec).PDT.METAUC2Abs(subj) = type2auc(correct(conf(:,:)>0),abs(conf(conf(:,:)>0)));
            else
                EXP(expec).PDT.METAUCAbs(subj) = NaN;
                EXP(expec).PDT.METAUC2Abs(subj) = NaN;
            end
            
            [EXP(expec).PDT.METACOG(subj),EXP(expec).PDT.METACrit(subj)]...
                = type2dprime(correct, conf, met_thresh); % type2roc(correct,abs(conf),4);
            [EXP(expec).PDT.METPrs(subj),EXP(expec).PDT.METCritPrs(subj)]...
                = type2dprime(correct(conf(:,:)<0),conf(conf(:,:)<0), met_thresh); % type2roc(correct(conf(:,:)<0),abs(conf(conf(:,:)<0)),4);
            [EXP(expec).PDT.METAbs(subj),EXP(expec).PDT.METCritAbs(subj)]...
                = type2dprime(correct(conf(:,:)>0),conf(conf(:,:)>0), met_thresh); % type2roc(correct(conf(:,:)>0),abs(conf(conf(:,:)>0)),4);
            
            conf = abs(conf);
            EXP(expec).PDT.CONFID(subj) = mean(conf);
            EXP(expec).PDT.CONFID_correct(subj) = mean(conf(correct==1)); EXP(expec).PDT.CONFID_cn(subj) = size(conf(correct==1),1);
            EXP(expec).PDT.CONFID_incorrect(subj) = mean(conf(correct==0)); EXP(expec).PDT.CONFID_in(subj) = size(conf(correct==0),1);
            
            pres=decision(signal<0); abse=decision(signal>0);
            
            % Present "it's present!"
            pp=pres(pres<0); EXP(expec).PDT.CONFsd_pp(subj) = abs(mean(pp));
            % Present "it's absent!"
            pa=pres(pres>0); EXP(expec).PDT.CONFsd_pa(subj) = abs(mean(pa));
            % Absent "it's present!"
            ap=abse(abse<0); EXP(expec).PDT.CONFsd_ap(subj) = abs(mean(ap));
            % absent "it's absent!"
            aa=abse(abse>0); EXP(expec).PDT.CONFsd_aa(subj) = abs(mean(aa));
            
            EXP(expec).CDT.Accuracy(subj)=mean(mean(dual_tasks.cent_accuracy(:,:,subj)));
            
            signal = reshape(dual_tasks.cent_signal(:,:,subj),[],1);
            decision = reshape(dual_tasks.cent_confid(:,:,subj),[],1);
            [~,~,EXP(expec).CDT.OBJPER(subj)] = type1auc(signal,decision);
            [~,~,EXP(expec).CDT.OBJPER(subj)...
                    ,EXP(expec).CDT.dPRIME(subj),EXP(expec).CDT.CRIT(subj)...
                    ] = type1auc(signal,decision);
            
            correct = reshape(dual_tasks.cent_accuracy(:,:,subj),[],1);
            conf = reshape(dual_tasks.cent_confid(:,:,subj),[],1);
            
            if ~isempty(abs(conf(conf(:,:)<0)))
                EXP(expec).CDT.METAUCPrs(subj) = type2roc(correct(conf(:,:)<0),abs(conf(conf(:,:)<0)),4);
            else
                EXP(expec).CDT.METAUCPrs(subj) = NaN;
            end
            
            if ~isempty(abs(conf(conf(:,:)>0)))
                EXP(expec).CDT.METAUCAbs(subj) = type2roc(correct(conf(:,:)>0),abs(conf(conf(:,:)>0)),4);
            else
                EXP(expec).CDT.METAUCAbs(subj) = NaN;
            end
            
            [EXP(expec).CDT.METACOG(subj),EXP(expec).CDT.METACrit(subj)]...
                = type2dprime(correct,conf, met_thresh); % type2roc(correct,abs(conf),4);
            [EXP(expec).CDT.METPrs(subj),EXP(expec).CDT.METCritPrs(subj)]...
                = type2dprime(correct(conf(:,:)<0),conf(conf(:,:)<0), met_thresh); % type2roc(correct(conf(:,:)<0),abs(conf(conf(:,:)<0)),4);
            [EXP(expec).CDT.METAbs(subj),EXP(expec).CDT.METCritAbs(subj)]...
                = type2dprime(correct(conf(:,:)>0),conf(conf(:,:)>0), met_thresh); % type2roc(correct(conf(:,:)>0),abs(conf(conf(:,:)>0)),4);
            
            conf = abs(conf);
            EXP(expec).CDT.CONFID(subj) = mean(conf);
            EXP(expec).CDT.CONFID_correct(subj) = mean(conf(correct==1)); EXP(expec).CDT.CONFID_cn(subj) = size(conf(correct==1),1);
            EXP(expec).CDT.CONFID_incorrect(subj) = mean(conf(correct==0)); EXP(expec).CDT.CONFID_in(subj) = size(conf(correct==0),1);
            
        end
    end
    GROUP(gr).EXP = EXP;
end

%% PLOT ACCURACY

ymits = [.5 1];

% 9 plots, x-axes are group (CONTROL, FMD, ORGANIC), y-axes are expectations (HIGH, MED, LOW)
figure

subplot(3,1,1)
gr = 1;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.Accuracy;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.Accuracy;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.Accuracy;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.Accuracy;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.Accuracy;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.Accuracy;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

% tempWSE=within_subject_error(...
%     [GROUP(gr).EXP(3).PST.Accuracy;GROUP(gr).EXP(2).PST.Accuracy;GROUP(gr).EXP(1).PST.Accuracy;...
%     GROUP(gr).EXP(3).PDT.Accuracy;GROUP(gr).EXP(2).PDT.Accuracy;GROUP(gr).EXP(1).PDT.Accuracy;...
%     ]',[2,3]);
% PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

ylim(ymits)
xticklabels({'25%' '50%' '75%'})
% xlabel('Expectations')
box off
set(gca,'Tickdir','out')
words = sprintf('CONTROL (n=%s)',mat2str(subnum)); ylabel(words);
title('Task Accuracy')

legend({'full attention' 'diverted attention'})
legend('boxoff')

subplot(3,1,2)
gr = 2;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.Accuracy;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.Accuracy;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.Accuracy;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.Accuracy;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.Accuracy;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.Accuracy;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

% tempWSE=within_subject_error(...
%     [GROUP(gr).EXP(3).PST.Accuracy;GROUP(gr).EXP(2).PST.Accuracy;GROUP(gr).EXP(1).PST.Accuracy;...
%     GROUP(gr).EXP(3).PDT.Accuracy;GROUP(gr).EXP(2).PDT.Accuracy;GROUP(gr).EXP(1).PDT.Accuracy;...
%     ]',[2,3]);
% PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

ylim(ymits)
xticklabels({'25%' '50%' '75%'})
% xlabel('Expectations')
box off
set(gca,'Tickdir','out')
words = sprintf('FMD (n=%s)',mat2str(subnum)); ylabel(words);

subplot(3,1,3)
gr = 3;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.Accuracy;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.Accuracy;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.Accuracy;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.Accuracy;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.Accuracy;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.Accuracy;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

% tempWSE=within_subject_error(...
%     [GROUP(gr).EXP(3).PST.Accuracy;GROUP(gr).EXP(2).PST.Accuracy;GROUP(gr).EXP(1).PST.Accuracy;...
%     GROUP(gr).EXP(3).PDT.Accuracy;GROUP(gr).EXP(2).PDT.Accuracy;GROUP(gr).EXP(1).PDT.Accuracy;...
%     ]',[2,3]);
% PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)
% sigstar({[1,3]},.005)
ylim(ymits)
xticklabels({'25%' '50%' '75%'})
xlabel('Expectations of Target Presence')
box off
set(gca,'Tickdir','out')
words = sprintf('ORGANIC (n=%s)',mat2str(subnum)); ylabel(words);

%% ACCURACY STATISTICS

EXP = GROUP(1).EXP; subnum = length(EXP(1).PDT.Accuracy);
all_subjs = [EXP(1).PDT.Accuracy,EXP(1).PST.Accuracy;...
    EXP(2).PDT.Accuracy,EXP(2).PST.Accuracy;...
    EXP(3).PDT.Accuracy,EXP(3).PST.Accuracy]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

EXP = GROUP(2).EXP; subnum = length(EXP(1).PDT.Accuracy);
all_subjs = [EXP(1).PDT.Accuracy,EXP(1).PST.Accuracy;...
    EXP(2).PDT.Accuracy,EXP(2).PST.Accuracy;...
    EXP(3).PDT.Accuracy,EXP(3).PST.Accuracy]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

EXP = GROUP(3).EXP; subnum = length(EXP(1).PDT.Accuracy);
all_subjs = [EXP(1).PDT.Accuracy,EXP(1).PST.Accuracy;...
    EXP(2).PDT.Accuracy,EXP(2).PST.Accuracy;...
    EXP(3).PDT.Accuracy,EXP(3).PST.Accuracy]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

%% PLOT MEAN CONFIDENCE

ymits = [.5 4.5];

% 9 plots, x-axes are group (CONTROL, FMD, ORGANIC), y-axes are expectations (HIGH, MED, LOW)
figure

subplot(3,1,1)
gr = 1;

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.CONFID;
PST_mean(1) = nanmean(measure); PST_SEM(1) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(2).PST.CONFID;
PST_mean(2) = nanmean(measure); PST_SEM(2) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(1).PST.CONFID;
PST_mean(3) = nanmean(measure); PST_SEM(3) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));

% Measures ['lo' 'med' 'hi'] PDT CORRECT
measure = GROUP(gr).EXP(3).PDT.CONFID;
PDT_mean(1) = nanmean(measure); PDT_SEM(1) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(2).PDT.CONFID;
PDT_mean(2) = nanmean(measure); PDT_SEM(2) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(1).PDT.CONFID;
PDT_mean(3) = nanmean(measure); PDT_SEM(3) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

ylim(ymits)
xticklabels({'25%' '50%' '75%'})

box off
set(gca,'Tickdir','out')
words = sprintf('CONTROL (n=%s)',mat2str(subnum)); ylabel(words);
title('Mean Confidence')

legend({'full attention' 'diverted attention'})
legend('boxoff')

subplot(3,1,2)
gr = 2;

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.CONFID;
PST_mean(1) = nanmean(measure); PST_SEM(1) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(2).PST.CONFID;
PST_mean(2) = nanmean(measure); PST_SEM(2) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(1).PST.CONFID;
PST_mean(3) = nanmean(measure); PST_SEM(3) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));

% Measures ['lo' 'med' 'hi'] PDT CORRECT
measure = GROUP(gr).EXP(3).PDT.CONFID;
PDT_mean(1) = nanmean(measure); PDT_SEM(1) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(2).PDT.CONFID;
PDT_mean(2) = nanmean(measure); PDT_SEM(2) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(1).PDT.CONFID;
PDT_mean(3) = nanmean(measure); PDT_SEM(3) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

ylim(ymits)
xticklabels({'25%' '50%' '75%'})
box off
set(gca,'Tickdir','out')
words = sprintf('FMD (n=%s)',mat2str(subnum)); ylabel(words);

subplot(3,1,3)
gr = 3;

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.CONFID;
PST_mean(1) = nanmean(measure); PST_SEM(1) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(2).PST.CONFID;
PST_mean(2) = nanmean(measure); PST_SEM(2) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(1).PST.CONFID;
PST_mean(3) = nanmean(measure); PST_SEM(3) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));

% Measures ['lo' 'med' 'hi'] PDT CORRECT
measure = GROUP(gr).EXP(3).PDT.CONFID;
PDT_mean(1) = nanmean(measure); PDT_SEM(1) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(2).PDT.CONFID;
PDT_mean(2) = nanmean(measure); PDT_SEM(2) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(1).PDT.CONFID;
PDT_mean(3) = nanmean(measure); PDT_SEM(3) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

ylim(ymits)
xticklabels({'25%' '50%' '75%'})
xlabel('Expectations of Target Presence')
box off
set(gca,'Tickdir','out')
words = sprintf('ORGANIC (n=%s)',mat2str(subnum)); ylabel(words);

%% CONFIDENCE STATISTICS

EXP = GROUP(1).EXP; subnum = length(EXP(1).PDT.CONFID);
all_subjs = [EXP(1).PDT.CONFID,EXP(1).PST.CONFID;...
    EXP(2).PDT.CONFID,EXP(2).PST.CONFID;...
    EXP(3).PDT.CONFID,EXP(3).PST.CONFID]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

EXP = GROUP(2).EXP; subnum = length(EXP(1).PDT.CONFID);
all_subjs = [EXP(1).PDT.CONFID,EXP(1).PST.CONFID;...
    EXP(2).PDT.CONFID,EXP(2).PST.CONFID;...
    EXP(3).PDT.CONFID,EXP(3).PST.CONFID]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

EXP = GROUP(3).EXP; subnum = length(EXP(1).PDT.CONFID);
all_subjs = [EXP(1).PDT.CONFID,EXP(1).PST.CONFID;...
    EXP(2).PDT.CONFID,EXP(2).PST.CONFID;...
    EXP(3).PDT.CONFID,EXP(3).PST.CONFID]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

%% PLOT C/I CONFIDENCE

ymits = [.5 4.5];

% 9 plots, x-axes are group (CONTROL, FMD, ORGANIC), y-axes are expectations (HIGH, MED, LOW)
figure

subplot(3,1,1)
gr = 1;

% Measures ['lo' 'med' 'hi'] PST CORRECT
measure = GROUP(gr).EXP(3).PST.CONFID_correct;
PSTc_mean(1) = nanmean(measure); PSTc_SEM(1) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(2).PST.CONFID_correct;
PSTc_mean(2) = nanmean(measure); PSTc_SEM(2) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(1).PST.CONFID_correct;
PSTc_mean(3) = nanmean(measure); PSTc_SEM(3) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));

% Measures ['lo' 'med' 'hi'] PST INCORRECT
measure = GROUP(gr).EXP(3).PST.CONFID_incorrect;
PSTi_mean(1) = nanmean(measure); PSTi_SEM(1) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(2).PST.CONFID_incorrect;
PSTi_mean(2) = nanmean(measure); PSTi_SEM(2) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(1).PST.CONFID_incorrect;
PSTi_mean(3) = nanmean(measure); PSTi_SEM(3) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));

% Measures ['lo' 'med' 'hi'] PDT CORRECT
measure = GROUP(gr).EXP(3).PDT.CONFID_correct;
PDTc_mean(1) = nanmean(measure); PDTc_SEM(1) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(2).PDT.CONFID_correct;
PDTc_mean(2) = nanmean(measure); PDTc_SEM(2) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(1).PDT.CONFID_correct;
PDTc_mean(3) = nanmean(measure); PDTc_SEM(3) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));

% Measures ['lo' 'med' 'hi'] PDT INCORRECT
measure = GROUP(gr).EXP(3).PDT.CONFID_incorrect;
PDTi_mean(1) = nanmean(measure); PDTi_SEM(1) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(2).PDT.CONFID_incorrect;
PDTi_mean(2) = nanmean(measure); PDTi_SEM(2) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(1).PDT.CONFID_incorrect;
PDTi_mean(3) = nanmean(measure); PDTi_SEM(3) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));


err = [PSTc_SEM;PSTi_SEM;PDTc_SEM;PDTi_SEM]';
colormap(gray);
h = bar([PSTc_mean;PSTi_mean;PDTc_mean;PDTi_mean]');
hold on
errbar(h(1).XData-.27,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData-.09,h(2).YData,err(:,2),'r','LineWidth',2)
errbar(h(3).XData+.09,h(3).YData,err(:,3),'r','LineWidth',2)
errbar(h(4).XData+.27,h(4).YData,err(:,4),'r','LineWidth',2)

ylim(ymits)
xticklabels({'25%' '50%' '75%'})

box off
set(gca,'Tickdir','out')
words = sprintf('CONTROL (n=%s)',mat2str(subnum)); ylabel(words);
title('Task Confidence (Correct vs. Incorrect)')

legend({'full correct' 'full incorrect' 'diverted correct' 'diverted incorrect'})
legend('boxoff')

subplot(3,1,2)
gr = 2;

% Measures ['lo' 'med' 'hi'] PST CORRECT
measure = GROUP(gr).EXP(3).PST.CONFID_correct;
PSTc_mean(1) = nanmean(measure); PSTc_SEM(1) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(2).PST.CONFID_correct;
PSTc_mean(2) = nanmean(measure); PSTc_SEM(2) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(1).PST.CONFID_correct;
PSTc_mean(3) = nanmean(measure); PSTc_SEM(3) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));

% Measures ['lo' 'med' 'hi'] PST INCORRECT
measure = GROUP(gr).EXP(3).PST.CONFID_incorrect;
PSTi_mean(1) = nanmean(measure); PSTi_SEM(1) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(2).PST.CONFID_incorrect;
PSTi_mean(2) = nanmean(measure); PSTi_SEM(2) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(1).PST.CONFID_incorrect;
PSTi_mean(3) = nanmean(measure); PSTi_SEM(3) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));

% Measures ['lo' 'med' 'hi'] PDT CORRECT
measure = GROUP(gr).EXP(3).PDT.CONFID_correct;
PDTc_mean(1) = nanmean(measure); PDTc_SEM(1) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(2).PDT.CONFID_correct;
PDTc_mean(2) = nanmean(measure); PDTc_SEM(2) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(1).PDT.CONFID_correct;
PDTc_mean(3) = nanmean(measure); PDTc_SEM(3) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));

% Measures ['lo' 'med' 'hi'] PDT INCORRECT
measure = GROUP(gr).EXP(3).PDT.CONFID_incorrect;
PDTi_mean(1) = nanmean(measure); PDTi_SEM(1) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(2).PDT.CONFID_incorrect;
PDTi_mean(2) = nanmean(measure); PDTi_SEM(2) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(1).PDT.CONFID_incorrect;
PDTi_mean(3) = nanmean(measure); PDTi_SEM(3) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));

err = [PSTc_SEM;PSTi_SEM;PDTc_SEM;PDTi_SEM]';
colormap(gray);
h = bar([PSTc_mean;PSTi_mean;PDTc_mean;PDTi_mean]');
hold on
errbar(h(1).XData-.27,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData-.09,h(2).YData,err(:,2),'r','LineWidth',2)
errbar(h(3).XData+.09,h(3).YData,err(:,3),'r','LineWidth',2)
errbar(h(4).XData+.27,h(4).YData,err(:,4),'r','LineWidth',2)

ylim(ymits)
xticklabels({'25%' '50%' '75%'})
box off
set(gca,'Tickdir','out')
words = sprintf('FMD (n=%s)',mat2str(subnum)); ylabel(words);

subplot(3,1,3)
gr = 3;

% Measures ['lo' 'med' 'hi'] PST CORRECT
measure = GROUP(gr).EXP(3).PST.CONFID_correct;
PSTc_mean(1) = nanmean(measure); PSTc_SEM(1) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(2).PST.CONFID_correct;
PSTc_mean(2) = nanmean(measure); PSTc_SEM(2) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(1).PST.CONFID_correct;
PSTc_mean(3) = nanmean(measure); PSTc_SEM(3) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));

% Measures ['lo' 'med' 'hi'] PST INCORRECT
measure = GROUP(gr).EXP(3).PST.CONFID_incorrect;
PSTi_mean(1) = nanmean(measure); PSTi_SEM(1) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(2).PST.CONFID_incorrect;
PSTi_mean(2) = nanmean(measure); PSTi_SEM(2) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(1).PST.CONFID_incorrect;
PSTi_mean(3) = nanmean(measure); PSTi_SEM(3) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));

% Measures ['lo' 'med' 'hi'] PDT CORRECT
measure = GROUP(gr).EXP(3).PDT.CONFID_correct;
PDTc_mean(1) = nanmean(measure); PDTc_SEM(1) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(2).PDT.CONFID_correct;
PDTc_mean(2) = nanmean(measure); PDTc_SEM(2) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(1).PDT.CONFID_correct;
PDTc_mean(3) = nanmean(measure); PDTc_SEM(3) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));

% Measures ['lo' 'med' 'hi'] PDT INCORRECT
measure = GROUP(gr).EXP(3).PDT.CONFID_incorrect;
PDTi_mean(1) = nanmean(measure); PDTi_SEM(1) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(2).PDT.CONFID_incorrect;
PDTi_mean(2) = nanmean(measure); PDTi_SEM(2) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));
measure = GROUP(gr).EXP(1).PDT.CONFID_incorrect;
PDTi_mean(3) = nanmean(measure); PDTi_SEM(3) = nanstd(measure)/sqrt(length(measure)-sum(isnan(measure)));

err = [PSTc_SEM;PSTi_SEM;PDTc_SEM;PDTi_SEM]';
colormap(gray);
h = bar([PSTc_mean;PSTi_mean;PDTc_mean;PDTi_mean]');
hold on
errbar(h(1).XData-.27,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData-.09,h(2).YData,err(:,2),'r','LineWidth',2)
errbar(h(3).XData+.09,h(3).YData,err(:,3),'r','LineWidth',2)
errbar(h(4).XData+.27,h(4).YData,err(:,4),'r','LineWidth',2)
% sigstar({[1,3]},.005)
ylim(ymits)
xticklabels({'25%' '50%' '75%'})
xlabel('Expectations of Target Presence')
box off
set(gca,'Tickdir','out')
words = sprintf('ORGANIC (n=%s)',mat2str(subnum)); ylabel(words);

%% PLOT dPRIME

% 9 plots, x-axes are group (CONTROL, FMD, ORGANIC), y-axes are expectations (HIGH, MED, LOW)
figure

subplot(3,1,1)
gr = 1;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.dPRIME;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.dPRIME;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.dPRIME;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.dPRIME;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.dPRIME;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.dPRIME;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

tempWSE=within_subject_error(...
    [GROUP(gr).EXP(3).PST.dPRIME;GROUP(gr).EXP(2).PST.dPRIME;GROUP(gr).EXP(1).PST.dPRIME;...
    GROUP(gr).EXP(3).PDT.dPRIME;GROUP(gr).EXP(2).PDT.dPRIME;GROUP(gr).EXP(1).PDT.dPRIME;...
    ]',[2,3]);
PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
% xlabel('Expectations')
box off
set(gca,'Tickdir','out')
words = sprintf('CONTROL (n=%s)',mat2str(subnum)); ylabel(words);
title('Objective Performance')

legend({'full attention' 'diverted attention'})
legend('boxoff')

subplot(3,1,2)
gr = 2;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.dPRIME;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.dPRIME;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.dPRIME;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.dPRIME;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.dPRIME;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.dPRIME;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

tempWSE=within_subject_error(...
    [GROUP(gr).EXP(3).PST.dPRIME;GROUP(gr).EXP(2).PST.dPRIME;GROUP(gr).EXP(1).PST.dPRIME;...
    GROUP(gr).EXP(3).PDT.dPRIME;GROUP(gr).EXP(2).PDT.dPRIME;GROUP(gr).EXP(1).PDT.dPRIME;...
    ]',[2,3]);
PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
% xlabel('Expectations')
box off
set(gca,'Tickdir','out')
words = sprintf('FMD (n=%s)',mat2str(subnum)); ylabel(words);

subplot(3,1,3)
gr = 3;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.dPRIME;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.dPRIME;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.dPRIME;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.dPRIME;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.dPRIME;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.dPRIME;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

tempWSE=within_subject_error(...
    [GROUP(gr).EXP(3).PST.dPRIME;GROUP(gr).EXP(2).PST.dPRIME;GROUP(gr).EXP(1).PST.dPRIME;...
    GROUP(gr).EXP(3).PDT.dPRIME;GROUP(gr).EXP(2).PDT.dPRIME;GROUP(gr).EXP(1).PDT.dPRIME;...
    ]',[2,3]);
PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)
sigstar({[1,3]},.005)
ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
xlabel('Expectations of Target Presence')
box off
set(gca,'Tickdir','out')
words = sprintf('ORGANIC (n=%s)',mat2str(subnum)); ylabel(words);

%% dPRIME STATISTICS

EXP = GROUP(1).EXP; subnum = length(EXP(1).PDT.dPRIME);
all_subjs = [EXP(1).PDT.dPRIME,EXP(1).PST.dPRIME;...
    EXP(2).PDT.dPRIME,EXP(2).PST.dPRIME;...
    EXP(3).PDT.dPRIME,EXP(3).PST.dPRIME]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

EXP = GROUP(2).EXP; subnum = length(EXP(1).PDT.dPRIME);
all_subjs = [EXP(1).PDT.dPRIME,EXP(1).PST.dPRIME;...
    EXP(2).PDT.dPRIME,EXP(2).PST.dPRIME;...
    EXP(3).PDT.dPRIME,EXP(3).PST.dPRIME]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

EXP = GROUP(3).EXP; subnum = length(EXP(1).PDT.dPRIME);
all_subjs = [EXP(1).PDT.dPRIME,EXP(1).PST.dPRIME;...
    EXP(2).PDT.dPRIME,EXP(2).PST.dPRIME;...
    EXP(3).PDT.dPRIME,EXP(3).PST.dPRIME]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

%% CRITERION Plot 'em

figure

subplot(3,1,1)
gr = 1;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.CRIT;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.CRIT;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.CRIT;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.CRIT;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.CRIT;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.CRIT;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

tempWSE=within_subject_error(...
    [GROUP(gr).EXP(3).PST.CRIT;GROUP(gr).EXP(2).PST.CRIT;GROUP(gr).EXP(1).PST.CRIT;...
    GROUP(gr).EXP(3).PDT.CRIT;GROUP(gr).EXP(2).PDT.CRIT;GROUP(gr).EXP(1).PDT.CRIT;...
    ]',[2,3]);
PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
sigstar({[1,3]},.0003)
% xlabel('Expectations')
box off
set(gca,'Tickdir','out')
words = sprintf('CONTROL (n=%s)',mat2str(subnum)); ylabel(words);
title('Type 1 criterion')

legend({'full' 'diverted'})
legend('boxoff')

subplot(3,1,2)
gr = 2;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.CRIT;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.CRIT;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.CRIT;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.CRIT;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.CRIT;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.CRIT;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

tempWSE=within_subject_error(...
    [GROUP(gr).EXP(3).PST.CRIT;GROUP(gr).EXP(2).PST.CRIT;GROUP(gr).EXP(1).PST.CRIT;...
    GROUP(gr).EXP(3).PDT.CRIT;GROUP(gr).EXP(2).PDT.CRIT;GROUP(gr).EXP(1).PDT.CRIT;...
    ]',[2,3]);
PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
sigstar({[1,3]},.0049)
% xlabel('Expectations')
box off
set(gca,'Tickdir','out')
words = sprintf('FMD (n=%s)',mat2str(subnum)); ylabel(words);

subplot(3,1,3)
gr = 3;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.CRIT;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.CRIT;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.CRIT;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.CRIT;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.CRIT;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.CRIT;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

tempWSE=within_subject_error(...
    [GROUP(gr).EXP(3).PST.CRIT;GROUP(gr).EXP(2).PST.CRIT;GROUP(gr).EXP(1).PST.CRIT;...
    GROUP(gr).EXP(3).PDT.CRIT;GROUP(gr).EXP(2).PDT.CRIT;GROUP(gr).EXP(1).PDT.CRIT;...
    ]',[2,3]);
PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% sigstar({[1,3]},.0258)
% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
xlabel('Expectations of Target Presence')
box off
set(gca,'Tickdir','out')
words = sprintf('ORGANIC (n=%s)',mat2str(subnum)); ylabel(words);

%% CRITERION STATISTICS

EXP = GROUP(1).EXP; subnum = length(EXP(1).PDT.CRIT);
all_subjs = [EXP(1).PDT.CRIT,EXP(1).PST.CRIT;...
    EXP(2).PDT.CRIT,EXP(2).PST.CRIT;...
    EXP(3).PDT.CRIT,EXP(3).PST.CRIT]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

EXP = GROUP(2).EXP; subnum = length(EXP(1).PDT.CRIT);
all_subjs = [EXP(1).PDT.CRIT,EXP(1).PST.CRIT;...
    EXP(2).PDT.CRIT,EXP(2).PST.CRIT;...
    EXP(3).PDT.CRIT,EXP(3).PST.CRIT]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

EXP = GROUP(3).EXP; subnum = length(EXP(1).PDT.CRIT);
all_subjs = [EXP(1).PDT.CRIT,EXP(1).PST.CRIT;...
    EXP(2).PDT.CRIT,EXP(2).PST.CRIT;...
    EXP(3).PDT.CRIT,EXP(3).PST.CRIT]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

%% BETWEEN GROUPS ANALYSIS

for gr = 1:3
    EXP = GROUP(gr).EXP; subnum = length(EXP(1).PDT.CRIT);
    THIS(gr).all_subjs = [EXP(1).PDT.CRIT,EXP(1).PST.CRIT;...
        EXP(2).PDT.CRIT,EXP(2).PST.CRIT;...
        EXP(3).PDT.CRIT,EXP(3).PST.CRIT]';
    
    THIS(gr).g1 = repmat([repmat({'25%'},subnum*2,1);repmat({'50%'},subnum*2,1);repmat({'75%'},subnum*2,1)],1,1); % EXPECTATIONS
    
    THIS(gr).g2 = repmat([repmat({'D'},subnum,1);repmat({'F'},subnum,1)],3,1); % ATTENTION
    
    THIS(gr).g3 = repmat(mat2str(gr),subnum,1); % TREATMENT
    
end

testing = [reshape(THIS(1).all_subjs,[],1);reshape(THIS(2).all_subjs,[],1);reshape(THIS(3).all_subjs,[],1)];
g1 = [THIS(1).g1; THIS(2).g1; THIS(3).g1];
g2 = [THIS(1).g2; THIS(2).g2; THIS(3).g2];
re3 = [repmat(THIS(1).g3,6,1); repmat(THIS(2).g3,6,1); repmat(THIS(3).g3,6,1)];

tested = 'confmd';
if strcmp(tested,'all')
    g3 = [THIS(1).g3; THIS(2).g3; THIS(3).g3];
    
    t = table(g3,...
        [GROUP(1).EXP(1).PST.CRIT,GROUP(2).EXP(1).PST.CRIT,GROUP(3).EXP(1).PST.CRIT]',...
        [GROUP(1).EXP(2).PST.CRIT,GROUP(2).EXP(2).PST.CRIT,GROUP(3).EXP(2).PST.CRIT]',...
        [GROUP(1).EXP(3).PST.CRIT,GROUP(2).EXP(3).PST.CRIT,GROUP(3).EXP(3).PST.CRIT]',...
        [GROUP(1).EXP(1).PDT.CRIT,GROUP(2).EXP(1).PDT.CRIT,GROUP(3).EXP(1).PDT.CRIT]',...
        [GROUP(1).EXP(2).PDT.CRIT,GROUP(2).EXP(2).PDT.CRIT,GROUP(3).EXP(2).PDT.CRIT]',...
        [GROUP(1).EXP(3).PDT.CRIT,GROUP(2).EXP(3).PDT.CRIT,GROUP(3).EXP(3).PDT.CRIT]',...
        'VariableNames',{'group','FL','FM','FH','DL','DM','DH'});
    
elseif strcmp(tested,'confmd')
    g3 = [THIS(1).g3; THIS(2).g3];
    % Control vs. FMD comparison
    t = table(g3,...
        [GROUP(1).EXP(1).PST.CRIT,GROUP(2).EXP(1).PST.CRIT]',...
        [GROUP(1).EXP(2).PST.CRIT,GROUP(2).EXP(2).PST.CRIT]',...
        [GROUP(1).EXP(3).PST.CRIT,GROUP(2).EXP(3).PST.CRIT]',...
        [GROUP(1).EXP(1).PDT.CRIT,GROUP(2).EXP(1).PDT.CRIT]',...
        [GROUP(1).EXP(2).PDT.CRIT,GROUP(2).EXP(2).PDT.CRIT]',...
        [GROUP(1).EXP(3).PDT.CRIT,GROUP(2).EXP(3).PDT.CRIT]',...
        'VariableNames',{'group','FL','FM','FH','DL','DM','DH'});
end

within = table([1 1 1 2 2 2]',[1 2 3 1 2 3]','VariableNames',{'attention','expectations'});
rm = fitrm(t,'FL-DH ~ group','WithinDesign',within);

ranova(rm,'WithinModel','attention*expectations')

% fectsiz = reshape(testing,[],1);
[p,tbl,stat] = anovan(testing,{g1 g2 re3},'model','full',...
    'varnames',{'Expectations','Attention','Treatment'})
c = multcompare(stat,'dimension',[2 3])

%% Conversion to TABLE for JASP analysis

JASP.criterion = [[GROUP(1).EXP(3).PST.CRIT'; GROUP(2).EXP(3).PST.CRIT';GROUP(3).EXP(3).PST.CRIT'],...
[GROUP(1).EXP(2).PST.CRIT'; GROUP(2).EXP(2).PST.CRIT';GROUP(3).EXP(2).PST.CRIT'],...
[GROUP(1).EXP(1).PST.CRIT'; GROUP(2).EXP(1).PST.CRIT';GROUP(3).EXP(1).PST.CRIT'],...
[GROUP(1).EXP(3).PDT.CRIT'; GROUP(2).EXP(3).PDT.CRIT';GROUP(3).EXP(3).PDT.CRIT'],...
[GROUP(1).EXP(2).PDT.CRIT'; GROUP(2).EXP(2).PDT.CRIT';GROUP(3).EXP(2).PDT.CRIT'],...
[GROUP(1).EXP(1).PDT.CRIT'; GROUP(2).EXP(1).PDT.CRIT';GROUP(3).EXP(1).PDT.CRIT']];

JASP.dPrime = [[GROUP(1).EXP(3).PST.dPRIME'; GROUP(2).EXP(3).PST.dPRIME';GROUP(3).EXP(3).PST.dPRIME'],...
[GROUP(1).EXP(2).PST.dPRIME'; GROUP(2).EXP(2).PST.dPRIME';GROUP(3).EXP(2).PST.dPRIME'],...
[GROUP(1).EXP(1).PST.dPRIME'; GROUP(2).EXP(1).PST.dPRIME';GROUP(3).EXP(1).PST.dPRIME'],...
[GROUP(1).EXP(3).PDT.dPRIME'; GROUP(2).EXP(3).PDT.dPRIME';GROUP(3).EXP(3).PDT.dPRIME'],...
[GROUP(1).EXP(2).PDT.dPRIME'; GROUP(2).EXP(2).PDT.dPRIME';GROUP(3).EXP(2).PDT.dPRIME'],...
[GROUP(1).EXP(1).PDT.dPRIME'; GROUP(2).EXP(1).PDT.dPRIME';GROUP(3).EXP(1).PDT.dPRIME']];

JASP.MetPresent = [[GROUP(1).EXP(3).PST.METPrs'; GROUP(2).EXP(3).PST.METPrs';GROUP(3).EXP(3).PST.METPrs'],...
[GROUP(1).EXP(2).PST.METPrs'; GROUP(2).EXP(2).PST.METPrs';GROUP(3).EXP(2).PST.METPrs'],...
[GROUP(1).EXP(1).PST.METPrs'; GROUP(2).EXP(1).PST.METPrs';GROUP(3).EXP(1).PST.METPrs'],...
[GROUP(1).EXP(3).PDT.METPrs'; GROUP(2).EXP(3).PDT.METPrs';GROUP(3).EXP(3).PDT.METPrs'],...
[GROUP(1).EXP(2).PDT.METPrs'; GROUP(2).EXP(2).PDT.METPrs';GROUP(3).EXP(2).PDT.METPrs'],...
[GROUP(1).EXP(1).PDT.METPrs'; GROUP(2).EXP(1).PDT.METPrs';GROUP(3).EXP(1).PDT.METPrs']];

JASP.MetAbsent = [[GROUP(1).EXP(3).PST.METAbs'; GROUP(2).EXP(3).PST.METAbs';GROUP(3).EXP(3).PST.METAbs'],...
[GROUP(1).EXP(2).PST.METAbs'; GROUP(2).EXP(2).PST.METAbs';GROUP(3).EXP(2).PST.METAbs'],...
[GROUP(1).EXP(1).PST.METAbs'; GROUP(2).EXP(1).PST.METAbs';GROUP(3).EXP(1).PST.METAbs'],...
[GROUP(1).EXP(3).PDT.METAbs'; GROUP(2).EXP(3).PDT.METAbs';GROUP(3).EXP(3).PDT.METAbs'],...
[GROUP(1).EXP(2).PDT.METAbs'; GROUP(2).EXP(2).PDT.METAbs';GROUP(3).EXP(2).PDT.METAbs'],...
[GROUP(1).EXP(1).PDT.METAbs'; GROUP(2).EXP(1).PDT.METAbs';GROUP(3).EXP(1).PDT.METAbs']];

JASP.MetCritPresent = [[GROUP(1).EXP(3).PST.METCritPrs'; GROUP(2).EXP(3).PST.METCritPrs';GROUP(3).EXP(3).PST.METCritPrs'],...
[GROUP(1).EXP(2).PST.METCritPrs'; GROUP(2).EXP(2).PST.METCritPrs';GROUP(3).EXP(2).PST.METCritPrs'],...
[GROUP(1).EXP(1).PST.METCritPrs'; GROUP(2).EXP(1).PST.METCritPrs';GROUP(3).EXP(1).PST.METCritPrs'],...
[GROUP(1).EXP(3).PDT.METCritPrs'; GROUP(2).EXP(3).PDT.METCritPrs';GROUP(3).EXP(3).PDT.METCritPrs'],...
[GROUP(1).EXP(2).PDT.METCritPrs'; GROUP(2).EXP(2).PDT.METCritPrs';GROUP(3).EXP(2).PDT.METCritPrs'],...
[GROUP(1).EXP(1).PDT.METCritPrs'; GROUP(2).EXP(1).PDT.METCritPrs';GROUP(3).EXP(1).PDT.METCritPrs']];

JASP.MetCritAbsent = [[GROUP(1).EXP(3).PST.METCritAbs'; GROUP(2).EXP(3).PST.METCritAbs';GROUP(3).EXP(3).PST.METCritAbs'],...
[GROUP(1).EXP(2).PST.METCritAbs'; GROUP(2).EXP(2).PST.METCritAbs';GROUP(3).EXP(2).PST.METCritAbs'],...
[GROUP(1).EXP(1).PST.METCritAbs'; GROUP(2).EXP(1).PST.METCritAbs';GROUP(3).EXP(1).PST.METCritAbs'],...
[GROUP(1).EXP(3).PDT.METCritAbs'; GROUP(2).EXP(3).PDT.METCritAbs';GROUP(3).EXP(3).PDT.METCritAbs'],...
[GROUP(1).EXP(2).PDT.METCritAbs'; GROUP(2).EXP(2).PDT.METCritAbs';GROUP(3).EXP(2).PDT.METCritAbs'],...
[GROUP(1).EXP(1).PDT.METCritAbs'; GROUP(2).EXP(1).PDT.METCritAbs';GROUP(3).EXP(1).PDT.METCritAbs']];

JASP.confidc = [[GROUP(1).EXP(3).PST.CONFID_correct'; GROUP(2).EXP(3).PST.CONFID_correct';GROUP(3).EXP(3).PST.CONFID_correct'],...
[GROUP(1).EXP(2).PST.CONFID_correct'; GROUP(2).EXP(2).PST.CONFID_correct';GROUP(3).EXP(2).PST.CONFID_correct'],...
[GROUP(1).EXP(1).PST.CONFID_correct'; GROUP(2).EXP(1).PST.CONFID_correct';GROUP(3).EXP(1).PST.CONFID_correct'],...
[GROUP(1).EXP(3).PDT.CONFID_correct'; GROUP(2).EXP(3).PDT.CONFID_correct';GROUP(3).EXP(3).PDT.CONFID_correct'],...
[GROUP(1).EXP(2).PDT.CONFID_correct'; GROUP(2).EXP(2).PDT.CONFID_correct';GROUP(3).EXP(2).PDT.CONFID_correct'],...
[GROUP(1).EXP(1).PDT.CONFID_correct'; GROUP(2).EXP(1).PDT.CONFID_correct';GROUP(3).EXP(1).PDT.CONFID_correct']];

JASP.confidi = [[GROUP(1).EXP(3).PST.CONFID_incorrect'; GROUP(2).EXP(3).PST.CONFID_incorrect';GROUP(3).EXP(3).PST.CONFID_incorrect'],...
[GROUP(1).EXP(2).PST.CONFID_incorrect'; GROUP(2).EXP(2).PST.CONFID_incorrect';GROUP(3).EXP(2).PST.CONFID_incorrect'],...
[GROUP(1).EXP(1).PST.CONFID_incorrect'; GROUP(2).EXP(1).PST.CONFID_incorrect';GROUP(3).EXP(1).PST.CONFID_incorrect'],...
[GROUP(1).EXP(3).PDT.CONFID_incorrect'; GROUP(2).EXP(3).PDT.CONFID_incorrect';GROUP(3).EXP(3).PDT.CONFID_incorrect'],...
[GROUP(1).EXP(2).PDT.CONFID_incorrect'; GROUP(2).EXP(2).PDT.CONFID_incorrect';GROUP(3).EXP(2).PDT.CONFID_incorrect'],...
[GROUP(1).EXP(1).PDT.CONFID_incorrect'; GROUP(2).EXP(1).PDT.CONFID_incorrect';GROUP(3).EXP(1).PDT.CONFID_incorrect']];

% Present "it's present"
JASP.conf_pp = [[GROUP(1).EXP(3).PST.CONFsd_pp'; GROUP(2).EXP(3).PST.CONFsd_pp';GROUP(3).EXP(3).PST.CONFsd_pp'],...
[GROUP(1).EXP(2).PST.CONFsd_pp'; GROUP(2).EXP(2).PST.CONFsd_pp';GROUP(3).EXP(2).PST.CONFsd_pp'],...
[GROUP(1).EXP(1).PST.CONFsd_pp'; GROUP(2).EXP(1).PST.CONFsd_pp';GROUP(3).EXP(1).PST.CONFsd_pp'],...
[GROUP(1).EXP(3).PDT.CONFsd_pp'; GROUP(2).EXP(3).PDT.CONFsd_pp';GROUP(3).EXP(3).PDT.CONFsd_pp'],...
[GROUP(1).EXP(2).PDT.CONFsd_pp'; GROUP(2).EXP(2).PDT.CONFsd_pp';GROUP(3).EXP(2).PDT.CONFsd_pp'],...
[GROUP(1).EXP(1).PDT.CONFsd_pp'; GROUP(2).EXP(1).PDT.CONFsd_pp';GROUP(3).EXP(1).PDT.CONFsd_pp']];

% Present "it's absent"
JASP.conf_pa = [[GROUP(1).EXP(3).PST.CONFsd_pa'; GROUP(2).EXP(3).PST.CONFsd_pa';GROUP(3).EXP(3).PST.CONFsd_pa'],...
[GROUP(1).EXP(2).PST.CONFsd_pa'; GROUP(2).EXP(2).PST.CONFsd_pa';GROUP(3).EXP(2).PST.CONFsd_pa'],...
[GROUP(1).EXP(1).PST.CONFsd_pa'; GROUP(2).EXP(1).PST.CONFsd_pa';GROUP(3).EXP(1).PST.CONFsd_pa'],...
[GROUP(1).EXP(3).PDT.CONFsd_pa'; GROUP(2).EXP(3).PDT.CONFsd_pa';GROUP(3).EXP(3).PDT.CONFsd_pa'],...
[GROUP(1).EXP(2).PDT.CONFsd_pa'; GROUP(2).EXP(2).PDT.CONFsd_pa';GROUP(3).EXP(2).PDT.CONFsd_pa'],...
[GROUP(1).EXP(1).PDT.CONFsd_pa'; GROUP(2).EXP(1).PDT.CONFsd_pa';GROUP(3).EXP(1).PDT.CONFsd_pa']];

% Absent "it's absent"
JASP.conf_aa = [[GROUP(1).EXP(3).PST.CONFsd_aa'; GROUP(2).EXP(3).PST.CONFsd_aa';GROUP(3).EXP(3).PST.CONFsd_aa'],...
[GROUP(1).EXP(2).PST.CONFsd_aa'; GROUP(2).EXP(2).PST.CONFsd_aa';GROUP(3).EXP(2).PST.CONFsd_aa'],...
[GROUP(1).EXP(1).PST.CONFsd_aa'; GROUP(2).EXP(1).PST.CONFsd_aa';GROUP(3).EXP(1).PST.CONFsd_aa'],...
[GROUP(1).EXP(3).PDT.CONFsd_aa'; GROUP(2).EXP(3).PDT.CONFsd_aa';GROUP(3).EXP(3).PDT.CONFsd_aa'],...
[GROUP(1).EXP(2).PDT.CONFsd_aa'; GROUP(2).EXP(2).PDT.CONFsd_aa';GROUP(3).EXP(2).PDT.CONFsd_aa'],...
[GROUP(1).EXP(1).PDT.CONFsd_aa'; GROUP(2).EXP(1).PDT.CONFsd_aa';GROUP(3).EXP(1).PDT.CONFsd_aa']];

% Absent "it's present"
JASP.conf_ap = [[GROUP(1).EXP(3).PST.CONFsd_ap'; GROUP(2).EXP(3).PST.CONFsd_ap';GROUP(3).EXP(3).PST.CONFsd_ap'],...
[GROUP(1).EXP(2).PST.CONFsd_ap'; GROUP(2).EXP(2).PST.CONFsd_ap';GROUP(3).EXP(2).PST.CONFsd_ap'],...
[GROUP(1).EXP(1).PST.CONFsd_ap'; GROUP(2).EXP(1).PST.CONFsd_ap';GROUP(3).EXP(1).PST.CONFsd_ap'],...
[GROUP(1).EXP(3).PDT.CONFsd_ap'; GROUP(2).EXP(3).PDT.CONFsd_ap';GROUP(3).EXP(3).PDT.CONFsd_ap'],...
[GROUP(1).EXP(2).PDT.CONFsd_ap'; GROUP(2).EXP(2).PDT.CONFsd_ap';GROUP(3).EXP(2).PDT.CONFsd_ap'],...
[GROUP(1).EXP(1).PDT.CONFsd_ap'; GROUP(2).EXP(1).PDT.CONFsd_ap';GROUP(3).EXP(1).PDT.CONFsd_ap']];

everything = [JASP.criterion,JASP.dPrime,JASP.MetPresent,JASP.MetAbsent,...
    JASP.MetCritPresent,JASP.MetCritAbsent,JASP.confidc,JASP.confidi];

confsd = [JASP.conf_pp,JASP.conf_pa,JASP.conf_aa,JASP.conf_ap];

% csvwrite('../../data/JASP/criterion.csv',JASP.criterion);
% csvwrite('../../data/JASP/dPrime.csv',JASP.dPrime);
% csvwrite('../../data/JASP/MetPres.csv',JASP.MetPresent);
% csvwrite('../../data/JASP/MetAbs.csv',JASP.MetAbsent);
% csvwrite('../../data/JASP/MetCritPres.csv',JASP.MetCritPresent);
% csvwrite('../../data/JASP/MetCritAbs.csv',JASP.MetCritAbsent);
% csvwrite('../../data/JASP/confid.csv',[JASP.confidc, JASP.confidi]);
csvwrite('../../data/JASP/confsd.csv',confsd);
csvwrite('../../data/JASP/everything.csv',everything);

%% 3-WAY CRITERION ANALYSIS

% Need to subset Present judgment criterion vs. Absent

% for gr = 1:3
%     EXP = GROUP(gr).EXP; subnum = length(EXP(1).PDT.CRIT);
%     all_subjs = [EXP(1).PDT.CRIT,EXP(1).PST.CRIT;...
%         EXP(2).PDT.CRIT,EXP(2).PST.CRIT;...
%         EXP(3).PDT.CRIT,EXP(3).PST.CRIT;...
%         EXP(1).PDT.CRIT,EXP(1).PST.CRIT;...
%         EXP(2).PDT.CRIT,EXP(2).PST.CRIT;...
%         EXP(3).PDT.CRIT,EXP(3).PST.CRIT]';
%     
%     g1 = repmat([repmat({'25%'},subnum*2,1);repmat({'50%'},subnum*2,1);repmat({'75%'},subnum*2,1)],2,1); % EXPECTATIONS
%     
%     g2 = repmat([repmat({'D'},subnum,1);repmat({'F'},subnum,1)],6,1); % ATTENTION
%     
%     g3 = [repmat({'A'},subnum*6,1);repmat({'P'},subnum*6,1)]; % REPORT
%     
%     fectsiz = reshape(all_subjs,[],1);
%     [p,tbl,stat] = anovan(fectsiz,{g1 g2 g3},'model','full',...
%         'varnames',{'Expectations','Attention','Report'})
% end

%% REPORT PRESENT METACOG Plot 'em

figure

subplot(3,1,1)
gr = 1;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.METPrs;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.METPrs;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.METPrs;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.METPrs;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.METPrs;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.METPrs;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

tempWSE=within_subject_error(...
    [GROUP(gr).EXP(3).PST.METPrs;GROUP(gr).EXP(2).PST.METPrs;GROUP(gr).EXP(1).PST.METPrs;...
    GROUP(gr).EXP(3).PDT.METPrs;GROUP(gr).EXP(2).PDT.METPrs;GROUP(gr).EXP(1).PDT.METPrs;...
    ]',[2,6]);
PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
% sigstar({[1,3]})
% xlabel('Expectations')
box off
set(gca,'Tickdir','out')
words = sprintf('CONTROL (n=%s)',mat2str(subnum)); ylabel(words);
title('Report Present Type 2 Dprime')

legend({'full' 'diverted'})
legend('boxoff')

subplot(3,1,2)
gr = 2;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.METPrs;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.METPrs;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.METPrs;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.METPrs;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.METPrs;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.METPrs;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

tempWSE=within_subject_error(...
    [GROUP(gr).EXP(3).PST.METPrs;GROUP(gr).EXP(2).PST.METPrs;GROUP(gr).EXP(1).PST.METPrs;...
    GROUP(gr).EXP(3).PDT.METPrs;GROUP(gr).EXP(2).PDT.METPrs;GROUP(gr).EXP(1).PDT.METPrs;...
    ]',[2,6]);
PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
% sigstar({[1,3]})
% xlabel('Expectations')
box off
set(gca,'Tickdir','out')
words = sprintf('FMD (n=%s)',mat2str(subnum)); ylabel(words);

subplot(3,1,3)
gr = 3;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.METPrs;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.METPrs;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.METPrs;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.METPrs;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.METPrs;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.METPrs;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

tempWSE=within_subject_error(...
    [GROUP(gr).EXP(3).PST.METPrs;GROUP(gr).EXP(2).PST.METPrs;GROUP(gr).EXP(1).PST.METPrs;...
    GROUP(gr).EXP(3).PDT.METPrs;GROUP(gr).EXP(2).PDT.METPrs;GROUP(gr).EXP(1).PDT.METPrs;...
    ]',[2,6]);
PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% sigstar({[1,3]})
% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
xlabel('Expectations of Target Presence')
box off
set(gca,'Tickdir','out')
words = sprintf('ORGANIC (n=%s)',mat2str(subnum)); ylabel(words);

%% REPORT PRESENT METACOG STATISTICS

EXP = GROUP(1).EXP; subnum = length(EXP(1).PDT.METPrs);
all_subjs = [EXP(1).PDT.METPrs,EXP(1).PST.METPrs;...
    EXP(2).PDT.METPrs,EXP(2).PST.METPrs;...
    EXP(3).PDT.METPrs,EXP(3).PST.METPrs]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Signcificance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

EXP = GROUP(2).EXP; subnum = length(EXP(1).PDT.METPrs);
all_subjs = [EXP(1).PDT.METPrs,EXP(1).PST.METPrs;...
    EXP(2).PDT.METPrs,EXP(2).PST.METPrs;...
    EXP(3).PDT.METPrs,EXP(3).PST.METPrs]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

EXP = GROUP(3).EXP; subnum = length(EXP(1).PDT.METPrs);
all_subjs = [EXP(1).PDT.METPrs,EXP(1).PST.METPrs;...
    EXP(2).PDT.METPrs,EXP(2).PST.METPrs;...
    EXP(3).PDT.METPrs,EXP(3).PST.METPrs]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

%% REPORT ABSENT METACOG Plot 'em

figure

subplot(3,1,1)
gr = 1;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.METAbs;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.METAbs;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.METAbs;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.METAbs;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.METAbs;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.METAbs;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

tempWSE=within_subject_error(...
    [GROUP(gr).EXP(3).PST.METAbs;GROUP(gr).EXP(2).PST.METAbs;GROUP(gr).EXP(1).PST.METAbs;...
    GROUP(gr).EXP(3).PDT.METAbs;GROUP(gr).EXP(2).PDT.METAbs;GROUP(gr).EXP(1).PDT.METAbs;...
    ]',[2,6]);
PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
% sigstar({[1,3]})
% xlabel('Expectations')
box off
set(gca,'Tickdir','out')
words = sprintf('CONTROL (n=%s)',mat2str(subnum)); ylabel(words);
title('Report Absent Type 2 Dprime')

legend({'full' 'diverted'})
legend('boxoff')

subplot(3,1,2)
gr = 2;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.METAbs;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.METAbs;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.METAbs;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.METAbs;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.METAbs;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.METAbs;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

tempWSE=within_subject_error(...
    [GROUP(gr).EXP(3).PST.METAbs;GROUP(gr).EXP(2).PST.METAbs;GROUP(gr).EXP(1).PST.METAbs;...
    GROUP(gr).EXP(3).PDT.METAbs;GROUP(gr).EXP(2).PDT.METAbs;GROUP(gr).EXP(1).PDT.METAbs;...
    ]',[2,6]);
PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
% sigstar({[1,3]})
% xlabel('Expectations')
box off
set(gca,'Tickdir','out')
words = sprintf('FMD (n=%s)',mat2str(subnum)); ylabel(words);

subplot(3,1,3)
gr = 3;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.METAbs;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.METAbs;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.METAbs;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.METAbs;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.METAbs;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.METAbs;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

tempWSE=within_subject_error(...
    [GROUP(gr).EXP(3).PST.METAbs;GROUP(gr).EXP(2).PST.METAbs;GROUP(gr).EXP(1).PST.METAbs;...
    GROUP(gr).EXP(3).PDT.METAbs;GROUP(gr).EXP(2).PDT.METAbs;GROUP(gr).EXP(1).PDT.METAbs;...
    ]',[2,6]);
PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% sigstar({[1,3]})
% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
xlabel('Expectations of Target Presence')
box off
set(gca,'Tickdir','out')
words = sprintf('ORGANIC (n=%s)',mat2str(subnum)); ylabel(words);

%% REPORT ABSENT METACOG STATISTICS

EXP = GROUP(1).EXP; subnum = length(EXP(1).PDT.METAbs);
all_subjs = [EXP(1).PDT.METAbs,EXP(1).PST.METAbs;...
    EXP(2).PDT.METAbs,EXP(2).PST.METAbs;...
    EXP(3).PDT.METAbs,EXP(3).PST.METAbs]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

EXP = GROUP(2).EXP; subnum = length(EXP(1).PDT.METAbs);
all_subjs = [EXP(1).PDT.METAbs,EXP(1).PST.METAbs;...
    EXP(2).PDT.METAbs,EXP(2).PST.METAbs;...
    EXP(3).PDT.METAbs,EXP(3).PST.METAbs]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

EXP = GROUP(3).EXP; subnum = length(EXP(1).PDT.METAbs);
all_subjs = [EXP(1).PDT.METAbs,EXP(1).PST.METAbs;...
    EXP(2).PDT.METAbs,EXP(2).PST.METAbs;...
    EXP(3).PDT.METAbs,EXP(3).PST.METAbs]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

%% 3-WAY METACOG STATISTICS

for gr = 1:3
    EXP = GROUP(gr).EXP; subnum = length(EXP(1).PDT.METAbs);
    all_subjs = [EXP(1).PDT.METAbs,EXP(1).PST.METAbs;...
        EXP(2).PDT.METAbs,EXP(2).PST.METAbs;...
        EXP(3).PDT.METAbs,EXP(3).PST.METAbs;...
        EXP(1).PDT.METPrs,EXP(1).PST.METPrs;...
        EXP(2).PDT.METPrs,EXP(2).PST.METPrs;...
        EXP(3).PDT.METPrs,EXP(3).PST.METPrs]';
    
    g1 = repmat([repmat({'25%'},subnum*2,1);repmat({'50%'},subnum*2,1);repmat({'75%'},subnum*2,1)],2,1); % EXPECTATIONS
    
    g2 = repmat([repmat({'D'},subnum,1);repmat({'F'},subnum,1)],6,1); % ATTENTION
    
    g3 = [repmat({'A'},subnum*6,1);repmat({'P'},subnum*6,1)]; % REPORT
    
    fectsiz = reshape(all_subjs,[],1);
    [p,tbl,stat] = anovan(fectsiz,{g1 g2 g3},'model','full',...
        'varnames',{'Expectations','Attention','Report'})
end

% Only effect of report is seen for Type 2 D' (lower metacog Report Absent)
% Ratings are very uniform for FMD (attention & exp effects p>.25)
% Trend in Controls for effect of expectations

%% REPORT PRESENT METACOG (Type 2 AUC: Fleming & Lau, 2014)

figure

subplot(3,1,1)
gr = 1;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.METAUCPrs;
PST_mean(1) = nanmean(measure); PST_SEM(1) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.METAUCPrs;
PST_mean(2) = nanmean(measure); PST_SEM(2) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.METAUCPrs;
PST_mean(3) = nanmean(measure); PST_SEM(3) = nanstd(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.METAUCPrs;
PDT_mean(1) = nanmean(measure); PDT_SEM(1) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.METAUCPrs;
PDT_mean(2) = nanmean(measure); PDT_SEM(2) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.METAUCPrs;
PDT_mean(3) = nanmean(measure); PDT_SEM(3) = nanstd(measure)/sqrt(subnum);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
% sigstar({[1,3]})
plot([0 4],[.5,.5],'k:');
xlim([0.5 3.5])
box off
set(gca,'Tickdir','out')
words = sprintf('CONTROL (n=%s)',mat2str(subnum)); ylabel(words);
title('Report Present Type 2 AUC')

legend({'full' 'diverted'})
legend('boxoff')

subplot(3,1,2)
gr = 2;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.METAUCPrs;
PST_mean(1) = nanmean(measure); PST_SEM(1) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.METAUCPrs;
PST_mean(2) = nanmean(measure); PST_SEM(2) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.METAUCPrs;
PST_mean(3) = nanmean(measure); PST_SEM(3) = nanstd(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.METAUCPrs;
PDT_mean(1) = nanmean(measure); PDT_SEM(1) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.METAUCPrs;
PDT_mean(2) = nanmean(measure); PDT_SEM(2) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.METAUCPrs;
PDT_mean(3) = nanmean(measure); PDT_SEM(3) = nanstd(measure)/sqrt(subnum);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
% sigstar({[1,3]})
% xlabel('Expectations')
plot([0 4],[.5,.5],'k:');
xlim([0.5 3.5])
box off
set(gca,'Tickdir','out')
words = sprintf('FMD (n=%s)',mat2str(subnum)); ylabel(words);

subplot(3,1,3)
gr = 3;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.METAUCPrs;
PST_mean(1) = nanmean(measure); PST_SEM(1) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.METAUCPrs;
PST_mean(2) = nanmean(measure); PST_SEM(2) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.METAUCPrs;
PST_mean(3) = nanmean(measure); PST_SEM(3) = nanstd(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.METAUCPrs;
PDT_mean(1) = nanmean(measure); PDT_SEM(1) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.METAUCPrs;
PDT_mean(2) = nanmean(measure); PDT_SEM(2) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.METAUCPrs;
PDT_mean(3) = nanmean(measure); PDT_SEM(3) = nanstd(measure)/sqrt(subnum);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% sigstar({[1,3]})
% ylim([0.5,3])

xticklabels({'25%' '50%' '75%'})
xlabel('Expectations of Target Presence')
plot([0 4],[.5,.5],'k:');
xlim([0.5 3.5])
box off
set(gca,'Tickdir','out')
words = sprintf('ORGANIC (n=%s)',mat2str(subnum)); ylabel(words);

%% REPORT PRESENT METACOG AUC STATISTICS

EXP = GROUP(1).EXP; subnum = length(EXP(1).PDT.METAUCPrs);
all_subjs = [EXP(1).PDT.METAUCPrs,EXP(1).PST.METAUCPrs;...
    EXP(2).PDT.METAUCPrs,EXP(2).PST.METAUCPrs;...
    EXP(3).PDT.METAUCPrs,EXP(3).PST.METAUCPrs]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Signcificance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

EXP = GROUP(2).EXP; subnum = length(EXP(1).PDT.METAUCPrs);
all_subjs = [EXP(1).PDT.METAUCPrs,EXP(1).PST.METAUCPrs;...
    EXP(2).PDT.METAUCPrs,EXP(2).PST.METAUCPrs;...
    EXP(3).PDT.METAUCPrs,EXP(3).PST.METAUCPrs]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

EXP = GROUP(3).EXP; subnum = length(EXP(1).PDT.METAUCPrs);
all_subjs = [EXP(1).PDT.METAUCPrs,EXP(1).PST.METAUCPrs;...
    EXP(2).PDT.METAUCPrs,EXP(2).PST.METAUCPrs;...
    EXP(3).PDT.METAUCPrs,EXP(3).PST.METAUCPrs]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

%% REPORT ABSENT METACOG Type 2 AUC

figure

subplot(3,1,1)
gr = 1;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.METAUCAbs;
PST_mean(1) = nanmean(measure); PST_SEM(1) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.METAUCAbs;
PST_mean(2) = nanmean(measure); PST_SEM(2) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.METAUCAbs;
PST_mean(3) = nanmean(measure); PST_SEM(3) = nanstd(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.METAUCAbs;
PDT_mean(1) = nanmean(measure); PDT_SEM(1) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.METAUCAbs;
PDT_mean(2) = nanmean(measure); PDT_SEM(2) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.METAUCAbs;
PDT_mean(3) = nanmean(measure); PDT_SEM(3) = nanstd(measure)/sqrt(subnum);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
% sigstar({[1,3]})
plot([0 4],[.5,.5],'k:');
xlim([0.5 3.5])
box off
set(gca,'Tickdir','out')
words = sprintf('CONTROL (n=%s)',mat2str(subnum)); ylabel(words);
title('Report Absent Type 2 AUC')

legend({'full attention' 'diverted attention'})
legend('boxoff')

subplot(3,1,2)
gr = 2;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.METAUCAbs;
PST_mean(1) = nanmean(measure); PST_SEM(1) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.METAUCAbs;
PST_mean(2) = nanmean(measure); PST_SEM(2) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.METAUCAbs;
PST_mean(3) = nanmean(measure); PST_SEM(3) = nanstd(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.METAUCAbs;
PDT_mean(1) = nanmean(measure); PDT_SEM(1) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.METAUCAbs;
PDT_mean(2) = nanmean(measure); PDT_SEM(2) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.METAUCAbs;
PDT_mean(3) = nanmean(measure); PDT_SEM(3) = nanstd(measure)/sqrt(subnum);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
% sigstar({[1,3]})
% xlabel('Expectations')
plot([0 4],[.5,.5],'k:');
xlim([0.5 3.5])
box off
set(gca,'Tickdir','out')
words = sprintf('FMD (n=%s)',mat2str(subnum)); ylabel(words);

subplot(3,1,3)
gr = 3;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.METAUCAbs;
PST_mean(1) = nanmean(measure); PST_SEM(1) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.METAUCAbs;
PST_mean(2) = nanmean(measure); PST_SEM(2) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.METAUCAbs;
PST_mean(3) = nanmean(measure); PST_SEM(3) = nanstd(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.METAUCAbs;
PDT_mean(1) = nanmean(measure); PDT_SEM(1) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.METAUCAbs;
PDT_mean(2) = nanmean(measure); PDT_SEM(2) = nanstd(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.METAUCAbs;
PDT_mean(3) = nanmean(measure); PDT_SEM(3) = nanstd(measure)/sqrt(subnum);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% sigstar({[1,3]})
% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
plot([0 4],[.5,.5],'k:');
xlim([0.5 3.5])
xlabel('Expectations of Target Presence')
box off
set(gca,'Tickdir','out')
words = sprintf('ORGANIC (n=%s)',mat2str(subnum)); ylabel(words);

%% REPORT ABSENT METACOG AUC STATISTICS

EXP = GROUP(1).EXP; subnum = length(EXP(1).PDT.METAUCAbs);
all_subjs = [EXP(1).PDT.METAUCAbs,EXP(1).PST.METAUCAbs;...
    EXP(2).PDT.METAUCAbs,EXP(2).PST.METAUCAbs;...
    EXP(3).PDT.METAUCAbs,EXP(3).PST.METAUCAbs]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Signcificance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

EXP = GROUP(2).EXP; subnum = length(EXP(1).PDT.METAUCAbs);
all_subjs = [EXP(1).PDT.METAUCAbs,EXP(1).PST.METAUCAbs;...
    EXP(2).PDT.METAUCAbs,EXP(2).PST.METAUCAbs;...
    EXP(3).PDT.METAUCAbs,EXP(3).PST.METAUCAbs]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking

EXP = GROUP(3).EXP; subnum = length(EXP(1).PDT.METAUCAbs);
all_subjs = [EXP(1).PDT.METAUCAbs,EXP(1).PST.METAUCAbs;...
    EXP(2).PDT.METAUCAbs,EXP(2).PST.METAUCAbs;...
    EXP(3).PDT.METAUCAbs,EXP(3).PST.METAUCAbs]';

expec = reshape(repmat(1:3,subnum*2,1),[],1); % Expectation levels
attens = reshape(repmat(1:2,subnum,3),[],1); % Attention levels
fectsiz = reshape(all_subjs,[],1);

mes2way(fectsiz,[expec attens],'partialeta2'); % Significance testing
% [p,tbl,stat] = anova2(all_subjs,subnum); % Double-checking


%% 3-WAY METACOG AUCSTATISTICS

for gr = 1:3
    EXP = GROUP(gr).EXP; subnum = length(EXP(1).PDT.METAUCAbs);
    all_subjs = [EXP(1).PDT.METAUCAbs,EXP(1).PST.METAUCAbs;...
        EXP(2).PDT.METAUCAbs,EXP(2).PST.METAUCAbs;...
        EXP(3).PDT.METAUCAbs,EXP(3).PST.METAUCAbs;...
        EXP(1).PDT.METAUCPrs,EXP(1).PST.METAUCPrs;...
        EXP(2).PDT.METAUCPrs,EXP(2).PST.METAUCPrs;...
        EXP(3).PDT.METAUCPrs,EXP(3).PST.METAUCPrs]';
    
    g1 = repmat([repmat({'25%'},subnum*2,1);repmat({'50%'},subnum*2,1);repmat({'75%'},subnum*2,1)],2,1); % EXPECTATIONS
    
    g2 = repmat([repmat({'D'},subnum,1);repmat({'F'},subnum,1)],6,1); % ATTENTION
    
    g3 = [repmat({'A'},subnum*6,1);repmat({'P'},subnum*6,1)]; % REPORT
    
    fectsiz = reshape(all_subjs,[],1);
    [p,tbl,stat] = anovan(fectsiz,{g1 g2 g3},'model','full',...
        'varnames',{'Expectations','Attention','Report'})
end

% Only effect of report is seen for Type 2 D' (lower metacog Report Absent)
% Ratings are very uniform for FMD (attention & exp effects p>.25)
% Trend in Controls for effect of expectations

%% REPORT PRESENT METACOG CRITERION Plot 'em

figure

subplot(3,1,1)
gr = 1;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.METCritPrs;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.METCritPrs;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.METCritPrs;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.METCritPrs;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.METCritPrs;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.METCritPrs;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

tempWSE=within_subject_error(...
    [GROUP(gr).EXP(3).PST.METCritPrs;GROUP(gr).EXP(2).PST.METCritPrs;GROUP(gr).EXP(1).PST.METCritPrs;...
    GROUP(gr).EXP(3).PDT.METCritPrs;GROUP(gr).EXP(2).PDT.METCritPrs;GROUP(gr).EXP(1).PDT.METCritPrs;...
    ]',[2,3]);
PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
% sigstar({[1,3]},.0012)
% xlabel('Expectations')
box off
set(gca,'Tickdir','out')
words = sprintf('CONTROL (n=%s)',mat2str(subnum)); ylabel(words);
title('Report Present Type 2 Criterion')

legend({'full' 'diverted'})
legend('boxoff')

subplot(3,1,2)
gr = 2;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.METCritPrs;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.METCritPrs;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.METCritPrs;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.METCritPrs;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.METCritPrs;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.METCritPrs;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

tempWSE=within_subject_error(...
    [GROUP(gr).EXP(3).PST.METCritPrs;GROUP(gr).EXP(2).PST.METCritPrs;GROUP(gr).EXP(1).PST.METCritPrs;...
    GROUP(gr).EXP(3).PDT.METCritPrs;GROUP(gr).EXP(2).PDT.METCritPrs;GROUP(gr).EXP(1).PDT.METCritPrs;...
    ]',[2,3]);
PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
% sigstar({[1,3]},.0239)
% xlabel('Expectations')
box off
set(gca,'Tickdir','out')
words = sprintf('FMD (n=%s)',mat2str(subnum)); ylabel(words);

subplot(3,1,3)
gr = 3;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.METCritPrs;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.METCritPrs;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.METCritPrs;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.METCritPrs;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.METCritPrs;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.METCritPrs;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

tempWSE=within_subject_error(...
    [GROUP(gr).EXP(3).PST.METCritPrs;GROUP(gr).EXP(2).PST.METCritPrs;GROUP(gr).EXP(1).PST.METCritPrs;...
    GROUP(gr).EXP(3).PDT.METCritPrs;GROUP(gr).EXP(2).PDT.METCritPrs;GROUP(gr).EXP(1).PDT.METCritPrs;...
    ]',[2,3]);
PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% sigstar({[1,3]},.0258)
% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
xlabel('Expectations of Target Presence')
box off
set(gca,'Tickdir','out')
words = sprintf('ORGANIC (n=%s)',mat2str(subnum)); ylabel(words);

%% REPORT ABSENT METACOG CRITERION Plot 'em

figure

subplot(3,1,1)
gr = 1;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.METCritAbs;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.METCritAbs;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.METCritAbs;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.METCritAbs;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.METCritAbs;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.METCritAbs;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

tempWSE=within_subject_error(...
    [GROUP(gr).EXP(3).PST.METCritAbs;GROUP(gr).EXP(2).PST.METCritAbs;GROUP(gr).EXP(1).PST.METCritAbs;...
    GROUP(gr).EXP(3).PDT.METCritAbs;GROUP(gr).EXP(2).PDT.METCritAbs;GROUP(gr).EXP(1).PDT.METCritAbs;...
    ]',[2,3]);
PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
% sigstar({[1,3]},.0012)
% xlabel('Expectations')
box off
set(gca,'Tickdir','out')
words = sprintf('CONTROL (n=%s)',mat2str(subnum)); ylabel(words);
title('Report Absent Type 2 Criterion')

legend({'full' 'diverted'})
legend('boxoff')

subplot(3,1,2)
gr = 2;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.METCritAbs;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.METCritAbs;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.METCritAbs;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.METCritAbs;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.METCritAbs;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.METCritAbs;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

tempWSE=within_subject_error(...
    [GROUP(gr).EXP(3).PST.METCritAbs;GROUP(gr).EXP(2).PST.METCritAbs;GROUP(gr).EXP(1).PST.METCritAbs;...
    GROUP(gr).EXP(3).PDT.METCritAbs;GROUP(gr).EXP(2).PDT.METCritAbs;GROUP(gr).EXP(1).PDT.METCritAbs;...
    ]',[2,3]);
PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
% sigstar({[1,3]},.0239)
% xlabel('Expectations')
box off
set(gca,'Tickdir','out')
words = sprintf('FMD (n=%s)',mat2str(subnum)); ylabel(words);

subplot(3,1,3)
gr = 3;
subnum = length(GROUP(gr).COND(1).subjs);

% Measures ['lo' 'med' 'hi']
measure = GROUP(gr).EXP(3).PST.METCritAbs;
PST_mean(1) = mean(measure); PST_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PST.METCritAbs;
PST_mean(2) = mean(measure); PST_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PST.METCritAbs;
PST_mean(3) = mean(measure); PST_SEM(3) = std(measure)/sqrt(subnum);

measure = GROUP(gr).EXP(3).PDT.METCritAbs;
PDT_mean(1) = mean(measure); PDT_SEM(1) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(2).PDT.METCritAbs;
PDT_mean(2) = mean(measure); PDT_SEM(2) = std(measure)/sqrt(subnum);
measure = GROUP(gr).EXP(1).PDT.METCritAbs;
PDT_mean(3) = mean(measure); PDT_SEM(3) = std(measure)/sqrt(subnum);

tempWSE=within_subject_error(...
    [GROUP(gr).EXP(3).PST.METCritAbs;GROUP(gr).EXP(2).PST.METCritAbs;GROUP(gr).EXP(1).PST.METCritAbs;...
    GROUP(gr).EXP(3).PDT.METCritAbs;GROUP(gr).EXP(2).PDT.METCritAbs;GROUP(gr).EXP(1).PDT.METCritAbs;...
    ]',[2,3]);
PST_SEM=tempWSE(1:3); PDT_SEM=tempWSE(4:6);

err = [PST_SEM;PDT_SEM]';
colormap(gray);
h = bar([PST_mean;PDT_mean]');
hold on
errbar(h(1).XData-.14,h(1).YData,err(:,1),'r','LineWidth',2)
errbar(h(2).XData+.14,h(2).YData,err(:,2),'r','LineWidth',2)

% sigstar({[1,3]},.0258)
% ylim([0.5,3])
xticklabels({'25%' '50%' '75%'})
xlabel('Expectations of Target Presence')
box off
set(gca,'Tickdir','out')
words = sprintf('ORGANIC (n=%s)',mat2str(subnum)); ylabel(words);

%% 3-WAY CRITERION ANALYSIS

for gr = 1:3
    EXP = GROUP(gr).EXP; subnum = length(EXP(1).PDT.METCritAbs);
    all_subjs = [EXP(1).PDT.METCritAbs,EXP(1).PST.METCritAbs;...
        EXP(2).PDT.METCritAbs,EXP(2).PST.METCritAbs;...
        EXP(3).PDT.METCritAbs,EXP(3).PST.METCritAbs;...
        EXP(1).PDT.METCritPrs,EXP(1).PST.METCritPrs;...
        EXP(2).PDT.METCritPrs,EXP(2).PST.METCritPrs;...
        EXP(3).PDT.METCritPrs,EXP(3).PST.METCritPrs]';
    
    g1 = repmat([repmat({'25%'},subnum*2,1);repmat({'50%'},subnum*2,1);repmat({'75%'},subnum*2,1)],2,1); % EXPECTATIONS
    
    g2 = repmat([repmat({'D'},subnum,1);repmat({'F'},subnum,1)],6,1); % ATTENTION
    
    g3 = [repmat({'A'},subnum*6,1);repmat({'P'},subnum*6,1)]; % REPORT
    
    fectsiz = reshape(all_subjs,[],1);
    [p,tbl,stat] = anovan(fectsiz,{g1 g2 g3},'model','full',...
        'varnames',{'Expectations','Attention','Report'})
end

% Controls: sig (p=.045) interaction between attention & report on type 2 criterion, else p>.25
% FMD: Report main effect ns (p=.11), else p>.25
% Organic: Attention main effect ns (p=.089), else p>.25

%% OLD STUFF

figure

subplot(1,2,1)

y = [ones(1,subnum)*.2;pCON.full_tr'];
e = [zeros(1,subnum);pCON.full_SEM_tr'];
c = cool(subnum);
er = errorbar(y,e,'o-');
for subject = 1:subnum
    set(er(subject),'color',c(subject,:))
end
box off
set(gca,'TickDir','out','XTickLabel',0:8)
xlim([1,9])
ylim([0,.65])
ylabel('Gabor Contrast (higher reflects greater visibility)');
xlabel('Training Block');
title('Gabor Contrast under Full Attention');
% legend(subjs,'Location','BestOutside');

subplot(1,2,2)
y = pCON.divert_tr';
e = pCON.divert_SEM_tr';
c = cool(subnum);
er = errorbar(y,e,'o-');
for subject = 1:subnum
    set(er(subject),'color',c(subject,:))
end
box off
set(gca,'TickDir','out','XTick',1:9)
xlim([0,8])
ylim([0,.65])
% ylabel('Gabor Contrast (higher reflects greater visibility)');
xlabel('Training Block');
title('Gabor Contrast under Diverted Attention');
% legend(subjs,'Location','BestInside');

%%
% Contrast of peripheral gabor during training under full and diverted attention. Each
% line represents a different participant with errorbars reflecting
% standard error within each block.
%
% These are the first 96 trials subjects respond to in full/diverted attention
%
% In the leftmost plot, their task is to ignore the central letter display and
% attempt to discriminate whether the peripheral gabor or "grill" appeared.
% The gabor appears on 50% of trials.
%
% Contrast is initially set at 0.2 and, unbeknownst to the participant, is
% thresholded over time such that their performance is held at roughly 80%
% accuracy. Typically, this requires a lowering of contrast but we are
% already witnessing a broad range of contrasts for this threshold level.
%
% The rightmost plot represents the 96 trials worth of training under the dual-task
% "diverted attention" condition. Subjects are required to make two
% responses to this task:

%%
% # Did the peripheral gabor appear? (50% probability)
% # Did the letter 'T' appear in the central letter stimulus?

%%
% Contrast starts a small proportion higher than in the full attention
% condition and is thresholded such that performance on the peripheral task
% remains at roughly 80% accuracy despite the concurrent task.
%
% Judging from Sherman et al. (2015) we expect to see an effect of
% attention under dual-task conditions. This should be reflected by
% significantly higher contrast under diverted attention than full
% attention.

%% COMPARISON OF FULL vs. DIVERTED ATTENTION FOLLOWING TRAINING
%

% NB. within_subject_error.m is available from: github.com/julian-matthews/stats-tools

figure

subplot(1,2,1)
y = [mean(pCON.full_tr(:,end)), mean(pCON.divert_tr(:,end))];
% Attention is a fixed factor, consider this in the future
e = within_subject_error([pCON.full_tr(:,end),pCON.divert_tr(:,end)]);
subjerr = [(nanstd(pCON.full_tr(:,end))/sqrt(subnum)),...
    (nanstd(pCON.divert_tr(:,end))/sqrt(subnum))];
bar([1 2],y,'w')
hold on
errbar([1 2],y,e,'k','LineWidth',4) % Within-subject error
errbar([1.1,2.1],y,subjerr,'r','LineWidth',4); % SEM
box off
colormap('white')
set(gca,'TickDir','out','XTickLabel',{'Full','Diverted'},'XTick',[1 2])
xlim([0.25,2.75])
ylim([0,.25])
ylabel('Gabor Contrast');
xlabel('Attention Condition');
title('Gabor Contrast Comparison (post-training)');

subplot(1,2,2)

y = [mean(pCON.full(:,end)), mean(pCON.divert(:,end))];
% e = within_subject_error(pCON.full(:,end),pCON.divert(:,end));
e = within_subject_error([pCON.full(:,end),pCON.divert(:,end)]);
subjerr = [(nanstd(pCON.full(:,end))/sqrt(subnum)),...
    (nanstd(pCON.divert(:,end))/sqrt(subnum))];
bar([1 2],y,'w')
hold on
errbar([1 2],y,e,'k','LineWidth',4) % Within-subject error
errbar([1.1,2.1],y,subjerr,'r','LineWidth',4); % SEM
box off
colormap('white')
set(gca,'TickDir','out','XTickLabel',{'Full','Diverted'},'XTick',[1 2])
xlim([0.25,2.75])
ylim([0,.25])
% ylabel('Gabor Contrast');
xlabel('Attention Condition');
title('Gabor Contrast Comparison (post-experiment)');

%%
% Comparison of gabor contrast between attention conditions post-training and experiment.
% Error bars reflect standard error of the mean within-subjects (in black)
% and between-subjects (in red).
%
% We appear to see an effect of attention in contrast thresold
% emerging post-training but this is clarified by the end of the experiment
%
% Sherman et al. (2015) finished adjusting contrast at this point assuming
% training would no longer effect performance. I'm skeptical of this so in
% our version we performed 6 additional contrast thresholds during 50%
% expectation blocks (the same expectation level as training).
%
% Let's compare contrast levels at the end of the experiment...
%
% Comparison of attention condition gabor contrast post-experiment.
% Error bars reflect standard error of the mean within-subjects (in black)
% and between-subjects (in red).
%
% We see a stronger difference between the attention conditions which may
% stand up to (within-subject) significance testing. This reflects some
% degree of training on the tasks, I think we were right to continue
% thresholding.

% CONTRAST THRESHOLDS DURING EXPERIMENT
% Here we examine how contrast shifted for each subject during the
% experiment

figure
subplot(1,2,1)
y = pCON.full';
e = pCON.full_SEM';
c = cool(subnum);
er = errorbar(y,e,'o-');
for subject = 1:subnum
    set(er(subject),'color',c(subject,:))
end
box off
set(gca,'TickDir','out','XTick',1:6)
xlim([0.5,6.5])
ylim([0,.65])
ylabel('Gabor Contrast');
xlabel('Experimental Block');
title('Full Attention');
% legend(subjs,'Location','BestOutside');

subplot(1,2,2)
y = pCON.divert';
e = pCON.divert_SEM';
c = cool(subnum);
er = errorbar(y,e,'o-');
for subject = 1:subnum
    set(er(subject),'color',c(subject,:))
end
box off
set(gca,'TickDir','out','XTick',1:6)
xlim([0.5,6.5])
ylim([0,.65])
% ylabel('Gabor Contrast');
xlabel('Experimental Block');
title('Diverted Attention');
% legend(subjs,'Location','BestOutside');



%%
% For the sake of completion, here are the individual subject plots of
% contrast threshold over the experimental trials. Error bars reflect
% standard error within each block of trials.
%
% Contrast was adjusted in the respective full or diverted attention
% condition when expectations of gabor presence were set at 50%. Six blocks
% worth of trials were performed for each hence the change in x-axis.
%
% One of the noticeable details is that the two subjects that especially
% struggled with the task (high contrast level) seemed to improve under the
% full attention condition. Their contrast drops progressively as the
% experiment continues.

% PERFORMANCE

% Subject x block dimensions
clear peri_OP;
subnum = length(subjs);
for condition = 1:length(conditions)
    switch COND(condition).task
        case 'lo_full'
            errors = size(COND(condition).pCON,1);
            peri_OP.full_L = (reshape(mean(COND(condition).pCON),[],subnum))';
            peri_OP.full_LSEM = (reshape(std(COND(condition).pCON)/sqrt(errors),[],subnum))';
        case 'lo_divert'
            errors = size(COND(condition).pCON,1);
            peri_OP.divert_L = (reshape(mean(COND(condition).pCON),[],subnum))';
            peri_OP.divert_LSEM = (reshape(std(COND(condition).pCON)/sqrt(errors),[],subnum))';
        case 'med_full'
            errors = size(COND(condition).pCON,1);
            peri_OP.full_M = (reshape(mean(COND(condition).pCON),[],subnum))';
            peri_OP.full_MSEM = (reshape(std(COND(condition).pCON)/sqrt(errors),[],subnum))';
        case 'med_divert'
            errors = size(COND(condition).pCON,1);
            peri_OP.divert_M = (reshape(mean(COND(condition).pCON),[],subnum))';
            peri_OP.divert_MSEM = (reshape(std(COND(condition).pCON)/sqrt(errors),[],subnum))';
    end
end
