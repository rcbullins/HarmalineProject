function [] = Behavioral_Metrics(animals,BASEPATH,exper_conditions)
% PURPOSE
%   For the harmaline dataset, compare behavioral metrics accuracy and path
%   length for control sessions vs harmaline sessions.
%   4 conditions: Control (baseline vs stim)
%                 Harmaline (baseline vs stim)
% INPUT
%   animals
%       struct of animal names
%   BASEPATH
%       path to harmaline project code folder
%   exper_conditions
%       struct with experimental conditions (control or harmaline)
%   USER
%       user name in pathway, used in animal directory
% DEPENDENDECIES
%   Runs at the end of Behavioral Comprehensive Script
%       - calls upon 2 mat files created here
%             - PathLengths.mat
%             - baselineAccuracy.mat  (accuracy on baseline trials)
%             - pertAccuracy.mat      (accuracy on stimulation trials)
%             - BLIsolateAccuracy.mat (accuracy for only failures and
%                                      eventual success)
%             - pertIsolateAccuracy.mat (accuracy for only failures and
%                                      eventual success)
% OUTPUT
%   - Makes bar plot and line plot of accuracy over conditions for each
%     animal
%   - Makes bar plot for eventual success/ eventual success + failures
%   - Makes bar plot of accuracy over conditions, harmaline vs control only
% NOTES: may need to change code in behavioral comprehensive script to not
% consider grooming in accuracy measures
%% Load Directory

Directory_Animals;
%% Compare accuracy harmaline vs control (baseline vs stim trials)
% Loop through subjects
for isub = 1:length(animals)
    % Identify subject
    SUB = animals{isub};
    % Define path to store figures
    COMPARISON_FIGS = [BASEPATH 'Figures/' SUB '/Behavior/HarmVsControl/'];
    figure;
    % Initiate i = 0, this is the exper session num, reset each condition
    i = 0;
    % Initiate pathLengthMat to store all path lengths in to bar plot later
    % NOTE: to make bar plots, must have in same struct
    pathLengthMat = {};
    % Loop through experimental conditions (control vs harmaline)
    for iexper = 1:length(exper_conditions)
        EXPER_COND = exper_conditions{iexper};
        ExperSessions = eval(sprintf('%s_%sBehaviorVideos',SUB,EXPER_COND));
        % Loop through sessions
        for isession = 1:length(ExperSessions)
            EXPER_SESSION = ExperSessions{isession};
            if isempty(EXPER_SESSION)
                continue;
            end
            % Load path length for session
            load([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_PathLengths.mat'],'totalPathLength');
            % Add path length to struct
            pathLengthMat{1+(2*i)} = totalPathLength.nbase;
            pathLengthMat{2+(2*i)} = totalPathLength.npert;
            % Load accuracy for session and add to accuracy matrix
            load([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_baselineAccuracy.mat'],'y_BL_acc');
            AccMat(1+(2*i),:) = y_BL_acc*100;
            load([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_pertAccuracy.mat'],'y_Pert_acc');
            AccMat(2+(2*i),:) = y_Pert_acc*100;
            %load subset of data accuracy (eventual success/eventual
            %success + failures)
            load([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_BLIsolateAccuracy.mat'],'y_BL_acc_isolate');
            load([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_pertIsolateAccuracy.mat'],'y_Pert_acc_isolate');
            AccIsolateMat(1+(2*i),:) = y_BL_acc_isolate*100;
            AccIsolateMat(2+(2*i),:) = y_Pert_acc_isolate*100;
            % Update index of experimental session
            i = i +1;
        end % experiment session
    end % experiment condition

    %% Plot Bar chart accuracy comparing harmaline vs control (baseline vs stim)
    h = bar(AccMat);
    title([ SUB ': Accuracy']);
    ylabel('Accuracy (%)');
    xticklabels({'Control Baseline';'Control Stim';'Harmaline Baseline';'Harmaline Stim'});
    l = cell(1,4);
    l{1}='First Reach Success'; l{2}='Eventual Success'; l{3}='No Success'; l{4}='No Reach';
    legend(h,l);
    savefig([COMPARISON_FIGS SUB '_BarplotAccuracy.fig']);
    %% Plot line chart accuracy
    figure;
    plot(1:length(AccMat(:,1)),AccMat(:,1));
    hold on;
    plot(1:length(AccMat(:,2)),AccMat(:,2));
    plot(1:length(AccMat(:,3)),AccMat(:,3));
    plot(1:length(AccMat(:,4)),AccMat(:,4));
    %plot(1:length(AccMat(:,5)),AccMat(:,5));
    title([ SUB ': Accuracy']);
    ylabel('Accuracy (%)');
    xticks([1 2 3 4]);
    xticklabels({'Control Baseline';'Control Stim';'Harmaline Baseline';'Harmaline Stim'});
    legend('First Reach Success','Eventual Success','No Success','No Reach');
    savefig([COMPARISON_FIGS SUB '_LinePlotAccuracy.fig']);
    %% Plot path length
    figure();
%     A1 =  pathLengthMat{1};
%     A2 =  pathLengthMat{2};
%     A3 =  pathLengthMat{3};
%     A4 =  pathLengthMat{4};
%     G = [ones(size(A1))  2*ones(size(A3)) 3*ones(size(A4)) 4*ones(size(A4))];
%     X = [A1, A2, A3, A4];
%     boxplot(X,G,'notch','on','colors',[0 0 0],'symbol','','labels',{'data1','data2','data3','data4'});
%     title([ SUB ':Path Length']);
%     bar2plot = [mean(A1) mean(A3)];
%     bar(bar2plot);
    %% Plot accuracy for eventual success and failures
    bar(AccIsolateMat);
    title([ SUB ': Eventual Accuracy']);
    ylabel('Eventual Success/(eventual + failures)(%)');
    xticklabels({'Control Baseline';'Control Stim';'Harmaline Baseline';'Harmaline Stim'});
    %% Updated Plot (Baseline:Control vs Harmaline) 
    % x axis is Eventual Success, no Reach, etc
    % Define control
    control_acc = AccMat(1,1:4);
    % eventual success should be eventual/eventual + failures
    control_acc(1,2) = AccIsolateMat(1);
    % Define harmaline
    harm_acc = AccMat(3,1:4);
    % eventual success shouldbe eventual/eventual + failures
    harm_acc(1,2) = AccIsolateMat(3);
    % m groups by n bars (m x n)
    % so (M = Success, evetual, no success, no reach)
    % n = 2 (control vs harmaline)
    % Control
    compMatBar(:,1) = control_acc(1,:);
    % harmaline
    compMatBar(:,2) = harm_acc(1,:);
   figure;
    h = bar(compMatBar) %, 'FaceColor','flat');
    hold on;
%     h(1).CData(1,:) = [0 0 1];
%     h(1).CData(2,:) = [0 0 1];
%     h(1).CData(3,:) = [0 0 1];
%     h(1).CData(4,:) = [0 0 1];
%     h(2).CData(1,:) = [1 0 0];
%     h(2).CData(2,:) = [1 0 0];
%     h(2).CData(3,:) = [1 0 0];
%     h(2).CData(4,:) = [1 0 0];    
%     h(1).FaceAlpha = .6;
%     h(2).FaceAlpha = .6;
%              
%                 l = cell(1,2);
%                 l{1}='Control'; l{2}='Harmaline';
if strcmp(SUB,'M340')
    title('Reaching during harmaline (10mg/kg)');
elseif strcmp(SUB,'M341')
    title('Reaching during harmaline (20mg/kg)');
end
ylim([0 100]);
    xticklabels({'First Reach Success';'Eventual Success';'No Success';'No Reach'});
    legend('Control','Harmaline');
     ylabel('Accuracy (%)');
    clear('AccMat');
%% Updated plot 2 (no success not included: first reach vs eventual vs no reach)
% Take out no success bar 
    % Control and harmaline
    figure;
    compMatBar(3,:) = [];
    h = bar(compMatBar) %, 'FaceColor','flat');
    hold on;
if strcmp(SUB,'M340')
    title('Reaching during harmaline (10mg/kg)');
elseif strcmp(SUB,'M341')
    title('Reaching during harmaline (20mg/kg)');
end
ylim([0 100]);
    xticklabels({'First Reach Success';'Eventual Success';'No Reach'});
    legend('Control','Harmaline');
     ylabel('Accuracy (%)');
clear('compMatBar');
end % subject
end

