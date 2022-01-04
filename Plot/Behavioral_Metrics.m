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
    l = cell(1,5);
    l{1}='Success'; l{2}='Eventual Success'; l{3}='No Success'; l{4}='No Reach'; l{5}='Grooming';
    legend(h,l);
    savefig([COMPARISON_FIGS SUB '_BarplotAccuracy.fig']);
    %% Plot line chart accuracy
    figure;
    plot(1:length(AccMat(:,1)),AccMat(:,1));
    hold on;
    plot(1:length(AccMat(:,2)),AccMat(:,2));
    plot(1:length(AccMat(:,3)),AccMat(:,3));
    plot(1:length(AccMat(:,4)),AccMat(:,4));
    plot(1:length(AccMat(:,5)),AccMat(:,5));
    title([ SUB ': Accuracy']);
    ylabel('Accuracy (%)');
    xticks([1 2 3 4]);
    xticklabels({'Control Baseline';'Control Stim';'Harmaline Baseline';'Harmaline Stim'});
    legend('Success','Eventual Success','No Success','No Reach','Grooming');
    savefig([COMPARISON_FIGS SUB '_LinePlotAccuracy.fig']);
    clear('AccMat');
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
end % subject
end

