<<<<<<< HEAD
<<<<<<< HEAD
% Behavioral metrics script (Accuracy, path length)
% Find accuracy for control vs harmaline
%   4 conditions: Control (baseline vs stim)
%                 Harmaline (baseline vs stim)
%   Makes bar plot and line plot of accuracy over conditions for each
%   animal
%% Add code paths
BASEPATH = 'C:/Users/rcbul/OneDrive - University of North Carolina at Chapel Hill/Hantman_Lab/Harmaline_Project/';
CODE_REAGAN = [BASEPATH 'Code/reagan_code/'];
CODE_CALIB = [BASEPATH 'Code/britton_code/calibration_files/camera_calibration_Jay_7-28-16/Calib_Results_stereo.mat'];
CODE_TRKR = [BASEPATH 'Code/britton_code/tracker'];
CODE_BRITTON = [BASEPATH 'Code/britton_code/code'];
CODE_BRITTON_PLOT = [BASEPATH 'Code/britton_code/other_code'];
addpath(genpath(CODE_REAGAN));
addpath(CODE_CALIB);
addpath(CODE_TRKR);
addpath(genpath(CODE_BRITTON));
addpath(genpath(CODE_BRITTON_PLOT));

Directory_Animals;
%%
%% Prep which experimental condition to analyze
exper_conditions = {'control','harm'};
%% Run
for isub = 1:length(animals)
    SUB = animals{isub};
    COMPARISON_FIGS = [BASEPATH 'Figures/' SUB '/Behavior/HarmVsControl/'];
    figure;
    i = 0;
    pathLengthMat = {};
    for iexper = 1:length(exper_conditions)
        EXPER_COND = exper_conditions{iexper};
        ExperSessions = eval(sprintf('%s_%sBehaviorVideos',SUB,EXPER_COND));
        for isession = 1:length(ExperSessions)
            EXPER_SESSION = ExperSessions{isession};
            if isempty(EXPER_SESSION)
                continue;
            end
            % Load path length
            load([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_PathLengths.mat']);
            pathLengthMat{1+(2*i)} = totalPathLength.nbase;
            pathLengthMat{2+(2*i)} = totalPathLength.npert;
%             % load accuracy
            load([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_baselineAccuracy.mat']);
            AccMat(1+(2*i),:) = y_BL_acc*100;
            load([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_pertAccuracy.mat']);
            AccMat(2+(2*i),:) = y_Pert_acc*100;
           
            AccIsolateMat(1+(2*i),:) = y_BL_acc_isolate*100;
            AccIsolateMat(2+(2*i),:) = y_Pert_acc_isolate*100;
          
             i = i +1;
        end % experiment session
    end % experiment condition
    % MAKE BAR CHART Accuracy
    h = bar(AccMat);
    title([ SUB ': Accuracy']);
    ylabel('Accuracy (%)');
    xticklabels({'Control Baseline';'Control Stim';'Harmaline Baseline';'Harmaline Stim'});
    l = cell(1,5);
    l{1}='Success'; l{2}='Eventual Success'; l{3}='No Success'; l{4}='No Reach'; l{5}='Grooming';
    legend(h,l);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/HarmVsControl/' SUB '_BarplotAccuracy.fig']);
    % MAKE LINE CHART Accuracy
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
    savefig([BASEPATH 'Figures/' SUB '/Behavior/HarmVsControl/' SUB '_LinePlotAccuracy.fig']);
    clear('AccMat');
    % MAKE PATH LENGTH
    figure();
    A1 =  pathLengthMat{1};
    A2 =  pathLengthMat{2};
    A3 =  pathLengthMat{3};
    A4 =  pathLengthMat{4};
    G = [ones(size(A1))  2*ones(size(A3)) 3*ones(size(A4)) 4*ones(size(A4))];
    X = [A1, A2, A3, A4];
    boxplot(X,G,'notch','on','colors',[0 0 0],'symbol','','labels',{'data1','data2','data3','data4'});
    title([ SUB ':Path Length']);
    bar2plot = [mean(A1) mean(A3)];
    bar(bar2plot);
  % MAKE ACCURACY PERCENTAGE FOR eventual success and failures
    h = bar(AccIsolateMat);
    title([ SUB ': Eventual Accuracy']);
    ylabel('Eventual Success/(eventual + failures)(%)');
    xticklabels({'Control Baseline';'Control Stim';'Harmaline Baseline';'Harmaline Stim'});
end

=======
% Behavioral metrics script (Accuracy, path length)
% Find accuracy for control vs harmaline
%   4 conditions: Control (baseline vs stim)
%                 Harmaline (baseline vs stim)
%   Makes bar plot and line plot of accuracy over conditions for each
%   animal
%% Add code paths
BASEPATH = 'C:/Users/rcbul/OneDrive - University of North Carolina at Chapel Hill/Hantman_Lab/Harmaline_Project/';
CODE_REAGAN = [BASEPATH 'Code/reagan_code/'];
CODE_CALIB = [BASEPATH 'Code/britton_code/calibration_files/camera_calibration_Jay_7-28-16/Calib_Results_stereo.mat'];
CODE_TRKR = [BASEPATH 'Code/britton_code/tracker'];
CODE_BRITTON = [BASEPATH 'Code/britton_code/code'];
CODE_BRITTON_PLOT = [BASEPATH 'Code/britton_code/other_code'];
addpath(genpath(CODE_REAGAN));
addpath(CODE_CALIB);
addpath(CODE_TRKR);
addpath(genpath(CODE_BRITTON));
addpath(genpath(CODE_BRITTON_PLOT));

Directory_Animals;
%%
%% Prep which experimental condition to analyze
exper_conditions = {'control','harm'};
%% Run
for isub = 1:length(animals)
    SUB = animals{isub};
    COMPARISON_FIGS = [BASEPATH 'Figures/' SUB '/Behavior/HarmVsControl/'];
    figure;
    i = 0;
    pathLengthMat = {};
    for iexper = 1:length(exper_conditions)
        EXPER_COND = exper_conditions{iexper};
        ExperSessions = eval(sprintf('%s_%sBehaviorVideos',SUB,EXPER_COND));
        for isession = 1:length(ExperSessions)
            EXPER_SESSION = ExperSessions{isession};
            if isempty(EXPER_SESSION)
                continue;
            end
            % Load path length
            load([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_PathLengths.mat']);
            pathLengthMat{1+(2*i)} = totalPathLength.nbase;
            pathLengthMat{2+(2*i)} = totalPathLength.npert;
%             % load accuracy
            load([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_baselineAccuracy.mat']);
            AccMat(1+(2*i),:) = y_BL_acc*100;
            load([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_pertAccuracy.mat']);
            AccMat(2+(2*i),:) = y_Pert_acc*100;
            i = i +1;
            
        end % experiment session
    end % experiment condition
    % MAKE BAR CHART Accuracy
    h = bar(AccMat);
    title([ SUB ': Accuracy']);
    ylabel('Accuracy (%)');
    xticklabels({'Control Baseline';'Control Stim';'Harmaline Baseline';'Harmaline Stim'});
    l = cell(1,5);
    l{1}='Success'; l{2}='Eventual Success'; l{3}='No Success'; l{4}='No Reach'; l{5}='Grooming';
    legend(h,l);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/HarmVsControl/' SUB '_BarplotAccuracy.fig']);
    % MAKE LINE CHART Accuracy
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
    savefig([BASEPATH 'Figures/' SUB '/Behavior/HarmVsControl/' SUB '_LinePlotAccuracy.fig']);
    clear('AccMat');
    % MAKE PATH LENGTH
    figure();
    A1 =  pathLengthMat{1};
    A2 =  pathLengthMat{2};
    A3 =  pathLengthMat{3};
    A4 =  pathLengthMat{4};
    G = [ones(size(A1))  2*ones(size(A3)) 3*ones(size(A4)) 4*ones(size(A4))];
    X = [A1, A2, A3, A4];
    boxplot(X,G,'notch','on','colors',[0 0 0],'symbol','','labels',{'data1','data2','data3','data4'});
    title([ SUB ':Path Length']);
    bar2plot = [mean(A1) mean(A3)];
    bar(bar2plot);
  
    
end

>>>>>>> da7b9ddee33fff3453e6f5e0546d20d35e92172a
=======
% Behavioral metrics script (Accuracy, path length)
% Find accuracy for control vs harmaline
%   4 conditions: Control (baseline vs stim)
%                 Harmaline (baseline vs stim)
%   Makes bar plot and line plot of accuracy over conditions for each
%   animal
%% Add code paths
BASEPATH = 'C:/Users/rcbul/OneDrive - University of North Carolina at Chapel Hill/Hantman_Lab/Harmaline_Project/';
CODE_REAGAN = [BASEPATH 'Code/reagan_code/'];
CODE_CALIB = [BASEPATH 'Code/britton_code/calibration_files/camera_calibration_Jay_7-28-16/Calib_Results_stereo.mat'];
CODE_TRKR = [BASEPATH 'Code/britton_code/tracker'];
CODE_BRITTON = [BASEPATH 'Code/britton_code/code'];
CODE_BRITTON_PLOT = [BASEPATH 'Code/britton_code/other_code'];
addpath(genpath(CODE_REAGAN));
addpath(CODE_CALIB);
addpath(CODE_TRKR);
addpath(genpath(CODE_BRITTON));
addpath(genpath(CODE_BRITTON_PLOT));

Directory_Animals;
%%
%% Prep which experimental condition to analyze
exper_conditions = {'control','harm'};
%% Run
for isub = 1:length(animals)
    SUB = animals{isub};
    COMPARISON_FIGS = [BASEPATH 'Figures/' SUB '/Behavior/HarmVsControl/'];
    figure;
    i = 0;
    pathLengthMat = {};
    for iexper = 1:length(exper_conditions)
        EXPER_COND = exper_conditions{iexper};
        ExperSessions = eval(sprintf('%s_%sBehaviorVideos',SUB,EXPER_COND));
        for isession = 1:length(ExperSessions)
            EXPER_SESSION = ExperSessions{isession};
            if isempty(EXPER_SESSION)
                continue;
            end
            % Load path length
            load([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_PathLengths.mat']);
            pathLengthMat{1+(2*i)} = totalPathLength.nbase;
            pathLengthMat{2+(2*i)} = totalPathLength.npert;
%             % load accuracy
            load([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_baselineAccuracy.mat']);
            AccMat(1+(2*i),:) = y_BL_acc*100;
            load([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_pertAccuracy.mat']);
            AccMat(2+(2*i),:) = y_Pert_acc*100;
            i = i +1;
            
        end % experiment session
    end % experiment condition
    % MAKE BAR CHART Accuracy
    h = bar(AccMat);
    title([ SUB ': Accuracy']);
    ylabel('Accuracy (%)');
    xticklabels({'Control Baseline';'Control Stim';'Harmaline Baseline';'Harmaline Stim'});
    l = cell(1,5);
    l{1}='Success'; l{2}='Eventual Success'; l{3}='No Success'; l{4}='No Reach'; l{5}='Grooming';
    legend(h,l);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/HarmVsControl/' SUB '_BarplotAccuracy.fig']);
    % MAKE LINE CHART Accuracy
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
    savefig([BASEPATH 'Figures/' SUB '/Behavior/HarmVsControl/' SUB '_LinePlotAccuracy.fig']);
    clear('AccMat');
    % MAKE PATH LENGTH
    figure();
    A1 =  pathLengthMat{1};
    A2 =  pathLengthMat{2};
    A3 =  pathLengthMat{3};
    A4 =  pathLengthMat{4};
    G = [ones(size(A1))  2*ones(size(A3)) 3*ones(size(A4)) 4*ones(size(A4))];
    X = [A1, A2, A3, A4];
    boxplot(X,G,'notch','on','colors',[0 0 0],'symbol','','labels',{'data1','data2','data3','data4'});
    title([ SUB ':Path Length']);
    bar2plot = [mean(A1) mean(A3)];
    bar(bar2plot);
  
    
end

>>>>>>> a39e1c0843197dbc2ab958057a9f139f35e0535d
