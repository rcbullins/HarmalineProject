% Behavior script - all behavior
clear all;
close ALL;
clc;
%% Specify what to plot
plotSampleTrajectories = 0;
plot1DOverlayTrajectories = 0;
plot3DOverlayTrajectories = 0;
plotAccuracy = 0;
score = 'all'; %1, 0, 2, -1, 'all'
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
%% Prep which experimental condition to analyze
exper_conditions = {'control';'harm'};
%% Score
if score == 1
    SCORE = 'idealSuccess';
elseif score == 0
    SCORE = 'noSuccess';
elseif score == -1
    SCORE = 'noReach';
elseif score == 2
    SCORE = 'eventualSuccess';
elseif score == 'all'
    score = [1 0 2];
    SCORE = 'allScores'
end

%% Initate CSV to hold metrics in
behavioralFile = [BASEPATH 'Data_Analyzed\behaviorMetrics.csv'];
if exist(behavioralFile,'file')==2
    delete(behavioralFile);
end
% Open the CSV file
allBehavFID = fopen(behavioralFile,'w');
% Make labels on columns: Subject, sustained, selective, accuracy, and
% reaction time
fprintf(allBehavFID,'subject,session,condition,trialIdx,trialScore,startFrame,endFrame, reactionTime, pathLength\n');
%fprintf(allBehavFID,'%s,sus,sel,%0.5f,%0.5f\n',subject,selSus_acc,selSus_rt);
%%
for isub = 1:length(animals)
    % Find how many sessions for animal in experimental condition of choice
    % (control or harmaline)
    SUB = animals{isub};
    for iexper = 1:length(exper_conditions)
        EXPER_COND = exper_conditions{iexper};
        ExperSessions = eval(sprintf('%s_%sBehaviorVideos',SUB,EXPER_COND));
        % Add path for figures
        BEHAV_FIG = [BASEPATH 'Figures\' SUB '\Behavior\'];
        for isession = 1:length(ExperSessions)
            EXPER_SESSION = ExperSessions{isession};
            if isempty(EXPER_SESSION)
                continue;
            end
            TRK = [BASEPATH 'Data/' SUB 'necab1_Chr2/' EXPER_SESSION 'movies/trk'];
            addpath(TRK);
            %% Load Data
            nframes=2500; %variable! %number of frames to analyze. Going to depend on how many frames in movie
            [traj conf] = get_traj_3D_v2(TRK,CODE_CALIB,nframes);
            % traj (trials x coordinate xyz x frames); % x is
            % dim1,vertical
            % conf is the confidence value assigned by the classifier
            
            %% Set trials
            % Define Trials by baseline, stim, washout M340_20210811_controlTrials.nbase
            trialIdxs = eval(sprintf('%s_%s_%sTrials',SUB,EXPER_SESSION,EXPER_COND));
            base.trialIdxs = trialIdxs.nbase;
            pert.trialIdxs = trialIdxs.npert;
            wash.trialIdxs = trialIdxs.nwash;
            %% Assign trials based on trial type
            % traj(trials x coordinate xyz x frames);
            base.traj = traj(base.trialIdxs,:,:);
            pert.traj = traj(pert.trialIdxs,:,:);
            wash.traj = traj(wash.trialIdxs,:,:);
            %% Set Graph Defaults
            SetGraphDefaults;
            %% Get num of success, no success, etc for each baseline,stim, and wash
            % baseline accuracy
            numBL.idealSuccess = length(find(trialIdxs.trialScore(base.trialIdxs) == 1));
            numBL.eventualSuccess = length(find(trialIdxs.trialScore(base.trialIdxs) == 2));
            numBL.noSuccess = length(find(trialIdxs.trialScore(base.trialIdxs) == 0));
            numBL.noReach = length(find(trialIdxs.trialScore(base.trialIdxs) == -1));
            numBL.grooming = length(find(trialIdxs.trialScore(base.trialIdxs) == 'g'));
            % stim accuracy
            numPert.idealSuccess = length(find(trialIdxs.trialScore(pert.trialIdxs)==1));
            numPert.eventualSuccess = length(find(trialIdxs.trialScore(pert.trialIdxs) == 2));
            numPert.noSuccess = length(find(trialIdxs.trialScore(pert.trialIdxs) == 0));
            numPert.noReach = length(find(trialIdxs.trialScore(pert.trialIdxs) == -1));
            numPert.grooming = length(find(trialIdxs.trialScore(pert.trialIdxs) == 'g'));
            % washout baseline
            numWash.idealSuccess = length(find(trialIdxs.trialScore(wash.trialIdxs)==1));
            numWash.eventualSuccess = length(find(trialIdxs.trialScore(wash.trialIdxs) == 2));
            numWash.noSuccess = length(find(trialIdxs.trialScore(wash.trialIdxs) == 0));
            numWash.noReach = length(find(trialIdxs.trialScore(wash.trialIdxs) == -1));
            numWash.grooming = length(find(trialIdxs.trialScore(wash.trialIdxs) == 'g'));
            %% Plot accuracy for baseline and save mat
            numTotalBL = length(base.trialIdxs);
            y_BL = [numBL.idealSuccess; numBL.eventualSuccess ;numBL.noSuccess ;numBL.noReach; numBL.grooming]
            y_BL_acc = y_BL./numTotalBL;
            save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_baselineAccuracy.mat'],'y_BL_acc');
            if plotAccuracy
                h = bar(y_BL_acc*100,'FaceColor','flat');
                colorScheme = colorGradient([1,0,.5],[.5,0,1],5);
                h.CData(1,:) = colorScheme(1,:);
                h.CData(2,:) = colorScheme(2,:);
                h.CData(3,:) = colorScheme(3,:);
                h.CData(4,:) = colorScheme(4,:);
                h.CData(5,:) = colorScheme(5,:);
                l = cell(1,5);
                l{1}='Success'; l{2}='Eventual Success'; l{3}='No Success'; l{4}='No Reach'; l{5}='Grooming';
                set(gca,'xticklabel', l);
                ylabel('Accuracy (%)');
                title([SUB ' ' EXPER_SESSION ' ' EXPER_COND ' :BL Accuracy']);
            end
            %% Plot stim accruacy and save mat
            numTotalPert = length(pert.trialIdxs);
            y_Pert = [numPert.idealSuccess; numPert.eventualSuccess ;numPert.noSuccess ;numPert.noReach; numPert.grooming]
            y_Pert_acc = y_Pert./numTotalPert;
            save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_pertAccuracy.mat'],'y_Pert_acc');
            if plotAccuracy
                figure;
                h = bar(y_Pert_acc*100,'FaceColor','flat');
                colorScheme = colorGradient([1,0,.5],[.5,0,1],5);
                h.CData(1,:) = colorScheme(1,:);
                h.CData(2,:) = colorScheme(2,:);
                h.CData(3,:) = colorScheme(3,:);
                h.CData(4,:) = colorScheme(4,:);
                h.CData(5,:) = colorScheme(5,:);
                l = cell(1,5);
                l{1}='Success'; l{2}='Eventual Success'; l{3}='No Success'; l{4}='No Reach'; l{5}='Grooming';
                set(gca,'xticklabel', l);
                ylabel('Accuracy (%)');
                title([SUB ' ' EXPER_SESSION ' ' EXPER_COND ' : Stim Accuracy']);
                %            savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajector/' SUB '_' EXPER_SESSION '_' EXPER_COND '_zDim_washout.fig']);
            end
            %% Plot Trajectory examples
            if plotSampleTrajectories
                %% Plot Sample 1D Trajectories in Z
                % INPUT PARSERS AVAILABLE: coordinate_dim, startFrame, endFrame
                plot1DTrajectories(base.traj,'endFrame',nframes);
                sgtitle([SUB ' ' EXPER_COND ': handpaths baseline trials ZDim']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_zDim_baseline.fig']);
                plot1DTrajectories(pert.traj, 'endFrame',nframes);
                sgtitle([SUB ' ' EXPER_COND ': handpaths stimulation trials ZDim']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_zDim_stim.fig']);
                plot1DTrajectories(wash.traj, 'endFrame',nframes);
                sgtitle([SUB ' ' EXPER_COND ': handpaths washout trials ZDim']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_zDim_washout.fig']);
                %% Plot Sample 1D Trajectories in y
                % INPUT PARSERS AVAILABLE: coordinate_dim, startFrame, endFrame
                plot1DTrajectories(base.traj,'endFrame',nframes,'coordinate_dim',2);
                sgtitle([SUB ' ' EXPER_COND ': handpaths baseline trials YDim']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_yDim_baseline.fig']);
                plot1DTrajectories(pert.traj, 'endFrame',nframes,'coordinate_dim',2);
                sgtitle([SUB ' ' EXPER_COND ': handpaths stimulation trials YDim']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_yDim_stim.fig']);
                plot1DTrajectories(wash.traj, 'endFrame',nframes,'coordinate_dim',2);
                sgtitle([SUB ' ' EXPER_COND ': handpaths washout trials YDim']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_yDim_washout.fig']);
                
                %% Plot 3D Trajectories
                plot3DTrajectories(base.traj, 'endFrame',nframes);
                sgtitle([SUB ' ' EXPER_COND ': 3D handpaths baseline trials']);
                savefig([BASEPATH 'Figures\' SUB '\Behavior\3DTrajectories\' SUB '_' EXPER_SESSION '_' EXPER_COND '_baseline.fig']);
                plot3DTrajectories(pert.traj, 'endFrame',nframes);
                sgtitle([SUB ' ' EXPER_COND ': 3D handpaths stimulation trials']);
                savefig([BASEPATH 'Figures\' SUB '\Behavior\3DTrajectories\' SUB '_' EXPER_SESSION '_' EXPER_COND '_stim.fig']);
                plot3DTrajectories(wash.traj, 'endFrame',nframes);
                sgtitle([SUB ' ' EXPER_COND ': 3D handpaths washout trials']);
                savefig([BASEPATH 'Figures\' SUB '\Behavior\3DTrajectories\' SUB '_' EXPER_SESSION '_' EXPER_COND '_washout.fig']);
            end
            %% Find start frame of movement
            
            % Only for trials with a score of 1 during baseline and stimulation
            [movStrt.nbase] = getMovStartFrame(score,traj,trialIdxs.trialScore,trialIdxs.nbase);
            [movStrt.npert] = getMovStartFrame(score,traj,trialIdxs.trialScore,trialIdxs.npert);
            [movStrt.nwash] = getMovStartFrame(score,traj,trialIdxs.trialScore,trialIdxs.nwash);
            %% Find end frame of reach movement
            % only for trials with a score of 1 suring baseline and stimulation
            [movEnd.nbase,movIdx.nbase] = getMovEndFrame(score,traj,trialIdxs.trialScore,trialIdxs.nbase, movStrt.nbase);
            [movEnd.npert,movIdx.npert] = getMovEndFrame(score,traj,trialIdxs.trialScore,trialIdxs.npert, movStrt.npert);
            [movEnd.nwash,movIdx.nwash] = getMovEndFrame(score,traj,trialIdxs.trialScore,trialIdxs.nwash, movStrt.nwash);
            %% Plot Overlay 1D - only reach
            if plot1DOverlayTrajectories
                dimPlot = 3; % if change add dimPlot input to plotOverlay1D function
                % Plot overlay for baseline and stim trials
                plotOverlay1D(movStrt.nbase, movEnd.nbase,movIdx.nbase,traj);
                title([SUB ' ' EXPER_COND ': Overlay ' SCORE ' baseline Trials (Dim ' num2str(dimPlot) ')']);
                savefig([BASEPATH 'Figures\' SUB '\Behavior\1DTrajectories\' SCORE '\' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE '_overlayReachBaseline_dim' num2str(dimPlot) '.fig']);
                
                plotOverlay1D(movStrt.npert, movEnd.npert,movIdx.npert,traj);
                title([SUB ' ' EXPER_COND ': Overlay' SCORE ' stim Trials (Dim ' num2str(dimPlot) ')']);
                savefig([BASEPATH 'Figures\' SUB '\Behavior\1DTrajectories\' SCORE '\' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE '_overlayReachStim_dim' num2str(dimPlot) '.fig']);
                
                % PRE and POST baseline
                movStrtBL = [movStrt.nbase movStrt.nwash];
                movEndBL = [movEnd.nbase movEnd.nwash];
                movStrtIdxs = [movIdx.nbase movIdx.nwash]';
                plotOverlay1D(movStrtBL, movEndBL,movStrtIdxs,traj);
                title([SUB ' ' EXPER_COND ': Overlay ' SCORE 'pre and post stim (Dim ' num2str(dimPlot) ')']);
                savefig([BASEPATH 'Figures\' SUB '\Behavior\1DTrajectories\' SCORE '\' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE '_overlayReachAllBL_dim' num2str(dimPlot) '.fig']);
                %% Plot Overlay 1D - whole movement
                % Plot baseline pre stim trials
                plotOverlay1D(movStrt.nbase, movEnd.nbase,movIdx.nbase, traj,'PlotEntireReach',1);
                title([SUB ' ' EXPER_COND ': Overlay ' SCORE ' baseline Trials (Dim ' num2str(dimPlot) ')']);
                savefig([BASEPATH 'Figures\' SUB '\Behavior\1DTrajectories\' SCORE '\' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE '_overlayBaseline_dim' num2str(dimPlot) '.fig']);
                % Plot stim
                plotOverlay1D(movStrt.npert, movEnd.npert,movIdx.npert,traj,'PlotEntireReach',1);
                title([SUB ' ' EXPER_COND ': Overlay ' SCORE ' stim Trials (Dim ' num2str(dimPlot) ')']);
                savefig([BASEPATH 'Figures\' SUB '\Behavior\1DTrajectories\' SCORE '\' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE '_overlayStim_dim' num2str(dimPlot) '.fig']);
                % Plot baseline pre stim and post stim
                movStrtBL = [movStrt.nbase movStrt.nwash];
                movEndBL = [movEnd.nbase movEnd.nwash];
                movStrtIdxs = [movIdx.nbase movIdx.nwash]';
                plotOverlay1D(movStrtBL, movEndBL,movStrtIdxs,traj,'PlotEntireReach',1);
                title([SUB ' ' EXPER_COND ': Overlay ' SCORE ' pre and post stim (Dim ' num2str(dimPlot) ')']);
                savefig([BASEPATH 'Figures\' SUB '\Behavior\1DTrajectories\' SCORE '\' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE 'overlayAllBL_dim' num2str(dimPlot) '.fig']);
            end
            %% Plot Overlay 3D - reach movement
            save([BASEPATH 'Data_Analyzed\' SUB '\Behavior\' SCORE '\' SUB '_' EXPER_SESSION '_' EXPER_COND '_Overlay3DVariables.mat'],'movStrt', 'movEnd','movIdx', 'traj');
            if plot3DOverlayTrajectories
                plotOverlay3D(movStrt.nbase, movEnd.nbase,movIdx.nbase, traj);
                title([SUB ' ' EXPER_COND ': Overlay ' SCORE ' baseline Trials']);
                savefig([BASEPATH 'Figures\' SUB '\Behavior\3DTrajectories\' SCORE '\' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE '_overlayReachBaseline.fig']);
                % Plot post and pre stim baseline
                movStrtBL = [movStrt.nbase movStrt.nwash];
                movEndBL = [movEnd.nbase movEnd.nwash];
                movStrtIdxs = [movIdx.nbase movIdx.nwash]';
                plotOverlay3D(movStrtBL, movEndBL,movStrtIdxs,traj);
                title([SUB ' ' EXPER_COND ': Overlay ' SCORE ' pre and post stim']);
                savefig([BASEPATH 'Figures\' SUB '\Behavior\3DTrajectories\' SCORE '\' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE 'overlayReachAllBL.fig']);
            end
            %% Get path length
            % diff(Traj(start:stop))
            %   movStrt.nbase, movEnd.nbase,movIdx.nbase,traj
            % get path lengths for baseline trials
            for itrial = 1:length(movIdx.nbase)
                totalLength = 0;
                frameStart = movStrt.nbase(1,itrial);
                frameEnd = movEnd.nbase(1,itrial);
                totalFramesTraj = frameEnd-frameStart;
                for iframe = 1:totalFramesTraj
                  iframeTraj = sqrt(sum((traj(movIdx.nbase(itrial),:,frameStart-1+iframe)) .^ 2));
                  totalLength = totalLength + iframeTraj;
                end
                totalPathLength.nbase(itrial) = totalLength;
            end
            % for stimulation trials
            for itrial = 1:length(movIdx.npert)
                totalLength = 0;
                frameStart = movStrt.npert(1,itrial);
                frameEnd = movEnd.npert(1,itrial);
                totalFramesTraj = frameEnd-frameStart;
                for iframe = 1:totalFramesTraj
                  iframeTraj = sqrt(sum((traj(movIdx.npert(itrial),:,frameStart-1+iframe)) .^ 2));
                  totalLength = totalLength + iframeTraj;
                end
                totalPathLength.npert(itrial) = totalLength;
            end
            %for washout baseline trials trials
            for itrial = 1:length(movIdx.nwash)
                totalLength = 0;
                frameStart = movStrt.nwash(1,itrial);
                frameEnd = movEnd.nwash(1,itrial);
                totalFramesTraj = frameEnd-frameStart;
                for iframe = 1:totalFramesTraj
                  iframeTraj = sqrt(sum((traj(movIdx.nwash(itrial),:,frameStart-1+iframe)) .^ 2));
                  totalLength = totalLength + iframeTraj;
                end
                totalPathLength.nwash(itrial) = totalLength;
            end  
            save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_PathLengths.mat'],'totalPathLength');
            end %session
        end %exper conditions
        
        % fprintf(allBehavFid,'%s,%s,%s,%i,%i,%i,%i,%i\n',SUB,EXPER_SESSION,EXPER_COND,trialIdx,trialScore,startFrame,endFrame,pathLength);
    end % animal
    % Compare 3D trajectories, harm vs control
    Behavioral_Comparison_Conditions;
    
    % Compare accuracy across conditions
    Behavioral_Metrics;
