% Behavioral_Comparison_Conditions
% PURPOSE
%   Compare 3D trajectories and endpoints for reach and grab movements for
%   control subjects vs harmaline subjects. Want to see differences/similarities in
%   behavior for these two groups.
% INPUTS
%   animals
%       struct of animal names
%   BASEPATH
%       path to harmaline project code folder
%   exper_conditions
%       struct with experimental conditions (control or harmaline)
%   SCORE
%       string label for what types of trials to look at
%   USER
%       directory user name
% OUTPUTS
%       3D trajectories (control vs harmaline)
%       3D endpoints (control vs harmaline)
%       3D endpoints aligned, success vs failure trial
% DEPENDENCIES
%       Runs at the end of Behavioral Comprehensive Script
%           - calls upon 2 mat files created here
%               - overlay3Dvariables.mat
%               - pelletTraj.mat
% HISTORY
%   11.23.2021 Reagan Bullins
clear;
%% Variables
block = 'nbase';% 'npert', 'nwash'
score = 'all'; % Options: 1, 0, 2, -1, 'all'
               % Code: (1) one grab and success
               %       (0) grab and failure
               %       (2) multiple reaches and eventual success
               %      (-1) no reach attempts
               %   ('all') all scores where some attempt was made
%% Add code paths
USER = 'bullinsr';
RAWDATA_BASEPATH = 'D:/rbullins/'; % Computer at lab only
BASEPATH = ['C:/Users/' USER '/OneDrive - University of North Carolina at Chapel Hill/Hantman_Lab/Harmaline_Project/'];
CODE_REAGAN = [BASEPATH 'Code/reagan_code/'];
CODE_PROCESS_EVENTS = [RAWDATA_BASEPATH 'Code/process_events/'];
addpath(genpath(CODE_REAGAN));
addpath(genpath(CODE_PROCESS_EVENTS));

%% Run Directory with animal
Directory_Animals;
%% Experimental conditions 
exper_conditions = {'control';'harm'};
%% Score
scoreLabel = num2str(score);
if strcmp(scoreLabel, '1')
    SCORE = 'idealSuccess';
elseif strcmp(scoreLabel, '0')
    SCORE = 'noSuccess';
elseif strcmp(scoreLabel, '-1')
    SCORE = 'noReach';
elseif strcmp(scoreLabel, '2')
    SCORE = 'eventualSuccess';
elseif strcmp(scoreLabel, 'all')
    score = [1 0 2];
    SCORE = 'allScores';
end
%% Plot endpoints only of movement (o) and position of pellet (.)
% Also calculates the x,y,z position of the mean of all pellet locations,
% and saves to a matrix called upon later for aligment
for isub = 1:length(animals)
    SUB = animals{isub};
    COMPARISON_FIGS = [BASEPATH 'Figures/' SUB '/Behavior/HarmVsControl/'];
    figure;
    % Make matrix to store average pellet location for each condition
    pelletPos = zeros(length(exper_conditions),3);
    for iexper = 1:length(exper_conditions)
        EXPER_COND = exper_conditions{iexper};
        ExperSessions = eval(sprintf('%s_%sBehaviorVideos',SUB,EXPER_COND));
        % Add path for analyzed mat files
        ANALYZED_MAT = [BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SCORE '/'];
        for isession = 1:length(ExperSessions)
            EXPER_SESSION = ExperSessions{isession};
            if isempty(EXPER_SESSION)
                continue;
            end
            % Load in variables
            MAT_FILE = [ANALYZED_MAT SUB '_' EXPER_SESSION '_' EXPER_COND '_Overlay3DVariables.mat'];
            MAT_PELLET = [ANALYZED_MAT SUB '_' EXPER_SESSION '_' EXPER_COND '_PelletTraj.mat'];
            load(MAT_FILE,'movIdx','movEnd','movStrt','traj');
            load(MAT_PELLET,'pellet_traj');

            movIdx = eval(sprintf('movIdx.%s',block));
            movEnd = eval(sprintf('movEnd.%s',block));%movEnd.nbase;
            movStart = eval(sprintf('movStrt.%s',block));%movStrt.nbase;
            movTraj = traj;

            % Define color for plotting based on condition (harmaline or
            % control)
            if strcmp(EXPER_COND, 'harm')
                colorMap = [1 0 0];
                colorMapPellet = [139/255,0,0];
            else
                colorMap = [0 0 1];
                colorMapPellet = [0,0,139/255];
            end
            allPelFrames = [];
            pelGrabFrame =350;%avg_frames_turn;
            % Loop through trials
            for itrial = 1:length(movIdx)
                frameStart = movStart(1,itrial);
                frameEnd = movEnd(1,itrial);
                if frameEnd > length(movTraj(movIdx(itrial),1,:))
                    frameEnd = length(movTraj(movIdx(itrial),1,:));
                end
                % Identify trajectory of movement and pellet
                thisTraj = squeeze(movTraj(movIdx(itrial),:,frameStart:frameEnd));
                thisPelletTraj = squeeze(pellet_traj(movIdx(itrial),:,:)); %3rd was : (1/6/2022)
   
                %                 thisTraj_Smooth = smoothdata(thisTraj,'gaussian',5);
                %                 hold on;
                %                 plot3(thisTraj_Smooth(1,end),thisTraj_Smooth(2,end),thisTraj_Smooth(3,end),'o','Color',colorMap(1,:));
                %             end
                % Plot trajectory endpoint (o) and pellet endpoint (.)
                plot3(thisTraj(1,end),thisTraj(2,end),thisTraj(3,end),'o','Color',colorMap(1,:));
                hold on;
                plot3(thisPelletTraj(1,pelGrabFrame),...round(size(thisPelletTraj,2)/2))
                    thisPelletTraj(2,pelGrabFrame),...
                    thisPelletTraj(3,pelGrabFrame),'.','Color',colorMapPellet(1,:));
            end % trials
            xlabel('x')
            ylabel('y')
            hold on;
            % Find pellet coordinate, x y z (time and velocity of pellet)
            % Find majority x position
            %%% WORKING HERE - have pellet frames which correspond to
            %%% index, now need to get the x position at that index and
            %%% frame into a vector
            x_position = squeeze(pellet_traj(movIdx,1,pelGrabFrame));figure;
            % Bin data
            x_val = histogram(x_position,70);
            % Find bin with most values in it
            idx_x = find(x_val.Values == max(x_val.Values));
            % If there are more than one bin that is the max, take the
            % average of the two (round to be a whole number)
            if length(idx_x) > 1
                idx_x = round(mean(idx_x));
            end
            x_maj = x_val.BinEdges(idx_x);close;
            % Find majority y position
            y_position = squeeze(pellet_traj(movIdx,2,pelGrabFrame));figure;
            y_val = histogram(y_position,70);
            idx_y = find(y_val.Values == max(y_val.Values));
            if length(idx_y) > 1
                idx_y = round(mean(idx_y));
            end
            y_maj = y_val.BinEdges(idx_y);close;
            % Find majority z position
            z_position = squeeze(pellet_traj(movIdx,3,pelGrabFrame));figure;
            z_val = histogram(z_position,70);
            idx_z = find(z_val.Values == max(z_val.Values));
            if length(idx_z) > 1
                idx_z = round(mean(idx_z));
            end
            z_maj = z_val.BinEdges(idx_z);close;
            % Define average pellet position and plot as a black circle
            pelletPos(iexper,1) = x_maj;
            pelletPos(iexper,2) = y_maj;
            pelletPos(iexper,3) = z_maj;
            plot3(x_maj, y_maj, z_maj,'o','MarkerFaceColor','k','MarkerEdgeColor','k');
        end % session
    end % experimental conditions
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title({[SUB ': Harmaline (red) vs Control (blue) ' SCORE ' Reach'],['Block: ' block]});
    savefig([COMPARISON_FIGS SUB '_3DOverlay_' SCORE 'Trials_' block '_Raw_Endpoints.fig']);
    % misalignment between the pellets in different experimental
    % conditions,but the same animal
    add2control = pelletPos(2,:)-pelletPos(1,:);
    % Save alignment configuration for the next piece of code
    save([ANALYZED_MAT SUB '_controlVsHarm_pelletAlignment.mat'],'add2control');
end
%% Aligned Data, Plot endpoints of harmaline(red) vs control(blue) for Success (o) vs failures (x)
% (failures is anything but ideal success)
% Loop through each subject
for isub = 1:length(animals)
    SUB = animals{isub};
    COMPARISON_FIGS = [BASEPATH 'Figures/' SUB '/Behavior/HarmVsControl/'];
    figure;
    for iexper = 1:length(exper_conditions)
        EXPER_COND = exper_conditions{iexper};
        ExperSessions = eval(sprintf('%s_%sBehaviorVideos',SUB,EXPER_COND));
        % Add path for analyzed mat files
        ANALYZED_MAT = [BASEPATH 'Data_Analyzed\' SUB '\Behavior\' SCORE '\'];
        for isession = 1:length(ExperSessions)
            EXPER_SESSION = ExperSessions{isession};
            trialIdxs = eval(sprintf('%s_%s_%sTrials',SUB,EXPER_SESSION,EXPER_COND));
            if isempty(EXPER_SESSION)
                continue;
            end
            %load in variables and plot on figure
            MAT_FILE = [ANALYZED_MAT SUB '_' EXPER_SESSION '_' EXPER_COND '_Overlay3DVariables.mat'];
            MAT_PELLET = [ANALYZED_MAT SUB '_' EXPER_SESSION '_' EXPER_COND '_PelletTraj.mat'];
            load(MAT_FILE,'movIdx','movEnd','movStrt','traj');
            load(MAT_PELLET,'pellet_traj');

            movIdx = eval(sprintf('movIdx.%s',block));
            movEnd = eval(sprintf('movEnd.%s',block));%movEnd.nbase;
            movStart = eval(sprintf('movStrt.%s',block));%movStrt.nbase;
            movTraj = traj;
            if strcmp(EXPER_COND, 'harm')
                colorMap = [1 0 0];
                colorMapPellet = [139/255,0,0];
            else
                colorMap = [0 0 1];
                colorMapPellet = [0,0,139/255];
            end
            % load in pellet alignment variable
            load([ANALYZED_MAT SUB '_controlVsHarm_pelletAlignment.mat'],'add2control'); %add2control
            % For each trial, plot end point and determine if success or
            % failure

            allTrial_End_x = [];
            allTrial_End_y = [];
            allTrial_End_z = [];
            pelGrabFrame =350;
            for itrial = 1:length(movIdx)
                frameStart = movStart(1,itrial);
                frameEnd = movEnd(1,itrial);
                if frameEnd > length(movTraj(movIdx(itrial),1,:))
                    frameEnd = length(movTraj(movIdx(itrial),1,:));
                end
               thisPelletTraj = squeeze(pellet_traj(movIdx(itrial),:,:));
                thisTraj = squeeze(movTraj(movIdx(itrial),:,frameStart:frameEnd));
                % Determine if trial is success or failure, plot shape as
                % such
                % If success plot o
                if trialIdxs.trialScore(movIdx(itrial)) == 1
                    plotShape = 'o';
                    % If failure plot x
                else
                    plotShape = 'x';
                end
                % Determine if harmaline or control condition
                % If control condition add alignment factor

                if strcmp(EXPER_COND,'control')
                    plot3(thisTraj(1,end)+add2control(1),...
                        thisTraj(2,end)+add2control(2),...
                        thisTraj(3,end)+add2control(3),plotShape,'Color',colorMap(1,:));
                    hold on;
                    plot3(thisPelletTraj(1,pelGrabFrame)+add2control(1),...round(size(thisPelletTraj,2)/2))
                        thisPelletTraj(2,pelGrabFrame)+add2control(2),...
                        thisPelletTraj(3,pelGrabFrame)+add2control(3),'.','Color',colorMapPellet(1,:));
                    % Gather all ends, to get mean later on
                    allTrial_End_x = [allTrial_End_x thisTraj(1,end)+add2control(1)];
                    allTrial_End_y = [allTrial_End_y thisTraj(2,end)+add2control(2)];
                    allTrial_End_z = [allTrial_End_z thisTraj(3,end)+add2control(3)];
                    % If harmaline condition, plot as is
                else
                    plot3(thisTraj(1,end),thisTraj(2,end),thisTraj(3,end),plotShape,'Color',colorMap(1,:));
                    hold on;
                    plot3(thisPelletTraj(1,pelGrabFrame),...round(size(thisPelletTraj,2)/2))
                        thisPelletTraj(2,pelGrabFrame),...
                        thisPelletTraj(3,pelGrabFrame),'.','Color',colorMapPellet(1,:));
                    % Gather all ends, to get mean later on
                    allTrial_End_x = [allTrial_End_x thisTraj(1,end)];
                    allTrial_End_y = [allTrial_End_y thisTraj(2,end)];
                    allTrial_End_z = [allTrial_End_z thisTraj(3,end)];
                end
            end % trials
            hold on;
            % find successes and plot mean
            success_idx = find(trialIdxs.trialScore(movIdx) == 1);
            if success_idx >=1
            plot3(mean(allTrial_End_x(1,success_idx)),...
                mean(allTrial_End_y(1,success_idx)),...
                mean(allTrial_End_z(1,success_idx)),...
                'o','LineWidth',2.0,'Color',colorMapPellet(1,:),'MarkerSize',10);
            end
            % Find failures and plot mean
            failure_idx = find(trialIdxs.trialScore(movIdx) ~= 1);
            if failure_idx >=1 
            plot3(mean(allTrial_End_x(1,failure_idx)),...
                mean(allTrial_End_y(1,failure_idx)),...
                mean(allTrial_End_z(1,failure_idx)),...
                'x','LineWidth',2.0,'Color',colorMapPellet(1,:),'MarkerSize',10);
            end
        end % experimental session
    end % experimental conditions
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title({[SUB ': Harmaline (red) vs Control (blue) ' SCORE ' Reach - Aligned'],['Block: ' block]});
    savefig([COMPARISON_FIGS SUB '_3DOverlay_' SCORE 'Trials_' block '_Raw_EndpointsAligned.fig']);
end % subjects
%% PART 2 Aligned Data, Plot endpoints of harmaline(red) vs control(blue) for Success (o) vs failures (x)
% (failures is anything but ideal success)
% Loop through each subject
for isub = 1:length(animals)
    SUB = animals{isub};
    COMPARISON_FIGS = [BASEPATH 'Figures/' SUB '/Behavior/HarmVsControl/'];
    figure;
    for iexper = 1:length(exper_conditions)
        EXPER_COND = exper_conditions{iexper};
        ExperSessions = eval(sprintf('%s_%sBehaviorVideos',SUB,EXPER_COND));
        % Add path for analyzed mat files
        ANALYZED_MAT = [BASEPATH 'Data_Analyzed\' SUB '\Behavior\' SCORE '\'];
        for isession = 1:length(ExperSessions)
            EXPER_SESSION = ExperSessions{isession};
            trialIdxs = eval(sprintf('%s_%s_%sTrials',SUB,EXPER_SESSION,EXPER_COND));
            if isempty(EXPER_SESSION)
                continue;
            end
            %load in variables and plot on figure
            MAT_FILE = [ANALYZED_MAT SUB '_' EXPER_SESSION '_' EXPER_COND '_Overlay3DVariables.mat'];
            MAT_PELLET = [ANALYZED_MAT SUB '_' EXPER_SESSION '_' EXPER_COND '_PelletTraj.mat'];
            load(MAT_FILE,'movIdx','movEnd','movStrt','traj');
            load(MAT_PELLET,'pellet_traj');

            movIdx = eval(sprintf('movIdx.%s',block));
            movEnd = eval(sprintf('movEnd.%s',block));%movEnd.nbase;
            movStart = eval(sprintf('movStrt.%s',block));%movStrt.nbase;
            movTraj = traj;
            if strcmp(EXPER_COND, 'harm')
                colorMap = [1 0 0];
                colorMapPellet = [139/255,0,0];
            else
                colorMap = [0 0 1];
                colorMapPellet = [0,0,139/255];
            end
            % load in pellet alignment variable
            load([ANALYZED_MAT SUB '_controlVsHarm_pelletAlignment.mat'],'add2control'); %add2control
            % For each trial, plot end point and determine if success or
            % failure

            allTrial_End_x = [];
            allTrial_End_y = [];
            allTrial_End_z = [];
            pelGrabFrame =350;
            for itrial = 1:length(movIdx)
                frameStart = movStart(1,itrial);
                frameEnd = movEnd(1,itrial);
                if frameEnd > length(movTraj(movIdx(itrial),1,:))
                    frameEnd = length(movTraj(movIdx(itrial),1,:));
                end
               thisPelletTraj = squeeze(pellet_traj(movIdx(itrial),:,:));
                thisTraj = squeeze(movTraj(movIdx(itrial),:,frameStart:frameEnd));
                % Determine if trial is success or failure, plot shape as
                % such
                % If success plot o
                if trialIdxs.trialScore(movIdx(itrial)) == 1
                    plotShape = 'o';
                    % If failure plot x
                else
                    plotShape = 'x';
                end
                % Determine if harmaline or control condition
                % If control condition add alignment factor

                if strcmp(EXPER_COND,'control')
                    plot3(thisTraj(1,end)+add2control(1),...
                        thisTraj(2,end)+add2control(2),...
                        thisTraj(3,end)+add2control(3),plotShape,'Color',colorMap(1,:));
                    hold on;
                    plot3(thisPelletTraj(1,pelGrabFrame)+add2control(1),...round(size(thisPelletTraj,2)/2))
                        thisPelletTraj(2,pelGrabFrame)+add2control(2),...
                        thisPelletTraj(3,pelGrabFrame)+add2control(3),'.','Color',colorMapPellet(1,:));
                   plot3(thisTraj(1,:)+ add2control(1),...
                        thisTraj(2,:)+ add2control(2),...
                        thisTraj(3,:)+ add2control(3),...
                        'Color',colorMap(1,:));
                else
                    plot3(thisTraj(1,end),thisTraj(2,end),thisTraj(3,end),plotShape,'Color',colorMap(1,:));
                    hold on;
                    plot3(thisPelletTraj(1,pelGrabFrame),...round(size(thisPelletTraj,2)/2))
                        thisPelletTraj(2,pelGrabFrame),...
                        thisPelletTraj(3,pelGrabFrame),'.','Color',colorMapPellet(1,:));
                    plot3(thisTraj(1,:),thisTraj(2,:),thisTraj(3,:),'Color',colorMap(1,:));
                end
            end % trials
            hold on;
        end % experimental session
    end % experimental conditions
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title({[SUB ': Harmaline (red) vs Control (blue) ' SCORE ' Reach - Aligned'],['Block: ' block]});
    savefig([COMPARISON_FIGS SUB '_3DOverlay_' SCORE 'Trials_' block '_Raw_TrajectoriesAligned.fig']);
end % subjects


% Pellet locating code old part 
% % back up
%             % Loop through trials
%             for itrial = 1:length(movIdx)
%                 frameStart = movStart(1,itrial);
%                 frameEnd = movEnd(1,itrial);
%                 if frameEnd > length(movTraj(movIdx(itrial),1,:))
%                     frameEnd = length(movTraj(movIdx(itrial),1,:));
%                 end
%                 % Identify trajectory of movement and pellet
%                 thisTraj = squeeze(movTraj(movIdx(itrial),:,frameStart:frameEnd));
%                 thisPelletTraj = squeeze(pellet_traj(movIdx(itrial),:,:));
%                 % find where pellet has velocity of zero for more than 150
%                 % time points in x direction
% % %                 PelletTrajTemp_x = smoothdata(thisPelletTraj(1,:),'gaussian',30);
% % %                 velocityPellet_x = diff(PelletTrajTemp_x);
% % %                 velocityLow_x = find(velocityPellet_x < .01 & velocityPellet_x > -.01);
% % %                 velocityLowAndOvertime_x = velocityLow_x(find(velocityLow_x > 350));
% % %                 %find where pellet has velocity of zero for more than 150
% % %                 %in y direction
% % %                 PelletTrajTemp_y = smoothdata(thisPelletTraj(2,:),'gaussian',30);
% % %                 velocityPellet_y = diff(PelletTrajTemp_y);
% % %                 velocityLow_y = find(velocityPellet_y < .01 & velocityPellet_y > -.01);
% % %                 velocityLowAndOvertime_y = velocityLow_y(find(velocityLow_y > 350));
% % %                 % find positions where x and y are both zero
% % %                 if velocityLowAndOvertime_x(30) <= velocityLowAndOvertime_x(1) + 150
% % %                     if velocityLowAndOvertime_y(30) <= velocityLowAndOvertime_y(1) + 150
% % %                         idx_pull = velocityLowAndOvertime_x(ismember(velocityLowAndOvertime_x,velocityLowAndOvertime_y));
% % %                         pelGrabFrame = idx_pull(1);
% % %                         allPelFrames = [allPelFrames pelGrabFrame];
% % %                     end
% % %                 else
% % %                     disp('Work pellet finding code');
% % % 
% % %                 end
% 
%                 % find only velocities above 350 time points
%                 %                 thisTraj_Smooth = smoothdata(thisTraj,'gaussian',5);
%                 %                 hold on;
%                 %                 plot3(thisTraj_Smooth(1,end),thisTraj_Smooth(2,end),thisTraj_Smooth(3,end),'o','Color',colorMap(1,:));
%                 %             end
%                 % Plot trajectory endpoint (o) and pellet endpoint (.)
%                 plot3(thisTraj(1,end),thisTraj(2,end),thisTraj(3,end),'o','Color',colorMap(1,:));
%                 hold on;
%                 plot3(thisPelletTraj(1,pelGrabFrame),...round(size(thisPelletTraj,2)/2))
%                     thisPelletTraj(2,pelGrabFrame),...
%                     thisPelletTraj(3,pelGrabFrame),'.','Color',colorMapPellet(1,:));
%             end % trials
%             xlabel('x')
%             ylabel('y')
%             hold on;
%             % Find pellet coordinate, x y z (time and velocity of pellet)
%             % Find majority x position
%             %%% WORKING HERE - have pellet frames which correspond to
%             %%% index, now need to get the x position at that index and
%             %%% frame into a vector
%             x_position = squeeze(pellet_traj(movIdx,1,allPelFrames));figure;
%             % Bin data
%             x_val = histogram(x_position,70);
%             % Find bin with most values in it
%             idx_x = find(x_val.Values == max(x_val.Values));
%             % If there are more than one bin that is the max, take the
%             % average of the two (round to be a whole number)
%             if idx_x > 1
%                 idx_x = round(mean(idx_x));
%             end
%             x_maj = x_val.BinEdges(idx_x);close;
%             % Find majority y position
%             y_position = squeeze(pellet_traj(movIdx,2,frameStart:frameEnd));figure;
%             y_val = histogram(y_position,70);
%             idx_y = find(y_val.Values == max(y_val.Values));
%             if idx_y > 1
%                 idx_y = round(mean(idx_y));
%             end
%             y_maj = y_val.BinEdges(idx_y);close;
%             % Find majority z position
%             z_position = squeeze(pellet_traj(movIdx,3,frameStart:frameEnd));figure;
%             z_val = histogram(z_position,70);
%             idx_z = find(z_val.Values == max(z_val.Values));
%             if idx_z > 1
%                 idx_z = round(mean(idx_z));
%             end
%             z_maj = z_val.BinEdges(idx_z);close;
%             % Define average pellet position and plot as a black circle
%             pelletPos(iexper,1) = x_maj;
%             pelletPos(iexper,2) = y_maj;
%             pelletPos(iexper,3) = z_maj;
%             plot3(x_maj, y_maj, z_maj,'o','MarkerFaceColor','k','MarkerEdgeColor','k');
%         end % session
%     end % experimental conditions
%     xlabel('X');
%     ylabel('Y');
%     zlabel('Z');
%     title([SUB ': Harmaline (red) vs Control (blue) ' SCORE ' Reach']);
%     savefig([COMPARISON_FIGS SUB '_3DOverlay_' SCORE 'Trials_Raw_Endpoints.fig']);
%     % misalignment between the pellets in different experimental
%     % conditions,but the same animal
%     add2control = pelletPos(2,:)-pelletPos(1,:);
%     % Save alignment configuration for the next piece of code

% %%%%%OTHER OPTION
% f_in = [RAWDATA_BASEPATH 'Data/' SUB 'necab1_Chr2/' EXPER_SESSION 'ephys/' SUB '_' EXPER_SESSION '_1000_g0/' SUB '_' EXPER_SESSION '_1000_g0_t0.nidq.bin'];
%             f_out = [BASEPATH 'Data_Analyzed/' SUB '/' SUB '_' EXPER_SESSION '_tableTurnData.mat'];
%             if exist(f_out) == 2
%                 load(f_out);
%                 tableTurnData = tableTurnData(3,:);
%             else
%                 foo = strfind(f_in,'/');
%                 dir_out = f_in(1:(foo(end)-1));
%                 thresh = [.001 .002 .001 .002 .002];
%                 [~, tableTurnData] = get_event_ind_v2(f_in,5,[1:5],thresh,0);
%                 save([BASEPATH 'Data_Analyzed/' SUB '/' SUB '_' EXPER_SESSION '_tableTurnData.mat'],'tableTurnData','-v7.3');
%                 tableTurnData = tableTurnData(3,:);
%             end
%             % Find when the table stops turning
%             diff_vector = diff(tableTurnData);
%             turnTableStop_Idx = find(diff_vector < -.003);
%             turnTableStart_Idx = find(diff_vector > .003);
%             avg_ms_turn = mean((turnTableStop_Idx - turnTableStart_Idx)./16.949);
%             avg_frames_turn = round(avg_ms_turn/2);