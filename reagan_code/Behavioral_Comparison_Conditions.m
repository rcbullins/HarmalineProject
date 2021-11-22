%Behavioral_Comparison_Conditions
%%%%% PART 2: COMPARISONS OF CONTROL AND HAMRALINE %%%%%%%
% 3D Overlays of harmaline vs controls (success trials only)
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
            if isempty(EXPER_SESSION)
                continue;
            end
            % Overlay 3D Figure
            %load in variables and plot on figure
            MAT_FILE = [ANALYZED_MAT SUB '_' EXPER_SESSION '_' EXPER_COND '_Overlay3DVariables.mat'];
            MAT_PELLET = [ANALYZED_MAT SUB '_' EXPER_SESSION '_' EXPER_COND '_PelletTraj.mat'];
            load(MAT_FILE);
            load(MAT_PELLET);
            movIdx = movIdx.nbase;
            movEnd = movEnd.nbase;
            movStart = movStrt.nbase;
            movTraj = traj;
            if strcmp(EXPER_COND, 'harm')
                colorMap = [1 0 0];
                colorMapPellet = [139/255,0,0];
            else
                colorMap = [0 0 1];
                colorMapPellet = [0 0 139/255];
            end

            for itrial = 1:length(movIdx)
                frameStart = movStart(1,itrial);

                frameEnd = movEnd(1,itrial);

                if frameEnd > length(movTraj(movIdx(itrial),1,:))
                    frameEnd = length(movTraj(movIdx(itrial),1,:));
                end

                thisTraj = squeeze(movTraj(movIdx(itrial),:,frameStart:frameEnd));
                thisPelletTraj = squeeze(pellet_traj(movIdx(itrial),:,frameStart:frameEnd));

                %thisTraj_Smooth = smoothdata(thisTraj,'gaussian',5);
                %                 plot3(thisTraj_Smooth(1,:),thisTraj_Smooth(2,:),thisTraj_Smooth(3,:),'Color',colorMap(1,:));
                %                 hold on;
                %                 plot3(thisTraj_Smooth(1,end),thisTraj_Smooth(2,end),thisTraj_Smooth(3,end),'o','Color',colorMap(1,:));
                %             end

                plot3(thisTraj(1,:),thisTraj(2,:),thisTraj(3,:),'Color',colorMap(1,:));
                hold on;
                plot3(thisTraj(1,end),thisTraj(2,end),thisTraj(3,end),'o','Color',colorMap(1,:));
                plot3(thisPelletTraj(1,end),...round(size(thisPelletTraj,2)/2))
                    thisPelletTraj(2,end),...
                    thisPelletTraj(3,end),'o','Color',colorMapPellet(1,:));
            end
            hold on;

        end % session
    end % experimental conditions
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title([SUB ': Harmaline (red) vs Control (blue) ' SCORE ' Reach']);
    savefig([COMPARISON_FIGS SUB '_3DOverlay_' SCORE 'Trials_Raw.fig']);
end %animal
%% Same thing with just end points
for isub = 1:length(animals)
    SUB = animals{isub};
    COMPARISON_FIGS = [BASEPATH 'Figures/' SUB '/Behavior/HarmVsControl/'];
    figure;
    pelletPos = zeros(length(exper_conditions),3);
    for iexper = 1:length(exper_conditions)
        EXPER_COND = exper_conditions{iexper};
        ExperSessions = eval(sprintf('%s_%sBehaviorVideos',SUB,EXPER_COND));
        % Make mat file to store average pellet location in per condition
            % yeah so if there is more than one session for each condition,
            % you might need to change this part :)
            
        % Add path for analyzed mat files
        ANALYZED_MAT = [BASEPATH 'Data_Analyzed\' SUB '\Behavior\' SCORE '\'];
        for isession = 1:length(ExperSessions)
            EXPER_SESSION = ExperSessions{isession};
            if isempty(EXPER_SESSION)
                continue;
            end
            % Overlay 3D Figure
            %load in variables and plot on figure
            MAT_FILE = [ANALYZED_MAT SUB '_' EXPER_SESSION '_' EXPER_COND '_Overlay3DVariables.mat'];
            load(MAT_FILE);
            clear thisPellet
            MAT_PELLET = [ANALYZED_MAT SUB '_' EXPER_SESSION '_' EXPER_COND '_PelletTraj.mat'];
            load(MAT_PELLET);

            movIdx = movIdx.nbase;
            movEnd = movEnd.nbase;
            movStart = movStrt.nbase;
            movTraj = traj;
            if strcmp(EXPER_COND, 'harm')
                colorMap = [1 0 0];
                colorMapPellet = [139/255,0,0];
            else
                colorMap = [0 0 1];
                colorMapPellet = [0,0,139/255];
            end

            for itrial = 1:length(movIdx)
                frameStart = movStart(1,itrial);

                frameEnd = movEnd(1,itrial);

                if frameEnd > length(movTraj(movIdx(itrial),1,:))
                    frameEnd = length(movTraj(movIdx(itrial),1,:));
                end

                thisTraj = squeeze(movTraj(movIdx(itrial),:,frameStart:frameEnd));

                thisPelletTraj = squeeze(pellet_traj(movIdx(itrial),:,frameStart:frameEnd));
                %                 thisTraj_Smooth = smoothdata(thisTraj,'gaussian',5);
                %                 hold on;
                %                 plot3(thisTraj_Smooth(1,end),thisTraj_Smooth(2,end),thisTraj_Smooth(3,end),'o','Color',colorMap(1,:));
                %             end
                
                plot3(thisTraj(1,end),thisTraj(2,end),thisTraj(3,end),'o','Color',colorMap(1,:));
                hold on;
                plot3(thisPelletTraj(1,end),...round(size(thisPelletTraj,2)/2))
                    thisPelletTraj(2,end),...
                    thisPelletTraj(3,end),'o','Color',colorMapPellet(1,:));
            end
            hold on;
            % Find pellet coordinate, x y z (majority of pellet position
            % found)
            x_position = squeeze(pellet_traj(movIdx,1,frameStart:frameEnd));figure;
            x_val = histogram(x_position,70);
            idx_x = find(x_val.Values == max(x_val.Values));
            if idx_x > 1
                idx_x = round(mean(idx_x));
            end
            x_maj = x_val.BinEdges(idx_x);close;
            y_position = squeeze(pellet_traj(movIdx,2,frameStart:frameEnd));figure;
            y_val = histogram(y_position,70);
            idx_y = find(y_val.Values == max(y_val.Values));
            if idx_y > 1
                idx_y = round(mean(idx_y));
            end
            y_maj = y_val.BinEdges(idx_y);close;
            z_position = squeeze(pellet_traj(movIdx,3,frameStart:frameEnd));figure;
            z_val = histogram(z_position,70);
            idx_z = find(z_val.Values == max(z_val.Values));
            if idx_z > 1
                idx_z = round(mean(idx_z));
            end
            z_maj = z_val.BinEdges(idx_z);close;
            pelletPos(iexper,1) = x_maj;
            pelletPos(iexper,2) = y_maj;
            pelletPos(iexper,3) = z_maj;
            plot3(x_maj, y_maj, z_maj,'o','MarkerFaceColor','k','MarkerEdgeColor','k');
            
        end % session
    end % experimental conditions
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title([SUB ': Harmaline (red) vs Control (blue) ' SCORE ' Reach']);
    savefig([COMPARISON_FIGS SUB '_3DOverlay_' SCORE 'Trials_Raw_Endpoints.fig']);

    % misalignment between the pellets in different experimental
    % conditions,but the same animal
    add2control = pelletPos(2,:)-pelletPos(1,:);
    save([ANALYZED_MAT SUB '_controlVsHarm_pelletAlignment.mat'],'add2control');
end



%% Same thing with but with data aligned to pellets (endpoints only)
% hamraline in red, control in blue
% success in o, failures in x (failures is anything but ideal success)
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
            % Overlay 3D Figure
            %load in variables and plot on figure
            MAT_FILE = [ANALYZED_MAT SUB '_' EXPER_SESSION '_' EXPER_COND '_Overlay3DVariables.mat'];
            load(MAT_FILE);
            clear thisPelle
            MAT_PELLET = [ANALYZED_MAT SUB '_' EXPER_SESSION '_' EXPER_COND '_PelletTraj.mat'];
            load(MAT_PELLET);

            movIdx = movIdx.nbase;
            movEnd = movEnd.nbase;
            movStart = movStrt.nbase;
            movTraj = traj;
            if strcmp(EXPER_COND, 'harm')
                colorMap = [1 0 0];
                colorMapPellet = [139/255,0,0];
            else
                colorMap = [0 0 1];
                colorMapPellet = [0,0,139/255];
            end
            % load in pellet alignment variable
            load([ANALYZED_MAT SUB '_controlVsHarm_pelletAlignment.mat']); %add2control
            % Plot
            for itrial = 1:length(movIdx)
                frameStart = movStart(1,itrial);
                frameEnd = movEnd(1,itrial);
                if frameEnd > length(movTraj(movIdx(itrial),1,:))
                    frameEnd = length(movTraj(movIdx(itrial),1,:));
                end
                thisPelletTraj = squeeze(pellet_traj(movIdx(itrial),:,frameStart:frameEnd));            
                thisTraj = squeeze(movTraj(movIdx(itrial),:,frameStart:frameEnd));
                % if success it will be a 0
                if trialIdxs.trialScore(movIdx(itrial)) == 1
                plotShape = 'o';
                else
                    plotShape = 'x';
                end
                % if control condition add alignment factor
                if strcmp(EXPER_COND,'control')
                % with correction add2control
                plot3(thisTraj(1,end)+add2control(1),...
                    thisTraj(2,end)+add2control(1),...
                    thisTraj(3,end)+add2control(3),plotShape,'Color',colorMap(1,:));
                hold on;
                plot3(thisPelletTraj(1,end)+add2control(1),...round(size(thisPelletTraj,2)/2))
                    thisPelletTraj(2,end)+add2control(2),...
                    thisPelletTraj(3,end)+add2control(3),plotShape,'Color',colorMapPellet(1,:));
                else
                % with correction add2control
                plot3(thisTraj(1,end),thisTraj(2,end),thisTraj(3,end),plotShape,'Color',colorMap(1,:));
                hold on;
                plot3(thisPelletTraj(1,end),...round(size(thisPelletTraj,2)/2))
                    thisPelletTraj(2,end),...
                    thisPelletTraj(3,end),plotShape,'Color',colorMapPellet(1,:));
                end
            end
            hold on;
            
        end % session
    end % experimental conditions
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title([SUB ': Harmaline (red) vs Control (blue) ' SCORE ' Reach - Aligned']);
    savefig([COMPARISON_FIGS SUB '_3DOverlay_' SCORE 'Trials_Raw_EndpointsAligned.fig']);
end
