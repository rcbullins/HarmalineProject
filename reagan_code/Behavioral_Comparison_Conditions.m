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

        end % session
    end % experimental conditions
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title([SUB ': Harmaline (red) vs Control (blue) ' SCORE ' Reach']);
    savefig([COMPARISON_FIGS SUB '_3DOverlay_' SCORE 'Trials_Raw_Endpoints.fig']);

    %% Find corrected for pellet difference (diff sessions have slightly different x and y)
   % [n x] = hist(squeeze(pellet_traj(movIdx),:,frameStart:frameEnd),100)

end