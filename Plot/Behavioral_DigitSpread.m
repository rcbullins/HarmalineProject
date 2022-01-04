% Behavioral Digit Spread


%  numBL.idealSuccess = length(find(trialIdxs.trialScore(base.trialIdxs) == 1));
%% Load Directory
USER = 'bullinsr';
RAWDATA_BASEPATH = 'D:/rbullins/'; % Computer at lab only
BASEPATH = ['C:/Users/' USER '/OneDrive - University of North Carolina at Chapel Hill/Hantman_Lab/Harmaline_Project/'];

Directory_Animals;
%% Compare digit spread harmaline vs control (baseline vs stim trials)
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
    digitSpreadMat = {};
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
            % Load grab distance for digits (5,1) for success, eventual
            % succes, etc
            load([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_DigitSpread.mat'],'digitSpread','movIdx','trialIdxs');
            % make mat with accuracy too
            %baseline
            idxBL.idealSuccess    = find(trialIdxs.trialScore(movIdx.nbase) == 1);
            idxBL.eventualSuccess = find(trialIdxs.trialScore(movIdx.nbase) == 2);
            idxBL.noSuccess       = find(trialIdxs.trialScore(movIdx.nbase) == 0);
            idxBL.noReach         = find(trialIdxs.trialScore(movIdx.nbase) == -1);
            idxBL.grooming        = find(trialIdxs.trialScore(movIdx.nbase) == 'g');
            % stim
            idxPert.idealSuccess    = find(trialIdxs.trialScore(movIdx.npert) == 1);
            idxPert.eventualSuccess = find(trialIdxs.trialScore(movIdx.npert) == 2);
            idxPert.noSuccess       = find(trialIdxs.trialScore(movIdx.npert) == 0);
            idxPert.noReach         = find(trialIdxs.trialScore(movIdx.npert) == -1);
            idxPert.grooming        = find(trialIdxs.trialScore(movIdx.npert) == 'g');
            % washout
            idxWash.idealSuccess    = find(trialIdxs.trialScore(movIdx.nwash) == 1);
            idxWash.eventualSuccess = find(trialIdxs.trialScore(movIdx.nwash) == 2);
            idxWash.noSuccess       = find(trialIdxs.trialScore(movIdx.nwash) == 0);
            idxWash.noReach         = find(trialIdxs.trialScore(movIdx.nwash) == -1);
            idxWash.grooming        = find(trialIdxs.trialScore(movIdx.nwash) == 'g');
            % find digit spread at this points
            % baseline
            digSpread_BL.idealSuccess = digitSpread.nbase(movIdx.nbase(ismember(movIdx.nbase, idxBL.idealSuccess)));
            digSpread_BL.eventualSuccess = digitSpread.nbase(movIdx.nbase(ismember(movIdx.nbase, idxBL.eventualSuccess)));
            digSpread_BL.noSuccess = digitSpread.nbase(movIdx.nbase(ismember(movIdx.nbase, idxBL.noSuccess)));
            digSpread_BL.noReach = digitSpread.nbase(movIdx.nbase(ismember(movIdx.nbase, idxBL.noReach)));
            digSpread_BL.grooming = digitSpread.nbase(movIdx.nbase(ismember(movIdx.nbase, idxBL.grooming)));
            % stim
            digSpread_Pert.idealSuccess = digitSpread.npert(movIdx.npert(ismember(movIdx.npert, idxPert.idealSuccess)));
            digSpread_Pert.eventualSuccess = digitSpread.npert(movIdx.npert(ismember(movIdx.npert, idxPert.eventualSuccess)));
            digSpread_Pert.noSuccess = digitSpread.npert(movIdx.npert(ismember(movIdx.npert, idxPert.noSuccess)));
            digSpread_Pert.noReach = digitSpread.npert(movIdx.npert(ismember(movIdx.npert, idxPert.noReach)));
            digSpread_Pert.grooming = digitSpread.npert(movIdx.npert(ismember(movIdx.npert, idxPert.grooming)));
            % washout
            digSpread_Wash.idealSuccess = digitSpread.nwash(movIdx.nwash(ismember(movIdx.nwash, idxWash.idealSuccess)));
            digSpread_Wash.eventualSuccess = digitSpread.nwash(movIdx.nwash(ismember(movIdx.nwash, idxWash.eventualSuccess)));
            digSpread_Wash.noSuccess = digitSpread.nwash(movIdx.nwash(ismember(movIdx.nwash, idxWash.noSuccess)));
            digSpread_Wash.noReach = digitSpread.nwash(movIdx.nwash(ismember(movIdx.nwash, idxWash.noReach)));
            digSpread_Wash.grooming = digitSpread.nwash(movIdx.nwash(ismember(movIdx.nwash, idxWash.grooming)));
            % Digit distance mean 
            digitSpreadSuccess_Mean(1+(2*i),:) = [mean(digSpread_BL.idealSuccess); mean(digSpread_BL.eventualSuccess); ...
                mean(digSpread_BL.noSuccess); mean(digSpread_BL.noReach); mean(digSpread_BL.grooming)];
            digitSpreadSuccess_Mean(2+(2*i),:) = [mean(digSpread_Pert.idealSuccess); mean(digSpread_Pert.eventualSuccess); ...
                mean(digSpread_Pert.noSuccess); mean(digSpread_Pert.noReach); mean(digSpread_Pert.grooming)];
            % Digit distance sem
            digitSpreadSuccess_sem(1+(2*i),:) = [std(digSpread_BL.idealSuccess)/sqrt(length(digSpread_BL.idealSuccess));...
                                                 std(digSpread_BL.eventualSuccess)/sqrt(length(digSpread_BL.eventualSuccess)); ...
                                                 std(digSpread_BL.noSuccess)/sqrt(length(digSpread_BL.noSuccess));...
                                                 std(digSpread_BL.noReach)/sqrt(length(digSpread_BL.noReach)); ...
                                                 std(digSpread_BL.grooming)/sqrt(length(digSpread_BL.grooming))];
            digitSpreadSuccess_sem(2+(2*i),:) = [std(digSpread_Pert.idealSuccess)/sqrt(length(digSpread_Pert.idealSuccess));...
                                                 std(digSpread_Pert.eventualSuccess)/sqrt(length(digSpread_Pert.eventualSuccess)); ...
                                                 std(digSpread_Pert.noSuccess)/sqrt(length(digSpread_Pert.noSuccess));...
                                                 std(digSpread_Pert.noReach)/sqrt(length(digSpread_Pert.noReach)); ...
                                                 std(digSpread_Pert.grooming)/sqrt(length(digSpread_Pert.grooming))];
            % Digit Scatter Mat
            digitSpreadSuccess_Scatter(1+(10*i),:) = digSpread_BL.idealSuccess;
            digitSpreadSuccess_Scatter(2+(10*i),:) = digSpread_BL.eventualSuccess;
            digitSpreadSuccess_Scatter(3+(10*i),:) = digSpread_BL.noSuccess;
            digitSpreadSuccess_Scatter(4+(10*i),:) = digSpread_BL.noReach;
            digitSpreadSuccess_Scatter(5+(10*i),:) = digSpread_BL.grooming;
            digitSpreadSuccess_Scatter(6+(10*i),:) = digSpread_Pert.idealSuccess;
            digitSpreadSuccess_Scatter(7+(10*i),:) = digSpread_Pert.eventualSuccess;
            digitSpreadSuccess_Scatter(8+(10*i),:) = digSpread_Pert.noSuccess;
            digitSpreadSuccess_Scatter(9+(10*i),:) = digSpread_Pert.noReach;
            digitSpreadSuccess_Scatter(10+(10*i),:) = digSpread_Pert.grooming;
            % Add digit distance to struct
            digitSpreadMat{1+(2*i)} = digitSpread.nbase;
            digitSpreadMat{2+(2*i)} = digitSpread.npert;
            % Update index of experimental session
            i = i +1;
        end % experiment session
    end % experiment condition
    digitSpreadSuccess_Mean(isnan(digitSpreadSuccess_Mean)) = 0;
    %% Plot Bar chart digit spread comparing harmaline vs control (baseline vs stim)
    % digit spread mean
    digSpread_mean(1,1) = mean(digitSpreadMat{1});
    digSpread_mean(1,2) = mean(digitSpreadMat{2});
    digSpread_mean(1,3) = mean(digitSpreadMat{3});
    digSpread_mean(1,4) = mean(digitSpreadMat{4});
    % digit spread sem
    digSpread_sem(1,1) = std(digitSpreadMat{1}) / sqrt(length(digitSpreadMat{1}));
    digSpread_sem(1,2) = std(digitSpreadMat{2}) / sqrt(length(digitSpreadMat{2}));
    digSpread_sem(1,3) = std(digitSpreadMat{3}) / sqrt(length(digitSpreadMat{3}));
    digSpread_sem(1,4) = std(digitSpreadMat{4}) / sqrt(length(digitSpreadMat{4}));
    % make the figure
    figure;
    h = bar(digSpread_mean,'facealpha',.4, 'edgecolor', 'none', 'FaceColor','flat');
    hold on;
    h.CData(1,:) = [0 .5 1]; %blue (baseline)
    h.CData(2,:) = [1 0 0];  %red (stimulation)
    h.CData(3,:) = [0 .5 1]; %blue
    h.CData(4,:) = [1 0 0];  %red
    % Add data point for each participant
    scatter(1*ones(1,length(digitSpreadMat{1})), digitSpreadMat{1}, 'jitter', 'on', 'jitterAmount', 0.3, 'MarkerFaceColor', [0 .5 1], 'MarkerEdgeColor', 'none');
    scatter(2*ones(1,length(digitSpreadMat{2})), digitSpreadMat{2}, 'jitter', 'on', 'jitterAmount', 0.3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
    scatter(3*ones(1,length(digitSpreadMat{3})), digitSpreadMat{3}, 'jitter', 'on', 'jitterAmount', 0.3, 'MarkerFaceColor', [0 .5 1], 'MarkerEdgeColor', 'none');
    scatter(4*ones(1,length(digitSpreadMat{4})), digitSpreadMat{4}, 'jitter', 'on', 'jitterAmount', 0.3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
    % finish labeling
    xticklabels({'Control Baseline';'Control Stim';'Harmaline Baseline';'Harmaline Stim'});
    title([ SUB ': Digit Spread at Grab']);
    x =[];
    errorbar(x',digSpread_mean,digSpread_sem,'k','LineWidth', 1.5,'linestyle','none','HandleVisibility','off');
    ylabel('Digit Spread (pixels)');
    %% Figure two: fig one in more detail - with success level too
    % get breakdown of accuracy
    figure;
    h = bar(digitSpreadSuccess_Mean);
    hold on;
    title([ SUB ': Digit Spread at Grab']);
    xticklabels({'Control Baseline';'Control Stim';'Harmaline Baseline';'Harmaline Stim'});
    l = cell(1,5);
    l{1}='Success'; l{2}='Eventual Success'; l{3}='No Success'; l{4}='No Reach'; l{5}='Grooming';
    legend(h,l); 
    errorbar(h(1).XEndPoints,digitSpreadSuccess_Mean,digitSpreadSuccess_sem,'k','LineWidth', 1.5,'linestyle','none','HandleVisibility','off');
    % add scatter to bars
    % add control baseline
    scatter(h(1).XEndPoints(1)*ones(length(digitSpreadSuccess_Scatter(1,:))), digitSpreadSuccess_Scatter(1,:),'jitter', 'on', 'jitterAmount', 0.3, 'MarkerFaceColor',h(1).FaceColor, 'MarkerEdgeColor', 'none');
    scatter(h(2).XEndPoints(1)*ones(length(digitSpreadSuccess_Scatter(2,:))), digitSpreadSuccess_Scatter(2,:),'jitter', 'on', 'jitterAmount', 0.3, 'MarkerFaceColor',h(2).FaceColor, 'MarkerEdgeColor', 'none');
    scatter(h(3).XEndPoints(1)*ones(length(digitSpreadSuccess_Scatter(3,:))), digitSpreadSuccess_Scatter(3,:),'jitter', 'on', 'jitterAmount', 0.3,'MarkerFaceColor', h(3).FaceColor, 'MarkerEdgeColor', 'none');
    scatter(h(4).XEndPoints(1)*ones(length(digitSpreadSuccess_Scatter(4,:))), digitSpreadSuccess_Scatter(4,:),'jitter', 'on', 'jitterAmount', 0.3,'MarkerFaceColor', h(4).FaceColor, 'MarkerEdgeColor', 'none');
    scatter(h(5).XEndPoints(1)*ones(length(digitSpreadSuccess_Scatter(5,:))), digitSpreadSuccess_Scatter(5,:),'jitter', 'on', 'jitterAmount', 0.3,'MarkerFaceColor', h(5).FaceColor, 'MarkerEdgeColor', 'none');
    % control stim
    scatter(h(1).XEndPoints(2)*ones(length(digitSpreadSuccess_Scatter(6,:))), digitSpreadSuccess_Scatter(6,:),'jitter', 'on', 'jitterAmount', 0.3,'MarkerFaceColor', h(1).FaceColor, 'MarkerEdgeColor', 'none');
    scatter(h(2).XEndPoints(2)*ones(length(digitSpreadSuccess_Scatter(7,:))), digitSpreadSuccess_Scatter(7,:),'jitter', 'on', 'jitterAmount', 0.3,'MarkerFaceColor', h(2).FaceColor, 'MarkerEdgeColor', 'none');
    scatter(h(3).XEndPoints(2)*ones(length(digitSpreadSuccess_Scatter(8,:))), digitSpreadSuccess_Scatter(8,:),'jitter', 'on', 'jitterAmount', 0.3,'MarkerFaceColor', h(3).FaceColor, 'MarkerEdgeColor', 'none');
    scatter(h(4).XEndPoints(2)*ones(length(digitSpreadSuccess_Scatter(9,:))), digitSpreadSuccess_Scatter(9,:),'jitter', 'on', 'jitterAmount', 0.3,'MarkerFaceColor', h(4).FaceColor, 'MarkerEdgeColor', 'none');
    scatter(h(5).XEndPoints(2)*ones(length(digitSpreadSuccess_Scatter(10,:))), digitSpreadSuccess_Scatter(10,:),'jitter', 'on', 'jitterAmount', 0.3,'MarkerFaceColor', h(5).FaceColor, 'MarkerEdgeColor', 'none');
    % harm baseline
    scatter(h(1).XEndPoints(3)*ones(length(digitSpreadSuccess_Scatter(11,:))), digitSpreadSuccess_Scatter(11,:),'jitter', 'on', 'jitterAmount', 0.3,'MarkerFaceColor', h(1).FaceColor, 'MarkerEdgeColor', 'none');
    scatter(h(2).XEndPoints(3)*ones(length(digitSpreadSuccess_Scatter(12,:))), digitSpreadSuccess_Scatter(12,:),'jitter', 'on', 'jitterAmount', 0.3,'MarkerFaceColor', h(2).FaceColor, 'MarkerEdgeColor', 'none');
    scatter(h(3).XEndPoints(3)*ones(length(digitSpreadSuccess_Scatter(13,:))), digitSpreadSuccess_Scatter(13,:),'jitter', 'on', 'jitterAmount', 0.3,'MarkerFaceColor', h(3).FaceColor, 'MarkerEdgeColor', 'none');
    scatter(h(4).XEndPoints(3)*ones(length(digitSpreadSuccess_Scatter(14,:))), digitSpreadSuccess_Scatter(14,:),'jitter', 'on', 'jitterAmount', 0.3,'MarkerFaceColor', h(4).FaceColor, 'MarkerEdgeColor', 'none');
    scatter(h(5).XEndPoints(3)*ones(length(digitSpreadSuccess_Scatter(15,:))), digitSpreadSuccess_Scatter(15,:),'jitter', 'on', 'jitterAmount', 0.3,'MarkerFaceColor', h(5).FaceColor, 'MarkerEdgeColor', 'none');
    % harm stim
    scatter(h(1).XEndPoints(4)*ones(length(digitSpreadSuccess_Scatter(16,:))), digitSpreadSuccess_Scatter(16,:),'jitter', 'on', 'jitterAmount', 0.3,'MarkerFaceColor', h(1).FaceColor, 'MarkerEdgeColor', 'none');
    scatter(h(2).XEndPoints(4)*ones(length(digitSpreadSuccess_Scatter(17,:))), digitSpreadSuccess_Scatter(17,:),'jitter', 'on', 'jitterAmount', 0.3,'MarkerFaceColor', h(2).FaceColor, 'MarkerEdgeColor', 'none');
    scatter(h(3).XEndPoints(4)*ones(length(digitSpreadSuccess_Scatter(18,:))), digitSpreadSuccess_Scatter(18,:),'jitter', 'on', 'jitterAmount', 0.3,'MarkerFaceColor', h(3).FaceColor, 'MarkerEdgeColor', 'none');
    scatter(h(4).XEndPoints(4)*ones(length(digitSpreadSuccess_Scatter(19,:))), digitSpreadSuccess_Scatter(19,:),'jitter', 'on', 'jitterAmount', 0.3,'MarkerFaceColor', h(4).FaceColor, 'MarkerEdgeColor', 'none');
    scatter(h(5).XEndPoints(4)*ones(length(digitSpreadSuccess_Scatter(20,:))), digitSpreadSuccess_Scatter(20,:),'jitter', 'on', 'jitterAmount', 0.3,'MarkerFaceColor', h(5).FaceColor, 'MarkerEdgeColor', 'none');

end % subject