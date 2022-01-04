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
            
            % Add digit distance to struct
            digitSpreadMat{1+(2*i)} = digitSpread.nbase;
            digitSpreadMat{2+(2*i)} = digitSpread.npert;
            % Update index of experimental session
            i = i +1;
        end % experiment session
    end % experiment condition

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
    h.CData(1,:) = [0 .5 1]; %blue (control)
    h.CData(2,:) = [1 0 0];  %red (harm)
    h.CData(3,:) = [0 .5 1]; %blue
    h.CData(4,:) = [1 0 0];  %red
    % Add data point for each participant
    scatter(1*ones(1,length(digitSpreadMat{1})), digitSpreadMat{1}, 'jitter', 'on', 'jitterAmount', 0.3, 'MarkerFaceColor', [0 .5 1], 'MarkerEdgeColor', 'none');
    scatter(2*ones(1,length(digitSpreadMat{2})), digitSpreadMat{2}, 'jitter', 'on', 'jitterAmount', 0.3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
    scatter(3*ones(1,length(digitSpreadMat{3})), digitSpreadMat{3}, 'jitter', 'on', 'jitterAmount', 0.3, 'MarkerFaceColor', [0 .5 1], 'MarkerEdgeColor', 'none');
    scatter(4*ones(1,length(digitSpreadMat{4})), digitSpreadMat{4}, 'jitter', 'on', 'jitterAmount', 0.3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
   % Add error bars
    for plotIdx = 1:4
        plot([plotIdx plotIdx],[digSpread_mean(1,plotIdx) - digSpread_sem(1,plotIdx) ...
            digSpread_mean(1,plotIdx) + digSpread_sem(1,plotIdx)],'k');
    end
    % finish labeling 
    xticklabels({'Control Baseline';'Control Stim';'Harmaline Baseline';'Harmaline Stim'});
    title([ SUB ': Digit Spread']);
    ylabel('Digit Spread (pixels)');
    %% Figure two: break up by success level


    %     xticklabels({'Control Baseline';'Control Stim';'Harmaline Baseline';'Harmaline Stim'});
    %     l = cell(1,5);
    %     l{1}='Success'; l{2}='Eventual Success'; l{3}='No Success'; l{4}='No Reach'; l{5}='Grooming';
    %     legend(h,l);


end % subject