function [movEndFrames, score_total_idx] = getMovEndFrame(trialScoreThr,traj, trialScores,trialBlock, movStartFrames,varargin)
% PURPOSE
%   get start frame indexes of move time intiation.
% INPUTS
%   trialScoreThr (which trials you want to look at)
%       1 = success and one reach (ideal scenario)
%       2 = multiple reaches and success
%       0 = reach and grab without success
%       -1 = nothing
%       g = grooming
%   trialScores
%       vecotr of the trial score for each trial
%   trialBlock
%       trialIdxs.nbase (baseline trials only)
%       trialIdxs.npert (stimulation trials only)
%       trialIdxs.nwash (post stim trials only)
%   traj
%       trajectory output from get_traj_3D_v2
%   movStartFrames (when movement frames start)
%   plotSanity (0 or 1) default = 1
% OUTPUS
%   movEndFrames
%       vector of all end frames for all trials when reach movement stops
%   block_score_idx
%       corresponding indexes in data
%   sanity check plot
%       plots sample trajectories with identified end time
% EXAMPLE
%   getMovEndFrame(1,traj, trialIdxs.trialScore,trialIdxs.nbase, movStrt.nbase) will find
%   the end frames for all trials with a score of 1 during baseline
% HISTORY
%   Reagan Bullins 11.12.2021
%% Input Parser
p = inputParser;
addParameter(p,'plotSanity',0,@isnumeric);
parse(p,varargin{:});
plotSanity = p.Results.plotSanity;
%%
idx_specified = [nan];
for itrialScore = 1:length(trialScoreThr)
    thisScore_idx = find(trialScores == trialScoreThr(itrialScore));
    idx_specified = [idx_specified thisScore_idx];
end
if length(trialScoreThr) > 1
    idx_specified(1) = [];
end
% Find where the block condition is met (baseline vs stimulation)
[~,block_score_idx,~] = intersect(trialBlock,idx_specified);
score_total_idx = trialBlock(block_score_idx);
movEndFrames = NaN(1,length(score_total_idx));
if plotSanity
    figure;
end
for itrial = 1:length(score_total_idx)
    thisTraj = traj(score_total_idx(itrial),:,:);
    z = squeeze(traj(score_total_idx(itrial),3,:));
    z_smooth = smoothdata(z,'movmedian',60);
    z_deriv = diff(z_smooth);
    z_deriv = smoothdata(z_deriv,'movmedian',60);
    max_z_deriv = max(z_deriv);
    thr = max_z_deriv*.1;
    overthr_idx = find(z_deriv <= thr);
    attenuated_idx = overthr_idx(find(overthr_idx > movStartFrames(1,itrial)));
    %             overbaseThr_idx = find(z < -30);
    %             [~,thr_idx,~] = intersect(overthr_idx,overbaseThr_idx);
    movEndFrames(1,itrial) = attenuated_idx(1);
    if plotSanity && itrial <=16
        subplot(4,4,itrial)
        plot(z, 'Color',[.5 0 .5]);
          hold on;
         xline(movEndFrames(1,itrial));
    end
end
end