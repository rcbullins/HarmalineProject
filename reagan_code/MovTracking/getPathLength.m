function [totalPathLength] = getPathLength(traj, movIdx,movStrt, movEnd)
% PURPOSE
%   Find the total path length from movement initiation to end of movement.
% INPUT
%   traj (trajectories of all trials)
%   movIdx (index of trials you want to find the pathlength for)
%   movStrt (start frames of movements)
%   movEnd (end frames of movements)
% OUTPUT
%   totalPathLength (vector of path length per trial in pixel space)
% HISTORY
%   11.23.2021 Reagan Bullins
% EXAMPLE
%         [totalPathLength.nbase] =
%         getPathLength(traj, movIdx.nbase,movStrt.nbase, movEnd.nbase)
%% Get path lengths for defined trials by movIdx
for itrial = 1:length(movIdx)
    % Initiate total length variable
    totalLength = 0;
    % Identify start and end frame of reach
    frameStart = movStrt(1,itrial);
    frameEnd = movEnd(1,itrial);
    % Find total frames within reach movement
    totalFramesTraj = frameEnd-frameStart;
    % For each frame, find the distance between x,y,z for the
    % previous frame and add it to totalPathLength
    for iframe = 1:totalFramesTraj
        iframeTraj = sqrt(sum((traj(movIdx(itrial),:,frameStart-1+iframe)) .^ 2));
        totalLength = totalLength + iframeTraj;
    end
    totalPathLength(itrial) = totalLength;
end
end