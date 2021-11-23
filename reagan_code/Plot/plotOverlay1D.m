function [] = plotOverlay1D(movStartFrames, movEndFrames,movIdx,traj, varargin)
% PURPOSE
%   Plot overlay of one deminsional movements for given start frames
% INPUT
%   movStartFrames (vector of start movement times)
%       output from getMovStartFrame
%   movEndFrames (vecotr of end frames)
%   traj
%       output from get_traj_3D_v2
%   dimPlot (default = 3, z dimension) what dimension to plot in
%   numFramesPlot (default = 1500) how many frames to plot
%   plotEntireReach (default = 0) 1 will plot only reach to grab
% OUTPUT
%   plot of overlayed movements for all specified start frames
% HISTORY
%   Reagan Bullins 11.12.2021

%% Input Parser
p = inputParser;
addParameter(p,'dimPlot',3,@isnumeric);
addParameter(p,'numFramesPlot',1500,@isnumeric);
addParameter(p,'plotEntireReach',0,@isnumeric);
parse(p,varargin{:});
dimPlot = p.Results.dimPlot;
numFramesPlot = p.Results.numFramesPlot;
plotEntireReach = p.Results.plotEntireReach;
%% Plotting
figure;
for itrial = 1:length(movStartFrames)
    frameStart = movStartFrames(1,itrial);
    if plotEntireReach
        frameEnd = frameStart + numFramesPlot;
    else
        frameEnd = movEndFrames(1,itrial);
    end
    if frameEnd > length(traj(movIdx(itrial),dimPlot,:))
        frameEnd = length(traj(movIdx(itrial),dimPlot,:));
    end
    thisTraj = squeeze(traj(movIdx(itrial),dimPlot,frameStart:frameEnd));
    plot(thisTraj);
    hold on;
end
end