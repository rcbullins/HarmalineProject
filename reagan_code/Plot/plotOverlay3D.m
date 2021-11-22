<<<<<<< HEAD
<<<<<<< HEAD
function [] = plotOverlay3D(movStartFrames, movEndFrames,movIdx,traj, varargin)
% PURPOSE
%   Plot overlay of 3 deminsional movements for given start frames to end
%   frames
% INPUT
%   movStartFrames (vector of start movement times)
%       output from getMovStartFrame
%   movEndFrames (vector of end frames)
%   movIdx (movement indexes corresponding to frames)
%   traj
%       output from get_traj_3D_v2
%   numFramesPlot (default = 1500) how many frames to plot
%   plotEntireReach (default = 0) 1 will plot only reach to grab
% OUTPUT
%   plot of overlayed movements for all specified start frames
% HISTORY
%   Reagan Bullins 11.12.2021

%% Input Parser
p = inputParser;
addParameter(p,'numFramesPlot',1500,@isnumeric);
addParameter(p,'plotEntireReach',0,@isnumeric);
parse(p,varargin{:});
numFramesPlot = p.Results.numFramesPlot;
plotEntireReach = p.Results.plotEntireReach;
%%
color_map=colorGradient([0,255,255],[255,0,255],length(movStartFrames));
figure;
for itrial = 1:length(movStartFrames)
    frameStart = movStartFrames(1,itrial);
    if plotEntireReach
        frameEnd = frameStart + numFramesPlot;
    else
        frameEnd = movEndFrames(1,itrial);
    end
    if frameEnd > length(traj(movIdx(itrial),1,:))
        frameEnd = length(traj(movIdx(itrial),1,:));
    end
    
    thisTraj = squeeze(traj(movIdx(itrial),:,frameStart:frameEnd));
    plot3(thisTraj(1,:),thisTraj(2,:),thisTraj(3,:),'Color',color_map(itrial,:));
    hold on;
    %plot3(thisTraj(1,1),thisTraj(2,1),thisTraj(3,1),'ks','MarkerSize',20,'MarkerFaceColor',color_map(i,:)) 
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
end
=======
function [] = plotOverlay3D(movStartFrames, movEndFrames,movIdx,traj, varargin)
% PURPOSE
%   Plot overlay of 3 deminsional movements for given start frames to end
%   frames
% INPUT
%   movStartFrames (vector of start movement times)
%       output from getMovStartFrame
%   movEndFrames (vector of end frames)
%   movIdx (movement indexes corresponding to frames)
%   traj
%       output from get_traj_3D_v2
%   numFramesPlot (default = 1500) how many frames to plot
%   plotEntireReach (default = 0) 1 will plot only reach to grab
% OUTPUT
%   plot of overlayed movements for all specified start frames
% HISTORY
%   Reagan Bullins 11.12.2021

%% Input Parser
p = inputParser;
addParameter(p,'numFramesPlot',1500,@isnumeric);
addParameter(p,'plotEntireReach',0,@isnumeric);
parse(p,varargin{:});
numFramesPlot = p.Results.numFramesPlot;
plotEntireReach = p.Results.plotEntireReach;
%%
color_map=colorGradient([0,255,255],[255,0,255],length(movStartFrames));
figure;
for itrial = 1:length(movStartFrames)
    frameStart = movStartFrames(1,itrial);
    if plotEntireReach
        frameEnd = frameStart + numFramesPlot;
    else
        frameEnd = movEndFrames(1,itrial);
    end
    if frameEnd > length(traj(movIdx(itrial),1,:))
        frameEnd = length(traj(movIdx(itrial),1,:));
    end
    
    thisTraj = squeeze(traj(movIdx(itrial),:,frameStart:frameEnd));
    plot3(thisTraj(1,:),thisTraj(2,:),thisTraj(3,:),'Color',color_map(itrial,:));
    hold on;
    %plot3(thisTraj(1,1),thisTraj(2,1),thisTraj(3,1),'ks','MarkerSize',20,'MarkerFaceColor',color_map(i,:)) 
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
end
>>>>>>> da7b9ddee33fff3453e6f5e0546d20d35e92172a
=======
function [] = plotOverlay3D(movStartFrames, movEndFrames,movIdx,traj, varargin)
% PURPOSE
%   Plot overlay of 3 deminsional movements for given start frames to end
%   frames
% INPUT
%   movStartFrames (vector of start movement times)
%       output from getMovStartFrame
%   movEndFrames (vector of end frames)
%   movIdx (movement indexes corresponding to frames)
%   traj
%       output from get_traj_3D_v2
%   numFramesPlot (default = 1500) how many frames to plot
%   plotEntireReach (default = 0) 1 will plot only reach to grab
% OUTPUT
%   plot of overlayed movements for all specified start frames
% HISTORY
%   Reagan Bullins 11.12.2021

%% Input Parser
p = inputParser;
addParameter(p,'numFramesPlot',1500,@isnumeric);
addParameter(p,'plotEntireReach',0,@isnumeric);
parse(p,varargin{:});
numFramesPlot = p.Results.numFramesPlot;
plotEntireReach = p.Results.plotEntireReach;
%%
color_map=colorGradient([0,255,255],[255,0,255],length(movStartFrames));
figure;
for itrial = 1:length(movStartFrames)
    frameStart = movStartFrames(1,itrial);
    if plotEntireReach
        frameEnd = frameStart + numFramesPlot;
    else
        frameEnd = movEndFrames(1,itrial);
    end
    if frameEnd > length(traj(movIdx(itrial),1,:))
        frameEnd = length(traj(movIdx(itrial),1,:));
    end
    
    thisTraj = squeeze(traj(movIdx(itrial),:,frameStart:frameEnd));
    plot3(thisTraj(1,:),thisTraj(2,:),thisTraj(3,:),'Color',color_map(itrial,:));
    hold on;
    %plot3(thisTraj(1,1),thisTraj(2,1),thisTraj(3,1),'ks','MarkerSize',20,'MarkerFaceColor',color_map(i,:)) 
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
end
>>>>>>> a39e1c0843197dbc2ab958057a9f139f35e0535d
end