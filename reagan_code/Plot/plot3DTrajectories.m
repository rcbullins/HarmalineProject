function [] = plot3DTrajectories(traj, varargin)
% PURPOSE
%   Plot example trials of the data in 3D space.
% INPUT
%   trajectory infromation (output from get_traj_3D_v2)
%   startFrame (default = 1)
%   endFrame   (default = 1500)
% OUTPUT
%   subplot containing sample trials from the dataset of 3D trajectory.
% DEPENDENCIES
%   bitton code
% HISTORY
%   Adapted from Kevin: kinematic_analysis_sandbox_v2.m
%   Reagan Bullins: 11.11.2021
%% Input Parsers
p = inputParser;
addParameter(p,'startFrame',1,@isnumeric);
addParameter(p,'endFrame',1500,@isnumeric);
parse(p,varargin{:});
startFrame = p.Results.startFrame;
endFrame = p.Results.endFrame;
%% plot 3d trajectories
color_map=colorGradient([0,255,255],[255,0,255],16);
figure
hold on
for i=1:16
    subplot(4,4,i)
    hold on;
    tmp_traj=squeeze(traj(i,:,startFrame:endFrame));
    size(tmp_traj)
    plot3(tmp_traj(1,:),tmp_traj(2,:),tmp_traj(3,:),'Color',color_map(i,:))
    hold on;
    plot3(tmp_traj(1,1),tmp_traj(2,1),tmp_traj(3,1),'ks','MarkerSize',20,'MarkerFaceColor',color_map(i,:)) %plot starting point
end
end