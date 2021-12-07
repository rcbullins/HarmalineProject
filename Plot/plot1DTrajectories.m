function [] = plot1DTrajectories(traj, varargin)
%% Behavioral - Motor reaching kinematics analysis script
% PURPOSE
%   Plot 1D trajectories of animal movement
% INPUT
%   trajectory infromation (output from get_traj_3D_v2)
%   coordinate_dim (default = 3 for z axis)
%        example use - plot1DTrajectories(traj,'coordinate_dim',2)
%   startFrame (default = 1)
%   endFrame   (default = 1500)
% OUTPUT
%   sample trial trajectory information (default in z dimension)
% DEPENDENCIES
%   bitton code
% HISTORY
%   Adapted from Kevin: kinematic_analysis_sandbox_v2.m
%   Reagan Bullins: 11.11.2021
%% Input parsers
p = inputParser;
addParameter(p,'coordinate_dim',3,@isnumeric);
addParameter(p,'startFrame',1,@isnumeric);
addParameter(p,'endFrame',1500,@isnumeric);
parse(p,varargin{:});
coordinate_dim    = p.Results.coordinate_dim;
startFrame = p.Results.startFrame;
endFrame = p.Results.endFrame;
%% Plot individual trials on the specified axis
color_map=colorGradient([0,255,255],[255,0,255],16);

%plot first 16 baseline control trials
figure;
hold on;
for i=1:16
    subplot(4,4,i)
    hold on
    tmp_traj=squeeze(traj(i,coordinate_dim,startFrame:endFrame));
    plot((startFrame:endFrame)*2,tmp_traj,'Color',color_map(i,:))
end

end