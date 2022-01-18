clear ALL;
close ALL;

%% load trk data and calibration file
calib_file='C:\Users\Hantman Lab\Documents\Kevin\setup\camera_calibration_Jay_7-28-16\Calib_Results_stereo.mat';
trk_files='D:\harmaline_data\M340necab1_Chr2\20210813movies\trk';
foo = strfind(trk_files,'\');
dir_ctx = trk_files(1:(foo(end)-1));

nframes=3500; %number of frames to analyze. Going to depend on how many frames in movie
[traj conf] = get_traj_3D_v2(trk_files,calib_file,nframes);
% traj contains the hand 3d position in space. traj is a 3 dimension matrix with trial, coordinate and time as the
% first, second and third index respectively

% conf is the confidence value assigned by the classifier


%% data sets and trials where magic happens
% Set file names and directories.
nbase=1:31;  
npert=32:81;  
nwash=82:165; 

%% assign trials based on trial type
base_reach_traj=traj(nbase,:,:);
pert_reach_traj=traj(npert,:,:);
wash_reach_traj=traj(nwash,:,:);
%% plot individual trials on the z axis
n_coord=3; %which coordinate to look at, 1=ant/post position, 2=medial/lateral, 3=height

%frames to plot
start_fr=1; 
end_fr=1500;
%set up plotting parameters
nbins=floor(length(nbase)/5);
color_map=colorGradient([0,255,255],[255,0,255],16);

%plot first 16 baseline control trials
figure
hold on
for i=1:16
    subplot(4,4,i)
    hold on   
    tmp_traj=squeeze(base_reach_traj(i,n_coord,start_fr:end_fr));
    plot((start_fr:end_fr)*2,tmp_traj,'Color',color_map(i,:))
end
title('Handpaths during baseline trials')


%plot first 16 stimulation trials
figure
hold on
for i=1:16
    subplot(4,4,i)
    hold on  
    tmp_traj=squeeze(pert_reach_traj(i,n_coord,start_fr:end_fr));
    plot((start_fr:end_fr)*2,tmp_traj,'Color',color_map(i,:))
end
title('Handpaths during stim trials')    



%plot first 16 washout trials
figure
hold on
for i=1:16
    subplot(4,4,i)
    hold on   
    tmp_traj=squeeze(wash_reach_traj(i,n_coord,start_fr:end_fr));
    plot((start_fr:end_fr)*2,tmp_traj,'Color',color_map(i,:))

end
title('Handpaths during washout trials')


%% plot 3d trajectories
%plot first 16 baseline control trials
color_map=colorGradient([0,255,255],[255,0,255],16);
figure
hold on
for i=1:16
    subplot(4,4,i)
    hold on  
    tmp_traj=squeeze(base_reach_traj(i,:,start_fr:end_fr));
    plot3(tmp_traj(1,:),tmp_traj(2,:),tmp_traj(3,:),'Color',color_map(i,:))
    plot3(tmp_traj(1,1),tmp_traj(2,1),tmp_traj(3,1),'ks','MarkerSize',20,'MarkerFaceColor',color_map(i,:)) %plot starting point
end
title('Handpaths during baseline trials')


%plot first 16 stimulation trials
figure
hold on
count=1;
for i=1:16 
    subplot(4,4,i)
    hold on  
    tmp_traj=squeeze(pert_reach_traj(i,:,start_fr:end_fr));
    plot3(tmp_traj(1,:),tmp_traj(2,:),tmp_traj(3,:),'Color',color_map(i,:))
    plot3(tmp_traj(1,1),tmp_traj(2,1),tmp_traj(3,1),'ks','MarkerSize',20,'MarkerFaceColor',color_map(i,:)) %plot starting point
end
title('Handpaths during stim trials')




%plot first 16 washout trials
figure
hold on
for i=1:16
    subplot(4,4,i)
    hold on  
    tmp_traj=squeeze(wash_reach_traj(i,:,start_fr:end_fr));
    plot3(tmp_traj(1,:),tmp_traj(2,:),tmp_traj(3,:),'Color',color_map(i,:))
    plot3(tmp_traj(1,1),tmp_traj(2,1),tmp_traj(3,1),'ks','MarkerSize',20,'MarkerFaceColor',color_map(i,:)) %plot starting point
end
title('Handpaths during washout trials')
