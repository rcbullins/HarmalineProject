%% Load opto behavior data.

% These files are the outputs of DCN_hand_traj_neural.m and
% FIGURE_hand_traj_neural.m.

load /Volumes/hantmanlab/britton/data/reach_cerebellum/processed_DCN_BPN_stim/behav_stim_PN_DCN
load /Volumes/hantmanlab/britton/data/reach_cortex_BPN_perturbation/processed/behav_stim_PN_CTX

% REMOVE REDUNDANT DOUBLE-RECORDING SESSIONS!
foo = [];
for i = 1:length(behav_stim_PN_DCN.labels_dset)
    bar = strcmp(behav_stim_PN_CTX.labels_dset,behav_stim_PN_DCN.labels_dset{i});
    if sum(bar) > 0
        foo(size(foo,1)+1,:) = [find(bar) i];
    end
end
i2keep_dcn = setdiff(1:length(behav_stim_PN_DCN.labels_dset),foo(:,2));

type_all = [behav_stim_PN_CTX.type_all behav_stim_PN_DCN.type_all(i2keep_dcn)];
trials_control = [behav_stim_PN_CTX.trials_control behav_stim_PN_DCN.trials_control(i2keep_dcn)];
trials_laser = [behav_stim_PN_CTX.trials_laser behav_stim_PN_DCN.trials_laser(i2keep_dcn)];
trials_all = [behav_stim_PN_CTX.trials_all behav_stim_PN_DCN.trials_all(i2keep_dcn)];
trials_control_ALL = [behav_stim_PN_CTX.trials_control_ALL behav_stim_PN_DCN.trials_control_ALL(i2keep_dcn)];
trials_laser_ALL = [behav_stim_PN_CTX.trials_laser_ALL behav_stim_PN_DCN.trials_laser_ALL(i2keep_dcn)];
labels_dset = [behav_stim_PN_CTX.labels_dset' behav_stim_PN_DCN.labels_dset(i2keep_dcn)'];
succ_type = [behav_stim_PN_CTX.succ_type behav_stim_PN_DCN.succ_type(i2keep_dcn)];
t2lift_all = [behav_stim_PN_CTX.t2lift_all behav_stim_PN_DCN.t2lift_all(i2keep_dcn)];
t2grab_all = [behav_stim_PN_CTX.t2grab_all behav_stim_PN_DCN.t2grab_all(i2keep_dcn)];
t2handopen_all = [behav_stim_PN_CTX.t2handopen_all behav_stim_PN_DCN.t2handopen_all(i2keep_dcn)];
video_root = [behav_stim_PN_CTX.video_root; behav_stim_PN_DCN.video_root(i2keep_dcn)'];
t2cue = [behav_stim_PN_CTX.t_first_frame_to_cue behav_stim_PN_DCN.t_first_frame_to_cue(i2keep_dcn)];

clear behav_stim_PN_CTX behav_stim_PN_DCN

%% Load finger tracking data and set up parameters.
calib_file = [repmat({'/Volumes/hantmanlab/britton/reports/setups/camera_calibration/camera_calibration_probes_6_24_16/Calib_Results_stereo.mat'},32,1); ...
    repmat({'/Volumes/hantmanlab/britton/reports/setups/camera_calibration/camera_calibration_corner_5-22-19/Calib_Results_stereo.mat'},6,1)];
n_frames = 1994;
w_frame = [384*ones(1,32) 448*ones(1,6)];
dset_ex = 1;
g_traj = flt(1,'gaussian');

load([video_root{dset_ex} '/tracked_fingers/tracked_comb.mat']);
load(calib_file{dset_ex});

%% Pull out trials and triangulate.
n_trials = size(pred_locs,1)/n_frames;
for i = 1:n_trials
    
    i_tmp = (1:n_frames) + (i-1)*n_frames;
    
    for j = 1:4
        % Stereo triangulation.
        xy_s = squeeze(pred_locs(i_tmp,j,:))';
        xy_f = squeeze(pred_locs(i_tmp,j+4,:))';
        xy_f(2,:) = xy_f(2,:) - w_frame(dset_ex);% Subtract x offset for front cam, because combined movie used.
        foo = stereo_triangulation(xy_s,xy_f,...
                om,T,fc_left,cc_left,kc_left,alpha_c_left,fc_right,cc_right,kc_right,alpha_c_right);
        foo = foo([1 3 2],:); % Permute dims.
        foo(3,:) = -foo(3,:); % Invert dim.
        for k = 1:3
            traj_fingers(i,j,k,:) = conv(squeeze(foo(k,:)),g_traj,'same');
            traj_fingers(i,j,k,1:5) = foo(k,6);
            traj_fingers(i,j,k,(end-4):end) = foo(k,end-5);
        end
        
        % Confidence score (min of two views for each object).
        conf_fingers(i,j,:) = min(pred_conf(i_tmp,j+[0 4]),[],2);
    end
        
end

% Insert NaN into low-confidence points.
traj_fingers(conf_fingers < .9) = NaN;

%% Plot cue-aligned position for every 5th trial.
%{
my_col = jet(4);
i2plot = 300:500;
t2plot = 2*i2plot;
figure('position',[0 0 1000 1000]);
for i = 1:5:size(traj_fingers,1)
    for j = 1:4
        plot3(squeeze(traj_fingers(i,j,1,i2plot)),squeeze(traj_fingers(i,j,2,i2plot)),...
            squeeze(traj_fingers(i,j,3,i2plot)),'color',my_col(j,:)); hold on
    end
end
xlim(5 + [-35 35])
ylim(135 + [-35 35])
zlim(2.5 + [-35 35]);
set(gca,'xtick',[],'ytick',[],'ztick',[]);
set_fig_white
%}

%% Get traj peri lift.
i_peri_lift = -125:500;
t_peri_lift = 2*i_peri_lift;

% Control.
traj_lift_c = NaN(length(trials_control_ALL{dset_ex}),size(traj_fingers,2),size(traj_fingers,3),length(t_peri_lift));
for i = 1:length(trials_control_ALL{dset_ex})
    ii = trials_control_ALL{dset_ex}(i);
    i_tmp = round((t2lift_all{dset_ex}(ii) + t2cue(dset_ex))/2);
    if i_tmp > 200 && i_tmp < 1000
        traj_lift_c(i,:,:,:) = traj_fingers(ii,:,:,i_tmp+i_peri_lift);
    end
end

% Laser.
traj_lift_l = NaN(length(trials_laser_ALL{dset_ex}),size(traj_fingers,2),size(traj_fingers,3),length(t_peri_lift));
for i = 1:length(trials_laser_ALL{dset_ex})
    ii = trials_laser_ALL{dset_ex}(i);
    i_tmp = round((t2lift_all{dset_ex}(ii) + t2cue(dset_ex))/2);
    if i_tmp > 200 && i_tmp < 1000
        traj_lift_l(i,:,:,:) = traj_fingers(ii,:,:,i_tmp+i_peri_lift);
    end
end

%% Plot position vs time (not so useful, since they're close).
col_c = [.4:.2:1; .4:.2:1; 0:.2:.6]';
col_l = [0:.2:.6; .4:.2:1; .4:.2:1]';
figure('position',[0 0 2000 900])
for i = 1:3
    subplot(3,2,2*i-1)
    for j = 1:4
        plot(t_peri_lift,squeeze(traj_lift_c(:,j,i,:)),'color',col_c(j,:)); hold on
    end
    xlim(t_peri_lift([1 end]))
    foo = ylim;
    foo = mean(foo)+[-50 50];
    ylim(mean(foo)+[-50 50])
    set(gca,'tickdir','out','xtick',-250:250:1000,'ytick',foo(1):25:foo(2),'yticklabel',[])
end
for i = 1:3
    subplot(3,2,2*i)
    for j = 1:4
        plot(t_peri_lift,squeeze(traj_lift_l(:,j,i,:)),'color',col_l(j,:)); hold on
    end
    xlim(t_peri_lift([1 end]))
    foo = ylim;
    foo = mean(foo)+[-50 50];
    ylim(mean(foo)+[-50 50])
    set(gca,'tickdir','out','xtick',-250:250:1000,'ytick',foo(1):25:foo(2),'yticklabel',[])
end

%% Plot position 3D.
i2plot_3D = find(t_peri_lift >= -50 & t_peri_lift <= 200);
figure('position',[0 0 1000 1000])
for i = 1:size(traj_lift_c,1)
    for j = 1:4
        plot3(squeeze(traj_lift_c(i,j,1,i2plot_3D)),squeeze(traj_lift_c(i,j,2,i2plot_3D)),...
            squeeze(traj_lift_c(i,j,3,i2plot_3D)),'color',col_c(j,:)); hold on
    end
end
for i = 1:size(traj_lift_l,1)
    for j = 1:4
        plot3(squeeze(traj_lift_l(i,j,1,i2plot_3D)),squeeze(traj_lift_l(i,j,2,i2plot_3D)),...
            squeeze(traj_lift_l(i,j,3,i2plot_3D)),'color',col_l(j,:)); hold on
    end
end
xlim(5 + [-35 35])
ylim(135 + [-35 35])
zlim(2.5 + [-35 35]);
set(gca,'xtick',[],'ytick',[],'ztick',[]);
set_fig_white
view([46.2 12.6]);

%% Plot mean trajectories for control and laser.
mean_c = squeeze(nanmean(traj_lift_c,1));
mean_l = squeeze(nanmean(traj_lift_l,1));
t_curves = 0:50:200;
for i = 1:length(t_curves)
    i_curves(i) = find(t_peri_lift == t_curves(i));
end
figure('position',[0 0 1000 1000])
%whitebg([.8 .8 .8])
for j = 1:4
    plot3(squeeze(mean_c(j,1,i2plot_3D)),squeeze(mean_c(j,2,i2plot_3D)),...
        squeeze(mean_c(j,3,i2plot_3D)),'color',col_c(j,:),'linewidth',1.5); hold on
end
for j = 1:4
	plot3(squeeze(mean_l(j,1,i2plot_3D)),squeeze(mean_l(j,2,i2plot_3D)),...
        squeeze(mean_l(j,3,i2plot_3D)),'color',col_l(j,:),'linewidth',1.5); hold on
end
xlim(10 + [-30 30])
ylim(136 + [-30 30])
zlim(0 + [-30 30]);
set(gca,'xtick',[],'ytick',[],'ztick',[]);
set_fig_white
view([46.2 12.6]);

%% Grasp aperture (width) and hand angle (coronal plane).
diff_tmp_c = squeeze(traj_lift_c(:,1,:,:)) - squeeze(traj_lift_c(:,4,:,:));
w_hand_c = sqrt(squeeze(sum(diff_tmp_c.^2,2)));
a_hand_c = pi+atan2(squeeze(diff_tmp_c(:,3,:)),squeeze(diff_tmp_c(:,2,:)));
a_hand_c(a_hand_c > pi) = a_hand_c(a_hand_c > pi) - 2*pi;

diff_tmp_l = squeeze(traj_lift_l(:,1,:,:)) - squeeze(traj_lift_l(:,4,:,:));
w_hand_l = sqrt(squeeze(sum(diff_tmp_l.^2,2)));
a_hand_l = pi+atan2(squeeze(diff_tmp_l(:,3,:)),squeeze(diff_tmp_l(:,2,:)));
a_hand_l(a_hand_l > pi) = a_hand_l(a_hand_l > pi) - 2*pi;

figure('position',[0 0 800 340])
plot(t_peri_lift,a_hand_c,'color',[.7 .7 0]); hold on
plot(t_peri_lift,a_hand_l,'color',[0 .7 .7]); hold on
xlim(t_peri_lift([1 end]))
ylim([-pi pi])
set(gca,'tickdir','out','xtick',-250:250:1000,'ytick',-pi:pi/2:pi,'yticklabel',[]);
set_fig_white

figure('position',[0 0 800 340])
plot(t_peri_lift,w_hand_c,'color',[.7 .7 0]); hold on
plot(t_peri_lift,w_hand_l,'color',[0 .7 .7]); hold on
xlim(t_peri_lift([1 end]))
ylim([0 15])
set(gca,'tickdir','out','xtick',-250:250:1000,'ytick',0:5:15,'yticklabel',[]);
set_fig_white