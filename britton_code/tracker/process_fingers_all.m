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

%% Set parameters.
calib_file = [repmat({'/Volumes/hantmanlab/britton/reports/setups/camera_calibration/camera_calibration_probes_6_24_16/Calib_Results_stereo.mat'},32,1); ...
    repmat({'/Volumes/hantmanlab/britton/reports/setups/camera_calibration/camera_calibration_corner_5-22-19/Calib_Results_stereo.mat'},length(video_root)-32,1)];
n_frames = 1994*ones(1,length(video_root));
w_frame = [384*ones(1,32) 448*ones(1,length(video_root)-32)];
g_traj = flt(1,'gaussian');

%% Iterate over sessions.
dset2proc = 1:length(video_root);
%dset2proc = setdiff(1:32,[11 20 22]);
 for dd = 1:length(dset2proc)
dset = dset2proc(dd);

load([video_root{dset} '/tracked_fingers/tracked_comb.mat']);
load(calib_file{dset});

%% Pull out trials and triangulate.
%n_trials = size(pred_locs,1)/n_frames(dset);
load([video_root{dset} '/tracked_fingers/job_info.mat']);
traj_fingers = NaN(length(t2lift_all{dset}),4,3,n_frames(dset));
conf_fingers = NaN(length(t2lift_all{dset}),4,n_frames(dset));
for ii = 1:length(trial)
    
    i = trial(ii);
    
    i_tmp = (1:n_frames(dset)) + (ii-1)*n_frames(dset);
    
	% Confidence score (min of two views for each object).
    for j = 1:4
        conf_fingers(i,j,:) = min(pred_conf(i_tmp,j+[0 4]),[],2);
    end
    i2del = squeeze(min(conf_fingers(i,[1 4],:),[],2)) < .99;

    for j = 1:4
        % Stereo triangulation.
        xy_s = squeeze(pred_locs(i_tmp,j,:))';
        xy_f = squeeze(pred_locs(i_tmp,j+4,:))';
        xy_f(2,:) = xy_f(2,:) - w_frame(dset);% Subtract x offset for front cam, because combined movie used.
        foo = stereo_triangulation(xy_s,xy_f,...
                om,T,fc_left,cc_left,kc_left,alpha_c_left,fc_right,cc_right,kc_right,alpha_c_right);
        foo = foo([1 3 2],:); % Permute dims.
        foo(3,:) = -foo(3,:); % Invert dim.
        foo(:,i2del) = NaN; % Delete low-confidence points.
        for k = 1:3
            traj_fingers(i,j,k,:) = conv(squeeze(foo(k,:)),g_traj,'same');
            traj_fingers(i,j,k,1:5) = foo(k,6);
            traj_fingers(i,j,k,(end-4):end) = foo(k,end-5);
        end

    end    
    
end

%% Get traj peri lift.
i_peri_lift = -125:500;
t_peri_lift = 2*i_peri_lift;

% Control.
traj_lift_c = NaN(length(trials_control_ALL{dset}),size(traj_fingers,2),size(traj_fingers,3),length(t_peri_lift));
for i = 1:length(trials_control_ALL{dset})
    ii = trials_control_ALL{dset}(i);
    i_tmp = round((t2lift_all{dset}(ii) + t2cue(dset))/2) + i_peri_lift;
    i_insert = ~(i_tmp < 1 | i_tmp > n_frames(dset) | isnan(i_tmp));
    
    traj_lift_c(i,:,:,i_insert) = traj_fingers(ii,:,:,i_tmp(i_insert));
end

% Laser.
traj_lift_l = NaN(length(trials_laser_ALL{dset}),size(traj_fingers,2),size(traj_fingers,3),length(t_peri_lift));
for i = 1:length(trials_laser_ALL{dset})
    ii = trials_laser_ALL{dset}(i);
	i_tmp = round((t2lift_all{dset}(ii) + t2cue(dset))/2) + i_peri_lift;
    i_insert = ~(i_tmp < 1 | i_tmp > n_frames(dset) | isnan(i_tmp));
    
    traj_lift_l(i,:,:,i_insert) = traj_fingers(ii,:,:,i_tmp(i_insert));
end


%% Grasp aperture (width) and hand angle (coronal plane).
diff_tmp_c = squeeze(traj_lift_c(:,1,:,:)) - squeeze(traj_lift_c(:,4,:,:));
w_hand_c = sqrt(squeeze(sum(diff_tmp_c.^2,2)));
a_hand_c = pi+atan2(squeeze(diff_tmp_c(:,3,:)),squeeze(diff_tmp_c(:,2,:)));
a_hand_c(a_hand_c > pi) = a_hand_c(a_hand_c > pi) - 2*pi;

diff_tmp_l = squeeze(traj_lift_l(:,1,:,:)) - squeeze(traj_lift_l(:,4,:,:));
w_hand_l = sqrt(squeeze(sum(diff_tmp_l.^2,2)));
a_hand_l = pi+atan2(squeeze(diff_tmp_l(:,3,:)),squeeze(diff_tmp_l(:,2,:)));
a_hand_l(a_hand_l > pi) = a_hand_l(a_hand_l > pi) - 2*pi;

%% Package into output structure.
out(dd).traj_lift_c = traj_lift_c;
out(dd).traj_lift_l = traj_lift_l;
out(dd).a_hand_c = a_hand_c;
out(dd).a_hand_l = a_hand_l;
out(dd).w_hand_c = w_hand_c;
out(dd).w_hand_l = w_hand_l;
out(dd).conf_fingers = conf_fingers;
out(dd).dset = dset;

display(['Finished dataset ' num2str(dset)]);

end

%% Save.
save('/Volumes/hantmanlab/britton/data/reach_cortex_BPN_perturbation/processed/fingers_opto_ALL','out');