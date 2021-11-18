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

load /Volumes/hantmanlab/britton/data/reach_cortex_BPN_perturbation/processed/fingers_opto_ALL

dset_ex = 1;
i_peri_lift = -125:500;
t_peri_lift = 2*i_peri_lift;

%% Plot position vs time?
%col_c = [.4:.2:1; .4:.2:1; 0:.2:.6]';
%col_l = [0:.2:.6; .4:.2:1; .4:.2:1]';

col_c = [.35 .35 0; .5 .5 0; .65 .65 0; .8 .8 0];
col_l = [0 .35 .35; 0 .5 .5; 0 .65 .65; 0 .8 .8];

%{
figure('position',[0 0 2560 1440])
i = 1; % Which spatial dim. to plot.
for kk = 1:length(out)
    
    subplot(5,9,kk)

    for j = 1:4
        plot(t_peri_lift,squeeze(out(kk).traj_lift_c(:,j,i,:)),'color',col_c(j,:)); hold on
    end
    xlim(t_peri_lift([1 end]))
    foo = ylim;
    foo = mean(foo)+[-50 50];
    ylim(mean(foo)+[-50 50])
    set(gca,'tickdir','out','xtick',-250:250:1000,'ytick',foo(1):25:foo(2),'yticklabel',[])

    for j = 1:4
        plot(t_peri_lift,squeeze(out(kk).traj_lift_l(:,j,i,:)),'color',col_l(j,:)); hold on
    end
    xlim(t_peri_lift([1 end]))
    foo = ylim;
    foo = mean(foo)+[-50 50];
    ylim(mean(foo)+[-50 50])
    set(gca,'tickdir','out','xtick',-250:250:1000,'ytick',foo(1):25:foo(2),'yticklabel',[])
    
    title(labels_dset{out(kk).dset});
    
end
%}
%% Plot position 3D.
%{
i2plot_3D = find(t_peri_lift >= -50 & t_peri_lift <= 200);
figure('position',[0 0 1000 1000])
for i = 1:size(out(dset_ex).traj_lift_c,1)
    for j = 1:4
        plot3(squeeze(out(dset_ex).traj_lift_c(i,j,1,i2plot_3D)),squeeze(out(dset_ex).traj_lift_c(i,j,2,i2plot_3D)),...
            squeeze(out(dset_ex).traj_lift_c(i,j,3,i2plot_3D)),'color',col_c(j,:)); hold on
    end
end
for i = 1:size(out(dset_ex).traj_lift_l,1)
    for j = 1:4
        plot3(squeeze(out(dset_ex).traj_lift_l(i,j,1,i2plot_3D)),squeeze(out(dset_ex).traj_lift_l(i,j,2,i2plot_3D)),...
            squeeze(out(dset_ex).traj_lift_l(i,j,3,i2plot_3D)),'color',col_l(j,:)); hold on
    end
end
xlim(5 + [-35 35])
ylim(135 + [-35 35])
zlim(2.5 + [-35 35]);
set(gca,'xtick',[],'ytick',[],'ztick',[]);
set_fig_white
view([46.2 12.6]);
%}

%% Plot mean trajectories for control and laser.
%{
figure('position',[0 0 2560 1440])
i2plot_3D = find(t_peri_lift >= -50 & t_peri_lift <= 200);
for kk = 1:length(out)
    subplot(5,9,kk);
    mean_c = squeeze(nanmean(out(kk).traj_lift_c,1));
    mean_l = squeeze(nanmean(out(kk).traj_lift_l,1));
    t_curves = 0:50:200;
    for i = 1:length(t_curves)
        i_curves(i) = find(t_peri_lift == t_curves(i));
    end
    
    whitebg([.8 .8 .8])
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
    title(labels_dset{out(kk).dset})
end
%}
%% Grasp aperture (width) and hand angle (coronal plane).

figure; whitebg('w'); close

figure('position',[0 0 2560 1440])

for kk = 1:length(out)
    
    subplot(5,9,kk);
    
    plot(t_peri_lift,out(kk).a_hand_c,'color',[.7 .7 0]); hold on
    plot(t_peri_lift,out(kk).a_hand_l,'color',[0 .7 .7]); hold on
    xlim(t_peri_lift([1 end]))
    ylim([-pi pi])
    set(gca,'tickdir','out','xtick',-250:250:1000,'ytick',-pi:pi/2:pi,'yticklabel',[]);

end
set_fig_white


figure('position',[0 0 2560 1440])

for kk = 1:length(out)
    
    subplot(5,9,kk);
    
    plot(t_peri_lift,out(kk).w_hand_c,'color',[.7 .7 0]); hold on
    plot(t_peri_lift,out(kk).w_hand_l,'color',[0 .7 .7]); hold on
    xlim(t_peri_lift([1 end]))
    ylim([0 15])
    set(gca,'tickdir','out','xtick',-250:250:1000,'ytick',0:5:15,'yticklabel',[]);

end
set_fig_white

%% Hand angle at grab. Larger values -> more supinated.
%{
for i = 1:length(out)
    a_hand_c{i} = NaN(1,length(trials_control_ALL{i}));
    for j = 1:length(trials_control_ALL{i})
        foo = t2grab_all{i}(trials_control_ALL{i}(j)) -...
            t2lift_all{i}(trials_control_ALL{i}(j));
        if foo <= t_peri_lift(end) && foo >= 30
            bar = find_closest(foo,t_peri_lift);
            a_hand_c{i}(j) = out(i).a_hand_c(j,bar);
        end
    end
    
    a_hand_l{i} = NaN(1,length(trials_laser_ALL{i}));
    for j = 1:length(trials_laser_ALL{i})
        foo = t2grab_all{i}(trials_laser_ALL{i}(j)) -...
            t2lift_all{i}(trials_laser_ALL{i}(j));
        if foo <= t_peri_lift(end) && foo >= 30
            bar = find_closest(foo,t_peri_lift);
            a_hand_l{i}(j) = out(i).a_hand_l(j,bar);
        end
    end
    
end
for i = 1:length(a_hand_c)
    a_hand_mean_c(i) = circ_mean(a_hand_c{i}(~isnan(a_hand_c{i}))');
    a_hand_mean_l(i) = circ_mean(a_hand_l{i}(~isnan(a_hand_l{i}))');
end
figure('position',[0 0 600 538])
plot(a_hand_mean_c,a_hand_mean_l,'ok'); hold on;
mylim = [0 2];
plot(mylim,mylim,'-r');
xlim(mylim);
ylim(mylim);
xlabel('Average hand angle at grab, control');
ylabel('Average hand angle at grab, laser');
%}

%% Skeleton or surface plots of hand trajectory for an example dataset.
i2plot_3D = find(t_peri_lift >= 0 & t_peri_lift <= 200);
i_highlight = 1:length(i2plot_3D);%find(ismember(t_peri_lift(i2plot_3D),0:50:200));

trials_cand_l = find(squeeze(sum(sum(sum(isnan(out(dset_ex).traj_lift_l(:,:,:,i2plot_3D)),4),3),2)) == 0);
trials_cand_c = find(squeeze(sum(sum(sum(isnan(out(dset_ex).traj_lift_c(:,:,:,i2plot_3D)),4),3),2)) == 0);

trials_3_c = trials_cand_c(1);
trials_3_l = trials_cand_l(5);
%{
% Try skeleton ...
figure('position',[0 0 1200 1200])
whitebg('w')
for i = 1:length(trials_3_c)
    
    % Try just as skeleton - might look better.
	xx = squeeze(out(dset_ex).traj_lift_c(trials_3_c(i),:,1,i2plot_3D))';
    yy = squeeze(out(dset_ex).traj_lift_c(trials_3_c(i),:,2,i2plot_3D))';
    zz = squeeze(out(dset_ex).traj_lift_c(trials_3_c(i),:,3,i2plot_3D))';
    
    for j = 1:4
        plot3(xx(:,j),yy(:,j),zz(:,j),'color',col_c(j,:)); hold on
    end
    for j = 1:length(i_highlight)
        for k = 1:3
            plot3(xx(i_highlight(j),k:k+1),yy(i_highlight(j),k:k+1),zz(i_highlight(j),k:k+1),...
                'color',col_c(k,:)); hold on
        end
    end
    
end

for i = 1:length(trials_3_l)
    
    % Try just as skeleton - might look better.
	xx = squeeze(out(dset_ex).traj_lift_l(trials_3_l(i),:,1,i2plot_3D))';
    yy = squeeze(out(dset_ex).traj_lift_l(trials_3_l(i),:,2,i2plot_3D))';
    zz = squeeze(out(dset_ex).traj_lift_l(trials_3_l(i),:,3,i2plot_3D))';
    
    for j = 1:4
        plot3(xx(:,j),yy(:,j),zz(:,j),'color',col_l(j,:)); hold on
    end
    for j = 1:length(i_highlight)
        for k = 1:3
            plot3(xx(i_highlight(j),k:k+1),yy(i_highlight(j),k:k+1),zz(i_highlight(j),k:k+1),...
                'color',col_l(k,:)); hold on
        end
    end
    
end


xlim(10 + [-30 30])
ylim(136 + [-30 30])
zlim(0 + [-30 30]);
set(gca,'xtick',[],'ytick',[],'ztick',[]);
set_fig_white
view([46.2 12.6]);
%}
%% Try surface?
i_highlight = 1:length(i2plot_3D);%find(ismember(t_peri_lift(i2plot_3D),0:50:200));
figure('position',[0 0 1200 1200])
whitebg('w')

% Control.
for i = 1:length(trials_3_c)
    
	xx = squeeze(out(dset_ex).traj_lift_c(trials_3_c(i),:,1,i2plot_3D))';
    yy = squeeze(out(dset_ex).traj_lift_c(trials_3_c(i),:,2,i2plot_3D))';
    zz = squeeze(out(dset_ex).traj_lift_c(trials_3_c(i),:,3,i2plot_3D))';
    
    
    px = [];
    py = [];
    pz = [];
    for j = 1:size(xx,1)-1
        for k = 1:3
            
            px = [px;...
                [xx(j,k) xx(j+1,k) xx(j,k+1)]; ...
                [xx(j+1,k) xx(j,k+1) xx(j+1,k+1)]];
            
            py = [py;...
                [yy(j,k) yy(j+1,k) yy(j,k+1)]; ...
                [yy(j+1,k) yy(j,k+1) yy(j+1,k+1)]];
            
            pz = [pz;...
                [zz(j,k) zz(j+1,k) zz(j,k+1)]; ...
                [zz(j+1,k) zz(j,k+1) zz(j+1,k+1)]];

        end
    end
    
    for j = 1:size(px,1)
        s_tmp = fill3(px(j,:),py(j,:),pz(j,:),[.8 .8 0]); hold on
        set(s_tmp,'edgecolor','none','facealpha',.7)
    end
    plot3(xx(:,1),yy(:,1),zz(:,1),'color',[.2 .2 0]); hold on
    plot3(xx(:,4),yy(:,4),zz(:,4),'color',[.2 .2 0]); hold on
        
    
    for j = 1:length(i_highlight)
        %{
        for k = 1:3
            plot3(xx(i_highlight(j),k:k+1),yy(i_highlight(j),k:k+1),zz(i_highlight(j),k:k+1),...
                'color',col_c(k,:)); hold on
        end
        %}
        plot3(xx(i_highlight(j),:),yy(i_highlight(j),:),...
            zz(i_highlight(j),:),'color',[.2 .2 0]); hold on
    end
    
end

% Laser.
for i = 1:length(trials_3_l)
    
	xx = squeeze(out(dset_ex).traj_lift_l(trials_3_l(i),:,1,i2plot_3D))';
    yy = squeeze(out(dset_ex).traj_lift_l(trials_3_l(i),:,2,i2plot_3D))';
    zz = squeeze(out(dset_ex).traj_lift_l(trials_3_l(i),:,3,i2plot_3D))';
    
    
    px = [];
    py = [];
    pz = [];
    for j = 1:size(xx,1)-1
        for k = 1:3
            
            px = [px;...
                [xx(j,k) xx(j+1,k) xx(j,k+1)]; ...
                [xx(j+1,k) xx(j,k+1) xx(j+1,k+1)]];
            
            py = [py;...
                [yy(j,k) yy(j+1,k) yy(j,k+1)]; ...
                [yy(j+1,k) yy(j,k+1) yy(j+1,k+1)]];
            
            pz = [pz;...
                [zz(j,k) zz(j+1,k) zz(j,k+1)]; ...
                [zz(j+1,k) zz(j,k+1) zz(j+1,k+1)]];

        end
    end
    
    for j = 1:size(px,1)
        s_tmp = fill3(px(j,:),py(j,:),pz(j,:),[0 .8 .8]); hold on
        set(s_tmp,'edgecolor','none','facealpha',.7)
    end
    plot3(xx(:,1),yy(:,1),zz(:,1),'color',[0 .2 .2]); hold on
    plot3(xx(:,4),yy(:,4),zz(:,4),'color',[0 .2 .2]); hold on
        
    
    for j = 1:length(i_highlight)
        %{
        for k = 1:3
            plot3(xx(i_highlight(j),k:k+1),yy(i_highlight(j),k:k+1),zz(i_highlight(j),k:k+1),...
                'color',col_l(k,:)); hold on
        end
        %}
        plot3(xx(i_highlight(j),:),yy(i_highlight(j),:),...
            zz(i_highlight(j),:),'color',[0 .2 .2]); hold on
    end
    
end

xlim(10 + [-30 30])
ylim(136 + [-30 30])
zlim(0 + [-30 30]);
set(gca,'xtick',[],'ytick',[],'ztick',[]);
set_fig_white
view([-378.4 38.2]);
set(gcf,'renderer','opengl')
my_light = light('Position',[-0.4 0.2 0.9],'Style','infinite');