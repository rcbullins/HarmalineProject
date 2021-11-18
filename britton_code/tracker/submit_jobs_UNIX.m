video_root =    {...
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M270Slc17a7_Chr2_BPN/20180813'         
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M270Slc17a7_Chr2_BPN/20180815'         
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M270Slc17a7_Chr2_BPN/20180816'         
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M271Slc17a7_Chr2_BPN/20180813'         
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M271Slc17a7_Chr2_BPN/20180815'         
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M271Slc17a7_Chr2_BPN/20180816'         
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M238Slc17a7_Chr2/20170824'             
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M238Slc17a7_Chr2/20170825'             
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M238Slc17a7_Chr2/20170828'             
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M240Slc17a7_Chr2/20170919'             
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M240Slc17a7_Chr2/20170920'             
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M240Slc17a7_Chr2/20170921'             
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M240Slc17a7_Chr2/20170922'             
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M240Slc17a7_Chr2/20170925'             
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M241Slc17a7_Chr2/20171003'             
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M241Slc17a7_Chr2/20171004'             
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M241Slc17a7_Chr2/20171005'             
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M241Slc17a7_Chr2/20171006'             
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M241Slc17a7_Chr2/20171009'             
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M241Slc17a7_Chr2/20171010'             
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M281Slc17a7_Chr2_Bpn/20181203'         
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M281Slc17a7_Chr2_Bpn/20181205'         
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M281Slc17a7_Chr2_Bpn/20181207'         
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M282Slc17a7_Chr2_Bpn/20181126'         
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M282Slc17a7_Chr2_Bpn/20181128'         
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M282Slc17a7_Chr2_Bpn/20181203'         
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M282Slc17a7_Chr2_Bpn/20181205'         
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M282Slc17a7_Chr2_Bpn/20181206'         
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M282Slc17a7_Chr2_Bpn/20181207'         
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M282Slc17a7_Chr2_Bpn/20181208'         
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M282Slc17a7_Chr2_Bpn/20181210'         
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M282Slc17a7_Chr2_Bpn/20181212'         
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M302Slc17a7_Chr2_Bpn_DCN/20191209probe'
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M302Slc17a7_Chr2_Bpn_DCN/20191210probe'
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M303Slc17a7_Chr2_Bpn_DCN/20191211probe'
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M303Slc17a7_Chr2_Bpn_DCN/20191212probe'
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M304Slc17a7_Chr2_Bpn_DCN/20191216probe'
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M304Slc17a7_Chr2_Bpn_DCN/20191217probe'
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M288Slc17a7_Chr2_Bpn_DCN/20190619'     
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M291Slc17a7_Chr2_Bpn_DCN/20190626'     
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M291Slc17a7_Chr2_Bpn_DCN/20190627'     
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M292Slc17a7_Chr2_Bpn_DCN/20190702_G4'  
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M292Slc17a7_Chr2_Bpn_DCN/20190703_G1'  
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M295Slc17a7_Chr2_Bpn_DCN/20190909'     
    '/Volumes/hantmanlab/from_tier2/Jay/videos/M296Slc17a7_Chr2_Bpn_DCN/20190925'     };

%dirs2track = 1:32; % Experiments with TALL ROI!!!
dirs2track = 33:length(video_root); % Experiments with SHORT ROI!!!
%dirs2track = [11 20 22]; % Failed sessions; short vids now removed.
for i = 1:length(dirs2track)
    dir_tmp = ['/groups/hantman' video_root{dirs2track(i)}(9:end)];
    load([dir_tmp '/tracked_fingers/job_info.mat'])
    tObj.trackListFile(jsons{1},outfiles{1});
    clear jsons outfiles
    display(['Submitted dataset ' num2str(dirs2track(i))]);
end