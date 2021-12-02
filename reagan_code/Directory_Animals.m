
% Prerun Script, set up animal list and behavior videos
%% List animals
animals = {'M340';'M341'};

%% Videos/Days
% List control videos for each animal in same order (if more than one, need
% to rewrite code a bit)
M340_controlBehaviorVideos = {'20210811'}; %is 18 a control?
M341_controlBehaviorVideos = {'20210818'};
% List harmaline videos for each animal in same order
M340_harmBehaviorVideos = {'20210813'};
M341_harmBehaviorVideos = {'20210819'};

%% Path to excel file info
m340_excelFile = 'D:\rbullins\Data\M340necab1_Chr2\M340Reachresults.xlsx';
m341_excelFile = 'D:\rbullins\Data\M341necab1_Chr2\M341Reachresults.xlsx';
%% Control trials - selected by hand from excel sheets     
% m340-20210811
M340_20210811_controlTrials.trialScore =  xlsread(m340_excelFile,'B133:EZ133');
M340_20210811_controlTrials.nbase = 1:38;
M340_20210811_controlTrials.npert = [39:81 83:89];
M340_20210811_controlTrials.nwash = 90:155;
% m341-20210818
M341_20210818_controlTrials.trialScore = xlsread(m341_excelFile,'B177:EM177');
M341_20210818_controlTrials.nbase = 1:30;
M341_20210818_controlTrials.npert = [31:55 57:58 60:71 73:86 88:90];
M341_20210818_controlTrials.nwash = 91:142;
%% Harmaline Trials
% m340-20210813
M340_20210813_harmTrials.trialScore = xlsread(m340_excelFile,'B147:FJ147');
M340_20210813_harmTrials.nbase = 1:31;
M340_20210813_harmTrials.npert = [32:34 36:81];
M340_20210813_harmTrials.nwash = 82:165;
% 341 -20210819
M341_20210819_harmTrials.trialScore = xlsread(m341_excelFile,'B205:KT205');
M341_20210819_harmTrials.nbase = 1:305;
M341_20210819_harmTrials.npert = [];
M341_20210819_harmTrials.nwash = [];

