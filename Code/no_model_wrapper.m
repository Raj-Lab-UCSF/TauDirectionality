%% 0. Loading data
clear; clc; rng(0);
matdir = '/Users/justintorok/Documents/MATLAB/TauDirectionality_ChrisPaper_Project/TauDirectionality/raw_data_mouse';
load([matdir filesep 'eNDM_mousedata.mat'],'Networks','tpts','netlabs','data426','seed426','labnms');
C_ret = Networks.ret; C_ant = Networks.ant; C_nd = Networks.nd;
taudata_all = struct;
tpts_all = struct;
seed_all = struct;
labnms_all = labnms;
for i = 1:length(labnms)
    taudata_all.(labnms{i}) = data426.(labnms{i});
    tpts_all.(labnms{i}) = tpts.(labnms{i});
    seed_all.(labnms{i}) = seed426.(labnms{i});
end
clearvars -except C_ret C_ant C_nd taudata_all tpts_all labnms_all seed_all netlabs matdir

load([matdir filesep 'KaufmanDiamond_datasets_dat&seed.mat'],'tpts','data426','seed426','labnms');
labnms_all = [labnms_all.' labnms];
for i = 1:length(labnms)
    taudata_all.(labnms{i}) = data426.(labnms{i});
    tpts_all.(labnms{i}) = tpts.(labnms{i});
    seed_all.(labnms{i}) = seed426.(labnms{i});
end
clearvars -except C_ret C_ant C_nd taudata_all tpts_all labnms_all seed_all netlabs matdir

%% 1. Graph metrics
