clear all;
close all;
addpath(genpath('W:\6_SEEG_Bandit\1_Analysis_banditOnline\2_ANALYSIS_VKF\cbm-master\codes'));
fdata = load('alldatabandit.mat');
data = fdata.alldata;
% 2nd input: a cell input containing function handle to models
models = {@fit_kf,@fit_vkf};
% note that by handle, I mean @ before the name of the function
% 3rd input: another cell input containing file-address to files saved by cbm_lap


v=6.25;



%kf with loss/win sigma/omega
prior_kf = struct('mean',zeros(3,1),'variance',v);
fname_kf= 'lap_kf_v6p25.mat';
cbm_lap(data,@fit_kf, prior_kf,fname_kf);

%vkf
prior_vkf = struct('mean',zeros(4,1),'variance',v);
fname_vkf = 'lap_vkf_v6p25.mat';
cbm_lap(data,@fit_vkf, prior_vkf,fname_vkf);