% Toolbox Name: main_SPICE_DMRSI_processing
%
% Description:
%   The toolbox provides a collection of functions and tools for analyzing, 
%   processing, and visualizing dynamic deuterium magnnatic resonance spectroscopic
%   imaging (DMRSI) signals.
%
% Key Features:
%   - B0 field correction
%   - Spectral profiles and temporal bases derivation
%   - DMRSI denoising and decomposition
%
% Target Audience:
%   This toolbox is intended for researchers, engineers, and students working in the field 
%   of deuterium MRSI of human brain and rodent brain, and related areas.
%
% Requirements:
%   - MATLAB R2014a or later
%   - Signal Processing Toolbox (required for some functions)
%   - SPICE Toolbox
%
% Author:
%   Liang's Lab @ UIUC
%
% Data:
%   2025-02-22
%
% License:
%   This toolbox is distributed under the MIT License. See LICENSE.txt for details.
%
% Version:
%   v1.0.0 (Feb 2025)
%   - Initial release

clear all;
restoredefaultpath;
set(0,'defaultfigureColor','w');
defaultFigPos = [100,100,800,600];
set(0,'defaultfigureposition',defaultFigPos);

%% Set Path
type = 'WholeBrain';
dataID = 'HC031';                              %%%% 'TP006'                       or 'HC031'
dataName = 'fid_combined5D_42measWSVD.mat';  %%%% 'fid_combined5D_34measWSVD_2' or 'fid_combined5D_42measWSVD'
dataVariableName = 'fid_combined';            %%%% 'fid_combined2'               or'fid_combined'
% data path
homePath     = './';
rawDataPath  = fullfile(homePath,'raw_data','whole_brain','WhiteningSVD', dataID, filesep);
processPath  = fullfile(homePath,'processed_data','whole_brain','WhiteningSVD', dataID, filesep);
SaveDataPath = fullfile(processPath,'data',filesep);
if ~exist(SaveDataPath,'dir')
    mkdir(SaveDataPath);
end
SaveFigsPath = fullfile(processPath,'figs',filesep);
if ~exist(SaveFigsPath,'dir')
    mkdir(SaveFigsPath);
end
% library path
libPath      = fullfile(homePath,'support',filesep);
addpath(genpath(libPath));
% pretrained profile path
pretrainedDataPath         = fullfile(homePath,'pretrained_data',filesep);
pretrainedDataName         = 'Group_Temporal_Subspace.mat';
pretrainedDataVariableName = 'Casorati';
disp('Path setup complete');

%% load data
% initial frame
ini_frame = 1:4;
% throw bad data
num_throw = 0;
% usable frames
use_frame = 1:20;
% load data
sxt_dyn_highres   = myLoadVars([rawDataPath,dataName], dataVariableName);
sxt_dyn_highres   = permute(sxt_dyn_highres,[2,3,4,1,5]);
sxt_dyn_highres   = sxt_dyn_highres(:,:,:,1+num_throw:end,use_frame);
sxt_dyn_highres   = conj(sxt_dyn_highres);
[L1_h,L2_h,L3_h,M,Nt_h] = size(sxt_dyn_highres);
Nt = length(use_frame);

%% Set Parameters
% sequence params
dt            = 8.5143e-4;
% acquisition parameters: 
% global constants
gyro          = 6.536;   % gyromagnetic ratio
B0            = 7;       % field length
f_H2          = gyro*B0; % frequency at targeting scanner [MHz]
% chemical shift parameters
shift_f0      = 8.0270;
shift_d2o     = 4.8;
shift_glu     = 3.773;
shift_glx2    = 2.306;
shift_lac     = 1.4;
shift_meta    = {shift_d2o, shift_glu, shift_glx2, shift_lac};
% generate ppm vector
ppm_vec       = gen_ppm_vec(M,dt,f_H2*1e6,shift_f0);

%% ------------------------ Step 0: Preparation
%% prepare the QM physics-based bases
meta       = {'D2O','Glc','Glx','Lac'};
num_meta   = length(shift_meta);
freq_meta  = zeros(num_meta,1);
for i = 1:num_meta
    freq_meta(i) = (shift_meta{i} - shift_f0)*f_H2;
end
tvec_D2    = (0:M-1).'*dt;
bases      = exp(1i*2*pi*freq_meta*tvec_D2.').';

%% brainMask
brainMask = normRange(sumOfSqr(mean(sxt_dyn_highres(:,:,:,:,1:end),5)))>0.1;
figure;montagesc(normRange(sumOfSqr(mean(sxt_dyn_highres(:,:,:,:,1:end),5))));axis off;
figure;montagesc(brainMask);axis off;

%% check data (high-res)
vxl  = [10,10,7];

figure;plot(ppm_vec,abs(F1_t2f(vec(mean(sxt_dyn_highres(vxl(1),vxl(2),vxl(3),:,:),5)))));
MakeFigPretty(gcf,2,2); xlim([-5,20]); set(gca,'XDir','reverse');
mySaveFigPng(gcf,fullfile(SaveFigsPath,'data_checking'),'highres_avg_spectrum');

figure;plot(ppm_vec,abs(F1_t2f(squeeze(sxt_dyn_highres(vxl(1),vxl(2),vxl(3),:,:)))));
MakeFigPretty(gcf,2,2); xlim([-5,20]); set(gca,'XDir','reverse');
mySaveFigPng(gcf,fullfile(SaveFigsPath,'data_checking'),'highres_indv_spectrum');

%% collect pipeline information
DMRSI_info.SaveDataPath = SaveDataPath;
DMRSI_info.SaveFigsPath = SaveFigsPath;
DMRSI_info.shift_meta   = shift_meta;
DMRSI_info.brainMask    = brainMask;
DMRSI_info.ini_frame    = ini_frame;
DMRSI_info.shift_f0     = shift_f0;
DMRSI_info.num_meta     = num_meta;
DMRSI_info.tvec_D2      = tvec_D2;
DMRSI_info.ppm_vec      = ppm_vec;
DMRSI_info.bases        = bases;
DMRSI_info.meta         = meta;
DMRSI_info.f_H2         = f_H2;
DMRSI_info.dt           = dt;

%% collect pretrained file information
PretrainedInfo.pretrainedDataPath         = pretrainedDataPath;
PretrainedInfo.pretrainedDataName         = pretrainedDataName;
PretrainedInfo.pretrainedDataVariableName = pretrainedDataVariableName;

%% ------------------------ Step 1: field correction
opt = struct('flag_disp', true, 'vxl', vxl);
sxt_dyn_highres_corr = SPICE_dynDMRSI_field_correction(sxt_dyn_highres, DMRSI_info, opt);
mySaveVars([SaveDataPath,'sxt_zp_corr.mat'],'sxt_dyn_highres_corr','sxt_dyn_highres');
% rename the original dmrsi data
sxt_orig = sxt_dyn_highres_corr; clear sxt_dyn_highres_corr;

%% ------------------------ Step 2: spectral profiles and temporal bases derivation
%% obtain spectral profiles
opt = struct('flag_disp',true,'Ns_tmp',[2 2 1 1]);
sxt_quant_profile = SPICE_dynDMRSI_spectProfile_derive(sxt_orig, DMRSI_info, opt);
mySaveVars([SaveDataPath,'SpectralTemporalBases.mat'],'sxt_quant_profile');

%% obtain temporal bases (data-driven)
opt = struct('flag_disp',false,'type',type);
[Vt_temporal, ~] = SPICE_dynDMRSI_temporalBase_derive(sxt_orig,sxt_quant_profile,DMRSI_info,PretrainedInfo,opt);
mySaveVars([SaveDataPath,'SpectralTemporalBases.mat'],'Vt_temporal');

%% ------------------------ Step 3: dynamic DMRSI denoising and decomposition
opt = struct('flag_disp',true,'flag_discrepancy_principle',false,...
             'Rank_temporal',[4,4,2,2], 'Ns_tmp',[2,2], ...
             'ker_size',[3,3,3], 'localV_rank',5);
[sxt_denoised, sxt_decomp, sxt_decomp_orig] = SPICE_dynDMRSI_denoising_decomposition(sxt_orig,sxt_quant_profile,Vt_temporal,opt);

%% ------------------------ Step 4: Save the processed results

mySaveVars(fullfile(SaveDataPath,'processedData_XinCMBD_WSVD.mat'),...
           'sxt_orig',...
           'sxt_denoised',...
           'sxt_decomp',...
           'sxt_decomp_orig');




