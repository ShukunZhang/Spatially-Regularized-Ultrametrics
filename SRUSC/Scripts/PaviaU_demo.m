 addpath(genpath('../../SRUSC'));
clear all;
close all;

load('paviaU.mat');
X=paviaU([300:340],[140:189],:);
load('paviaU_gt.mat');
LabelsGT=double(paviaU_gt([300:340],[140:189]));
X=reshape(X,41*50,103);

%% Compute Ultrametric SC with LLPD 
% The Kmeans and SC with Euc distance will be computed among our method

SetDefaultParameters
SRUSCopts.KNN = 600 ; %Number of nearest neighbors in underlying graph
SRUSCopts.UseFixedNumScales = 0; %Forces a fixed number of scales in LLPD approximation
SRUSCopts.LogRatio=1.1;
DenoisingOpts.Method='Cutoff'; %Look at figure to denoise, or set a cutoff?
DenoisingOpts.Cutoff = 750;
SpectralOpts.SigmaScaling = 'Manual';
SpectralOpts.SigmaValues = 400;
SpectralOpts.NumEig = 10;
ComparisonOpts.RunEucSC = 0;
ComparisonOpts.EucSCSigmaScaling = 'Manual'; 
ComparisonOpts.EucSCSigmaValues = 400;
ComparisonOpts.RunKmeans = 0;


SpatialReg.UseReg = 1;
SpatialReg.Width = 41; %41
SpatialReg.Height = 50; %50
SpatialReg.r=30;

MajorV.Use = 1;
MajorV.Radius = 3;
MajorV.VoteRate = 0.7;
MajorV.Width = 41;
MajorV.Height = 50;

 
  GeneralScript_SRUSC

%figure; imagesc(reshape(Labels_SRUSC_FullData,41,50));