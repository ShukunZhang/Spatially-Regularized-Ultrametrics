addpath(genpath('../../SRUSC'));

clear all;
close all;

load('SalinasA_gt.mat');
LabelsGT=salinasA_gt;
X=load('SalinasA.mat');
X=X.salinasA;
X=reshape(X,83*86,224);

%% Compute Ultrametric SC with SRUSC 
% The Kmeans and SC with Euc distance will be computed among our method
SetDefaultParameters

SRUSCopts.KNN = 1500; %Number of nearest neighbors in underlying graph
SRUSCopts.UseFixedNumScales = 0; %Forces a fixed number of scales in SRUSC approximation
SRUSCopts.LogRatio=1.1;

DenoisingOpts.Method='Cutoff'; %Look at figure to denoise, or set a cutoff?
DenoisingOpts.Cutoff = 290;
ComparisonOpts.RunEucSC = 0;
SpectralOpts.SigmaScaling = 'Manual';
SpectralOpts.SigmaValues = 200;
SpectralOpts.NumEig = 20;

SpatialReg.UseReg = 1;
SpatialReg.Width = 83;
SpatialReg.Height = 86;
SpatialReg.r=65;

MajorV.Use = 1;
MajorV.Radius = 6;
MajorV.VoteRate = 0.7;
MajorV.Width = 83;
MajorV.Height = 86;

 
 GeneralScript_SRUSC

figure; imagesc(reshape(Labels_SRUSC_FullData,83,86));
