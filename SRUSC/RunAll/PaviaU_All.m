
addpath(genpath('../../SRUSC'));
clear all;
close all;

load('paviaU.mat');
X=paviaU([300:340],[140:189],:);
load('paviaU_gt.mat');
LabelsGT=double(paviaU_gt([300:340],[140:189]));
X=reshape(X,41*50,103);
Yuse=UniqueGT(LabelsGT(LabelsGT>0));
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
ComparisonOpts.RunEucSC = 1;
ComparisonOpts.RunKmeans = 0;
ComparisonOpts.EucSCSigmaScaling = 'Automatic'; 


SpatialReg.UseReg = 0;
SpatialReg.Width = 41; %41
SpatialReg.Height = 50; %50
SpatialReg.r=30;

MajorV.Use = 0;
MajorV.Radius = 3;
MajorV.VoteRate = 0.7;
MajorV.Width = 41;
MajorV.Height = 50;

GeneralScript_SRUSC;

%% Kmeans
tic;
[Labels_KmeansFull,~]=kmeans(X,K,'Replicates', SpectralOpts.NumReplicates);
Time_Kmeans=toc;
[OA_KmeansFull,AA_KmeansFull,Kappa_KmeansFull]=GetAccuracies(Labels_KmeansFull(LabelsGT>0),Yuse,K);


%% Do PCA dimension reduction first, then compute Kmeans
tic;
V=GetPC(X);
Labels_PCA=kmeans(X*V(:,1:10),K,'Replicates', SpectralOpts.NumReplicates)';
Labels_PCA=reshape(Labels_PCA,41*50,1);
Time_PCA=toc;
[OA_PCA,AA_PCA,Kappa_PCA]=GetAccuracies(Labels_PCA(LabelsGT>0),Yuse,K);


%% Compute Gaussian Mixture Model
tic;
gmfit=fitgmdist(X,K, 'CovarianceType', 'full',  'RegularizationValue', 1e-10);
Labels_GMM= cluster(gmfit,X);
Time_GMM=toc;
[OA_GMM,AA_GMM,Kappa_GMM]=GetAccuracies(Labels_GMM(LabelsGT>0),Yuse,K);


%% Compute Hierarchical NMF
tic;
Labels_NMF=hierclust2nmf((X-min(min(X)))',K);
Time_NMF=toc;
 [OA_NMF,AA_NMF,KappaNMF]=GetAccuracies(Labels_NMF(LabelsGT>0),Yuse,K);


%% Compute Diffusion Learning
DataInput.datasource='HSI';
DataInput.HSIname='Pavia';
[X_d,Y,Ytuple,K_GT,d,DataParameters,PeakOpts.Window_Diffusion]=GenerateData(DataInput);
GT=DataParameters.GT;
M=DataParameters.M;
N=DataParameters.N;

% Set Parameters

%Data name
PeakOpts.Data=DataInput.HSIname;

%Should we plot everything or not?  Should we save things?
UserPlot=0;
SavePlots=0;
PeakOpts.UserPlot=UserPlot;

%Detect number of eigenvectors to use automatically
PeakOpts.DiffusionOpts.K='automatic'; 

%How many nearest neighbors to use for diffusion distance
PeakOpts.DiffusionOpts.kNN=100;

%Force probability of self-loop to exceed .5.
PeakOpts.DiffusionOpts.LazyWalk=0;

%Set epsilon for weight matrix; use same one for spectral clustering
PeakOpts.DiffusionOpts.epsilon=1;

% Set time to evaluate diffusion distance
PeakOpts.DiffusionTime=30;

%Kmeans parameters
NumReplicates=10;
MaxIterations=100;

%Set how many cores we want to learn.  If K_Learned is different from
%K_GT, it will not be possible to compare easily and immediately against
%the ground truth, though other patterns may be studied in this way.
K_Learned=K_GT;

%How many nearest neighbors to use for KDE
DensityNN=20;

%Spatial window size for labeling regime
SpatialWindowSize=3;

tic;
if isempty(GT)
    [I,J]=find(sum(DataParameters.RawHSI,3)~=0);
else
   [I,J]=find(GT>-1);
end
Time_temp=toc;

% Diffusion distances
tic;
PeakOpts.ModeDetection='Diffusion';
[CentersDiffusion, G_dif, DistStructDiffusion] = DensityPeaksEstimation(X_d, K_Learned, DensityNN, PeakOpts, M, N);
Labels_DL=FSFclustering(X,K,DistStructDiffusion,CentersDiffusion);
Labels_DL=reshape(Labels_DL,41*50,1);
Time_DL=Time_temp+toc;

% Euclidean distances
tic;
PeakOpts.ModeDetection='Euclidean';
[CentersEuclidean, G_euc, DistStructEuclidean] = DensityPeaksEstimation(X_d, K_Learned, DensityNN, PeakOpts, M, N);
Labels_FSFDPC=FSFclustering(X,K_GT,DistStructEuclidean,CentersEuclidean);
Labels_FSFDPC=reshape(Labels_FSFDPC,41*50,1);
Time_FSFDPC=Time_temp+toc;

% Diffusion distances with spatial regularization
%Compute OA, AA, Kappa after aligning
%Diffusion Learning
[OA_DL,AA_DL,Kappa_DL]=GetAccuracies(Labels_DL(LabelsGT>0),Yuse,K);
[OA_FSFDPC,AA_FSFDPC,Kappa_FSFDPC]=GetAccuracies(Labels_FSFDPC(LabelsGT>0),Yuse,K_Learned);


%% Compute LCMR

load('paviaU.mat');
[RD_hsi]=paviaU([300:340],[140:189],:);
load('paviaU_gt.mat');
labels=double(paviaU_gt([300:340],[140:189]));
labels(find(labels==6))=1;
labels(find(labels==7))=2;
labels(find(labels==9))=3;
sz = size(RD_hsi);
no_classes = 3;
wnd_sz =10;
K_n = 100;
tic;
train_number = ones(1,no_classes)*5;
[lcmrfea_all] =  fun_LCMR_all(RD_hsi,wnd_sz,K_n);
Time_LCMR=toc;

tic;
for flag = 1:5
        [train_SL,test_SL,test_number]= GenerateSample(labels,train_number,no_classes);
        train_id = train_SL(1,:);
        train_label = train_SL(2,:);
        test_id = test_SL(1,:);
        test_label = test_SL(2,:);
        train_cov = lcmrfea_all(:,:,train_id);
        test_cov = lcmrfea_all;
        KMatrix_Train = logmkernel(train_cov, train_cov);
        KMatrix_Test = logmkernel(train_cov, test_cov);
        Ktrain = [(1:size(KMatrix_Train,1))',KMatrix_Train];     
        model = svmtrain(train_label', Ktrain, '-t 4');  
        Ktest = [(1:size(KMatrix_Test,2))', KMatrix_Test'];  
        tmp = ones(1,size(KMatrix_Test,2));
        [predict_label, accuracy, P1] = svmpredict(tmp',Ktest,model); 
        
end

Labels_LCMR = predict_label;
Time_LCMR=Time_LCMR+toc;
[OA_LCMR,AA_LCMR,KappaLCMR]=GetAccuracies(Labels_LCMR(LabelsGT>0),Yuse,3);

