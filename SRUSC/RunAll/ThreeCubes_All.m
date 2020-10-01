close all;
clear all;
num=20;
X=zeros(num^3,3);

for i=1:num
    for j=1:num
        for k=1:num
            X((i-1)*num^2+(j-1)*num+k,1) = i/num;
            X((i-1)*num^2+(j-1)*num+k,2) = j/num;
            X((i-1)*num^2+(j-1)*num+k,3) = k/num;
        end
    end
end
X=vertcat(X,X,X);
scatter3(X(:,1),X(:,2),X(:,3),[],'filled');

X([num^3+1:2*(num^3)],3) = X([num^3+1:2*(num^3)],3)+1.5;
X([2*num^3+1:3* num^3],3) = X([2*num^3+1:3* num^3],3)+3;

scatter3(X(:,1),X(:,2),X(:,3),[],'filled');
X=horzcat(X,zeros(3*num^3,196));
 
[Q,R] = qr(randn(199));
X=X*Q;
X=horzcat(X,ones(3*num^3,1));
for i=1:3
     for j=1:num^3
         X((i-1)*num^3 + j,200) = i/10; 
     end
 end


C_1=randi([num^3/2 num^3/2+100],10,1);
C_2=randi([2*num^3+num^3/2 2*num^3+num^3/2+100],10,1);
  
temp =  X(C_1,:);
X(C_1,:) = X(C_2,:);
X(C_2,:) = temp;

addpath(genpath('../../SRUSC'));
load('TC_GT_100x240.mat');

SetDefaultParameters 
DenoisingOpts.Method='None';
 
 %using spatial 
SRUSCopts.KNN = 8100; 
SpectralOpts.SigmaScaling = 'Manual';
SpectralOpts.NumEig = 10;
SpectralOpts.SigmaValues = linspace(1,10,5);



SpatialReg.UseReg = 1;
SpatialReg.Width = 240;
SpatialReg.Height = 100;

SpatialReg.r=50;
MajorV.Use = 0;

 
ComparisonOpts.RunEucSC = 0;
ComparisonOpts.Kmeans = 0;
ComparisonOpts.EucSCSigmaScaling = 'Automatic'; 
% 
% 
% 
  GeneralScript_SRUSC

%% Kmeans
tic;
[Labels_KmeansFull,~]=kmeans(X,K,'Replicates', SpectralOpts.NumReplicates);
Time_Kmeans=toc;
[OA_KmeansFull,AA_KmeansFull,Kappa_KmeansFull]=GetAccuracies(Labels_KmeansFull(LabelsGT>0),LabelsGT(LabelsGT>0),K);


 %% Do PCA dimension reduction first, then compute Kmeans
tic;
V=GetPC(X);
Labels_PCA=kmeans(X*V(:,1:2),K,'Replicates', SpectralOpts.NumReplicates)';
Labels_PCA=reshape(Labels_PCA,24000,1);
Time_PCA=toc;
[OA_PCA,AA_PCA,Kappa_PCA]=GetAccuracies(Labels_PCA(LabelsGT>0),reshape(LabelsGT,24000,1),K);


%% Compute Gaussian Mixture Model
tic;
gmfit=fitgmdist(X,K, 'CovarianceType', 'full',  'RegularizationValue', 1e-10);
Labels_GMM= cluster(gmfit,X);
Time_GMM=toc;
[OA_GMM,AA_GMM,Kappa_GMM]=GetAccuracies(Labels_GMM(LabelsGT>0),reshape(LabelsGT,24000,1),K);


%% Compute Hierarchical NMF
tic;
Labels_NMF=hierclust2nmf((X-min(min(X)))',K);
Time_NMF=toc;
[OA_NMF,AA_NMF,KappaNMF]=GetAccuracies(Labels_NMF(LabelsGT>0),reshape(LabelsGT,24000,1),K);



%% Compute Diffusion Learning
Y=LabelsGT;
Yuse=UniqueGT(LabelsGT(LabelsGT>0));
GT=Y;
Y=reshape(Y,size(Y,1)*size(Y,2),size(Y,3));
Y=Y';

Y=Y(Y>0);
Y=UniqueGT(Y);
    
K_GT=max(Y);
    
Y=LabelsGT';  
%X=permute(X,[2,1,3]); 
M=100;
N=240;
PeakOpts.Window_Diffusion = 20;


%%
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
K_Learned=K;

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
[CentersDiffusion, G_dif, DistStructDiffusion] = DensityPeaksEstimation(X, K_Learned, DensityNN, PeakOpts, M, N);
Labels_DL=FSFclustering(X,K,DistStructDiffusion,CentersDiffusion);
Labels_DL=reshape(Labels_DL,24000,1);
Time_DL=Time_temp+toc;
% 

% 
% % Euclidean distances
tic;
PeakOpts.ModeDetection='Euclidean';
[CentersEuclidean, G_euc, DistStructEuclidean] = DensityPeaksEstimation(X, K_Learned, DensityNN, PeakOpts, M, N);
Labels_FSFDPC=FSFclustering(X,K_GT,DistStructEuclidean,CentersEuclidean);
Labels_FSFDPC=reshape(Labels_FSFDPC,24000,1);
Time_FSFDPC=Time_temp+toc;
% 
% % Diffusion distances with spatial regularization
% %Compute OA, AA, Kappa after aligning
% %Diffusion Learning
[OA_DL,AA_DL,Kappa_DL]=GetAccuracies(Labels_DL(LabelsGT>0),Yuse,K);
[OA_FSFDPC,AA_FSFDPC,Kappa_FSFDPC]=GetAccuracies(Labels_FSFDPC(LabelsGT>0),Yuse,K_Learned);

% 
%

%% Compute LCMR

[RD_hsi]= reshape(X,240,100,200);

load('ThreeCubes_gt.mat');
labels=TCGT;
sz = size(RD_hsi);
no_classes = 3;
wnd_sz =10;
K_n = 100;
tic;
train_number = ones(1,no_classes)*5;
[lcmrfea_all] =  fun_LCMR_all(RD_hsi,wnd_sz,K_n);
     Time_LCMR=toc;   

%%
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
[OA_LCMR,AA_LCMR,KappaLCMR]=GetAccuracies(Labels_LCMR(LabelsGT>0),Yuse,3);
 
