clear all;
close all;
n=4900; % number of points that you want

 radius = 1.7; % radius of the circle
 angle1 = 2*pi*rand(n,1);
 angle2 = 2*pi*rand(n,1);
 A=rand(n,2);
 B=rand(n,2);
 C=rand(n,2);
 D=rand(n,2);

 for i=1:99
    center1 = [1,3];
    r1 = radius*sqrt(rand(n,1));
    X_1 = r1.*cos(angle1)+ center1(1);
    Y_1 = r1.*sin(angle1)+ center1(2);

    center2 = [1,5];
    r2 = radius*sqrt(rand(n,1));
    X_2 = r2.*cos(angle1)+ center2(1);
    Y_2 = r2.*sin(angle1)+ center2(2);
 
    center3 = [1,7];
    r3 = radius*sqrt(rand(n,1));
    X_3 = r3.*cos(angle1)+ center3(1);
    Y_3 = r3.*sin(angle1)+ center3(2);
    
    
    center4 = [5,5];
    r4 = radius*sqrt(rand(n,1));
    X_4 = r4.*cos(angle2)+ center4(1);
    Y_4 = r4.*sin(angle2)+ center4(2);
    
    
    A=horzcat(A,X_1,Y_1);
    B=horzcat(B,X_2,Y_2);
    C=horzcat(C,X_3,Y_3);
    D=horzcat(D,X_4,Y_4);
    
end

X=vertcat(A,B,C,D);
addpath(genpath('../../SRUSC'));
 
load('FS_GT_140x140.mat');
%if you want to change the number of points in each sphere, you need to 
%construct the ground truth by yourself



 SetDefaultParameters
 
 DenoisingOpts.Method='None';
 LLPDopts.KNN = 5000; %Number of nearest neighbors in underlying graph
 SpectralOpts.SigmaScaling = 'Manual';
 SpectralOpts.SigmaValues = 8;
 
 SpatialReg.UseReg = 0;
 SpatialReg.Width = 140;
 SpatialReg.Height = 140;

 SpatialReg.r=65;
 MajorV.Use = 0;

 
 ComparisonOpts.RunEucSC = 1;
 ComparisonOpts.Kmeans = 0;
 ComparisonOpts.EucSCSigmaScaling = 'Manual';
 ComparisonOpts.EucSCSigmaValues = 8;

 
 GeneralScript_SRUSC

%% Kmeans
tic;
[Labels_KmeansFull,~]=kmeans(X,K,'Replicates', SpectralOpts.NumReplicates);
Time_Kmeans=toc;
[OA_KmeansFull,AA_KmeansFull,Kappa_KmeansFull]=GetAccuracies(Labels_KmeansFull(LabelsGT>0),LabelsGT(LabelsGT>0),K);

imagesc(reshape(Labels_KmeansFull,140,140));
%% Do PCA dimension reduction first, then compute Kmeans
tic;
V=GetPC(X);
Labels_PCA=kmeans(X*V(:,1:2),K,'Replicates', SpectralOpts.NumReplicates)';
Labels_PCA=reshape(Labels_PCA,19600,1);
Time_PCA=toc;
[OA_PCA,AA_PCA,Kappa_PCA]=GetAccuracies(Labels_PCA(LabelsGT>0),LabelsGT(LabelsGT>0),K);
imagesc(reshape(Labels_PCA,140,140));

%% Compute Gaussian Mixture Model
tic;
gmfit=fitgmdist(X,2, 'CovarianceType', 'full',  'RegularizationValue', 1e-10);
Labels_GMM= cluster(gmfit,X);
Time_GMM=toc;
[OA_GMM,AA_GMM,Kappa_GMM]=GetAccuracies(Labels_GMM(LabelsGT>0),LabelsGT(LabelsGT>0),K);
imagesc(reshape(Labels_GMM,140,140));

%% Compute Hierarchical NMF
tic;
Labels_NMF=hierclust2nmf((X-min(min(X)))',K);
 [OA_NMF,AA_NMF,KappaNMF]=GetAccuracies(Labels_NMF(LabelsGT>0),LabelsGT(LabelsGT>0),K);
 Time_NMF=toc;

imagesc(reshape(Labels_NMF,140,140));
%% 
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
M=140;
N=140;
PeakOpts.Window_Diffusion = 20;


%%
%Should we plot everything or not?  Should we save things?
UserPlot=0;
SavePlots=0;
PeakOpts.UserPlot=UserPlot;

%Detect number of eigenvectors to use automatically
PeakOpts.DiffusionOpts.K='automatic'; 

%How many nearest neighbors to use for diffusion distance
PeakOpts.DiffusionOpts.kNN=1000;

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
Labels_DL=reshape(Labels_DL,19600,1);
Time_DL=Time_temp+toc;
[OA_DL,AA_DL,Kappa_DL]=GetAccuracies(Labels_DL(LabelsGT>0),Yuse,K);
imagesc(reshape(Labels_DL,140,140));


%% Euclidean distances
tic;
 PeakOpts.ModeDetection='Euclidean';
 [CentersEuclidean, G_euc, DistStructEuclidean] = DensityPeaksEstimation(X, K_Learned, DensityNN, PeakOpts, M, N);
 Labels_FSFDPC=FSFclustering(X,K_GT,DistStructEuclidean,CentersEuclidean);
 Labels_FSFDPC=reshape(Labels_FSFDPC,19600,1);
 [OA_FSFDPC,AA_FSFDPC,Kappa_FSFDPC]=GetAccuracies(Labels_FSFDPC(LabelsGT>0),Yuse,K_Learned);
 Time_FSFDPC=Time_temp+toc;
 imagesc(reshape(Labels_FSFDPC,140,140));
% 
% 
% %Compute OA, AA, Kappa after aligning
% %Diffusion Learning
 
%%
tic;
Y=reshape(X,140,140,200);
gt=LabelsGT;
row=size(Y,1);
col=size(Y,2);
nband=size(Y,3);
Y=reshape(Y,row*col,nband);
lambda=1e-8;
alpha=0.001;
rho=1;
thr=[2*10^-4,2*10^-4];
maxIter=100;
n=10;

C=L2LRC(Y,lambda,alpha,rho,maxIter,nband);
CKSym=C.*C/max(max(C));
groups = SC(CKSym,n);
CC=zeros(row*col,n);
CC1=zeros(row*col,n);
for i=1:n
    R=find(groups(:)==i);
    s=zeros(row*col,1);
    for j=1:length(R)
        s=s+Y(:,R(j));
    end
    CC1(:,i)=s./length(R);
    dist=[];
    for k=1:length(R)
        dist(k)=norm((Y(:,R(k))-CC1(:,i)),2);
    end
    center=min(dist);
    R1=find(dist(:)==center);
    CC(:,i)=Y(:,R(R1(1)));
end
datasubset=reshape(CC,row,col,n);
band=reshape(datasubset,row*col,n);
[Labels_LLR,~]=kmeans(band,2,'Replicates', 10);
Time_LLR=toc;
imagesc(reshape(Labels_LLR,140,140));
axis off;
axis equal;
%%
[OA_LLR,AA_LLR,Kappa_LLR]=GetAccuracies(Labels_LLR ,Yuse,2);
L=figure; imagesc(reshape(Labels_LLR,140,140));
%% Compute LCMR
Yuse=UniqueGT(LabelsGT(LabelsGT>0));
[RD_hsi]= reshape(X,140,140,200);


labels=LabelsGT;
sz = size(RD_hsi);
no_classes = 2;
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
[OA_LCMR,AA_LCMR,KappaLCMR]=GetAccuracies(Labels_LCMR(LabelsGT>0),Yuse,2);
