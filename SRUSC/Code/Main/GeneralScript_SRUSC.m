% Before running script, data matrix X should be stored in workspace and
% the structures SRUSCopts, DenoisingOpts, SpectralOpts, and ComparisonOpts
% should be defined (this can be done by running SetDefaultParameters.m).
% See the README for more details and documentation of options. If a ground truth file is
% available for accuracy quantification, it should be stored as a column 
% vector named LabelsGT. The script denoises the data and then runs SRUSC
% spectral clustering; number of clusters K and scale sigma are
% automatically selected. The script can run K-means and Euclidean spectral
% clustering for comparison if desired.

profile off;
profile on;

if size(X,2)==2
    figure
    scatter(X(:,1),X(:,2),'filled')
    title('Original Data','fontsize',16)
end

%%  Compute eigenvectors of SRUSC Laplacian, while denoising
tic;
if SpatialReg.UseReg
    [EigVals_SRUSC,EigVecs_SRUSC,Idx_Retain,Sigma_SRUSC] = FastEigensolverDenoisingS(X,SRUSCopts,SpectralOpts,DenoisingOpts,SpatialReg);
    Time_USC=toc;
else
    [EigVals_SRUSC,EigVecs_SRUSC,Idx_Retain,Sigma_SRUSC] = FastEigensolverDenoisingS(X,SRUSCopts,SpectralOpts,DenoisingOpts,SpatialReg);
     Time_USC=toc;
end

if ComparisonOpts.RunEucSC
    tic;
    [EigVals_Euc,EigVecs_Euc, Sigma_Euc] = EuclideanEigensolver(X(Idx_Retain,:),SRUSCopts,ComparisonOpts);
    Time_SC=toc;
end 

%% Relabel the ground truth

if exist('LabelsGT')
    LabelsGT=UniqueGT(LabelsGT); 
    if(ismember([0],LabelsGT))
        K = length(unique(LabelsGT(Idx_Retain)))-1; % 0 is not a class!
    else
        K=length(unique(LabelsGT(Idx_Retain)));
    end
end

%%  Compare K-means, Euclidean spectral clustering, and SRUSC spectral clustering on the denoised data

[K_EstimateSRUSC,MaxGapScaleSRUSC,SizeMaxGapSRUSC,SigmaIndexMaxGapSRUSC,NewIndexSRUSC]=ComputeEigengaps(EigVals_SRUSC);
Sigma_EstimateSRUSC = Sigma_SRUSC(SigmaIndexMaxGapSRUSC);

if ComparisonOpts.RunEucSC ==1
    tic;
    [K_EstimateEuc,MaxGapScaleEuc,SizeMaxGapEuc,SigmaIndexMaxGapEuc,NewIndexEuc]=ComputeEigengaps(EigVals_Euc);
    Sigma_EstimateEuc = Sigma_Euc(SigmaIndexMaxGapEuc);
    Time_SC=Time_SC+toc;
end

if exist('LabelsGT')==0
    K_SRUSC = K_EstimateSRUSC;   
    SigmaIndexSRUSC = SigmaIndexMaxGapSRUSC;
    if ComparisonOpts.RunEucSC
        K_Euc = K_EstimateEuc;
        SigmaIndexEuc = SigmaIndexMaxGapEuc;
    end
elseif exist('LabelsGT')
    tic;
    K_SRUSC = K;
    try
        SigmaIndexSRUSC = find(EigVals_SRUSC(K_EstimateSRUSC+1,:)-EigVals_SRUSC(K_EstimateSRUSC,:)==max(EigVals_SRUSC(K_EstimateSRUSC+1,:)-EigVals_SRUSC(K_EstimateSRUSC,:)));
    catch 
        keyboard
    end
    Time_USC=Time_USC + toc;
        if ComparisonOpts.RunEucSC
            tic;
        K_Euc = K;
        SigmaIndexEuc = find(EigVals_Euc(K_EstimateEuc+1,:)-EigVals_Euc(K_EstimateEuc,:)==max(EigVals_Euc(K_EstimateEuc+1,:)-EigVals_Euc(K_EstimateEuc,:)));
        Time_SC=Time_SC+toc;
        end
    
end

tic;
if SpectralOpts.RowNormalization==0
    Labels_SRUSC=kmeans(real(EigVecs_SRUSC(:,1:K_SRUSC,SigmaIndexSRUSC)),K_SRUSC,'Replicates',SpectralOpts.NumReplicates);
elseif SpectralOpts.RowNormalization==1
    Labels_SRUSC=kmeans(normr(real(EigVecs_SRUSC(:,1:K_SRUSC,SigmaIndexSRUSC))),K_SRUSC,'Replicates',SpectralOpts.NumReplicates);
end
Labels_SRUSC_FullData = zeros(size(X,1),1);
Labels_SRUSC_FullData(Idx_Retain) = Labels_SRUSC;
Time_USC=Time_USC+toc;


tic;
if ComparisonOpts.RunEucSC==1
   if SpectralOpts.RowNormalization==0 
       Labels_EuclideanSC=kmeans(real(EigVecs_Euc(:,1:K_Euc,SigmaIndexEuc)),K_Euc,'Replicates',ComparisonOpts.EucSCNumReplicates); 
   elseif SpectralOpts.RowNormalization==1
       Labels_EuclideanSC=kmeans(normr(real(EigVecs_Euc(:,1:K_Euc,SigmaIndexEuc))),K_Euc,'Replicates',ComparisonOpts.EucSCNumReplicates);
   end
   Labels_EuclideanSC_FullData = zeros(size(X,1),1);
   Labels_EuclideanSC_FullData(Idx_Retain) = Labels_EuclideanSC; 
   Time_SC=Time_SC+toc;
   if size(X,2)==2
       figure
       scatter(X(Idx_Retain,1),X(Idx_Retain,2),[],Labels_EuclideanSC,'filled')
       title('Data Labeled by Euclidean SC','fontsize',16)
   end
end

if ComparisonOpts.RunKmeans==1
    [Labels_Kmeans,Centroids_Kmeans]=kmeans(X(Idx_Retain,:),K_SRUSC,'Replicates',SpectralOpts.NumReplicates);
    Labels_Kmeans_FullData = zeros(size(X,1),1);
    Labels_Kmeans_FullData(Idx_Retain) = Labels_Kmeans;
    if size(X,2)==2
        figure
        scatter(X(Idx_Retain,1),X(Idx_Retain,2),[],Labels_Kmeans,'filled')
        title('Data Labeled by K-means','fontsize',16)
    end    
end
%%
tic;
if MajorV.Use ==1
    Labels_SRUSC_FullData_NS=Labels_SRUSC_FullData;
    Labels_SRUSC_FullData=MajorVote(X,Labels_SRUSC_FullData,MajorV);
    Time_USC=Time_USC+toc;
    if ComparisonOpts.RunEucSC==1
        tic;    
        Labels_EuclideanSC_FullData=MajorVote(X,Labels_EuclideanSC_FullData,MajorV);
        Time_SC=Time_SC+toc;
    end
end

    
%% Get accuracies if Labels are available

if exist('LabelsGT')
    
    % Define confusion matrix:
     M = confusionmat(LabelsGT(Idx_Retain),Labels_SRUSC);
     ClusterCounts = M(any(M,2),any(M,1));

    LabelsGT_Retained=LabelsGT(Idx_Retain);
    [OA_SRUSC_SC,AA_SRUSC_SC,Kappa_SRUSC_SC]=GetAccuracies(Labels_SRUSC(LabelsGT_Retained>0),LabelsGT_Retained(LabelsGT_Retained>0),K_SRUSC);

    if ComparisonOpts.RunKmeans==1
        [OA_Kmeans,AA_Kmeans,Kappa_Kmeans]=GetAccuracies(Labels_Kmeans(LabelsGT_Retained>0),LabelsGT_Retained(LabelsGT_Retained>0),K_SRUSC);
    end

    if ComparisonOpts.RunEucSC==1
        [OA_EuclideanSC,AA_EuclideanSC,Kappa_EuclideanSC]=GetAccuracies(Labels_EuclideanSC(LabelsGT_Retained>0),LabelsGT_Retained(LabelsGT_Retained>0),K_Euc);
    end
    
end

%% Plot Multiscale Eigenvalues SRUSC:

if length(Sigma_SRUSC) > 1
    if isempty(NewIndexSRUSC)
        figure
        plot(Sigma_SRUSC,real(EigVals_SRUSC(1,:)),'linewidth',2)
        hold on
        xlim([min(Sigma_SRUSC) max(Sigma_SRUSC)])
        ylim([-.1 1.1])
        for i=2:SpectralOpts.NumEig
            plot(Sigma_SRUSC,real(EigVals_SRUSC(i,:)),'linewidth',2);
        end
    else
        figure
        plot(Sigma_SRUSC(NewIndexSRUSC),real(EigVals_SRUSC(1,NewIndexSRUSC)),'linewidth',2)
        hold on
        xlim([min(Sigma_SRUSC(NewIndexSRUSC)) max(Sigma_SRUSC(NewIndexSRUSC))])
        ylim([-.1 1.1])
        for i=2:SpectralOpts.NumEig
            plot(Sigma_SRUSC(NewIndexSRUSC),real(EigVals_SRUSC(i,NewIndexSRUSC)),'linewidth',2);
        end
    end
    title('Multiscale Eigenvalues for SRUSC','fontsize',16)
end

% Plot Multiscale Eigenvalues Euc SC:
if ComparisonOpts.RunEucSC == 1
    if length(Sigma_Euc)>1
        if isempty(NewIndexEuc)
            figure
            plot(Sigma_Euc,EigVals_Euc(1,:),'linewidth',2)
            hold on
            xlim([min(Sigma_Euc) max(Sigma_Euc)])
            ylim([-.1 1.1])
            for i=2:ComparisonOpts.EucSCNumEig
                plot(Sigma_Euc,real(EigVals_Euc(i,:)),'linewidth',2);
            end
        else
            figure
            plot(Sigma_Euc(NewIndexEuc),real(EigVals_Euc(1,NewIndexEuc)),'linewidth',2)
            hold on
            xlim([min(Sigma_Euc(NewIndexEuc)) max(Sigma_Euc(NewIndexEuc))])
            ylim([-.1 1.1])
            for i=2:ComparisonOpts.EucSCNumEig
                plot(Sigma_Euc(NewIndexEuc),real(EigVals_Euc(i,NewIndexEuc)),'linewidth',2);
            end
        end
        title('Multiscale Eigenvalues for Euc SC','fontsize',16)
    end
end
