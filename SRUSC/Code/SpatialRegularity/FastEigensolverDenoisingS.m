function [EigVals,FullEigVecs,IdxRetain,Sigma] = FastEigensolverDenoisingS(X,SRUSCopts,SpectralOpts,DenoisingOpts,SpatialReg)

% Input: 
%
% -X: nxD data matrix
% -SRUSCopts: Structure of options for how to approximate LLPD
% -SpectralOpts: Structure of options for LLPD spectral clustering
% -DenoisingOpts: Structure of options for how to denoise with Knn LLPD
% distances

%
% Output:
%
% -EigVals: Eigenvalues of LLPD graph Laplacian
% -EigVecs: Eigenvectors of LLPD graph Laplacian
% -Sigma: Range of sigma scaling values used. 

%% Compute the connected components on the non-denoised data, then denoise the Data


[SortedCCmatrixInitial, ThresholdsInitial] = ConnectedComponents(X, SRUSCopts);
Tree = ConstructTree(SortedCCmatrixInitial, ThresholdsInitial, DenoisingOpts.KNN);

p = size(Tree.CC_array,1); %Number of CC's at finest scale
n=size(X,1);
ydata = zeros(1, n);
index = 0;

for i=1:p
    ydata(index+1:index+Tree.SS(i)) = Tree.NNDistances(i);
    index = index + Tree.SS(i);
end

if strcmp(DenoisingOpts.Method, 'None')
    cutoff=1.01*max(Tree.NNDistances(1:p));
    CCsToKeepIndex = find(Tree.NNDistances(1:p) < cutoff);
    RowsToKeepIndex = find(ismember(Tree.PointPartitionAtFinestScale(:,1),CCsToKeepIndex)==1);
    IdxRetain = Tree.PointPartitionAtFinestScale(RowsToKeepIndex,2);
    denoisedX = X(IdxRetain,:);

elseif strcmp(DenoisingOpts.Method, 'Cutoff')
    cutoff=DenoisingOpts.Cutoff;
    CCsToKeepIndex = find(Tree.NNDistances(1:p) < cutoff);
    RowsToKeepIndex = find(ismember(Tree.PointPartitionAtFinestScale(:,1),CCsToKeepIndex)==1);
    IdxRetain = Tree.PointPartitionAtFinestScale(RowsToKeepIndex,2);
    denoisedX = X(IdxRetain,:);
    
    figure
    scatter(1:n, sort(ydata),'filled')
    hold on
    scatter(1:n, cutoff*ones(1,n), 'x')
    xlim([0 n]);
    title('Plot of SRUSC to NN and Denoising Cutoff','fontsize',16)
    xlabel('Point Index', 'fontsize',14)
    ylabel(['Path Distance to K=' num2str(DenoisingOpts.KNN) ' th Nearest Neighbor'], 'fontsize',14)
    
elseif strcmp(DenoisingOpts.Method,'ExamineFigure')
    close all;
    figure
    scatter(1:n, sort(ydata),'filled')
    xlim([0 n]);
    title('Plot of LLPD to NN','fontsize',16)
    xlabel('Point Index', 'fontsize',14)
    ylabel(['Path Distance to K=' num2str(DenoisingOpts.KNN) ' th Nearest Neighbor'], 'fontsize',14)
    Prompt='At what Knn-LLPD value to threshold for denoising? \n';
    
    cutoff=input(Prompt);
    CCsToKeepIndex = find(Tree.NNDistances(1:p) < cutoff);
    if isempty(CCsToKeepIndex)
        error('Error: all points cannot be removed during denoising. Please increase cutoff.')
        return
    end
    RowsToKeepIndex = find(ismember(Tree.PointPartitionAtFinestScale(:,1),CCsToKeepIndex)==1);
    IdxRetain = Tree.PointPartitionAtFinestScale(RowsToKeepIndex,2);
    denoisedX = X(IdxRetain,:);
    
elseif strcmp(DenoisingOpts.Method, 'MultipleCutoffs')
    cutoff_ranges = DenoisingOpts.multiple_cutoffs;
    m = length(cutoff_ranges);
    CCsToKeepIndex = [];
    for i=1:p
        throw_away = 0; %don't remove the CC unless it falls into one of the specified ranges
        for j=0:(m-3)/2
            if Tree.NNDistances(i)>cutoff_ranges(2*j+1) && Tree.NNDistances(i)<cutoff_ranges(2*j+2)
                throw_away = 1;
            end
        end
        if Tree.NNDistances(i) > cutoff_ranges(end)
            throw_away = 1;
        end
        if throw_away == 0
        CCsToKeepIndex = vertcat(CCsToKeepIndex,i);
        end
    end
    RowsToKeepIndex = find(ismember(Tree.PointPartitionAtFinestScale(:,1),CCsToKeepIndex)==1);
    IdxRetain = Tree.PointPartitionAtFinestScale(RowsToKeepIndex,2);
    denoisedX = X(IdxRetain,:);
end

%% Recompute Connected Components Matrix on the De-noised Data

% try different sigmas, LLPD different
    NumEig = SpectralOpts.NumEig;
    if strcmp(DenoisingOpts.Method, 'None')
        sortedCCmatrix = SortedCCmatrixInitial;
        th = ThresholdsInitial;
    else
        [sortedCCmatrix, th] = ConnectedComponents(denoisedX, SRUSCopts);
    end

    LLPD = knnLLPD(sortedCCmatrix,th,SRUSCopts.KNN);
    
    S=ComputeSpatial_Fast(SpatialReg);
    S=S(IdxRetain,IdxRetain);

    if strcmp(SpectralOpts.SigmaScaling,'Automatic')
        min_th_idx = ceil(length(th)/5); %largest th index to use as a sigma value
        max_th_idx = ceil(length(th)/2); %largest th index to use as a sigma value
        Sigma = th(min_th_idx:max_th_idx);
    elseif strcmp(SpectralOpts.SigmaScaling,'Manual')
        Sigma = SpectralOpts.SigmaValues;
    end
     
    flag = zeros(size(Sigma));
    EigVals=zeros(NumEig,length(Sigma));
    EigVecs=zeros(size(denoisedX,1),NumEig,length(Sigma));
    
    LLPD=full(LLPD);
    LLPD(LLPD==0)=Inf; 

    
    for i=1:length(Sigma)
        W=exp(-LLPD.^2/Sigma(i).^2);
    

        W_SS=sparse(W.*S);
     
 
        D=zeros(size(W));
        for j=1:size(W_SS,1)      
                D(j,j)=sum(W_SS(j,:));   
        end
        I=eye(size(W,1));
        D_inv=D^(-1/2);
        W_SS_1=D_inv * W_SS;
        W_SS = W_SS_1 * D_inv;
        L_ss= sparse(I - W_SS);
     
 
        [V_temp, D_temp, flag(i)] = eigs(L_ss,NumEig,'sr');
        [EigVals(:,i), I] = sort(diag(D_temp), 'ascend');
         EigVecs(:,:,i) = V_temp(:,I); 
        
    end
    
    GoodSigma_idx = find(flag==0);
    Sigma = Sigma(GoodSigma_idx);
    EigVals = EigVals(:,GoodSigma_idx);
    EigVecs = EigVecs(:,:,GoodSigma_idx);
    FullEigVecs=EigVecs;
   
end