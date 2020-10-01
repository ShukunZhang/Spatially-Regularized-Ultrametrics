function [EigVals,FullEigVecs,IdxRetain,Sigma] = FastEigensolverDenoising(X,LLPDopts,SpectralOpts,DenoisingOpts)

% Input: 
%
% -X: nxD data matrix
% -LLPDopts: Structure of options for how to approximate LLPD
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


[SortedCCmatrixInitial, ThresholdsInitial] = ConnectedComponents(X, LLPDopts);
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
    title('Plot of LLPD to NN and Denoising Cutoff','fontsize',16)
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

[sortedCCmatrix, th] = ConnectedComponents(denoisedX, LLPDopts);

%% Build a tree 

Tree = ConstructTree(sortedCCmatrix, th, DenoisingOpts.KNN);
p = size(Tree.CC_array,1);

if strcmp(SpectralOpts.SigmaScaling,'Automatic')
    min_th_idx = ceil(length(th)/5); %largest th index to use as a sigma value
    max_th_idx = ceil(length(th)/2); %largest th index to use as a sigma value
    Sigma = th(min_th_idx:max_th_idx);
elseif strcmp(SpectralOpts.SigmaScaling,'Manual')
    Sigma = SpectralOpts.SigmaValues;
end

NumEig = SpectralOpts.NumEig;
if length(Sigma)==1
    MVM = @(x) FastMultiplication(Tree,Sigma,x,SpectralOpts);
    %[V, D] = eigs(MVM,p,p,opts);
    [V, D] = eigs(MVM,p,NumEig,'sr','SubspaceDimension',2*NumEig,'issym',1);
    [EigVals, I] = sort(diag(D), 'ascend');
    EigVecs = V(:,I); %Note: V is p by p
    FullEigVecs = zeros(size(denoisedX,1),NumEig);
    for s=1:p
        CCindex = find(Tree.PointPartitionAtFinestScale(:,1) == s);
        Point_index_of_CC = Tree.PointPartitionAtFinestScale(CCindex,2);
        FullEigVecs(Point_index_of_CC,:) = ones(length(Point_index_of_CC),1)*EigVecs(s,:);
    end
end


if length(Sigma)>1
    flag = zeros(size(Sigma));
    EigValsPath=zeros(NumEig,length(Sigma));
    EigVecsPath=zeros(size(denoisedX,1),NumEig,length(Sigma));
    for i=1:length(Sigma)
        MVM = @(x) FastMultiplication(Tree,Sigma(i),x,SpectralOpts);
        [V_temp, D_temp, flag(i)] = eigs(MVM,p,NumEig,'sr','SubspaceDimension',2*NumEig,'issym',1,'maxit',100);
        [EigVals(:,i), I] = sort(diag(D_temp), 'ascend');
        EigVecs(:,:,i) = V_temp(:,I); %Note: V is p by p
    end
    FullEigVecs = zeros(size(denoisedX,1),NumEig,length(Sigma));
    for s=1:p
        CCindex = find(Tree.PointPartitionAtFinestScale(:,1) == s);
        Point_index_of_CC = Tree.PointPartitionAtFinestScale(CCindex,2);
        for j=1:NumEig
            for i=1:length(Sigma)
                FullEigVecs(Point_index_of_CC,j,i) = EigVecs(s,j,i);
            end
        end
    end
    %Discard any sigma values where all eigenvalues did not converge
    GoodSigma_idx = find(flag==0);
    Sigma = Sigma(GoodSigma_idx);
    EigVals = EigVals(:,GoodSigma_idx);
    FullEigVecs = FullEigVecs(:,:,GoodSigma_idx);

end

    

end