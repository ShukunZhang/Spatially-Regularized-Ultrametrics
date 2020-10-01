% Compute the Euclidean eigendecomposition.  This is done without the use of the learned
% blackbox eigensolver available for LLPD.  

function [EigVals,EigVecs, Sigma] = EuclideanEigensolver(X,LLPDopts,ComparisonOpts)

% Input: 
%
% X: nxD data matrix
% LLPDopts: Structure of options for how to compute adjacency matrix
% ComparisonOpts: Structure of options for how to run Euclidean spectral
% clustering as a comparison method
%
% Output:
%
% EigVals: Eigenvalues of Euclidean graph Laplacian
% EigVecs: Eigenvectors of Euclidean graph Laplacian
% Sigma: Range of sigma scaling values used.  

%% Compute Euclidean distance matrix D

if strcmp(ComparisonOpts.EucSCSparsity,'KNN')

    D=ComputeAdjacency(X,LLPDopts);
    lin_idx_D = find(D>0);
    
    % Define Sigma scales
    if strcmp(ComparisonOpts.EucSCSigmaScaling,'Automatic')
        tmin = prctile(D(lin_idx_D),20); %Take 20th percentile as smallest scale
        tmax = max(max(D(lin_idx_D)));
        num_scales = ComparisonOpts.EucSCNumScales;
        Sigma = exp(linspace(log(tmin),log(tmax),num_scales));
    elseif strcmp(ComparisonOpts.EucSCSigmaScaling,'Manual')
        Sigma = ComparisonOpts.EucSCSigmaValues;
    end

    
elseif strcmp(ComparisonOpts.EucSCSparsity,'None')
    
    n=size(X,1);
    D=zeros(n,n);
    for i=1:n
        for j=i+1:n
            D(i,j) = norm(X(i,:)-X(j,:));
            D(j,i) = D(i,j);
        end
    end
    
    % Define Sigma scales
    if strcmp(ComparisonOpts.EucSCSigmaScaling,'Automatic')
        tmin = prctile(D(D>0),5); %Take 5th percentile as smallest scale
        tmax = prctile(D(D>0),50);
        num_scales = ComparisonOpts.EucSCNumScales;
        Sigma = exp(linspace(log(tmin),log(tmax),num_scales));
    elseif strcmp(ComparisonOpts.EucSCSigmaScaling,'Manual')
        Sigma = ComparisonOpts.EucSCSigmaValues;
    end
        
end

%% Compute Euclidean graph Laplacian and spectral decomposition.

NumEig = ComparisonOpts.EucSCNumEig;
if length(Sigma)==1
    if strcmp(ComparisonOpts.EucSCSparsity,'KNN')
        W = sparse(size(X,1), size(X,1));
        W(lin_idx_D) = exp(-D(lin_idx_D).^2./Sigma.^2);
    elseif strcmp(ComparisonOpts.EucSCSparsity,'None')
        W=exp(-D.^2./Sigma.^2);
    end
    Deg=sum(W,1);
    L=eye(size(W))-diag(Deg.^(-.5))*W*diag(Deg.^(-.5));
    [V, D] = eigs(L,NumEig,'sr');
    [EigVals, I] = sort(diag(D), 'ascend');
    EigVecs = V(:,I);
end

%% If Sigma is a vector, compute a spectral decomposition for each scale in Sigma

if length(Sigma)>1
    flag = zeros(size(Sigma));
    for i=1:length(Sigma)
        if strcmp(ComparisonOpts.EucSCSparsity,'KNN')
            W = sparse(size(X,1), size(X,1));
            W(lin_idx_D) = exp(-D(lin_idx_D).^2./Sigma(i).^2);
        elseif strcmp(ComparisonOpts.EucSCSparsity,'None')
            W=exp(-D.^2./Sigma(i).^2);
        end
        Deg=sum(W,1);
        L=eye(size(W))-diag(Deg.^(-.5))*W*diag(Deg.^(-.5));
        [V_temp, D_temp, flag(i)] = eigs(L,NumEig,'sr');
        [EigVals(:,i), I] = sort(diag(D_temp), 'ascend');
        EigVecs(:,:,i) = V_temp(:,I); 
    end
    % Discard sigma values where all eigenvalues failed to converge
    GoodSigma_idx = find(flag==0);
    Sigma = Sigma(GoodSigma_idx);
    EigVals = EigVals(:,GoodSigma_idx);
    EigVecs = EigVecs(:,:,GoodSigma_idx);
end
    


end