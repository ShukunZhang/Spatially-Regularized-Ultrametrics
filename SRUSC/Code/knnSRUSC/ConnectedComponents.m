%Computes a sequence of scales and calculates the connected components and at each scale

function [sortedCCmatrix, th] = ConnectedComponents(X, LLPDopts)
    
% Input: 
% 
% -X: Data
% -LLPDopts: Structure of options for computing LLPD
%
% Output:
%
% -sortedCCmatrix: a matrix containing the connected components at each scale
% and containing point index as the final column. This matrix is sorted
% from rightmost column to left for efficient path distance querying.
% -th: vector of thresholds corresponding to different scales



%% Make graph on X

Gknn= ComputeAdjacency(X, LLPDopts);
n=size(X,1);
EdgesSorted=sort(Gknn(Gknn>0));
EdgesSorted=full(EdgesSorted);

%% Determine thresholds

if isempty(LLPDopts)
    LLPDOpts.ThresholdMethod = 'LogScaling';
    LLPDOpts.UseFixedNumScales = 1;
    LLPDOpts.NumScales = 20;
end

if strcmp(LLPDopts.ThresholdMethod, 'LogScaling')
    
    tmin = prctile(EdgesSorted(EdgesSorted>0),1); %Take 1st percentile as smallest scale
    tmax = 1.01*max(max(EdgesSorted));
    if LLPDopts.UseFixedNumScales == 1
        num_scales = LLPDopts.NumScales;
        th = exp(linspace(log(tmin),log(tmax),num_scales));
    else
        num_scales = log(tmax/tmin)/log(LLPDopts.LogRatio);
        th = tmin.*(LLPDopts.LogRatio).^(1:ceil(num_scales));
    end
    th=[0,th]; %Necessary to ensure one CC at largest scale

    
elseif strcmp(LLPDopts.ThresholdMethod, 'PercentileScaling')
    
    if LLPDopts.UseFixedNumScales == 1
        PercentileSize = 100/LLPDopts.NumScales;
    else
        PercentileSize = LLPDopts.PercentileSize;
    end
    th=prctile(EdgesSorted(EdgesSorted>0),0:PercentileSize:100);
    th=[th,1.01*max(max(EdgesSorted))];
    
end

%% Threshold graph

r=length(th)-1;
[I,J]=find(Gknn>0);
Gbin=sparse(I,J,discretize(Gknn(Gknn>0),th,'IncludedEdge','right'));

E{1}=Gknn.*(Gbin==1);
G{1} = E{1};
[~,CC{1}]=AlternateConnComp(G{1});

for i=2:r
    E{i}=Gknn.*(Gbin==i);
    G{i} = G{i-1}+E{i}; 
    [S{i}, CC{i}] = AlternateConnComp(G{i});
end

%% Create Sorted Matrix of Component Indices

m=length(CC); %number of scales
CCmatrix = zeros(n,m);

for j=1:m
    CCmatrix(:,j)=CC{j};
end

point_index = (1:n)';
CCmatrix = [CCmatrix point_index];
sortedCCmatrix = sortrows(CCmatrix, m:-1:1); 

% Note: last column contains the point index so that when the rows are 
% sorted the original index can be recovered

end