%Function that computes the weighted adjacency matrix of the graph formed by
%connecting each point of X with its k nearest neighbors (ensures it is
%symmetric; edge is present if either x is among y's knn or vice versa);
%weights are the Euclidean distances

function A = ComputeAdjacency(X,LLPDopts)

% Input:
%
% -X: data matrix
% -LLPDopts: A structure containing options for computing adjacency
% matrix, including:
%   -LLPDopts.FullyConnected: either 0 or 1;
%   if the knn graph is not connected and LLPDopts.FullyConnected=1,the algorithm will randomly
%   sample points from the connected components and add an approximately
%   minimal edge connecting the CC's
%   -LLPDopts.NumberofSamplePoints : specify how many random points to
%   sample from each connected component; required if LLPDopts.FullyConnected=1
%
% Output:
%
% -A: (sparse) adjacency matrix

k = LLPDopts.KNN;
n=size(X,1);
Base=ones(n,k);

[IDX,D] = knnsearch(X,X,'K',k);

for i=1:n
    Base(i,:)=i*Base(i,:);
end

D=D';
D=double(D(:));
IDX=IDX';
IDX=IDX(:);
Base=Base';
Base=Base(:);

try
    A=sparse(Base,IDX,D,n,n);
catch
    keyboard
end

try
    A = max(A, A');
catch
    keyboard
end

if LLPDopts.FullyConnected==1
    m = LLPDopts.NumberofSamplePoints;
    % Add in edges connected the CC's:
    [~, CC] = AlternateConnComp(A);
    NumberCC = size(unique(CC),2);
    RandomSampleCC = zeros(NumberCC,m);
    for i=1:NumberCC
        CCindex{i} = find(CC==i);
        if length(CCindex{i}) == 1
            RandomSampleCC(i,:) = CCindex{i}*ones(1,m);
        else
           RandomSampleCC(i,:) = randsample(CCindex{i},m,true); %for sampling with replacement
        end
    end

I=[];
J=[];
MinDistance=[];

for i=1:NumberCC
    for j=i+1:NumberCC
        SampleDistances = zeros(m,m);
        for s=1:m
           for q=1:m
               SampleDistances(s,q) = norm(X(RandomSampleCC(i,s),:)-X(RandomSampleCC(j,q),:));
           end
        end
        MinDistance(end+1) = min(min(SampleDistances));
        [s0, q0] = find(SampleDistances==MinDistance(end));
        I(end+1) = RandomSampleCC(i,s0(1));
        J(end+1) = RandomSampleCC(j,q0(1));
    end
end

A=A+sparse(I,J,MinDistance,n,n);
A=max(A,A');

end

end
