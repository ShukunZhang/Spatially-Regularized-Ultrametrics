% Align clusters in Clustering2 to maximize overall accuracy with
% Clustering1.

function NewLabels=AlignClustersHungarian(Clustering1,Clustering2,K)

% Input:
%
% -Clustering1: one assignment of labels of a dataset
% -Clustering2: another assignment of labels of the same dataset as
% clustering1
% -K: number of clusters; should be the same for Clustering1 and
% Clustering2
%
% Output:
%
% -NewLabels: new labels of Clustering2, designed to maximize correspondence with
% clustering 1
%
    
NewLabels=zeros(size(Clustering2));

%% If number of input classes is different from number of labels, do nothing.

if K~=length(unique(Clustering1))
    error('Class Mismatch');
else
    
    Overlap=zeros(K,K);
    
    try
        for i=1:K
            for j=1:K
                Overlap(i,j)=sum((Clustering1==i).*(Clustering2==j));
            end
        end
    catch
        
        keyboard
    end
    
    Overlap=Overlap';
            
end

%% Run the Hungarian algorithm

Cost=repmat(max(Overlap')',[1,K])-Overlap;
Assignment = munkres(Cost);

for i=1:K
    NewLabels(Clustering2==i)=Assignment(i);
end


