% Compute accuracy statistics of a set of labels with respect to the group
% truth.  

function [OA, AA, CA, kappa, LabelsAligned] = GetAccuracies(Labels,GT,K)

% Input:
%
% -Labels: Learned labels
% -GT: Ground truth labels 
% -K: Number clusters
%
% Output: 
%
% -OA: Overall accuracy
% -AA: Average classwise accuracy
% -CA: Individual class accuracies
% -kappa: Cohen's kappa
% -LabelsAligned: Labels aligned with GT, as performed by AlignClusters

%% Align clusters

LabelsAligned=AlignClustersHungarian(GT,Labels,K);

%% Compute classwise accuracies

for k=1:K
    temp=find(GT==k);
    ClassIdxs{k}=squeeze(temp);
    ClassSize(k)=squeeze(length(ClassIdxs{k}));
    
    Assignments{k}=LabelsAligned(ClassIdxs{k});
    NumCorrect(k)=length(find(Assignments{k}==k));
    CA(k)=NumCorrect(k)/ClassSize(k);
    
    for r=1:K
        M(k,r)=length(find(Assignments{k}==r));
    end
    
end

%% Compute statistics

AA=mean(CA);
OA=sum(NumCorrect)/sum(ClassSize);
p=sum(M,2)'*sum(M)'/(sum(sum(M)))^2;
kappa=(OA-p)/(1-p);

end

