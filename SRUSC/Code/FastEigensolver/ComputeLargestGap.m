%Compute the largest gap in the eigenvalues in order to estimate the number
%of clusters K.

function [AllCandidates,AllSizeMaxGapScores,LargestMaxGapEstimate,SigmaIndexMaxGapEstimate]=ComputeLargestGap(MaxGap,SizeMaxGap,NewIndex)

% Input: 
%
% -MaxGapScale: Maximum eigengap index at each scale
% -SizeMaxGap:  Maximum eigengap value at each scale
% -NewIndex:  Scales where the largest gap is not between the first and second
% 
% Output:
% 
% -AllCandidates: Across all scales, the indices of the maximum eigengap
% -AllSizeMaxGapScores: Across all scales, the values of the maximum eigengap
% -LargestMaxGapEstimate: The index of the maximum eigengap across
% AllCandidates
% -SigmaIndexMaxGapEstimate: Sigma value corresponding to LargestMaxGapEstimate


SigmaIndexMaxGapEstimate = NewIndex(find(SizeMaxGap(NewIndex)==max(SizeMaxGap(NewIndex))));
LargestMaxGapEstimate = MaxGap(SigmaIndexMaxGapEstimate); 

AllCandidates = unique(MaxGap(NewIndex));
AllSizeMaxGapScores = zeros(1, length(AllCandidates));
for i=1:length(AllCandidates)
    AllSizeMaxGapScores(i) = max(SizeMaxGap(find(MaxGap==AllCandidates(i))));
end

end