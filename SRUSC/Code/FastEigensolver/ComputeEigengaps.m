% Compute the maximum eigengap across a series of scales, i.e. compute
% K_hat = argmax_(k,sigma) [lambda_k+1(sigma)-lambda_k(sigma)]

function [Khat,MaxGapScale,SizeMaxGap,SigmaIndexMaxGap,NewIndex]=ComputeEigengaps(EigVals)

% Input: 
%
% -EigVals: NumEig x NumScales matrix of real numbers (since L is symmetric)
% 
% Output:
% 
% -Khat: argmax_(k,sigma) [lambda_k+1(sigma)-lambda_k(sigma)], i.e. maximum
% of MaxGapScale across scales
% -MaxGapScale: Maximum eigengap index at each scale
% -SizeMaxGap:  Maximum eigengap value at each scale
% -SigmaIndexMaxGap: Sigma index corresponding to maximum gap
% -NewIndex:  Scales where the largest gap is not between the first and second
% eigenvalue

NumEig = size(EigVals,1); %Number of eigenvalue curves
NumScales = size(EigVals,2); %Number of scales

%% Compute gaps
Gaps = zeros(NumEig-1, NumScales);
MaxGapScale = zeros(1, NumScales);
SizeMaxGap = zeros(1, NumScales);
for i=1:NumScales
    Gaps(:,i) = diff(EigVals(:,i));
    MaxGapScale(i) = find(Gaps(:,i)==max(Gaps(:,i))); % index of max gap at scale sigma(i)
    SizeMaxGap(i) = EigVals(MaxGapScale(i)+1,i) - EigVals(MaxGapScale(i),i);
end

%% Throw away scales where largest gap is between lambda_1 and lambda_2

NewIndex=find(MaxGapScale>1);

%% Compute Khat

SigmaIndexMaxGap = NewIndex(SizeMaxGap(NewIndex)==max(SizeMaxGap(NewIndex)));
Khat = MaxGapScale(SigmaIndexMaxGap); 

if isempty(Khat)
    Khat=0;
end

end
