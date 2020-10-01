% Forces the ground truth classes to run from 1 to K=number classes.

function GT_New=UniqueGT(GT)

%Input:
%
% -GT: Original ground truth data
%
%Output:
%
% -New, modified ground truth data with simplified class labels

GT_Classes=unique(GT(GT>0));
for k=1:length(GT_Classes)
    GT(find(GT==GT_Classes(k)))=k;
end

GT_New=GT;

end