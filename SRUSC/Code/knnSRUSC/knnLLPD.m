% Computes LLPD from the connected components of a multiscale sequence of 
% adjacency graphs obtained by throwing away edges which exceed a sequence 
% of scales; LLPD is estimated by determining at which 
% scale two points are in the same connected component. 

function LLPD = knnLLPD(sortedCCmatrix,th,k2)


% Input:
%
% -sortedCCmatrix: An n by m matrix giving 
% the connected components of the sequence of thresholded graphs (m is the 
% number of scales used in the path distance approximation, determined by LLPDopts)
% -th: sequence of scales corresponding to sortedCCmatrix
% -k2: Number of nearest path distance neighbors to return for each point/row of sortedCCmatrix
% 
% Note: Given data matrix X, sortedCCmatrix and th are created using the ConnectedComponents function. 
%
% Output:
%
% -LLPD: An n by n sparse matrix giving the approximate path distance of each 
% point to it's k2 nearest neighbors 


%% Determine the point index and path distance of each point's k2 nearest neighbors
n = size(sortedCCmatrix,1);
m = size(sortedCCmatrix,2)-1; %(Since last column contains point indices, not CC's)
knnindex = zeros(n,k2);
knndistances = zeros(n,k2);

for s=1:n %Loop over all points in the order they appear in the last column of sortedCCmatrix
    
    knnindex(sortedCCmatrix(s,m+1),1)=sortedCCmatrix(s,m+1); %Each point is it's first nearest neighbor
    knndistances(sortedCCmatrix(s,m+1),1)=th(2); %Define the path distance of a point to itself as the finest scale (note: th(1)=0)
    
    % Starting in row s, we search both up and down to find a point in
    % the same connected component. When one is found, it's index and
    % path distance are recorded. When one is not find, we move left in
    % sortedCCmatrix to search for points that are in the same
    % connected component at the next larger scale.
    
    iup = s; %will tell us the index of a nearest neighbor when it is found
    idown = s; %will tell us the index of a nearest neighbor when it is found
    counter = 2; %keeps track of how many nearest neighbors you have found
    for j=1:m %defines the approximate path distance
        % Move up until you find a point not in the same CC then move
        % left and repeat
        while (iup>1 && sortedCCmatrix(iup,j)==sortedCCmatrix(iup-1,j) && counter<k2+1)
            iup = iup-1;
            knnindex(sortedCCmatrix(s,m+1),counter)=sortedCCmatrix(iup,m+1);
            knndistances(sortedCCmatrix(s,m+1),counter)=th(j+1);
            counter=counter+1;
        end
        % Move down until you find a point not in the same CC then move
        % left and repeat
        while (idown<n && sortedCCmatrix(idown,j)==sortedCCmatrix(idown+1,j) && counter<k2+1)
            idown = idown+1;
            knnindex(sortedCCmatrix(s,m+1),counter)=sortedCCmatrix(idown,m+1);
            knndistances(sortedCCmatrix(s,m+1),counter)=th(j+1);
            counter=counter+1;
        end
    end
    % If we have not found k2 NN b/c we have disconnected components at
    % the largest scale, simply move up/down to add additional
    % neighbors at infinite distance
    while (iup>1 && counter<k2+1)
        iup = iup-1;
        knnindex(sortedCCmatrix(s,m+1),counter)=sortedCCmatrix(iup,m+1);
        knndistances(sortedCCmatrix(s,m+1),counter)=Inf;
        counter=counter+1;
    end
    while (idown<n && counter<k2+1)
        idown = idown+1;
        knnindex(sortedCCmatrix(s,m+1),counter)=sortedCCmatrix(idown,m+1);
        knndistances(sortedCCmatrix(s,m+1),counter)=Inf;
        counter=counter+1;
    end
end

% Construct sparse LLPD matrix:

Base=ones(n,k2);

for i=1:n
    Base(i,:)=i*Base(i,:);
end

knndistances=knndistances';
knnindex=knnindex';
Base=Base';

LLPD = sparse(Base(:),knnindex(:),knndistances(:));
LLPD = max(LLPD, LLPD');


end
