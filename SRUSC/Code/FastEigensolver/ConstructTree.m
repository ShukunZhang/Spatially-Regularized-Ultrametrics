% Function for constructing tree for fast matrix-vector multiplication


function [Tree,p] = ConstructTree(SortedCCmatrix, th, k3)

% Input:
%
% SortedCCmatrix: matrix giving which connected component each point is in
% at each scale (last column gives point order)
% th: sequence of scales which were used to create sortedCCmatrix 
% k3: number of nearest neighbors to be used for denoising
%
%
% Output:
% Tree: structure containing all information about the tree
% p: number of CC's at the finest scale

%% Store Point Partition at Finest Scale (first columns gives CC, last column which points are in it)

PointPartitionAtFinestScale = [SortedCCmatrix(:,1) SortedCCmatrix(:,end)];

%% Remove last columns which contains point permutations (when using output of knnLLPD)

SortedCCmatrix = SortedCCmatrix(:,1:end-1); 

%% Make all CC labels distinct

m=size(SortedCCmatrix,2); %number of scales

for j=2:m
    SortedCCmatrix(:,j)=max(SortedCCmatrix(:,j-1))+SortedCCmatrix(:,j);
end

UniqueSortedCCmatrix = unique(SortedCCmatrix,'rows','stable');
p=size(UniqueSortedCCmatrix,1); %Number of CC's at the finest scale

%% Create tree structure:

Nodes = sort(unique(UniqueSortedCCmatrix));

%% Find the scale associated with each node

NodeScales = zeros(size(Nodes));

for j=1:m
    NodeScales(UniqueSortedCCmatrix(:,j)) = th(j+1);
end

%% Find parents of all nodes

Parents = zeros(size(Nodes));

for j=1:m-1
    for i=1:size(UniqueSortedCCmatrix,1)
        try 
        Parents( UniqueSortedCCmatrix(i,j)) = UniqueSortedCCmatrix(i,j+1);
        NodeScales( UniqueSortedCCmatrix(i,j)) = th(j+1);
        catch 
            keyboard
        end
    end
end



%% Find children of all nodes (using cell array):

Children = cell(size(Nodes));

for j=2:m
    for i=1:p
        if isempty(Children{ UniqueSortedCCmatrix(i,j) })
            Children{ UniqueSortedCCmatrix(i,j) }(end+1) = UniqueSortedCCmatrix(i,j-1);
        elseif Children{ UniqueSortedCCmatrix(i,j) }(end) ~= UniqueSortedCCmatrix(i,j-1)
            Children{ UniqueSortedCCmatrix(i,j) }(end+1) = UniqueSortedCCmatrix(i,j-1);
        end
    end
end

%% Find sample size associated with each node

SS=zeros(p,1);

for i=1:size(SortedCCmatrix(:,1))
    SS(SortedCCmatrix(i,1))=SS(SortedCCmatrix(i,1))+1;
end

%% Find sample sizes of nodes at larger scales by adding the SS's of their children

for i=p+1:length(Nodes)
    SS(i) = sum(SS( Children{i}));
end

%% Create SS array for fast multiplication:

full_SS_array = SS(UniqueSortedCCmatrix);

change_SS_array = [full_SS_array(:,1) diff(full_SS_array,1,2)];

%% Compute distance to k3^rd nearest neighbor for each CC

NNDistances = zeros(size(Nodes));

try
    for i=1:p
        ancestor = Parents(i);
        NN = SS(ancestor);
        while NN < k3
            ancestor = Parents(ancestor);
            NN = SS(ancestor);
        end
        NNDistances(i) = NodeScales(ancestor);
    end
catch
    keyboard
end
%% Store all needed information:
Tree.Scales = full(th(2:length(th)));
Tree.PointPartitionAtFinestScale = PointPartitionAtFinestScale; %First column contains CC index, second column point index
Tree.CC_array = UniqueSortedCCmatrix;
Tree.Nodes = Nodes;
Tree.Parents = Parents;
Tree.Children = Children;
Tree.SS = SS; % the sum of SS's associated with each node
Tree.SS_array = change_SS_array;
Tree.NNDistances = NNDistances;
end
        