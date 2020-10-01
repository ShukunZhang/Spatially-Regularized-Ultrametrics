% Perform fast blackbox multiplication Lx, where the structure of L is
% encoded in the tree structure Tree.

function Y = FastMultiplication(Tree, sigma, x, SpectralOpts)

% Input: 
% 
% -Tree: tree structure encoding the multiscale hierarchy of the LLPD 
% distances
% -sigma: scaling parameter for L
% -x: nx1 vector
% -Spectralopts: Structure of options for spectral decomposition
%
% Output: 
% -Y = L(x)=Lx, multiplication by normalized Laplacian, either 
% the random walk Laplacian L=I-D^(-1)W or the symmetric normalized Laplacian
% L=I-D^(-1/2)WD^(-1/2).  
% 
% Remark: Both x and Y are ordered so that x(i), Y(i) corresponds to the value on
% the i^th CC (e.g. row 5 corresponds to the values of CC 5)


m=length(Tree.Scales); %Number of scales
p=size(Tree.CC_array,1); %Number of CC's at the finest scale

%Create the kernel values (evaluate W(scale(i),sigma))
W = zeros(1,m);
for i=1:m
    W(i) = exp(-Tree.Scales(i)^2/(2*sigma^2));
end

%% Store the degree of each component:

CCDeg(Tree.CC_array(1:p,1),1) = (Tree.SS_array*W');

%% Compute Lx, depending on type of Laplacian

if strcmp(SpectralOpts.Laplacian, 'RandomWalk_Laplacian')
    
    %First we need to construct a function which maps each node of the tree to the sum of x-values associated with the CC's represented by the node
    %(same as in the construction of SS in ConstructTree)

    sum_x_val = zeros( size(Tree.Nodes) );
    % Associate the values of x with the finest scale nodes:
    sum_x_val(1:p) = x.*Tree.SS(1:p);
    % Find sum of x_values of nodes at larger scales by adding the sum_x_val's of their
    % children
    for i=p+1:length(Tree.Nodes)
        sum_x_val(i) = sum( sum_x_val( Tree.Children{i} ) );
    end

    %Create the sum_x_val_array for fast access/multiplication
    full_sum_x_val_array = sum_x_val(Tree.CC_array);

    change_sum_x_val_array = [full_sum_x_val_array(:,1) diff(full_sum_x_val_array,1,2)];

    %Create Y=Wx

    Wx = zeros(p,1);
    for i=1:p
        Wx(Tree.CC_array(i,1)) = sum(change_sum_x_val_array(i,:).*W);
    end

    %Finally compute Y = (I - D^(-1)W)x = x - D^(-1)(Wx)
    Y = x - Wx./CCDeg;
    
elseif strcmp(SpectralOpts.Laplacian, 'Symmetric_Laplacian')
    
    %First we need to compute D^{-1/2)x, i.e. degree normalize the vector x
    degnx = zeros(size(x));
    for i=1:p
        degnx(i) = x(i)/sqrt(CCDeg(i));
    end
    
    %Now we need to construct a function which maps each node of the tree to the sum of degnx-values associated with the CC's represented by the node
    %(same as in the construction of SS in ConstructTree)

    sum_degnx_val = zeros( size(Tree.Nodes) );
    % Associate the values of x with the finest scale nodes:
    sum_degnx_val(1:p) = degnx.*Tree.SS(1:p);
    % Find sum of degnx_values of nodes at larger scales by adding the sum_x_val's of their
    % children
    for i=p+1:length(Tree.Nodes)
        sum_degnx_val(i) = sum( sum_degnx_val( Tree.Children{i} ) );
    end

    %Create the sum_degnx_val_array for fast access/multiplication
    full_sum_degnx_val_array = sum_degnx_val(Tree.CC_array);

    change_sum_degnx_val_array = [full_sum_degnx_val_array(:,1) diff(full_sum_degnx_val_array,1,2)];


    %Create Y=W(D^{-1/2}x)=W*degnx

    Wdegnx = zeros(p,1);
    for i=1:p
        Wdegnx(Tree.CC_array(i,1)) = sum(change_sum_degnx_val_array(i,:).*W);
    end

    %Finally compute Y = (I - D^(-1/2)WD^{-1/2})x = x - D^(-1/2)(Wdegx)
    Y = x - Wdegnx./sqrt(CCDeg);
    
end


end