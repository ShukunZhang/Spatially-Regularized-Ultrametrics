function [S,C] = AlternateConnComp(A)
  % Replacement for graphconncomp.m.  A is an n x n adjacency matrix, corresponding
  % to a nearest neighbor graog.  The function identifies the S
  % connected components C. This is done ala (Pothen, A. and Fan, C.J., 
  % 1990. Computing the block triangular form of a sparse matrix. ACM 
  % Transactions on Mathematical Software (TOMS), 16(4), pp.303-324).  This
  % is based on code from the gptoolbox by Alec Jacobson.
  %
  %
  % Input:
  %   -A:  n x n adjacency matrix
  % Outputs:
  %   -S:  scalar number of connected components
  %   -C:  matrix of connected component labels
  
  [p,~,r] = dmperm(A+speye(size(A)));
  S = numel(r)-1;
  C = cumsum(full(sparse(1,r(1:end-1),1,1,size(A,1))));
  C(p) = C;
end