function S=krons(A,B)
% S=KRONS(A,B)
% Kronecker sum of matrices A and B
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if issparse(A) || issparse(B)
    S=kron(A,speye(size(B)))+kron(speye(size(A)),B);
else
    S=kron(A,eye(size(B)))+kron(eye(size(A)),B);
end
end