function C=cellzeros(c,d,m,n)
% C=CELLZEROS(c,n)
% Creates a cell array of size cxc, each cell contains a nxn zero matrix
%
% C=CELLZEROS(c,d,m,n)
% Creates a cell array of size cxd, each cell contains a mxn zero matrix
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if nargin==2
    m=d;
    n=d;
    d=1;
end
if nargin==3
    n=m;
end
C = cell(c,d);
for i=1:c
    for j=1:d
        C{i,j}=zeros(m,n);
    end
end
end