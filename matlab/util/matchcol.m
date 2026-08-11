function pos = matchcol(matrix, col)
% pos = MATCHCOL(M, c)
% Return position of column c in matrix M if unique, -1 otherwise
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if all(matrix(:,end) == col)
    pos = size(matrix,2);
else
    pos = find(all(bsxfun(@eq,matrix,col),1),1);
    if isempty(pos)
        pos = -1;
    end
end
end