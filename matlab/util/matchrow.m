function pos = matchrow(matrix, row)
% pos = matchrow(M, r)
% Return position of row r in matrix M if unique, -1 otherwise
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if size(matrix,2) ~= size(row,2)
    %Incompatible matrix and row sizes
    pos = -2;
    return
end
if all(matrix(end,:) == row)
    pos = size(matrix,1);
else
    pos = find(all(bsxfun(@eq,matrix,row),2),1);
    if isempty(pos)
        pos = -1;
    end
end
end