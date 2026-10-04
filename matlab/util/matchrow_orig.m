function pos = matchrow_old(matrix, row)
% pos = matchrow(M, r)
% Finds position of row r in matrix M if unique
% Returns -1 otherwise
%
% Copyright (c) 2015-2022, Imperial College London
% All rights reserved.
I = size(matrix,1);
%L = 1:ceil(I/8):I;
L = [1,I];
for l = 1:(numel(L)-1)
    pos=L(l):L(l+1);
    for col=1:size(matrix,2)
        pos = pos(matrix(pos,col)==row(col));
        if length(pos)==1 && all(all(matrix(pos,:)==row,1))
            return
        end
    end
end
pos=-1;
end
