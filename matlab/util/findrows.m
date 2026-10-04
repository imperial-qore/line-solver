function idxs = findrows(M, r)
% idxs = FINDROWS(M, r)
% Find index of rows in matrix M identical to row vector r
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
I = size(M,1);
idxs=1:I;
for col=1:size(M,2)
    idxs = idxs(M(idxs,col)==r(col));
    if length(idxs)==1 && all(all(M(idxs,:)==r,1))
        return
    end
end
end