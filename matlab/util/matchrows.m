function I = matchrows(matrix, rows)
% I = MATCHROWS(M, RS)
% Run matchrow(M,R) for every row R in the input matrix RS
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

I = zeros(size(rows,1),1);
for i=1:size(rows,1)
    I(i) = matchrow(matrix,rows(i,:));
end
end
