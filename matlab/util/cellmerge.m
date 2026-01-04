function d = cellmerge(c)
% D=CELLMERGE(C)
% Vertically stacks the matrix elements of a cell array C.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
d = c{1};
for i=2:length(c)
    d(end+1:end+length(c{i}),:)=c{i};
end
end