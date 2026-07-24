function idx = findstring(myCell, myStr)
% [I] = FINDSTRING(A, B) returns the row index I
% corresponding to the string B in a cell of strings A.
% It returns -1 if the string B is not in A.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if iscell(myStr)
    idx = {};
    for c=1:length(myStr)
        idx{c} = findstring(myCell, myStr{c});
    end
else
    idx = find(strcmp(myStr,myCell));
    if isempty(idx)
        idx = -1;
    end
end
end