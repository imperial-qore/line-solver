function S=cellsum(C)
% S=CELLSUM(C)
% Returns sum of non-empty elements in cell array C
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
S = [];
for i=1:length(C)
    if ~isempty(C{i})
        if isempty(S)
            S = C{i};
        else
            S = S + C{i};
        end
    end
end
end