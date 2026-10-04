function S=cellsum(C)
% S=CELLSUM(C)
% Returns sum of non-empty elements in cell array C
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
S = [];
coder.varsize('S');
% numel, NOT length: length(C) is max(size(C)), so on a 2-D cell with both
% dimensions >= 2 it would walk only the first max(R,C) entries in
% column-major order and silently return a partial sum. Every caller today
% passes a 1-D cell (visits/nodevisits are cell(1,nchains), Dfilt is
% cell(1,A), the LN and SSA accumulators are cell(1,n)), for which numel and
% length agree, so this is a no-op guard rather than a behaviour change. The
% one 2-D caller that did exist, cellsum(P) over the K x K class-pair cell in
% @MNetwork/link.m, was a real defect and has been removed.
for i=1:numel(C)
    if ~isempty(C{i})
        if isempty(S)
            S = C{i};
        else
            S = S + C{i};
        end
    end
end
end