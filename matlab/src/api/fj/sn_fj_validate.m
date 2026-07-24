function sn_fj_validate(sn)
% SN_FJ_VALIDATE(SN)
%
% Validate that a fork-join model is within the supported feature set of
% the native CTMC/SSA fork-join implementation (v1): closed classes only,
% non-nested fork-join pairs, standard join strategy (wait for all
% siblings), one task per output link. Unsupported models raise a
% line_error with a specific message; further structural checks (class
% switching on a branch, nesting, sibling traps) are performed during the
% branch discovery in ModelAdapter.fjtag.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = sn.nclasses;
forkIndexes = find(sn.nodetype == NodeType.Fork)';
joinIndexes = find(sn.nodetype == NodeType.Join)';
Vnodes = cellsum(sn.nodevisits);

for f=forkIndexes
    j = find(sn.fj(f,:));
    if isempty(j)
        line_error(mfilename,'Fork nodes without a matched Join are not supported by the native CTMC/SSA fork-join implementation.');
    end
    if length(j)>1
        line_error(mfilename,'Multiple Join nodes per Fork are not supported by the native CTMC/SSA fork-join implementation.');
    end
    if isfield(sn.nodeparam{f},'fanOut') && any(sn.nodeparam{f}.fanOut ~= round(sn.nodeparam{f}.fanOut))
        line_error(mfilename,'Non-integer tasksPerLink is not supported by the native CTMC/SSA fork-join implementation.');
    end
    for r=find(Vnodes(f,:)>0)
        c = find(sn.chains(:,r), 1);
        if isinf(sn.njobs(r)) || (~isempty(c) && any(isinf(sn.njobs(sn.chains(c,:)))))
            line_error(mfilename,'Open classes routed through a Fork are not supported by the native CTMC/SSA fork-join implementation.');
        end
    end
end

for j=joinIndexes
    if ~any(sn.fj(:,j))
        line_error(mfilename,'Join nodes without a matched Fork are not supported by the native CTMC/SSA fork-join implementation.');
    end
    if isfield(sn.nodeparam{j},'joinStrategy')
        for r=1:K
            if length(sn.nodeparam{j}.joinStrategy) >= r && ~isempty(sn.nodeparam{j}.joinStrategy{r}) && sn.nodeparam{j}.joinStrategy{r} ~= JoinStrategy.STD
                line_error(mfilename,'Only JoinStrategy.STD is supported by the native CTMC/SSA fork-join implementation.');
            end
        end
    end
end

end
