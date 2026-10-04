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
    % TASKS PER LINK. The tag-augmented construction carries an integer weight
    % per branch, so a link whose count is an integer is served exactly, whether
    % or not it differs from its siblings'. What it cannot carry is a count that
    % is not fixed at build time: a DISTRIBUTION has to be drawn per firing, and
    % the auxiliary class capacity would have to be its support's maximum with a
    % different tag occupancy per draw.
    if isfield(sn.nodeparam{f},'fanOutDist') && ~all(cellfun(@isempty, sn.nodeparam{f}.fanOutDist(:)))
        line_error(mfilename,'A random tasks-per-link distribution (Fork.setTasksPerLinkDistribution) is not supported by the native CTMC/SSA fork-join implementation; use SolverJMT or SolverLDES, which draw the degree at the fork epoch.');
    end
    if isfield(sn.nodeparam{f},'fanOutLink') && ~isempty(sn.nodeparam{f}.fanOutLink)
        taken = sn.nodeparam{f}.fanOutProb > 0;
        vals = sn.nodeparam{f}.fanOutLink(taken);
        if any(vals ~= round(vals))
            line_error(mfilename,'Non-integer tasksPerLink is not supported by the native CTMC/SSA fork-join implementation.');
        end
    elseif isfield(sn.nodeparam{f},'fanOut') && any(sn.nodeparam{f}.fanOut ~= round(sn.nodeparam{f}.fanOut))
        line_error(mfilename,'Non-integer tasksPerLink is not supported by the native CTMC/SSA fork-join implementation.');
    end

    % BRANCH PROBABILITIES. A branch that may decline makes the SET of siblings
    % random, so the firing has 2^B outcomes and the Join's required count is a
    % function of which subset fired. The tag construction records no such
    % per-firing state -- `fj.required{r}` is fixed when the state space is
    % built -- so this path would emit every branch anyway and answer with the
    % certain-fork number. It is refused by name instead: silently returning the
    % model's answer for a DIFFERENT model is the one outcome worth avoiding.
    if isfield(sn.nodeparam{f},'fanOutProb') && ~isempty(sn.nodeparam{f}.fanOutProb)
        pb = sn.nodeparam{f}.fanOutProb;
        taken = pb > 0;
        if any(pb(taken) < 1)
            line_error(mfilename,'A branch activation probability below one (Fork.setBranchProbability) is not supported by the native CTMC/SSA fork-join implementation: the sibling SET would be random and the tag construction fixes it when the state space is built. Use SolverJMT or SolverLDES, which draw the activation at the fork epoch, or SolverMVA, whose MMT transform sees the expected degree.');
        end
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
    % JOIN STRATEGY. PARTIAL is served: `fjsn.nodeparam{j}.fj.required{r}` is a
    % per-branch count, so a quorum is that count lowered rather than a
    % different mechanism. What it needs is a quorum that is reachable -- a
    % quorum above the siblings the fork emits would never fire.
    if isfield(sn.nodeparam{j},'joinStrategy')
        for r=1:K
            if length(sn.nodeparam{j}.joinStrategy) < r || isempty(sn.nodeparam{j}.joinStrategy{r})
                continue
            end
            st = sn.nodeparam{j}.joinStrategy{r};
            if st ~= JoinStrategy.STD && st ~= JoinStrategy.PARTIAL
                line_error(mfilename,'Only JoinStrategy.STD and JoinStrategy.PARTIAL are supported by the native CTMC/SSA fork-join implementation.');
            end
            if st == JoinStrategy.PARTIAL && isfield(sn.nodeparam{j},'fanIn') && ...
                    length(sn.nodeparam{j}.fanIn) >= r && ~isempty(sn.nodeparam{j}.fanIn{r})
                fq = sum(sn.nodeparam{j}.fanIn{r});
                fsrc = find(sn.fj(:,j), 1);
                if ~isempty(fsrc) && isfield(sn.nodeparam{fsrc},'fanOutLink') && ...
                        ~isempty(sn.nodeparam{fsrc}.fanOutLink)
                    emitted = sum(sn.nodeparam{fsrc}.fanOutLink(:,r) .* (sn.nodeparam{fsrc}.fanOutProb(:,r) == 1));
                    if fq > emitted
                        line_error(mfilename,'A partial Join asks for more siblings than its Fork is certain to emit, so it could never fire.');
                    end
                end
            end
        end
    end
end

end
