function members = fj_branch_members(lqn, joinaidx)
% FJ_BRANCH_MEMBERS Activities belonging to each branch of an AND-join.
%
% MEMBERS = FJ_BRANCH_MEMBERS(LQN, JOINAIDX) returns a cell array with one entry
% per branch feeding the AND-join activity JOINAIDX. Each entry lists the global
% activity indices lying on that branch, from its head (the activity spawned by
% the AND-fork) to its tail (the immediate predecessor of the join).
%
% A branch is recovered by walking backwards from each immediate predecessor of
% the join until an activity marked POST_AND is reached, that activity being the
% branch head spawned by the fork. Branches between a fork and its join are
% disjoint paths, so the walk is unambiguous.

members = {};
preds = find(lqn.graph(:, joinaidx) > 0)';
ashift = lqn.ashift;
nacts = lqn.nacts;

for tail = preds
    if tail <= ashift || tail > ashift + nacts
        continue;  % not an activity
    end
    chain = tail;
    cur = tail;
    guard = 0;
    while guard < nacts
        guard = guard + 1;
        % The branch head is the activity the fork spawned.
        if full(lqn.actposttype(cur)) == ActivityPrecedenceType.POST_AND
            break;
        end
        prevs = find(lqn.graph(:, cur) > 0)';
        prevs = prevs(prevs > ashift & prevs <= ashift + nacts);
        if numel(prevs) ~= 1
            break;  % a merge or the start of the graph: stop here
        end
        cur = prevs(1);
        chain(end+1) = cur; %#ok<AGROW>
    end
    members{end+1} = chain; %#ok<AGROW>
end
end
