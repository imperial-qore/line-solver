function [alpha, T] = lqn_ph_serial_law(wf)
% [ALPHA, T] = LQN_PH_SERIAL_LAW(WF)
%
% Composed law of a workflow in which the branches of an AND fork are SERIAL
% rather than concurrent, that is, the total work the branches request rather
% than the elapsed time until the last of them finishes.
%
% This is the law of the PROCESSOR demand of an LQN entry. Two branches of an
% AND fork are two activity threads of the same task instance: they overlap in
% time, so the entry response time is the maximum of the branches, but they run
% on ONE processor, so the demand they place on it is the sum. Composing the
% host law with Workflow.toPH would charge the processor the maximum and let the
% layer report a utilization below the true one, which no amount of iterating
% recovers.
%
% Every other node keeps its own composition rule: an OR fork is a mixture, a
% loop is a geometric compound, so the correlation within a branch survives.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tree = wf.getSPTree();
if isempty(tree)
    line_error(mfilename, sprintf(['Workflow %s is not series-parallel, so it has no exact ' ...
        'phase-type reduction.'], wf.getName()));
end
[alpha, T] = composeSerialized(wf, tree, tree.root);
end

% ------------------------------------------------------------------------
function [alpha, T] = composeSerialized(wf, tree, k)
kids = tree.kids{k};
switch tree.type{k}
    case 'leaf'
        [alpha, T] = wf.activities{tree.act(k)}.getPHRepresentation();
    case {'serial','par'}
        [alpha, T] = composeSerialized(wf, tree, kids(1));
        for i = 2:numel(kids)
            [a2, T2] = composeSerialized(wf, tree, kids(i));
            [alpha, T] = Workflow.composeSerial(alpha, T, a2, T2);
        end
    case 'or'
        alphas = cell(1, numel(kids));
        Ts = cell(1, numel(kids));
        for i = 1:numel(kids)
            [alphas{i}, Ts{i}] = composeSerialized(wf, tree, kids(i));
        end
        [alpha, T] = Workflow.composeMixture(alphas, Ts, tree.probs{k});
    case 'loop'
        [a1, T1] = composeSerialized(wf, tree, kids(1));
        [alpha, T] = Workflow.composeLoopGeometric(a1, T1, tree.count(k));
    otherwise
        line_error(mfilename, sprintf('Unknown series-parallel node type "%s".', tree.type{k}));
end
end
