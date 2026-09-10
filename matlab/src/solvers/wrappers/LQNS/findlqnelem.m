function idx = findlqnelem(lqn, name, elemType)
% IDX = FINDLQNELEM(LQN, NAME, ELEMTYPE)
% Index of the layered element called NAME whose lqn.type is ELEMTYPE, or -1
% when the struct carries no such element. ELEMTYPE is a
% LayeredNetworkElement constant (HOST, TASK, ENTRY, ACTIVITY, CALL).
%
% THE KIND IS PART OF THE KEY, and has to be. A LINE-generated layered model
% routinely gives a processor, its task and that task's entry the SAME name,
% and lqn.names holds all three, so findstring(lqn.names,name) returns EVERY
% match and an assignment indexed by it writes one element's result into all of
% them: the processor row ends up carrying the task's numbers and then the
% entry's, so one .lqxo yields three different wrong answers. The result file
% states which kind each row describes -- it is the tag being read -- so the
% ambiguity does not have to exist. On a model whose names are unique this
% agrees with findstring element for element.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

match = strcmp(name, lqn.names(:).') & (lqn.type(:).' == elemType);
idx = find(match);
if isempty(idx)
    idx = -1;
else
    % First declaration wins, as find(strcmp(...))(1) would.
    idx = idx(1);
end
end
