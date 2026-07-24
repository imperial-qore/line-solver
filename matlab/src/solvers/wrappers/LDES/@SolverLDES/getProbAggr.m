function ProbAggr = getProbAggr(self, node, state)
% PROBAGGR = GETPROBAGGR(NODE, STATE)
% Aggregated state probability at a node. For LDES the sample paths are already
% per-class, so this equals getProb(). Fully JSON-mediated. Mirrors the
% Python-native getProbAggr().
if nargin < 3
    state = [];
end
ProbAggr = self.getProb(node, state);
end
