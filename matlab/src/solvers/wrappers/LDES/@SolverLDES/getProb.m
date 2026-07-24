function Prob = getProb(self, node, state)
% PROB = GETPROB(NODE, STATE)
% Steady-state marginal state probability at a node, estimated as the
% time-weighted fraction of the simulated trajectory spent in STATE. Fully
% JSON-mediated (via sample()). Mirrors the Python-native getProb().
%
% If STATE is omitted, the model's current state for the node is used. Returns
% 0 if the state is not observed during the simulation.

if nargin < 3
    state = [];
end
if GlobalConstants.DummyMode
    Prob = NaN;
    return
end

sampleResult = self.sample(node, 0);
if isempty(sampleResult) || ~isstruct(sampleResult) || isempty(sampleResult.t)
    Prob = 0;
    return
end

t = sampleResult.t(:);
stateMatrix = sampleResult.state;
nt = numel(t);
if nt < 2 || isempty(stateMatrix)
    Prob = 0;
    return
end

sn = self.model.getStruct;
if isa(node, 'Node')
    nodeIdx = node.index;
else
    nodeIdx = node;
end
if isempty(state)
    isf = sn.nodeToStateful(nodeIdx);
    if iscell(sn.state) && numel(sn.state) >= isf && ~isempty(sn.state{isf})
        state = sn.state{isf};
    else
        state = zeros(1, size(stateMatrix, 2));
    end
end
target = state(:).';
L = min(numel(target), size(stateMatrix, 2));
target = target(1:L);

total = t(end) - t(1);
if total <= 0
    Prob = 0;
    return
end

timeInState = 0;
for ti = 1:nt-1
    dt = t(ti+1) - t(ti);
    if all(abs(stateMatrix(ti, 1:L) - target) < 1e-10)
        timeInState = timeInState + dt;
    end
end
Prob = timeInState / total;
end
