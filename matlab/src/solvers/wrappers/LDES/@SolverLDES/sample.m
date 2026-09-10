function sampleResult = sample(self, node, numEvents)
% SAMPLERESULT = SAMPLE(NODE, NUMEVENTS)
% Generate a sample path (state trajectory) for a stateful node via a fully
% JSON-mediated transient LDES simulation. Mirrors the Python-native sample().
%
% Returns a struct with fields:
%   .handle     - the node handle passed in
%   .t          - (nTimePoints x 1) time vector
%   .state      - (nTimePoints x nclasses) per-class queue-length trajectory
%   .isaggregate- false
%   .nodeIndex  - node index
%   .numEvents  - number of events (horizon) used

if nargin < 3 || isempty(numEvents)
    numEvents = 0;
end

if GlobalConstants.DummyMode
    sampleResult = [];
    return
end

sn = self.model.getStruct;
if isa(node, 'Node')
    nodeIdx = node.index;
else
    nodeIdx = node;
end

res = self.runTransientJson(numEvents);

sampleResult = struct();
sampleResult.handle = node;
sampleResult.nodeIndex = nodeIdx;
sampleResult.isaggregate = false;
sampleResult.numEvents = numEvents;

if isempty(res.t) || isempty(res.QNt)
    sampleResult.t = [];
    sampleResult.state = [];
    return
end

% The engine's QNt is STATION-major (see runTransientJson), so the row is found
% by station index; sn.nodeToStateful indexed a different space and returned
% another station's path, or an empty one, whenever the two differ.
ist = sn.nodeToStation(nodeIdx);
nt = numel(res.t);
R = sn.nclasses;
state = zeros(nt, R);
if ist >= 1 && ist <= numel(res.QNt)
    classCells = res.QNt{ist};
    for k = 1:min(R, numel(classCells))
        cd = classCells{k};
        if ~isempty(cd)
            n = min(nt, size(cd, 1));
            state(1:n, k) = cd(1:n, 1);
        end
    end
end

sampleResult.t = res.t;
sampleResult.state = state;
end
