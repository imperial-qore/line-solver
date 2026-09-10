function tranSysState = sampleSys(self, numEvents)
% TRANSYSSTATE = SAMPLESYS(NUMEVENTS)
% System-wide sample path over all stateful nodes via a fully JSON-mediated
% transient LDES simulation. Mirrors the Python-native sampleSys().
%
% Returns a struct with fields:
%   .handle      - column of stateful-node handles
%   .t           - (nTimePoints x 1) time vector
%   .state       - 1 x nstateful cell; {isf} is (nTimePoints x nclasses)
%   .isaggregate - false
%   .numEvents   - horizon used

if nargin < 2 || isempty(numEvents)
    numEvents = 0;
end

if GlobalConstants.DummyMode
    tranSysState = [];
    return
end

res = self.runTransientJson(numEvents);

tranSysState = struct();
tranSysState.handle = self.model.getStatefulNodes';
tranSysState.isaggregate = false;
tranSysState.numEvents = numEvents;

if isempty(res.t) || isempty(res.QNt)
    tranSysState.t = [];
    tranSysState.state = {};
    return
end

sn = self.model.getStruct;
nt = numel(res.t);
R = sn.nclasses;
% One entry per STATEFUL node, as the contract says, each filled from that
% node's STATION row: the engine's QNt is station-major and the two index spaces
% coincide only when every stateful node is a station. A stateful node that is
% not one -- a Router, a Cache -- holds no queue-length series and keeps zeros.
nstateful = sn.nstateful;
statefulToStation = zeros(1, nstateful);
for ind = 1:sn.nnodes
    isfInd = sn.nodeToStateful(ind);
    if isfInd >= 1
        statefulToStation(isfInd) = sn.nodeToStation(ind);
    end
end
states = cell(1, nstateful);
for isf = 1:nstateful
    nodeState = zeros(nt, R);
    ist = statefulToStation(isf);
    if ist >= 1 && ist <= numel(res.QNt)
        classCells = res.QNt{ist};
        for k = 1:min(R, numel(classCells))
            cd = classCells{k};
            if ~isempty(cd)
                n = min(nt, size(cd, 1));
                nodeState(1:n, k) = cd(1:n, 1);
            end
        end
    end
    states{isf} = nodeState;
end

tranSysState.t = res.t;
tranSysState.state = states;
end
