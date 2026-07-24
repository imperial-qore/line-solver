function ProbSys = getProbSys(self)
% PROBSYS = GETPROBSYS()
% Joint steady-state probability of the current system state, estimated as the
% time-weighted fraction of the system-wide simulated trajectory in that joint
% state. Fully JSON-mediated (via sampleSys()). Mirrors Python-native getProbSys().

if GlobalConstants.DummyMode
    ProbSys = NaN;
    return
end

sysResult = self.sampleSys(0);
if isempty(sysResult) || ~isstruct(sysResult) || isempty(sysResult.t)
    ProbSys = 0;
    return
end

t = sysResult.t(:);
states = sysResult.state;
nt = numel(t);
if nt < 2 || isempty(states)
    ProbSys = 0;
    return
end

sn = self.model.getStruct;
nstateful = numel(states);
targets = cell(1, nstateful);
for isf = 1:nstateful
    if iscell(sn.state) && numel(sn.state) >= isf && ~isempty(sn.state{isf})
        targets{isf} = sn.state{isf}(:).';
    else
        targets{isf} = zeros(1, size(states{isf}, 2));
    end
end

total = t(end) - t(1);
if total <= 0
    ProbSys = 0;
    return
end

timeInState = 0;
for ti = 1:nt-1
    dt = t(ti+1) - t(ti);
    allMatch = true;
    for isf = 1:nstateful
        L = min(numel(targets{isf}), size(states{isf}, 2));
        if ~all(abs(states{isf}(ti, 1:L) - targets{isf}(1:L)) < 1e-10)
            allMatch = false;
            break
        end
    end
    if allMatch
        timeInState = timeInState + dt;
    end
end
ProbSys = timeInState / total;
end
