function varargout = getAvgOrbitTable(self,varargin)
% ORBITTABLE = GETAVGORBITTABLE()
%
% Table of the mean orbit length of every retrial station-class pair, with the
% station population and the in-service population it decomposes into.
%
% Reported as a separate table rather than as an extra column of getAvgTable so
% that the average table keeps its shape for models without retrials.
%
% See also NetworkSolver.getAvgOrbit, Queue.setOrbit
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
% The result recorder captures the returned table together with the solver
% that produced it, so cross-codebase parity is asserted against the values a
% solver RETURNED rather than the text it printed. Off unless a run asked for
% it (LineResultRecorder.enable), and then it costs one appdata lookup here.
% The wrapper exists so that recording happens on EVERY exit path, including
% the early returns inside the implementation below.
[scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
[varargout{1:max(nargout,1)}] = getAvgOrbitTable_impl(self,varargin{:});
LineResultRecorder.capture(scope, self, 'orbit', varargout{1});
end

function OrbitTable = getAvgOrbitTable_impl(self)
% GETAVGORBITTABLE_IMPL Implementation of GETAVGORBITTABLE; see the wrapper above.

ON = self.getAvgOrbit();
[QN,~,~,TN] = self.getAvg();
sn = self.model.getStruct();

Station = {};
JobClass = {};
QLen = [];
InService = [];
Orbit = [];
for ist = 1:size(ON,1)
    for r = 1:size(ON,2)
        if size(sn.retrialProc,1) < ist || size(sn.retrialProc,2) < r ...
                || isempty(sn.retrialProc{ist,r})
            continue
        end
        Station{end+1,1} = sn.nodenames{sn.stationToNode(ist)}; %#ok<AGROW>
        JobClass{end+1,1} = sn.classnames{r}; %#ok<AGROW>
        QLen(end+1,1) = QN(ist,r); %#ok<AGROW>
        rate_ir = sn.rates(ist,r);
        if isfinite(rate_ir) && rate_ir > 0
            InService(end+1,1) = TN(ist,r) / rate_ir; %#ok<AGROW>
        else
            InService(end+1,1) = 0; %#ok<AGROW>
        end
        Orbit(end+1,1) = ON(ist,r); %#ok<AGROW>
    end
end

OrbitTable = Table(Station, JobClass, QLen, InService, Orbit);
OrbitTable = IndexedTable(OrbitTable);
end
