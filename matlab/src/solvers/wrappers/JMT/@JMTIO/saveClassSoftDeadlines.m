function [simDoc, node] = saveClassSoftDeadlines(self, simDoc, node, ind)
% [SIMDOC, NODE] = SAVECLASSSOFTDEADLINES(SIMDOC, NODE, NODEIDX)
%
% Writes the per-class due dates an EDD or EDF buffer is ordered by, as the
% <classSoftDeadlines> child of <node> that SIMmodeldefinition.xsd puts ahead
% of the sections. jmt.engine.simEngine.SimLoader parses its element children
% positionally into a double[] and hands it to Queue.setSoftDeadlines, and
% Queue.process then stamps each arriving job with softDeadlines(classId) plus
% the current time. The ORDER therefore has to be the order saveClasses.m
% emits <userClass> in, since that is what fixes JobClass.getId, which is why
% the same getExportableClasses filter is applied here.
%
% ONLY EDD AND EDF STATIONS GET THE ELEMENT. No other put strategy reads the
% array, and writing it everywhere would also start feeding the station-level
% Tardiness measure, which is a separate question from buffer ordering.
%
% The deadline is RELATIVE to arrival AT THE STATION, since the engine adds
% the current time on arrival. That is the same reading as the reference
% SolverLDES, whose Solver_ssj sets absoluteDeadline to
% ssjSim.time() + classDeadline(r) at each arriveAtQueue, so the two engines
% agree on what sn.classdeadline means rather than only on its name.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;
if ~sn.isstation(ind)
    return
end
sched = sn.sched(sn.nodeToStation(ind));
if sched ~= SchedStrategy.EDD && sched ~= SchedStrategy.EDF
    return
end

% The same predicate SolverJMT.supportsModelMethod asks through
% jmtMethodRefusal, so the gate that decides whether to OFFER a JSIM method
% and this writer cannot answer differently.
reason = jmtDeadlineRefusal(sn);
if ~isempty(reason)
    line_error(mfilename, reason);
end

exportClasses = self.getExportableClasses();
deadlinesNode = simDoc.createElement('classSoftDeadlines');
for r = 1:sn.nclasses
    if ~exportClasses(r)
        continue;
    end
    if isfinite(sn.classdeadline(r))
        deadline = sn.classdeadline(r);
    else
        % jmtDeadlineRefusal has already refused every class this station
        % SERVES without a deadline, so reaching here means class r is not
        % served at this station and the engine never reads the slot. 0.0 is
        % the value saveClasses.m already writes for an absent deadline.
        deadline = 0;
    end
    deadlineNode = simDoc.createElement('softDeadline');
    deadlineNode.appendChild(simDoc.createTextNode(sprintf('%.12f', deadline)));
    deadlinesNode.appendChild(deadlineNode);
end
node.appendChild(deadlinesNode);
end
