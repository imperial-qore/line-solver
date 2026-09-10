function reason = jmtDeadlineRefusal(sn)
% REASON = JMTDEADLINEREFUSAL(SN)
%
% The due date an EDD or EDF station needs and LINE may not carry, as a
% sentence; '' when every such station can be exported.
%
% WHY THE RULE EXISTS. JMT 1.2.x serves EDD and EDF for real -- the classes
% jmt.engine.NetStrategies.QueuePutStrategies.EDDStrategy and .EDFStrategy are
% in JMT.jar, the SIMmodeldefinition.xsd carries <classSoftDeadlines> under
% <node>, and jmt.engine.simEngine.SimLoader feeds it to
% Queue.setSoftDeadlines. Queue.process then sets each arriving job's
% currentStationSoftDeadline to softDeadlines(classId) + now, and EDDStrategy
% orders the buffer by it. But the field is initialised to -1.0 in Job, and
% EDDStrategy.put THROWS on that value rather than degrading to FCFS:
%
%   java.lang.IllegalArgumentException: Attempting to schedule job with no
%   soft deadline   (EDDStrategy.java:19)
%
% measured on 2026-09-05 by running jmt.commandline.Jmt on a two-class M/M/1
% whose Queue put strategy was EDDStrategy and whose node carried no
% <classSoftDeadlines>. So a class without a due date is not an approximation
% here, it is an aborted run, and the model is refused by name instead.
%
% The <userClass softDeadline="..."> attribute saveClasses.m already writes
% does NOT satisfy this: SimLoader passes it to JobClass.setSoftDeadline,
% which feeds the SYSTEM tardiness and earliness measures only. The buffer
% ordering reads the per-node array and nothing else.
%
% WHAT COUNTS AS SERVED. A class whose rate at the station is NaN is not
% served there and can never enter its buffer, so it needs no due date; that
% is the same test jmtStationCapRefusal applies. saveClassSoftDeadlines fills
% those slots with 0.0, which is never read.
%
% ONE PREDICATE, TWO CALLERS. saveClassSoftDeadlines raises it while writing
% the JSIM document, and jmtMethodRefusal returns it so findSolver and
% SolverAUTO never offer a JSIM method on a model the writer will refuse.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

reason = '';
for ist = 1:sn.nstations
    sc = sn.sched(ist);
    if sc ~= SchedStrategy.EDD && sc ~= SchedStrategy.EDF
        continue
    end
    for r = 1:sn.nclasses
        if isnan(sn.rates(ist, r))
            continue % class r is not served here, so it never enters this buffer
        end
        if isfinite(sn.classdeadline(r))
            continue
        end
        reason = sprintf(['Station %s is scheduled %s but class %s carries no deadline. ', ...
            'JMT orders an EDD or EDF buffer by a per-station due date it reads from the ', ...
            '<classSoftDeadlines> element, and a job that reaches such a buffer without one ', ...
            'aborts the run with "Attempting to schedule job with no soft deadline" rather ', ...
            'than being served FCFS. Give the class a finite deadline, e.g. ', ...
            'OpenClass(model, ''%s'', prio, 5.0), or use SolverLDES, which sorts a class of ', ...
            'infinite deadline last instead of refusing it.'], ...
            sn.nodenames{sn.stationToNode(ist)}, upper(SchedStrategy.toText(sc)), ...
            sn.classnames{r}, sn.classnames{r});
        return
    end
end
end
