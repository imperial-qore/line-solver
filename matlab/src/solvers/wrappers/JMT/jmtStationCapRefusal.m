function reason = jmtStationCapRefusal(sn, ist)
% REASON = JMTSTATIONCAPREFUSAL(SN, IST)
%
% The binding station capacity JSIM cannot express, as a sentence; '' when the
% buffer at station IST is exportable. Refused on two counts.
%
% (1) THE RULE IS ONE JMT CANNOT READ. Its queue section recognizes exactly
% four dropStrategies strings -- 'drop', 'BAS blocking', 'waiting queue',
% 'retrial' (a lookupswitch on String.hashCode in
% jmt/engine/NodeSections/Queue.class; the Storage section of a Place is even
% narrower and drops 'retrial'). An unrecognized value falls through the
% default arm with NO flag set, so BBS, RSRD and retrial-with-limit are not
% approximated, they are IGNORED: the capacity stops being enforced and JMT
% returns the unconstrained answer. Those three are refused by name.
%
% (2) THE RULE IS WAITQ AND A CLOSED CLASS CAN REACH THE LIMIT, for the same
% reason JMTCLASSCAPASSERT in saveRegions.m refuses the per-class one: JMT
% cannot hold a blocked closed job at its upstream station. Note this is the
% case where NO blocking rule is declared. That the limit CAN be reached is the
% caller's to establish and is not retested here: both callers reach this
% function only for a capacity strictly below the reachable population, which
% is the one thing that makes a buffer a buffer. A model that does declare BAS
% is exported as JMT "BAS blocking", which is the same queueing model, under
% either declaration form -- see jmtIsBasDestination.
%
% That entry advised expressing the limit as the STATION capacity instead.
% Measured on 2026-08-19, that advice was wrong, and neither of the two
% strategies that a WAITQ station maps onto reproduces the UNDECLARED case:
%
%   waiting queue  does not enforce <size> at all. On a closed 3-queue tandem,
%                  N=6, Exp(1) FCFS, cap 2 at Q2, JMT returned the
%                  UNCONSTRAINED [2.03 1.99 1.98], X = 0.750, against the exact
%                  [3.6090 0.9711 1.4199], X = 0.6522.
%   BAS blocking   enforces it, but completes the service BEFORE blocking, so
%                  the blocked job moves the instant room frees. That is a
%                  different queueing model, not a rounding: same fixture,
%                  [2.871 1.373 1.756], X = 0.7126.
%
% With no rule declared LINE instead disables the upstream departure while the
% destination is full, which for exponential service is repetitive service (RS)
% and is what SolverCTMC, SolverSSA and SolverLDES all agree on. So THAT model
% is refused rather than exported as either of the two things JMT can say. See
% BUG-81. A declared-BAS model is a different model and is exported, not
% refused: blocking after service is precisely what JMT's "BAS blocking" does.
%
% ONE PREDICATE, TWO CALLERS. saveBufferCapacity raises it while writing the
% JSIM document, and jmtMethodRefusal returns it so that findSolver and
% SolverAUTO never offer jmt.jsim on a model the writer will refuse. It used to
% be a local function of the writer, which is why the gate could not see it.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

reason = '';
if isinf(sn.cap(ist))
    return
end
for r = 1:sn.nclasses
    if isnan(sn.rates(ist, r))
        continue % class r is not served here
    end
    dr = sn.droprule(ist, r);
    switch dr
        case {DropStrategy.BBS, DropStrategy.RSRD, DropStrategy.RETRIAL_WITH_LIMIT}
            % Unmappable for EITHER class type, so this test precedes the open-class skip
            reason = sprintf(['Station %s applies drop strategy "%s" to class %s and carries a finite ', ...
                'capacity %d it can reach. JMT''s queue section reads only "drop", "waiting queue", "BAS blocking" ', ...
                'and "retrial"; it does not approximate anything else, it ignores it, so the capacity would stop ', ...
                'being enforced and the run would return the unconstrained answer. Use SolverCTMC, SolverSSA or ', ...
                'SolverLDES.'], sn.nodenames{sn.stationToNode(ist)}, DropStrategy.toText(dr), sn.classnames{r}, ...
                sn.cap(ist));
            return
    end
    if isinf(sn.njobs(r))
        continue % open class: JMT loses its arrivals, as LINE does
    end
    if dr ~= DropStrategy.WAITQ
        continue % a mappable declared blocking rule is exported as itself
    end
    if jmtIsBasDestination(sn, ist, r)
        continue % BAS declared on the UPSTREAM station: jmtDropStrategyText
                 % moves it onto this one, which is where JMT reads it
    end
    reason = sprintf(['Station %s carries a finite capacity %d that binds for the closed class %s. ', ...
        'LINE blocks a closed job that finds no room -- the upstream departure is disabled and the job stays where ', ...
        'it is -- and no JMT drop strategy reproduces that: "waiting queue" does not enforce the size at all, and ', ...
        '"BAS blocking" completes the service before blocking, which is a different queueing model. Use SolverCTMC, ', ...
        'SolverSSA or SolverLDES, or declare DropStrategy.BAS if blocking after service is the model you want, which ', ...
        'SolverJMT does export.'], sn.nodenames{sn.stationToNode(ist)}, sn.cap(ist), sn.classnames{r});
    return
end
end
