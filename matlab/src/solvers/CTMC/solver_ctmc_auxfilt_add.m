function DfiltAux = solver_ctmc_auxfilt_add(DfiltAux, sn, node, s, ns, w, starttag, preempttag, irow)
% DFILTAUX = SOLVER_CTMC_AUXFILT_ADD(DFILTAUX, SN, NODE, S, NS, W, STARTTAG, PREEMPTTAG, IROW)
%
% Accumulate the START/PREEMPT annotation of one successor row into the
% derived filtrations of the station behind NODE. W is the same weight the
% caller added to Dfilt{a}(S,NS), so the filtration integrates rate * count
% and pi*F*e is a rate of starts (or of preemptions) per unit time.
%
% IROW selects the annotated successor row; a node that is not a station, or
% an event that tagged nothing, contributes nothing.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if w == 0 || ~sn.isstation(node)
    return
end
ist = sn.nodeToStation(node);
if ~isempty(starttag) && irow <= size(starttag,1)
    for r = find(starttag(irow,:) ~= 0)
        DfiltAux.start{ist,r}(s,ns) = DfiltAux.start{ist,r}(s,ns) + w * starttag(irow,r);
    end
end
if ~isempty(preempttag) && irow <= size(preempttag,1)
    for r = find(preempttag(irow,:) ~= 0)
        DfiltAux.preempt{ist,r}(s,ns) = DfiltAux.preempt{ist,r}(s,ns) + w * preempttag(irow,r);
    end
end
end
