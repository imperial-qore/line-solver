function tf = sn_has_immfeed(sn)
% TF = SN_HAS_IMMFEED(SN)
% Whether immediate feedback is EFFECTIVE anywhere in the model.
%
% A declaration alone is not enough. Queue.setImmediateFeedback marks a station
% and JobClass.setImmediateFeedback marks a class -- and the class spelling
% marks EVERY station, since sn.immfeed is the OR of the two -- so a model with
% no self-loop at all can carry a full sn.immfeed matrix while the feature
% changes nothing. Reading the raw matrix made every solver that consults it
% warn, or refuse, on a plain M/M/1 that merely mentioned the flag.
%
% Immediate feedback is effective at (station i, class r) when sn.immfeed(i,r)
% holds AND the routing table has a self-loop INTO (i,r) from some class s at
% the same station, which is the only way a job can come back to the server it
% just left. A class switch on the way round is folded into sn.rt by
% refreshRouting, so the incoming class s need not be r.
%
% Solvers that handle immediate feedback look at the SYNCHRONIZATION instead;
% see State.immfeedSyncMask, which applies the same test per sync.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tf = false;
if ~isfield(sn,'immfeed') || isempty(sn.immfeed) || ~any(sn.immfeed(:))
    return
end
% Without a routing table there is nothing to qualify the declaration with, so
% report it as declared rather than silently dropping it.
if ~isfield(sn,'rt') || isempty(sn.rt)
    tf = true;
    return
end
R = sn.nclasses;
for ist = 1:min(sn.nstations, size(sn.immfeed,1))
    if ~any(sn.immfeed(ist,:))
        continue
    end
    ind = sn.stationToNode(ist);
    if ind < 1 || ind > sn.nnodes || ~sn.isstateful(ind)
        continue
    end
    isf = sn.nodeToStateful(ind);
    for r = 1:min(R, size(sn.immfeed,2))
        if ~sn.immfeed(ist,r)
            continue
        end
        col = (isf-1)*R + r;
        if col > size(sn.rt,2)
            continue
        end
        rows = (isf-1)*R + (1:R);
        rows = rows(rows <= size(sn.rt,1));
        if any(sn.rt(rows, col) > 0)
            tf = true;
            return
        end
    end
end
end
