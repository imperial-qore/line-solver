function trips = pollingBlocks(pinfo, srvclass, nbuf, R)
% TRIPS = POLLINGBLOCKS(PINFO, SRVCLASS, NBUF, R)
%
% Every controller configuration compatible with a service facility holding a
% class-SRVCLASS job (0 when empty) and buffers holding NBUF. Returns one
% [pos, swk, ctr] triple per configuration, in full regardless of which columns
% State.pollingInfo materializes (project with State.pollingProject), and zero
% rows when the combination is unoccupiable.
%
% The three tangible configurations of State.pollingInfo map as follows.
%
% SERVING(p): the server stands at the buffer of the job it is serving, so pos
% is pinned to srvclass. The visit budget is free within the bounds its
% discipline can have left it in:
%   GATED        ctr counts the jobs admitted at the polling instant that have
%                not completed, the one in service included, so ctr >= 1; and
%                the ctr-1 still uncompleted ones are all waiting in the
%                buffer, so ctr-1 <= nbuf(p).
%   KLIMITED     1 <= ctr <= K, the one in service counted.
%   DECREMENTING ctr is the population the visit is driving the class down to;
%                it started at one below the level found and the visit is still
%                running, so the current population nbuf(p)+1 exceeds it.
%
% SWITCHING(q): the facility must be empty, and q must be a buffer whose
% switchover takes time, since a zero-time one is never dwelt in. Every phase
% of that switchover is reachable. Jobs may wait meanwhile: this is exactly the
% configuration that makes a polling station non-work-conserving.
%
% PARKED: the facility and the whole station are empty and no switchover takes
% time. With work waiting and only immediate switchovers the server would have
% reached it in zero time, so an idle facility and a non-empty station is
% unoccupiable and yields no rows at all. That pruning matters even when no
% column is materialized: leaving those rows in the state space would make the
% generator reducible.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

trips = zeros(0,3);
if srvclass > 0
    if ~pinfo.polled(srvclass)
        return % a buffer outside the cyclic order can hold no job
    end
    switch pinfo.ptype
        case PollingType.EXHAUSTIVE
            ctrset = 0;
        case PollingType.GATED
            ctrset = 1:(nbuf(srvclass)+1);
        case PollingType.KLIMITED
            ctrset = 1:pinfo.pk;
        case PollingType.DECREMENTING
            ctrset = 0:nbuf(srvclass);
    end
    for ctr = ctrset
        trips(end+1,:) = [srvclass, 0, ctr]; %#ok<AGROW>
    end
else
    for q = find(pinfo.hasSw)
        for swk = 1:pinfo.Ksw(q)
            trips(end+1,:) = [q, swk, 0]; %#ok<AGROW>
        end
    end
    if ~any(pinfo.hasSw) && sum(nbuf) == 0
        trips(end+1,:) = [find(pinfo.polled,1), 0, 0]; % parked, pos canonical
    end
end
end
