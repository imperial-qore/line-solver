function [outspace, outrate, outprob] = afterEventStationReply(sn, ind, ist, class, K, Ks, S, pie, space_buf, space_srv, space_var)
% [OUTSPACE, OUTRATE, OUTPROB] = AFTEREVENTSTATIONREPLY(SN, IND, IST, CLASS, K, KS, S, PIE, SPACE_BUF, SPACE_SRV, SPACE_VAR)
%
% Arrival of a REPLY signal class at the FCFS station that is holding a server
% for the matching synchronous call. The reply completes the call:
%
%   1. one held server of the calling class is released (the reply block
%      counter is decremented);
%   2. if jobs are waiting, the head of line takes the released server, since
%      those jobs queued before the reply arrived (FCFS buffer is right-aligned,
%      head of line = rightmost occupied slot);
%   3. the reply itself then joins like an ordinary arrival: into a server if
%      one is still free, otherwise at the tail of the buffer.
%
% Unlike a NEGATIVE or CATASTROPHE signal the reply is NOT annihilated: it is a
% job that carries the call result onward, so it is served here (typically
% Immediate) and routed on by the ordinary routing matrix. This mirrors LDES,
% where a REPLY arrival frees the blocked server and then continues routing.
%
% The event is passive: the rate is set by the active departure at the callee.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

R = sn.nclasses;
rinfo = State.replyBlockInfo(sn, ind);

% The calling class this reply releases: the class whose expected reply is
% CLASS and which holds a block here. sn.syncreply is stored 0-based.
callclass = 0;
for r = rinfo.classes
    if sn.syncreply(r) + 1 == class
        callclass = r;
        break
    end
end

outspace = [];
outrate = [];
outprob = [];
for row = 1:size(space_srv,1)
    buf = space_buf(row,:);
    srv = space_srv(row,:);
    var = space_var(row,:);

    if callclass > 0 && var(rinfo.slot(callclass)) > 0
        var(rinfo.slot(callclass)) = var(rinfo.slot(callclass)) - 1;
        % PASS-THROUGH: the released server is taken by the REPLY itself, never
        % by a waiting job. The reply is the released work of a call this
        % station already paid for, so it does not queue behind the residents
        % (queueing it gave QLen 0.125 and residence 0.0996 against 0 in LDES,
        % and stole capacity, costing 5% of throughput). Its service is
        % typically Immediate, so the server is handed back at once and the
        % ordinary FCFS departure path then promotes the head of line -- which
        % also keeps sum(srv) <= S, unlike admitting the reply on top of a
        % promoted job.
    end
    [sp, pr] = sub_joinReply(sn, ind, ist, class, K, Ks, S, pie, buf, srv, var);
    outspace = [outspace; sp]; %#ok<AGROW>
    outprob = [outprob; pr]; %#ok<AGROW>
end
outrate = -1*ones(size(outspace,1),1); % passive action
end

function [sp, pr] = sub_joinReply(sn, ind, ist, class, K, Ks, S, pie, buf, srv, var)
% The reply job joins the station: into a free server (enumerating its entry
% phase) or, if all remaining servers are busy, at the tail of the buffer.
[~, nb] = State.replyBlocked(sn, ind, var);
Seff = S(ist) - nb;
sp = [];
pr = [];
if sum(srv) < Seff
    pentry = pie{ist}{class};
    for kentry = 1:K(class)
        if pentry(kentry) <= 0
            continue
        end
        srv_k = srv;
        srv_k(Ks(class)+kentry) = srv_k(Ks(class)+kentry) + 1;
        sp = [sp; buf, srv_k, var]; %#ok<AGROW>
        pr = [pr; pentry(kentry)]; %#ok<AGROW>
    end
    return
end
% All available servers busy: queue at the tail. The FCFS buffer is
% right-aligned with the newest job leftmost, so the tail is the last empty
% slot; grow the buffer by one column when it is full, as the ordinary arrival
% branch of State.afterEventStation does.
emptypos = find(buf == 0, 1, 'last');
if isempty(emptypos)
    buf = [0, buf];
    emptypos = 1;
end
buf(emptypos) = class;
sp = [buf, srv, var];
pr = 1;
end
