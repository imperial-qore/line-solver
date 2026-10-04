function block = pollingInit(sn, ind, nbuf, srvclass)
% BLOCK = POLLINGINIT(SN, IND, NBUF, SRVCLASS)
%
% Canonical initial value of the polling controller columns for a station whose
% buffers hold NBUF and whose service facility holds a class-SRVCLASS job (0
% when empty). Returns a 1-by-pinfo.width row, empty when the node is not a
% polling station or carries no materialized controller column.
%
% The row is always one that State.pollingBlocks also enumerates, so that the
% initial state is a member of the generated state space.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

block = [];
pinfo = State.pollingInfo(sn, ind);
if isempty(pinfo)
    return
end
R = sn.nclasses;

if srvclass > 0
    % A visit to srvclass is under way. Take it to have started with the whole
    % class present, which is the state reached by a server that has just
    % arrived at a buffer holding nbuf(srvclass)+1 jobs and pulled one of them
    % into service.
    pos = srvclass;
    swk = 0;
    ctr = State.pollingBudget(pinfo, nbuf(srvclass)+1);
else
    % The facility is empty: let the server walk from a canonical position and
    % settle wherever the discipline puts it. With work waiting this cannot
    % return a visit, because initialization always fills the single server of
    % a non-empty station (see initDefault), so only a switchover or a park is
    % reachable here.
    [q, mode] = State.pollingNext(pinfo, 1, nbuf, R, true);
    pos = q;
    ctr = 0;
    switch mode
        case 2
            swk = find(pinfo.swpie{q} > 0, 1, 'first');
        otherwise
            swk = 0;
    end
end

block = State.pollingProject(pinfo, [pos, swk, ctr]);
end
