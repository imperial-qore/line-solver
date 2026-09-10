function [b, nb] = replyBlocked(sn, ind, space_var)
% [B, NB] = REPLYBLOCKED(SN, IND, SPACE_VAR)
%
% Servers held at node IND by jobs that made a synchronous call and are waiting
% for their REPLY signal. B is (rows x nclasses) per-class counts, one row per
% row of SPACE_VAR, and NB is the (rows x 1) total. Returns zeros when the node
% carries no synchronous-call block, so callers can subtract NB from the server
% count unconditionally.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

R = sn.nclasses;
nrows = max(1, size(space_var,1));
b = zeros(nrows,R);
nb = zeros(nrows,1);
if isempty(space_var) || size(sn.nvars,2) < 2*R+1+R
    return
end
rinfo = State.replyBlockInfo(sn, ind);
if rinfo.width == 0
    return
end
for r = rinfo.classes
    if rinfo.slot(r) <= size(space_var,2)
        b(:,r) = space_var(:, rinfo.slot(r));
    end
end
nb = sum(b,2);
end
