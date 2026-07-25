function rinfo = replyBlockInfo(sn, ind)
% RINFO = REPLYBLOCKINFO(SN, IND)
%
% Layout of the synchronous-call (REPLY signal) blocked-server block that node
% IND carries in its local-variable state, as a struct with fields
%
%   classes : calling classes that can hold a server at this node
%   off     : offset of the block inside the node's local-variable vector
%   slot    : 1 x nclasses index into the local-variable vector, 0 when the
%             class holds no slot here
%   width   : number of columns in the block
%
% A job of a class R with SN.SYNCREPLY(R) >= 0 makes a synchronous call: it
% leaves this station for the callee but KEEPS its server, which stays held
% until the matching REPLY signal class comes back here. The held servers are
% not derivable from the marginal state (the job is at the callee, not here),
% so they are counted per class in this block. LDES keys the same information
% by job id (Solver_ssj.pendingReplyMap); a CTMC has no job identity, so it
% carries counts.
%
% The block trails the modulation, routing and node blocks of SN.NVARS
% (columns 1..R, R+1..2R and 2R+1), occupying columns 2R+1+r. Appending keeps
% every existing nvars reader valid, and the columns stay zero-width for models
% without reply signals, so no other model changes state width.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

R = sn.nclasses;
rinfo.classes = [];
rinfo.slot = zeros(1,R);
rinfo.width = 0;
rinfo.off = 0;
if size(sn.nvars,2) < 2*R+1+R
    return
end
rinfo.off = sum(sn.nvars(ind,1:(2*R+1)));
pos = rinfo.off;
for r=1:R
    if sn.nvars(ind, 2*R+1+r) > 0
        pos = pos + 1;
        rinfo.slot(r) = pos;
        rinfo.classes(end+1) = r; %#ok<AGROW>
        rinfo.width = rinfo.width + 1;
    end
end
end
