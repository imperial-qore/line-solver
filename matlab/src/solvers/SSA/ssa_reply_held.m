function held = ssa_reply_held(sn, tranSysState)
% HELD = SSA_REPLY_HELD(SN, TRANSYSSTATE)
% Time-averaged servers held for a pending REPLY signal, per (station, calling class).
%
% A job of a class r with SN.SYNCREPLY(r) >= 0 leaves its caller for the callee
% but KEEPS its server there until the REPLY returns. The hold is counted in the
% reply block at the tail of the caller's local state (State.replyBlockInfo), not
% in its buffer or servers, so State.toMarginal does not see it. SolverCTMC adds
% the same block mean to QLen and Util (solver_ctmc_analyzer); this is its
% sample-path twin, read from the recorded trajectory: TRANSYSSTATE{1} is the
% cumulative time and TRANSYSSTATE{1+isf} the state of stateful node isf held
% over each step. HELD is all-zero for a model without synchronous calls.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;
K = sn.nclasses;
held = zeros(M,K);
if ~isfield(sn,'replyblock') || isempty(sn.replyblock) || ~any(sn.replyblock(:) > 0) ...
        || ~iscell(tranSysState) || numel(tranSysState) < 2 || isempty(tranSysState{1})
    return
end
cumt = tranSysState{1}(:);
dt = [cumt(1); diff(cumt)];
T = sum(dt);
if T <= 0
    return
end
for ist = 1:M
    ind = sn.stationToNode(ist);
    if size(sn.replyblock,1) < ind || ~any(sn.replyblock(ind,:) > 0)
        continue
    end
    rinfo = State.replyBlockInfo(sn, ind);
    if rinfo.width == 0
        continue
    end
    isf = sn.nodeToStateful(ind);
    X = tranSysState{1+isf};
    % the block trails the row; the buffer may have grown on the left only
    bcols = (size(X,2)-rinfo.width+1):size(X,2);
    bmean = (dt' * X(:,bcols)) / T;
    pos = 0;
    for r = rinfo.classes
        pos = pos + 1;
        held(ist,r) = bmean(pos);
    end
end
end
