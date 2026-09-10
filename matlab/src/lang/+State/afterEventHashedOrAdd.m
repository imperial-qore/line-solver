function [outhash, outrate, outprob, sn] =  afterEventHashedOrAdd(sn, ind, inhash, event, class)
% [OUTHASH, OUTRATE, OUTPROB, QN] =  AFTEREVENTHASHEDORADD(QN, IND, INHASH, EVENT, CLASS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if inhash == 0
    outhash = -1;
    outrate = 0;
    return
end
% ind: node index
%ist = sn.nodeToStation(ind);
isf = sn.nodeToStateful(ind);

inspace = sn.space{isf}(inhash,:);
isSimulation = true; % allow state vector to grow, e.g. for FCFS buffers
[outspace, outrate, outprob] =  State.afterEvent(sn, ind, inspace, event, class,isSimulation);
if isempty(outspace)
    outhash = -1;
    outrate = 0;
    return
else
    [outhash, sn] = State.getHashOrAdd(sn, ind, outspace);
end
end
