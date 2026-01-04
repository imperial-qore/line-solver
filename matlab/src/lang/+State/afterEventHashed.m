function [outhash, outrate, outprob] =  afterEventHashed(sn, ind, inhash, event, class)
% [OUTHASH, OUTRATE, OUTPROB] =  AFTEREVENTHASHED(QN, IND, INHASH, EVENT, CLASS)

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
isSimulation = false;

[outspace, outrate, outprob] =  State.afterEvent(sn, ind, inspace, event, class, isSimulation);
if isempty(outspace)
    outhash = -1;
    outrate = 0;
    return
else
    outhash = State.getHash(sn, ind, outspace);
end
end
