function [outhash, outrate, outprob, outstart, outpreempt] =  afterEventHashed(sn, ind, inhash, event, class)
% [OUTHASH, OUTRATE, OUTPROB, OUTSTART, OUTPREEMPT] =  AFTEREVENTHASHED(QN, IND, INHASH, EVENT, CLASS)
%
% OUTSTART and OUTPREEMPT are the START/PREEMPT annotations of each successor,
% carried alongside the hashed successor because both CTMC and SSA reach the
% state machine through here. See State.afterEvent.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

outprob = [];
outstart = zeros(0, sn.nclasses);
outpreempt = zeros(0, sn.nclasses);

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

[outspace, outrate, outprob, ~, outstart, outpreempt] =  State.afterEvent(sn, ind, inspace, event, class, isSimulation);
if isempty(outspace)
    outhash = -1;
    outrate = 0;
    outstart = zeros(0, sn.nclasses);
    outpreempt = zeros(0, sn.nclasses);
    return
else
    outhash = State.getHash(sn, ind, outspace);
end
end
