function [outhash, outrate, outprob, outstart, outpreempt] =  afterEventHashed(sn, ind, inhash, event, class)
% [OUTHASH, OUTRATE, OUTPROB, OUTSTART, OUTPREEMPT] =  AFTEREVENTHASHED(QN, IND, INHASH, EVENT, CLASS)
%
% OUTSTART and OUTPREEMPT are the START/PREEMPT annotations of each successor,
% carried alongside the hashed successor because both CTMC and SSA reach the
% state machine through here. See State.afterEvent.
%
% THE INPUT IS A HASH, SO THE INPUT MUST BE A STATE OF THE CHAIN. The departure
% half of an immediate-feedback self-loop leaves the server idle with jobs still
% waiting, which is NOT a state of the chain -- the enumerated space has no row
% for it, because between events a work-conserving station never idles with a
% queue -- so that pair must compose on the raw vector instead. See
% State.afterEventHashedFromSpace and State.immfeedSyncMask.

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
