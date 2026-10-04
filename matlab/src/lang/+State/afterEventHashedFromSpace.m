function [outhash, outrate, outprob, outstart, outpreempt] = afterEventHashedFromSpace(sn, ind, inspace, event, class)
% [OUTHASH, OUTRATE, OUTPROB, OUTSTART, OUTPREEMPT] = AFTEREVENTHASHEDFROMSPACE(SN, IND, INSPACE, EVENT, CLASS)
%
% State.afterEventHashed with a RAW input state instead of a hash: the
% successors are still hashed, but the predecessor need not be a state of the
% chain.
%
% This exists for the IMMEDIATE-FEEDBACK self-loop (sn.immfeed). Its departure
% half runs with noPromote and so leaves the server idle with jobs still
% waiting, which has no row in the enumerated space -- between events a
% work-conserving station never idles with a queue, so State.spaceGenerator
% never emits it. Hashing that intermediate returns -1 and the whole arc is
% dropped, which is worse than ignoring the flag: the generator loses the
% transitions instead of merely re-queueing the job. The two halves therefore
% compose on the raw vector, and only the state the PAIR ends in is hashed.
%
% See State.immfeedSyncMask for which synchronizations these are.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

outprob = [];
outstart = zeros(0, sn.nclasses);
outpreempt = zeros(0, sn.nclasses);

isSimulation = false;
[outspace, outrate, outprob, ~, outstart, outpreempt] = State.afterEvent(sn, ind, inspace, event, class, isSimulation);
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
