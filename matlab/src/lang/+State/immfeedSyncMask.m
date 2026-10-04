function mask = immfeedSyncMask(sn, sync)
% MASK = IMMFEEDSYNCMASK(SN, SYNC)
%
% Logical vector, one entry per synchronization, marking the ones that are an
% IMMEDIATE-FEEDBACK SELF-LOOP: a departure whose passive half is an arrival
% back at the same station, in a class the station holds its server for
% (sn.immfeed, set by Queue.setImmediateFeedback or JobClass.setImmediateFeedback).
%
% The departure half of such a synchronization must run with NOPROMOTE, so the
% vacated server is not given to a waiting job and the fed-back arrival seizes
% it instead. Every solver that walks sn.sync must consult this, or it answers
% a model in which the self-looping job re-queues behind the waiting jobs.
%
% The class read is the PASSIVE one: after a class switch the job comes back as
% the class it was switched INTO, and refreshSync folds the switch into sn.rt,
% so a class-switching self-loop is still one synchronization with
% passive.node == active.node.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

A = length(sync);
mask = false(A,1);
if ~isfield(sn,'immfeed') || isempty(sn.immfeed) || ~any(sn.immfeed(:))
    return
end
for a = 1:A
    node_a = sync{a}.active{1}.node;
    node_p = sync{a}.passive{1}.node;
    if sync{a}.active{1}.event ~= EventType.DEP || node_p ~= node_a
        continue
    end
    if node_a < 1 || node_a > sn.nnodes || ~sn.isstation(node_a)
        continue
    end
    ist = sn.nodeToStation(node_a);
    class_p = sync{a}.passive{1}.class;
    if ist >= 1 && ist <= size(sn.immfeed,1) && class_p >= 1 && class_p <= size(sn.immfeed,2)
        mask(a) = logical(sn.immfeed(ist, class_p));
    end
end
end
