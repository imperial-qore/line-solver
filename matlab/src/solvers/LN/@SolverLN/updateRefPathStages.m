function updateRefPathStages(self, it)
% UPDATEREFPATHSTAGES(SELF, IT)
%
% Recompute the service time of every reference-path stage from the current
% iterate. There are TWO kinds, and they partition the path.
%
% A HOP stage stands for an entry of a task that has no activity graph in the
% layer it appears in: the customer descending from the reference task passes
% through it, and what it must be charged there is
%
%   stage(e) = residt(e)                        the entry's own residence
%            - sum over the re-expanded calls of e of
%                 callresidt(c) / entryvisits(e) the descents the chain now
%                                                walks explicitly, converted
%                                                from the task scale to the
%                                                entry scale
%            + W(task of e)                     the wait its threads queue for
%
% A GATE stage stands for the descent into a MEMBER, a caller whose subgraph IS
% present in this layer. Every residence that member owns is therefore already
% explicit and nothing is charged or subtracted, leaving only
%
%   stage(gate) = W(callee task)                the wait its threads queue for
%
% which is the one part of a member's cost the explicit subgraph cannot show,
% because a task's own thread pool is a station of its own task layer.
%
% W is taken in the POPULATION domain from the task's own layer,
% max(0, Qtot - m*U)/X, where both operands are O(m). The time-domain form
% callservt(c)/callmean(c) - servt(e) must not be used: both of ITS operands
% carry the whole descent to the bottom of the tree and are nearly equal when
% the pool is uncontended, which is the cancellation updateLayersPH forbids
% after it produced RespT 1.4e47 on lqn_sockshop.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(self.refpath_classes_updmap)
    return
end
lqn = self.lqn;
ensemble = self.ensemble;
idxhash = self.idxhash;
rows = self.refpath_classes_updmap;
nkeys = size(self.refpath_keys,1);

terms = zeros(size(rows,1),1);
for r = 1:size(rows,1)
    kind = rows(r,4);
    elem = rows(r,5);
    switch kind
        case 1 % the entry's own residence, per entry invocation
            terms(r) = self.residt(elem);
        case 2 % a descent the chain re-expands, converted to the entry scale
            eidx_from = rows(r,7);
            v = 1;
            if eidx_from >= 1 && self.entryvisits(eidx_from) > GlobalConstants.FineTol
                v = self.entryvisits(eidx_from);
            end
            terms(r) = self.callresidt(elem) / v;
        case 3 % the thread-acquisition wait of the entry's task
            terms(r) = threadWait(self, lqn, elem);
    end
    if ~isfinite(terms(r))
        terms(r) = 0;
    end
    terms(r) = terms(r) * rows(r,6);
end

stage = accumarray(self.refpath_group, terms, [nkeys, 1]);
stage = max(stage, 0);

% Under-relaxation on the same schedule as the think times, including the
% crash-recovery escape: a previous value an order of magnitude above the raw
% one is a divergent iterate, not a smoothing target.
omega = self.relax_omega;
for k = 1:nkeys
    if omega < 1.0 && it > 1 && ~isnan(self.refpath_stage_prev(k))
        prevS = self.refpath_stage_prev(k);
        if prevS > 10 * stage(k) && stage(k) > GlobalConstants.FineTol
            prevS = stage(k);
        end
        stage(k) = omega * stage(k) + (1 - omega) * prevS;
    end
    self.refpath_stage_prev(k) = stage(k);
    e = idxhash(self.refpath_keys(k,1));
    node = ensemble{e}.nodes{self.refpath_keys(k,2)};
    cls = ensemble{e}.classes{self.refpath_keys(k,3)};
    if stage(k) > GlobalConstants.FineTol
        node.setService(cls, Exp.fitMean(stage(k)));
    else
        % Exp.fitMean(0) is a rate of Inf, which sn_refresh_visits reads as a
        % disabled (station, class) pair and drops from the chain
        node.setService(cls, Immediate.getInstance());
    end
end
end

function w = threadWait(self, lqn, tidx)
% Mean time a request queues for a thread of task TIDX, read from TIDX's own
% layer in the population domain. An infinite or unbounded pool never queues.
w = 0;
if tidx < 1 || tidx > lqn.nidx
    return
end
if lqn.sched(tidx) == SchedStrategy.INF || ~isfinite(lqn.maxmult(tidx))
    return
end
[e, k] = self.layerOf(tidx);
if isnan(e) || isnan(k) || size(self.results,1) < 1 || isempty(self.results{end,e})
    return
end
Qtot = sum(self.results{end,e}.QN(k,:), 2);
Btot = full(lqn.maxmult(tidx)) * self.util(tidx);
Xtot = self.tput(tidx) / max(1, lqn.repl(tidx));
if Xtot > GlobalConstants.FineTol
    w = max(GlobalConstants.Zero, Qtot - Btot) / Xtot;
end
end
