function phComposeEntryLaws(self)
% PHCOMPOSEENTRYLAWS() Recompose the entry service laws of method 'srvn.ph'
%
% Sets every leaf of each entry workflow from the current fixed-point iterate --
% an activity leaf to its processor residence time, a call leaf to the geometric
% compound of the call response law -- and recomposes the entry law by the
% series-parallel reduction. Only the path from a changed leaf to the root is
% recomposed; a subtree whose leaves did not move is reused from the cache.
%
% An activity leaf keeps its shape and its order, because only its mean moves
% (setActivityDemandMean rescales the law in time). A call leaf is refitted from
% two moments, so its order can change and its ancestors are resized: reuse is at
% subtree granularity, not at matrix-entry granularity.
%
% The composed mean is NOT the sum of the leaf means when the graph forks: the
% branches of an AND fork overlap, and the entry finishes with the last of them.
% The ratio of the two, the overlap factor, is what the caller-side aggregates
% are scaled by, so that the pieces of a cycle still add up to the cycle.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

lqn = self.lqn;
ph = self.ph;

%% Entry laws, and the overlap between the branches of their AND forks
overlap = ones(lqn.nidx,1);
setupshare = zeros(lqn.nidx,1);
for e = 1:lqn.nentries
    eidx = lqn.eshift + e;
    if isempty(ph.wf{eidx})
        continue
    end
    wf = ph.wf{eidx};
    ex = ph.execs{eidx};
    entrysum = 0;
    procsum = 0;
    for aidx = lqn.actsof{eidx}
        m = self.residt(aidx) + lqn_act_thinktime(lqn, aidx);
        procsum = procsum + ex(aidx) * m;
        wf.setActivityDemandMean(lqn.names{aidx}, max(m, GlobalConstants.FineTol));
        for cidx = lqn.callsof{aidx}
            if lqn.calltype(cidx) ~= CallType.SYNC
                continue
            end
            wf.setActivityDemand(lqn.callhashnames{cidx}, callBurstLaw(self, lqn, ph, cidx));
            m = m + self.callservt(cidx);
        end
        entrysum = entrysum + ex(aidx) * m;
    end
    law = wf.refreshPH();
    alpha = reshape(law.getInitProb(),1,[]);
    T = law.getSubgenerator();
    [m1, scv] = lqn_ph_moments(alpha, T);
    % All activities of an entry run on ONE processor, so the branches of an
    % AND fork cannot overlap the processor residence they request: the
    % composed maximum is a lower bound on the entry service time only above
    % that total. Where it falls below, the law is rescaled in time to it,
    % which keeps its shape, its SCV and its order.
    if procsum > m1 + GlobalConstants.FineTol
        T = T * (m1 / procsum);
        m1 = procsum;
    end
    % A SetupTask powers a thread down when it goes idle, so a request may find
    % it off and pay a cold start before the entry runs at all. The setup is not
    % part of the activity graph and never enters the series-parallel reduction:
    % it is prefixed to the composed law afterwards, as the mixture
    %   p * (setup THEN entry) + (1-p) * entry
    % which is again phase-type. See phSetupProb for p.
    p = phSetupProb(self, lqn, eidx);
    if p > GlobalConstants.FineTol
        [alpha_s, T_s] = phSetupLaw(lqn, lqn.parent(eidx));
        if ~isempty(T_s)
            [alpha_c, T_c] = Workflow.composeSerial(alpha_s, T_s, alpha, T);
            [alpha, T] = Workflow.composeMixture({alpha_c, alpha}, {T_c, T}, [p, 1-p]);
            [m1, scv] = lqn_ph_moments(alpha, T);
            % The share of the entry law that is cold start and not work. The
            % surrogate-delay closure measures a thread's cycle in WORK, so it
            % must not read a station utilization that this has inflated --
            % see updateThinkTimesPH.
            setupshare(eidx) = p * phSetupMean(lqn, 'setuptime', lqn.parent(eidx)) / max(m1, GlobalConstants.FineTol);
        end
    end
    ph.entryalpha{eidx} = alpha;
    ph.entryT{eidx} = T;
    ph.entrymean(eidx) = m1;
    ph.entryscv(eidx) = scv;
    if entrysum > GlobalConstants.FineTol
        overlap(eidx) = min(1, m1 / entrysum);
    end
end
ph.overlap = overlap;
% Per task, the share-weighted fraction of its station service that is cold
% start rather than work.
ph.setupshare = zeros(lqn.nidx,1);
for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    if self.ignore(tidx)
        continue
    end
    for eidx = lqn.entriesof{tidx}
        ph.setupshare(tidx) = ph.setupshare(tidx) + ph.share(eidx) * setupshare(eidx);
    end
end

%% Expected number of calls per invocation, and the caller-side aggregates
ph.ncalls = zeros(lqn.nidx, lqn.nidx);
ph.calltime = zeros(lqn.nidx, lqn.nidx); % [caller task, called task] time per invocation
ph.procresid = zeros(lqn.nidx, 1);
ph.actthinkt = zeros(lqn.nidx, 1);
ph.calltotal = zeros(lqn.nidx, 1);

for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    if self.ignore(tidx)
        continue
    end
    for eidx = lqn.entriesof{tidx}
        w = ph.share(eidx);
        if w <= 0 || isempty(ph.execs{eidx})
            continue
        end
        ex = ph.execs{eidx};
        r = overlap(eidx);
        for aidx = lqn.actsof{eidx}
            ph.procresid(tidx) = ph.procresid(tidx) + w * r * ex(aidx) * self.residt(aidx);
            ph.actthinkt(tidx) = ph.actthinkt(tidx) + w * r * ex(aidx) * lqn_act_thinktime(lqn, aidx);
            for cidx = lqn.callsof{aidx}
                if lqn.calltype(cidx) ~= CallType.SYNC
                    continue
                end
                tgte = lqn.callpair(cidx,2);
                tgtt = lqn.parent(tgte);
                % the COUNT of calls does not change with the overlap, only the
                % time the caller is held by them
                ph.ncalls(tidx, tgte) = ph.ncalls(tidx, tgte) + w * ex(aidx) * lqn.callproc_mean(cidx);
                ph.calltime(tidx, tgtt) = ph.calltime(tidx, tgtt) + w * r * ex(aidx) * self.callservt(cidx);
                ph.calltotal(tidx) = ph.calltotal(tidx) + w * r * ex(aidx) * self.callservt(cidx);
            end
        end
    end
end

self.ph = ph;
end

% ------------------------------------------------------------------------
function law = callBurstLaw(self, lqn, ph, cidx)
% Law of the total time one execution of the issuing activity spends in call
% CIDX: the geometric compound, of mean LQN.CALLPROC_MEAN, of the response law
% of the called entry. The response law is fitted to the response time reported
% by the callee's layer and to the SCV of the callee's own composed law, so no
% extra solver output is needed.
m = lqn.callproc_mean(cidx);
eidx = lqn.callpair(cidx,2);
if m <= GlobalConstants.FineTol
    law = Immediate.getInstance();
    return
end
R = self.callservt(cidx) / m;
scv = ph.entryscv(eidx);
if ~isfinite(scv) || scv <= GlobalConstants.FineTol
    scv = 1.0;
end
base = APH.fitMeanAndSCV(max(R, GlobalConstants.FineTol), scv);
[alpha, T] = Workflow.composeLoopGeometric(base.getInitProb(), base.getSubgenerator(), m);
if Workflow.isAcyclicGenerator(T)
    law = APH(alpha, T);
else
    law = PH(alpha, T);
end
end

function p = phSetupProb(self, lqn, eidx)
% Probability that a request for entry EIDX finds its task's thread powered off.
%
% A thread is released at a reply and starts a delay-off countdown D of mean d;
% it powers off if D expires before the next request arrives, and a request that
% arrives first cancels the countdown and pays nothing. With the idle interval I
% seen by one thread,
%
%   p = P(D < I),  and for exponential D,  p = E[I] / (E[I] + d).
%
% The estimate of E[I] lives in lqn_setup_charge, which method 'default' calls
% too, so the two methods charge the same cold start.
p = 0;
tidx = lqn.parent(eidx);
if ~isfield(lqn,'hassetup') || isempty(lqn.hassetup) || tidx > numel(lqn.hassetup) ...
        || ~full(lqn.hassetup(tidx))
    return
end
d = phSetupMean(lqn, 'delayofftime', tidx);
s = phSetupMean(lqn, 'setuptime', tidx);
if ~(d > GlobalConstants.FineTol) || ~(s > GlobalConstants.FineTol)
    return
end
% ONE closure for both methods: lqn_setup_charge returns p*s, so p is that over s.
% It also answers p = 1 during construction, before the first solve has sized
% tput or util: nothing has arrived, so the thread is down when the first does.
p = min(1, max(0, lqn_setup_charge(self, lqn, tidx) / s));
end

function [alpha, T] = phSetupLaw(lqn, tidx)
% Phase-type law of task TIDX's setup time, empty when it declares none.
alpha = [];
T = [];
if ~isfield(lqn,'setuptime') || isempty(lqn.setuptime) || tidx > numel(lqn.setuptime)
    return
end
proc = lqn.setuptime{tidx};
if isempty(proc) || ~isa(proc,'Distribution')
    return
end
m = proc.getMean();
if ~isfinite(m) || m <= GlobalConstants.FineTol
    return
end
scv = proc.getSCV();
if ~isfinite(scv) || scv <= GlobalConstants.FineTol
    scv = 1.0;
end
law = APH.fitMeanAndSCV(m, scv);
alpha = reshape(law.getInitProb(),1,[]);
T = law.getSubgenerator();
end

function m = phSetupMean(lqn, fieldname, tidx)
% Mean of task TIDX's setup or delay-off time, 0 when it declares none.
m = 0;
if ~isfield(lqn, fieldname) || isempty(lqn.(fieldname)) || tidx > numel(lqn.(fieldname))
    return
end
proc = lqn.(fieldname){tidx};
if isempty(proc) || ~isa(proc,'Distribution')
    return
end
m = proc.getMean();
if ~isfinite(m), m = 0; end
end
