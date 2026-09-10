function updateMetricsPH(self, it)
% UPDATEMETRICSPH(IT) Reconstruct the LQN metrics of a PH encoding
%
% A layer of these methods reports one row per caller task, not one per entry,
% activity and call, so the per-element quantities the rest of SolverLN reads --
% servt, residt, callservt, callresidt, tput -- are recovered analytically from
% the series-parallel weights of the entry workflows.
%
% Serves 'srvn.ph' and 'flat.ph' alike: what a station reports about a caller is
% the same measurement whether the station sits alone in a layer of its own or
% beside every other server, so the reconstruction reads ph.layer{idx} without
% asking which layering built it.
%
% The split is conservative by construction. A station reports a residence time
% R per visit against a service law of mean S, so the queueing inflation R/S is
% attributed to every leaf of that visit in proportion to its own mean: the
% pieces sum back to R exactly.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

lqn = self.lqn;
ph = self.ph;
nidx = lqn.nidx;

self.servt = zeros(nidx,1);
self.residt = zeros(nidx,1);
self.callservt = zeros(lqn.ncalls,1);
self.callresidt = zeros(lqn.ncalls,1);

%% Host layers: the queueing inflation of the processor demand
inflNum = zeros(nidx,1);
inflDen = zeros(nidx,1);
taskTput = zeros(nidx,1);
openTput = zeros(nidx,1);
for hidx = 1:lqn.nhosts
    e = self.idxhash(hidx);
    if isnan(e) || isempty(ph.layer{hidx})
        continue
    end
    L = ph.layer{hidx};
    res = self.results{end,e};
    qs = L.qstations;
    npop = phLayerPop(self, L, hidx);
    for c = L.callers
        k = L.classOfCaller(c);
        X = sum(res.TN(qs,k));
        R = phResidence(sum(res.QN(qs,k)), X, res.RN(qs(1),k));
        f = inflationOf(R, L.svcmeanByClass(k), npop);
        if ~isfinite(X) || X < 0
            X = 0;
        end
        % TOTAL over the replicas. The processor layer of a replicated element
        % models ONE representative replica, so X is one replica's rate and the
        % element's own rate is REPL times it: LDES puts the replicated task of
        % the two-replica probe at Tput 8.30737, exactly its caller's rate,
        % because the caller calls it once per invocation. The matching per
        % replica quantity is ph.xdemand, which the think-time closure divides
        % down for the same reason.
        taskTput(c) = taskTput(c) + lqn.repl(c) * X;
        for eidx = lqn.entriesof{c}
            w = max(ph.share(eidx),0) * X;
            inflNum(eidx) = inflNum(eidx) + w * f;
            inflDen(eidx) = inflDen(eidx) + w;
        end
    end
    for r = 1:size(L.openArrivals,1)
        tag = L.openArrivals(r,2);
        if tag <= 0
            continue % an async call is served in the task layer, not here
        end
        k = L.openArrivals(r,1);
        eidx = tag;
        X = sum(res.TN(qs,k));
        if ~isfinite(X) || X <= 0
            continue
        end
        f = inflationOf(phResidence(sum(res.QN(qs,k)), X, res.RN(qs(1),k)), ...
            L.svcmeanByClass(k), npop);
        inflNum(eidx) = inflNum(eidx) + X * f;
        inflDen(eidx) = inflDen(eidx) + X;
        openTput(eidx) = openTput(eidx) + X;
        taskTput(lqn.parent(eidx)) = taskTput(lqn.parent(eidx)) + X;
    end
end

for e = 1:lqn.nentries
    eidx = lqn.eshift + e;
    f = 1;
    if inflDen(eidx) > GlobalConstants.FineTol
        f = inflNum(eidx) / inflDen(eidx);
    end
    if ~isfinite(f) || f < 1
        f = 1; % a residence time cannot fall below the demand it contains
    end
    for aidx = lqn.actsof{eidx}
        self.residt(aidx) = f * lqn.hostdem_mean(aidx);
    end
end

%% Task layers: the response time of every call
relw = zeros(nidx,1);
for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    e = self.idxhash(tidx);
    if isnan(e) || isempty(ph.layer{tidx})
        continue
    end
    L = ph.layer{tidx};
    res = self.results{end,e};
    qs = L.qstations;
    npop = phLayerPop(self, L, tidx);
    for c = L.callers
        k = L.classOfCaller(c);
        X = sum(res.TN(qs,k));
        if ~isfinite(X) || X < 0
            X = 0;
        end
        g = inflationOf(phResidence(sum(res.QN(qs,k)), X, res.RN(qs(1),k)), ...
            L.svcmeanByClass(k), npop);
        for cidx = syncCallsBetween(lqn, c, tidx)
            eidx = lqn.callpair(cidx,2);
            self.callservt(cidx) = lqn.callproc_mean(cidx) * g * ph.entrymean(eidx);
            self.callresidt(cidx) = self.callservt(cidx);
        end
        for eidx = lqn.entriesof{tidx}
            relw(eidx) = relw(eidx) + X * ph.ncalls(c, eidx);
        end
    end
    for r = 1:size(L.openArrivals,1)
        tag = L.openArrivals(r,2);
        if tag >= 0
            continue
        end
        k = L.openArrivals(r,1);
        cidx = -tag;
        eidx = lqn.callpair(cidx,2);
        X = sum(res.TN(qs,k));
        R = res.RN(qs(1),k);
        if isfinite(R) && R > 0
            self.callservt(cidx) = R * lqn.callproc_mean(cidx);
            self.callresidt(cidx) = self.callservt(cidx);
        end
        if isfinite(X) && X > 0
            relw(eidx) = relw(eidx) + X;
        end
    end
end

%% Entry shares and throughputs
for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    entries = lqn.entriesof{tidx};
    if isempty(entries)
        continue
    end
    % How the requests SPLIT over the entries is a flow-balance question, and is
    % answered at the task layer: a caller class reaches that server once per
    % invocation of the caller, carrying its whole call burst in its service
    % law, so the station rate counts caller cycles and the per-entry rate is
    % that rate times the calls the caller makes.
    w = relw(entries) + openTput(entries);
    tot = sum(w);
    if tot > GlobalConstants.FineTol
        ph.share(entries) = w / tot;
    else
        ph.share(entries) = ones(1,numel(entries)) / numel(entries);
    end
    % HOW MANY requests the task completes is a different question, and the
    % flow-balance total does not answer it: that total is what the callers
    % DEMAND, not what the task's threads can deliver. A thread cycles through
    % its host demand AND then through the task think time, and only the
    % processor layer of the task carries both, so the rate is read there.
    % Reading the demand instead lets a think-time-throttled task report a rate
    % its threads cannot sustain: lqn_basic's T3 has 25 threads and a think time
    % of 4, hence at most 25/(4+0.02) = 6.219 requests per second, whatever its
    % callers ask for.
    if taskTput(tidx) > GlobalConstants.FineTol
        self.tput(tidx) = taskTput(tidx);
    else
        % no processor layer of its own: fall back on the demand
        self.tput(tidx) = tot;
    end
    self.tput(entries) = self.tput(tidx) * ph.share(entries);
    % The DEMAND is kept apart because it, and not the rate just reported, is
    % what closes the surrogate delay: normalising the think time by a rate the
    % same think time produced makes the processor layer self-referential and it
    % settles wherever it started -- see updateThinkTimesPH.
    %
    % PER REPLICA, because the thread count it is paired with there is per
    % replica: `njobs` counts the threads of ONE replica. Closing a per-replica
    % thread cycle with the rate summed over all replicas makes the think time
    % REPL times too small, the processor layer that much too fast, and the
    % replicated task reported 5.86443 where LDES measures 8.30737.
    nrep = max(1, lqn.repl(tidx));
    if tot > GlobalConstants.FineTol
        ph.xdemand(tidx) = tot / nrep;
    else
        ph.xdemand(tidx) = self.tput(tidx) / nrep;
    end
end

%% Recovery, under-relaxation, and the derived per-element quantities
omega = self.relax_omega;
for aidx = (lqn.ashift+1):(lqn.ashift+lqn.nacts)
    v = self.residt(aidx);
    if (isinf(v) || isnan(v)) && it > 1 && ~isnan(self.residt_prev(aidx))
        v = self.residt_prev(aidx);
    end
    if omega < 1.0 && it > 1 && ~isnan(self.residt_prev(aidx))
        v = omega * v + (1 - omega) * self.residt_prev(aidx);
    end
    self.residt(aidx) = v;
    self.residt_prev(aidx) = v;
end
for cidx = 1:lqn.ncalls
    v = self.callservt(cidx);
    if (isinf(v) || isnan(v))
        if it > 1 && isfinite(self.callservt_prev(cidx))
            v = self.callservt_prev(cidx);
        else
            v = 0;
        end
    end
    if omega < 1.0 && it > 1 && ~isnan(self.callservt_prev(cidx))
        v = omega * v + (1 - omega) * self.callservt_prev(cidx);
    end
    self.callservt(cidx) = v;
    self.callresidt(cidx) = v;
    self.callservt_prev(cidx) = v;
    self.callresidt_prev(cidx) = v;
    if v > 0
        self.callservtproc{cidx} = Exp.fitMean(v);
    end
end

self.ph = ph;

% Recompose the entry laws from the iterate just computed. The entry service
% time is then the mean of the COMPOSED law and not the sum of the parts: the
% branches of an AND fork overlap, so an entry that forks finishes with the last
% of its branches and is not charged their sum.
self.phComposeEntryLaws();
ph = self.ph;

for e = 1:lqn.nentries
    eidx = lqn.eshift + e;
    if isempty(ph.execs{eidx})
        continue
    end
    ex = ph.execs{eidx};
    for aidx = lqn.actsof{eidx}
        sa = self.residt(aidx) + lqn_act_thinktime(lqn, aidx);
        for cidx = lqn.callsof{aidx}
            if lqn.calltype(cidx) == CallType.SYNC
                sa = sa + self.callservt(cidx);
            end
        end
        self.servt(aidx) = sa;
        self.servt_prev(aidx) = sa;
        self.tput(aidx) = self.tput(eidx) * ex(aidx);
        self.tput_prev(aidx) = self.tput(aidx);
        % Exp admits rate 0, but a never-called activity has no arrivals: Disabled says so, as the python twin does.
        if self.tput(aidx) > 0
            self.tputproc{aidx} = Exp.fitRate(self.tput(aidx));
        else
            self.tputproc{aidx} = Disabled.getInstance();
        end
        if sa > 0
            self.servtproc{aidx} = Exp.fitMean(sa);
        end
    end
    self.servt(eidx) = ph.entrymean(eidx);
    self.residt(eidx) = ph.entrymean(eidx);
    if self.servt(eidx) > 0
        self.servtproc{eidx} = Exp.fitMean(self.servt(eidx));
    end
end
end

% ------------------------------------------------------------------------
function n = phLayerPop(self, L, idx)
% Total closed population of the MODEL the server sits in, i.e. how many jobs a
% job can queue behind. Under 'flat.ph' that is every caller of the single
% network and not only the callers of this one station, which is why it is taken
% from the layer record rather than recomputed from the callers here.
if isfield(L,'npop') && isfinite(L.npop) && L.npop >= 1
    n = L.npop;
    return
end
n = 0;
for c = L.callers
    v = self.njobs(c, idx);
    if isfinite(v) && v > 0
        n = n + v;
    end
end
if n < 1
    n = 1;
end
end

% ------------------------------------------------------------------------
function R = phResidence(Q, X, RN)
% Residence time per visit, by Little from the queue length rather than from the
% reported RN.
%
% A layer that saturates can come back from AMVA with an RN that no closed model
% can produce -- lqn_sockshop reached RN 1.74e40 at a station whose 1000 jobs and
% service mean 8.23 bound it by 8226 -- and a reconstruction that trusts it feeds
% the impossible value straight back into the call response times. Q is bounded
% by the layer population by construction, so Q/X cannot run away. This is the
% same substitution updateMetricsDefault makes, for the same reason.
R = RN;
if isfinite(Q) && Q >= 0 && isfinite(X) && X > GlobalConstants.FineTol
    R = Q / X;
end
end

% ------------------------------------------------------------------------
function f = inflationOf(R, S, npop)
% Ratio of a residence time to the mean of the law it was measured against, i.e.
% how much the queue at that station stretched the demand.
%
% Bounded above by the layer population: a job can wait behind at most every
% other job in a closed layer, so its residence cannot exceed N service times.
f = 1;
if S > GlobalConstants.FineTol && isfinite(R) && R > 0
    f = R / S;
end
if ~isfinite(f) || f < 1
    f = 1;
end
if nargin >= 3 && isfinite(npop) && npop >= 1 && f > npop
    f = npop;
end
end

% ------------------------------------------------------------------------
function cidxs = syncCallsBetween(lqn, c, tidx)
% Synchronous calls issued by task C to an entry of task TIDX
cidxs = [];
for cidx = 1:lqn.ncalls
    if lqn.calltype(cidx) ~= CallType.SYNC
        continue
    end
    if lqn.parent(lqn.callpair(cidx,1)) == c && lqn.parent(lqn.callpair(cidx,2)) == tidx
        cidxs(end+1) = cidx; %#ok<AGROW>
    end
end
end
