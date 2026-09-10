function [QN,UN,RN,TN,AN,WN] = getEnsembleAvgPH(self)
% [QN,UN,RN,TN,AN,WN] = GETENSEMBLEAVGPH(SELF)
%
% LQN-level results of the PH encodings, 'srvn.ph' and 'flat.ph' alike. The
% layers report per caller task, so every entry, activity and call figure is
% rebuilt from the converged fixed point rather than read off a class row, in
% the same layout getEnsembleAvg returns: QN carries the entry and task
% utilizations, UN the processor utilizations, RN the response times and WN the
% residence times. Nothing here reads the ensemble, so the two layerings share
% it verbatim.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Solver console: see getEnsembleAvg. The guard must live until this
% function returns.
consoleGuard = LineConsole.beginRun(self, self.options); %#ok<NASGU>
LineConsole.loop(['solving the layered fixed point over %d layers ' ...
    '(phase-type encoding)'], self.nlayers);
lnRuntime = iterate(self);

lqn = self.lqn;
nidx = lqn.nidx;
QN = nan(nidx,1);
UN = nan(nidx,1);
RN = nan(nidx,1);
TN = nan(nidx,1);
AN = nan(nidx,1);
WN = nan(nidx,1);
PN = nan(nidx,1);   % processor utilization
UT = nan(nidx,1);   % task and entry utilization

for a = 1:lqn.nacts
    aidx = lqn.ashift + a;
    tidx = lqn.parent(aidx);
    if self.ignore(tidx)
        continue
    end
    hidx = lqn.parent(tidx);
    TN(aidx) = self.tput(aidx);
    RN(aidx) = self.servt(aidx);
    UT(aidx) = self.tput(aidx) * self.servt(aidx);
    % LINE scales the utilization of a queueing station into [0,1] whatever its
    % multiplicity, and reports a mean number of busy servers at an infinite
    % server: the processor share of an activity follows the same convention
    PN(aidx) = self.tput(aidx) * lqn.hostdem_mean(aidx) / hostServers(lqn, hidx);
    if isnan(PN(hidx)), PN(hidx) = 0; end
    PN(hidx) = PN(hidx) + PN(aidx);
end

for e = 1:lqn.nentries
    eidx = lqn.eshift + e;
    tidx = lqn.parent(eidx);
    if self.ignore(tidx)
        continue
    end
    TN(eidx) = self.tput(eidx);
    RN(eidx) = self.servt(eidx);
    UT(eidx) = self.tput(eidx) * self.servt(eidx);
    acts = lqn.actsof{eidx};
    if ~isempty(acts)
        PN(eidx) = sum(PN(acts));
    end
    % ResidT is reported per visit to the TASK, not per execution of the
    % activity: an activity of this entry runs EXECS times per invocation, and
    % the entry takes SHARE of the task's invocations. RespT stays per
    % execution. This is the normalization updateMetricsDefault applies through
    % the task/entry throughput ratio.
    if ~isempty(self.ph.execs{eidx})
        ex = self.ph.execs{eidx};
        w = self.ph.share(eidx);
        for aidx = acts
            WN(aidx) = w * ex(aidx) * self.residt(aidx);
        end
    end
    if isnan(UT(tidx)), UT(tidx) = 0; end
    UT(tidx) = UT(tidx) + UT(eidx);
end

for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    if self.ignore(tidx)
        continue
    end
    TN(tidx) = self.tput(tidx);
    acts = lqn.actsof{tidx};
    if ~isempty(acts)
        PN(tidx) = sum(PN(acts));
        WN(tidx) = sum(WN(acts(~isnan(WN(acts)))));
    end
end

for hidx = 1:lqn.nhosts
    TN(hidx) = NaN; % kept NaN for consistency with LQNS
end

% Idle, not undefined -- the same rule getEnsembleAvg applies, and for the
% same reason: an unreachable element reports zero for the measures its kind
% HAS and NaN for the ones it never has, so that the table's NaN mask survives
% a disconnected component. Reported columns here are QLen=UT, Util=PN,
% RespT=RN, ResidT=WN, ArvR=AN, Tput=TN; the pre-swap QN and UN are discarded
% below and are not written.
for idx = find(self.ignore)'
    PN(idx) = 0;   % every kind reports a utilization
    AN(idx) = NaN; % nothing reports an arrival rate on an LQN
    switch lqn.type(idx)
        case LayeredNetworkElement.PROCESSOR
            UT(idx) = NaN; RN(idx) = NaN; WN(idx) = NaN; TN(idx) = NaN;
        case LayeredNetworkElement.TASK
            UT(idx) = 0;   RN(idx) = NaN; WN(idx) = 0;   TN(idx) = 0;
        case LayeredNetworkElement.ENTRY
            UT(idx) = 0;   RN(idx) = 0;   WN(idx) = NaN; TN(idx) = 0;
        case LayeredNetworkElement.ACTIVITY
            UT(idx) = 0;   RN(idx) = 0;   WN(idx) = 0;   TN(idx) = 0;
    end
end

QN = UT;
UN = PN;

% Closing banner, the line every NetworkSolver prints. SolverLN is an
% EnsembleSolver and never reaches NetworkSolver.setAvgResults, so it had none.
self.reportCompletion(lnRuntime);
end

% ------------------------------------------------------------------------
function c = hostServers(lqn, hidx)
% Divisor that scales a processor utilization into [0,1]. An infinite server
% reports a mean number of busy servers instead, so it divides by one.
c = 1;
if lqn.sched(hidx) == SchedStrategy.INF
    return
end
m = lqn.maxmult(hidx);
if isfinite(m) && m > 0
    c = m;
end
end
