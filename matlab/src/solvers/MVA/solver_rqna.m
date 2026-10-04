function [Q,U,R,T,C,X,lG,totiter] = solver_rqna(sn, options)
% [Q,U,R,T,C,X,lG,totiter] = SOLVER_RQNA(SN, OPTIONS)
%
% Robust Queueing Network Analyzer (RQNA) based on indices of dispersion.
% Approximates the steady-state performance of a single-class open queueing
% network of single-server FCFS queues with Markovian routing and general
% (non-renewal) external arrival and (non-exponential) service processes.
%
% Reference: W. Whitt and W. You (2018), "A Robust Queueing Network Analyzer
% Based on Indices of Dispersion", INFORMS J. on Computing. This routine
% implements Algorithm 1 (general framework): traffic-rate equations, limiting
% variability equations, time-dependent IDC equations with default correction
% terms alpha (eq. 34) and beta (eqs. 38-39), and the robust-queueing (RQ)
% workload approximation (eq. 13) with the derived performance measures.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% One predicate for the gate and the run: SolverMVA.supportsModelMethod asks
% the same question before the report offers 'rqna', so the sentence a caller
% reads here is the sentence that kept the row off the report.
[rqnaOk, rqnaReason] = SolverMVA.supportsSingleClassOpen(sn, 'rqna');
if ~rqnaOk
    line_error(mfilename, rqnaReason);
end
if any(isfinite(sn.njobs))
    line_error(mfilename, 'RQNA supports open networks only (no closed classes).');
end

M = sn.nstations;
K = sn.nclasses;   % == 1

Q = zeros(M,K); U = zeros(M,K); R = zeros(M,K); T = zeros(M,K);
X = zeros(1,K); lG = NaN; totiter = 1;

% ----- identify source and queueing stations -----
isSource = false(M,1);
schedInf = false(M,1);
for i = 1:M
    nd = sn.stationToNode(i);
    isSource(i) = (sn.nodetype(nd) == NodeType.Source);
    schedInf(i) = (sn.sched(i) == SchedStrategy.INF);
end
srcList = find(isSource);
qstat = find(~isSource);          % queueing (and delay) stations
nq = numel(qstat);

% single-class station-to-station routing matrix
rtS = full(sn.rt);                % (nstateful x nstateful) for single class

% ----- external arrival process (from the source) -----
if isempty(srcList)
    line_error(mfilename, 'RQNA requires an open network with a Source station.');
end
src = srcList(1);
arvMAP = sn.proc{src,1}{1};
lambda_src = map_lambda(arvMAP);
c2_src = map_idc(arvMAP);
arvIdc = @(t) map_count_idc(arvMAP, t);

% ----- per-queue data -----
mu = zeros(nq,1); cs2 = zeros(nq,1);
lambda0 = zeros(nq,1); c2a0 = zeros(nq,1);
svcMAP = cell(nq,1);
qsplit = zeros(nq,1);             % source split prob into each queue
P = zeros(nq,nq);
for a = 1:nq
    ia = qstat(a);
    mu(a)  = sn.rates(ia,1);
    cs2(a) = sn.scv(ia,1);
    svcMAP{a} = sn.proc{ia,1}{1};
    qsplit(a) = rtS(src, ia);
    lambda0(a) = lambda_src * qsplit(a);
    for b = 1:nq
        ib = qstat(b);
        P(a,b) = rtS(ia, ib);
    end
end

% external arrival IDC seen by each queue = split of the source process (eq. 31)
% I_{a,0,b}(t) = q_b I_src(t) + (1-q_b);  c2_{a,0,b} = q_b c2_src + (1-q_b)
c2a0 = qsplit .* c2_src + (1 - qsplit);
a0IdcFun = @(t) qsplit .* arvIdc(t) + (1 - qsplit);

% service IDC vector-valued function I_{s,a}(t)
sIdcFun = @(t) local_svc_idc(svcMAP, t);

% ----- traffic variability equations (limiting + time-dependent) -----
corrections = struct();
if isfield(options,'config') && isstruct(options.config)
    if isfield(options.config,'rqna_alpha'), corrections.alpha = options.config.rqna_alpha; end
    if isfield(options.config,'rqna_beta'),  corrections.beta  = options.config.rqna_beta;  end
end
ctx = npfqn_traffic_idc(lambda0, P, c2a0, a0IdcFun, mu, cs2, sIdcFun, corrections);
lambda = ctx.lambda;
rho = ctx.rho;

% near-immediate feedback elimination is applied by default (Whitt-You
% Algorithm 2 / Section 4.2); set options.config.rqna_feedback_elim=false for
% the plain RQNA (Algorithm 1) without elimination.
doElim = true;
if isfield(options,'config') && isstruct(options.config) && isfield(options.config,'rqna_feedback_elim')
    doElim = logical(options.config.rqna_feedback_elim);
end

% ----- per-queue robust-queueing workload and performance -----
for a = 1:nq
    ia = qstat(a);
    T(ia,1) = lambda(a);
    if lambda(a) <= 0
        continue;
    end
    if schedInf(ia)
        % infinite-server (delay) station: no waiting
        U(ia,1) = lambda(a)/mu(a);
        Q(ia,1) = lambda(a)/mu(a);
        R(ia,1) = 1/mu(a);
        continue;
    end
    phat = 0;
    if doElim
        phat = local_phat(P, rho, a);
    end
    if phat > options.tol
        % ----- near-immediate feedback elimination at station a -----
        Ra = local_elim_response(a, P, rho, lambda, mu, cs2, svcMAP, ...
            lambda0, arvMAP, qsplit, corrections);
        R(ia,1) = Ra;
    else
        IaFun_a = @(x) local_component(ctx, x, a);
        [~, Wa] = qsys_gig1_rq(rho(a), mu(a), cs2(a), IaFun_a);
        R(ia,1) = Wa + 1/mu(a);   % per-visit response time = waiting + service
    end
    U(ia,1) = rho(a);
    Q(ia,1) = lambda(a) * R(ia,1);% mean number in system (Little, incl. service)
end

% source station bookkeeping
T(src,1) = lambda_src;

C = sum(R,1);
X(1) = lambda_src;
Q(isnan(Q)) = 0; U(isnan(U)) = 0; R(isnan(R)) = 0; C(isnan(C)) = 0;
end

function phat = local_phat(P, rho, a)
% Near-immediate feedback probability at station a: the probability that a
% customer departing a returns to a before visiting any station with strictly
% higher traffic intensity (Whitt-You [52] eq. 3.8/3.9, H={a}).
%
% Delegated to npfqn_feedback_elim so that the solver and the API function
% cannot drift apart: they answer the same question, and a private copy of the
% rule here is how the two came to differ on ties in the first place.
res = npfqn_feedback_elim(P, rho);
phat = res.feedbackProb(a);
end

function Ra = local_elim_response(a, P, rho, lambda, mu, cs2, svcMAP, ...
        lambda0, arvMAP, qsplit, corrections)
% Near-immediate feedback elimination at station a (Whitt-You Algorithm 2 /
% Corollary 4.2). Reduced network: collapse only the equal/lower-rho cloud Hc
% to instantaneous switches (censoring, eq. 3.8), retaining the strictly
% higher-rho stations as real queues so that the global traffic rates and the
% arrival variability from the upstream bottleneck(s) are preserved. The
% near-immediate self-return at a (self-loop phat in the censored network) is
% removed by immediate-feedback elimination (Section 4.1): remove the self-loop,
% renormalize a's remaining routing by 1/(1-phat) and replace a's service by the
% geometric-sum service. The IDC equations are re-solved on the reduced network,
% then RQ is applied at a with the per-visit adjustment W = (1-phat) Wtilde.
Hc = setdiff(find(rho <= rho(a) + 1e-9), a);   % equal/lower-rho cloud (switches)
Hi = setdiff(find(rho >  rho(a) + 1e-9), a);   % strictly higher-rho stations
R  = [a, reshape(Hi,1,[])];                     % retained stations, a first
m  = numel(R);
if isempty(Hc)
    Fhc = []; G = zeros(0,m); Pred = P(R,R);
else
    Fhc = inv(eye(numel(Hc)) - P(Hc,Hc));       % fundamental matrix of the cloud
    G = Fhc * P(Hc,R);                          % first-passage Hc-station -> R
    Pred = P(R,R) + P(R,Hc) * Fhc * P(Hc,R);    % censored routing among R
end
phat = Pred(1,1);                               % near-immediate return prob at a
phat = min(max(phat,0), 1 - 1e-9);

% immediate-feedback elimination at a (position 1 in R)
if phat > 0
    Pred(1,:) = Pred(1,:) / (1 - phat);
end
Pred(1,1) = 0;

% first-passage external arrival rate and IDC into each retained station
arvIdc = @(t) map_count_idc(arvMAP, t);
c2fun_i = @(i,t) qsplit(i)*arvIdc(t) + (1 - qsplit(i));
lam0R = zeros(m,1);
for rr = 1:m
    lam0R(rr) = lambda0(R(rr));
    for ii = 1:numel(Hc)
        lam0R(rr) = lam0R(rr) + lambda0(Hc(ii)) * G(ii,rr);
    end
end
c2a0R = zeros(m,1);
for rr = 1:m
    if lam0R(rr) <= 0, continue; end
    s = lambda0(R(rr)) * (qsplit(R(rr))*map_idc(arvMAP) + (1 - qsplit(R(rr))));
    for ii = 1:numel(Hc)
        g = G(ii,rr);
        ci = qsplit(Hc(ii))*map_idc(arvMAP) + (1 - qsplit(Hc(ii)));
        s = s + lambda0(Hc(ii)) * g * (g*ci + (1-g));
    end
    c2a0R(rr) = s / lam0R(rr);
end
a0IdcR = @(t) local_fp_ext_idc_R(t, R, Hc, G, lambda0, c2fun_i, lam0R);

% service data on R; station a gets the geometric-sum (folded) service
muR = mu(R); cs2R = cs2(R);
svcR = cell(m,1);
for rr = 1:m, svcR{rr} = svcMAP{R(rr)}; end
svcR{1}  = local_geom_map(svcMAP{a}, phat);
muR(1)   = (1-phat) * mu(a);
cs2R(1)  = phat + (1-phat) * cs2(a);
sIdcR = @(t) local_svc_idc(svcR, t);

% re-solve the IDC equations on the reduced network and apply RQ at a
ctxR = npfqn_traffic_idc(lam0R, Pred, c2a0R, a0IdcR, muR, cs2R, sIdcR, corrections);
IaFunA = @(x) local_component(ctxR, x, 1);
[~, Wt] = qsys_gig1_rq(rho(a), muR(1), cs2R(1), IaFunA);
Ra = (1-phat) * Wt + 1/mu(a);   % per-visit adjustment (mean visits 1/(1-phat))
end

function Ir = local_fp_ext_idc_R(t, R, Hc, G, lambda0, c2fun_i, lam0R)
% First-passage external arrival IDC into each retained station of the reduced
% network (superposition of the direct external and the cloud-entering splits).
m = numel(R);
Ir = ones(m,1);
for rr = 1:m
    if lam0R(rr) <= 0, continue; end
    num = lambda0(R(rr)) * c2fun_i(R(rr), t);
    for ii = 1:numel(Hc)
        g = G(ii,rr);
        num = num + lambda0(Hc(ii)) * g * (g*c2fun_i(Hc(ii),t) + (1-g));
    end
    Ir(rr) = num / lam0R(rr);
end
end

function out = local_geom_map(map, p)
% Geometric random sum of i.i.d. PH service times with success prob (1-p):
% PH(alpha, T) -> PH(alpha, T + p*t0*alpha), t0 = -T*e. Yields the head-of-line
% immediate-feedback-eliminated service (Whitt-You Section 4.1).
D0 = map{1};
e  = ones(size(D0,1),1);
t0 = -D0 * e;
al = map_pie(map);
D0m = D0 + p * (t0 * al);
D1m = (1-p) * (t0 * al);
out = {D0m, D1m};
end

function Is = local_svc_idc(svcMAP, t)
% t is either a scalar time (same for all queues) or a length-nq vector giving
% the per-queue evaluation time (service IDC evaluated at rho_a * t).
nq = numel(svcMAP);
Is = zeros(nq,1);
if isscalar(t)
    tt = repmat(t, nq, 1);
else
    tt = t(:);
end
for a = 1:nq
    Is(a) = map_count_idc(svcMAP{a}, tt(a));
end
end

function ia = local_component(ctx, x, a)
v = ctx.IaFun(x);
ia = v(a);
end
