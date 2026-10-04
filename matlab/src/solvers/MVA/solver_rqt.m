function [Q,U,R,T,C,X,lG,totiter] = solver_rqt(sn, options)
% [Q,U,R,T,C,X,lG,totiter] = SOLVER_RQT(SN, OPTIONS)
%
% Robust Queueing Network Analyzer (RQNA) of Robust Queueing Theory. Estimates
% the steady-state performance of a single-class open network of FCFS queues
% with Markovian routing by replacing the stochastic primitives with polyhedral
% uncertainty sets and taking a worst-case view of each node in isolation.
%
% The algorithm is Section 7.2 of the reference: (1) the external streams get
% Gamma_a = sigma_a, (2) the effective arrival process at each node follows from
% the network characterization of Theorem 10 (NPFQN_TRAFFIC_RQT), (3) the
% service variability parameter follows from the adaptation of Section 7.1
% (QSYS_GIGK_RQT_GAMMA), and (4) the system time at each node is the worst-case
% bound of Theorem 3 (QSYS_GIGK_RQT). Step 3 of the published algorithm, the
% path enumeration, is not needed here: LINE aggregates per-node system times
% into per-class response times through the visit ratios.
%
% The adaptation of step (3) is regressed against simulation in heavy traffic,
% so accuracy degrades at low utilization: on M/M/1 the error is about 8% at
% rho=0.9 but over 50% at rho=0.5. Use 'qna' or an exact method for lightly
% loaded models.
%
% Tail coefficients default to alpha=2, the finite-variance regime; set
% options.config.rqt_alpha_a and options.config.rqt_alpha_s (scalars, or vectors
% over stations) in (1,2] to declare heavy-tailed streams, which LINE's
% distributions cannot express since they carry finite moments.
% options.config.rqt_regime selects the Table 1 adaptation regime
% ('independent', 'normal', 'pareto') and options.config.rqt_exact replaces the
% closed-form bound by the exact worst case over the uncertainty sets.
%
% Reference: C. Bandi, D. Bertsimas, N. Youssef (2015). Robust Queueing Theory.
% Operations Research 63(3), 676-700.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% One predicate for the gate and the run: SolverMVA.supportsModelMethod asks
% the same question before the report offers 'rqt', so the sentence a caller
% reads here is the sentence that kept the row off the report.
[rqtOk, rqtReason] = SolverMVA.supportsSingleClassOpen(sn, 'rqt');
if ~rqtOk
    line_error(mfilename, rqtReason);
end
if any(isfinite(sn.njobs))
    line_error(mfilename, 'RQT supports open networks only (no closed classes).');
end

M = sn.nstations;
K = sn.nclasses;   % == 1

Q = zeros(M,K); U = zeros(M,K); R = zeros(M,K); T = zeros(M,K);
X = zeros(1,K); lG = NaN; totiter = 1;

% ----- configuration -----
regime = 'independent';
useExact = false;
alpha_a_cfg = []; alpha_s_cfg = [];
if isfield(options,'config') && isstruct(options.config)
    if isfield(options.config,'rqt_regime'), regime = options.config.rqt_regime; end
    if isfield(options.config,'rqt_exact'), useExact = logical(options.config.rqt_exact); end
    if isfield(options.config,'rqt_alpha_a'), alpha_a_cfg = options.config.rqt_alpha_a; end
    if isfield(options.config,'rqt_alpha_s'), alpha_s_cfg = options.config.rqt_alpha_s; end
end

% ----- identify source and queueing stations -----
isSource = false(M,1);
schedInf = false(M,1);
for i = 1:M
    nd = sn.stationToNode(i);
    isSource(i) = (sn.nodetype(nd) == NodeType.Source);
    schedInf(i) = (sn.sched(i) == SchedStrategy.INF);
end
srcList = find(isSource);
if isempty(srcList)
    line_error(mfilename, 'RQT requires an open network with a Source station.');
end
src = srcList(1);
qstat = find(~isSource);
nq = numel(qstat);

rtS = full(sn.rt);   % station-to-station routing for a single class

% ----- external arrival process -----
arvMAP = sn.proc{src,1}{1};
lambda_src = map_lambda(arvMAP);
sigma_a_src = sqrt(map_scv(arvMAP)) / lambda_src;
alpha_a_src = 2;
if ~isempty(alpha_a_cfg)
    alpha_a_src = alpha_a_cfg(1);
end

% ----- per-node primitives -----
mu = zeros(nq,1); sigma_s = zeros(nq,1); nserv = ones(nq,1);
alpha_s = repmat(2, nq, 1);
lambda0 = zeros(nq,1); Gamma0 = zeros(nq,1); alpha0 = repmat(alpha_a_src, nq, 1);
F = zeros(nq,nq);
for a = 1:nq
    ia = qstat(a);
    mu(a) = sn.rates(ia,1);
    sigma_s(a) = sqrt(sn.scv(ia,1)) / mu(a);
    if isfinite(sn.nservers(ia)) && sn.nservers(ia) > 0
        nserv(a) = sn.nservers(ia);
    end
    if ~isempty(alpha_s_cfg)
        alpha_s(a) = alpha_s_cfg(min(a,numel(alpha_s_cfg)));
    end
    % the source stream reaches node a thinned by q, Theorem 6
    q = rtS(src, ia);
    lambda0(a) = lambda_src * q;
    if q > 0
        Gamma0(a) = sigma_a_src * (1/q)^(1/alpha_a_src);
    end
    for b = 1:nq
        F(a,b) = rtS(ia, qstat(b));
    end
end

% ----- effective arrival processes, Theorem 10 -----
[lambda, Gamma_a, alpha_a] = npfqn_traffic_rqt(lambda0, Gamma0, alpha0, F);

% ----- per-node worst-case analysis -----
for a = 1:nq
    ia = qstat(a);
    T(ia,1) = lambda(a);
    if lambda(a) <= 0
        continue
    end
    if schedInf(ia)
        U(ia,1) = lambda(a)/mu(a);
        Q(ia,1) = lambda(a)/mu(a);
        R(ia,1) = 1/mu(a);
        continue
    end
    rho = lambda(a) / (nserv(a)*mu(a));
    if rho >= 1
        line_warning(mfilename, 'Station %d is unstable (rho=%f), RQT returns an infinite system time.', ia, rho);
    end
    Gamma_s = qsys_gigk_rqt_gamma(rho, mu(a), Gamma_a(a), sigma_s(a), nserv(a), alpha_a(a), regime);
    [Wa, ~, Sworst] = qsys_gigk_rqt(lambda(a), mu(a), Gamma_a(a), Gamma_s, nserv(a), alpha_a(a), alpha_s(a));
    if useExact
        Wa = Sworst;
    end
    R(ia,1) = Wa;
    U(ia,1) = rho;
    Q(ia,1) = lambda(a) * Wa;   % Little's law, number in system
end

T(src,1) = lambda_src;
C = sum(R,1);
X(1) = lambda_src;
Q(isnan(Q)) = 0; U(isnan(U)) = 0; R(isnan(R)) = 0; C(isnan(C)) = 0;
end
