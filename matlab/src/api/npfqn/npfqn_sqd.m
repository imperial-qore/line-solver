function [X,Q,U,R] = npfqn_sqd(sn, N, calibrationMode, serverBlockingTime, neighborMode, v1Policy, initialV1)
% [X,Q,U,R] = NPFQN_BAS(SN, N, CALIBRATIONMODE, SERVERBLOCKINGTIME, NEIGHBORMODE, V1POLICY, INITIALV1)
%
% Smith Queue Decomposition (SQD): approximate MVA for closed Blocking-After-Service (BAS) networks.
%
% Solves a finite-buffer closed queueing network under Blocking-After-Service
% (manufacturing/transfer) blocking directly from its NetworkStruct. The method is an
% AMVA-style population recursion in which each finite-capacity station is described by
% a load-dependent effective service rate calibrated from an M/M/1/K blocking
% probability; downstream blocking is propagated through the effective routing between
% service stations. Delay (INF/EXT) stations are treated as infinite-capacity pure-delay
% nodes. Single-chain (chain-aggregated) demands only.
%
% Returns per-station throughput X, queue length Q, utilization U, residence time R.
%
% Originally contributed as SolverDBT by Avinash Bommareddy (Imperial College London
% FYP, 2026); refactored here into an sn-based API function.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

INITIAL_V1 = 692.192;
CALIBRATION_EPSILON = 0.05;

M = sn.nstations;
Rcl = sn.nclasses;

if nargin < 2 || isempty(N), N = sn.nclosedjobs; end
if nargin < 3 || isempty(calibrationMode), calibrationMode = 0; end
if nargin < 4 || isempty(serverBlockingTime), serverBlockingTime = true; end
if nargin < 5 || isempty(neighborMode), neighborMode = 'downstream'; end
if nargin < 6 || isempty(v1Policy), v1Policy = 'compound'; end
if nargin < 7, initialV1 = []; end

[~,STchain,Vchain] = sn_get_demands_chain(sn);

V = Vchain(:,1);
ST = STchain(:,1);

isDelay = false(M,1);
cap = zeros(M,1);
for i = 1:M
    s = sn.sched(i);
    isDelay(i) = (s == SchedStrategy.INF || s == SchedStrategy.EXT);
    if isDelay(i)
        cap(i) = Inf;
    else
        c = sn.cap(i);
        if isinf(c) || c > 1e14
            cap(i) = Inf;
        else
            cap(i) = c;
        end
    end
end

pEff = computeEffectiveRouting(sn, M, Rcl, isDelay);

V1 = zeros(M,1);
v1init = zeros(M,1);
L_buf = zeros(M,1);
L_svr = zeros(M,1);
for i = 1:M
    if isempty(initialV1)
        v1init(i) = INITIAL_V1;
    else
        v1init(i) = initialV1(i);
    end
    V1(i) = v1init(i);
end

X = 0.0;
W_buf = zeros(M,1);
W_svr = zeros(M,1);

ownserver = strcmpi(neighborMode, 'ownserver');
fresh = strcmpi(v1Policy, 'fresh');

for pop = 1:N

    % wait times
    for i = 1:M
        if isDelay(i)
            W_buf(i) = 0.0;
            W_svr(i) = ST(i);
            continue;
        end

        if ownserver
            if ~isinf(cap(i))
                pBlockDown = mm1kBlocking(cap(i), X * V(i) * ST(i));
            else
                pBlockDown = 0.0;
            end
        else
            pBlockDown = 0.0;
            for j = 1:M
                if ~isDelay(j) && pEff(i,j) > 0 && ~isinf(cap(j))
                    rho_j = X * V(j) * ST(j);
                    pBlockDown = pBlockDown + pEff(i,j) * mm1kBlocking(cap(j), rho_j);
                end
            end
        end

        if ownserver
            n = L_svr(i);
        else
            n = 0.0;
            for j = 1:M
                n = n + pEff(i,j) * L_svr(j);
            end
        end

        if n < 1e-10
            mu_n = V1(i);
        else
            [beta, gamma] = computeBetaGamma(cap(i), pBlockDown, calibrationMode, CALIBRATION_EPSILON);
            base = max(0.0, (n - 1.0) / beta);
            expArg = base ^ gamma;
            mu_n = n * V1(i) * exp(-expArg);   % Eq.13
        end
        if mu_n < 1e-10, mu_n = 1e-10; end

        W_buf(i) = (1.0 / mu_n) * (1.0 + n);   % Eq.18
        W_svr(i) = ST(i) * (1.0 + L_svr(i));   % Eq.17

        if serverBlockingTime
            bt = 0.0;
            for j = 1:M
                if ~isDelay(j) && pEff(i,j) > 0 && ~isinf(cap(j))
                    rho_j = X * V(j) * ST(j);
                    pBj = mm1kBlocking(cap(j), rho_j);
                    denom = ST(i) + ST(j);
                    if denom > 1e-15
                        theta = ST(j) / denom;
                    else
                        theta = 0.0;
                    end
                    bt = bt + pEff(i,j) * pBj * ST(j) * theta;
                end
            end
            W_svr(i) = W_svr(i) + bt;
        end
    end

    % throughput
    sumVW = sum(V .* (W_buf + W_svr));
    if sumVW > 1e-15
        X = pop / sumVW;
    else
        X = 0.0;
    end

    % queue lengths
    L_buf = X .* V .* W_buf;
    L_svr = X .* V .* W_svr;

    % adjust V1
    if pop < N
        for i = 1:M
            if isDelay(i), continue; end
            if ownserver
                if ~isinf(cap(i))
                    pBlock = mm1kBlocking(cap(i), X * V(i) * ST(i));
                else
                    pBlock = 0.0;
                end
            else
                pBlock = 0.0;
                for j = 1:M
                    if ~isDelay(j) && pEff(i,j) > 0 && ~isinf(cap(j))
                        rho_j = X * V(j) * ST(j);
                        pBlock = pBlock + pEff(i,j) * mm1kBlocking(cap(j), rho_j);
                    end
                end
            end
            if fresh
                V1(i) = v1init(i) * (1.0 - pBlock);
            else
                V1(i) = V1(i) * (1.0 - pBlock);
            end
            if V1(i) < 1e-10, V1(i) = 1e-10; end
        end
    end
end

XN = zeros(M,1);
QN = zeros(M,1);
UN = zeros(M,1);
WN = zeros(M,1);
for i = 1:M
    T_i = X * V(i);
    Q_i = L_buf(i) + L_svr(i);
    W_tot = W_buf(i) + W_svr(i);
    if T_i > 1e-15
        R_i = Q_i / T_i;
    else
        R_i = W_tot;
    end
    U_i = min(1.0, T_i * ST(i));
    XN(i) = T_i;
    QN(i) = Q_i;
    UN(i) = U_i;
    WN(i) = R_i;
end

X = XN;
Q = QN;
U = UN;
R = WN;

end

function [beta, gamma] = computeBetaGamma(K, pBlockDown, mode, CALIBRATION_EPSILON)
% Calibrate (beta, gamma) of the load-dependent effective service rate.
if K <= 2 || isinf(K)
    beta = K; gamma = 1.0; return;
end

a = 2.0;
b = K;

switch mode
    case 2   % blocking-aware
        pK = max(1e-6, min(pBlockDown, 1.0 - 1e-6));
        Va = 1.0 - pK * (a - 1.0) / (b - 1.0);
        Vb = 1.0 - pK;
    case 1   % fixed heuristic
        Va = (b - a) / b;
        Vb = CALIBRATION_EPSILON;
    otherwise % mode 0: base
        beta = K; gamma = 1.0; return;
end

Va = min(Va, 0.999);
Va = max(Va, 0.01);
Vb = max(Vb, 1e-6);
if Vb >= Va
    beta = K; gamma = 1.0; return;
end

lnVa = log(Va);
lnVb = log(Vb);
gamma = log(lnVa / lnVb) / log((a - 1.0) / (b - 1.0));
gamma = max(0.5, min(gamma, 10.0));
beta = (a - 1.0) / (-lnVa) ^ (1.0 / gamma);

if ~isfinite(gamma) || ~isfinite(beta) || beta <= 0
    beta = K; gamma = 1.0; return;
end
end

function p = computeEffectiveRouting(sn, M, Rcl, isDelay)
% Build station-to-station effective routing, collapsing pass-through delay nodes.
p = zeros(M, M);
for i = 1:M
    if isDelay(i), continue; end
    sf_i = sn.stationToStateful(i);
    for j = 1:M
        sf_j = sn.stationToStateful(j);
        p_ij = sn.rt((sf_i - 1) * Rcl + 1, (sf_j - 1) * Rcl + 1);
        if p_ij <= 0, continue; end
        if ~isDelay(j)
            p(i,j) = p(i,j) + p_ij;
        else
            for k = 1:M
                if ~isDelay(k)
                    sf_k = sn.stationToStateful(k);
                    p_jk = sn.rt((sf_j - 1) * Rcl + 1, (sf_k - 1) * Rcl + 1);
                    if p_jk > 0
                        p(i,k) = p(i,k) + p_ij * p_jk;
                    end
                end
            end
        end
    end
end
end

function pb = mm1kBlocking(K, rho)
% Steady-state blocking probability of an M/M/1/K queue at load rho.
if rho <= 1e-15
    pb = 0.0;
elseif abs(rho - 1.0) < 1e-9
    pb = 1.0 / (K + 1.0);
else
    pb = (1.0 - rho) * rho ^ K / (1.0 - rho ^ (K + 1));
end
end
