function [Qt, Ut, Tt] = solver_mam_transient_qbd(sn, options)
% SOLVER_MAM_TRANSIENT_QBD Transient analysis of a single-class open queue via
% the Laplace-domain transient QBD method plus numerical inverse Laplace.
%
% Supports single-server MAP/MAP/1 (infinite buffer) and MAP/MAP/1/N (finite
% buffer), where arrival and service are read uniformly as (D0,D1) MAPs from
% sn.proc; this subsumes M/M/1, M/PH/1 and correlated-MAP arrival/service that
% the libQBD/expm fast path (solver_mam_ldqbd_transient) cannot represent.
%
% The transient level-to-level transform V(s,0,m) is computed by
% mam_transient2_open (infinite) or mam_transient2 (finite) and inverted with
% the CME-based matrix inverse Laplace transform. Metrics returned in the
% [metric, time] layout used by getTranAvg.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;
K = sn.nclasses;
if K ~= 1
    line_error(mfilename, 'Transient QBD method requires a single-class model.');
end
N = sn.njobs';
if ~isinf(N)
    line_error(mfilename, 'Transient QBD method requires an open model.');
end

sourceIdx = find(sn.sched == SchedStrategy.EXT);
queueIdx  = find(sn.sched == SchedStrategy.FCFS);
if numel(sourceIdx) ~= 1 || numel(queueIdx) ~= 1
    line_error(mfilename, 'Transient QBD method requires exactly one Source and one FCFS Queue.');
end
if sn.nservers(queueIdx) ~= 1
    line_error(mfilename, 'Transient QBD (Laplace) method supports single-server queues only.');
end

%% Arrival and service MAPs (uniform (D0,D1) representation)
arrProc = sn.proc{sourceIdx}{1};
svcProc = sn.proc{queueIdx}{1};
Da0 = arrProc{1}; Da1 = arrProc{2};
Ds0 = svcProc{1}; Ds1 = svcProc{2};
na = size(Da0, 1);
ns = size(Ds0, 1);
Ina = eye(na); Ins = eye(ns);

% Repeating-level blocks (levels >= 1): matches qbd_rg convention.
Lrep = krons(Da0, Ds0);          % local, both phases evolve
Frep = kron(Da1, Ins);           % arrival, level up
Brep = kron(Ina, Ds1);           % service completion, level down

% Level-0 (empty) block: service phase frozen, only arrival evolves.
Lv0  = kron(Da0, Ins);
F0   = kron(Da1, Ins);
B0   = kron(Ina, Ds1);

% Initial distribution: empty system, stationary arrival and service phases.
pi_arr = map_prob({Da0, Da1});
pi_svc = map_prob({Ds0, Ds1});
pi0 = kron(pi_arr, pi_svc);      % 1 x (na*ns), level-0 phase distribution

% Service-completion rate vector over (arrival,service) phases at a busy level.
sExit = Ds1 * ones(ns, 1);
wDep  = kron(ones(na, 1), sExit);   % departure-rate weight, (na*ns) x 1
wOne  = ones(na * ns, 1);

bufCap = sn.cap(queueIdx);
isFinite = ~isinf(bufCap);

%% Time grid and ILT budget
T_start = options.timespan(1);
T_end   = options.timespan(2);
dur = T_end - T_start;
nTimePoints = min(101, max(11, round(dur * 10)));
times = linspace(T_start, T_end, nTimePoints)';
% ILT singular at t=0, fill empty-system IC directly; see _kb/06-solver-catalog.md for rationale
posMask = times > 0;
tpos = times(posMask);
if isfield(options, 'iter_max') && ~isempty(options.iter_max) && options.iter_max > 1
    maxFnEvals = min(1000, max(11, round(options.iter_max)));
else
    maxFnEvals = 100;
end

if isFinite
    Ncap = bufCap;
    % Closed piecewise QBD: T=[0,Ncap], K=1 regime.
    % Lv{1}=level 0, L{1}=interior repeating, Lv{2}=top level Ncap.
    LvTop = kron(Da0 + Da1, Ins) + kron(Ina, Ds0);   % arrivals blocked (lost) at full buffer
    Bc = {Brep};
    Lc = {Lrep};
    Fc = {Frep};
    Lvc = {Lv0, LvTop};
    Tc = [0, Ncap];
    % Level dimensions: level 0 and levels 1..Ncap all na*ns here.
    lvsz = na * ns;

    % EN(s), U(s), Dep(s) as scalar Laplace functions inverted per metric.
    EN = zeros(size(times)); DEP = zeros(size(times)); Uval = zeros(size(times));
    EN(posMask)  = matlab_ilt(@(s) laplace_EN_finite(s, Bc, Lc, Fc, Lvc, Tc, Ncap, pi0), tpos, maxFnEvals);
    P0 = matlab_ilt(@(s) laplace_P0(s, Bc, Lc, Fc, Lvc, Tc, pi0, wOne, false), tpos, maxFnEvals);
    DEP(posMask) = matlab_ilt(@(s) laplace_DEP_finite(s, Bc, Lc, Fc, Lvc, Tc, Ncap, pi0, wDep), tpos, maxFnEvals);
    Uval(posMask) = 1 - P0;
else
    % Open piecewise QBD: T=[0,1], K=2 (regime 1 = level 0 boundary; regime 2 repeats).
    Bo = {B0, Brep};
    Lo = {[], Lrep};
    Fo = {F0, Frep};
    Lvo = {Lv0, Lrep};
    To = [0, 1];

    % Homogeneous open case: V(s,0,m) = V(s,0,1) * R^(m-1) for m>=1, so the
    % level sums have exact closed forms (no truncation):
    %   E[N](s)  = pi0 * V(s,0,1) * (I-R)^-2 * 1
    %   Tput(s)  = pi0 * V(s,0,1) * (I-R)^-1 * wDep
    %   P0(s)    = pi0 * V(s,0,0) * 1
    EN = zeros(size(times)); DEP = zeros(size(times)); Uval = zeros(size(times));
    EN(posMask)  = matlab_ilt(@(s) laplace_EN_open(s, Bo, Lo, Fo, Lvo, To, pi0), tpos, maxFnEvals);
    P0 = matlab_ilt(@(s) laplace_P0(s, Bo, Lo, Fo, Lvo, To, pi0, wOne, true), tpos, maxFnEvals);
    DEP(posMask) = matlab_ilt(@(s) laplace_DEP_open(s, Bo, Lo, Fo, Lvo, To, pi0, wDep), tpos, maxFnEvals);
    Uval(posMask) = 1 - P0;
end

%% Package results
Qt = cell(M, K); Ut = cell(M, K); Tt = cell(M, K);
Qt{queueIdx, 1} = [EN, times];
Ut{queueIdx, 1} = [Uval, times];
Tt{queueIdx, 1} = [DEP, times];
end

% ------------------------------------------------------------------------
function val = laplace_P0(s, B, L, F, Lv, T, pi0, wOne, isOpen)
% Laplace transform of P(N(t)=0): pi0 * V(s,0,0) * 1.
if isOpen
    V = mam_transient2_open(B, L, F, Lv, T, 0, 0, s);
else
    V = mam_transient2(B, L, F, Lv, T, 0, 0, s);
end
val = pi0 * V * wOne;
end

% ------------------------------------------------------------------------
function val = laplace_EN_finite(s, B, L, F, Lv, T, Ncap, pi0)
% sum_{m=0}^{Ncap} m * pi0 * V(s,0,m) * 1.
acc = 0;
for m = 1:Ncap
    V = mam_transient2(B, L, F, Lv, T, 0, m, s);
    acc = acc + m * (pi0 * V * ones(size(V, 2), 1));
end
val = acc;
end

% ------------------------------------------------------------------------
function val = laplace_DEP_finite(s, B, L, F, Lv, T, Ncap, pi0, wDep)
% sum_{m=1}^{Ncap} pi0 * V(s,0,m) * wDep (departure rate).
acc = 0;
for m = 1:Ncap
    V = mam_transient2(B, L, F, Lv, T, 0, m, s);
    acc = acc + pi0 * V * wDep;
end
val = acc;
end

% ------------------------------------------------------------------------
function val = laplace_EN_open(s, B, L, F, Lv, T, pi0)
% E[N](s) = pi0 * V(s,0,1) * (I-R)^-2 * 1 (exact geometric first moment).
K = length(T);
IK = eye(size(L{K}, 1));
[~, R] = qbd_fundmat(B{K}, L{K} - s*IK, F{K}, 'GR');
V1 = mam_transient2_open(B, L, F, Lv, T, 0, 1, s);
ImR = IK - R;
val = pi0 * V1 * (ImR \ (ImR \ ones(size(IK, 1), 1)));
end

% ------------------------------------------------------------------------
function val = laplace_DEP_open(s, B, L, F, Lv, T, pi0, wDep)
% Tput(s) = pi0 * V(s,0,1) * (I-R)^-1 * wDep (exact geometric zeroth moment).
K = length(T);
IK = eye(size(L{K}, 1));
[~, R] = qbd_fundmat(B{K}, L{K} - s*IK, F{K}, 'GR');
V1 = mam_transient2_open(B, L, F, Lv, T, 0, 1, s);
val = pi0 * V1 * ((IK - R) \ wDep);
end
