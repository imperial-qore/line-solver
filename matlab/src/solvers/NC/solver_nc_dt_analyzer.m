function [Q,U,R,T,C,X,lG,runtime,iter,method] = solver_nc_dt_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER,METHOD] = SOLVER_NC_DT_ANALYZER(SN, OPTIONS)
%
% Exact normalizing-constant analysis of a discrete-time (slotted) queueing
% model, selected by OPTIONS.CONFIG.SLOTTED and classified by NC_IS_DT_MODEL.
% Two families are covered, both from Daduna (2001):
%
%   chapter 2  a Bernoulli server fed by a Bernoulli arrival stream, with an
%              unbounded buffer (theorem 2.3, corollary 2.7), a finite buffer
%              (corollary 2.8) or a load-dependent service probability
%              (example 2.10); evaluated by DQSYS_BERNOULLI1;
%   chapter 3  a closed cycle of Bernoulli servers (theorem 3.2, corollary
%              3.4); evaluated by DPFQN_NC when the service probabilities are
%              state independent and by DPFQN_NCLD otherwise.
%
% Every metric is expressed on the slot lattice: a rate is a per-slot
% probability and a time is a number of slots. OPTIONS.CONFIG.SLOTLENGTH
% rescales both to model time units, dividing throughputs and multiplying
% times, so that a caller who declared a slot of length d reads back the same
% units it built the model in.
%
% On the cycle route the per-class split is proportional to the per-class
% population. Service in the cycle is type independent and FCFS forbids
% overtaking, so the cyclic order of the jobs is frozen; the marginal law of
% the queue lengths therefore carries no class information, and the long-run
% share of station j held by chain g is its population share N_g/N. That is
% the sense in which section 3.2 of the reference calls the multichain case a
% direct adaptation of the unichain one.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
Tstart = tic;
iter = 1;

dt = nc_is_dt_model(sn, options);
switch dt.kind
    case 'bernoulli1'
        method = 'dt.bernoulli1';
        [Q,U,R,T,C,X,lG] = dt_single(sn, dt);
    case 'cycle'
        if size(dt.serviceProb, 2) >= 1 && all(all(abs(dt.serviceProb - dt.serviceProb(:,1)) <= GlobalConstants.FineTol))
            method = 'dt.cycle';
        else
            method = 'dt.cycleld';
        end
        [Q,U,R,T,C,X,lG] = dt_cycle(sn, dt, method);
    otherwise
        line_error(mfilename, 'options.config.slotted was requested but the model is not a discrete-time product-form model: %s', dt.reason);
end

% Slot length: metrics above are per slot, the model is in time units.
d = 1;
if isfield(options, 'config') && isfield(options.config, 'slotlength') && ~isempty(options.config.slotlength)
    d = options.config.slotlength;
end
if ~isnumeric(d) || ~isscalar(d) || ~isreal(d) || ~isfinite(d) || d <= 0
    line_error(mfilename, 'options.config.slotlength must be a positive finite scalar');
end
if d ~= 1
    T = T / d;  X = X / d;  R = R * d;  C = C * d;
end

runtime = toc(Tstart);
end

% ==========================================================================
function [Q,U,R,T,C,X,lG] = dt_single(sn, dt)
% Chapter 2 route: one Bernoulli server, one open class.
M = sn.nstations;
K = sn.nclasses;
Q = zeros(M,K); U = zeros(M,K); R = zeros(M,K); T = zeros(M,K);
C = zeros(1,K); X = zeros(1,K);

res = dqsys_bernoulli1(dt.arrivalProb, dt.serviceProb, dt.capacity);
ist = dt.station;

Q(ist,1) = res.meanQueueLength;
U(ist,1) = res.utilization;
T(ist,1) = res.throughput;
R(ist,1) = res.meanSojournTime;
X(1) = res.throughput;
C(1) = res.meanSojournTime;

% The Source row carries the offered stream, as on the continuous-time route.
src = dt.source;
T(src,1) = dt.arrivalProb;
R(src,1) = 0;
lG = log(res.normConst);
end

% ==========================================================================
function [Q,U,R,T,C,X,lG] = dt_cycle(sn, dt, method)
% Chapter 3 route: closed cycle of Bernoulli servers.
M = sn.nstations;
K = sn.nclasses;
N = dt.population;
order = dt.order;
P = dt.serviceProb;

Qs = zeros(1,M); Us = zeros(1,M); Ts = zeros(1,M);
if strcmp(method, 'dt.cycle')
    % State independent: propositions 3.18 and 3.19 end to end, with the
    % index of corollary 3.20(a) corrected (see DPFQN_NC).
    p = P(:,1)';
    q = 1 - p;
    [lG, G, G1] = dpfqn_nc(p, N);
    xtput = G1(N+1) / G;
    for k = 1:M
        tail = 0;
        for n = 1:N
            tail = tail + (q(k)/p(k))^n / q(k) * G1(N-n+1+1) / G;
        end
        Qs(order(k)) = tail;                      % E[X_j] = sum_{n>=1} P(X_j>=n)
        Us(order(k)) = xtput / p(k);              % P(X_j >= 1)
        Ts(order(k)) = xtput;
    end
else
    % State dependent: theorem 3.2 through the convolution of DPFQN_NCLD.
    [lG, G, W, Gc] = dpfqn_ncld(P, N);
    for k = 1:M
        marg = W(k,:) .* Gc(k, N+1:-1:1) / G(N+1);
        Qs(order(k)) = sum(marg .* (0:N));
        Us(order(k)) = 1 - marg(1);
        Ts(order(k)) = sum(marg(2:end) .* P(k,1:N));
    end
end

% Per-class split by population share; see the header for why this is exact.
Nr = sn.njobs(:)';
Nr(~isfinite(Nr)) = 0;
share = zeros(1,K);
if N > 0
    share = Nr / N;
end
Q = zeros(M,K); U = zeros(M,K); R = zeros(M,K); T = zeros(M,K);
C = zeros(1,K); X = zeros(1,K);
for ist = 1:M
    for r = 1:K
        Q(ist,r) = Qs(ist) * share(r);
        U(ist,r) = Us(ist) * share(r);
        T(ist,r) = Ts(ist) * share(r);
        if T(ist,r) > 0
            R(ist,r) = Q(ist,r) / T(ist,r);
        end
    end
end
for r = 1:K
    X(r) = 0;
    ref = sn.refstat(r);
    if ref >= 1 && ref <= M
        X(r) = T(ref,r);
    end
    if X(r) > 0
        C(r) = Nr(r) / X(r);
    end
end
end
