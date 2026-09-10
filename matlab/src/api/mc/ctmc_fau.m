function [pit,info] = ctmc_fau(pi0,Q,t,epsilon,delta,maxsteps)
% [PIT,INFO]=CTMC_FAU(PI0,Q,T,EPSILON,DELTA,MAXSTEPS)
%
% Transient distribution of a CTMC at time T by fast adaptive uniformization.
%
% ADAPTIVE UNIFORMIZATION (van Moorsel and Sanders, 1994). Ordinary
% uniformization fixes one rate q >= max_i |q_ii| over the WHOLE state space
% and mixes the powers of P = I + Q/q against a Poisson(q*t) law. Its cost is
% therefore set by the fastest state anywhere, including states that carry no
% probability at time t. Adaptive uniformization instead picks a rate per step
% from the states the iterate actually occupies,
%
%   Lambda_n >= max{ |q_ii| : i in supp(u^(n)) },  u^(n+1) = u^(n)(I + Q/Lambda_n),
%
% which keeps every entry of u^(n+1) nonnegative. The subordinating process is
% then no longer Poisson but the pure birth process N(t) with rates
% Lambda_0, Lambda_1, ..., and
%
%   pi(t) = sum_{n>=0} P{N(t)=n} u^(n).
%
% FAST ADAPTIVE UNIFORMIZATION (Mateescu, Wolf, Didier and Henzinger, 2010)
% adds the second half: an entry of u^(n) below DELTA is dropped rather than
% propagated, so the support tracks the states of non-negligible occupancy
% instead of the reachable set. Dropping only removes nonnegative
% contributions, so PIT is a componentwise LOWER BOUND on the exact
% distribution. Nothing is renormalized anywhere, so the error is not
% estimated but MEASURED: the three approximations (the birth index truncated
% at K, the Poisson window of the weight computation, the DELTA threshold)
% each remove mass and none puts any back, whence
%
%   0 <= pi(t) - PIT  componentwise, and
%   |pi(t) - PIT|_1 = sum(PI0) - sum(PIT) = INFO.errorBound.
%
% This is what replaces a blind population cutoff on an open model by a
% numerical one whose error is a returned quantity. It is a transient method:
% it produces no stationary distribution.
%
% THE BIRTH WEIGHTS ARE COMPUTED EXACTLY, not quadratured. The rates
% Lambda_0..Lambda_K generate a bidiagonal generator B on the birth index
% 0..K plus one absorbing overflow index, and P{N(t)=n} is the transient
% distribution of that scalar chain, obtained by uniformizing it at
% Lstar = max_n Lambda_n and applying the shipped Fox-Glynn weights. Every
% entry of the uniformized bidiagonal kernel lies in [0,1], so there is no
% cancellation, and the mass reaching the overflow index IS the truncation
% error INFO.weightTail. The alternative found in most implementations,
% integrating the convolution
% b_n(s) = Lambda_{n-1} int_0^s b_{n-1}(v) exp(-Lambda_n (s-v)) dv
% on a time grid, carries a quadrature error that the reported bound would
% then have to absorb.
%
% THE SWEEP RUNS TWICE, and this is the one real cost of the exact weights.
% b_n(t) needs the rates up to n, which are not known before the sweep ends,
% while u^(n) is needed after them; storing every iterate would cost
% K * |support| doubles. The sweep is therefore replayed, and since it is
% deterministic the second pass reproduces the first rate for rate and drop
% for drop. Two sparse sweeps still beat ordinary uniformization whenever
% Lstar is well below max_i |q_ii| or the occupied support is well below n.
%
% STOPPING IS CERTIFIED, not heuristic. S_{K+1} = sum_{m<=K} Exp(Lambda_m)
% dominates an Erlang(K+1, Lstar) stochastically, so
% P{N(t) > K} = P{S_{K+1} <= t} <= P{Poisson(Lstar*t) >= K+1}, bounded above
% by the Chernoff exponent of that Poisson tail. The sweep stops at the first
% K meeting EPSILON, which also shows that this method never takes more steps
% than uniformization at the largest rate it visited.
%
% -- Input
% PI0      : initial distribution, 1xn
% Q        : infinitesimal generator, nxn (sparse is used as given)
% T        : time horizon, T >= 0
% EPSILON  : birth-process truncation tolerance (default 1e-6)
% DELTA    : occupancy threshold below which a state is dropped (default 1e-12)
% MAXSTEPS : cap on birth steps; nonpositive for the default cap of 1e6
%
% -- Output
% PIT  : 1xn defective distribution at time T, a lower bound on pi(T)
% INFO : struct with fields
%        steps        number of birth steps K+1 actually taken
%        lambdaMin    smallest adaptive rate used
%        lambdaMax    largest adaptive rate used, the Lstar above
%        uniformRate  max_i |q_ii|, the rate ordinary uniformization would use
%        weightTail   mass reaching the overflow index, i.e. P{N(T) > K}
%        weightWindow Poisson mass outside the Fox-Glynn window of the weights
%        droppedMass  probability removed by the DELTA threshold
%        errorBound   sum(PI0)-sum(PIT), which IS |pi(T)-PIT|_1
%        supportMax   largest occupied support over the sweep
%        supportFinal support at the last step
%        truncated    true if MAXSTEPS stopped the sweep
%        absorbed     true if the support emptied or became absorbing

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if nargin<4 || isempty(epsilon)
    epsilon = 1e-6;
end
if nargin<5 || isempty(delta)
    delta = 1e-12;
end
if nargin<6 || isempty(maxsteps)
    maxsteps = -1;
end
if epsilon<=0
    epsilon = 1e-6;
end
if delta<0
    delta = 0;
end
if maxsteps<=0
    maxsteps = 1e6;
end
n = size(Q,1);
if size(Q,2)~=n
    line_error(mfilename,'Q must be square.');
end
pi0 = full(pi0(:)).';
if numel(pi0)~=n
    line_error(mfilename,'PI0 and Q have inconsistent sizes.');
end
if t<0
    line_error(mfilename,'T must be nonnegative.');
end

d = full(-diag(Q)).';                 % exit rates, nonnegative on a generator
info = struct('steps',0,'lambdaMin',0,'lambdaMax',0,'uniformRate',max([d,0]), ...
    'weightTail',0,'weightWindow',0,'droppedMass',0,'errorBound',0, ...
    'supportMax',0,'supportFinal',0,'truncated',false,'absorbed',false);
if t==0 || n==0
    pit = pi0;
    info.steps = 1;
    info.supportMax = nnz(pi0);
    info.supportFinal = info.supportMax;
    return
end
Qt = Q.';                             % columns of Qt are rows of Q

% Pass one: the adaptive rate sequence, and where it stops.
[lambda,truncated,absorbed] = fau_rates(pi0,Qt,d,t,delta,maxsteps,epsilon);

% The birth-process weights of that rate sequence, exactly.
[b,wtail,wwin] = fau_weights(lambda,t,epsilon);

% Pass two: the same sweep again, accumulating sum_n b_n u^(n).
[pit,dropped,suppmax,suppfin] = fau_accumulate(pi0,Qt,d,delta,lambda,b);

info.steps = numel(lambda);
info.lambdaMin = min(lambda);
info.lambdaMax = max(lambda);
info.weightTail = wtail;
info.weightWindow = wwin;
info.droppedMass = dropped;
% Every approximation removes nonnegative mass and none is put back, so the
% error is not estimated but measured: what is missing from PIT is exactly
% what separates it from pi(T).
info.errorBound = sum(pi0) - sum(pit);
info.supportMax = suppmax;
info.supportFinal = suppfin;
info.truncated = truncated;
info.absorbed = absorbed;
end

function [lambda,truncated,absorbed] = fau_rates(pi0,Qt,d,t,delta,maxsteps,epsilon)
% Sweep the iterate to collect the adaptive rates Lambda_0..Lambda_K, stopping
% when the Poisson-dominance bound on P{N(t)>K} falls to EPSILON.
lambda = zeros(1,0);
truncated = false;
absorbed = false;
u = pi0;
act = find(u>0);
lstar = 0;
while true
    if isempty(act)
        absorbed = true;
        break
    end
    L = max(d(act));
    lambda(end+1) = L; %#ok<AGROW>
    if L<=0
        % Every occupied state is absorbing: the birth process stops here and
        % the remaining weight falls entirely on this iterate.
        absorbed = true;
        break
    end
    lstar = max(lstar,L);
    if fau_tailbound(lstar,t,numel(lambda))<=epsilon
        break
    end
    if numel(lambda)>=maxsteps
        truncated = true;
        break
    end
    [u,act,~] = fau_step(u,act,Qt,L,delta);
end
end

function [pit,dropped,suppmax,suppfin] = fau_accumulate(pi0,Qt,d,delta,lambda,b)
% Replay the sweep of FAU_RATES, accumulating sum_n b_n u^(n). The arithmetic
% is identical, so the rates and the drops reproduce those of the first pass.
nsteps = numel(lambda);
pit = zeros(1,numel(pi0));
dropped = 0;
u = pi0;
act = find(u>0);
suppmax = numel(act);
suppfin = numel(act);
for m=1:nsteps
    if isempty(act)
        break
    end
    suppmax = max(suppmax,numel(act));
    suppfin = numel(act);
    pit(act) = pit(act) + b(m)*u(act);
    if m<nsteps
        L = max(d(act));
        if L<=0
            break
        end
        [u,act,dropStep] = fau_step(u,act,Qt,L,delta);
        dropped = dropped + dropStep;
    end
end
end

function [u,act,dropStep] = fau_step(u,act,Qt,L,delta)
% One adaptive uniformization step u <- u(I + Q/L), touching only the rows of
% Q in the current support, followed by the FAU drop rule. States with a zero
% exit rate are absorbing: their row of Q is empty, so they hold their mass
% and stay in the support.
cs = Qt(:,act)*u(act).';
idx = find(cs);
dropStep = 0;
if ~isempty(idx)
    idx = reshape(idx,1,[]);
    inflow = full(cs(idx));
    vals = u(idx) + reshape(inflow,1,[])/L;
    small = vals<delta;
    if any(small)
        dropStep = sum(max(vals(small),0));
        vals(small) = 0;
    end
    u(idx) = vals;
    % A row of Q that cancels exactly leaves its state untouched and out of
    % IDX, so the surviving part of the old support is carried over too.
    act = unique([act(u(act)>0), idx(vals>0)]);
end
end

function [b,wtail,wwin] = fau_weights(lambda,t,tol)
% Transient distribution of the pure birth process with rates LAMBDA at time
% T, i.e. b_n = P{N(t)=n} for n=0..K, plus the mass WTAIL that reached the
% absorbing overflow index K+1 and therefore bounds P{N(t)>K}. The chain is
% uniformized at Lstar=max(LAMBDA) and mixed against Fox-Glynn Poisson
% weights, so the kernel entries 1-Lambda_n/Lstar and Lambda_n/Lstar are
% probabilities and nothing cancels. The weights are taken UNNORMALIZED, so
% the Poisson mass outside the Fox-Glynn window WWIN is missing from B rather
% than redistributed over it: B is then a sub-distribution, every term of the
% mixture is an underestimate, and the error stays measurable as missing mass.
k1 = numel(lambda);
lstar = max(lambda);
if lstar<=0 || t<=0
    b = zeros(1,k1);
    b(1) = 1;
    wtail = 0;
    wwin = 0;
    return
end
[left,right,w] = foxglynn_weights(lstar*t,tol,-1,false);
wwin = max(1-sum(w),0);
a = 1 - lambda/lstar;
c = lambda/lstar;
v = zeros(1,k1+1);
v(1) = 1;
b = zeros(1,k1+1);
for k=0:right
    if k>=left
        b = b + w(k-left+1)*v;
    end
    if k<right
        v = [v(1:k1).*a, v(k1+1)] + [0, v(1:k1).*c];
    end
end
wtail = b(k1+1);
b = b(1:k1);
end

function bound = fau_tailbound(lstar,t,k)
% Upper bound on P{N(t)>=k} for the birth process, through the stochastic
% domination of its k-th jump epoch by an Erlang(k,LSTAR): the bound is the
% Poisson(LSTAR*t) upper tail P{X>=k}, taken at its Chernoff exponent
% lambda*h(k/lambda) with h(u)=u*log(u)-u+1. The exponent bounds the upper
% tail only above the mean, so below it the bound is left vacuous.
lambda = lstar*t;
if k<=lambda
    bound = 1;
    return
end
bound = exp(-(lambda - k + k*log(k/lambda)));
end
