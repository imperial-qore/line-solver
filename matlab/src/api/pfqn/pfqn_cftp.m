%{
%{
 % @file pfqn_cftp.m
 % @brief Perfect/approximate stationary state sampler for closed
 %        single-class multiserver product-form networks.
%}
%}

%{
%{
 % @brief Draws states exactly distributed according to the product-form
 %        stationary distribution of a closed single-class Jackson network
 %        with multiple servers, using monotone Coupling From The Past
 %        (Propp-Wilson) as proposed by Kijima and Matsui (WSC 2005,
 %        "Approximate/Perfect Samplers for Closed Jackson Networks").
 % @fn pfqn_cftp(L, N, S, nsamples, method)
 % @param L Service demand (visit ratio / service rate) vector, one entry
 %          per station, L(i) = theta_i/mu_i.
 % @param N Total closed population K (scalar, single class).
 % @param S Number of servers per station (default: 1; use Inf for a
 %          delay/infinite-server station).
 % @param nsamples Number of independent samples to draw (default: 1).
 % @param method 'cftp' for exact/perfect sampling (default) or 'approx'
 %               for the rapidly-mixing approximate sampler M_A.
 % @return Q Empirical mean queue length per station (1 x M).
 % @return X Sampled states, one per row (nsamples x M), each row sums to N.
 % @return T Per-sample coalescence horizon ('cftp') or mixing steps used
 %           ('approx'), returned as (nsamples x 1) for diagnostics.
%}
%}
function [Q,X,T] = pfqn_cftp(L,N,S,nsamples,method)
% [Q,X,T] = PFQN_CFTP(L,N,S,NSAMPLES,METHOD)
%
% Exact (perfect) stationary state sampling for closed single-class
% multiserver product-form networks via monotone Coupling From The Past.
%
% Reference: S. Kijima and T. Matsui, "Approximate/Perfect Samplers for
% Closed Jackson Networks", Proc. Winter Simulation Conference, 2005.
%
% Input:
% L        - demands (stations x 1), L(i) = theta_i/mu_i
% N        - total population K (scalar)
% S        - servers per station (stations x 1), Inf for infinite server
% NSAMPLES - number of independent samples (default 1)
% METHOD   - 'cftp' (exact, default) or 'approx' (rapidly mixing M_A)
%
% Output:
% Q - empirical mean queue length (1 x stations)
% X - sampled states (nsamples x stations), each row sums to N
% T - coalescence horizon / mixing steps per sample (nsamples x 1)

L = L(:).';                       % row vector of demands
M = numel(L);
if nargin < 3 || isempty(S)
    S = ones(1,M);
end
S = S(:).';
if nargin < 4 || isempty(nsamples)
    nsamples = 1;
end
if nargin < 5 || isempty(method)
    method = 'cftp';
end
K = round(N);

if M < 2
    error('pfqn_cftp: at least two stations are required.');
end
if any(L <= 0)
    error('pfqn_cftp: all demands L must be strictly positive.');
end

% Precompute cumulative log service factors:
%   logfac(i,m+1) = sum_{t=1}^m log(min(t,S(i))),  m = 0..K
% so that log alpha_i(m) = m*log(L(i)) - logfac(i,m+1).
logfac = zeros(M,K+1);
for i = 1:M
    acc = 0;
    for m = 1:K
        acc = acc + log(min(m,S(i)));
        logfac(i,m+1) = acc;
    end
end
logL = log(L);

X = zeros(nsamples,M);
T = zeros(nsamples,1);
for smp = 1:nsamples
    switch lower(method)
        case 'cftp'
            [X(smp,:),T(smp)] = draw_cftp(logL,logfac,K,M);
        case 'approx'
            [X(smp,:),T(smp)] = draw_approx(logL,logfac,K,M);
        otherwise
            error('pfqn_cftp: unknown method ''%s''.',method);
    end
end
Q = mean(X,1);
end

% ---- one exact draw via monotone CFTP ------------------------------------
function [x,horizon] = draw_cftp(logL,logfac,K,M)
% top state x_U = (K,0,...,0), bottom state x_L = (0,...,0,K)
u = [];                            % u(t): uniform for step t before present
Tback = 1;
while true
    old = numel(u);
    % prepend randomness for the newly exposed (older) steps T .. 2T-1;
    % older steps sit at higher indices and are applied first
    u(old+1:Tback) = rand(1,Tback-old); %#ok<AGROW>
    xU = [K, zeros(1,M-1)];
    xL = [zeros(1,M-1), K];
    for t = Tback:-1:1             % oldest step first
        xU = monotone_update(xU,u(t),logL,logfac,M);
        xL = monotone_update(xL,u(t),logL,logfac,M);
    end
    if isequal(xU,xL)
        x = xU;
        horizon = Tback;
        return;
    end
    Tback = 2*Tback;
end
end

% ---- monotone update on a consecutive pair -------------------------------
function x = monotone_update(x,u,logL,logfac,M)
% single uniform u in [0,1) encodes: pair index (integer part) and the
% split Lambda (fractional part), per Kijima-Matsui section 4.1
lam = 1 + u*(M-1);                 % in [1,M)
j = floor(lam);                    % pair (j, j+1), 1 <= j <= M-1
if j > M-1, j = M-1; end
Lambda = lam - j;                  % uniform in [0,1)
k = x(j) + x(j+1);
l = split_index(logL,logfac,j,j+1,k,Lambda);
x(j)   = l;
x(j+1) = k - l;
end

% ---- inverse-CDF split: smallest l with Lambda <= g^k_{ij}(l) ------------
function l = split_index(logL,logfac,i,j,k,Lambda)
% w(s) propto alpha_i(s) * alpha_j(k-s), s = 0..k
s = 0:k;
lw = (s.*logL(i) - logfac(i,s+1)) + ((k-s).*logL(j) - logfac(j,k-s+1));
lw = lw - max(lw);
w  = exp(lw);
cdf = cumsum(w);
cdf = cdf / cdf(end);
l = find(Lambda <= cdf,1,'first') - 1;   % s index (0-based) -> l
if isempty(l), l = k; end
end

% ---- approximate rapidly-mixing sampler M_A ------------------------------
function [x,steps] = draw_approx(logL,logfac,K,M,eps)
if nargin < 5, eps = 1e-2; end
steps = ceil(M*(M-1)/2 * log(K/eps));     % mixing-time bound (Theorem 1)
x = mnrnd_start(K,M);
for t = 1:steps
    p = randperm(M,2);            % distinct pair, not necessarily adjacent
    i = p(1); j = p(2);
    k = x(i) + x(j);
    l = split_index(logL,logfac,i,j,k,rand);
    x(i) = l;
    x(j) = k - l;
end
end

function x = mnrnd_start(K,M)
% arbitrary feasible start: all customers at station 1
x = [K, zeros(1,M-1)];
end
