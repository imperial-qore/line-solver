%{
%{
 % @file pfqn_sdr.m
 % @brief Exact product form of a multiclass network with state-dependent routing.
%}
%}

%{
%{
 % @brief Exact product form of a multiclass network with state-dependent routing.
 % @fn pfqn_sdr(S, xi, N, sdr, alpha)
 % @param S Mean service times per station and chain.
 % @param xi Relative visit counts per station and chain.
 % @param N Chain population vector.
 % @param sdr State-dependent routing structure.
 % @param alpha Optional load-dependent rate scalings.
 % @return Q Mean queue lengths per station and chain.
 % @return X Mean throughputs per station and chain.
 % @return U Mean utilizations per station and chain.
 % @return R Mean response times per station and chain.
 % @return G Normalizing constant.
 % @return lG Logarithm of the normalizing constant.
 % @return prob Stationary probability of every enumerated state.
 % @return states Enumerated state list.
%}
%}
function [Q, X, U, R, G, lG, prob, states] = pfqn_sdr(S, xi, N, sdr, alpha)
% [Q, X, U, R, G, LG, PROB, STATES] = PFQN_SDR(S, XI, N, SDR, ALPHA)
%
% Exact evaluation of the product-form joint probability distribution of a
% multiclass closed queueing network with state-dependent routing, Krzesinski
% (1987), "Multiclass Queueing Networks with State-Dependent Routing",
% Performance Evaluation 7:125-143, eq. (16):
%
%   P(n) = G^-1 prod_i f_i(n_i)
%                prod_t [Omega_{t-1,t}(v_t)/Omega_tt(v_t)]
%                prod_{b in A_t - A_{t+1}} Delta_tb(m_b)
%
% with f_i(n_i) = [n_i!/beta_i(n_i)] prod_j gamma_ij^{n_ij}/n_ij!,
% gamma_ij = xi_ij/mu_ij, beta_i(n) = alpha_i(n) beta_i(n-1), beta_i(0) = 1,
% and the cumulative coefficients
%
%   Omega_{t-1,t}(v) = omega_{t-1,t}(v-1) Omega_{t-1,t}(v-1),  Omega(0) = 1,
%   Omega_tt(v)      = omega_tt(v-1) Omega_tt(v-1),            Omega(0) = 1,
%   Delta_tb(m)      = delta_tb(m-1) Delta_tb(m-1),            Delta(0) = 1.
%
% This form is general in the branch topology: a branch may hold several
% interconnected centers. Only the MVA and convolution algorithm of the paper's
% Section 4, implemented in PFQN_SDRMVA, is restricted to single-center
% branches. The normalizing constant here is obtained by summing the
% unnormalized weights over the whole reachable state space, which is exact for
% any branch topology but grows combinatorially with the populations.
%
% Inputs:
%   S      MxJ mean service times 1/mu_ij of a chain j customer at center i
%   XI     MxJ coefficients xi_ij of Section 3.2. For the entry and departure
%          centers of Q(V,V) and of every branch these are all equal and are
%          NOT relative visit counts; use PFQN_SDRVISITS to obtain them
%   N      1xJ chain populations
%   SDR    state-dependent routing structure, see PFQN_SDRCOEFF
%   ALPHA  optional Mxmax(N) matrix of load-dependent rate scalings, with
%          ALPHA(i,k) = alpha_i(k) the rate multiplier when k customers are
%          present at center i. Defaults to 1, a fixed-rate center. Use
%          ALPHA(i,k) = k for an infinite-server center and min(k,c) for a
%          c-server center
%
% S and XI are required separately rather than as their product: under SDR the
% xi_ij are not visit ratios, so the per-center throughputs cannot be recovered
% from the demands alone.
%
% Outputs Q, X, U and R are MxJ. X holds the per-center chain throughputs
%   T_ij = sum_n P(n) alpha_i(n_i) (n_ij/n_i) / S_ij,
% U_ij = T_ij S_ij is the mean number of chain j customers in service, and
% R_ij = Q_ij / T_ij is the mean response time at the center.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = size(S,1);
J = size(S,2);
N = N(:)';
if numel(N) ~= J
    line_error(mfilename,'The population vector must have one entry per chain.');
end
if any(size(xi) ~= [M,J])
    line_error(mfilename,'S and xi must have the same size.');
end
if nargin < 5 || isempty(alpha)
    alpha = ones(M, max(1,sum(N)));
end
if size(alpha,2) < sum(N)
    alpha = [alpha, ones(M, sum(N)-size(alpha,2))];
end

c = pfqn_sdrcoeff(sdr);
if max([c.entry, c.departure, cell2mat(c.branch(2:end))]) > M
    line_error(mfilename,'The SDR structure references a station index beyond the number of centers.');
end

Ntot = sum(N);
gamma = xi .* S;

% log beta_i(n), n = 0..Ntot
logbeta = zeros(M, Ntot+1);
for i = 1:M
    for k = 1:Ntot
        logbeta(i,k+1) = logbeta(i,k) + log(alpha(i,k));
    end
end

% log Delta_tb(m) and log Omega(v), with a flag marking the population beyond
% which the cumulative product hits a nonpositive factor and the weight is zero
[logDelta, okDelta] = sub_cumlog(c.C(c.level(2:end)), c.d(sub2ind(size(c.d), c.level(2:end), 2:c.B)), Ntot);
logOmTT = zeros(c.T, Ntot+1); okOmTT = true(c.T, Ntot+1);
logOmPrev = zeros(c.T, Ntot+1); okOmPrev = true(c.T, Ntot+1);
for t = 1:c.T
    [logOmTT(t,:), okOmTT(t,:)] = sub_cumlog(c.C(t), c.Dtt(t), Ntot);
    if t > 1
        [logOmPrev(t,:), okOmPrev(t,:)] = sub_cumlog(c.C(t-1), c.Dprev(t), Ntot);
    end
end

states = sub_states(M, N);
ns = size(states,1);
logw = -Inf(ns,1);
for s = 1:ns
    nmat = reshape(states(s,:), M, J);
    ni = sum(nmat,2)';
    lw = 0;
    ok = true;
    for i = 1:M
        lw = lw + gammaln(ni(i)+1) - logbeta(i,ni(i)+1);
        for j = 1:J
            if nmat(i,j) > 0
                if gamma(i,j) <= 0
                    ok = false; break
                end
                lw = lw + nmat(i,j)*log(gamma(i,j)) - gammaln(nmat(i,j)+1);
            end
        end
        if ~ok, break; end
    end
    if ~ok, continue; end
    m = zeros(1,c.B);
    for b = 2:c.B
        m(b) = sum(ni(c.branch{b}));
    end
    for b = 2:c.B
        if ~okDelta(b-1, m(b)+1), ok = false; break; end
        lw = lw + logDelta(b-1, m(b)+1);
    end
    if ~ok, continue; end
    for t = 1:c.T
        v = sum(m(c.inA{t}));
        if ~okOmTT(t, v+1), ok = false; break; end
        lw = lw - logOmTT(t, v+1);
        if t > 1
            if ~okOmPrev(t, v+1), ok = false; break; end
            lw = lw + logOmPrev(t, v+1);
        end
    end
    if ~ok, continue; end
    logw(s) = lw;
end

lmax = max(logw);
if isinf(lmax)
    line_error(mfilename,'The SDR network has no reachable state at the given populations: the routing coefficients forbid every state.');
end
w = exp(logw - lmax);
Gs = sum(w);
prob = w / Gs;
lG = lmax + log(Gs);
G = exp(lG);

Q = zeros(M,J);
X = zeros(M,J);
for s = 1:ns
    if prob(s) == 0, continue; end
    nmat = reshape(states(s,:), M, J);
    ni = sum(nmat,2);
    Q = Q + prob(s)*nmat;
    for i = 1:M
        if ni(i) > 0
            X(i,:) = X(i,:) + prob(s)*alpha(i,ni(i))*(nmat(i,:)./ni(i))./S(i,:);
        end
    end
end
U = X .* S;
R = zeros(M,J);
nz = X > 0;
R(nz) = Q(nz) ./ X(nz);
end

function [lc, ok] = sub_cumlog(Cv, dv, nmax)
% [LC, OK] = SUB_CUMLOG(CV, DV, NMAX)
% Cumulative log-product prod_{k=0}^{n-1} (CV k + DV) for n = 0..NMAX, one row
% per entry of CV. OK(.,n+1) is false once a nonpositive factor has been met,
% which is where the SDR population bound closes the branch or subnetwork.
p = numel(Cv);
lc = zeros(p, nmax+1);
ok = true(p, nmax+1);
for q = 1:p
    for n = 1:nmax
        f = Cv(q)*(n-1) + dv(q);
        if ~ok(q,n) || f <= 0
            ok(q,n+1) = false;
            lc(q,n+1) = -Inf;
        else
            lc(q,n+1) = lc(q,n) + log(f);
        end
    end
end
end

function states = sub_states(M, N)
% STATES = SUB_STATES(M, N)
% Every MxJ population matrix with column sums N, flattened row by row into
% the columns of STATES in column-major (station-major within chain) order.
J = numel(N);
percl = cell(1,J);
for j = 1:J
    percl{j} = sub_compositions(N(j), M);
end
states = percl{1};
for j = 2:J
    a = states; b = percl{j};
    na = size(a,1); nb = size(b,1);
    states = [repmat(a, nb, 1), reshape(repmat(b(:)', na, 1), na*nb, M)];
end
end

function C = sub_compositions(n, m)
% C = SUB_COMPOSITIONS(N, M)
% All nonnegative integer M-vectors summing to N, one per row.
if m == 1
    C = n;
    return
end
C = [];
for k = 0:n
    tail = sub_compositions(n-k, m-1);
    C = [C; k*ones(size(tail,1),1), tail]; %#ok<AGROW>
end
end
