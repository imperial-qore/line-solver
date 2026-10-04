%{
%{
 % @file pfqn_sdrprob.m
 % @brief Krzesinski state-dependent routing probabilities.
%}
%}

%{
%{
 % @brief Krzesinski state-dependent routing probabilities.
 % @fn pfqn_sdrprob(sdr, n)
 % @param sdr State-dependent routing structure.
 % @param n Per-station total population vector.
 % @return P Routing probabilities from the entry center to each branch entry.
 % @return Ped Probability of routing from the entry center to the departure center.
%}
%}
function [P, Ped] = pfqn_sdrprob(sdr, n)
% [P, PED] = PFQN_SDRPROB(SDR, N)
%
% State-dependent routing (SDR) probabilities of Krzesinski (1987),
% "Multiclass Queueing Networks with State-Dependent Routing", Performance
% Evaluation 7:125-143, eq. (10):
%
%   P_{e,e(b)}(N) = delta_tb(m_b) * prod_{s=1}^{t} omega_{s-1,s}(v_s)/omega_ss(v_s)
%
% for the branch b at level t = level(b), and zero whenever omega_ss(v_s) = 0
% for any s <= t. Here m_b is the total population of branch b, v_s the total
% population of the subnetwork Q(V_s,V_s), and
%
%   delta_tb(m) = C_t m + d_tb,    omega_tt(v) = C_t v + D_tt,
%   omega_{t-1,t}(v) = C_{t-1} v + D_{t-1,t},   omega_{0,1}(v) = 1.
%
% N is the 1xM vector of total station populations. P is 1xB with P(b) the
% probability of proceeding from the entry center e of Q(V,V) to the entry
% center e(b) of branch b; P(1) is zero because branch index 1 denotes the
% complement M-V, which is not reached from e. PED = 1 - sum(P) is the
% probability of proceeding directly to the departure center d, that is, of
% being denied entry into Q(V,V) and returned to e (Sec. 2.5).
%
% These probabilities are chain independent: they are functions of the total
% branch and subnetwork populations, not of the per-chain populations. The
% chain-dependent form of eq. (1) has no published product form (the paper
% defers it to an unpublished IBM report) and is refused elsewhere.
%
% A branch population m_b with delta_tb(m_b) < 0 lies beyond the bound that
% SDR itself enforces and so is unreachable; the probability returned there is
% zero, consistent with eq. (10) never routing a customer into such a state.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

c = pfqn_sdrcoeff(sdr);
n = n(:)';

m = zeros(1,c.B);
for b = 2:c.B
    m(b) = sum(n(c.branch{b}));
end
v = zeros(1,c.T);
for t = 1:c.T
    v(t) = sum(m(c.inA{t}));
end

% omega_ss(v_s) and omega_{s-1,s}(v_s) at the current state
om = zeros(1,c.T);
omprev = ones(1,c.T);
for s = 1:c.T
    om(s) = c.C(s)*v(s) + c.Dtt(s);
    if s > 1
        omprev(s) = c.C(s-1)*v(s) + c.Dprev(s);
    end
end

P = zeros(1,c.B);
for b = 2:c.B
    t = c.level(b);
    if any(om(1:t) <= 0)
        continue % eq. (10): the branch is closed to new arrivals
    end
    delta = c.C(t)*m(b) + c.d(t,b);
    if delta <= 0
        continue
    end
    P(b) = delta * prod(omprev(1:t) ./ om(1:t));
end
Ped = 1 - sum(P);
end
