%{
%{
 % @file pfqn_xia.m
 % @brief Xia's asymptotic approximation of the load-dependent normalizing constant.
%}
%}

%{
%{
 % @brief Xia's asymptotic approximation of the load-dependent normalizing constant.
 % @fn pfqn_xia(L, N, s)
 % @param L Service demand vector (M x 1).
 % @param N Closed population (scalar).
 % @param s Server counts (M x 1).
 % @return lGasy Logarithm of the approximate normalizing constant.
%}
%}
function lGasy = pfqn_xia(L, N, s)
% LGASY = PFQN_XIA(L, N, S)
%
% Xia's asymptotic approximation of the normalizing constant of a
% load-dependent (multiserver) single-class closed network.
%
% The demands are first rescaled so that the largest per-server utilization
% rho_i = L_i/s_i is one. The stations that attain it are the BOTTLENECK SET B;
% they saturate and contribute the M/M/s saturated term, while every other
% station contributes its finite-capacity Erlang-like partial sum
%
%   F(u,k) = sum_{j<k} u^j/j! + (u^k/k!)/(1 - u/k),
%
% the closed form of the geometric tail beyond the k-th server. The result is
%
%   log G ~ -log((|B|-1)!) - N log(c) + sum_{b in B} [ s_b log L_b - log(s_b!) ]
%                                     + sum_{k not in B} log F(L_k, s_k),
%
% c being the rescaling factor. The leading behaviour in N enters ONLY through
% -N log(c): this is the large-population limit, so the approximation does not
% resolve the O(1) corrections a finite population carries.
%
% A non-bottleneck station with u > k gives a NEGATIVE F, whose logarithm is
% complex here and NaN in the reference; NaN is returned so the refusal reads
% the same in every codebase. Only an infinite F (u == k exactly) is dropped:
% suppressing a negative term would quietly return a plausible number for a
% model the expansion does not cover. The condition cannot arise when every
% station has one server.
%
% See also PFQN_GLD, PFQN_NCLD, PFQN_PANACEALD.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

L = L(:);
s = s(:);
M = numel(L);
if M == 0
    line_error(mfilename, 'pfqn_xia requires at least one station.');
end
if numel(s) ~= M
    line_error(mfilename, 'pfqn_xia: L and s disagree on the station count.');
end
if any(L <= 0)
    line_error(mfilename, 'pfqn_xia requires positive demands.');
end
if any(s <= 0)
    line_error(mfilename, 'pfqn_xia requires positive server counts.');
end

rho = L ./ s;
scalefactor = 1 / max(rho);
Ls = L * scalefactor;
rs = rho * scalefactor;

bnkset = find(rs == max(rs));
nbnkset = setdiff((1:M)', bnkset);
B = numel(bnkset);

lGasy = -factln(B - 1) - N * log(scalefactor);
for b = bnkset(:)'
    lGasy = lGasy + s(b) * log(Ls(b)) - factln(s(b));
end
for k = nbnkset(:)'
    f = xia_F(Ls(k), s(k));
    if ~isfinite(f)
        continue
    end
    if f < 0
        lGasy = NaN;   % log of a negative F poisons the whole constant
        return
    end
    lGasy = lGasy + log(f);
end
end

function lf = factln(n)
lf = gammaln(1 + n);
end

% ---- F(u,k) = sum_{j<k} u^j/j! + (u^k/k!)/(1 - u/k) -----------------------
function ret = xia_F(u, k)
ret = 0;
j = 0;
while j < k
    ret = ret + xia_pow_over_fact(u,j);      % u^j/j! through logs
    j = j + 1;
end
ret = ret + xia_pow_over_fact(u,k) / (1 - u / k);
end

function v = xia_pow_over_fact(u,j)
% u^j/j! without forming either half: the quotient is bounded by exp(u) but both
% u^j and j! leave the double range for j >~ 171.
if j == 0
    v = 1;
elseif u == 0
    v = 0;
else
    v = exp(j*log(u) - gammaln(j+1));
end
end
