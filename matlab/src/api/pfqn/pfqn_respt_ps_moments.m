function [W, W2, out] = pfqn_respt_ps_moments(S, N, Z, method)
% [W,W2,OUT] = PFQN_RESPT_PS_MOMENTS(S, N, Z, METHOD)
%
% Sojourn-time moments at the processor-sharing station of the closed
% terminal-driven system of Mitra and Morrison (1983): a bank of terminals in
% series with a single processor-sharing CPU, with class-dependent exponential
% think times (mean Z(r)) and class-dependent exponential service times (mean
% S(r)), and N(r) jobs of class r cycling between the two.
%
% Two routes to the moments are implemented, both from that paper:
%
%   'exact'      solves the linear system c'[A - q_J I] = -pi'B of Proposition
%                3 on the state space {n : 0 <= n <= K}, K being the population
%                vector with the tagged class decremented by one. The moments
%                are then E[W_J] = sum_n c(n) and (q_J/2) E[W_J^2] =
%                sum_n (n'1+1) c(n). Exact to solver precision, at the cost of
%                a linear solve of dimension prod_r (K(r)+1).
%
%   'asymptotic' evaluates the two leading terms of the asymptotic expansion in
%                inverse powers of the large parameter Nexp = max_r Z(r)/S(r),
%                E[W_J^2] ~ c0 + c1/Nexp, of Proposition 6. The cost is a
%                linear system of dimension R, the number of classes, and is
%                therefore independent of the populations. Note that the
%                expansion parameter is the think-to-service ratio and NOT the
%                population, so a model with short think times is expanded in a
%                small parameter no matter how many jobs it holds.
%
%   'auto'       (default) takes the exact route when the state space has at
%                most 4096 states and the asymptotic route otherwise.
%
% The asymptotic route requires the normal-usage condition alpha > 0, where
% alpha = 1 - sum_r lambda_r/q_r with lambda_r = K(r)/Z(r) and q_r = 1/S(r), is
% the unutilized fraction of the CPU in the corresponding open system. Where it
% fails and the exact route is not affordable, the entry of W and W2 is NaN and
% OUT.method records 'unavailable'; asking for 'asymptotic' explicitly in that
% regime is an error rather than a blank.
%
% Input:
%   S      - per-class mean service times at the PS station (1,R), positive
%   N      - per-class populations (1,R), non-negative integers
%   Z      - per-class mean think times (1,R), positive where N > 0
%   method - 'auto' (default), 'exact' or 'asymptotic'
%
% Output:
%   W   - per-class mean sojourn times at the PS station (1,R)
%   W2  - per-class second moments of the sojourn time (1,R)
%   out - struct with fields method (1,R cell), c0, c1 (1,R, asymptotic route
%         only), alpha (1,R), nstates (1,R) and expansionParam (scalar Nexp)
%
% A class with N(r) = 0 has no sojourn time and its entries are NaN.
%
% Reference: D. Mitra, J. A. Morrison, "Asymptotic Expansions of Moments of
% the Waiting Time in Closed and Open Processor-Sharing Systems with Multiple
% Job Classes", Adv. Appl. Prob. 15(4), 1983, Propositions 3 and 6.
%
% See also: qsys_mm1_ps (the open counterpart, exact in closed form).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4 || isempty(method)
    method = 'auto';
end
method = lower(strtrim(char(method)));
if ~any(strcmp(method, {'auto', 'exact', 'asymptotic'}))
    line_error(mfilename, 'method must be one of auto, exact, asymptotic');
end

S = S(:)'; N = N(:)'; Z = Z(:)';
R = numel(S);
if numel(N) ~= R || numel(Z) ~= R
    line_error(mfilename, 'S, N and Z must have the same number of classes');
end
if any(~isfinite(S)) || any(S <= 0)
    line_error(mfilename, 'S must be finite and positive');
end
if any(~isfinite(N)) || any(N < 0) || any(N ~= round(N))
    line_error(mfilename, 'N must contain non-negative integers');
end
act = find(N > 0);
if any(~isfinite(Z(act))) || any(Z(act) <= 0)
    line_error(mfilename, 'Z must be finite and positive for every populated class');
end

AUTO_MAX = 4096;        % state-space size below which auto goes exact
EXACT_MAX = 65536;      % hard bound on an explicitly requested exact solve

W = nan(1, R); W2 = nan(1, R);
out = struct('method', {repmat({'none'}, 1, R)}, 'c0', nan(1, R), ...
    'c1', nan(1, R), 'alpha', nan(1, R), 'nstates', nan(1, R), ...
    'expansionParam', NaN);
if isempty(act)
    return
end

qa = 1 ./ S(act);
pa = 1 ./ Z(act);
out.expansionParam = max(qa ./ pa);

for J = act
    jj = find(act == J);
    K = N(act);
    K(jj) = K(jj) - 1;
    ns = prod(K + 1);
    lambda = pa .* K;
    alpha = 1 - sum(lambda ./ qa);
    out.alpha(J) = alpha;
    out.nstates(J) = ns;
    useExact = strcmp(method, 'exact') || (strcmp(method, 'auto') && ns <= AUTO_MAX);
    if useExact
        if ns > EXACT_MAX
            line_error(mfilename, sprintf(['the exact route needs a linear solve of dimension %d, ' ...
                'above the bound of %d; use method = ''asymptotic'''], ns, EXACT_MAX));
        end
        [W(J), W2(J)] = exactMoments(pa, qa, K, jj);
        out.method{J} = 'exact';
        continue
    end
    if alpha <= 0
        if strcmp(method, 'asymptotic')
            line_error(mfilename, sprintf(['the asymptotic expansion needs normal usage alpha > 0, ' ...
                'but class %d gives alpha = %.6f'], J, alpha));
        end
        out.method{J} = 'unavailable';
        continue
    end
    [W(J), W2(J), out.c0(J), out.c1(J)] = asymptoticMoments(pa, qa, K, jj);
    out.method{J} = 'asymptotic';
end
end

% =========================================================================
function [W, W2] = exactMoments(p, q, K, J)
% Proposition 3: the moments follow from c, the solution of c'[A - q_J I] =
% -pi'B, with A the generator-like operator of equation (26) and B the diagonal
% operator B(n,n) = n'1+1.
R = numel(K);
dims = K + 1;
ns = prod(dims);
stride = cumprod([1, dims(1:end-1)]);
lin = (0:ns-1)';
states = zeros(ns, R);
res = lin;
for j = 1:R
    states(:, j) = mod(res, dims(j));
    res = floor(res ./ dims(j));
end
tot = sum(states, 2);

% stationary law (15), in logs so that large populations do not overflow
r = p ./ q;
logpi = gammaln(tot + 1);
for j = 1:R
    nj = states(:, j);
    logpi = logpi + gammaln(K(j) + 1) - gammaln(nj + 1) - gammaln(K(j) - nj + 1);
    if r(j) > 0
        logpi = logpi + nj * log(r(j));
    else
        logpi(nj > 0) = -Inf;
    end
end
logpi = logpi - max(logpi);
pin = exp(logpi);
pin = pin / sum(pin);

rows = []; cols = []; vals = [];
diagv = zeros(ns, 1);
for j = 1:R
    nj = states(:, j);
    dn = nj >= 1;
    if any(dn)
        rows = [rows; lin(dn) - stride(j) + 1]; %#ok<AGROW>
        cols = [cols; lin(dn) + 1];             %#ok<AGROW>
        vals = [vals; p(j) * (K(j) - nj(dn) + 1) .* tot(dn)]; %#ok<AGROW>
    end
    up = nj <= K(j) - 1;
    if any(up)
        rows = [rows; lin(up) + stride(j) + 1]; %#ok<AGROW>
        cols = [cols; lin(up) + 1];             %#ok<AGROW>
        vals = [vals; (nj(up) + 1) * q(j)];     %#ok<AGROW>
    end
    diagv = diagv - (p(j) * (K(j) - nj) .* (tot + 1) + nj * q(j));
end
rows = [rows; lin + 1]; cols = [cols; lin + 1]; vals = [vals; diagv];
A = sparse(rows, cols, vals, ns, ns);

c = (A - q(J) * speye(ns))' \ (-(tot + 1) .* pin);
W = sum(c);
W2 = (2 / q(J)) * sum((tot + 1) .* c);
end

% =========================================================================
function [W, W2, c0, c1] = asymptoticMoments(p, q, K, J)
% Proposition 6: the two leading terms of the expansion in 1/Nexp. Equation
% numbers below are those of Mitra and Morrison (1983).
lambda = p .* K;
alpha = 1 - sum(lambda ./ q);
Nexp = max(q ./ p);                              % (50)
Gam = Nexp * p ./ q;                             % (51)
beta = K / Nexp;                                 % (51)
qJ = q(J);

den = 1 - sum(lambda ./ (q + qJ));
F10 = (-1 / (alpha^2 * qJ)) * (1 - sum(lambda .* (q - qJ) ./ (q .* (q + qJ)))) / den;   % (110)
c0 = -2 / qJ * F10;

bg2 = sum(beta .* Gam.^2);
f1 = lambda ./ (q + qJ) .* (F10 - 2 ./ (alpha^2 * q));                                  % (113iii)
S2j = 6 / alpha^4 * (alpha * beta .* Gam.^2 + 2 * bg2 * beta .* Gam);                   % (113i)
S2js = 3 / alpha^3 * ((beta .* Gam)' * (beta .* Gam));                                  % (113ii)

R = numel(K);
Amat = eye(R);
rhs = zeros(R, 1);
for j = 1:R
    for s = 1:R
        d = q(j) + q(s) + qJ;
        Amat(j, s) = Amat(j, s) - lambda(j) / d;
        Amat(j, j) = Amat(j, j) - lambda(s) / d;
        rhs(j) = rhs(j) + S2js(j, s) / d;
    end
    rhs(j) = rhs(j) - f1(j);
end
F2 = (Amat \ rhs)';                                                                     % (112)

f10 = -3 / (alpha^3 * qJ) * bg2;                                                        % (98)
F20 = (sum((2 * Gam .* q .* F2 + S2j) ./ (q + qJ)) - f10) / den;                         % (111)
c1 = -2 / qJ * F20 + c0 / alpha^2 * bg2;                                                % (114ii)

W = 1 / (alpha * qJ) * (1 - 2 / Nexp * bg2 / alpha^2);                                  % (68)
W2 = c0 + c1 / Nexp;                                                                    % (114i)
end
