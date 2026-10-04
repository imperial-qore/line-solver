function [g, Xi, Phi] = pfqn_oi_insvc(oirate, N, options)
% [G, XI, PHI] = PFQN_OI_INSVC(OIRATE, N, OPTIONS)
%
% Conditional mean number of IN-SERVICE jobs per class at an order-independent
% (OI) station, as a function of the per-class count vector n. This is the
% quantity underlying the LINE utilization convention at OI stations,
%   U_r = E[sir_r] / c,
% with c the number of servers and sir_r the number of class-r jobs receiving a
% strictly positive service rate.
%
% In an OI station the state is the ordered list c = (c_1,...,c_n) of job
% classes (position 1 = head) and the job in position p is served at the rank
% rate increment
%   Delta_p(c) = mu(c_1..c_p) - mu(c_1..c_{p-1}),
% so that the total rate telescopes to mu(c). Position p is IN SERVICE when
% Delta_p(c) > 0, and
%   sir_r(c) = #{p : c_p = r, Delta_p(c) > 0}.
% Note that sir_r counts JOBS, not servers: a single job served concurrently by
% several compatible servers counts once. This matches the definition used by
% the exact CTMC solver (State.toMarginal, PAS branch) and by LDES.
%
% Because mu is permutation-invariant (a function of the count vector), the
% unnormalized weight of an ordering c of the multiset n factorizes over its
% prefixes as w(c) = prod_{p=1}^{|n|} 1/mu(n(c_1..c_p)), and the OI balance
% function is Phi(n) = sum_{orderings c of n} w(c), which obeys the standard
% balanced-fairness recursion (condition on the tail element c_{|n|}):
%
%   Phi(0) = 1,   Phi(n) = (1/mu(n)) sum_{r: n_r>0} Phi(n - e_r).
%
% Conditioning the same way on the tail element and using
% sir_r(c) = sir_r(c_1..c_{|n|-1}) + [c_{|n|} = r] * 1{mu(n) > mu(n - e_r)}
% gives the companion recursion for the sir-weighted balance
% Xi_r(n) = sum_{orderings c of n} w(c) sir_r(c):
%
%   Xi_r(0) = 0,
%   Xi_r(n) = (1/mu(n)) [ sum_{s: n_s>0} Xi_r(n - e_s)
%                         + 1{n_r > 0} 1{mu(n) > mu(n - e_r)} Phi(n - e_r) ].
%
% Given n, every ordering carries the same class-weight factor prod_r x_r^{n_r},
% so the conditional law of the ordering is w(c)/Phi(n) and
%
%   E[sir_r | n] = Xi_r(n) / Phi(n) =: g_r(n),
%
% a function of the count vector alone. The station's in-service mean then
% follows from the count marginal pM by E[sir_r] = sum_n pM(n) g_r(n), or, in
% normalizing-constant form, from the functional-server identity of
% PFQN_OI_FNC applied to f(n) = g_r(n) (note g_r(0) = 0, as required).
%
% Parameters:
%   oirate  - function handle mu(n) returning the OI total service rate for the
%             per-class count vector n (1 x R). mu(0) is taken as 0.
%   N       - (1 x R) closed population vector, finite.
%   options - solver options (optional, currently unused).
%
% Returns:
%   g   - (prod(N+1) x R) table, column-major over the lattice 0 <= n <= N, with
%         g(1 + sum(n .* stride), r) = E[sir_r | n].
%   Xi  - (prod(N+1) x R) table with the sir-weighted balance Xi_r(n).
%   Phi - (prod(N+1) x 1) table with the OI balance function Phi(n).
%
% See also PFQN_OI_FNC, PFQN_NCOI, PFQN_MVAOI, PFQN_MVAOI_MARG.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3
    options = struct(); %#ok<NASGU>
end
if ~isa(oirate, 'function_handle')
    line_error(mfilename, 'oirate must be a function handle mu(n).');
end
N = round(N(:)');
if any(~isfinite(N)) || any(N < 0)
    line_error(mfilename, 'pfqn_oi_insvc requires finite nonnegative populations.');
end
R = numel(N);
shp = N + 1;
stride = ones(1, R);
for d = 2:R
    stride(d) = stride(d-1) * shp(d-1);
end
total = prod(shp);

% Decode the lattice once and tabulate the rank rate mu(n).
subs = zeros(total, R);
for i = 1:total
    li = i - 1;
    for d = 1:R
        subs(i, d) = mod(li, shp(d));
        li = floor(li / shp(d));
    end
end
muv = zeros(total, 1);
for i = 1:total
    if sum(subs(i, :)) > 0
        muv(i) = oirate(subs(i, :));
    end
end

Phi = zeros(total, 1);
Xi = zeros(total, R);
for i = 1:total
    n = subs(i, :);
    if sum(n) == 0
        Phi(i) = 1;
        continue
    end
    mun = muv(i);
    if mun <= 0
        % Unreachable state (no server can serve this composition): zero weight.
        continue
    end
    sPhi = 0;
    sXi = zeros(1, R);
    for s = 1:R
        if n(s) > 0
            j = i - stride(s);
            sPhi = sPhi + Phi(j);
            sXi = sXi + Xi(j, :);
        end
    end
    Phi(i) = sPhi / mun;
    for r = 1:R
        acc = sXi(r);
        if n(r) > 0
            j = i - stride(r);
            if mun > muv(j)
                acc = acc + Phi(j);       % the tail class-r job is in service
            end
        end
        Xi(i, r) = acc / mun;
    end
end

g = zeros(total, R);
for i = 1:total
    if Phi(i) > 0
        g(i, :) = Xi(i, :) / Phi(i);
    end
end
end
