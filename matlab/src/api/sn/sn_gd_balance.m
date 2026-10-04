function [viol, nworst] = sn_gd_balance(phi, cutoffs)
% [VIOL, NWORST] = SN_GD_BALANCE(PHI, CUTOFFS)
%
% Worst relative violation of the Whittle balance property by a globally
% state-dependent rate scaling PHI (the handle declared through
% setGlobalDependence). For every state n of the lattice 0..CUTOFFS and every
% pair of stations (s,t) populated in n, the property requires
%
%   phi_s(n) phi_t(n-e_s) = phi_t(n) phi_s(n-e_t).
%
% When it holds, the chain is reversible with pi(n) ~ Phi(n) prod rho^n for the
% balance function Phi implied by phi, and the stationary law is insensitive to
% the service-time distribution beyond its mean. When it fails, the model is
% still solvable by SolverCTMC but has no product form and is sensitive.
%
% PHI is evaluated on an (nstations x 1) population column, i.e. the
% single-class reading of the (nstations x nclasses) contract, and must return a
% scalar or an (nstations x 1) column. CUTOFFS is a scalar (same bound at every
% station) or an (nstations x 1) vector.
%
% VIOL is the worst relative violation, NWORST the state attaining it.
%
% Reference: P. Whittle, "Partial balance and insensitivity", J. Appl. Prob.
% 22(1), 1985; T. Bonald, A. Proutiere, "Insensitivity in processor-sharing
% networks", Perf. Eval. 49, 2002.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~isa(phi,'function_handle')
    line_error(mfilename, 'PHI must be a function handle.');
end
cutoffs = cutoffs(:);
S = numel(cutoffs);
if S == 1
    line_error(mfilename, 'CUTOFFS must have one entry per station (at least two stations are needed for a balance pair).');
end

viol = 0;
nworst = zeros(S,1);
base = cutoffs + 1;
for idx = 0:prod(base)-1
    n = zeros(S,1);
    rem = idx;
    for s = 1:S
        n(s) = mod(rem, base(s));
        rem = floor(rem / base(s));
    end
    for s = 1:S
        if n(s) == 0
            continue
        end
        for t = s+1:S
            if n(t) == 0
                continue
            end
            es = zeros(S,1); es(s) = 1;
            et = zeros(S,1); et(t) = 1;
            xn = evalphi(phi, n, S);
            xs = evalphi(phi, n-es, S);
            xt = evalphi(phi, n-et, S);
            lhs = xn(s) * xs(t);
            rhs = xn(t) * xt(s);
            scale = max(abs([lhs, rhs]));
            if scale > 0
                v = abs(lhs-rhs) / scale;
                if v > viol
                    viol = v;
                    nworst = n;
                end
            end
        end
    end
end
end

function v = evalphi(phi, n, S)
v = phi(n);
if isscalar(v)
    v = v * ones(S,1);
else
    v = v(:);
end
if numel(v) ~= S
    line_error(mfilename, sprintf('PHI must return a scalar or a vector of length %d.', S));
end
end
