function nevals = pfqn_cub_evals(M, order, Z)
% NEVALS = PFQN_CUB_EVALS(M, ORDER, Z)
%
% Number of integrand evaluations pfqn_cub performs on an M-station model at
% the given cubature order. The Grundmann-Moeller rule of degree ORDER on the
% (M-1)-simplex evaluates sum_{d=0..ORDER} nchoosek(M-1+2d, M-1) points, and a
% non-zero think time makes pfqn_cub repeat the whole rule at each of its
% v-quadrature steps. pfqn_nc prices CUB against this before selecting it.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

steps = 1e4; % must match the v-integration in pfqn_cub

n = M - 1;
nodes = 0;
for d = 0:order
    nodes = nodes + nchoosek(n + 2*d, n);
end

if nargin >= 3 && ~isempty(Z) && sum(Z(:)) >= GlobalConstants.FineTol
    nevals = nodes * steps;
else
    nevals = nodes;
end
end
