function pi = ag_birth_death_solve(Q)
% BIRTH_DEATH_SOLVE Solve equilibrium of a birth-death (tridiagonal) CTMC
%
% For a birth-death chain with birth rate lambda_n = Q(n, n+1) and
% death rate mu_n = Q(n, n-1), the equilibrium is computed using the
% recursion pi(n) = pi(n-1) * lambda(n-1) / mu(n).
%
% This is numerically stable and avoids the ill-conditioned linear system
% that plagues null-space methods for large state spaces.

n = size(Q, 1);
if n <= 1
    pi = 1;
    return;
end

pi = zeros(1, n);
pi(1) = 1.0;

for i = 2:n
    birth_rate = Q(i-1, i);
    death_rate = Q(i, i-1);
    if death_rate > 0
        pi(i) = pi(i-1) * birth_rate / death_rate;
    else
        pi(i) = 0;
    end
end

total = sum(pi);
if total > 0
    pi = pi / total;
else
    pi = ones(1, n) / n;
end

end
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
