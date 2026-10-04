function pi = ag_qbd_finite_solve(Q, m, nlev)
% QBD_FINITE_SOLVE Stationary vector of a finite block-tridiagonal generator.
%
% Linear level reduction: censor the chain level by level from the top,
%   C(nlev-1) = B(nlev-1),  C(n) = B(n) + F(n) * (-C(n+1))^-1 * D(n+1),
% with B, F and D the diagonal, up and down blocks. C(0) is the generator of
% the chain censored on level 0, so pi_0 is its stationary vector and the rest
% follows from pi_(n+1) = pi_n F(n) (-C(n+1))^-1. This is the block form of
% BIRTH_DEATH_SOLVE and reduces to it entry for entry when m == 1.

if nlev <= 1
    pi = ctmc_solve(Q);
    return;
end

C = cell(1, nlev);
C{nlev} = Q(ag_blk(nlev - 1, m), ag_blk(nlev - 1, m));
for n = (nlev - 2):-1:0
    F = Q(ag_blk(n, m), ag_blk(n + 1, m));
    D = Q(ag_blk(n + 1, m), ag_blk(n, m));
    C{n + 1} = Q(ag_blk(n, m), ag_blk(n, m)) + F * ((-C{n + 2}) \ D);
end

pi = zeros(1, nlev * m);
pi(ag_blk(0, m)) = ag_stat_vector(C{1});
for n = 0:(nlev - 2)
    F = Q(ag_blk(n, m), ag_blk(n + 1, m));
    pi(ag_blk(n + 1, m)) = (pi(ag_blk(n, m)) * F) / (-C{n + 2});
end

total = sum(pi);
if total > 0
    pi = pi / total;
else
    pi = ones(1, nlev * m) / (nlev * m);
end
end
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
