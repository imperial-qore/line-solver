function pi = ag_solve_component(Qk, mph, nlev, lvl)
% AG_SOLVE_COMPONENT Stationary vector of one isolated component.
%
% PI = AG_SOLVE_COMPONENT(QK, MPH, NLEV, LVL)
%
% Shared by every execution backend and by the cluster fallback, so a remote
% agent and a local one are solved by the same code.
%
% A component with a single phase per level is the birth-death chain the
% analyzer has always built, and the ratio recursion is both exact and stable
% there; a phase-expanded component is block tridiagonal instead, and the
% matrix analogue of that recursion (linear level reduction) keeps the same
% stability at the 100-level truncation, where a null-space solve is already
% ill-conditioned. Anything that reaches beyond the neighbouring level -- a
% catastrophe, a batch removal -- is neither, and falls back to ctmc_solve.
if mph == 1
    if ag_is_tridiagonal(Qk)
        pi = ag_birth_death_solve(Qk);
        return;
    end
elseif ag_is_block_tridiagonal(Qk, lvl)
    pi = ag_qbd_finite_solve(Qk, mph, nlev);
    return;
end
pi = ctmc_solve(Qk);
end
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
