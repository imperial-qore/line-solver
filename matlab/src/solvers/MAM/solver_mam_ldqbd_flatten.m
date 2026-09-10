function [Q, levelOf] = solver_mam_ldqbd_flatten(ld)
% SOLVER_MAM_LDQBD_FLATTEN Assemble a flat generator from LD-QBD blocks.
%
% [Q, LEVELOF] = SOLVER_MAM_LDQBD_FLATTEN(LD) takes the block-tridiagonal
% representation LD (fields Q0/Q1/Q2/Nlev as produced by solver_mam_ldqbd) and
% returns the dense infinitesimal generator Q over the level/phase state space,
% together with LEVELOF(s) = the queue level of flat state s.
%
% Level 0 is a single (empty-queue) state; levels 1..Nlev each carry nPhases
% phases (PH service) or a single state (exponential). Used by the SolverENV
% state-vector analyzer to propagate an entry distribution over the LD-QBD.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Q0 = ld.Q0; Q1 = ld.Q1; Q2 = ld.Q2;
Nlev = ld.Nlev;

% Per-level state-space size and offsets.
levelSize = zeros(1, Nlev + 1);
for n = 0:Nlev
    levelSize(n+1) = size(Q1{n+1}, 1);
end
levelStart = [0, cumsum(levelSize)];   % levelStart(n+1) = offset before level n
dim = levelStart(end);

Q = zeros(dim, dim);
levelOf = zeros(1, dim);
for n = 0:Nlev
    rows = levelStart(n+1) + (1:levelSize(n+1));
    levelOf(rows) = n;
    Q(rows, rows) = Q1{n+1};                       % within-level
    if n < Nlev
        cols = levelStart(n+2) + (1:levelSize(n+2));
        Q(rows, cols) = Q0{n+1};                   % upward (arrival)
    end
    if n >= 1
        cols = levelStart(n) + (1:levelSize(n));
        Q(rows, cols) = Q2{n};                     % downward (departure)
    end
end
end
