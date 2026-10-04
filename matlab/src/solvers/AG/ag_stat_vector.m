function v = ag_stat_vector(C)
% Stationary vector of a generator C, allowing a reducible one.
%
% Replaces the first balance equation by the normalization, which is the
% equation it is redundant with (the columns of a generator sum to zero), and
% solves the resulting square system. Unlike a null-space solve this stays well
% posed when the chain is reducible with ONE closed class, which the level-0
% chain of a phase-expanded component routinely is: a phase-type restarts in the
% support of alpha, so every service phase outside that support is unreachable
% once the queue has emptied at least once.
m = size(C, 1);
Amat = C;
Amat(:, 1) = 1;
rhs = zeros(1, m);
rhs(1) = 1;
v = rhs / Amat;
end
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
