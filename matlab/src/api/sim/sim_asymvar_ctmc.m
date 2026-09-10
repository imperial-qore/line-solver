function result = sim_asymvar_ctmc(A, f, pi)
% SIM_ASYMVAR_CTMC Asymptotic variance of a reward on a CTMC.
%
% RESULT = SIM_ASYMVAR_CTMC(A, F) returns sigma^2 = 2 sum_x pi(x) g(x) d(x),
% where g = F - E_pi[F] and d solves A d = -g: the DEVIATION vector, the
% accumulated future excess reward started from each state.
%
% This is the general form of the quantity SIM_ASYMVAR_MM1 gives in closed form
% for M/M/1, and it is what run-length planning needs for any model LINE can
% build a generator for. A alone is singular: a constant may be added to d
% without changing the variance, so one normalization is needed to pin it.
%
% Returns a struct with fields mean, variance, asymptoticVariance and deviation.
%
% Reference: W. Whitt (1989). Planning queueing simulations. Management Science
% 35(11), 1341-1366.
%
% See also SIM_RUNLENGTH, SIM_ASYMVAR_MM1, CTMC_SOLVE.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

n = size(A, 1);
if size(A, 2) ~= n
    line_error(mfilename, 'The generator must be square.');
end
if max(abs(sum(A, 2))) > 1e-8
    line_error(mfilename, 'The generator rows must sum to zero.');
end
f = f(:);
if numel(f) ~= n
    line_error(mfilename, 'One reward per state is required.');
end
if nargin < 3 || isempty(pi)
    M = [A.'; ones(1, n)];
    M(n, :) = [];
    M = [M; ones(1, n)];
    b = zeros(n, 1);
    b(end) = 1;
    pi = M \ b;
end
pi = pi(:).';
mean_ = pi*f;
g = f - mean_;
% A d = -g pins d only up to a constant, so one equation of A is redundant and
% one normalization replaces it. WHICH equation is dropped matters: the rows of A
% are related by pi A = 0, so a row whose pi is tiny is only nominally redundant
% and dropping it loses real information -- on a queue truncated where pi has
% underflowed, that alone puts sigma^2 out by orders of magnitude. Dropping the
% row with the LARGEST pi is the well-conditioned choice.
[~, drop] = max(pi);
keep = setdiff(1:n, drop);
d = [A(keep,:); pi] \ [-g(keep); 0];
result.mean = mean_;
result.variance = pi*(g.^2);
result.asymptoticVariance = 2*(pi*(g.*d));
result.deviation = d;
end
