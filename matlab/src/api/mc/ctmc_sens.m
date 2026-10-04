function dpi = ctmc_sens(Q, dQ, pi)
% DPI = CTMC_SENS(Q, DQ, PI)
%
% Sensitivity of the steady-state distribution of a CTMC to a scalar
% parameter theta, given the generator Q, its derivative DQ = dQ/dtheta, and
% the steady-state vector PI.
%
% Differentiating the balance equations pi*Q = 0 and pi*e = 1 with respect to
% theta gives the linear system
%
%   (dpi/dtheta) * Q = -pi * (dQ/dtheta),   sum_i dpi_i/dtheta = 0,
%
% i.e. Trivedi and Bobbio (2017), Eq. (9.81). The system has the same
% coefficient matrix as the steady-state solve itself, so obtaining a
% sensitivity costs one extra solve against a matrix that is already
% assembled. The normalization replaces one column of the singular Q, exactly
% as in the steady-state solve.
%
% @param Q Generator matrix (n x n)
% @param dQ Derivative of the generator with respect to theta (n x n)
% @param pi Steady-state distribution (1 x n); computed if omitted
% @return dpi Derivative of the steady-state distribution (1 x n)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

n = size(Q, 1);
if nargin < 3 || isempty(pi)
    pi = ctmc_solve(Q);
end
pi = pi(:)';

if any(size(dQ) ~= size(Q))
    line_error(mfilename, 'dQ must have the same size as Q');
end

% Right-hand side of Eq. (9.81)
b = -pi * dQ;

% Solve dpi * Q = b subject to sum(dpi) = 0. Transpose to column form and
% replace the last equation by the normalization, mirroring ctmc_solve.
A = Q';
A(n, :) = ones(1, n);
b = b(:);
b(n) = 0;

dpi = (A \ b)';
end
