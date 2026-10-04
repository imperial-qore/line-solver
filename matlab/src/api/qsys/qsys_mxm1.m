function [W, Wq, U, Q] = qsys_mxm1(lambda_batch, mu, E_X, E_X2_or_Var_X, varargin)
% QSYS_MXM1 - MX/M/1 queue with batch arrivals
%
% [W, Wq, U, Q] = QSYS_MXM1(LAMBDA_BATCH, MU, E_X, E_X2)
% [W, Wq, U, Q] = QSYS_MXM1(LAMBDA_BATCH, MU, BATCH_SIZES, PMF)
%
% Analytical solution for MX/M/1 queueing system with batch Markovian arrivals
% and exponential service.
%
% Input formats:
%   1. Moment-based: qsys_mxm1(lambda_batch, mu, E_X, E_X2)
%      - lambda_batch: Batch arrival rate
%      - mu: Service rate
%      - E_X: Mean batch size
%      - E_X2: Second moment of batch size
%
%   2. PMF-based: qsys_mxm1(lambda_batch, mu, batch_sizes, pmf)
%      - lambda_batch: Batch arrival rate
%      - mu: Service rate
%      - batch_sizes: Array of batch sizes (e.g., [1, 2, 4, 8])
%      - pmf: Probability mass function for batch sizes
%
%   3. Variance-based: qsys_mxm1(lambda_batch, mu, E_X, Var_X, 'variance')
%      - lambda_batch: Batch arrival rate
%      - mu: Service rate
%      - E_X: Mean batch size
%      - Var_X: Variance of batch size
%      - 'variance': Flag to indicate variance mode
%
% Outputs:
%   W - Mean time in system
%   Wq - Mean waiting time in queue
%   U - Server utilization
%   Q - Mean queue length (including service)
%
% The formula accounts for both queueing delay and internal batch delay:
%   Wq = rho/(mu*(1-rho)) + (E[X²] - E[X])/(2*mu*E[X]*(1-rho))
%        ^M/M/1 term       ^Internal batch delay
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Determine which input format is used
if nargin >= 5 && ischar(varargin{1}) && strcmpi(varargin{1}, 'variance')
    % Format 3: Variance-based
    E_X = E_X;
    Var_X = E_X2_or_Var_X;
    E_X2 = Var_X + E_X^2;
elseif nargin == 4 && isvector(E_X) && length(E_X) > 1
    % Format 2: PMF-based (E_X is actually batch_sizes, E_X2_or_Var_X is pmf)
    batch_sizes = E_X;
    pmf = E_X2_or_Var_X;

    % Validate inputs
    if length(batch_sizes) ~= length(pmf)
        line_error(mfilename, 'Batch sizes and PMF must have the same length');
    end

    % Normalize PMF
    pmf = pmf / sum(pmf);

    % Compute moments
    E_X = sum(batch_sizes .* pmf);
    E_X2 = sum(batch_sizes.^2 .* pmf);
else
    % Format 1: Moment-based (default)
    E_X2 = E_X2_or_Var_X;
end

% Compute effective job arrival rate
lambda = lambda_batch * E_X;

% Compute utilization
rho = lambda / mu;

% Check stability condition
if rho >= 1
    line_error(mfilename, sprintf('System is unstable: rho = %.6f >= 1', rho));
end

% Mean waiting time in queue includes:
% 1. M/M/1 queueing delay: rho/(mu*(1-rho))
% 2. Internal batch delay: (E[X²] - E[X])/(2*mu*E[X]*(1-rho))
Wq = rho / (mu * (1 - rho)) + (E_X2 - E_X) / (2 * mu * E_X * (1 - rho));

% Mean time in system
W = Wq + 1/mu;

% Server utilization
U = rho;

% Mean queue length (Little's Law)
Q = lambda * W;

end
