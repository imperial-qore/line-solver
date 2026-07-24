function result = qsys_phm1(alpha, T, mu)
% QSYS_PHM1 Exact PH/M/1 (GI/M/1 with phase-type inter-arrivals).
%
% RESULT = QSYS_PHM1(ALPHA, T, MU) solves the GI/M/1 sigma-root
%   sigma = psi_A(mu*(1 - sigma)),     psi_A(s) = ALPHA*(s*I - T)\(-T*1)
% and returns time-average performance metrics.
%
% Inputs:
%   ALPHA - PH entry probability row vector (1 x k)
%   T     - PH sub-generator (k x k)
%   MU    - Exponential service rate (positive scalar)
%
% RESULT struct fields:
%   meanQueueLength  - L = rho/(1 - sigma)
%   meanWaitingQueue - Lq = rho*sigma/(1 - sigma)
%   meanWaitingTime  - Wq = Lq/lambda
%   meanSojournTime  - W = Wq + 1/MU
%   utilization      - rho = lambda/MU
%   sigma            - GI/M/1 root in (0, 1)
%   analyzer         - identifier string
%
% See also qsys_gm1, qsys_dmc

if mu <= 0
    line_error(mfilename, 'Service rate MU must be positive');
end
alpha = alpha(:)';     % row vector
k = size(T, 1);
if size(T, 2) ~= k
    line_error(mfilename, 'T must be square');
end
if numel(alpha) ~= k
    line_error(mfilename, 'ALPHA length must match T dimension');
end

ones_k = ones(k, 1);
t_vec = -T * ones_k;          % column

% Mean inter-arrival = alpha * (-T)\1
mean_ia = alpha * ((-T) \ ones_k);
if mean_ia <= 0
    line_error(mfilename, sprintf('Non-positive mean inter-arrival: %g', mean_ia));
end
lambda = 1.0 / mean_ia;
rho = lambda / mu;
if rho >= 1 - 1e-12
    line_error(mfilename, sprintf('Load rho=%g must be strictly less than 1', rho));
end

I = eye(k);
lst = @(s) alpha * ((s*I - T) \ t_vec);
f   = @(sigma) sigma - lst(mu*(1 - sigma));

try
    sigma = fzero(f, [1e-12, 1 - 1e-12]);
catch
    % Fall back to fixed-point iteration if bracket fails
    sigma = 0.5;
    for it = 1:2000
        sNew = lst(mu*(1 - sigma));
        if abs(sNew - sigma) < 1e-13
            sigma = sNew;
            break;
        end
        sigma = sNew;
    end
end

L = rho / (1 - sigma);
Lq = rho * sigma / (1 - sigma);
Wq = Lq / lambda;
W = Wq + 1.0 / mu;

result = struct();
result.meanQueueLength = L;
result.meanWaitingQueue = Lq;
result.meanWaitingTime = Wq;
result.meanSojournTime = W;
result.utilization = rho;
result.sigma = sigma;
result.analyzer = 'PHM1:fzero';

end
