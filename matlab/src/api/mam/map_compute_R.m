function R = map_compute_R(C, D, mu)
% MAP_COMPUTE_R Compute rate matrix R for MAP/M/1 queue
%
% R = MAP_COMPUTE_R(C, D, mu) computes the matrix R, which is the minimal
% nonnegative solution of the matrix equation:
%   D + R(C - mu*I) + mu*R^2 = O
%
% This matrix is used in the analysis of MAP/M/1 queues based on
% quasi-birth-death processes.
%
% Input:
%   C   - M x M matrix governing MAP transitions without arrivals
%   D   - M x M matrix governing MAP transitions with arrivals
%   mu  - Service rate (scalar)
%
% Output:
%   R   - M x M rate matrix (minimal nonnegative solution)
%
% Reference:
%   Masuyama, H., & Takine, T. (2003). Sojourn time distribution in a
%   MAP/M/1 processor-sharing queue. Operations Research Letters, 31(6),
%   406-412.
%
% See also: MAP_M1PS_SOJOURN, MAP_M1PS_H_RECURSIVE

% Copyright (c) 2012-2025, Imperial College London
% All rights reserved.

M = size(C, 1);
I = eye(M);

% Use iterative method to solve for R: R = -D(C - mu*I + mu*R)^{-1}
% Start with R = O
R = zeros(M);
max_iter = 1000;
tol = 1e-10;

for iter = 1:max_iter
    R_old = R;
    % R_{k+1} = -D * inv(C - mu*I + mu*R_k)
    R = -D / (C - mu*I + mu*R);

    % Check convergence
    if norm(R - R_old, inf) < tol
        break;
    end
end

if iter == max_iter
    warning('MAP_COMPUTE_R:NoConvergence', ...
        'R computation did not converge within %d iterations', max_iter);
end

% Ensure non-negativity (numerical errors might produce small negative values)
R(R < 0) = 0;

end
