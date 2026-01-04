function sample = me_sample(ME, n, xs)
% SAMPLE = ME_SAMPLE(ME, N, XS) - Generate random samples from ME distribution
%
% Generate samples from a Matrix Exponential (ME) distribution using
% inverse CDF interpolation method.
%
% Input:
%   ME: ME distribution as either:
%       - ME object
%       - Process cell array {D0, D1}
%   N: Number of samples to generate (default: 1)
%   XS: Optional pre-computed grid for CDF evaluation
%       If not provided, auto-generates grid based on mean and variance
%
% Output:
%   SAMPLE: Column vector of N samples from the ME distribution
%
% Algorithm:
%   Uses inverse CDF interpolation with binary search:
%   1. Compute CDF on a fine grid (default: 1000 points)
%   2. For each sample, generate u ~ Uniform(0,1)
%   3. Use interpolation to find x such that CDF(x) = u
%
% Examples:
%   me = ME([0.3, 0.7], [-2, 1; 0.5, -1.5]);
%   samples = me_sample(me, 10000);
%
%   % Or use process representation directly
%   ME_proc = {[-2, 1; 0.5, -1.5], [0.4, 0.6; 0.3, 0.7]};
%   samples = me_sample(ME_proc, 10000);
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Handle input arguments
if nargin < 2
    n = 1;
end

% Extract process representation
if isa(ME, 'ME')
    ME_proc = ME.getProcess();
elseif iscell(ME)
    ME_proc = ME;
else
    error('ME must be either an ME object or a cell array {D0, D1}');
end

% Auto-generate grid if not provided
if nargin < 3 || isempty(xs)
    mean_val = map_mean(ME_proc);
    var_val = map_var(ME_proc);
    std_val = sqrt(var_val);

    % Create grid from 0 to mean + 10*sigma with 1000 points
    xs = linspace(0, mean_val + 10*std_val, 1000);
end

% Compute CDF at grid points
Fxs = map_cdf(ME_proc, xs);

% Ensure CDF is strictly increasing for interpolation
% (handle numerical precision issues)
for i = 2:length(Fxs)
    if Fxs(i) <= Fxs(i-1)
        Fxs(i) = Fxs(i-1) + eps;
    end
end

% Generate samples via inverse CDF interpolation
sample = zeros(n, 1);
for i = 1:n
    u = rand();

    % Handle edge cases
    if u <= Fxs(1)
        sample(i) = xs(1);
    elseif u >= Fxs(end)
        % Extrapolate beyond grid if needed
        % Use exponential tail approximation
        sample(i) = xs(end) + log(1 / (1 - u + eps));
    else
        % Linear interpolation between grid points
        sample(i) = interp1(Fxs, xs, u, 'linear');
    end
end

end
