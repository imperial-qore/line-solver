function MAP = map_bernstein(f, n)
% MAP = MAP_BERNSTEIN(F, N) - Convert distribution to MAP via Bernstein approximation
%
% This function implements Bernstein polynomial approximation to convert any
% continuous distribution (specified by its PDF) to a Markovian Arrival Process
% (MAP) representation.
%
% Input:
%   f: PDF function handle @(x) - probability density function
%   n: Number of phases for the approximation (default: 20)
%
% Output:
%   MAP: MAP representation {D0, D1}
%        Note: Caller must rescale to target mean using map_scale
%
% Example:
%   % Convert a gamma distribution with shape=2, scale=1
%   pdf_func = @(x) gampdf(x, 2, 1);
%   MAP = map_bernstein(pdf_func, 20);
%   MAP = map_scale(MAP, targetMean);  % Rescale to desired mean
%
% Reference:
%   Bernstein polynomial approximation for phase-type distributions
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    n = 20;
end

% Bernstein approximation of order n
c = 0;
for i = 1:n
    xi = -log(i/n);
    fi = f(xi);
    if isfinite(fi) && fi > 0
        c = c + fi / i;
    end
end

% Handle case where normalizing constant is invalid
if c <= 0 || ~isfinite(c)
    % Fallback to Erlang-n with unit mean (caller will rescale)
    MAP = map_erlang(1, n);
    return;
end

% Build subgenerator T
T = diag(-[1:n]) + diag([1:(n-1)], 1);

% Build initial probability vector alpha
alpha = zeros(1, n);
for i = 1:n
    xi = -log(i/n);
    fi = f(xi);
    if isfinite(fi) && fi > 0
        alpha(i) = fi / (i * c);
    end
end

% Normalize alpha to ensure it sums to 1
if sum(alpha) > 0
    alpha = alpha / sum(alpha);
else
    alpha(1) = 1; % Default to starting in first phase
end

% Conversion to MAP
P = repmat(alpha, n, 1);
MAP = {T, -T * P};

end
