function result = qsys_maxima_twomoment(n, mean_, cs2, varargin)
% QSYS_MAXIMA_TWOMOMENT Two-moment approximation for the maximum of n variables.
%
% RESULT = QSYS_MAXIMA_TWOMOMENT(N, MEAN, CS2) approximates the mean of the
% maximum of N iid non-negative variables with the given mean and squared
% coefficient of variation.
%
% THE SHAPE OF THE ANSWER. For a law with an exponential-like tail the maximum of
% N samples grows like c~^2 (log N + ...): doubling N adds a constant, it does
% not scale the answer. What the two moments buy is the SLOPE c~^2 of that
% logarithm and an offset eta:
%
%   x_n(q) = c~^2 [log(n eta) - log log(1/q)],  E[M_n] = c~^2 [log(n eta) + gamma]
%
% with, for CS2 >= 1, c~^2 = CS2 and eta = (CS2+1)/(2 CS2^2) from the H2
% representative, and for CS2 < 1 the shifted-exponential representative
% c~^2 = sqrt(CS2), eta = exp((1-sqrt(CS2))/sqrt(CS2)).
%
% WHEN NOT TO USE IT. The extreme-value form needs N past a threshold
% n* ~ CS2/q, because with a highly variable law most of the N samples come from
% the short component and only about N p of them can contend for the maximum.
% Measured against exact maxima, the closed form is within a few percent for
% N >= 100 at CS2 = 4 and 16, and useless at N = 10 for CS2 = 16 -- which is
% exactly what n* predicts. RESULT.reliable reports the test.
%
% AND WHEN TWO MOMENTS ARE NOT ENOUGH. Below CS2 = 1 the maximum is genuinely
% family-dependent: an Erlang and a shifted exponential with the same two moments
% have maxima that differ by tens of percent and diverge as N grows, because
% their tails decay at different rates. Measured on Erlang-4, both the closed
% form and the fitted shifted exponential are 15-23% high, and they agree with
% each other, so the gap is the model's, not the arithmetic's.
%
% Options: 'q' (a quantile level in (0,1); the mean is returned without it),
% 'exactFitted' (also compute the maximum exactly from the fitted representative
% by integrating 1-F^n, which is the paper's other recommendation; default true).
%
% Returns a struct with fields value, slope, eta, threshold, reliable, family and
% exactFittedValue.
%
% Reference: C. Crow, D. Goldberg, W. Whitt (2007). Two-moment approximations
% for maxima. Operations Research 55(3), 532-548.
%
% See also FJ_RMAX, FJ_XMAX_ERLANG.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

options = struct('q', [], 'exactfitted', true);
for i = 1:2:numel(varargin)
    if i+1 > numel(varargin)
        line_error(mfilename, sprintf('option %s has no value', char(varargin{i})));
    end
    name = lower(char(varargin{i}));
    if ~isfield(options, name)
        line_error(mfilename, sprintf('unknown option %s', char(varargin{i})));
    end
    options.(name) = varargin{i+1};
end
if n < 1
    line_error(mfilename, 'At least one sample is required.');
end
if mean_ <= 0
    line_error(mfilename, 'The mean must be positive.');
end
if cs2 <= 0
    line_error(mfilename, 'The squared coefficient of variation must be positive.');
end
q = options.q;
if ~isempty(q) && (q <= 0 || q >= 1)
    line_error(mfilename, 'The quantile level must lie in (0,1).');
end

EULER = 0.5772156649015329;
if cs2 >= 1
    ct = cs2;
    eta = (cs2 + 1)/(2*cs2^2);
    family = 'H2';
else
    ct = sqrt(cs2);
    eta = exp((1 - sqrt(cs2))/sqrt(cs2));
    family = 'shifted exponential';
end
if isempty(q)
    inner = log(n*eta) + EULER;
    qq = 0.5;
else
    inner = log(n*eta) - log(log(1/q));
    qq = q;
end
result.value = mean_*ct*inner;
result.slope = mean_*ct;
result.eta = eta;
result.threshold = cs2/qq;                  % eq. (4.22)
result.reliable = n >= result.threshold;
result.family = family;

if options.exactfitted
    % The paper's other recommendation: fit the representative law, then compute
    % the maximum exactly from F^n rather than from its tail.
    if cs2 >= 1
        p1 = 0.5*(1 + sqrt((cs2-1)/(cs2+1)));       % H2, balanced means
        l1 = 2*p1/mean_;
        l2 = 2*(1-p1)/mean_;
        ccdf = @(t) p1*exp(-l1*t) + (1-p1)*exp(-l2*t);
        hi = 40*mean_*max(cs2, 1);
    else
        d = mean_*(1 - sqrt(cs2));
        m = mean_*sqrt(cs2);
        ccdf = @(t) (t <= d) + (t > d).*exp(-(t - d)/m);
        hi = d + 40*m;
    end
    grid = linspace(0, hi, 200001);
    cdfn = (1 - ccdf(grid)).^n;
    if isempty(q)
        result.exactFittedValue = trapz(grid, 1 - cdfn);
    else
        idx = find(cdfn >= q, 1, 'first');
        if isempty(idx)
            idx = numel(grid);
        end
        result.exactFittedValue = grid(idx);
    end
end
end
