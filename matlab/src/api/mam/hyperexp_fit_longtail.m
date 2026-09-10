function result = hyperexp_fit_longtail(ccdf, varargin)
% HYPEREXP_FIT_LONGTAIL Fit a hyperexponential to a long-tail distribution.
%
% RESULT = HYPEREXP_FIT_LONGTAIL(CCDF) fits a mixture of exponentials to the
% distribution whose complementary cdf is the function handle CCDF, recursively
% over time scales.
%
% WHY MOMENTS ARE THE WRONG HANDLE. A Pareto law with tail index below 2 has
% infinite variance, so no two- or three-moment fit exists at all; and even when
% the moments are finite, matching them says nothing about the several ORDERS OF
% MAGNITUDE of time scale over which a long-tail distribution actually acts.
% This procedure matches the CCDF ITSELF at points spread across those decades.
%
% THE RECURSION. Order the components so that lambda_1 < ... < lambda_k. In the
% far tail only the slowest component survives, so (p_1,lambda_1) can be fitted
% there alone, from the ccdf at c_1 and b*c_1:
%
%   lambda_1 = ln(F^c(c_1)/F^c(b c_1))/((b-1)c_1),  p_1 = F^c(c_1)exp(lambda_1 c_1).
%
% Subtract that component from the ccdf and repeat one decade lower, and so on
% (eqs. 4.6-4.11). The last component takes whatever probability is left,
% p_k = 1 - sum_{j<k} p_j, and its rate follows from the ccdf at c_k
% (eqs. 4.12-4.14). This is Prony's method applied to a ccdf.
%
% Options:
%   'k', K        - number of components; [] takes one per decade between the
%                   0.9 quantile and the 1e-6 quantile, retrying with fewer if
%                   the recursion runs out of probability
%   'c1', C       - the largest fitting argument, default the 1e-6 quantile
%   'b', B        - the within-scale spacing, 1 < B < c_i/c_{i+1}
%   'decade', D   - the ratio between successive fitting arguments
%   'points', VEC - explicit decreasing fitting arguments, overriding c1/decade
%
% The default pair (b,decade) = (1.5,4) is not the paper's illustrative (2,10):
% the algorithm is exact AT the fitting arguments and free between them, and
% measured on a Weibull(0.3) the tighter grid cuts the worst between-point error
% from about 54% to 12%, at the cost of more components.
%
% Returns a struct with fields p, lambda, points, mean, targetMean (the original
% mean over the covered range), coverage, maxRelError (at the fitting arguments)
% and maxRelErrorGrid (on a log grid across the coverage). The last component
% matches only at c_k, its weight being fixed by the total probability, so the
% error at b*c_k is not zero by construction.
%
% Example:
%   res = hyperexp_fit_longtail(@(t) (1+t).^-1.5);   % Pareto, infinite variance
%   he = HyperExp(res.p, res.lambda);
%
% Reference: A. Feldmann, W. Whitt (1998). Fitting mixtures of exponentials to
% long-tail distributions to analyze network performance models. Performance
% Evaluation 31, 245-279, Section 4.
%
% See also HYPEREXP, MAP_FIT.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

options = struct('k', [], 'c1', [], 'b', 1.5, 'decade', 4, 'points', []);
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
b = options.b;
decade = options.decade;
if b <= 1
    line_error(mfilename, 'The spacing b must exceed 1.');
end
if decade <= b
    line_error(mfilename, ['The decade ratio must exceed the spacing b, or the fitting ' ...
        'arguments would interleave.']);
end

if isempty(options.k) && isempty(options.points)
    % One component per decade between the body and the tail: that is what the
    % spacing buys, and asking for more components than decades is exactly what
    % breaks the recursion.
    if isempty(options.c1)
        top = hyperexp_fit_longtail_quantile(ccdf, 1e-6);
    else
        top = options.c1;
    end
    body = hyperexp_fit_longtail_quantile(ccdf, 0.9);
    if body <= 0 || top <= body
        line_error(mfilename, 'the ccdf gives no usable range of time scales to fit over');
    end
    k0 = max(2, round(log(top/body)/log(decade)) + 1);
    % The recursion needs each component to dominate at its own scale. Near the
    % body of a law with a lot of mass there (a Pareto, say) that fails and the
    % remaining probability runs out; back off one component at a time until it
    % holds. Only the AUTOMATIC count retries: an explicit k that cannot be
    % fitted is an error the caller asked for.
    for kk = k0:-1:2
        try
            result = hyperexp_fit_longtail_core(ccdf, kk, top, b, decade, []);
            return
        catch
            continue
        end
    end
    line_error(mfilename, sprintf(['no component count from %d down to 2 admits the recursion; ' ...
        'the ccdf may not be long-tailed enough for this scheme'], k0));
end

result = hyperexp_fit_longtail_core(ccdf, options.k, options.c1, b, decade, options.points);
end

function result = hyperexp_fit_longtail_core(ccdf, k, c1, b, decade, points)
% The recursion at a fixed component count.
if ~isempty(points)
    cs = points(:).';
    k = numel(cs);
    if any(diff(cs) >= 0)
        line_error(mfilename, 'The fitting arguments must be strictly decreasing.');
    end
else
    k = round(k);
    if k < 1
        line_error(mfilename, 'At least one exponential component is required.');
    end
    if isempty(c1)
        c1 = hyperexp_fit_longtail_quantile(ccdf, 1e-6);
    end
    cs = c1 * decade.^-(0:k-1);
end

p = zeros(1, k);
lam = zeros(1, k);
for i = 1:k
    ci = cs(i);
    % Eqs. (4.6)-(4.7): what the already-fitted, slower components leave.
    residC = ccdf(ci) - sum(p(1:i-1) .* exp(-lam(1:i-1)*ci));
    residBC = ccdf(b*ci) - sum(p(1:i-1) .* exp(-lam(1:i-1)*b*ci));
    if i < k
        if residC <= 0 || residBC <= 0 || residC <= residBC
            line_error(mfilename, sprintf(['the residual ccdf is not positive and decreasing at ' ...
                'fitting argument %g. The recursion needs the arguments well separated, ' ...
                'c_i/c_(i+1) >> b, so that only the slowest surviving component matters at each ' ...
                'scale; widen decade, lower k, or move c1 further into the tail'], ci));
        end
        lam(i) = log(residC/residBC)/((b-1)*ci);      % eq. (4.10)
        p(i) = residC * exp(lam(i)*ci);               % eq. (4.11)
    else
        % Eqs. (4.12)-(4.14): the last component takes the rest of the mass.
        p(i) = 1 - sum(p(1:i-1));
        if p(i) <= 0
            line_error(mfilename, ['the fitted components already carry all the probability, so ' ...
                'the last one has none left; lower k or move c1 further into the tail']);
        end
        if residC <= 0
            line_error(mfilename, ['the residual ccdf has gone non-positive at the last fitting ' ...
                'argument; lower k or move c1 further into the tail']);
        end
        lam(i) = log(p(i)/residC)/ci;                 % eq. (4.14)
    end
    if lam(i) <= 0
        line_error(mfilename, sprintf(['a non-positive rate came out of the fit at argument %g; ' ...
            'the ccdf is not decaying fast enough there for this many components'], ci));
    end
end

fitted = @(t) sum(p .* exp(-lam*t));
hi = cs(1)*b;
% The target mean, for the caller to compare against: the fit is constrained
% only on [c_k, b c_1], and a mean lives wherever the body is, so a k too small
% to reach the body shows up here and nowhere else.
grid = linspace(0, hi, 20001);
result.p = p;
result.lambda = lam;
result.points = cs;
result.mean = sum(p ./ lam);
result.targetMean = trapz(grid, arrayfun(ccdf, grid));
result.coverage = [cs(end), hi];
errs = zeros(1, 2*k);
for i = 1:k
    for j = 1:2
        t = cs(i)*(1 + (j-1)*(b-1));
        target = ccdf(t);
        if target > 0
            errs((i-1)*2+j) = abs(fitted(t) - target)/target;
        end
    end
end
result.maxRelError = max(errs);
% The fit is exact at the fitting arguments by construction, so this is what
% says whether it also holds BETWEEN them.
ts = exp(linspace(log(cs(end)), log(hi), 200));
worst = 0;
for i = 1:numel(ts)
    target = ccdf(ts(i));
    if target > 1e-300
        worst = max(worst, abs(fitted(ts(i)) - target)/target);
    end
end
result.maxRelErrorGrid = worst;
end

function t = hyperexp_fit_longtail_quantile(ccdf, prob)
% Smallest t with F^c(t) <= prob, by doubling then bisection.
hi = 1;
while ccdf(hi) > prob
    hi = 2*hi;
    if hi > 1e15
        line_error(mfilename, 'the ccdf does not decay, so there is no tail to fit');
    end
end
lo = 0;
for i = 1:200
    mid = (lo + hi)/2;
    if ccdf(mid) > prob
        lo = mid;
    else
        hi = mid;
    end
end
t = (lo + hi)/2;
end
