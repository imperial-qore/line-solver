function sample = me_sample(ME, n, xs)
% SAMPLE = ME_SAMPLE(ME, N, XS) - Generate random samples from ME distribution
%
% Generate N independent samples from a Matrix Exponential (ME) distribution
% by numerical inversion of its exact CDF.
%
% Input:
%   ME: ME distribution as either:
%       - ME object
%       - Process cell array {D0, D1}
%   N: Number of samples to generate (default: 1)
%   XS: Optional pre-computed grid for CDF evaluation. When supplied, the
%       grid is used as given and no horizon extension is performed, so the
%       caller is responsible for covering the tail.
%
% Output:
%   SAMPLE: Column vector of N samples from the ME distribution
%
% Algorithm:
%   Inversion of F(t) = 1 - alpha*expm(A*t)*e, which is valid for every ME
%   representation. The CTMC walk of map_sample is not applicable here: it
%   presumes a phase-type reading of the representation, which fails when
%   alpha has negative entries or A has negative off-diagonal entries.
%     1. The horizon is doubled until the survival function falls below
%        TAILTOL, so the tabulated range covers all but a negligible mass.
%        A fixed horizon of mean + 10*sigma is not enough for a heavy tail:
%        on an MMPP with SCV 4.21 it leaves 5e-04 of the mass untabulated.
%     2. The CDF is tabulated on that horizon and forced nondecreasing to
%        absorb roundoff.
%     3. Each variate is located by binary search, then polished by Newton
%        steps on the exact CDF and density, so the result is not limited by
%        the linear interpolation of the table.
%     4. Variates beyond the last tabulated CDF value are placed by
%        exponential extrapolation using the dominant (least negative)
%        eigenvalue of D0, which governs the decay of the tail. Extrapolating
%        with a unit rate instead, as this function previously did, inflates
%        the tail: the same MMPP then samples an SCV of 5.27 against an exact
%        4.21. Clamping to the grid endpoint truncates it and biases the mean
%        the other way.
%
% Note:
%   Passing a MAP or a RAP samples its stationary marginal independently,
%   which is a deliberate renewal approximation: the autocorrelation is
%   discarded. Use rap_sample to retain it.
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

TAILTOL = 1e-12;    % survival mass left beyond the horizon
MAXDBL = 40;        % cap on horizon doublings
GRIDPTS = 1000;
NEWTON = 3;

% Auto-generate grid if not provided
if nargin < 3 || isempty(xs)
    mean_val = map_mean(ME_proc);
    var_val = map_var(ME_proc);
    std_val = sqrt(max(var_val, 0));

    horizon = mean_val + 10*std_val;
    if ~isfinite(horizon) || horizon <= 0
        horizon = 1;
    end
    % Extend until the tabulated range covers all but TAILTOL of the mass.
    for k = 1:MAXDBL
        if (1 - map_cdf(ME_proc, horizon)) < TAILTOL
            break;
        end
        horizon = 2*horizon;
    end

    xs = linspace(0, horizon, GRIDPTS);
end

% Compute CDF at grid points
Fxs = map_cdf(ME_proc, xs);
Fxs = Fxs(:)';
xs = xs(:)';

% Force the tabulated CDF to be nondecreasing (roundoff can break monotonicity
% near the tail, where consecutive values differ by less than eps).
for i = 2:length(Fxs)
    if Fxs(i) < Fxs(i-1)
        Fxs(i) = Fxs(i-1);
    end
end

% Dominant eigenvalue of D0 sets the decay rate of the tail.
eta = max(real(eig(full(ME_proc{1}))));
if ~(eta < 0) || ~isfinite(eta)
    if isfinite(mean_val) && mean_val > 0
        eta = -1/mean_val;
    else
        eta = -1;
    end
end

xEnd = xs(end);
sEnd = max(1 - Fxs(end), 0);

sample = zeros(n, 1);
for i = 1:n
    u = rand();

    if u <= Fxs(1)
        sample(i) = xs(1);
        continue;
    end
    if u >= Fxs(end)
        % Exponential tail: S(x) ~ S(xEnd)*exp(eta*(x-xEnd)).
        tailProb = 1 - u;
        if sEnd <= 0 || tailProb <= 0
            sample(i) = xEnd;
        else
            sample(i) = xEnd + log(sEnd/tailProb)/(-eta);
        end
        continue;
    end

    % Binary search for the bracketing interval.
    lo = 1;
    hi = length(Fxs);
    while hi - lo > 1
        mid = floor((lo + hi)/2);
        if Fxs(mid) <= u
            lo = mid;
        else
            hi = mid;
        end
    end

    den = Fxs(lo+1) - Fxs(lo);
    if den > 0
        x = xs(lo) + (u - Fxs(lo))/den*(xs(lo+1) - xs(lo));
    else
        x = xs(lo);
    end

    % Newton polish on the exact CDF and density, kept inside the bracket.
    for k = 1:NEWTON
        fx = map_pdf(ME_proc, x);
        fx = fx(1);
        if ~(fx > 0)
            break;
        end
        err = map_cdf(ME_proc, x) - u;
        if abs(err) < 1e-14
            break;
        end
        xn = x - err/fx;
        if xn <= xs(lo) || xn >= xs(lo+1)
            break;
        end
        if abs(xn - x) <= 1e-14*max(1, abs(x))
            x = xn;
            break;
        end
        x = xn;
    end

    sample(i) = x;
end

end
