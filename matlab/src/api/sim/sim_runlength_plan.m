function plan = sim_runlength_plan(means, ciHalfWidth, samplesUsed, varargin)
% PLAN = SIM_RUNLENGTH_PLAN(MEANS, CIHALFWIDTH, SAMPLESUSED)
%
% How long a simulation run should have been, from the one it already did.
%
% A batch-means half-width H at confidence 1-alpha over a run of N samples pins
% the ASYMPTOTIC variance of the estimator,
%
%   sigma^2 = (H/z)^2 N,   z = Phi^-1((1+confidence)/2),
%
% and that is the quantity a run length is planned from -- NOT the stationary
% variance, which on M/M/1 differs from it by a factor blowing up like
% (1-rho)^-2. SIM_RUNLENGTH then turns it into the sample count that reaches a
% requested RELATIVE precision.
%
% MEANS and CIHALFWIDTH are matrices of the same shape, one entry per
% (station, class); an entry with a non-positive mean or half-width is left NaN,
% since there is nothing to plan from there.
%
% Options: 'relprecision' (default 0.05), 'confidence' (default 0.95, and it
% must be the level the half-widths were computed at).
%
% Returns a struct with fields relprecision, confidence, samplesUsed,
% asymptoticVariance and requiredSamples.
%
% Reference: W. Whitt (1989). Planning queueing simulations. Management Science
% 35(11), 1341-1366.
%
% See also SIM_RUNLENGTH, SIM_ASYMVAR_CTMC, SIM_ASYMVAR_MM1.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

options = struct('relprecision', 0.05, 'confidence', 0.95);
for i = 1:2:numel(varargin)
    name = lower(char(varargin{i}));
    if ~isfield(options, name)
        line_error(mfilename, sprintf('unknown option %s', char(varargin{i})));
    end
    options.(name) = varargin{i+1};
end
if samplesUsed <= 0
    line_error(mfilename, 'The number of samples already used must be positive.');
end

z = sqrt(2)*erfinv(options.confidence);
[M, K] = size(means);
plan = struct('relprecision', options.relprecision, 'confidence', options.confidence, ...
    'samplesUsed', samplesUsed, 'asymptoticVariance', nan(M,K), 'requiredSamples', nan(M,K));
for i = 1:M
    for r = 1:K
        h = ciHalfWidth(i,r);
        m = means(i,r);
        if ~isfinite(h) || h <= 0 || ~isfinite(m) || m <= 0
            continue
        end
        asymVar = (h/z)^2 * samplesUsed;
        plan.asymptoticVariance(i,r) = asymVar;
        res = sim_runlength(m, asymVar, 'relprecision', options.relprecision, ...
            'confidence', options.confidence);
        plan.requiredSamples(i,r) = res.requiredRunLength;
    end
end
end
