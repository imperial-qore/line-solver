function result = sim_runlength(mean_, asymVar, varargin)
% SIM_RUNLENGTH Run length for a steady-state estimate of a given precision.
%
% RESULT = SIM_RUNLENGTH(MEAN, ASYMVAR) returns the simulated time needed to
% estimate a steady-state mean to within 5% at 95% confidence, given the
% ASYMPTOTIC VARIANCE of the estimator.
%
% A time average over [0,t] has standard error sqrt(sigma^2/t), so a two-sided
% interval of half-width z sqrt(sigma^2/t) reaches relative precision eps when
%
%   t* = (z/eps)^2 sigma^2 / mean^2.
%
% THE POINT OF THE FORMULA is that everything expensive sits in sigma^2/mean^2,
% the squared coefficient of variation of the TIME AVERAGE rather than of the
% process. Halving the tolerance quadruples the run.
%
% Options: 'relPrecision' (default 0.05), 'confidence' (default 0.95),
% 'runLength' (an actual run length, to report the precision it buys).
%
% Returns a struct with fields requiredRunLength, z and, when runLength is
% given, halfWidth and achievedRelPrecision.
%
% Reference: W. Whitt (1989). Planning queueing simulations. Management Science
% 35(11), 1341-1366.
%
% See also SIM_ASYMVAR_MM1, SIM_ASYMVAR_CTMC.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

options = struct('relprecision', 0.05, 'confidence', 0.95, 'runlength', []);
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

if mean_ == 0
    line_error(mfilename, 'A relative precision is meaningless for a zero mean.');
end
if asymVar < 0
    line_error(mfilename, 'The asymptotic variance cannot be negative.');
end
if options.relprecision <= 0
    line_error(mfilename, 'The relative precision must be positive.');
end
z = sim_runlength_z(options.confidence);
result.requiredRunLength = (z/options.relprecision)^2 * asymVar/mean_^2;
result.z = z;
if ~isempty(options.runlength) && options.runlength > 0
    result.halfWidth = z*sqrt(asymVar/options.runlength);
    result.achievedRelPrecision = result.halfWidth/abs(mean_);
end
end

function z = sim_runlength_z(confidence)
% Two-sided normal quantile, by bisection on erfc so no toolbox is needed.
if confidence <= 0 || confidence >= 1
    line_error(mfilename, 'The confidence must lie in (0,1).');
end
target = 1 - confidence;
lo = 0;
hi = 40;
for i = 1:200
    mid = (lo + hi)/2;
    if erfc(mid/sqrt(2)) > target
        lo = mid;
    else
        hi = mid;
    end
end
z = (lo + hi)/2;
end
