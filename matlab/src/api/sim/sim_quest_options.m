function opt = sim_quest_options(options)
% OA_QUEST_OPTIONS Validate and complete the option struct of the QUEST procedures.
%
% OPT = OA_QUEST_OPTIONS(OPTIONS) fills the fields OA_FQUEST and OA_FIRQUEST
% recognize with the published FQUEST defaults and rejects unknown or
% inadmissible ones. See OA_FQUEST for what each field controls.
%
% The defaults b0 = 50, m0 = 500, s = [32 24 16 10], beta = 0.30, eta = 0.2 and
% theta = 2.3 are the ones the article reports after its own experimentation:
% b0 = 50 gives the warmup randomness test enough power, 32 batches suffice to
% estimate the variance parameter while fewer than 10 make the interval
% unreliable, and the decaying warmup significance keeps the batch size from
% growing so far that truncation eats a short sample. With these values the
% fourth warmup iteration runs at beta*exp(-0.2*3^2.3) = 0.025.
%
% See also OA_FQUEST, OA_FIRQUEST
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 1 || isempty(options)
    options = struct();
end
if ~isstruct(options) || ~isscalar(options)
    line_error(mfilename, 'options must be a scalar struct');
end

opt = struct('b0', 50, 'm0', 500, 's', [32 24 16 10], 'beta', 0.30, ...
    'eta', 0.2, 'theta', 2.3, 'weight', sqrt(12), 'force', true);

known = fieldnames(opt);
given = fieldnames(options);
for i = 1:numel(given)
    if ~any(strcmp(known, given{i}))
        line_error(mfilename, 'unknown option ''%s''', given{i});
    end
    opt.(given{i}) = options.(given{i});
end

if ~isscalar(opt.b0) || opt.b0 ~= floor(opt.b0) || opt.b0 < 3
    line_error(mfilename, 'b0 must be an integer >= 3');
end
if ~isscalar(opt.m0) || opt.m0 ~= floor(opt.m0) || opt.m0 < 1
    line_error(mfilename, 'm0 must be a positive integer');
end
opt.s = opt.s(:)';
if isempty(opt.s) || any(opt.s ~= floor(opt.s)) || any(opt.s < 1)
    line_error(mfilename, 's must be a nonempty vector of positive integers');
end
if any(diff(opt.s) >= 0)
    line_error(mfilename, 's must be strictly decreasing');
end
if ~isscalar(opt.beta) || ~isreal(opt.beta) || opt.beta <= 0 || opt.beta >= 1
    line_error(mfilename, 'beta must be a real scalar in (0,1)');
end
if ~isscalar(opt.eta) || ~isreal(opt.eta) || opt.eta < 0
    line_error(mfilename, 'eta must be a nonnegative real scalar');
end
if ~isscalar(opt.theta) || ~isreal(opt.theta) || opt.theta <= 0
    line_error(mfilename, 'theta must be a positive real scalar');
end
if ~isscalar(opt.weight) || ~isreal(opt.weight) || opt.weight == 0
    line_error(mfilename, 'weight must be a nonzero real scalar');
end
if ~isscalar(opt.force) || ~(islogical(opt.force) || isnumeric(opt.force))
    line_error(mfilename, 'force must be a logical scalar');
end
opt.force = logical(opt.force);
end
