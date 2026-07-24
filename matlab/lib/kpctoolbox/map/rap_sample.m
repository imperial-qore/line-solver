function sample = rap_sample(RAP, n)
% SAMPLE = RAP_SAMPLE(RAP, N) - Generate a correlated sample path from a RAP
%
% Generate N successive inter-event times along one sample path of a
% Rational Arrival Process. A RAP has no underlying Markov chain over the
% phases to walk, so the sampler carries the conditional phase vector V
% across events instead of a phase index.
%
% Input:
%   RAP: RAP distribution as either:
%        - RAP object
%        - Process cell array {H0, H1}
%   N: Number of samples to generate (default: 1)
%
% Output:
%   SAMPLE: Column vector of N inter-event times along a single sample path
%
% Algorithm:
%   Let V be a row vector with V*e = 1, initialized to the arrival-embedded
%   equilibrium vector map_pie. For each draw:
%     1. Draw U ~ Uniform(0,1) and solve the conditional survival equation
%        S(X) = V*expm(H0*X)*e = 1-U for X, by bracket expansion followed by
%        safeguarded Newton iteration. S is monotone decreasing from 1 to 0,
%        and its derivative is S'(X) = V*expm(H0*X)*H0*e, so the conditional
%        density is -S'(X) = V*expm(H0*X)*H1*e.
%     2. Update V <- V*expm(H0*X)*H1 / (V*expm(H0*X)*H1*e).
%
%   Step 2 is what reproduces the autocorrelation of the process. Sampling
%   from the ME marginal alone (me_sample) yields independent inter-event
%   times with the correct marginal but zero autocorrelation, and the MAP
%   CTMC walk (map_sample) is invalid here because a general RAP has
%   negative off-diagonal entries, for which the walk has no probabilistic
%   interpretation. A MAP is the special case in which H0 and H1 are
%   nonnegative, and on that input this sampler agrees with map_sample.
%
% Examples:
%   rap = RAP([-2, 1; 0.5, -1.5], [0.5, 0.5; 0.5, 0.5]);
%   samples = rap_sample(rap, 10000);
%
%   % Or use process representation directly
%   RAP_proc = {[-2, 1; 0.5, -1.5], [0.5, 0.5; 0.5, 0.5]};
%   samples = rap_sample(RAP_proc, 10000);
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    n = 1;
end

% Extract process representation
if isa(RAP, 'RAP')
    RAP_proc = RAP.getProcess();
elseif iscell(RAP)
    RAP_proc = RAP;
else
    error('RAP must be either a RAP object or a cell array {H0, H1}');
end

H0 = RAP_proc{1};
H1 = RAP_proc{2};
nphases = size(H0, 1);
e = ones(nphases, 1);

% Conditional phase vector, normalized so that v*e = 1. The arrival-embedded
% equilibrium vector is the stationary choice for the first inter-event time.
v = map_pie(RAP_proc);
v = v / (v * e);

% Scale for the initial bracket and for the Newton safeguard.
meanval = map_mean(RAP_proc);
if ~isfinite(meanval) || meanval <= 0
    meanval = 1;
end

sample = zeros(n, 1);
for i = 1:n
    target = 1 - rand();  % target survival level in (0,1]

    x = invert_survival(v, H0, e, target, meanval);
    sample(i) = x;

    % Advance the conditional vector across the event.
    vnext = v * expm(H0 * x) * H1;
    mass = vnext * e;
    if mass <= 0 || ~isfinite(mass)
        % The conditional vector has lost its normalization to roundoff,
        % which can only happen at survival levels far into the tail. Restart
        % from the embedded equilibrium rather than propagate a meaningless
        % vector.
        v = map_pie(RAP_proc);
        v = v / (v * e);
    else
        v = vnext / mass;
    end
end

end

function x = invert_survival(v, H0, e, target, scale)
% Solve v*expm(H0*x)*e = target for x, with target in (0,1].

surv = @(t) v * expm(H0 * t) * e;

% Bracket the root: survival is 1 at t = 0 and decreases to 0.
lo = 0;
hi = scale;
smax = 200;
k = 0;
while surv(hi) > target && k < smax
    lo = hi;
    hi = hi * 2;
    k = k + 1;
end
if k >= smax
    x = hi;
    return;
end

% Bisection to a tight bracket, then Newton polish. Bisection alone is used
% first because the conditional density of a RAP need not be monotone, so an
% unguarded Newton step from an arbitrary start can leave the bracket.
for k = 1:60
    mid = 0.5 * (lo + hi);
    if surv(mid) > target
        lo = mid;
    else
        hi = mid;
    end
    if (hi - lo) <= 1e-12 * max(1, hi)
        break;
    end
end

x = 0.5 * (lo + hi);
for k = 1:3
    ex = expm(H0 * x);
    s = v * ex * e;
    d = v * ex * H0 * e;   % derivative of the survival function, negative
    if ~(d < 0)
        break;
    end
    xn = x - (s - target) / d;
    if xn <= lo || xn >= hi
        break;
    end
    if abs(xn - x) <= 1e-14 * max(1, abs(x))
        x = xn;
        break;
    end
    x = xn;
end

end
