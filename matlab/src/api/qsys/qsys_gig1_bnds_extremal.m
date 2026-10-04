function result = qsys_gig1_bnds_extremal(lambda, mu, ca, cs, varargin)
% QSYS_GIG1_BNDS_EXTREMAL Extremal two-moment bounds for the GI/GI/1 queue.
%
% RESULT = QSYS_GIG1_BNDS_EXTREMAL(LAMBDA, MU, CA, CS) returns the tightest known
% interval for the mean steady-state waiting time of a GI/GI/1 queue known ONLY
% through the first two moments of its interarrival and service times, together
% with the classical bounds it improves on.
%
% WHAT THE INTERVAL MEANS. Two moments do not determine E[W]; they determine a
% SET of possible values, and the width of that set is the honest uncertainty in
% any two-moment approximation. The extremal distributions are the ones that
% attain its ends:
%   - the lower end is attained by deterministic interarrival times and a
%     three-point service law concentrated on multiples of that interval, and
%     has the closed form rho((1+cs^2)rho - 1)^+ / (2(1-rho)) (eq. 2.12);
%   - the upper end is attained (asymptotically) by TWO-POINT laws: an
%     interarrival law with an atom at 0, and a service law whose upper atom
%     runs off to infinity while its probability vanishes. Making an
%     interarrival time larger only empties the queue once, but making a service
%     time larger delays every customer behind it, which is why the two ends
%     look so different.
%
% HOW THE UPPER END IS COMPUTED. Chen and Whitt reduce that limit to a
% D(1/p)/RS(D(rho),p)/1 model with p = 1/(1+ca^2) and RS a geometric random sum,
% then evaluate its mean waiting time by Spitzer's identity with the negative
% binomial pmf (their Algorithm 1). The sum is truncated in both indices, so this
% bound is a numerical limit, not a formula. The closed-form companion (eq. 3.4)
% uses the D/M/1 root delta = exp(-(1-delta)/rho) and is within about 1% of it.
%
% Options:
%   'K', K   - truncation of the negative binomial value, default 4000
%   'N', N   - truncation of the random-walk length, default 2000
%   'skipTight', TF - skip the O(K*N) tight bound and return the closed forms
%                     only, default false
%
% Returns a struct whose waiting-time fields are TIMES IN QUEUE (add 1/MU for a
% response time), all for the given LAMBDA:
%   trafficIntensity - rho = lambda/mu
%   lowerBound       - the tight lower bound, eq. (2.12)
%   upperBound       - the conjectured tight upper bound, eq. (3.2) by Algorithm 1
%   upperBoundClosed - the closed-form upper bound, eq. (3.4)
%   upperBoundDaley  - Daley's bound, eq. (2.7)
%   upperBoundKingman- Kingman's bound, eq. (2.6)
%   heavyTraffic     - the heavy-traffic approximation, eq. (2.9)
%   delta            - the D/M/1 root behind upperBoundClosed
%   relativeWidth    - (upperBound-lowerBound)/upperBound, the fraction of the
%                      answer that two moments genuinely leave undetermined
%
% Example:
%   res = qsys_gig1_bnds_extremal(0.5, 1, 2, 2);   % rho = 0.5, ca^2 = cs^2 = 4
%   [res.lowerBound res.upperBound res.upperBoundKingman]   % 0.750 3.470 5.000
%
% Reference: Y. Chen, W. Whitt (2020). Algorithms for the upper bound mean
% waiting time in the GI/GI/1 queue. Queueing Systems 94, 327-356. The lower
% bound is classical, restated there as eq. (2.12); Kingman (1962) and Daley
% (1977) are the two established upper bounds.
%
% See also QSYS_GIG1_LBND, QSYS_GIG1_UBND_KINGMAN, QSYS_GIG1_RQ.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

options = struct('k', 4000, 'n', 2000, 'skiptight', false);
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

if lambda <= 0 || mu <= 0
    line_error(mfilename, 'The arrival and service rates must be positive.');
end
rho = lambda / mu;
if rho >= 1
    line_error(mfilename, 'The bounds require a stable queue, rho < 1.');
end
ca2 = ca^2;
cs2 = cs^2;

% The reference sets E[U] = 1; every waiting time below therefore carries the
% factor 1/lambda, which is that time unit expressed in the caller's units.
scale = 1/lambda;

result.trafficIntensity = rho;
result.lowerBound = scale * rho*max((1+cs2)*rho - 1, 0)/(2*(1-rho));            % (2.12)
result.upperBoundKingman = scale * rho^2*(ca2/rho^2 + cs2)/(2*(1-rho));         % (2.6)
result.upperBoundDaley = scale * rho^2*((2-rho)*ca2/rho + cs2)/(2*(1-rho));     % (2.7)
result.heavyTraffic = scale * rho^2*(ca2 + cs2)/(2*(1-rho));                    % (2.9)

delta = qsys_gig1_bnds_extremal_delta(rho);                                     % (3.5)
result.delta = delta;
result.upperBoundClosed = scale * (2*(1-rho)*rho/(1-delta)*ca2 + rho^2*cs2)/(2*(1-rho)); % (3.4)

if options.skiptight
    result.upperBound = result.upperBoundClosed;
    result.tightComputed = false;
else
    result.upperBound = scale * qsys_gig1_bnds_extremal_tight(rho, ca2, cs2, options.k, options.n);
    result.tightComputed = true;
end
if result.upperBound > 0
    result.relativeWidth = (result.upperBound - result.lowerBound)/result.upperBound;
else
    result.relativeWidth = 0;
end
end

function delta = qsys_gig1_bnds_extremal_delta(rho)
% The D/M/1 root of eq. (3.5), delta = exp(-(1-delta)/rho), in (0,1).
%
% g(delta) = delta - exp(-(1-delta)/rho) is negative at 0 and positive just below
% 1, where the second root delta = 1 sits, so bisection on [0,1) finds the one
% that is wanted without ever landing on the trivial root.
lo = 0;
hi = 1 - 1e-15;
for i = 1:200
    mid = (lo + hi)/2;
    if mid - exp(-(1-mid)/rho) < 0
        lo = mid;
    else
        hi = mid;
    end
end
delta = (lo + hi)/2;
end

function EW = qsys_gig1_bnds_extremal_tight(rho, ca2, cs2, K, N)
% Algorithm 1 of the reference: the mean waiting time of the extremal model,
%
%   E[W(F0,Gu*)] = rho*ca2 + rho^2 cs2/(2(1-rho)) + E[W(D(1/p),RS(D(rho),p))],
%
% the last term by Spitzer's identity sum_n E[Sn^+]/n with
% Sn = rho(NB(n,1-p)+n) - n/p, evaluated from the negative binomial pmf, which
% is formed in LOG space so that neither the factorials nor the products
% overflow at the thousands of terms the truncation needs.
p = 1/(1 + ca2);
EW = rho*ca2 + rho^2*cs2/(2*(1-rho));
n = (1:N);
lognfac = gammaln(n);
total = 0;
for k = 1:K
    % log P(NB(n,1-p)=k) = lgamma(n+k)-lgamma(k+1)-lgamma(n)+n log p+k log(1-p).
    % DIVERGENCE from Algorithm 1, which steps the pmf by the ratio
    % P(n+1)/P(n) = ((n+k)/n)p: the same numbers, but a K*N scalar recursion is
    % slow here, while that recursion cannot be vectorized as a cumulative
    % product without overflowing before its (1-p)^k factor tames it. The JAR
    % and the C++ port keep the recursion, where loops are cheap.
    logp = gammaln(n + k) - gammaln(k + 1) - lognfac + n*log(p) + k*log1p(-p);
    step = max((n + k)*rho - n/p, 0);
    nz = step > 0;
    if any(nz)
        total = total + sum(exp(logp(nz)) .* step(nz) ./ n(nz));
    end
end
EW = EW + total;
end
