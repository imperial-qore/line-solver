function result = qsys_mmk_qed_staffing(lambda, mu, target, varargin)
% QSYS_MMK_QED_STAFFING Square-root staffing of the M/M/s queue.
%
% RESULT = QSYS_MMK_QED_STAFFING(LAMBDA, MU, TARGET) returns the smallest number
% of servers whose Halfin-Whitt delay probability is at most TARGET, together
% with the quality-of-service parameter behind it.
%
% THE RULE. Invert alpha(beta) = TARGET for the server slack beta, then staff
%
%   s = ceil( a + beta sqrt(a) ),      a = lambda/mu,
%
% the SQUARE-ROOT STAFFING rule: the base a erlangs of work plus a safety
% cushion that grows only as the square root of the load. Doubling the load
% needs only sqrt(2) times the cushion, which is why large service systems can
% be both highly utilized and responsive, and why the QED regime is the one
% large systems are actually run in.
%
% RESULT = QSYS_MMK_QED_STAFFING(..., 'criterion', C) chooses what TARGET means:
%   'delay'      - P(W > 0) <= TARGET (the default)
%   'meanwait'   - E[W] <= TARGET, in time units
%   'servicelevel' - P(W <= TARGET.deadline) >= TARGET.level, TARGET being a
%                  struct with those two fields. In the QED regime the wait is
%                  exponential with rate s*mu - lambda given that it is
%                  positive, so P(W > t) = alpha(beta) exp(-(s*mu-lambda) t).
%
% RESULT = QSYS_MMK_QED_STAFFING(..., 'exact', true) then walks s up or down
% until the EXACT Erlang C measure meets the target, starting from the
% square-root answer. That costs an O(s) Erlang C evaluation per step and is
% what to use when the answer must be defensible rather than asymptotic.
%
% Returns a struct with fields:
%   numServers     - the recommended s
%   beta           - the server slack achieved, (s-a)/sqrt(s)
%   betaTarget     - the slack the target asks for, before rounding s up
%   offeredLoad    - a = lambda/mu
%   probDelay      - the QED delay probability at the recommended s
%   meanWait       - the QED mean wait at the recommended s
%   serviceLevel   - P(W <= deadline) at the recommended s, for the
%                    'servicelevel' criterion
%   exactUsed      - whether the exact Erlang C refinement was applied
%
% Example:
%   res = qsys_mmk_qed_staffing(1000, 1, 0.2);            % 20% delayed at most
%   res = qsys_mmk_qed_staffing(1000, 1, struct('deadline',0.02,'level',0.8), ...
%                               'criterion','servicelevel');
%
% Reference: S. Halfin, W. Whitt (1981). Heavy-traffic limits for queues with
% many exponential servers. Operations Research 29(3), 567-588. The staffing
% form is the standard reading of that limit; see also W. Whitt (2007), What you
% should know about queueing models to set staffing requirements in service
% systems, Naval Research Logistics 54(5), 476-484.
%
% See also QSYS_MMK_QED, QSYS_MMK_QED_ALPHA, QSYS_MMK.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

options = struct('criterion', 'delay', 'exact', false, 'maxservers', 1e7);
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

if lambda <= 0
    line_error(mfilename, 'The arrival rate lambda must be positive.');
end
if mu <= 0
    line_error(mfilename, 'The service rate mu must be positive.');
end
a = lambda / mu;
criterion = lower(options.criterion);

switch criterion
    case 'delay'
        if ~(isnumeric(target) && isscalar(target)) || target <= 0 || target >= 1
            line_error(mfilename, 'For the delay criterion the target must be in (0,1).');
        end
        betaTarget = qsys_mmk_qed_invalpha(target);
    case 'meanwait'
        if ~(isnumeric(target) && isscalar(target)) || target <= 0
            line_error(mfilename, 'For the meanwait criterion the target must be positive.');
        end
        % E[W] = alpha(beta)/(s*mu-lambda) and s*mu-lambda = mu*beta*sqrt(s), so
        % with s ~ a + beta*sqrt(a) the target reads alpha(beta)/(mu beta sqrt(a)) = T.
        % The residual is written target - E[W] so that it increases in beta, as
        % the bisection below requires.
        betaTarget = qsys_mmk_qed_solve(@(b) target - qsys_mmk_qed_alpha(b)/(mu*b*sqrt(a)));
    case 'servicelevel'
        if ~isstruct(target) || ~isfield(target, 'deadline') || ~isfield(target, 'level')
            line_error(mfilename, ['For the servicelevel criterion the target must be a struct ' ...
                'with fields deadline and level.']);
        end
        if target.level <= 0 || target.level >= 1 || target.deadline <= 0
            line_error(mfilename, 'The service level must be in (0,1) and the deadline positive.');
        end
        % P(W > t) = alpha(beta) exp(-(s*mu-lambda) t), s*mu-lambda = mu beta sqrt(s).
        betaTarget = qsys_mmk_qed_solve(@(b) ...
            (1 - qsys_mmk_qed_alpha(b)*exp(-mu*b*sqrt(a + b*sqrt(a))*target.deadline)) - target.level);
    otherwise
        line_error(mfilename, sprintf('unknown criterion %s', criterion));
end

s = max(1, ceil(a + betaTarget*sqrt(a)));
if s*mu <= lambda
    s = floor(a) + 1;
end

exactUsed = false;
if options.exact
    [s, exactUsed] = qsys_mmk_qed_refine(lambda, mu, s, target, criterion, options.maxservers);
end

qed = qsys_mmk_qed(lambda, mu, s);
result.numServers = s;
result.beta = qed.beta;
result.betaTarget = betaTarget;
result.offeredLoad = a;
result.probDelay = qed.probDelay;
result.meanWait = qed.meanWait;
result.exactUsed = exactUsed;
if strcmp(criterion, 'servicelevel')
    result.serviceLevel = 1 - qed.probDelay*exp(-(s*mu - lambda)*target.deadline);
end
end

function beta = qsys_mmk_qed_invalpha(alphaTarget)
% Invert the strictly decreasing alpha(beta) by bisection.
beta = qsys_mmk_qed_solve(@(b) alphaTarget - qsys_mmk_qed_alpha(b));
end

function b = qsys_mmk_qed_solve(f)
% Bisection for a root of f on (0, hi], f being monotone increasing in beta over
% the range every caller uses. The bracket grows until the sign changes.
lo = 1e-9;
hi = 1;
flo = f(lo);
while f(hi) < 0
    hi = 2*hi;
    if hi > 1e6
        line_error(mfilename, 'no server slack meets the target; the target is unattainable');
    end
end
if flo > 0
    b = lo;
    return
end
for i = 1:200
    mid = (lo + hi)/2;
    if f(mid) < 0
        lo = mid;
    else
        hi = mid;
    end
end
b = (lo + hi)/2;
end

function [s, used] = qsys_mmk_qed_refine(lambda, mu, s0, target, criterion, maxServers)
% Walk s until the EXACT Erlang C measure meets the target. The exact measure is
% monotone in s, so a single direction of travel suffices.
used = true;
s = max(1, s0);
while ~qsys_mmk_qed_meets(lambda, mu, s, target, criterion)
    s = s + 1;
    if s > maxServers
        line_error(mfilename, 'the exact refinement passed maxServers without meeting the target');
    end
end
while s > 1 && qsys_mmk_qed_meets(lambda, mu, s-1, target, criterion)
    s = s - 1;
end
end

function ok = qsys_mmk_qed_meets(lambda, mu, s, target, criterion)
% The exact M/M/s measure against the target.
if s*mu <= lambda
    ok = false;
    return
end
C = qsys_mmk_qed_erlangc(s, lambda, mu);
Wq = C / (s*mu - lambda);
switch criterion
    case 'delay'
        ok = C <= target;
    case 'meanwait'
        ok = Wq <= target;
    case 'servicelevel'
        ok = (1 - C*exp(-(s*mu - lambda)*target.deadline)) >= target.level;
    otherwise
        ok = false;
end
end

function C = qsys_mmk_qed_erlangc(s, lambda, mu)
% Erlang C by the recursion B_j = a B_{j-1}/(j + a B_{j-1}) on the Erlang B
% blocking probability, which never forms a^j/j! and so never overflows. That
% matters here: this function is called at the s the staffing rule proposes,
% which is routinely in the thousands, where the factorial form is already Inf.
a = lambda / mu;
b = 1;
for j = 1:s
    b = a*b / (j + a*b);
end
rho = a / s;
if rho >= 1
    C = 1;
else
    C = b / (1 - rho*(1 - b));
end
end
