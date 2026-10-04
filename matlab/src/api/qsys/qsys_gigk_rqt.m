function [W,rhohat,Sworst] = qsys_gigk_rqt(lambda,mu,Gamma_a,Gamma_s,k,alpha_a,alpha_s)
% [W,RHOHAT,SWORST]=QSYS_GIGK_RQT(LAMBDA,MU,GAMMA_A,GAMMA_S,K,ALPHA_A,ALPHA_S)
%
% Robust Queueing Theory (RQT) worst-case system time of a G/G/k FCFS queue.
% The arrival and service processes are not described by distributions but by
% the polyhedral uncertainty sets
%   U^a = { T : (sum_{i=k+1}^n T_i - (n-k)/lambda)/(n-k)^(1/alpha_a) >= -Gamma_a }
%   U^s = { X : (sum_{i=k}^n X_i - (n-k+1)/mu)/(n-k+1)^(1/alpha_s) <= Gamma_s }
% whose shape follows the (generalized) central limit theorem: alpha=2 is the
% finite-variance regime, alpha in (1,2) the heavy-tailed one. The performance
% analysis is then a worst-case optimization rather than an expectation.
%
% The returned W is the closed-form bound of Theorem 3 (Theorem 8 when the two
% tail coefficients differ, with alphabar = min(alpha_a,alpha_s)),
%
%   W <= (alphabar-1)/alphabar^(alphabar/(alphabar-1))
%        * lambda^(1/(alphabar-1)) (Gamma_a + Gamma_s/k^(1/alphabar))^(alphabar/(alphabar-1))
%          / (1-rho)^(1/(alphabar-1))  +  k/lambda,
%
% which for k=1 reduces to Theorem 2 and, at alphabar=2, to the Kingman-like
% form (lambda/4)(Gamma_a+Gamma_s)^2/(1-rho) + 1/lambda. SWORST returns instead
% the exact worst case over the uncertainty sets, eq. (45), i.e. the supremum
% over the integer x = nu-j+1 >= 1 of
%   x/mu + Gamma_s x^(1/alpha_s) - k(x-1)/lambda + Gamma_a (k(x-1))^(1/alpha_a),
% a one-dimensional problem. W >= SWORST by construction, and the two agree
% closely in heavy traffic. The arrival deviation ADDS to the worst case, since
% the adversary shortens the interarrival times: the sign printed in eq. (12) is
% easily misread as a subtraction of the whole arrival bracket, and reading it
% that way puts SWORST an order of magnitude below W.
%
% Note that W is a SYSTEM time (waiting plus service), not a waiting time, and
% that its additive term is k/lambda rather than the mean service time 1/mu.
%
% Inputs:
%   LAMBDA  - Arrival rate
%   MU      - Service rate of each server
%   GAMMA_A - Variability parameter of the arrival uncertainty set
%   GAMMA_S - Variability parameter of the service uncertainty set
%   K       - Number of servers (default 1)
%   ALPHA_A - Arrival tail coefficient in (1,2] (default 2)
%   ALPHA_S - Service tail coefficient in (1,2] (default 2)
%
% Returns:
%   W       - Closed-form bound on the system time (Theorem 3 / Theorem 8)
%   RHOHAT  - Modified utilization (so that M/M/1 formulas still hold)
%   SWORST  - Exact worst-case system time over the uncertainty sets, eq. (45)
%
% Reference: C. Bandi, D. Bertsimas, N. Youssef (2015). Robust Queueing Theory.
% Operations Research 63(3), 676-700.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5 || isempty(k), k = 1; end
if nargin < 6 || isempty(alpha_a), alpha_a = 2; end
if nargin < 7 || isempty(alpha_s), alpha_s = 2; end

if alpha_a <= 1 || alpha_a > 2 || alpha_s <= 1 || alpha_s > 2
    line_error(mfilename, 'RQT tail coefficients must lie in (1,2].');
end

rho = lambda / (k*mu);
if lambda <= 0
    W = 1/mu; rhohat = 0; Sworst = 1/mu;
    return
end
if rho >= 1
    W = Inf; rhohat = 1; Sworst = Inf;
    return
end

% Theorem 8 collapses to Theorem 3 when the two tails agree
ab = min(alpha_a, alpha_s);
beta = Gamma_a + Gamma_s / k^(1/ab);
if beta <= 0
    % a nonpositive effective variability leaves only the deterministic term
    W = k/lambda;
else
    W = (ab-1)/ab^(ab/(ab-1)) * lambda^(1/(ab-1)) * beta^(ab/(ab-1)) ...
        / (1-rho)^(1/(ab-1)) + k/lambda;
end
rhohat = W*lambda/(1+W*lambda); % so that M/M/1 formulas still hold

if nargout < 3
    return
end

% Exact worst case, eq. (45): sup over the integer x = nu-j+1 >= 1. The
% objective grows like x(1/mu - k/lambda) < 0, so the supremum is attained.
    function g = obj(x)
        y = max(x - 1, 0);
        g = x/mu + Gamma_s*x^(1/alpha_s) - k*y/lambda + Gamma_a*(k*y)^(1/alpha_a);
    end
% the continuous maximizer of the bounding problem, eq. (16), sizes the scan
if beta > 0
    xstar = (lambda*beta/(ab*(1-rho)))^(ab/(ab-1));
else
    xstar = 1;
end
xhi = max(4, ceil(4*xstar));
xs = unique(round(logspace(0, log10(xhi), 400)));
gs = arrayfun(@obj, xs);
[Sworst, imax] = max(gs);
% refine on the continuous relaxation, then round back onto the integer lattice
lo = xs(max(1,imax-1)); hi = xs(min(numel(xs),imax+1));
if hi > lo
    xc = fminbnd(@(x) -obj(x), lo, hi, optimset('TolX',1e-8));
    for x = [floor(xc), ceil(xc)]
        if x >= 1
            Sworst = max(Sworst, obj(x));
        end
    end
end
end
