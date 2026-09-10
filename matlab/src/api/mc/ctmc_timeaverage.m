function [piTimeAvg,piExit,kmax]=ctmc_timeaverage(pi0,Q,t,tol,maxiter)
% [PITIMEAVG,PIEXIT,KMAX]=CTMC_TIMEAVERAGE(PI0,Q,T,TOL,MAXITER)
%
% Time-averaged transient distribution of a CTMC with generator Q over [0,T],
% starting from the (arbitrary) initial distribution PI0, via Jensen's
% uniformization. Companion of CTMC_UNIFORMIZATION, which returns only the
% endpoint PI0*exp(Q*T); this function additionally returns the time average
%
%   PITIMEAVG = PI0 * (1/T) * \int_0^T exp(Q*tau) d(tau)
%
% as well as the endpoint PIEXIT = PI0*exp(Q*T) (computed from the same series).
% Both are obtained without forming any dense matrix exponential.
%
% Uniformization: with q = 1.1*max|diag(Q)| and P = I + Q/q (row-stochastic),
%   PI0*exp(Q*T)            = sum_j w_j(qT) * (PI0*P^j)
%   PI0*\int_0^T exp(Q*tau) = (1/q) * sum_j (1 - W_j(qT)) * (PI0*P^j)
% where w_j and W_j are the Poisson(qT) PMF and CDF. The time average divides
% the integral by T (equivalently the integral sum by q*T).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<4 % tol
    tol = 1e-12;
end
if nargin<5 % maxiter
    maxiter = -1;
end

q  = 1.1*max(abs(diag(Q)));
% Split the horizon into equal segments with q*tSeg below the underflow
% bound (exp(-745)==0): the integral over [0,t] is the sum of segment
% integrals, each started from the previous segment's endpoint
MAXQT = 500;
if q*t > MAXQT
    nSeg = ceil(q*t/MAXQT);
    tSeg = t/nSeg;
    piCur = pi0;
    integral = zeros(1,size(Q,1));
    for seg=1:nSeg
        [avgSeg,piCur,kmax] = ctmc_timeaverage(piCur,Q,tSeg,tol,maxiter);
        integral = integral + tSeg*avgSeg(:)';
    end
    piTimeAvg = integral/t;
    piExit = piCur;
    return
end
if maxiter<=0
    % The Poisson(q*t) mass concentrates around q*t with spread O(sqrt(q*t));
    % a fixed cap silently truncates the series for large horizons
    maxiter = max(100, ceil(q*t+10*sqrt(q*t)+20));
end
Qs = speye(size(Q)) + sparse(Q)/q;
qt = q*t;

% Number of Poisson terms needed (right-tail below tol), as in ctmc_uniformization
k = 0; s = 1; r = 1; iter = 0; kmax = 1;
while iter < maxiter
    iter = iter + 1;
    k = k + 1;
    r = r*qt/k;
    s = s + r;
    if (1 - exp(-qt)*s) <= tol
        kmax = k;
        break;
    end
    % Best-effort truncation depth if the loop exhausts maxiter
    kmax = k;
end

% Accumulate endpoint and integral over the shared PI0*P^j sequence
w  = exp(-qt);          % Poisson PMF  w_0
W  = w;                 % Poisson CDF  W_0
P  = pi0;               % PI0*P^0
piExit   = w * P;
piIntSum = max(1 - W, 0) * P;
for j = 1:kmax
    P = P * Qs;         % PI0*P^j
    w = w * qt / j;     % w_j
    W = W + w;          % W_j
    piExit   = piExit   + w * P;
    piIntSum = piIntSum + max(1 - W, 0) * P;
end

piTimeAvg = piIntSum / qt;   % (1/q)*sum / t  =  sum/(q*t)
end
