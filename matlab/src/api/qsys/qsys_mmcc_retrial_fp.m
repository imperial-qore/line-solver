function [blocProb, r, niter] = qsys_mmcc_retrial_fp(lambda, mu, c, tol, maxiter)
% [blocProb, r, niter] = QSYS_MMCC_RETRIAL_FP(LAMBDA, MU, C, TOL, MAXITER)
%
% Fixed-point approximation for M/M/c/c retrial queues.
%
% Customers arrive at rate LAMBDA to a system with C servers, each with
% service rate MU. Blocked customers join an orbit and retry. Under the
% assumption that the retrial rate is small relative to the service rate,
% the total arrival flow (fresh + retrial) is approximated by a Poisson
% process with rate LAMBDA + r, where r satisfies the fixed-point equation:
%
%   r = (lambda + r) * B(lambda/mu + r/mu, c)
%
% and B(a, c) is the Erlang-B blocking probability for offered load a and
% c servers.
%
% INPUT:
%   lambda  : arrival rate
%   mu      : service rate per server
%   c       : number of servers (= capacity, no waiting room)
%   tol     : convergence tolerance (default: 1e-10)
%   maxiter : maximum iterations (default: 10000)
%
% OUTPUT:
%   blocProb : blocking probability (fraction of arrivals lost or retried)
%   r        : additional arrival rate due to retrials
%   niter    : number of iterations to converge
%
% REFERENCE:
%   Cohen (1957), fixed-point approximation for M/M/c/c retrial queues.
%   Phung-Duc, "Retrial Queueing Models: A Survey on Theory and
%   Applications", 2019, Eq. (1).

if nargin < 4 || isempty(tol)
    tol = 1e-10;
end
if nargin < 5 || isempty(maxiter)
    maxiter = 10000;
end

r = 0;
niter = 0;
for iter = 1:maxiter
    niter = iter;
    a = (lambda + r) / mu;  % offered load
    b = ErlangB(a, c);
    r_new = (lambda + r) * b;
    if abs(r_new - r) < tol
        r = r_new;
        break;
    end
    r = r_new;
end

blocProb = ErlangB((lambda + r) / mu, c);
end

function B = ErlangB(a, c)
% Erlang-B formula using the recursive method (numerically stable).
% B(a, c) = blocking probability for offered load a and c servers.
B = 1.0;
for i = 1:c
    B = a * B / (i + a * B);
end
end
