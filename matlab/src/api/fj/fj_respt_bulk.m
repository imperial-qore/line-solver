%{ @file fj_respt_bulk.m
 %  @brief Centralized splitting analysed as an M[K]/M/c bulk arrival system
 %
 %  @author LINE Development Team
%}

%{
 % @brief Centralized splitting analysed as an M[K]/M/c bulk arrival system
 %
 % @details
 % Under centralized splitting a request forks into K tasks that are held in a
 % single central queue and served by c identical servers, so the same server
 % may serve several tasks of the same request. That is an M[K]/M/c queue with
 % fixed batch size K. Its stationary distribution has no product form, and is
 % obtained here by solving the balance equations of the level chain truncated
 % at nmax, which is exact up to the tail mass that truncation discards.
 %
 % Two response times are reported. The task response time follows from
 % Little's law on the mean number in system. The request response time is the
 % completion of the LAST of the K tasks of a tagged request: by PASTA the
 % batch finds n tasks in system, its last task is the (n+K)-th in line, and
 % under first come first served with c exponential servers it starts service
 % after max(0, n+K-c) departures, each an exponential of rate c*mu, whence
 %
 %   E[R_request] = sum_n p_n * [ max(0, n+K-c)/(c*mu) + 1/mu ].
 %
 % This is the centralized counterpart of the distributed splitting fork-join
 % system, and lower bounds it because no task is bound to a particular server.
 %
 % @par Syntax:
 % @code
 % [Rreq, Rtask] = fj_respt_bulk(K, lambda, mu, c)
 % [Rreq, Rtask, Q, p] = fj_respt_bulk(K, lambda, mu, c, nmax)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Batch size, that is the number of tasks per request
 % <tr><td>lambda<td>Arrival rate of requests (batches)
 % <tr><td>mu<td>Task service rate at each server
 % <tr><td>c<td>Number of servers
 % <tr><td>nmax<td>Truncation level of the task-count chain (optional)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Rreq<td>Mean request response time, that is the last of the K tasks
 % <tr><td>Rtask<td>Mean task response time
 % <tr><td>Q<td>Mean number of tasks in system
 % <tr><td>p<td>Stationary distribution of the number of tasks in system
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Section 6.1
 % on page 17:33.
 %
 % Original: R. Nelson, D. Towsley, A. N. Tantawi, "Performance Analysis of
 % Parallel Processing Systems", IEEE Trans. Software Eng. 14(4), 1988.
%}
function [Rreq, Rtask, Q, p] = fj_respt_bulk(K, lambda, mu, c, nmax)

if K < 1 || K ~= round(K)
    line_error(mfilename, 'K must be a positive integer. Got K=%g.', K);
end
if c < 1 || c ~= round(c)
    line_error(mfilename, 'c must be a positive integer. Got c=%g.', c);
end
if lambda <= 0 || mu <= 0
    line_error(mfilename, 'lambda and mu must be positive. Got lambda=%g, mu=%g.', lambda, mu);
end

rho = lambda * K / (c * mu);
if rho >= 1
    line_error(mfilename, 'System is unstable: rho = lambda*K/(c*mu) = %.4f >= 1.', rho);
end

if nargin < 5 || isempty(nmax)
    % Deep enough that the geometric tail of the truncated chain is negligible
    nmax = max(200, ceil(K + 40 * c / (1 - rho)));
end

% Generator of the task-count chain: batch arrivals of K, service rate min(n,c)*mu
n = 0:nmax;
ns = numel(n);
Qgen = zeros(ns, ns);
for i = 1:ns
    ni = n(i);
    if ni > 0
        srv = min(ni, c) * mu;
        Qgen(i, i - 1) = Qgen(i, i - 1) + srv;
        Qgen(i, i) = Qgen(i, i) - srv;
    end
    j = ni + K;
    if j <= nmax
        Qgen(i, j + 1) = Qgen(i, j + 1) + lambda;
        Qgen(i, i) = Qgen(i, i) - lambda;
    end
end

% Stationary distribution: p*Qgen = 0 with the normalization replacing a column
A = [Qgen'; ones(1, ns)];
b = [zeros(ns, 1); 1];
p = (A \ b)';
p = max(p, 0);
p = p / sum(p);

Q = sum(n .* p);
Rtask = Q / (lambda * K);

% Last of the K tasks of a tagged batch: it is the (n+K)-th in line on arrival
wait = max(0, n + K - c) / (c * mu);
Rreq = sum(p .* (wait + 1 / mu));

end
