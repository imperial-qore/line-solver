%{ @file fj_tsm_capacity.m
 %  @brief Saturation throughput of the team service model
 %
 %  @author LINE Development Team
%}

%{
 % @brief Saturation throughput of the team service model
 %
 % @details
 % In the team service model a class-k job seizes r(k) of the s servers at
 % once, holds them for a mean x(k), and releases them all together. The
 % apparent saturation rate is the one that drives the server utilization to
 % one,
 %
 %   Lambda_max = s / sum_k f(k)*r(k)*x(k),
 %
 % but it is attainable only when the scheduler can pack jobs into execution
 % states that leave no server idle. The attainable capacity is the largest
 % arrival rate for which some mixture p over the feasible execution states
 % balances every class,
 %
 %   maximise Lambda subject to  sum_j p_j * n(j,k) / x(k) = Lambda * f(k),
 %                               sum_j p_j = 1,  p >= 0,
 %
 % over the execution states j, each a multiset of jobs with total server
 % demand at most s. n(j,k) counts the class-k jobs in state j. The linear
 % program is solved by a dense simplex on its standard form; the optimum
 % equals Lambda_max exactly when every state carrying positive probability is
 % full capacity.
 %
 % For the two-server two-class special case with r = (1,2), the strict first
 % come first served discipline cannot pack at all and reaches only
 %
 %   lambda_FCFS = 2*mu1*mu2 / (f1^2*mu2 + 2*f2^2*mu1 + 2*f1*f2*(mu1+mu2)),
 %
 % which is returned whenever the arguments describe that case.
 %
 % @par Syntax:
 % @code
 % Lmax = fj_tsm_capacity(s, f, r, x)
 % [Lmax, Llp, Lfcfs, p] = fj_tsm_capacity(s, f, r, x)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>s<td>Number of servers
 % <tr><td>f<td>Vector of class frequencies in the arrival stream, summing to one
 % <tr><td>r<td>Vector of per-class server requirements (positive integers, r(k) <= s)
 % <tr><td>x<td>Vector of per-class mean service times (positive)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Lmax<td>Apparent saturation rate that drives the utilization to one
 % <tr><td>Llp<td>Attainable capacity from the linear program, Llp <= Lmax
 % <tr><td>Lfcfs<td>Strict first come first served capacity, NaN outside the two-server two-class case
 % <tr><td>p<td>Execution-state probabilities attaining Llp, as a struct with fields states and prob
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Section 6.3,
 % Eqs. (57)-(58) and the linear programming formulation on page 17:36.
 %
 % Original: A. Thomasian, "Performance Evaluation of Centralized Databases
 % with Static Locking", IEEE Trans. Software Eng. SE-11(4), 1985;
 % K. Omahen, L. Schrage, "A Queueing Analysis of a Multiprocessor System with
 % Shared Memory", Symp. Computer-Communication Networks, 1972.
%}
function [Lmax, Llp, Lfcfs, p] = fj_tsm_capacity(s, f, r, x)

f = f(:)';
r = r(:)';
x = x(:)';
K = numel(f);

if s < 1 || s ~= round(s)
    line_error(mfilename, 's must be a positive integer. Got s=%g.', s);
end
if numel(r) ~= K || numel(x) ~= K
    line_error(mfilename, 'f, r and x must have the same length. Got %d, %d, %d.', K, numel(r), numel(x));
end
if any(f < 0) || abs(sum(f) - 1) > 1e-9
    line_error(mfilename, 'The class frequencies must be non-negative and sum to one. Got %g.', sum(f));
end
if any(r < 1) || any(r ~= round(r)) || any(r > s)
    line_error(mfilename, 'The server requirements must be integers in 1..%d.', s);
end
if any(x <= 0)
    line_error(mfilename, 'The mean service times must be positive.');
end

% Apparent saturation rate, Eq. (57)
Lmax = s / sum(f .* r .* x);

% Feasible execution states: multisets of jobs whose server demand fits in s
states = enumerate_states(r, s, K);
nstates = size(states, 1);

% Linear program: variables [p_1..p_nstates, Lambda], maximise Lambda
%   sum_j p_j * n(j,k)/x(k) - Lambda*f(k) = 0   for every class k with f(k) > 0
%   sum_j p_j = 1
active = find(f > 0);
Aeq = zeros(numel(active) + 1, nstates + 1);
beq = zeros(numel(active) + 1, 1);
for idx = 1:numel(active)
    k = active(idx);
    Aeq(idx, 1:nstates) = states(:, k)' / x(k);
    Aeq(idx, nstates + 1) = -f(k);
end
Aeq(end, 1:nstates) = 1;
beq(end) = 1;

cost = zeros(1, nstates + 1);
cost(nstates + 1) = -1;   % simplex minimises, so minimise -Lambda

[sol, status] = fj_simplex(Aeq, beq, cost);
if status ~= 0
    line_error(mfilename, 'The capacity linear program did not solve (status %d).', status);
end
Llp = sol(nstates + 1);
p = struct('states', states, 'prob', sol(1:nstates)');

% Strict first come first served capacity of the two-server two-class case
Lfcfs = NaN;
if s == 2 && K == 2 && all(sort(r) == [1 2])
    if r(1) == 1
        f1 = f(1); f2 = f(2); mu1 = 1 / x(1); mu2 = 1 / x(2);
    else
        f1 = f(2); f2 = f(1); mu1 = 1 / x(2); mu2 = 1 / x(1);
    end
    Lfcfs = 2 * mu1 * mu2 / (f1^2 * mu2 + 2 * f2^2 * mu1 + 2 * f1 * f2 * (mu1 + mu2));
end

end

function states = enumerate_states(r, s, K)
% Every multiset of jobs whose total server demand is at most s, excluding the
% empty state which can carry no completions
states = zeros(0, K);
stack = zeros(1, K);
states = recurse(states, stack, 1, s, r, K);
end

function states = recurse(states, stack, k, left, r, K)
if k > K
    if any(stack > 0)
        states = [states; stack];
    end
    return
end
nmax = floor(left / r(k));
for n = 0:nmax
    stack(k) = n;
    states = recurse(states, stack, k + 1, left - n * r(k), r, K);
end
end
