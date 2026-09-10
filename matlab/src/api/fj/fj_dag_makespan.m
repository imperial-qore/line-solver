%{ @file fj_dag_makespan.m
 %  @brief Makespan of a task system with precedence constraints
 %
 %  @author LINE Development Team
%}

%{
 % @brief Makespan of a task system with precedence constraints
 %
 % @details
 % A task system is a set of n tasks whose precedence relations form a directed
 % acyclic graph. Because the precedence relation is acyclic, so is the chain
 % whose state is the SET of completed tasks, and the makespan can be swept
 % level by level instead of solved as a linear system.
 %
 % In a state whose completed set is S, the tasks eligible to run are those all
 % of whose predecessors lie in S; a task i among the k of them completes at
 % rate rate(i,k), so the state is held for M(S) = 1/T(S) with T(S) the sum of
 % those rates, and moves to S union {i} with probability b(S,i) = rate(i,k)/T(S).
 % Making the rate depend on how many tasks run concurrently is what couples the
 % task system to the queueing network underneath it: two tasks sharing a
 % processor each run slower than either would alone.
 %
 % The weighted delay to reach a state and the probability of reaching it obey
 %
 %   p(R) = sum_{S -> R} p(S)*b(S,R),
 %   D(R) = M(R)*p(R) + sum_{S -> R} b(S,R)*D(S),
 %
 % started at p(empty) = 1, and the makespan is D at the fully completed state.
 % The initiation, completion and execution times of each task follow from the
 % same sweep, D being accumulated over the transitions that start and that
 % finish the task.
 %
 % @par Syntax:
 % @code
 % C = fj_dag_makespan(pred, rate)
 % [C, I, Cend, E] = fj_dag_makespan(pred, rate)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>pred<td>n by n logical matrix, pred(i,j) true when task i precedes task j
 % <tr><td>rate<td>Completion rates: an n-vector of constant rates, or an n by n matrix whose entry (i,k) is the rate of task i while k tasks run concurrently
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>C<td>Mean makespan, the time until every task has completed
 % <tr><td>I<td>Vector of mean initiation times, one per task
 % <tr><td>Cend<td>Vector of mean completion times, one per task
 % <tr><td>E<td>Vector of mean execution times, Cend - I
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Eq. (66) and
 % Section 7.4 on page 17:49.
 %
 % Original: A. Thomasian, P. F. Bay, "Analytic Queueing Network Models for
 % Parallel Processing of Task Systems", IEEE Trans. Computers C-35(12), 1986.
%}
function [C, I, Cend, E] = fj_dag_makespan(pred, rate)

pred = logical(pred);
n = size(pred, 1);

if size(pred, 2) ~= n
    line_error(mfilename, 'pred must be square. Got %dx%d.', n, size(pred, 2));
end
if n < 1
    line_error(mfilename, 'At least one task is required.');
end
if n > 20
    line_error(mfilename, 'The completed-set sweep enumerates 2^n states; n=%d is too large.', n);
end

% Rates as an n by n table indexed by (task, concurrency level)
if isvector(rate)
    rate = rate(:);
    if numel(rate) ~= n
        line_error(mfilename, 'The rate vector must have one entry per task. Got %d for %d tasks.', numel(rate), n);
    end
    rate = repmat(rate, 1, n);
elseif ~isequal(size(rate), [n n])
    line_error(mfilename, 'The rate matrix must be %dx%d. Got %dx%d.', n, n, size(rate, 1), size(rate, 2));
end
if any(rate(:) <= 0)
    line_error(mfilename, 'All completion rates must be positive.');
end

% Reject a cycle in the precedence relation by attempting a topological order
indeg = sum(pred, 1);
remaining = n;
seen = false(1, n);
for pass = 1:n
    idx = find(~seen & indeg == 0, 1, 'first');
    if isempty(idx)
        break
    end
    seen(idx) = true;
    indeg = indeg - double(pred(idx, :));
    indeg(idx) = Inf;
    remaining = remaining - 1;
end
if remaining > 0
    line_error(mfilename, 'The precedence relation contains a cycle.');
end

% Predecessor masks: bit j-1 of predmask(j) marks a task that must precede j
predmask = zeros(1, n);
for j = 1:n
    for i = 1:n
        if pred(i, j)
            predmask(j) = bitset(predmask(j), i);
        end
    end
end

nmask = 2^n;
p = zeros(1, nmask);
D = zeros(1, nmask);
p(1) = 1;

% Eligibility masks of every completed set, precomputed once
eligmask = zeros(1, nmask);
closed = false(1, nmask);
for mask = 0:(nmask - 1)
    ok = true;
    em = 0;
    for i = 1:n
        if bitget(mask, i)
            % A completed task must have all of its predecessors completed
            if bitand(predmask(i), mask) ~= predmask(i)
                ok = false;
                break
            end
        elseif bitand(predmask(i), mask) == predmask(i)
            em = bitset(em, i);
        end
    end
    closed(mask + 1) = ok;
    if ok
        eligmask(mask + 1) = em;
    end
end

% Per-task accumulators over the transitions that start and that finish a task
Istart = zeros(1, n);
Cfin = zeros(1, n);

for mask = 0:(nmask - 1)
    if ~closed(mask + 1)
        continue
    end
    elig = find(bitget(eligmask(mask + 1), 1:n) == 1);
    k = numel(elig);
    if k == 0
        continue
    end
    Ttot = sum(rate(elig, k));
    % Holding time of this state, weighted by the probability of reaching it
    D(mask + 1) = D(mask + 1) + p(mask + 1) / Ttot;
    for i = elig
        b = rate(i, k) / Ttot;
        nxt = bitset(mask, i);
        contrib = b * D(mask + 1);
        p(nxt + 1) = p(nxt + 1) + p(mask + 1) * b;
        D(nxt + 1) = D(nxt + 1) + contrib;
        % Task i completes on this transition
        Cfin(i) = Cfin(i) + contrib;
        % Tasks that first become eligible on this transition start on it
        fresh = bitand(eligmask(nxt + 1), bitcmp_n(eligmask(mask + 1), n));
        for j = find(bitget(fresh, 1:n) == 1)
            Istart(j) = Istart(j) + contrib;
        end
    end
end

C = D(nmask);
Cend = Cfin;
I = Istart;
E = Cend - I;

end

function y = bitcmp_n(x, n)
% Complement of x within n bits, avoiding any dependence on the integer class
y = bitxor(x, 2^n - 1);
end
