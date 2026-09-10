%{ @file fj_simplex.m
 %  @brief Two-phase dense simplex for small equality-form linear programs
 %
 %  @author LINE Development Team
%}

%{
 % @brief Two-phase dense simplex for small equality-form linear programs
 %
 % @details
 % Solves min c'z subject to A z = b and z >= 0 with a dense two-phase simplex
 % using Bland's rule, which terminates without cycling at the cost of a slower
 % pivot sequence. The programs raised by the fork-join capacity analyses have
 % at most a few hundred columns, so the dense tableau is the right trade and
 % keeps the four language ports free of any linear programming dependency.
 %
 % Phase one drives a full set of artificial variables out of the basis; if
 % their sum cannot be driven to zero the program is infeasible. Phase two then
 % optimises the true objective on the remaining basis.
 %
 % @par Syntax:
 % @code
 % [z, status] = fj_simplex(A, b, c)
 % [z, status, obj] = fj_simplex(A, b, c)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>A<td>Equality constraint matrix, m by n
 % <tr><td>b<td>Right hand side, m by 1
 % <tr><td>c<td>Cost row, 1 by n, minimised
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>z<td>Optimal point, 1 by n
 % <tr><td>status<td>0 optimal, 1 infeasible, 2 unbounded, 3 iteration limit
 % <tr><td>obj<td>Optimal objective value
 % </table>
%}
function [z, status, obj] = fj_simplex(A, b, c)

[m, n] = size(A);
b = b(:);
c = c(:)';

tol = 1e-10;
maxit = 200 * (m + n) + 1000;

% Normalise the right hand side to be non-negative
for i = 1:m
    if b(i) < 0
        A(i, :) = -A(i, :);
        b(i) = -b(i);
    end
end

% Phase one tableau: original columns, artificials, right hand side
T = [A, eye(m), b];
basis = n + (1:m);
% Reduced costs of minimising the sum of the artificials
obj1 = [-sum(A, 1), zeros(1, m), -sum(b)];

[T, basis, obj1, flag] = pivot_loop(T, basis, obj1, maxit, tol);
if flag == 2
    z = zeros(1, n); status = 2; obj = -Inf; return
end
if -obj1(end) > 1e-7 * max(1, norm(b, 1))
    z = zeros(1, n); status = 1; obj = NaN; return
end

% Drive any remaining artificial out of the basis, dropping redundant rows
row = 1;
while row <= size(T, 1)
    if basis(row) > n
        piv = find(abs(T(row, 1:n)) > tol, 1, 'first');
        if isempty(piv)
            T(row, :) = [];
            basis(row) = [];
            continue
        end
        T = do_pivot(T, row, piv);
        basis(row) = piv;
    end
    row = row + 1;
end

% Phase two on the original columns only
T = T(:, [1:n, n + m + 1]);
m2 = size(T, 1);
obj2 = [c, 0];
for i = 1:m2
    obj2 = obj2 - obj2(basis(i)) * T(i, :);
end

[T, basis, obj2, flag] = pivot_loop(T, basis, obj2, maxit, tol);
if flag == 2
    z = zeros(1, n); status = 2; obj = -Inf; return
end
if flag == 3
    z = zeros(1, n); status = 3; obj = NaN; return
end

z = zeros(1, n);
for i = 1:numel(basis)
    z(basis(i)) = T(i, end);
end
obj = -obj2(end);
status = 0;

end

function [T, basis, obj, flag] = pivot_loop(T, basis, obj, maxit, tol)
% Bland's rule: smallest index with a negative reduced cost enters
flag = 0;
ncol = size(T, 2) - 1;
for it = 1:maxit
    enter = 0;
    for j = 1:ncol
        if obj(j) < -tol
            enter = j;
            break
        end
    end
    if enter == 0
        return
    end
    leave = 0;
    best = Inf;
    for i = 1:size(T, 1)
        if T(i, enter) > tol
            ratio = T(i, end) / T(i, enter);
            if ratio < best - tol || (abs(ratio - best) <= tol && (leave == 0 || basis(i) < basis(leave)))
                best = ratio;
                leave = i;
            end
        end
    end
    if leave == 0
        flag = 2;
        return
    end
    T = do_pivot(T, leave, enter);
    obj = obj - obj(enter) * T(leave, :);
    basis(leave) = enter;
end
flag = 3;
end

function T = do_pivot(T, row, col)
T(row, :) = T(row, :) / T(row, col);
for i = 1:size(T, 1)
    if i ~= row && T(i, col) ~= 0
        T(i, :) = T(i, :) - T(i, col) * T(row, :);
    end
end
end
