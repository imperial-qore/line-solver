%{ @file fj_lst_max_het.m
 %  @brief Laplace-Stieltjes transform of the maximum of heterogeneous exponentials
 %
 %  @author LINE Development Team
%}

%{
 % @brief Laplace-Stieltjes transform of the maximum of heterogeneous exponentials
 %
 % @details
 % Evaluates L*_K(s) = E[exp(-s Y)] for Y = max(X_1,...,X_K) with independent
 % X_i ~ Exp(lambda_i), using the Harrison and Zertal recurrence over the
 % sub-collections obtained by deleting one rate at a time:
 %
 %   ( s + sum_{j=1..m} lambda_j ) L*_m(lambda, s)
 %       = sum_{j=1..m} lambda_j * L*_{m-1}(lambda \ j, s),   1 <= m <= K,
 %
 % anchored at L*_0 = 1 because the maximum of an empty collection is zero.
 % The recurrence is evaluated bottom-up over the 2^K sub-collections, each
 % identified by a bit mask, so every value is computed once.
 %
 % @par Syntax:
 % @code
 % L = fj_lst_max_het(lambda, s)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lambda<td>Vector of K positive exponential rates
 % <tr><td>s<td>Transform argument, scalar or vector, s >= 0
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>L<td>Transform value, same shape as s
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Eq. (29).
 %
 % Original: P. G. Harrison, S. Zertal, "Queueing Models of RAID Systems with
 % Maxima of Waiting Times", Performance Evaluation 64(7-8), 2007.
%}
function L = fj_lst_max_het(lambda, s)

lambda = lambda(:)';
K = numel(lambda);

if K < 1
    line_error(mfilename, 'lambda must contain at least one rate.');
end
if any(lambda <= 0)
    line_error(mfilename, 'All exponential rates must be positive.');
end
if any(s < 0)
    line_error(mfilename, 'The transform argument s must be non-negative.');
end
if K > 22
    line_error(mfilename, 'The recurrence enumerates 2^K sub-collections; K=%d is too large.', K);
end

L = zeros(size(s));
for idx = 1:numel(s)
    L(idx) = lst_at(lambda, K, s(idx));
end

end

function val = lst_at(lambda, K, s)
% Bottom-up sweep over the 2^K masks; mask 0 is the empty collection
nmask = 2^K;
tab = zeros(1, nmask);
tab(1) = 1;
for mask = 1:(nmask - 1)
    bits = bitget(mask, 1:K);
    members = find(bits == 1);
    tot = sum(lambda(members));
    acc = 0;
    for j = members
        acc = acc + lambda(j) * tab(bitxor(mask, bitshift(1, j - 1)) + 1);
    end
    tab(mask + 1) = acc / (s + tot);
end
val = tab(nmask);
end
