%{ @file fj_order_stat.m
 %  @brief CDF and expected value of k-th order statistic
 %
 %  @author LINE Development Team
%}

%{
 % @brief CDF and expected value of k-th order statistic
 %
 % @details
 % Computes the CDF and expected value of the k-th order statistic
 % (k-th smallest) of K i.i.d. random variables.
 %
 % For maximum (k=K):
 %   F_{Y_K}(y) = [F_X(y)]^K
 %
 % For k-th order statistic (k-th smallest):
 %   F_{Y_k}(y) = sum_{j=k}^{K} C(K,j) * F_X(y)^j * (1-F_X(y))^{K-j}
 %
 % Expected value of maximum:
 %   E[Y_K] = integral_0^inf [1 - F_X(y)^K] dy
 %
 % @par Syntax:
 % @code
 % F_Yk = fj_order_stat(y, k, K, F_X)
 % [F_Yk, E_Yk] = fj_order_stat(y, k, K, F_X)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>y<td>Value(s) at which to evaluate CDF
 % <tr><td>k<td>Order of the statistic (1 = minimum, K = maximum)
 % <tr><td>K<td>Total number of random variables
 % <tr><td>F_X<td>CDF function handle: F_X(y) returns CDF value
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>F_Yk<td>CDF of k-th order statistic at y
 % <tr><td>E_Yk<td>Expected value of k-th order statistic (if requested)
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Eq. (18) and (19) on page 17:15.
%}
function [F_Yk, E_Yk] = fj_order_stat(y, k, K, F_X)

if k < 1 || k > K
    line_error(mfilename, 'k must satisfy 1 <= k <= K. Got k=%d, K=%d.', k, K);
end

% Evaluate the base CDF at y
F_y = F_X(y);

if k == K
    % Maximum (K-th order statistic): F_{Y_K}(y) = [F_X(y)]^K
    F_Yk = F_y.^K;
else
    % General k-th order statistic:
    % F_{Y_k}(y) = sum_{j=k}^{K} C(K,j) * F_X(y)^j * (1-F_X(y))^{K-j}
    F_Yk = zeros(size(y));
    for j = k:K
        F_Yk = F_Yk + nchoosek(K, j) .* (F_y.^j) .* ((1 - F_y).^(K - j));
    end
end

% Compute expected value if requested
if nargout > 1
    if k == K
        % Expected value of maximum: E[Y_K] = integral_0^inf [1 - F_X(y)^K] dy
        integrand = @(t) 1 - F_X(t).^K;
    elseif k == 1
        % Expected value of minimum: E[Y_1] = integral_0^inf [1 - F_X(y)]^K dy
        integrand = @(t) (1 - F_X(t)).^K;
    else
        % General case: use the formula involving order statistic pdf
        % E[Y_k] = integral_0^inf y * f_{Y_k}(y) dy
        % where f_{Y_k}(y) = K * C(K-1,k-1) * f_X(y) * F_X(y)^{k-1} * (1-F_X(y))^{K-k}

        % Use numerical differentiation for pdf (approximate)
        eps_val = 1e-8;
        f_X = @(t) (F_X(t + eps_val) - F_X(t - eps_val)) / (2 * eps_val);

        coeff = K * nchoosek(K - 1, k - 1);
        integrand = @(t) t .* coeff .* f_X(t) .* (F_X(t).^(k - 1)) .* ((1 - F_X(t)).^(K - k));
    end

    % Numerical integration
    % Find reasonable upper limit (where CDF is very close to 1)
    upper_limit = find_upper_limit(F_X, 0.999999);
    E_Yk = integral(integrand, 0, upper_limit, 'RelTol', 1e-8, 'AbsTol', 1e-10);
end

end

function u = find_upper_limit(F_X, target)
% Find value u such that F_X(u) >= target
u = 1;
max_iter = 100;
for iter = 1:max_iter
    if F_X(u) >= target
        break;
    end
    u = u * 2;
end
end
