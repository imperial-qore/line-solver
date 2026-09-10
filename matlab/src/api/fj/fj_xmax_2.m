%{ @file fj_xmax_2.m
 %  @brief Expected maximum of 2 exponential random variables
 %
 %  @author LINE Development Team
%}

%{
 % @brief Expected maximum of 2 exponential random variables
 %
 % @details
 % Computes the expected value of the maximum of 2 independent
 % exponential random variables with rates lambda1 and lambda2.
 %
 % Y_2^max = 1/lambda1 + 1/lambda2 - 1/(lambda1 + lambda2)
 %
 % This is the exact closed-form solution for K=2 case.
 %
 % For identical rates (lambda1 = lambda2 = lambda):
 %   Y_2^max = 2/lambda - 1/(2*lambda) = 3/(2*lambda) = 1.5/lambda = H_2/lambda
 %
 % which matches the general formula X_K^max = H_K/mu for exponentials.
 %
 % @par Syntax:
 % @code
 % Xmax = fj_xmax_2(lambda1, lambda2)
 % Xmax = fj_xmax_2(lambda)  % Same rate for both
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lambda1<td>Rate of first exponential
 % <tr><td>lambda2<td>Rate of second exponential (optional, default=lambda1)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Xmax<td>Expected maximum Y_2^max
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Eq. (27) on page 17:17.
%}
function Xmax = fj_xmax_2(lambda1, lambda2)

if nargin < 2
    lambda2 = lambda1;  % Default to identical rates
end

if lambda1 <= 0 || lambda2 <= 0
    line_error(mfilename, 'Rates must be positive. Got lambda1=%.4f, lambda2=%.4f.', lambda1, lambda2);
end

% Y_2^max = 1/lambda1 + 1/lambda2 - 1/(lambda1 + lambda2)
Xmax = 1/lambda1 + 1/lambda2 - 1/(lambda1 + lambda2);

end
