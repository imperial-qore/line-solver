%{ @file cache_prob_fpi.m
 %  @brief Computes cache hit probabilities using fixed-point iteration
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes cache hit probabilities using FPI method
 %
 % @details
 % This function computes cache hit probability distribution using the
 % fixed-point iteration (FPI) method.
 %
 % @par Syntax:
 % @code
 % prob = cache_prob_fpi(gamma, m)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item popularity probabilities
 % <tr><td>m<td>Cache capacity vector
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>prob<td>Cache hit probability distribution
 % </table>
%}
function prob=cache_prob_fpi(gamma,m)
% FPI method
[n,h]=size(gamma);
xi = cache_xi_fp(gamma,m);
prob = zeros(n,h);
for i=1:n
    prob(i,1) = 1/(1+gamma(i,:)*xi(:));
    prob(i,2:(1+h)) = gamma(i,:)*xi(:) ./ (1+gamma(i,:)*xi(:));
end
end
