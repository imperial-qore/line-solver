%{ @file dtmc_rand.m
 %  @brief Generates a random stochastic transition matrix
 %
 %  @author LINE Development Team
%}

%{
 % @brief Generates a random stochastic transition matrix
 %
 % @details
 % Creates a random stochastic matrix of specified size.
 %
 % @par Syntax:
 % @code
 % P = dtmc_rand(n)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>n<td>Size of the matrix (n x n)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>P<td>Random stochastic transition matrix
 % </table>
%}
function P=dtmc_rand(n)
P=ctmc_randomization(ctmc_rand(n));
end