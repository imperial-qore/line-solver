%{ @file ctmc_makeinfgen.m
 %  @brief Normalizes a matrix to be a valid infinitesimal generator
 %
 %  @author LINE Development Team
%}

%{
 % @brief Normalizes a matrix to be a valid infinitesimal generator
 %
 % @details
 % Sets diagonal elements such that row sums are zero.
 %
 % @par Syntax:
 % @code
 % Q = ctmc_makeinfgen(Q)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Q<td>Input matrix
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Q<td>Valid infinitesimal generator matrix where rows sum to zero
 % </table>
%}
function Q=ctmc_makeinfgen(Q)
A=Q-diag(diag(Q)); Q=A-diag(sum(A,2));
end
