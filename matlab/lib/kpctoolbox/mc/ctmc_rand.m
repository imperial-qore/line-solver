%{ @file ctmc_rand.m
 %  @brief Generates a random infinitesimal generator matrix
 %
 %  @author LINE Development Team
%}

%{
 % @brief Generates a random infinitesimal generator matrix
 %
 % @details
 % Creates a random infinitesimal generator of specified size.
 %
 % @par Syntax:
 % @code
 % Q = ctmc_rand(n)
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
 % <tr><td>Q<td>Random infinitesimal generator
 % </table>
%}
function Q=ctmc_rand(n)
    Q=ctmc_makeinfgen(rand(n));
end
