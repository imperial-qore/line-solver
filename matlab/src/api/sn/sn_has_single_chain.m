%{ @file sn_has_single_chain.m
 %  @brief Checks if the network has a single chain
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network has a single chain
 %
 % @details
 % Returns true if the network has exactly one chain.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_single_chain(sn)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>bool<td>True if number of chains equals one
 % </table>
%}
function bool = sn_has_single_chain(sn)

bool = sn.nchains == 1;
end