%{ @file sn_has_multi_chain.m
 %  @brief Checks if the network has multiple chains
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network has multiple chains
 %
 % @details
 % Returns true if the network has more than one chain.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_multi_chain(sn)
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
 % <tr><td>bool<td>True if number of chains is greater than one
 % </table>
%}
function bool = sn_has_multi_chain(sn)

bool = sn.nchains > 1;
end