%{ @file sn_has_fork_join.m
 %  @brief Checks if the network contains fork-join nodes
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network contains fork-join nodes
 %
 % @details
 % Returns true if the network contains any fork or join nodes.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_fork_join(sn)
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
 % <tr><td>bool<td>True if the network contains fork or join nodes
 % </table>
%}
function bool = sn_has_fork_join(sn)
bool = any(sn.fj(:) > 0);
end