%{ @file sn_has_siro.m
 %  @brief Checks if the network uses SIRO scheduling
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network uses SIRO scheduling
 %
 % @details
 % Returns true if any station uses Service In Random Order scheduling.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_siro(sn)
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
 % <tr><td>bool<td>True if any station uses SIRO scheduling
 % </table>
%}
function bool = sn_has_siro(sn)

bool = any(sn.sched==SchedStrategy.SIRO);
end
