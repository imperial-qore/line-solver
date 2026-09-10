%{ @file sn_has_lept.m
 %  @brief Checks if the network uses LEPT scheduling
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network uses LEPT scheduling
 %
 % @details
 % Returns true if any station uses Longest Expected Processing Time scheduling.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_lept(sn)
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
 % <tr><td>bool<td>True if any station uses LEPT scheduling
 % </table>
%}
function bool = sn_has_lept(sn)

bool = any(sn.sched==SchedStrategy.LEPT);
end