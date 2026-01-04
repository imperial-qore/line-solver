%{ @file sn_has_sept.m
 %  @brief Checks if the network uses SEPT scheduling
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network uses SEPT scheduling
 %
 % @details
 % Returns true if any station uses Shortest Expected Processing Time scheduling.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_sept(sn)
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
 % <tr><td>bool<td>True if any station uses SEPT scheduling
 % </table>
%}
function bool = sn_has_sept(sn)

bool = any(sn.sched==SchedStrategy.SEPT);
end