%{ @file sn_has_sjf.m
 %  @brief Checks if the network uses SJF scheduling
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network uses SJF scheduling
 %
 % @details
 % Returns true if any station in the network uses Shortest Job First scheduling.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_sjf(sn)
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
 % <tr><td>bool<td>True if any station uses SJF scheduling
 % </table>
%}
function bool = sn_has_sjf(sn)

bool = any(sn.sched==SchedStrategy.SJF);
end