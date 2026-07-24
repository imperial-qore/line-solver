%{ @file sn_has_inf.m
 %  @brief Checks if the network has infinite server stations
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network has infinite server stations
 %
 % @details
 % Returns true if any station in the network has an infinite number of servers.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_inf(sn)
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
 % <tr><td>bool<td>True if any station has infinite servers
 % </table>
%}
function bool = sn_has_inf(sn)

bool = any(sn.sched==SchedStrategy.INF);
end