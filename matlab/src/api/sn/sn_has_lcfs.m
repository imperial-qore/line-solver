%{ @file sn_has_lcfs.m
 %  @brief Checks if the network uses LCFS scheduling
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network uses LCFS scheduling
 %
 % @details
 % Returns true if any station in the network uses Last-Come First-Served scheduling.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_lcfs(sn)
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
 % <tr><td>bool<td>True if any station uses LCFS scheduling
 % </table>
%}
function bool = sn_has_lcfs(sn)

bool = any(sn.sched==SchedStrategy.LCFS);
end