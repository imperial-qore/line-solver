%{ @file sn_has_fcfs.m
 %  @brief Checks if the network has FCFS scheduling
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network has FCFS scheduling
 %
 % @details
 % Returns true if any station in the network uses First-Come First-Served scheduling.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_fcfs(sn)
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
 % <tr><td>bool<td>True if any station uses FCFS scheduling
 % </table>
%}
function bool = sn_has_fcfs(sn)

bool = any(sn.sched==SchedStrategy.FCFS);
end