%{ @file sn_has_lcfs_pi.m
 %  @brief Checks if the network uses LCFS with preemptive-identical scheduling
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network uses LCFS with preemptive-identical scheduling
 %
 % @details
 % Returns true if any station uses LCFS preemptive-identical scheduling.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_lcfs_pi(sn)
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
 % <tr><td>bool<td>True if any station uses LCFS-PI scheduling
 % </table>
%}
function bool = sn_has_lcfs_pi(sn)

bool = any(sn.sched==SchedStrategy.LCFSPI);
end