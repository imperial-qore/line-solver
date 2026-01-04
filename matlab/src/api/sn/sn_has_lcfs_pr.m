%{ @file sn_has_lcfs_pr.m
 %  @brief Checks if the network uses LCFS with preemptive-resume scheduling
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network uses LCFS with preemptive-resume scheduling
 %
 % @details
 % Returns true if any station uses LCFS preemptive-resume scheduling.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_lcfs_pr(sn)
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
 % <tr><td>bool<td>True if any station uses LCFS-PR scheduling
 % </table>
%}
function bool = sn_has_lcfs_pr(sn)

bool = any(sn.sched==SchedStrategy.LCFSPR);
end