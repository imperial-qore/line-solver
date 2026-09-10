%{ @file sn_has_gps_prio.m
 %  @brief Checks if the network uses GPS with priorities
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network uses GPS with priorities
 %
 % @details
 % Returns true if any station in the network uses Generalized Processor Sharing with priorities.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_gps_prio(sn)
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
 % <tr><td>bool<td>True if any station uses GPS with priorities
 % </table>
%}
function bool = sn_has_gps_prio(sn)

bool = any(sn.sched==SchedStrategy.GPSPRIO);
end