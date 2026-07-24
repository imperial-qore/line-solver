%{ @file sn_has_ps_prio.m
 %  @brief Checks if the network uses PS with priorities
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network uses PS with priorities
 %
 % @details
 % Returns true if any station uses Processor Sharing with priorities.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_ps_prio(sn)
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
 % <tr><td>bool<td>True if any station uses PS with priorities
 % </table>
%}
function bool = sn_has_ps_prio(sn)

bool = any(sn.sched==SchedStrategy.PSPRIO);
end