%{ @file sn_has_dps.m
 %  @brief Checks if the network uses Discriminatory Processor Sharing
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network uses Discriminatory Processor Sharing
 %
 % @details
 % Returns true if any station in the network uses DPS scheduling.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_dps(sn)
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
 % <tr><td>bool<td>True if any station uses DPS scheduling
 % </table>
%}
function bool = sn_has_dps(sn)

bool = any(sn.sched==SchedStrategy.DPS);
end