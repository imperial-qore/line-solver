%{ @file sn_has_gps.m
 %  @brief Checks if the network uses GPS scheduling
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network uses GPS scheduling
 %
 % @details
 % Returns true if any station in the network uses Generalized Processor Sharing.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_gps(sn)
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
 % <tr><td>bool<td>True if any station uses GPS scheduling
 % </table>
%}
function bool = sn_has_gps(sn)

bool = any(sn.sched==SchedStrategy.GPS);
end