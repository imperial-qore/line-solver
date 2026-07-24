%{ @file sn_has_ps.m
 %  @brief Checks if the network uses PS scheduling
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network uses PS scheduling
 %
 % @details
 % Returns true if any station in the network uses Processor Sharing scheduling.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_ps(sn)
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
 % <tr><td>bool<td>True if any station uses PS scheduling
 % </table>
%}
function bool = sn_has_ps(sn)

bool = any(sn.sched==SchedStrategy.PS);
end