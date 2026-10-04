%{ @file sn_has_lps.m
 %  @brief Checks if the network uses LPS scheduling
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network uses LPS scheduling
 %
 % @details
 % Returns true if any station in the network uses Least Progress Scheduling.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_lps(sn)
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
 % <tr><td>bool<td>True if any station uses LPS scheduling
 % </table>
%}
function bool = sn_has_lps(sn)

bool = any(sn.sched==SchedStrategy.LPS);
end
