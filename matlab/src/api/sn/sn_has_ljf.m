%{ @file sn_has_ljf.m
 %  @brief Checks if the network uses LJF scheduling
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network uses LJF scheduling
 %
 % @details
 % Returns true if any station in the network uses Longest Job First scheduling.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_ljf(sn)
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
 % <tr><td>bool<td>True if any station uses LJF scheduling
 % </table>
%}
function bool = sn_has_ljf(sn)

bool = any(sn.sched==SchedStrategy.LJF);
end