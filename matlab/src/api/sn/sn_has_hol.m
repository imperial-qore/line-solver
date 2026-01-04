%{ @file sn_has_hol.m
 %  @brief Checks if the network uses HOL scheduling
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network uses HOL scheduling
 %
 % @details
 % Returns true if any station in the network uses Head-Of-Line scheduling.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_hol(sn)
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
 % <tr><td>bool<td>True if any station uses HOL scheduling
 % </table>
%}
function bool = sn_has_hol(sn)

bool = any(sn.sched==SchedStrategy.HOL);
end