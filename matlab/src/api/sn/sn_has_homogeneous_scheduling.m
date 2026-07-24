%{ @file sn_has_homogeneous_scheduling.m
 %  @brief Checks if all stations use the same scheduling policy
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if all stations use the same scheduling policy
 %
 % @details
 % Returns true if all stations in the network use identical scheduling policies.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_homogeneous_scheduling(sn)
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
 % <tr><td>bool<td>True if all stations use the same scheduling policy
 % </table>
%}
function bool = sn_has_homogeneous_scheduling(sn, strategy)

bool = length(findstring(sn.sched,strategy)) == sn.nstations;
end