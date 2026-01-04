%{ @file sn_has_priorities.m
 %  @brief Checks if the network uses priority scheduling
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network uses priority scheduling
 %
 % @details
 % Returns true if any station uses priority-based scheduling.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_priorities(sn)
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
 % <tr><td>bool<td>True if any station uses priority scheduling
 % </table>
%}
function bool = sn_has_priorities(sn)
bool = any(sn.classprio(:) > 0);
end