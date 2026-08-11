%{ @file sn_has_closed_classes.m
 %  @brief Checks if the network contains closed classes
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network contains closed classes
 %
 % @details
 % Returns true if the network contains any closed class (finite population).
 %
 % @par Syntax:
 % @code
 % bool = sn_has_closed_classes(sn)
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
 % <tr><td>bool<td>True if any class has a finite population
 % </table>
%}
function bool = sn_has_closed_classes(sn)

bool = any(isfinite(sn.njobs));
end