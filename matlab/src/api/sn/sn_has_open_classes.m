%{ @file sn_has_open_classes.m
 %  @brief Checks if the network contains open classes
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network contains open classes
 %
 % @details
 % Returns true if the network contains any open class (infinite population).
 %
 % @par Syntax:
 % @code
 % bool = sn_has_open_classes(sn)
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
 % <tr><td>bool<td>True if any class has infinite population
 % </table>
%}
function bool = sn_has_open_classes(sn)

bool = any(isinf(sn.njobs));
end