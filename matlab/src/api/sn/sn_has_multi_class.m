%{ @file sn_has_multi_class.m
 %  @brief Checks if the network has multiple classes
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network has multiple classes
 %
 % @details
 % Returns true if the network has more than one job class.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_multi_class(sn)
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
 % <tr><td>bool<td>True if number of classes is greater than one
 % </table>
%}
function bool = sn_has_multi_class(sn)

bool = sn.nclasses > 1;
end