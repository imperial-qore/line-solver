%{ @file sn_has_single_class.m
 %  @brief Checks if the network has a single class
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network has a single class
 %
 % @details
 % Returns true if the network has exactly one job class.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_single_class(sn)
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
 % <tr><td>bool<td>True if number of classes equals one
 % </table>
%}
function bool = sn_has_single_class(sn)

bool = sn.nclasses == 1;
end