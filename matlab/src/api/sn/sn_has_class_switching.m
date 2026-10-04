%{ @file sn_has_class_switching.m
 %  @brief Checks if the network has class switching
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network has class switching
 %
 % @details
 % Returns true if the network has class switching, indicated by the number
 % of classes differing from the number of chains.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_class_switching(sn)
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
 % <tr><td>bool<td>True if number of classes differs from number of chains
 % </table>
%}
function bool = sn_has_class_switching(sn)

bool = sn.nclasses ~= sn.nchains;
end