%{ @file sn_is_mixed_model.m
 %  @brief Checks if the network is a mixed model
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network is a mixed model
 %
 % @details
 % Returns true if the network contains both open and closed classes.
 %
 % @par Syntax:
 % @code
 % bool = sn_is_mixed_model(sn)
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
 % <tr><td>bool<td>True if the network contains both open and closed classes
 % </table>
%}
function bool = sn_is_mixed_model(sn)
bool = sn_has_mixed_classes(sn);
end
