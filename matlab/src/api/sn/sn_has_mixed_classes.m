%{ @file sn_has_mixed_classes.m
 %  @brief Checks if the network has both open and closed classes
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network has both open and closed classes
 %
 % @details
 % Returns true if the network contains at least one open class and one closed class.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_mixed_classes(sn)
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
 % <tr><td>bool<td>True if the network has both open and closed classes
 % </table>
%}
function bool = sn_has_mixed_classes(sn)

bool = sn_has_closed_classes(sn) && sn_has_open_classes(sn);
end
    
