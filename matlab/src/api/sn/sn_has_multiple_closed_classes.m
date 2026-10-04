%{ @file sn_has_multiple_closed_classes.m
 %  @brief Checks if the network has multiple closed classes
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network has multiple closed classes
 %
 % @details
 % Returns true if the network has more than one closed (finite population) class.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_multiple_closed_classes(sn)
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
 % <tr><td>bool<td>True if there are multiple closed classes
 % </table>
%}
function bool = sn_has_multiple_closed_classes(sn)

bool = nnz(isfinite(sn.njobs))>1;
end