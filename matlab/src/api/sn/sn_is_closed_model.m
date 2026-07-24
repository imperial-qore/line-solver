%{ @file sn_is_closed_model.m
 %  @brief Checks if the network is a closed model
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network is a closed model
 %
 % @details
 % Returns true if all classes in the network have finite population.
 %
 % @par Syntax:
 % @code
 % bool = sn_is_closed_model(sn)
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
 % <tr><td>bool<td>True if all classes have a finite population
 % </table>
%}
function bool = sn_is_closed_model(sn)
bool = all(isfinite(sn.njobs));
end