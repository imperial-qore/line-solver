%{ @file sn_is_open_model.m
 %  @brief Checks if the network is an open model
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network is an open model
 %
 % @details
 % Returns true if all classes in the network have infinite population.
 %
 % @par Syntax:
 % @code
 % bool = sn_is_open_model(sn)
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
 % <tr><td>bool<td>True if all classes have infinite population
 % </table>
%}
function bool = sn_is_open_model(sn)
bool = all(isinf(sn.njobs));
end