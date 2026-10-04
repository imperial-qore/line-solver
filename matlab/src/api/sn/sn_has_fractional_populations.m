%{ @file sn_has_fractional_populations.m
 %  @brief Checks if the network has fractional population values
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network has fractional population values
 %
 % @details
 % Returns true if any class has a non-integer population value.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_fractional_populations(sn)
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
 % <tr><td>bool<td>True if any class has fractional population
 % </table>
%}
function bool = sn_has_fractional_populations(sn)
% Checks if the network has closed classes with non-integer populations
bool = any(sn.njobs ~= round(sn.njobs));
end