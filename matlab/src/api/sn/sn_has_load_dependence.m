%{ @file sn_has_load_dependence.m
 %  @brief Checks if the network has load-dependent service rates
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network has load-dependent service rates
 %
 % @details
 % Returns true if any station has service rates that depend on queue length.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_load_dependence(sn)
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
 % <tr><td>bool<td>True if any station has load-dependent service
 % </table>
%}
function bool = sn_has_load_dependence(sn)

bool = size(sn.lldscaling,2)>0;
end