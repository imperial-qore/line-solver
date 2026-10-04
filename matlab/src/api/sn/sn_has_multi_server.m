%{ @file sn_has_multi_server.m
 %  @brief Checks if the network has multi-server stations
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network has multi-server stations
 %
 % @details
 % Returns true if any station has more than one server.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_multi_server(sn)
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
 % <tr><td>bool<td>True if any station has multiple servers
 % </table>
%}
function bool = sn_has_multi_server(sn)

bool = any(sn.nservers(isfinite(sn.nservers)) > 1);
end