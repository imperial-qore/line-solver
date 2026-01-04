%{ @file sn_set_servers.m
 %  @brief Sets the number of servers at a station
 %
 %  @author LINE Development Team
%}

%{
 % @brief Sets the number of servers at a station
 %
 % @details
 % Directly modifies the server count in NetworkStruct.
 %
 % @par Syntax:
 % @code
 % sn = sn_set_servers(sn, stationIdx, nServers)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>stationIdx<td>Station index (1-based)
 % <tr><td>nServers<td>Number of servers (positive, or Inf)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Modified network structure
 % </table>
%}
function sn = sn_set_servers(sn, stationIdx, nServers)

sn.nservers(stationIdx) = nServers;

end
