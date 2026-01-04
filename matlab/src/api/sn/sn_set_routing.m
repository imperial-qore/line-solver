%{ @file sn_set_routing.m
 %  @brief Sets the routing matrix
 %
 %  @author LINE Development Team
%}

%{
 % @brief Sets the routing matrix for stateful nodes
 %
 % @details
 % Replaces the rt matrix in NetworkStruct.
 %
 % @par Syntax:
 % @code
 % sn = sn_set_routing(sn, rt)
 % sn = sn_set_routing(sn, rt, autoRefresh)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>rt<td>New routing matrix (nstateful*nclasses x nstateful*nclasses)
 % <tr><td>autoRefresh<td>If true, refresh visit ratios (default false)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Modified network structure
 % </table>
%}
function sn = sn_set_routing(sn, rt, autoRefresh)

if nargin < 3 || isempty(autoRefresh)
    autoRefresh = false;
end

% Update rt matrix
sn.rt = rt;

% Auto-refresh visit ratios if requested
if autoRefresh
    sn = sn_refresh_visits(sn, sn.chains, sn.rt, sn.rtnodes);
end

end
