%{ @file sn_set_routing_prob.m
 %  @brief Sets a routing probability between two stateful node-class pairs
 %
 %  @author LINE Development Team
%}

%{
 % @brief Sets a routing probability
 %
 % @details
 % Updates a single entry in the rt matrix.
 %
 % @par Syntax:
 % @code
 % sn = sn_set_routing_prob(sn, fromStateful, fromClass, toStateful, toClass, prob)
 % sn = sn_set_routing_prob(sn, fromStateful, fromClass, toStateful, toClass, prob, autoRefresh)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>fromStateful<td>Source stateful node index (1-based)
 % <tr><td>fromClass<td>Source class index (1-based)
 % <tr><td>toStateful<td>Destination stateful node index (1-based)
 % <tr><td>toClass<td>Destination class index (1-based)
 % <tr><td>prob<td>Routing probability [0, 1]
 % <tr><td>autoRefresh<td>If true, refresh visit ratios (default false)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Modified network structure
 % </table>
%}
function sn = sn_set_routing_prob(sn, fromStateful, fromClass, toStateful, toClass, prob, autoRefresh)

if nargin < 7 || isempty(autoRefresh)
    autoRefresh = false;
end

K = sn.nclasses;

% Calculate indices in rt matrix (1-based)
fromIdx = (fromStateful - 1) * K + fromClass;
toIdx = (toStateful - 1) * K + toClass;

% Update rt matrix
sn.rt(fromIdx, toIdx) = prob;

% Auto-refresh visit ratios if requested
if autoRefresh
    [~, ~, sn] = sn_refresh_visits(sn, sn.chains, sn.rt, sn.rtnodes);
end

end
