%{ @file sn_get_state_aggr.m
 %  @brief Returns the aggregated initial state of the network
 %
 %  @author LINE Development Team
%}

%{
 % @brief Returns the aggregated initial state of the network
 %
 % @details
 % This function extracts and aggregates the initial state of all stateful
 % nodes in the network.
 %
 % @par Syntax:
 % @code
 % initialStateAggr = sn_get_state_aggr(sn)
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
 % <tr><td>initialStateAggr<td>Aggregated initial state cell array
 % </table>
%}
function [initialStateAggr] = sn_get_state_aggr(sn)

initialState = sn.state;
initialStateAggr = cell(size(initialState));
for isf=1:length(initialStateAggr)
    ind = sn.statefulToNode(isf);
    [~,initialStateAggr{isf}] = State.toMarginalAggr(sn, ind, initialState{isf});
end
end