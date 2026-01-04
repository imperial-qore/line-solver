%{ @file sn_has_sd_routing.m
 %  @brief Checks if the network has state-dependent routing strategies
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network has state-dependent routing strategies
 %
 % @details
 % Returns true if the network contains routing strategies that are state-dependent
 % and therefore violate the product-form assumption. State-dependent routing
 % strategies include Round-Robin, Weighted Round-Robin, Join Shortest Queue,
 % Power of K Choices, and Reinforcement Learning.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_sd_routing(sn)
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
 % <tr><td>bool<td>True if the network has state-dependent routing
 % </table>
%}
function bool = sn_has_sd_routing(sn)

% Product-form requires state-independent (Markovian) routing.
% PROB and RAND are product-form compatible.
% RROBIN, WRROBIN, JSQ, KCHOICES, RL are state-dependent and violate product-form.

if isempty(sn.routing)
    bool = false;
    return;
end

% Check for any state-dependent routing strategy
bool = any(sn.routing(:) == RoutingStrategy.RROBIN | ...
           sn.routing(:) == RoutingStrategy.WRROBIN | ...
           sn.routing(:) == RoutingStrategy.JSQ | ...
           sn.routing(:) == RoutingStrategy.KCHOICES | ...
           sn.routing(:) == RoutingStrategy.RL);
end
