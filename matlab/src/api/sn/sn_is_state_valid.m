%{ @file sn_is_state_valid.m
 %  @brief Validates a network state
 %
 %  @author LINE Development Team
%}

%{
 % @brief Validates a network state
 %
 % @details
 % Checks if a given state vector represents a valid state for the network.
 %
 % @par Syntax:
 % @code
 % bool = sn_is_state_valid(sn)
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
 % <tr><td>bool<td>True if the state is valid for this network
 % </table>
%}
function [isvalid] = sn_is_state_valid(sn)
% [ISVALID] = ISSTATEVALID()
nir = [];
sir = [];
for ist=1:sn.nstations
    isf = sn.stationToStateful(ist);
    if size(sn.state{isf},1)>1
        line_warning(mfilename,sprintf('isStateValid will ignore some states of station %d, define a unique initial state to address this problem.\n',ist));
        sn.state{isf} = sn.state{isf}(1,:);
    end
    [~, nir(ist,:), sir(ist,:), ~] = State.toMarginal(sn, sn.stationToNode(ist), sn.state{isf});
end
isvalid = State.isValid(sn, nir, sir);
end