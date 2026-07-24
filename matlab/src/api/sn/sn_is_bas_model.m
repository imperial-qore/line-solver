%{ @file sn_is_bas_model.m
 %  @brief Checks if the network is a closed single-class Blocking-After-Service model
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network is a closed single-class BAS model
 %
 % @details
 % Returns true for a closed, single-class network with Blocking-After-Service
 % (BAS) finite-buffer blocking, which solver_sqd handles but exact/AMVA MVA
 % does not. The Smith queue-decomposition approximation models a single
 % circulating population, so an open class or more than one class disqualifies
 % the model.
 %
 % @par Syntax:
 % @code
 % bool = sn_is_bas_model(sn)
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
 % <tr><td>bool<td>True if the model is a closed single-class BAS model
 % </table>
%}
function bool = sn_is_bas_model(sn)
bool = false;
if sn.nclasses ~= 1 || sn.nclosedjobs <= 0 || isempty(sn.droprule)
    return;
end
if any(isinf(sn.njobs))
    return; % open class present
end
bool = any(sn.droprule(:) == DropStrategy.BAS);
end
