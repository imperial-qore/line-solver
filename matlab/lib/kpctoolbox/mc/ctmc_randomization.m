%{ @file ctmc_randomization.m
 %  @brief Uniformization (randomization) of a continuous-time Markov chain
 %
 %  @author LINE Development Team
%}

%{
 % @brief Uniformization (randomization) of a continuous-time Markov chain
 %
 % @details
 % Applies uniformization to transform a CTMC into a DTMC.
 %
 % @par Syntax:
 % @code
 % [P, q] = ctmc_randomization(Q)
 % [P, q] = ctmc_randomization(Q, q)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Q<td>Infinitesimal generator matrix
 % <tr><td>q<td>(Optional) Uniformization rate. Default: 1.05*max(|Q|)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>P<td>Uniformized discrete-time stochastic matrix
 % <tr><td>q<td>The rate used for uniformization
 % </table>
%}
function [P, q] = ctmc_randomization(Q, q)
    if nargin == 1
        % Deterministic default: an unseeded rand made P, and everything
        % downstream of it, irreproducible run to run. 1.05*max|Q| is the
        % rate ctmc_courtois already derives, and every quantity taken
        % from P is invariant to the choice of q above max|q_ii|.
        q = 1.05 * max(max(abs(Q)));
    end
    if issparse(Q)
        P = Q / q + speye(size(Q));
    else
        P = Q / q + eye(size(Q));
    end
    P = dtmc_makestochastic(P);
end
