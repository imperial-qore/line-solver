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
 % <tr><td>q<td>(Optional) Uniformization rate. Default: max(|diag(Q)|) + rand
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
        q = (max(max(abs(Q)))) + rand;
    end
    if issparse(Q)
        P = Q / q + speye(size(Q));
    else
        P = Q / q + eye(size(Q));
    end
    P = dtmc_makestochastic(P);
end
