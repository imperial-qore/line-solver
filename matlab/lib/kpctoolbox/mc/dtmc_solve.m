%{ @file dtmc_solve.m
 %  @brief Equilibrium distribution of the discrete-time Markov chain
 %
 %  @author LINE Development Team
%}

%{
 % @brief Equilibrium distribution of the discrete-time Markov chain
 %
 % @details
 % Calculates the equilibrium distribution of a discrete-time Markov chain given its transition matrix.
 %
 % @par Syntax:
 % @code
 % PROB = dtmc_solve(P)
 % PROB = dtmc_solve(P, options)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>P<td>Stochastic transition matrix of the discrete-time Markov chain
 % <tr><td>options<td>(Optional) Solver options
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>PROB<td>Equilibrium distribution vector
 % </table>
 %
 % @par Examples:
 % @code
 % PROB = dtmc_solve([0.5,0.5;0.2,0.8]);
 % @endcode
%}
function PROB=dtmc_solve(P,options)

if nargin<2
    if issparse(P)
        PROB=ctmc_solve(P-speye(size(P)));
    else
        PROB=ctmc_solve(P-eye(size(P)));
    end
else
    if issparse(P)
        PROB=ctmc_solve(P-speye(size(P)),options);
    else
        PROB=ctmc_solve(P-eye(size(P)),options);
    end
end
end