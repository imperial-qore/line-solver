%{ @file dtmc_timereverse.m
 %  @brief Computes the time-reversed transition matrix of a DTMC
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes the time-reversed transition matrix of a DTMC
 %
 % @details
 % Computes the transition matrix of the time-reversed process.
 %
 % @par Syntax:
 % @code
 % Prev = dtmc_timereverse(P)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>P<td>Stochastic transition matrix of the original process
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Prev<td>Stochastic transition matrix of the time-reversed process
 % </table>
%}
function Prev=dtmc_timereverse(P)
K=length(P);
Prev=P;
pie=dtmc_solve(P);
for i=1:K
    for j=1:K
    Prev(i,j)=P(i,j)*pie(i)/pie(j);
    end
end
Prev=Prev';
end