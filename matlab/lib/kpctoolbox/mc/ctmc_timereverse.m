%{ @file ctmc_timereverse.m
 %  @brief Computes the time-reversed generator of a CTMC
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes the time-reversed generator of a CTMC
 %
 % @details
 % Computes the infinitesimal generator of the time-reversed process.
 %
 % @par Syntax:
 % @code
 % Qrev = ctmc_timereverse(Q)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Q<td>Infinitesimal generator matrix of the original process
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Qrev<td>Infinitesimal generator matrix of the time-reversed process
 % </table>
%}
function Qrev=ctmc_timereverse(Q)
K=length(Q);
Qrev=Q;
pie=ctmc_solve(Q);
for i=1:K
    for j=1:K
    Qrev(i,j)=Q(i,j)*pie(i)/pie(j);
    end
end
Qrev=Qrev';
end