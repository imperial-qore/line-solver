%{ @file dmap_lambda.m
 %  @brief Event rate per slot of a discrete-time batch arrival stream
 %
 %  @author LINE Development Team
%}

%{
 % @brief Event rate per slot of a discrete-time batch arrival stream
 %
 % @details
 % Returns pi * sum_k k*A_k * e with pi the stationary vector of the phase
 % chain sum_k A_k. This counts EVENTS per slot, not slots with at least one
 % event, so it is the quantity that Little's law consumes downstream. For a
 % plain DMAP {D0,D1} it reduces to pi*D1*e.
 %
 % @par Syntax:
 % @code
 % lambda = dmap_lambda(A)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>A<td>Cell {A_0, A_1, ...} of batch matrices
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lambda<td>Mean number of events per slot
 % </table>
%}
function lambda = dmap_lambda(A)

m = size(A{1}, 1);
P = zeros(m, m);
for k = 1:length(A)
    P = P + A{k};
end
pi = dtmc_solve(P);
W = zeros(m, m);
for k = 2:length(A)
    W = W + (k-1) * A{k};
end
lambda = pi * W * ones(m, 1);

end
