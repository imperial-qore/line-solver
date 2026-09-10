%{ @file dph_to_dmap.m
 %  @brief Renewal DMAP of a discrete phase-type law
 %
 %  @author LINE Development Team
%}

%{
 % @brief Renewal DMAP of a discrete phase-type law
 %
 % @details
 % Maps (alpha, A) to {D0, D1} with D0 = A and D1 = a*alpha, a = e - A*e. The
 % resulting DMAP renews the phase at every event, so its interevent times are
 % i.i.d. copies of the DPH and D0+D1 is stochastic by construction.
 %
 % @par Syntax:
 % @code
 % DMAP = dph_to_dmap(alpha, A)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>alpha<td>1xm initial phase probability row vector
 % <tr><td>A<td>mxm substochastic transient matrix
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>DMAP<td>Cell {D0, D1}
 % </table>
%}
function DMAP = dph_to_dmap(alpha, A)

m = size(A, 1);
a = ones(m, 1) - A * ones(m, 1);
DMAP = {A, a * alpha(:)'};

end
