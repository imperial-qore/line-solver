%{ @file dmap_to_dph.m
 %  @brief Discrete phase-type law underlying a renewal DMAP
 %
 %  @author LINE Development Team
%}

%{
 % @brief Discrete phase-type law underlying a renewal DMAP
 %
 % @details
 % Inverts dph_to_dmap: with D1 = a*alpha of rank one, alpha is recovered by
 % normalizing any nonzero row of D1 and A is D0. Errors on a DMAP that does
 % not renew, since no discrete phase-type law then describes its interevent
 % times.
 %
 % @par Syntax:
 % @code
 % [alpha, A] = dmap_to_dph(DMAP)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>DMAP<td>Cell {D0, D1} of a renewal DMAP
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>alpha<td>1xm initial phase probability row vector
 % <tr><td>A<td>mxm substochastic transient matrix
 % </table>
%}
function [alpha, A] = dmap_to_dph(DMAP)

if ~dmap_is_renewal(DMAP)
    line_error(mfilename, 'The DMAP does not renew at events, so it has no discrete phase-type form.');
end

A = DMAP{1};
D1 = DMAP{2};
rowMass = sum(D1, 2);
[~, pivotRow] = max(rowMass);
if rowMass(pivotRow) <= GlobalConstants.Zero
    line_error(mfilename, 'The DMAP has no events: D1 is the zero matrix.');
end
alpha = D1(pivotRow, :) / rowMass(pivotRow);

end
