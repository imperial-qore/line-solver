%{ @file fes_map_interp_edge.m
 %  @brief Endpoint slope of the monotone cubic Hermite interpolant
 %
 %  @author LINE Development Team
%}

%{
 % @brief Noncentered three-point endpoint slope with monotonicity clamps
 %
 % @details
 % Implements the endpoint rule of de Boor used by shape-preserving cubic
 % Hermite interpolation: the three-point estimate is set to zero when it
 % disagrees in sign with the adjacent secant, and is clipped to three
 % times that secant when the two secants disagree in sign and the estimate
 % is too steep.
 %
 % @par Syntax:
 % @code
 % d = fes_map_interp_edge(h1, h2, del1, del2)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>h1<td>Width of the interval adjacent to the endpoint
 % <tr><td>h2<td>Width of the next interval
 % <tr><td>del1<td>Secant slope adjacent to the endpoint
 % <tr><td>del2<td>Secant slope of the next interval
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>d<td>Endpoint derivative
 % </table>
 %
 % @see fes_map_interp
%}
function d = fes_map_interp_edge(h1, h2, del1, del2)

d = ((2*h1 + h2)*del1 - h1*del2)/(h1 + h2);
if sign(d) ~= sign(del1)
    d = 0;
elseif (sign(del1) ~= sign(del2)) && (abs(d) > abs(3*del1))
    d = 3*del1;
end
end
