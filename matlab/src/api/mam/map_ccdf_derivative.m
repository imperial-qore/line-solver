%{ @file map_ccdf_derivative.m
 %  @brief Derivative at zero of a MAP CCDF
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes derivative at zero of MAP complementary CDF
 %
 % @details
 % This function computes the derivative at 0 of a MAP's Complementary
 % Cumulative Distribution Function (CCDF). Based on A. Horvath et al.
 % "A Joint Moments Based Analysis of Networks of MAP/MAP/1 Queues".
 %
 % @par Syntax:
 % @code
 % nu = map_ccdf_derivative(MAP, i)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>MAP<td>Markovian Arrival Process {D0, D1}
 % <tr><td>i<td>Order of the derivative
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>nu<td>Calculated derivative value
 % </table>
%}
function nu = map_ccdf_derivative(MAP, i)
% derivative at 0 of a MAP CCDF
% A. Horvath et al. A Joint Moments Based Analysis of Networks of
% MAP/MAP/1 Queues

n = length(MAP{1});

pie = map_pie(MAP);
nu = pie * MAP{1}^i *ones(n,1);
end
