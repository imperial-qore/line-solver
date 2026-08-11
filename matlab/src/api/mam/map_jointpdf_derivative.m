%{ @file map_jointpdf_derivative.m
 %  @brief Partial derivative at zero of a MAP's joint PDF
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes partial derivative at zero of MAP joint PDF
 %
 % @details
 % This function computes the partial derivative at 0 of a MAP's joint
 % probability density function. Based on A. Horvath et al.
 % "A Joint Moments Based Analysis of Networks of MAP/MAP/1 Queues".
 %
 % @par Syntax:
 % @code
 % gamma = map_jointpdf_derivative(MAP, iset)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>MAP<td>Markovian Arrival Process {D0, D1}
 % <tr><td>iset<td>Vector of indices defining the partial derivative
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Calculated derivative value
 % </table>
%}
function gamma = map_jointpdf_derivative(MAP, iset)
% partial derivative at 0 of a MAP's joint PDF
% A. Horvath et al. A Joint Moments Based Analysis of Networks of
% MAP/MAP/1 Queues

n = length(MAP{1});

gamma = map_pie(MAP);
for j=iset(:)'
    gamma = gamma * MAP{1}^j * MAP{2};
end
gamma = gamma * ones(n,1);
end
