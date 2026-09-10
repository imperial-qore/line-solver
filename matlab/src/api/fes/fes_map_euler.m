%{ @file fes_map_euler.m
 %  @brief Quadrature of v*(-T0)^-1 by uniformized Euler integration
 %
 %  @author LINE Development Team
%}

%{
 % @brief Integrates v*int_0^inf exp(T0 t) dt by the trapezoid rule
 %
 % @details
 % Evaluates the product v*(-T0)^-1 without any factorization, as done in
 % Section 5.2.2 of Casale, Mi, Cherkasova and Smirni, IEEE Trans. Soft.
 % Eng. 37(5), 2011. The propagated vector uses the Euler approximation
 % exp(T0 dt) ~ I + T0 dt, so only sparse vector-matrix products are
 % performed and the sparsity of T0 is preserved throughout. The
 % integration stops when the propagated vector has lost the fraction tol
 % of its initial mass.
 %
 % @par Syntax:
 % @code
 % y = fes_map_euler(v, T0, dt, tol, iter_max)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>v<td>Row vector to be multiplied by (-T0)^-1
 % <tr><td>T0<td>Hidden transitions of the MAP, a stable matrix
 % <tr><td>dt<td>Integration step, below 1/max(abs(diag(T0)))
 % <tr><td>tol<td>Relative mass left when the integration stops
 % <tr><td>iter_max<td>Maximum number of integration steps
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>y<td>Row vector approximating v*(-T0)^-1
 % </table>
 %
 % @see fes_map_moments
%}
function y = fes_map_euler(v, T0, dt, tol, iter_max)

y = zeros(size(v));
z = v;
nrm0 = norm(v,1);
converged = false;

for it = 1:iter_max
    znext = z + dt*(z*T0);
    y = y + dt*(z + znext)/2;
    z = znext;
    if norm(z,1) <= tol*nrm0
        converged = true;
        break
    end
end

if ~converged
    line_warning(mfilename,'The Euler quadrature did not converge in %d steps, %.3e of the mass is left.\n', iter_max, norm(z,1)/nrm0);
end
end
