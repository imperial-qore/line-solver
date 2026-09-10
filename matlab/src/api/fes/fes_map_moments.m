%{ @file fes_map_moments.m
 %  @brief Moments and index of dispersion of an inter-departure MAP
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes the first three moments, the lag-1 joint moment and the
 % index of dispersion of a MAP without forming (-T0)^-1
 %
 % @details
 % Evaluates equations (4), (5) and (7) of Casale, Mi, Cherkasova and
 % Smirni, IEEE Trans. Soft. Eng. 37(5), 2011, on the block bidiagonal MAP
 % returned by fes_map_interdeparture. The inverse (-T0)^-1 is dense even
 % when T0 is sparse, so it is never formed: the moments are obtained by
 % the vector recursion v_{k+1} = v_k (-T0)^-1, each step being a sparse
 % linear solve. Method 'euler' replaces the solve by the quadrature
 %
 %   v (-T0)^-1 = v * int_0^inf exp(T0 t) dt
 %
 % of Section 5.2.2, integrated by the trapezoid rule with the Euler
 % approximation exp(T0 dt) ~ I + T0 dt and step dt below the inverse of
 % the largest diagonal element in absolute value, as in the uniformization
 % method. Method 'ssolve' is the default because it is exact and faster;
 % 'euler' reproduces the reference implementation of the paper.
 %
 % @par Syntax:
 % @code
 % [e1,e2,e3,e11,idc] = fes_map_moments(T0,T1)
 % [e1,e2,e3,e11,idc] = fes_map_moments(T0,T1,options)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>T0<td>Hidden transitions of the MAP
 % <tr><td>T1<td>Marked transitions of the MAP
 % <tr><td>options<td>(Optional) struct with fields method ('ssolve' or
 %                    'euler'), tol (default 1e-12), step_safety (default
 %                    0.1) and iter_max (default 1e6)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>e1<td>Mean inter-departure time
 % <tr><td>e2<td>Second moment of the inter-departure times
 % <tr><td>e3<td>Third moment of the inter-departure times
 % <tr><td>e11<td>Joint moment E[X_k X_{k+1}] of consecutive samples
 % <tr><td>idc<td>Asymptotic index of dispersion
 % </table>
 %
 % @see fes_map_interdeparture, map2_fit_idc, map_idc
%}
function [e1,e2,e3,e11,idc] = fes_map_moments(T0,T1,options)

if nargin < 3
    options = struct();
end
if ~isfield(options,'method') || isempty(options.method)
    options.method = 'ssolve';
end
if ~isfield(options,'tol')
    options.tol = 1e-12;
end
if ~isfield(options,'step_safety')
    options.step_safety = 0.1;
end
if ~isfield(options,'iter_max')
    options.iter_max = 1e6;
end

dim = size(T0,1);
e = ones(dim,1);
Q = T0 + T1;

phi = ctmc_solve(Q);
phi = reshape(phi, 1, dim);
pie = phi*T1;
lambda = sum(pie);
pie = pie/lambda;

switch lower(options.method)
    case 'ssolve'
        solve = @(v) v/(-T0);
    case 'euler'
        dt = options.step_safety/max(abs(diag(T0)));
        solve = @(v) fes_map_euler(v, T0, dt, options.tol, options.iter_max);
    otherwise
        line_error(mfilename, sprintf('Unknown method %s, use ssolve or euler.', options.method));
end

v1 = solve(pie);      e1 = v1*e;
v2 = solve(v1);       e2 = 2*(v2*e);
v3 = solve(v2);       e3 = 6*(v3*e);
v4 = solve(v2*T1);    e11 = v4*e;

% equation (7), with pie*inv(Q+e*phi) obtained from the rank-one update
% y*Q = pie-phi under the normalization y*e = 1
A = Q;
A(:,dim) = 1;
rhs = pie - phi;
rhs(dim) = 1;
y = rhs/A;
idc = 1 + 2*(lambda - y*T1*e);
end
