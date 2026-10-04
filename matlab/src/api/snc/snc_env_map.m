%{ @file snc_env_map.m
 %  @brief MGF arrival envelope of a MAP/MMPP flow with unit-size jobs
 %
 %  @author LINE Development Team
%}

%{
 % @brief MGF arrival envelope of a MAP/MMPP flow with unit-size jobs
 %
 % @details
 % For a Markovian arrival process (D0,D1) counting N(0,t) unit-work jobs,
 %
 %   E[exp(theta*N(0,t))] = pi*expm((D0+D1*exp(theta))*t)*1.
 %
 % Let lambda* be the eigenvalue of maximal real part of A(theta)=D0+D1*e^theta
 % and v>0 its right Perron eigenvector. Bounding 1 <= v/min(v) entrywise and
 % using the nonnegativity of expm(A*t) off the diagonal gives
 %
 %   rho(theta) = lambda*/theta,   sigma(theta) = log(max(v)/min(v))/theta,
 %
 % which is the standard exponential-form envelope of a Markov-modulated
 % source. The burst term is what the modulating chain contributes: it is 0 for
 % a one-phase MAP (Poisson) and grows with the phase disparity of an MMPP.
 %
 % @par Syntax:
 % @code
 % [sigma,rho] = snc_env_map(D0,D1,theta)
 % arv = @(theta) snc_env_map(D0,D1,theta);
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>D0<td>Hidden-transition generator block of the MAP
 % <tr><td>D1<td>Arrival-transition block of the MAP
 % <tr><td>theta<td>Chernoff parameter, theta > 0
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sigma<td>Burst term contributed by the modulating chain
 % <tr><td>rho<td>Rate term of the envelope, jobs per slot
 % </table>
 %
 % @par Reference:
 % C.-S. Chang, "Performance Guarantees in Communication Networks", Springer,
 % 2000, Ch. 7 (effective bandwidth of a Markov-modulated source).
%}
function [sigma, rho] = snc_env_map(D0, D1, theta)

if nargin < 3
    line_error(mfilename, 'Usage: [sigma,rho] = snc_env_map(D0,D1,theta).');
end
if theta <= 0
    line_error(mfilename, 'theta must be positive. Got %g.', theta);
end
if any(size(D0) ~= size(D1)) || size(D0,1) ~= size(D0,2)
    line_error(mfilename, 'D0 and D1 must be square and of equal size.');
end

A = full(D0) + full(D1) * exp(theta);
[V, L] = eig(A);
[~, imax] = max(real(diag(L)));
lstar = real(L(imax, imax));
v = real(V(:, imax));
if max(v) < 0
    v = -v; % eigenvector sign is arbitrary, take the positive representative
end
if min(v) <= 0
    line_error(mfilename, 'MAP is not irreducible: the Perron eigenvector is not positive.');
end

rho = lstar / theta;
sigma = log(max(v) / min(v)) / theta;
end
