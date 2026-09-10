%{ @file fj_serialization.m
 %  @brief Blocking probability and pseudoserver delay of serialization phases
 %
 %  @author LINE Development Team
%}

%{
 % @brief Blocking probability and pseudoserver delay of serialization phases
 %
 % @details
 % A serialization phase is a stretch of a job's execution that is protected by
 % an exclusive lock, so at most one of the M circulating jobs may occupy it.
 % The queueing network is no longer product form, and the delay in entering
 % phase s is represented by a pseudoserver that is bypassed when the phase is
 % free. A job entering phase s is blocked when at least one of the other M-1
 % jobs is inside it, and treating those jobs as independently placed in
 % proportion to the residence times gives
 %
 %   P_s(M) = 1 - [ 1 - R_s(M)/R(M) ]^(M-1),   R(M) = sum_s R_s(M),
 %
 % where the sum runs over the nonserialized phase and every serialization
 % phase. The delay charged at the pseudoserver is alpha*R_s(M), with alpha
 % depending on the residence time distribution and on where within the phase
 % the blocked job arrives; alpha = 1/2 is the value for an arrival uniform in
 % a phase of low utilization, which is the regime in which the approximation
 % is stated.
 %
 % @par Syntax:
 % @code
 % P = fj_serialization(Rs, R0, M)
 % [P, delay, Rtot] = fj_serialization(Rs, R0, M, alpha)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Rs<td>Vector of mean residence times inside each serialization phase
 % <tr><td>R0<td>Mean residence time in the nonserialized phase
 % <tr><td>M<td>Number of circulating jobs (positive integer)
 % <tr><td>alpha<td>Fraction of the phase charged to a blocked job (optional, default 0.5)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>P<td>Vector of probabilities of being blocked on entering each phase
 % <tr><td>delay<td>Vector of expected pseudoserver delays, P(s)*alpha*Rs(s)
 % <tr><td>Rtot<td>Mean cycle time including the serialization delays
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Section 7.1
 % on page 17:43.
 %
 % Original: A. Thomasian, "Queueing Network Models to Estimate Serialization
 % Delays in Computer Systems", Performance, 1983.
%}
function [P, delay, Rtot] = fj_serialization(Rs, R0, M, alpha)

if nargin < 4 || isempty(alpha)
    alpha = 0.5;
end

Rs = Rs(:)';
S = numel(Rs);

if S < 1
    line_error(mfilename, 'At least one serialization phase must be supplied.');
end
if any(Rs < 0)
    line_error(mfilename, 'The residence times inside the serialization phases must be non-negative.');
end
if R0 < 0
    line_error(mfilename, 'The nonserialized residence time must be non-negative. Got R0=%g.', R0);
end
if M < 1 || M ~= round(M)
    line_error(mfilename, 'M must be a positive integer. Got M=%g.', M);
end
if alpha < 0 || alpha > 1
    line_error(mfilename, 'alpha must lie in [0,1]. Got alpha=%g.', alpha);
end

% Mean residence time over the nonserialized phase and every serialized phase
R = R0 + sum(Rs);
if R <= 0
    line_error(mfilename, 'The total residence time vanished; every phase has zero demand.');
end

P = 1 - (1 - Rs / R).^(M - 1);
delay = P .* (alpha * Rs);
Rtot = R + sum(delay);

end
