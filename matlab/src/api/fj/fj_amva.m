%{ @file fj_amva.m
 %  @brief Mean value analysis of a closed network of fork-join subnetworks
 %
 %  @author LINE Development Team
%}

%{
 % @brief Mean value analysis of a closed network of fork-join subnetworks
 %
 % @details
 % Population-by-population recursion for a single-class closed queueing
 % network whose stations are P(n)-way fork-join subnetworks with per-visit
 % demand D(n). The response time of a subnetwork inflates the arrival-instant
 % queue length by the harmonic number of its fork degree, which is the
 % Varki bound on the residence time of a parallel subsystem:
 %
 %   R_n(m) = D_n * [ H_{P_n} + Q_n(m-1) ]
 %   X(m)   = m / (Z + sum_n R_n(m))
 %   Q_n(m) = X(m) * R_n(m)
 %
 % started from Q_n(0) = 0. Setting P(n) = 1 for every n recovers the exact
 % single-class mean value analysis of Reiser and Lavenberg, because H_1 = 1.
 % For P(n) > 1 the recursion is an approximation whose per-subnetwork
 % residence time is an upper bound in the sense of Varki.
 %
 % @par Syntax:
 % @code
 % [R, Q, X, U] = fj_amva(D, P, M)
 % [R, Q, X, U] = fj_amva(D, P, M, Z)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>D<td>Vector of per-visit service demands, one per subnetwork
 % <tr><td>P<td>Vector of fork degrees, one per subnetwork (P(n) >= 1)
 % <tr><td>M<td>Number of circulating jobs (positive integer)
 % <tr><td>Z<td>Think time (optional, default 0)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>R<td>Vector of mean residence times, one per subnetwork
 % <tr><td>Q<td>Vector of mean queue lengths, one per subnetwork
 % <tr><td>X<td>System throughput
 % <tr><td>U<td>Vector of utilizations of the busiest queue of each subnetwork
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Eqs. (68)-(69).
 %
 % Original: E. Varki, "Mean Value Technique for Closed Fork-Join Networks",
 % ACM SIGMETRICS, 1999; G. Casale, R. R. Muntz, G. Serazzi, "Geometric
 % Bounds", IEEE Trans. Computers 57(6), 2008.
%}
function [R, Q, X, U] = fj_amva(D, P, M, Z)

if nargin < 4 || isempty(Z)
    Z = 0;
end

D = D(:)';
P = P(:)';

if numel(D) ~= numel(P)
    line_error(mfilename, 'D and P must have the same number of elements. Got %d and %d.', numel(D), numel(P));
end
if any(D < 0)
    line_error(mfilename, 'Service demands must be non-negative.');
end
if any(P < 1) || any(P ~= round(P))
    line_error(mfilename, 'Fork degrees P must be positive integers.');
end
if M < 1 || M ~= round(M)
    line_error(mfilename, 'M must be a positive integer. Got M=%g.', M);
end
if Z < 0
    line_error(mfilename, 'Think time Z must be non-negative.');
end

N = numel(D);
H = zeros(1, N);
for n = 1:N
    H(n) = fj_harmonic(P(n));
end

Q = zeros(1, N);
R = zeros(1, N);
X = 0;
for m = 1:M
    R = D .* (H + Q);
    Rtot = sum(R);
    if Rtot <= 0
        line_error(mfilename, 'Total residence time vanished at population %d; all demands are zero.', m);
    end
    X = m / (Z + Rtot);
    Q = X * R;
end

% Each subnetwork holds P(n) queues sharing the demand equally, so the
% per-queue utilization is the subnetwork demand divided by the fork degree
U = X * D ./ P;

end
