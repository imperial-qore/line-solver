%{ @file fj_qgb.m
 %  @brief Geometric bound on the queue length of a fork-join subnetwork
 %
 %  @author LINE Development Team
%}

%{
 % @brief Geometric bound on the queue length of a fork-join subnetwork
 %
 % @details
 % Non-iterative geometric bound on the mean queue length of each fork-join
 % subnetwork of a closed queueing network. Subnetwork n consists of P(n)
 % parallel queues traversed as a P(n)-way fork-join request, and carries a
 % per-visit service demand D(n).
 %
 %   y_n(M) = D_n * M / (Z + sum_j D_j * H_{P_j} + Dmax * M)
 %   Q_n(M) = H_{P_n} * [ y_n/(1-y_n) - y_n^(M+1)/(1-y_n) ]
 %
 % with Dmax = max_j D_j and H_k the k-th harmonic number. The harmonic
 % weights are what distinguishes this from the ordinary geometric bound of
 % pfqn_qzgblow: a P-way fork-join subnetwork inflates its own demand by H_P
 % in the denominator and its queue length by H_P in the numerator. Setting
 % P(n) = 1 for every n recovers pfqn_qzgblow exactly.
 %
 % @par Syntax:
 % @code
 % [Q, y] = fj_qgb(D, P, M)
 % [Q, y] = fj_qgb(D, P, M, Z)
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
 % <tr><td>Q<td>Vector of bounded mean queue lengths, one per subnetwork
 % <tr><td>y<td>Vector of geometric ratios y_n(M)
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Eq. (70).
 %
 % Original: G. Casale, R. R. Muntz, G. Serazzi, "Geometric Bounds: A
 % Noniterative Analysis Technique for Closed Queueing Networks", IEEE Trans.
 % Computers 57(6), 2008.
%}
function [Q, y] = fj_qgb(D, P, M, Z)

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

% Harmonic-weighted total demand and the heaviest subnetwork demand
Dtot = sum(D .* H);
Dmax = max(D);

y = zeros(1, N);
Q = zeros(1, N);
for n = 1:N
    y(n) = D(n) * M / (Z + Dtot + Dmax * M);
    if y(n) < 1
        Q(n) = H(n) * (y(n) / (1 - y(n)) - y(n)^(M + 1) / (1 - y(n)));
    else
        % Degenerate ratio: the bound collapses onto the full population
        Q(n) = M;
    end
end

end
