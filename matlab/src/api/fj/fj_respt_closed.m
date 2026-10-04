%{ @file fj_respt_closed.m
 %  @brief Varki bound on the residence time of a closed fork-join subnetwork
 %
 %  @author LINE Development Team
%}

%{
 % @brief Varki bound on the residence time of a closed fork-join subnetwork
 %
 % @details
 % Bounds the mean residence time of a K-way fork-join subnetwork of
 % homogeneous exponential servers with mean service time x. Writing A for the
 % mean number of jobs an arriving job finds at the subnetwork,
 %
 %   R_{P_K}(M) <= x * [ H_K + A ].
 %
 % In a closed network the arrival theorem supplies A = Q(M-1); when the
 % network consists of the parallel subsystem alone, every one of the other
 % M-1 jobs is necessarily inside it, so A = M-1 and
 %
 %   R_{P_K}(M) <= x * [ H_K + M - 1 ],
 %
 % which holds with equality for K = 2 and is therefore exact there.
 %
 % @par Syntax:
 % @code
 % R = fj_respt_closed(K, x, M)
 % [R, exact] = fj_respt_closed(K, x, M, A)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of parallel servers (positive integer)
 % <tr><td>x<td>Mean service time of each server
 % <tr><td>M<td>Number of jobs circulating in the closed network
 % <tr><td>A<td>Mean queue length seen on arrival (optional, default M-1)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>R<td>Upper bound on the mean residence time of the subnetwork
 % <tr><td>exact<td>True when the bound is known to be tight, i.e. K = 2 and A = M-1
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Eq. (67)
 % and Theorem 4.1 on page 17:47.
 %
 % Original: E. Varki, "Response Time Analysis of Parallel Computer and
 % Storage Systems", IEEE Trans. Parallel Distrib. Syst. 12(11), 2001.
%}
function [R, exact] = fj_respt_closed(K, x, M, A)

if K < 1 || K ~= round(K)
    line_error(mfilename, 'K must be a positive integer. Got K=%g.', K);
end
if x <= 0
    line_error(mfilename, 'Mean service time x must be positive. Got x=%g.', x);
end
if M < 1 || M ~= round(M)
    line_error(mfilename, 'M must be a positive integer. Got M=%g.', M);
end

isolated = (nargin < 4) || isempty(A);
if isolated
    % Closed network made of the parallel subsystem alone: every other job is
    % inside it, so the arrival theorem gives A = M-1 (Theorem 4.1)
    A = M - 1;
end
if A < 0
    line_error(mfilename, 'Arrival-instant queue length A must be non-negative. Got A=%g.', A);
end

H_K = fj_harmonic(K);
R = x * (H_K + A);

% Theorem 4.1 holds with equality at K = 2 for the isolated parallel subsystem
exact = isolated && (K == 2);

end
