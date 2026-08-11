%{ @file fj_sm_tput.m
 %  @brief Maximum throughput for Split-Merge queueing system
 %
 %  @author LINE Development Team
%}

%{
 % @brief Maximum throughput for Split-Merge queueing system
 %
 % @details
 % Computes the maximum throughput for a K-way Split-Merge (SM) queueing
 % system with exponential service times.
 %
 % In an SM system, all tasks of a request must complete before the next
 % request can be issued. The maximum throughput is:
 %
 %   lambda_K^SM = 1 / X_K^max = mu / H_K
 %
 % where X_K^max = H_K / mu is the expected maximum service time.
 %
 % For comparison, the maximum throughput of a K-way F/J system is
 % lambda_K^FJ = mu, which exceeds lambda_K^SM by a factor H_K.
 %
 % @par Syntax:
 % @code
 % lambda_max = fj_sm_tput(K, mu)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of parallel servers (positive integer)
 % <tr><td>mu<td>Service rate at each server
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lambda_max<td>Maximum throughput lambda_K^SM = mu / H_K
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Page 17:2.
%}
function lambda_max = fj_sm_tput(K, mu)

if K < 1
    line_error(mfilename, 'K must be a positive integer. Got K=%d.', K);
end

if mu <= 0
    line_error(mfilename, 'Service rate mu must be positive. Got mu=%.4f.', mu);
end

% Compute harmonic number H_K
H_K = fj_harmonic(K);

% Maximum throughput: lambda_K^SM = mu / H_K = 1 / X_K^max
lambda_max = mu / H_K;

end
