%{ @file fj_rmax.m
 %  @brief Compute maximum response time R_K^max(rho) for K M/M/1 queues
 %
 %  @author LINE Development Team
%}

%{
 % @brief Compute maximum response time R_K^max(rho) for K M/M/1 queues
 %
 % @details
 % Computes the expected value of the maximum of K response times from
 % independent M/M/1 queues. This serves as an upper bound to the K-way
 % Fork-Join response time R_K^{F/J}(rho).
 %
 % R_K^max(rho) = H_K * R(rho) = H_K / (mu - lambda)
 %
 % where H_K is the K-th Harmonic number and R(rho) is the M/M/1 mean
 % response time.
 %
 % This formula holds because for M/M/1 queues, the response time
 % distribution is exponential with rate (mu - lambda), and the expected
 % maximum of K i.i.d. exponential random variables with rate theta is
 % H_K / theta.
 %
 % @par Syntax:
 % @code
 % Rmax = fj_rmax(K, lambda, mu)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of parallel servers (positive integer)
 % <tr><td>lambda<td>Arrival rate
 % <tr><td>mu<td>Service rate (mu > lambda for stability)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Rmax<td>Expected maximum response time R_K^max(rho) = H_K * R(rho)
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Page 17:8, Eq. (1), and surrounding text.
%}
function Rmax = fj_rmax(K, lambda, mu)

if K < 1
    line_error(mfilename, 'K must be a positive integer. Got K=%d.', K);
end

% Compute harmonic number H_K
H_K = fj_harmonic(K);

% Compute M/M/1 mean response time
R_rho = qsys_mm1(lambda, mu);

% Maximum response time: R_K^max(rho) = H_K * R(rho)
Rmax = H_K * R_rho;

end
