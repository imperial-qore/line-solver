%{ @file fj_respt_nt.m
 %  @brief Nelson-Tantawi approximation for K-way F/J response time
 %
 %  @author LINE Development Team
%}

%{
 % @brief Nelson-Tantawi approximation for K-way F/J response time
 %
 % @details
 % Computes an approximate mean response time for a K-way Fork-Join
 % queueing system with Poisson arrivals and exponential service times.
 % Valid for 2 <= K <= 32.
 %
 % R_K^{F/J}(rho) ≈ [H_K/H_2 + (1 - H_K/H_2) * 4*rho/11] * (1.5 - rho/8) / (mu - lambda)
 %
 % The approximation is based on the exact 2-way solution and a scaling
 % approximation that upper and lower bounds increase at the same rate.
 % The coefficient alpha(rho) ≈ 4*rho/11 was obtained from simulation.
 %
 % @par Syntax:
 % @code
 % R = fj_respt_nt(K, lambda, mu)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of parallel servers (2 <= K <= 32)
 % <tr><td>lambda<td>Arrival rate
 % <tr><td>mu<td>Service rate (mu > lambda for stability)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>R<td>Approximate K-way F/J response time R_K^{F/J}(rho)
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Table I, Eq. (3) on page 17:10.
 %
 % Original: R. Nelson and A.N. Tantawi, "Approximate Analysis of Fork/Join
 % Synchronization in Parallel Queues", IEEE Trans. Computers, 37(6), 1988.
%}
function R = fj_respt_nt(K, lambda, mu)

if K < 2
    line_error(mfilename, 'Nelson-Tantawi approximation requires K >= 2. Got K=%d.', K);
end

if K > 32
    line_warning(mfilename, 'Nelson-Tantawi approximation was validated for K <= 32. Using K=%d may have reduced accuracy.', K);
end

% Compute utilization
rho = lambda / mu;

if rho >= 1
    line_error(mfilename, 'System is unstable: rho = lambda/mu = %.4f >= 1. Require lambda < mu.', rho);
end

% Compute harmonic numbers
H_K = fj_harmonic(K);
H_2 = fj_harmonic(2);  % H_2 = 1.5

% Nelson-Tantawi approximation (Eq. 3):
% R_K^{F/J}(rho) ≈ [H_K/H_2 + (1 - H_K/H_2) * 4*rho/11] * (1.5 - rho/8) * (mu - lambda)^{-1}
S_K = H_K / H_2 + (1 - H_K / H_2) * (4 * rho / 11);
R_2_factor = 1.5 - rho / 8;  % This is (H_2 - rho/8)
R_rho = 1 / (mu - lambda);

R = S_K * R_2_factor * R_rho;

end
