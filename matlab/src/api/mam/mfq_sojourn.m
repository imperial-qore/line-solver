%{ @file mfq_sojourn.m
 %  @brief Sojourn-time distribution of a Markov-modulated fluid queue
 %
 %  @author LINE Development Team
%}

%{
 % @brief Sojourn-time distribution of a fluid queue as an ME/PH representation.
 %
 % @details
 % Thin LINE wrapper around FluidQueueSTD. Returns the distribution of the time
 % a fluid drop spends in a Markov-modulated fluid queue with input rate matrix
 % Rin and output (service) rate matrix Rout, modulated by generator Q, as a
 % matrix-exponential (ME) or phase-type (PH) representation (alpha, A).
 %
 % @par Syntax:
 % @code
 % [alpha,A] = mfq_sojourn(Q,Rin,Rout)
 % [alpha,A] = mfq_sojourn(Q,Rin,Rout,Q0,transToPH)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Q<td>(N,N) generator of the background Markov chain
 % <tr><td>Rin<td>(N,N) diagonal input fluid-rate matrix
 % <tr><td>Rout<td>(N,N) diagonal output (service) fluid-rate matrix
 % <tr><td>Q0<td>(Optional) level-0 generator (default: Q)
 % <tr><td>transToPH<td>(Optional) true to return a PH (else ME) representation
 % </table>
 %
 % @par Returns:
 % alpha, A: ME/PH representation of the sojourn-time distribution.
%}
function [alpha,A] = mfq_sojourn(Q,Rin,Rout,Q0,transToPH)
if ~exist('Q0','var'), Q0 = []; end
if ~exist('transToPH','var'), transToPH = false; end
[alpha,A] = FluidQueueSTD(Q,Rin,Rout,Q0,transToPH);
end
