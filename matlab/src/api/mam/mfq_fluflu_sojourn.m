%{ @file mfq_fluflu_sojourn.m
 %  @brief Sojourn-time distribution of a fluid/fluid queue
 %
 %  @author LINE Development Team
%}

%{
 % @brief Sojourn-time distribution of a fluid queue with fluid-modulated service.
 %
 % @details
 % Thin LINE wrapper around FluFluSTD. Returns the sojourn-time distribution of
 % a fluid queue in which both the arrival and the service processes are
 % Markov-modulated fluid flows, as a matrix-exponential (ME) or phase-type
 % (PH) representation (alpha, A).
 %
 % @par Syntax:
 % @code
 % [alpha,A] = mfq_fluflu_sojourn(Qin,Rin,Qout,Rout,srv0stop)
 % [alpha,A] = mfq_fluflu_sojourn(Qin,Rin,Qout,Rout,srv0stop,transToPH)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Qin<td>(Na,Na) generator of the arrival-modulating chain
 % <tr><td>Rin<td>(Na,Na) diagonal arrival fluid-rate matrix
 % <tr><td>Qout<td>(Ns,Ns) generator of the service-modulating chain
 % <tr><td>Rout<td>(Ns,Ns) diagonal service fluid-rate matrix
 % <tr><td>srv0stop<td>true if service stops when the server fluid level hits zero
 % <tr><td>transToPH<td>(Optional) true to return a PH (else ME) representation
 % </table>
 %
 % @par Returns:
 % alpha, A: ME/PH representation of the sojourn-time distribution.
%}
function [alpha,A] = mfq_fluflu_sojourn(Qin,Rin,Qout,Rout,srv0stop,transToPH)
if ~exist('transToPH','var'), transToPH = false; end
[alpha,A] = FluFluSTD(Qin,Rin,Qout,Rout,srv0stop,transToPH);
end
