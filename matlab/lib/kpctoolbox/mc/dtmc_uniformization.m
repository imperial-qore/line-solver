%{ @file dtmc_uniformization.m
 %  @brief Computes transient probabilities for a DTMC using uniformization
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes transient probabilities for a DTMC using uniformization
 %
 % @details
 % Applies the uniformization method to compute transient state probabilities for a DTMC.
 %
 % @par Syntax:
 % @code
 % [pi, kmax] = dtmc_uniformization(pi0, P)
 % [pi, kmax] = dtmc_uniformization(pi0, P, t, tol, maxiter)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>pi0<td>Initial probability distribution vector
 % <tr><td>P<td>Transition probability matrix
 % <tr><td>t<td>(Optional) Time point for transient analysis. Default: 1e4
 % <tr><td>tol<td>(Optional) Error tolerance. Default: 1e-12
 % <tr><td>maxiter<td>(Optional) Maximum number of iterations. Default: 100
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>pi<td>Probability distribution vector at time t
 % <tr><td>kmax<td>Number of iterations performed
 % </table>
%}
function [pi,kmax]=dtmc_uniformization(pi0,P,t,tol,maxiter)
if ~exist('t','var')
    t = 1e4;
end
if ~exist('tol','var')
    tol = 1e-12;
end
if ~exist('maxiter','var')
    maxiter = 100;
end
[pi,kmax] = ctmc_uniformization(pi0,ctmc_makeinfgen(P),t,tol,maxiter);
end
