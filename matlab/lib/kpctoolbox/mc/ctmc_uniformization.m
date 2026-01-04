%{ @file ctmc_uniformization.m
 %  @brief Computes transient probabilities using uniformization
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes transient probabilities using uniformization
 %
 % @details
 % Applies the uniformization method to compute transient state probabilities.
 %
 % @par Syntax:
 % @code
 % [pi, kmax] = ctmc_uniformization(pi0, Q, t)
 % [pi, kmax] = ctmc_uniformization(pi0, Q, t, tol, maxiter)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>pi0<td>Initial probability distribution vector
 % <tr><td>Q<td>Infinitesimal generator matrix
 % <tr><td>t<td>Time point for transient analysis
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
function [pi,kmax]=ctmc_uniformization(pi0,Q,t,tol,maxiter)
if ~exist('tol','var')
    tol = 1e-12;
end
if ~exist('maxiter','var')
    maxiter = 100;
end
q=1.1*max(abs(diag(Q)));
Qs=speye(size(Q))+sparse(Q)/q;
k=0;
s=1;
r=1;
iter=0;
kmax=1;
while iter<maxiter
    iter=iter+1;
    k=k+1;
    r=r*(q*t)/k;
    s=s+r;
    if (1-exp(-q*t)*s)<=tol
        kmax=k;
        break;
    end
end

pi=pi0*(exp(-q*t));
P=pi0;
ri=exp(-q*t);
for j=1:kmax
    P=P*Qs;
    ri=ri*(q*t/j);
    pi=pi+ri*P;
end
end
