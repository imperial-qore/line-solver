%{
%{
 % @file pfqn_gflinearizer.m
 % @brief Generalized fixed-point Linearizer with uniform scaling exponent.
%}
%}

%{
%{
 % @brief Generalized fixed-point Linearizer with uniform scaling exponent.
 % @fn pfqn_gflinearizer(L, N, Z, type, tol, maxiter, alpha)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param type Scheduling strategy type per station.
 % @param tol Convergence tolerance; 'cn' or NaN selects the Chandy-Neuse (1982) population-scaled termination test, see pfqn_cntol.
 % @param maxiter Maximum number of iterations.
 % @param alpha Uniform scaling exponent for all classes.
 % @param QN0 (M x R) queue lengths that warm-start the Bard-Schweitzer initialization; empty for the default cold start.
 % @return Q Mean queue lengths.
 % @return U Utilization.
 % @return W Waiting times.
 % @return C Cycle times.
 % @return X System throughput.
 % @return totiter Total iterations performed.
%}
%}
function [Q,U,W,C,X,totiter] = pfqn_gflinearizer(L,N,Z,type,tol,maxiter,alpha,QN0)
if nargin<8
    QN0 = [];
end
R = size(L,2);
[Q,U,W,C,X,totiter] = pfqn_egflinearizer(L,N,Z,type,tol,maxiter,alpha*ones(1,R),QN0);
end


