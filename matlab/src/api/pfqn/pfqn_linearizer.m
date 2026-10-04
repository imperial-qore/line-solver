%{
%{
 % @file pfqn_linearizer.m
 % @brief Linearizer approximation for single-server stations.
%}
%}

%{
%{
 % @brief Linearizer approximation for single-server stations.
 % @fn pfqn_linearizer(L, N, Z, type, tol, maxiter)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param type Scheduling strategy type per station.
 % @param tol Convergence tolerance; 'cn' or NaN selects the Chandy-Neuse (1982) population-scaled termination test, see pfqn_cntol.
 % @param maxiter Maximum number of iterations.
 % @param QN0 (M x R) queue lengths that warm-start the Bard-Schweitzer initialization; empty for the default cold start.
 % @return Q Mean queue lengths.
 % @return U Utilization.
 % @return W Waiting times.
 % @return C Cycle times.
 % @return X System throughput.
 % @return totiter Total iterations performed.
%}
%}
function [Q,U,W,C,X,totiter] = pfqn_linearizer(L,N,Z,type,tol,maxiter,QN0)
% Single-server version of linearizer
if nargin<7
    QN0 = [];
end
alpha = ones(size(N));
[Q,U,W,C,X,totiter] = pfqn_egflinearizer(L,N,Z,type,tol,maxiter,alpha,QN0);
end
