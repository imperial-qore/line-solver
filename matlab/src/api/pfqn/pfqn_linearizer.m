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
 % @param tol Convergence tolerance.
 % @param maxiter Maximum number of iterations.
 % @return Q Mean queue lengths.
 % @return U Utilization.
 % @return W Waiting times.
 % @return C Cycle times.
 % @return X System throughput.
 % @return totiter Total iterations performed.
%}
%}
function [Q,U,W,C,X,totiter] = pfqn_linearizer(L,N,Z,type,tol,maxiter)
% Single-server version of linearizer
alpha = ones(size(N));
[Q,U,W,C,X,totiter] = pfqn_egflinearizer(L,N,Z,type,tol,maxiter,alpha);
end
