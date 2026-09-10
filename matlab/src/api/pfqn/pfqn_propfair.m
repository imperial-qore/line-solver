%{
%{
 % @file pfqn_propfair.m
 % @brief Proportionally fair allocation approximation (Walton 2009).
%}
%}

%{
%{
 % @brief Proportionally fair allocation approximation (Walton 2009).
 % @fn pfqn_propfair(L, N, Z)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @return G Normalizing constant estimate.
 % @return lG Logarithm of normalizing constant.
 % @return Xasy Asymptotic throughput vector.
%}
%}
function [G,lG,Xasy] = pfqn_propfair(L,N,Z)
% [G,LOG] = PFQN_PROPFAIR(L,N,Z)

% Proportionally fair allocation
%
% Estimate the normalizing constant using a convex optimization program
% that is asymptotically exact in models with single-server PS queues only.
% The underlying optimization program is convex.
% The script implements a heuristic to estimate the solution in the
% presence of delay stations.
%
% Schweitzer, P. J. (1979). Approximate analysis of multiclass closed networks of
% queues. In Proceedings of the International Conference on Stochastic Control
% and Optimization. Free Univ., Amsterdam.
%
% Walton, Proportional fairness and its relationship with multi-class
% queueing networks, 2009. 

[M,R]=size(L);
optimopt = optimoptions(@fmincon,'MaxFunctionEvaluations',1e6,'Display','none');
obj = @(X) -sum((N-X.*Z).*log(X+1e-6));
[Xasy] = fmincon(@(x) obj(x), zeros(1,R), L, ones(M,1), [],[], zeros(1,R), [],[],optimopt);
lG = 0;
for r=1:R
    if Z(r)>0
        lG = lG + Xasy(r)*Z(r).*log(Z(r));
    end
end
lG = sum((N-Xasy.*Z).*log(1./Xasy)) - sum(factln(Xasy.*Z));
G = exp(lG);
end
