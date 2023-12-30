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
x0 = N ./ (sum(L,1) + (sum(N)-1)*max(L,[],1) + Z);
[Xopt] = fmincon(@(x) obj(x), x0, L, ones(M,1), [],[], zeros(1,R), [],[],optimopt);
Xasy = Xopt .* (sum(L,1) + (sum(N)-1)*max(L,[],1)) ./ (sum(L,1) + (sum(N)-1)*max(L,[],1) + Z); 
lG = -sum(N.*log(Xasy));
G = exp(lG);
end
