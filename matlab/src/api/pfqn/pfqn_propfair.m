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
obj = @(X) -sum(N.*log(X+1e-6));
x0 = N ./ (N.*sum(L,1)+Z);
[Xopt] = fmincon(@(x) obj(x), x0, L, ones(M,1), [],[], zeros(1,R), [],[],optimopt);
Xasy = Xopt .* N ./ (N + Z.*Xopt); % add think time. If embedded in the optimization program it goes worse.
lG = -sum(N.*log(Xasy));
G = exp(lG);
end
