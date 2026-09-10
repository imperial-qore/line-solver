%{
%{
 % @file pfqn_lekt.m
 % @brief The common corrected asymptotic expansion (LE-KT), computed on the cheaper side.
%}
%}

%{
%{
 % @brief The common corrected asymptotic expansion (LE-KT), computed on the cheaper side.
 % @fn pfqn_lekt(L, N, Z)
 % @param L Service demand matrix (MxR).
 % @param N Population vector (1xR).
 % @param Z Think time vector (1xR, default: zeros).
 % @return Gn Estimated normalizing constant.
 % @return lGn Logarithm of normalizing constant.
 % @return route 'kt' or 'le', the side the value was computed on.
%}
%}
function [Gn,lGn,route]=pfqn_lekt(L,N,Z)
% [GN,LGN,ROUTE]=PFQN_LEKT(L,N,Z)
%
% PFQN_LEKT The corrected logistic expansion and the corrected Knessl-Tier
% expansion are ONE estimator, evaluated in M-1 and in R dimensions; this
% routine computes it on whichever side is cheaper.
%
% With a think time, pfqn_ble (LE plus M units of 1-log(2*pi)/2, one per
% Laplaced station direction) and pfqn_bkt (KT minus the Stirling remainder
% s(N_r) of every Laplaced class direction) are the same function of (L,N,Z):
% their stationary points are one point in dual coordinates, xi_r =
% N_r/(Z_r + v u'L_r) being the class throughputs of the LE fixed point and
% v u_k = 1/(1-U_k) the M/M/1 factor of the KT saddle, and Sylvester's identity
% exchanges the R x R Hessian determinant for the M x M one, after which every
% 2*pi cancels on both sides. They agree to the accuracy of the two saddle-point
% solvers (~1e-7 nats; 1e-14 with polished saddles). Without a think time the
% LE branch integrates the radius exactly as Gamma(N+M) while KT Laplaces it, so
% the two differ by the constant (1-log(2*pi)/2) - r(N+M), with r the Stirling
% remainder of a Gamma direction; the common estimator is defined as the KT
% value, which is the better of the two on interior modes, and the LE side here
% carries M*(1-log(2*pi)/2) - r(N+M) rather than pfqn_ble's (M-1)*(1-log(2*pi)/2).
%
% ROUTE. The KT side solves an R-dimensional convex problem and an R x R
% determinant, the LE side an M-dimensional fixed point and an (M-1) x (M-1)
% one, so the KT side is taken when R <= M. It is also taken whenever a class
% self-loops (a single nonzero demand and no think time), since pfqn_kt extracts
% that class exactly where the logistic expansion only approximates it.
%
% Input:
% L : MxR demand matrix. L(i,r) is the demand of class-r at queue i
% N : 1xR population vector. N(r) is the number of jobs in class r
% Z : 1xR think time vector. Z(r) is the total think time of class r
%
% Output:
% Gn   : estimated normalizing constant
% lGn  : logarithm of Gn. If Gn exceeds the floating-point range, only lGn
%        will be correctly estimated.
% route: 'kt' or 'le'
%
% References:
% G. Casale. Accelerating performance inference over closed systems by
% asymptotic methods. ACM SIGMETRICS 2017.
% C. Knessl, C. Tier. Asymptotic expansions for large closed queueing networks
% with multiple job classes. IEEE Trans. Computers, 41(4):480-488, 1992.

if nargin<3 || isempty(Z)
    Z = zeros(1,numel(N));
end
[M,R] = size(L);

selfloop = false;
if R > 1
    selfloop = any(sum(L>0,1)==1 & Z(:)'==0);
end
if R <= M || selfloop
    route = 'kt';
    [Gn,lGn] = pfqn_bkt(L,N,Z);
    return
end

route = 'le';
[Gn,lGn] = pfqn_ble(L,N,Z);
if isempty(L) || isempty(N) || sum(N)==0 || sum(L(:))<1e-4
    return % pfqn_ble's degenerate branch: the delay term is exact, nothing to shift
end
if sum(Z(:)) < GlobalConstants.Zero
    % the Z=0 branch of pfqn_ble counts M-1 directions; the common estimator
    % carries M kappa - r(N+M), see the header
    eta = sum(N) + M;
    kappa = 1 - log(2*pi)/2;
    r = gammaln(eta) - (eta-0.5)*log(eta) + eta - 0.5*log(2*pi);
    lGn = lGn + kappa - r;
    Gn = exp(lGn);
end
end
