%{
%{
 % @file pfqn_ncldmx.m
 % @brief Normalizing constant for mixed open/closed networks with limited load dependence.
%}
%}

%{
%{
 % @brief Normalizing constant for mixed open/closed networks with limited load dependence.
 %
 % The closed-conditional normalizing constant of a mixed limited load-dependent
 % (LLD) network equals a purely closed load-dependent normalizing constant in
 % which every queueing station i carries the Bruell-Balbo-Afshari effective
 % capacity rate mu_i^eff(n) = 1/EC_i(n), where EC is returned by
 % pfqn_ldmx_ec and folds the open classes into the closed subnetwork. The
 % open classes contribute the separable prefactor lGopen = sum_i log E_i(0),
 % which reduces to -sum_i log(1-rho_i) in the load-independent limit.
 %
 % Mean closed metrics follow from the standard normalizing-constant ratios,
 % e.g. X_r = G(N-e_r)/G(N), matching the exact pfqn_mvaldmx solver.
 %
 % @fn pfqn_ncldmx(lambda, D, N, Z, mu, S, varargin)
 % @param lambda Arrival rate vector (0 on closed classes).
 % @param D Service demand matrix (MxR).
 % @param N Population vector (Inf on open classes).
 % @param Z Think time vector (closed classes).
 % @param mu Load-dependent rate matrix (Mx>=sum(N_closed)).
 % @param S Number of servers per station (currently informational).
 % @param varargin Optional solver parameters forwarded to pfqn_ncld.
 % @return lG Logarithm of the closed-conditional normalizing constant.
 % @return G Closed-conditional normalizing constant (exp(lG)).
 % @return lGopen Logarithm of the open-class normalizing prefactor sum_i log E_i(0).
%}
%}
function [lG,G,lGopen] = pfqn_ncldmx(lambda,D,N,Z,mu,S,varargin)
% [LG,G,LGOPEN] = PFQN_NCLDMX(LAMBDA,D,N,Z,MU,S,VARARGIN)

[M,R] = size(D);
if nargin<5 || isempty(mu)
    mu = ones(M,max(1,sum(N(isfinite(N)))));
end
if nargin<6 || isempty(S)
    S = ones(M,1); %#ok<NASGU> % kept for signature parity with pfqn_mvaldmx
end
if nargin<4 || isempty(Z)
    Z = zeros(1,R);
end

openClasses = find(isinf(N));
closedClasses = setdiff(1:R, openClasses);

if any(N(intersect(find(lambda),closedClasses))>0)
    line_error(mfilename,'Arrival rate cannot be specified on closed classes.');
end

Nc = N(closedClasses);
Kc = sum(Nc);

% open-class normalizing prefactor sum_i log E_i(0); needs the effective
% capacity terms even when there are no closed jobs.
mup = mu;
if size(mup,2) < max(1,Kc)
    mup = [mup, repmat(mup(:,end), 1, max(1,Kc)-size(mup,2))];
end
mup = [mup, mup(:,end)]; % one extra column as in pfqn_mvaldmx
lambdao = zeros(1,R); lambdao(openClasses) = lambda(openClasses);
[EC,E] = pfqn_ldmx_ec(lambdao, D, mup);
lGopen = sum(log(E(1:M,1)));

% closed-conditional normalizing constant
if Kc==0
    lG = 0;
    G = 1;
    return
end
Dc = D(:,closedClasses);
Zc = Z(closedClasses);
muEff = 1 ./ EC(:,1:Kc);
[lG,G] = pfqn_ncld(Dc, Nc, Zc, muEff, varargin{:});
end
