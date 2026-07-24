%{
%{
 % @file pfqn_gld.m
 % @brief Exact normalizing constant for load-dependent queueing networks.
%}
%}

%{
%{
 % @brief Exact normalizing constant for load-dependent queueing networks.
 % @fn pfqn_gld(L, N, mu, options)
 % @param L Service demand matrix (MxR).
 % @param N Population vector (1xR).
 % @param mu Load-dependent rate matrix (Mx sum(N)).
 % @param options Solver options.
 % @return G Normalizing constant.
 % @return lG Logarithm of normalizing constant.
%}
%}
function [G,lG]=pfqn_gld(L,N,mu,options)
% [G,LG]=PFQN_GLD(L,N,MU,OPTIONS)

% G=pfqn_gld(L,N,mu)
% mu: MxN matrix of load-dependent rates
[M,R]=size(L);
lambda = zeros(1,R);
if M==1       
    lG = factln(sum(N)) - sum(factln(N)) + N(L>0)*log(L(L>0))' - sum(log(mu(1,1:sum(N))));
    G = exp(lG);
    return
end

if R==1
    [lG,G] = pfqn_gldsingle(L,N,mu);
    return
end

if isempty(L)
    G = 0; lG = -Inf; return
end

if nargin==2
    mu=ones(M,sum(N));
end

if nargin<4
    options = SolverNC.defaultOptions;
end

isLoadDep = false;
isInfServer = [];
for ist=1:M
    if min(mu(ist,1:sum(N))) == 1 & max(mu(ist,1:sum(N))) == 1
        isInfServer(ist) = false;
        continue; % this is a LI station
    elseif all(mu(ist,1:sum(N)) == 1:sum(N))
        isInfServer(ist) = true;
        continue; % this is a infinite server station
    else
        isInfServer(ist) = false;
        isLoadDep = true;
    end
end

if ~isLoadDep
    % if load-independent model then use faster pfqn_gmva solver
    Lli = L(find(~isInfServer),:);
    if isempty(Lli)
        Lli = 0*N;
    end
    Zli = L(find(isInfServer),:);
    if isempty(Zli)
        Zli = 0*N;
    end
    options.method='exact';    
    lG = pfqn_nc(lambda,Lli, N, sum(Zli,1), options);
    G = exp(lG);
    return
end

G=0;
if M==0 
	G=0; 
	lG=log(G); 
	return; 
end
if sum(N==zeros(1,R))==R 
	G=1; 
	lG=log(G); 
	return; 
end

if R==1
    G=pfqn_gldsingle(L,N,mu);
    lG=log(G);
    return
end

G=G + pfqn_gld(L(1:(M-1),:),N,mu(1:(M-1),:),options);
for r=1:R
    if N(r)>0
        if R>1
            N_1 = oner(N,r);
        else
            N_1 = N-1;
        end
        G = G + (L(M,r)/mu(M,1))*pfqn_gld(L,N_1,pfqn_mushift(mu,M),options);
    end
end
lG=log(G);
return
end
