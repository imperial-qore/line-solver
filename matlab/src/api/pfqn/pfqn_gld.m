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

% SYMBOLIC INPUT. Every branch below that inspects the VALUES of L or mu, the
% zero-demand guard, the L>0 mask and the load-independence scan, is a
% comparison that a sym cannot resolve: it raises "Unable to prove the
% statement" or "Conversion to logical from sym is not possible" rather than
% returning a logical. Each is therefore replaced by a numeric test on N, which
% is always concrete, or skipped in favour of the recursion, which is +, * and /
% throughout and so is field-agnostic. The single-class kernel is routed to
% pfqn_gldsingle rather than pfqn_lldsingle for the same reason: the LLD
% threshold is found by COMPARING rates, and that is undecidable on a sym.
isSym = isa(L,'sym') || isa(mu,'sym');

% The mu default has to precede the M==1 and R==1 branches below, both of which
% READ mu: left further down, next to the options default, a two-argument call
% reached them with mu undefined and died on "Not enough input arguments"
% instead of running the load-independent model the default describes.
if nargin<3 || isempty(mu)
    mu=ones(M,sum(N));
end

if M==1
    % A CLASS WITH JOBS AND NO DEMAND AT THE ONLY STATION MAKES THE CONSTANT
    % ZERO. Its factor is L_r^N_r = 0, so the whole product vanishes; dropping
    % the class from the sum below instead answers with the constant of a
    % DIFFERENT model, the one without it. The recursion at the foot of this
    % file reaches this base case with the full population every time it peels
    % a station, so the error surfaces on any load-dependent model carrying a
    % zero demand. See _kb/07-cross-language-parity.md.
    % abs(), not L>0: the modulus is the right notion of "no demand" for the
    % complex demands infradius_h hands this routine, and it is the same test
    % on a real one
    if isSym
        % abs(L)==0 and L>0 are undecidable here, so the classes are selected by
        % the numeric test N>0 instead: one with no jobs contributes nothing
        % whatever its demand is, and one with jobs is assumed to have a
        % nonzero symbolic demand, there being no way to prove otherwise.
        act = N>0;
        % The multinomial is formed EXACTLY, as a ratio of factorials, and the
        % product is taken directly rather than through exp(factln(...)):
        % factln returns a DOUBLE, which sym would carry as a rational
        % approximation of a float, leaving exp(6243314768165359/45035996...)
        % where the integer 4 belongs and making the expression unusable.
        G = factorial(sym(sum(N))) / prod(factorial(sym(N))) ...
            * prod(L(1,act).^N(act)) / prod(mu(1,1:sum(N)));
        lG = log(G);
        return
    end
    if any(N>0 & abs(L)==0)
        lG = -Inf; G = 0;
        return
    end
    lG = factln(sum(N)) - sum(factln(N)) + N(L>0)*log(L(L>0))' - sum(log(mu(1,1:sum(N))));
    G = exp(lG);
    return
end

if R==1
    if isSym
        [lG,G] = pfqn_gldsingle(L,N,mu);
    else
        [lG,G] = pfqn_lldsingle(L,N,mu);
    end
    return
end

if isempty(L)
    G = 0; lG = -Inf; return
end

if nargin<4
    options = SolverNC.defaultOptions;
end

isLoadDep = false;
isInfServer = [];
if isSym
    % min(row)==1 and all(row==1:sum(N)) are symbolic comparisons with no
    % logical value, so the load-independent shortcut is skipped: the recursion
    % at the foot of this file compares nothing and handles the model as given
    isLoadDep = true;
    isInfServer = false(1,M);
else
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
