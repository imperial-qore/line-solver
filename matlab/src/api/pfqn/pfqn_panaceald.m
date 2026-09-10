%{
%{
 % @file pfqn_panaceald.m
 % @brief PANACEA asymptotic expansion for load-dependent closed networks.
%}
%}

%{
%{
 % @brief PANACEA asymptotic expansion for load-dependent closed networks.
 % @fn pfqn_panaceald(L, N, Z, mu, terms)
 % @param L Service demand matrix (MxR).
 % @param N Population vector (1xR).
 % @param Z Think time matrix (DxR), summed over rows.
 % @param mu Load-dependent rate matrix (Mx sum(N)).
 % @param terms Number of terms in the normal-usage asymptotic series
 %        (1, 2, or 3; default 3), as in pfqn_panacea.
 % @return Gn Normalizing constant.
 % @return lGn Logarithm of normalizing constant.
%}
%}
function [Gn,lGn]=pfqn_panaceald(L,N,Z,mu,terms)
% [GN,LGN]=PFQN_PANACEALD(L,N,Z,MU,TERMS)

% Mitra-McKenna (JACM 33(3):568-592, 1986) load-dependent PANACEA: the
% expansion coefficients A_n are linear combinations of partition functions
% of a pseudonetwork whose load dependence is the phi(n) transform of the
% original {f(n)}. See _kb/03-api-layer.md (pfqn/ family, panaceald).
[M,R]=size(L);
if nargin<3 || isempty(Z)
    Z = zeros(1,R);
end
if nargin<5 || isempty(terms)
    terms = 3;
end
if ~isscalar(terms) || terms<1 || terms>3 || terms~=round(terms)
    line_error(mfilename,'The terms parameter must be 1, 2, or 3 (higher-order coefficients are not implemented).');
end
Gn = NaN; lGn = NaN;
Nt = round(sum(N)); % total population, also the truncation point of {f(n)}
if nargin<4 || isempty(mu)
    mu = ones(M,max(1,Nt));
end
if Nt==0
    Gn = 1; lGn = 0;
    return
end
if size(mu,2) < Nt
    mu = [mu, repmat(mu(:,end),1,Nt-size(mu,2))];
end

% Type-3 (infinite-server) rows are absent from the pseudonetwork and enter
% only through rho_j0; solver_ncld encodes them as mu(i,n)=n rows of L.
isIS = false(M,1);
for i=1:M
    if all(abs(mu(i,1:Nt) - (1:Nt)) < GlobalConstants.FineTol)
        isIS(i) = true;
    end
end
Ztot = sum(Z,1) + sum(L(isIS,:),1);
Lq = L(~isIS,:);
muq = mu(~isIS,1:Nt);
Mq = size(Lq,1);

if any(N>0 & Ztot<=0)
    % no IS center on the route of a populated class: the expansion parameter
    % rho_j0 is undefined and PANACEA does not apply
    return
end
if Mq==0
    lGn = -sum(factln(N)) + sum(xlogy(N,Ztot));
    Gn = exp(lGn);
    return
end
if any(muq(:)<=0) || any(~isfinite(muq(:)))
    return
end

r = zeros(Mq,R);
for j=1:R
    if Ztot(j)>0
        r(:,j) = Lq(:,j)/Ztot(j);
    end
end
lambda = r*N(:);        % offered load per queueing center, sum_j K_j e_ji/rho_j0
muK = muq(:,Nt);        % saturation rates, c_i = 1/mu_i(Ntot)
alpha = 1 - lambda./muK;
if min(alpha)<=0
    % model is not in normal usage: the {phi(n)} series diverges
    return
end

% log-partial products log prod_{k=1}^{s} mu_i(k), s=0..Nt
lPi = [zeros(Mq,1), cumsum(log(muq),2)];

nmax = 2*(terms-1);     % largest pseudonetwork per-class population needed
lpsi = zeros(Mq,nmax+1);
for i=1:Mq
    for n=0:nmax
        lpsi(i,n+1) = logpsi(n,lambda(i),lPi(i,:),muK(i),alpha(i),Nt);
    end
end

% load dependence of the pseudonetwork centers: psi_i(n) = psi_i(0) n! /
% prod_{k=1}^{n} mups_i(k)
mups = ones(Mq,max(1,nmax));
for i=1:Mq
    for n=1:nmax
        mups(i,n) = exp(log(n) + lpsi(i,n) - lpsi(i,n+1));
    end
end

% Expansion coefficients (5.4). The large parameter N cancels identically
% between beta_j=K_j/N, Gamma=N*r and the 1/N^n scaling, so the demands are
% taken as r and beta as N.
A = zeros(1,3);
A(1) = 1;
if terms>=2
    for j=1:R
        m = zeros(1,R); m(j)=2;
        A(2) = A(2) - N(j) * pseudonet(r,m,mups);
    end
end
if terms>=3
    for j=1:R
        m = zeros(1,R); m(j)=3;
        A(3) = A(3) + 2 * N(j) * pseudonet(r,m,mups);
        m = zeros(1,R); m(j)=4;
        A(3) = A(3) + 3 * N(j)^2 * pseudonet(r,m,mups);
        for k=1:R
            if k~=j
                m = zeros(1,R); m(j)=2; m(k)=2;
                A(3) = A(3) + 0.5 * N(j) * N(k) * pseudonet(r,m,mups);
            end
        end
    end
end
I = sum(A(1:terms));
if I<=0
    return
end

lGn = -sum(factln(N)) + sum(xlogy(N,Ztot)) + sum(lpsi(:,1)) + log(I);
Gn = exp(lGn);
if ~isfinite(lGn)
    Gn = NaN;
    lGn = NaN;
end
end

function lp = logpsi(n,lambda,lPirow,muK,alpha,K)
% LP = LOGPSI(N,LAMBDA,LPIROW,MUK,ALPHA,K)

% log of psi(n) = sum_{s>=n} [s!/(s-n)!] lambda^(s-n) / prod_{k=1}^{s} mu(k),
% the mu-free part of the phi(n) transform in eq. (3.7)-(3.8a). The series is
% split into the exact head s<=K and a geometric tail summed in closed form
% via the Vandermonde identity, all terms positive.
t = -Inf(1, max(0,K-n+1) + n+1);
c = 0;
for s=n:K
    c = c + 1;
    t(c) = factln(s) - factln(s-n) + xlogy(s-n,lambda) - lPirow(s+1);
end
T = max(n,K+1);
for i=0:n
    c = c + 1;
    t(c) = factln(n) + factln(T) - factln(n-i) - factln(T-n+i) ...
        + xlogy(T+i-n,lambda) + (K-T-i)*log(muK) ...
        - (i+1)*log(alpha) - lPirow(K+1);
end
t = t(1:c);
tmax = max(t);
if ~isfinite(tmax)
    lp = tmax;
else
    lp = tmax + log(sum(exp(t-tmax)));
end
end

function G = pseudonet(gam,k,mups)
% G = PSEUDONET(GAM,K,MUPS)

% Partition function of the pseudonetwork at population k, normalized so that
% G(0)=1. Populations are at most 2*(terms-1), so a direct load-dependent
% convolution over the mixed-radix population lattice is used.
nz = find(k>0);
gam = gam(:,nz);
k = k(nz);
[Mq,Rp] = size(gam);
sizes = k+1;
npop = prod(sizes);
sterm = zeros(Mq,npop);
for i=1:Mq
    for jdx=1:npop
        m = idx2vec(jdx,sizes,Rp);
        sm = sum(m);
        v = factln(sm);
        for rr=1:Rp
            if m(rr)>0
                if gam(i,rr)<=0
                    v = -Inf;
                    break
                end
                v = v + m(rr)*log(gam(i,rr)) - factln(m(rr));
            end
        end
        if isfinite(v) && sm>0
            v = v - sum(log(mups(i,1:sm)));
        end
        sterm(i,jdx) = exp(v);
    end
end
g = zeros(Mq+1,npop);
g(1,1) = 1;
for i=1:Mq
    for idx=1:npop
        n = idx2vec(idx,sizes,Rp);
        acc = 0;
        for jdx=1:npop
            m = idx2vec(jdx,sizes,Rp);
            if all(m<=n)
                acc = acc + g(i,vec2idx(n-m,sizes,Rp)) * sterm(i,jdx);
            end
        end
        g(i+1,idx) = acc;
    end
end
G = g(Mq+1,npop);
end

function v = idx2vec(idx,sizes,Rp)
% V = IDX2VEC(IDX,SIZES,RP)

v = zeros(1,Rp);
t = idx-1;
for rr=1:Rp
    v(rr) = mod(t,sizes(rr));
    t = floor(t/sizes(rr));
end
end

function idx = vec2idx(v,sizes,Rp)
% IDX = VEC2IDX(V,SIZES,RP)

idx = 1;
mult = 1;
for rr=1:Rp
    idx = idx + mult*v(rr);
    mult = mult*sizes(rr);
end
end

function y = xlogy(e,x)
% Y = XLOGY(E,X)

% e.*log(x) with the convention 0*log(0)=0
y = zeros(size(e));
nz = e~=0;
y(nz) = e(nz).*log(x(nz));
end
