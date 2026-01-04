%{
%{
 % @file pfqn_rd.m
 % @brief Reduced Decomposition (RD) method for load-dependent networks.
%}
%}

%{
%{
 % @brief Reduced Decomposition (RD) method for load-dependent networks.
 % @fn pfqn_rd(L, N, Z, mu, options)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param mu Load-dependent rate matrix.
 % @param options Solver options.
 % @return lGN Logarithm of normalizing constant.
 % @return Cgamma Gamma correction factor.
%}
%}
function [lGN,Cgamma] = pfqn_rd(L,N,Z,mu,options)
[M,R]=size(L);
lambda=zeros(1,R);
if sum(N)<0
    lGN=-Inf;
    return
end
if nargin<5
    options=SolverNC.defaultOptions;
end
% L
% N
% Z
% mu
for ist=1:M
    if all(mu(ist,:)==mu(ist,1)) % LI station
        L(ist,:) = L(ist,:) / mu(ist,1);
        mu(ist,:) = 1;
        %isLI(i) = true;
    end
end
if sum(N)==0
    lGN=0;
    return
end
gamma = ones(M,ceil(sum(N)));
mu = mu(:,1:sum(N));
%mu(mu==0)=Inf;
mu(isnan(mu))=Inf;
s = zeros(M,1);
for ist=1:M
    if isfinite(mu(ist,end))
        s(ist) = min(find(abs(mu(ist,:)-mu(ist,end))<options.tol));
        if s(ist)==0
            s(ist) = sum(N);
        end
    else
        s(ist) = sum(N);
    end
end
isDelay = false(M,1);
isLI = false(M,1);
y = L;

for ist=1:M
    if isinf(mu(ist,s(ist)))
        lastfinite=max(find(isfinite(mu(ist,:))));
        s(ist) = lastfinite;
    end
    y(ist,:) = y(ist,:) / mu(ist,s(ist));
end
for ist=1:M
    gamma(ist,:) = mu(ist,:)/mu(ist,s(ist));
    if max(abs(mu(ist,:)-(1:sum(N)))) < options.tol
        %isDelay(i) = true;
    end
end

% eliminating the delays seems to produce problems
% Z = sum([Z; L(isDelay,:)],1);
% L(isDelay,:)=[];
% mu(isDelay,:)=[];
% gamma(isDelay,:)=[];
% y(isDelay,:)=[];
% isLI(isDelay) = [];
% M = M - sum(isDelay);

beta = ones(M,ceil(sum(N)));
for ist=1:M
    beta(ist,1) = gamma(ist,1) / (1-gamma(ist,1)) ;
    for j=2:sum(N)
        beta(ist,j) = (1-gamma(ist,j-1)) * (gamma(ist,j) / (1-gamma(ist,j)));
    end
end
beta(isnan(beta))=Inf;

if (all(beta==Inf))
    options.method='default';
    lGN = pfqn_nc(lambda,L,N,Z,options);
    return
else
    Cgamma=0;
    sld = s(s>1);
    vmax = min(sum(sld-1),ceil(sum(N)));
    Y = pfqn_mva(y,N,0*N); % single class model
    rhoN = y*Y';
    for vtot=1:vmax
        lEN(vtot+1) = real(pfqn_gldsingle(rhoN,vtot,beta));
    end
    
    for vtot=0:vmax
        EN = exp(lEN(vtot+1));
        Cgamma = Cgamma + ((sum(N)-max(0,max(vtot-1)))/sum(N)) * EN;
    end
    options.method='default';
    lGN = pfqn_nc(lambda,y,N,Z,options);
    lGN = lGN + log(Cgamma);
end
end
