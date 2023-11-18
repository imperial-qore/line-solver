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
for i=1:M
    if all(mu(i,:)==mu(i,1)) % LI station
        L(i,:) = L(i,:) / mu(i,1);
        mu(i,:) = 1;
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
for i=1:M
    if isfinite(mu(i,end))
        s(i) = min(find(abs(mu(i,:)-mu(i,end))<options.tol));
        if s(i)==0
            s(i) = sum(N);
        end
    else
        s(i) = sum(N);
    end
end
isDelay = false(M,1);
isLI = false(M,1);
y = L;

for i=1:M
    if isinf(mu(i,s(i)))
        lastfinite=max(find(isfinite(mu(i,:))));
        s(i) = lastfinite;
    end
    y(i,:) = y(i,:) / mu(i,s(i));
end
for i=1:M
    gamma(i,:) = mu(i,:)/mu(i,s(i));
    if max(abs(mu(i,:)-(1:sum(N)))) < options.tol
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
for i=1:M
    beta(i,1) = gamma(i,1) / (1-gamma(i,1)) ;
    for j=2:sum(N)
        beta(i,j) = (1-gamma(i,j-1)) * (gamma(i,j) / (1-gamma(i,j)));
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
