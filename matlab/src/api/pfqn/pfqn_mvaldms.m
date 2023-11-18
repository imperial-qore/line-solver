function [XN,QN,UN,CN,lGN] = pfqn_mvaldms(lambda,D,N,Z,S)
% [XN,QN,UN,CN] = PFQN_MVALDMS(LAMBDA,D,N,Z,S)
% Wrapper for pfqn_mvaldmx that adjusts utilizations to account for
% multiservers
[M,R] = size(D);
Nct = sum(N(isfinite(N)));
mu = ones(M,Nct);
for i=1:M
    mu(i,:) = min(1:Nct,S(i));
end
if isempty(Z)
    Z = zeros(1,R);
end
[XN,QN,~,CN,lGN] = pfqn_mvaldmx(lambda,D,N,Z,mu,S);

openClasses = find(isinf(N));
closedClasses = setdiff(1:length(N), openClasses);
UN = zeros(M,R);
for r=closedClasses
    for i=1:M
        UN(i,r) = XN(r) * D(i,r)/S(i);
    end
end

for r=openClasses
    for i=1:M
        UN(i,r) = lambda(r) * D(i,r)/S(i);
    end
end

end
