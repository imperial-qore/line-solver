%{
%{
 % @file pfqn_mvaldms.m
 % @brief Load-dependent MVA for multiserver mixed networks (wrapper for pfqn_mvaldmx).
%}
%}

%{
%{
 % @brief Load-dependent MVA for multiserver mixed networks (wrapper for pfqn_mvaldmx).
 % @fn pfqn_mvaldms(lambda, D, N, Z, S)
 % @param lambda Arrival rate vector.
 % @param D Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param S Number of servers per station.
 % @return XN System throughput.
 % @return QN Mean queue lengths.
 % @return UN Utilization (adjusted for multiservers).
 % @return CN Cycle times.
 % @return lGN Logarithm of normalizing constant.
%}
%}
function [XN,QN,UN,CN,lGN] = pfqn_mvaldms(lambda,D,N,Z,S)
% [XN,QN,UN,CN] = PFQN_MVALDMS(LAMBDA,D,N,Z,S)
% Wrapper for pfqn_mvaldmx that adjusts utilizations to account for
% multiservers
[M,R] = size(D);
Nct = sum(N(isfinite(N)));
mu = ones(M,Nct);
for ist=1:M
    mu(ist,:) = min(1:Nct,S(ist));
end
if isempty(Z)
    Z = zeros(1,R);
end
[XN,QN,~,CN,lGN] = pfqn_mvaldmx(lambda,D,N,Z,mu,S);

openClasses = find(isinf(N));
closedClasses = setdiff(1:length(N), openClasses);
UN = zeros(M,R);
for r=closedClasses
    for ist=1:M
        UN(ist,r) = XN(r) * D(ist,r)/S(ist);
    end
end

for r=openClasses
    for ist=1:M
        UN(ist,r) = lambda(r) * D(ist,r)/S(ist);
    end
end

end
