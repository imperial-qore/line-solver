%{
%{
 % @file pfqn_mvamx.m
 % @brief Exact MVA for mixed open/closed single-server networks.
%}
%}

%{
%{
 % @brief Exact MVA for mixed open/closed single-server networks.
 % @fn pfqn_mvamx(lambda, D, N, Z, mi)
 % @param lambda Arrival rate vector.
 % @param D Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param mi Queue replication factors (default: ones).
 % @return XN System throughput.
 % @return QN Mean queue lengths.
 % @return UN Utilization.
 % @return CN Cycle times.
 % @return lGN Logarithm of normalizing constant.
%}
%}
function [XN,QN,UN,CN,lGN] = pfqn_mvamx(lambda,D,N,Z, mi)
% [XN,QN,UN,CN,LGN] = PFQN_MVAMX(LAMBDA,D,N,Z, MI)

if any(N(lambda>0)>0 & isfinite(N(lambda>0)))
    line_error(mfilename,'Arrival rate cannot be specified on closed classes.');
end
[M,R] = size(D);
if nargin<5 %~exist('mi','var')
    mi = ones(M,1);
end
openClasses = find(isinf(N));
closedClasses = setdiff(1:length(N), openClasses);

XN = zeros(1,R);
UN = zeros(M,R);
CN = zeros(M,R);
QN = zeros(M,R);
for r=openClasses
    for ist=1:M
        UN(ist,r) = lambda(r)*D(ist,r);
    end
    XN(r) = lambda(r);
end

UNt = sum(UN,2);

if isempty(Z)
    Z = zeros(1,R);
end
Dc = D(:,closedClasses) ./ (1-repmat(UNt,1,length((closedClasses))));
[XNc,QNc,~,CNc,lGN] = pfqn_mva(Dc,N(closedClasses),Z(closedClasses),mi);
XN(closedClasses) = XNc;
QN(:,closedClasses) = QNc;
CN(:,closedClasses) = CNc;
for ist = 1:M
    for r=closedClasses
        UN(ist,r) = XN(r)*D(ist,r);
    end
end
for ist = 1:M
    for r=openClasses
        if isempty(QNc)
            CN(ist,r) = D(ist,r) / (1-UNt(ist));
        else
            CN(ist,r) = D(ist,r) * (1+sum(QNc(ist,:))) / (1-UNt(ist));
        end
        QN(ist,r) = CN(ist,r) * XN(r);
    end
end
end
