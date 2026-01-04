%{
%{
 % @file pfqn_linearizermx.m
 % @brief Linearizer for mixed open/closed queueing networks.
%}
%}

%{
%{
 % @brief Linearizer for mixed open/closed queueing networks.
 % @fn pfqn_linearizermx(lambda, L, N, Z, nservers, type, tol, maxiter, method)
 % @param lambda Arrival rate vector (inf for closed classes).
 % @param L Service demand matrix.
 % @param N Population vector (inf for open classes).
 % @param Z Think time vector.
 % @param nservers Number of servers per station.
 % @param type Scheduling strategy per station.
 % @param tol Convergence tolerance (default: 1e-8).
 % @param maxiter Maximum iterations (default: 1000).
 % @param method Linearizer variant ('lin', 'gflin', 'egflin', default: 'egflin').
 % @return QN Mean queue lengths.
 % @return UN Utilization.
 % @return WN Waiting times.
 % @return CN Cycle times.
 % @return XN System throughput.
 % @return totiter Total iterations.
%}
%}
function [QN,UN,WN,CN,XN,totiter] = pfqn_linearizermx(lambda,L,N,Z,nservers,type,tol,maxiter,method)
% function [Q,U,W,C,X,totiter] = PFQN_LINEARIZERMX(lambda,L,N,Z,nservers,type,tol,maxiter,method)

if nargin<4
    maxiter = 1000;
end
if nargin<5
    tol = 1e-8;
end
if nargin<9
    method = 'egflin';
end

if any(N(lambda>0)>0 & isfinite(N(lambda>0)))
    line_error(mfilename,'Arrival rate cannot be specified on closed classes.');
end

[M,R] = size(L);
if nargin<6 || isempty(type)
    type = ones(M,1)*SchedStrategy.PS;    
end

lambda(isnan(lambda))= 0;
L(isnan(L))=0;
Z(isnan(Z))=0;
openClasses = find(isinf(N));
closedClasses = setdiff(1:length(N), openClasses);

XN = zeros(1,R);
UN = zeros(M,R);
WN = zeros(M,R);
QN = zeros(M,R);
CN = zeros(1,R);
for r=openClasses
    for ist=1:M
        UN(ist,r) = lambda(r)*L(ist,r);
    end
    XN(r) = lambda(r);
end

UNt = sum(UN,2);

if isempty(Z)
    Z = zeros(1,R);
else
    Z = sum(Z,1);
end
Dc = L(:,closedClasses) ./ (1-repmat(UNt,1,length((closedClasses))));

if max(nservers)==1
    switch method
        case {'default','lin'}
            [QNc,UNc,WNc,CNc,XNc,totiter] = pfqn_linearizer(Dc,N(closedClasses),Z(closedClasses),type,tol,maxiter);
        case 'gflin'
            linAlpha = 2.0;
            [QNc,UNc,WNc,CNc,XNc,totiter] = pfqn_gflinearizer(Dc,N(closedClasses),Z(closedClasses),type,tol,maxiter,linAlpha);
        case 'egflin'
            alphaM = zeros(1,R);
            for r=closedClasses
                alphaM(r) = 0.6 + 1.4 * exp(-8 * exp(-0.8 * N(r)));
            end
            [QNc,UNc,WNc,CNc,XNc,totiter] = pfqn_egflinearizer(Dc,N(closedClasses),Z(closedClasses),type,tol,maxiter,alphaM);
        otherwise
            [QNc,UNc,WNc,CNc,XNc,totiter] = pfqn_linearizer(Dc,N(closedClasses),Z(closedClasses),type,tol,maxiter);
    end
else
    [QNc,UNc,WNc,CNc,XNc,totiter] = pfqn_linearizerms(Dc,N(closedClasses),Z(closedClasses),nservers,type,tol,maxiter);
end

XN(closedClasses) = XNc;
QN(:,closedClasses) = QNc;
WN(:,closedClasses) = WNc;
UN(:,closedClasses) = UNc;
CN(closedClasses) = CNc;

for ist = 1:M
    for r=closedClasses
        UN(ist,r) = XN(r)*L(ist,r);
    end
end
for ist = 1:M
    for r=openClasses
        if isempty(QNc)
            WN(ist,r) = L(ist,r) / (1-UNt(ist));
        else
            WN(ist,r) = L(ist,r) * (1+sum(QNc(ist,:))) / (1-UNt(ist));
        end
        QN(ist,r) = WN(ist,r) * XN(r);
    end
end
CN(openClasses) = sum(WN(:,openClasses),1);
end

