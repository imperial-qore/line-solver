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
    for i=1:M
        UN(i,r) = lambda(r)*D(i,r);
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
for i = 1:M
    for r=closedClasses
        UN(i,r) = XN(r)*D(i,r);
    end
end
for i = 1:M
    for r=openClasses
        if isempty(QNc)
            CN(i,r) = D(i,r) / (1-UNt(i));
        else
            CN(i,r) = D(i,r) * (1+sum(QNc(i,:))) / (1-UNt(i));
        end
        QN(i,r) = CN(i,r) * XN(r);
    end
end
end
