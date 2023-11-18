function [XN,QN,UN,pqueue,R,eta,G,A_1,A0,A1,U,MAPs]=qbd_mapmap1(MAPa,MAPs,util)
% [XN,QN,UN,PQUEUE,R,ETA]=QBD_MAPMAP1(MAPA,MAPS,UTIL)

%[XN,QN,UN,pqueue,R]=qbd_mapmap1(MAPa,MAPs,util)
na = length(MAPa{1});
ns = length(MAPs{1});

if nargin>=3%exist('util','var')
    MAPs = map_scale(MAPs,util/map_lambda(MAPa));
end
util = map_lambda(MAPa) / map_lambda(MAPs);
A1 = kron(MAPa{2},eye(ns));
A0 = krons(MAPa{1},MAPs{1});
A_1 = kron(eye(na),MAPs{2});
A0bar = kron(MAPa{1},eye(ns));

[G,R,U] = QBD_CR(A_1,A0,A1);

alpha = max(-diag(A0));
A0dt = A0/alpha+eye(size(A0));
A_1dt = A_1/alpha;
A1dt = A1/alpha;
if nargout>5
    [eta] = QBD_Caudal(A_1dt,A0dt,A1dt);
end
%%warning off;
pqueue = QBD_pi(A_1,A0bar,R); %,'MaxNumComp',1e5);
if sum(sum(pqueue(2:end,:)))<util*0.99
    %pqueue = QBD_pi(A_1,A0bar,R,'MaxNumComp',2e4);
end

if na == 1 && ns == 1
    UN = 1 - pqueue(1);
    %QN=(0:(size(pqueue,1)-1))*pqueue;
    QN=pqueue(2,:)*inv(eye(size(R))-R)^2*ones(length(R),1);
else
    UN= 1 - sum(pqueue(1,:));
    %QN=(0:(size(pqueue,1)-1))*sum(pqueue,2);
    QN=pqueue(2,:)*inv(eye(size(R))-R)^2*ones(length(R),1);
end
XN=map_lambda(MAPa);
end
