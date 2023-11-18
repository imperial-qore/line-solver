function [XN,QN,UN,pqueue,R,eta,G,B,L,F]=qbd_raprap1(RAPa,RAPs,util)
% [XN,QN,UN,PQUEUE,R,ETA]=QBD_RAPRAP1(RAPA,RAPS,UTIL)

%[XN,QN,UN,pqueue,R]=qbd_raprap1(RAPa,RAPs,util)
na = length(RAPa{1});
ns = length(RAPs{1});

if nargin>=3 %exist('util','var')
    RAPs = map_scale(RAPs,util/map_lambda(RAPa));
end
util = map_lambda(RAPa) / map_lambda(RAPs);

[QN,pqueue,R,G,B,L,F] = Q_RAP_RAP_1(RAPa{1},RAPa{2},RAPs{1},RAPs{2});
eta = max(abs(eigs(R,1)));

if na == 1 && ns == 1
    UN = 1 - pqueue(1);
else
    UN= 1 - sum(pqueue(1,:));
end
XN=map_lambda(RAPa);
end
