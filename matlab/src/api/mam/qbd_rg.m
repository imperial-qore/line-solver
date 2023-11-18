function [R,G,B,L,F,U]=qbd_rg(MAPa,MAPs,util)
% [XN,QN,UN,PQUEUE,R,ETA]=QBD_MAPMAP1(MAPA,MAPS,UTIL)

%[XN,QN,UN,pqueue,R]=qbd_mapmap1(MAPa,MAPs,util)
na = length(MAPa{1});
ns = length(MAPs{1});

if nargin>=3%exist('util','var')
    MAPs = map_scale(MAPs,util/map_lambda(MAPa));
end
util = map_lambda(MAPa) / map_lambda(MAPs);
F = kron(MAPa{2},eye(ns));
L = krons(MAPa{1},MAPs{1});
B = kron(eye(na),MAPs{2});
A0bar = kron(MAPa{1},eye(ns));

[G,R,U] = QBD_CR(B,L,F);
end
