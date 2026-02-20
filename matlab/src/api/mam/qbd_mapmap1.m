%{ @file qbd_mapmap1.m
 %  @brief Solves a MAP/MAP/1 queue using QBD methods
 %
 %  @author LINE Development Team
%}

%{
 % @brief Analyzes MAP/MAP/1 queue using Quasi-Birth-Death process
 %
 % @details
 % This function solves a MAP/MAP/1 queue using QBD (Quasi-Birth-Death)
 % methods, computing throughput, queue length, utilization, and other
 % performance metrics.
 %
 % @par Syntax:
 % @code
 % [XN, QN, UN, pqueue, R, eta, G, A_1, A0, A1, U, MAPs] = qbd_mapmap1(MAPa, MAPs)
 % [XN, QN, UN, pqueue, R, eta, G, A_1, A0, A1, U, MAPs] = qbd_mapmap1(MAPa, MAPs, util)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>MAPa<td>Arrival process (MAP)
 % <tr><td>MAPs<td>Service process (MAP)
 % <tr><td>util<td>(Optional) Target utilization to scale service rate
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>XN<td>System throughput
 % <tr><td>QN<td>Mean queue length
 % <tr><td>UN<td>Utilization
 % <tr><td>pqueue<td>Queue length distribution
 % <tr><td>R<td>Rate matrix R
 % <tr><td>eta<td>Caudal characteristic
 % <tr><td>G<td>Rate matrix G
 % <tr><td>A_1<td>Downward transition block
 % <tr><td>A0<td>Local transition block
 % <tr><td>A1<td>Upward transition block
 % <tr><td>U<td>Matrix U
 % <tr><td>MAPs<td>Scaled service process
 % </table>
%}
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

pqueue = QBD_pi(A_1,A0bar,R,'MaxNumComp',1e2);
pqueue = reshape(pqueue',length(R),length(pqueue)/length(R))';
if sum(sum(pqueue(2:end,:)))<util*0.99
    pqueue = QBD_pi(A_1,A0bar,R,'MaxNumComp',2e4);
end

if na == 1 && ns == 1
    UN = 1 - pqueue(1);
    QN=(0:(size(pqueue,1)-1))*pqueue;
    %QN=pqueue(2,:)*inv(eye(size(R))-R)^2*ones(length(R),1);
else
    UN= 1 - sum(pqueue(1,:));
    QN=(0:(size(pqueue,1)-1))*sum(pqueue,2);
    %QN=pqueue(2,:)*inv(eye(size(R))-R)^2*ones(length(R),1);
end
XN=map_lambda(MAPa);
end
