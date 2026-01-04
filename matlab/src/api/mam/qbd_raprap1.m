%{ @file qbd_raprap1.m
 %  @brief Solves a RAP/RAP/1 queue using QBD methods
 %
 %  @author LINE Development Team
%}

%{
 % @brief Analyzes RAP/RAP/1 queue using Quasi-Birth-Death process
 %
 % @details
 % This function solves a RAP/RAP/1 queue (Rational Arrival Process) using
 % QBD methods, computing throughput, queue length, utilization, and other
 % performance metrics.
 %
 % @par Syntax:
 % @code
 % [XN, QN, UN, pqueue, R, eta, G, B, L, F] = qbd_raprap1(RAPa, RAPs)
 % [XN, QN, UN, pqueue, R, eta, G, B, L, F] = qbd_raprap1(RAPa, RAPs, util)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>RAPa<td>Arrival process (RAP)
 % <tr><td>RAPs<td>Service process (RAP)
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
 % <tr><td>B<td>Backward transition block
 % <tr><td>L<td>Local transition block
 % <tr><td>F<td>Forward transition block
 % </table>
%}
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
