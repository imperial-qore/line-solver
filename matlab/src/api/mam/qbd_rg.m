%{ @file qbd_rg.m
 %  @brief Computes QBD matrices R and G for MAP/MAP/1 queue
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes QBD transition matrices for MAP/MAP/1 queue
 %
 % @details
 % This function computes the QBD (Quasi-Birth-Death) matrices R and G
 % along with transition blocks B, L, F for a MAP/MAP/1 queue.
 %
 % @par Syntax:
 % @code
 % [R, G, B, L, F, U] = qbd_rg(MAPa, MAPs)
 % [R, G, B, L, F, U] = qbd_rg(MAPa, MAPs, util)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>MAPa<td>Markovian Arrival Process for arrivals
 % <tr><td>MAPs<td>Markovian Arrival Process for service
 % <tr><td>util<td>(Optional) Target utilization for scaling service rate
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>R<td>Rate matrix R (forward transitions)
 % <tr><td>G<td>Rate matrix G (backward transitions)
 % <tr><td>B<td>Backward transition matrix
 % <tr><td>L<td>Local transition matrix
 % <tr><td>F<td>Forward transition matrix
 % <tr><td>U<td>Utilization
 % </table>
%}
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
