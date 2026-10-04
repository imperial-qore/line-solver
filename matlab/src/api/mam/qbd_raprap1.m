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
 % The two RAPs are INDEPENDENT of each other, so the QBD phase space is the
 % product of the two phase spaces and the blocks factor as Kronecker
 % products. This function is the thin product-space wrapper; the analysis
 % itself is the block-level core qbd_rap.m, which solves an arbitrary QBD
 % with RAP components. A model whose arrival process and sequence of
 % service times share a phase space, and are therefore cross-correlated,
 % has no such product structure and must call qbd_rap directly.
 %
 % @par References:
 % N. G. Bean and B. F. Nielsen, "Quasi-Birth-and-Death Processes with
 % Rational Arrival Process Components", Stochastic Models, 26(3), 2010,
 % pp. 309-334. The analysis rests on the prediction-process interpretation
 % of a RAP due to Asmussen and Bladt, which is what allows a QBD argument
 % to be carried over to matrices that are not nonnegative. The same
 % prediction process underlies the conditional-vector RAP sampler in
 % rap_sample.m.
 %
 % @par Phase ordering:
 % The QBD phase is the pair (arrival phase, service phase) laid out with
 % the ARRIVAL phase major and the service phase minor, i.e. the phase index
 % is (a-1)*ns + s. That is the ordering produced by kron(RAPa, eye(ns)) and
 % kron(eye(na), RAPs), and pqueue is indexed by it downstream in
 % solver_mam_basic, so the two Kronecker factors must not be swapped.
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

na = length(RAPa{1});
ns = length(RAPs{1});

if nargin>=3 %exist('util','var')
    RAPs = map_scale(RAPs,util/map_lambda(RAPa));
end

% see _kb/03-api-layer.md (QBD boundary conventions) for rationale
F = kron(RAPa{2},eye(ns));           % arrivals, level up
L = kron(RAPa{1},eye(ns)) + kron(eye(na),RAPs{1});
B = kron(eye(na),RAPs{2});           % service completions, level down
B1 = kron(RAPa{1},eye(ns));

% Theorem 7 of Bean and Nielsen (2010), in qbd_rap.m: G from the quadratic
% matrix equation, U = L + F*G, R = F*inv(-U), and pi0 from the boundary
% equation pi0*(B1 + R*B) = 0 normalised by pi0*inv(I-R)*e = 1.
[~,~,R,G,~,eta,~,pi0] = qbd_rap(F, L, B, F, B1, 0);

% see _kb/03-api-layer.md (QBD boundary conventions) for rationale
maxNumComp = 100;
pqueue = pi0;
sumpi = sum(pi0);
numit = 1;
while sumpi < 1-1e-10 && numit < 1+maxNumComp
    pqueue(numit+1,1:(na*ns)) = pqueue(numit,:)*R;
    numit = numit+1;
    sumpi = sumpi + sum(pqueue(numit,:));
end

numLevels = size(pqueue,1);
levelProb = sum(pqueue,2);
QN = (0:(numLevels-1))*levelProb;

if na == 1 && ns == 1
    UN = 1 - pqueue(1);
else
    UN= 1 - sum(pqueue(1,:));
end
XN=map_lambda(RAPa);
end
