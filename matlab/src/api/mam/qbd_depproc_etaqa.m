%{ @file qbd_depproc_etaqa.m
 %  @brief Constructs MAP departure process for MAP/MAP/1-FCFS via ETAQA truncation
 %
 %  @author LINE Development Team
%}

%{
 % @brief Builds a MAP approximation of the departure process for a MAP/MAP/1 FCFS queue
 %
 % @details
 % This function constructs a finite-state MAP {D0, D1} representation of the
 % departure process from a MAP/MAP/1 queue with FCFS discipline, using ETAQA
 % truncation at QBD level n. The resulting MAP captures both the SCV (squared
 % coefficient of variation) and autocorrelation structure of inter-departure
 % times.
 %
 % The QBD process is formed with forward matrix F = kron(MAPa{2}, I_ns),
 % local matrix L = krons(MAPa{1}, MAPs{1}), and backward matrix
 % B = kron(I_na, MAPs{2}). Levels 0..n-1 are represented explicitly, and
 % level n uses a tail approximation via the G matrix.
 %
 % @par Syntax:
 % @code
 % D = qbd_depproc_etaqa(MAPa, MAPs, n)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>MAPa<td>Arrival process in MAP format {D0, D1}
 % <tr><td>MAPs<td>Service process in MAP format {D0, D1}
 % <tr><td>n<td>Truncation level (number of QBD levels to represent explicitly)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>D<td>Departure process in MAP format {D0, D1}
 % </table>
%}
function [D]=qbd_depproc_etaqa(MAPa,MAPs,n)

na = length(MAPa{1});
ns = length(MAPs{1});
lvlsz = ns*na;

F = kron(MAPa{2},eye(ns));
L = krons(MAPa{1},MAPs{1});
B = kron(eye(na),MAPs{2});
L0 = kron(MAPa{1},eye(ns));

[~,R,~] = QBD_CR(B,L,F);
G = inv(-L-R*B)*B;
Lhat = F+L;
Bbar = B+F*G;
Bhat = F*G;

vn1 = ones(1,n-1);
vn = zeros(1,n); vn(end)=1;

D0=kron(diag([1,vn1]),L)+kron(diag(vn1,1),F);
D0=[zeros(size(L0,1),size(D0,2)); D0];
D0 = [zeros(size(D0,1),size(B,2)),D0];
D0(1:size(L0,1),1:(size(L0,2)+size(F,2)))=[L0,F];
D0(((n-1)*lvlsz+1):n*lvlsz,((n-1)*lvlsz+1):n*lvlsz)=Lhat;

D1=kron(diag(vn,-1),Bbar)+kron(diag([0,vn]),Bhat);
D1(1:n*lvlsz,1:n*lvlsz)=kron(diag(vn1,-1),B);
D={D0,D1};
end
