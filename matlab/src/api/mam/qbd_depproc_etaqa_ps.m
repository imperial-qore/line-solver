%{ @file qbd_depproc_etaqa_ps.m
 %  @brief Constructs MAP departure process for MAP/MAP/1-PS via ETAQA truncation
 %
 %  @author LINE Development Team
%}

%{
 % @brief Builds a MAP approximation of the departure process for a MAP/MAP/1 PS queue
 %
 % @details
 % This function constructs a finite-state MAP {D0, D1} representation of the
 % departure process from a MAP/MAP/1 queue with Processor Sharing (PS)
 % discipline, using ETAQA truncation at QBD level n.
 %
 % Compared to the FCFS variant (qbd_depproc_etaqa), the PS discipline splits
 % service completions at level j into a departure component B*(1/j) and an
 % internal transition component B*(1-1/j), reflecting the rate-dependent
 % sharing of the server among j jobs. The tail approximation at level n uses
 % Bbar and Bhat matrices weighted by 1/n and (n-1)/n respectively.
 %
 % @par Syntax:
 % @code
 % D = qbd_depproc_etaqa_ps(MAPa, MAPs, n)
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
function [D]=qbd_depproc_etaqa_ps(MAPa,MAPs,n)

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

D0 = D0 + kron(diag(vn,-1),Bbar)+kron(diag([0,vn]),Bhat)*(n-1)/n;
D1=kron(diag(vn,-1),Bbar)+kron(diag([0,vn]),Bhat)/n;
for j=1:(n-1)
    D0((j*lvlsz+1):(j+1)*lvlsz, ((j-1)*lvlsz+1):j*lvlsz)=B*(1-1/j);
    D1((j*lvlsz+1):(j+1)*lvlsz, ((j-1)*lvlsz+1):j*lvlsz)=B*(1/j);
end
D={D0,D1};
end
