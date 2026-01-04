%{ @file qbd_setupdelayoff.m
 %  @brief Queue-length analysis for system with setup delay and turn-off
 %
 %  @author LINE Development Team
%}

%{
 % @brief Analyzes queue length for system with setup and turn-off phases
 %
 % @details
 % This function performs queue-length analysis for a queueing system with
 % setup delay and turn-off periods using QBD methods.
 %
 % @par Syntax:
 % @code
 % QN = qbd_setupdelayoff(lambda, mu, alpharate, alphascv, betarate, betascv)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lambda<td>Arrival rate
 % <tr><td>mu<td>Service rate
 % <tr><td>alpharate<td>Rate of setup delay phase
 % <tr><td>alphascv<td>Squared coefficient of variation for setup delay
 % <tr><td>betarate<td>Rate of turn-off phase
 % <tr><td>betascv<td>Squared coefficient of variation for turn-off period
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>QN<td>Average queue length
 % </table>
%}
function QN = qbd_setupdelayoff(lambda, mu, alpharate, alphascv, betarate, betascv)

alpha = APH.fitMeanAndSCV(1/alpharate, alphascv).getProcess;
na = length(alpha{1});

beta = APH.fitMeanAndSCV(1/betarate, betascv).getProcess;
nb = length(beta{1});
n = na+nb;

F = zeros(n); % forward transitions
B = zeros(n);  % backward transitions
for i=1:na
    F(i,i) = lambda;
end
for i=1:nb
    F(na+i,na+1) = lambda;
end
F(na+1,na+1) = lambda;
B(na+1,na+1) = mu;

L = zeros(n);  % local transitions
for i=1:na
    L(i,i) = alpha{1}(i,i) -lambda;
    if i<na
        L(i,(i+1):na) = alpha{1}(i,(i+1):na);
    else
        L(na,na+1) =  -alpha{1}(end,end);
    end
end
L(na+1,na+1) = -mu -lambda;
for i=2:nb
    L(na+i,na+i) = -lambda;
end
%[B,L,F]

L0 = zeros(2);  % initial block
L0(2,1) =  beta{2}(1,1);
L0(1,1) = -lambda;
L0(2,2) = -beta{2}(1,1) -lambda;

L0 = zeros(na+nb);  % local transitions
for i=1:na
    L0(i,i) = -lambda;
end
for i=1:nb
    L0(na+i,na+i) = beta{1}(i,i) - lambda;
    if i==nb
        L0(na+i,1) = -beta{1}(i,i);
    else
        L0(na+i,na+i+1) = -beta{1}(i,i);
    end
end
%[L0,F]

[~,R,~] = QBD_CR(B,L,F);
pn = QBD_pi(B,L0,R);

QN = 0; % queue-lengths
j = n+1;
ni = 0;
while 1
    ni = ni + 1;
    QN = QN + ni*sum(pn(j:(j+n)));
    j = j +n;
    if j+n>length(pn)
        break;
    end
end
end
