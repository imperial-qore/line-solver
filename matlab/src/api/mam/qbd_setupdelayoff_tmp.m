function [QN] = qbd_setupdelayoff(lambda, mu, alpharate, alphascv, betarate, betascv)

alpha = APH.fitMeanAndSCV(1/alpharate, alphascv).getRepres;
na = length(alpha{1});

betascv = 1.0;
beta = APH.fitMeanAndSCV(1/betarate, betascv).getRepres;
nb = length(beta{1})
n = na+nb;

F = diag(lambda*ones(1,n)); % forward transitions
B = zeros(n);  % backward transitions
for i=1:nb
    B(na+i,na+i) = mu;
end

L = zeros(2);  % local transitions
for i=1:na
    L(i,i) = alpha{1}(i,i) -lambda;
    if i<na
        L(i,(i+1):na) = alpha{1}(i,(i+1):na);
    else
        L(na,na+1) =  alpha{2}(end,end);
    end
end
for i=1:nb
    L(na+i,na+i) = -mu -lambda;
end

L0 = zeros(2);  % initial block
L0(2,1) =  beta{2}(1,1);
L0(1,1) = -lambda;
L0(2,2) = -beta{2}(1,1) -lambda;

L0 = zeros(2);  % local transitions
for i=1:na
    L0(i,i) = alpha{1}(i,i) -lambda;
    if i<na
        L0(i,(i+1):na) = alpha{1}(i,(i+1):na);
    else
        L0(na,na+1) =  alpha{2}(end,end);
    end
end
for i=1:nb
    L0(na+i,na+i-1) = beta{2}(i,i);
    L0(na+i,na+i) = -beta{2}(i,i) -lambda;
end

[~,R,~] = QBD_CR(B,L,F);
pn = QBD_pi(B,L0,R);
pn = pn(1:n:end);
n


QN = 0; % queue-lengths
for n=1:(length(pn)-1)
    QN = QN + n*pn(n+1);
end
end