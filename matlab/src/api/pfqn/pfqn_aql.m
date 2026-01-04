%{
%{
 % @file pfqn_aql.m
 % @brief Approximate Queue Length (AQL) algorithm for product-form networks.
%}
%}

function [XN,QN,UN,RN,numIters,AN]=pfqn_aql(L,N,Z,TOL,MAXITER,QN0)
%{
%{
 % @brief Approximate Queue Length (AQL) algorithm for product-form networks.
 % @fn pfqn_aql(L, N, Z, TOL, MAXITER, QN0)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param TOL Tolerance for convergence.
 % @param MAXITER Maximum number of iterations.
 % @param QN0 Initial guess for queue lengths.
 % @return XN System throughput.
 % @return QN Mean queue lengths.
 % @return UN Utilization.
 % @return RN Residence times.
 % @return numIters Number of iterations.
 % @return AN Average arrival rate at nodes.
%}
%}
[M,K]=size(L);
if nargin<3
    Z=zeros(1,K);
end
if nargin<6 || isempty(QN0) %~exist('QN','var')
    QN0 = repmat(N,M,1)/M;
else
    QN0 = QN0+eps; % 0 gives problems
end
if nargin<4
    TOL = 1e-7;
end
if nargin<5
    MAXITER=1e3;
end
Q=cell(1,K+1);
R=cell(1,K+1);
X=cell(1,K+1);
gamma = zeros(M,K);

for t=0:K
    n=oner(N,t);
    for k=1:M
        Q{t+1}(k,1)=QN0(k);
    end
end
it = 0;
while 1
    Q_olditer = Q;
    it = it + 1;
    for t=0:K
        n=oner(N,t);
        for k=1:M
            for s=1:K
                R{t+1}(k,s) = L(k,s)*(1+(sum(n)-1)*(Q{t+1}(k)/sum(n)-gamma(k,s)));
            end
        end
        for s=1:K
            X{t+1}(s) = n(s)/(Z(s)+sum(R{t+1}(:,s)));
        end
        for k=1:M
            Q{t+1}(k) = X{t+1}(:)'*R{t+1}(k,:)';
        end
    end % for t
    for k=1:M
        for s=1:K
            gamma(k,s) = (Q{0+1}(k)/sum(N)) - (Q{s+1}(k)/(sum(N)-1));
        end
    end
    if max(abs((Q_olditer{1}(:)-Q{1}(:))./Q{1}(:))) < TOL || it == MAXITER
        numIters=it;
        break
    end
end
XN = X{1};
RN = R{1};
for k=1:M
    for s=1:K
        UN(k,s) = XN(s)*L(k,s);
        QN(k,s) = UN(k,s)*(1+Q{s+1}(k));
        AN(k,s) = Q{s+1}(k);
    end
end
end