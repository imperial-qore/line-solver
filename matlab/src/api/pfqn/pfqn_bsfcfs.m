%{
%{
 % @file pfqn_bsfcfs.m
 % @brief Bard-Schweitzer approximate MVA for FCFS scheduling with weighted priorities.
%}
%}

%{
%{
 % @brief Bard-Schweitzer approximate MVA for FCFS scheduling with weighted priorities.
 % @fn pfqn_bsfcfs(L, N, Z, tol, maxiter, QN, weight)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector (default: zeros).
 % @param tol Convergence tolerance (default: 1e-6).
 % @param maxiter Maximum number of iterations (default: 1000).
 % @param QN Initial queue length matrix (default: uniform distribution).
 % @param weight Weight matrix for relative priorities (default: ones).
 % @return XN System throughput.
 % @return QN Mean queue lengths.
 % @return UN Utilization.
 % @return RN Residence times.
 % @return it Number of iterations performed.
%}
%}
function [XN,QN,UN,RN,it]=pfqn_bsfcfs(L,N,Z,tol,maxiter,QN,weight)
% [XN,QN,UN,RN]=PFQN_BSFCFS(L,N,Z,TOL,MAXITER,QN,WEIGHT)

if nargin<3%~exist('Z','var')
    Z=0*N;
end
if nargin<4%~exist('tol','var')
    tol = 1e-6;
end
if nargin<5%~exist('maxiter','var')
    maxiter = 1000;
end
[M,R]=size(L);
CN=zeros(M,R);
if nargin<6 %~exist('QN','var')
    QN=repmat(N,M,1)/M;
else
    QN = QN+eps; % 0 gives problems
end
if nargin<7 %~exist('weight','var')
    weight = ones(M,R);
end
XN=zeros(1,R);
UN=zeros(M,R);
relprio=zeros(M,R);
for it=1:maxiter
    QN_1 = QN;
    for ist=1:M
        for r=1:R
            relprio(ist,r) = (QN(ist,r)*weight(ist,r));
        end
    end
    for r=1:R
        for ist=1:M
            CN(ist,r) = L(ist,r);
            for s=1:R
                if s~=r
                    % FCFS approximation
                    CN(ist,r) = CN(ist,r) + L(ist,s)*QN(ist,s)*relprio(ist,s)/relprio(ist,r);
                else
                    CN(ist,r) = CN(ist,r) + L(ist,r)*QN(ist,r)*(N(r)-1)/N(r)*relprio(ist,s)/relprio(ist,r);
                end
            end
        end
        XN(r) = N(r)/(Z(r)+sum(CN(:,r)));
    end
    for r=1:R
        for ist=1:M
            QN(ist,r) = XN(r)*CN(ist,r);
        end
    end
    for r=1:R
        for ist=1:M
            UN(ist,r) = XN(r)*L(ist,r);
        end
    end
    if max(abs(1-QN./QN_1)) < tol
        break
    end
end
RN = QN ./ repmat(XN,M,1);
end
