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
    for i=1:M
        for r=1:R
            relprio(i,r) = (QN(i,r)*weight(i,r));
        end
    end
    for r=1:R
        for i=1:M
            CN(i,r) = L(i,r);
            for s=1:R
                if s~=r
                    % FCFS approximation
                    CN(i,r) = CN(i,r) + L(i,s)*QN(i,s)*relprio(i,s)/relprio(i,r);
                else
                    CN(i,r) = CN(i,r) + L(i,r)*QN(i,r)*(N(r)-1)/N(r)*relprio(i,s)/relprio(i,r);
                end
            end
        end
        XN(r) = N(r)/(Z(r)+sum(CN(:,r)));
    end
    for r=1:R
        for i=1:M
            QN(i,r) = XN(r)*CN(i,r);
        end
    end
    for r=1:R
        for i=1:M
            UN(i,r) = XN(r)*L(i,r);
        end
    end
    if max(abs(1-QN./QN_1)) < tol
        break
    end
end
RN = QN ./ repmat(XN,M,1);
end
