function [XN,QN,UN,RN,it]=pfqn_bs(L,N,Z,tol,maxiter,QN0,type)
% [XN,QN,UN,RN]=PFQN_BS(L,N,Z,TOL,MAXITER,QN)

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
if nargin<6 || isempty(QN0) %~exist('QN','var')
    QN = repmat(N,M,1)/M;
else
    QN = QN0;
end
if nargin<7
    type = SchedStrategy.ID_PS * ones(M,1);
end

XN=zeros(1,R);
UN=zeros(M,R);
for it=1:maxiter
    QN_1 = QN;
    for r=1:R
        for i=1:M
            CN(i,r) = L(i,r);
            if L(i,r) == 0
                % 0 service demand at this station => this class does not visit the current node
                continue;
            end
            for s=1:R
                if s~=r
                    if type(i) == SchedStrategy.ID_FCFS
                        CN(i,r) = CN(i,r) + L(i,s)*QN(i,s);
                    else
                        CN(i,r) = CN(i,r) + L(i,r)*QN(i,s);
                    end
                else
                    CN(i,r) = CN(i,r) + L(i,r)*QN(i,r)*(N(r)-1)/N(r);
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
