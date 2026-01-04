%{
%{
 % @file pfqn_bs.m
 % @brief Bard-Schweitzer Approximate Mean Value Analysis (MVA).
%}
%}

function [XN,QN,UN,RN,it]=pfqn_bs(L,N,Z,tol,maxiter,QN0,type)
%{
%{
 % @brief Bard-Schweitzer Approximate Mean Value Analysis (MVA).
 % @fn pfqn_bs(L, N, Z, tol, maxiter, QN0, type)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param tol Tolerance for convergence.
 % @param maxiter Maximum number of iterations.
 % @param QN0 Initial guess for queue lengths.
 % @param type Scheduling strategy type (default: PS).
 % @return XN System throughput.
 % @return QN Mean queue lengths.
 % @return UN Utilization.
 % @return RN Residence times.
 % @return it Number of iterations performed.
%}
%}
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
    type = SchedStrategy.PS * ones(M,1);
end

XN=zeros(1,R);
UN=zeros(M,R);
for it=1:maxiter
    QN_1 = QN;
    for r=1:R
        for ist=1:M
            CN(ist,r) = L(ist,r);
            if L(ist,r) == 0
                % 0 service demand at this station => this class does not visit the current node
                continue;
            end
            for s=1:R
                if s~=r
                    if type(ist) == SchedStrategy.FCFS
                        CN(ist,r) = CN(ist,r) + L(ist,s)*QN(ist,s);
                    else
                        CN(ist,r) = CN(ist,r) + L(ist,r)*QN(ist,s);
                    end
                else
                    CN(ist,r) = CN(ist,r) + L(ist,r)*QN(ist,r)*(N(r)-1)/N(r);
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
