%{
%{
 % @file pfqn_lcp.m
 % @brief Bard Large Customer Population (LCP) approximate MVA.
%}
%}

function [XN,QN,UN,RN,it]=pfqn_lcp(L,N,Z,tol,maxiter,QN0,type)
%{
%{
 % @brief Bard Large Customer Population (LCP) approximate MVA.
 %
 % Y. Bard, "Some extensions to multiclass queueing network analysis", in
 % Performance of Computer Systems, North-Holland, 1979. The first
 % approximate MVA algorithm: it estimates the arrival-instant queue length
 % by the time-averaged one WITHOUT removing the arriving customer,
 %
 %   A_k^(c)(N) = Q_k(N - 1_c) ~= Q_k(N) = sum_s Q_ks(N),
 %
 % since with a large population one customer less cannot change the mean
 % queue lengths appreciably. Setting the Bard-Schweitzer proportional term
 % Q_kc(N)/N_c to zero recovers this algorithm, so LCP is uniformly more
 % pessimistic than pfqn_bs and is inaccurate at small populations.
 %
 % @fn pfqn_lcp(L, N, Z, tol, maxiter, QN0, type)
 % @param L Service demand matrix (stations x classes).
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

if nargin<3
    Z=0*N;
end
if nargin<4
    tol = 1e-6;
end
if nargin<5
    maxiter = 1000;
end

[M,R]=size(L);
CN=zeros(M,R);
if nargin<6 || isempty(QN0)
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
        if N(r) == 0
            % Empty class: contributes no jobs anywhere, see pfqn_bs
            XN(r) = 0;
            CN(:,r) = 0;
            QN(:,r) = 0;
            UN(:,r) = 0;
            continue;
        end
        for ist=1:M
            CN(ist,r) = L(ist,r);
            if L(ist,r) == 0
                continue;
            end
            for s=1:R
                % LCP eq (2.8): the arriving customer is NOT removed, so the
                % class-r term carries no (N(r)-1)/N(r) factor
                if type(ist) == SchedStrategy.FCFS && s~=r
                    CN(ist,r) = CN(ist,r) + L(ist,s)*QN(ist,s);
                else
                    CN(ist,r) = CN(ist,r) + L(ist,r)*QN(ist,s);
                end
            end
        end
        XN(r) = N(r)/(Z(r)+sum(CN(:,r)));
    end
    for r=1:R
        for ist=1:M
            QN(ist,r) = XN(r)*CN(ist,r);
            UN(ist,r) = XN(r)*L(ist,r);
        end
    end
    nz = N > 0;
    if isempty(find(nz,1)) || max(max(abs(1-QN(:,nz)./QN_1(:,nz)))) < tol
        break
    end
end
RN = QN ./ repmat(XN,M,1);
RN(:,N==0) = 0;
end
