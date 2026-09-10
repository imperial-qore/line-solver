%{
%{
 % @file pfqn_chow.m
 % @brief Chow Second Approximation (SA) approximate MVA.
%}
%}

function [XN,QN,UN,RN,it]=pfqn_chow(L,N,Z,tol,maxiter,QN0,type,variant)
%{
%{
 % @brief Chow Second Approximation (SA) approximate MVA.
 %
 % W.-M. Chow, "Approximations for large scale closed queueing networks",
 % Perform. Eval. 3(1), 1983. The arrival-instant queue length is written
 % exactly as
 %
 %   A_k^(c)(N) = Q_k(N - 1_c) = Q_k(N) (1 + theta_ck),
 %   theta_ck   = [Q_k(N - 1_c) - Q_k(N)] / Q_k(N),
 %
 % and the theta-terms are estimated ONCE, off the Bard LCP solution, before
 % the fixed point over (2.2)-(2.6) is run. Two estimators are given:
 %
 %   'backward'  theta_ck = [Qhat_k(N - 1_c) - Qhat_k(N)] / Qhat_k(N)      (2.16)
 %   'forward'   theta_ck = [Qhat_k(N) - Qhat_k(N + 1_c)] / Qhat_k(N + 1_c) (2.17)
 %
 % Chow reports the forward form to be the more accurate of the two, so it is
 % the default here. Setting every theta to zero recovers pfqn_lcp.
 %
 % @fn pfqn_chow(L, N, Z, tol, maxiter, QN0, type, variant)
 % @param L Service demand matrix (stations x classes).
 % @param N Population vector.
 % @param Z Think time vector.
 % @param tol Tolerance for convergence.
 % @param maxiter Maximum number of iterations.
 % @param QN0 Initial guess for queue lengths.
 % @param type Scheduling strategy type (default: PS).
 % @param variant 'forward' (eq. 2.17, default) or 'backward' (eq. 2.16).
 % @return XN System throughput.
 % @return QN Mean queue lengths.
 % @return UN Utilization.
 % @return RN Residence times.
 % @return it Number of iterations performed.
%}
%}

if nargin<3 || isempty(Z)
    Z=0*N;
end
if nargin<4 || isempty(tol)
    tol = 1e-6;
end
if nargin<5 || isempty(maxiter)
    maxiter = 1000;
end
[M,R]=size(L);
if nargin<6
    QN0 = [];
end
if nargin<7 || isempty(type)
    type = SchedStrategy.PS * ones(M,1);
end
if nargin<8 || isempty(variant)
    variant = 'forward';
end

%% theta-terms from the LCP solution
[~,Qlcp] = pfqn_lcp(L,N,Z,tol,maxiter,QN0,type);
Qtot = sum(Qlcp,2);   % Qhat_k(N)
theta = zeros(M,R);
for r=1:R
    if N(r) == 0
        continue
    end
    switch variant
        case 'backward'
            Nr = oner(N,r);
            [~,Qr] = pfqn_lcp(L,Nr,Z,tol,maxiter,QN0,type);
            base = Qtot;
            delta = sum(Qr,2) - Qtot;
        otherwise % 'forward', eq. (2.17)
            Np = N; Np(r) = Np(r) + 1;
            [~,Qp] = pfqn_lcp(L,Np,Z,tol,maxiter,QN0,type);
            base = sum(Qp,2);
            delta = Qtot - base;
    end
    for ist=1:M
        if base(ist) > 0
            theta(ist,r) = delta(ist)/base(ist);
        end
    end
end

%% fixed point with A_k^(c) = Q_k (1 + theta_ck)
CN=zeros(M,R);
if isempty(QN0)
    QN = repmat(N,M,1)/M;
else
    QN = QN0;
end
XN=zeros(1,R);
UN=zeros(M,R);
for it=1:maxiter
    QN_1 = QN;
    for r=1:R
        if N(r) == 0
            XN(r) = 0; CN(:,r) = 0; QN(:,r) = 0; UN(:,r) = 0;
            continue;
        end
        for ist=1:M
            CN(ist,r) = L(ist,r);
            if L(ist,r) == 0
                continue;
            end
            for s=1:R
                if type(ist) == SchedStrategy.FCFS && s~=r
                    CN(ist,r) = CN(ist,r) + L(ist,s)*QN(ist,s)*(1+theta(ist,r));
                else
                    CN(ist,r) = CN(ist,r) + L(ist,r)*QN(ist,s)*(1+theta(ist,r));
                end
            end
            % a theta below -1 would make the arrival-instant queue negative
            CN(ist,r) = max(CN(ist,r), L(ist,r));
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
