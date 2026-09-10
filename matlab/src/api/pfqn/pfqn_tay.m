%{
%{
 % @file pfqn_tay.m
 % @brief Tay's arrival-instant approximate mean value analysis.
%}
%}

function [XN,QN,UN,RN,it,QNarr]=pfqn_tay(L,N,Z,tol,maxiter,QN0)
%{
%{
 % @brief Approximate MVA for closed multiclass product-form networks in which
 %        the arrival-instant queue lengths are estimated from the THROUGHPUT
 %        ELASTICITIES rather than from a population-shift heuristic (Tay 1987;
 %        presented as eqs. 4.8.2-1..3 of the Schweitzer-Serazzi-Broglia
 %        survey, where it is benchmarked against exact, Linearizer and
 %        Bard-Schweitzer on Tay's Example 4).
 %
 %        Let E_mkc = (D_mk/X_c) dX_c/dD_mk be the elasticity of the class-c
 %        throughput with respect to the class-k demand at station m. Tay shows
 %        that the elasticities satisfy the R linear equations
 %          E_mkj sum_t B_tj Q_jt (1+Q_jt) =
 %              -[(delta_jk + Q_jm) B_mk Q_km + sum_{c~=j} E_mkc sum_t B_tc Q_jt Q_ct]
 %        with B_ir = 1/(1 + D_ir X_r/N_r), and that the arrival-instant queue
 %        length is then simply
 %          Q_km^(r) = Q_km + E_mkr,
 %        which closes the MVA recursion R_rm = D_rm (1 + sum_k Q_km^(r)).
 %        One R x R solve per (station, class) pair per iteration, so the cost
 %        is O(M R (R^3 + M R^2)) per sweep: more than Bard-Schweitzer, less
 %        than Linearizer's R+1 auxiliary networks.
 %
 %        Delay stations enter through Z only. They are "AS" servers in the
 %        survey's notation (d_t = 0), so they contribute Z_j X_j to the
 %        denominator of the elasticity equations but nothing to its numerator.
 %
 %        Reference: Y. C. Tay, R. Suri, "Error bounds for performance
 %        prediction in queuing networks", ACM TOCS 3(4), 1985; Y. C. Tay,
 %        "An approach to analyzing the behavior of some queueing networks",
 %        Operations Research 40(S2), 1992; P. J. Schweitzer, G. Serazzi,
 %        M. Broglia, "A survey of bottleneck analysis in closed queueing
 %        networks", Sec. 4.8.2.
 % @fn pfqn_tay(L, N, Z, tol, maxiter, QN0)
 % @param L Service demand matrix (M x R).
 % @param N Population vector (1 x R).
 % @param Z Think time vector (1 x R). Default: zeros.
 % @param tol Convergence tolerance on the queue lengths. Default: 1e-6.
 % @param maxiter Maximum number of iterations. Default: 1000.
 % @param QN0 Initial guess for the queue lengths (M x R). Default: uniform.
 % @return XN Per-class throughputs (1 x R).
 % @return QN Mean queue lengths (M x R).
 % @return UN Utilizations (M x R).
 % @return RN Residence times (M x R).
 % @return it Number of iterations performed.
 % @return QNarr Arrival-instant queue lengths (M x R x R): QNarr(m,k,r) is the
 %        class-k queue length at station m as seen by an arriving class-r job.
 %        These are the auxiliary quantities the method is tabulated on and they
 %        are NOT the queue lengths of the model re-solved at N - e_r, which is
 %        the same object only for an exact solution.
%}
%}
if nargin<3 || isempty(Z)
    Z = 0*N;
end
if nargin<4 || isempty(tol)
    tol = 1e-6;
end
if nargin<5 || isempty(maxiter)
    maxiter = 1000;
end

[M,R] = size(L);
N = N(:)';
if size(Z,1) > 1
    Z = sum(Z,1);      % several delay stations aggregate into one think time
end
Z = Z(:)';
if numel(Z) < R, Z = sum(Z)*ones(1,R); end

XN = zeros(1,R);
QN = zeros(M,R);
UN = zeros(M,R);
RN = zeros(M,R);
QNarr = zeros(M,R,R);
it = 0;

% Empty classes contribute no jobs anywhere and make the elasticity system
% singular (their denominator is identically zero); solve the model without
% them and re-expand, as in pfqn_bs.
act = find(N > 0);
if isempty(act)
    return
end
if numel(act) < R
    [Xa,Qa,Ua,Ra,it,Qarra] = pfqn_tay(L(:,act),N(act),Z(act),tol,maxiter);
    XN(act) = Xa; QN(:,act) = Qa; UN(:,act) = Ua; RN(:,act) = Ra;
    QNarr(:,act,act) = Qarra;
    return
end

if nargin<6 || isempty(QN0)
    QN = repmat(N,M,1)/M;
else
    QN = QN0;
end
XN = N./(Z + sum(L,1).*(1+sum(QN,1)));

E = zeros(R,1);
for it=1:maxiter
    QN_1 = QN;

    B = 1./(1 + L.*repmat(XN./N,M,1));

    % Denominators of the elasticity system, one per class; the delay term
    % Z_j X_j is the AS-server contribution (d_t = 0 leaves B = 1).
    den = zeros(1,R);
    for j = 1:R
        den(j) = sum(B(:,j).*QN(:,j).*(1+QN(:,j))) + Z(j)*XN(j);
    end

    % Off-diagonal coupling C(j,c) = sum_t B_tc Q_jt Q_ct, class-pair symmetric
    % in the station sum but not in the division by den(j) applied below.
    C = zeros(R,R);
    for j = 1:R
        for c = 1:R
            C(j,c) = sum(B(:,c).*QN(:,j).*QN(:,c));
        end
    end

    Qarr = zeros(M,R,R);       % Qarr(m,k,r) = Q_km as seen by an arriving r
    for m = 1:M
        for k = 1:R
            A = eye(R);
            b = zeros(R,1);
            for j = 1:R
                for c = [1:j-1, j+1:R]
                    A(j,c) = C(j,c)/den(j);
                end
                b(j) = -((j==k) + QN(m,j))*B(m,k)*QN(m,k)/den(j);
            end
            E = A\b;
            Qarr(m,k,:) = QN(m,k) + E;
        end
    end

    for r = 1:R
        RN(:,r) = L(:,r).*(1 + sum(Qarr(:,:,r),2));
    end
    XN = N./(Z + sum(RN,1));
    QN = RN.*repmat(XN,M,1);

    if max(abs(QN(:)-QN_1(:))) < tol
        break
    end
end
UN = L.*repmat(XN,M,1);
QNarr = Qarr;
end
