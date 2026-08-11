%{
%{
 % @file pfqn_sens_mva.m
 % @brief Exact per-station queue-length variances and covariances for closed
 %        product-form queueing networks, computed by an MVA-like moment
 %        recursion that does not require the full sensitivity Jacobian.
%}
%}

function mom = pfqn_sens_mva(L,N,Z,mi)
%{
%{
 % @brief Exact second moments (variances and per-station covariances) of the
 %        queue lengths of a closed product-form (BCMP) queueing network. The
 %        moments are obtained by an MVA-type recursion evaluated on the same
 %        population lattice as pfqn_mva, so no derivative of the model is ever
 %        formed and the cost is O(M*R^2) per lattice point rather than the
 %        O(M^2*R^2) of the differentiated-MVA kernel used by pfqn_sens.
 %
 %        The recursion is obtained by differentiating the Reiser-Lavenberg MVA
 %        equation Q(j,v|N) = X(v|N) * L(j,v) * (mi(j) + Qtot(j|N-e_v)) with
 %        respect to the visit ratio theta(i,k) of class k at station i and
 %        rescaling. Writing W(k,i;v,j|N) = Cov[n(i,k),n(j,v)] at population N,
 %
 %          W(k,i;v,j|N) = Q(j,v|N) * ( Q(i,k|N-e_v) - Q(i,k|N) )
 %                       + [i==j & k==v] * Q(j,v|N)
 %                       + X(v|N) * L(j,v) * sum_t W(k,i;t,j|N-e_v)
 %
 %        with W(.|0) = 0. This routine evaluates the same-station case i==j,
 %        which is self-contained: the inner sum then only involves same-station
 %        terms, so a single scalar Ssum(j,k|N) = sum_t W(k,j;t,j|N) carried
 %        along the lattice closes the recursion. The cross-station case i~=j is
 %        not self-contained (it couples every station pair) and costs as much as
 %        the full Jacobian, so it is left to pfqn_sens.
 %
 %        Setting i==j and mi==1 reproduces Corollary 1 of the reference below,
 %        i.e. its equations (2.9a) for the variance and (2.10) for the
 %        covariance; the station multiplicity mi cancels identically because
 %        X(v|N)*L(j,v)*(mi(j)+Qtot(j|N-e_v)) = Q(j,v|N) is the MVA equation for
 %        any mi. The equivalent statement for an infinite-server station,
 %        equation (2.9b) of the reference, is recovered automatically because
 %        LINE folds the delay into the think time Z, which enters only through
 %        X(v|N) and carries no queue-length moment of its own.
 %
 %        Reference: E. de Souza e Silva and R. R. Muntz, "Simple Relationships
 %        Among Moments of Queue Lengths in Product Form Queueing Networks",
 %        IEEE Trans. Computers 37(9):1125-1129, 1988 (Theorems 1-2 and
 %        Corollary 1). The underlying identity Cov[n(i,k),n(j,v)] =
 %        theta(i,k) * dQ(j,v)/dtheta(i,k) is the k=2 case of Theorem 1 of
 %        I. F. Akyildiz and J. C. Strelen, "Moment Analysis for Load-Dependent
 %        Mixed Product Form Queueing Networks", IEEE Trans. Communications
 %        39(6):828-832, 1991.
 %
 % @fn pfqn_sens_mva(L, N, Z, mi)
 % @param L  Service demand matrix (M x R), L(i,r) = visits_ir / rate_ir.
 % @param N  Population vector (1 x R).
 % @param Z  Think time vector (1 x R). Default: zeros.
 % @param mi (Optional) Server multiplicity vector (1 x M). Default: ones.
 % @return mom A struct with the base measures and their second moments:
 %   .X (1 x R), .Q (M x R), .U (M x R), .R (M x R)  base MVA measures,
 %       identical to pfqn_mva(L,N,Z,mi).
 %   .QCov (M x R x R)  QCov(i,r,s) = Cov[n(i,r),n(i,s)], the queue-length
 %       covariance of classes r and s at station i. Symmetric in (r,s).
 %   .QVar (M x R)      QVar(i,r) = QCov(i,r,r) = Var[n(i,r)].
 %   .QTotVar (M x 1)   QTotVar(i) = Var[sum_r n(i,r)], the variance of the
 %       total queue length at station i, i.e. sum_{r,s} QCov(i,r,s). This is
 %       Theorem 3 of the reference, obtained here without a capacity
 %       derivative.
 %   .QCovAsym (scalar)  max |W(r,s) - W(s,r)| over the covariance entries
 %       before symmetrization. The two triangles come from differentiating two
 %       different classes' MVA equations, so this is an independent residual of
 %       the recursion and should sit at roundoff; a large value signals a bug.
 %
 % Notes:
 %  - Restricted to closed populations. Mixed and load-dependent models are
 %    handled by pfqn_sens_mvaldmx.
 %  - For a station of multiplicity mi(i)>1, which LINE treats as mi(i)
 %    identical replicas sharing the demand row L(i,:), the moments returned are
 %    those of the aggregate queue length over the replicas.
%}
%}
[M,R] = size(L);
N = ceil(N(:)');
if nargin < 3 || isempty(Z)
    Z = zeros(1,R);
end
Z = Z(:)';
if nargin < 4 || isempty(mi)
    mi = ones(1,M);
end
mi = mi(:)';
if length(N) ~= R
    line_error(mfilename,'demand matrix and population vector have different number of classes');
end
if any(isinf(N))
    line_error(mfilename,'pfqn_sens_mva requires a closed population; use pfqn_sens_mvaldmx for mixed models');
end

X = zeros(1,R); Q = zeros(M,R); U = zeros(M,R); C = zeros(M,R);
QCov = zeros(M,R,R);

if ~any(N > 0)
    mom = pack(X,Q,U,C,QCov);
    return;
end

% population-lattice odometer, identical to pfqn_mva and to the sens_mva kernel
% of pfqn_sens so that the base measures agree entry by entry
prods = zeros(1,R-1);
for w = 1:R-1
    prods(w) = prod(ones(1,R-(w+1)+1) + N(w+1:R));
end
firstnonempty = R;
while N(firstnonempty) == 0
    firstnonempty = firstnonempty - 1;
end
totpop = prod(N+1);
ctr = totpop;
Qtot = zeros(totpop,M);     % Qtot(m,i)   = sum_r Q(i,r) at population m
Qcls = zeros(totpop,M,R);   % Qcls(m,i,r) = Q(i,r) at population m
Xall = zeros(totpop,R);     % Xall(m,r)   = X(r) at population m
Ssum = zeros(totpop,M,R);   % Ssum(m,j,k) = sum_t Cov[n(j,k),n(j,t)] at pop m
currentpop = 2;
n = zeros(1,R);
n(firstnonempty) = 1;
rows = ones(1,R);           % rows(s) = lattice index of n - e_s

while ctr
    % ---- mean value analysis step at population n -----------------------
    s = 1;
    while s <= R
        pos = 0;
        if n(s) > 0
            n(s) = n(s) - 1;
            pos = n(R);
            w = 1;
            while w <= R-1
                pos = pos + n(w)*prods(w);
                w = w + 1;
            end
            n(s) = n(s) + 1;
        end
        % when n(s)==0 the index collapses to the empty population, whose
        % stored moments are zero; X(s) is then zero and every term that reads
        % rows(s) is annihilated, so no guard is needed
        row = 1 + pos;
        rows(s) = row;
        CNtot = 0;
        for i = 1:M
            C(i,s) = L(i,s) * (mi(i) + Qtot(row,i));
            CNtot = CNtot + C(i,s);
        end
        den = Z(s) + CNtot;
        X(s) = n(s) / den;
        Xall(currentpop,s) = X(s);
        for i = 1:M
            Q(i,s) = X(s) * C(i,s);
            Qcls(currentpop,i,s) = Q(i,s);
            Qtot(currentpop,i) = Qtot(currentpop,i) + Q(i,s);
        end
        s = s + 1;
    end

    % ---- moment step at population n ------------------------------------
    % W(k,j;t,j|n) = Q(j,t|n)*(Q(j,k|n-e_t) - Q(j,k|n)) + [k==t]*Q(j,t|n)
    %              + X(t|n)*L(j,t)*Ssum(j,k|n-e_t)
    for j = 1:M
        for k = 1:R
            Qjk = Qcls(currentpop,j,k);
            sk = 0;
            for t = 1:R
                Qjt = Qcls(currentpop,j,t);
                wkt = Qjt * (Qcls(rows(t),j,k) - Qjk);
                if k == t
                    wkt = wkt + Qjt;
                end
                wkt = wkt + Xall(currentpop,t) * L(j,t) * Ssum(rows(t),j,k);
                QCov(j,k,t) = wkt;
                sk = sk + wkt;
            end
            Ssum(currentpop,j,k) = sk;
        end
    end

    % ---- odometer advance ------------------------------------------------
    s = R;
    while (s>0 && n(s)==N(s)) || s>firstnonempty
        s = s - 1;
    end
    if s == 0
        break;
    end
    n(s) = n(s) + 1;
    s = s + 1;
    while s <= R
        n(s) = 0;
        s = s + 1;
    end
    ctr = ctr - 1;
    currentpop = currentpop + 1;
end

% utilization
for i = 1:M
    for r = 1:R
        U(i,r) = X(r) * L(i,r);
    end
end

mom = pack(X,Q,U,C,QCov);
end

% =========================================================================
function mom = pack(X,Q,U,C,QCov)
[M,R] = size(Q);
mom.X = X; mom.Q = Q; mom.U = U; mom.R = C;
% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
QCovRaw = QCov;
QCov = (QCov + permute(QCov,[1 3 2])) / 2;
mom.QCovAsym = max(max(max(abs(QCovRaw - permute(QCovRaw,[1 3 2])))));
if isempty(mom.QCovAsym)
    mom.QCovAsym = 0;
end
QVar = zeros(M,R);
QTotVar = zeros(M,1);
for i = 1:M
    for r = 1:R
        QVar(i,r) = QCov(i,r,r);
    end
    QTotVar(i) = sum(sum(QCov(i,:,:)));
end
mom.QCov = QCov;
mom.QVar = QVar;
mom.QTotVar = QTotVar;
end
