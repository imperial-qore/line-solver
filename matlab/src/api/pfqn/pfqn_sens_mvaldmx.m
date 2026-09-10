%{
%{
 % @file pfqn_sens_mvaldmx.m
 % @brief Exact queue-length variances and covariances for mixed open/closed
 %        product-form queueing networks with limited load dependence.
%}
%}

function mom = pfqn_sens_mvaldmx(lambda,D,N,Z,mu,S)
%{
%{
 % @brief Exact second moments (variances and covariances) of the queue lengths
 %        of a mixed open/closed product-form queueing network with limited
 %        load-dependent service rates. This is the load-dependent and mixed
 %        counterpart of pfqn_sens_mva, which is restricted to closed
 %        load-independent models.
 %
 %        The method is the moment analysis of Akyildiz and Strelen. Their
 %        Theorem 1, equation (11), states that multiplying a queue-length
 %        moment by one further factor Q_jT costs one derivative with respect to
 %        a parameter y_j that scales the service demands s_ir of the classes
 %        r in T at station j:
 %
 %          E[Q_jT^k ...] = d/dy_j E[Q_jT^(k-1) ...]|_{y_j=1}
 %                          + nbar_jT E[Q_jT^(k-1) ...]
 %
 %        Taking k=2 and T={s} gives the second moment, hence
 %
 %          Cov[n(i,r),n(j,s)] = d nbar(i,r) / dy_(j,s) |_{y=1}
 %
 %        which is evaluated here by forward-mode differentiation of the mixed
 %        load-dependent MVA of Bruell-Balbo-Afshari, i.e. of exactly the
 %        recursion implemented by pfqn_mvaldmx. The differentiated equations
 %        are (13) for the residence times, (15)-(17) for the conditional
 %        marginal probabilities, (18) for the throughputs, (19) and (24)-(31)
 %        for the effective capacities (delegated to pfqn_sens_ldmx_ec), (32) for
 %        the closed-class queue lengths and (33) for the open-class ones.
 %
 %        Because a demand-scaling parameter perturbs the whole network through
 %        the closed-class throughputs, the derivatives must be propagated for
 %        every parameter, so the cross-station covariances come out at no extra
 %        cost and are returned in .QCovFull. This is unlike pfqn_sens_mva, whose
 %        cheaper same-station recursion cannot reach them.
 %
 %        Reference: I. F. Akyildiz and J. C. Strelen, "Moment Analysis for
 %        Load-Dependent Mixed Product Form Queueing Networks", IEEE Trans.
 %        Communications 39(6):828-832, 1991. The closed load-independent case
 %        reduces to E. de Souza e Silva and R. R. Muntz, IEEE Trans. Computers
 %        37(9):1125-1129, 1988, which pfqn_sens_mva implements directly.
 %
 % @fn pfqn_sens_mvaldmx(lambda, D, N, Z, mu, S)
 % @param lambda Arrival rate vector (1 x R). Must be zero on closed classes.
 % @param D Service demand matrix (M x R).
 % @param N Population vector (1 x R). Inf entries denote open classes.
 % @param Z Think time vector (1 x R).
 % @param mu Load-dependent rate matrix (M x sum(N)), limited load dependence.
 % @param S Number of servers per station (M x 1). Accepted for signature
 %        compatibility with pfqn_mvaldmx, which likewise does not read it: the
 %        multiserver behaviour is carried entirely by the rates mu.
 % @return mom A struct with the base measures and their second moments:
 %   .X (1 x R), .Q (M x R), .U (M x R), .R (M x R)  base measures, identical
 %       to pfqn_mvaldmx(lambda,D,N,Z,mu,S).
 %   .QCov (M x R x R)      QCov(i,r,s) = Cov[n(i,r),n(i,s)], same-station.
 %   .QCovFull (M x R x M x R)  QCovFull(i,r,j,s) = Cov[n(i,r),n(j,s)].
 %   .QVar (M x R)          QVar(i,r) = Var[n(i,r)].
 %   .QTotVar (M x 1)       QTotVar(i) = Var[sum_r n(i,r)].
 %   .QCovAsym (scalar)     max |QCovFull(i,r,j,s) - QCovFull(j,s,i,r)| before
 %       symmetrization. The two entries are produced by differentiating two
 %       different classes' equations, so this residual is an independent check
 %       of the recursion and should sit at roundoff.
 %
 % Notes:
 %  - Open classes are supported: an open class contributes to the load Lo(i)
 %    that drives the effective capacities, and equation (21) supplies the
 %    corresponding dLo/dy.
 %  - The moments of an open class are those of its queue length at a station,
 %    which is finite even though its population is infinite.
%}
%}
if nargin<5
    mu=ones(size(D,1),sum(N(isfinite(N))));
    S=ones(size(D,1),1); %#ok<NASGU>
end
if nargin<6
    S=ones(size(D,1),1); %#ok<NASGU>
end
if size(mu,2) < sum(N(isfinite(N)))
    line_error(mfilename,'PFQN_SENS_MVALDMX requires to specify the load-dependent rates with one job more than the maximum closed population.');
end
if any(N(find(lambda))>0 & isfinite(N(find(lambda)))) %#ok<*FNDSB>
    line_error(mfilename,'Arrival rate cannot be specified on closed classes.');
end
[M,R] = size(D);
lambda = lambda(:)';
N = N(:)';
Z = Z(:)';
openClasses = find(isinf(N));
closedClasses = setdiff(1:length(N), openClasses);
C = length(closedClasses);
if C == 0
    line_error(mfilename,'pfqn_sens_mvaldmx requires at least one closed class; use the open-class formulas directly otherwise.');
end

mu = [mu, mu(:,size(mu,2))]; % up to sum(N)+1, limited load dependence
[EC,E,Eprime,~,dEC_dLo] = pfqn_sens_ldmx_ec(lambda,D,mu);

Dc = D(:,closedClasses);
Nc = N(closedClasses);
Zc = Z(closedClasses);
NCtot = sum(Nc);

% ---- parameter list: y(j,r) multiplies the demand D(j,r) -----------------
P = M*R;
pidx = zeros(M,R);
pj = zeros(1,P); pr = zeros(1,P);
p = 0;
for j = 1:M
    for r = 1:R
        p = p + 1;
        pidx(j,r) = p; pj(p) = j; pr(p) = r;
    end
end
% eq. (21): Lo(i) = sum_o lambda(o)*D(i,o), so only an open-class parameter at
% station i perturbs Lo(i), and no parameter perturbs Lo at another station
dLo = zeros(M,P);
for p = 1:P
    if isinf(N(pr(p)))
        dLo(pj(p),p) = lambda(pr(p)) * D(pj(p),pr(p));
    end
end

% ---- population recursion ------------------------------------------------
prods = zeros(1,C);
for r = 1:C
    prods(r) = prod(Nc(1:r-1)+1);
end
NT = prod(1+Nc);
Pc  = zeros(M,1+NCtot,NT);
dPc = zeros(M,1+NCtot,NT,P);
x   = zeros(C,NT);
dx  = zeros(C,NT,P);
w   = zeros(M,C,NT);
dw  = zeros(M,C,NT,P);

nvec = pprod(Nc);
for ist = 1:M
    Pc(ist, 1+0, hashpop(nvec,Nc,C,prods)) = 1.0;   % eq. (16)
end

while nvec >= 0
    hnvec = hashpop(nvec,Nc,C,prods);
    nc = sum(nvec);

    % ---- residence times, eq. (12) and its derivative eq. (13) ----------
    for ist = 1:M
        for c = 1:C
            if nvec(c) > 0
                hnvec_c = hashpop(oner(nvec,c),Nc,C,prods);
                cls = closedClasses(c);
                acc = 0;
                dacc = zeros(1,P);
                for n = 1:nc
                    Pprev = Pc(ist, 1+(n-1), hnvec_c);
                    acc = acc + n * EC(ist,n) * Pprev;
                    for q = 1:P
                        dacc(q) = dacc(q) + n * ( dEC_dLo(ist,n)*dLo(ist,q)*Pprev ...
                                                + EC(ist,n)*dPc(ist,1+(n-1),hnvec_c,q) );
                    end
                end
                w(ist,c,hnvec) = Dc(ist,c) * acc;
                for q = 1:P
                    dwq = Dc(ist,c) * dacc(q);
                    if pj(q) == ist && pr(q) == cls
                        dwq = dwq + Dc(ist,c) * acc;   % d(D*y)/dy = D
                    end
                    dw(ist,c,hnvec,q) = dwq;
                end
            end
        end
    end

    % ---- throughputs, eq. (18) ------------------------------------------
    for c = 1:C
        den = Zc(c) + sum(w(1:M,c,hnvec));
        x(c,hnvec) = nvec(c) / den;
        if nvec(c) > 0
            for q = 1:P
                sdw = 0;
                for ist = 1:M
                    sdw = sdw + dw(ist,c,hnvec,q);
                end
                dx(c,hnvec,q) = -nvec(c) / den^2 * sdw;
            end
        end
    end

    % ---- conditional marginal probabilities, eq. (14)-(15) --------------
    for ist = 1:M
        for n = 1:nc
            for c = 1:C
                if nvec(c) > 0
                    hnvec_c = hashpop(oner(nvec,c),Nc,C,prods);
                    cls = closedClasses(c);
                    Pprev = Pc(ist, 1+(n-1), hnvec_c);
                    Pc(ist, 1+n, hnvec) = Pc(ist, 1+n, hnvec) ...
                        + Dc(ist,c) * EC(ist,n) * x(c,hnvec) * Pprev;
                    for q = 1:P
                        dt = Dc(ist,c) * ( dEC_dLo(ist,n)*dLo(ist,q)*x(c,hnvec)*Pprev ...
                                         + EC(ist,n)*dx(c,hnvec,q)*Pprev ...
                                         + EC(ist,n)*x(c,hnvec)*dPc(ist,1+(n-1),hnvec_c,q) );
                        if pj(q) == ist && pr(q) == cls
                            dt = dt + Dc(ist,c) * EC(ist,n) * x(c,hnvec) * Pprev;
                        end
                        dPc(ist,1+n,hnvec,q) = dPc(ist,1+n,hnvec,q) + dt;
                    end
                end
            end
        end
        % eq. (17). The primal keeps pfqn_mvaldmx's max(eps,.) floor so that the
        % base measures agree entry by entry; the derivative is the exact
        % -sum of the derivatives, since the floor is a numerical guard and not
        % part of the model.
        Pc(ist, 1+0, hnvec) = max(eps, 1-sum(Pc(ist, 1+(1:nc), hnvec)));
        for q = 1:P
            dPc(ist,1+0,hnvec,q) = -sum(dPc(ist,1+(1:nc),hnvec,q));
        end
    end

    nvec = pprod(nvec, Nc);
end

% ---- measures and their derivatives at the full population ---------------
hnvec = hashpop(Nc,Nc,C,prods);
XN = zeros(1,R); QN = zeros(M,R); UN = zeros(M,R); CN = zeros(M,R);
dQN = zeros(M,R,P);

% closed classes, eq. (32)
for c = 1:C
    cls = closedClasses(c);
    XN(cls) = x(c,hnvec);
    hnvec_c = hashpop(oner(Nc,c),Nc,C,prods);
    for ist = 1:M
        CN(ist,cls) = w(ist,c,hnvec);
        QN(ist,cls) = XN(cls) * CN(ist,cls);
        for q = 1:P
            dQN(ist,cls,q) = dx(c,hnvec,q)*w(ist,c,hnvec) + x(c,hnvec)*dw(ist,c,hnvec,q);
        end
        uacc = 0;
        for n = 1:NCtot
            uacc = uacc + Dc(ist,c) * x(c,hnvec) * Eprime(ist,1+n-1) / E(ist,1+n-1) ...
                   * Pc(ist, 1+n-1, hnvec_c);
        end
        UN(ist,cls) = uacc;
    end
end

% open classes, eq. (33)
for ridx = 1:length(openClasses)
    r = openClasses(ridx);
    XN(r) = lambda(r);
    for ist = 1:M
        acc = 0;
        dacc = zeros(1,P);
        for n = 0:NCtot
            Pn = Pc(ist, 1+n, hnvec);
            acc = acc + (n+1) * EC(ist,n+1) * Pn;
            for q = 1:P
                dacc(q) = dacc(q) + (n+1) * ( dEC_dLo(ist,n+1)*dLo(ist,q)*Pn ...
                                            + EC(ist,n+1)*dPc(ist,1+n,hnvec,q) );
            end
        end
        QN(ist,r) = lambda(r) * D(ist,r) * acc;
        CN(ist,r) = QN(ist,r) / lambda(r);
        for q = 1:P
            dq = lambda(r) * D(ist,r) * dacc(q);
            if pj(q) == ist && pr(q) == r
                dq = dq + lambda(r) * D(ist,r) * acc;
            end
            dQN(ist,r,q) = dq;
        end
        uacc = 0;
        for n = 0:NCtot
            uacc = uacc + lambda(r) * Eprime(ist,1+n+1) / E(ist,1+n+1) * Pc(ist, 1+n, hnvec);
        end
        UN(ist,r) = uacc;
    end
end

% ---- moments -------------------------------------------------------------
% Cov[n(i,r),n(j,s)] = d nbar(i,r) / dy_(j,s)
QCovFull = zeros(M,R,M,R);
for i = 1:M
    for r = 1:R
        for j = 1:M
            for s = 1:R
                QCovFull(i,r,j,s) = dQN(i,r,pidx(j,s));
            end
        end
    end
end
QCovRaw = QCovFull;
QCovFull = (QCovFull + permute(QCovFull,[3 4 1 2])) / 2;
mom.QCovAsym = max(abs(QCovRaw(:) - reshape(permute(QCovRaw,[3 4 1 2]),[],1)));
if isempty(mom.QCovAsym)
    mom.QCovAsym = 0;
end

QCov = zeros(M,R,R);
QVar = zeros(M,R);
QTotVar = zeros(M,1);
for i = 1:M
    for r = 1:R
        for s = 1:R
            QCov(i,r,s) = QCovFull(i,r,i,s);
        end
        QVar(i,r) = QCov(i,r,r);
    end
    QTotVar(i) = sum(sum(QCov(i,:,:)));
end

mom.X = XN; mom.Q = QN; mom.U = UN; mom.R = CN;
mom.QCov = QCov;
mom.QCovFull = QCovFull;
mom.QVar = QVar;
mom.QTotVar = QTotVar;
end
