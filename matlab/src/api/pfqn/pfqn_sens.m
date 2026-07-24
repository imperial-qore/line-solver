%{
%{
 % @file pfqn_sens.m
 % @brief Exact analytic performance sensitivities for closed product-form
 %        queueing networks. Dispatches to a CoMoM-backed kernel for the
 %        single-station repairman model and to differentiated MVA otherwise.
%}
%}

function sens = pfqn_sens(L,N,Z,mi)
%{
%{
 % @brief Exact derivatives of the mean performance measures {X,Q,U,R} of a
 %        closed product-form (BCMP) queueing network with respect to the
 %        service demands L(i,r) and the think times Z(r). The derivatives are
 %        analytic (exact to machine precision), not finite differences.
 %
 %        Two exact kernels are dispatched transparently (local subfunctions):
 %          - sens_comom : Class-Oriented Method of Moments. Selected for the
 %            repairman model (a single single-server queue, M=1, plus a
 %            think-time delay with sum(Z)>0 and strictly positive demands).
 %            Polynomial in the number of classes R via normalizing-constant
 %            moment relations (queue-length covariances / replica moments).
 %          - sens_mva : forward-mode differentiation of the exact
 %            Reiser-Lavenberg MVA recursion. Used for every other model.
 %        Both kernels return the identical struct layout, so the choice is an
 %        internal optimization invisible to callers.
 %
 %        Reference: Z. Liu and P. Nain, "Sensitivity Results in Open, Closed
 %        and Mixed Product-Form Queueing Networks", INRIA RR-1144, 1989;
 %        X.-R. Cao and D.-J. Ma, "Performance sensitivity formulae, algorithms
 %        and estimates for closed queueing networks with exponential servers",
 %        Performance Evaluation 26:181-199, 1996; G. Casale, "CoMoM: Efficient
 %        Class-Oriented Evaluation of Multiclass Performance Models", IEEE TSE
 %        2011.
 %
 % @fn pfqn_sens(L, N, Z, mi)
 % @param L  Service demand matrix (M x R), L(i,r) = visits_ir / rate_ir.
 % @param N  Population vector (1 x R).
 % @param Z  Think time vector (1 x R). Default: zeros.
 % @param mi (Optional) Station residence multiplicity (1 x M). Default: ones.
 % @return sens A struct with base measures and their Jacobians:
 %   .X (1 x R), .Q (M x R), .U (M x R), .R (M x R)  base MVA measures
 %       (X system throughput per class, Q mean queue length, U utilization,
 %        R residence time per visit-chain, i.e. CN of pfqn_mva).
 %   .params  1 x P struct array describing each differentiation parameter,
 %            fields .type ('L' or 'Z'), .station (i, 0 for Z), .class (r).
 %   .dX (R x P), .dQ (M x R x P), .dU (M x R x P), .dR (M x R x P)
 %       derivatives of each base measure w.r.t. parameter p. For a 'L'
 %       parameter at (i,r) the derivative is d(.)/dL(i,r); for a 'Z'
 %       parameter at class r it is d(.)/dZ(r).
 %   .QCov (M x R x M x R)  QCov(i,r,j,s) = Cov[n(i,r),n(j,s)], the exact
 %       queue-length covariance, a by-product of the Jacobian.
 %   .QVar (M x R)          QVar(i,r) = Var[n(i,r)].
 %   .QTotVar (M x 1)       QTotVar(i) = Var[sum_r n(i,r)].
 %   .QCovAsym (scalar)     roundoff-level residual of the moment recursion,
 %       see pfqn_sens_mva.
 %
 % Notes:
 %  - Mirrors pfqn_mva(L,N,Z,mi) exactly for the base measures (single-server
 %    or residence-multiplicity mi stations plus an infinite-server delay Z).
 %  - Derivatives w.r.t. a service rate mu(i,r) follow by the chain rule
 %    d(.)/dmu(i,r) = -(L(i,r)/mu(i,r)) * d(.)/dL(i,r).
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

% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
populated = N > 0;
useComom = (M == 1) && all(mi == 1) && any(populated) && ...
           all(Z(populated) > GlobalConstants.FineTol) && ...
           all(L(1, populated) > GlobalConstants.FineTol);
if useComom
    sens = sens_comom(L,N,Z);
else
    sens = sens_mva(L,N,Z,mi);
end

% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
mom = pfqn_sens_mva(L,N,Z,mi);
[sens.QCov, sens.QVar] = queue_cov(L,sens,mom);
sens.QTotVar = mom.QTotVar;
sens.QCovAsym = mom.QCovAsym;
end

% =========================================================================
function [QCov,QVar] = queue_cov(L,sens,mom)
% QCov(i,r,j,s) = Cov[n_{i,r},n_{j,s}] = D_{j,s} dQ_{i,r}/dD_{j,s}.
[M,R] = size(L);
pL = zeros(M,R);
for p = 1:numel(sens.params)
    pr = sens.params(p);
    if pr.type == 'L'
        pL(pr.station,pr.class) = p;
    end
end
QCov = zeros(M,R,M,R);
for i = 1:M
    for r = 1:R
        for j = 1:M
            for s = 1:R
                if i == j
                    QCov(i,r,j,s) = mom.QCov(i,r,s);
                else
                    p = pL(j,s);
                    if p > 0
                        QCov(i,r,j,s) = L(j,s) * sens.dQ(i,r,p);
                    end
                end
            end
        end
    end
end
QVar = mom.QVar;
end

% =========================================================================
% CoMoM-backed kernel (M=1 repairman model)
% =========================================================================
function sens = sens_comom(L,N,Z)
% Exact derivatives of {X,Q,U,R} for a single single-server queue plus a
% think-time delay, using pfqn_comomrm to evaluate normalizing constants of the
% base model and of its 2- and 3-replica extensions. The queue-length moment
% relation
%     D_{1,s} dQ_{1,r}/dD_{1,s} = Cov[n_{1,r},n_{1,s}]
%                               = Q_{1,r}(delta_{rs}+2 Q^{+1}_{1,s}(N-1_r)-Q_{1,s})
% yields the demand Jacobian; think-time derivatives use the exact identity
% d log G_m(n)/dZ_s = G_m(n-1_s)/G_m(n) (valid for any Z_s>=0). Q^{+1}_{1,s}(n)
% is the class-s queue at one of two identical replicas of the station:
%     Q^{+1}_{1,s}(n) = D_{1,s} G_3(n-1_s)/G_2(n),
% obtained from CoMoM's replication factor m.
[M,R] = size(L); %#ok<ASGLU>
N = ceil(N(:)');
Z = Z(:)';
D = L(1,:);

% parameter list: L(1,r) first, then Z(r) (mirrors sens_mva ordering)
P  = 2*R;
pL = 1:R;
pZ = R + (1:R);
paramType = cell(1,P); paramStation = zeros(1,P); paramClass = zeros(1,P);
for r = 1:R
    paramType{pL(r)} = 'L'; paramStation(pL(r)) = 1; paramClass(pL(r)) = r;
    paramType{pZ(r)} = 'Z'; paramStation(pZ(r)) = 0; paramClass(pZ(r)) = r;
end

X = zeros(1,R); Q = zeros(1,R); U = zeros(1,R); C = zeros(1,R);
dX = zeros(R,P); dQ = zeros(1,R,P); dU = zeros(1,R,P); dC = zeros(1,R,P);

if ~any(N > 0)
    sens = pack(X,Q,U,C,dX,dQ,dU,dC,paramType,paramStation,paramClass);
    return;
end

cache = containers.Map('KeyType','char','ValueType','double');
    function lg = lgm(m,n)
        % memoized log normalizing constant of the m-replica model; zero
        % population classes are stripped (they leave the NC unchanged) so
        % pfqn_comomrm does not index a stale class count after its sanitize.
        n = round(n);
        nz = n>0;
        if ~any(nz), lg = 0; return; end
        key = [num2str(m) '|' sprintf('%d,',n)];
        if isKey(cache,key)
            lg = cache(key); return;
        end
        lg = pfqn_comomrm(D(nz),n(nz),Z(nz),m,GlobalConstants.FineTol);
        cache(key) = lg;
    end
    function e = ei(s)
        e = zeros(1,R); e(s) = 1;
    end
    function q = qmean(n,s)
        % mean class-s queue at population n (single station, m=1 model)
        if n(s) < 1, q = 0; return; end
        q = D(s) * exp(lgm(2,n-ei(s)) - lgm(1,n));
    end
    function x = xput(m,n,s)
        % class-s throughput X_s^{(m)}(n) = G_m(n-1_s)/G_m(n) in m-replica model
        if n(s) < 1, x = 0; return; end
        x = exp(lgm(m,n-ei(s)) - lgm(m,n));
    end
    function q = qplus(n,s)
        % Q^{+1}_{1,s}(n): class-s queue at one replica of the doubled station
        if n(s) < 1, q = 0; return; end
        q = D(s) * exp(lgm(3,n-ei(s)) - lgm(2,n));
    end

% base measures
for r = 1:R
    if N(r) < 1, continue; end
    X(r) = xput(1,N,r);
end
for s = 1:R
    Q(s) = qmean(N,s);
end
for r = 1:R
    U(r) = X(r) * D(r);
    if X(r) > 0, C(r) = Q(r) / X(r); end
end

% Jacobian
for r = 1:R
    if N(r) < 1, continue; end   % empty class: X=Q=0, all derivatives 0
    Nr = N - ei(r);
    Xr = X(r); Qr = Q(r);
    for s = 1:R
        % L(1,s) parameter
        p = pL(s);
        Vrs  = Qr * ((r==s) + 2*qplus(Nr,s) - Q(s));  % Cov[n_r,n_s]
        dQ_L = Vrs / D(s);
        dX_L = Xr * (qmean(Nr,s) - Q(s)) / D(s);
        dU_L = dX_L * D(r) + Xr * (r==s);
        dQ(1,r,p) = dQ_L;
        dX(r,p)   = dX_L;
        dU(1,r,p) = dU_L;
        if Xr > 0
            dC(1,r,p) = (dQ_L*Xr - Qr*dX_L) / Xr^2;
        end

        % Z(s) parameter: d log G_m(n)/dZ_s = G_m(n-1_s)/G_m(n) (exact, any Z_s>=0)
        p = pZ(s);
        dX_Z = Xr * (xput(1,Nr,s) - X(s));
        dQ_Z = Qr * (xput(2,Nr,s) - X(s));
        dU_Z = dX_Z * D(r);
        dQ(1,r,p) = dQ_Z;
        dX(r,p)   = dX_Z;
        dU(1,r,p) = dU_Z;
        if Xr > 0
            dC(1,r,p) = (dQ_Z*Xr - Qr*dX_Z) / Xr^2;
        end
    end
end

sens = pack(X,Q,U,C,dX,dQ,dU,dC,paramType,paramStation,paramClass);
end

% =========================================================================
% Differentiated-MVA kernel (general model)
% =========================================================================
function sens = sens_mva(L,N,Z,mi)
% Forward-mode differentiation of the exact Reiser-Lavenberg MVA recursion.
[M,R] = size(L);
N = ceil(N(:)');
Z = Z(:)';
mi = mi(:)';

% parameter list: all L(i,r), then all Z(r)
P = M*R + R;
paramType = cell(1,P);
paramStation = zeros(1,P);
paramClass = zeros(1,P);
pL = zeros(M,R);
pZ = zeros(1,R);
p = 0;
for i = 1:M
    for r = 1:R
        p = p + 1;
        paramType{p} = 'L'; paramStation(p) = i; paramClass(p) = r;
        pL(i,r) = p;
    end
end
for r = 1:R
    p = p + 1;
    paramType{p} = 'Z'; paramStation(p) = 0; paramClass(p) = r;
    pZ(r) = p;
end

X = zeros(1,R); Q = zeros(M,R); U = zeros(M,R); C = zeros(M,R);
dX = zeros(R,P); dQ = zeros(M,R,P); dU = zeros(M,R,P); dC = zeros(M,R,P);

if ~any(N > 0)
    sens = pack(X,Q,U,C,dX,dQ,dU,dC,paramType,paramStation,paramClass);
    return;
end

% population-lattice odometer, identical to pfqn_mva
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
Qtot  = zeros(totpop,M);
Qtotd = zeros(totpop,M,P);
currentpop = 2;
n = zeros(1,R);
n(firstnonempty) = 1;
C_d_is = zeros(M,P);

while ctr
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
        row = 1 + pos;
        CNtot = 0;
        CNtotd = zeros(1,P);
        for i = 1:M
            base = mi(i) + Qtot(row,i);
            C(i,s) = L(i,s) * base;
            Cd = L(i,s) * reshape(Qtotd(row,i,:),1,P);
            Cd(pL(i,s)) = Cd(pL(i,s)) + base;
            C_d_is(i,:) = Cd;
            CNtot = CNtot + C(i,s);
            CNtotd = CNtotd + Cd;
        end
        den = Z(s) + CNtot;
        X(s) = n(s) / den;
        Xd = -n(s) * CNtotd / den^2;
        Xd(pZ(s)) = Xd(pZ(s)) - n(s) / den^2;
        dX(s,:) = Xd;
        for i = 1:M
            Q(i,s) = X(s) * C(i,s);
            Qd = Xd * C(i,s) + X(s) * C_d_is(i,:);
            dQ(i,s,:) = reshape(Qd,1,1,P);
            dC(i,s,:) = reshape(C_d_is(i,:),1,1,P);
            Qtot(currentpop,i)  = Qtot(currentpop,i) + Q(i,s);
            Qtotd(currentpop,i,:) = Qtotd(currentpop,i,:) + reshape(Qd,1,1,P);
        end
        s = s + 1;
    end
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

% utilization and its derivatives
for i = 1:M
    for r = 1:R
        U(i,r) = X(r) * L(i,r);
        Ud = dX(r,:) * L(i,r);
        Ud(pL(i,r)) = Ud(pL(i,r)) + X(r);
        dU(i,r,:) = reshape(Ud,1,1,P);
    end
end

sens = pack(X,Q,U,C,dX,dQ,dU,dC,paramType,paramStation,paramClass);
end

% =========================================================================
function sens = pack(X,Q,U,C,dX,dQ,dU,dC,paramType,paramStation,paramClass)
P = numel(paramType);
params = struct('type',cell(1,P),'station',cell(1,P),'class',cell(1,P));
for p = 1:P
    params(p).type = paramType{p};
    params(p).station = paramStation(p);
    params(p).class = paramClass(p);
end
sens.X = X; sens.Q = Q; sens.U = U; sens.R = C;
sens.params = params;
sens.dX = dX; sens.dQ = dQ; sens.dU = dU; sens.dR = dC;
end
