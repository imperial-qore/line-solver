function pfqn_sens_mvaldmx_validate()
%{
%{
 % @file pfqn_sens_mvaldmx_validate.m
 % @brief Validation harness for pfqn_sens_mvaldmx, the mixed load-dependent
 %        moment analysis of Akyildiz and Strelen (1991).
 %
 %        Five independent references, chosen so that every channel of the
 %        derivation is exercised by something that does not share its code:
 %
 %          A. pfqn_mvaldmx for the base measures X, Q, U, R. The primal must be
 %             reproduced entry by entry, otherwise the derivative is of the
 %             wrong function.
 %          B. central finite differences of pfqn_mvaldmx with respect to the
 %             demand-scaling parameter y(j,s). This checks the differentiated
 %             recursion itself, including the load-dependent channel
 %             dEC/dLo and the open-class channel dLo/dy of eq. (21), but does
 %             not check the identity Cov = d nbar / dy.
 %          C. pfqn_sens_mva in the closed load-independent limit. This checks
 %             the identity against the independently validated de Souza e Silva
 %             and Muntz recursion.
 %          D. brute-force enumeration of the product-form equilibrium
 %             distribution. Closed load-dependent models are enumerated exactly;
 %             mixed models are enumerated with the open populations truncated,
 %             which converges geometrically and is therefore checked at a
 %             looser tolerance. This is the only check that closes the loop on
 %             the identity in the mixed load-dependent case.
 %          E. symmetry of QCovFull. Cov[n(i,r),n(j,s)] and Cov[n(j,s),n(i,r)]
 %             are computed by differentiating two different classes'
 %             equations, so their agreement is a nontrivial structural check.
%}
%}
rng(1);
tolMva  = 1e-12;
tolFd   = 1e-6;
tolMom  = 1e-9;
tolBrt  = 1e-8;    % exact enumeration, closed load-dependent
tolBrtT = 5e-5;    % truncated enumeration, mixed
tolSym  = 1e-8;

errMva = 0; errFd = 0; errMom = 0; errBrt = 0; errBrtT = 0; errSym = 0;
nFd = 0; nMom = 0; nBrt = 0; nBrtT = 0;

% =====================================================================
% A/B/E on random mixed load-dependent models
% =====================================================================
for trial = 1:12
    M = randi([1 2]);
    Ropen = randi([0 1]);
    D = 0.2 + 0.6*rand(M,1+Ropen);
    R = 1 + Ropen;
    N = zeros(1,R);
    N(1) = randi([1 3]);              % closed class
    lambda = zeros(1,R);
    if Ropen == 1
        N(2) = Inf;                   % open class
        lambda(2) = 0.05 + 0.15*rand; % keep the station loads modest
    end
    Z = zeros(1,R);
    Z(1) = 0.5*rand;
    NCtot = sum(N(isfinite(N)));
    % limited load dependence: rates grow up to level b then saturate
    b = randi([1 3]);
    mu = zeros(M,max(NCtot,1));
    for i = 1:M
        for n = 1:size(mu,2)
            mu(i,n) = min(n,b) * (0.8 + 0.4*rand);
        end
    end
    % keep the geometric tail of the limited load dependence stable
    Lo = zeros(M,1);
    for i = 1:M
        Lo(i) = lambda*D(i,:)';
    end
    if any(Lo ./ mu(:,end) > 0.6)
        continue;
    end

    mom = pfqn_sens_mvaldmx(lambda,D,N,Z,mu,ones(M,1));

    % ---- A. base measures ------------------------------------------
    [XN,QN,UN,CN] = pfqn_mvaldmx(lambda,D,N,Z,mu,ones(M,1));
    errMva = max([errMva, relerr(mom.X,XN), relerr(mom.Q,QN), ...
                  relerr(mom.U,UN), relerr(mom.R,CN)]);

    % ---- E. symmetry -----------------------------------------------
    errSym = max(errSym, mom.QCovAsym);

    % ---- B. finite differences -------------------------------------
    h = 1e-6;
    for j = 1:M
        for s = 1:R
            if D(j,s) <= 0, continue; end
            Dp = D; Dp(j,s) = D(j,s)*(1+h);
            Dm = D; Dm(j,s) = D(j,s)*(1-h);
            [~,QNp] = pfqn_mvaldmx(lambda,Dp,N,Z,mu,ones(M,1));
            [~,QNm] = pfqn_mvaldmx(lambda,Dm,N,Z,mu,ones(M,1));
            fd = (QNp - QNm) / (2*h);   % d nbar / dy at y=1
            an = zeros(M,R);
            for i = 1:M
                for r = 1:R
                    an(i,r) = mom.QCovFull(i,r,j,s);
                end
            end
            errFd = max(errFd, relerr(an,fd));
            nFd = nFd + 1;
        end
    end
end

% =====================================================================
% C. closed load-independent limit against pfqn_sens_mva
% =====================================================================
for trial = 1:12
    M = randi([1 3]);
    R = randi([1 2]);
    D = 0.2 + rand(M,R);
    N = randi([1 3],1,R);
    Z = 0.4*rand(1,R);
    lambda = zeros(1,R);
    mu = ones(M,sum(N));
    mom = pfqn_sens_mvaldmx(lambda,D,N,Z,mu,ones(M,1));
    ref = pfqn_sens_mva(D,N,Z);
    errMom = max([errMom, relerr(mom.Q,ref.Q), relerr(mom.QCov,ref.QCov), ...
                  relerr(mom.QVar,ref.QVar), relerr(mom.QTotVar,ref.QTotVar)]);
    nMom = nMom + 1;
end

% =====================================================================
% D. brute-force product form
% =====================================================================
% D1. closed load-dependent, exact enumeration
for trial = 1:10
    M = 2; R = 1;
    D = 0.3 + 0.5*rand(M,R);
    N = randi([2 4],1,R);
    Z = 0.3*rand(1,R);
    lambda = 0;
    b = randi([2 3]);
    mu = zeros(M,sum(N));
    for i = 1:M
        for n = 1:size(mu,2)
            mu(i,n) = min(n,b) * (0.8 + 0.4*rand);
        end
    end
    mom = pfqn_sens_mvaldmx(lambda,D,N,Z,mu,ones(M,1));
    [Qb,QCovb] = brute_ldmx(lambda,D,N,Z,mu,0);
    errBrt = max([errBrt, relerr(mom.Q,Qb), relerr(mom.QCov,QCovb)]);
    nBrt = nBrt + 1;
end

% D2. mixed load-dependent, truncated enumeration
for trial = 1:6
    M = 2; R = 2;
    D = 0.3 + 0.4*rand(M,R);
    N = [randi([1 2]), Inf];
    Z = [0.3*rand, 0];
    lambda = [0, 0.05 + 0.1*rand];
    b = randi([1 2]);
    mu = zeros(M,sum(N(isfinite(N))));
    for i = 1:M
        for n = 1:size(mu,2)
            mu(i,n) = min(n,b) * (1.0 + 0.3*rand);
        end
    end
    Lo = zeros(M,1);
    for i = 1:M
        Lo(i) = lambda*D(i,:)';
    end
    if any(Lo ./ mu(:,end) > 0.4)
        continue;
    end
    mom = pfqn_sens_mvaldmx(lambda,D,N,Z,mu,ones(M,1));
    [Qb,QCovb] = brute_ldmx(lambda,D,N,Z,mu,60);
    errBrtT = max([errBrtT, relerr(mom.Q,Qb), relerr(mom.QCov,QCovb)]);
    nBrtT = nBrtT + 1;
end

fprintf('\n=== pfqn_sens_mvaldmx validation (max relative error) ===\n');
fprintf('  A. pfqn_mvaldmx base measures            : %.3e  (tol %.1e)\n', errMva, tolMva);
fprintf('  B. finite differences      (%3d params)  : %.3e  (tol %.1e)\n', nFd, errFd, tolFd);
fprintf('  C. pfqn_sens_mva closed LI  (%3d models)  : %.3e  (tol %.1e)\n', nMom, errMom, tolMom);
fprintf('  D1. brute force closed LD  (%3d models)  : %.3e  (tol %.1e)\n', nBrt, errBrt, tolBrt);
fprintf('  D2. brute force mixed LD   (%3d models)  : %.3e  (tol %.1e)\n', nBrtT, errBrtT, tolBrtT);
fprintf('  E. QCovFull symmetry (raw)               : %.3e  (tol %.1e)\n', errSym, tolSym);

ok = errMva <= tolMva && errFd <= tolFd && errMom <= tolMom && ...
     errBrt <= tolBrt && errBrtT <= tolBrtT && errSym <= tolSym;
if ~ok
    error('pfqn_sens_mvaldmx_validate:mismatch','one or more checks exceeded tolerance');
end
fprintf('  ALL CHECKS PASSED\n');
end

% =========================================================================
function e = relerr(a,b)
a = a(:); b = b(:);
d = abs(a-b);
scale = max(1, max(abs(a),abs(b)));
e = max(d ./ scale);
if isempty(e)
    e = 0;
end
end

% =========================================================================
function [Q,QCov] = brute_ldmx(lambda,D,N,Z,mu,Kopen)
% Exact moments by enumerating the mixed load-dependent product form
%   p(n) ~ prod_i [ n_i! prod_r a(i,r)^n(i,r)/n(i,r)! prod_{j=1}^{n_i} 1/mu(i,j) ]
%          * prod_{closed c} Z(c)^n(0,c)/n(0,c)!
% with a(i,r) = D(i,r) for a closed class and a(i,r) = lambda(r)*D(i,r) for an
% open class, n_i the total population at station i, and the closed classes
% constrained to sum to N. Open classes are truncated at Kopen jobs per station.
[M,R] = size(D);
openClasses = find(isinf(N));
closedClasses = setdiff(1:R, openClasses);

a = zeros(M,R);
for r = 1:R
    if isinf(N(r))
        a(:,r) = lambda(r) * D(:,r);
    else
        a(:,r) = D(:,r);
    end
end

% per-class allocation lists
alloc = cell(1,R);
for r = 1:R
    if isinf(N(r))
        alloc{r} = compositions_leq(Kopen,M);   % open: 0..Kopen total, free
    else
        alloc{r} = compositions_leq(N(r),M);    % closed: remainder in the delay
    end
end

idx = ones(1,R);
tot = 1;
for r = 1:R
    tot = tot * size(alloc{r},1);
end
W = zeros(tot,1);
NIR = zeros(tot,M*R);
k = 0;
while true
    k = k + 1;
    nir = zeros(M,R);
    for r = 1:R
        nir(:,r) = alloc{r}(idx(r),:)';
    end
    lw = 0;
    ok = true;
    for i = 1:M
        ni = sum(nir(i,:));
        lw = lw + gammaln(ni+1);
        for j = 1:ni
            lw = lw - log(mu_at(mu,i,j));
        end
        for r = 1:R
            if nir(i,r) > 0
                if a(i,r) <= 0
                    ok = false; break;
                end
                lw = lw + nir(i,r)*log(a(i,r)) - gammaln(nir(i,r)+1);
            end
        end
        if ~ok, break; end
    end
    if ok
        for ci = 1:numel(closedClasses)
            c = closedClasses(ci);
            n0c = N(c) - sum(nir(:,c));
            if n0c > 0
                if Z(c) <= 0
                    ok = false; break;
                end
                lw = lw + n0c*log(Z(c)) - gammaln(n0c+1);
            end
        end
    end
    if ok
        W(k) = exp(lw);
    end
    NIR(k,:) = reshape(nir',1,[]);
    r = R;
    while r >= 1
        idx(r) = idx(r) + 1;
        if idx(r) <= size(alloc{r},1)
            break;
        end
        idx(r) = 1;
        r = r - 1;
    end
    if r == 0
        break;
    end
end
W = W / sum(W);

Q = zeros(M,R);
for k = 1:size(NIR,1)
    Q = Q + W(k)*reshape(NIR(k,:),R,M)';
end
QCov = zeros(M,R,R);
for i = 1:M
    for r = 1:R
        for s = 1:R
            m2 = 0;
            for k = 1:size(NIR,1)
                nir = reshape(NIR(k,:),R,M)';
                m2 = m2 + W(k)*nir(i,r)*nir(i,s);
            end
            QCov(i,r,s) = m2 - Q(i,r)*Q(i,s);
        end
    end
end
end

% =========================================================================
function v = mu_at(mu,i,j)
% Limited load dependence: the rate saturates at its last tabulated value.
if j <= size(mu,2)
    v = mu(i,j);
else
    v = mu(i,end);
end
end

% =========================================================================
function C = compositions_leq(n,M)
% All nonnegative integer vectors of length M summing to at most n.
if M == 1
    C = (0:n)';
    return;
end
C = zeros(0,M);
for first = 0:n
    sub = compositions_leq(n-first,M-1);
    C = [C; [repmat(first,size(sub,1),1), sub]]; %#ok<AGROW>
end
end
