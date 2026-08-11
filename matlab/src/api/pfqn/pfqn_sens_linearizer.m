%{
%{
 % @file pfqn_sens_linearizer.m
 % @brief Approximate higher moments of the queue lengths of a closed
 %        product-form queueing network, by differentiating the Linearizer
 %        fixed point. Polynomial in the population, unlike the exact
 %        pfqn_sens_mom.
%}
%}

function res = pfqn_sens_linearizer(L,N,Z,tol,maxiter)
%{
%{
 % @brief Approximate moments E[Q_i], Var[Q_i], Cov[Q_i,Q_j], E[Q_i^2] and
 %        E[Q_i^3] of the per-station total queue lengths of a closed
 %        product-form queueing network, by the LINEARIZER-2 / LINEARIZER-3
 %        algorithms of the reference (Section 5).
 %
 %        Motivation. The exact moment analysis of pfqn_sens_mom evaluates the
 %        MVA recursion on the whole population lattice, so it costs
 %        O(prod(N+1)) and is unusable once the populations are large. The
 %        Linearizer replaces that lattice by a fixed point over a handful of
 %        populations, and the reference observes that the same trick applies to
 %        the derivatives: differentiate the Linearizer equations, append the
 %        differentiated equations to the originals, and iterate all of them
 %        together. This routine does that, carrying both the first and the
 %        second derivative, so it returns everything (3.2) needs, including the
 %        third moment. Carrying only the first derivative is the reference's
 %        LINEARIZER-2; carrying the second as well is its LINEARIZER-3.
 %
 %        The approximation. CORE (equations (5.1)-(5.2)) estimates the queue
 %        lengths at population n - 1_l from those at n by
 %
 %          v_i(l)                 = m_i^(n)(l) / n(l)
 %          m_i^(n-1_l')(l)        = (n - 1_l')_l * ( v_i(l) + delta_i(l',l) )
 %
 %        and substitutes them into the exact MVA equations. Setting the delta
 %        terms to zero gives Bard-Schweitzer; Linearizer instead estimates them
 %        from (5.3), delta_i^(N)(l',l) = v_i^(N-1_l')(l) - v_i^(N)(l), by
 %        running CORE at each of the N - 1_l populations, and holds them fixed
 %        across populations (the heuristic (5.4)). Differentiating (5.1)-(5.3)
 %        gives (5.5)-(5.8), which are carried alongside.
 %
 %        Accuracy. The reference reports, over 51 networks including 34 stress
 %        cases, relative errors below 2.1% on E[Q], 4.1% on E[Q^2] and 6.2% on
 %        E[Q^3]. pfqn_sens_linearizer_validate measures the error against the
 %        exact pfqn_sens_mom on models small enough for both, and asserts
 %        bands of that order rather than machine precision: this routine is an
 %        approximation and is expected to disagree with the exact answer.
 %
 %        Reference: J. C. Strelen, "Moment Analysis for Closed Queuing Networks
 %        and its Linearizer", Performance Evaluation 11:127-142, 1990,
 %        Section 5, equations (5.1)-(5.8), the CORE-2 and LINEARIZER-2
 %        algorithms. The delta definition follows the standard Chandy-Neuse
 %        Linearizer, v_i^(N-1_l')(l) - v_i^(N)(l), which is what LINE's
 %        pfqn_linearizer implements.
 %
 % @fn pfqn_sens_linearizer(L, N, Z, tol, maxiter)
 % @param L  Service demand matrix (M x R), L(i,r) = visits_ir / rate_ir.
 % @param N  Population vector (1 x R).
 % @param Z  Think time vector (1 x R). Default: zeros.
 % @param tol (Optional) Convergence tolerance of the CORE fixed point on the
 %        mean queue lengths. Default: the test of the reference,
 %        1/(4000+16*sum(n)), which is also applied to the variances at 1e-3.
 % @param maxiter (Optional) Maximum CORE iterations. Default: 200.
 % @return res A struct with:
 %   .X (1 x R), .Q (M x R), .U (M x R), .W (M x R)  approximate base measures.
 %   .m (M x 1)      approximate E[Q_i], the total queue length at station i.
 %   .dm (M x M)     dm(i,h) = x_h dm_i/dx_h, the scaled first derivative.
 %   .d2m (M x 1)    d2m(i) = x_i^2 d^2m_i/dx_i^2.
 %   .Var (M x 1), .Cov (M x M), .M2 (M x 1), .M3 (M x 1), .Skew (M x 1)
 %       the moments of (3.2), formed exactly as in pfqn_sens_mom but from the
 %       approximate derivatives.
 %   .CovAsym (scalar)  raw asymmetry of Cov before symmetrization. Unlike the
 %       exact routines, this is NOT expected to sit at roundoff: the Linearizer
 %       fixed point does not enforce the symmetry that the product form
 %       guarantees, so this is a useful measure of the approximation error.
 %   .iter (scalar)  total CORE iterations performed.
 %
 % Notes:
 %  - Single-server stations plus an optional delay Z, matching pfqn_linearizer.
 %  - The exact counterpart is pfqn_sens_mom; the per-class exact second moments
 %    are in pfqn_sens_mva.
%}
%}
[M,R] = size(L);
N = ceil(N(:)');
if nargin < 3 || isempty(Z)
    Z = zeros(1,R);
end
Z = Z(:)';
if nargin < 4 || isempty(tol)
    tol = [];
end
if nargin < 5 || isempty(maxiter)
    maxiter = 200;
end
if length(N) ~= R
    line_error(mfilename,'demand matrix and population vector have different number of classes');
end
if any(isinf(N))
    line_error(mfilename,'pfqn_sens_linearizer requires a closed population');
end

if ~any(N > 0)
    res = pack(zeros(1,R),zeros(M,R),zeros(M,R),zeros(M,R),zeros(M,1), ...
               zeros(M,M),zeros(M,1),0);
    return;
end

% ---- Linearizer state ----------------------------------------------------
% pops: index 1 = N, index 1+l = N - e_l
npops = 1 + R;
pv = zeros(npops,R);
pv(1,:) = N;
for l = 1:R
    pv(1+l,:) = N;
    if N(l) > 0
        pv(1+l,l) = N(l) - 1;
    end
end

% m{p}(i,l) = estimate of m_i^{(pv(p,:))}(l), and its derivatives w.r.t. y_h
mE = cell(npops,1); dmE = cell(npops,1); d2mE = cell(npops,1);
for p = 1:npops
    mE{p} = zeros(M,R);
    for l = 1:R
        mE{p}(:,l) = pv(p,l)/M;      % initialization of the reference
    end
    dmE{p} = zeros(M,R,M);
    d2mE{p} = zeros(M,R,M);
end
% delta(i,l',l), indexed by the removed class l' and the class l
delta = zeros(M,R,R);
ddelta = zeros(M,R,R,M);
d2delta = zeros(M,R,R,M);

totiter = 0;
Xf = zeros(1,R); Qf = zeros(M,R); Wf = zeros(M,R);
dmf = zeros(M,M); d2mf = zeros(M,1);

for outer = 1:3
    % ---- Step 1: CORE at the full population --------------------------
    [mE{1}, dmE{1}, d2mE{1}, Xf, Wf, it] = core2(L,Z,N,delta,ddelta,d2delta, ...
                                                 mE{1}, dmE{1}, d2mE{1}, tol, maxiter);
    totiter = totiter + it;

    if outer < 3
        % ---- Step 2: CORE at each of the N - e_l populations ------------
        for l = 1:R
            if N(l) == 0
                continue;
            end
            [mE{1+l}, dmE{1+l}, d2mE{1+l}, ~, ~, it] = core2(L,Z,pv(1+l,:), ...
                delta,ddelta,d2delta, mE{1+l}, dmE{1+l}, d2mE{1+l}, tol, maxiter);
            totiter = totiter + it;
        end

        % ---- Step 3: refresh delta from (5.1) and (5.3) ----------------
        [vN, dvN, d2vN] = fractions(mE{1}, dmE{1}, d2mE{1}, N);
        for lp = 1:R
            if N(lp) == 0
                continue;
            end
            [vL, dvL, d2vL] = fractions(mE{1+lp}, dmE{1+lp}, d2mE{1+lp}, pv(1+lp,:));
            for i = 1:M
                for l = 1:R
                    delta(i,lp,l) = vL(i,l) - vN(i,l);
                    for h = 1:M
                        ddelta(i,lp,l,h) = dvL(i,l,h) - dvN(i,l,h);
                        d2delta(i,lp,l,h) = d2vL(i,l,h) - d2vN(i,l,h);
                    end
                end
            end
        end
    end
end

% ---- final measures ------------------------------------------------------
Qf = mE{1};
m = sum(Qf,2);
for i = 1:M
    for h = 1:M
        dmf(i,h) = sum(dmE{1}(i,:,h));
    end
    d2mf(i) = sum(d2mE{1}(i,:,i));
end
U = zeros(M,R);
for i = 1:M
    for r = 1:R
        U(i,r) = Xf(r) * L(i,r);
    end
end

res = pack(Xf,Qf,U,Wf,m,dmf,d2mf,totiter);
end

% =========================================================================
function [v, dv, d2v] = fractions(m, dm, d2m, n)
% (5.1) and (5.5): v_i(l) = m_i(l)/n(l), and likewise for the derivatives.
[M,R] = size(m);
v = zeros(M,R); dv = zeros(M,R,M); d2v = zeros(M,R,M);
for l = 1:R
    if n(l) <= 0
        continue;
    end
    v(:,l) = m(:,l) / n(l);
    dv(:,l,:) = dm(:,l,:) / n(l);
    d2v(:,l,:) = d2m(:,l,:) / n(l);
end
end

% =========================================================================
function [m, dm, d2m, lam, w, it] = core2(L,Z,n,delta,ddelta,d2delta,m,dm,d2m,tol,maxiter)
% CORE-2 of the reference, extended to second derivatives. Iterates (5.1),
% (5.2) and the MVA equations (1.4)-(1.5) together with their first and second
% derivatives (5.5)-(5.6) and (3.4), until the mean queue lengths and the
% variances both stop moving.
[M,R] = size(L);
nc = sum(n);
if isempty(tol)
    tolm = 1/(4000 + 16*nc);     % termination test of the reference
else
    tolm = tol;
end
tolv = 1e-3;

lam = zeros(1,R); w = zeros(M,R);
varprev = zeros(M,1);
it = 0;
for iter = 1:maxiter
    it = iter;
    mprev = m;

    % ---- (5.1)-(5.2) and (5.5)-(5.6): queue lengths one job down --------
    % maux(i,l,l2) estimates m_i^{(n - e_l)}(l2)
    [v, dv, d2v] = fractions(m, dm, d2m, n);
    mtot = zeros(M,R);                 % mtot(i,l) = sum_l2 m_i^{(n-e_l)}(l2)
    dmtot = zeros(M,R,M);
    d2mtot = zeros(M,R,M);
    for l = 1:R
        if n(l) <= 0
            continue;
        end
        for i = 1:M
            acc = 0; dacc = zeros(1,M); d2acc = zeros(1,M);
            for l2 = 1:R
                cnt = n(l2) - (l2 == l);       % (n - 1_l)_{l2}
                if cnt <= 0
                    continue;
                end
                acc = acc + cnt * (v(i,l2) + delta(i,l,l2));
                for h = 1:M
                    dacc(h) = dacc(h) + cnt * (dv(i,l2,h) + ddelta(i,l,l2,h));
                    d2acc(h) = d2acc(h) + cnt * (d2v(i,l2,h) + d2delta(i,l,l2,h));
                end
            end
            mtot(i,l) = acc;
            for h = 1:M
                dmtot(i,l,h) = dacc(h);
                d2mtot(i,l,h) = d2acc(h);
            end
        end
    end

    % ---- MVA (1.4)-(1.5) and its derivatives (3.4) ----------------------
    % w_i(l) = y_i * L(i,l) * (1 + mtot(i,l))
    w = zeros(M,R); dw = zeros(M,R,M); d2w = zeros(M,R,M);
    for l = 1:R
        if n(l) <= 0
            continue;
        end
        for i = 1:M
            A = 1 + mtot(i,l);
            w(i,l) = L(i,l) * A;
            for h = 1:M
                dA = dmtot(i,l,h);
                d2A = d2mtot(i,l,h);
                if i == h
                    dw(i,l,h) = L(i,l) * (A + dA);
                    d2w(i,l,h) = L(i,l) * (2*dA + d2A);
                else
                    dw(i,l,h) = L(i,l) * dA;
                    d2w(i,l,h) = L(i,l) * d2A;
                end
            end
        end
    end
    lam = zeros(1,R); dlam = zeros(R,M); d2lam = zeros(R,M);
    for l = 1:R
        if n(l) <= 0
            continue;
        end
        den = Z(l) + sum(w(:,l));
        lam(l) = n(l) / den;
        for h = 1:M
            dden = sum(dw(:,l,h));
            d2den = sum(d2w(:,l,h));
            dlam(l,h) = -n(l) * dden / den^2;
            d2lam(l,h) = -n(l) * d2den / den^2 + 2*n(l) * dden^2 / den^3;
        end
    end
    m = zeros(M,R); dm = zeros(M,R,M); d2m = zeros(M,R,M);
    for l = 1:R
        if n(l) <= 0
            continue;
        end
        for i = 1:M
            m(i,l) = lam(l) * w(i,l);
            for h = 1:M
                dm(i,l,h) = dlam(l,h)*w(i,l) + lam(l)*dw(i,l,h);
                d2m(i,l,h) = d2lam(l,h)*w(i,l) + 2*dlam(l,h)*dw(i,l,h) ...
                             + lam(l)*d2w(i,l,h);
            end
        end
    end

    % ---- termination test of the reference ------------------------------
    dev = 0;
    for l = 1:R
        if n(l) <= 0
            continue;
        end
        dev = max(dev, max(abs(m(:,l) - mprev(:,l))) / n(l));
    end
    varnow = zeros(M,1);
    for i = 1:M
        varnow(i) = sum(dm(i,:,i));
    end
    sv = sum(varnow);
    if sv > 0
        vdev = max(abs(varnow - varprev)) / sv;
    else
        vdev = 0;
    end
    varprev = varnow;
    if dev <= tolm && vdev <= tolv
        break;
    end
end
end

% =========================================================================
function res = pack(X,Q,U,W,m,dm,d2m,it)
M = size(Q,1);
res.X = X; res.Q = Q; res.U = U; res.W = W;
res.m = m; res.dm = dm; res.d2m = d2m; res.iter = it;

% Cov(i,j) = x_j dm_i/dx_j. The product form makes this symmetric, but the
% Linearizer fixed point does not enforce that, so the raw asymmetry is a
% genuine error indicator here rather than a roundoff residual.
Cov = dm;
res.CovAsym = max(max(abs(Cov - Cov.')));
if isempty(res.CovAsym)
    res.CovAsym = 0;
end
Cov = (Cov + Cov.') / 2;
res.Cov = Cov;

Var = zeros(M,1); M2 = zeros(M,1); M3 = zeros(M,1); Skew = zeros(M,1);
for i = 1:M
    Var(i) = dm(i,i);
    M2(i) = dm(i,i) + m(i)^2;
    M3(i) = d2m(i) + (1 + 3*m(i))*dm(i,i) + m(i)^3;
    mu3 = M3(i) - 3*m(i)*M2(i) + 2*m(i)^3;
    if Var(i) > 0
        Skew(i) = mu3 / Var(i)^1.5;
    else
        Skew(i) = NaN;
    end
end
res.Var = Var; res.M2 = M2; res.M3 = M3; res.Skew = Skew;
end
