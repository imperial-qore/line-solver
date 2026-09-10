%{
%{
 % @file pfqn_sens_respt.m
 % @brief Exact moments of the sojourn time of a job at FCFS multiserver
 %        centers of a closed product-form queueing network.
%}
%}

function res = pfqn_sens_respt(S,V,N,Z,b,tmax)
%{
%{
 % @brief Exact raw moments E[W_(i,l)^t], t = 1..tmax, of the sojourn time of a
 %        class-l job at an FCFS b-server center i of a closed product-form
 %        queueing network, together with the variance of that sojourn time.
 %
 %        This is Theorem 4.1 of the reference. Its mechanism is the arrival
 %        theorem of Lavenberg-Reiser and Sevcik-Mitrani: a class-l job arriving
 %        at center i finds j jobs already there with probability
 %        p_i(j, N - 1_l). Conditioning the sojourn time on j and inverting the
 %        Laplace transform of the conditional density gives
 %
 %          E[W_(i,l)^t] = t!/mu^t + sum_{tau=0..t} a_(t,tau)(0) E[Qt_i^tau]
 %                         - sum_{j=0..b-1} p_i(j,N-1_l) sum_{tau=0..t} a_(t,tau)(0) j^tau
 %
 %        where mu = 1/S(i) is the rate of each of the b servers, Qt_i is the
 %        total queue length at center i at population N - 1_l (so its moments
 %        are those of pfqn_sens_mom evaluated one job down in class l), and the
 %        coefficients a_(t,tau)(0) depend only on b and mu, not on the network
 %        (Remark 4.3 of the reference). The moments E[Qt_i^tau] up to tau = 3
 %        need the second derivative of the MVA recursion, so this routine
 %        carries a second-order forward-mode pass exactly as pfqn_sens_mom
 %        does, but over the b-server recursion (4.1)-(4.2) rather than the
 %        single-server one.
 %
 %        For b = 1 the coefficients a_(t,0)(0) vanish identically and the
 %        double-sum correction disappears, so no marginal probabilities are
 %        needed (Remark 4.2 of the reference); the routine still evaluates the
 %        general expression, which reduces to that case on its own.
 %
 %        Only FCFS centers are covered. The reference is explicit that the
 %        sojourn-time distribution at PS and LCFS centers is in general not
 %        known, so no analogue exists there. FCFS in a BCMP network further
 %        requires the service time to be exponential and class-independent,
 %        which is why this routine takes a per-station service time S(i) and a
 %        separate visit-ratio matrix V rather than a demand matrix: the sojourn
 %        time is per visit, so the per-visit rate mu = 1/S(i) must be known and
 %        cannot be recovered from the demand L(i,l) = S(i)*V(i,l) alone.
 %
 %        Reference: J. C. Strelen, "Moment Analysis for Closed Queuing Networks
 %        and its Linearizer", Performance Evaluation 11:127-142, 1990,
 %        Theorem 4.1 with equations (4.1)-(4.5) and Remarks 4.2-4.3.
 %
 % @fn pfqn_sens_respt(S, V, N, Z, b, tmax)
 % @param S  Service time at each station (M x 1), common to all classes.
 % @param V  Visit ratio matrix (M x R). The demand is L(i,r) = S(i)*V(i,r).
 % @param N  Population vector (1 x R).
 % @param Z  Think time vector (1 x R). Default: zeros.
 % @param b  Number of servers at each station (M x 1). Default: ones.
 % @param tmax Highest sojourn-time moment to return, 1..3. Default: 3. The
 %        coefficients a_(t,tau)(0) are tabulated in the reference up to t = 3.
 % @return res A struct with:
 %   .X (1 x R), .Q (M x R), .U (M x R)  base measures at population N.
 %   .W (M x R)        W(i,l) = E[W_(i,l)], the mean sojourn time per visit of
 %       a class-l job at station i. Zero where class l does not visit i.
 %   .WM (M x R x tmax) WM(i,l,t) = E[W_(i,l)^t].
 %   .WVar (M x R)     Var[W_(i,l)] = E[W^2] - E[W]^2. Requires tmax >= 2.
 %   .WSkew (M x R)    skewness of W_(i,l). Requires tmax >= 3; NaN if the
 %       variance is zero.
 %   .m (M x 1)        E[Q_i] at population N, the total queue length.
 %   .Var (M x 1)      Var[Q_i] at population N.
 %   .p (M x max(b))   p(i,1+j) = P[Q_i = j] at population N, for j = 0..b_i-1.
 %       These are the only marginal probabilities the b-server recursion needs,
 %       so the matrix is RAGGED: row i is meaningful only up to column b_i and
 %       is zero-padded out to max(b). A padded entry is not P[Q_i = j]; it is
 %       simply not computed. Read row i as p(i,1:b(i)).
 %   .Wresid (M x R)   the residence time w_i(l) of the MVA recursion. The
 %       identity W(i,l) = w_i(l)/V(i,l) is an independent check of the t = 1
 %       case of (4.5) and is asserted by pfqn_sens_respt_validate.
 %
 % Notes:
 %  - Restricted to closed populations.
 %  - Load-dependent rates are not covered here; the b-server dependence is the
 %    only state dependence, and it is carried exactly by (4.1)-(4.2).
%}
%}
[M,R] = size(V);
S = S(:);
N = ceil(N(:)');
if nargin < 4 || isempty(Z)
    Z = zeros(1,R);
end
Z = Z(:)';
if nargin < 5 || isempty(b)
    b = ones(M,1);
end
b = round(b(:));
if nargin < 6 || isempty(tmax)
    tmax = 3;
end
if tmax < 1 || tmax > 3
    line_error(mfilename,'tmax must be 1, 2 or 3: the coefficients a_{t,tau}(0) are tabulated in the reference up to order three.');
end
if length(N) ~= R
    line_error(mfilename,'visit matrix and population vector have different number of classes');
end
if any(isinf(N))
    line_error(mfilename,'pfqn_sens_respt requires a closed population');
end
if any(b < 1)
    line_error(mfilename,'the number of servers must be at least one at every station');
end
if any(S <= 0)
    line_error(mfilename,'every FCFS station must have a strictly positive service time');
end

rho = zeros(M,R);
for i = 1:M
    for l = 1:R
        rho(i,l) = S(i) * V(i,l);
    end
end
mu = 1 ./ S;
bmax = max(b);

X = zeros(1,R); Q = zeros(M,R); U = zeros(M,R);
m = zeros(M,1);
W = zeros(M,R); WM = zeros(M,R,tmax); Wresid = zeros(M,R);
pN = zeros(M,bmax);

if ~any(N > 0)
    res = pack(X,Q,U,m,zeros(M,1),pN,W,WM,Wresid,tmax);
    return;
end

% ---- population lattice --------------------------------------------------
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

% state carried along the lattice: the mean total queue length at each station,
% the marginal probabilities p_i(j) for j = 0..bmax-1, and the first and second
% derivatives of both with respect to y_h, a scaling of station h's service
% time. At y = 1 the y-derivatives are the scaled x-derivatives of (3.2).
Mrow  = zeros(totpop,M);
D1m   = zeros(totpop,M,M);
D2m   = zeros(totpop,M,M);
Prow  = zeros(totpop,M,bmax);
D1p   = zeros(totpop,M,bmax,M);
D2p   = zeros(totpop,M,bmax,M);
Prow(1,:,1) = 1;   % empty population: every station holds zero jobs

currentpop = 2;
n = zeros(1,R);
n(firstnonempty) = 1;
rows = ones(1,R);

while ctr
    hnvec = currentpop;
    % ---- residence times, eq. (4.1), and their first two derivatives ----
    wv = zeros(M,R); d1w = zeros(M,R,M); d2w = zeros(M,R,M);
    for s = 1:R
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
        rows(s) = row;
        if n(s) == 0
            continue;   % w and every derivative stay zero, as does X(s)
        end
        for i = 1:M
            % bracket = 1 + m_i(n-e_s) + sum_{j=0}^{b_i-2} (b_i-1-j) p_i(j,n-e_s)
            brk = 1 + Mrow(row,i);
            for j = 0:(b(i)-2)
                brk = brk + (b(i)-1-j) * Prow(row,i,1+j);
            end
            wv(i,s) = (rho(i,s)/b(i)) * brk;
            for h = 1:M
                dbrk = D1m(row,i,h);
                d2brk = D2m(row,i,h);
                for j = 0:(b(i)-2)
                    dbrk = dbrk + (b(i)-1-j) * D1p(row,i,1+j,h);
                    d2brk = d2brk + (b(i)-1-j) * D2p(row,i,1+j,h);
                end
                % w = y_i * (rho/b) * brk
                if i == h
                    d1w(i,s,h) = (rho(i,s)/b(i)) * (brk + dbrk);
                    d2w(i,s,h) = (rho(i,s)/b(i)) * (2*dbrk + d2brk);
                else
                    d1w(i,s,h) = (rho(i,s)/b(i)) * dbrk;
                    d2w(i,s,h) = (rho(i,s)/b(i)) * d2brk;
                end
            end
        end
    end

    % ---- throughputs and their derivatives ------------------------------
    lam = zeros(1,R); d1lam = zeros(R,M); d2lam = zeros(R,M);
    for s = 1:R
        if n(s) == 0
            continue;
        end
        den = Z(s) + sum(wv(:,s));
        lam(s) = n(s) / den;
        for h = 1:M
            dden = sum(d1w(:,s,h));
            d2den = sum(d2w(:,s,h));
            d1lam(s,h) = -n(s) * dden / den^2;
            d2lam(s,h) = -n(s) * d2den / den^2 + 2*n(s) * dden^2 / den^3;
        end
    end

    % ---- mean queue lengths ---------------------------------------------
    for i = 1:M
        acc = 0;
        for s = 1:R
            if n(s) == 0, continue; end
            acc = acc + lam(s) * wv(i,s);
        end
        Mrow(hnvec,i) = acc;
        for h = 1:M
            d1acc = 0; d2acc = 0;
            for s = 1:R
                if n(s) == 0, continue; end
                d1acc = d1acc + d1lam(s,h)*wv(i,s) + lam(s)*d1w(i,s,h);
                d2acc = d2acc + d2lam(s,h)*wv(i,s) + 2*d1lam(s,h)*d1w(i,s,h) ...
                        + lam(s)*d2w(i,s,h);
            end
            D1m(hnvec,i,h) = d1acc;
            D2m(hnvec,i,h) = d2acc;
        end
    end

    % ---- marginal probabilities, eq. (4.2), and their derivatives --------
    nc = sum(n);
    for i = 1:M
        % p_i(j,n) = (1/j) sum_l lam(l) * rho_i(l)*y_i * p_i(j-1, n-e_l)
        for j = 1:(b(i)-1)
            if j > nc
                Prow(hnvec,i,1+j) = 0;
                continue;
            end
            acc = 0;
            for l = 1:R
                if n(l) == 0, continue; end
                acc = acc + lam(l) * rho(i,l) * Prow(rows(l),i,1+(j-1));
            end
            Prow(hnvec,i,1+j) = acc / j;
            for h = 1:M
                d1acc = 0; d2acc = 0;
                for l = 1:R
                    if n(l) == 0, continue; end
                    pprev = Prow(rows(l),i,1+(j-1));
                    d1prev = D1p(rows(l),i,1+(j-1),h);
                    d2prev = D2p(rows(l),i,1+(j-1),h);
                    % g = rho * u * v with u = lam, v = y_i * pprev
                    if i == h
                        v1 = pprev + d1prev;
                        v2 = 2*d1prev + d2prev;
                    else
                        v1 = d1prev;
                        v2 = d2prev;
                    end
                    d1acc = d1acc + rho(i,l) * (d1lam(l,h)*pprev + lam(l)*v1);
                    d2acc = d2acc + rho(i,l) * (d2lam(l,h)*pprev + 2*d1lam(l,h)*v1 ...
                                                + lam(l)*v2);
                end
                D1p(hnvec,i,1+j,h) = d1acc / j;
                D2p(hnvec,i,1+j,h) = d2acc / j;
            end
        end
        % u_i = sum_l lam(l) * rho_i(l)*y_i  (mean number of busy servers)
        ui = 0; d1ui = zeros(1,M); d2ui = zeros(1,M);
        for l = 1:R
            if n(l) == 0, continue; end
            ui = ui + lam(l) * rho(i,l);
            for h = 1:M
                if i == h
                    d1ui(h) = d1ui(h) + rho(i,l) * (d1lam(l,h) + lam(l));
                    d2ui(h) = d2ui(h) + rho(i,l) * (d2lam(l,h) + 2*d1lam(l,h));
                else
                    d1ui(h) = d1ui(h) + rho(i,l) * d1lam(l,h);
                    d2ui(h) = d2ui(h) + rho(i,l) * d2lam(l,h);
                end
            end
        end
        % p_i(0,n) = 1 - (1/b)(u_i + sum_{j=1}^{b-1} (b-j) p_i(j,n))
        acc0 = ui;
        for j = 1:(b(i)-1)
            acc0 = acc0 + (b(i)-j) * Prow(hnvec,i,1+j);
        end
        Prow(hnvec,i,1) = 1 - acc0/b(i);
        for h = 1:M
            d1acc0 = d1ui(h); d2acc0 = d2ui(h);
            for j = 1:(b(i)-1)
                d1acc0 = d1acc0 + (b(i)-j) * D1p(hnvec,i,1+j,h);
                d2acc0 = d2acc0 + (b(i)-j) * D2p(hnvec,i,1+j,h);
            end
            D1p(hnvec,i,1,h) = -d1acc0/b(i);
            D2p(hnvec,i,1,h) = -d2acc0/b(i);
        end
    end

    % keep the measures of the last (full) population
    X = lam;
    for i = 1:M
        for s = 1:R
            Wresid(i,s) = wv(i,s);
            Q(i,s) = lam(s) * wv(i,s);
            U(i,s) = lam(s) * rho(i,s);
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

lastrow = currentpop;
for i = 1:M
    m(i) = Mrow(lastrow,i);
    for j = 0:(bmax-1)
        pN(i,1+j) = Prow(lastrow,i,1+j);
    end
end
Var = zeros(M,1);
for i = 1:M
    Var(i) = D1m(lastrow,i,i);
end

% ---- sojourn-time moments, eq. (4.5) -------------------------------------
% index of N - e_l on the lattice
rowsN = ones(1,R);
for l = 1:R
    if N(l) > 0
        nn = N; nn(l) = nn(l) - 1;
        pos = nn(R);
        for w = 1:R-1
            pos = pos + nn(w)*prods(w);
        end
        rowsN(l) = 1 + pos;
    end
end

for i = 1:M
    for l = 1:R
        if N(l) == 0 || V(i,l) <= 0
            continue;
        end
        rl = rowsN(l);
        % moments of the queue length seen by an arriving class-l job, i.e. of
        % the total queue at station i at population N - e_l, from (3.2)
        mt = Mrow(rl,i);
        d1t = D1m(rl,i,i);
        d2t = D2m(rl,i,i);
        EQ = zeros(1,4);         % EQ(1+tau) = E[Qt_i^tau]
        EQ(1) = 1;
        EQ(2) = mt;
        EQ(3) = d1t + mt^2;
        EQ(4) = d2t + (1 + 3*mt)*d1t + mt^3;
        acoef = strelen_a(b(i), mu(i), tmax);
        for t = 1:tmax
            val = factorial(t) / mu(i)^t;
            for tau = 0:t
                val = val + acoef(t,1+tau) * EQ(1+tau);
            end
            % correction over the states in which a server is idle
            for j = 0:(b(i)-1)
                inner = 0;
                for tau = 0:t
                    inner = inner + acoef(t,1+tau) * jpow(j,tau);
                end
                val = val - Prow(rl,i,1+j) * inner;
            end
            WM(i,l,t) = val;
        end
        W(i,l) = WM(i,l,1);
    end
end

res = pack(X,Q,U,m,Var,pN,W,WM,Wresid,tmax);
end

% =========================================================================
function v = jpow(j,tau)
% j^tau with the convention j^0 = 1, so that 0^0 = 1 as the reference states.
if tau == 0
    v = 1;
else
    v = j^tau;
end
end

% =========================================================================
function a = strelen_a(b,mu,tmax)
% Coefficients a_{t,tau}(0) of Remark 4.3 of the reference. They depend only on
% the number of servers b and on the per-server rate mu, not on the network.
a = zeros(3,4);
a(1,1) = (1-b)/(b*mu);
a(1,2) = 1/(b*mu);
if tmax >= 2
    a(2,1) = (2 - b - b^2)/(b^2*mu^2);
    a(2,2) = 3/(b^2*mu^2);
    a(2,3) = 1/(b^2*mu^2);
end
if tmax >= 3
    a(3,1) = (6 - 5*b + 3*b^2 - 4*b^3)/(b^3*mu^3);
    a(3,2) = (11 - 3*b + 3*b^2)/(b^3*mu^3);
    a(3,3) = 6/(b^3*mu^3);
    a(3,4) = 1/(b^3*mu^3);
end
end

% =========================================================================
function res = pack(X,Q,U,m,Var,p,W,WM,Wresid,tmax)
[M,R] = size(Q);
res.X = X; res.Q = Q; res.U = U;
res.m = m; res.Var = Var; res.p = p;
res.W = W; res.WM = WM; res.Wresid = Wresid;
WVar = zeros(M,R); WSkew = zeros(M,R);
if tmax >= 2
    for i = 1:M
        for l = 1:R
            WVar(i,l) = WM(i,l,2) - WM(i,l,1)^2;
        end
    end
end
if tmax >= 3
    for i = 1:M
        for l = 1:R
            mu3 = WM(i,l,3) - 3*WM(i,l,1)*WM(i,l,2) + 2*WM(i,l,1)^3;
            if WVar(i,l) > 0
                WSkew(i,l) = mu3 / WVar(i,l)^1.5;
            else
                WSkew(i,l) = NaN;
            end
        end
    end
end
res.WVar = WVar; res.WSkew = WSkew;
end
