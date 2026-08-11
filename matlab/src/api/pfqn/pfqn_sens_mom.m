%{
%{
 % @file pfqn_sens_mom.m
 % @brief Exact higher moments (up to order three) of the per-station total
 %        queue lengths of a closed product-form queueing network, by
 %        second-order differentiation of the MVA recursion.
%}
%}

function mom = pfqn_sens_mom(L,N,Z,mi,groups)
%{
%{
 % @brief Exact moments E[Q_i], E[Q_i^2], E[Q_i^3] and the covariances
 %        Cov[Q_i,Q_j] of the TOTAL queue lengths Q_i = sum_r n(i,r) of a closed
 %        product-form (BCMP) queueing network.
 %
 %        The method is the moment analysis of Strelen. Its Theorem 3.1 states
 %        that one further factor Q_i in a moment costs one differentiation with
 %        respect to x_i, the reciprocal of the capacity of station i, i.e. a
 %        parameter that scales the service times of ALL classes at station i:
 %
 %          E[Q_i^j] = m_i E[Q_i^(j-1)] + x_i d/dx_i E[Q_i^(j-1)],  Q_i^0 = 1
 %
 %        Iterating from E[Q_i^0] = 1 gives, with m_i = E[Q_i] (equation (3.2)):
 %
 %          Var[Q_i]    = x_i dm_i/dx_i
 %          Cov[Q_i,Q_j]= x_j dm_i/dx_j = x_i dm_j/dx_i
 %          E[Q_i^2]    = x_i dm_i/dx_i + m_i^2
 %          E[Q_i^3]    = x_i^2 d^2m_i/dx_i^2 + (x_i + 3 x_i m_i) dm_i/dx_i + m_i^3
 %
 %        so the third moment requires the SECOND derivative of the MVA
 %        recursion, which is what this routine adds over pfqn_sens_mva and
 %        pfqn_sens (both first order only). The derivatives are obtained by
 %        second-order forward-mode differentiation of the Reiser-Lavenberg
 %        recursion, i.e. by carrying, for each parameter, the value together
 %        with its first and second derivative along the population lattice
 %        (Theorem 3.2 for one class, Theorem 3.5 for several).
 %
 %        The parameter need not scale a whole column. Theorem 1 of Akyildiz and
 %        Strelen states the same recursion for a parameter that scales the
 %        service times of an arbitrary class subset T at station i, and the
 %        moments it generates are then those of Q_(i,T) = sum_(r in T) n(i,r).
 %        GROUPS supplies that subset structure: it partitions the classes, and
 %        the routine reports the moments of each group's queue length at each
 %        station. The three useful settings are
 %          groups = ones(1,R)  the whole column: per-station TOTALS (default,
 %                              this is Strelen's x_i);
 %          groups = 1:R        one class per group: PER-CLASS moments, so that
 %                              even the third moment is per class;
 %          groups = chain(r)   one group per chain: PER-CHAIN moments.
 %        Strelen states only the first; the generalization is Akyildiz and
 %        Strelen's T. pfqn_sens_mom_validate checks the per-class setting
 %        against brute force and against pfqn_sens_mva's second moments, which
 %        it must reproduce exactly.
 %
 %        Reference: J. C. Strelen, "Moment Analysis for Closed Queuing Networks
 %        and its Linearizer", Performance Evaluation 11:127-142, 1990,
 %        Theorems 2.1, 3.1, 3.2, 3.5 and equation (3.2).
 %
 % @fn pfqn_sens_mom(L, N, Z, mi)
 % @param L  Service demand matrix (M x R), L(i,r) = visits_ir / rate_ir.
 % @param N  Population vector (1 x R).
 % @param Z  Think time vector (1 x R). Default: zeros.
 % @param mi (Optional) Server multiplicity vector (1 x M). Default: ones.
 % @param groups (Optional) Class-to-group map (1 x R), a partition of the
 %        classes into G = max(groups) groups labelled 1..G. The moments
 %        returned are those of each group's queue length at each station.
 %        Default: ones(1,R), i.e. one group holding every class, which is the
 %        per-station total.
 % @return mom A struct with:
 %   .X (1 x R), .Q (M x R), .U (M x R), .R (M x R)  base MVA measures,
 %       identical to pfqn_mva(L,N,Z,mi).
 %   .m (M x G)      m(i,g) = E[Q_(i,g)], the mean queue length of group g at
 %       station i. With the default groups this is (M x 1), the station total.
 %   .d2m (M x G)    d2m(i,g) = the scaled pure second derivative with respect
 %       to the parameter of (i,g). Only that entry is needed by (3.2); the
 %       mixed second derivatives are not required for moments of a single
 %       Q_(i,g) and would cost an extra factor M*G to carry.
 %   .Var (M x G)    Var[Q_(i,g)].
 %   .M2 (M x G)     E[Q_(i,g)^2].
 %   .M3 (M x G)     E[Q_(i,g)^3].
 %   .Skew (M x G)   skewness of Q_(i,g). NaN where Var is zero.
 %   .Cov, .dm       Cov((i,g),(j,g')) = Cov[Q_(i,g),Q_(j,g')] and the scaled
 %       first derivative it comes from. Shaped (M x G x M x G) in general, but
 %       COLLAPSED to (M x M) in the default single-group case, where the group
 %       index carries no information and an (M x 1 x M x 1) array would only be
 %       awkward to index.
 %   .CovAsym (scalar)  raw asymmetry of Cov before symmetrization. The two
 %       triangles are distinct expressions that must agree, so this is a live
 %       residual of the recursion; expect roundoff.
 %
 % Notes:
 %  - Restricted to closed populations, as is the moment analysis of the
 %    reference. Mixed and load-dependent second moments are in
 %    pfqn_sens_mvaldmx.
 %  - Moments of the sojourn times at FCFS centers are built on top of these
 %    queue-length moments by pfqn_sens_respt, following Theorem 4.1.
 %  - The exact recursion costs O(prod(N+1)) lattice points; pfqn_sens_linearizer
 %    approximates the same quantities in polynomial time.
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
if nargin < 5 || isempty(groups)
    groups = ones(1,R);
end
groups = round(groups(:)');
if length(N) ~= R
    line_error(mfilename,'demand matrix and population vector have different number of classes');
end
if any(isinf(N))
    line_error(mfilename,'pfqn_sens_mom requires a closed population');
end
if length(groups) ~= R || any(groups < 1)
    line_error(mfilename,'groups must be a (1 x R) vector of group labels starting at 1');
end
G = max(groups);
if ~isequal(unique(groups), 1:G)
    line_error(mfilename,'groups must label the classes consecutively from 1 to max(groups), with no empty group');
end

X = zeros(1,R); Q = zeros(M,R); U = zeros(M,R); C = zeros(M,R);
m = zeros(M,G); dm = zeros(M,G,M,G); d2m = zeros(M,G);

if ~any(N > 0)
    mom = pack(X,Q,U,C,m,dm,d2m,G);
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
% Parameters are indexed by (station h, group g): y_(h,g) scales L(h,r) for
% every class r in group g. At y = 1 the y-derivatives are exactly the scaled
% x-derivatives that (3.2) asks for, since a pure rescaling x -> x*y gives
% d/dy = x d/dx and d2/dy2 = x^2 d2/dx2.
%   Qtot(row,i)   = sum_r Q(i,r) at population row; the recursion needs only the
%                   station total, whatever the grouping.
%   D1(row,i,p)   = d/dy_p Qtot(i),  D2(row,i,p) = d^2/dy_p^2 Qtot(i)
% The per-group queue lengths and their derivatives are accumulated for the
% CURRENT population only and overwritten each step, so at the end of the walk
% they hold the values at N, exactly as X, Q and C do.
P = M*G;
pidx = zeros(M,G);
for h = 1:M
    for g = 1:G
        pidx(h,g) = (h-1)*G + g;
    end
end
Qtot = zeros(totpop,M);
D1 = zeros(totpop,M,P);
D2 = zeros(totpop,M,P);
Qg = zeros(M,G); D1Qg = zeros(M,G,P); D2Qg = zeros(M,G,P);
currentpop = 2;
n = zeros(1,R);
n(firstnonempty) = 1;

Cs = zeros(M,1); dCs = zeros(M,P); d2Cs = zeros(M,P);

while ctr
    % the group accumulators describe one population only
    Qg(:) = 0; D1Qg(:) = 0; D2Qg(:) = 0;
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

        % ---- residence times and their first two derivatives -----------
        % C(i,s) = L(i,s)*y_(i,g(s))*(mi(i)+Qtot(i|n-e_s)) =: L(i,s)*y_(i,g(s))*A
        % The parameter (h,g) touches C(i,s) directly only when i == h AND class
        % s belongs to group g; otherwise it acts only through A.
        %   d/dy_p   C = L(i,s)*( [i==h & g(s)==g]*A + dA/dy_p )
        %   d2/dy_p2 C = L(i,s)*( 2*[i==h & g(s)==g]*dA/dy_p + d2A/dy_p2 )
        gs = groups(s);
        CNtot = 0;
        dCNtot = zeros(1,P);
        d2CNtot = zeros(1,P);
        for i = 1:M
            A = mi(i) + Qtot(row,i);
            Cs(i) = L(i,s) * A;
            C(i,s) = Cs(i);
            CNtot = CNtot + Cs(i);
            for p = 1:P
                dA = D1(row,i,p);
                d2A = D2(row,i,p);
                if p == pidx(i,gs)
                    dCs(i,p) = L(i,s) * (A + dA);
                    d2Cs(i,p) = L(i,s) * (2*dA + d2A);
                else
                    dCs(i,p) = L(i,s) * dA;
                    d2Cs(i,p) = L(i,s) * d2A;
                end
                dCNtot(p) = dCNtot(p) + dCs(i,p);
                d2CNtot(p) = d2CNtot(p) + d2Cs(i,p);
            end
        end

        % ---- throughput and its first two derivatives -------------------
        % X(s) = n(s)/den,  den = Z(s) + sum_i C(i,s)
        %   dX   = -n(s)*dden/den^2
        %   d2X  = -n(s)*d2den/den^2 + 2*n(s)*dden^2/den^3
        den = Z(s) + CNtot;
        X(s) = n(s) / den;
        dX = zeros(1,P); d2X = zeros(1,P);
        for p = 1:P
            dX(p) = -n(s) * dCNtot(p) / den^2;
            d2X(p) = -n(s) * d2CNtot(p) / den^2 + 2*n(s) * dCNtot(p)^2 / den^3;
        end

        % ---- queue lengths ---------------------------------------------
        % Q = X*C,  dQ = dX*C + X*dC,  d2Q = d2X*C + 2*dX*dC + X*d2C
        for i = 1:M
            Q(i,s) = X(s) * Cs(i);
            Qtot(currentpop,i) = Qtot(currentpop,i) + Q(i,s);
            Qg(i,gs) = Qg(i,gs) + Q(i,s);
            for p = 1:P
                dQ = dX(p)*Cs(i) + X(s)*dCs(i,p);
                d2Q = d2X(p)*Cs(i) + 2*dX(p)*dCs(i,p) + X(s)*d2Cs(i,p);
                D1(currentpop,i,p) = D1(currentpop,i,p) + dQ;
                D2(currentpop,i,p) = D2(currentpop,i,p) + d2Q;
                D1Qg(i,gs,p) = D1Qg(i,gs,p) + dQ;
                D2Qg(i,gs,p) = D2Qg(i,gs,p) + d2Q;
            end
        end
        s = s + 1;
    end

    % ---- odometer advance ---------------------------------------------
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

% moments at the full population: Qg, D1Qg and D2Qg were overwritten on every
% population sweep, so they now hold the values at N.
for i = 1:M
    for g = 1:G
        m(i,g) = Qg(i,g);
        for j = 1:M
            for g2 = 1:G
                dm(i,g,j,g2) = D1Qg(i,g,pidx(j,g2));
            end
        end
        d2m(i,g) = D2Qg(i,g,pidx(i,g));
    end
end

mom = pack(X,Q,U,C,m,dm,d2m,G);
end

% =========================================================================
function mom = pack(X,Q,U,C,m,dm,d2m,G)
M = size(Q,1);
mom.X = X; mom.Q = Q; mom.U = U; mom.R = C;
mom.m = m; mom.d2m = d2m;

% Cov((i,g),(j,g')) = the scaled derivative of m_(i,g) w.r.t. the parameter of
% (j,g'), and the transposed entry, are distinct expressions for the same
% quantity; report the raw disagreement, then symmetrize.
flat = reshape(dm, M*G, M*G);
mom.CovAsym = max(max(abs(flat - flat.')));
if isempty(mom.CovAsym)
    mom.CovAsym = 0;
end
flat = (flat + flat.') / 2;
if G == 1
    % the group index carries no information here; an (M x 1 x M x 1) array
    % would only be awkward to index
    mom.Cov = reshape(flat, M, M);
    mom.dm = reshape(dm, M, M);
else
    mom.Cov = reshape(flat, M, G, M, G);
    mom.dm = dm;
end

Var = zeros(M,G); M2 = zeros(M,G); M3 = zeros(M,G); Skew = zeros(M,G);
for i = 1:M
    for g = 1:G
        d1 = dm(i,g,i,g);
        Var(i,g) = d1;                                            % (3.2)
        M2(i,g) = d1 + m(i,g)^2;                                  % (3.2)
        M3(i,g) = d2m(i,g) + (1 + 3*m(i,g))*d1 + m(i,g)^3;        % (3.2)
        % third central moment mu3 = E[Q^3] - 3 m E[Q^2] + 2 m^3
        mu3 = M3(i,g) - 3*m(i,g)*M2(i,g) + 2*m(i,g)^3;
        if Var(i,g) > 0
            Skew(i,g) = mu3 / Var(i,g)^1.5;
        else
            Skew(i,g) = NaN;
        end
    end
end
mom.Var = Var; mom.M2 = M2; mom.M3 = M3; mom.Skew = Skew;
end
