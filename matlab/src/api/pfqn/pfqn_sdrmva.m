%{
%{
 % @file pfqn_sdrmva.m
 % @brief MVA and convolution solution of a network with state-dependent routing.
%}
%}

%{
%{
 % @brief MVA and convolution solution of a network with state-dependent routing.
 % @fn pfqn_sdrmva(S, xi, N, sdr, alpha)
 % @param S Mean service times per station and chain.
 % @param xi Relative visit counts per station and chain.
 % @param N Chain population vector.
 % @param sdr State-dependent routing structure.
 % @param alpha Optional load-dependent rate scalings.
 % @return Q Mean queue lengths per station and chain.
 % @return X Mean throughputs per station and chain.
 % @return U Mean utilizations per station and chain.
 % @return R Mean response times per station and chain.
 % @return lG Logarithm of the normalizing constant.
%}
%}
function [Q, X, U, R, lG] = pfqn_sdrmva(S, xi, N, sdr, alpha)
% [Q, X, U, R, LG] = PFQN_SDRMVA(S, XI, N, SDR, ALPHA)
%
% Mean value analysis and convolution of a closed multiclass network with the
% state-dependent routing of Krzesinski (1987), "Multiclass Queueing Networks
% with State-Dependent Routing", Performance Evaluation 7:125-143, Section 4.
%
% The solution proceeds level by level, outermost subnetwork last:
%   Sec. 4.2.1  a modified MVA of the centers of V_t - V_{t+1} in isolation,
%               whose arrival theorem carries the SDR admission coefficient
%               delta_ti(n) = d_ti - n and the center's own queue-length
%               distribution;
%   Sec. 4.2.2  convolution of that level with the already solved inner
%               subnetwork Q(V,V_{t+1}), weighted by Omega_{t-1,t}/Omega_tt;
%   Sec. 4.2.3  the same convolution re-normalizes the inner queue lengths;
%   Sec. 4.3    a final convolution against the complement M-V.
%
% Same signature and same outputs as PFQN_SDR, which evaluates the product form
% (16) exactly by state enumeration, so the two are directly comparable. This
% routine costs O(J T M (V_1...V_J)^2) rather than the size of the state space,
% and is the one to use when the populations make enumeration impractical.
%
% RESTRICTIONS, both from the paper itself:
%   - every SDR branch must be a SINGLE center. The paper's Section 4 is stated
%     that way and defers the general case to an unpublished technical report;
%     PFQN_SDR has no such restriction.
%   - every C_t must be negative. Section 2.5 assumes it, and Section 4 is
%     written throughout for C_t = -1. A structure with C_t < 0 but not -1 is
%     rescaled internally, which leaves both the routing probabilities of
%     eq. (10) and the weights of eq. (16) unchanged: replacing (C_t, d_tb) by
%     (C_t/k_t, d_tb/k_t) scales delta_tb and omega_tt by 1/k_t and
%     omega_{t-1,t} by 1/k_{t-1}, and the factors telescope away.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = size(S,1);
J = size(S,2);
N = round(N(:)');
if nargin < 5 || isempty(alpha)
    alpha = ones(M, max(1,sum(N)));
end
if size(alpha,2) < sum(N)
    alpha = [alpha, ones(M, sum(N)-size(alpha,2))];
end

c0 = pfqn_sdrcoeff(sdr);
for b = 2:c0.B
    if numel(c0.branch{b}) ~= 1
        line_error(mfilename,['pfqn_sdrmva requires every SDR branch to hold a single center: ', ...
            'the MVA and convolution of Krzesinski (1987) Section 4 is stated that way and its general ', ...
            'case is in an unpublished technical report. Use pfqn_sdr, which evaluates eq. (16) exactly ', ...
            'for any branch topology.']);
    end
end
if any(c0.C >= 0)
    line_error(mfilename,['pfqn_sdrmva requires every C_t to be negative: Section 2.5 assumes it and ', ...
        'Section 4 is written for C_t = -1. A nonnegative C_t leaves the branch populations unbounded ', ...
        'and the queue-length recursion of Section 4.2.1 without a top.']);
end

% Rescale each level to C_t = -1, which leaves eqs. (10) and (16) unchanged
sdr1 = sdr;
for t = 1:c0.T
    k = -c0.C(t);
    sdr1.C(t) = -1;
    sdr1.d(t,:) = sdr.d(t,:) / k;
end
c = pfqn_sdrcoeff(sdr1);

gamma = xi .* S;
Ntot = sum(N);
[latt, lidx] = sub_lattice(N);
nl = size(latt,1);

inV = false(1,M);
dvec = nan(1,M);          % delta_i(n) = dvec(i) - n; NaN marks a center outside V
lvl = zeros(1,M);
for b = 2:c.B
    i = c.branch{b};
    inV(i) = true;
    lvl(i) = c.level(b);
    dvec(i) = c.d(c.level(b), b);
end
mv = find(~inV);

% Sec. 4.1: the complement M-V is an ordinary BCMP subnetwork, so its MVA is
% the same recursion with delta = 1
[Qc, Tc, gc] = sub_set_mva(mv, nan(1,numel(mv)), gamma, alpha, N, latt, lidx, Ntot);

% Sec. 4.2: outermost level last
Gin = zeros(1,nl); Gin(lidx(sub_key(zeros(1,J), N))) = 1;  % G(L, V_{T+1}) = [L == 0]
Qin = zeros(M, J, nl);
Tin = zeros(J, nl);
Ain = zeros(M, nl);       % E[P_{e,e(i)} | population of the subnetwork], see below
for t = c.T:-1:1
    St = find(lvl == t);
    [Qt, Tt, gt] = sub_set_mva(St, dvec(St), gamma, alpha, N, latt, lidx, Ntot);
    Gnew = zeros(1,nl);
    Qnew = zeros(M, J, nl);
    Tnew = zeros(J, nl);
    Anew = zeros(M, nl);
    inner = find(inV & lvl > t);
    for v = 1:nl
        V = latt(v,:);
        vv = sum(V);
        omr = sub_omega_cum(c, t, vv);
        if omr == 0
            continue % the subnetwork is closed at this population
        end
        sub = sub_sublattice(V);
        anum = zeros(M,1);
        for q = 1:size(sub,1)
            L = sub(q,:);
            iL = lidx(sub_key(L, N));
            iVL = lidx(sub_key(V - L, N));
            pb = omr * gt(iVL) * Gin(iL);
            if pb == 0, continue; end
            Gnew(v) = Gnew(v) + pb;
            Qnew(St,:,v) = Qnew(St,:,v) + Qt(St,:,iVL) * pb;
            Tnew(:,v) = Tnew(:,v) + Tt(:,iVL) * pb;
            if ~isempty(inner)
                Qnew(inner,:,v) = Qnew(inner,:,v) + Qin(inner,:,iL) * pb;
                anum(inner) = anum(inner) + Ain(inner,iL) * pb;
            end
        end
        if Gnew(v) > 0
            Qnew(:,:,v) = Qnew(:,:,v) / Gnew(v);
            Tnew(:,v) = Tnew(:,v) / Gnew(v);
            % The SDR probability of eq. (10) carries, besides delta_ti, the
            % single-step ratios omega_{s-1,s}(v_s)/omega_ss(v_s) for every level
            % s down to the center's own. Conditioning on the population of
            % Q(V,V_t) fixes v_t, so that ratio comes out of the expectation and
            % the rest recurses through the same convolution as the queue
            % lengths. Section 4.2.3 prints only [d_1i - Q_i], which drops these
            % ratios; that is exact only when they are deterministic.
            st = sub_omega_step(c, t, vv);
            for i = St
                Anew(i,v) = st * (dvec(i) - sum(Qnew(i,:,v)));
            end
            for i = inner
                Anew(i,v) = st * anum(i) / Gnew(v);
            end
        end
    end
    Gin = Gnew; Qin = Qnew; Tin = Tnew; Ain = Anew;
end

% Sec. 4.3: convolve Q(V,V) against the complement, for every population on the
% lattice, because the per-center throughputs below read Q at N - 1_j
Qall = zeros(M, J, nl);
Tall = zeros(J, nl);
Aall = zeros(M, nl);
Gall = zeros(1, nl);
for v = 1:nl
    Np = latt(v,:);
    sub = sub_sublattice(Np);
    for q = 1:size(sub,1)
        V = sub(q,:);
        iV = lidx(sub_key(V, N));
        iC = lidx(sub_key(Np - V, N));
        pb = gc(iC) * Gin(iV);
        if pb == 0, continue; end
        Gall(v) = Gall(v) + pb;
        Qall(mv,:,v) = Qall(mv,:,v) + Qc(mv,:,iC) * pb;
        Qall(inV,:,v) = Qall(inV,:,v) + Qin(inV,:,iV) * pb;
        Tall(:,v) = Tall(:,v) + Tc(:,iC) * pb;
        Aall(inV,v) = Aall(inV,v) + Ain(inV,iV) * pb;
    end
    if Gall(v) > 0
        Qall(:,:,v) = Qall(:,:,v) / Gall(v);
        Tall(:,v) = Tall(:,v) / Gall(v);
        Aall(:,v) = Aall(:,v) / Gall(v);
    end
end

vN = lidx(sub_key(N, N));
Q = Qall(:,:,vN);
lG = log(Gall(vN));

% Per-center throughputs. Outside Q(V,V) the visit ratio is a constant, inside
% it the admission coefficient of the center's own level enters (Sec. 4.2.3,
% which prints d_1i; the level-indexed d_{t(i),i} is what eq. (10) puts in front
% of branch i and is what reproduces the exact product form)
% Shifting one chain j customer off center i in the weight of eq. (16) gives
%   alpha_i(n_i)(n_ij/n_i) w(n) = gamma_ij P_{e,e(i)}(n - e_ij) w(n - e_ij),
% so summing over states yields the exact identity
%   T_ij = xi_ij T_j(N,M) E_{N-1_j}[P_{e,e(i)}],
% the arrival-theorem statement that a departing customer sees the network at
% N - 1_j. Outside Q(V,V) the SDR probability is replaced by the constant visit
% ratio.
X = zeros(M,J);
for j = 1:J
    if N(j) == 0, continue; end
    Nm = N; Nm(j) = Nm(j) - 1;
    vm = lidx(sub_key(Nm, N));
    for i = 1:M
        if inV(i)
            X(i,j) = xi(i,j) * Tall(j,vN) * Aall(i,vm);
        else
            X(i,j) = xi(i,j) * Tall(j,vN);
        end
    end
end
U = X .* S;
R = zeros(M,J);
nz = X > 0;
R(nz) = Q(nz) ./ X(nz);
end

function [Qs, Ts, gs] = sub_set_mva(cidx, dvec, gamma, alpha, N, latt, lidx, Ntot)
% [QS, TS, GS] = SUB_SET_MVA(CIDX, DVEC, GAMMA, ALPHA, N, LATT, LIDX, NTOT)
% Section 4.2.1: MVA of the centers CIDX in isolation. DVEC(k) is the SDR
% admission coefficient of center CIDX(k), so delta_k(n) = DVEC(k) - n; NaN
% means delta = 1, which is the ordinary BCMP arrival theorem and is what the
% complement M-V uses.
%
% The factor [d_ti - Q_i(V - 1_j)] that the paper places in both W_ij and the
% denominator of T_j cancels between them, so it is not formed here: it can
% vanish at the population bound, where it would divide by zero without
% changing any queue length or throughput. It reappears only where it is a
% quantity in its own right, in the per-center throughputs of the caller.
M = size(gamma,1);
J = size(gamma,2);
nc = numel(cidx);
nl = size(latt,1);
Qs = zeros(M, J, nl);
Ts = zeros(J, nl);
gs = zeros(1, nl);
Ps = zeros(max(nc,1), Ntot+1, nl);
z = lidx(sub_key(zeros(1,J), N));
gs(z) = 1;
if nc == 0
    % An empty set holds no customers: the only feasible population is zero
    for v = 1:nl
        if sum(latt(v,:)) == 0, gs(v) = 1; else, gs(v) = 0; end
    end
    return
end
Ps(:,1,z) = 1;

[~, ord] = sort(sum(latt,2));
for oi = 1:numel(ord)
    v = ord(oi);
    V = latt(v,:);
    vv = sum(V);
    if vv == 0, continue; end
    A = zeros(nc, J);
    vm = zeros(1,J);
    for j = 1:J
        if V(j) == 0, continue; end
        Vm = V; Vm(j) = Vm(j) - 1;
        vm(j) = lidx(sub_key(Vm, N));
        for k = 1:nc
            acc = 0;
            for n = 1:vv
                df = sub_delta(dvec(k), n-1);
                if df <= 0, break; end
                acc = acc + n * (df / alpha(cidx(k), n)) * Ps(k, n, vm(j));
            end
            A(k,j) = acc;
        end
    end
    for j = 1:J
        if V(j) == 0, continue; end
        den = 0;
        for k = 1:nc
            den = den + gamma(cidx(k), j) * A(k,j);
        end
        if den > 0
            Ts(j,v) = V(j) / den;
        end
    end
    for k = 1:nc
        for j = 1:J
            if V(j) == 0, continue; end
            Qs(cidx(k), j, v) = gamma(cidx(k), j) * Ts(j,v) * A(k,j);
        end
    end
    for k = 1:nc
        tot = 0;
        for n = 1:vv
            df = sub_delta(dvec(k), n-1);
            if df <= 0, break; end
            acc = 0;
            for j = 1:J
                if V(j) == 0, continue; end
                acc = acc + gamma(cidx(k), j) * Ts(j,v) * Ps(k, n, vm(j));
            end
            Ps(k, n+1, v) = (df / alpha(cidx(k), n)) * acc;
            tot = tot + Ps(k, n+1, v);
        end
        Ps(k, 1, v) = 1 - tot;
    end
    for j = 1:J
        if V(j) > 0 && Ts(j,v) > 0
            gs(v) = gs(vm(j)) / Ts(j,v);
            break
        end
    end
end
end

function d = sub_delta(dk, n)
% D = SUB_DELTA(DK, N)
% delta_i(n) = DK - n after the rescaling to C_t = -1, or 1 outside Q(V,V).
if isnan(dk)
    d = 1;
else
    d = dk - n;
end
end

function r = sub_omega_cum(c, t, v)
% R = SUB_OMEGA_CUM(C, T, V)
% Omega_{t-1,t}(v)/Omega_tt(v), the cumulative ratio of eq. (16) at level T.
% Omega_{0,1} is one, so level 1 contributes only the reciprocal of Omega_11.
num = 1; den = 1;
for k = 0:(v-1)
    f = -k + c.Dtt(t);
    if f <= 0, r = 0; return; end
    den = den * f;
    if t > 1
        g = -k + c.Dprev(t);
        if g <= 0, r = 0; return; end
        num = num * g;
    end
end
r = num / den;
end

function r = sub_omega_step(c, t, v)
% R = SUB_OMEGA_STEP(C, T, V)
% omega_{t-1,t}(v)/omega_tt(v), the single-step ratio of eq. (10). Distinct from
% SUB_OMEGA_CUM, which is the cumulative Omega ratio of eq. (16): the routing
% probability carries the lowercase omega, the normalizing constant the
% uppercase one.
den = -v + c.Dtt(t);
if den <= 0
    r = 0;
    return
end
if t > 1
    num = -v + c.Dprev(t);
    if num <= 0
        r = 0;
        return
    end
else
    num = 1; % omega_{0,1} is one
end
r = num / den;
end

function [latt, lidx] = sub_lattice(N)
% [LATT, LIDX] = SUB_LATTICE(N)
% Every population vector V with 0 <= V <= N, and the mixed-radix map from a
% vector key to its row.
J = numel(N);
tot = prod(N+1);
latt = zeros(tot, J);
for r = 0:(tot-1)
    rem = r;
    for j = 1:J
        latt(r+1, j) = mod(rem, N(j)+1);
        rem = floor(rem / (N(j)+1));
    end
end
lidx = 1:tot;
end

function k = sub_key(V, N)
% K = SUB_KEY(V, N)
% Mixed-radix row of the population vector V in the lattice of SUB_LATTICE.
k = 0; mul = 1;
for j = 1:numel(N)
    k = k + V(j) * mul;
    mul = mul * (N(j)+1);
end
k = k + 1;
end

function sub = sub_sublattice(V)
% SUB = SUB_SUBLATTICE(V)
% Every L with 0 <= L <= V, one per row.
J = numel(V);
tot = prod(V+1);
sub = zeros(tot, J);
for r = 0:(tot-1)
    rem = r;
    for j = 1:J
        sub(r+1, j) = mod(rem, V(j)+1);
        rem = floor(rem / (V(j)+1));
    end
end
end
