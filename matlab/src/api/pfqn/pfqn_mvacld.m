%{
%{
 % @file pfqn_mvacld.m
 % @brief MVAC (Mean Value Analysis by Chain) for load-dependent networks.
%}
%}

%{
%{
 % @brief MVAC (Mean Value Analysis by Chain) for load-dependent networks.
 % @fn pfqn_mvacld(L, N, Z, mu)
 % @param L Service demand matrix of the queue-length dependent centers (M x R).
 % @param N Population vector (1 x R).
 % @param Z Think time matrix of the infinite-server centers (1 x R or Mz x R).
 % @param mu Load-dependent rate matrix (M x Nt), mu(j,n) = rate with n jobs.
 % @return XN Per-class throughput (1 x R).
 % @return QN Per-class mean queue-length at the centers (M x R).
 % @return UN Utilization of each center (M x 1), i.e. 1-P_j(0).
 % @return CN Per-class cycle time exclusive of think time (1 x R).
 % @return pij Marginal queue-length probabilities (M x (Nt+1)).
%}
%}
function [XN,QN,UN,CN,pij] = pfqn_mvacld(L,N,Z,mu)
% [XN,QN,UN,CN,PIJ] = PFQN_MVACLD(L,N,Z,MU)
%
% Exact mean value analysis by chain (MVAC) of a closed multichain product-form
% queueing network that may contain queue-length dependent (QLD) service
% centers. This is the extension of Section V of Conway, de Souza e Silva and
% Lavenberg, "Mean Value Analysis by Chain of Product Form Queueing Networks",
% IEEE Trans. Computers 38(3):432-442, 1989; PFQN_MVAC implements Sections II-IV,
% which cover single-server fixed-rate (SSFR) and infinite-server (IS) centers
% only.
%
% Where PFQN_MVAC propagates the MEAN queue-lengths L^k_j(v) through eq. (7) and
% closes the recursion with the arrival-theorem identity (10), the QLD extension
% propagates the MARGINAL queue-length DISTRIBUTIONS P^k_j(n,v) instead. That is
% forced by load dependence -- the rate seen by a job depends on the whole
% occupancy, so a mean no longer suffices -- but it also SIMPLIFIES the
% recursion: eq. (21)-(25) read level k-1 only at the shifted vectors v + 1_i,
% so the basic step sweeps v in I_k alone, where PFQN_MVAC must sweep the larger
% I_k u ... u I_K. The extra information comes almost for free, and PIJ is
% returned as a first-class output.
%
% Notation follows PFQN_MVAC: j and i index centers (j = 1,...,J1 are the QLD
% centers of L, j = J1+1,...,J the IS centers of Z), k and l index the
% single-customer chains, a_jk = theta_jk T_jk is the demand of chain k at
% center j, K = sum(N) and I_k = {v : sum_j v_j = K - k} with v_j the number of
% self-looping single-customer (SCSL) chains pinned at center j. P^k_j(n,v) is
% the probability of n customers at center j -- EXCLUDING the v_j SCSL customers
% there -- in the network with normalizing constant G_k(v). Writing tau_k(v,i)
% for the throughput of an SCSL chain that replaces chain k at center i, the
% recursion of Section V reads
%   tau_k(v,i)  = T_ik^-1 sum_{n=0}^{k-1} P^{k-1}_i(n,v+1_i) mu_i(n+v_i+1)
%                                                          / (n+v_i+1),      (21)
%   tau_k(v,i)  = T_ik^-1,                                       i IS,       (22)
%   L^k_{jk}(v) = theta_jk tau_k(v,j)^-1 / sum_m theta_mk tau_k(v,m)^-1,     (23)
%   lambda^k_k(v) = tau_k(v,j(k)) L^k_{j(k)k}(v),                            (24)
%   P^k_j(n,v)  = L^k_{jk}(v) P^{k-1}_j(n-1,v+1_j)
%                 + sum_{m ~= j} L^k_{mk}(v) P^{k-1}_j(n,v+1_m),             (25)
% with P^0_j(0,v) = 1 and P^{k-1}_j(n,.) = 0 for n < 0 or n > k-1. Eq. (21) is
% just "the mean rate at which an SCSL chain is served": given n other customers
% the processor-sharing rate share is mu_i(n+v_i+1)/(n+v_i+1), averaged over the
% distribution of those others. The queueing discipline may be assumed PS with
% no loss of generality, since product-form measures do not depend on it.
%
% This implementation writes (23)-(24) in the reference-station-free form
%   c_i(k,v)      = sum_{n=0}^{k-1} P^{k-1}_i(n,v+1_i) mu_i(n+v_i+1)/(n+v_i+1),
%   L^k_{jk}(v)   = (a_jk/c_j) / sum_m (a_mk/c_m),
%   lambda^k_k(v) = 1 / sum_m (a_mk/c_m),
% which follows from theta_jk tau_k(v,j)^-1 = a_jk/c_j and theta_{j(k)k} = 1, so
% only the demands a_jk are needed and the visit ratios never appear separately.
% For an IS center c_i = 1 identically, which is exactly (22). Eq. (25) is
% self-normalizing, sum_n P^k_j(n,v) = sum_m L^k_{mk}(v) = 1, so no normalizing
% constant is formed and the recursion involves only positive quantities: unlike
% the classic load-dependent MVA of PFQN_MVALD it cannot produce negative
% probabilities and needs no stabilization.
%
% Parts 2 and 3 of the algorithm are unchanged from PFQN_MVAC, since eq. (6)
% holds verbatim in the presence of QLD centers: part 2 resolves the chains that
% visit at least one IS center, and the chains that visit no IS center are
% resolved by re-executing part 1 with their label interchanged with K.
%
% Parameters:
%   L  - (M x R) service demand matrix of the QLD centers.
%   N  - (1 x R) closed population vector, finite and nonnegative.
%   Z  - (1 x R) or (Mz x R) demand matrix of the IS centers, one row per
%        center. Defaults to zeros(1,R).
%   mu - (M x Nt) load-dependent rates, Nt >= sum(N); mu(j,n) is the total
%        service rate of center j with n jobs present. mu(j,:) = 1 is a
%        single-server fixed-rate queue, mu(j,n) = min(n,c) a c-server queue,
%        mu(j,n) = n an infinite server. Defaults to ones(M,sum(N)), i.e. all
%        centers SSFR, in which case the results agree with PFQN_MVAC.
%
% Returns:
%   XN  - (1 x R) per-class throughput at the reference station.
%   QN  - (M x R) per-class mean queue-length at the QLD centers.
%   UN  - (M x 1) utilization of each center, 1 - P_j(0). PER-STATION, not
%         per-class, as in PFQN_MVALD and PFQN_DAC: for a load-dependent center
%         the per-class product XN(r)*L(j,r) of PFQN_MVAC is NOT the utilization.
%   CN  - (1 x R) per-class cycle time exclusive of think time, N(r)/XN(r)-Z(r),
%         as in PFQN_MVALD and PFQN_DAC. NOT the (M x R) per-station residence
%         time of PFQN_MVAC: the whole LD family reports a cycle time here.
%   pij - (M x (sum(N)+1)) marginal queue-length probabilities,
%         pij(j,n+1) = P(n jobs at center j).
%
% See also PFQN_MVAC, PFQN_MVALD, PFQN_DAC, PFQN_GLD.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[M,R] = size(L);
if nargin < 3 || isempty(Z)
    Z = zeros(1,R);
end
if size(Z,2) ~= R
    line_error(mfilename,'the think time matrix and the demand matrix have a different number of classes');
end
N = round(N(:)');
if length(N) ~= R
    line_error(mfilename,'the demand matrix and the population vector have a different number of classes');
end
if any(N < 0) || any(isinf(N))
    line_error(mfilename,'the population vector must be finite and nonnegative');
end
K = sum(N); % number of single-customer chains in the original network
if nargin < 4 || isempty(mu)
    mu = ones(M,max(K,1));
end
if size(mu,1) ~= M
    line_error(mfilename,'the rate matrix and the demand matrix have a different number of centers');
end
if K > 0 && size(mu,2) < K
    line_error(mfilename,'the rate matrix must supply a rate for every population up to sum(N)');
end

XN = zeros(1,R);
QN = zeros(M,R);
UN = zeros(M,1);
CN = zeros(1,R);
pij = zeros(M,K+1);
pij(:,1) = 1; % an unvisited center holds no jobs
if K == 0
    return
end

% Discard the centers that no chain visits: they carry no customers and would
% only inflate the multiplicity vector v.
ldIdx = find(any(L > 0, 2))';
isIdx = find(any(Z > 0, 2))';
J1 = length(ldIdx);
J = J1 + length(isIdx);
if J == 0
    line_error(mfilename,'all service demands are zero, the throughput is unbounded');
end
% A(j,r) = a_jr, centers ordered QLD first, then IS, as required by (21)-(22)
A = [L(ldIdx,:); Z(isIdx,:)];
MU = mu(ldIdx,1:K); % (J1 x K) rates of the QLD centers
if any(MU(:) <= 0)
    line_error(mfilename,'the service rates must be strictly positive for every population up to sum(N)');
end

%% Partition the chains into subsets of identical single-customer chains.
% see _kb/03-api-layer.md (pfqn_mvacld -- MVAC with queue-dependent centers)
posr = find(N > 0);
[Adist,~,grpOfClass] = unique(A(:,posr)','rows','stable');
Dall = size(Adist,1);
visitsIS = false(1,Dall);
for g = 1:Dall
    visitsIS(g) = any(Adist(g,(J1+1):J) > 0);
end
% see _kb/03-api-layer.md (pfqn_mvacld -- MVAC with queue-dependent centers)
gorder = [find(visitsIS), find(~visitsIS)];
D = Dall;
S = sum(~visitsIS);
Agrp = Adist(gorder,:)'; % (J x D) demands of the representative of each subset
mult = zeros(1,D);
for g = 1:D
    mult(g) = sum(N(posr(grpOfClass == gorder(g))));
end

% Chain labels: the D representatives are labelled K-D+1,...,K, the remaining
% K-D identical chains fill the labels 1,...,K-D in any order.
chainGroup = zeros(1,K);
chainGroup((K-D+1):K) = 1:D;
p = 0;
for g = 1:D
    for c = 1:(mult(g)-1)
        p = p + 1;
        chainGroup(p) = g;
    end
end
a = Agrp(:,chainGroup); % (J x K) relative utilizations a_jk

%% Enumerate the multiplicity vectors.
% see _kb/03-api-layer.md (pfqn_mvacld -- MVAC with queue-dependent centers)
Vlist = [];
cnt = zeros(1,K+1);
off = zeros(1,K+1);
for t = 0:K
    Vt = multichoose(J,t);
    off(t+1) = size(Vlist,1);
    cnt(t+1) = size(Vt,1);
    Vlist = [Vlist; Vt]; %#ok<AGROW>
end
nv = size(Vlist,1);
% succ(vi,i) is the global index of v + 1_i, or 0 when out of range (sum(v)==K)
succ = zeros(nv,J);
for j = 1:J
    Vplus = Vlist;
    Vplus(:,j) = Vplus(:,j) + 1;
    [tf,loc] = ismember(Vplus,Vlist,'rows');
    succ(:,j) = loc .* tf;
end

%% MVAC recursion.
Pall = cell(1,K+1); % Pall{k+1}(j,vloc,n+1) = P^k_j(n,v), v the vloc-th of I_k
Ljkall = cell(1,K+1); % Ljkall{k+1}(j,vloc) = L^k_{jk}(v)
lamall = cell(1,K+1); % lamall{k+1}(vloc) = lambda^k_k(v)
Pall{1} = ones(J1,cnt(K+1),1); % P^0_j(0,v) = 1, I_0 is the block of sum K

lamChain = zeros(1,K); % per-chain throughput, indexed by the ORIGINAL chain label
Lchain = zeros(J,K); % per-chain mean queue-length, same indexing

% First execution of the basic step, k0 = 1.
[Pall,Ljkall,lamall] = mvacld_part1(1,K,J,J1,a,MU,Vlist,off,cnt,succ,Pall,Ljkall,lamall);
lamChain(K) = lamall{K+1}(1);
Lchain(:,K) = Ljkall{K+1}(:,1);
% The marginals of the ORIGINAL network are read at k = K and v = 0; the label
% interchanges below overwrite this level, so capture them now.
pij(ldIdx,:) = reshape(Pall{K+1}(:,1,:),J1,K+1);

% Part 2: measures of the chains that visit at least one IS center. Eq. (6) is
% unchanged by load dependence.
lmaxK = min(K-1,K-S);
if D >= 2 && lmaxK >= K-D+1
    L2prev = cell(1,K);
    for k = (K-D+2):K
        t = K-k;
        L2cur = cell(1,K);
        for l = (K-D+1):min(k-1,K-S)
            acc = zeros(J,cnt(t+1));
            for vloc = 1:cnt(t+1)
                vi = off(t+1) + vloc;
                for j = 1:J
                    sloc = succ(vi,j) - off(t+2);
                    if l == k-1
                        % base case of (6): L^{k-1}_{i,k-1} comes from (23)
                        prev = Ljkall{k}(:,sloc);
                    else
                        prev = L2prev{l}(:,sloc);
                    end
                    acc(:,vloc) = acc(:,vloc) + Ljkall{k+1}(j,vloc) * prev;
                end
            end
            L2cur{l} = acc;
        end
        if k == K
            for l = (K-D+1):lmaxK
                Lchain(:,l) = L2cur{l}(:,1);
                % Little's law at an IS center visited by chain l
                jIS = find(a((J1+1):J,l) > 0,1) + J1;
                lamChain(l) = Lchain(jIS,l) / a(jIS,l);
            end
        end
        L2prev = L2cur;
    end
end

% see _kb/03-api-layer.md (pfqn_mvacld -- MVAC with queue-dependent centers)
perm = 1:K; % perm(k) is the original label of the chain now labelled k
for l = 1:(S-1)
    perm([K-l,K]) = perm([K,K-l]);
    a(:,[K-l,K]) = a(:,[K,K-l]);
    [Pall,Ljkall,lamall] = mvacld_part1(K-l,K,J,J1,a,MU,Vlist,off,cnt,succ,Pall,Ljkall,lamall);
    lamChain(perm(K)) = lamall{K+1}(1);
    Lchain(:,perm(K)) = Ljkall{K+1}(:,1);
end

%% Expand the per-chain measures back to the per-class measures.
% see _kb/03-api-layer.md (pfqn_mvacld -- MVAC with queue-dependent centers)
for g = 1:D
    kg = K-D+g; % original label of the representative of subset g
    for r = posr(grpOfClass == gorder(g))
        XN(r) = N(r) * lamChain(kg);
        QN(ldIdx,r) = N(r) * Lchain(1:J1,kg);
    end
end
% Utilization of a load-dependent center is 1-P_j(0), NOT XN(r)*L(j,r).
UN = 1 - pij(:,1);
% Cycle time exclusive of think time, as in PFQN_MVALD and PFQN_DAC. An empty
% class has zero throughput, hence no cycle time.
CN = zeros(1,R);
CN(posr) = N(posr) ./ XN(posr) - sum(Z(:,posr),1);
end

function [Pall,Ljkall,lamall] = mvacld_part1(k0,K,J,J1,a,MU,Vlist,off,cnt,succ,Pall,Ljkall,lamall)
% Part 1 of the basic step for networks with QLD centers: evaluate (21)-(25) for
% k = k0,...,K over v in I_k, the block of multiplicity vectors of sum K-k.
for k = k0:K
    t = K - k;
    nvk = cnt(t+1);
    Pp = Pall{k}; % P^{k-1}, over I_{k-1}, third dimension n = 0,...,k-1
    Lk = zeros(J,nvk);
    lam = zeros(1,nvk);
    Pk = zeros(J1,nvk,k+1);
    for vloc = 1:nvk
        vi = off(t+1) + vloc;
        v = Vlist(vi,:);
        sloc = succ(vi,:) - off(t+2); % local index of v + 1_i within I_{k-1}
        % see _kb/03-api-layer.md (pfqn_mvacld -- MVAC with queue-dependent centers)
        c = ones(1,J);
        for i = 1:J1
            ci = 0;
            for n = 0:(k-1)
                ci = ci + Pp(i,sloc(i),n+1) * MU(i,n+v(i)+1) / (n+v(i)+1);
            end
            c(i) = ci;
        end
        % (23)-(24) in reference-station-free form: theta_jk/tau_k(v,j) = a_jk/c_j
        w = zeros(J,1);
        for j = 1:J
            if a(j,k) > 0
                w(j) = a(j,k) / c(j);
            end
        end
        sw = sum(w);
        lam(vloc) = 1 / sw;
        Lk(:,vloc) = w / sw;
        % (25): condition on the center holding the single chain-k customer
        for j = 1:J1
            for n = 0:k
                s = 0;
                if n >= 1
                    % chain k is at j, so n-1 of the k-1 others are there too
                    s = Lk(j,vloc) * Pp(j,sloc(j),n);
                end
                if n <= k-1
                    % chain k is elsewhere, so all n are from the k-1 others
                    for m = 1:J
                        if m ~= j
                            s = s + Lk(m,vloc) * Pp(j,sloc(m),n+1);
                        end
                    end
                end
                Pk(j,vloc,n+1) = s;
            end
        end
    end
    Pall{k+1} = Pk;
    Ljkall{k+1} = Lk;
    lamall{k+1} = lam;
end
end
