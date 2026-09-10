%{
%{
 % @file pfqn_mvac.m
 % @brief MVAC (Mean Value Analysis by Chain) for product-form queueing networks.
%}
%}

%{
%{
 % @brief MVAC (Mean Value Analysis by Chain) for product-form queueing networks.
 % @fn pfqn_mvac(L, N, Z)
 % @param L Service demand matrix of the single-server fixed-rate queues (M x R).
 % @param N Population vector (1 x R).
 % @param Z Think time matrix of the infinite-server centers (1 x R or Mz x R).
 % @return XN Per-class throughput (1 x R).
 % @return QN Per-class mean queue-length at the queues (M x R).
 % @return UN Per-class utilization at the queues (M x R).
 % @return CN Per-class residence time at the queues (M x R).
%}
%}
function [XN,QN,UN,CN] = pfqn_mvac(L,N,Z)
% [XN,QN,UN,CN] = PFQN_MVAC(L,N,Z)
%
% Exact mean value analysis by chain (MVAC) of a closed multichain product-form
% queueing network composed of single-server fixed-rate (SSFR) queues and
% infinite-server (IS) centers, as given in Conway, de Souza e Silva and
% Lavenberg, "Mean Value Analysis by Chain of Product Form Queueing Networks",
% IEEE Trans. Computers 38(3):432-442, 1989.
%
% Unlike the classic MVA recursion of PFQN_MVA, which recurs on the population
% vector and therefore costs O(prod(N+1)), MVAC recurs on the *chains*: each
% chain is reduced to single-customer chains and the removed chains are replaced
% by self-looping single-customer (SCSL) chains pinned at a service center. The
% multiplicity vector v = (v_1,...,v_J), where v_j is the number of SCSL chains
% at center j, indexes the recursion in place of the population vector. MVAC is
% therefore attractive for networks with few centers and many chains, where its
% cost grows only polynomially in the number of distinct chains, and it is the
% mean-value counterpart of the RECAL normalizing-constant recursion of
% PFQN_RECAL. Since no normalizing constant is formed, MVAC does not suffer the
% floating-point underflow/overflow that complicates RECAL and convolution.
%
% Throughout, j and i index service centers (j = 1,...,J1 are SSFR and
% j = J1+1,...,J are IS), k and l index single-customer chains, a_jk = theta_jk
% T_jk is the relative utilization of chain k at center j (i.e., its service
% demand), a_k = sum_j a_jk, and I_k = {v : sum_j v_j = K - k} with K the total
% number of single-customer chains. Writing L^k_j(v) for the mean number of
% customers at center j (SCSL customers excluded), L^k_{jl}(v) for the mean
% number of chain-l customers at center j, and lambda^k_k(v) for the throughput
% of chain k, all for the network with normalizing constant G_k(v), the
% recursion of Section II of the paper reads
%   lambda^k_k(v) = 1 / (a_k + sum_{j=1}^{J1} a_jk (L^{k-1}_j(v) + v_j)),   (10)
%   L^k_{jk}(v)   = lambda^k_k(v) a_jk (1 + L^{k-1}_j(v) + v_j),  j SSFR,   (9a)
%   L^k_{jk}(v)   = lambda^k_k(v) a_jk,                           j IS,     (9b)
%   L^k_i(v)      = sum_j L^k_{jk}(v) L^{k-1}_i(v + 1_j) + L^k_{ik}(v),     (7)
%   L^k_{il}(v)   = sum_j L^k_{jk}(v) L^{k-1}_{il}(v + 1_j),  l = 1,...,k-1 (6)
% with L^0_j(v) = 0. Equation (10) is the arrival-theorem closure obtained by
% summing (9a)-(9b) over all centers, since chain k holds a single customer.
% The measures of the original network are read off at k = K and v = 0.
%
% Part 1 of the basic step evaluates (10), (9) and (7) and yields the measures
% of chain K; part 2 evaluates (6) and yields those of the chains that visit at
% least one IS center, whose throughput then follows from Little's law at that
% center. Chains that visit only SSFR centers require instead a re-execution of
% part 1 with their label interchanged with K, which is cheap because the levels
% below the interchanged label are unaffected and are reused from the first
% execution. Classes with N_r > 1, and classes with identical demand columns,
% collapse into a single subset of identical single-customer chains: only one
% representative per subset is analyzed and its per-chain measures are scaled by
% the class population, so the cost depends on the number D of *distinct* chains
% rather than on K.
%
% Parameters:
%   L - (M x R) service demand matrix of the SSFR queues.
%   N - (1 x R) closed population vector, finite and nonnegative.
%   Z - (1 x R) or (Mz x R) service demand matrix of the IS centers; each row is
%       one IS center. Defaults to zeros(1,R).
%
% Returns:
%   XN - (1 x R) per-class throughput at the reference station.
%   QN - (M x R) per-class mean queue-length at the SSFR queues.
%   UN - (M x R) per-class utilization, XN(r) * L(i,r).
%   CN - (M x R) per-class residence time, QN(i,r) / XN(r).
%
% See also PFQN_MVA, PFQN_RECAL, PFQN_CONV.
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

XN = zeros(1,R);
QN = zeros(M,R);
UN = zeros(M,R);
CN = zeros(M,R);

K = sum(N); % number of single-customer chains in the original network
if K == 0
    return
end

% Discard the centers that no chain visits: they carry no customers and would
% only inflate the multiplicity vector v.
ssfrIdx = find(any(L > 0, 2))';
isIdx = find(any(Z > 0, 2))';
J1 = length(ssfrIdx);
J = J1 + length(isIdx);
if J == 0
    line_error(mfilename,'all service demands are zero, the throughput is unbounded');
end
% A(j,r) = a_jr, centers ordered SSFR first as required by (9a)-(9b) and (10)
A = [L(ssfrIdx,:); Z(isIdx,:)];

%% Partition the chains into subsets of identical single-customer chains.
% Two chains are identical when they have the same demand at every center, so
% the subsets are the distinct columns of A restricted to the populated classes.
posr = find(N > 0);
[Adist,~,grpOfClass] = unique(A(:,posr)','rows','stable');
Dall = size(Adist,1);
visitsIS = false(1,Dall);
for g = 1:Dall
    visitsIS(g) = any(Adist(g,(J1+1):J) > 0);
end
% Chains that visit at least one IS center are labelled K-D+1,...,K-S and are
% resolved by part 2; the S chains that visit only SSFR centers are labelled
% K-S+1,...,K and are resolved by re-executions of part 1.
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
ak = sum(a,1); % (1 x K) a_k = sum_j a_jk, over ALL centers

%% Enumerate the multiplicity vectors.
% see _kb/03-api-layer.md (pfqn_mvacld -- MVAC with queue-dependent centers)
Vlist = [];
for t = 0:K
    Vlist = [Vlist; multichoose(J,t)]; %#ok<AGROW>
end
nv = size(Vlist,1);
vsum = sum(Vlist,2);
% succ(vi,j) is the index of v + 1_j, or 0 when out of range (sum(v) == K)
succ = zeros(nv,J);
for j = 1:J
    Vplus = Vlist;
    Vplus(:,j) = Vplus(:,j) + 1;
    [tf,loc] = ismember(Vplus,Vlist,'rows');
    succ(:,j) = loc .* tf;
end
z0 = find(vsum == 0,1); % index of v = 0

%% MVAC recursion.
Lall = cell(1,K+1); % Lall{k+1}(i,vi) = L^k_i(v)
Ljkall = cell(1,K+1); % Ljkall{k+1}(j,vi) = L^k_{jk}(v)
lamall = cell(1,K+1); % lamall{k+1}(vi) = lambda^k_k(v)
Lall{1} = zeros(J,nv); % L^0_i(v) = 0

lamChain = zeros(1,K); % per-chain throughput, indexed by the ORIGINAL chain label
Lchain = zeros(J,K); % per-chain mean queue-length, same indexing

% First execution of the basic step, k0 = 1.
[Lall,Ljkall,lamall] = mvac_part1(1,K,J,J1,a,ak,Vlist,vsum,succ,Lall,Ljkall,lamall);
lamChain(K) = lamall{K+1}(z0);
Lchain(:,K) = Ljkall{K+1}(:,z0);

% Part 2: measures of the chains that visit at least one IS center.
lmaxK = min(K-1,K-S);
if D >= 2 && lmaxK >= K-D+1
    L2prev = cell(1,K);
    for k = (K-D+2):K
        idxI = find(vsum == K-k); % v in I_k
        L2cur = cell(1,K);
        for l = (K-D+1):min(k-1,K-S)
            acc = zeros(length(idxI),J);
            for j = 1:J
                if l == k-1
                    % base case of (6): L^{k-1}_{i,k-1} comes from (9)
                    prev = Ljkall{k}(:,succ(idxI,j))';
                else
                    prev = L2prev{l}(:,succ(idxI,j))';
                end
                acc = acc + Ljkall{k+1}(j,idxI)' .* prev;
            end
            L2cur{l} = zeros(J,nv);
            L2cur{l}(:,idxI) = acc';
        end
        if k == K
            for l = (K-D+1):lmaxK
                Lchain(:,l) = L2cur{l}(:,z0);
                % Little's law at an IS center visited by chain l
                jIS = find(a((J1+1):J,l) > 0,1) + J1;
                lamChain(l) = Lchain(jIS,l) / a(jIS,l);
            end
        end
        L2prev = L2cur;
    end
end

% Re-executions of part 1: measures of the chains that visit only SSFR centers.
% Interchanging the labels of chains K-l and K leaves chains 1,...,K-l-1 in
% place, so the levels up to K-l-1 computed by the first execution stay valid
% and the basic step can restart at k0 = K-l.
perm = 1:K; % perm(k) is the original label of the chain now labelled k
for l = 1:(S-1)
    perm([K-l,K]) = perm([K,K-l]);
    a(:,[K-l,K]) = a(:,[K,K-l]);
    ak([K-l,K]) = ak([K,K-l]);
    [Lall,Ljkall,lamall] = mvac_part1(K-l,K,J,J1,a,ak,Vlist,vsum,succ,Lall,Ljkall,lamall);
    lamChain(perm(K)) = lamall{K+1}(z0);
    Lchain(:,perm(K)) = Ljkall{K+1}(:,z0);
end

%% Expand the per-chain measures back to the per-class measures.
% Every chain of a subset has the same measures, hence a class with N_r chains
% in that subset carries N_r times the measures of the subset representative.
for g = 1:D
    kg = K-D+g; % original label of the representative of subset g
    for r = posr(grpOfClass == gorder(g))
        XN(r) = N(r) * lamChain(kg);
        QN(ssfrIdx,r) = N(r) * Lchain(1:J1,kg);
        UN(ssfrIdx,r) = XN(r) * L(ssfrIdx,r);
        CN(ssfrIdx,r) = QN(ssfrIdx,r) / XN(r);
    end
end
% An empty class has zero throughput and queue-length, so its residence time is
% reported as the bare service demand, as in PFQN_MVA.
CN(:,N == 0) = L(:,N == 0);
end

function [Lall,Ljkall,lamall] = mvac_part1(k0,K,J,J1,a,ak,Vlist,vsum,succ,Lall,Ljkall,lamall)
% Part 1 of the basic step: evaluate (10), (9a)-(9b) and (7) for k = k0,...,K
% over v in V_k = I_k u ... u I_K.
for k = k0:K
    idxV = find(vsum <= K-k);
    Lp = Lall{k}; % L^{k-1}
    Vv = Vlist(idxV,:);
    Lpv = Lp(:,idxV)'; % L^{k-1}_j(v)
    ajk = a(:,k)'; % a_jk
    % (10): the chain-k customer is at some center, so sum_j L^k_{jk}(v) = 1
    den = ak(k) + sum((Lpv(:,1:J1) + Vv(:,1:J1)) .* ajk(1:J1), 2);
    lam = 1 ./ den;
    Ljkv = zeros(length(idxV),J);
    Ljkv(:,1:J1) = (lam .* (1 + Lpv(:,1:J1) + Vv(:,1:J1))) .* ajk(1:J1); % (9a)
    Ljkv(:,(J1+1):J) = lam .* ajk((J1+1):J); % (9b)
    % (7)
    Lkv = Ljkv;
    for j = 1:J
        Lkv = Lkv + Ljkv(:,j) .* Lp(:,succ(idxV,j))';
    end
    Lk = zeros(J,length(vsum));
    Lk(:,idxV) = Lkv';
    Ljk = zeros(J,length(vsum));
    Ljk(:,idxV) = Ljkv';
    lamv = zeros(1,length(vsum));
    lamv(idxV) = lam;
    Lall{k+1} = Lk;
    Ljkall{k+1} = Ljk;
    lamall{k+1} = lamv;
end
end
