%{
%{
 % @file pfqn_dac.m
 % @brief DAC (Distribution Analysis by Chain) method for joint queue-length distributions.
%}
%}

%{
%{
 % @brief DAC (Distribution Analysis by Chain) method for joint queue-length distributions.
 % @fn pfqn_dac(L, N, Z, mu)
 % @param L Service demand matrix (MxR).
 % @param N Population vector (1xR).
 % @param Z Think time vector (1xR, default: zeros).
 % @param mu Load-dependent rate matrix (MxNt, default: ones).
 % @return Pjoint Joint queue-length probabilities over the aggregate state space.
 % @return states Aggregate states, one per row of Pjoint.
 % @return XN Chain throughputs.
 % @return QN Mean queue lengths.
 % @return UN Station utilizations.
 % @return CN Cycle times.
 % @return pi Marginal queue-length probabilities.
%}
%}
function [Pjoint,states,XN,QN,UN,CN,pi] = pfqn_dac(L,N,Z,mu)
% [PJOINT,STATES,XN,QN,UN,CN,PI]=PFQN_DAC(L,N,Z,MU)
%
% Distribution Analysis by Chain (DAC) for closed product-form queueing
% networks with single-server fixed-rate, infinite-server and queue-dependent
% service centers. Unlike MVA, RECAL or MVAC, the recursion returns the whole
% set of joint queue-length probabilities, which are required e.g. in
% availability modeling.
%
% The recursion proceeds chain by chain over a related network in which every
% chain holds a single customer, a transformation that leaves the aggregate
% queue-length distribution unchanged. Given the distribution of a network
% with k-1 such chains, adding one customer of a chain with demands r gives
%
%   c_j       = sum_{n=1..k} (n/mu_j(n)) * P_j^{k-1}(n-1)
%   lambda_k  = 1 / sum_j r_j c_j
%   P^k(n)    = lambda_k * sum_j r_j (n_j/mu_j(n_j)) * P^{k-1}(n-e_j)
%
% where lambda_k is the throughput of the customer being added and 1/c_j is
% the throughput of a chain visiting center j only. The recursion conserves
% probability mass by construction, hence it is numerically stable.
%
% Inputs:
%   L  (MxR)  service demand of chain r at station j
%   N  (1xR)  number of customers of chain r
%   Z  (1xR)  think time of chain r. If sum(Z)>0 an extra infinite-server
%             station is appended, so that STATES has M+1 columns and its
%             last column holds the think-station population.
%   mu (MxNt) service rate of station j with n customers, Nt=sum(N).
%             Defaults to ones(M,Nt), i.e. single-server fixed rate. Use
%             mu(j,:)=1:Nt for infinite server and mu(j,n)=min(n,c) for a
%             c-server station.
%
% Outputs:
%   Pjoint (SxR) probability of the aggregate state in the corresponding row
%          of STATES, S=nchoosek(Nt+J-1,J-1) with J the number of centers
%   states (SxJ) aggregate states, states(s,j) = customers at center j
%   XN     (1xR) throughput of chain r
%   QN     (MxR) mean number of chain-r customers at station j
%   UN     (Mx1) utilization of station j, i.e. 1-P_j(0)
%   CN     (1xR) cycle time of chain r, exclusive of think time
%   pi     (Mx(Nt+1)) pi(j,n+1) = marginal probability of n jobs at station j
%
% Example (availability model, de Souza e Silva 1987, Section 3):
%   L = [5,0; 0,10; 2,1]; N = [1,3]; mu = [1,1,1,1; 1,2,2,2; 1,1,1,1];
%   [P,S] = pfqn_dac(L,N,[0,0],mu);
%   AV = sum(P(S(:,1)==1 & S(:,2)>=1)); % system available
%
% References:
% E. de Souza e Silva, "Distribution Analysis of Product Form Queueing
% Networks", UCLA Computer Science Department, CSD-870023, April 1987.

[M,R] = size(L);
if nargin<3 || isempty(Z)
    Z = zeros(1,R);
end
N = N(:)';
Z = Z(:)';
Nt = sum(N);
if nargin<4 || isempty(mu)
    mu = ones(M,max(Nt,1));
end

if any(N<0)
    line_error(mfilename,'Population vector must be non-negative.');
end
if size(mu,1)~=M
    line_error(mfilename,'The mu matrix must have one row per station.');
end
if Nt>0 && size(mu,2)<Nt
    line_error(mfilename,'The mu matrix must have at least sum(N) columns.');
end

% A non-zero think time is modelled as an appended infinite-server center.
hasZ = sum(Z)>0;
if hasZ
    Lx = [L; Z];
    mux = [mu(:,1:max(Nt,1)); 1:max(Nt,1)];
else
    Lx = L;
    mux = mu(:,1:max(Nt,1));
end
J = size(Lx,1);

% Trivial population
if Nt==0
    states = zeros(1,J);
    Pjoint = 1;
    XN = zeros(1,R);
    QN = zeros(M,R);
    UN = zeros(M,1);
    CN = zeros(1,R);
    pi = zeros(M,1);
    pi(:,1) = 1;
    return
end

active = find(N>0);
D = numel(active);
for idx=1:D
    if all(Lx(:,active(idx))<=0)
        line_error(mfilename,'Chain %d has null demand at every center.',active(idx));
    end
end

% Enumerate the aggregate state space level by level and precompute, for each
% state at level k, the index of that state with one more job at center j.
[lvstates,succ] = dac_lattice(J,Nt);

% Chains are ordered so that the last D added are distinct, which lets all
% per-chain measures reuse the common prefix of the recursion.
prefix = zeros(1,0);
for idx=1:D
    prefix = [prefix, repmat(active(idx),1,N(active(idx))-1)]; %#ok<AGROW>
end
tail = active;

% Prefix of the recursion, shared by every per-chain run.
p = 1;
k = 0;
for idx=1:numel(prefix)
    [p,~,~] = dac_step(p,Lx(:,prefix(idx)),mux,k,lvstates,succ);
    k = k+1;
end

% Base run over the tail, saving the intermediate distributions S_0..S_{D-1}.
Sp = cell(1,D);
Sp{1} = p;
pb = p;
kb = k;
for idx=1:D
    [pb,lam,Lq] = dac_step(pb,Lx(:,tail(idx)),mux,kb,lvstates,succ);
    kb = kb+1;
    if idx<D
        Sp{idx+1} = pb;
    end
end
Pjoint = pb;
states = lvstates{Nt+1};

XN = zeros(1,R);
QN = zeros(M,R);
% The base run already places chain tail(D) last.
r = tail(D);
XN(r) = N(r)*lam;
QN(:,r) = N(r)*Lq(1:M);

% Re-run the tail with chain tail(idx) moved last, restarting from S_{idx}.
for idx=1:D-1
    pc = Sp{idx};
    kc = k+idx-1;
    order = [tail(idx+1:D), tail(idx)];
    for t=1:numel(order)
        [pc,lam,Lq] = dac_step(pc,Lx(:,order(t)),mux,kc,lvstates,succ);
        kc = kc+1;
    end
    r = tail(idx);
    XN(r) = N(r)*lam;
    QN(:,r) = N(r)*Lq(1:M);
end

% Marginal queue-length probabilities at the full population.
pi = zeros(M,Nt+1);
for j=1:M
    pi(j,:) = accumarray(states(:,j)+1,Pjoint,[Nt+1,1])';
end
UN = 1-pi(:,1);
CN = zeros(1,R);
CN(active) = N(active)./XN(active) - Z(active);
end

% Enumerate the aggregate state space and the successor index map.
function [lvstates,succ] = dac_lattice(J,Nt)
lvstates = cell(1,Nt+1);
succ = cell(1,Nt+1);
for k=0:Nt
    lvstates{k+1} = dac_compositions(J,k);
end
% Binomial table for ranking compositions in lexicographic order.
C = zeros(Nt+J,J+1);
for a=0:(Nt+J-1)
    for b=0:min(a,J)
        if b==0
            C(a+1,b+1) = 1;
        else
            C(a+1,b+1) = C(a,b) + C(a,b+1);
        end
    end
end
for k=0:Nt-1
    S = lvstates{k+1};
    nS = size(S,1);
    ix = zeros(nS,J);
    for i=1:nS
        for j=1:J
            t = S(i,:);
            t(j) = t(j)+1;
            ix(i,j) = dac_rank(t,k+1,J,C);
        end
    end
    succ{k+1} = ix;
end
end

% All J-part compositions of k, in lexicographic order.
function A = dac_compositions(J,k)
if J==1
    A = k;
    return
end
A = zeros(0,J);
for v=0:k
    B = dac_compositions(J-1,k-v);
    A = [A; repmat(v,size(B,1),1), B]; %#ok<AGROW>
end
end

% Lexicographic rank of composition n of k into J parts.
function idx = dac_rank(n,k,J,C)
idx = 1;
rem = k;
for j=1:J-1
    for v=0:n(j)-1
        s = rem-v;
        idx = idx + C(s+J-j,J-j);
    end
    rem = rem-n(j);
end
end

% One step of the recursion: add a single customer with demands r to a network
% holding k customers, returning the distribution at level k+1, the throughput
% of the added customer and its per-center presence probabilities.
function [pn,lam,Lq] = dac_step(p,r,mu,k,lvstates,succ)
J = numel(r);
S = lvstates{k+1};
% Marginal queue lengths of the network with k customers.
marg = zeros(J,k+1);
for j=1:J
    marg(j,:) = accumarray(S(:,j)+1,p,[k+1,1])';
end
c = zeros(J,1);
for j=1:J
    for n=1:k+1
        c(j) = c(j) + (n/mu(j,n))*marg(j,n);
    end
end
lam = 1/(r(:)'*c);
Lq = lam*r(:).*c;
ix = succ{k+1};
pn = zeros(size(lvstates{k+2},1),1);
for j=1:J
    if r(j)<=0
        continue
    end
    nj = S(:,j)+1;
    w = lam*r(j)*(nj./mu(j,nj)').*p;
    pn = pn + accumarray(ix(:,j),w,[numel(pn),1]);
end
end
