%{
%{
 % @file pfqn_ncldmx.m
 % @brief Normalizing constant and mean measures for mixed open/closed networks with limited load dependence.
%}
%}

%{
%{
 % @brief Normalizing constant and mean measures for mixed open/closed networks with limited load dependence.
 %
 % The closed-conditional normalizing constant of a mixed limited load-dependent
 % (LLD) network equals a purely closed load-dependent normalizing constant in
 % which every queueing station i carries the Bruell-Balbo-Afshari effective
 % capacity rate mu_i^eff(n) = 1/EC_i(n), where EC is returned by
 % pfqn_ldmx_ec and folds the open classes into the closed subnetwork. The
 % open classes contribute the separable prefactor lGopen = sum_i log E_i(0),
 % which reduces to -sum_i log(1-rho_i) in the load-independent limit.
 %
 % Mean measures follow from that identification without ever enumerating the
 % closed population lattice, which is what makes this the normalizing-constant
 % counterpart of pfqn_mvaldmx rather than a rename of it:
 %
 %   - closed throughputs are the ratios X_r = G(N-e_r)/G(N);
 %   - closed queue lengths are the conditional normalizing-constant recursion
 %     of the load-dependent closed network (pfqn_mushift / pfqn_fnc), applied
 %     to the effective-capacity rates;
 %   - open queue lengths are the Bruell-Balbo-Afshari sum
 %     Q_ir = lambda_r D_ir sum_n (n+1) EC_i(n+1) P_i(n) with its SATURATED TAIL
 %     FOLDED ONTO THE CLOSED MEAN. EC_i(n) is constant for n >= b_i, the level
 %     where the rate row stops growing, so writing EC_i(n) = EC_i^inf + delta_i(n)
 %     with delta_i(n)=0 for n >= b_i leaves
 %       Q_ir = lambda_r D_ir [ EC_i^inf (Q_i^closed + 1)
 %                              + sum_{n=0}^{b_i-2} (n+1) delta_i(n+1) P_i(n) ]
 %     using sum_n P_i(n)=1 and sum_n n P_i(n)=Q_i^closed. Only the first b_i-1
 %     marginal probabilities survive, and b_i is the number of servers, not the
 %     population: a single-server station needs none at all, and the formula
 %     collapses to the classical lambda_r D_ir (1+Q_i^closed)/(1-rho_i).
 %
 % The marginals that remain are themselves normalizing-constant ratios,
 % P_i(n) = sum_{|k|=n} F_i(k) G_{-i}(N-k) / G(N), evaluated as in
 % solver_nc_margaggr.
 %
 % Reproduces pfqn_mvaldmx to machine precision on multiserver, nonlinear-mu,
 % think-time and multi-chain models; see test_pfqn_ncldmx.
 %
 % @fn pfqn_ncldmx(lambda, D, N, Z, mu, S, varargin)
 % @param lambda Arrival rate vector (0 on closed classes).
 % @param D Service demand matrix (MxR).
 % @param N Population vector (Inf on open classes).
 % @param Z Think time vector (closed classes).
 % @param mu Load-dependent rate matrix (Mx>=sum(N_closed)).
 % @param S Number of servers per station (currently informational).
 % @param varargin Optional solver parameters forwarded to pfqn_ncld.
 % @return lG Logarithm of the closed-conditional normalizing constant.
 % @return G Closed-conditional normalizing constant (exp(lG)).
 % @return lGopen Logarithm of the open-class normalizing prefactor sum_i log E_i(0).
 % @return XN Throughputs (1xR): G(N-e_r)/G(N) on closed classes, lambda_r on open ones.
 % @return QN Mean queue lengths (MxR).
%}
%}
function [lG,G,lGopen,XN,QN] = pfqn_ncldmx(lambda,D,N,Z,mu,S,varargin)
% [LG,G,LGOPEN,XN,QN] = PFQN_NCLDMX(LAMBDA,D,N,Z,MU,S,VARARGIN)

[M,R] = size(D);
if nargin<5 || isempty(mu)
    mu = ones(M,max(1,sum(N(isfinite(N)))));
end
if nargin<6 || isempty(S)
    S = ones(M,1); %#ok<NASGU> % kept for signature parity with pfqn_mvaldmx
end
if nargin<4 || isempty(Z)
    Z = zeros(1,R);
end

openClasses = find(isinf(N));
closedClasses = setdiff(1:R, openClasses);

if any(N(intersect(find(lambda),closedClasses))>0)
    line_error(mfilename,'Arrival rate cannot be specified on closed classes.');
end

Nc = N(closedClasses);
Kc = sum(Nc);

% open-class normalizing prefactor sum_i log E_i(0); needs the effective
% capacity terms even when there are no closed jobs.
mup = mu;
if size(mup,2) < max(1,Kc)
    mup = [mup, repmat(mup(:,end), 1, max(1,Kc)-size(mup,2))];
end
mup = [mup, mup(:,end)]; % one extra column as in pfqn_mvaldmx
lambdao = zeros(1,R); lambdao(openClasses) = lambda(openClasses);
[EC,E] = pfqn_ldmx_ec(lambdao, D, mup);
lGopen = sum(log(E(1:M,1)));

% closed-conditional normalizing constant
Dc = D(:,closedClasses);
Zc = Z(closedClasses);
muEff = 1 ./ EC(:,1:max(1,Kc));
if Kc==0
    lG = 0;
    G = 1;
else
    [lG,G] = pfqn_ncld(Dc, Nc, Zc, muEff, varargin{:});
end

if nargout <= 3
    return
end

%% mean measures
XN = zeros(1,R);
QN = zeros(M,R);
XN(openClasses) = lambda(openClasses);

Cc = numel(closedClasses);
lGr = zeros(1,Cc);
if Kc > 0
    % closed throughputs: X_r = G(N-e_r)/G(N)
    for rc=1:Cc
        if Nc(rc) <= 0
            continue
        end
        lGr(rc) = ncld_local(Dc, oner(Nc,rc), Zc, muEff, varargin);
        XN(closedClasses(rc)) = exp(lGr(rc) - lG);
    end
    % closed queue lengths: conditional normalizing-constant recursion of the
    % load-dependent closed network, on the effective-capacity rates
    for ist=1:M
        if ~any(Dc(ist,:) > 0)
            continue
        end
        muhat = pfqn_mushift(muEff, ist);
        [muhat_f, cshift] = pfqn_fnc(muhat(ist,:));
        Dminus = Dc; Dminus(ist,:) = [];
        muminus = muEff; muminus(ist,:) = [];
        for rc=1:Cc
            if Nc(rc) <= 0 || Dc(ist,rc) <= 0
                continue
            end
            Ncr = oner(Nc,rc);
            lGhat = ncld_local(Dc, Ncr, Zc, muhat, varargin);
            lGhatf = ncld_local([Dc; Dc(ist,:)], Ncr, Zc, [muhat; muhat_f], varargin);
            lGminus = ncld_local(Dminus, Ncr, Zc, muminus, varargin);
            CQ = (exp(lGhatf - lGhat) - 1) + cshift(1)*(exp(lGminus - lGhat) - 1);
            ldDemand = log(Dc(ist,rc)) + lGhat - log(muEff(ist,1)) - lGr(rc);
            QN(ist,closedClasses(rc)) = exp(ldDemand) * XN(closedClasses(rc)) * (1 + CQ);
        end
    end
end

% open queue lengths, with the saturated tail of EC folded onto the closed mean
if ~isempty(openClasses)
    if Cc > 0
        Qtot = sum(QN(:,closedClasses),2);
    else
        Qtot = zeros(M,1);
    end
    for ist=1:M
        b = lld_level(mup(ist,:));
        ECinf = EC(ist, min(b, size(EC,2)));
        acc = ECinf * (Qtot(ist) + 1);
        if b >= 2
            if Kc > 0
                Dminus = Dc; Dminus(ist,:) = [];
                muminus = muEff; muminus(ist,:) = [];
            end
            for n=0:(b-2)
                delta = EC(ist,n+1) - ECinf;
                if delta == 0
                    continue
                end
                % WITH NO CLOSED POPULATION THE MARGINAL IS DEGENERATE, not
                % absent: P_i(0)=1 and P_i(n)=0 above it, so only the n=0 term
                % survives and acc collapses to EC_i(1), the exact open
                % load-dependent mean. Skipping the loop instead left acc at
                % EC_i^inf, i.e. read a c-server station as if every arrival
                % found it saturated -- an M/M/3 at lambda=1.5 came back with
                % mean 1 against the exact 1.7368.
                if Kc == 0
                    Pn = double(n == 0);
                else
                    Pn = marginal_local(n, ist, Dc, Nc, Zc, muEff, Dminus, muminus, lG, varargin);
                end
                acc = acc + (n+1) * delta * Pn;
            end
        end
        for r=openClasses
            QN(ist,r) = lambda(r) * D(ist,r) * acc;
        end
    end
end
end

function lG = ncld_local(L, N, Z, mu, opts)
% pfqn_ncld, with the empty-station residual network handled explicitly: a
% network reduced to its think times alone has G(N) = prod_r Z_r^N_r / N_r!,
% and no G at all when a class has jobs but neither demand nor think time.
if size(L,1) == 0
    Zt = sum(Z,1);
    if all(N <= 0)
        lG = 0;
        return
    end
    if all(Zt <= 0)
        lG = -Inf;
        return
    end
    lG = 0;
    for r=1:numel(N)
        if N(r) <= 0
            continue
        end
        if Zt(r) <= 0
            lG = -Inf;
            return
        end
        lG = lG + N(r)*log(Zt(r)) - factln(N(r));
    end
    return
end
lG = pfqn_ncld(L, N, Z, mu, opts{:});
end

function P = marginal_local(n, ist, Dc, Nc, Zc, muEff, Dminus, muminus, lG, opts)
% P_i(n) = sum_{|k|=n, k<=Nc} F_i(k) G_{-i}(Nc-k) / G(Nc), where F_i(k) is the
% one-station constant of station i at population k (solver_nc_margaggr's lF_i).
ks = compositions_local(n, Nc);
P = 0;
for t=1:size(ks,1)
    k = ks(t,:);
    if n == 0
        lF = 0;
    else
        lF = ncld_local(Dc(ist,:), k, 0*k, muEff(ist,:), opts);
    end
    lGbar = ncld_local(Dminus, Nc-k, Zc, muminus, opts);
    P = P + exp(lF + lGbar - lG);
end
end

function ks = compositions_local(n, Nc)
% Non-negative integer vectors k with sum(k)==n and k<=Nc, one per row.
Cc = numel(Nc);
if Cc == 0
    if n == 0
        ks = zeros(1,0);
    else
        ks = zeros(0,0);
    end
    return
end
if Cc == 1
    if n <= Nc(1)
        ks = n;
    else
        ks = zeros(0,1);
    end
    return
end
ks = zeros(0,Cc);
for v=0:min(n,Nc(1))
    sub = compositions_local(n-v, Nc(2:end));
    if ~isempty(sub)
        ks = [ks; repmat(v,size(sub,1),1), sub]; %#ok<AGROW>
    end
end
end

function b = lld_level(murow)
% First column of the trailing constant run of a limited load-dependence row,
% i.e. the level b with mu(n)=mu(b) for every n>=b. This is the level
% pfqn_ldmx_ec infers, and hence the one past which EC is constant.
b = numel(murow);
if b == 0
    b = 1;
    return
end
while b > 1 && murow(b-1) == murow(b)
    b = b - 1;
end
end
