function pfqn_sens_mva_validate()
%{
%{
 % @file pfqn_sens_mva_validate.m
 % @brief Validation harness for pfqn_sens_mva. Checks the MVA-like moment
 %        recursion of de Souza e Silva and Muntz (1988), Corollary 1, against
 %        three independent references:
 %
 %          A. brute-force enumeration of the closed product-form equilibrium
 %             distribution (ground truth, mi==1);
 %          B. the differentiated-MVA Jacobian of pfqn_sens, via the identity
 %             Cov[n(i,r),n(i,s)] = L(i,s) * dQ(i,r)/dL(i,s) (covers mi>1);
 %          C. pfqn_mva for the base measures X, Q, U, R.
 %
 %        It also reports the raw asymmetry of QCov before symmetrization: the
 %        recursion computes W(k,j;t,j) and W(t,j;k,j) by numerically distinct
 %        expressions, so their agreement is a nontrivial check of the formula.
%}
%}
rng(0);
tolBrute = 1e-9;
tolSens  = 1e-9;
tolMva   = 1e-10;
tolSym   = 1e-9;

errBrute = 0; errSens = 0; errMva = 0; errSym = 0;
nBrute = 0; nSens = 0;

for trial = 1:40
    M = randi([1 3]);
    R = randi([1 3]);
    L = 0.2 + rand(M,R);
    N = randi([0 3],1,R);
    if ~any(N > 0)
        N(1) = 2;
    end
    if mod(trial,2) == 0
        Z = 0.3 + rand(1,R);
    else
        Z = zeros(1,R);
    end
    % exercise a zero-demand column now and then: class r never visits station i
    if mod(trial,5) == 0 && M > 1
        L(1,1) = 0;
    end

    mom = pfqn_sens_mva(L,N,Z);

    % ---- C. base measures against pfqn_mva -----------------------------
    [XN,QN,UN,CN] = pfqn_mva(L,N,Z);
    errMva = max([errMva, relerr(mom.X,XN), relerr(mom.Q,QN), ...
                  relerr(mom.U,UN), relerr(mom.R,CN)]);

    % ---- A. brute force -------------------------------------------------
    if prod(N+1) <= 64 && M <= 3
        [Qb,QCovb] = brute_moments(L,N,Z);
        errBrute = max([errBrute, relerr(mom.Q,Qb), relerr(mom.QCov,QCovb)]);
        nBrute = nBrute + 1;
    end

    % ---- B. pfqn_sens Jacobian, and the raw asymmetry -------------------
    for micase = 1:2
        if micase == 1
            mi = ones(1,M);
        else
            mi = randi([1 3],1,M);
        end
        momi = pfqn_sens_mva(L,N,Z,mi);
        sens = pfqn_sens(L,N,Z,mi);
        % Read the reference off the raw Jacobian, NOT off sens.QCov: pfqn_sens
        % now sources its same-station blocks from pfqn_sens_mva, so comparing
        % against sens.QCov would compare the recursion with itself.
        pL = zeros(M,R);
        for p = 1:numel(sens.params)
            if sens.params(p).type == 'L'
                pL(sens.params(p).station,sens.params(p).class) = p;
            end
        end
        CovRef = zeros(M,R,R);
        for i = 1:M
            for r = 1:R
                for s = 1:R
                    if pL(i,s) > 0
                        CovRef(i,r,s) = L(i,s) * sens.dQ(i,r,pL(i,s));
                    end
                end
            end
        end
        errSens = max(errSens, relerr(momi.QCov,CovRef));
        errSym = max(errSym, momi.QCovAsym);
        % Theorem 3: variance of the total queue length at a station
        TotRef = zeros(M,1);
        for i = 1:M
            TotRef(i) = sum(sum(CovRef(i,:,:)));
        end
        errSens = max(errSens, relerr(momi.QTotVar,TotRef));
        nSens = nSens + 1;
    end
end

fprintf('\n=== pfqn_sens_mva validation (max relative error) ===\n');
fprintf('  A. brute-force product form (%d models) : %.3e  (tol %.1e)\n', nBrute, errBrute, tolBrute);
fprintf('  B. pfqn_sens Jacobian      (%d models) : %.3e  (tol %.1e)\n', nSens, errSens, tolSens);
fprintf('  C. pfqn_mva base measures              : %.3e  (tol %.1e)\n', errMva, tolMva);
fprintf('  D. QCov symmetry (raw, pre-symmetrize) : %.3e  (tol %.1e)\n', errSym, tolSym);

ok = errBrute <= tolBrute && errSens <= tolSens && errMva <= tolMva && errSym <= tolSym;
if ~ok
    error('pfqn_sens_mva_validate:mismatch','one or more checks exceeded tolerance');
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
function [Q,QCov] = brute_moments(L,N,Z)
% Exact moments by enumerating the closed product-form equilibrium
% distribution. Stations 1..M are single-server fixed-rate centers; the think
% time Z is an infinite-server station indexed 0 and carries no moment.
%   p(n) ~ prod_i [ n_i! prod_r L(i,r)^n(i,r)/n(i,r)! ] * prod_r Z(r)^n(0,r)/n(0,r)!
[M,R] = size(L);
states = enumerate_states(N,M);   % each row: [n(1,1..R) n(2,1..R) ... n(M,1..R)]
K = size(states,1);
w = zeros(K,1);
for k = 1:K
    nir = reshape(states(k,:),R,M)';   % M x R
    lw = 0;
    ok = true;
    for i = 1:M
        ni = sum(nir(i,:));
        lw = lw + gammaln(ni+1);
        for r = 1:R
            if nir(i,r) > 0
                if L(i,r) <= 0
                    ok = false; break;
                end
                lw = lw + nir(i,r)*log(L(i,r)) - gammaln(nir(i,r)+1);
            end
        end
        if ~ok, break; end
    end
    if ok
        for r = 1:R
            n0r = N(r) - sum(nir(:,r));
            if n0r > 0
                if Z(r) <= 0
                    ok = false; break;
                end
                lw = lw + n0r*log(Z(r)) - gammaln(n0r+1);
            end
        end
    end
    if ok
        w(k) = exp(lw);
    else
        w(k) = 0;
    end
end
w = w / sum(w);

Q = zeros(M,R);
QCov = zeros(M,R,R);
for k = 1:K
    nir = reshape(states(k,:),R,M)';
    Q = Q + w(k)*nir;
end
for i = 1:M
    for r = 1:R
        for s = 1:R
            m2 = 0;
            for k = 1:K
                nir = reshape(states(k,:),R,M)';
                m2 = m2 + w(k)*nir(i,r)*nir(i,s);
            end
            QCov(i,r,s) = m2 - Q(i,r)*Q(i,s);
        end
    end
end
end

% =========================================================================
function states = enumerate_states(N,M)
% All allocations of N(r) class-r jobs over M stations (the remainder sits in
% the delay). Returns a matrix whose rows are [n(1,:) n(2,:) ... n(M,:)] with
% the inner index running over stations for each class, flattened class-major.
R = numel(N);
per = cell(1,R);
for r = 1:R
    per{r} = compositions_leq(N(r),M);   % rows: n(1..M,r) with sum <= N(r)
end
states = zeros(1,R*M);
states = states([],:);
idx = ones(1,R);
while true
    row = zeros(M,R);
    for r = 1:R
        row(:,r) = per{r}(idx(r),:)';
    end
    states(end+1,:) = reshape(row',1,[]); %#ok<AGROW>
    r = R;
    while r >= 1
        idx(r) = idx(r) + 1;
        if idx(r) <= size(per{r},1)
            break;
        end
        idx(r) = 1;
        r = r - 1;
    end
    if r == 0
        break;
    end
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
