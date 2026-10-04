function pfqn_sens_mom_validate()
%{
%{
 % @file pfqn_sens_mom_validate.m
 % @brief Validation harness for pfqn_sens_mom, the higher-moment analysis of
 %        Strelen (1990).
 %
 %        Five references:
 %          A. brute-force enumeration of the closed product-form distribution,
 %             which is ground truth for m, Var, Cov, E[Q^2] and E[Q^3].
 %          B. pfqn_sens_mva. Summing its per-class covariance matrix at station
 %             i over all class pairs must give Var[Q_i], since
 %             Var[sum_r n(i,r)] = sum_{r,s} Cov[n(i,r),n(i,s)]. This ties the
 %             per-station-total moments of Strelen to the finer per-class
 %             moments of de Souza e Silva and Muntz.
 %          C. pfqn_mva for the base measures.
 %          D. Cov symmetry: x_j dm_i/dx_j and x_i dm_j/dx_i are computed by
 %             different derivative tracks and must agree.
 %          E. the published table of Example 3.4 of the reference (the
 %             Kobayashi central-server model), which pins the second
 %             derivative against numbers the author printed rather than
 %             against our own code.
 %
 %        Reference: J. C. Strelen, "Moment Analysis for Closed Queuing Networks
 %        and its Linearizer", Performance Evaluation 11:127-142, 1990.
%}
%}
rng(3);
tolBrute = 1e-9;
tolMva   = 1e-10;
tolTot   = 1e-9;
tolSym   = 1e-9;
tolPaper = 5e-5;   % the paper prints 5 significant digits

errBrute = 0; errMva = 0; errTot = 0; errSym = 0; errGrp = 0;
nBrute = 0; nGrp = 0;

for trial = 1:40
    M = randi([1 3]);
    R = randi([1 2]);
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
    if mod(trial,3) == 0
        mi = randi([1 3],1,M);
    else
        mi = ones(1,M);
    end

    mom = pfqn_sens_mom(L,N,Z,mi);

    % ---- C. base measures ------------------------------------------
    [XN,QN,UN,CN] = pfqn_mva(L,N,Z,mi);
    errMva = max([errMva, relerr(mom.X,XN), relerr(mom.Q,QN), ...
                  relerr(mom.U,UN), relerr(mom.R,CN)]);

    % ---- D. symmetry -----------------------------------------------
    errSym = max(errSym, mom.CovAsym);

    % ---- B. total variance against the per-class covariances -------
    ref = pfqn_sens_mva(L,N,Z,mi);
    errTot = max(errTot, relerr(mom.Var, ref.QTotVar));

    % ---- A. brute force --------------------------------------------
    if prod(N+1) <= 32 && M <= 3 && all(mi == 1)
        [mb,Varb,Covb,M2b,M3b] = brute_totals(L,N,Z);
        errBrute = max([errBrute, relerr(mom.m,mb), relerr(mom.Var,Varb), ...
                        relerr(mom.Cov,Covb), relerr(mom.M2,M2b), ...
                        relerr(mom.M3,M3b)]);
        nBrute = nBrute + 1;

        % ---- F. the per-class grouping ------------------------------
        % groups = 1:R scales one class at a time, which is Akyildiz and
        % Strelen's Theorem 1 with T = {r}. It must reproduce the per-class
        % moments of the brute-force distribution, INCLUDING the third, and its
        % second moments must equal pfqn_sens_mva's exactly.
        momc = pfqn_sens_mom(L,N,Z,mi,1:R);
        [mc,Varc,M2c,M3c] = brute_perclass(L,N,Z);
        errGrp = max([errGrp, relerr(momc.m,mc), relerr(momc.Var,Varc), ...
                      relerr(momc.M2,M2c), relerr(momc.M3,M3c)]);
        errGrp = max(errGrp, relerr(momc.Var, ref.QVar));
        nGrp = nGrp + 1;

        % ---- G. an intermediate grouping ----------------------------
        % Two classes in one group must give the moments of their sum, and a
        % grouping that puts every class in one group must reproduce the
        % default (station totals). Both are checked against brute force.
        if R == 2
            momg = pfqn_sens_mom(L,N,Z,mi,[1 1]);
            errGrp = max([errGrp, relerr(momg.m, mom.m), ...
                          relerr(momg.Var, mom.Var), relerr(momg.M3, mom.M3)]);
        end
    end
end

% =====================================================================
% E. Example 3.4 of the reference: Kobayashi central-server model.
%    12 type-1 queues, one class. Queues 1-9: x=0.0215, e=9.333;
%    queues 10,11: x=0.104, e=10.5; queue 12: x=0.019, e=105.
%    The paper prints E(Q_i) and sigma^2_{Q_i} for n=3, 2, 1.
% =====================================================================
xs = [repmat(0.0215,1,9), 0.104, 0.104, 0.019];
es = [repmat(9.333,1,9),  10.5,  10.5,  105];
Lk = (xs .* es)';          % demands, 12 x 1
paperM   = [0.07606, 0.05835, 0.03353;   % rows: queues 1-9, 10-11, 12
            0.53316, 0.36327, 0.18246;
            1.24917, 0.74835, 0.33334];
paperVar = [0.07893, 0.05873, 0.03240;
            0.57689, 0.34341, 0.14917;
            1.02546, 0.56250, 0.22222];
errPaper = 0;
for col = 1:3
    nJobs = 4 - col;   % col 1 -> n=3, col 2 -> n=2, col 3 -> n=1
    mk = pfqn_sens_mom(Lk, nJobs, 0);
    got = [mk.m(1), mk.m(10), mk.m(12)];
    gotV = [mk.Var(1), mk.Var(10), mk.Var(12)];
    errPaper = max([errPaper, relerr(got(:), paperM(:,col)), ...
                    relerr(gotV(:), paperVar(:,col))]);
    % the nine identical queues must be identical, and so must 10 and 11
    errPaper = max([errPaper, relerr(mk.m(1:9), repmat(mk.m(1),9,1)), ...
                    relerr(mk.m(10), mk.m(11))]);
end

fprintf('\n=== pfqn_sens_mom validation (max relative error) ===\n');
fprintf('  A. brute-force product form (%d models) : %.3e  (tol %.1e)\n', nBrute, errBrute, tolBrute);
fprintf('  B. Var vs pfqn_sens_mva QTotVar         : %.3e  (tol %.1e)\n', errTot, tolTot);
fprintf('  C. pfqn_mva base measures               : %.3e  (tol %.1e)\n', errMva, tolMva);
fprintf('  D. Cov symmetry (raw, pre-symmetrize)   : %.3e  (tol %.1e)\n', errSym, tolSym);
fprintf('  E. Strelen Example 3.4 published table  : %.3e  (tol %.1e)\n', errPaper, tolPaper);
fprintf('  F. per-class grouping vs brute force + pfqn_sens_mva (%d models): %.3e  (tol %.1e)\n', nGrp, errGrp, tolBrute);

ok = errBrute <= tolBrute && errTot <= tolTot && errMva <= tolMva && ...
     errSym <= tolSym && errPaper <= tolPaper && errGrp <= tolBrute;
if ~ok
    error('pfqn_sens_mom_validate:mismatch','one or more checks exceeded tolerance');
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
function [m,Var,Cov,M2,M3] = brute_totals(L,N,Z)
% Moments of the per-station total queue lengths by enumerating the closed
% product-form equilibrium distribution.
[M,R] = size(L);
states = enumerate_states(N,M);
K = size(states,1);
w = zeros(K,1);
for k = 1:K
    nir = reshape(states(k,:),R,M)';
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
    end
end
w = w / sum(w);

tot = zeros(K,M);
for k = 1:K
    nir = reshape(states(k,:),R,M)';
    tot(k,:) = sum(nir,2)';
end
m = (w' * tot)';
M2 = (w' * (tot.^2))';
M3 = (w' * (tot.^3))';
Var = M2 - m.^2;
Cov = zeros(M,M);
for i = 1:M
    for j = 1:M
        Cov(i,j) = sum(w .* tot(:,i) .* tot(:,j)) - m(i)*m(j);
    end
end
end

% =========================================================================
function [m,Var,M2,M3] = brute_perclass(L,N,Z)
% Per-class moments of n(i,r) by enumerating the closed product form.
[M,R] = size(L);
states = enumerate_states(N,M);
K = size(states,1);
w = zeros(K,1);
for k = 1:K
    nir = reshape(states(k,:),R,M)';
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
    end
end
w = w / sum(w);
m = zeros(M,R); M2 = zeros(M,R); M3 = zeros(M,R);
for k = 1:K
    nir = reshape(states(k,:),R,M)';
    m = m + w(k)*nir;
    M2 = M2 + w(k)*nir.^2;
    M3 = M3 + w(k)*nir.^3;
end
Var = M2 - m.^2;
end

% =========================================================================
function states = enumerate_states(N,M)
R = numel(N);
per = cell(1,R);
for r = 1:R
    per{r} = compositions_leq(N(r),M);
end
states = zeros(0,R*M);
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
