function pfqn_sens_respt_validate()
%{
%{
 % @file pfqn_sens_respt_validate.m
 % @brief Validation harness for pfqn_sens_respt, the FCFS sojourn-time moment
 %        analysis of Strelen (1990), Theorem 4.1.
 %
 %        Five references:
 %          A. brute-force enumeration. This is the strongest check because it
 %             shares none of Theorem 4.1's algebra. The equilibrium product
 %             form is enumerated to get the exact arrival-theorem marginals
 %             p_i(j, N-1_l); the sojourn time conditioned on finding j jobs is
 %             known in closed form (Exp(mu) if j < b, otherwise an
 %             Erlang(j-b+1, b*mu) queueing delay plus an Exp(mu) service), so
 %             its moments are formed directly and mixed over j. Neither the
 %             coefficients a_(t,tau)(0) nor the recursion (4.2) enter, so the
 %             agreement tests both.
 %          B. the published table of Example 3.4 (continued) of the reference,
 %             which prints E(W_i) and sigma^2_(W_i) for the Kobayashi model.
 %          C. the internal identity W(i,l) = w_i(l)/V(i,l): the t = 1 case of
 %             (4.5) must reproduce the MVA residence time divided by the visit
 %             ratio, which is a completely different expression.
 %          D. pfqn_mva for the base measures in the single-server case.
 %          E. the single-job network, where an arriving job always finds an
 %             empty station, so W is exactly Exp(mu) and every moment is known
 %             in closed form.
 %
 %        Reference: J. C. Strelen, "Moment Analysis for Closed Queuing Networks
 %        and its Linearizer", Performance Evaluation 11:127-142, 1990.
%}
%}
rng(5);
tolBrute = 1e-9;
tolPaper = 5e-4;   % the paper prints 5 decimals on small sojourn times
tolIdent = 1e-10;
tolMva   = 1e-10;
tolExp   = 1e-12;

errBrute = 0; errIdent = 0; errMva = 0; errExp = 0;
nBrute = 0;

% =====================================================================
% A/C/D. random models, single- and multi-server
% =====================================================================
for trial = 1:36
    M = randi([1 3]);
    R = randi([1 2]);
    S = 0.2 + rand(M,1);
    V = 0.3 + rand(M,R);
    if mod(trial,4) == 0 && M > 1
        V(1,1) = 0;      % a class that skips a station
    end
    N = randi([1 3],1,R);
    if mod(trial,2) == 0
        Z = 0.3 + rand(1,R);
    else
        Z = zeros(1,R);
    end
    if mod(trial,3) == 0
        b = randi([1 3],M,1);
    else
        b = ones(M,1);
    end

    res = pfqn_sens_respt(S,V,N,Z,b,3);

    % ---- C. W = w/V ------------------------------------------------
    for i = 1:M
        for l = 1:R
            if V(i,l) > 0 && N(l) > 0
                errIdent = max(errIdent, relerr(res.W(i,l), res.Wresid(i,l)/V(i,l)));
            end
        end
    end

    % ---- D. base measures, single-server only ----------------------
    if all(b == 1)
        L = zeros(M,R);
        for i = 1:M
            L(i,:) = S(i) * V(i,:);
        end
        [XN,QN,UN] = pfqn_mva(L,N,Z);
        errMva = max([errMva, relerr(res.X,XN), relerr(res.Q,QN), relerr(res.U,UN)]);
    end

    % ---- A. brute force --------------------------------------------
    if prod(N+1) <= 24 && M <= 3
        [Wb, pb] = brute_respt(S,V,N,Z,b,3);
        errBrute = max(errBrute, relerr(res.WM, Wb));
        % .p is ragged: station i only defines j = 0..b_i-1, the range the
        % b-server recursion needs, and the rest of the row is zero padding out
        % to max(b). Comparing the padding against the true marginal would be
        % comparing against something the algorithm never claims to compute.
        for i = 1:M
            errBrute = max(errBrute, relerr(res.p(i,1:b(i)), pb(i,1:b(i))));
        end
        nBrute = nBrute + 1;
    end
end

% =====================================================================
% E. one job: the arriving job always finds the station empty, so W~Exp(mu)
% =====================================================================
S1 = [0.4; 0.25];
V1 = [1; 2];
r1 = pfqn_sens_respt(S1,V1,1,0.7,[1;1],3);
for i = 1:2
    mu = 1/S1(i);
    errExp = max([errExp, relerr(r1.WM(i,1,1), 1/mu), ...
                  relerr(r1.WM(i,1,2), 2/mu^2), relerr(r1.WM(i,1,3), 6/mu^3), ...
                  relerr(r1.WVar(i,1), 1/mu^2)]);
end

% =====================================================================
% B. Example 3.4 (continued): Kobayashi central-server model, n = 3
%    The paper prints E(W_i) and sigma^2_(W_i).
% =====================================================================
xs = [repmat(0.0215,1,9), 0.104, 0.104, 0.019]';
es = [repmat(9.333,1,9),  10.5,  10.5,  105]';
rk = pfqn_sens_respt(xs, es, 3, 0, ones(12,1), 3);
paperW   = [0.02275; 0.14178; 0.03322];
paperWV  = [0.00052; 0.01846; 0.00083];
gotW  = [rk.W(1,1);    rk.W(10,1);    rk.W(12,1)];
gotWV = [rk.WVar(1,1); rk.WVar(10,1); rk.WVar(12,1)];
errPaper = max(relerr(gotW,paperW), relerr(gotWV,paperWV));

fprintf('\n=== pfqn_sens_respt validation (max relative error) ===\n');
fprintf('  A. brute force, arrival theorem (%d models) : %.3e  (tol %.1e)\n', nBrute, errBrute, tolBrute);
fprintf('  B. Strelen Example 3.4 published sojourn    : %.3e  (tol %.1e)\n', errPaper, tolPaper);
fprintf('  C. identity W(i,l) = w_i(l)/V(i,l)          : %.3e  (tol %.1e)\n', errIdent, tolIdent);
fprintf('  D. pfqn_mva base measures (b=1)             : %.3e  (tol %.1e)\n', errMva, tolMva);
fprintf('  E. single job, W ~ Exp(mu) exactly          : %.3e  (tol %.1e)\n', errExp, tolExp);
fprintf('     (paper E(W_12)=%.5f got %.5f ; sigma2=%.5f got %.5f)\n', ...
        paperW(3), gotW(3), paperWV(3), gotWV(3));

ok = errBrute <= tolBrute && errPaper <= tolPaper && errIdent <= tolIdent && ...
     errMva <= tolMva && errExp <= tolExp;
if ~ok
    error('pfqn_sens_respt_validate:mismatch','one or more checks exceeded tolerance');
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
function [WM, pN] = brute_respt(S,V,N,Z,b,tmax)
% Sojourn-time moments from first principles: enumerate the product form to get
% the exact arrival-theorem marginals p_i(j, N-e_l), then mix the conditional
% sojourn-time moments over j. Uses none of Theorem 4.1.
[M,R] = size(V);
bmax = max(b);
WM = zeros(M,R,tmax);
for l = 1:R
    if N(l) == 0
        continue;
    end
    Nl = N; Nl(l) = Nl(l) - 1;
    pj = brute_marginals(S,V,Nl,Z,b);      % pj(i,1+j) at population N - e_l
    for i = 1:M
        if V(i,l) <= 0
            continue;
        end
        mu = 1/S(i);
        for t = 1:tmax
            acc = 0;
            for j = 0:(size(pj,2)-1)
                acc = acc + pj(i,1+j) * cond_moment(j,b(i),mu,t);
            end
            WM(i,l,t) = acc;
        end
    end
end
pAll = brute_marginals(S,V,N,Z,b);
pN = zeros(M,bmax);
for i = 1:M
    for j = 0:(bmax-1)
        pN(i,1+j) = pAll(i,1+j);
    end
end
end

% =========================================================================
function v = cond_moment(j,b,mu,t)
% E[(W|j)^t] where a job arriving to find j jobs at an FCFS b-server station
% waits an Erlang(max(0,j-b+1), b*mu) and is then served for an Exp(mu).
k = max(0, j - b + 1);
v = 0;
for s = 0:t
    % E[X^s] with X ~ Exp(mu)
    EX = factorial(s) / mu^s;
    p = t - s;
    % E[Y^p] with Y ~ Erlang(k, b*mu), and Y = 0 when k = 0
    if k == 0
        if p == 0
            EY = 1;
        else
            EY = 0;
        end
    else
        theta = b*mu;
        EY = 1;
        for a = 0:(p-1)
            EY = EY * (k + a);
        end
        EY = EY / theta^p;
    end
    v = v + nchoosek(t,s) * EX * EY;
end
end

% =========================================================================
function pj = brute_marginals(S,V,N,Z,b)
% P[Q_i = j] for every station, by enumerating the closed product form of a
% network of FCFS b-server stations:
%   f_i(q_i) = q_i! prod_l (a(i,l)^q_il / q_il!) prod_{j=1}^{q_i} 1/min(j,b_i)
% with a(i,l) = S(i)*V(i,l), plus the delay term for the think times.
[M,R] = size(V);
a = zeros(M,R);
for i = 1:M
    a(i,:) = S(i) * V(i,:);
end
states = enumerate_states(N,M);
K = size(states,1);
w = zeros(K,1);
tot = zeros(K,M);
for k = 1:K
    nir = reshape(states(k,:),R,M)';
    lw = 0;
    ok = true;
    for i = 1:M
        ni = sum(nir(i,:));
        lw = lw + gammaln(ni+1);
        for j = 1:ni
            lw = lw - log(min(j,b(i)));
        end
        for r = 1:R
            if nir(i,r) > 0
                if a(i,r) <= 0
                    ok = false; break;
                end
                lw = lw + nir(i,r)*log(a(i,r)) - gammaln(nir(i,r)+1);
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
    tot(k,:) = sum(nir,2)';
end
w = w / sum(w);
maxj = max(sum(N),1);
pj = zeros(M,1+maxj);
for k = 1:K
    for i = 1:M
        pj(i,1+tot(k,i)) = pj(i,1+tot(k,i)) + w(k);
    end
end
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
