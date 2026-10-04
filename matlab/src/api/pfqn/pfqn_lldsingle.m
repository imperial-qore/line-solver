%{
%{
 % @file pfqn_lldsingle.m
 % @brief Exact normalizing constant for single-class limited load-dependent models.
%}
%}

%{
%{
 % @brief Exact normalizing constant for single-class limited load-dependent models.
 %
 % Same recursion, same arithmetic and same result as pfqn_gldsingle, but with
 % the auxiliary rate-offset axis truncated at the LIMITED LOAD-DEPENDENCE
 % threshold instead of at the population. Unrolling the recursion of
 % pfqn_gldsingle shows that its third index is an offset into station m's rate
 % function,
 %
 %    g(m,n,t) = sum_{j=0..n} prod_{i=0..j-1} L_m/alpha_m(t+i) * g(m-1,n-j,1)
 %
 % so once t >= s_m, where s_m is the population past which alpha_m stays
 % constant, every factor is alpha_m(s_m) and the product collapses to
 % (L_m/alpha_m(s_m))^j. Hence
 %
 %    g(m,n,t) = g(m,n,s_m)   for all t >= s_m
 %
 % and the N-s_m upper slices that pfqn_gldsingle computes are duplicates of
 % one another. Capping the offset at s_m and reading g(m,n-1,min(t+1,s_m))
 % keeps every value it needs.
 %
 % COST. O(N * sum_k s_k) time against O(M N^2) for pfqn_gldsingle, and
 % O(N * max_k s_k) space against O(M N^2), the levels being rolled. On a
 % multiserver model, where s_k is the server count, this is LINEAR in the
 % population rather than quadratic. The two agree to the last bit, since the
 % arithmetic performed is a subset of pfqn_gldsingle's: the log-domain branch
 % remains a sum of nonnegative terms and loses no digits to cancellation,
 % unlike the closed form of pfqn_explicit_ld, which reaches the same
 % asymptotics through an alternating sum.
 %
 % There is no gain on a station whose rates never settle, an infinite server
 % alpha(n)=n being the usual case: it gets s_k = N and costs what it costs in
 % pfqn_gldsingle. The saving is over the OTHER stations, so a model carrying
 % one delay among M queues drops from O(M N^2) to O(N^2 + N sum_k s_k).
 %
 % The threshold is detected per station rather than declared, so an arbitrary
 % rate matrix is accepted and simply yields s_k = N, at which point this is
 % pfqn_gldsingle with its slices rolled. A MISSED tie only costs time; a FALSE
 % tie would be a wrong answer, hence the strict tolerance, which follows
 % pfqn_explicit_ld's scan.
 %
 % @fn pfqn_lldsingle(L, N, mu, options)
 % @param L Service demand vector (Mx1).
 % @param N Population (scalar).
 % @param mu Load-dependent rate matrix (MxN), alpha_i(j) = mu(i,j).
 % @param options Solver options.
 % @return lG Logarithm of normalizing constant.
 % @return G Normalizing constant.
 % @return s Detected per-station thresholds (Mx1), alpha_i(n)=alpha_i(s_i) for n>=s_i.
%}
%}
function [lG,G,s]=pfqn_lldsingle(L,N,mu,options)
% G=PFQN_LLDSINGLE(L,N,MU)

if nargin<4
    options = [];
end

[M,R]=size(L);
if R>1
    line_error(mfilename,'multiclass model detected. pfqn_lldsingle is for single class models.');
end
Nscal = N(1); % codegen: ensure scalar loop bound
if Nscal<=0
    % the empty product, as pfqn_gldsingle returns from its unexecuted loops
    lG = 0; G = 1; s = ones(M,1);
    return
end

% ---- s_k: the smallest offset past which the rate row is constant ----
% Equality is tested first so that an infinite rate, which pfqn_gldsingle
% admits and zeroes through log(mu)=+Inf, ties with itself instead of
% producing Inf-Inf. The tolerance is then confined to FINITE pairs: at
% tail=Inf the bound eps*max(abs(tail),1) is itself Inf and abs(prev-Inf)<=Inf
% would tie every finite rate to it, collapsing the row on a false tie.
s = Nscal*ones(M,1);
for m=1:M
    tail = mu(m,Nscal);
    for n=Nscal:-1:2
        if mu(m,n-1)==tail || (isfinite(tail) && isfinite(mu(m,n-1)) ...
                && abs(mu(m,n-1)-tail) <= eps*max(abs(tail),1))
            s(m) = n-1;
        else
            break
        end
    end
end
s = max(s,1);

% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
useLog = isreal(L) && isreal(mu) && all(L(:)>=0) && all(mu(:)>0);

if useLog
    lL = log(L);       % -Inf where the demand is zero
    lmu = log(mu);     % +Inf where the rate is infinite, zeroing the term
    lgprev = -Inf(Nscal+1,1);
    lgprev(0 +1) = 0;  % g(0,0,1)=1, and g(0,n,1)=0 for n>=1: no station holds them
    for m=1:M
        sm = s(m);
        cur = -Inf(Nscal+1, sm);
        for tm=1:sm
            cur(0 +1, tm) = 0; % log(1): zero jobs
        end
        for n=1:Nscal
            % offsets above Nscal-n+1 are never read back, exactly as in
            % pfqn_gldsingle, so the triangle is kept
            for tm=1:min(sm,Nscal-n+1)
                a = lgprev(n +1);
                b = lL(m) + cur(n-1 +1, min(tm+1,sm)) - lmu(m,tm);
                % pairwise log-sum-exp of a and b, stable when either is -Inf
                if a > b
                    if b == -Inf
                        cur(n +1, tm) = a;
                    else
                        cur(n +1, tm) = a + log1p(exp(b-a));
                    end
                else
                    if a == -Inf
                        cur(n +1, tm) = b;
                    else
                        cur(n +1, tm) = b + log1p(exp(a-b));
                    end
                end
            end
        end
        lgprev = cur(:,1);
    end
    lG = lgprev(Nscal +1);
    G = exp(lG);
else
    gprev = zeros(Nscal+1,1);
    gprev(0 +1) = 1;
    for m=1:M
        sm = s(m);
        cur = zeros(Nscal+1, sm);
        for tm=1:sm
            cur(0 +1, tm) = 1;
        end
        for n=1:Nscal
            for tm=1:min(sm,Nscal-n+1)
                cur(n +1, tm) = gprev(n +1) + L(m)*cur(n-1 +1, min(tm+1,sm))/mu(m,tm);
            end
        end
        gprev = cur(:,1);
    end
    G = gprev(Nscal +1);
    lG = log(G);
end
end
