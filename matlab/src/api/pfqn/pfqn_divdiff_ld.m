%{
%{
 % @file pfqn_divdiff_ld.m
 % @brief Exact load-dependent normalizing constant by divided differences of a single-class constant.
%}
%}

%{
%{
 % @brief Exact load-dependent normalizing constant by divided differences of a single-class constant.
 %
 % Evaluates the multiclass load-dependent normalizing constant as the
 % N-th order finite difference, in the class populations, of a SINGLE-CLASS
 % load-dependent constant:
 %
 %    G(N) = sum_{0<=n<=N} (-1)^(|N|-|n|) / (N_1! ... N_R!)
 %                         * prod_r nchoosek(N_r,n_r) * Gld_n(|N|)
 %
 % where |N| = sum_r N_r, |n| = sum_r n_r and Gld_n(|N|) is the normalizing
 % constant of the single-class model over the same M stations, with the same
 % load-dependent rates alpha_i(.) = mu(i,.), total population |N|, and
 % aggregated demands rho_i(n) = sum_r n_r * L(i,r). Since
 % nchoosek(N_r,n_r)/N_r! = 1/(n_r! (N_r-n_r)!), the coefficient is evaluated
 % here in the equivalent factorial form, which needs no binomial.
 %
 % The identity holds for arbitrary rate functions alpha_i(.), hence it covers
 % load-independent (alpha=1), infinite-server (alpha(j)=j), multiserver
 % (alpha(j)=min(j,s_i)) and limited load-dependent stations alike. Proof:
 % expand rho_i(n)^k_i multinomially in Gld_n(|N|); the R-fold difference
 % operator annihilates every monomial whose degree in n_r is below N_r, and
 % since the total degree is |N| the only surviving monomial is prod_r n_r^N_r,
 % whose coefficient is the multiclass constant times prod_r N_r!.
 %
 % Cost is prod_r (N_r+1) evaluations of pfqn_gldsingle, i.e. O(M |N|^2)
 % time each. The O(1) space of the theoretical statement is attained only
 % where the single-class constant itself has a closed form (e.g. Gordon's
 % formula in the multiserver case); this implementation uses the standard
 % recursion instead, so the space is that of pfqn_gldsingle.
 %
 % NUMERICS. The sum alternates in sign and its terms are much larger than
 % the result, so it is evaluated as a signed log-sum-exp: this removes the
 % floating-point RANGE problem but not the cancellation. The third output
 % reports the decimal digits lost to cancellation and a warning is raised
 % once the loss exceeds what double precision carries; multiprecision
 % arithmetic is needed beyond that point.
 %
 % @fn pfqn_divdiff_ld(L, N, Z, mu, options)
 % @param L Service demand matrix (MxR).
 % @param N Population vector (1xR).
 % @param Z Think time vector (1xR or DxR), folded in as an infinite server.
 % @param mu Load-dependent rate matrix (Mx sum(N)), alpha_i(j) = mu(i,j).
 % @param options Solver options.
 % @return lG Logarithm of the normalizing constant.
 % @return G Normalizing constant.
 % @return lossDigits Decimal digits lost to cancellation in the alternating sum.
%}
%}
function [lG,G,lossDigits] = pfqn_divdiff_ld(L,N,Z,mu,options)
N = N(:)';
R = length(N);
lossDigits = 0;
if sum(N)<0
    lG = -Inf; G = 0;
    return
end
if sum(N)==0
    lG = 0; G = 1;
    return
end
if nargin<3
    Z = zeros(1,R);
end
if nargin<5 || isempty(options)
    options = SolverNC.defaultOptions;
end
Nt = sum(N);
if nargin<4 || isempty(mu)
    mu = ones(size(L,1),Nt);
end
if size(mu,1) ~= size(L,1)
    line_error(mfilename,'the load-dependent rate matrix must have one row per station of L.');
end
if size(mu,2) < Nt
    line_error(mfilename,'the load-dependent rate matrix must have at least sum(N) columns.');
end
mu = mu(:,1:Nt); % trim so that every rate used downstream is positive
if ~isempty(Z) && sum(Z(:))>0
    L = [L; sum(Z,1)];
    mu(end+1,1:Nt) = 1:Nt; % the delay is an infinite server station
end
if isempty(L) || size(L,1)==0
    lG = -Inf; G = 0;
    return
end
if any(L(:)<0)
    line_error(mfilename,'the demand matrix must be nonnegative.');
end
if any(mu(:)<=0)
    line_error(mfilename,'the load-dependent rates must be strictly positive.');
end
if R==1
    % in a single class the difference operator is the identity, so the
    % recursion is used directly and no cancellation is incurred
    lG = pfqn_gldsingle(L,Nt,mu,options);
    G = exp(lG);
    return
end

% ---- alternating sum over the sub-populations 0 <= n <= N ----
nterms = prod(N+1);
lterm = -Inf(nterms,1);
sterm = zeros(nterms,1);
idx = 0;
n = pprod(N);
while n>=0
    idx = idx + 1;
    if sum(n)>0
        % n = 0 leaves every induced demand at zero, so its single-class
        % constant vanishes at sum(N)>0 and the term is left out of the sum
        lgld = pfqn_gldsingle(L*n(:), Nt, mu, options);
        lterm(idx) = lgld - sum(factln(n)) - sum(factln(N-n));
        sterm(idx) = (-1)^(Nt-sum(n));
    end
    n = pprod(n,N);
end

[lG,sgn,lossDigits] = signedlogsumexp(lterm,sterm);
if sgn==0
    lG = -Inf; G = 0;
    return
elseif sgn<0
    line_warning(mfilename,'The alternating sum returned a negative value, double precision is exhausted by cancellation (%.1f digits lost). Use an exact method such as ''exact'' or ''comomld''.\n',lossDigits);
    lG = NaN; G = NaN;
    return
end
G = exp(lG);
if lossDigits > 15
    line_warning(mfilename,'Cancellation in the alternating sum has consumed about %.1f decimal digits, more than double precision carries. The result is unreliable, multiprecision arithmetic is required.\n',lossDigits);
end
end

% Signed log-sum-exp of S = sum_i sterm(i)*exp(lterm(i)). Returns log|S|,
% sign(S) and the decimal digits lost to cancellation.
function [lS,sgnS,lossDigits] = signedlogsumexp(lterm,sterm)
keep = isfinite(lterm) & sterm~=0;
if ~any(keep)
    lS = -Inf; sgnS = 0; lossDigits = 0;
    return
end
lterm = lterm(keep);
sterm = sterm(keep);
a = max(lterm);
s = sum(sterm .* exp(lterm - a));
sgnS = sign(s);
if s==0
    lS = -Inf; lossDigits = Inf;
    return
end
lS = a + log(abs(s));
% max(exp(lterm-a)) is 1, so -log10|s| is the shortfall of the sum against its
% largest term. Every one of the n terms carries a rounding error of order
% eps*max_term, so the accumulated absolute error is n*eps*max_term and the
% digits actually lost are that shortfall PLUS log10(n). Dropping the count
% understates the loss by log10(n) and lets a wrong answer past the guard.
lossDigits = max(0,log10(numel(sterm)/abs(s)));
end
