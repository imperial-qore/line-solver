%{
%{
 % @file pfqn_explicit_ld.m
 % @brief Explicit closed-form normalizing constant of a multiclass limited load-dependent network.
%}
%}

%{
%{
 % @brief Explicit closed-form normalizing constant of a multiclass limited load-dependent network.
 %
 % Load-dependent counterpart of pfqn_explicit. It evaluates the same
 % divided-difference form of Casale, "Accelerating Performance Inference over
 % Closed Systems by Asymptotic Methods", ACM SIGMETRICS 2017, Corollary 3.2,
 %
 %    G(N) = sum_{0<=t<=N} (-1)^(|N|-|t|)/(N_1!...N_R!) prod_r nchoosek(N_r,t_r) h_t(|N|)
 %
 % but substitutes for the single-class constant h_t(|N|) the LIMITED
 % LOAD-DEPENDENT closed form of Casale, Harrison and Ong, "Facilitating
 % Load-Dependent Queueing Analysis Through Factorization", Perform. Eval. 2021,
 % Theorem 1, Eq. (8),
 %
 %    h_theta(N) = sum_{0<=v<s} g_sigma(N-|v|) prod_k phi_k(v_k)
 %
 %    phi_k(v_k) = theta_k^v_k / prod_{t=1..v_k} alpha_k(t) * (1 - alpha_k(v_k)/alpha_k(s_k))
 %
 % at the induced demands theta_k(t) = sum_r t_r L(k,r). Here alpha_k(.) = mu(k,.)
 % is the load-dependent scaling of station k, s_k the population past which it
 % stays constant, sigma_k = theta_k/alpha_k(s_k) the SCALED demands, and
 % g_sigma the FIXED-RATE single-class constant at those scaled demands, which
 % is exactly what pfqn_explicit evaluates in closed form (Eqs. 15 and 16). The
 % result is therefore explicit throughout, with no recursion over population:
 % pfqn_divdiff_ld is the same outer sum carried over pfqn_gldsingle's
 % O(M|N|^2) recursion instead.
 %
 % Two conventions of Theorem 1 are not those of the equilibrium distribution
 % and are easy to get wrong. alpha_k(0) is taken as ZERO inside the bracket of
 % phi_k, so that phi_k(0) = 1, even though the state probabilities use
 % alpha_k(0) = 1; and g_sigma(n) = 0 for n < 0, which caps the outer sum at
 % |v| <= |N|. With alpha_k(n) = min(n,s_k) the expression collapses to Gordon's
 % multi-server formula, Oper. Res. 38(5), 1990, Eq. (29), but unlike that one it
 % needs neither a multi-server shape nor distinct scaled demands.
 %
 % LIMITED LOAD DEPENDENCE. Theorem 1 holds for any s_k with
 % alpha_k(n) = alpha_k(s_k) for all n >= s_k, and a LARGER s_k is always
 % admissible, so s_k is detected here as the smallest index whose value the
 % tail of mu(k,:) repeats to within tol. A station whose rates never settle
 % (an infinite server, mu(k,n) = n) gets s_k = |N|, which is still exact:
 % populations above |N| do not occur, so redefining alpha_k there changes
 % nothing. It is merely expensive, since the inner sum costs prod_k s_k terms,
 % capped by |v| <= |N|. Think time is not admissible: a delay would have to
 % enter g_sigma, whose closed form covers queues only.
 %
 % NUMERICS. Both sums alternate in sign with terms far larger than the result,
 % so they are evaluated as signed log-sum-exps: this removes the floating-point
 % RANGE problem but not the cancellation. phi_k is sign-definite when alpha_k
 % increases, as a multi-server station does, and changes sign where alpha_k
 % decreases, so a decreasing rate function costs digits in the inner sum too.
 % The fourth output reports the decimal digits lost and a warning is raised
 % once the loss exceeds what double precision carries.
 %
 % SINGLE CLASS. At R=1 the divided difference is the identity, since h_theta(N)
 % is homogeneous of degree N in theta exactly as in the fixed-rate case, so the
 % outer sum is skipped and Theorem 1 is evaluated once at theta = L.
 %
 % @fn pfqn_explicit_ld(L, N, mu, tol, method, maxloss)
 % @param L Service demand matrix (MxR).
 % @param N Population vector (1xR).
 % @param mu Load-dependent rate matrix (Mx sum(N)), alpha_i(j) = mu(i,j); default all ones.
 % @param tol Relative tolerance declaring two scaled demands redundant, and the rate tail constant (default: eps).
 % @param method 'auto' (default), 'distinct' to force Eq. (15), 'repeated' to force Eq. (16).
 % @param maxloss Cancellation budget in decimal digits. Finite values turn the
 %        warnings into a silent REFUSAL (lG=NaN) once the budget is exceeded,
 %        for callers that hold a fallback; default Inf keeps the warnings.
 % @return lG Logarithm of the normalizing constant.
 % @return G Normalizing constant.
 % @return method Expression actually used for g_sigma, 'distinct' (Eq. 15) or 'repeated' (Eq. 16).
 % @return lossDigits Decimal digits lost to cancellation.
%}
%}
function [lG,G,method,lossDigits] = pfqn_explicit_ld(L,N,mu,tol,method,maxloss)
N = N(:)';
R = length(N);
lossDigits = 0;
if nargin<4 || isempty(tol)
    tol = eps; % machine precision, the tolerance is relative to max(sigma)
end
if nargin<6 || isempty(maxloss)
    maxloss = Inf; % warn rather than refuse, the caller has nowhere else to go
end
if nargin<5 || isempty(method)
    method = 'auto';
end
if ~any(strcmp(method,{'auto','distinct','repeated'}))
    line_error(mfilename,'Unrecognized method, use ''auto'', ''distinct'' (Eq. 15) or ''repeated'' (Eq. 16).');
end
if sum(N)<0
    lG = -Inf; G = 0; method = 'distinct';
    return
end
if sum(N)==0
    lG = 0; G = 1; method = 'distinct';
    return
end
if isempty(L)
    lG = -Inf; G = 0; method = 'distinct';
    return
end
if size(L,2) ~= R
    line_error(mfilename,'the demand matrix must have one column per class of N.');
end
if any(L(:)<0)
    line_error(mfilename,'the demand matrix must be nonnegative.');
end
[M,~] = size(L);
Nt = sum(N);
if nargin<3 || isempty(mu)
    mu = ones(M,Nt);
end
if size(mu,1) ~= M
    line_error(mfilename,'the load-dependent rate matrix must have one row per station of L.');
end
if size(mu,2) < Nt
    line_error(mfilename,'the load-dependent rate matrix must have at least sum(N) columns.');
end
mu = mu(:,1:Nt);
if any(mu(:)<=0)
    line_error(mfilename,'the load-dependent rates must be strictly positive.');
end

% ---- s_k: the smallest index whose value the tail of the rate row repeats ----
% Any larger s_k also satisfies alpha_k(n)=alpha_k(s_k) for n>=s_k, so a missed
% tie only adds terms; a false tie would be a wrong answer, hence the strict tol.
s = ones(M,1);
for i=1:M
    s(i) = Nt;
    tail = mu(i,Nt);
    for n=Nt:-1:2
        if abs(mu(i,n-1)-tail) <= tol*max(abs(tail),1)
            s(i) = n-1;
        else
            break
        end
    end
end
alphaS = mu(sub2ind(size(mu),(1:M)',s));

% ---- per-station phi tables, in the log domain, indexed by v_k = 0..s_k-1 ----
lcum = cell(M,1); % sum_{t=1..v} log alpha_k(t)
lbr = cell(M,1);  % log|1 - alpha_k(v)/alpha_k(s_k)|, with alpha_k(0) := 0
sbr = cell(M,1);
for i=1:M
    lcum{i} = [0; cumsum(log(mu(i,1:(s(i)-1))'))];
    br = 1 - [0; mu(i,1:(s(i)-1))']/alphaS(i);
    lbr{i} = -Inf(s(i),1);
    lbr{i}(br~=0) = log(abs(br(br~=0)));
    sbr{i} = sign(br);
end
vcap = min(s-1,Nt); % g_sigma vanishes below zero population, Eq. (8) caps |v|<=|N|

% ---- redundancy scan: are the SCALED induced demands pairwise distinct? ----
% The scan MUST form sigma exactly as hlld does, (L*t)./alphaS and not
% (L./alphaS)*t: the two orderings differ in the last ulp, so an exact tie can
% clear an eps-relative gap under one and not the other, and Eq. (15) would then
% divide by that ulp. Measured on L=[0 1.3;0.9 0.7], mu=min(n,2), N=[2 3]: at
% t=[2 3] both scaled demands are 1.95, the scaled-first ordering reports a
% 4.4e-16 gap and misses the tie, the demand-first ordering reports 2.2e-16 and
% catches it.
isRedundant = false;
if R==1
    % the scaled demands at t are t*sigma, so both the tie structure and the
    % relative tolerance are those of sigma itself, at every t at once
    th = sort(L(:,1)./alphaS);
    scale = th(end);
    isRedundant = scale>0 && any(diff(th) <= tol*scale);
else
    t = pprod(N);
    while t>=0
        if sum(t)>0
            th = sort((L*t(:))./alphaS);
            scale = th(end);
            % scale==0 leaves every scaled demand at zero, so the term takes no
            % part in the sum
            if scale>0 && any(diff(th) <= tol*scale)
                isRedundant = true;
                break
            end
        end
        t = pprod(t,N);
    end
end
switch method
    case 'auto'
        if isRedundant
            method = 'repeated';
        else
            method = 'distinct';
        end
    case 'distinct'
        if isRedundant
            line_error(mfilename,'Eq. (15) requires pairwise distinct scaled demands, but two of them agree to within tol. Use ''auto'' or ''repeated''.');
        end
end

if R==1
    % ---- single class: the divided difference is the identity, evaluate h ----
    [lG,sgn,lossDigits] = hlld(L(:,1),M,Nt,alphaS,vcap,lcum,lbr,sbr,method,tol);
else
    % ---- outer divided-difference sum over 0 <= t <= N ----
    nterms = prod(N+1);
    lterm = -Inf(nterms,1);
    sterm = zeros(nterms,1);
    innerLoss = 0;
    idx = 0;
    t = pprod(N);
    while t>=0
        idx = idx + 1;
        if sum(t)>0
            th = L*t(:);
            if max(th)>0
                [lh,sh,dl] = hlld(th,M,Nt,alphaS,vcap,lcum,lbr,sbr,method,tol);
                innerLoss = max(innerLoss,dl);
                if sh~=0
                    lterm(idx) = lh - sum(factln(t)) - sum(factln(N-t));
                    sterm(idx) = sh*(-1)^(Nt-sum(t));
                end
            end
        end
        t = pprod(t,N);
    end

    [lG,sgn,lossDigits] = explicit_signedlogsumexp(lterm,sterm);
    lossDigits = max(lossDigits,innerLoss);
end
% A caller that named a cancellation budget has a fallback and wants a verdict,
% not a warning: refuse quietly. lossDigits is Inf when the sum vanished
% identically, which is a total loss rather than a legitimate G=0.
if isfinite(maxloss) && (sgn<0 || lossDigits>maxloss)
    lG = NaN; G = NaN;
    return
end
if sgn==0
    lG = -Inf; G = 0;
    return
elseif sgn<0
    line_warning(mfilename,'The explicit expression returned a negative value, double precision is exhausted by cancellation (%.1f digits lost). Multiprecision arithmetic is required.\n',lossDigits);
    lG = NaN; G = NaN;
    return
end
G = exp(lG);
if lossDigits > 15
    line_warning(mfilename,'Cancellation has consumed about %.1f decimal digits, more than double precision carries. The result is unreliable, multiprecision arithmetic is required.\n',lossDigits);
end
end

% Theorem 1 of Casale-Harrison-Ong (2021), Eq. (8): the single-class limited
% load-dependent constant at induced demands th and total population Nt, as the
% finite sum over 0 <= v < s of the fixed-rate constant at the scaled demands
% th/alpha(s), one population level lower for every job held back by v.
function [lh,sh,lossDigits] = hlld(th,M,Nt,alphaS,vcap,lcum,lbr,sbr,method,tol)
sigma = th ./ alphaS;
lth = -Inf(M,1);
lth(th>0) = log(th(th>0));
nterms = prod(vcap+1);
lterm = -Inf(nterms,1);
sterm = zeros(nterms,1);
lossDigits = 0;
idx = 0;
v = pprod(vcap');
while v>=0
    idx = idx + 1;
    nv = sum(v);
    if nv<=Nt
        lval = 0;
        sval = 1;
        for i=1:M
            vi = v(i);
            if vi>0
                if ~isfinite(lth(i))
                    sval = 0; % theta_k = 0 kills every v_k>0, and 0^0=1 keeps v_k=0
                    break
                end
                % kept inside the guard because 0*(-Inf) is NaN, not 0
                lval = lval + vi*lth(i);
            end
            lval = lval - lcum{i}(vi+1) + lbr{i}(vi+1);
            sval = sval * sbr{i}(vi+1);
        end
        if sval~=0 && isfinite(lval)
            if Nt-nv == 0
                % g_sigma(0) = 1 by definition. Reading it off the partial
                % fraction instead would spend digits on an alternating sum
                % whose value is known exactly.
                lg = 0; sg = 1; dl = 0;
            elseif strcmp(method,'distinct')
                [lg,sg,dl] = explicit_gdistinct(sigma,Nt-nv,M);
            else
                [lg,sg,dl] = explicit_grepeated(sigma,Nt-nv,M,tol);
            end
            lossDigits = max(lossDigits,dl);
            if sg~=0
                lterm(idx) = lval + lg;
                sterm(idx) = sval * sg;
            end
        end
    end
    v = pprod(v,vcap');
end
[lh,sh,dl] = explicit_signedlogsumexp(lterm,sterm);
lossDigits = max(lossDigits,dl);
end
