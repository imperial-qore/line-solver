%{
%{
 % @file pfqn_explicit.m
 % @brief Explicit closed-form normalizing constant of a multiclass closed network.
%}
%}

%{
%{
 % @brief Explicit closed-form normalizing constant of a multiclass closed network.
 %
 % Evaluates the two explicit expressions of Casale, "Accelerating Performance
 % Inference over Closed Systems by Asymptotic Methods", ACM SIGMETRICS 2017,
 % Eqs. (15) and (16). Both instantiate the divided-difference form of
 % Corollary 3.2,
 %
 %    G(N) = sum_{0<=t<=N} (-1)^(|N|-|t|)/(N_1!...N_R!) prod_r nchoosek(N_r,t_r) g_t(|N|)
 %
 % by substituting a closed form for the single-class constant g_t(|N|) at the
 % induced demands theta_k(t) = sum_r t_r L(k,r):
 %
 %  Eq. (15), induced demands PAIRWISE DISTINCT. Gordon's partial-fraction
 %  formula g_t(|N|) = sum_k theta_k^(|N|+K-1) / prod_{i~=k}(theta_k-theta_i),
 %  with the convention 0/0 = 0. It is O(K) per term of the outer sum.
 %
 %  Eq. (16), induced demands REDUNDANT (repeated). With K' distinct induced
 %  demands theta_j of multiplicity m_j, the general partial-fraction
 %  expansion is used instead,
 %
 %    g_t(|N|) = sum_j (-1)^(m_j-1) theta_j^(|N|+K-m_j)
 %               * sum_{r>=0, |r|=m_j-1} (-1)^r_j nchoosek(|N|+r_j,r_j)
 %                 prod_{k~=j} nchoosek(m_k+r_k-1,r_k) theta_k^r_k / (theta_j-theta_k)^(m_k+r_k)
 %
 %  which reduces to Eq. (15) when every m_j is one. It costs
 %  sum_j nchoosek(K'+m_j-2,m_j-1) evaluations per term of the outer sum.
 %
 % The choice between the two is automatic: the induced demands are scanned
 % over every t of the outer sum and Eq. (16) is used as soon as two of them
 % are closer than tol relative to the largest induced demand at that t,
 % Eq. (15) otherwise. tol defaults to machine precision (eps).
 %
 % Only single-server load-independent queues are admissible: infinite servers
 % need the integral form of Corollary 3.4 and load-dependent rates need
 % pfqn_explicit_ld, which keeps this closed form as its inner kernel, or
 % pfqn_divdiff_ld, which generalizes the outer sum over a recursion instead.
 %
 % NUMERICS. Both expressions alternate in sign with terms far larger than the
 % result, so they are evaluated as signed log-sum-exps: this removes the
 % floating-point RANGE problem but not the cancellation, which is what makes
 % multiprecision arithmetic necessary on all but small models. The fourth
 % output reports the decimal digits lost, and a warning is raised once the
 % loss exceeds what double precision carries.
 %
 % SINGLE CLASS. At R=1 the multiclass constant IS the single-class constant
 % at demands L, so the outer sum is skipped: g_t(N) = t^N g_1(N) and
 % sum_t (-1)^(N-t) t^N/(t!(N-t)!) = S(N,N) = 1. Running the difference anyway
 % would add N alternating terms, and their cancellation, to a closed form that
 % carries none of them. What is left is O(K^2) work at any population, which is
 % why the single-class route is the cheap one on large N.
 %
 % @fn pfqn_explicit(L, N, tol, method, maxloss)
 % @param L Service demand matrix (KxR) of single-server load-independent queues.
 % @param N Population vector (1xR).
 % @param tol Relative tolerance declaring two induced demands redundant (default: eps).
 % @param method 'auto' (default), 'distinct' to force Eq. (15), 'repeated' to force Eq. (16).
 % @param maxloss Cancellation budget in decimal digits. Finite values turn the
 %        warnings into a silent REFUSAL (lG=NaN) once the budget is exceeded,
 %        for callers that hold a fallback; default Inf keeps the warnings.
 % @return lG Logarithm of the normalizing constant.
 % @return G Normalizing constant.
 % @return method Expression actually used, 'distinct' (Eq. 15) or 'repeated' (Eq. 16).
 % @return lossDigits Decimal digits lost to cancellation.
%}
%}
function [lG,G,method,lossDigits] = pfqn_explicit(L,N,tol,method,maxloss)
N = N(:)';
R = length(N);
lossDigits = 0;
if nargin<3 || isempty(tol)
    tol = eps; % machine precision, the tolerance is relative to max(theta)
end
if nargin<5 || isempty(maxloss)
    maxloss = Inf; % warn rather than refuse, the caller has nowhere else to go
end
if nargin<4 || isempty(method)
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
[K,~] = size(L);
Nt = sum(N);

% ---- redundancy scan: are the induced demands pairwise distinct at every t? ----
isRedundant = false;
if R==1
    % the induced demands at t are t*L, so both the tie structure and the
    % relative tolerance are those of L itself, at every t at once
    th = sort(L(:,1));
    scale = th(end);
    isRedundant = scale>0 && any(diff(th) <= tol*scale);
else
    t = pprod(N);
    while t>=0
        if sum(t)>0
            th = sort(L*t(:));
            scale = th(end);
            % scale==0 leaves every induced demand at zero, so g_t(|N|)=0 at
            % sum(N)>0 and the term takes no part in the sum
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
            line_error(mfilename,'Eq. (15) requires pairwise distinct induced demands, but two of them agree to within tol. Use ''auto'' or ''repeated''.');
        end
end

if R==1
    % ---- single class: the divided difference is the identity, evaluate g ----
    % The multiclass constant at R=1 is the single-class constant at demands L,
    % so there is nothing for the difference to extract; see the header note.
    if strcmp(method,'distinct')
        [lG,sgn,lossDigits] = explicit_gdistinct(L(:,1),Nt,K);
    else
        [lG,sgn,lossDigits] = explicit_grepeated(L(:,1),Nt,K,tol);
    end
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
                if strcmp(method,'distinct')
                    [lg,sg,dl] = explicit_gdistinct(th,Nt,K);
                else
                    [lg,sg,dl] = explicit_grepeated(th,Nt,K,tol);
                end
                innerLoss = max(innerLoss,dl);
                if sg~=0
                    lterm(idx) = lg - sum(factln(t)) - sum(factln(N-t));
                    sterm(idx) = sg*(-1)^(Nt-sum(t));
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
