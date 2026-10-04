%{
%{
 % @file pfqn_rgfmc.m
 % @brief Multiclass Recursion by Generating Functions (RGF): the normalizing
 %        constant by iterated residues, with think times.
%}
%}

function [G,lG] = pfqn_rgfmc(L,N,Z,options)
%{
%{
 % @brief Exact normalizing constant of a closed MULTICLASS product-form
 %        network by eliminating one class at a time by residues, finishing in
 %        the single-class convolution of pfqn_rgf.
 %
 %        P. G. Harrison, S. Coury, "On the asymptotic behaviour of closed
 %        multiclass queueing networks", Performance Evaluation 47:131-138,
 %        2002, Thm 1, expresses the generating function of a q-class network
 %        in terms of those of (q-1)-class networks; P. G. Harrison, T. T. Lee,
 %        "A new recursive algorithm for computing generating functions in
 %        closed multi-class queueing networks", IEEE MASCOTS 2004, eqs. (4)-(5),
 %        turns it into the RGF algorithm, bottoming out in a collection of
 %        single-class normalizing constants memoised by load vector (Sec. 3.4).
 %
 % @par Think times are not in either paper
 %        Both write the generating function as the RATIONAL
 %          H_q(M,X;z) = prod_i (1 - rho_i z)^-m_i,
 %        with every node a load-independent single server. An infinite server
 %        multiplies this by the ENTIRE exp(sum_r Z_r z_r), and that breaks the
 %        step Thm 1 rests on: G_n(z') = -sum_i r_i holds only because the
 %        residues of n(z)/d(z) sum to zero when deg d >= deg n + 2 (Bertozzi
 %        and McKenna, SIAM Review 35(2):239-268, 1993, fact (IV), p. 246), and
 %        an exponential numerator does not decay at infinity.
 %
 %        The delay is therefore carried by their own repair, eqs. (3.19)-(3.21)
 %        of the same paper: only the first k_r+1 Taylor coefficients of
 %        exp(Z_r z_r) can reach the coefficient of z_r^k_r, so replacing the
 %        exponential by that polynomial is EXACT, not an approximation, and
 %        leaves a rational integrand the residue calculus handles unchanged.
 %        The price is that the eliminated class's population re-enters the term
 %        count, which is exactly the population-insensitivity Harrison-Lee
 %        Sec. 4 advertises; the class kept for the base case pays nothing.
 %
 % @par Degeneracy
 %        Thm 1 assumes rho_iq ~= rho_lq and its Conclusion leaves the tied case
 %        open. Two affine forms name the SAME pole only when they are
 %        PROPORTIONAL, so fusing proportional forms into one factor of summed
 %        multiplicity disposes of the degeneracy with nothing else changed; a
 %        tie in the eliminated class alone merely leaves a form with no
 %        constant term, which the recursion carries.
 %
 % @par Conditioning
 %        The elimination is exact in exact arithmetic but is an ALTERNATING sum
 %        over residues, so near-coincident loads over an eliminated class
 %        destroy significance. The worst cancellation ratio is tracked and the
 %        routine REFUSES past options.maxcancel rather than returning a
 %        confidently wrong lG. Everything else runs in the log domain, so no
 %        Poisson weight or binomial is ever formed as a naive ratio.
 %
 % @fn pfqn_rgfmc(L, N, Z, options)
 % @param L Service demand matrix (M x R).
 % @param N Population vector (1 x R), nonnegative integers.
 % @param Z Think time vector (1 x R). Default: zeros.
 % @param options Struct with optional fields:
 %          .tol       relative tolerance for calling two affine forms
 %                     proportional, hence one pole. Default: 1e-12.
 %          .maxterms  cap on residue terms carried between eliminations.
 %                     Default: 1e6.
 %          .maxcancel nats of cancellation tolerated before refusing.
 %                     Default: 15.
 % @return G Normalizing constant.
 % @return lG Logarithm of the normalizing constant.
%}
%}
if nargin < 3 || isempty(Z), Z = zeros(1,size(L,2)); end
if nargin < 4 || isempty(options), options = struct(); end
if ~isfield(options,'tol'), options.tol = 1e-12; end
if ~isfield(options,'maxterms'), options.maxterms = 1e6; end
if ~isfield(options,'maxcancel'), options.maxcancel = 15; end

N = N(:)'; Z = Z(:)';
R = size(L,2);
if numel(N) ~= R || numel(Z) ~= R
    line_error(mfilename,'pfqn_rgfmc requires numel(N) and numel(Z) to match the number of columns of L.');
end
if any(L(:) < 0) || any(Z < 0) || any(N < 0)
    line_error(mfilename,'pfqn_rgfmc requires nonnegative L, N and Z.');
end
if any(abs(N-round(N)) > 0)
    line_error(mfilename,'pfqn_rgfmc requires integer populations.');
end
keepr = N > 0;
L = L(:,keepr); N = N(keepr); Z = Z(keepr); R = sum(keepr);
if R == 0, G = 1; lG = 0; return; end
L = L(any(L>0,2),:);
M = size(L,1);
if M == 0
    lG = sum(N.*log(Z) - gammaln(N+1));
    G = exp(lG);
    return
end
if R == 1
    [G,lG] = pfqn_rgf(L(:,1),N(1),Z(1));
    return
end

% Base class = smallest population: that population is the degree the
% sign-indefinite base series is carried to, so it drives the cancellation.
[~,ord] = sort(N,'ascend');
L = L(:,ord); N = N(ord); Z = Z(ord);
cs = max([max(L,[],1); Z], [], 1);
cs(cs <= 0) = 1;
L = L ./ repmat(cs,M,1); Z = Z ./ cs;
lGscale = sum(N .* log(cs));

terms = {struct('lc',0,'sc',1,'F',[ones(M,1), -L],'m',ones(M,1))};
cond = 0;
% F column p+1 carries class p, so eliminating class p means column p+1.
for col = R:-1:2
    [terms,cond] = rgf_step(terms, col, N(col), Z(col), options, cond);
    if isempty(terms), G = 0; lG = -Inf; return; end
end
[lG,sG,cond] = rgf_base(terms, N(1), Z(1), options, cond);
if sG == 0, G = 0; lG = -Inf; return; end
if sG < 0 || cond > options.maxcancel
    line_error(mfilename,sprintf(['pfqn_rgfmc: the residue sum cancelled %.1f nats, past the %.1f ' ...
        'allowed, so lG carries no significant digits. The eliminated classes have ' ...
        'near-coincident loads over the stations. Use method ''ca'' for the exact convolution.'], ...
        cond, options.maxcancel));
end
lG = lG + lGscale;
G = exp(lG);
end

% ------------------------------------------------------------------ helpers

function [ls,sg,cond] = rgf_slogsum(lv,sv,cond)
% Signed log-domain sum, tracking the worst cancellation ratio seen.
k = isfinite(lv) & sv~=0;
if ~any(k), ls = -Inf; sg = 0; return; end
lv = lv(k); sv = sv(k);
mx = max(lv); e = exp(lv-mx); tot = sum(sv.*e);
if tot == 0, ls = -Inf; sg = 0; return; end
ls = mx + log(abs(tot));
if tot > 0, sg = 1; else, sg = -1; end
if numel(lv) > 1
    cond = max(cond, mx + log(sum(e)) - ls);
end
end

function [lc,sc,cond] = rgf_slogconv(lu,su,lv,sv,cond)
% Signed log-domain linear convolution truncated at the common length.
n = numel(lu);
lc = -inf(1,n); sc = zeros(1,n);
for k = 1:n
    [lc(k),sc(k),cond] = rgf_slogsum(lu(1:k)+lv(k:-1:1), su(1:k).*sv(k:-1:1), cond);
end
end

function v = rgf_lbinom(n,r)
% log C(n,r); never a factorial quotient.
v = gammaln(n+1) - gammaln(r+1) - gammaln(n-r+1);
end

function C = rgf_compositions(total,parts)
% Nonnegative integer rows of length PARTS summing to TOTAL.
if parts == 0
    if total == 0, C = zeros(1,0); else, C = zeros(0,0); end
    return
end
if parts == 1, C = total; return; end
C = [];
for first = 0:total
    sub = rgf_compositions(total-first, parts-1);
    C = [C; [repmat(first,size(sub,1),1), sub]]; %#ok<AGROW>
end
end

function [F,m,lc,sc] = rgf_merge(F,m,tol)
% Fuse PROPORTIONAL affine forms into one factor of summed multiplicity: two
% forms name the same pole exactly when they are proportional, and that is the
% degeneracy Harrison-Coury Thm 1 excludes.
T = size(F,1); used = false(T,1);
outF = []; outm = []; lc = 0; sc = 1;
for i = 1:T
    if used(i), continue; end
    used(i) = true;
    Fi = F(i,:); mi = m(i);
    [~,pi] = max(abs(Fi));
    if Fi(pi) == 0
        line_error(mfilename,'pfqn_rgfmc met an identically zero factor.');
    end
    for j = i+1:T
        if used(j), continue; end
        Fj = F(j,:); nj = max(abs(Fj));
        if nj <= 0, continue; end
        r = Fj(pi)/Fi(pi);
        if r == 0, continue; end
        if all(abs(Fj - r*Fi) <= tol*nj)
            lc = lc - m(j)*log(abs(r));
            if r < 0, sc = sc * (-1)^m(j); end
            mi = mi + m(j); used(j) = true;
        end
    end
    outF = [outF; Fi]; outm = [outm; mi]; %#ok<AGROW>
end
F = outF; m = outm;
end

function [out,cond] = rgf_step(terms, col, kr, Zr, options, cond)
% Eliminate the class in column COL of F. Harrison-Coury Thm 1 as a partial
% fraction: poles of z_col at A_j/B_j of order m_j, the residue there naming a
% network with one class fewer.
tol = options.tol;
out = {};
for it = 1:numel(terms)
    t = terms{it};
    [F,m,dlc,dsc] = rgf_merge(t.F, t.m, tol);
    lc = t.lc + dlc; sc = t.sc * dsc;
    A = F(:,1:col); B = -F(:,col+1);
    scale = max(abs(F),[],2); scale(scale==0) = 1;
    isMono = all(abs(A) <= repmat(tol*scale,1,col), 2);
    shift = 0;
    if any(isMono)
        lc = lc - sum(m(isMono).*log(abs(B(isMono))));
        sc = sc * prod(sign(-B(isMono)).^m(isMono));
        shift = sum(m(isMono));
        keep = ~isMono;
        A = A(keep,:); B = B(keep); m = m(keep);
    end
    S = find(B ~= 0); P = find(B == 0);
    Ntot = kr + shift;
    if Zr > 0, plist = 0:Ntot; else, plist = 0; end
    for p = plist
        % Bertozzi-McKenna truncation, in logs: the naive Z^p/p! overflows.
        if p > 0
            lcz = lc + p*log(Zr) - gammaln(p+1);
        else
            lcz = lc;
        end
        n = Ntot - p;
        if isempty(S)
            if n == 0
                out{end+1} = struct('lc',lcz,'sc',sc,'F',A(P,:),'m',m(P)); %#ok<AGROW>
            end
            continue
        end
        for jj = 1:numel(S)
            j = S(jj);
            oth = S(S ~= j); no = numel(oth);
            Bj = B(j); Aj = A(j,:);
            if no > 0
                Cjl = (A(oth,:)*Bj - B(oth)*Aj)/Bj;
            else
                Cjl = zeros(0,col);
            end
            for k = 0:(m(j)-1)
                lbase = lcz - k*log(abs(Bj)) + rgf_lbinom(n+m(j)-k-1, n) + n*log(abs(Bj));
                sbase = sc * sign(-Bj)^k * sign(Bj)^n;
                comps = rgf_compositions(k, no);
                for cc = 1:size(comps,1)
                    jl = comps(cc,:);
                    lt = lbase; st = sbase;
                    for a = 1:no
                        if jl(a) > 0
                            lt = lt + rgf_lbinom(m(oth(a))+jl(a)-1, jl(a)) + jl(a)*log(abs(B(oth(a))));
                            st = st * sign(B(oth(a)))^jl(a);
                        end
                    end
                    newF = [Aj; Cjl; A(P,:)];
                    newm = [n+m(j)-k; m(oth)+jl(:); m(P)];
                    out{end+1} = struct('lc',lt,'sc',st,'F',newF,'m',newm); %#ok<AGROW>
                end
            end
        end
    end
    if numel(out) > options.maxterms
        line_error(mfilename,sprintf(['pfqn_rgfmc exceeded maxterms=%g. The residue term count grows ' ...
            'as C(S+M-1,M-1) per further elimination; use method ''ca''.'], options.maxterms));
    end
end
end

function [lg,sg,cond] = rgf_base(terms, k1, Z1, options, cond)
% Single-class base case, memoised on (loads, multiplicities) per Harrison-Lee
% Sec. 3.4, with loads of either sign: an eliminated class leaves pole
% differences that are not sign-definite.
tol = options.tol;
lv = []; sv = [];
keys = {}; vals = {};
for it = 1:numel(terms)
    t = terms{it};
    [F,m,dlc,dsc] = rgf_merge(t.F, t.m, tol);
    lc = t.lc + dlc; sc = t.sc * dsc;
    shift = 0; loads = []; mults = [];
    for a = 1:size(F,1)
        a0 = F(a,1); a1 = F(a,2); ma = m(a);
        sca = max(abs(a0),abs(a1));
        if sca == 0
            line_error(mfilename,'pfqn_rgfmc met an identically zero factor.');
        end
        if abs(a0) <= tol*sca
            lc = lc - ma*log(abs(a1));
            sc = sc * sign(a1)^ma;
            shift = shift + ma;
        else
            lc = lc - ma*log(abs(a0));
            sc = sc * sign(a0)^ma;
            if abs(a1) > tol*sca
                loads(end+1) = -a1/a0; mults(end+1) = ma; %#ok<AGROW>
            end
        end
    end
    Ntot = k1 + shift;
    key = sprintf('%d|%.17g|%s', Ntot, Z1, mat2str([loads(:).'; mults(:).'],17));
    hit = find(strcmp(keys,key),1);
    if isempty(hit)
        [lg1,sg1,cond] = rgf_base_kernel(loads, mults, Ntot, Z1, cond);
        keys{end+1} = key; vals{end+1} = [lg1 sg1]; %#ok<AGROW>
    else
        lg1 = vals{hit}(1); sg1 = vals{hit}(2);
    end
    if sg1 ~= 0
        lv(end+1) = lc + lg1; sv(end+1) = sc * sg1; %#ok<AGROW>
    end
end
[lg,sg,cond] = rgf_slogsum(lv, sv, cond);
end

function [lg,sg,cond] = rgf_base_kernel(loads, mults, N, Z, cond)
% [z^N] exp(Z z) prod_t (1 - p_t z)^-m_t, signed and in the log domain: this is
% Coury-Harrison (1997) Property 1 with the loads allowed to be negative.
kk = 0:N;
lgv = -inf(1,N+1); sgv = zeros(1,N+1);
lgv(1) = 0; sgv(1) = 1;
if Z > 0
    [lgv,sgv,cond] = rgf_slogconv(lgv, sgv, kk*log(Z) - gammaln(kk+1), ones(1,N+1), cond);
end
for a = 1:numel(loads)
    p = loads(a); mm = mults(a);
    if p == 0, continue; end
    if mm == 1
        lr = kk*log(abs(p));
    else
        lr = gammaln(kk+mm) - gammaln(kk+1) - gammaln(mm) + kk*log(abs(p));
    end
    if p > 0, sr = ones(1,N+1); else, sr = (-1).^kk; end
    [lgv,sgv,cond] = rgf_slogconv(lgv, sgv, lr, sr, cond);
end
lg = lgv(N+1); sg = sgv(N+1);
end
