%{
%{
 % @file pfqn_gerasimov.m
 % @brief Gerasimov's residue (closed-form) normalizing constant of a
 %        multiclass closed product-form network, generalized to R classes.
%}
%}

function [G,lG] = pfqn_gerasimov(L,N,Z,options)
%{
%{
 % @brief Exact normalizing constant of a closed multiclass product-form
 %        network obtained by ITERATED RESIDUES of its rational generating
 %        function, one class at a time.
 %
 %        A. I. Gerasimov, "On Normalizing Constants in Multiclass Queueing
 %        Networks", Operations Research 43(4):704-711, 1995, evaluates
 %          G(N_1,...,N_R) = (2 pi i)^-R int_G1 ... int_GR
 %                             prod_s z_s^(N_s-1) prod_i (1-sum_s x_is/z_s)^-1
 %        by residues, and gives the resulting CLOSED FORM only for R = 1
 %        (Thm 1-2) and R = 2 (Thm 3 for simple poles, Thm 4 for multiple
 %        ones), stating that "for three or more classes of customers, the
 %        normalizing constants can be found by numerical methods". This
 %        routine implements the residue elimination itself, so the closed
 %        form is produced for ANY R; at R = 2 it reproduces Thm 3/4 term by
 %        term (see the correspondence below).
 %
 % @par The recursion
 %        Write the same object as a coefficient of the u_s = 1/z_s series,
 %          G(N) = [prod_s u_s^(N_s)] exp(sum_s Z_s u_s)
 %                                     prod_i (1 - sum_s x_is u_s)^-1,
 %        and eliminate one class at a time. Every factor is AFFINE in u,
 %        f(u) = c_0 + sum_s c_s u_s, so when u_r is singled out it reads
 %        f = A - B u_r with A affine in the surviving variables. Partial
 %        fractions in u_r give, for a pole of order m_j at u_r = A_j/B_j,
 %          [u_r^n] prod_j (A_j - B_j u_r)^-m_j
 %            = sum_j sum_{k=0}^{m_j-1} (-B_j)^-k C(n+m_j-k-1,n) B_j^n
 %              A_j^-(n+m_j-k) [t^k] prod_{l~=j} (C_jl - B_l t)^-m_l,
 %          C_jl = (A_l B_j - B_l A_j)/B_j,
 %        and C_jl, A_j are again AFFINE. The class-r elimination therefore
 %        maps a sum of products of affine powers into another one with one
 %        variable fewer: the structure is closed under the residue step, and
 %        R-1 steps leave a univariate coefficient extraction.
 %
 % @par Why this is Gerasimov's two-class formula
 %        At R = 2 the first (and only) step has all m_j = 1, so it returns
 %        one term per station i,
 %          x_i2^(N_2+M-1) / prod_{k~=i}(x_i2-x_k2)
 %             * (1-x_i1 u_1)^-(N_2+1) prod_{k~=i} (1 - z_1ik u_1)^-1,
 %          z_1ik = (x_k1 x_i2 - x_i1 x_k2)/(x_i2 - x_k2),
 %        which is exactly the paper's outer factor and its set Omega_i =
 %        {x_i1, z_1ij}. The pole of order N_2+1 at x_i1 is the paper's
 %        tau_ik = taubar_ik + N_2, and coincidences inside Omega_i (the
 %        xi_i < M case of Thm 4) are the multiple poles the step already
 %        handles. The paper's own hypotheses are relaxed in three places:
 %        x_i2 = x_k2 (his z_1ik is then undefined) leaves an affine form
 %        with no constant term, x_i2 = 0 leaves a factor with no pole in
 %        u_2, and identical station rows merge into one factor of doubled
 %        multiplicity. None of the three is a special case here.
 %
 % @par Cost
 %        Let M be the number of stations and order the populations
 %        N_(1) <= ... <= N_(R). The first elimination turns the single input
 %        term into M, and every later one multiplies the count by
 %        C(S+M-1,M-1) + M-1, where S is the total population already
 %        eliminated: a pole of order S+1 has to be differentiated against
 %        the M-1 remaining ones. The innermost extraction then convolves M
 %        series of length N_(1). Hence
 %          R = 1   O(M N), Buzen's own cost;
 %          R = 2   O(M^2 N_(1)^2), INDEPENDENT OF N_(2);
 %          R >= 3  the same times prod_{r=3}^{R} C(N_(r)+M-1, M-1).
 %        The R = 2 line is the reason to reach for this method. A population
 %        removed by residues enters only as a pole ORDER, i.e. through
 %        binomial coefficients, so it costs nothing at all: on a 4-station
 %        two-class model at N = [6, 20000] this returns lG in 0.4 ms where
 %        pfqn_ca needs 1.8 s, to the same 1.3e-16. For R >= 3 the term count
 %        is polynomial in the populations of degree (M-1)(R-2) and
 %        exponential in R, which is why the paper stops at two classes and
 %        why options.maxterms exists.
 %
 % @par Conditioning
 %        The sum is an alternating one over residues, exactly as the paper
 %        writes it, and two decisions keep it usable. Near-coincident poles
 %        are merged under a RELATIVE tolerance, so they are one multiple
 %        pole rather than two nearly cancelling simple ones. And the class
 %        left for the innermost extraction is the one with the SMALLEST
 %        population, because that population is the degree the final,
 %        sign-indefinite series is carried to, whereas the eliminated ones
 %        cancel nothing; basing on N = 100 instead of N = 6 in one
 %        four-station model cost 39 nats of lG. What remains, measured on
 %        372 random models against pfqn_ca: median 2.0e-16, p90 4.6e-15, p99
 %        1.1e-11, worst 1.5e-09. On an ill-conditioned demand matrix pfqn_ca
 %        or pfqn_nc are still the safer routes to the same number.
 %
 % @fn pfqn_gerasimov(L, N, Z, options)
 % @param L Service demand matrix (M x R), L(i,r) = demand of class r at
 %          queueing station i.
 % @param N Population vector (1 x R), nonnegative integers.
 % @param Z Think time vector (1 x R). Default: zeros. A delay contributes
 %          the entire factor exp(sum_s Z_s u_s), which is handled exactly by
 %          convolving its Poisson coefficients into each elimination.
 % @param options Struct with optional fields:
 %          .tol      relative tolerance for declaring two affine forms
 %                    proportional, hence one pole rather than two.
 %                    Default: 1e-12.
 %          .maxterms cap on the number of residue terms carried between
 %                    eliminations. Default: 2e5. Exceeding it is an error,
 %                    not a truncation: a truncated residue sum is not a
 %                    bound or an approximation of G, it is a wrong number.
 % @return G Normalizing constant.
 % @return lG Logarithm of the normalizing constant.
%}
%}

if nargin < 3 || isempty(Z)
    Z = zeros(1,size(L,2));
end
if nargin < 4 || isempty(options)
    options = struct();
end
if ~isfield(options,'tol') || isempty(options.tol)
    options.tol = 1e-12;
end
if ~isfield(options,'maxterms') || isempty(options.maxterms)
    options.maxterms = 2e5;
end

R = size(L,2);
N = N(:)';
Z = Z(:)';
if numel(N) ~= R
    line_error(mfilename,'pfqn_gerasimov requires numel(N) to match the number of columns of L.');
end
if numel(Z) ~= R
    line_error(mfilename,'pfqn_gerasimov requires numel(Z) to match the number of columns of L.');
end
if any(L(:) < 0) || any(Z < 0) || any(N < 0)
    line_error(mfilename,'pfqn_gerasimov requires nonnegative L, N and Z.');
end
if any(N ~= round(N))
    line_error(mfilename,'pfqn_gerasimov requires integer populations.');
end
N = round(N);

% A class with no jobs is eliminated by evaluating the generating function at
% u_r = 0, i.e. by deleting its column outright.
keepr = N > 0;
L = L(:,keepr); N = N(keepr); Z = Z(keepr);
R = numel(N);
if R == 0
    G = 1; lG = 0; return
end

% A station with no demand at all contributes the factor 1.
L = L(any(L > 0, 2),:);

% Per-class scaling. The residue coefficients carry x_ir^(N_r+M-1), which in
% double overflows well before G itself does: at x = 4 and N_r = 400 the factor
% alone is 1e240 while G is finite. Dividing column r by c_r divides G by
% exactly c_r^N_r (substitute u_r -> u_r/c_r in the generating function), so the
% scaling is exact and is undone in the log domain at the end.
cs = max([max(L,[],1); Z], [], 1);
cs(cs <= 0) = 1;
L = L ./ repmat(cs, size(L,1), 1);
Z = Z ./ cs;
lGscale = sum(N .* log(cs));

% Class order, which decides both the cost and the accuracy.
%  - The class left for the innermost extraction sets the CONDITIONING. Its
%    population is the degree the final series is carried to, and the poles of
%    the reduced problem have arbitrary sign, so that series cancels; the
%    populations eliminated by residues enter only as pole ORDERS, through
%    binomial coefficients, and cancel nothing. Basing on N = 100 rather than on
%    N = 6 in one 4-station model cost 39 nats of lG. The SMALLEST population
%    therefore goes to the base.
%  - Eliminating class r leaves a pole of order N_r+1 that every LATER
%    elimination has to differentiate, so the remaining classes are eliminated
%    smallest-first to keep the multiplicities low for as long as possible.
% Eliminations run from index R down to 2, so indices 2..R hold the remaining
% populations in DECREASING order and index 1 holds the smallest.
[~,asc] = sort(N,'ascend');
ord = [asc(1), fliplr(asc(2:end))];
L = L(:,ord); N = N(ord); Z = Z(ord);

M = size(L,1);
terms = cell(1,1);
terms{1} = struct('c',1,'F',[ones(M,1), -L],'m',ones(M,1));

for r = R:-1:2
    terms = geras_step(terms, r, N(r), Z(r), options);
    if isempty(terms)
        G = 0; lG = -Inf; return
    end
end
Gs = geras_base(terms, N(1), Z(1), options);
if Gs <= 0
    G = 0; lG = -Inf; return
end
lG = log(Gs) + lGscale;
G = exp(lG);
end

% ------------------------------------------------------------------------
% One residue elimination: integrate out u_r and return the surviving sum of
% products of affine powers, each form narrowed from r+1 to r columns.
% ------------------------------------------------------------------------
function out = geras_step(terms, r, Nr, Zr, options)
out = cell(1,0);
for it = 1:numel(terms)
    [F,m,c] = geras_merge(terms{it}.F, terms{it}.m, terms{it}.c, options.tol);
    A = F(:,1:r);           % affine part in u_1..u_(r-1), column 1 = constant
    B = -F(:,r+1);          % f_j = A_j - B_j u_r
    scale = max(abs(F),[],2);
    isMono = all(abs(A) <= options.tol*repmat(scale,1,r), 2);
    if any(isMono & abs(B) <= options.tol*scale)
        line_error(mfilename,'pfqn_gerasimov met an identically zero factor, which cannot happen after merging proportional ones.');
    end
    % A factor -B u_r carries no finite pole: it only shifts the exponent.
    shift = 0;
    if any(isMono)
        c = c * prod((-B(isMono)).^(-m(isMono)));
        shift = sum(m(isMono));
        A = A(~isMono,:); B = B(~isMono); m = m(~isMono);
    end
    S = find(B ~= 0);       % factors carrying a pole in u_r
    P = find(B == 0);       % factors free of u_r, carried through unchanged
    Ntot = Nr + shift;
    if Zr > 0
        plist = 0:Ntot;     % Poisson coefficients of exp(Z_r u_r)
    else
        plist = 0;
    end
    for p = plist
        % Poisson weight Z_r^p/p! through logs. The naive ratio overflows for
        % p >~ 171 in double, which is reachable whenever the ELIMINATED class
        % carries think time and a sizeable population; the quotient itself is
        % bounded by exp(Z_r), so only the two halves overflow, not the answer.
        if p > 0
            cz = c * exp(p*log(Zr) - gammaln(p+1));
        else
            cz = c;
        end
        Neff = Ntot - p;
        if isempty(S)
            if Neff == 0
                out{end+1} = struct('c',cz,'F',A,'m',m); %#ok<AGROW>
            end
            continue
        end
        for jj = 1:numel(S)
            j = S(jj);
            oth = S(S ~= j);
            oth = oth(:);
            no = numel(oth);
            if no > 0
                Cjl = (A(oth,:)*B(j) - B(oth)*A(j,:))/B(j);
            else
                Cjl = zeros(0,r);
            end
            for k = 0:(m(j)-1)
                comps = geras_compositions(k, no);
                for ci = 1:size(comps,1)
                    nk = comps(ci,:)';
                    coef = cz * (-B(j))^(-k) * geras_binom(Neff+m(j)-k-1, Neff) * B(j)^Neff;
                    if no > 0
                        for l = 1:no
                            coef = coef * geras_binom(m(oth(l))+nk(l)-1, nk(l)) * B(oth(l))^nk(l);
                        end
                    end
                    if coef == 0
                        continue
                    end
                    out{end+1} = struct('c',coef, ...
                        'F',[A(j,:); Cjl; A(P,:)], ...
                        'm',[Neff+m(j)-k; m(oth)+nk; m(P)]); %#ok<AGROW>
                end
            end
        end
    end
    if numel(out) > options.maxterms
        line_error(mfilename,sprintf('pfqn_gerasimov exceeded options.maxterms (%d) while eliminating class %d; the residue expansion of this model is too large. Use pfqn_ca or pfqn_nc.',options.maxterms,r));
    end
end
end

% ------------------------------------------------------------------------
% Last class: a univariate coefficient extraction. Summing the residues here
% too would repeat the step above, but convolving the series of each factor
% returns the same number without expanding the multiple poles, so it is both
% cheaper and better conditioned.
% ------------------------------------------------------------------------
function G = geras_base(terms, N1, Z1, options)
G = 0;
for it = 1:numel(terms)
    [F,m,c] = geras_merge(terms{it}.F, terms{it}.m, terms{it}.c, options.tol);
    A = F(:,1);
    B = -F(:,2);
    scale = max(abs(F),[],2);
    isMono = abs(A) <= options.tol*scale;
    if any(isMono & abs(B) <= options.tol*scale)
        line_error(mfilename,'pfqn_gerasimov met an identically zero factor at the innermost coefficient extraction.');
    end
    shift = 0;
    if any(isMono)
        c = c * prod((-B(isMono)).^(-m(isMono)));
        shift = sum(m(isMono));
        A = A(~isMono); B = B(~isMono); m = m(~isMono);
    end
    Ntot = N1 + shift;
    s = zeros(1,Ntot+1); s(1) = 1;
    if Z1 > 0
        pois = exp((0:Ntot)*log(Z1) - gammaln((0:Ntot)+1));
        s = geras_conv(s, pois, Ntot);
    end
    for j = 1:numel(A)
        c = c * A(j)^(-m(j));
        if B(j) == 0
            continue
        end
        ratio = B(j)/A(j);
        seq = zeros(1,Ntot+1);
        for n = 0:Ntot
            seq(n+1) = geras_binom(m(j)+n-1, n) * ratio^n;
        end
        s = geras_conv(s, seq, Ntot);
    end
    G = G + c * s(Ntot+1);
end
end

% Truncated convolution of two coefficient sequences, kept to degree n.
function y = geras_conv(a, b, n)
y = conv(a, b);
y = y(1:(n+1));
end

% Merge proportional affine forms: f_k = lambda f_j means one pole of order
% m_j+m_k, not two nearby simple ones, and lambda^-m_k moves into the scalar.
function [F,m,c] = geras_merge(F, m, c, tol)
nf = size(F,1);
if nf <= 1
    return
end
keep = true(nf,1);
for j = 1:nf
    if ~keep(j)
        continue
    end
    [~,pj] = max(abs(F(j,:)));
    for k = (j+1):nf
        if ~keep(k)
            continue
        end
        if abs(F(j,pj)) == 0
            continue
        end
        lambda = F(k,pj)/F(j,pj);
        if norm(F(k,:) - lambda*F(j,:), Inf) <= tol*max(norm(F(k,:),Inf), norm(F(j,:),Inf)) && lambda ~= 0
            c = c * lambda^(-m(k));
            m(j) = m(j) + m(k);
            keep(k) = false;
        end
    end
end
F = F(keep,:); m = m(keep);
end

% All k-tuples of nonnegative integers summing to n.
function C = geras_compositions(n, k)
if k == 0
    if n == 0
        C = zeros(1,0);
    else
        C = zeros(0,0);
    end
    return
end
if k == 1
    C = n;
    return
end
C = zeros(0,k);
for a = 0:n
    sub = geras_compositions(n-a, k-1);
    if ~isempty(sub)
        C = [C; repmat(a,size(sub,1),1), sub]; %#ok<AGROW>
    end
end
end

% Binomial coefficient, exact while the value stays representable.
function b = geras_binom(n, k)
if k < 0 || n < 0 || k > n
    b = 0;
    return
end
k = min(k, n-k);
b = 1;
for i = 1:k
    b = b * (n-k+i)/i;
end
if b < 2^53
    b = round(b);
end
end
