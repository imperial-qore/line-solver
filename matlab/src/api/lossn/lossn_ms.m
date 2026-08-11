%{ @file lossn_ms.m
 %  @brief Exact analysis of loss networks by the Manjunath-Sikdar transform
 %
 %  @author LINE Development Team
%}

%{
 % @brief Exact normalization constant and blocking probabilities of a loss network
 %
 % @details
 % Calls on route r arrive according to a Poisson process of rate nu_r with
 % unit mean holding time, so nu_r is the offered load. A call is admitted
 % only if it leaves every link within capacity,
 %
 %      sum_r A(j,r) n_r <= C(j),   j = 1 ... J,
 %
 % and is lost otherwise. The admissible set is coordinate convex, so by
 % Kelly's truncation theorem the stationary distribution is the truncated
 % product form
 %
 %      p(n) = (1/g(C)) prod_r nu_r^{n_r} / n_r!,   A n <= C,
 %
 % and every metric follows from ratios of the normalization constant
 %
 %      g(C) = sum_{A n <= C} prod_r nu_r^{n_r} / n_r!
 %      E[n_r]  = nu_r g(C - A e_r) / g(C)
 %      Loss_r  = 1 - g(C - A e_r) / g(C)
 %
 % because a class r call is blocked exactly when the state cannot absorb one
 % more unit of its own requirement vector.
 %
 % g(C) is evaluated exactly by the transform technique of Manjunath and
 % Sikdar. Writing the indicator of each constraint as a contour integral
 % turns the sum into a J-fold integral over the unit circle whose integrand
 % factorizes into the z-transforms of the per-route terms,
 %
 %      g(C) = oint ... oint prod_j [ (z_j^{C_j+1}-1) / (z_j^{C_j+1}(z_j-1)) ]
 %                           prod_r Fcal_r(z_1^{A_1r} ... z_J^{A_Jr}) dz_1 ... dz_J.
 %
 % Inside the unit circle the only pole in z_j is at the origin, of order
 % C_j+1, so each integration is a residue and hence a Taylor coefficient.
 % The routine therefore builds the generating function as a multivariate
 % power series truncated at degree C_j in z_j, one shift-and-accumulate
 % convolution per route, and each '<=' constraint is discharged by summing
 % the coefficients of degrees 0 ... C_j along that dimension. Truncation is
 % exact because A is nonnegative: a monomial above degree C_j can never
 % contribute to an extracted coefficient.
 %
 % Contour integrations are interleaved with the product rather than deferred
 % to the end: variable z_j is created when the first route with A(j,r) ~= 0
 % is multiplied in and integrated out immediately after the last one, so the
 % peak memory is the product of (C_j+1) over the simultaneously live links,
 % an induced width, rather than over all J links.
 %
 % Unlike lossn_erlangfp this is exact rather than a reduced-load
 % approximation, and unlike lossn_mci it carries no sampling error, which
 % matters for rare blocking: a loss probability of 1e-4 recovered from a
 % simulated or sampled throughput is dominated by the estimator variance.
 % The cost is prod_j (C_j+1) memory over the live links, so it is exact but
 % not unconditionally cheap; use lossn_mci when that product is prohibitive.
 %
 % A and C must be integer valued, since the residue argument counts whole
 % units of capacity. Each row is divided by the greatest common divisor of
 % its entries and of C(j), which is exact and shrinks the corresponding
 % dimension. Routes appearing in no constraint never block and contribute a
 % factor exp(nu_r) to g(C).
 %
 % @par Syntax:
 % @code
 % [QLen, Loss, lG, niter] = lossn_ms(nu, A, C)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>nu<td>Offered load of route (class) r (1xR vector)
 % <tr><td>A<td>Capacity requirement of link j for route r (JxR nonnegative integers)
 % <tr><td>C<td>Available capacity of link j (Jx1 nonnegative integers)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>QLen<td>Mean carried load E[n_r] for route r (1xR)
 % <tr><td>Loss<td>Blocking probability for route r (1xR)
 % <tr><td>lG<td>Log of the exact normalization constant g(C)
 % <tr><td>niter<td>Number of iterations, always 1 (direct method)
 % </table>
 %
 % @par References:
 % D. Manjunath and B. Sikdar, Integral Expressions for the Numerical
 % Evaluation of Product Form Expressions Over Irregular Multidimensional
 % Integer Spaces.
%}
function [QLen, Loss, lG, niter] = lossn_ms(nu, A, C)
nu = nu(:)';
C  = C(:);
R = numel(nu);
J = numel(C);
niter = 1;

if size(A,1) ~= J || size(A,2) ~= R
    line_error(mfilename, sprintf('A must be %dx%d (J x R).', J, R));
end
if any(A(:) < 0) || any(abs(A(:) - round(A(:))) > 1e-9)
    line_error(mfilename, 'A must contain nonnegative integers; use lossn_mci or lossn_erlangfp otherwise.');
end
if any(C < 0) || any(abs(C - round(C)) > 1e-9)
    line_error(mfilename, 'C must contain nonnegative integers; use lossn_mci or lossn_erlangfp otherwise.');
end
if any(nu < 0)
    line_error(mfilename, 'nu must be nonnegative.');
end
A = round(A);
C = round(C);

% Drop rows that constrain nothing, and reduce each remaining row by the gcd
% of its entries together with its right hand side. Both are exact and the
% second shrinks the truncation degree, which is the dominant cost.
keep = any(A ~= 0, 2);
A = A(keep, :);
C = C(keep);
J = numel(C);
for j = 1:J
    g = C(j);
    for r = find(A(j,:) ~= 0)
        g = gcd(g, A(j,r));
    end
    if g > 1
        A(j,:) = A(j,:) / g;
        C(j) = floor(C(j) / g);
    end
end

% Routes absent from every constraint never block and factor out of g(C).
free = true(1, R);
if J > 0
    free = ~any(A ~= 0, 1);
end
lGfree = sum(nu(free));            % log prod_r exp(nu_r) over free routes

QLen = zeros(1, R);
Loss = zeros(1, R);
QLen(free) = nu(free);

if J == 0 || all(free)
    lG = lGfree;
    return
end

% Per route truncation implied by the constraints, and the scaled terms
% f_r(n) = nu_r^n / n!. Each sequence is divided by its largest entry so that
% moderate loads do not overflow; the scale cancels in every ratio below and
% is added back in the log of the normalization constant.
nmax = zeros(1, R);
f = cell(1, R);
logscale = 0;
for r = 1:R
    if free(r)
        nmax(r) = 0;
        f{r} = 1;
        continue
    end
    v = inf;
    for j = find(A(:,r)' > 0)
        v = min(v, floor(C(j) / A(j,r)));
    end
    nmax(r) = v;
    n = 0:v;
    fr = exp(n * log(max(nu(r), realmin)) - gammaln(n + 1));
    if nu(r) == 0
        fr = zeros(1, v + 1);
        fr(1) = 1;
    end
    s = max(fr);
    if s > 0 && isfinite(s)
        fr = fr / s;
        logscale = logscale + log(s);
    end
    f{r} = fr;
end

G = lossn_ms_series(f, A, C, nmax);
if G <= 0
    line_error(mfilename, 'Empty admissible set: no state satisfies A n <= C.');
end
lG = log(G) + logscale + lGfree;

for r = 1:R
    if free(r)
        continue
    end
    Cr = C - A(:,r);
    if any(Cr < 0)
        Loss(r) = 1;               % a single class r call already overflows
        QLen(r) = 0;
        continue
    end
    Gr = lossn_ms_series(f, A, Cr, nmax);
    ratio = Gr / G;
    QLen(r) = nu(r) * ratio;
    Loss(r) = 1 - ratio;
end
end

% ------------------------------------------------------------------------

function G = lossn_ms_series(f, A, C, nmaxFull)
% Coefficient-domain evaluation of the J-fold contour integral. The series is
% held as a truncated multivariate polynomial over the currently live links;
% link j is created at its first route and integrated out after its last.

R = numel(f);
J = numel(C);
Crow = C(:)';

first = zeros(1, J);
last = zeros(1, J);
for j = 1:J
    idx = find(A(j,:) ~= 0);
    first(j) = idx(1);
    last(j) = idx(end);
end

curdim = ones(1, J);
Aser = 1;

for r = 1:R
    for j = find(first == r)
        Aser = lossn_ms_expand(Aser, curdim, j, C(j) + 1);
        curdim(j) = C(j) + 1;
    end

    s = A(:,r)';
    if any(s ~= 0)
        nmax = min(nmaxFull(r), numel(f{r}) - 1);
        for j = find(s > 0)
            nmax = min(nmax, floor(C(j) / s(j)));
        end
        [SUB, stride] = lossn_ms_index(curdim);
        Anew = zeros(size(Aser));
        for n = 0:nmax
            c = f{r}(n + 1);
            if c == 0
                continue
            end
            if n == 0
                Anew = Anew + c * Aser;
            else
                shift = s * n;
                ok = all(SUB + shift <= Crow, 2);
                if ~any(ok)
                    break                  % shifts only grow with n
                end
                tgt = 1 + (SUB(ok,:) + shift) * stride(:);
                Anew(tgt) = Anew(tgt) + c * Aser(ok);
            end
        end
        Aser = Anew;
    else
        Aser = sum(f{r}(1:nmaxFull(r) + 1)) * Aser;
    end

    for j = find(last == r)
        Aser = lossn_ms_extract(Aser, curdim, j);
        curdim(j) = 1;
    end
end

G = Aser;
end

function A = lossn_ms_expand(A, curdim, j, newdim)
% Create link variable z_j, keeping the existing content at degree zero.
pre = prod(curdim(1:j-1));
post = prod(curdim(j+1:end));
T = zeros(pre, newdim, post);
T(:,1,:) = reshape(A, pre, 1, post);
A = T(:);
end

function A = lossn_ms_extract(A, curdim, j)
% Contour integration in z_j for a '<=' constraint: the multiplier
% (z^{C+1}-1)/(z-1) turns the residue into the partial sum of the
% coefficients of degrees 0 ... C_j, which is the sum along that dimension.
pre = prod(curdim(1:j-1));
dj = curdim(j);
post = prod(curdim(j+1:end));
T = reshape(A, pre, dj, post);
A = reshape(sum(T, 2), [], 1);
end

function [SUB, stride] = lossn_ms_index(curdim)
q = numel(curdim);
P = prod(curdim);
SUB = zeros(P, q);
rep = 1;
for k = 1:q
    SUB(:,k) = repmat(kron((0:curdim(k)-1)', ones(rep,1)), P / (rep * curdim(k)), 1);
    rep = rep * curdim(k);
end
stride = cumprod([1, curdim(1:end-1)]);
end
