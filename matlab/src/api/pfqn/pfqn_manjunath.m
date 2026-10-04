%{ @file pfqn_manjunath.m
 %  @brief Exact normalizing constant of a product-form network over a state
 %         space cut by linear integer constraints (Manjunath-Sikdar)
 %
 %  @author LINE Development Team
%}

%{
 % @brief Exact normalizing constant of a closed multiclass product-form network
 %        whose state space carries arbitrary linear integer constraints
 %
 % @details
 % This is the queueing-network half of the transform technique of Manjunath and
 % Sikdar, of which lossn_manjunath is the loss-network half.
 % The two solve the same problem -- sum a product form over an irregular
 % integer state space -- but from opposite ends of the paper: lossn_manjunath
 % implements Section 2.2, a set of '<=' constraints over Poisson terms
 % nu^n/n!, while this routine implements Section 3 together with Section 5.3,
 % a MIXED set of '=', '<=' and '>' constraints over the BCMP terms, where the
 % population constraint of a closed network is itself one of the equalities.
 %
 % @par The model
 % M queueing stations (FCFS, PS or LCFS, rows of L) and Mz delay stations
 % (rows of Z) serve R closed classes with populations N. Writing n_ir for the
 % number of class r jobs at station i and n_i = sum_r n_ir, the BCMP product
 % form of Baskett-Chandy-Muntz-Palacios is
 %
 %   p(n) = (1/G) prod_{i queueing} n_i! prod_r L_ir^{n_ir}/n_ir!
 %                prod_{i delay}          prod_r Z_ir^{n_ir}/n_ir!
 %
 % and G is its sum over the admissible set. Every state obeys the R population
 % equalities sum_i n_ir = N_r; on top of those the caller may impose any number
 % of further rows
 %
 %   sum_{i,r} A(j, i + S(r-1)) n_ir  {=, <=, >}  b(j),   S = M + Mz,
 %
 % i.e. A acts on n(:), the (M+Mz)-by-R occupancy matrix read column by column
 % with the queueing stations first. With no extra rows the routine returns
 % exactly the normalizing constant of pfqn_ca, which is the parity oracle used
 % by the tests; with extra rows it answers a question no other routine in the
 % pfqn family can, because the convolution and MVA recursions are built around
 % the population constraint alone and have nowhere to carry a second one.
 %
 % @par Why the generating function is a product, and why the n_i! disappears
 % Marking class r by z_r and constraint row j by y_j, and abbreviating the
 % monomial that one class r job at station i contributes as
 %
 %   u_ir = z_r prod_j y_j^{A(j, i + S(r-1))},
 %
 % the sum over the occupancies of a single QUEUEING station is, by the
 % multinomial theorem,
 %
 %   sum_{n_i.} n_i! prod_r (L_ir u_ir)^{n_ir}/n_ir!
 %     = sum_k (sum_r L_ir u_ir)^k = 1 / (1 - sum_r L_ir u_ir),
 %
 % so the n_i! that couples the classes at a queueing station is exactly what
 % turns the station's factor from an exponential into a geometric one. The
 % paper reaches the same place through the Euler integral n! = int_0^inf
 % e^{-t} t^n dt (Eqns 16-18), which is that geometric series evaluated; the
 % closed form is used here because there is then no quadrature to discretize.
 % A DELAY station has no n_i! and keeps its exponential,
 % prod_r exp(Z_ir u_ir). Hence
 %
 %   F(z,y) = prod_{i=1}^{M} 1/(1 - sum_r L_ir u_ir)
 %            prod_{k=1}^{Mz} prod_r exp(Z_kr u_kr)
 %
 % and G is read off F as a coefficient: degree exactly N_r in z_r for every
 % class, and for row j the degree dictated by its sense -- exactly b_j for '=',
 % the sum of degrees 0 ... b_j for '<=' (the multiplier (y^{b+1}-1)/(y-1) of
 % the paper's Eqn 5, whose residue is that partial sum), and the complement of
 % the latter for '>' (Eqn 6).
 %
 % @par Why this is a coefficient computation and not a quadrature
 % The contour integrals of Eqn 9 all have their only pole at the origin, of
 % order one more than the right-hand side, so each is a residue and hence a
 % Taylor coefficient. The routine therefore never integrates: it carries F as a
 % multivariate power series truncated at degree N_r in z_r and b_j in y_j.
 % Truncation is exact because A is nonnegative -- no monomial above the cut can
 % be brought back down by a later factor.
 %
 % Each queueing station is applied by SOLVING (1 - sum_r L_ir u_ir) x = ser
 % rather than by expanding the geometric series, which is what keeps the cost
 % at one pass. Every monomial of the operator carries z_r to a strictly higher
 % power, so sweeping the lattice in increasing total class degree lets each
 % coefficient read only coefficients already final: a Gauss-Seidel sweep whose
 % result is the exact solve, not an iterate. A delay station has no such
 % recurrence and is convolved with exp term by term, which is where the extra
 % factor of the population in the cost comes from.
 %
 % @par The elimination order is the memory bound
 % Variable y_j is created when the first station its row touches is multiplied
 % in and is discharged immediately after the last one, so peak memory is
 % prod_r (N_r+1) times the product of (b_j+1) over the SIMULTANEOUSLY LIVE
 % rows, not over all rows. A row that constrains one station therefore costs
 % essentially nothing. The class dimensions are live throughout, so
 % prod_r (N_r+1) is a floor on the cost -- the same lattice pfqn_ca walks.
 %
 % @par Scope
 % Load-dependent and multiserver stations are NOT covered: their per-station
 % term is not geometric, and while the paper admits an arbitrary f_i(n_i) in
 % the single-class case (Section 2), the multiclass n_i! coupling used above
 % then breaks. Use pfqn_gld or pfqn_conwayms for those. A and b must be integer
 % valued and A nonnegative, since the residue argument counts whole units; a
 % fractional entry is refused rather than rounded.
 %
 % @par The per-class decomposition, and where the blocked jobs sit
 % Asking for the fourth output STATS returns the whole solution of the
 % constrained network, not just its normalizing constant. It requires the ONE
 % configuration in which the truncated product form is the EXACT stationary
 % law: a single queueing station inside the region and a SINGLE DELAY STATION
 % OUTSIDE IT. Anything else is refused by name rather than answered wrongly.
 %
 % Why that configuration is not a convenience: with one queueing station the
 % state is the queue occupancy alone (the delay holds the complement) and every
 % transition moves one job of one class by one unit, so the chain is a
 % multidimensional birth-death process. That process is reversible, and Kelly's
 % truncation theorem then applies verbatim -- restricting it to the
 % coordinate-convex set A n <= b and renormalizing gives exactly the truncated
 % product form. Add a second queueing station and the delay -> q1 -> q2 -> delay
 % cycle destroys reversibility; truncation no longer preserves the product form,
 % measured at 131% relative error on the stationary law of a 2-class, N = [2 2]
 % instance. G and lG remain correct as a sum over the admissible set in every
 % configuration; only the metrics are withheld.
 %
 % Everything follows from two ratios of normalizing constants, both taken in the
 % log domain so the internal rescaling cancels without being reconstructed:
 %
 %   X_r = G(N - e_r ; b - A(:,qcol_r)) / G(N ; b)
 %   P(n_qr = k) = G(N ; b, with the added row n_qr = k) / G(N ; b)
 %
 % The first is the loss network's g(C - A e_r) in another guise: removing one
 % class r job from the queue leaves a state whose admission rule is shifted by
 % that job's own requirement column. The second is what an '=' row is for.
 %
 % A refused admission is a DELETED transition, so a blocked job never leaves the
 % delay, and since the think time is exponential a held job is indistinguishable
 % from one still thinking. The delay population carries both, and Little's law
 % separates them:
 %
 %   delay_r   = N_r - Q_r          (everything not at the queue)
 %   think_r   = X_r Z_r            (genuinely thinking)
 %   blocked_r = delay_r - think_r  (held at the delay by the constraint)
 %
 % This is NOT the WAITQ rule of SolverSSA/SolverCTMC/JMT, which moves a refused
 % job out of the delay into a per-region FIFO counted at no station. That is a
 % different chain, and on the reference instance below its class throughputs
 % differ by 15%.
 %
 % @par Syntax:
 % @code
 % [G, lG] = pfqn_manjunath(L, N)
 % [G, lG] = pfqn_manjunath(L, N, Z)
 % [G, lG, peak] = pfqn_manjunath(L, N, Z, A, b, sense)
 % [G, lG, peak, stats] = pfqn_manjunath(L, N, Z, A, b, sense)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>L<td>Service demand of class r at queueing station i (MxR)
 % <tr><td>N<td>Population of class r (1xR nonnegative integers)
 % <tr><td>Z<td>Think time of class r at delay station k (MzxR), default zeros
 % <tr><td>A<td>Extra constraint coefficients on n(:) (Jx((M+Mz)*R), nonnegative integers), default empty
 % <tr><td>b<td>Extra constraint right-hand sides (Jx1 integers), default empty
 % <tr><td>sense<td>Row senses, char vector of 'E' (=), 'L' (<=), 'G' (>), default all 'L'
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>G<td>Normalizing constant
 % <tr><td>lG<td>Logarithm of the normalizing constant
 % <tr><td>peak<td>Peak number of live series coefficients, the realised cost
 % <tr><td>stats<td>Per-class decomposition (1xR fields Q, X, U, think, blocked, delay). Needs one queueing station and one delay station outside the region
 % </table>
 %
 % @par stats fields:
 % <table>
 % <tr><th>Field<th>Description
 % <tr><td>Q<td>Mean class r jobs at the queueing station
 % <tr><td>X<td>Class r cycle throughput
 % <tr><td>U<td>Class r utilization of the queueing station (X_r times its demand)
 % <tr><td>think<td>Class r jobs genuinely thinking, X_r Z_r
 % <tr><td>blocked<td>Class r jobs held at the delay by the constraint
 % <tr><td>delay<td>Class r jobs at the delay, think + blocked
 % </table>
 %
 % @par References:
 % D. Manjunath and B. Sikdar, Integral Expressions for the Numerical
 % Evaluation of Product Form Expressions Over Irregular Multidimensional
 % Integer Spaces. Sections 3 and 5.3.
%}
function [G, lG, peak, stats] = pfqn_manjunath(L, N, Z, A, b, sense)
N = N(:)';
R = numel(N);
if nargin < 1 || isempty(L)
    L = zeros(0, R);
end
if size(L,2) ~= R
    line_error(mfilename, sprintf('L must have %d columns, one per class.', R));
end
if nargin < 3 || isempty(Z)
    Z = zeros(0, R);
end
if size(Z,2) ~= R
    line_error(mfilename, sprintf('Z must have %d columns, one per class.', R));
end
M = size(L,1);
Mz = size(Z,1);
S = M + Mz;
if nargin < 4
    A = [];
end
if nargin < 5
    b = [];
end
if nargin < 6
    sense = [];
end
b = b(:);
J = numel(b);
if isempty(A)
    A = zeros(J, S*R);
end
if J == 0
    A = zeros(0, S*R);
end
if size(A,1) ~= J || size(A,2) ~= S*R
    line_error(mfilename, sprintf('A must be %dx%d (J x (M+Mz)*R), acting on n(:).', J, S*R));
end
if isempty(sense)
    sense = repmat('L', 1, J);
end
sense = upper(char(sense(:)'));
if numel(sense) ~= J
    line_error(mfilename, 'sense must have one character per row of b.');
end
if any(~ismember(sense, 'ELG'))
    line_error(mfilename, 'sense must contain only ''E'' (=), ''L'' (<=) or ''G'' (>).');
end
if any(N < 0)
    G = 0; lG = -Inf; peak = 0;
    if nargout >= 4, stats = pfqn_manjunath_emptystats(R); end
    return
end
if any(abs(N - round(N)) > 1e-9)
    line_error(mfilename, 'N must contain nonnegative integers.');
end
N = round(N);
if any(A(:) < 0) || any(abs(A(:) - round(A(:))) > 1e-9)
    line_error(mfilename, 'A must contain nonnegative integers; the residue argument counts whole units.');
end
if any(abs(b - round(b)) > 1e-9)
    line_error(mfilename, 'b must contain integers; the residue argument counts whole units.');
end
A = round(A);
b = round(b);

% A row of zeros constrains nothing, so it is decided here rather than carried
% as a one-coefficient dimension: 0 = b, 0 <= b and 0 > b are each settled by
% the sign of b alone. Same for a negative right-hand side, which no
% nonnegative combination can meet ('E', 'L') or can fail to beat ('G').
keep = true(J,1);
for j = 1:J
    trivial = ~any(A(j,:) ~= 0);
    switch sense(j)
        case 'E'
            if b(j) < 0 || (trivial && b(j) ~= 0)
                G = 0; lG = -Inf; peak = 0;
                if nargout >= 4, stats = pfqn_manjunath_emptystats(R); end
                return
            end
            if trivial, keep(j) = false; end
        case 'L'
            if b(j) < 0
                G = 0; lG = -Inf; peak = 0;
                if nargout >= 4, stats = pfqn_manjunath_emptystats(R); end
                return
            end
            if trivial, keep(j) = false; end
        case 'G'
            if b(j) < 0
                keep(j) = false;            % 0 > negative always holds
            elseif trivial
                G = 0; lG = -Inf; peak = 0; % 0 > nonnegative never holds
                if nargout >= 4, stats = pfqn_manjunath_emptystats(R); end
                return
            end
    end
end
A = A(keep,:);
b = b(keep);
sense = sense(keep);
J = numel(b);

if S == 0
    % No station: the only state is empty, admissible when every class is too.
    if all(N == 0)
        G = 1; lG = 0; peak = 1;
    else
        G = 0; lG = -Inf; peak = 0;
    end
    if nargout >= 4, stats = pfqn_manjunath_emptystats(R); end
    return
end

% Every monomial that survives the extraction has total degree sum(N) in the
% demands, so a common rescaling of L and Z moves lG by a known amount and
% nothing else. The exponent is chosen as in pfqn_ca, from the largest single
% term the network can produce, so that the series is centred near 1.
Nt = sum(N);
if Nt == 0
    cscale = 1;
else
    lGest = -Inf;
    for i = 1:M
        t = 0; ok = true;
        for r = 1:R
            if N(r) > 0
                if L(i,r) > 0
                    t = t + N(r)*log(L(i,r));
                else
                    ok = false; break
                end
            end
        end
        if ok, lGest = max(lGest, t); end
    end
    if Mz > 0
        t = 0; ok = true;
        for r = 1:R
            if N(r) > 0
                if sum(Z(:,r)) > 0
                    t = t + N(r)*log(sum(Z(:,r))) - factln(N(r));
                else
                    ok = false; break
                end
            end
        end
        if ok, lGest = max(lGest, t); end
    end
    if ~isfinite(lGest)
        cscale = 1;
    else
        cscale = pow2(round(lGest/(Nt*log(2))));
    end
end
Ls = L / cscale;
Zs = Z / cscale;

% A '>' row is the complement of a '<=' row at the same right-hand side, which
% is how the paper discharges it (Eqn 6). With several such rows the product of
% the complements expands by inclusion-exclusion, so the series is evaluated
% once per subset of them, with the subset's rows re-entered as '<=' and the
% rest dropped. Exact, and the only place the cost is exponential -- in the
% number of '>' rows, which is normally zero.
gt = find(sense == 'G');
K = numel(gt);
base = setdiff(1:J, gt);
Gs = 0;
peak = 0;
for mask = 0:(2^K - 1)
    sel = base;
    for t = 1:K
        if bitand(mask, bitshift(1, t-1))
            sel(end+1) = gt(t); %#ok<AGROW>
        end
    end
    sel = sort(sel);
    subsense = sense(sel);
    subsense(subsense == 'G') = 'L';
    [g, pk] = pfqn_manjunath_series(Ls, Zs, N, A(sel,:), b(sel), subsense);
    sgn = 1 - 2*mod(sum(bitget(mask, 1:max(K,1))), 2);
    Gs = Gs + sgn*g;
    peak = max(peak, pk);
end

if Gs <= 0
    % Either the admissible set is empty or the '>' complements cancelled it.
    G = 0;
    lG = -Inf;
    if nargout >= 4
        stats = pfqn_manjunath_emptystats(R);
    end
    return
end
lG = log(Gs) + Nt*log(cscale);
G = exp(lG);

if nargout >= 4
    stats = pfqn_manjunath_stats(L, N, Z, A, b, sense, lG, M, Mz, S, R);
end
end

% ------------------------------------------------------------------------

function stats = pfqn_manjunath_emptystats(R)
stats = struct('Q', zeros(1,R), 'X', zeros(1,R), 'U', zeros(1,R), ...
    'think', zeros(1,R), 'blocked', zeros(1,R), 'delay', zeros(1,R));
end

function stats = pfqn_manjunath_stats(L, N, Z, A, b, sense, lG, M, Mz, S, R)
% Per-class decomposition of the constrained closed network, for the ONE
% configuration in which the truncated product form is the exact stationary law:
% a single queueing station inside the finite capacity region and a single delay
% station outside it.
%
% WHY THE CONFIGURATION IS NOT A CONVENIENCE. With one queueing station the
% state is the queue occupancy alone (the delay holds the complement), every
% transition moves ONE job of ONE class by one unit, and the chain is a
% multidimensional birth-death process. That process is reversible, so Kelly's
% truncation theorem applies verbatim: restricting it to the coordinate-convex
% set A n <= b and renormalizing gives exactly the truncated product form. Add a
% second queueing station and the delay -> q1 -> q2 -> delay cycle destroys
% reversibility; truncation then does NOT preserve the product form, measured at
% 131% relative error on the stationary law of a 2-class, N=[2 2] instance. The
% metrics are therefore refused rather than reported wrong.
%
% WHERE THE BLOCKED JOBS SIT. Nowhere special: a refused admission is a DELETED
% transition, so the job never leaves the delay, and because the think time is
% exponential a held job is indistinguishable from one still thinking. The delay
% population therefore carries both, and Little's law splits them: X_r*Z_r are
% genuinely thinking and the remainder is held. This is NOT the WAITQ rule of
% SolverSSA/SolverCTMC/JMT, which moves a refused job out of the delay into a
% per-region FIFO counted at no station, and which is a different chain.

if Mz ~= 1
    line_error(mfilename, sprintf(['the per-class decomposition needs exactly one delay ' ...
        'station, got %d. Pass Z as a 1xR row of think times.'], Mz));
end
if M ~= 1
    line_error(mfilename, sprintf(['the per-class decomposition needs exactly one queueing ' ...
        'station, got %d. With two or more the delay->q1->q2->delay cycle makes the chain ' ...
        'irreversible, Kelly truncation no longer holds, and the truncated product form is ' ...
        'not the stationary law (measured at 131%% error). G and lG are still returned and ' ...
        'still correct as a sum over the admissible set.'], M));
end
% The delay must sit OUTSIDE the region: its columns are S*r, r = 1..R.
for r = 1:R
    if any(A(:, S*r) ~= 0)
        line_error(mfilename, sprintf(['constraint row(s) reference the delay station in ' ...
            'class %d (column %d). The delay must lie OUTSIDE the finite capacity region, ' ...
            'because the decomposition charges every held job to it.'], r, S*r));
    end
end

% Ratios of normalizing constants are taken in the LOG domain, so the internal
% power-of-two rescaling cancels without ever being reconstructed.
Q = zeros(1,R); X = zeros(1,R);
eqrow = zeros(1, S*R);
for r = 1:R
    qcol = 1 + S*(r-1);                 % column of (queueing station, class r)

    % Throughput. One class r job removed from the queue leaves a state of
    % population N - e_r whose admission rule is shifted by that job's own
    % requirement column, exactly as the loss network's g(C - A e_r):
    %   X_r = G(N - e_r ; b - A(:,qcol)) / G(N ; b).
    Nr = N; Nr(r) = Nr(r) - 1;
    if Nr(r) >= 0
        [~, lGr] = pfqn_manjunath(L, Nr, Z, A, b - A(:,qcol), sense);
        X(r) = exp(lGr - lG);
    end

    % Mean queue length, from the marginal law. An '=' row is discharged by
    % picking a single coefficient, so each call returns the mass of exactly
    % that occupancy.
    for k = 1:N(r)
        row = eqrow; row(qcol) = 1;
        [~, lGk] = pfqn_manjunath(L, N, Z, [A; row], [b; k], [sense 'E']);
        Q(r) = Q(r) + k*exp(lGk - lG);
    end
end

U = X .* L(1,:);            % single server, one visit: U_r = X_r * demand
think = X .* Z(1,:);        % Little's law at the delay, genuinely thinking
delay = N - Q;              % everything not at the queue is at the delay
blocked = delay - think;    % the remainder is held there by the constraint
stats = struct('Q', Q, 'X', X, 'U', U, 'think', think, 'blocked', blocked, ...
    'delay', delay);
end

% ------------------------------------------------------------------------

function [G, peak] = pfqn_manjunath_series(L, Z, N, A, b, sense)
% Coefficient-domain evaluation of the multiple contour integral of Eqn 9 for a
% set of '=' and '<=' rows. The series is held as a truncated multivariate
% polynomial whose first R axes are the class markers z_r, of extent N_r+1
% throughout, and whose remaining J axes are the row markers y_j, of extent 1
% while row j is not live and b_j+1 while it is.

M = size(L,1);
Mz = size(Z,1);
S = M + Mz;
R = numel(N);
J = numel(b);

% Row j is created at the first station it touches and discharged after the
% last, so only an induced width of rows is ever live. A row that reached here
% touches at least one station, the trivial ones having been decided already.
first = zeros(1,J);
last = zeros(1,J);
for j = 1:J
    touched = false(1,S);
    for i = 1:S
        for r = 1:R
            if A(j, i + S*(r-1)) ~= 0
                touched(i) = true;
            end
        end
    end
    idx = find(touched);
    first(j) = idx(1);
    last(j) = idx(end);
end

dims = [N + 1, ones(1,J)];
ser = zeros(prod(dims), 1);
ser(1) = 1;
peak = numel(ser);
[SUB, stride, lev] = pfqn_manjunath_index(dims, R);

for i = 1:S
    for j = find(first == i)
        newdim = b(j) + 1;
        ser = pfqn_manjunath_expand(ser, dims, R+j, newdim);
        dims(R+j) = newdim;
        peak = max(peak, numel(ser));
        [SUB, stride, lev] = pfqn_manjunath_index(dims, R);
    end

    if i <= M
        % Queueing station: solve (1 - sum_r L_ir u_ir) x = ser in place. Every
        % monomial of the operator raises the total class degree by one, so a
        % sweep in increasing total class degree reads only final coefficients
        % and the sweep IS the solve.
        coef = L(i,:);
        [off, ok] = pfqn_manjunath_shifts(A, SUB, stride, i, S, R, J, dims);
        for ell = 1:max(lev)
            at = find(lev == ell);
            for r = 1:R
                if coef(r) == 0 || isempty(off{r})
                    continue
                end
                sel = at(ok{r}(at));
                if isempty(sel)
                    continue
                end
                ser(sel) = ser(sel) + coef(r)*ser(sel - off{r});
            end
        end
    else
        % Delay station: no n_i! coupling, so the factor is a product of
        % exponentials, one per class, each convolved in term by term. There is
        % no first-order recurrence to exploit here, which is why the delay
        % costs a factor of the population that the queueing station does not.
        coef = Z(i-M,:);
        [off, ok] = pfqn_manjunath_shifts(A, SUB, stride, i, S, R, J, dims);
        for r = 1:R
            if coef(r) == 0 || isempty(off{r}) || N(r) == 0
                continue
            end
            nxt = ser;
            term = ser;
            for n = 1:N(r)
                shifted = zeros(size(term));
                sel = find(ok{r});
                shifted(sel) = term(sel - off{r});
                term = (coef(r)/n)*shifted;
                if ~any(term)
                    break
                end
                nxt = nxt + term;
            end
            ser = nxt;
        end
    end

    for j = find(last == i)
        ser = pfqn_manjunath_reduce(ser, dims, R+j, sense(j), b(j));
        dims(R+j) = 1;
        [SUB, stride, lev] = pfqn_manjunath_index(dims, R);
    end
end

if numel(ser) ~= prod(N+1)
    line_error(mfilename, 'a constraint row was never discharged; the elimination order is inconsistent.');
end
% The closed network's own equalities: degree exactly N_r in every class.
G = ser(1 + sum(N .* cumprod([1, N(1:end-1)+1])));
end

% ------------------------------------------------------------------------

function [off, ok] = pfqn_manjunath_shifts(A, SUB, stride, i, S, R, J, dims)
% Flat offset and in-range mask of the monomial one class r job at station i
% contributes: z_r gains one degree and y_j gains A(j, i + S(r-1)). A row that
% is not live at this station has a zero entry here by construction of
% first/last, so a dead axis is never shifted.
off = cell(1,R);
ok = cell(1,R);
for r = 1:R
    delta = zeros(1, R+J);
    delta(r) = 1;
    for j = 1:J
        delta(R+j) = A(j, i + S*(r-1));
    end
    if any(delta > dims - 1)
        off{r} = [];                 % a single job already breaks the cut
        ok{r} = false(size(SUB,1),1);
        continue
    end
    off{r} = delta*stride(:);
    ok{r} = all(SUB >= delta, 2);
end
end

function [SUB, stride, lev] = pfqn_manjunath_index(dims, R)
% Subscripts, column-major strides and total class degree of every lattice
% point, rebuilt whenever an axis is created or discharged.
q = numel(dims);
P = prod(dims);
SUB = zeros(P, q);
rep = 1;
for k = 1:q
    SUB(:,k) = repmat(kron((0:dims(k)-1)', ones(rep,1)), P/(rep*dims(k)), 1);
    rep = rep*dims(k);
end
stride = cumprod([1, dims(1:end-1)]);
lev = sum(SUB(:,1:R), 2);
end

function A = pfqn_manjunath_expand(A, dims, k, newdim)
% Create marker k, keeping the existing content at degree zero: nothing
% multiplied in so far carries any power of it.
pre = prod(dims(1:k-1));
post = prod(dims(k+1:end));
T = zeros(pre, newdim, post);
T(:,1,:) = reshape(A, pre, 1, post);
A = T(:);
end

function A = pfqn_manjunath_reduce(A, dims, k, s, rhs)
% Discharge marker k. The multiplier (y^{b+1}-1)/(y-1) of a '<=' row turns its
% residue into the partial sum of the coefficients of degrees 0 ... b, and the
% multiplier 1/y^{b+1} of an '=' row picks the single coefficient of degree b.
pre = prod(dims(1:k-1));
dk = dims(k);
post = prod(dims(k+1:end));
T = reshape(A, pre, dk, post);
if s == 'E'
    A = reshape(T(:, rhs+1, :), [], 1);
else
    A = reshape(sum(T, 2), [], 1);
end
end
