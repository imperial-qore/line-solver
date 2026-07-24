%{ @file qbd_rap.m
 %  @brief Solves a general QBD process with Rational Arrival Process components
 %
 %  @author LINE Development Team
%}

%{
 % @brief Equilibrium analysis of a Quasi-Birth-and-Death process with
 % Rational Arrival Process (RAP) components
 %
 % @details
 % The process is specified directly by its level-independent blocks
 % (A0,A1,A2) and its boundary blocks (B0,B1), where A0 drives level
 % increases, A2 drives level decreases and A1 the within-level evolution.
 % Unlike a Markovian QBD the blocks need not be nonnegative: they are only
 % required to be conservative, (A0+A1+A2)*e = 0, and to define a genuine
 % RAP through the prediction-process interpretation. This makes qbd_rap
 % strictly more general than qbd_raprap1, which builds a product-space QBD
 % from two INDEPENDENT RAPs; here the arrival process and the sequence of
 % service times may be driven from a shared phase space and therefore be
 % cross-correlated.
 %
 % @par References:
 % N. G. Bean and B. F. Nielsen, "Quasi-Birth-and-Death Processes with
 % Rational Arrival Process Components", Stochastic Models, 26(3), 2010,
 % pp. 309-334 (DTU technical report IMM-2007-20). The equilibrium
 % construction below is their Theorem 7 and the stability test is their
 % Corollary 8. The argument rests on the prediction-process interpretation
 % of a RAP due to Asmussen and Bladt, which is what allows a QBD argument
 % to be carried over to matrices that are not nonnegative; the same
 % prediction process underlies the conditional-vector RAP sampler in
 % rap_sample.m.
 %
 % @par Algorithm (Theorem 7):
 % 1. Solve A0*G^2 + A1*G + A2 = 0 for G.
 % 2. U = A1 + A0*G.
 % 3. R = A0*inv(-U).
 % 4. Find the row vector pihat0 with pihat0*(B1 + R*A2) = 0, pihat0*e = 1.
 % 5. pi0 = K*pihat0 with K chosen so that pi0*inv(I-R)*e = 1.
 % 6. pi_n = pi0*R^n, and the marginal level probability is pi_n*e.
 % The process is positive recurrent iff Sp(R) < 1 and step 4 has a solution.
 %
 % @par Computation of G:
 % The blocks are not nonnegative, so the probabilistic iterations used for
 % Markovian QBDs (logarithmic reduction, cyclic reduction) carry no
 % convergence guarantee here, and the paper explicitly leaves the general
 % case open ("The issue of justifying algorithms for the evaluation of the
 % matrix G for such processes has not been undertaken", Section 6). Two
 % paths are therefore taken:
 % - If A2 has rank one, A2 = u*v, then G = e*v/(v*e) solves the equation
 %   exactly. This is immediate from conservativity: G is idempotent and
 %   (A0+A1)*e = -A2*e = -u*(v*e), so (A0+A1)*e*v/(v*e) = -u*v = -A2. This
 %   is the case covered by the paper's example, and the residual is checked.
 % - Otherwise the quadratic matrix equation is solved numerically by
 %   natural functional iteration G <- inv(-A1)*(A2 + A0*G^2) used as a warm
 %   start, followed by Newton's method on the Sylvester-form Jacobian,
 %   (A0*G+A1)*H + A0*H*G = -(A0*G^2 + A1*G + A2), solved through its
 %   Kronecker expansion. If the residual does not reach roundoff level, or
 %   the iterate does not satisfy the stochastic-analogue constraint G*e = e,
 %   an error is raised rather than returning an unconverged G.
 %
 % @par Syntax:
 % @code
 % [levelProb, QN, R, G, U, spr, pqueue, pi0] = qbd_rap(A0, A1, A2)
 % [levelProb, QN, R, G, U, spr, pqueue, pi0] = qbd_rap(A0, A1, A2, B0, B1)
 % [levelProb, QN, R, G, U, spr, pqueue, pi0] = qbd_rap(A0, A1, A2, B0, B1, numLevels)
 % @endcode
 %
 % This is the block-level core of the RAP QBD family. qbd_raprap1 is the
 % thin wrapper over it that builds the product-space blocks of two
 % INDEPENDENT RAPs; callers with a coupled model, in which arrivals and
 % services share a phase space, must call qbd_rap directly because no
 % product form exists to factor out.
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>A0<td>Level-up block (m x m)
 % <tr><td>A1<td>Local block (m x m)
 % <tr><td>A2<td>Level-down block (m x m)
 % <tr><td>B0<td>(Optional) boundary level-up block, default A0
 % <tr><td>B1<td>(Optional) boundary local block, default A1
 % <tr><td>numLevels<td>(Optional) highest level reported, default 20
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>levelProb<td>Row vector of marginal level probabilities, levels 0..numLevels
 % <tr><td>QN<td>Mean queue length (exact, pi0*R*inv(I-R)^2*e)
 % <tr><td>R<td>Rate matrix R
 % <tr><td>G<td>Matrix G solving A0*G^2 + A1*G + A2 = 0
 % <tr><td>U<td>Matrix U = A1 + A0*G
 % <tr><td>spr<td>Spectral radius Sp(R)
 % <tr><td>pqueue<td>(numLevels+1) x m matrix of the vectors pi_n
 % <tr><td>pi0<td>Level-0 vector pi_0, the boundary vector of Theorem 7
 % </table>
 %
 % @par Example:
 % @code
 % g  = 0.25;
 % A1 = [-1 0 0; -2/3 -1 1; 2/3 -1 -1];
 % Da = [14/5 -9/10 -9/10; 26/15 -8/15 -8/15; 58/15 -19/15 -19/15];
 % Ds = [1; 2/3; 4/3]*[3 -1 -1];
 % [levelProb, QN] = qbd_rap(g*Da, A1, (1-g)*Ds, g*Da, g*A1, 8);
 % % levelProb reproduces Table 1 of Bean and Nielsen (2010)
 % @endcode
%}
function [levelProb, QN, R, G, U, spr, pqueue, pi0] = qbd_rap(A0, A1, A2, B0, B1, numLevels)
% [LEVELPROB,QN,R,G,U,SPR,PQUEUE,PI0]=QBD_RAP(A0,A1,A2,B0,B1,NUMLEVELS)

if nargin < 3
    line_error(mfilename, 'qbd_rap requires at least the blocks A0, A1 and A2.');
end
if nargin < 4 || isempty(B0)
    B0 = A0;
end
if nargin < 5 || isempty(B1)
    B1 = A1;
end
if nargin < 6 || isempty(numLevels)
    numLevels = 20;
end

m = size(A1,1);
if size(A1,2) ~= m || any(size(A0) ~= [m m]) || any(size(A2) ~= [m m]) ...
        || any(size(B0) ~= [m m]) || any(size(B1) ~= [m m])
    line_error(mfilename, 'All QBD blocks must be square and of the same order.');
end
if numLevels < 0 || numLevels ~= round(numLevels)
    line_error(mfilename, 'numLevels must be a nonnegative integer.');
end

e = ones(m,1);
I = eye(m);
blockScale = max([norm(A0,'fro'), norm(A1,'fro'), norm(A2,'fro'), 1]);

% Conservativity of the repeating portion, (A0+A1+A2)*e = 0. This is the
% RAP analogue of the generator row-sum condition and every step below
% relies on it.
if norm((A0+A1+A2)*e, inf) > 1e-8*blockScale
    line_error(mfilename, sprintf(['The repeating blocks are not conservative: ' ...
        '||(A0+A1+A2)*e||_inf = %g. A QBD with RAP components requires ' ...
        '(A0+A1+A2)*e = 0.'], norm((A0+A1+A2)*e, inf)));
end
if norm((B0+B1)*e, inf) > 1e-8*blockScale
    line_error(mfilename, sprintf(['The boundary blocks are not conservative: ' ...
        '||(B0+B1)*e||_inf = %g. A QBD with RAP components requires ' ...
        '(B0+B1)*e = 0 at level 0.'], norm((B0+B1)*e, inf)));
end

% Step 1: matrix G.
G = qbd_rap_g(A0, A1, A2, blockScale);

% Steps 2 and 3: U and R.
U = A1 + A0*G;
if rcond(-U) < eps
    line_error(mfilename, 'The matrix U = A1 + A0*G is singular, R = A0*inv(-U) does not exist.');
end
R = A0/(-U);

% Corollary 8(i): positive recurrence.
spr = max(abs(eig(R)));
% The threshold carries a 1e-12 margin: at the null-recurrent boundary
% Sp(R) equals 1 in exact arithmetic but rounds to either side, and the
% three codebases must agree on rejecting it. Any model within 1e-12 of the
% boundary has an unbounded queue regardless.
if spr >= 1 - 1e-12
    line_error(mfilename, sprintf(['The process is not positive recurrent: ' ...
        'Sp(R) = %.15g >= 1 (Corollary 8 of Bean and Nielsen, 2010).'], spr));
end

% Step 4: boundary vector, pihat0*(B1 + R*A2) = 0 normalised to pihat0*e = 1.
% The left null space is extracted from the SVD of V.' with a relative
% tolerance: V is only singular up to the accuracy with which R was
% computed, so the absolute default tolerance of null() is too strict here.
V = B1 + R*A2;
[~, Sv, W] = svd(V.');
sv = diag(Sv);
nullTol = 1e-8*max(sv(1),1);
if sv(end) > nullTol
    line_error(mfilename, sprintf(['The boundary equation x*(B1 + R*A2) = 0 has no ' ...
        'nontrivial solution (smallest singular value %g against tolerance %g), so ' ...
        'the process is not positive recurrent (Corollary 8(ii) of Bean and ' ...
        'Nielsen, 2010).'], sv(end), nullTol));
end
if m > 1 && sv(end-1) <= nullTol
    line_error(mfilename, ['The boundary equation x*(B1 + R*A2) = 0 has a solution ' ...
        'space of dimension greater than one, the equilibrium vector is not unique.']);
end
pihat0 = W(:,end).';
if abs(pihat0*e) < 1e-12*norm(pihat0,inf)
    line_error(mfilename, 'The boundary vector cannot be normalised, x*e = 0.');
end
pihat0 = pihat0/(pihat0*e);

% Step 5: level-0 vector.
K = 1/(pihat0*((I-R)\e));
pi0 = K*pihat0;

% Consistency of the supplied boundary up-block: the level-0 balance
% equation pi0*B0 + pi1*A1 + pi2*A2 = 0 must hold with pi_n = pi0*R^n.
bal = pi0*B0 + pi0*R*A1 + pi0*R*R*A2;
if norm(bal, inf) > 1e-8*blockScale*max(norm(pi0,inf),1)
    line_error(mfilename, sprintf(['The boundary block B0 is inconsistent with the ' ...
        'repeating blocks: ||pi0*B0 + pi1*A1 + pi2*A2||_inf = %g. The level-0 ' ...
        'balance equation of Theorem 7 requires pi0*(B0-A0) = 0.'], norm(bal,inf)));
end

% Step 6: level vectors and marginal level distribution.
pqueue = zeros(numLevels+1, m);
pin = pi0;
for n = 0:numLevels
    pqueue(n+1,:) = pin;
    pin = pin*R;
end
levelProb = (pqueue*e).';

% Exact mean queue length, sum_n n*pi0*R^n*e = pi0*R*inv(I-R)^2*e.
QN = pi0*R*((I-R)\((I-R)\e));
end

% see _kb/03-api-layer.md (QBD boundary conventions) for rationale
function G = qbd_rap_g(A0, A1, A2, blockScale)
m = size(A1,1);
e = ones(m,1);
resTol = 1e-10*blockScale;

% Rank-one A2 = u*v admits the closed form G = e*v/(v*e), which is the case
% deliberately chosen in the example of Bean and Nielsen (2010) precisely so
% that G is available a priori.
sv = svd(A2);
if numel(sv) > 1 && sv(1) > 0 && sv(2) <= 1e-10*sv(1)
    [~,~,W] = svd(A2);
    v = W(:,1).';
    if abs(v*e) < 1e-12*norm(v,inf)
        line_error(mfilename, ['A2 has rank one but its right factor v satisfies ' ...
            'v*e = 0, so the closed form G = e*v/(v*e) is undefined.']);
    end
    G = e*(v/(v*e));
    res = norm(A0*G*G + A1*G + A2, 'fro');
    if res > resTol
        line_error(mfilename, sprintf(['The rank-one closed form for G leaves a ' ...
            'residual ||A0*G^2 + A1*G + A2||_F = %g, which is above the roundoff ' ...
            'level %g.'], res, resTol));
    end
    return
end

% General case. Natural functional iteration G <- inv(-A1)*(A2 + A0*G^2)
% gives a warm start; it is the standard Markovian iteration but has no
% convergence guarantee for blocks that are not nonnegative.
if rcond(-A1) < eps
    line_error(mfilename, 'The local block A1 is singular, the iteration for G cannot be started.');
end
G = zeros(m);
for it = 1:200
    Gnew = (-A1)\(A2 + A0*G*G);
    if norm(Gnew-G, 'fro') <= 1e-14*max(norm(G,'fro'),1)
        G = Gnew;
        break
    end
    G = Gnew;
    if ~all(isfinite(G(:)))
        break
    end
end
if ~all(isfinite(G(:)))
    G = zeros(m);
end

% Newton's method on F(G) = A0*G^2 + A1*G + A2. The derivative in the
% direction H is (A0*G+A1)*H + A0*H*G, a Sylvester-type operator solved here
% through its Kronecker expansion (I kron (A0*G+A1) + G' kron A0)*vec(H).
Im = eye(m);
for it = 1:100
    res = A0*G*G + A1*G + A2;
    if norm(res, 'fro') <= resTol
        break
    end
    J = kron(Im, A0*G + A1) + kron(G.', A0);
    if rcond(J) < eps
        break
    end
    H = reshape(-(J\res(:)), m, m);
    G = G + H;
    if ~all(isfinite(G(:)))
        break
    end
end

res = Inf;
if all(isfinite(G(:)))
    res = norm(A0*G*G + A1*G + A2, 'fro');
end
if ~(res <= resTol) || norm(G*e - e, inf) > 1e-8
    line_error(mfilename, sprintf(['Could not compute the matrix G for this QBD with ' ...
        'RAP components: residual ||A0*G^2 + A1*G + A2||_F = %g against a tolerance ' ...
        'of %g, and ||G*e-e||_inf = %g. The blocks are not nonnegative, so neither ' ...
        'the functional iteration nor Newton''s method is guaranteed to converge, and ' ...
        'the justification of algorithms for G in this setting is left as an open ' ...
        'problem in Section 6 of N. G. Bean and B. F. Nielsen, "Quasi-Birth-and-Death ' ...
        'Processes with Rational Arrival Process Components", Stochastic Models, 26(3), ' ...
        '2010, pp. 309-334. Supply a model with a rank-one A2, for which G is available ' ...
        'in closed form.'], res, resTol, norm(G*e-e,inf)));
end
end
