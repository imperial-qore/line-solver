function pf = spn_pf(model, options)
% PF = SPN_PF(MODEL)
% PF = SPN_PF(MODEL, OPTIONS)
% Decide whether a stochastic Petri net has a product-form equilibrium
% distribution and, when it has, return the per-level factors g_l that
% MDD_REC and SPN_METRICS take as input.
%
% THIS IS THE PART THE MDD-REC PAPER DECLARES OUT OF SCOPE (FGCS Sec. 3.2).
% Every other function in api/spn/ receives the g_l already formed; this one
% derives them from the net, which is what lets a solver reach them.
%
% -- The theory, in one paragraph
% Write I(t), O(t) for the input and output vectors of mode t and lambda_t for
% its rate constant. Henderson-Taylor and Coleman-Henderson-Taylor show that a
% net whose firing rate has the form
%
%   r_t(m) = lambda_t psi(m - I(t)) / psi(m),     m >= I(t)
%
% has invariant measure pi(m) = psi(m) prod_l y_l^{m_l} whenever the positive
% vector y satisfies COMPLEX BALANCE: reading the distinct vectors appearing as
% some I(t) or O(t) as the COMPLEXES of the net, the flow into every complex
% must equal the flow out of it,
%
%   sum_{t : O(t)=v} lambda_t y^{I(t)} = ( sum_{t : I(t)=v} lambda_t ) y^v.
%
% Two choices of psi are realisable in LINE's own rate law, and they are the
% two this function tests for:
%
%   psi = 1              r_t = lambda_t, the rate of a SINGLE-SERVER mode.
%                        pi(m) = prod_l y_l^{m_l}, so g_l(k) = y_l^k.
%   psi = prod_l 1/m_l!  r_t = lambda_t prod_l m_l!/(m_l-I_l)!, MASS ACTION,
%                        reached in LINE through Transition.setFiringRateDependence
%                        or, for a mode drawing one token from one place, by
%                        infinite-server semantics.
%                        pi(m) = prod_l y_l^{m_l}/m_l!, so g_l(k) = y_l^k/k!.
%
% Which one holds is not guessed from the model API: the effective rate LINE
% would use, lambda_t min(enabling degree, servers) g(m), is EVALUATED at every
% reachable marking and compared against both laws. A net that matches neither
% under one common psi is refused by name, never approximated.
%
% -- Solving for y
% Complex balance reads A_lambda Psi(y) = 0 with A_lambda the Laplacian of the
% weighted digraph on complexes and Psi(y)_v = y^v. That Laplacian is the
% TRANSPOSED GENERATOR of a Markov chain that hops from complex to complex at
% the rate of the mode joining them, so its kernel on one linkage class is that
% chain's stationary distribution and CTMC_SOLVE returns it -- strictly
% positive exactly when the class is strongly connected, which is weak
% reversibility. With that positive vector kappa in hand y follows from the
% LINEAR system in x = log y
%
%   (v - v0) x = log kappa_v - log kappa_v0,   v, v0 in the same linkage class,
%
% solved in minimum norm. Feinberg's Deficiency Zero Theorem says this system
% is consistent for every choice of rate constants when the net is weakly
% reversible and its deficiency c - l - s is zero, which is why those two
% numbers are reported; but consistency is CHECKED rather than assumed, so a
% net of positive deficiency whose particular rates still admit a
% complex-balanced point is accepted on the evidence.
%
% -- The gauge, and why the minimum-norm solution is the canonical one
% Complex balance fixes y only up to y -> y .* exp(u) for any u orthogonal to
% the stoichiometric subspace S. Such a shift multiplies pi(m) by exp(u'm),
% which is CONSTANT on one compatibility class, so every reported measure is
% invariant under it -- but the normalising constant G itself is not, it scales
% by that constant. A gauge must therefore be FIXED, or the four codebases
% would report four different G on the same net. The one fixed here is
% x in the row space of the constraint matrix, i.e. the minimum-norm solution,
% and it is reached in a form that is unique whichever least-squares primitive
% a codebase carries: solve (rows*rows')w = rhs and set x = rows'*w. Any two
% solutions w of that system give the SAME rows'*w, so the answer does not
% depend on how the rank-deficient solve breaks its tie.
%
% -- Input
% MODEL   : a Network holding Places and Transitions
% OPTIONS : struct, fields
%             bound   - per-place-level token bound, passed to SPN_MDD
%             tol     - relative tolerance of the rate-law and complex-balance
%                       checks (default 1e-9)
%             verbose - print the certificate (default false)
% -- Output
% PF      : struct with fields
%             g          - 1 x L cell, g{l}(k+1) = g_l(k), ready for MDD_REC
%             y          - 1 x L positive vector solving complex balance
%             kind       - 'geometric' or 'massaction', the psi that was found
%             complexes  - C x L matrix of the distinct complexes
%             deficiency - c - l - s
%             linkage    - number of linkage classes
%             srank      - rank of the stoichiometric subspace
%             weaklyreversible - true when every linkage class is strongly connected
%             residual   - relative complex-balance residual at y
%             mdds, info - the reachable set and metadata SPN_MDD returned
%
% -- Reference
% J. L. Coleman, W. Henderson, P. G. Taylor, "Product form equilibrium
% distributions and a convolution algorithm for stochastic Petri nets",
% Performance Evaluation 26(3), 1996.
% M. Feinberg, "Complex balancing in general kinetic systems", Arch. Rational
% Mech. Anal. 49, 1972.
% D. F. Anderson, G. Craciun, T. G. Kurtz, "Product-form stationary
% distributions for deficiency zero chemical reaction networks", Bull. Math.
% Biol. 72, 2010.
%
% See also MDD_REC, SPN_METRICS, SPN_MDD, SPN_CONV.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(options), options = struct(); end
if ~isfield(options, 'tol'),     options.tol = 1e-9; end
if ~isfield(options, 'verbose'), options.verbose = false; end
if ~isfield(options, 'bound'),   options.bound = []; end

mddopt = struct('descriptor', false, 'bound', options.bound);
[mdds, ~, info] = spn_mdd(model, mddopt);

L = info.nplacelevels;
md = info.modes;
E = numel(md);
if E == 0
    line_error(mfilename, 'the net has no timed mode');
end

% ---- a queueing place holds an embedded server, not a token container
sn = model.getStruct();
for pp = 1:numel(info.places)
    ist = sn.nodeToStation(info.places(pp));
    if ist >= 1 && sn.sched(ist) ~= SchedStrategy.INF
        line_error(mfilename, sprintf(['place %s is a QUEUEING place (scheduling %s): its ' ...
            'embedded service is state that the marking does not carry, so the net is not the ' ...
            'token-container Petri net this product form is written for'], ...
            sn.nodenames{info.places(pp)}, SchedStrategy.toText(sn.sched(ist))));
    end
end

% ---- the rate constants and the structural vectors
lambda = zeros(1, E);
Iv = zeros(E, L); Ov = zeros(E, L);
for e = 1:E
    lambda(e) = md(e).D1(1);
    Iv(e, :) = md(e).enab;
    Ov(e, :) = md(e).fire;
    if any(isfinite(md(e).inhib))
        line_error(mfilename, sprintf(['mode %d of node %d has an inhibitor arc. An inhibitor ' ...
            'zeroes the firing rate on markings that still satisfy m >= I(t), so the rate is ' ...
            'not lambda*psi(m-I)/psi(m) on any psi and the net has no product form of this ' ...
            'kind'], md(e).mode, md(e).trans));
    end
    if md(e).srv ~= 1 && ~any(md(e).enab > 0)
        line_error(mfilename, sprintf(['mode %d of node %d has %g servers but consumes from no ' ...
            'place, so its enabling degree is unbounded and its firing rate undefined'], ...
            md(e).mode, md(e).trans, md(e).srv));
    end
    if ~(lambda(e) > 0)
        line_error(mfilename, sprintf('mode %d of node %d has a non-positive firing rate', ...
            md(e).mode, md(e).trans));
    end
end

% ---- which psi does LINE's own rate law follow on this net?
states = info.mdd.enumerate();
kind = i_ratelaw(states, md, lambda, info, options.tol);

% ---- complexes and the weighted digraph on them
[C, src, dst] = i_complexes(Iv, Ov);
c = size(C, 1);

% ---- Laplacian of the complex graph: A(j,i) = rate of the arc i -> j
A = zeros(c, c);
for e = 1:E
    if src(e) == dst(e), continue; end        % a mode that moves nothing
    A(dst(e), src(e)) = A(dst(e), src(e)) + lambda(e);
    A(src(e), src(e)) = A(src(e), src(e)) - lambda(e);
end

% ---- linkage classes, weak reversibility, deficiency
[lclass, nlink] = i_linkage(c, src, dst);
wr = i_weaklyreversible(c, src, dst, lclass, nlink);
srank = rank(Ov - Iv, 1e-9);
deficiency = c - nlink - srank;

% ---- kappa: the positive kernel of the Laplacian on each linkage class
kappa = zeros(c, 1);
for b = 1:nlink
    idx = find(lclass == b);
    if numel(idx) == 1
        kappa(idx) = 1; continue
    end
    % A(idx,idx) is the transposed generator of the complex-hopping chain, so
    % its kernel is that chain's stationary law.
    v = ctmc_solve(A(idx, idx)')';
    if any(v <= 0)
        line_error(mfilename, sprintf(['linkage class %d of the complex graph carries no flow ' ...
            'through complex %d, so the net admits no positive complex-balanced point. A weakly ' ...
            'reversible net has a strictly positive balance flow on every linkage class; this ' ...
            'one is %s'], b, idx(find(v <= 0, 1)), i_wrtext(wr)));
    end
    kappa(idx) = v / max(v);
end

% ---- x = log y from the linear system on each linkage class
rows = zeros(0, L); rhs = zeros(0, 1);
for b = 1:nlink
    idx = find(lclass == b);
    v0 = idx(1);
    for jj = 2:numel(idx)
        v = idx(jj);
        rows(end + 1, :) = C(v, :) - C(v0, :);                      %#ok<AGROW>
        rhs(end + 1, 1) = log(kappa(v)) - log(kappa(v0));           %#ok<AGROW>
    end
end
if isempty(rows)
    x = zeros(L, 1);
else
    % minimum norm through the row space; see the gauge note in the header
    x = rows' * (pinv(rows * rows') * rhs);
    res = norm(rows * x - rhs, inf);
    scale = max(1, norm(rhs, inf));
    if res > options.tol * scale
        line_error(mfilename, sprintf(['the complex-balance equations are inconsistent (residual ' ...
            '%.3e): this net has no product form of the tested kind at these rates. Its ' ...
            'deficiency is %d and it is %s; the Deficiency Zero Theorem guarantees a solution ' ...
            'only at deficiency 0 with weak reversibility'], res, deficiency, i_wrtext(wr)));
    end
end
y = exp(x(:)');

% ---- verify complex balance itself, which is what makes pi stationary
Psi = ones(c, 1);
for v = 1:c, Psi(v) = prod(y .^ C(v, :)); end
resb = norm(A * Psi, inf);
scaleb = max(norm(abs(A) * Psi, inf), realmin);
if resb > options.tol * scaleb
    line_error(mfilename, sprintf(['complex balance fails at the computed point (relative ' ...
        'residual %.3e), so the product form would not be stationary'], resb / scaleb));
end

% ---- the per-level factors, tabulated over the reachable domain
g = cell(1, L);
for l = 1:L
    d = mdds.domain(l);
    k = 0:(d - 1);
    switch kind
        case 'geometric'
            g{l} = y(l) .^ k;
        case 'massaction'
            g{l} = (y(l) .^ k) ./ gamma(k + 1);
    end
end

pf = struct('g', {g}, 'y', y, 'kind', kind, 'complexes', C, ...
    'deficiency', deficiency, 'linkage', nlink, 'srank', srank, ...
    'weaklyreversible', wr, 'residual', resb / scaleb, 'mdds', mdds, 'info', info);

if options.verbose
    line_printf('\nSPN product form: %s, %d complexes, %d linkage classes, rank %d, deficiency %d, %s\n', ...
        kind, c, nlink, srank, deficiency, i_wrtext(wr));
    line_printf('  y = %s (complex-balance residual %.2e)\n', mat2str(y, 6), resb / scaleb);
end
end

% ------------------------------------------------------------------------
function kind = i_ratelaw(states, md, lambda, info, tol)
% Evaluate the rate LINE would actually use at every reachable marking and
% report which psi reproduces it. Both candidates are tried on every mode; a
% net matching neither, or matching different ones on different modes, has no
% product form of this family and is refused by name.
E = numel(md);
L = info.nplacelevels;
okgeo = true; okma = true;
for s = 1:size(states, 1)
    m = states(s, 1:L);
    Marc = i_arcmatrix(m, info);
    for e = 1:E
        if any(m < md(e).enab), continue; end
        actual = lambda(e) * i_servers(m, md(e));
        if ~isempty(md(e).dep)
            actual = actual * md(e).dep(Marc);
        end
        okgeo = okgeo && i_close(actual, lambda(e), tol);
        okma  = okma  && i_close(actual, lambda(e) * i_massaction(m, md(e).enab), tol);
        if ~okgeo && ~okma
            line_error('spn_pf', sprintf(['mode %d of node %d fires at rate %g in a reachable ' ...
                'marking, which is neither its rate constant (single-server, psi = 1) nor its ' ...
                'mass-action rate %g (psi = prod 1/m!). LINE''s rate law on this mode is ' ...
                'lambda*min(enabling degree, servers)*g(m), and no psi puts that in the form ' ...
                'lambda*psi(m-I)/psi(m)'], md(e).mode, md(e).trans, actual, ...
                lambda(e) * i_massaction(m, md(e).enab)));
        end
    end
end
if okgeo
    kind = 'geometric';
else
    kind = 'massaction';
end
end

% ------------------------------------------------------------------------
function n = i_servers(m, mde)
% min(enabling degree, servers): the number of sets of tokens firing at once.
deg = Inf;
for l = find(mde.enab > 0)
    deg = min(deg, floor(m(l) / mde.enab(l)));
end
if isinf(deg), deg = 1; end                   % consumes nothing: always one set
n = min(deg, mde.srv);
end

% ------------------------------------------------------------------------
function r = i_massaction(m, enab)
% prod_l m_l!/(m_l - I_l)!, the number of ordered ways to pick the input tokens.
r = 1;
for l = find(enab > 0)
    for j = 0:(enab(l) - 1)
        r = r * (m(l) - j);
    end
end
end

% ------------------------------------------------------------------------
function tf = i_close(a, b, tol)
tf = abs(a - b) <= tol * max(1, max(abs(a), abs(b)));
end

% ------------------------------------------------------------------------
function M = i_arcmatrix(m, info)
% Place-major level vector -> the (nnodes x nclasses) marking matrix that a
% Transition.setFiringRateDependence handle is written against.
R = info.nclasses;
M = zeros(info.nnodes, R);
for pp = 1:numel(info.places)
    for k = 1:R
        M(info.places(pp), k) = m((pp - 1) * R + k);
    end
end
end

% ------------------------------------------------------------------------
function [C, src, dst] = i_complexes(Iv, Ov)
% The distinct input and output vectors, and the arc each mode draws.
E = size(Iv, 1);
[C, ~, ic] = unique([Iv; Ov], 'rows', 'stable');
src = ic(1:E)';
dst = ic(E + (1:E))';
end

% ------------------------------------------------------------------------
function [lclass, nlink] = i_linkage(c, src, dst)
% Connected components of the UNDIRECTED complex graph.
lclass = zeros(1, c);
nlink = 0;
for v = 1:c
    if lclass(v) > 0, continue; end
    nlink = nlink + 1;
    stack = v; lclass(v) = nlink;
    while ~isempty(stack)
        u = stack(end); stack(end) = [];
        nb = [dst(src == u), src(dst == u)];
        for w = nb
            if lclass(w) == 0
                lclass(w) = nlink; stack(end + 1) = w; %#ok<AGROW>
            end
        end
    end
end
end

% ------------------------------------------------------------------------
function tf = i_weaklyreversible(c, src, dst, lclass, nlink)
% Every linkage class strongly connected in the DIRECTED complex graph.
tf = true;
for b = 1:nlink
    idx = find(lclass == b);
    reach = i_reach(idx(1), src, dst, c);
    back = i_reach(idx(1), dst, src, c);
    if ~all(reach(idx)) || ~all(back(idx))
        tf = false; return
    end
end
end

% ------------------------------------------------------------------------
function seen = i_reach(v0, from, to, c)
seen = false(1, c); seen(v0) = true; stack = v0;
while ~isempty(stack)
    u = stack(end); stack(end) = [];
    for w = to(from == u)
        if ~seen(w), seen(w) = true; stack(end + 1) = w; end %#ok<AGROW>
    end
end
end

% ------------------------------------------------------------------------
function s = i_wrtext(wr)
if wr, s = 'weakly reversible'; else, s = 'not weakly reversible'; end
end
