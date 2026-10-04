function th = fluid_petri_theta(terms, x, s2)
% TH = FLUID_PETRI_THETA(TERMS, X, S2)
%
% The closed enabling term of every transition mode, and its derivative.
%
% A mode fires at a rate proportional to the number of BINDINGS the marking
% supports, min over its input arcs of m_a/w_a, capped at its server count and
% zeroed by any inhibitor arc that has reached its threshold. Under the Gaussian
% marginal that whole expression is
%
%   theta_j = ( prod_b Phi((thr_b - m_b)/sd_b) ) * E[ min_a(m_a/w_a), c_j ]
%
% where the many-argument min is closed by FLUID_MINMULTI_CLOSURE and the
% inhibitor INDICATOR is closed by the normal CDF. Both collapse to their
% first-order form at zero variance -- Phi becomes 1{m_b < thr_b} and the min
% closure becomes min() -- so the mean-field limit is one code path, not two.
%
% THE INHIBITOR GATE IS WHY A PETRI NET NEEDS A SMOOTHED CLOSURE AT ALL, quite
% apart from the accuracy argument. 1{m_b < thr_b} is a step, and a Newton
% solver has no derivative to descend on a step: the drift residual would be
% piecewise constant in m_b and the iterate would either sit still or chatter
% across the threshold. Phi((thr_b - m_b)/sd_b) is the same function smoothed by
% the marginal the closure already carries, so the Jacobian is defined
% everywhere the variance is positive.
%
% THE VARIANCES ARE UNKNOWNS, NOT FUNCTIONS OF X. S2 holds one entry per
% TERMS.covPairs row, pinned by its own consistency row in the DAE, so the
% derivative below is with respect to the MEANS only -- exactly how
% FLUID_DRIFT_JACOBIAN differentiates the queueing closure.
%
% Parameters:
%   terms - FLUID_PETRI_TERMS output
%   x     - state vector [marking; phase distributions]
%   s2    - (npair x 1) the Sigma entries the closure reads
%
% Returns:
%   th - struct with fields
%          theta   - (nmodes x 1) closed enabling term
%          dslot   - (nmodes x 1) cell of coordinate indices
%          dval    - (nmodes x 1) cell of dTHETA/dX at those coordinates
%          dep     - (nmodes x 1) marking-dependent firing multiplier
%          depslot - (nmodes x 1) cell of coordinate indices
%          depval  - (nmodes x 1) cell of dDEP/dX at those coordinates
%
% See also FLUID_MINMULTI_CLOSURE, FLUID_PETRI_RATES, FLUID_PETRI_JACOBIAN.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

nmod = numel(terms.modes);
th = struct();
th.theta = zeros(nmod,1);
th.dslot = cell(nmod,1);
th.dval = cell(nmod,1);
th.dep = ones(nmod,1);
th.depslot = cell(nmod,1);
th.depval = cell(nmod,1);

if isempty(s2)
    s2 = zeros(terms.npair,1);
end

for j = 1:nmod
    md = terms.modes(j);
    A = numel(md.arcSlot);
    mu = zeros(A,1);
    Sarg = zeros(A);
    for a = 1:A
        mu(a) = x(md.arcSlot(a)) / md.arcW(a);
    end
    if md.closable
        for a = 1:A
            for b = 1:A
                Sarg(a,b) = i_sig(terms, s2, md.arcSlot(a), md.arcSlot(b)) / (md.arcW(a)*md.arcW(b));
            end
        end
    end
    [hmin, gmin] = fluid_minmulti_closure(mu, Sarg, md.c);

    % the inhibitor gates, and the product rule over them
    nb = numel(md.inhSlot);
    gate = ones(nb,1);
    dgate = zeros(nb,1);
    for b = 1:nb
        mb = x(md.inhSlot(b));
        thr = md.inhThr(b);
        vb = i_sig(terms, s2, md.inhSlot(b), md.inhSlot(b));
        if vb > 0
            sd = sqrt(vb);
            zb = (thr - mb)/sd;
            gate(b) = 0.5*erfc(-zb/sqrt(2));
            dgate(b) = -exp(-0.5*zb^2)/(sqrt(2*pi)*sd);
        else
            gate(b) = double(mb < thr - GlobalConstants.FineTol*max(1,thr));
            dgate(b) = 0;
        end
    end
    ginh = prod(gate);

    slots = zeros(0,1); vals = zeros(0,1);
    for a = 1:A
        slots(end+1,1) = md.arcSlot(a); %#ok<AGROW>
        vals(end+1,1) = ginh * gmin(a) / md.arcW(a); %#ok<AGROW>
    end
    for b = 1:nb
        if gate(b) ~= 0
            others = ginh / gate(b);
        else
            others = prod(gate([1:b-1, b+1:nb]));
        end
        slots(end+1,1) = md.inhSlot(b); %#ok<AGROW>
        vals(end+1,1) = hmin * others * dgate(b); %#ok<AGROW>
    end
    % an arc and an inhibitor arc may share a coordinate, so accumulate
    [uslots, ~, ic] = unique(slots);
    uvals = accumarray(ic, vals, [numel(uslots), 1]);

    th.theta(j) = ginh * hmin;
    th.dslot{j} = uslots;
    th.dval{j} = uvals;

    if ~isempty(md.dep)
        [g, gslot, gval] = i_dep(terms, md, x);
        th.dep(j) = g;
        th.depslot{j} = gslot;
        th.depval{j} = gval;
    end
end
end

% -------------------------------------------------------------------------
function v = i_sig(terms, s2, a, b)
% One entry of the closure covariance, zero where the drift never reads it.
t = terms.pairIndex(a,b);
if t == 0
    v = 0;
else
    v = s2(t);
end
end

% -------------------------------------------------------------------------
function [g, slots, vals] = i_dep(terms, md, x)
% The marking-dependent firing multiplier and its gradient.
%
% G is a user function of the (nnodes x nclasses) marking matrix, so it is
% evaluated at the MEAN marking -- a first-order closure of g, which is the
% same order at which SolverCTMC evaluates it per state and the only one
% available without the distribution of the marking. Its gradient has no
% analytic form, so it is taken by central differences over the marking
% coordinates.
mm = zeros(terms.I, terms.K);
for s = 1:terms.nm
    mm(terms.coordNode(s), terms.coordClass(s)) = x(s);
end
g = double(md.dep(mm));
slots = (1:terms.nm)';
vals = zeros(terms.nm,1);
for s = 1:terms.nm
    h = max(1e-6*abs(x(s)), 1e-6);
    mp = mm; mp(terms.coordNode(s), terms.coordClass(s)) = mm(terms.coordNode(s), terms.coordClass(s)) + h;
    mn = mm; mn(terms.coordNode(s), terms.coordClass(s)) = max(0, mm(terms.coordNode(s), terms.coordClass(s)) - h);
    hh = mp(terms.coordNode(s), terms.coordClass(s)) - mn(terms.coordNode(s), terms.coordClass(s));
    if hh <= 0
        continue
    end
    vals(s) = (double(md.dep(mp)) - double(md.dep(mn)))/hh;
end
keepv = vals ~= 0;
slots = slots(keepv);
vals = vals(keepv);
end
