function imm = fluid_petri_immediate(terms, x, imm)
% IMM = FLUID_PETRI_IMMEDIATE(TERMS, X)
% IMM = FLUID_PETRI_IMMEDIATE(TERMS, X, IMM)
%
% The active set of the IMMEDIATE transitions, and the equations that pin their
% flows.
%
% AN IMMEDIATE TRANSITION HAS NO RATE, so it cannot be a column of the drift
% with a rate factor like every other event. Its fluid limit is a singular
% perturbation: the transition fires infinitely fast, so a marking that enables
% it cannot persist, and what survives in the limit is a FLOW. That flow is an
% algebraic unknown of the DAE, pinned by the constraint that the transition's
% binding input place holds no mass:
%
%   phi_j >= 0,   x_b = 0 for the coordinate b that binds mode j
%
% which is the complementarity condition of the limit, and is the same object
% the capacity caps already solve as an active set -- an inequality has no
% residual for a Newton solver, so the loop iterates over WHICH constraints
% bind and each pass is an equality-constrained solve.
%
% TWO MODES DRAINING ONE PLACE NEED ONE MORE EQUATION THAN THE PIN GIVES, and
% the extra equation is the GSPN's own conflict rule: among the enabled
% immediate modes of highest firing priority the branch is taken in proportion
% to the firing weights, so
%
%   phi_j * weight_l = phi_l * weight_j
%
% for the members of one conflict group, and phi = 0 for a member of lower
% priority. That is exactly what the exact engines do at a vanishing marking
% (SOLVER_SSA_NRM's SPNCOLLAPSE, STATE.AFTERGLOBALEVENT's immediate branch),
% carried over to a continuous flow.
%
% THE COUNT IS SQUARE BY CONSTRUCTION. Each active mode is assigned exactly one
% binding coordinate; the pinned coordinates are the image of that assignment,
% so V pins plus (F - V) ratio rows is F equations for F flows. A coordinate
% that happens to be empty without being anybody's binding coordinate is not
% pinned -- its own drift row determines it -- and if the solve then pushes it
% negative, the caller REBINDS a mode to it and solves again, which is the
% active-set move.
%
% Parameters:
%   terms - FLUID_PETRI_TERMS output
%   x     - state to read the assignment off (the seed, or the last iterate)
%   imm   - a previous assignment to refresh; omitted to initialise
%
% Returns:
%   imm - struct with fields
%           n       - number of immediate modes
%           active  - (n x 1) logical, whether the mode carries a flow unknown
%           bind    - (n x 1) marking coordinate pinned for that mode, 0 if none
%           pins    - the distinct pinned coordinates
%           rows    - (nrow x 1) struct array describing each equation:
%                       .kind 'pin' | 'ratio' | 'zero'
%                       .a, .b, .wa, .wb as the kind requires
%
% See also SOLVER_FLUID_PETRI, FLUID_PETRI_TERMS.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

n = numel(terms.immIdx);
if nargin < 3 || isempty(imm)
    imm = struct('n', n, 'active', true(n,1), 'bind', zeros(n,1), 'pins', [], 'rows', []);
    % An inhibited mode never fires, so it neither carries a flow nor empties a
    % place: leaving it active would pin a coordinate the model does not empty.
    % An empty input place is NOT a reason to deactivate -- that is the normal
    % state of an enabled immediate mode, and its flow is then whatever the
    % inflow is.
    for k = 1:n
        j = terms.immIdx(k);
        md = terms.modes(j);
        if isempty(md.arcSlot)
            line_error(mfilename, sprintf(['Immediate mode %s has no enabling arc, so nothing bounds its ' ...
                'firing flow and the net has no fluid limit. Give it an input place, or make it timed.'], md.label));
        end
        if i_inhibited(terms, md, x)
            imm.active(k) = false;
        end
    end
end

% ---- the assignment: each active mode binds the input arc it is shortest of
for k = 1:n
    if ~imm.active(k)
        imm.bind(k) = 0;
        continue
    end
    if imm.bind(k) > 0 && any(terms.modes(terms.immIdx(k)).arcSlot == imm.bind(k))
        continue % a binding the caller set explicitly is kept
    end
    md = terms.modes(terms.immIdx(k));
    lev = x(md.arcSlot(:)) ./ md.arcW(:);
    [~, a] = min(lev);
    imm.bind(k) = md.arcSlot(a);
end

% ---- the equations
imm.pins = unique(imm.bind(imm.active & imm.bind > 0));
rows = struct('kind',{},'a',{},'b',{},'wa',{},'wb',{});
for p = imm.pins(:)'
    rows(end+1) = struct('kind','pin','a',p,'b',0,'wa',0,'wb',0); %#ok<AGROW>
    grp = find(imm.active & imm.bind == p);
    if numel(grp) <= 1
        continue
    end
    prio = zeros(numel(grp),1);
    wgt = zeros(numel(grp),1);
    for t = 1:numel(grp)
        prio(t) = terms.modes(terms.immIdx(grp(t))).prio;
        wgt(t) = terms.modes(terms.immIdx(grp(t))).weight;
    end
    top = grp(prio == max(prio));
    low = grp(prio < max(prio));
    tw = wgt(prio == max(prio));
    for t = 2:numel(top)
        rows(end+1) = struct('kind','ratio','a',top(1),'b',top(t), ...
            'wa',tw(1),'wb',tw(t)); %#ok<AGROW>
    end
    for t = 1:numel(low)
        rows(end+1) = struct('kind','zero','a',low(t),'b',0,'wa',0,'wb',0); %#ok<AGROW>
    end
end
for k = find(~imm.active(:)).'
    rows(end+1) = struct('kind','zero','a',k,'b',0,'wa',0,'wb',0); %#ok<AGROW>
end
imm.rows = rows;
end

% -------------------------------------------------------------------------
function tf = i_inhibited(terms, md, x) %#ok<INUSL>
% True when an inhibitor arc of this mode has reached its threshold, which
% disables the mode outright whatever its input places hold.
%
% A HARD TEST ON THE MEAN, AND A KNOWN WRONG ANSWER WHEN THE MEAN SITS ON THE
% THRESHOLD. A timed mode closes the same indicator as Phi((thr-m)/sd) in
% FLUID_PETRI_THETA; this path does not, so a mode whose inhibitor place hovers
% at its threshold is switched off outright and its whole downstream branch
% carries zero flow (spn_open_sevenplaces: T5 dead, P7 at 0 against 0.207 under
% LDES). The smoothed gate ALONE would not fix it -- an immediate mode's input
% is pinned at zero, and a gate strictly inside (0,1) breaks that premise, since
% the input accumulates while the mode is inhibited. See
% _kb/06-solver-catalog.md for the measurement and for what a real fix costs.
tf = false;
for b = 1:numel(md.inhSlot)
    if x(md.inhSlot(b)) >= md.inhThr(b)
        tf = true;
        return
    end
end
end
