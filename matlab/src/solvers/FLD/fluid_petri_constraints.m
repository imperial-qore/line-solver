function con = fluid_petri_constraints(sn, terms)
% CON = FLUID_PETRI_CONSTRAINTS(SN, TERMS)
%
% Every finite place capacity as a linear constraint on the fluid marking.
%
%     Arow * x <= b
%
% with one row per finite total capacity (SETCAPACITY on a Place) and one per
% finite per-class capacity (SETCLASSCAPACITY). This is the same family of rows
% the queueing DAE builds in FLUID_CAPACITY_CONSTRAINTS, restricted to what a
% Petri net can declare: a place has no scheduling, no region and no drop rule
% of its own, so there is exactly one gate and it is the same one the exact
% engines apply.
%
% THE GATE IS A LOSS ON THE DEPOSIT. LINE loses the tokens a firing would push
% past a place's capacity -- the NRM's APPLYPLACECAPS clamps the just-deposited
% slots after every firing, which is JMT's and the CTMC's semantics for a
% bounded place (an M/M/1/1 place at rho = 0.5 holds 1/3 of a token, not 1). The
% fluid analogue scales the DEPOSIT leg of every event that adds mass to the
% capped place, and leaves the removal leg alone: the firing still happens, and
% only the mass that does not fit is lost.
%
% CONSERVATION AND LOSS CANNOT BOTH HOLD. A row that binds destroys mass, so any
% conserved quantity supported on the capped coordinates stops being conserved.
% SOLVER_FLUID_PETRI therefore drops those conservation rows for as long as the
% cap is active, which is the honest statement of what the model declares: a
% closed net whose place drops tokens is not closed while the cap binds. In
% practice this only arises where it should -- an open net bounded to make it
% ergodic (the M/M/1/K place) has no conserved quantity to lose in the first
% place, since its arrivals already broke them.
%
% Parameters:
%   sn    - NetworkStruct
%   terms - FLUID_PETRI_TERMS output
%
% Returns:
%   con - struct with fields
%           A     - (ncon x nstate) constraint rows
%           b     - (ncon x 1) right-hand sides
%           label - one description per row
%           cover - (ncon x nstate) logical, the coordinates the row caps, i.e.
%                   the deposits its multiplier throttles
%
% See also SOLVER_FLUID_PETRI, FLUID_CAPACITY_CONSTRAINTS.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

A = zeros(0, terms.nstate);
b = zeros(0,1);
label = {};

for ind = terms.places
    ist = sn.nodeToStation(ind);
    slots = zeros(0,1);
    for k = 1:terms.K
        if terms.pidx(ind,k) > 0
            slots(end+1,1) = terms.pidx(ind,k); %#ok<AGROW>
        end
    end
    if isempty(slots)
        continue
    end
    if ist <= numel(sn.cap) && isfinite(sn.cap(ist))
        row = zeros(1, terms.nstate);
        row(slots) = 1;
        A(end+1,:) = row; %#ok<AGROW>
        b(end+1,1) = sn.cap(ist); %#ok<AGROW>
        label{end+1} = sprintf('capacity %g of place %s', sn.cap(ist), sn.nodenames{ind}); %#ok<AGROW>
    end
    for k = 1:terms.K
        s = terms.pidx(ind,k);
        if s == 0 || ist > size(sn.classcap,1) || ~isfinite(sn.classcap(ist,k))
            continue
        end
        row = zeros(1, terms.nstate);
        row(s) = 1;
        A(end+1,:) = row; %#ok<AGROW>
        b(end+1,1) = sn.classcap(ist,k); %#ok<AGROW>
        label{end+1} = sprintf('class-%d capacity %g of place %s', k, sn.classcap(ist,k), sn.nodenames{ind}); %#ok<AGROW>
    end
end

% TWO ROWS THAT SAY THE SAME THING ARE A SINGULAR NEWTON SYSTEM, not a
% redundancy the least squares absorbs, so an exact duplicate is pruned and the
% tighter bound survives -- the same pruning FLUID_CAPACITY_CONSTRAINTS does.
keep = true(size(A,1),1);
for c = 1:size(A,1)
    if ~keep(c)
        continue
    end
    for d = c+1:size(A,1)
        if keep(d) && isequal(A(c,:), A(d,:))
            if b(d) < b(c)
                keep(c) = false;
                break
            end
            keep(d) = false;
        end
    end
end
A = A(keep,:); b = b(keep); label = label(keep);

con = struct('A', A, 'b', b, 'label', {label}, 'cover', A > 0);
end
