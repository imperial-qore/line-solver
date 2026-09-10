function cons = fluid_petri_conservation(terms)
% CONS = FLUID_PETRI_CONSERVATION(TERMS)
%
% The conserved quantities of a fluid Petri net, as equations.
%
% Everywhere else in SolverFLD conservation is a CONSEQUENCE of the drift and
% therefore holds only to integrator tolerance. The DAE form states it, and for
% a Petri net the right statement is not "one row per closed chain" -- a net has
% no chains -- but the left null space of the jump matrix:
%
%   u' D = 0   =>   u' x is constant along every trajectory
%
% On the marking coordinates those u are exactly the net's P-INVARIANTS (an
% S-invariant is a non-negative left null vector of the incidence matrix, which
% is what D restricted to the places is). On a mode's phase block the all-ones
% vector is one of them, and it is the statement that the phase coordinates are
% a DISTRIBUTION: sum_h z(j,h) = 1. Both come out of the same null space, so the
% phase normalisation needs no separate row and no separate mechanism.
%
% AN OPEN NET LOSES THE ROWS ITS ARRIVALS BREAK, automatically and for the right
% reason: the arrival columns are part of D, so a u that an arrival moves is not
% in the null space and never becomes a constraint.
%
% THE BASIS IS RATIONAL, not orthonormal. D is integral (arc multiplicities and
% unit phase moves), so NULL(.,'r') returns exact rational rows through RREF,
% which keeps each constraint readable as a conservation statement about named
% places instead of an arbitrary orthogonal mixture of them.
%
% Parameters:
%   terms - FLUID_PETRI_TERMS output
%
% Returns:
%   cons - struct with fields
%            C     - (ncons x nstate) conserved directions
%            N     - (ncons x 1) their value at the initial marking
%            leak  - max|C*D|, zero by construction; a positive value means the
%                    basis and the jump matrix disagree
%            label - one description per row
%
% See also SOLVER_FLUID_PETRI, FLUID_PETRI_TERMS, SPN_SINVARIANTS.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

D = terms.D;
if isempty(D) || size(D,2) == 0
    cons = struct('C', zeros(0,terms.nstate), 'N', zeros(0,1), 'leak', 0, 'label', {{}});
    return
end

Z = null(D.', 'r');
if isempty(Z)
    Z = zeros(terms.nstate, 0);
end
C = Z.';

% Clear the numerical dust RREF leaves on a rational basis and normalise each
% row so that its smallest nonzero entry is one, which makes an invariant read
% as "these places hold this many tokens" rather than as a scaled copy of it.
C(abs(C) < 1e-12) = 0;
for c = 1:size(C,1)
    nz = abs(C(c,:));
    nz = nz(nz > 0);
    if ~isempty(nz)
        C(c,:) = C(c,:) / min(nz);
    end
end
C = C(any(C ~= 0, 2), :);

N = C * terms.x0;
leak = 0;
if ~isempty(C)
    leak = max(max(abs(C * D)));
end

label = cell(size(C,1),1);
for c = 1:size(C,1)
    parts = {};
    for s = find(C(c,:) ~= 0)
        parts{end+1} = i_coordname(terms, s, C(c,s)); %#ok<AGROW>
    end
    label{c} = strjoin(parts, ' + ');
end

cons = struct('C', C, 'N', N, 'leak', leak, 'label', {label});
end

% -------------------------------------------------------------------------
function s = i_coordname(terms, idx, w)
if idx <= terms.nm
    nm = sprintf('%s(class %d)', terms.namesNode{terms.coordNode(idx)}, terms.coordClass(idx));
else
    nm = 'phase';
    for j = 1:numel(terms.modes)
        if ~isempty(terms.modes(j).zblk) && any(terms.modes(j).zblk == idx)
            nm = sprintf('%s phase %d', terms.modes(j).label, find(terms.modes(j).zblk == idx));
            break
        end
    end
end
if w == 1
    s = nm;
else
    s = sprintf('%g*%s', w, nm);
end
end
