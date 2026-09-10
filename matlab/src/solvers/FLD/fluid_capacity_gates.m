function gates = fluid_capacity_gates(sn, terms, con)
% GATES = FLUID_CAPACITY_GATES(SN, TERMS, CON)
%
% Per cap and per event: is this event an admission the cap throttles, and WHERE
% DOES THE STOPPED MASS GO.
%
% An admission is an event that pushes the constrained quantity UP, read off
% A*D rather than off the topology so that a per-class or memory-weighted row
% picks out its own admissions with no extra code.
%
% THE THREE ANSWERS ARE THE MODEL, AND THEY ARE NOT INTERCHANGEABLE. LINE's own
% semantics decides which one a cap gets, and the choice is visible in the
% answer -- a job held upstream is still counted at that station, a lost job is
% counted nowhere and breaks flow balance across the cap on purpose, a staged job
% is counted in neither and is reported as blocked mass:
%
%   STAGED  a finite capacity region under a waiting queue (CON.STAGED). The job
%           COMPLETES upstream service and waits outside the region, which is
%           what JMT and LDES simulate; the upstream station empties exactly as
%           it would with no region. See FLUID_CAPACITY_STAGING.
%   HELD    a station buffer reached by a CLOSED class. State.arrivalIsLost
%           refuses to lose a closed job -- population conservation is a defining
%           invariant -- and returns an empty successor instead, which DISABLES
%           the upstream departure until room frees. The job is therefore still
%           at the upstream station, in service as far as that station's own
%           metrics are concerned, so the fluid analogue scales the WHOLE event:
%           both the removal upstream and the arrival.
%   LOSS    a station buffer reached by an OPEN class. The same predicate loses
%           it: the external stream is memoryless, so a job that finds the buffer
%           full simply never enters, the arrival event still fires (ArvR counts
%           the offered job) and only the carried flow is admitted. The fluid
%           analogue scales the ARRIVAL leg alone and destroys the difference.
%
% Parameters:
%   sn    - NetworkStruct, read for njobs (open vs closed decides held vs loss)
%   terms - the event representation from FLUID_MOMENT_TERMS
%   con   - the constraint set from FLUID_CAPACITY_CONSTRAINTS
%
% Returns:
%   gates - struct with fields
%             gate (ncon x nevents) does this row throttle this event
%             held (ncon x nevents) scale the WHOLE event: the job stays upstream
%             loss (ncon x nevents) scale the ARRIVAL leg: the rest is destroyed
%             Dn, DnExt, Dp (nstate x nevents) the jump matrix split in three, so
%                    that an admission can remove mass upstream at one rate and
%                    deliver it at another. DN carries the removal at real
%                    stations, DNEXT the removal from the EXT source pool -- which
%                    a LOST arrival must be returned to, since that coordinate is
%                    a normalisation and its drift row a real equation -- and DP
%                    the arrival. Empty where no cap exists, which keeps the
%                    uncapped drift a single matrix product.
%
% See also SOLVER_FLUID_DAE, FLUID_CAPACITY_CONSTRAINTS.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

nevents = numel(terms.rateBase);
ncon = numel(con.b);
gates = struct('gate', false(ncon,nevents), 'held', false(ncon,nevents), ...
    'loss', false(ncon,nevents), 'Dn', [], 'DnExt', [], 'Dp', []);
if ncon == 0
    return
end
tol = sqrt(GlobalConstants.Zero);
gates.Dn = min(terms.D, 0);
gates.Dp = max(terms.D, 0);
gates.DnExt = zeros(size(gates.Dn));
for i = 1:terms.M
    if terms.isExt(i)
        blk = terms.stationBlock{i};
        if ~isempty(blk)
            gates.DnExt(blk,:) = gates.Dn(blk,:);
            gates.Dn(blk,:) = 0;
        end
    end
end

for c = 1:ncon
    delta = con.A(c,:) * terms.D;
    for e = 1:nevents
        if delta(e) <= tol
            continue
        end
        gates.gate(c,e) = true;
        if con.staged(c)
            continue
        end
        % the class is read off the coordinate the mass LANDS on, inside the
        % capped station: a class switch on entry would otherwise ask the class
        % the job is leaving behind whether it may be lost
        landing = find(gates.Dp(:,e) > tol & con.A(c,:)' > 0, 1);
        isopen = false;
        if ~isempty(landing)
            k = con.coordClass(landing);
            if k >= 1 && numel(sn.njobs) >= k
                isopen = ~isfinite(sn.njobs(k));
            end
        end
        gates.loss(c,e) = isopen;
        gates.held(c,e) = ~isopen;
    end
    if ~any(gates.gate(c,:))
        line_error(mfilename, sprintf(['No event increases %s, so the cap can never be approached and ' ...
            'there is no admission flow for the constraint to throttle. This is a malformed limit ' ...
            'rather than a solvable one.'], con.label{c}));
    end
end
end
