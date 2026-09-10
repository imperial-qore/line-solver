function stg = fluid_capacity_staging(terms, con)
% STG = FLUID_CAPACITY_STAGING(TERMS, CON)
%
% The waiting room outside a capped region, as fluid coordinates.
%
% WHY THE CONSTRAINT ALONE IS NOT ENOUGH, FOR A REGION. Throttling the admission
% events of a region does hold its population at the cap, but it holds it by
% slowing the UPSTREAM STATION'S COMPLETIONS -- an admission event is that
% station finishing a job -- so the blocked mass piles up at a station it has
% already finished being served by. Where that station is a delay the error is
% visible as a broken Little's law: measured against LDES the delay reported 12
% jobs at a throughput of 1.6 and a think time of 1, a factor of 7.5 out, while
% the region population and the throughput were both right. A waiting queue means
% the job COMPLETES upstream service and then waits; it is somewhere else, and
% the model needs somewhere else to put it.
%
% A STATION BUFFER GETS NO ROOM, and that is not an omission. There the job
% genuinely does stay where it was: LINE disables the upstream departure
% (State.arrivalIsLost) rather than moving the job out, so the blocked mass is
% still at the upstream station and still counted there. Giving it a room would
% move mass the reference keeps in place. See FLUID_CAPACITY_GATES.
%
% So each capped region gains one coordinate per class, and every admission into
% it is split in two:
%
%   upstream -> staging    at the nominal rate, untouched, so the upstream
%                          station empties exactly as it would with no region
%   staging  -> region     at theta * s, the throttled leg, with theta the
%                          algebraic unknown the constraint pins
%
% The two jumps sum to the original one, so nothing about the event set changes
% except where the mass rests in between.
%
% ONE ROOM PER (REGION, CLASS) AND ONE MULTIPLIER PER ROW, which is what lets two
% caps of one region bind at once -- the case a single per-region throttle had to
% refuse. A room gated by several active rows drains at the HARMONIC composition
% of their rates, 1/theta = sum 1/theta_f, because the waits a job serves in turn
% add; with one row that is exactly theta_f, so the single-cap answer is
% unchanged. GATEDBY is what records which rows gate which room.
%
% Parameters:
%   terms - event representation from FLUID_MOMENT_TERMS
%   con   - constraint set from FLUID_CAPACITY_CONSTRAINTS, carrying per-region
%           membership and which rows stage
%
% Returns:
%   stg - struct with fields
%           n         number of staging coordinates
%           region    (n x 1) region each belongs to
%           class     (n x 1) class each carries
%           idx       (nregions x K) coordinate index, 0 where there is none
%           adm       (nevents x 1) is this event an admission into a region
%           admRegion (nevents x 1) which region, 0 if not an admission
%           admStage  (nevents x 1) which staging coordinate, 0 if not
%           gatedBy   (ncon x n) which rows gate which room
%
% See also FLUID_CAPACITY_CONSTRAINTS, FLUID_CAPACITY_EXTEND, SOLVER_FLUID_DAE.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = terms.K;
nevents = numel(terms.rateBase);
ncon = numel(con.b);
tol = sqrt(GlobalConstants.Zero);

stg = struct('n', 0, 'region', zeros(0,1), 'class', zeros(0,1), ...
    'idx', zeros(max(con.nregions,1),K), 'adm', false(nevents,1), 'admRegion', zeros(nevents,1), ...
    'admStage', zeros(nevents,1), 'gatedBy', false(ncon,0));

if con.nregions == 0 || ncon == 0 || ~any(con.staged)
    return
end

F = con.nregions;
D = terms.D;
Dp = max(D, 0);
stg.idx = zeros(F, K);

stagedRegions = unique(con.region(con.staged));
stagedRegions = stagedRegions(stagedRegions > 0);

region = zeros(0,1); klass = zeros(0,1);
for f = 1:F
    if ~any(stagedRegions == f)
        continue
    end
    memberRow = con.member(f,:);
    if ~any(memberRow)
        continue
    end
    % net change of this region's population per event: positive means the event
    % brings mass in from outside, which is what the waiting queue feeds
    delta = memberRow * D;
    for e = 1:nevents
        if delta(e) <= tol
            continue
        end
        % the class is read off the coordinate the mass lands on, INSIDE the
        % region, not off the upstream one: a class switch on entry would
        % otherwise stage the job under the class it is leaving behind
        landing = find(Dp(:,e) > tol & memberRow(:));
        if isempty(landing)
            continue
        end
        c = con.coordClass(landing(1));
        if c < 1 || c > K
            continue
        end
        if stg.idx(f,c) == 0
            region(end+1,1) = f; %#ok<AGROW>
            klass(end+1,1) = c;  %#ok<AGROW>
            stg.idx(f,c) = numel(region);
        end
        stg.adm(e) = true;
        stg.admRegion(e) = f;
        stg.admStage(e) = stg.idx(f,c);
    end
end

stg.n = numel(region);
stg.region = region;
stg.class = klass;

% WHICH ROWS GATE WHICH ROOM. One drain rate per room and one equality per row,
% and the two are not in bijection: a region-global cap gates every room of its
% region, a per-class cap only the room of its class.
stg.gatedBy = false(ncon, stg.n);
for c = 1:ncon
    if ~con.staged(c) || con.region(c) < 1
        continue
    end
    for j = 1:stg.n
        if stg.region(j) ~= con.region(c)
            continue
        end
        sel = con.member(con.region(c),:) & (con.coordClass == stg.class(j));
        if any(sel) && max(con.A(c,sel)) > 0
            stg.gatedBy(c,j) = true;
        end
    end
end
end
