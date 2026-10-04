function [rup, rin, ds, drain] = fluid_capacity_legs(terms, gates, stg, staged, active, x, sg, r, mult, stagedFlow)
% [RUP, RIN, DS, DRAIN] = FLUID_CAPACITY_LEGS(TERMS, GATES, STG, STAGED, ACTIVE, X, SG, R, MULT, STAGEDFLOW)
%
% The two legs of every event under the active caps, and the waiting rooms.
%
% Shared by the steady-state residual and the transient right-hand side so that
% the two solve the SAME model and not two spellings of it. What differs between
% them is only what a STAGED cap's multiplier means:
%
%   stagedFlow = false (steady state) -- a drain RATE theta. The room mass is then
%       pinned by the cap (s = inflow/theta at the fixed point), which is what
%       makes the algebraic system square, and the proportional split across a
%       region's rooms falls out of the common rate.
%   stagedFlow = true (transient) -- the admitted FLOW itself. A rate cannot start
%       the constrained phase at all: at the instant the region fills the room is
%       EMPTY, so theta*s is zero however large theta is, and holding the cap needs
%       a finite admitted flow immediately. The flow is split across the room
%       masses, or across their inflows while the rooms are still empty -- and the
%       two rules AGREE at a fixed point, where s_j is proportional to inflow_j, so
%       the transient and the steady state describe one model.
%
% A held or lost cap composes as a PRODUCT of fractions either way, which is what
% independent blocking gives and what keeps every active cap present in the
% Jacobian: a product has a live derivative in each factor. Several staged caps
% gating one room compose HARMONICALLY, 1/theta = sum 1/theta_f, because the waits
% a job serves in turn add; with one cap that is exactly theta_f.
%
% Parameters:
%   terms      - the event representation from FLUID_MOMENT_TERMS
%   gates      - which cap throttles which event, from FLUID_CAPACITY_GATES
%   stg        - the waiting rooms, from FLUID_CAPACITY_STAGING
%   staged     - CON.STAGED, true where the blocked job waits in a room
%   active     - indices of the caps that currently bind
%   x, sg      - the state and the room populations
%   r          - the nominal rate of every event
%   mult       - one multiplier per ACTIVE cap
%   stagedFlow - false for the steady state, true for the transient
%
% Returns:
%   rup   - the rate each event FIRES at
%   rin   - the rate mass LANDS at
%   ds    - the derivative of each room
%   drain - each room's total outflow
%
% See also SOLVER_FLUID_DAE, FLUID_CAPACITY_GATES, FLUID_CAPACITY_STAGING.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
nevents = numel(terms.rateBase);
nstg = stg.n;
nact = numel(active);
whole = ones(nevents,1);
entry = ones(nevents,1);
invTheta = zeros(nstg,1);
flowOf = zeros(nstg,1);
throttled = false(nstg,1);
roomsOf = cell(nact,1);
multOf = zeros(nact,1);
for k = 1:nact
    c = active(k);
    m = mult(k);
    if staged(c)
        rooms = find(stg.gatedBy(c,:));
        roomsOf{k} = rooms;
        multOf(k) = m;
        throttled(rooms) = true;
        if ~stagedFlow
            if m > GlobalConstants.Zero
                invTheta(rooms) = invTheta(rooms) + 1/m;
            else
                invTheta(rooms) = Inf;
            end
        end
    else
        heldMask = gates.held(c,:)';
        lossMask = gates.loss(c,:)';
        whole(heldMask) = whole(heldMask) * m;
        entry(lossMask) = entry(lossMask) * m;
    end
end

adm_split = false(nevents,1);
if nstg > 0
    live = find(stg.adm);
    for ii = 1:numel(live)
        e = live(ii);
        if throttled(stg.admStage(e))
            adm_split(e) = true;
        end
    end
end
idxSplit = find(adm_split);

% A staged event is NOT suppressed upstream even when a held cap also gates it:
% the job completes upstream service into the waiting room, so a station buffer
% inside a capped region holds the job in THAT ROOM rather than back at a station
% it has already left; the held fraction moves to the room's exit leg below.
% Suppressing both legs would park the same job in two places at once and the
% solve stalls -- a region of 6 with a buffer of 2 inside it stopped at residual
% 1.5 with its cap exceeded.
scale = whole;
scale(adm_split) = 1;
rup = r .* scale;
rin = rup .* entry;

inflow = zeros(nstg,1);
R = zeros(nstg,1);
for ii = 1:numel(idxSplit)
    e = idxSplit(ii);
    j = stg.admStage(e);
    inflow(j) = inflow(j) + rup(e);
    R(j) = R(j) + rup(e);
end

if stagedFlow
    for k = 1:nact
        rooms = roomsOf{k};
        if isempty(rooms)
            continue
        end
        mass = sum(sg(rooms));
        if mass > GlobalConstants.FineTol
            w = sg(rooms) / mass;
        else
            tot = sum(inflow(rooms));
            if tot > GlobalConstants.FineTol
                w = inflow(rooms) / tot;
            else
                w = ones(numel(rooms),1) / numel(rooms);
            end
        end
        flowOf(rooms) = flowOf(rooms) + multOf(k) * w(:);
    end
else
    % a room fed by a multiplier of zero drains at zero: the tightest cap wins in
    % the limit, which is what 1/Inf gives with no special case
    thetaOf = zeros(nstg,1);
    pos = invTheta > 0;
    thetaOf(pos) = 1 ./ invTheta(pos);
    flowOf = thetaOf .* sg;
end
drain = flowOf;

left = zeros(nstg,1);
for ii = 1:numel(idxSplit)
    e = idxSplit(ii);
    j = stg.admStage(e);
    if R(j) > GlobalConstants.FineTol
        q = flowOf(j) * rup(e) / R(j);
    else
        q = 0;
    end
    % the held fraction gates the room's EXIT: what it stops stays in the room,
    % which is what a waiting queue does with it. The lost fraction gates the
    % ARRIVAL: that mass leaves the room and is destroyed, so it is drained but
    % never delivered.
    rin(e) = q * whole(e) * entry(e);
    left(j) = left(j) + q * whole(e);
end
hot = R > GlobalConstants.FineTol;
drain(hot) = left(hot);

ds = inflow - drain;
% a waiting room with no cap above it must be EMPTY, not merely balanced
ds(~throttled) = sg(~throttled);
end
