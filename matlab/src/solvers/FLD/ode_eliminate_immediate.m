function [all_jumps_red, rateBase_red, eventIdx_red, state_map, Emap, absorb] = ...
    ode_eliminate_immediate(all_jumps, rateBase, eventIdx, sn, options) %#ok<INUSL>
% [ALL_JUMPS_RED, RATEBASE_RED, EVENTIDX_RED, STATE_MAP, EMAP, ABSORB] = ...
%     ODE_ELIMINATE_IMMEDIATE(ALL_JUMPS, RATEBASE, EVENTIDX, SN, OPTIONS)
%
% Stochastic complementation of the IMMEDIATE coordinates of the fluid ODE.
%
% A coordinate whose exit rate is GlobalConstants.Immediate (= 1/FineTol = 1e8)
% is not a fast coordinate, it is an INSTANTANEOUS one: the rate is LINE's
% stand-in for infinity, written by SolverLN for the branch of an activity that
% takes no time (an entry called with probability y < 1 carries a second PH
% phase at InfRate entered with probability 1-y). Integrating it numerically is
% meaningless work -- the mode relaxes 1e8 times faster than anything else in
% the model, so every integrator either crawls or steps over it -- and it is
% what made a two-station LN layer take 145 s in the JAR and 70 s in MATLAB for
% an answer identical, to every digit, to the one the reduced system gives in
% 0.35 s. See _kb/06-solver-catalog.md.
%
% THE REDUCTION IS EXACT, not an approximation. It is the ODE twin of
% CTMC_STOCHCOMP: the instantaneous coordinates F are absorbed into the timed
% ones S by the absorption probabilities of the embedded jump chain restricted
% to F, so the flow that would enter F is routed straight to where F would have
% sent it. On the layer above, the FCFS station's immediate phase folds into the
% delay -- Delay(0.3767) -> Queue(0.6105) with only 0.739012 of the departures
% entering the queue -- whose fluid limit is Q = [2.192891786795962,
% 1.807108213330636], the number the stiff integration spends 145 s reaching.
%
% WHY THIS IS A STRUCTURAL COMPOSITION AND NOT A GENERATOR ROUND TRIP. Every
% event of ODE_JUMPS_NEW is a single -1 at EVENTIDX and a single +1 at its
% destination, so a path through F composes to one event, -1 at the original
% source and +1 at the absorbing coordinate, that keeps the original source's
% GATING. Rebuilding the events from a reduced generator instead (the shape this
% function had before) loses that identity, and with it the event ORDER that
% FLUID_MOMENT_TERMS reads throughputs off -- which is the whole reason the
% moment-closure methods used to refuse the reduction outright. EMAP carries the
% identity across instead: EMAP(e,o) is the expected number of times the ORIGINAL
% event o fires per firing of the reduced event e, so a caller maps any
% per-event quantity with NEWATTR = EMAP * OLDATTR and gets an exact rate
% accounting. It is the identity when nothing is eliminated.
%
% A COMPOSED EVENT CAN BE A DEPARTURE AT TWO STATIONS AT ONCE, which is why a
% single evIsDeparture flag cannot survive the composition: a job that leaves
% the delay, passes through the queue's immediate phase and returns has
% completed at BOTH, and both throughputs must count it. EMAP gives it a row
% with weight on both original events, and the null jump it composes to (-1 and
% +1 on the same coordinate) correctly contributes nothing to the drift and
% nothing to the diffusion D*diag(r)*D'.
%
% Parameters:
%   all_jumps - [nstate x nevents] jump matrix of ODE_JUMPS_NEW
%   rateBase  - [nevents x 1] fixed part of each event rate
%   eventIdx  - [nevents x 1] source coordinate of each event, which is also
%               the coordinate whose occupancy gates it
%   sn        - NetworkStruct, unused, kept for the caller's signature
%   options   - solver options; options.config.immediate_tol overrides the
%               detection threshold
%
% Returns:
%   all_jumps_red - [nstate x nevents_red] jumps, in the ORIGINAL coordinate
%                   layout, with the eliminated rows identically zero
%   rateBase_red  - [nevents_red x 1] reduced rates
%   eventIdx_red  - [nevents_red x 1] source coordinates, original indexing
%   state_map     - the timed coordinates that survive, in increasing order
%   Emap          - [nevents_red x nevents] expected firings of each original
%                   event per firing of each reduced one
%   absorb        - [nstate x nstate] projector taking an initial condition to
%                   the reduced coordinates: identity on the timed rows, the
%                   absorption distribution on the immediate ones. Mass parked
%                   on an eliminated coordinate would otherwise be frozen there
%                   for the whole integration, because nothing moves it any more
%
% See also CTMC_STOCHCOMP, ODE_JUMPS_NEW, ODE_RATE_BASE, FLUID_MOMENT_TERMS.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

nstate = size(all_jumps, 1);
nevents = numel(rateBase);

% A rate at or above the threshold is the InfRate sentinel, not a fast rate the
% user wrote: the default sits just under GlobalConstants.Immediate so that only
% the sentinel qualifies. options.config.immediate_tol lowers it deliberately.
imm_tol = GlobalConstants.Immediate * (1 - 1e-2);
if isfield(options,'config') && isfield(options.config,'immediate_tol') ...
        && ~isempty(options.config.immediate_tol)
    imm_tol = options.config.immediate_tol;
end

identity_out = @() deal(all_jumps, rateBase(:), eventIdx(:), (1:nstate)', ...
    speye(nevents), speye(nstate));

imm_idx = find(rateBase(:) >= imm_tol);
if isempty(imm_idx)
    [all_jumps_red, rateBase_red, eventIdx_red, state_map, Emap, absorb] = identity_out();
    return
end

% The immediate coordinates are the SOURCES of the immediate events: it is the
% coordinate that empties instantaneously, not the event.
isImm = false(nstate,1);
isImm(eventIdx(imm_idx)) = true;

% Destination of each event. Every event of ODE_JUMPS_NEW is -1 at its source
% and +1 at one destination; a departure that re-enters its own coordinate
% cancels to an all-zero column, whose destination is that same coordinate.
dst = zeros(nevents,1);
src = eventIdx(:);
for e = 1:nevents
    d = find(all_jumps(:,e) > 0, 1);
    if isempty(d)
        dst(e) = src(e);
    else
        dst(e) = d;
    end
end

% A coordinate with no outflow at all cannot be complemented away, and one whose
% outflow is entirely a self-loop would make the fundamental matrix singular.
% Both are dropped from F rather than guessed at.
for f = find(isImm)'
    out_f = find(src == f);
    tot = sum(rateBase(out_f));
    if tot <= 0 || all(dst(out_f) == f)
        isImm(f) = false;
    end
end
if ~any(isImm)
    [all_jumps_red, rateBase_red, eventIdx_red, state_map, Emap, absorb] = identity_out();
    return
end

Fidx = find(isImm);
Sidx = find(~isImm);
nF = numel(Fidx);
if isempty(Sidx)
    % Nothing timed is left to absorb into; the reduced system would be empty.
    [all_jumps_red, rateBase_red, eventIdx_red, state_map, Emap, absorb] = identity_out();
    return
end
posF = zeros(nstate,1); posF(Fidx) = 1:nF;
posS = zeros(nstate,1); posS(Sidx) = 1:numel(Sidx);

% Branching of the embedded jump chain out of each immediate coordinate. The
% probabilities are the rate shares, so a coordinate carrying both an immediate
% and an ordinary exit gives the ordinary one its (vanishing) share rather than
% being special-cased.
PFF = sparse(nF, nF);
PFS = sparse(nF, numel(Sidx));
Cnt = sparse(nF, nevents); % per-entry firing probability of each original event
for a = 1:nF
    f = Fidx(a);
    out_f = find(src == f);
    tot = sum(rateBase(out_f));
    for e = out_f(:)'
        p = rateBase(e) / tot;
        Cnt(a,e) = Cnt(a,e) + p; %#ok<SPRIX>
        if isImm(dst(e))
            PFF(a, posF(dst(e))) = PFF(a, posF(dst(e))) + p; %#ok<SPRIX>
        else
            PFS(a, posS(dst(e))) = PFS(a, posS(dst(e))) + p; %#ok<SPRIX>
        end
    end
end

% Fundamental matrix of the instantaneous chain. (I-PFF) is invertible whenever
% every immediate coordinate reaches a timed one, which the drop above ensures
% for the self-loop case; a genuinely absorbing cycle of immediate coordinates
% is a modelling error and is left to the original system rather than solved.
Ifm = speye(nF) - PFF;
if rcond(full(Ifm)) < GlobalConstants.FineTol
    line_warning(mfilename, ['the immediate coordinates form a closed cycle, so they have no ' ...
        'absorption distribution; integrating the unreduced system instead']);
    [all_jumps_red, rateBase_red, eventIdx_red, state_map, Emap, absorb] = identity_out();
    return
end
Nfm = Ifm \ speye(nF);   % expected visits to each immediate coordinate
Aabs = Nfm * PFS;        % absorption distribution over the timed coordinates
ExpCnt = Nfm * Cnt;      % expected firings of each original event, per entry

% Compose the event list. An event sourced in F is dropped: its flow is already
% carried by whichever event feeds F.
keep = ~isImm(src);
nkeep = nnz(keep);
jumps_new = cell(1, nkeep + nnz(keep & isImm(dst)) * numel(Sidx));
rate_new = zeros(1, numel(jumps_new));
evidx_new = zeros(1, numel(jumps_new));
emap_i = []; emap_j = []; emap_v = [];
cnt = 0;
for e = find(keep)'
    if ~isImm(dst(e))
        cnt = cnt + 1;
        jumps_new{cnt} = all_jumps(:,e);
        rate_new(cnt) = rateBase(e);
        evidx_new(cnt) = src(e);
        emap_i(end+1) = cnt; emap_j(end+1) = e; emap_v(end+1) = 1; %#ok<AGROW>
        continue
    end
    % The event feeds an immediate coordinate: replace it by one event per
    % absorbing destination, keeping the original source and so the original
    % gating, since the rate of the composed flow IS the rate of the inflow.
    a = posF(dst(e));
    row = ExpCnt(a,:);
    for b = find(Aabs(a,:) > 0)
        s = Sidx(b);
        cnt = cnt + 1;
        jump = zeros(nstate,1);
        jump(src(e)) = jump(src(e)) - 1;
        jump(s) = jump(s) + 1;
        jumps_new{cnt} = jump;
        rate_new(cnt) = rateBase(e) * Aabs(a,b);
        evidx_new(cnt) = src(e);
        % Weighting every absorbing branch by the SAME unconditional expected
        % counts is what makes the rate accounting exact: the branch rates sum
        % back to rateBase(e), so the mapped total is rateBase(e)*row.
        emap_i(end+1) = cnt; emap_j(end+1) = e; emap_v(end+1) = 1; %#ok<AGROW>
        nz = find(row);
        for o = nz
            emap_i(end+1) = cnt; emap_j(end+1) = o; emap_v(end+1) = row(o); %#ok<AGROW>
        end
    end
end

all_jumps_red = zeros(nstate, cnt);
for e = 1:cnt
    all_jumps_red(:,e) = jumps_new{e};
end
rateBase_red = rate_new(1:cnt)';
eventIdx_red = evidx_new(1:cnt)';
Emap = sparse(emap_i, emap_j, emap_v, cnt, nevents);
state_map = Sidx;

absorb = speye(nstate);
absorb(Fidx,:) = 0;
for a = 1:nF
    for b = find(Aabs(a,:) > 0)
        absorb(Fidx(a), Sidx(b)) = Aabs(a,b);
    end
end

if isfield(options,'verbose') && options.verbose >= VerboseLevel.DEBUG
    line_printf(sprintf('\nEliminated %d immediate coordinates of %d, %d events of %d\n', ...
        nF, nstate, nevents - cnt, nevents));
end
end % ode_eliminate_immediate
