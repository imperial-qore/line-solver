function pinfo = pollingInfo(sn, ind)
% PINFO = POLLINGINFO(SN, IND)
%
% Derived description of the polling controller at the station of node IND.
% Returns [] when the node is not a polling station.
%
% The polling controller is stored in the trailing local-variable block of the
% station state, whose width is sn.nvars(ind,2*R+1). The block holds, in order,
% the columns [pos, swk, ctr]; each is materialized only when the discipline
% actually needs it (see the widths below), so that a polling station never
% carries state that its dynamics cannot distinguish.
%
%   pos  index of the buffer the server is currently at (serving) or heading to
%        (switching). It is materialized only when at least one switchover is
%        non-immediate: while a job is in service pos always equals the class of
%        that job, and while the station is empty and every switchover is
%        immediate the server position is unobservable (see PARKED below).
%   swk  0 when the server sits at pos, otherwise the phase of the switchover
%        PH into buffer pos. Materialized only when some switchover is
%        non-immediate.
%   ctr  the visit budget, materialized for every discipline except EXHAUSTIVE:
%          GATED        jobs of class pos admitted at the polling instant that
%                       have not completed yet (the job in service counts as
%                       one of them), so the visit ends when ctr reaches 0;
%          KLIMITED     services still permitted in this visit, the one in
%                       progress included;
%          DECREMENTING the target class-pos population: the visit ends once
%                       the population has dropped to ctr, i.e. one below the
%                       level found at the polling instant (semi-exhaustive).
%
% The three tangible controller configurations are therefore
%   SERVING(p)   swk=0, one class-p job in the service facility;
%   SWITCHING(p) swk>0, service facility empty;
%   PARKED       swk=0, service facility empty, station empty. Reachable only
%                when every switchover is immediate, in which case a server that
%                completes a full lap without finding work would otherwise cycle
%                in zero time forever. pos is then unobservable and canonical.
%
% Immediate() switchovers are NOT represented as states: Immediate.getProcess
% returns an exponential of rate GlobalConstants.Immediate, and taking that
% literally would put a ~1e8 rate in the generator (stiff, and a spurious state
% per buffer). They are instead folded into the enclosing transition by
% State.pollingNext, which walks the cyclic order until a tangible state.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

pinfo = [];
if ~sn.isstation(ind)
    return
end
ist = sn.nodeToStation(ind);
if sn.sched(ist) ~= SchedStrategy.POLLING
    return
end

% The description is a pure function of the model, but it is read once per
% state per synchronization when the generator is built, so it is memoized on
% nodeparam by refreshLocalVars rather than rebuilt (map_pie and all) per call.
if ~isempty(sn.nodeparam{ind}) && iscell(sn.nodeparam{ind}) ...
        && ~isempty(sn.nodeparam{ind}{1}) && isfield(sn.nodeparam{ind}{1},'pollinfo')
    pinfo = sn.nodeparam{ind}{1}.pollinfo;
    return
end

R = sn.nclasses;
pinfo.ptype = PollingType.EXHAUSTIVE;
pinfo.pk = 1;
pinfo.hasSw = false(1,R);
pinfo.Ksw = zeros(1,R);
pinfo.swD0 = cell(1,R);
pinfo.swD1 = cell(1,R);
pinfo.swpie = cell(1,R);

% Buffers the server actually visits. A class disabled at this station can
% never hold a job, and State.afterEvent short-circuits every event carrying
% it (K(class)==0), so it cannot be given a switchover to walk through: such a
% buffer is dropped from the cyclic order rather than polled forever.
pinfo.polled = true(1,R);
for r=1:R
    if isempty(sn.proc{ist}{r}) || any(any(isnan(sn.proc{ist}{r}{1})))
        pinfo.polled(r) = false;
    end
end
if ~any(pinfo.polled)
    line_error(mfilename, sprintf('Polling station %d has no class with an enabled service.', ist));
end

% see _kb/04-networkstruct.md (refreshGlobalSync/refreshSync) for rationale
swof = zeros(1,R);
polledlist = find(pinfo.polled);
for jj=1:length(polledlist)
    q = polledlist(jj);
    prevq = polledlist(mod(jj-2, length(polledlist)) + 1);
    swof(q) = prevq;
end

np = sn.nodeparam{ind};
for r=1:R
    if isempty(np) || length(np) < r || isempty(np{r})
        continue
    end
    if isfield(np{r},'pollingType') && ~isempty(np{r}.pollingType)
        % setPollingType writes the same discipline into every class buffer
        pinfo.ptype = PollingType.toId(np{r}.pollingType);
        if isfield(np{r},'pollingPar') && ~isempty(np{r}.pollingPar)
            pinfo.pk = round(np{r}.pollingPar(1));
        end
    end
end

% Switchover of the leg entering each polled buffer q, read off the buffer the
% server leaves to get there.
for q = polledlist
    src = swof(q);
    if isempty(np) || length(np) < src || isempty(np{src}) ...
            || ~isfield(np{src},'switchoverTime') || isempty(np{src}.switchoverTime)
        continue
    end
    procid = np{src}.switchoverProcId;
    proc = np{src}.switchoverTime;
    if iscell(proc) && ~isempty(proc) && iscell(proc{1})
        % from-to matrix form: polling only uses the per-buffer vector
        proc = proc{1};
        procid = procid(1);
    end
    if procid(1) == ProcessType.IMMEDIATE
        continue % zero-time switchover: folded, never a state
    end
    pinfo.hasSw(q) = true;
    pinfo.Ksw(q) = length(proc{1});
    pinfo.swD0{q} = proc{1};
    pinfo.swD1{q} = proc{2};
    pinfo.swpie{q} = map_pie(proc);
end

% Column widths. pos and swk exist only to encode SWITCHING(p); with no
% non-immediate switchover the server is either serving (pos = class in
% service) or parked (pos unobservable), so neither column carries
% information. ctr exists for every discipline that bounds a visit.
pinfo.wpos = double(any(pinfo.hasSw));
pinfo.wswk = double(any(pinfo.hasSw));
pinfo.wctr = double(pinfo.ptype ~= PollingType.EXHAUSTIVE);
pinfo.width = pinfo.wpos + pinfo.wswk + pinfo.wctr;

% Offset of the polling block inside space_var. The block is the 2*R+1 column
% of nvars, hence it trails the per-class Markov-modulation slots (1..R) and
% the routing slots (R+1..2*R).
pinfo.off = sum(sn.nvars(ind,1:2*R));

% Column indices inside space_var, 0 when the column is not materialized.
pinfo.ipos = 0;
pinfo.iswk = 0;
pinfo.ictr = 0;
c = pinfo.off;
if pinfo.wpos, c = c + 1; pinfo.ipos = c; end
if pinfo.wswk, c = c + 1; pinfo.iswk = c; end
if pinfo.wctr, c = c + 1; pinfo.ictr = c; end
end
