function checkServiceReachable(self)
% CHECKSERVICEREACHABLE  Refuse a class ROUTED TO a station that cannot serve it.
%
% SANITIZE disables the OUTGOING routing of a class a station cannot serve,
% which is what keeps it out of that station's visit ratios -- but nothing
% stopped the class being routed IN, and a class that arrives where it cannot be
% served is a flow sink: it enters and never leaves. The station-level guard next
% to it cannot see this, because it asks whether the station serves ANY class,
% not whether it serves the classes that reach it.
%
% One such model gave three different wrong answers, none flagged, on a closed
% cycle D <-> Q whose class C2 has no service at Q:
%
%     MVA   Q/C2 ArvR 1, Tput 0        flow not conserved
%     CTMC  drops class C2 entirely
%     SSA   D/C2 QLen 2e-06, Q rows absent
%
% AN ABSENT SLOT AND AN EXPLICIT Disabled() ARE TREATED ALIKE. sanitize fills an
% absent slot with Disabled(), so the two reach every solver as the identical
% struct and produce identical wrong numbers; the spelling cannot decide. What
% decides is whether anything ROUTES THE CLASS IN, which is what this tests, so
% the class-switching idiom -- each queue serving one class and marking the rest
% Disabled -- has no incoming flow there and is untouched.
%
% READS sn.rtnodes, WALKED FORWARD FROM THE FEED POINTS -- not sn.nodevisits,
% which this guard read until 2026-09-02 and which 94d5570f3 had made blind to
% the very case it exists for. That commit extended the `served` mask of
% SN_REFRESH_VISITS from the station chain to the NODE chain, and it had to: on a
% materialised LQN replica the unserved states close into a spurious cycle. But
% the mask zeroes exactly the (station, class) cell a flow sink shows up in. On
% Source -> Q -> Sink with class B unservable at Q, B's chain went from Q = 1 to
% Q = 0 and the guard fell silent, while the Sink still read 1 -- flow arriving
% downstream of a node it never visited. A MASKED VISIT VECTOR CANNOT ANSWER THIS
% QUESTION, because the mask IS the answer being looked for. Do not route this
% guard back through nodevisits or visits; both carry that mask.
%
% rtnodes on its own over-approximates -- it says where a class WOULD go if one
% existed -- and the WALK is what removes the slack. It starts only at (Source,
% class) pairs whose arrival is not Disabled, and at the reference station of
% each closed class with a positive population, so the Disabled-arrival row a
% class-switching Source carries is never entered. Three rules keep it honest:
% a Sink is ABSORBING (rtnodes wraps it back to the Source to close the kernel,
% and following that wrap re-enters every Source row, including the Disabled ones
% the seeding just excluded); a Source is expanded ONLY AS A SEED, for the same
% reason; and an unservable (station, class) is REACHED BUT NOT EXPANDED, since
% nothing leaves it -- that is the whole complaint -- so nothing downstream of it
% is evidence of anything.
%
% This SUBSUMES the fed-chain precondition the guard used to carry separately: a
% chain no job can enter has no seed, so its rows are never walked at all. That is
% strictly finer than the per-chain test it replaces, which admitted every class
% of a chain any one of whose classes was fed.
%
% A SYNCHRONOUS REPLY IS NOT SERVED BY THE STATION IT RETURNS TO: it releases the
% server that station held across the call, which is the whole content of
% setSyncReply. Its Disabled service there is the marker of the feature, not a
% flow sink, so the (station, reply class) pairs sn.replyblock marks are exempt.

if ~self.getChecks()
    return
end
sn = self.sn;
if isempty(sn) || isempty(sn.nclasses) || sn.nclasses < 1
    return
end
% Same exemption as the sanitize checks: a station of a cache, Petri-net or
% fork-join model legitimately carries no per-class service.
for ind = 1:length(self.nodes)
    nd = self.nodes{ind};
    if isa(nd,'Cache') || isa(nd,'Place') || isa(nd,'Transition') || isa(nd,'Fork') || isa(nd,'Join')
        return
    end
end
reached = local_reached(self, sn);
if isempty(reached)
    return
end
for ind = 1:length(self.nodes)
    nd = self.nodes{ind};
    if ~(isa(nd,'Queue') || isa(nd,'Delay'))
        continue
    end
    % A heterogeneous pool carries its service on the server types, so the
    % per-class process is legitimately disabled there.
    %
    % ISEMPTY IS NOT THE TEST. heteroServiceDistributions is a dictionary, and a
    % dictionary is a 1x1 object whose ISEMPTY is FALSE even with zero entries --
    % so `~isempty(...)` alone is true for EVERY queue and skipped every station
    % this guard exists to check. SANITIZE.m spells the same test with both
    % halves for exactly this reason; keep NUMENTRIES.
    if isprop(nd,'heteroServiceDistributions') && ~isempty(nd.heteroServiceDistributions) ...
            && numEntries(nd.heteroServiceDistributions) > 0
        continue
    end
    if ind > size(reached,1)
        continue
    end
    for r = 1:sn.nclasses
        if r > size(reached,2) || ~reached(ind,r)
            continue
        end
        svc = [];
        if r <= length(nd.server.serviceProcess) && ~isempty(nd.server.serviceProcess{r})
            svc = nd.server.serviceProcess{r}{end};
        end
        if ~isempty(svc) && ~isa(svc,'Disabled')
            continue
        end
        % A SYNCHRONOUS REPLY is not served by the station it returns to: it
        % releases the server that station held across the call, which is the
        % whole content of setSyncReply.
        if holdsReplyFor(sn, ind, r)
            continue
        end
        if isa(nd,'Delay')
            kind = 'Delay';
        else
            kind = 'Queue';
        end
        line_error(mfilename, sprintf(['%s ''%s'' has no service configured for job class ''%s'', but the class is routed to it. ' ...
            'Jobs would arrive and never leave. Call setService() for that class, route it elsewhere, ' ...
            'or disable this check with model.setChecks(false).'], ...
            kind, nd.getName(), self.classes{r}.getName()));
    end
end
end

function tf = holdsReplyFor(sn, ind, r)
% HOLDSREPLYFOR  True when node IND holds a server across a synchronous call
% whose reply class is R, i.e. R returns there to RELEASE a server rather than to
% be served by one. Both indices are 1-based.
%
% sn.syncreply is indexed by the CALLING class and holds the 0-BASED reply class
% (-1 where no reply is expected, which can never match a 1-based R);
% sn.replyblock marks the (node, calling class) pairs that hold a server.
tf = false;
if ~isfield(sn,'replyblock') || ~isfield(sn,'syncreply') || ...
        isempty(sn.replyblock) || isempty(sn.syncreply) || ind > size(sn.replyblock,1)
    return
end
for k = 1:min(numel(sn.syncreply), size(sn.replyblock,2))
    if sn.syncreply(k) + 1 == r && sn.replyblock(ind,k) ~= 0
        tf = true;
        return
    end
end
end

function reached = local_reached(self, sn)
% LOCAL_REACHED  The (node, class) pairs a job can actually ARRIVE at, as an
% NNODES x NCLASSES logical.
%
% A forward walk of sn.rtnodes from the feed points. See the header of
% CHECKSERVICEREACHABLE for why the evidence is the UNMASKED routing kernel
% rather than sn.nodevisits, and for the three rules -- absorbing Sink, Source
% expanded only as a seed, unservable pair reached but not expanded -- that keep
% the walk from over-approximating.
%
% Returns [] when rtnodes is unreadable or is not the expected (N*R) square,
% which leaves the caller checking nothing, exactly as before.
reached = [];
if ~isfield(sn,'rtnodes') || isempty(sn.rtnodes)
    return
end
N = length(self.nodes);
R = sn.nclasses;
rt = sn.rtnodes;
if N < 1 || R < 1 || size(rt,1) < N*R || size(rt,2) < N*R
    return
end
reached = false(N,R);
seed = false(N,R);
for ind = 1:N
    nd = self.nodes{ind};
    if ~isa(nd,'Source')
        continue
    end
    for r = 1:R
        if r <= numel(nd.input.sourceClasses) && ~isempty(nd.input.sourceClasses{r}) && ...
                ~nd.input.sourceClasses{r}{end}.isDisabled
            seed(ind,r) = true;
        end
    end
end
for r = 1:R
    if r > numel(sn.njobs) || ~isfinite(sn.njobs(r)) || sn.njobs(r) <= 0
        continue
    end
    if ~isfield(sn,'refstat') || r > numel(sn.refstat) || sn.refstat(r) < 1
        continue
    end
    ind = sn.stationToNode(sn.refstat(r));
    if ind >= 1 && ind <= N
        seed(ind,r) = true;
    end
end
[si, sr] = find(seed);
stack = [si(:), sr(:)];
for t = 1:size(stack,1)
    reached(stack(t,1), stack(t,2)) = true;
end
while ~isempty(stack)
    ind = stack(end,1);
    r = stack(end,2);
    stack(end,:) = [];
    nd = self.nodes{ind};
    if isa(nd,'Sink')
        continue
    end
    if isa(nd,'Source') && ~seed(ind,r)
        continue
    end
    if ~local_serves(self, sn, ind, r)
        continue
    end
    row = rt((ind-1)*R + r, 1:(N*R));
    for t = find(row > GlobalConstants.Zero)
        j = floor((t-1)/R) + 1;
        s = mod(t-1, R) + 1;
        if ~reached(j,s)
            reached(j,s) = true;
            stack(end+1,:) = [j, s]; %#ok<AGROW>
        end
    end
end
end

function tf = local_serves(self, sn, ind, r)
% LOCAL_SERVES  True when a job of class R can LEAVE node IND again.
%
% True for anything that is not a service station, and for a station that serves
% R, declares heterogeneous server types (its per-class slot is empty by
% construction) or holds a server across a synchronous call whose reply class is
% R. False only for the flow sink itself, which is what stops the walk -- the
% same three exemptions the guard applies, kept in one place so the walk and the
% verdict cannot drift apart.
tf = true;
nd = self.nodes{ind};
if ~(isa(nd,'Queue') || isa(nd,'Delay'))
    return
end
if isprop(nd,'heteroServiceDistributions') && ~isempty(nd.heteroServiceDistributions) ...
        && numEntries(nd.heteroServiceDistributions) > 0
    return
end
svc = [];
if r <= length(nd.server.serviceProcess) && ~isempty(nd.server.serviceProcess{r})
    svc = nd.server.serviceProcess{r}{end};
end
if ~isempty(svc) && ~isa(svc,'Disabled')
    return
end
tf = holdsReplyFor(sn, ind, r);
end
