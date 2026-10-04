%{ @file sn_to_qrf_blocking.m
 %  @brief Derives the QRF BAS blocking tables (f, MR, BB, MM, ZZ, MM1) from sn
 %
 %  @author LINE Development Team
%}

%{
 % @brief Builds the blocking configuration tables the QRF BAS bound needs
 %
 % @details
 % QRF_BAS describes a Blocking-After-Service network by a finite-capacity
 % queue F and an enumeration of the BLOCKING CONFIGURATIONS reachable behind
 % it. Everything in that enumeration is implied by the model, so this routine
 % derives it rather than asking the caller to hand-build it; the caller keeps
 % options.config.qrf_params as an explicit override.
 %
 % The tables, and the constraint that reads each one in QRF_BAS:
 %
 %   f          the ONE finite-capacity queue. The formulation carries a scalar
 %              f (ZERO4/ZERO7/ZERO8, THM30, THM3I, THM3L all index it), so a
 %              model with two binding buffers is refused here.
 %   F(i)       min(buffer size, N) for every queue; N where the buffer is
 %              unbounded, since no queue can hold more than the population.
 %   BB(m,i)    1 iff queue i is blocked in configuration m.
 %   ZZ(m)      nnz(BB(m,:)), the blocking depth of configuration m.
 %   MM(m,1)    head of the FIFO blocking order: the queue that takes the slot
 %              when f completes. QRF_BAS reads ONLY column 1 of MM.
 %   MM1(m,j)   index of the configuration reached from m when j becomes
 %              blocked. Read by THM3L alone, at depth ZM-1.
 %
 % THREE INVARIANTS, each a correctness condition rather than a convention:
 %
 %   1. Configuration 1 MUST be the empty one. ZERO4 iterates `m = 2:MR` and
 %      ZERO5/ZERO7/ZERO8 test `m >= 2` to mean "some queue is blocked".
 %   2. ZM = max(ZZ) MUST be the reachable maximum. QRF_BAS recomputes ZM from
 %      ZZ and closes the depth ladder there, so a truncated enumeration
 %      excises states the real chain visits and the polytope stops containing
 %      the true distribution -- the bound stops bounding. The size guard below
 %      therefore REFUSES; it never truncates.
 %   3. Blocking APPENDS at the tail: a queue that becomes blocked joins behind
 %      those already waiting, so MM1's successor is the configuration with j
 %      appended, and the head -- the queue MM(m,1) names -- never moves.
 %
 % THE ENUMERATION IS THE FULL ORDERED ONE, and it has to be. A (set, head)
 % collapse looks sound -- the LP reads configurations only through BB, ZZ,
 % MM(:,1) and MM1, and both objective and readout sum over m -- and it would
 % shrink MR from sum_z P(B,z) to 1 + sum_z C(B,z)*z. It was tried and it is
 % WRONG. Merging the depth-ZM configurations that share a set and a head makes
 % several THM3L rows, one per depth-(ZM-1) predecessor, reference the SAME
 % merged successor block. That is extra coupling the fine system does not
 % have, so the collapsed polytope is strictly SMALLER, not a projection of the
 % fine one, and it can cut off the true distribution. Measured on a 4-station
 % model with three feeders (B=3, ZM=3, MR 13 collapsed vs 16 full), the
 % collapse reported upper bounds of 0.681/0.979/0.768 where the full
 % enumeration gives 0.709/0.982/0.800: tighter, from a coarser state space,
 % which is the signature of a cut that is not valid. A coarser enumeration
 % that returns a tighter bound is reporting a number it cannot justify.
 %
 % So MR is factorial in the number of feeders B, and the size guard below is
 % what keeps that honest: it REFUSES an oversized instance rather than
 % trimming the enumeration, because trimming is the same unsound cut by
 % another name (invariant 2).
 %
 % WHO CAN BE BLOCKED is read from sn.isbasblocking / sn.isbasdestination, not
 % from sn.droprule: LINE accepts the BAS declaration on the upstream station
 % (cqn_bas_blocking.m) or on the full destination (the JMT/LDES form), and
 % reading droprule at the capped station sees only the second. That is the
 % BUG-83 distinction jmtIsBasDestination.m documents.
 %
 % @par Syntax:
 % @code
 % [blk, msg] = sn_to_qrf_blocking(sn, options)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>options<td>Solver options; reads config.qrf_maxvars
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>blk<td>Struct with fields f, F, MR, BB, MM, ZZ, MM1, ZM, blockers; empty when MSG is set
 % <tr><td>msg<td>Empty on success, otherwise the reason the tables cannot be derived
 % </table>
%}
function [blk, msg] = sn_to_qrf_blocking(sn, options)
blk = [];
msg = '';

if nargin < 2
    options = struct();
end

M = sn.nstations;
N = sum(sn.njobs);

[F, binding, msg] = sn_to_qrf_capacity(sn);
if ~isempty(msg)
    return
end

% ---- f: the one finite-capacity queue -----------------------------------
fcand = find(binding);
if numel(fcand) == 0
    % No binding buffer: the unblocked enumeration is the whole story. Callers
    % gate on sn_has_blocking before asking, so this is a defensive branch, not
    % a normal path.
    blk = qrf_empty_blocking(F, 1, M);
    return
end
if numel(fcand) > 1
    names = cell(1, numel(fcand));
    for c = 1:numel(fcand)
        names{c} = sn_station_name(sn, fcand(c));
    end
    msg = sprintf(['''qrf.bas'' models a single finite-capacity queue (its f is a scalar), but ' ...
        '%d stations have a binding buffer: %s. Use ''qrf.rsrd'', whose PBB constraint sums ' ...
        'over every full queue and therefore admits several, or cap only one station.'], ...
        numel(fcand), strjoin(names, ', '));
    return
end
f = fcand;

% ---- the queues that can be blocked behind f ----------------------------
blockers = qrf_blockers(sn, f, M);

% ---- reachable blocking depth -------------------------------------------
% Blocking needs f at capacity plus one held job per blocked queue, so the
% population caps the depth as tightly as the feeder count does.
ZM = min(numel(blockers), N - F(f));
if ZM < 0
    ZM = 0;
end

if ZM == 0
    % f cannot fill while any feeder still holds a job: no blocking state is
    % reachable and the single unblocked configuration is exact.
    blk = qrf_empty_blocking(F, f, M);
    blk.blockers = blockers;
    return
end

cfg = qrf_enumerate_permutations(blockers, ZM);

MR = numel(cfg);

% ---- size guard: refuse, never truncate (invariant 2) -------------------
Ktot = qrf_total_phases(sn, M);
nVars = MR * (N+1)^2 * Ktot^2 + Ktot;
maxVars = 5e5;
if isfield(options, 'config') && isstruct(options.config) && ...
        isfield(options.config, 'qrf_maxvars') && ~isempty(options.config.qrf_maxvars)
    maxVars = options.config.qrf_maxvars;
end
if nVars > maxVars
    msg = sprintf(['the QRF BAS linear program for this model would carry %.3g variables ' ...
        '(MR=%d blocking configurations, N=%d, %d service phases in total), above the ' ...
        'config.qrf_maxvars limit of %.3g. The enumeration cannot be truncated -- a depth ' ...
        'below the reachable maximum ZM=%d excises states the chain visits, and the result ' ...
        'would no longer bound. Reduce the population, the number of stations feeding %s, ' ...
        'or the phase counts; or raise config.qrf_maxvars deliberately.'], ...
        nVars, MR, N, Ktot, maxVars, ZM, sn_station_name(sn, f));
    return
end

% ---- fill BB, ZZ, MM, MM1 -----------------------------------------------
BB = zeros(MR, M);
ZZ = zeros(MR, 1);
MM = zeros(MR, max(2, numel(blockers)));
MM1 = zeros(MR, M);

for m = 1:MR
    seq = cfg{m};
    ZZ(m) = numel(seq);
    for z = 1:numel(seq)
        BB(m, seq(z)) = 1;
        MM(m, z) = seq(z); % only column 1 is read; the rest records the full order
    end
end

% MM1: successor under "j becomes blocked". Blocking appends at the tail, so
% the head is preserved (invariant 3) and the successor is unique on both
% enumerations.
keys = cell(MR, 1);
for m = 1:MR
    keys{m} = qrf_cfg_key(cfg{m});
end
for m = 1:MR
    if ZZ(m) >= ZM
        continue % no deeper configuration exists; THM3L reads MM1 only below ZM
    end
    for jj = 1:numel(blockers)
        j = blockers(jj);
        if BB(m, j) == 1
            continue
        end
        succ = [cfg{m}, j];
        sk = qrf_cfg_key(succ);
        for mp = 1:MR
            if strcmp(keys{mp}, sk)
                MM1(m, j) = mp;
                break
            end
        end
    end
end

blk = struct('f', f, 'F', F, 'MR', MR, 'BB', BB, 'MM', MM, 'ZZ', ZZ, 'MM1', MM1, ...
    'ZM', ZM, 'blockers', blockers);
end

function blk = qrf_empty_blocking(F, fidx, M)
% BLK = QRF_EMPTY_BLOCKING(F, FIDX, M)
%
% The one-configuration table for a model in which no blocking state is
% reachable. QRF_BAS reads it as a plain finite-buffer network.
blk = struct('f', fidx, 'F', F, 'MR', 1, 'BB', zeros(1, M), 'MM', zeros(1, 2), ...
    'ZZ', 0, 'MM1', zeros(1, M), 'ZM', 0, 'blockers', zeros(1, 0));
end

function blockers = qrf_blockers(sn, f, M)
% BLOCKERS = QRF_BLOCKERS(SN, F, M)
%
% Station indices that hold a completed job when F is full. A blocker must
% route into F, must not be F, and must not be an infinite server (which has a
% server per job and cannot be held). BAS itself is read from the BUG-83
% fields, with a structural fallback for an sn built without refreshLocalVars.
blockers = zeros(1, 0);
declared = false(1, M);
haveMarker = false;
if isfield(sn, 'isbasblocking') && ~isempty(sn.isbasblocking) && ...
        isfield(sn, 'stationToNode') && ~isempty(sn.stationToNode)
    haveMarker = true;
    for i = 1:M
        ind = sn.stationToNode(i);
        if ~isnan(ind) && ind >= 1 && numel(sn.isbasblocking) >= ind && sn.isbasblocking(ind) == 1
            declared(i) = true;
        end
    end
end
if ~haveMarker || ~any(declared)
    % Fallback: BAS declared on the upstream station or on the full
    % destination, read straight off sn.droprule. Same two forms
    % declaresBlockedMarker resolves, without needing refreshLocalVars.
    if isfield(sn, 'droprule') && ~isempty(sn.droprule)
        destBAS = any(sn.droprule(f, :) == DropStrategy.BAS);
        for i = 1:M
            if i == f
                continue
            end
            if destBAS || any(sn.droprule(i, :) == DropStrategy.BAS)
                declared(i) = true;
            end
        end
    end
end

for i = 1:M
    if i == f || ~declared(i)
        continue
    end
    if sn.sched(i) == SchedStrategy.INF
        continue
    end
    % Summed over class pairs so the test survives a multiclass sn, even though
    % the QRF gate upstream admits one class only.
    R = sn.nclasses;
    block = sn.rt((i-1)*R + (1:R), (f-1)*R + (1:R));
    if any(block(:) > 0)
        blockers(end+1) = i; %#ok<AGROW>
    end
end
blockers = sort(blockers);
end

function cfg = qrf_enumerate_permutations(blockers, ZM)
% CFG = QRF_ENUMERATE_PERMUTATIONS(BLOCKERS, ZM)
%
% Every ordered sequence of distinct blockers up to length ZM, as a cell of
% station vectors whose FIRST entry is the head. Deterministic order -- depth
% ascending, then subsets lexicographic by ascending station index, then the
% orders of each subset sorted -- so every codebase emits identical tables. The
% empty configuration is first (invariant 1).
cfg = {zeros(1,0)};
nb = numel(blockers);
for z = 1:ZM
    subsets = nchoosek(1:nb, z);
    for s = 1:size(subsets, 1)
        members = blockers(subsets(s, :));
        perms_z = perms(members);
        perms_z = sortrows(perms_z);
        for p = 1:size(perms_z, 1)
            cfg{end+1} = perms_z(p, :); %#ok<AGROW>
        end
    end
end
end

function key = qrf_cfg_key(seq)
% KEY = QRF_CFG_KEY(SEQ)
%
% Identity of a configuration: the whole blocking order, since that is what
% distinguishes configurations in the enumeration QRF_BAS is entitled to.
key = sprintf('%d,', seq);
end

function Ktot = qrf_total_phases(sn, M)
% KTOT = QRF_TOTAL_PHASES(SN, M)
%
% Total number of service phases across stations, which sizes the QRF variable
% space together with MR and the population.
Ktot = 0;
for i = 1:M
    if numel(sn.proc) >= i && ~isempty(sn.proc{i}) && ~isempty(sn.proc{i}{1})
        Ktot = Ktot + size(sn.proc{i}{1}{1}, 1);
    else
        Ktot = Ktot + 1;
    end
end
end

function name = sn_station_name(sn, ist)
% NAME = SN_STATION_NAME(SN, IST)
%
% Printable station name, falling back to the index when the struct carries no
% name table.
name = sprintf('station %d', ist);
if isfield(sn, 'nodenames') && ~isempty(sn.nodenames) && ...
        isfield(sn, 'stationToNode') && ~isempty(sn.stationToNode)
    ind = sn.stationToNode(ist);
    if ~isnan(ind) && ind >= 1 && numel(sn.nodenames) >= ind
        name = sn.nodenames{ind};
    end
end
end
