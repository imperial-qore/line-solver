function nodevisits = sn_fj_visits_spn(sn)
% SN_FJ_VISITS_SPN Compute fork-join node visit ratios via auxiliary SPN models
%
% NODEVISITS = SN_FJ_VISITS_SPN(SN) builds, for each class that passes
% through a fork-join pair, an auxiliary closed Stochastic Petri Net
% capturing the fork/join synchronization semantics. The SPN is solved
% with SolverCTMC and the throughput ratios give the per-node visit ratios.
%
% Population-preserving approach: the SPN uses population = B (number of
% fork branches). The pre-fork station consumes B tokens and produces 1
% per branch. The Join consumes 1 per branch and produces B. This ensures
% token conservation and correct CTMC analysis.
%
% The function handles:
%   - Multiple fork-join pairs (including serial fork-joins)
%   - Multiple classes
%
% Returns:
%   nodevisits - cell(1, nchains), each entry is (nnodes x nclasses)
%                with visit ratios normalized to reference station = 1
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

I = sn.nnodes;
K = sn.nclasses;
nchains = sn.nchains;
inchain = sn.inchain;
refstat = sn.refstat;

nodevisits = cell(1, nchains);
for c = 1:nchains
    nodevisits{c} = zeros(I, K);
end

if ~any(sn.fj(:))
    return
end

% For each chain, build and solve an SPN for each class in the chain
for c = 1:nchains
    for kidx = 1:length(inchain{c})
        r = inchain{c}(kidx);

        % Extract the single-class sub-routing from rtnodes
        P_r = zeros(I, I);
        for i = 1:I
            for j = 1:I
                P_r(i, j) = sn.rtnodes((i-1)*K + r, (j-1)*K + r);
            end
        end

        % Find nodes visited by this class (reachable from reference station)
        refnode = sn.stationToNode(refstat(r));
        visited = false(1, I);
        visited(refnode) = true;
        changed = true;
        while changed
            changed = false;
            for i = 1:I
                if visited(i)
                    for j = 1:I
                        if P_r(i, j) > 0 && ~visited(j)
                            visited(j) = true;
                            changed = true;
                        end
                    end
                end
            end
        end

        % Skip if this class doesn't pass through any fork
        has_fork = false;
        for i = 1:I
            if visited(i) && sn.nodetype(i) == NodeType.Fork
                has_fork = true;
                break;
            end
        end
        if ~has_fork
            for i = 1:I
                if visited(i)
                    nodevisits{c}(i, r) = 1;
                end
            end
            continue;
        end

        % Build auxiliary SPN model
        nodevisits{c}(:, r) = build_and_solve_spn(sn, P_r, visited, r, refnode);
    end

    % Normalize by reference station
    refnode_c = sn.stationToNode(refstat(inchain{c}(1)));
    for kidx = 1:length(inchain{c})
        r = inchain{c}(kidx);
        normVal = nodevisits{c}(refnode_c, r);
        if normVal > GlobalConstants.FineTol
            nodevisits{c}(:, r) = nodevisits{c}(:, r) / normVal;
        end
    end
end
end

function visits_r = build_and_solve_spn(sn, P_r, visited, r, refnode)
% BUILD_AND_SOLVE_SPN Construct and solve auxiliary SPN for one class
%
% Uses a population-preserving approach:
%   - Population = total leaf branches (B) across all fork-join pairs
%   - Pre-fork station: requires B tokens, produces 1 per leaf branch
%   - Branch stations: normal 1-in-1-out service
%   - Join: requires 1 from each branch done, produces B in post-join dest
%
% The throughput ratios are NOT uniform across the station Places, although an
% earlier version of this comment said so. The pre-fork transition consumes all
% B tokens at once and the Join returns them, so a station INSIDE a fork-join
% region fires once per B firings of the cycle: after the normalisation on the
% reference station, a station outside the region is 1, a station inside it is
% 1/B, and a Fork or a Join, which holds no Place, is 0. This solve is what
% defines the function; the JAR reproduces it through SolverCTMC, and the
% Python and C++ ports evaluate the same rule analytically.

I = sn.nnodes;
visits_r = zeros(I, 1);

visitedNodes = find(visited);
nVisited = length(visitedNodes);
if nVisited == 0
    return;
end

% Classify visited nodes
stationNodes = [];
forkNodes = [];
joinNodes = [];
for idx = 1:nVisited
    nd = visitedNodes(idx);
    if sn.nodetype(nd) == NodeType.Fork
        forkNodes(end+1) = nd; %#ok<AGROW>
    elseif sn.nodetype(nd) == NodeType.Join
        joinNodes(end+1) = nd; %#ok<AGROW>
    elseif sn.isstation(nd) && sn.nodetype(nd) ~= NodeType.Source && sn.nodetype(nd) ~= NodeType.Sink
        stationNodes(end+1) = nd; %#ok<AGROW>
    end
end

% Compute leaf count for each Join node (bottom-up, handles arbitrary nesting)
% join_leaves(jnd) = number of leaf tokens that flow through this Join
% (indexed by node index; NaN entries denote Joins not yet resolved)
join_leaves = nan(1, I);
for pass = 1:length(joinNodes) % iterate until stable
    for jidx = 1:length(joinNodes)
        jnd = joinNodes(jidx);
        srcNodes_j = find(P_r(:, jnd) > 0 & visited');
        lc = 0;
        for si = 1:length(srcNodes_j)
            srcNd = srcNodes_j(si);
            if sn.nodetype(srcNd) == NodeType.Join && ~isnan(join_leaves(srcNd))
                lc = lc + join_leaves(srcNd);
            else
                lc = lc + 1; % station or unresolved
            end
        end
        join_leaves(jnd) = lc;
    end
end

% Compute total leaf branches B (= population for conservation)
% For each fork, resolve to leaf stations; B = max leaf count across outermost forks
B = 0;
for idx = 1:length(forkNodes)
    fnd = forkNodes(idx);
    % Check if this Fork is an outermost fork (not a descendant of another fork)
    % An outermost fork has its upstream being a station, not another fork
    srcNodes = find(P_r(:, fnd) > 0 & visited');
    is_outermost = true;
    for si = 1:length(srcNodes)
        if sn.nodetype(srcNodes(si)) == NodeType.Join
            is_outermost = false; % serial fork - Join feeds this Fork
            break;
        end
    end
    if is_outermost
        [leafs, leafw] = resolve_fork_dests(sn, P_r, visited, fnd, r);
        % The EXPECTED sibling count, which is the leaf count exactly when every
        % link is certain and carries one task. Under a variable forking level it
        % is generally FRACTIONAL, and a fractional token population is not a net
        % anyone can enumerate: that is why the solve below is skipped exactly
        % there and the closed form answers instead.
        if isempty(leafs)
            line_error(mfilename, 'A Fork reaches no station on any branch, so the auxiliary net has nothing to synchronize.');
        end
        B = max(B, sum(leafw));
    end
end
if B == 0
    B = 1;
end
if B <= 0
    line_error(mfilename, 'A Fork emits no task in expectation, so its Join can never fire; at least one branch must be certain to emit at least one task.');
end

% THE NET BELOW IS THE PLAIN FORK'S NET. It puts B tokens in a closed SPN whose
% pre-fork transition emits ONE token per branch, so it is balanced only when
% every branch carries exactly one task: give it a fork whose links carry three
% each, and it deadlocks and reports zero throughput everywhere. Under a
% variable forking level B is the expected sibling count, which is generally not
% an integer either -- and a fractional token population is not a net anyone can
% enumerate. The rule that solve produces is available in closed form -- a
% station outside the fork-join region carries 1, one inside it carries 1/B, a
% Fork and a Join carry 0 -- and it is what answers whenever the fork varies.
% The plain fork keeps taking the solve, which is what holds the two in step.
isVariableFork = false;
for idx = 1:length(forkNodes)
    fnd = forkNodes(idx);
    if ~isfield(sn.nodeparam{fnd},'fanOutLink') || isempty(sn.nodeparam{fnd}.fanOutLink)
        continue
    end
    taken = sn.nodeparam{fnd}.fanOutProb > 0;
    if any(sn.nodeparam{fnd}.fanOutProb(taken) ~= 1) || ...
            any(sn.nodeparam{fnd}.fanOutLink(taken) ~= sn.nodeparam{fnd}.fanOut) || ...
            ~all(cellfun(@isempty, sn.nodeparam{fnd}.fanOutDist(:)))
        isVariableFork = true;
    end
end
if isVariableFork || abs(B - round(B)) > GlobalConstants.FineTol
    inRegion = false(1, I);
    for idx = 1:length(forkNodes)
        fnd = forkNodes(idx);
        srcNodes = find(P_r(:, fnd) > 0 & visited');
        if any(arrayfun(@(s) sn.nodetype(s) == NodeType.Join, srcNodes))
            continue % a serial stage, inside a region already fixed
        end
        frontier = fnd;
        while ~isempty(frontier)
            nd = frontier(1); frontier(1) = [];
            dests = find(P_r(nd, :) > 0 & visited);
            for di = 1:length(dests)
                j = dests(di);
                if sn.nodetype(j) == NodeType.Join || inRegion(j)
                    continue
                end
                inRegion(j) = true;
                frontier(end+1) = j; %#ok<AGROW>
            end
        end
    end
    for idx = 1:length(stationNodes)
        nd = stationNodes(idx);
        if inRegion(nd)
            visits_r(nd) = 1/B;
        else
            visits_r(nd) = 1;
        end
    end
    return
end

% Build SPN model
model = Network('fj_spn_aux');

% Create Places for station nodes
places = cell(1, I);
for idx = 1:length(stationNodes)
    nd = stationNodes(idx);
    places{nd} = Place(model, sprintf('P_%s', sn.nodenames{nd}));
end

% Create "done" Places for stations that route to a Join
pre_join = cell(1, I);
for idx = 1:length(stationNodes)
    nd = stationNodes(idx);
    dstNodes = find(P_r(nd, :) > 0 & visited);
    for di = 1:length(dstNodes)
        if sn.nodetype(dstNodes(di)) == NodeType.Join
            pre_join{nd} = Place(model, sprintf('P_%s_done', sn.nodenames{nd}));
            break;
        end
    end
end

% For inner Joins that route to outer Joins (nested fork-join)
for idx = 1:length(joinNodes)
    jnd = joinNodes(idx);
    dstNodes = find(P_r(jnd, :) > 0 & visited);
    for di = 1:length(dstNodes)
        if sn.nodetype(dstNodes(di)) == NodeType.Join
            pre_join{jnd} = Place(model, sprintf('P_%s_done', sn.nodenames{jnd}));
            break;
        end
    end
end

% Intermediate Places for Join -> Fork connections (serial fork-join)
inter_jf = cell(1, I);
for idx = 1:length(joinNodes)
    jnd = joinNodes(idx);
    dstNodes = find(P_r(jnd, :) > 0 & visited);
    for di = 1:length(dstNodes)
        if sn.nodetype(dstNodes(di)) == NodeType.Fork
            inter_jf{jnd} = Place(model, sprintf('P_%s_to_%s', sn.nodenames{jnd}, sn.nodenames{dstNodes(di)}));
            break;
        end
    end
end

% ClosedClass with B tokens at reference Place
jobclass = ClosedClass(model, 'Token', B, places{refnode});

transitions = {};
transInfo = {};

% Create timed service Transitions for each station
for idx = 1:length(stationNodes)
    nd = stationNodes(idx);
    ist = sn.nodeToStation(nd);

    % Service rate
    is_timed = false;
    rate = 0;
    if ist > 0 && ~isempty(sn.rates) && ist <= size(sn.rates, 1) && r <= size(sn.rates, 2) && sn.rates(ist, r) > 0 && ~isinf(sn.rates(ist, r))
        rate = sn.rates(ist, r);
        is_timed = true;
    end

    % Determine effective output Places by resolving through Fork/Join routing
    dstNodes = find(P_r(nd, :) > 0 & visited);
    out_places = {};
    enable_count = 1; % how many tokens to require/consume
    for di = 1:length(dstNodes)
        dstNd = dstNodes(di);
        if sn.nodetype(dstNd) == NodeType.Fork
            % Resolve Fork to leaf stations; require B tokens (population-preserving)
            forkDests = resolve_fork_dests(sn, P_r, visited, dstNd);
            enable_count = B; % require all B tokens
            for fi = 1:length(forkDests)
                if ~isempty(places{forkDests(fi)})
                    out_places{end+1} = places{forkDests(fi)}; %#ok<AGROW>
                end
            end
        elseif sn.nodetype(dstNd) == NodeType.Join
            if ~isempty(pre_join{nd})
                out_places{end+1} = pre_join{nd}; %#ok<AGROW>
            end
        elseif sn.isstation(dstNd)
            if ~isempty(places{dstNd})
                out_places{end+1} = places{dstNd}; %#ok<AGROW>
            end
        end
    end

    if isempty(out_places), continue; end

    tName = sprintf('T_svc_%s', sn.nodenames{nd});
    T = Transition(model, tName);
    mode = T.addMode('serve');

    if is_timed
        T.setDistribution(mode, Exp(rate));
    else
        T.setTimingStrategy(mode, TimingStrategy.IMMEDIATE);
        T.setFiringWeights(mode, 1.0);
    end

    T.setEnablingConditions(mode, jobclass, places{nd}, enable_count);
    T.setFiringOutcome(mode, jobclass, places{nd}, -enable_count);
    for oi = 1:length(out_places)
        T.setFiringOutcome(mode, jobclass, out_places{oi}, 1);
    end

    transitions{end+1} = T; %#ok<AGROW>
    transInfo{end+1} = struct('type', 'service', 'nd', nd, ...
        'inPlaces', {{places{nd}}}, 'outPlaces', {out_places}); %#ok<AGROW>
end

% Create Immediate Join Transitions
for idx = 1:length(joinNodes)
    jnd = joinNodes(idx);
    srcNodes = find(P_r(:, jnd) > 0 & visited');
    dstNodes = find(P_r(jnd, :) > 0 & visited);

    % Collect input "done" Places and their source nodes
    in_places = {};
    in_src_nodes = [];
    for si = 1:length(srcNodes)
        srcNd = srcNodes(si);
        if ~isempty(pre_join{srcNd})
            in_places{end+1} = pre_join{srcNd}; %#ok<AGROW>
            in_src_nodes(end+1) = srcNd; %#ok<AGROW>
        end
    end
    if isempty(in_places), continue; end

    % Compute enabling counts per input: 1 for stations, join_leaves for inner joins
    en_counts = ones(1, length(in_places));
    for ii = 1:length(in_src_nodes)
        srcNd = in_src_nodes(ii);
        if sn.nodetype(srcNd) == NodeType.Join && ~isnan(join_leaves(srcNd))
            en_counts(ii) = join_leaves(srcNd);
        end
    end

    % produce_count = join_leaves(jnd) = total leaf tokens through this Join
    if ~isnan(join_leaves(jnd))
        produce_count = join_leaves(jnd);
    else
        produce_count = sum(en_counts);
    end

    % Determine post-Join output Places
    out_places = {};
    for di = 1:length(dstNodes)
        dstNd = dstNodes(di);
        if sn.nodetype(dstNd) == NodeType.Fork
            % Serial fork-join
            if ~isempty(inter_jf{jnd})
                out_places{end+1} = inter_jf{jnd}; %#ok<AGROW>
            end
        elseif sn.nodetype(dstNd) == NodeType.Join
            % Nested: inner Join -> outer Join
            if ~isempty(pre_join{jnd})
                out_places{end+1} = pre_join{jnd}; %#ok<AGROW>
            end
        elseif sn.isstation(dstNd)
            if ~isempty(places{dstNd})
                out_places{end+1} = places{dstNd}; %#ok<AGROW>
            end
        end
    end
    if isempty(out_places), continue; end

    tName = sprintf('T_join_%s', sn.nodenames{jnd});
    T = Transition(model, tName);
    mode = T.addMode('sync');
    T.setTimingStrategy(mode, TimingStrategy.IMMEDIATE);
    T.setFiringWeights(mode, 1.0);

    for ii = 1:length(in_places)
        T.setEnablingConditions(mode, jobclass, in_places{ii}, en_counts(ii));
        T.setFiringOutcome(mode, jobclass, in_places{ii}, -en_counts(ii));
    end
    for oi = 1:length(out_places)
        T.setFiringOutcome(mode, jobclass, out_places{oi}, produce_count);
    end

    transitions{end+1} = T; %#ok<AGROW>
    transInfo{end+1} = struct('type', 'join', 'nd', jnd, ...
        'inPlaces', {in_places}, 'outPlaces', {out_places}); %#ok<AGROW>
end

% Create Immediate Fork Transitions for Join -> Fork connections (serial FJ)
for idx = 1:length(joinNodes)
    jnd = joinNodes(idx);
    if isempty(inter_jf{jnd}), continue; end

    dstNodes = find(P_r(jnd, :) > 0 & visited);
    for di = 1:length(dstNodes)
        dstNd = dstNodes(di);
        if sn.nodetype(dstNd) ~= NodeType.Fork, continue; end

        forkDests = resolve_fork_dests(sn, P_r, visited, dstNd);
        if isempty(forkDests), continue; end

        tName = sprintf('T_fork_%s_%s', sn.nodenames{jnd}, sn.nodenames{dstNd});
        T = Transition(model, tName);
        mode = T.addMode('fork');
        T.setTimingStrategy(mode, TimingStrategy.IMMEDIATE);
        T.setFiringWeights(mode, 1.0);

        T.setEnablingConditions(mode, jobclass, inter_jf{jnd}, B);
        T.setFiringOutcome(mode, jobclass, inter_jf{jnd}, -B);
        fork_out = {};
        for fi = 1:length(forkDests)
            if ~isempty(places{forkDests(fi)})
                T.setFiringOutcome(mode, jobclass, places{forkDests(fi)}, 1);
                fork_out{end+1} = places{forkDests(fi)}; %#ok<AGROW>
            end
        end

        transitions{end+1} = T; %#ok<AGROW>
        transInfo{end+1} = struct('type', 'fork', 'nd', jnd, ...
            'inPlaces', {{inter_jf{jnd}}}, 'outPlaces', {fork_out}); %#ok<AGROW>
    end
end

% Set up routing matrix
R = model.initRoutingMatrix();
for tidx = 1:length(transitions)
    T = transitions{tidx};
    info = transInfo{tidx};
    for ii = 1:length(info.inPlaces)
        R{1,1}(info.inPlaces{ii}, T) = 1;
    end
    for oi = 1:length(info.outPlaces)
        R{1,1}(T, info.outPlaces{oi}) = 1;
    end
end
model.link(R);

% Set initial state: B tokens at reference, 0 everywhere else
places{refnode}.setState(B);
for i = 1:I
    if ~isempty(places{i}) && i ~= refnode
        places{i}.setState(0);
    end
    if ~isempty(pre_join{i})
        pre_join{i}.setState(0);
    end
    if ~isempty(inter_jf{i})
        inter_jf{i}.setState(0);
    end
end

% Solve with CTMC
try
    solver = SolverCTMC(model);
    avg = solver.getAvgTable();

    % Extract throughput ratios as visit ratios
    stationNames = avg.Station;
    tputVals = avg.Tput;

    for idx = 1:nVisited
        nd = visitedNodes(idx);
        if ~isempty(places{nd})
            placeName = places{nd}.name;
            rowIdx = find(strcmp(string(stationNames), placeName));
            if ~isempty(rowIdx)
                visits_r(nd) = tputVals(rowIdx(1));
            end
        end
    end
catch ME
    line_warning(mfilename, sprintf('SPN CTMC solve failed for class %d: %s', r, ME.message));
end
end

function [stDests, stWeights] = resolve_fork_dests(sn, P_r, visited, forkNd, r, w)
% RESOLVE_FORK_DESTS Recursively resolve Fork destinations to station nodes
%
% The second output is the EXPECTED number of tasks each leaf receives per
% firing of the outermost fork: P(branch fires) times E[tasks on that link],
% multiplied down through any nesting. A plain fork gives every link exactly 1,
% so the weighted leaf count collapses to the leaf count it has always been.
if nargin < 5, r = 1; end
if nargin < 6, w = 1; end
stDests = [];
stWeights = [];
hasFan = isfield(sn.nodeparam{forkNd},'fanOutLink') && ~isempty(sn.nodeparam{forkNd}.fanOutLink);
branchDests = find(P_r(forkNd, :) > 0 & visited);
for bdi = 1:length(branchDests)
    bd = branchDests(bdi);
    wbd = w;
    if hasFan
        wbd = w * sn.nodeparam{forkNd}.fanOutProb(bd,r) * sn.nodeparam{forkNd}.fanOutLink(bd,r);
    end
    if sn.nodetype(bd) == NodeType.Fork
        [d2, w2] = resolve_fork_dests(sn, P_r, visited, bd, r, wbd);
        stDests = [stDests, d2]; %#ok<AGROW>
        stWeights = [stWeights, w2]; %#ok<AGROW>
    elseif sn.isstation(bd)
        stDests = [stDests, bd]; %#ok<AGROW>
        stWeights = [stWeights, wbd]; %#ok<AGROW>
    end
end
end
