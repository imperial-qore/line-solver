function buildLayers(self)
lqn = self.lqn;
LineConsole.step('building the layer submodels of %d tasks', lqn.ntasks);
self.ensemble = cell(lqn.ntasks,1);

self.servt_classes_updmap = cell(lqn.nhosts+lqn.ntasks,1);
self.call_classes_updmap = cell(lqn.nhosts+lqn.ntasks,1);
self.arvproc_classes_updmap = cell(lqn.nhosts+lqn.ntasks,1);
self.thinkt_classes_updmap = cell(lqn.nhosts+lqn.ntasks,1);
self.actthinkt_classes_updmap = cell(lqn.nhosts+lqn.ntasks,1);
self.route_prob_updmap = cell(lqn.nhosts+lqn.ntasks,1);
self.singleReplicaTasks = [];

% see _kb/06-solver-catalog.md (LN section) for the layering taxonomy
layering = 'srvn';
if isfield(self.options.config,'layering') && ~isempty(self.options.config.layering)
    layering = lower(self.options.config.layering);
end

% Method resolution. A method name carries both the LAYERING and the ENCODING:
% 'srvn.ph' replaces the routing encoding of the activity graph by a composed
% phase-type server law, 'srvn' is the alias that takes it where it can serve
% the model and 'srvn.cs' otherwise, and 'flat.cs' squashes every server into
% one submodel. The choice is made ONCE, here, and every later dispatch reads
% self.lnmethod rather than options.method, so the reconstruction can never
% disagree with the layers that were built.
% See _kb/06-solver-catalog.md (LN section).
requested = lqn_ln_method(self.options.method);
if any(strcmp(requested, {'flat.cs','flat.ph'}))
    % the method names the layering, so it sets it
    layering = 'flat';
end
if strcmp(requested, 'flat.ph')
    % the squashed layering with the composed law: ONE submodel holding every
    % server, and a caller visiting each of them once per invocation. The
    % feature gate is the srvn.ph one plus the refusals a single submodel
    % carries -- see buildLayersPH.
    self.ph = [];
    buildLayersPH(self, 'build', true);
    self.lnmethod = 'flat.ph';
    return
end
if any(strcmp(requested, {'srvn.ph','srvn'}))
    hard = strcmp(requested, 'srvn.ph');
    if ~strcmp(layering, 'srvn')
        if hard
            line_error(mfilename, sprintf(['method=''srvn.ph'' requires the srvn layering, ' ...
                'because it replaces each server by a submodel of its own; got ''%s''. ' ...
                'Use method=''srvn.cs'' for that layering.'], layering));
        end
        line_debug('LN: method=srvn cannot use srvn.ph under the %s layering', layering);
    else
        self.ph = [];
        if hard
            buildLayersPH(self);
            self.lnmethod = 'srvn.ph';
            return
        elseif buildLayersPH(self, 'probe')
            buildLayersPH(self);
            self.lnmethod = 'srvn.ph';
            return
        end
    end
end
% The label reports what was BUILT, so a model squashed through
% options.config.layering reads back as 'flat.cs' even when no method named it.
if strcmp(requested, 'moment3')
    self.lnmethod = 'moment3';
elseif any(strcmp(layering, {'flat','squashed'}))
    self.lnmethod = 'flat.cs';
else
    self.lnmethod = 'srvn.cs';
end

% The encoding rules, from the predicate SolverLN.supportsModelMethod asks, so
% the gate and this run speak one sentence: the series-parallel fork nesting,
% the routed call groups under srvn, and the squashing refusals under flat.
[okm, whym] = ln_method_refusal(lqn, self.lnmethod, layering);
if ~okm
    line_error(mfilename, whym);
end
assertCallGroupSolver(self, lqn);

switch layering
    case 'srvn'
        flatServers = [];
    case {'flat','squashed'}
        flatServers = flatServerSet(self, lqn);
    otherwise
        line_error(mfilename, sprintf('Unknown layering strategy ''%s'', use ''srvn'' or ''flat''.', layering));
end

if isempty(flatServers)
    %% build one subnetwork for every processor
    for hidx = 1:lqn.nhosts
        if ~self.ignore(hidx)
            callers = lqn.tasksof{hidx};
            self.buildLayersRecursive(hidx, callers, true);
        else
            self.ensemble{hidx} = [];
        end
    end

    %% build one subnetwork for every task
    for t = 1:lqn.ntasks
        tidx = lqn.tshift + t;
        if ~self.ignore(tidx) & ~lqn.isref(tidx) & ~(isempty(find(self.lqn.iscaller(tidx,:), 1)) & isempty(find(self.lqn.iscaller(:,tidx), 1)))  %#ok<OR2,AND2> % ignore isolated tasks and ref tasks
            % obtain the activity graph of each task that calls some entry in t
            [calling_idx, called_entries] = find(lqn.iscaller(:, lqn.entriesof{tidx})); %#ok<ASGLU>
            callers = intersect(lqn.tshift+(1:lqn.ntasks), unique(calling_idx)');
            if ~isempty(callers) % true if the server is a software task
                self.buildLayersRecursive(tidx, callers, false);
            else
                self.ensemble{tidx} = [];
            end
        else
            self.ensemble{tidx} = [];
        end
    end
else
    %% flat layering: a single subnetwork holding every processor and task
    flatCallers = lqn.tshift + find(~self.ignore(lqn.tshift+(1:lqn.ntasks)))';
    self.buildLayersRecursive(flatServers, flatCallers, false, true);
end

self.thinkt_classes_updmap = cell2mat(self.thinkt_classes_updmap);
self.actthinkt_classes_updmap = cell2mat(self.actthinkt_classes_updmap);
self.call_classes_updmap = cell2mat(self.call_classes_updmap);
self.servt_classes_updmap = cell2mat(self.servt_classes_updmap);
self.arvproc_classes_updmap = cell2mat(self.arvproc_classes_updmap);
self.route_prob_updmap = cell2mat(self.route_prob_updmap);

% we now calculate the new index of the models after removing the empty
% models associated to 'ref' tasks
emptymodels = cellfun(@isempty,self.ensemble);
self.ensemble(emptymodels) = [];
self.idxhash = [1:length(emptymodels)]' - cumsum(emptymodels);
self.idxhash(emptymodels) = NaN;

%% Classify layers as host (processor) or task for MOL iteration
self.hostLayerIndices = [];
self.taskLayerIndices = [];

if ~isempty(flatServers)
    % every server resolves to the single flat layer, which is at once the
    % host layer and the task layer
    self.idxhash = NaN(lqn.nhosts+lqn.ntasks,1);
    self.idxhash(flatServers) = 1;
    self.hostLayerIndices = 1;
    self.taskLayerIndices = 1;
    self.model.ensemble = self.ensemble;
    return
end

% Host layers: indices 1:nhosts (before idxhash remapping)
for hidx = 1:lqn.nhosts
    if ~isnan(self.idxhash(hidx))
        self.hostLayerIndices(end+1) = self.idxhash(hidx);
    end
end

% Task layers: indices tshift+1:tshift+ntasks
for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    if ~isnan(self.idxhash(tidx))
        self.taskLayerIndices(end+1) = self.idxhash(tidx);
    end
end

% Layers carrying an admission constraint need the region wait recovered in
% updateMetricsDefault -- see _kb/06-solver-catalog.md (LN section)
self.layerHasRegion = false(length(self.ensemble),1);
self.layerChains = cell(length(self.ensemble),1);
for e = 1:length(self.ensemble)
    self.layerHasRegion(e) = ~isempty(self.ensemble{e}.regions);
    if self.layerHasRegion(e)
        % layer structure is iteration-invariant, so cache the chain matrix
        lsn = self.ensemble{e}.getStruct();
        self.layerChains{e} = lsn.chains;
    end
end

self.model.ensemble = self.ensemble;
end

function servers = flatServerSet(self, lqn)
% Processors and called tasks that become stations of the flat layer. The
% features a single submodel cannot carry (a replica, a cache task, a setup
% task) are refused by ln_method_refusal before this runs, not dropped here.
servers = [];
for hidx = 1:lqn.nhosts
    if ~self.ignore(hidx) && ~isempty(lqn.tasksof{hidx})
        servers(end+1) = hidx; %#ok<AGROW>
    end
end
for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    if self.ignore(tidx) || lqn.isref(tidx)
        continue
    end
    if isempty(find(lqn.iscaller(tidx,:), 1)) && isempty(find(lqn.iscaller(:,tidx), 1))
        continue % isolated task
    end
    % ROW indices, and ANY over the entries to get them. iscaller(:, entriesof)
    % is a MATRIX as soon as the task has more than one entry, and a one-output
    % find on a matrix returns positions in the flattened column stack: a caller
    % of the task's k-th entry came back as its row plus (k-1)*nidx, which lands
    % outside the task range and reads as "nobody calls this task". A
    % multi-entry server was therefore left out of the flat layer, got idxhash
    % NaN, and crashed updateRoutingProbabilities on the routing update that
    % buildLayersRecursive had registered for it. Single-entry tasks were
    % unaffected, which is why this survived: for them the submatrix is a
    % column and the linear index IS the row index.
    calling_idx = find(any(lqn.iscaller(:, lqn.entriesof{tidx}), 2));
    if ~isempty(intersect(lqn.tshift+(1:lqn.ntasks), calling_idx'))
        servers(end+1) = tidx; %#ok<AGROW>
    end
end
if isempty(servers)
    line_error(mfilename, 'Flat layering found no server: the model has no processor with tasks.');
end
end

function assertCallGroupSolver(self, lqn)
% Rejects routed call groups under a layer solver that cannot dispatch them.
%
% The layering half of the rule (a group needs the squashed layering, where its
% targets share a submodel) is ln_method_refusal's, asked before this. What is
% left is a property of the SOLVER FACTORY, which no model gate can see: only a
% layer solver with state-dependent routing honours the strategy; MVA, NC and
% FLD would silently return the probabilistic split instead.
if ~isfield(lqn, 'callgroups') || isempty(lqn.callgroups)
    return
end
if ~layerSolverSupportsRoutedGroups(self)
    line_error(mfilename, ['Routed call groups need a layer solver with ', ...
        'state-dependent routing (CTMC or SSA); MVA, NC and FLD would silently ', ...
        'return the probabilistic split under a round-robin or JSQ label.']);
end
end

function tf = layerSolverSupportsRoutedGroups(self)
% True when the layer solver implements state-dependent routing. Matched on the
% factory's text because the handle carries no model to gate a featSet against;
% the python twin does the same. Both should move to a getFeatureSet() query on
% the first built layer, where a real Network exists.
tf = false;
if isempty(self.solverFactory)
    return
end
fname = func2str(self.solverFactory);
for tag = {'CTMC','SSA','LDES','JMT'}
    if contains(upper(fname), tag{1})
        tf = true;
        return
    end
end
end

