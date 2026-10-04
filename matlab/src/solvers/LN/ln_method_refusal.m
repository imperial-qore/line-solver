function [ok, reason] = ln_method_refusal(lqn, method, layering)
% [OK, REASON] = LN_METHOD_REFUSAL(LQN, METHOD, LAYERING)
% Whether a layered method can encode this model, asked as a predicate.
%
% ONE PREDICATE, TWO CALLERS. SolverLN.supportsModelMethod asks it, so a report
% never offers a pair the builders refuse on contact, and BUILDLAYERS /
% BUILDLAYERSPH ask it again on the run path and raise its REASON, so the gate
% and the run cannot say different things about one model. Until this predicate
% held every rule the builders carried their own copies: 'flat.cs' refused a
% replicated processor, a cache task or a setup task, every routing encoding
% refused a non-series-parallel fork graph and a routed call group under srvn,
% and the PH encodings refused a second phase, none of it visible to model.help.
%
% NONE OF THESE IS EXPRESSIBLE AS A FEATURE SET. A feature set says "I accept
% construct X", so it can refuse a model for HAVING one; but 'replicated
% processors need a submodel each' is a property of the SQUASHING rather than of
% the construct, and 'routed call groups state a dispatch order the composed law
% folds away' is a property of the ENCODING. Both are refusals about what the
% method does to the model, which is what a predicate is for.
%
% METHOD is any method name LQN_LN_METHOD accepts. The alias 'srvn' (and 'default')
% takes 'srvn.ph' where it can serve the model and 'srvn.cs' otherwise, and
% 'srvn.ph' never rescues a model 'srvn.cs' refuses, so the alias is gated as
% 'srvn.cs'. LAYERING is options.config.layering ('srvn' when omitted); a method
% that names a layering sets it, as in BUILDLAYERS.
%
% ONE RULE STAYS ON THE RUN PATH: whether an activity graph composes into one
% phase-type law under the PH encodings is found out by composing it
% (phInitLaws), which 'srvn' probes for and falls back from, so it is not a
% model property this predicate states.
%
% REASON is empty exactly when OK is true, and is the builder's own wording.

ok = true;
reason = '';
if nargin < 2 || isempty(method)
    return
end
if nargin < 3 || isempty(layering)
    layering = 'srvn';
end
if isempty(lqn) || ~isstruct(lqn)
    return
end
requested = lqn_ln_method(char(method));
layering = lower(char(layering));
if any(strcmp(requested, {'flat.cs','flat.ph'}))
    layering = 'flat'; % the method names the layering, so it sets it
end
isFlat = any(strcmp(layering, {'flat','squashed'}));
if strcmp(requested, 'srvn')
    requested = 'srvn.cs';
end
isPH = any(strcmp(requested, {'srvn.ph','flat.ph'}));
method = requested;
if strcmp(requested, 'srvn.ph') && isFlat
    % the same refusal buildLayers raises: the encoding needs a submodel per server
    ok = false;
    reason = sprintf(['method=''srvn.ph'' requires the srvn layering, because it replaces each ' ...
        'server by a submodel of its own; got ''%s''. Use method=''srvn.cs'' for that layering.'], layering);
    return
end

nelem = 0;
if isfield(lqn,'nhosts') && isfield(lqn,'ntasks')
    nelem = lqn.nhosts + lqn.ntasks;
end

% -- the squashing refusals, under either encoding ---------------------------
% Each carries PER-LAYER state that one submodel cannot hold, so they are
% properties of the flattening and not of the encoding: the srvn layering places
% the same stations across several networks and takes none of these.
if isFlat && nelem > 0
    if isfield(lqn,'repl') && any(full(lqn.repl(1:nelem)) > 1)
        ok = false;
        if isPH
            reason = ['method=''flat.ph'' does not support replicated processors or tasks, ' ...
                'whose replicas need a submodel each. Use method=''srvn.ph''.'];
        else
            reason = 'Flat layering does not support replicated processors or tasks, use the default ''srvn'' layering.';
        end
        return
    end
    if ~isPH && isfield(lqn,'iscache') && any(full(lqn.iscache(1:nelem)))
        ok = false;
        reason = 'Flat layering does not support cache tasks, use the default ''srvn'' layering.';
        return
    end
    if isfield(lqn,'hassetup') && any(full(lqn.hassetup(1:nelem)))
        ok = false;
        if isPH
            reason = ['method=''flat.ph'' does not support setup tasks, whose powered-down ' ...
                'threads are per-layer state. Use method=''srvn.ph''.'];
        else
            reason = 'Flat layering does not support setup tasks, use the default ''srvn'' layering.';
        end
        return
    end
end

if ~isPH
    % -- the routing encodings: srvn.cs, flat.cs and moment3 -----------------
    % buildLayersRecursive pairs a join with the most recent fork through a LIFO
    % stack of fork classes, so it can only represent series-parallel graphs.
    reason = seriesParallelForkRefusal(lqn);
    if ~isempty(reason)
        ok = false;
        return
    end
    % A routed call group states the order in which one caller visits several
    % callees. The srvn layering puts every callee in a submodel of its own and
    % replaces it, in the caller's submodel, by a surrogate delay, so the
    % callees are never co-resident and the order has nowhere to be expressed;
    % the squashed layering keeps them as stations of one model.
    if ~isFlat && isfield(lqn, 'callgroups') && ~isempty(lqn.callgroups)
        ok = false;
        reason = ['Call groups routed by a routing strategy require the ', ...
            'squashed layering; set options.config.layering=''flat''. Under srvn the ', ...
            'targets never share a submodel, so the dispatch order cannot be represented.'];
        return
    end
    return
end

% -- the composed-entry-law refusals, both PH encodings --------------------
% The composed entry law has no reply point, so a second phase has nowhere to
% go. lqn.actphase is written by LayeredNetwork.getStruct, so the rule is
% decidable from the model alone (SolverLN.hasPhase2 is this same test).
if isfield(lqn,'actphase') && any(full(lqn.actphase(:)) > 1)
    ok = false;
    reason = sprintf(['method=''%s'' does not support second-phase activities: the ' ...
        'composed entry law has no reply point. Use method=''default''.'], method);
    return
end
if isfield(lqn,'ncalls') && lqn.ncalls > 0 && isfield(lqn,'calltype') ...
        && any(lqn.calltype == CallType.FWD)
    ok = false;
    reason = sprintf(['method=''%s'' does not support forwarding calls, whose target is not ' ...
        'part of the caller''s activity graph. Use method=''default''.'], method);
    return
end
if isfield(lqn,'iscache') && any(full(lqn.iscache))
    ok = false;
    reason = sprintf('method=''%s'' does not support cache tasks. Use method=''default''.', method);
    return
end
% A SetupTask IS supported on a finite-multiplicity task: the setup is prefixed
% to the composed entry law as the mixture p*(setup THEN entry) + (1-p)*entry.
% An INF task is the exception, as in LDES: it holds no thread to power down.
if isfield(lqn,'hassetup') && any(full(lqn.hassetup)) && isfield(lqn,'sched') ...
        && isfield(lqn,'mult')
    for tidx_su = find(full(lqn.hassetup(:)))'
        if lqn.sched(tidx_su) == SchedStrategy.INF || ~isfinite(full(lqn.mult(tidx_su)))
            ok = false;
            reason = sprintf(['method=''%s'': task ''%s'' declares a setup time on an ' ...
                'infinite-server task, which holds no thread to power down; give it a ' ...
                'finite multiplicity.'], method, lqn.names{tidx_su});
            return
        end
    end
end
if isfield(lqn,'callgroups') && ~isempty(lqn.callgroups)
    % The group states the ORDER in which one caller visits several callees, and
    % the composed law folds every call into one visit, so the order has nowhere
    % to be expressed. Squashing does not recover it.
    ok = false;
    reason = sprintf(['method=''%s'' does not support routed call groups, whose dispatch ' ...
        'order is a routing property. Use method=''flat.cs''.'], method);
    return
end
if isfield(lqn,'lincon') && ~isempty(lqn.lincon) && any(~cellfun(@isempty, lqn.lincon(:,1)))
    ok = false;
    reason = sprintf(['method=''%s'' does not support admission constraints on a layer ' ...
        'station. Use method=''default''.'], method);
    return
end
% A queue-dependent service rate is a property of the layer STATION, and the
% composed law replaces that station by an entry law, so the scaling has nowhere
% to attach. Only buildLayersRecursive emits it.
depFields = {'lldscaling','cdscaling','jdscaling','pools'};
for kdep = 1:numel(depFields)
    fndep = depFields{kdep};
    if ~isfield(lqn, fndep) || isempty(lqn.(fndep))
        continue
    end
    sidxdep = find(~cellfun(@isempty, lqn.(fndep)(:)), 1);
    if ~isempty(sidxdep)
        ok = false;
        reason = sprintf(['method=''%s'' does not support queue-dependent service rates on ' ...
            'a layer station (''%s'' declares %s). Use method=''srvn.cs''.'], ...
            method, lqn.names{sidxdep}, fndep);
        return
    end
end
end

function reason = seriesParallelForkRefusal(lqn)
% Activity graphs whose AND forks and joins are not properly nested, which the
% routing encoding cannot represent; '' when every join closes one fork.
reason = '';
if ~isfield(lqn,'graph') || ~isfield(lqn,'actpretype') || ~isfield(lqn,'actposttype') ...
        || ~isfield(lqn,'nentries') || ~isfield(lqn,'nidx')
    return
end
ashift = lqn.nhosts + lqn.ntasks + lqn.nentries;
% kept as rows: a column mask would broadcast against the row index lists below
isPreAndAct = reshape(full(lqn.actpretype),1,[])==ActivityPrecedenceType.PRE_AND;
isPostAndAct = reshape(full(lqn.actposttype),1,[])==ActivityPrecedenceType.POST_AND;
acts = (ashift+1):lqn.nidx;
isFork = false(1,lqn.nidx);
for f = acts
    succ = find(lqn.graph(f,:));
    if any(isPostAndAct(succ(succ>ashift)))
        isFork(f) = true;
    end
end
for j = acts
    inputs = find(lqn.graph(:,j))';
    inputs = inputs(inputs>ashift & isPreAndAct(inputs));
    if numel(inputs) < 2
        continue
    end
    forks = arrayfun(@(a) enclosingFork(a, lqn, isFork, ashift), inputs);
    if numel(unique(forks)) > 1 || any(forks < 0)
        reason = sprintf(['Activity ''%s'' joins branches of different AND forks; ' ...
            'SolverLN supports only properly nested (series-parallel) fork-join graphs.'], lqn.hashnames{j});
        return
    end
end
end

function f = enclosingFork(a, lqn, isFork, ashift)
% Fork whose branch subtree contains activity A, -1 when A is outside every fork
seen = false(1,lqn.nidx);
queue = a;
while ~isempty(queue)
    cur = queue(1); queue(1) = [];
    if seen(cur), continue; end
    seen(cur) = true;
    preds = find(lqn.graph(:,cur))';
    preds = preds(preds>ashift);
    hit = preds(isFork(preds));
    if ~isempty(hit)
        f = hit(1);
        return
    end
    queue = [queue, preds]; %#ok<AGROW>
end
f = -1;
end
