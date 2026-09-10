function [nonfjmodel, fjclassmap, fjforkmap, fanout, prov] = mmt(model, forkLambda)
% build a model with fork-joins replaced by routers and delays and
% parallelism simulated by artificial classes
% forkLambda(s) is the arrival rate of artificial class s
% s = fjclassmap(r) for auxiliary class r gives the index s of the original class
% f = fjforkmap(r) for auxiliary class r gives the index of the associated fork node f
% fo = fanout(r) is the number of output jobs across all links for the (fork f,class s) pair modelled by auxiliary class r
%
% prov records the provenance of every service slot of nonfjmodel, so that a
% caller holding on to this transformation across an outer fixed point (see
% SolverMVA.runAnalyzer, driven by SolverLN) can re-feed it from the base model
% instead of rebuilding it -- the rebuild costs a full model.copy() per outer
% iteration, and the transformation depends only on the fork topology, which
% does not change. See ModelAdapter.refreshServicesFromBase. Fields:
%   prov.serviceSrc    (nslots x 3) rows [nodeIdx, nonfjClassIdx, baseClassIdx]
%                      -- base-derived: re-read from the base model on reuse.
%   prov.immediateSlots (nslots x 2) rows [nodeIdx, nonfjClassIdx]
%                      -- owned by this transformation (the joins it turns into
%                      zero-service delays). These must be RESTORED, not merely
%                      left alone: the fork loop overwrites each join with its
%                      current synchronisation delay on every pass, so a reused
%                      model still holds the previous outer iteration's value.
%   prov.forkLambdaInit, prov.baseModel -- inputs needed to reproduce the
%                      auxiliary-class source arrivals a cold call would set.
% Slots in neither map are owned by the forkLambda fixed point (the auxiliary
% arrivals) and are reset from forkLambdaInit. Reading a transformation-owned
% slot from the base model corrupts the transform silently; keeping a stale one
% silently warm-starts the fork loop.

%% this has been migrated to Java inside FJ.java as mmt

sn = model.getStruct;
if nargin < 2
    forkLambda = GlobalConstants.FineTol * ones(1,sn.nclasses);
end
fjclassmap = [];
fjforkmap = [];
fanout = [];
% Provenance of nonfjmodel's service slots (see the header). The original
% classes at the original stations carry over from the copy() below, so their
% provenance is the identity; the joins overwrite theirs with Immediate and are
% recorded separately.
prov = struct('serviceSrc', zeros(0,3), 'immediateSlots', zeros(0,2), ...
              'auxArrival', zeros(0,3), 'auxSource', [], ...
              'forkLambdaInit', forkLambda, 'baseModel', model);
for i = 1:sn.nnodes
    if sn.isstation(i)
        switch sn.nodetype(i)
            case {NodeType.Source, NodeType.Fork, NodeType.Join}
                % Source/Fork hold no per-class service; Join is overwritten
                % with Immediate below and so is transformation-owned.
            otherwise
                for r = 1:sn.nclasses
                    if ~isempty(model.nodes{i}.getService(model.classes{r}))
                        prov.serviceSrc(end+1,:) = [i, r, r]; %#ok<AGROW>
                    end
                end
        end
    end
end
% we create an equivalent model without fj stations
nonfjmodel = model.copy();
nonfjmodel.allowReplace = true;
P = nonfjmodel.getLinkedRoutingMatrix;
nonfjmodel.resetNetwork(true);
nonfjmodel.resetStruct();
if isempty(P)
    line_error(mfilename,'SolverMVA can process fork-join networks only if their routing topology has been generated using Network.link.');
end
Vnodes = cellsum(sn.nodevisits);
forkedClasses = {};
forkIndexes = find(sn.nodetype == NodeType.Fork)';
% replaces forks and joins with routers
fanout = [];
for f=forkIndexes
    for r=1:size(P,1)
        if length(model.nodes{f}.output.outputStrategy{r})>2
            origfanout(f,r) = length(model.nodes{f}.output.outputStrategy{r}{3});
            for s=1:size(P,2)
                P{r,s}(f,:) = P{r,s}(f,:) / origfanout(f,r);
            end
        else
            origfanout(f,r) = 0;
        end
    end
    % replace Join with a Router
    nonfjmodel.nodes{f} = Router(nonfjmodel, nonfjmodel.nodes{f}.name);
    % replace Fork with a StatelessClassSwitcher that doesn't
    % change the classes
    %nonfjmodel.nodes{f} = ClassSwitch(nonfjmodel, nonfjmodel.nodes{f}.name, eye(sn.nclasses));
    forkedClasses{f,1} = find(Vnodes(f,:)>0); %#ok<AGROW>
end
for j=find(sn.nodetype == NodeType.Join)'
    % replace Join with an Infinite Server
    nonfjmodel.nodes{j} = Delay(nonfjmodel, nonfjmodel.nodes{j}.name);
    nonfjmodel.stations{model.nodes{j}.stationIndex} = nonfjmodel.nodes{j};
    for c=1:length(nonfjmodel.classes)
        nonfjmodel.nodes{j}.setService(nonfjmodel.classes{c},Immediate());
        prov.immediateSlots(end+1,:) = [j, c]; %#ok<AGROW>
    end
end
nonfjmodel.stations={nonfjmodel.stations{1:model.getNumberOfStations}}'; % remove automatically added station and put it where the join was
% if they don't exist already, add source and sink
if nonfjmodel.hasOpenClasses
    source = nonfjmodel.getSource;
    sink = nonfjmodel.getSink;
else
    source = Source(nonfjmodel,'Source');
    sink = Sink(nonfjmodel,'Sink');
end
for r=1:size(P,1)
    for s=1:size(P,2)
        P{r,s}(length(nonfjmodel.nodes),length(nonfjmodel.nodes)) = 0;
        P{s,r}(length(nonfjmodel.nodes),length(nonfjmodel.nodes)) = 0;
    end
end
nonfjmodel.connections = zeros(length(nonfjmodel.nodes));
oclass = {};
all_aux_class_indices = []; % aux class indices for post-relink routing fix
for f=forkIndexes
    % find join associated to fork f
    joinIdx = find(sn.fj(f,:));
    if length(joinIdx)>1
        line_error(mfilename,'SolverMVA supports at present only a single join station per fork node.');
    end
    % find chains associated to classes forked by f
    forkedChains = find(sum(sn.chains(:,forkedClasses{f}),2));
    for fc=forkedChains'
        % create a new open class for each class in forkedChains
        oclass = {};
        inchain = find(sn.chains(fc,:)); inchain = inchain(:)';

        % see _kb/12-interfaces-and-docs.md (@ModelAdapter: JMT export helpers) for rationale
        refStationIdx = sn.refstat(inchain(1));
        refStatefulIdx = sn.stationToStateful(refStationIdx);
        refNodeIdx = sn.statefulToNode(refStatefulIdx);
        K = sn.nclasses;
        reachableFromRef = false(sn.nnodes, K);
        bfsQueue = zeros(0, 2); % [node, class] pairs
        for ci = inchain
            reachableFromRef(refNodeIdx, ci) = true;
            bfsQueue(end+1,:) = [refNodeIdx, ci]; %#ok<AGROW>
        end
        while ~isempty(bfsQueue)
            cn = bfsQueue(1,1); cc = bfsQueue(1,2);
            bfsQueue(1,:) = [];
            for destNode = 1:sn.nnodes
                for dci = inchain
                    srcIdx = (cn-1)*K + cc;
                    dstIdx = (destNode-1)*K + dci;
                    if ~reachableFromRef(destNode, dci) && sn.rtnodes(srcIdx, dstIdx) > 0
                        reachableFromRef(destNode, dci) = true;
                        bfsQueue(end+1,:) = [destNode, dci]; %#ok<AGROW>
                    end
                end
            end
        end

        for r=inchain
            oclass{end+1} = OpenClass(nonfjmodel,[nonfjmodel.classes{r}.name,'.',nonfjmodel.nodes{f}.name]); %#ok<AGROW>
            fjclassmap(oclass{end}.index) = nonfjmodel.classes{r}.index;
            fjforkmap(oclass{end}.index) = f;
            s = fjclassmap(oclass{end}.index); % auxiliary class index
            % fanout is the SIBLING count, links times tasksPerLink, which is what
            % the auxiliary open class's rate (fanout-1)*forkLambda must carry: one
            % sibling is the closed method name's own, the other fanout-1 are open
            % traffic. NetworkSolver.fjFixedPoint synchronises on the same count by
            % replicating each branch time tasksPerLink times.
            fanout(oclass{end}.index) = origfanout(f,r)*model.nodes{f}.output.tasksPerLink;
            all_aux_class_indices(end+1) = oclass{end}.index; %#ok<AGROW>
            disableAux = origfanout(f,r) == 0 || ~reachableFromRef(f, r);
            if disableAux
                source.setArrival(oclass{end},Disabled.getInstance);
            else
                source.setArrival(oclass{end},Exp(forkLambda(r)));
            end
            % Auxiliary arrivals are owned by the forkLambda fixed point, not by
            % the base model. Record what a cold call sets here (note the rate is
            % indexed by the ORIGINAL class r) so that reuse can restore it.
            prov.auxArrival(end+1,:) = [oclass{end}.index, r, disableAux]; %#ok<AGROW>
            prov.auxSource = source;
            % joins are now Delays, let us set their service time
            for i=1:sn.nnodes
                if sn.isstation(i)
                    switch sn.nodetype(i)
                        case NodeType.Join
                            nonfjmodel.nodes{i}.setService(oclass{end},Immediate());
                            prov.immediateSlots(end+1,:) = [i, oclass{end}.index]; %#ok<AGROW>
                        case {NodeType.Source, NodeType.Fork}
                            %no-op
                        otherwise
                            nonfjmodel.nodes{i}.setService(oclass{end},model.nodes{i}.getService(model.classes{r}).copy());
                            prov.serviceSrc(end+1,:) = [i, oclass{end}.index, r]; %#ok<AGROW>
                    end
                else
                end
            end
        end

        for r=inchain
            for s=inchain
                P{oclass{find(r==inchain,1)},oclass{find(s==inchain,1)}} = P{r,s};
            end
        end
        for r=inchain
            for s=inchain
                P{oclass{find(r==inchain,1)},oclass{find(s==inchain,1)}}(source,:) = 0.0;
                if ~isempty(joinIdx)
                    P{oclass{find(r==inchain,1)},oclass{find(s==inchain,1)}}(nonfjmodel.nodes{joinIdx},:) = 0.0;
                end
            end
            if origfanout(f,r) > 0
                P{oclass{find(r==inchain,1)},oclass{find(r==inchain,1)}}(source, nonfjmodel.nodes{f}) = 1.0;
                if ~isempty(joinIdx)
                    P{oclass{find(r==inchain,1)},oclass{find(r==inchain,1)}}(nonfjmodel.nodes{joinIdx},sink) = 1.0;
                end
            end
        end
        % see _kb/12-interfaces-and-docs.md (@ModelAdapter: JMT export helpers) for rationale
        pnnodes = size(P{inchain(1),inchain(1)}, 1);
        maxclass = max(inchain);
        % Source/Sink/Fork/Join rows are mmt infrastructure: their aux
        % routing is set explicitly above, so they are never cleared.
        infra_mask = false(1, pnnodes);
        infra_mask(f) = true;
        if ~isempty(joinIdx), infra_mask(joinIdx) = true; end
        for nd = 1:pnnodes
            if isa(nonfjmodel.nodes{nd}, 'Source') || isa(nonfjmodel.nodes{nd}, 'Sink')
                infra_mask(nd) = true;
            end
        end
        % Class-aware BFS: find all (node, class) pairs reachable from Fork
        visited = false(pnnodes, maxclass);
        bfs_q = zeros(0, 2); % [node, class] pairs
        % Seed: classes entering the fork
        for ri = 1:length(inchain)
            r = inchain(ri);
            for si = 1:length(inchain)
                s = inchain(si);
                if any(P{r,s}(f,:) > 0) && ~visited(f, r)
                    visited(f, r) = true;
                    bfs_q(end+1,:) = [f, r]; %#ok<AGROW>
                end
            end
        end
        while ~isempty(bfs_q)
            cn = bfs_q(1,1); cc = bfs_q(1,2);
            bfs_q(1,:) = [];
            for si = 1:length(inchain)
                s = inchain(si);
                for nd = 1:pnnodes
                    if P{cc,s}(cn, nd) > 0 && ~visited(nd, s)
                        visited(nd, s) = true;
                        % Continue BFS unless this is Join
                        if isempty(joinIdx) || nd ~= joinIdx
                            bfs_q(end+1,:) = [nd, s]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
        % see _kb/12-interfaces-and-docs.md (@ModelAdapter: JMT export helpers) for rationale
        for nd = 1:pnnodes
            if ~infra_mask(nd)
                for ri = 1:length(inchain)
                    if ~visited(nd, inchain(ri))
                        for si = 1:length(inchain)
                            P{oclass{ri}, oclass{si}}(nd, :) = 0.0;
                        end
                    end
                end
            end
        end
        % see _kb/12-interfaces-and-docs.md (@ModelAdapter: JMT export helpers) for rationale
        for ri = 1:length(inchain)
            if ~any(visited(:, inchain(ri)))
                P{oclass{ri}, oclass{ri}}(source, sink) = 1.0;
            end
        end
    end
end
% see _kb/12-interfaces-and-docs.md (@ModelAdapter: JMT export helpers) for rationale
for nd = 1:length(nonfjmodel.nodes)
    if ~isa(nonfjmodel.nodes{nd}, 'Station') && isprop(nonfjmodel.nodes{nd}, 'output') ...
            && ~isempty(nonfjmodel.nodes{nd}.output) && ismethod(nonfjmodel.nodes{nd}.output, 'initDispatcherJobClasses')
        nonfjmodel.nodes{nd}.output.initDispatcherJobClasses(nonfjmodel.classes);
    end
end
nonfjmodel.relink(P);
% see _kb/12-interfaces-and-docs.md (@ModelAdapter: JMT export helpers) for rationale
for nd = 1:length(nonfjmodel.nodes)
    os = nonfjmodel.nodes{nd}.output.outputStrategy;
    for ci = all_aux_class_indices
        if ci <= length(os)
            if isempty(os{ci}) || ~strcmp(os{ci}{2}, 'Probabilities')
                os{ci} = {nonfjmodel.classes{ci}.name, 'Probabilities', {}};
            end
        end
    end
    nonfjmodel.nodes{nd}.output.outputStrategy = os;
end
for f=forkIndexes
    for r=1:length(nonfjmodel.nodes{f}.output.outputStrategy)
        if strcmp(nonfjmodel.nodes{f}.output.outputStrategy{r}{2},RoutingStrategy.RAND)
            nonfjmodel.nodes{f}.output.outputStrategy{r}{1} = nonfjmodel.classes{r}.name;
            nonfjmodel.nodes{f}.output.outputStrategy{r}{2} = RoutingStrategy.DISABLED;
        end
    end
end
end
