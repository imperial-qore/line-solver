function [nonfjmodel, fjclassmap, fjforkmap, fanout] = mmt(model, forkLambda)
% build a model with fork-joins replaced by routers and delays and
% parallelism simulated by artificial classes
% forkLambda(s) is the arrival rate of artificial class s
% s = fjclassmap(r) for auxiliary class r gives the index s of the original class
% f = fjforkmap(r) for auxiliary class r gives the index of the associated fork node f
% fo = fanout(r) is the number of output jobs across all links for the (fork f,class s) pair modelled by auxiliary class r

%% this has been migrated to Java inside FJ.java as mmt

sn = model.getStruct;
if nargin < 2
    forkLambda = GlobalConstants.FineTol * ones(1,sn.nclasses);
end
fjclassmap = [];
fjforkmap = [];
fanout = [];
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
        for r=inchain
            oclass{end+1} = OpenClass(nonfjmodel,[nonfjmodel.classes{r}.name,'.',nonfjmodel.nodes{f}.name]); %#ok<AGROW>
            fjclassmap(oclass{end}.index) = nonfjmodel.classes{r}.index;
            fjforkmap(oclass{end}.index) = f;
            s = fjclassmap(oclass{end}.index); % auxiliary class index
            if model.nodes{f}.output.tasksPerLink > 1
                line_warning(mfilename, 'There are no synchronisation delays implemented in MMT for multiple tasks per link. Results may be inaccurate.');
            end
            fanout(oclass{end}.index) = origfanout(f,r)*model.nodes{f}.output.tasksPerLink;
            if sn.nodevisits{fc}(f,r) == 0
                source.setArrival(oclass{end},Disabled.getInstance);
            else
                source.setArrival(oclass{end},Exp(forkLambda(r)));
            end
            % joins are now Delays, let us set their service time
            for i=1:sn.nnodes
                if sn.isstation(i)
                    switch sn.nodetype(i)
                        case NodeType.Join
                            nonfjmodel.nodes{i}.setService(oclass{end},Immediate());
                        case {NodeType.Source, NodeType.Fork}
                            %no-op
                        otherwise
                            nonfjmodel.nodes{i}.setService(oclass{end},model.nodes{i}.getService(model.classes{r}).copy());
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
                P{oclass{find(r==inchain,1)},oclass{find(s==inchain,1)}}(nonfjmodel.nodes{joinIdx},:) = 0.0;
            end
            P{oclass{find(r==inchain,1)},oclass{find(r==inchain,1)}}(source, nonfjmodel.nodes{f}) = 1.0;
            P{oclass{find(r==inchain,1)},oclass{find(r==inchain,1)}}(nonfjmodel.nodes{joinIdx},sink) = 1.0;
        end
    end
end
nonfjmodel.relink(P);
for f=forkIndexes
    for r=1:length(nonfjmodel.nodes{f}.output.outputStrategy)
        if strcmp(nonfjmodel.nodes{f}.output.outputStrategy{r}{2},RoutingStrategy.RAND)
            nonfjmodel.nodes{f}.output.outputStrategy{r}{1} = nonfjmodel.classes{r}.name;
            nonfjmodel.nodes{f}.output.outputStrategy{r}{2} = RoutingStrategy.DISABLED;
        end
    end
end
end
