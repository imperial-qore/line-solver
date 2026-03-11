function lqn=QN2LQN(model)

lqn = LayeredNetwork(model.getName());
sn = model.getStruct;

PH = Host(lqn, model.getName(), Inf, SchedStrategy.INF); % pseudo host
for c=1:sn.nchains
    inchain = sn.inchain{c};    
    RT{c} = Task(lqn,['RefTask_',num2str(c)], sum(sn.njobs(inchain)), SchedStrategy.REF).on(PH); % reference task for chain c
    RE{c} = Entry(lqn,['Chain_',num2str(c)]).on(RT{c}); % entry on reference task for chain c
end

for i=1:sn.nnodes
    switch sn.nodetype(i)
        case {NodeType.Queue, NodeType.Delay}
            P{i} = Host(lqn, sn.nodenames{i}, sn.nservers(sn.nodeToStation(i)), SchedStrategy.fromId(sn.sched(sn.nodeToStation(i))));
            T{i} = Task(lqn,['T_',sn.nodenames{i}], Inf, SchedStrategy.INF).on(P{i});
            for r=1:sn.nclasses
                c = find(sn.chains(:,r)); % chain of class r
                if sn.visits{c}(i,r)>0
                    E{i,r} = Entry(lqn, ['E',num2str(i),'_',num2str(r)]).on(T{i});
                    A{i,r} = Activity(lqn,  ['Q',num2str(i),'_',num2str(r)], model.nodes{i}.getServiceProcess(model.classes{r})).on(T{i}).boundTo(E{i,r}).repliesTo(E{i,r});
                end
            end
        case {NodeType.ClassSwitch, NodeType.Router, NodeType.Logger}
            % no-op: passthrough routing nodes
        case {NodeType.Source, NodeType.Sink}
            % no-op: chain boundary nodes
        case {NodeType.Fork, NodeType.Join}
            % no-op: synchronization nodes, precedences added later
    end
end

PA = cell(sn.nchains, sn.nnodes, sn.nclasses);
boundToRE = cell(1,sn.nchains);
for i=1:sn.nnodes
    switch sn.nodetype(i)
        case {NodeType.ClassSwitch, NodeType.Router, NodeType.Logger}
            for r=1:sn.nclasses
                c = find(sn.chains(:,r)); % chain of class r
                if any(sn.rtnodes(:, ((i-1)*sn.nclasses + r))>0)
                    PA{c,i,r} = Activity(lqn,  ['CS','_',num2str(c),'_',num2str(i),'_',num2str(r)], Immediate()).on(RT{c}); % pseudo-activity in ref task
                end
            end
        case {NodeType.Fork, NodeType.Join}
            for r=1:sn.nclasses
                c = find(sn.chains(:,r)); % chain of class r
                if any(sn.rtnodes(:, ((i-1)*sn.nclasses + r))>0)
                    PA{c,i,r} = Activity(lqn,  ['FJ','_',num2str(c),'_',num2str(i),'_',num2str(r)], Immediate()).on(RT{c}); % pseudo-activity for fork/join
                end
            end
        case {NodeType.Queue, NodeType.Delay}
            for r=1:sn.nclasses
                c = find(sn.chains(:,r)); % chain of class r
                if sn.visits{c}(i,r)>0
                    inchain = sn.inchain{c};
                    if i == sn.refstat(inchain(1)) && r == inchain(1)
                        PA{c,i,r} = Activity(lqn,  ['A',num2str(i),'_',num2str(r)], Immediate()).on(RT{c}).boundTo(RE{c}).synchCall(E{i,r}); % pseudo-activity in ref task
                        boundToRE{c} = [i,r];
                    else
                        PA{c,i,r} = Activity(lqn,  ['A',num2str(i),'_',num2str(r)], Immediate()).on(RT{c}).synchCall(E{i,r}); % pseudo-activity in ref task
                    end
                    %RT{c}.addPrecedence(ActivityPrecedence.Serial(PN{c,i},PA{i,r}));
                end
            end
        case {NodeType.Source, NodeType.Sink}
            % no-op: chain boundary nodes
    end
end

usedInORFork = zeros(sn.nnodes,sn.nclasses);
routingNodeTypes = [NodeType.Queue, NodeType.Delay, NodeType.ClassSwitch, NodeType.Router, NodeType.Logger, NodeType.Fork, NodeType.Join];
for c=1:sn.nchains
    inchain = sn.inchain{c};
    %refstat = sn.refstat(inchain(1));
    for i=1:sn.nnodes
        if ~ismember(sn.nodetype(i), routingNodeTypes)
            continue;
        end
        for r=inchain
            orfork_prec = {};
            orfork_prob = [];
            for j=1:sn.nnodes
                if ~ismember(sn.nodetype(j), routingNodeTypes)
                    continue;
                end
                % Skip Join destinations - handled by AND-Join precedences below
                if sn.nodetype(j) == NodeType.Join
                    continue;
                end
                for s=inchain
                    pr = sn.rtnodes((i-1)*sn.nclasses + r, (j-1)*sn.nclasses + s);
                    if pr>0 && any(sn.rtnodes(:, ((i-1)*sn.nclasses + r))>0)
                        if ~isempty(boundToRE{c})
                            if boundToRE{c}(1)==j && boundToRE{c}(2)==s
                                if ~isempty(PA{c,i,r})
                                    orfork_prec{end+1} = Activity(lqn,  ['End_',num2str(c),'_',num2str(i),'_',num2str(r)], Immediate()).on(RT{c});
                                    orfork_prob(end+1)= pr;
                                end
                            else
                                orfork_prec{end+1} = PA{c,j,s};
                                orfork_prob(end+1) = pr;
                            end
                        end
                    end
                end
            end
            if ~isempty(orfork_prec)
                if ~isempty(PA{c,i,r})
                    if sn.nodetype(i) == NodeType.Fork
                        % Fork node: AND-Fork (all branches taken simultaneously)
                        RT{c}.addPrecedence(ActivityPrecedence.AndFork(PA{c,i,r}, orfork_prec));
                    else
                        RT{c}.addPrecedence(ActivityPrecedence.OrFork(PA{c,i,r}, orfork_prec, orfork_prob));
                    end
                    usedInORFork(i,r) = usedInORFork(i,r) + 1;
                end
            end
        end
    end
end

% AND-Join precedences for Join nodes
for c=1:sn.nchains
    inchain = sn.inchain{c};
    for j=1:sn.nnodes
        if sn.nodetype(j) ~= NodeType.Join
            continue;
        end
        for s=inchain
            if isempty(PA{c,j,s})
                continue;
            end
            join_pre = {};
            for i=1:sn.nnodes
                if ~ismember(sn.nodetype(i), routingNodeTypes)
                    continue;
                end
                for r=inchain
                    pr = sn.rtnodes((i-1)*sn.nclasses + r, (j-1)*sn.nclasses + s);
                    if pr > 0 && ~isempty(PA{c,i,r})
                        join_pre{end+1} = PA{c,i,r};
                    end
                end
            end
            if length(join_pre) > 1
                RT{c}.addPrecedence(ActivityPrecedence.AndJoin(join_pre, PA{c,j,s}));
            elseif length(join_pre) == 1
                RT{c}.addPrecedence(ActivityPrecedence({join_pre{1}.getName},{PA{c,j,s}.getName}));
            end
        end
    end
end

end