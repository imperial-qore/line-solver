function ln_host_task_priorities(lqn, model, idxSet)
% LN_HOST_TASK_PRIORITIES(LQN, MODEL, IDXSET) set task priorities on a host layer
%
% Give every class of the layer MODEL the priority of the task that owns it,
% for each processor in IDXSET that schedules by priority (HOL, 'pri', ...).
% lqns serves a LARGER task priority first while a LINE class priority of 0 is
% the highest, so the class priority is max(prio on the host) - prio(task). A
% layer solver with a priority arm (MVA, CTMC, ...) then applies it. Without
% this every class tied at priority 0 and a priority processor ran as FCFS.
if ~isfield(lqn,'prio') || isempty(lqn.prio)
    return
end
prioScheds = [SchedStrategy.HOL, SchedStrategy.FCFSPRIO, SchedStrategy.FCFSPRPRIO, ...
    SchedStrategy.FCFSPIPRIO, SchedStrategy.LCFSPRIO, SchedStrategy.LCFSPRPRIO, ...
    SchedStrategy.LCFSPIPRIO, SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO];
for s = idxSet(:)'
    if s > lqn.nhosts || ~ismember(lqn.sched(s), prioScheds)
        continue
    end
    hostTasks = lqn.tshift + find(lqn.parent(lqn.tshift+(1:lqn.ntasks)) == s);
    pmax = max([0; reshape(lqn.prio(hostTasks), [], 1)]);
    for k = 1:length(model.classes)
        owner = layerClassOwnerTask(lqn, model.classes{k});
        if owner > 0 && lqn.parent(owner) == s
            model.classes{k}.setPriority(round(pmax - lqn.prio(owner)));
        end
    end
end
end

function t = layerClassOwnerTask(lqn, cls)
% LQN task a layer class belongs to, from its [element type, index] attribute;
% a reference-path class ([-1, index]) is told apart by its name suffix
t = -1;
att = cls.attribute;
if ~isnumeric(att) || numel(att) < 2 || ~(att(2) >= 1)
    return
end
switch att(1)
    case LayeredNetworkElement.TASK
        t = att(2);
    case {LayeredNetworkElement.ENTRY, LayeredNetworkElement.ACTIVITY}
        t = lqn.parent(att(2));
    case LayeredNetworkElement.CALL
        t = lqn.parent(lqn.callpair(att(2),1));
    case -1
        if endsWith(cls.name, {'.RefPath','.RefRet'})
            t = att(2);
        elseif endsWith(cls.name, {'.RefHop','.RefHopRet'})
            t = lqn.parent(att(2));
        elseif endsWith(cls.name, {'.RefGate','.RefGate.Aux','.RefResume'})
            t = lqn.parent(lqn.callpair(att(2),1));
        end
end
end
