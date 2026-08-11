classdef rl_env < handle
    properties
        model;                  % queue model
        gamma;                  % discount factor
        idxOfQueueInNodes;      % model.nodes{i} are queues for i in idxOfQueueInNodes
        idxOfSourceInNodes;     % model.nodes{i} are sources for i in idxOfQueueInNodes
        stateSize;             % [s_i <= m] are state space to be considered
        actionSize;            % actionSize: number of queues
    end
    
    methods
        function obj=rl_env(model, idxOfQueueInNodes, idxOfSourceInNodes, stateSize, gamma)
            obj.model = model.copy;
            obj.idxOfQueueInNodes = idxOfQueueInNodes;                      % how to distinguish sources & queues in model
            obj.idxOfSourceInNodes = idxOfSourceInNodes;
            obj.stateSize = stateSize;
            obj.gamma = gamma;
            obj.actionSize = length(idxOfQueueInNodes);
        end
     
        function r = isInStateSpace(obj, nodes)
            r = true;
            for i = obj.idxOfQueueInNodes
                if sum(nodes{i}.state) > obj.stateSize
                    r = false;
                    break
                end
            end
        end

        function r = isInActionSpace(obj, nodes)
            r = true;
            for i = obj.idxOfQueueInNodes
                if sum(nodes{i}.state) > obj.stateSize - 1
                    r = false;
                    break
                end
            end
        end

        function [t, depNode]=sample(obj)                                   % how to identify the type of new event
            solver = SolverSSA(obj.model, 'verbose', false);
            sample = solver.sampleSysAggr(1);
            t = sample.t;
            if sample.event{1}.event==EventType.DEP
                depNode = sample.event{1}.node;
            elseif sample.event{2}.event==EventType.DEP
                depNode = sample.event{2}.node;
            end
        end

        function update(obj, newState)                                      % how to update model after an event
            for i = 1:length(obj.idxOfQueueInNodes)                        % newState = [#jobsInQueue_i]
                obj.model.nodes{obj.idxOfQueueInNodes(i)}.state = State.fromMarginal(obj.model, obj.model.nodes{obj.idxOfQueueInNodes(i)}, newState(i));
            end
        end

        function reset(obj)
            obj.model.reset();
            obj.model.initDefault();
        end
    end
end
