classdef rl_env_general < handle
    
    
    properties
        model;                                                              % network model
        gamma;                                                              % discount factor (e.g. = 0.95)
        idxOfQueueInNodes;                                                  % model.nodes{i} are queues for i in idxOfQueueInNodes
        nqueues;                                                            % number of queues
        idxOfActionNodes;                                                   % model.nodes{i} are nodes where actions are needed, for i in idxOfQueueInNodes
        stateSize;                                                          % Integer: maximal jobs in a queue where policy can be learned; exceeding will result in JSQ policy
        actionSpace;                                                        % {key: [idx1, ..., idxn]}; nodes{key} has a dispatch policy to node{idxi}
    end
    
    methods
        function obj = rl_env_general(model, idxOfQueueInNodes, idxOfActionNodes, stateSize, gamma) 
         
            obj.model = model.copy;
            obj.idxOfQueueInNodes = idxOfQueueInNodes;                      
            obj.nqueues = length(idxOfQueueInNodes);
            obj.idxOfActionNodes = idxOfActionNodes;                        
            obj.stateSize = stateSize;                                      
            obj.gamma = gamma;                                              

            obj.actionSpace = containers.Map('KeyType', 'int32', 'ValueType','any');
            for i = idxOfActionNodes
                obj.actionSpace(i) = find(model.connections(i,:)==1);
            end
        end
     
        function r = isInStateSpace(obj, state)
            if length(obj.idxOfQueueInNodes) ~= length(state)
                line_error(mfilename,' states size not match. ( state size=%i, required size=%i ) \n', length(state), length(obj.idxOfQueueInNodes));
            end

            r = true;
            for i = 1:length(obj.idxOfQueueInNodes)
                if state(i) > obj.stateSize
                    r = false;
                    break
                end
            end
        end

%         function r = isInStateSpace(obj, nodes)
%             r = true;
%             for i = obj.idxOfQueueInNodes
%                 if sum(nodes{i}.state) > obj.stateSize
%                     r = false;
%                     break
%                 end
%             end
%         end

        function r = isInActionSpace(obj, state)
            if length(obj.idxOfQueueInNodes) ~= length(state)
                line_error(mfilename,' states size not match. ( state size=%i, required size=%i ) \n', length(state), length(obj.idxOfQueueInNodes));
            end
            
            r = true;
            for i = 1:length(obj.idxOfQueueInNodes)
                if state(i) > obj.stateSize - 1
                    r = false;
                    break
                end
            end
        end

%         function r = isInActionSpace(obj, nodes)
%             r = true;
%             for i = obj.idxOfQueueInNodes
%                 if sum(nodes{i}.state) > obj.stateSize - 1
%                     r = false;
%                     break
%                 end
%             end
%         end

        function [dt, depNode, arvNode, sample]=sample(obj)                 
            solver = SolverSSA(obj.model, 'verbose', false);
            sample = solver.sampleSysAggr(1);
            dt = sample.t;
            switch sample.event{1}.event
                case EventType.DEP
                    depNode = sample.event{1}.node;
                case EventType.ARV
                    arvNode = sample.event{1}.node;
                otherwise
                    line_error(mfilename,' Not known event type in sample. ( Event1.event=%i ) \n', sample.event{1}.event);
            end
            switch sample.event{2}.event
                case EventType.DEP
                    depNode = sample.event{2}.node;
                case EventType.ARV
                    arvNode = sample.event{2}.node;
                otherwise
                    line_error(mfilename,' Not known event type in sample. ( Event2.event=%i ) \n', sample.event{2}.event);
            end
        end

%         function update(obj, newState)                                      
%             for i = 1:obj.nqueues                                           % newState = [#jobsInQueue_i]
%                 tmp = newState(i);
%                 q_idx = obj.idxOfQueueInNodes(i);
%                 obj.model.nodes{q_idx}.state = State.fromMarginal(obj.model, obj.model.nodes{q_idx}, tmp);
%             end
%         end

        function update(obj, sample)
            for i = 1:length(sample.event)
                [obj.model.nodes{sample.event{i}.node}.state, ~] = State.afterEvent(obj.model.getStruct, sample.event{i}.node, obj.model.nodes{sample.event{i}.node}.state, sample.event{i}.event, sample.event{1}.class, true); % isSimulation = true or false, both ok
            end
        end

        function reset(obj)
            obj.model.reset();
            obj.model.initDefault();
        end
    end
end

