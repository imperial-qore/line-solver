classdef rl_td_agent_general < handle                                            % class for TD learning and TD control

    properties
        v;                                                                  % value function
        vSize;                                                              % size of value function
        epsilon = 1;                                                        % explore-exploit rate
        eps_decay = 0.9999;                                                 % explore-exploit rate decay
        lr = 0.1;                                                           % learning rate
    end
    
    methods
        function obj = rl_td_agent_general(lr, eps, epsDecay)                   
            obj.lr = lr;
            obj.epsilon = eps;
            obj.eps_decay = epsDecay;
            obj.v = 0; 
            obj.vSize = 0;
        end
           
        function reset(obj, env)
            obj.v = 0; 
            obj.vSize = 0;
            env.reset();
        end
        
        function v = getValueFunction(obj)
            v = obj.v;
        end
      


        % TD learning for value function with heuristic routing strategy
        function v = solve_for_fixed_policy(obj, env, num_episodes)         % num_epsiodes = 10^4 ususally
            
            obj.reset(env);
            
            obj.v = zeros((zeros(1, env.nqueues)+env.stateSize + 1));       % value function
            obj.vSize = size(obj.v);
            
            t = 0;                                                          % time of current event
            c = 0;                                                          % incurred costs between the visits
            T = 0;                                                          % total discounted elapsed time
            C = 0;                                                          % total discounted costs
            x = zeros(1, env.nqueues);                                      % initial state
            n = zeros(1, env.nqueues);                                      % initial previous state
            
            
            j = 0;
            while j < num_episodes
                if mod(j, 1e3)==0
                    line_printf('running episode #%d \n',j);
                end
                       
                [dt, depNode, arvNode, sample] = env.sample();             
                t = dt + t;
                c = c + sum(x) * dt;
                
                if ismember(depNode, env.idxOfQueueInNodes)                 % Event involves departure from server i
                    depServer = find(env.idxOfQueueInNodes == depNode);
                    x(depServer) = x(depServer) - 1;
                end 

                if ismember(arvNode, env.idxOfQueueInNodes)                 % Event involves Arrival at server j
                    arvServer = find(env.idxOfQueueInNodes == arvNode);
                    x(arvServer) = x(arvServer) + 1;
                end
                
                env.update(sample);

                if env.isInStateSpace(x)
                    j = j + 1;
                    T = env.gamma * T + t;
                    C = env.gamma * C + c;
                    mean_cost_rate = C/T;
                    
                    prev_state = num2cell(n+1);                                                                                 % obj.get_state_from_loc(obj.vSize, n+1);
                    cur_state = num2cell(x+1);                                                                                  % obj.get_state_from_loc(obj.vSize, x+1);
                    obj.v(prev_state{:}) = (1-obj.lr)*obj.v(prev_state{:}) + obj.lr*(c - t*mean_cost_rate + obj.v(cur_state{:}));  %(1-obj.lr)*obj.v(prev_state) + obj.lr*(c - t*mean_cost_rate + obj.v(cur_state)); 
                    obj.v = obj.v - obj.v(1);

                    t = 0;
                    c = 0;
                    n = x;
                end
            
            end

            v = obj.v;
        end        



        % TD Control with Tabular value function
        function value_function = solve(obj, env, num_episodes)             % num_epsiodes = 10^4 ususally
            obj.reset(env);
            
            obj.v = zeros((zeros(1, env.nqueues)+env.stateSize + 1));       % value function
            obj.vSize = size(obj.v);

            t = 0;                                                          % time of current event
            c = 0;                                                          % incurred costs between the visits
            T = 0;                                                          % total discounted elapsed time
            C = 0;                                                          % total discounted costs
            x = zeros(1, env.nqueues);                                      % initial state
            n = zeros(1, env.nqueues);                                      % initial previous state

            eps = obj.epsilon;
            
            j = 0;
            while j < num_episodes
                if mod(j, 1e3)==0
                    line_printf('running episode #%d .\n',j);
                end
                
                
                eps = eps * obj.eps_decay;
                            
                [dt, depNode, arvNode, sample] = env.sample();                                                             
                t = dt + t;
                c = c + sum(x) * dt;
            
                if ismember(depNode, env.idxOfQueueInNodes)                 % Event involves departure from server i
                    depServer = find(env.idxOfQueueInNodes == depNode);
                    x(depServer) = max(0, x(depServer) - 1);
                end 

                if ismember(depNode, env.idxOfActionNodes) && env.isInActionSpace(x) % actions wanted at server i, and in action space

                    actions = env.actionSpace{depNode};                     % dep at node i, possible actions are [k,l,m]
                    
                   
                    % create an exploit-explore policy
                    next_values = obj.gen_next_values(env, x, actions);
                    policy = obj.createGreedyPolicy(next_values, eps, length(actions));
                    
                    arvNode = actions(sum(rand >= cumsum([0, policy]))); 
                    
                    % update sample with new arvNode
                    for i = 1:length(sample.event)
                        if sample.event{i}.event == EventType.ARV
                            sample.event{i}.node = arvNode;
                            break;
                        end
                    end

                end                                                        
                
                % update current state and model
                x(env.idxOfQueueInNodes == arvNode) = x(env.idxOfQueueInNodes == arvNode) + 1;
                env.update(sample); 
                
                       
                if env.isInStateSpace(x)                                    % in State Space, update state value
                    j = j + 1;
                    T = env.gamma * T + t;
                    C = env.gamma * C + c;
                    mean_cost_rate = C/T;
                       
                    prev_state = num2cell(n+1);
                    cur_state = num2cell(x+1);
                    obj.v(prev_state{:}) = (1-obj.lr)*obj.v(prev_state{:}) + obj.lr*(c - t*mean_cost_rate + obj.v(cur_state{:}));                % here "obj.v(cur_state) * env.gamma" ?
                    obj.v = obj.v - obj.v(1);
                    
                    t = 0;
                    c = 0;
                    n = x;
                end
            
            end
            value_function = obj.v;
        end


    

        % TD Control with HashMap value fn 
        function [X, Y]=solve_by_hashmap(obj, env, num_episodes)
            obj.reset(env);
                
            % Value function stored as parallel arrays keyed by state string:
            % pvKeys{i} is the state string, pvVals(i) its estimated value.
            pvKeys = {};                                                    % state-string keys
            pvVals = [];                                                    % corresponding values
            [pvKeys, pvVals] = pv_set(pvKeys, pvVals, num2str(zeros(1, env.nqueues)), 0);
            [pvKeys, pvVals] = pv_set(pvKeys, pvVals, 'external', 0);

            t = 0;                                                          % time of current event
            c = 0;                                                          % incurred costs between the visits
            T = 0;                                                          % total discounted elapsed time
            C = 0;                                                          % total discounted costs
            x = zeros(1, env.nqueues);                                      % initial state
            n = zeros(1, env.nqueues);                                      % initial previous state

            eps = obj.epsilon;

            j = 0;
            while j < num_episodes
                if mod(j, 1e3)==0
                    line_printf('running episode #%d .\n',j);
                end
                
                % if mod(j, 100) == 0
                eps = eps * obj.eps_decay;
                % end
                            
                [dt, depNode, arvNode, sample] = env.sample();                                  
                t = dt + t;
                c = c + sum(x) * dt;

                if ismember(depNode, env.idxOfQueueInNodes)                 % Event involves departure from server i
                    depServer = find(env.idxOfQueueInNodes == depNode);
                    x(depServer) = max(0, x(depServer) - 1);
                end
                

                if ismember(depNode, env.idxOfActionNodes) && env.isInActionSpace(x) % actions wanted at server i, and in action space

                    actions = env.actionSpace{depNode};                     % dep at node i, possible actions are [k,l,m]
                    
                    % create an exploit-explore policy
                    nextPointValues = zeros(1, length(actions));
                    for act_i = 1 : length(actions)
                        q_idx = find(env.idxOfQueueInNodes == actions(act_i));
                        tmp_next_state = x;
                        tmp_next_state(q_idx) = tmp_next_state(q_idx) + 1;
                        if pv_iskey(pvKeys, num2str(tmp_next_state))
                            nextPointValues(act_i) = pv_get(pvKeys, pvVals, num2str(tmp_next_state));
                        else
                            nextPointValues(act_i) = pv_get(pvKeys, pvVals, 'external');
                        end
                    end
                    policy = obj.createGreedyPolicy(nextPointValues, eps, length(actions));
                    
                    arvNode = actions(sum(rand >= cumsum([0, policy])));
                    
                    % update sample
                    for i = 1:length(sample.event)
                        if sample.event{i}.event == EventType.ARV
                            sample.event{i}.node = arvNode;
                            break;
                        end
                    end

                end
                
                % update current state and model
                x(env.idxOfQueueInNodes == arvNode) = x(env.idxOfQueueInNodes == arvNode) + 1;   
                env.update(sample); 
                
                if env.isInStateSpace(x)                                    % in State Space, update state value
                    j = j + 1;
                    T = env.gamma * T + t;
                    C = env.gamma * C + c;
                    mean_cost_rate = C/T;
                       
                    if ~pv_iskey(pvKeys, num2str(n))
                        [pvKeys, pvVals] = pv_set(pvKeys, pvVals, num2str(n), pv_get(pvKeys, pvVals, 'external'));
                    end

                    if pv_iskey(pvKeys, num2str(x))
                        [pvKeys, pvVals] = pv_set(pvKeys, pvVals, num2str(n), (1-obj.lr)*pv_get(pvKeys, pvVals, num2str(n)) + obj.lr*(c-t*mean_cost_rate + pv_get(pvKeys, pvVals, num2str(x))));
                    else
                        [pvKeys, pvVals] = pv_set(pvKeys, pvVals, num2str(n), (1-obj.lr)*pv_get(pvKeys, pvVals, num2str(n)) + obj.lr*(c-t*mean_cost_rate + pv_get(pvKeys, pvVals, 'external')));
                    end

                    if sum(n)==0
                        substractor = pv_get(pvKeys, pvVals, num2str(n));
                        pvVals = pvVals - substractor;
                    end

                    t = 0;
                    c = 0;
                    n = x;
                end
            
            end


            [pvKeys, pvVals] = pv_remove(pvKeys, pvVals, 'external');
            X = zeros(length(pvKeys), 1 + env.nqueues);
            Y = zeros(length(pvKeys), 1);
            for iterator = 1:length(pvKeys)
                X(iterator, :) = [1 str2num(pvKeys{iterator})]; %#ok<ST2NM>
                Y(iterator, :) = pvVals(iterator);
            end
            
        end


        
        % TD control using linear value fn approximator:  
        % v(q1,q2,...,qn) = w1*q1 + w2*q2 + ... + wn*qn (linear fn)
        function [X, Y, coeff]=solve_by_linear(obj, env, num_episodes)
            [X, Y] = obj.solve_by_hashmap(env, num_episodes);

            coeff = regress(Y, X);    
        end


        % TD control using quadratic value fn approximator:  
        % v(q1,q2,...,qn) = sum_{i,j} w_{ij} * q_i * q_j (quadratic fn)
        function [X, Y, coeff]=solve_by_quad(obj, env, num_episodes)
            [X, Y] = obj.solve_by_hashmap(env, num_episodes);

            sizeX = size(X);
            for i = 2:sizeX(2)
                for j = i:sizeX(2)
                    X(:,end+1) = X(:,i).* X(:,j);
                end
            end

            coeff = regress(Y, X);
        end


        function values=gen_next_values(obj, env, cur_state, actions)       % cur_state = x
            values = zeros(1, length(actions));
            for act_i = 1 : length(actions)
                q_idx = find(env.idxOfQueueInNodes == actions(act_i));
                tmp_loc = cur_state + 1;
                tmp_loc(q_idx) = tmp_loc(q_idx) + 1;
                tmp_idx = num2cell(tmp_loc);
                values(act_i) = obj.v(tmp_idx{:});
            end
        end
        
    end


    methods(Static)
        function policy = createGreedyPolicy(state_Q, epsilon, nA)
             policy = ones(1, nA) * epsilon / nA;
             argmin = find(state_Q == min(state_Q));
             policy(argmin) = policy(argmin) + (1-epsilon)/length(argmin);
        end

    end
end

% Local helpers implementing a string-keyed value map as parallel arrays
% (keys is a cell of char keys, vals a numeric vector of the same length).
function tf = pv_iskey(keys, k)
tf = any(strcmp(keys, k));
end

function v = pv_get(keys, vals, k)
i = find(strcmp(keys, k), 1);
v = vals(i);
end

function [keys, vals] = pv_set(keys, vals, k, v)
i = find(strcmp(keys, k), 1);
if isempty(i)
    keys{end+1} = k;
    vals(end+1) = v;
else
    vals(i) = v;
end
end

function [keys, vals] = pv_remove(keys, vals, k)
i = find(strcmp(keys, k), 1);
if ~isempty(i)
    keys(i) = [];
    vals(i) = [];
end
end

