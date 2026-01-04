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

                    actions = env.actionSpace(depNode);                     % dep at node i, possible actions are [k,l,m]
                    
                   
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
                
            pointValues = containers.Map;                                   % hashmap value fn
            pointValues(num2str(zeros(1, env.nqueues))) = 0;
            pointValues('external')=0;

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

                    actions = env.actionSpace(depNode);                     % dep at node i, possible actions are [k,l,m]
                    
                    % create an exploit-explore policy
                    nextPointValues = zeros(1, length(actions));
                    for act_i = 1 : length(actions)
                        q_idx = find(env.idxOfQueueInNodes == actions(act_i));
                        tmp_next_state = x;
                        tmp_next_state(q_idx) = tmp_next_state(q_idx) + 1;
                        if pointValues.isKey(num2str(tmp_next_state)) 
                            nextPointValues(act_i) = pointValues(num2str(tmp_next_state));
                        else
                            nextPointValues(act_i) = pointValues('external');
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
                       
                    if ~pointValues.isKey(num2str(n))
                        pointValues(num2str(n)) = pointValues('external');
                    end
                    
                    if pointValues.isKey(num2str(x))
                        pointValues(num2str(n)) = (1-obj.lr)*pointValues(num2str(n)) + obj.lr*(c-t*mean_cost_rate + pointValues(num2str(x)));
                    else
                        pointValues(num2str(n)) = (1-obj.lr)*pointValues(num2str(n)) + obj.lr*(c-t*mean_cost_rate + pointValues('external'));
                    end

                    if sum(n)==0
                        substractor = pointValues(num2str(n));
                        for k = keys(pointValues)
                            pointValues(k{1}) = pointValues(k{1}) - substractor;
                        end
                    end

                    t = 0;
                    c = 0;
                    n = x;
                end
            
            end


            pointValues.remove('external');
            X = zeros(pointValues.Count, 1 + env.nqueues);
            Y = zeros(pointValues.Count, 1);
            iterator = 1;
            for k = keys(pointValues)
                X(iterator, :) = [1 str2num(k{1})];
                Y(iterator, :) = pointValues(k{1});
                iterator = iterator + 1;
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

