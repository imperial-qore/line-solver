classdef rl_td_agent < handle
    properties
        v;                      % value function
        Q;                      % Q function
        vSize;                  % size of value function
        QSize;                  % size of Q function
        epsilon = 1;            % explore-exploit rate
        eps_decay = 0.99;       % explore-exploit rate decay
        lr = 0.05;              % learning rate
    end

    methods
        function obj=rl_td_agent(lr, eps, epsDecay)
            obj.lr = lr;
            obj.epsilon = eps;
            obj.eps_decay = epsDecay;
            obj.v = 0; 
            obj.vSize = 0;
            obj.Q = 0; 
            obj.QSize = 0;
        end

        function reset(obj, env)
            obj.v = 0; 
            obj.vSize = 0;
            obj.Q = 0; 
            obj.QSize = 0;
            env.reset();
        end
        
        function v = getValueFunction(obj)
            v = obj.v;
        end
        
        function Q = getQFunction(obj)
            Q = obj.Q;
        end

%         function p = getPolicy(obj)
%             p = zeros(size(obj.v));
%             
%         end

        function solve(obj, env)
            obj.reset(env);
            
            obj.v = zeros((zeros(1, env.actionSize)+env.stateSize + 5));                     % value function
            obj.Q = rand([(zeros(1, env.actionSize)+env.stateSize + 5), env.actionSize]);    % Q function
            obj.vSize = size(obj.v);
            obj.QSize = size(obj.Q);

            x = zeros(1, env.actionSize);       % initial state
            n = zeros(1, env.actionSize);       % initial previous state
            % t_prev = 0;                         % time of last event
            t = 0;                              % time of current event
            % dt = 0;                             % time period between two successive events
            c = 0;                              % incurred costs between the visits
            T = 0;                              % total discounted elapsed time
            C = 0;                              % total discounted costs

            num_episodes = 1e4;
            eps = obj.epsilon;
            j = 0;

            while j < num_episodes
                if mod(j, 1e3)==0
                    line_printf(mfilename,sprintf('running episode #%d .\n',j));
                end
                
                % if mod(j, 100) == 0
                eps = eps * obj.eps_decay;
                % end
                            
                % t_prev = t;
                [dt, depNode] = env.sample();                                       % how to successive sampling
                t = dt + t;
                % t = dt + t_prev;
                
                c = c + sum(x) * dt;
            
                if ismember(depNode, env.idxOfSourceInNodes)     % new job                           
                    if env.isInActionSpace(env.model.nodes)
                        % create an exploit-explore policy
                        next_locs = zeros(env.actionSize, env.actionSize) + x + 1 + eye(env.actionSize);
                        next_states = obj.get_state_from_locs(obj.vSize, next_locs);
                        policy = obj.createGreedyPolicy(obj.v(next_states), eps, env.actionSize);
                        
                        action = sum(rand >= cumsum([0, policy]));                                  
                    else
                        action = find(x==min(x));               % JSQ
                        if length(action)>1
                            action = randomsample(action, 1);
                        end  
                    end
            
                    x(action) = x(action) + 1;
                    env.update(x);
                    % env.model.nodes{action+1}.state = State.fromMarginal(model, model.stations{action+1}, sum(model.stations{action+1}.state)+1);
                    % State.afterEvent(model, queue1, queue1.space, sample.event{1}, 1)
                    % State.afterEvent(sn, ind, inspace, event, class, isSimulation)
                    
                elseif ismember(depNode, env.idxOfQueueInNodes) % dep from Queue, idx: Node{depNode}
                    x(env.idxOfQueueInNodes == depNode) = max(0, x(env.idxOfQueueInNodes == depNode) - 1);
                    env.update(x);
                    % State.afterEvent(sn, ind, inspace, event, class, isSimulation)
                end
                       
                if env.isInStateSpace(env.model.nodes)
                    j = j + 1;
                    T = env.gamma * T + t;
                    C = env.gamma * C + c;
                    mean_cost_rate = C/T;
                       
                    prev_state = obj.get_state_from_loc(obj.vSize, n+1);
                    cur_state = obj.get_state_from_loc(obj.vSize, x+1);
                    obj.v(prev_state) = (1-obj.lr)*obj.v(prev_state) + obj.lr*(c - t*mean_cost_rate + obj.v(cur_state)); % here "obj.v(cur_state) * env.gamma" ?
                    obj.v = obj.v - obj.v(1);
                    
                    t = 0;
                    c = 0;
                    n = x;
                end
            
            end
        end
        
        function s = get_state_from_locs(obj, objSize, locs)
            s = zeros(1, size(locs,1));
            for i=1:size(locs,1)
                s(i) = obj.get_state_from_loc(objSize, locs(i,:));
            end
        end
    end

    methods(Static)
        function policy = createGreedyPolicy(state_Q, epsilon, nA)
             policy = ones(1, nA) * epsilon / nA;
             argmin = find(state_Q-min(state_Q)<GlobalConstants.FineTol);
             policy(argmin) = policy(argmin) + (1-epsilon)/length(argmin);
        end

        function s = get_state_from_loc(objSize, loc)
            s = 0;
            if size(objSize,2) == size(loc, 2)
                for i=1:size(objSize,2)
                    if i==1
                        s = s + loc(i);
                    else
                        s = s + (loc(i)-1) * prod(objSize(1:(i-1)));
                    end
                end
            end
        end

    end

end


