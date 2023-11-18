%% M/M/4/infinity open network 
% Model Definition
clear; clc;

model = Network('Example2');
% Block 1: nodes
source = Source(model, 'Source');
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS);
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS);
queue3 = Queue(model, 'Queue3', SchedStrategy.FCFS);
queue4 = Queue(model, 'Queue4', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');
% Block 2: job classes
oclass = OpenClass(model, 'Class1');
source.setArrival(oclass, Exp(4));
queue1.setService(oclass, Exp(1));
queue2.setService(oclass, Exp(2));
queue3.setService(oclass, Exp(3));
queue4.setService(oclass, Exp(4));
% Block 3: topology
model.addLinks([source, queue1; source, queue2; source, queue3; source, queue4;...
                queue1, sink; queue2, sink; queue3, sink; queue4, sink]);
source.setProbRouting(oclass, queue1, 0.25);
source.setProbRouting(oclass, queue2, 0.25);
source.setProbRouting(oclass, queue3, 0.25);
source.setProbRouting(oclass, queue4, 0.25);

% Model View
% model.jsimgView
% model.printRoutingMatrix
[~,H] = model.getGraph();
plot(H,'EdgeLabel',H.Edges.Weight,'Layout','Layered')

%% RL control

% update model
model.reset();
source.setRouting(oclass, RoutingStrategy.JSQ);

% hyper-parameters
truncation = 6;
lrs = [0.1, 0.15, 0.2];
eps = 1.0;
eps_decays = [0.9, 0.99, 0.999];
episode = 2e4;
discount_factor = 0.95;
repeat_number = 5;

env = RLenvGeneral(model, [2,3,4,5], 1, truncation, discount_factor);
value_fn_control = containers.Map;
quad_approx_coeff = containers.Map;
quad_approx_X = containers.Map;
quad_approx_Y = containers.Map;

for i = 1:length(lrs)
    for j = 1:length(eps_decays)

        line_printf('running parameters: lr=%.2f, eps-decay=%.3f \n', lrs(i), eps_decays(j));
        keystr = strcat(num2str(i),num2str(j)); 

        td = TDAgentGeneral(lrs(i), eps, eps_decays(j));
        value_fn_control(keystr) = zeros(truncation+1, truncation+1);
        for repeat = 1:repeat_number
            temp_v = td.solve(env, episode);
            value_fn_control(keystr) = value_fn_control(keystr) + temp_v;
        end
        value_fn_control(keystr) = value_fn_control(keystr) / repeat_number;
        
        [tempX, tempY, tempCoeff] = td.solve_by_quad(env, episode);
        quad_approx_X(keystr) = tempX;
        quad_approx_Y(keystr) = tempY;
        quad_approx_coeff(keystr) = tempCoeff;
    
    end
end



%% RL performance 

% tabular
value_fn_control_best = value_fn_control("13");
model.reset();
source.setRouting(oclass, RoutingStrategy.RL, value_fn_control_best, {[1], 0});     % [1]: nodesNeedAction;  0: tabular value fn
SolverSSA(model).getAvgTable
SolverSSA(model).getAvgSysTable

% quad
coeff_best = quad_approx_coeff("13");
model.reset();
source.setRouting(oclass, RoutingStrategy.RL, coeff_best.', {[1], truncation+1});   % [1]: nodesNeedAction;  truncation+1: state size
SolverSSA(model).getAvgTable
SolverSSA(model).getAvgSysTable

%% heuristics performance

% equal prob
SolverSSA(model).getAvgTable
SolverSSA(model).getAvgSysTable


% JSQ
model.reset();
source.setRouting(oclass, RoutingStrategy.JSQ);
SolverSSA(model).getAvgTable
SolverSSA(model).getAvgSysTable

% Round Robin
model.reset();
source.setRouting(oclass, RoutingStrategy.RROBIN);
SolverSSA(model).getAvgTable
SolverSSA(model).getAvgSysTable


% Power-of-K
model.reset();
source.setRouting(oclass, RoutingStrategy.KCHOICES, 2);  % k = 2
SolverJMT(model).getAvgTable
SolverJMT(model).getAvgSysTable

