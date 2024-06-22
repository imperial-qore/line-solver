%% M/M/4/infinity open network with general topology
% Model Definition
clear; clc;

model = Network('Example3');
% Block 1: nodes
source = Source(model, 'Source');
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS);
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS);
queue3 = Queue(model, 'Queue3', SchedStrategy.FCFS);
queue4 = Queue(model, 'Queue4', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');
% Block 2: job classes
oclass = OpenClass(model, 'Class1');
source.setArrival(oclass, Exp(3));
queue1.setService(oclass, Exp(2));
queue2.setService(oclass, Exp(3));
queue3.setService(oclass, Exp(3));
queue4.setService(oclass, Exp(1));
% Block 3: topology
P = model.initRoutingMatrix;
P{oclass} = [0, 0.5, 0.5, 0, 0, 0;     % row: source
             0, 0, 0, 1, 0, 0;         % row: queue1
             0, 0, 0, .5, .5, 0;       % row: queue2
             0, 0, 0, 0, 0, 1;         % row: queue3
             0, 0, 0, 0, 0, 1;         % row: queue4
             0, 0, 0, 0, 0, 0];        % row: sink
model.link(P);


% Model View
% model.jsimgView
% model.printRoutingMatrix
[~,H] = model.getGraph();
plot(H,'EdgeLabel',H.Edges.Weight,'Layout','Layered')

%% heuristics performance

% Prob
SolverSSA(model).getAvgTable()
SolverSSA(model).getAvgSysTable()

% JSQ
model.reset();
source.setRouting(oclass, RoutingStrategy.JSQ);
queue2.setRouting(oclass, RoutingStrategy.JSQ);

SolverSSA(model).getAvgTable()
SolverSSA(model).getAvgSysTable()

% RR
model.reset();
source.setRouting(oclass, RoutingStrategy.RROBIN);
queue2.setRouting(oclass, RoutingStrategy.RROBIN);

SolverSSA(model).getAvgTable()
SolverSSA(model).getAvgSysTable()

%% RL control

% update model
model.reset();
source.setRouting(oclass, RoutingStrategy.JSQ);
queue2.setRouting(oclass, RoutingStrategy.JSQ);

% hyper-parameters
truncation = 10;
lrs = [0.1, 0.15, 0.2];
eps = 1.0;
eps_decays = [0.9, 0.99, 0.999];
episode = 3e4;
discount_factor = 0.95;
repeat_number = 3;

env = RLenvGeneral(model, [2,3,4,5], [1,3], truncation, discount_factor);
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
source.setRouting(oclass, RoutingStrategy.RL, value_fn_control_best, {[1,3], 0});   % [1,3]: nodesNeedAction;  0: tabular value fn
queue2.setRouting(oclass, RoutingStrategy.RL, value_fn_control_best, {[1,3], 0});

SolverSSA(model).getAvgTable
SolverSSA(model).getAvgSysTable

% quad
coeff_best = quad_approx_coeff("13");
model.reset();
source.setRouting(oclass, RoutingStrategy.RL, coeff_best.', {[1,3], truncation+1});   % [1,3]: nodesNeedAction;  truncation+1: state size
queue2.setRouting(oclass, RoutingStrategy.RL, coeff_best.', {[1,3], truncation+1});

SolverSSA(model).getAvgTable
SolverSSA(model).getAvgSysTable

%% runtime analysis for RL-learning and control

% % comment all previous codes except mode-definition
% 
% model.reset();
% source.setRouting(oclass, RoutingStrategy.JSQ);
% queue2.setRouting(oclass, RoutingStrategy.JSQ);
% 
% % hyper-parameters
% truncation = 10;
% lr = 0.1;
% eps = 1.0;
% eps_decay = 0.999;
% episode = 3e4;
% discount_factor = 0.95;
% 
% env = RLenvGeneral(model, [2,3, 4, 5], [1,3], truncation, discount_factor);
% td = TDAgentGeneral(lr, eps, eps_decay);
% 
% 
% % valueFnControl = td.solve(env, episode);
% [X, Y, Coeff] = td.solve_by_quad(env, episode);