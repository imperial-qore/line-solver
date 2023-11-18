% M/M/10/infinity open network with general topology
%% Model Definition
clear; clc;

model = Network('Example4');
% Block 1: nodes
source = Source(model, 'Source');
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS);
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS);
queue3 = Queue(model, 'Queue3', SchedStrategy.FCFS);
queue4 = Queue(model, 'Queue4', SchedStrategy.FCFS);
queue5 = Queue(model, 'Queue5', SchedStrategy.FCFS);
queue6 = Queue(model, 'Queue6', SchedStrategy.FCFS);
queue7 = Queue(model, 'Queue7', SchedStrategy.FCFS);
queue8 = Queue(model, 'Queue8', SchedStrategy.FCFS);
queue9 = Queue(model, 'Queue9', SchedStrategy.FCFS);
queue10 = Queue(model, 'Queue10', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');
% Block 2: job classes
oclass = OpenClass(model, 'Class1');
source.setArrival(oclass, Exp(3));
queue1.setService(oclass, Exp(4));
queue2.setService(oclass, Exp(3));
queue3.setService(oclass, Exp(3));
queue4.setService(oclass, Exp(1));
queue5.setService(oclass, Exp(3));
queue6.setService(oclass, Exp(2));
queue7.setService(oclass, Exp(3));
queue8.setService(oclass, Exp(4));
queue9.setService(oclass, Exp(1));
queue10.setService(oclass, Exp(3));
% Block 3: topology
model.addLinks([source, queue1; source, queue2; source, queue3; source, queue4; ...
                queue1, queue4; queue1, queue7; queue1, queue9; ...
                queue2, queue3; queue2, queue5; ...
                queue3, queue6; ...
                queue4, sink; ...
                queue5, queue6; queue5, queue7; ...
                queue6, queue8; ...
                queue7, queue10; ...
                queue8, sink; ...
                queue9, sink; ...
                queue10, sink]);

source.setProbRouting(oclass, queue1, 0.25);
source.setProbRouting(oclass, queue2, 0.25);
source.setProbRouting(oclass, queue3, 0.25);
source.setProbRouting(oclass, queue4, 0.25);

queue1.setProbRouting(oclass, queue4, 0.3);
queue1.setProbRouting(oclass, queue7, 0.4);
queue1.setProbRouting(oclass, queue9, 0.4);

queue2.setProbRouting(oclass, queue3, 0.5);
queue2.setProbRouting(oclass, queue5, 0.5);

queue5.setProbRouting(oclass, queue6, 0.5);
queue5.setProbRouting(oclass, queue7, 0.5);

% Model View
% model.jsimgView
% model.printRoutingMatrix
[~,H] = model.getGraph();
plot(H,'EdgeLabel',H.Edges.Weight,'Layout','Layered')

%% heuristic performance 

% Prob
SolverSSA(model).getAvgSysTable()

% JSQ
model.reset();
source.setRouting(oclass, RoutingStrategy.JSQ);
queue1.setRouting(oclass, RoutingStrategy.JSQ);
queue2.setRouting(oclass, RoutingStrategy.JSQ);
queue5.setRouting(oclass, RoutingStrategy.JSQ);

SolverSSA(model).getAvgSysTable()

%% RL control

% update model
model.reset();
source.setRouting(oclass, RoutingStrategy.JSQ);
queue1.setRouting(oclass, RoutingStrategy.JSQ);
queue2.setRouting(oclass, RoutingStrategy.JSQ);
queue5.setRouting(oclass, RoutingStrategy.JSQ);

% hyper-parameters
truncation = 10;
lrs = [0.1, 0.2];
eps = 1.0;
eps_decays = [0.9, 0.999];
episode = 5e4;
discount_factor = 0.95;

env_full = RLenvGeneral(model, [2,3,4,5,6,7,8,9,10,11], [1,2,3,6], truncation, discount_factor);
env_partial = RLenvGeneral(model, [2,3,4,5,6,7,8,9,10,11], 1, truncation, discount_factor);

quad_approx_coeff_full = containers.Map;
quad_approx_X_full = containers.Map;
quad_approx_Y_full = containers.Map;

quad_approx_coeff_partial = containers.Map;
quad_approx_X_partial = containers.Map;
quad_approx_Y_partial = containers.Map;

for i = 1:length(lrs)
    for j = 1:length(eps_decays)

        line_printf('running parameters: lr=%.2f, eps-decay=%.3f \n', lrs(i), eps_decays(j));
        keystr = strcat(num2str(i),num2str(j)); 

        td_full = TDAgentGeneral(lrs(i), eps, eps_decays(j));
        td_partial = TDAgentGeneral(lrs(i), eps, eps_decays(j));

        
        [tempX, tempY, tempCoeff] = td_full.solve_by_quad(env_full, episode);
        quad_approx_X_full(keystr) = tempX;
        quad_approx_Y_full(keystr) = tempY;
        quad_approx_coeff_full(keystr) = tempCoeff;

        [tempX, tempY, tempCoeff] = td_partial.solve_by_quad(env_partial, episode);
        quad_approx_X_partial(keystr) = tempX;
        quad_approx_Y_partial(keystr) = tempY;
        quad_approx_coeff_partial(keystr) = tempCoeff;
    
    end
end



%% RL performance

% full quad
coeff_best_full = quad_approx_coeff_full("12");
model.reset();
source.setRouting(oclass, RoutingStrategy.RL, coeff_best_full.', {[1,2,3,6], truncation+1});   % [1,2,3,6]: nodesNeedAction;  truncation+1: state size
queue1.setRouting(oclass, RoutingStrategy.RL, coeff_best_full.', {[1,2,3,6], truncation+1});
queue2.setRouting(oclass, RoutingStrategy.RL, coeff_best_full.', {[1,2,3,6], truncation+1});
queue5.setRouting(oclass, RoutingStrategy.RL, coeff_best_full.', {[1,2,3,6], truncation+1});

SolverSSA(model).getAvgSysTable

% partial quad
coeff_best_partial = quad_approx_coeff_partial("11");
model.reset();
source.setRouting(oclass, RoutingStrategy.RL, coeff_best_partial.', {[1], truncation+1});   % [1]: nodesNeedAction;  truncation+1: state size
queue1.setRouting(oclass, RoutingStrategy.RL, coeff_best_partial.', {[1], truncation+1});
queue2.setRouting(oclass, RoutingStrategy.RL, coeff_best_partial.', {[1], truncation+1});
queue5.setRouting(oclass, RoutingStrategy.RL, coeff_best_partial.', {[1], truncation+1});

SolverSSA(model).getAvgSysTable

%% runtime analysis for RL-learning and control
% comment all previous codes except mode-definition

% model.reset();
% source.setRouting(oclass, RoutingStrategy.JSQ);
% queue1.setRouting(oclass, RoutingStrategy.JSQ);
% queue2.setRouting(oclass, RoutingStrategy.JSQ);
% queue5.setRouting(oclass, RoutingStrategy.JSQ);
% 
% % hyper-parameters
% truncation = 10;
% lr = 0.1;
% eps = 1.0;
% eps_decay = 0.999;
% episode = 3e4;
% discount_factor = 0.95;
% 
% % env_full = RLenvGeneral(model, [2,3,4,5,6,7,8,9,10,11], [1,2,3,6], truncation, discount_factor);
% % td_full = TDAgentGeneral(lr, eps, eps_decay);
% % [X_full, Y_full, Coeff_full] = td_full.solve_by_quad(env_full, episode);
% 
% env_partial = RLenvGeneral(model, [2,3,4,5,6,7,8,9,10,11], 1, truncation, discount_factor);
% td_partial = TDAgentGeneral(lr, eps, eps_decay);
% [X_partial, Y_partial, Coeff_partial] = td_partial.solve_by_quad(env_partial, episode);
