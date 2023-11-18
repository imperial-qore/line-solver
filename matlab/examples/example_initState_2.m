clear node jobclass solver;

model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
jobclass{1} = ClosedClass(model, 'Class1', 3, node{2}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 2, node{2}, 0);
node{2}.setNumServers(3);

node{1}.setService(jobclass{1}, Exp(1));
node{1}.setService(jobclass{2}, Exp(1));
node{2}.setService(jobclass{1}, Exp(1.2));
node{2}.setService(jobclass{2}, Erlang.fitMeanAndSCV(1.0,0.5));

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

P = cell(K,K);
P{1,1} = [0.3,0.1; 0.2,0];
P{1,2} = [0.6,0; 0.8,0];
P{2,2} = [0,1; 0,0];
P{2,1} = [0,0; 1,0];

model.link(P);
[Qt,Ut,Tt] = model.getTranHandles();
options = Solver.defaultOptions;
options.verbose=1;
options.samples=1e4;
options.stiff=true;
options.timespan = [0,5];
%% This part illustrates the execution of different solvers
solver{1} = SolverCTMC(model,options);
%solver{end+1} = SolverJMT(model,options);
%solver{end+1} = SolverSSA(model,options);
solver{end+1} = SolverFluid(model,options);
%solver{end+1} = SolverMVA(model,options);
dashing = {'-','--'};

%%
model.initDefault;
disp('Prior 1: prior all on default initialization')
disp('Initial state is:')
state=model.getState();
[state{1}(1,:),state{2}(1,:)]
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    [QNt,UNt,TNt] = solver{s}.getTranAvg(Qt,Ut,Tt);
    subplot(3,1,1);
    plot(QNt{2,1}.t,QNt{2,1}.metric,dashing{s}); hold all
    solver{s}.reset();
end
title('Prior on default state');
ylabel('QLen- station 2, class 1');
ylim([0,5])
xlabel('Time t');
xlim(options.timespan)
legend('ctmc','fluid')
%%
model.initFromMarginal([0,0;4,1]);
disp('Prior 2: prior all on first found state with the same number of jobs')
disp('Initial state is:')
state=model.getState();
[state{1}(1,:),state{2}(1,:)]
for s=1:length(solver)
    solver{s}.reset();
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    [QNt_marg,UNt_marg,TNt_marg] = solver{s}.getTranAvg(Qt,Ut,Tt);
    subplot(3,1,2);
    plot(QNt_marg{2,1}.t,QNt_marg{2,1}.metric,dashing{s}); hold all
    solver{s}.reset();
end
title('Prior on first state with the same number of jobs');
ylabel('QLen- station 2, class 1');
ylim([0,5])
xlabel('Time t');
xlim(options.timespan)
legend('ctmc','fluid')
%%
disp('Prior 3: uniform prior over all states with the same number of jobs')
model.initFromMarginal([0,0;4,1]);
disp('Initial states are:')
state=model.getState();
[repmat(state{1}(1,:),size(state{2},1),1),state{2}]
prior = node{2}.getStatePrior;
prior = 0*prior; prior=ones(size(prior))/length(prior);
node{2}.setStatePrior(prior);
for s=1:length(solver)
    solver{s}.reset();
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    [QNt_unif,UNt_unif,TNt_unif] = solver{s}.getTranAvg(Qt,Ut,Tt);
    subplot(3,1,3);
    plot(QNt_unif{2,1}.t,QNt_unif{2,1}.metric,dashing{s}); hold all
end
title('Uniform prior on states with the same number of jobs');
ylabel('QLen- station 2, class 1');
ylim([0,5])
xlabel('Time t');
xlim(options.timespan)
legend('ctmc','fluid')
hold off
