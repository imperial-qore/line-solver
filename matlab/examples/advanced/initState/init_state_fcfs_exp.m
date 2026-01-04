clear node jobclass solver;

model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
jobclass{1} = ClosedClass(model, 'Class1', 5, node{2}, 0);

node{1}.setService(jobclass{1}, Exp(1));
node{2}.setService(jobclass{1}, Exp(0.7));

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

P = cell(K,K);
P{1} = circul(2);

model.link(P);
[Qt,Ut,Tt] = model.getTranHandles();
options = Solver.defaultOptions;
options.verbose=0;
options.samples=1e4;
options.stiff=true;
options.timespan = [0,40];

%% This part illustrates the execution of different solvers
solver={};
solver{end+1} = CTMC(model,options);
%solver{end+1} = JMT(model,options);
%solver{end+1} = SSA(model,options);
solver{end+1} = FLD(model,options);
%solver{end+1} = MVA(model,options);
dashing = {'-','+'};

%%
model.initDefault;
disp('Prior 1: prior all on default initialization')
disp('Initial state is:')
state=model.getState();
[state{1}(1,:),state{2}(1,:)]
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    [QNt,UNt,TNt] = solver{s}.getTranAvg(Qt,Ut,Tt);
    subplot(1,2,1);
    plot(QNt{2,1}.t,QNt{2,1}.metric,dashing{s}); hold on
    solver{s}.reset();
end
title('Prior on default state');
ylabel('Queue length - station 2, class 1');
ylim([3,5])
xlabel('Time t');
xlim(options.timespan)
legend('ctmc','fluid','Location','SouthEast')

%%
model.initFromMarginal([2;3]);
disp('Prior 2: prior all on first found state with given marginal')
disp('Initial state is:')
state=model.getState();
[state{1}(1,:),state{2}(1,:)]
for s=1:length(solver)
    solver{s}.reset();
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    [QNt_marg,UNt_marg,TNt_marg] = solver{s}.getTranAvg(Qt,Ut,Tt);
    subplot(1,2,2);
    plot(QNt_marg{2,1}.t,QNt_marg{2,1}.metric,dashing{s}); hold on
    solver{s}.reset();
end
title('Prior on state with 3 jobs in station 2');
ylabel('Queue length - station 2, class 1');
ylim([3,5])
xlabel('Time t');
xlim(options.timespan)
%legend('ctmc','fluid')
