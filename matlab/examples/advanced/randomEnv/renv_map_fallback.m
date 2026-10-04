% Random-environment fallback for MAP/MMPP models
%
% A solver that cannot consume a non-renewal (MAP/MMPP/MMAP) process does not
% reject the model any longer: NetworkSolver intercepts it, replaces every
% modulating chain by a set of random-environment stages in which the process
% is exponential with its phase-conditional intensity (map2renv), and solves
% the stages with the same solver through SolverENV.
%
% The interception is method-aware: a method that handles the process natively
% (MVA 'rqna', MAM, CTMC, FLD, SSA, JMT, LDES) runs unchanged.

clear;

%% Closed network with an MMPP2 server: MVA and NC have no native MAP support
model = Network('mmppClosed');
delay = Delay(model,'Think');
queue = Queue(model,'Q1',SchedStrategy.FCFS);
jobclass = ClosedClass(model,'C1',5,delay);
delay.setService(jobclass, Exp(1.0));
queue.setService(jobclass, MMPP2(1.0,10.0,0.2,0.3)); % slow phase 1, fast phase 10
model.link(Network.serialRouting(delay,queue));

%% The environment image built for the solver
[renv, info] = map2renv(model);
fprintf('Random-environment image: %d stages, phase orders %s, MMPP image: %d\n', ...
    info.nstages, mat2str(info.orders), info.isMMPP);
renv.getStageTable()

%% MVA through the fallback, CTMC natively (exact)
AvgTableMVA = MVA(model).getAvgTable()
AvgTableCTMC = CTMC(model).getAvgTable()

%% The two environment limits, forced
options = MVA.defaultOptions;
options.config.map_env_method = 'dec'; % quasi-stationary (slow environment)
AvgTableDec = MVA(model,options).getAvgTable()
options.config.map_env_method = 'avg'; % rate-averaged (fast environment)
AvgTableAvg = MVA(model,options).getAvgTable()

%% Opting out restores the feature rejection
options = MVA.defaultOptions;
options.config.map_env = 'off';
try
    MVA(model,options).getAvgTable();
catch ME
    fprintf('map_env=''off'': %s\n', ME.message);
end
