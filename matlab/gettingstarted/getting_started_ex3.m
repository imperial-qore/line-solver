% https://github.com/imperial-qore/line-solver/wiki/Getting-started#example-3-machine-interference-problem
model = Network('MRP');
%% Block 1: nodes
delay = Delay(model,'WorkingState');
queue = Queue(model, 'RepairQueue', SchedStrategy.FCFS);
queue.setNumberOfServers(2);
%% Block 2: classes
cclass = ClosedClass(model, 'Machines', 3, delay);
delay.setService(cclass, Exp(0.5));
queue.setService(cclass, Exp(4.0));
%% Block 3: topology
model.link(Network.serialRouting(delay,queue));
%% Block 4: solution
solver = SolverCTMC(model);
ctmcAvgTable = solver.getAvgTable

StateSpace = solver.getStateSpace()
InfGen = full(solver.getGenerator())

model.printInfGen(InfGen,StateSpace)

[StateSpace,nodeStateSpace] = solver.getStateSpace()
nodeStateSpace{delay}
nodeStateSpace{queue}
