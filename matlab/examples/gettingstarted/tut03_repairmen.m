% Example 3: Machine interference problem
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
solver = CTMC(model);
ctmcAvgTable = solver.avgTable()

StateSpace = solver.stateSpace()
InfGen = full(solver.generator())

CTMC.printInfGen(InfGen,StateSpace)

[StateSpace,nodeStateSpace] = solver.stateSpace()
nodeStateSpace{model.getNodeIndex(delay)}  % delay node state space
nodeStateSpace{model.getNodeIndex(queue)}  % queue node state space
