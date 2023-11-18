jmodel = JNetwork('MRP');
%% Block 1: nodes
delay = Delay(jmodel,'WorkingState');
queue = Queue(jmodel, 'RepairQueue', SchedStrategy.FCFS);
queue.setNumberOfServers(2);
%% Block 2: classes
cclass = ClosedClass(jmodel, 'Machines', 3, delay);
delay.setService(cclass, Exp(0.5));
queue.setService(cclass, Exp(4.0));
%% Block 3: topology
jmodel.link(Network.serialRouting(delay,queue));
%% Block 4: solution
%solver = SolverCTMC(model);
%ctmcAvgTable = solver.getAvgTable

%StateSpace = solver.getStateSpace()
%InfGen = full(solver.getGenerator())

%model.printInfGen(InfGen,StateSpace)

%[StateSpace,nodeStateSpace] = solver.getStateSpace()
%nodeStateSpace{delay}
%nodeStateSpace{queue}
