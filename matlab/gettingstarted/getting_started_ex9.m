% https://github.com/imperial-qore/line-solver/wiki/Getting-started#example-9-studying-a-departure-process
[model,source,queue,sink,oclass] = gallery_merl1;
%% Block 4: solution
solver = SolverCTMC(model,'cutoff',150,'seed',23000);

sa = solver.sampleSysAggr(1e5);
ind = model.getNodeIndex(queue);

filtEvent = cellfun(@(c) c.node == ind && c.event == EventType.ID_DEP, sa.event);
interDepTimes = diff(cellfun(@(c) c.t, {sa.event{filtEvent}}));

% estimated squared coeff. of variation of departures
SCVdEst = var(interDepTimes)/mean(interDepTimes)^2

util = solver.getAvgUtil(); 
util = util(queue);
avgWaitTime = solver.getAvgWaitT();  % Waiting time excluding service 
avgWaitTime = avgWaitTime(queue);
SCVa = source.getArrivalProcess(oclass).getSCV();
svcRate = queue.getServiceProcess(oclass).getRate();
SCVs = queue.getServiceProcess(oclass).getSCV();
% Marshall's exact formula
SCVd = SCVa + 2*util^2*SCVs - 2*util*(1-util)*svcRate*avgWaitTime 
