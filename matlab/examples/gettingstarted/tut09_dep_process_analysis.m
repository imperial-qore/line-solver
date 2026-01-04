% Example 9: Studying a departure process
[model,source,queue,sink,oclass] = gallery_merl1;
%% Block 4: solution
solver = CTMC(model,'cutoff',150);

sa = solver.sampleSysAggr(5e3);
ind = model.getNodeIndex(queue);

filtEvent = cellfun(@(c) c.node == ind && (isequal(c.event, EventType.DEP) || (isnumeric(c.event) && c.event == 2)), sa.event);
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

% Calculate relative error between simulated and theoretical SCV
relativeError = abs(SCVdEst - SCVd) / SCVd * 100;
fprintf('\n=== Departure Process Analysis Results ===\n');
fprintf('Simulated SCV of departures:   %.6f\n', SCVdEst);
fprintf('Theoretical SCV (Marshall):    %.6f\n', SCVd);
fprintf('Relative error:                %.2f%%\n', relativeError); 
