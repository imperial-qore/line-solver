clear node jobclass

N = 4; % number of jobs
c = 2; % number of servers
%%
model = Network('model');
node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
jobclass{1} = ClosedClass(model, 'Class1', N, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', N/2, node{1}, 0);
node{1}.setService(jobclass{1}, Exp.fitMean(1.0)); % mean = 1
node{1}.setService(jobclass{2}, Exp.fitMean(2.0)); % mean = 1
node{2}.setService(jobclass{1}, Exp.fitMean(1.5)); % mean = 1.5
node{2}.setService(jobclass{2}, Exp.fitMean(2.5)); % mean = 1.5
node{2}.setNumberOfServers(c);

P = model.initRoutingMatrix();
P{1,1} = model.serialRouting(node);
P{2,2} = model.serialRouting(node);
model.link(P);

msT=SolverMVA(model,'exact').getAvgTable
%%
ldmodel = Network('ldmodel');
node{1} = Delay(ldmodel, 'Delay');
node{2} = Queue(ldmodel, 'Queue1', SchedStrategy.PS);
jobclass{1} = ClosedClass(ldmodel, 'Class1', N, node{1}, 0);
jobclass{2} = ClosedClass(ldmodel, 'Class2', N/2, node{1}, 0);
node{1}.setService(jobclass{1}, Exp.fitMean(1.0)); % mean = 1
node{1}.setService(jobclass{2}, Exp.fitMean(2.0)); % mean = 1
node{2}.setService(jobclass{1}, Exp.fitMean(1.5)); % mean = 1.5
node{2}.setService(jobclass{2}, Exp.fitMean(2.5)); % mean = 1.5
node{2}.setLoadDependence(min(1:(N+N/2),c)); % multi-server with c servers

P = ldmodel.initRoutingMatrix();
P{1,1} = ldmodel.serialRouting(node);
P{2,2} = ldmodel.serialRouting(node);
ldmodel.link(P);

lldAvgTableCTMC=SolverCTMC(ldmodel).getAvgTable %exact

lldAvgTableNC=SolverNC(ldmodel).getAvgTable %exact
lldAvgTableRD=SolverNC(ldmodel,'method','rd').getAvgTable
lldAvgTableNRP=SolverNC(ldmodel,'method','nr.probit').getAvgTable
lldAvgTableNRL=SolverNC(ldmodel,'method','nr.logit').getAvgTable

lldAvgTableMVALD=SolverMVA(ldmodel,'method','exact').getAvgTable
lldAvgTableQD=SolverMVA(ldmodel,'method','qd').getAvgTable
