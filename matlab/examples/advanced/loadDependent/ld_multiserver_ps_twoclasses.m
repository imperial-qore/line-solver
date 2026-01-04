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

msT=MVA(model,'exact').getAvgTable
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

lldAvgTableCTMC=CTMC(ldmodel).getAvgTable %exact

lldAvgTableNC=NC(ldmodel).getAvgTable %exact
lldAvgTableRD=NC(ldmodel,'method','rd').getAvgTable
lldAvgTableNRP=NC(ldmodel,'method','nrp').getAvgTable
lldAvgTableNRL=NC(ldmodel,'method','nrl').getAvgTable

lldAvgTableMVALD=MVA(ldmodel,'method','exact').getAvgTable
lldAvgTableQD=MVA(ldmodel,'method','qd').getAvgTable

lldAvgTableJMT=JMT(ldmodel,'seed',23000,'samples',5000).getAvgTable
