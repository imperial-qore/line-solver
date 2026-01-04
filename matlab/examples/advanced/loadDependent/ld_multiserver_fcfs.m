clear node jobclass

N = 16; % number of jobs
c = 2; % number of servers

model = Network('model');
node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
jobclass{1} = ClosedClass(model, 'Class1', N, node{1}, 0);
node{1}.setService(jobclass{1}, Exp.fitMean(1.0)); % mean = 1
node{2}.setService(jobclass{1}, Exp.fitMean(1.5)); % mean = 1.5
node{2}.setNumberOfServers(c);

model.link(model.serialRouting(node));

msT=NC(model).getAvgTable

%% casted with scaling function that depends on the total queue
ldmodel = Network('model');
node{1} = Delay(ldmodel, 'Delay');
node{2} = Queue(ldmodel, 'Queue1', SchedStrategy.FCFS);
jobclass{1} = ClosedClass(ldmodel, 'Class1', N, node{1}, 0);
node{1}.setService(jobclass{1}, Exp.fitMean(1.0)); % mean = 1
node{2}.setService(jobclass{1}, Exp.fitMean(1.5)); % mean = 1.5
node{2}.setLoadDependence(min(1:N,c)); % multi-server with c servers

ldmodel.link(ldmodel.serialRouting(node));

lldAvgTableCTMC=CTMC(ldmodel).getAvgTable

lldAvgTableNC=NC(ldmodel).getAvgTable
lldAvgTableRD=NC(ldmodel,'method','rd').getAvgTable
lldAvgTableNRP=NC(ldmodel,'method','nrp').getAvgTable
lldAvgTableNRL=NC(ldmodel,'method','nrl').getAvgTable

lldAvgTableMVALD=MVA(ldmodel,'method','exact').getAvgTable
lldAvgTableQD=MVA(ldmodel,'method','qd').getAvgTable

lldAvgTableJMT=JMT(ldmodel,'seed',23000).getAvgTable

%% casted with scaling function that depends on the per-class queue population
cdmodel = Network('model');
node{1} = Delay(cdmodel, 'Delay');
node{2} = Queue(cdmodel, 'Queue1', SchedStrategy.FCFS);
jobclass{1} = ClosedClass(cdmodel, 'Class1', N, node{1}, 0);
node{1}.setService(jobclass{1}, Exp.fitMean(1.0)); % mean = 1
node{2}.setService(jobclass{1}, Exp.fitMean(1.5)); % mean = 1.5
node{2}.setClassDependence(@(ni) min(sum(ni),c)); % ni is a vector where ni(r) is the number of jobs in class r at station i

cdmodel.link(cdmodel.serialRouting(node));

cdAvgTableCTMC=CTMC(cdmodel).getAvgTable
cdAvgTableCD=MVA(cdmodel,'method','qd').getAvgTable
cdAvgTableJMT=JMT(cdmodel,'seed',23000).getAvgTable
