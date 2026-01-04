clear node jobclass

% class-dependent model
N = 16; % number of jobs
c = 2;
%%
cdmodel = Network('model');
node{1} = Delay(cdmodel, 'Delay');
node{2} = Queue(cdmodel, 'Queue1', SchedStrategy.PS);
jobclass{1} = ClosedClass(cdmodel, 'Class1', N, node{1}, 0);
jobclass{2} = ClosedClass(cdmodel, 'Class2', N/2, node{1}, 0);
node{1}.setService(jobclass{1}, Exp.fitMean(1.0)); % mean = 1
node{1}.setService(jobclass{2}, Exp.fitMean(2.0)); % mean = 1
node{2}.setService(jobclass{1}, Exp.fitMean(1.5)); % mean = 1.5
node{2}.setService(jobclass{2}, Exp.fitMean(2.5)); % mean = 1.5
node{2}.setClassDependence(@(ni) min(ni(1),c)); % multi-server only for class-1 jobs

P = cdmodel.initRoutingMatrix();
P{1,1} = cdmodel.serialRouting(node);
P{2,2} = cdmodel.serialRouting(node);
cdmodel.link(P);

cdAvgTableCTMC=CTMC(cdmodel).getAvgTable
cdAvgTableCD=MVA(cdmodel,'method','qd').getAvgTable
cdAvgTableJMT=JMT(cdmodel,'seed',23000).getAvgTable

model = cdmodel; % for test compatibility
