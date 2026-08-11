clear node jobclass

% joint-dependent (non-product-form) model
N = 16; % number of jobs
c = 2;
%%
jdmodel = Network('model');
node{1} = Delay(jdmodel, 'Delay');
node{2} = Queue(jdmodel, 'Queue1', SchedStrategy.PS);
jobclass{1} = ClosedClass(jdmodel, 'Class1', N, node{1}, 0);
jobclass{2} = ClosedClass(jdmodel, 'Class2', N/2, node{1}, 0);
node{1}.setService(jobclass{1}, Exp.fitMean(1.0)); % mean = 1
node{1}.setService(jobclass{2}, Exp.fitMean(2.0)); % mean = 1
node{2}.setService(jobclass{1}, Exp.fitMean(1.5)); % mean = 1.5
node{2}.setService(jobclass{2}, Exp.fitMean(2.5)); % mean = 1.5
node{2}.setJointDependence(@(ni) min(ni(1),c), c); % NON-product-form: rate reads the class-1 marginal only (eta_i), shared across classes; peak rate scaling = c (Util = T*S/c)

P = jdmodel.initRoutingMatrix();
P{1,1} = jdmodel.serialRouting(node);
P{2,2} = jdmodel.serialRouting(node);
jdmodel.link(P);

jdAvgTableCTMC=CTMC(jdmodel).getAvgTable
jdAvgTableJD=MVA(jdmodel,'method','qd').getAvgTable
% JMT is not solved here: the JSIM writer has no representation for the
% joint-dependence handle, so SolverJMT rejects the model rather than
% silently solving it unscaled (see SolverJMT.getFeatureSet).

model = jdmodel; % for test compatibility
