clear node jobclass

% class-dependent (product-form) model.
% Each class sees a service-rate scaling beta_{i,r} that depends ONLY on its
% own per-class population n_{i,r} at the station. This is the BCMP
% product-form case of QD-AMVA (Casale, Perez, Wang, IFIP PERFORMANCE 2015):
% D_{i,r}(n) = theta_{i,r} * beta_{i,r}(n_{i,r}). The handle returns a length-R
% vector [beta_{i,1}(n), beta_{i,2}(n)], of which the solver picks class r.
% Contrast ld_joint_dependence.m, where the scalar min(ni(1),c) reads a foreign
% class marginal and is therefore non-product-form (setJointDependence).
N = 16; % number of jobs
c = 2;
%%
cdmodel = Network('model');
node{1} = Delay(cdmodel, 'Delay');
node{2} = Queue(cdmodel, 'Queue1', SchedStrategy.PS);
jobclass{1} = ClosedClass(cdmodel, 'Class1', N, node{1}, 0);
jobclass{2} = ClosedClass(cdmodel, 'Class2', N/2, node{1}, 0);
node{1}.setService(jobclass{1}, Exp.fitMean(1.0)); % mean = 1
node{1}.setService(jobclass{2}, Exp.fitMean(2.0)); % mean = 2
node{2}.setService(jobclass{1}, Exp.fitMean(1.5)); % mean = 1.5
node{2}.setService(jobclass{2}, Exp.fitMean(2.5)); % mean = 2.5
% beta_{i,r}(n_{i,r}): class 1 scales up to c servers with its OWN count, class
% 2 is always single-server. Both entries read only their own marginal, so the
% demands satisfy the product-form recurrence. Peak rate scaling [c 1] per
% class normalizes Util = T*S/peak.
node{2}.setClassDependence(@(ni) [min(ni(1),c), 1], [c 1]);

P = cdmodel.initRoutingMatrix();
P{1,1} = cdmodel.serialRouting(node);
P{2,2} = cdmodel.serialRouting(node);
cdmodel.link(P);

cdAvgTableCTMC=CTMC(cdmodel,'exact').getAvgTable
cdAvgTableCD=MVA(cdmodel,'method','qd').getAvgTable
% JMT is not solved here: the JSIM writer has no representation for the
% class-dependence handle, so SolverJMT rejects the model rather than
% silently solving it unscaled (see SolverJMT.getFeatureSet).

model = cdmodel; % for test compatibility
