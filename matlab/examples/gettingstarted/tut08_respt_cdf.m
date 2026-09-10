% response-time-distribution-and-percentiles
model = Network('Model');

% Block 1: nodes
node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);

% Block 2: classes
jobclass{1} = ClosedClass(model, 'Class1', 5, node{1}, 0);
node{1}.setService(jobclass{1}, Exp(1.0));
node{2}.setService(jobclass{1}, Exp(0.5));

% Block 3: topology
model.link(Network.serialRouting(node{1},node{2}));
%%
% Block 4: solution
RDfluid = FLD(model,'seed',23000).getCdfRespT();
RDsim = JMT(model,'seed',23000,'samples',1e4).getCdfRespT();

% Plot results
if ~isempty(RDsim{2,1})
    semilogx(RDsim{2,1}(:,2),1-RDsim{2,1}(:,1),'r'); hold on;
    semilogx(RDfluid{2,1}(:,2),1-RDfluid{2,1}(:,1),'k--');
    legend('jmt-transient','fluid-steady','Location','Best');
    ylabel('Pr(T > t)'); xlabel('time t');
end
%%
% Compute CDF-derived scalar statistics
M = model.getNumberOfStations;
K = model.getNumberOfClasses;

AvgRespTfromCDFSim = zeros(M,K);
SqCoeffOfVariationRespTfromCDFSim = zeros(M,K);
for i=1:M
    for c=1:K
        if ~isempty(RDsim{i,c}) && size(RDsim{i,c},2) >= 2
            AvgRespTfromCDFSim(i,c) = diff(RDsim{i,c}(:,1))'*RDsim{i,c}(2:end,2);
            PowerMoment2_R(i,c) = diff(RDsim{i,c}(:,1))'*(RDsim{i,c}(2:end,2).^2);
            Variance_R(i,c) = PowerMoment2_R(i,c)-AvgRespTfromCDFSim(i,c)^2;
            SqCoeffOfVariationRespTfromCDFSim(i,c) = Variance_R(i,c)/AvgRespTfromCDFSim(i,c)^2;
        end
    end
end

AvgRespTfromCDFFluid = zeros(M,K);
SqCoeffOfVariationRespTfromCDFFluid = zeros(M,K);
for i=1:M
    for c=1:K
        if ~isempty(RDfluid{i,c}) && size(RDfluid{i,c},2) >= 2
            AvgRespTfromCDFFluid(i,c) = diff(RDfluid{i,c}(:,1))'*RDfluid{i,c}(2:end,2);
            PowerMoment2_R(i,c) = diff(RDfluid{i,c}(:,1))'*(RDfluid{i,c}(2:end,2).^2);
            Variance_R(i,c) = PowerMoment2_R(i,c)-AvgRespTfromCDFFluid(i,c)^2;
            SqCoeffOfVariationRespTfromCDFFluid(i,c) = Variance_R(i,c)/AvgRespTfromCDFFluid(i,c)^2;
        end
    end
end

AvgRespTfromCDFSim
AvgRespTfromCDFFluid
SqCoeffOfVariationRespTfromCDFSim
SqCoeffOfVariationRespTfromCDFFluid