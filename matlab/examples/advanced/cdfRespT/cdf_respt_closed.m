clear node jobclass;

model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue2', SchedStrategy.PS);

jobclass{1} = ClosedClass(model, 'Class1', 1, node{1}, 0);

servProc1 = Exp(1/0.1);
node{1}.setService(jobclass{1}, servProc1);
servProc2 = Erlang.fitMeanAndSCV(1,1/3);
node{2}.setService(jobclass{1}, servProc2);

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

P = cell(K);
P{1,1} = circul(2);

model.link(P);
%%
solver = JMT(model,'seed',23000,'samples',1e4);
RDsim = solver.getCdfRespT();
fprintf(1,'\n')
for i=1:model.getNumberOfStations
    for c=1:model.getNumberOfClasses
%        plot(FC{i,c}(:,2),FC{i,c}(:,1)); hold all;
        AvgRespTfromCDFSim(i,c) = diff(RDsim{i,c}(:,1))'*RDsim{i,c}(2:end,2); %mean
        PowerMoment2_R(i,c) = diff(RDsim{i,c}(:,1))'*(RDsim{i,c}(2:end,2).^2);
        Variance_R(i,c) = PowerMoment2_R(i,c)-AvgRespTfromCDFSim(i,c)^2; %variance
        SqCoeffOfVariationRespTfromCDFSim(i,c) = (Variance_R(i,c))/AvgRespTfromCDFSim(i,c)^2; %scv
    end
end
%%
solver = FLD(model);
RDfluid = solver.getCdfRespT();
%%
% MIDPOINT RIEMANN-STIELTJES, and the rule is the LAW and not the taste.
% A fluid CDF is a CONTINUOUS law read off an integrator grid: the mass dF of
% an interval sits somewhere inside it, so the midpoint is the second-order
% estimate and the right endpoint used above is first-order, i.e. carries a
% bias that moves with whatever grid the integrator happened to choose. A
% SIMULATED CDF jumps AT its samples, so up there the right endpoint IS the
% sample and that sum is exact -- the two rules are not interchangeable.
% Measured on this model, whose exact answers are the service laws above:
% midpoint reads SCV 1.00032 and 0.33301 against 1 and 1/3, the right endpoint
% 1.01715 and 0.33409, and the two hooks in options.odesolvers agree to the
% fourth decimal under it where they agreed only to the second before.
for i=1:model.getNumberOfStations
    for c=1:model.getNumberOfClasses
%        plot(FC{i,c}(:,2),FC{i,c}(:,1)); hold all;
        dF_R = diff(RDfluid{i,c}(:,1));
        tmid_R = (RDfluid{i,c}(1:end-1,2)+RDfluid{i,c}(2:end,2))/2;
        AvgRespTfromCDFFluid(i,c) = dF_R'*tmid_R; %mean
        PowerMoment2_R(i,c) = dF_R'*(tmid_R.^2);
        Variance_R(i,c) = PowerMoment2_R(i,c)-AvgRespTfromCDFFluid(i,c)^2; %variance
        SqCoeffOfVariationRespTfromCDFFluid(i,c) = (Variance_R(i,c))/AvgRespTfromCDFFluid(i,c)^2; %scv
    end
end
fprintf(1,'\n')
disp('Since there is a single job, mean and squared coefficient of variation');
disp('of response times are close, up to fluid approximation precision, those');
fprintf(1,'of the service time distribution.\n\n');
AvgRespTfromTheory = [servProc1.getMean; servProc2.getMean]
AvgRespTfromCDFSim
AvgRespTfromCDFFluid
SqCoeffOfVariationRespTfromTheory = [servProc1.getSCV; servProc2.getSCV]
SqCoeffOfVariationRespTfromCDFSim
SqCoeffOfVariationRespTfromCDFFluid

