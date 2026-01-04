clear node jobclass;

model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue2', SchedStrategy.PS);

jobclass{1} = ClosedClass(model, 'Class1', 1, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 0, node{1}, 0);
jobclass{3} = ClosedClass(model, 'Class3', 0, node{1}, 0);

jobclass{1}.completes = false;

node{1}.setService(jobclass{1}, Exp(1/1));
node{1}.setService(jobclass{2}, Exp(1/1));
node{1}.setService(jobclass{3}, Exp(1/1));

node{2}.setService(jobclass{1}, Exp(1/1));
node{2}.setService(jobclass{2}, Erlang(1/2,2));
node{2}.setService(jobclass{3}, Exp(1/0.01));

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

P = cell(K,K);

P{1,1} = [0,1; 0,0];
P{1,2} = [0,0; 1,0];
P{1,3} = zeros(M);

P{2,1} = [0,0; 1,0];
P{2,2} = [0,1; 0,0];
P{2,3} = zeros(M);

P{3,1} = zeros(M);
P{3,2} = zeros(M);
P{3,3} = circul(M);

model.link(P);
%%
options = FLD.defaultOptions;
options.method = 'statedep';
options.iter_max = 100;
solver = FLD(model, options);
AvgRespT = solver.getAvgRespT
solver = FLD(model, options);
RDfluid = solver.getCdfRespT();

%%
for i=1:model.getNumberOfStations
    for c=1:model.getNumberOfClasses
%        plot(FC{i,c}(:,2),FC{i,c}(:,1)); hold all;
        AvgRespTfromCDF(i,c) = diff(RDfluid{i,c}(:,1))'*RDfluid{i,c}(2:end,2); %mean
        PowerMoment2_R(i,c) = diff(RDfluid{i,c}(:,1))'*(RDfluid{i,c}(2:end,2).^2);
        Variance_R(i,c) = PowerMoment2_R(i,c)-AvgRespTfromCDF(i,c)^2; %variance
        SqCoeffOfVariationRespTfromCDF(i,c) = (Variance_R(i,c))/AvgRespTfromCDF(i,c)^2; %scv
    end
end
AvgRespTfromCDF
SqCoeffOfVariationRespTfromCDF