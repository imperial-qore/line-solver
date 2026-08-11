clear node jobclass;

model = Network('myModel');

%% Block 1: nodes
node{1} = Source(model, 'Source');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
node{3} = Queue(model, 'Queue2', SchedStrategy.FCFS);
node{4} = Sink(model, 'Sink');

%% Block 2: classes
jobclass{1} = OpenClass(model, 'Class1', 0);
jobclass{2} = OpenClass(model, 'Class2', 0);

node{1}.setArrival(jobclass{1}, Exp.fitMean(4.000000)); % (Source,Class1)
node{1}.setArrival(jobclass{2}, Exp.fitMean(4.000000)); % (Source,Class2)

node{2}.setService(jobclass{1}, Exp.fitMean(1.000000)); % (Queue1,Class1)
node{2}.setService(jobclass{2}, Exp.fitMean(1.000000)); % (Queue1,Class2)

node{3}.setService(jobclass{1}, Exp.fitMean(1.000000)); % (Queue2,Class1)
node{3}.setService(jobclass{2}, Exp.fitMean(1.000000)); % (Queue2,Class2)

%% Block 3: topology
P = model.initRoutingMatrix(); % initialize routing matrix 
P{1,1}(1,2) = 1;
P{1,2}(2,3) = 1;
P{2,1}(3,4) = 1;
P{2,2}(1,2) = 1;
P{2,1}(2,3) = 1;
P{1,2}(3,4) = 1;
model.link(P);

%%
options = FLD.defaultOptions;
options.iter_max = 300;
RDfluid = FLD(model,options).getCdfRespT()

%%
jmtoptions = JMT.defaultOptions;
jmtoptions.samples = 1e4;
jmtoptions.seed = 23000;
RDsim = JMT(model, jmtoptions).getTranCdfRespT();

%%
figure;
for i=2:model.getNumberOfStations
    subplot(model.getNumberOfStations-1,2,2*(i-2)+1)
    loglog(RDsim{i,1}(:,2),1-RDsim{i,1}(:,1),'r--')
    hold on;
    loglog(RDfluid{i,1}(:,2),1-RDfluid{i,1}(:,1),'b-.')
    legend('sim','fluid','Location','SouthWest');
    title(['RespT Tail: Node ',num2str(i),', Class ',num2str(1)]);

    subplot(model.getNumberOfStations-1,2,2*(i-2)+2)
    loglog(RDsim{i,2}(:,2),1-RDsim{i,2}(:,1),'r--')
    hold on;
    loglog(RDfluid{i,2}(:,2),1-RDfluid{i,2}(:,1),'b-.')
    legend('sim','fluid','Location','SouthWest');
    title(['RespT Tail: Node ',num2str(i),', Class ',num2str(2)]);
end

%%
for i=2:model.getNumberOfStations
    for c=1:model.getNumberOfClasses
        AvgRespTfromCDFfluid(i,c) = diff(RDfluid{i,c}(:,1))'*RDfluid{i,c}(2:end,2);
        AvgRespTfromCDFsim(i,c) = diff(RDsim{i,c}(:,1))'*RDsim{i,c}(2:end,2);
    end
end

AvgRespTfromCDFfluid
AvgRespTfromCDFsim