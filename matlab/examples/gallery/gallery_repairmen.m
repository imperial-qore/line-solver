function model = gallery_repairmen(nservers, seed)
model = Network('Finite repairmen CQN');
%% Block 1: nodes

M=1;

if nargin<=1
    nservers=1;
end

if nargin<3
    seed = 2300;
end
rng(seed);
node{1} = Queue(model, ['Queue',num2str(1)], SchedStrategy.PS);
node{1}.setNumberOfServers(nservers);
node{M+1} = Delay(model, 'Delay1');
%% Block 2: classes
jobclass = ClosedClass(model, 'Class1', round(rand*10*M+3), node{1}, 0);

node{1}.setService(jobclass, Exp.fitMean(rand()+1)); % (Queue 1,Class1)
node{M+1}.setService(jobclass, Exp.fitMean(2.000000)); % (Delay 1,Class1)

%% Block 3: topology
P = model.initRoutingMatrix(); % initialize routing matrix
P{jobclass,jobclass} = Network.serialRouting(node);
model.link(P);
end