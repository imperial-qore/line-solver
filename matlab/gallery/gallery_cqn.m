function model = gallery_cqn(M, useDelay, seed)
model = Network('Single-class CQN');
%% Block 1: nodes

if nargin==0
    M=2;
end
if nargin<2
    useDelay = false;
end
if nargin<3
    seed = 2300;
end
rng(seed);
for i=1:M
    node{i} = Queue(model, ['Queue ',num2str(i)], SchedStrategy.PS); 
end
if useDelay
    node{M+1} = DelayStation(model, 'Delay 1');
end
%% Block 2: classes
jobclass{1} = ClosedClass(model, 'Class1', round(rand*10*M+3), node{1}, 0);

for i=1:M
    node{i}.setService(jobclass{1}, Exp.fitMean(rand()+i)); % (Queue 1,Class1)
end
if useDelay
    node{M+1}.setService(jobclass{1}, Exp.fitMean(2.000000)); % (Delay 1,Class1)
end

%% Block 3: topology
P = model.initRoutingMatrix(); % initialize routing matrix 
P{jobclass{1},jobclass{1}} = Network.serialRouting(node);
model.link(P);
end