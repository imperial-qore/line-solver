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
    seed = 23000;
end
rng(seed);
for i=1:M
    station{i} = Queue(model, ['Queue',num2str(i)], SchedStrategy.PS); 
end
if useDelay
    station{M+1} = Delay(model, 'Delay1');
end
%% Block 2: classes
jobclass = ClosedClass(model, 'Class1', round(rand*M+2), station{1}, 0);

for i=1:M
    station{i}.setService(jobclass, Exp.fitMean(rand()+i)); % (Queue 1,Class1)
end

if useDelay
    station{M+1}.setService(jobclass, Exp.fitMean(2.0)); % (Delay 1,Class1)
end

%% Block 3: topology
P = model.initRoutingMatrix(); % initialize routing matrix 
P{jobclass,jobclass} = Network.serialRouting(station);
model.link(P);
end