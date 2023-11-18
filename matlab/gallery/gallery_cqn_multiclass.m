function model = gallery_cqn_multiclass(m,r,wantdelay)
if nargin == 0
    m = 1;
    r = 2;
    wantdelay = true;
end
model = Network('Multi-class CQN');
%% Block 1: nodes
for i=1:m
    node{i} = Queue(model, ['Queue ',num2str(i)], SchedStrategy.PS);
end
if wantdelay
    node{end+1} = DelayStation(model, 'Delay 1');
end
%% Block 2: classes
for s=1:r
    jobclass{s} = ClosedClass(model, ['Class',num2str(s)], 5, node{1}, 0);
end

for s=1:r
    for i=1:m
        node{i}.setService(jobclass{s}, Exp.fitMean(round(50*rand))); % (Queue 1,Class1)
    end
    if wantdelay
        node{end}.setService(jobclass{s}, Exp.fitMean(round(100*rand))); % (Delay 1,Class1)
    end
end

%% Block 3: topology
P = model.initRoutingMatrix(); % initialize routing matrix
for s=1:r
    P{jobclass{s},jobclass{s}} = Network.serialRouting(node);
end
model.link(P);
end