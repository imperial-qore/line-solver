function model = gallery_mm1_ps_feedback(p)
if ~exist('p','var')
    p = 1/3;
end
model = Network('M/M/1-PS-Feedback');
%% Block 1: nodes
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.PS);
sink = Sink(model, 'Sink');
%% Block 2: classes
oclass1 = OpenClass(model, 'Class1');
source.setArrival(oclass1, Exp.fitMean(1));
queue.setService(oclass1, Exp.fitMean(0.5));
%% Block 3: topology
P = model.initRoutingMatrix;
P{oclass1,oclass1}(source,queue)=1;
P{oclass1,oclass1}(queue,queue)=p;
P{oclass1,oclass1}(queue,sink)=1-p;
model.link(P);
end