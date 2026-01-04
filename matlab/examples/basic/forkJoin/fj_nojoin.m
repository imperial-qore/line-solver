clear solver AvgTable;
model = Network('model');

source = Source(model,'Source');
queue1 = Queue(model,'Queue1',SchedStrategy.PS);
queue2 = Queue(model,'Queue2',SchedStrategy.PS);
queue3 = Queue(model,'Queue3',SchedStrategy.PS);
fork = Fork(model,'Fork');
sink = Sink(model,'Sink');

jobclass1 = OpenClass(model, 'class1');

source.setArrival(jobclass1, Exp(0.5));
queue1.setService(jobclass1, Exp(1.0));
queue2.setService(jobclass1, Exp(2.0));
queue3.setService(jobclass1, Exp(3.0));

P = zeros(5);
P(source,fork) = 1;
P(fork,queue1) = 1.0;
P(fork,queue2) = 1.0;
P(fork,queue3) = 1.0;
P(queue1,sink) = 1.0;
P(queue2,sink) = 1.0;
P(queue3,sink) = 1.0;

model.link(P);

solver = {};
solver{end+1} = JMT(model,'seed',23000);
solver{end+1} = MVA(model);
%solver{end+1} = DES(model,'seed',23000);

AvgTable = {};
for s=1:length(solver)
    AvgTable{end+1} = solver{s}.getAvgTable;
    AvgTable{s}
end
