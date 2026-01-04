clear solver AvgTable;
model = Network('model');

source = Source(model,'Source');
queue1 = Queue(model,'Queue1',SchedStrategy.FCFS);
queue2 = Queue(model,'Queue2',SchedStrategy.FCFS);
fork1 = Fork(model,'Fork1');
join1 = Join(model,'Join1',fork1);
queue3 = Queue(model,'Queue3',SchedStrategy.FCFS);
queue4 = Queue(model,'Queue4',SchedStrategy.FCFS);
fork2 = Fork(model,'Fork2');
join2 = Join(model,'Join2',fork2);
sink = Sink(model,'Sink');

jobclass1 = OpenClass(model, 'class1');

source.setArrival(jobclass1, Exp(0.4));
queue1.setService(jobclass1, Exp(1.0));
queue2.setService(jobclass1, Exp(1.0));
queue3.setService(jobclass1, Exp(1.0));
queue4.setService(jobclass1, Exp(1.0));

P = cellzeros(2,2,6,6);
P{jobclass1,jobclass1}(source,fork1) = 1;
P{jobclass1,jobclass1}(fork1,queue1) = 1.0;
P{jobclass1,jobclass1}(fork1,queue2) = 1.0;
P{jobclass1,jobclass1}(queue1,join1) = 1.0;
P{jobclass1,jobclass1}(queue2,join1) = 1.0;
P{jobclass1,jobclass1}(join1,fork2) = 1.0;
P{jobclass1,jobclass1}(fork2,queue3) = 1.0;
P{jobclass1,jobclass1}(fork2,queue4) = 1.0;
P{jobclass1,jobclass1}(queue3,join2) = 1.0;
P{jobclass1,jobclass1}(queue4,join2) = 1.0;
P{jobclass1,jobclass1}(join2,sink) = 1.0;

model.link(P);

solver{1} = JMT(model,'seed',23000);
solver{end+1} = MVA(model);
%solver{end+1} = DES(model,'seed',23000);

for s=1:length(solver)
    AvgTable{s} = solver{s}.getAvgTable;
    AvgTable{s}
end
