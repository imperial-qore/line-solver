model = Network('model');

delay = Delay(model, 'Delay1');
queue1 = Queue(model,'Queue1',SchedStrategy.PS);
queue2 = Queue(model,'Queue2',SchedStrategy.PS);
queue3 = Queue(model,'Queue3',SchedStrategy.PS);
fork = Fork(model,'Fork');
join = Join(model,'Join', fork);

jobclass1 = ClosedClass(model, 'class1', 10, delay, 0);
jobclass2 = ClosedClass(model, 'class2', 10, delay, 0);

queue1.setService(jobclass1, Exp(1.5));
queue2.setService(jobclass1, Exp(1.1));
queue3.setService(jobclass1, Exp(2.5));
delay.setService(jobclass1, Exp(0.5));

queue1.setService(jobclass2, Exp(2.8));
queue2.setService(jobclass2, Exp(3));
queue3.setService(jobclass2, Exp(1.0));
delay.setService(jobclass2, Exp(0.8));

P = model.initRoutingMatrix;
P{jobclass1, jobclass1}(delay, fork) = 1.0;
P{jobclass1, jobclass1}(fork,queue1) = 1.0;
P{jobclass1, jobclass1}(fork,queue2) = 1.0;
P{jobclass1, jobclass1}(queue2, queue3) = 1.0;
P{jobclass1, jobclass1}(queue3, join) = 1.0;
P{jobclass1, jobclass1}(queue1,join) = 1.0;
P{jobclass1, jobclass1}(join,delay) = 1.0;

P{jobclass2, jobclass2}(delay, fork) = 1.0;
P{jobclass2, jobclass2}(fork,queue1) = 1.0;
P{jobclass2, jobclass2}(fork,queue2) = 1.0;
P{jobclass2, jobclass2}(queue2, queue3) = 1.0;
P{jobclass2, jobclass2}(queue3, join) = 1.0;
P{jobclass2, jobclass2}(queue1,join) = 1.0;
P{jobclass2, jobclass2}(join,delay) = 1.0;

model.link(P);

MVA(model).getAvgTable