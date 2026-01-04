% Example 4: Round robin load balancing
model = Network('RRLB');

source = Source(model, 'Source');
lb = Router(model, 'LB');
queue1 = Queue(model, 'Queue1', SchedStrategy.PS);
queue2 = Queue(model, 'Queue2', SchedStrategy.PS);
sink  = Sink(model, 'Sink');

oclass = OpenClass(model, 'Class1');
source.setArrival(oclass, Exp(1));
queue1.setService(oclass, Exp(2));
queue2.setService(oclass, Exp(2));

% Add links individually since array concatenation is not supported
model.addLink(source, lb);
model.addLink(lb, queue1);
model.addLink(lb, queue2);
model.addLink(queue1, sink);
model.addLink(queue2, sink);

lb.setRouting(oclass, RoutingStrategy.RAND);
jmtAvgTable = JMT(model,'seed',23000).avgTable()

lb.setRouting(oclass, RoutingStrategy.RROBIN);
model.reset();
jmtAvgTableRR = JMT(model,'seed',23000).avgTable()
