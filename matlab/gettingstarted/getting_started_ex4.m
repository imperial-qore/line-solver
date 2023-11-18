% https://github.com/imperial-qore/line-solver/wiki/Getting-started#example-4-round-robin-load-balancing
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

model.addLinks([source, lb; 
                lb,     queue1; 
                lb,     queue2; 
                queue1, sink; 
                queue2, sink]);
            
lb.setRouting(oclass, RoutingStrategy.RAND);
jmtAvgTable = SolverJMT(model,'seed',23000).getAvgTable

lb.setRouting(oclass, RoutingStrategy.RROBIN);
model.reset();
jmtAvgTableRR = SolverJMT(model,'seed',23000).getAvgTable
