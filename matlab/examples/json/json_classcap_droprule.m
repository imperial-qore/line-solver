function model = json_classcap_droprule()
% Open 2-class network with per-class capacity, drop rules, and load-dependence.
model = Network('ClassCap_DropRule');

source = Source(model, 'Source');
queue  = Queue(model, 'Queue', SchedStrategy.PS);
sink   = Sink(model, 'Sink');

queue.setNumberOfServers(2);
queue.setCapacity(20);

class1 = OpenClass(model, 'Class1');
class2 = OpenClass(model, 'Class2');

source.setArrival(class1, Exp(1));
source.setArrival(class2, Exp(0.5));
queue.setService(class1, Exp(3));
queue.setService(class2, Exp(2));

% Per-class capacity
queue.classCap(class1.index) = 8;
queue.classCap(class2.index) = 15;

% Per-class drop rules. WAITQ is not usable here: no solver honours "wait
% upstream" for an open class at a plain finite capacity (CTMC drops the
% arrival, JMT blocks at the source and ignores the cap), so refreshCapacity
% rejects that combination. BAS is the honoured blocking policy.
queue.setDropRule(class1, DropStrategy.DROP);
queue.setDropRule(class2, DropStrategy.BAS);

% Load-dependent scaling
queue.setLoadDependence([1.0, 0.9, 0.8, 0.7]);

P = model.initRoutingMatrix();
P{class1} = Network.serialRouting(source, queue, sink);
P{class2} = Network.serialRouting(source, queue, sink);
model.link(P);
end
