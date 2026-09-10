function model = json_signal_classes()
% Open network with G-network negative signal.
model = Network('SignalClasses');

source = Source(model, 'Source');
queue  = Queue(model, 'Queue', SchedStrategy.FCFS);
sink   = Sink(model, 'Sink');

class1  = OpenClass(model, 'Class1');
signal1 = OpenSignal(model, 'Signal1', SignalType.NEGATIVE);
signal1 = signal1.forJobClass(class1);
signal1.setRemovalPolicy(RemovalPolicy.RANDOM);

source.setArrival(class1, Exp(2));
source.setArrival(signal1, Exp(0.5));
queue.setService(class1, Exp(5));
queue.setService(signal1, Immediate());

P = model.initRoutingMatrix();
P{class1}(source, queue) = 1;
P{class1}(queue, sink) = 1;
P{signal1}(source, queue) = 1;
P{signal1}(queue, sink) = 1;
model.link(P);
end
