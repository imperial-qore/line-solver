function model = json_balking_retrial()
% Open 2-class network with balking, retrial, and patience.
model = Network('Balking_Retrial');

source = Source(model, 'Source');
queue  = Queue(model, 'Queue', SchedStrategy.FCFS);
sink   = Sink(model, 'Sink');

queue.setCapacity(15);

class1 = OpenClass(model, 'Class1');
class2 = OpenClass(model, 'Class2');

source.setArrival(class1, Exp(1));
source.setArrival(class2, Exp(0.5));
queue.setService(class1, Exp(2));
queue.setService(class2, Exp(3));

% Class1: balking based on queue length
queue.setBalking(class1, BalkingStrategy.QUEUE_LENGTH, ...
    {{5, 10, 0.3}, {11, Inf, 1.0}});

% Class1: patience (reneging)
queue.setPatience(class1, Exp(0.1));

% Class2: retrial with max attempts
queue.setRetrial(class2, Exp(0.5), 3);

P = model.initRoutingMatrix();
P{class1} = Network.serialRouting(source, queue, sink);
P{class2} = Network.serialRouting(source, queue, sink);
model.link(P);
end
