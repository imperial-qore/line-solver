function model = json_fcr_details()
% Open 2-class network with a Finite Capacity Region.
model = Network('FCR_Details');

source = Source(model, 'Source');
queue1 = Queue(model, 'Queue1', SchedStrategy.PS);
queue2 = Queue(model, 'Queue2', SchedStrategy.PS);
sink   = Sink(model, 'Sink');

class1 = OpenClass(model, 'Class1');
class2 = OpenClass(model, 'Class2');

source.setArrival(class1, Exp(1));
source.setArrival(class2, Exp(0.5));
queue1.setService(class1, Exp(3));
queue1.setService(class2, Exp(2));
queue2.setService(class1, Exp(4));
queue2.setService(class2, Exp(3));

P = model.initRoutingMatrix();
P{class1} = Network.serialRouting(source, queue1, queue2, sink);
P{class2} = Network.serialRouting(source, queue1, queue2, sink);
model.link(P);

% Finite Capacity Region
fcr = Region({queue1, queue2}, {class1, class2});
fcr.setGlobalMaxJobs(15);
fcr.setClassMaxJobs(class1, 8);
fcr.setClassMaxJobs(class2, 10);
fcr.setClassWeight(class1, 1.0);
fcr.setClassWeight(class2, 2.0);
fcr.setClassSize(class1, 1);
fcr.setClassSize(class2, 3);
end
