clear;
model = Network('myModel');
source = Source(model, 'Source');
queue = Queue(model, 'Queue1', SchedStrategy.FCFS);
delay = Delay(model, 'ThinkTime');
sink = Sink(model, 'Sink');
class1 = OpenClass(model, 'Class1');
class2 = ClosedClass(model, 'Class2', 5, queue);
class3 = ClosedClass(model, 'Class3', 0, queue, 0);

source.setArrival(class1, Exp.fitMean(10));
delay.setService(class1, Exp.fitMean(1));
delay.setService(class2, Exp.fitMean(2));
delay.setService(class3, Exp.fitMean(2));
queue.setService(class1, Exp.fitMean(1));
queue.setService(class2, Exp.fitMean(10));
queue.setService(class3, Exp.fitMean(2));

P = cellzeros(3,3,4);
P{class1,class1}(source,queue) = 1.0;
P{class1,class1}(queue,[queue,delay]) = [0.3,0.7];
P{class1,class1}(delay,sink) = 1.0;
P{class2,class3}(delay,queue) = 1.0; % closed class starts at delay
P{class3,class2}(queue,delay) = 1.0;
model.link(P);

SolverJMT(model).getAvgTable