from line_solver import *

if __name__ == "__main__":
    GlobalConstants.setVerbose(VerboseLevel.STD)

    model = Network("M/M/1 model")
    source = Source(model, "Source")
    queue = Queue(model, "Queue", SchedStrategy.FCFS)
    sink = Sink(model, "Sink")

    jobclass = OpenClass(model, "Class1")
    source.setArrival(jobclass, Exp(1.0))
    queue.setService(jobclass, Exp(2.0))

    model.addLink(source, queue)
    model.addLink(queue, sink)

    solver = JMT(model,verbose=True)
    table = solver.avg_table()





