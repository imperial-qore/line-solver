from line_solver import *

if __name__ == "__main__":
    GlobalConstants.setVerbose(VerboseLevel.STD)

    model = Network("M/M/1 model")
    source = Source(model, "mySource")
    queue = Queue(model, "myQueue", SchedStrategy.FCFS)
    sink = Sink(model, "mySink")

    # An M/M/1 queue with arrival rate 0.5 and service rate 1.0
    openclass = OpenClass(model, "Class1")
    source.setArrival(openclass, Exp(1.0))
    queue.setService(openclass, Exp(2.0))

    model.addLink(source, queue)
    model.addLink(queue, sink)

    solver = SolverJMT(model)
    table = solver.getAvgTable()  # pandas.DataFrame




