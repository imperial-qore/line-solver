from line_solver import *

if __name__ == "__main__":
    GlobalConstants.setVerbose(VerboseLevel.STD)

    model = Network("M/M/1 model")
    source = Source(model, "Source")
    queue = Queue(model, "Queue", SchedStrategy.FCFS)
    sink = Sink(model, "Sink")

    # An M/M/1 queue with arrival rate 1.0 and service rate 2.0
    jobclass = OpenClass(model, "Class1")
    source.setArrival(jobclass, Exp(1.0))
    queue.setService(jobclass, Exp(2.0))

    model.addLink(source, queue)
    model.addLink(queue, sink)

    solver = SolverJMT(model)
    table = solver.getAvgTable()  # pandas.DataFrame

    # uncomment the next lines to filter specific results
    #print(tget(table,"Queue"))
    #print(tget(table,"RespT","Queue"))
    #print(tget(table,"RespT","Queue","Class1"))



