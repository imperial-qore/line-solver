import unittest
from line_solver import *

class test_getting_started_ex1(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(test_getting_started_ex1, self).__init__(*args, **kwargs)

        model = Network("M/M/1 model")
        source = Source(model, "mySource")
        queue = Queue(model, "myQueue", SchedStrategy.FCFS)
        sink = Sink(model, "mySink")

        # An M/M/1 queue with arrival rate 0.5 and service rate 1.0
        openclass = OpenClass(model, "Class1")
        source.setArrival(openclass, Exp(0.5))
        queue.setService(openclass, Exp(1.0))

        model.addLink(source, queue)
        model.addLink(queue, sink)

        solver = SolverJMT(model, 'seed', 23000)
        self.table = solver.getAvgTable()  # pandas.DataFrame
        model.getStruct()

    def test_table(self):
        print(self.table)
        self.assertAlmostEqual(self.table.QLen[0], 0.0, 8)  # add assertion here
        self.assertAlmostEqual(self.table.QLen[1], 0.9441511643947065, 8)  # add assertion here
        self.assertAlmostEqual(self.table.Util[0], 0.0, 8)  # add assertion here
        self.assertAlmostEqual(self.table.Util[1], 0.5033610765632501, 8)  # add assertion here
        self.assertAlmostEqual(self.table.RespT[0], 0.0, 8)  # add assertion here
        self.assertAlmostEqual(self.table.RespT[1], 1.9701934041278817, 8)  # add assertion here
        self.assertAlmostEqual(self.table.ResidT[0], 0.0, 8)  # add assertion here
        self.assertAlmostEqual(self.table.ResidT[1], 1.9701934041278817, 8)  # add assertion here
        self.assertAlmostEqual(self.table.Tput[0], 0.5047874245996684, 8)  # add assertion here
        self.assertAlmostEqual(self.table.Tput[1], 0.5007972686846811, 8)  # add assertion here

if __name__ == '__main__':
    unittest.main()
