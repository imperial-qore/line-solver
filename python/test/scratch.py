from line_solver import *

if __name__ == "__main__":
    model = Network('LoadBalCQN')
    # Block 1: nodes
    delay = Delay(model, 'Think')
    queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.PS)
    # Block 2: classes
    cclass = ClosedClass(model, 'Job1', 16, delay)
    delay.setService(cclass, Exp(1))
    queue1.setService(cclass, Exp(0.75))
    queue2.setService(cclass, Exp(0.50))
    P = model.initRoutingMatrix()
    P.set(cclass, cclass, queue1, delay, 1.0)
    P.set(cclass, cclass, queue2, delay, 1.0)
    p = 0.60
    P.set(cclass, cclass, delay, queue1, p)
    P.set(cclass, cclass, delay, queue2, 1.0 - p)
    model.link(P)
    # solver = SolverMVA(model)
    # table = solver.getAvgTable()  # pandas.DataFrame
    # print(table)

    sn = model.getStruct()

    print(sn['inchain']) # this seems buggy

