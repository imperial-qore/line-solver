from line_solver import *
from line_solver import Erlang


# TODO: APH.fitCentral
# def gallery_aphm1():
#     model = Network('APH/M/1')
#     # Block 1: nodes
#     source = Source(model, 'mySource')
#     queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
#     sink = Sink(model, 'mySink')
#     # Block 2: classes
#     oclass = OpenClass(model, 'myClass')
#     source.setArrival(oclass, APH.fitCentral(1, 0.99, 1.999))
#     queue.setService(oclass, Exp(2))
#     # Block 3: topology
#     model.link(Network.serialRouting(source, queue, sink))
#     return model

# TODO: Coxian and Coxian.fitCentral
# def gallery_coxm1():
#     model = Network('APH/M/1')
#     # Block 1: nodes
#     source = Source(model, 'mySource')
#     queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
#     sink = Sink(model, 'mySink')
#     # Block 2: classes
#     oclass = OpenClass(model, 'myClass')
#     source.setArrival(oclass, Coxian.fitCentral(1, 0.99, 1.999))
#     queue.setService(oclass, Exp(2))
#     # Block 3: topology
#     model.link(Network.serialRouting(source, queue, sink))
#     return model

def gallery_detm1():
    model = Network('D/M/1')
    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Det(1))
    queue.setService(oclass, Exp(2))
    # Block 3: topology
    model.link(Network.serialRouting(source, queue, sink))
    return model


def gallery_erlm1():
    model = Network('Er/M/1')
    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Erlang.fitMeanAndOrder(1, 5))
    queue.setService(oclass, Exp(2))
    # Block 3: topology
    model.link(Network.serialRouting(source, queue, sink))
    return model


def gallery_erlm1ps():
    model = Network('Er/M/1-PS')
    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.PS)
    sink = Sink(model, 'mySink')
    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Erlang.fitMeanAndOrder(1, 5))
    queue.setService(oclass, Exp(2))
    # Block 3: topology
    model.link(Network.serialRouting(source, queue, sink))
    return model


def gallery_gamm1():
    model = Network('Gam/M/1')
    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Gamma.fitMeanAndSCV(1, 1 / 5))
    queue.setService(oclass, Exp(2))
    # Block 3: topology
    model.link(Network.serialRouting(source, queue, sink))
    return model


def gallery_hyperlk(k=2):
    model = Network('H/Er/k')
    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, HyperExp.fitMeanAndSCVBalanced(1.0 / 1.8, 4))
    queue.setService(oclass, Erlang.fitMeanAndSCV(1, 0.25))
    queue.setNumberOfServers(k)
    # Block 3: topology
    model.link(Network.serialRouting(source, queue, sink))
    return model


def gallery_hypm1():
    model = Network('H/M/1')
    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, HyperExp.fitMeanAndSCV(1, 64))
    queue.setService(oclass, Exp(2))
    # Block 3: topology
    model.link(Network.serialRouting(source, queue, sink))
    return model


def gallery_mhyp1():
    model = Network('M/H/1')
    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp(1))
    queue.setService(oclass, Coxian.fitMeanAndSCV(0.5, 4))
    # Block 3: topology
    model.link(Network.serialRouting(source, queue, sink))
    return model


def gallery_mm1():
    model = Network('M/M/1')
    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp(1))
    queue.setService(oclass, Exp(2))
    # Block 3: topology
    model.link(Network.serialRouting(source, queue, sink))
    return model


def gallery_mm1_linear(n=2, Umax=0.9):
    model = Network('M/M/1-Linear')

    # Block 1: nodes
    line = [Source(model, 'mySource')]
    for i in range(1, n + 1):
        line.append(Queue(model, 'Queue' + str(i), SchedStrategy.FCFS))
    line.append(Sink(model, 'mySink'))

    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    line[0].setArrival(oclass, Exp(1.0))

    if n == 2:  # linspace has a different behavior in np than matlab in this case
        means = np.linspace(Umax, Umax, 1)
    else:
        means = np.linspace(0.1, Umax, n // 2)

    if n % 2 == 0:
        means = np.concatenate([means, means[::-1]])
    else:
        means = np.concatenate([means, [Umax], means[::-1]])

    for i in range(1, n + 1):
        line[i].setService(oclass, Exp.fitMean(means[i - 1]))  # Replace with correct expression

    # Block 3: topology
    model.link(Network.serialRouting(line))
    return model


def gallery_mm1_tandem():
    return gallery_mm1_linear(2)


def gallery_mmk(k=2):
    model = Network('M/M/k')
    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp(1))
    queue.setService(oclass, Exp(2))
    queue.setNumberOfServers(k)
    # Block 3: topology
    model.link(Network.serialRouting(source, queue, sink))
    return model


def gallery_mpar1():
    model = Network('M/Par/1')
    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp(1))
    queue.setService(oclass, Pareto.fitMeanAndSCV(0.5, 64))
    # Block 3: topology
    model.link(Network.serialRouting(source, queue, sink))
    return model


def gallery_parm1():
    model = Network('Par/M/1')
    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Pareto.fitMeanAndSCV(1, 64))
    queue.setService(oclass, Exp(2))
    # Block 3: topology
    model.link(Network.serialRouting(source, queue, sink))
    return model


def gallery_um1():
    model = Network('U/M/1')
    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Uniform(1, 2))
    queue.setService(oclass, Exp(2))
    # Block 3: topology
    model.link(Network.serialRouting(source, queue, sink))
    return model
