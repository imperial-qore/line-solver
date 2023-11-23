from line_solver import *


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
    for i in range(1,n+1):
        line.append(Queue(model, 'Queue' + str(i), SchedStrategy.FCFS))
    line.append(Sink(model, 'mySink'))

    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    line[0].setArrival(oclass, Exp(1.0))

    if n == 2: # linspace has a different behavior in np than matlab in this case
        means = np.linspace(Umax, Umax, 1)
    else:
        means = np.linspace(0.1, Umax, n // 2)

    if n % 2 == 0:
        means = np.concatenate([means, means[::-1]])
    else:
        means = np.concatenate([means, [Umax], means[::-1]])

    for i in range(1, n+1):
        line[i].setService(oclass, Exp.fitMean(means[i - 1]))  # Replace with correct expression

    # Block 3: topology
    model.link(Network.serialRouting(line))
    return model

def gallery_mm1_tandem():
    return gallery_mm1_linear(2)
