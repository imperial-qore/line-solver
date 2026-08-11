# Use relative imports to avoid circular import issues
import os
import random
import numpy as np
from .lang import (
    Network, Source, Queue, Sink, Delay, OpenClass, ClosedClass,
    Fork, Join, ClassSwitch, Router, Cache
)
from .environment import Environment
from .layered import (
    LayeredNetwork, Processor, Task, Entry, Activity, ActivityPrecedence,
    CacheTask, ItemEntry
)
from .distributions import (
    APH, Cox2, Coxian, Det, Disabled, Erlang, Exp, Gamma, HyperExp,
    Immediate, MAP, Pareto, PH, Replayer, Uniform, DiscreteSampler
)
from .constants import SchedStrategy, RoutingStrategy
from .lang.base import ReplacementStrategy
from .gen import NetworkGenerator, LayeredNetworkGenerator, cyclic_graph


def gallery_aphm1():
    """
    Create an APH/M/1 queueing model.
    
    Models a single-server queue with Acyclic Phase-type (APH) arrivals
    and exponential service times. Demonstrates advanced arrival process modeling.
    
    Returns:
        Network: APH/M/1 queueing network model.
    """
    model = Network('APH/M/1')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    alpha = [1, 0]
    T = [[-2, 1.5], [0, -1]]
    e = [[0.5], [1]]
    source.setArrival(oclass, APH(alpha, T, e))
    queue.setService(oclass, Exp(2))
    model.link(Network.serial_routing(source, queue, sink))
    return model


def gallery_coxm1():
    """
    Create a Cox/M/1 queueing model.
    
    Models a single-server queue with Coxian arrivals (fitted to high variability)
    and exponential service times. Used for modeling bursty arrival processes.
    
    Returns:
        Network: Cox/M/1 queueing network model with SCV=4.0 arrivals.
    """
    model = Network('Cox/M/1')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Coxian.fitMeanAndSCV(1.0, 4.0))
    queue.setService(oclass, Exp(2))
    model.link(Network.serial_routing(source, queue, sink))
    return model

def gallery_detm1():
    """
    Create a D/M/1 queueing model.
    
    Models a single-server queue with deterministic (constant) arrivals
    and exponential service times. Classic model for studying the effect
    of deterministic arrivals on queueing performance.
    
    Returns:
        Network: D/M/1 queueing network model.
    """
    model = Network('D/M/1')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Det(1))
    queue.setService(oclass, Exp(2))
    model.link(Network.serial_routing(source, queue, sink))
    return model


def gallery_erlm1():
    """
    Create an Erlang/M/1 queueing model.
    
    Models a single-server queue with 5-phase Erlang arrivals and exponential
    service times. Demonstrates low-variability arrival processes with
    coefficient of variation < 1.
    
    Returns:
        Network: Er/M/1 queueing network model with 5-phase Erlang arrivals.
    """
    model = Network('Er/M/1')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Erlang.fitMeanAndOrder(1, 5))
    queue.setService(oclass, Exp(2))
    model.link(Network.serial_routing(source, queue, sink))
    return model


def gallery_erlm1ps():
    """
    Create an Erlang/M/1 queue with Processor Sharing.
    
    Models a single-server queue with 5-phase Erlang arrivals, exponential
    service times, and processor sharing scheduling. Demonstrates PS scheduling
    with controlled-variance arrivals.
    
    Returns:
        Network: Er/M/1-PS queueing network model.
    """
    model = Network('Er/M/1-PS')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.PS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Erlang.fitMeanAndOrder(1, 5))
    queue.setService(oclass, Exp(2))
    model.link(Network.serial_routing(source, queue, sink))
    return model


def gallery_gamm1():
    """
    Create a Gamma/M/1 queueing model.
    
    Models a single-server queue with Gamma-distributed arrivals and exponential
    service times. Uses Gamma distribution fitted to mean=1, SCV=0.2 for
    flexible arrival process modeling.
    
    Returns:
        Network: Gamma/M/1 queueing network model.
    """
    model = Network('Gam/M/1')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Gamma.fitMeanAndSCV(1, 1 / 5))
    queue.setService(oclass, Exp(2))
    model.link(Network.serial_routing(source, queue, sink))
    return model


def gallery_hyperlk(k=2):
    """
    Create a HyperExp/Erlang/k queueing model.
    
    Models a multi-server queue with high-variability hyper-exponential arrivals
    and low-variability Erlang service times. Demonstrates the interaction
    between high-variance arrivals and controlled-variance service.
    
    Args:
        k (int): Number of servers (default: 2).
        
    Returns:
        Network: H/Er/k queueing network model.
    """
    model = Network('H/Er/k')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, HyperExp.fitMeanAndSCVBalanced(1.0 / 1.8, 4))
    queue.setService(oclass, Erlang.fitMeanAndSCV(1, 0.25))
    queue.setNumberOfServers(k)
    model.link(Network.serial_routing(source, queue, sink))
    return model


def gallery_hypm1():
    """
    Create a HyperExp/M/1 queueing model.
    
    Models a single-server queue with extremely high-variability hyper-exponential
    arrivals (SCV=64) and exponential service times. Demonstrates modeling of
    very bursty arrival processes.
    
    Returns:
        Network: H/M/1 queueing network model with very high-variance arrivals.
    """
    model = Network('H/M/1')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, HyperExp.fitMeanAndSCV(1, 64))
    queue.setService(oclass, Exp(2))
    model.link(Network.serial_routing(source, queue, sink))
    return model


def gallery_mhyp1():
    model = Network('M/H/1')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp(1))
    queue.setService(oclass, Coxian.fitMeanAndSCV(0.5, 4))
    model.link(Network.serial_routing(source, queue, sink))
    return model


def gallery_merl1():
    model = Network('M/E/1')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp(1))
    queue.setService(oclass, Erlang.fitMeanAndOrder(0.5, 2))
    model.link(Network.serial_routing(source, queue, sink))
    return model, source, queue, sink, oclass


def gallery_mm1():
    model = Network('M/M/1')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp(1))
    queue.setService(oclass, Exp(2))
    model.link(Network.serial_routing(source, queue, sink))
    return model


def gallery_mm1_linear(n=2, Umax=0.9):
    """
    Create a linear tandem network of M/M/1 queues.
    
    Models a series of single-server queues in tandem, with utilizations
    that form a pattern (increasing then decreasing). Used for studying
    the behavior of jobs flowing through multiple service stages.
    
    Args:
        n (int): Number of queues in the tandem (default: 2).
        Umax (float): Maximum utilization level (default: 0.9).
        
    Returns:
        Network: Linear tandem network with n M/M/1 queues.
    """
    model = Network('M/M/1-Linear')

    line = [Source(model, 'mySource')]
    for i in range(1, n + 1):
        line.append(Queue(model, 'Queue' + str(i), SchedStrategy.FCFS))
    line.append(Sink(model, 'mySink'))

    oclass = OpenClass(model, 'myClass')
    line[0].setArrival(oclass, Exp(1.0))

    if n == 2:
        means = np.linspace(Umax, Umax, 1)
    else:
        means = np.linspace(0.1, Umax, n // 2)

    if n % 2 == 0:
        means = np.concatenate([means, means[::-1]])
    else:
        means = np.concatenate([means, [Umax], means[::-1]])

    for i in range(1, n + 1):
        line[i].setService(oclass, Exp.fitMean(means[i - 1]))

    model.link(Network.serial_routing(line))
    return model


def gallery_mm1_tandem():
    """
    Create a simple 2-queue M/M/1 tandem network.
    
    Convenience function that creates a 2-queue linear tandem network
    by calling gallery_mm1_linear(2). Represents the basic tandem
    queueing system.
    
    Returns:
        Network: 2-queue M/M/1 tandem network.
    """
    return gallery_mm1_linear(2)


def gallery_mmk(k=2):
    """
    Create an M/M/k multi-server queueing model.
    
    Models a multi-server queue with Poisson arrivals, exponential service times,
    and k identical servers. Demonstrates the performance benefits of
    multiple servers versus a single fast server.
    
    Args:
        k (int): Number of servers (default: 2).
        
    Returns:
        Network: M/M/k queueing network model.
    """
    model = Network('M/M/k')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp(1))
    queue.setService(oclass, Exp(2))
    queue.setNumberOfServers(k)
    model.link(Network.serial_routing(source, queue, sink))
    return model


def gallery_mpar1():
    """
    Create an M/Pareto/1 queueing model.
    
    Models a single-server queue with Poisson arrivals and heavy-tailed
    Pareto-distributed service times. Demonstrates modeling of service
    processes with very high variability and infinite variance.
    
    Returns:
        Network: M/Par/1 queueing network model with Pareto service times.
    """
    model = Network('M/Par/1')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp(1))
    queue.setService(oclass, Pareto.fitMeanAndSCV(0.5, 64))
    model.link(Network.serial_routing(source, queue, sink))
    return model


def gallery_parm1():
    """
    Create a Pareto/M/1 queueing model.
    
    Models a single-server queue with heavy-tailed Pareto arrivals and
    exponential service times. Demonstrates modeling of bursty arrival
    processes with power-law characteristics.
    
    Returns:
        Network: Par/M/1 queueing network model with Pareto arrivals.
    """
    model = Network('Par/M/1')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Pareto.fitMeanAndSCV(1, 64))
    queue.setService(oclass, Exp(2))
    model.link(Network.serial_routing(source, queue, sink))
    return model


def gallery_um1():
    """
    Create a Uniform/M/1 queueing model.
    
    Models a single-server queue with uniformly distributed arrivals
    and exponential service times. Demonstrates modeling with
    bounded inter-arrival times.
    
    Returns:
        Network: U/M/1 queueing network model with uniform arrivals.
    """
    model = Network('U/M/1')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Uniform(1, 2))
    queue.setService(oclass, Exp(2))
    model.link(Network.serial_routing(source, queue, sink))
    return model


def gallery_cqn(M=2, useDelay=False, seed=2300):
    """
    Create a closed queueing network (CQN) model.
    
    Models a closed network with fixed population, where jobs circulate
    between service stations. Can use either delay stations (infinite servers)
    or finite capacity queues depending on the useDelay parameter.
    
    Args:
        M (int): Number of stations in the network (default: 2).
        useDelay (bool): Whether to use delay stations (default: False).
        seed (int): Random seed for reproducible results (default: 2300).
        
    Returns:
        Network: Closed queueing network model.
    """
    model = Network('CQN')

    stations = []
    for i in range(M):
        station = Queue(model, f'Queue{i+1}', SchedStrategy.PS)
        stations.append(station)

    if useDelay:
        delay = Delay(model, 'Delay')
        stations.append(delay)

    refStation = stations[0] if not useDelay else stations[-1]
    jobclass = ClosedClass(model, 'Jobs', 20, refStation)

    np.random.seed(seed)
    for i, station in enumerate(stations):
        if isinstance(station, Queue):
            rate = 0.1 + 0.9 * np.random.random()
            station.setService(jobclass, Exp(rate))
        elif isinstance(station, Delay):
            station.setService(jobclass, Exp(0.1))

    if len(stations) == 1:
        model.link(Network.selfRouting(stations[0]))
    else:
        model.link(Network.serial_routing(stations))

    return model


def gallery_mm1_feedback(p=0.5):
    """
    Create an M/M/1 queue with probabilistic feedback.
    
    Models a single-server queue where jobs have probability p of returning
    to the queue after service completion, creating a feedback loop.
    This increases the effective service demand and response time.
    
    Args:
        p (float): Feedback probability (default: 0.5).
        
    Returns:
        Network: M/M/1 queueing network with probabilistic feedback.
    """
    model = Network('M/M/1-Feedback')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp(1))
    queue.setService(oclass, Exp(2))
    model.link(Network.serial_routing(source, queue, sink))
    return model


def gallery_mm1_prio():
    """
    Create an M/M/1 queue with priority classes.
    
    Models a single-server queue with two job classes having different
    priorities. High priority jobs are served before low priority jobs,
    demonstrating head-of-line priority scheduling.
    
    Returns:
        Network: M/M/1 queueing network with high and low priority classes.
    """
    model = Network('M[2]/M[2]/1')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.HOL)
    sink = Sink(model, 'mySink')

    oclass1 = OpenClass(model, 'myClass1', 1)
    source.setArrival(oclass1, Exp(1))
    queue.setService(oclass1, Exp(4))

    oclass2 = OpenClass(model, 'myClass2', 0)
    source.setArrival(oclass2, Exp(0.5))
    queue.setService(oclass2, Exp(4))

    P = model.init_routing_matrix()
    P[oclass1] = Network.serial_routing([source, queue, sink])
    P[oclass2] = Network.serial_routing([source, queue, sink])
    model.link(P)
    return model


def gallery_mm1_multiclass():
    """
    Create an M/M/1 queue with multiple job classes.
    
    Models a single-server queue with two different job classes arriving
    from separate sources, each with different arrival rates and service
    requirements. Demonstrates multi-class queueing behavior.
    
    Returns:
        Network: M/M/1 multi-class queueing network model.
    """
    model = Network('M/M/1-MultiClass')
    source1 = Source(model, 'Source1')
    source2 = Source(model, 'Source2')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')

    class1 = OpenClass(model, 'Class1')
    class2 = OpenClass(model, 'Class2')

    source1.setArrival(class1, Exp(0.4))
    source2.setArrival(class2, Exp(0.6))
    queue.setService(class1, Exp(1.5))
    queue.setService(class2, Exp(1.0))

    model.link(Network.serial_routing(source1, queue, sink))
    model.link(Network.serial_routing(source2, queue, sink))
    return model


def gallery_mapm1(map_arrival=None):
    """
    Create a MAP/M/1 queueing model with Markovian arrival process.
    
    Models a single-server queue with a Markovian Arrival Process (MAP)
    and exponential service times. MAP allows modeling of correlated
    arrivals and more complex arrival patterns than Poisson processes.
    
    Args:
        map_arrival: MAP arrival process (default: None, creates a standard MAP).
        
    Returns:
        Network: MAP/M/1 queueing network model.
    """
    model = Network('MAP/M/1')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')

    if map_arrival is None:
        D0 = [[-2, 0], [0, -1]]
        D1 = [[1.5, 0.5], [0.8, 0.2]]
        map_arrival = MAP(D0, D1)

    source.setArrival(oclass, map_arrival)
    queue.setService(oclass, Exp(2))
    model.link(Network.serial_routing(source, queue, sink))
    return model


def gallery_dm1():
    model = Network('D/M/1')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass1 = OpenClass(model, 'Class1')
    source.setArrival(oclass1, Det(1))
    queue.setService(oclass1, Exp(2))
    model.link(Network.serial_routing([source, queue, sink]))
    return model


def gallery_erldk(k=2):
    model = Network('Erl/D/k')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Erlang.fit_mean_and_order(1, 5))
    queue.setService(oclass, Det(2 / k))
    queue.setNumberOfServers(k)
    model.link(Network.serial_routing([source, queue, sink]))
    return model


def gallery_erlerl1(n=5):
    model = Network('Erl/Erl/1')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass1 = OpenClass(model, 'Class1')
    oclass2 = OpenClass(model, 'Class2')
    source.setArrival(oclass1, Erlang.fit_mean_and_order(1, n))
    source.setArrival(oclass2, Disabled())
    queue.setService(oclass1, Erlang.fit_mean_and_order(0.5, n))
    queue.setService(oclass2, Exp(3))
    P = model.init_routing_matrix()
    P.set(oclass1, oclass1, source, queue, 1)
    P.set(oclass2, oclass2, queue, sink, 1)
    model.link(P)
    return model


def gallery_erlerl1_reentrant():
    model = Network('Erl/Erl/1-Reentrant')
    n = 5
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass1 = OpenClass(model, 'Class1')
    oclass2 = OpenClass(model, 'Class2')
    source.setArrival(oclass1, Erlang.fit_mean_and_order(1, n))
    source.setArrival(oclass2, Disabled())
    queue.setService(oclass1, Erlang.fit_mean_and_order(0.5, n))
    queue.setService(oclass2, Exp(3))
    P = model.init_routing_matrix()
    P.set(oclass1, oclass1, source, queue, 1)
    P.set(oclass1, oclass2, queue, queue, 0.5)
    P.set(oclass2, oclass2, queue, sink, 0.5)
    model.link(P)
    return model


def gallery_erlm1_ps():
    model = Network('Er/M/1-PS')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.PS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Erlang.fit_mean_and_order(1, 5))
    queue.setService(oclass, Exp(2))
    model.link(Network.serial_routing([source, queue, sink]))
    return model


def gallery_erlm1_reentrant():
    model = Network('Er/M/1-Reentrant')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass1 = OpenClass(model, 'Class1')
    oclass2 = OpenClass(model, 'Class2')
    source.setArrival(oclass1, Erlang.fit_mean_and_order(1, 5))
    source.setArrival(oclass2, Disabled())
    queue.setService(oclass1, Exp(2))
    queue.setService(oclass2, Exp(3))
    P = model.init_routing_matrix()
    P.set(oclass1, oclass1, source, queue, 1)
    P.set(oclass1, oclass2, queue, queue, 1)
    P.set(oclass2, oclass2, queue, sink, 1)
    model.link(P)
    return model


def gallery_hyperl1_feedback():
    model = Network('Hyper/Erl/1-Feedback')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass1 = OpenClass(model, 'Class1')
    source.setArrival(oclass1, HyperExp.fit_mean_and_scv(1, 64))
    queue.setService(oclass1, Erlang.fit_mean_and_order(0.5, 5))
    P = model.init_routing_matrix()
    P.set(oclass1, oclass1, source, queue, 1)
    P.set(oclass1, oclass1, queue, queue, 0.9)
    P.set(oclass1, oclass1, queue, sink, 0.1)
    model.link(P)
    return model


def gallery_hyperl1_reentrant():
    model = Network('Hyper/Erl/1-Reentrant')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass1 = OpenClass(model, 'Class1')
    oclass2 = OpenClass(model, 'Class2')
    source.setArrival(oclass1, HyperExp.fit_mean_and_scv(1, 64))
    source.setArrival(oclass2, Disabled())
    queue.setService(oclass1, Erlang.fit_mean_and_order(0.5, 5))
    queue.setService(oclass2, Exp(3))
    P = model.init_routing_matrix()
    P.set(oclass1, oclass1, source, queue, 1)
    P.set(oclass1, oclass2, queue, queue, 1)
    P.set(oclass2, oclass2, queue, sink, 1)
    model.link(P)
    return model


def gallery_hyphyp1_linear(n=2, Umax=0.9):
    model = Network('Hyp/Hyp/1-Linear')
    line = [Source(model, 'mySource')]
    for i in range(n):
        line.append(Queue(model, f'Queue{i+1}', SchedStrategy.FCFS))
    line.append(Sink(model, 'mySink'))
    oclass = OpenClass(model, 'myClass')
    line[0].setArrival(oclass, HyperExp.fit_mean_and_scv(1, 2))
    means = [Umax] if n // 2 == 1 else list(np.linspace(0.1, Umax, n // 2))
    if n % 2 == 0:
        means = means + means[::-1]
    else:
        means = means + [Umax] + means[::-1]
    for i in range(n):
        line[i + 1].setService(oclass, HyperExp.fit_mean_and_scv(means[i], 1 + i + 1))
    model.link(Network.serial_routing(line))
    return model


def gallery_hyphyp1_reentrant():
    model = Network('Hyper/Hyper/1-Reentrant')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    queue.setNumberOfServers(2)
    sink = Sink(model, 'Sink')
    oclass1 = OpenClass(model, 'Class1')
    oclass2 = OpenClass(model, 'Class2')
    source.setArrival(oclass1, HyperExp.fit_mean_and_scv(1, 64))
    source.setArrival(oclass2, Disabled())
    queue.setService(oclass1, HyperExp.fit_mean_and_scv(0.5, 4))
    queue.setService(oclass2, Exp(3))
    P = model.init_routing_matrix()
    P.set(oclass1, oclass1, source, queue, 1)
    P.set(oclass1, oclass2, queue, queue, 1)
    P.set(oclass2, oclass2, queue, sink, 1)
    model.link(P)
    return model


def gallery_hyphyp1_tandem():
    return gallery_hyphyp1_linear(2)


def gallery_hypm1_reentrant():
    model = Network('Hyper/M/1-Reentrant')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass1 = OpenClass(model, 'Class1')
    oclass2 = OpenClass(model, 'Class2')
    source.setArrival(oclass1, HyperExp.fit_mean_and_scv(1, 4))
    source.setArrival(oclass2, Disabled())
    queue.setService(oclass1, Exp(2))
    queue.setService(oclass2, Exp(3))
    P = model.init_routing_matrix()
    P.set(oclass1, oclass1, source, queue, 1)
    P.set(oclass1, oclass2, queue, queue, 1)
    P.set(oclass2, oclass2, queue, sink, 1)
    model.link(P)
    return model


def gallery_lukumar_reentrant(sched_strategy='FCFS'):
    sched_strategy = sched_strategy.upper()
    if sched_strategy == 'HOL':
        sched1 = sched2 = SchedStrategy.HOL
    elif sched_strategy == 'PS':
        sched1 = sched2 = SchedStrategy.PS
    else:
        sched1 = sched2 = SchedStrategy.FCFS
    model = Network('Lu-Kumar-Reentrant')
    source = Source(model, 'Source')
    station1 = Queue(model, 'Station1', sched1)
    station2 = Queue(model, 'Station2', sched2)
    sink = Sink(model, 'Sink')
    class1 = OpenClass(model, 'Class1', 1)
    class2 = OpenClass(model, 'Class2', 0)
    class3 = OpenClass(model, 'Class3', 1)
    class4 = OpenClass(model, 'Class4', 0)
    arrival_rate = 0.08
    source.setArrival(class1, Exp(arrival_rate))
    source.setArrival(class2, Disabled())
    source.setArrival(class3, Exp(arrival_rate))
    source.setArrival(class4, Disabled())
    m1, m2, m3, m4 = 10.0, 1.0, 10.0, 1.0
    station1.setService(class1, Exp(1 / m1))
    station1.setService(class2, Disabled())
    station1.setService(class3, Disabled())
    station1.setService(class4, Exp(1 / m4))
    station2.setService(class1, Disabled())
    station2.setService(class2, Exp(1 / m2))
    station2.setService(class3, Exp(1 / m3))
    station2.setService(class4, Disabled())
    P = model.init_routing_matrix()
    P.set(class1, class1, source, station1, 1)
    P.set(class1, class2, station1, station2, 1)
    P.set(class2, class2, station2, sink, 1)
    P.set(class3, class3, source, station2, 1)
    P.set(class3, class4, station2, station1, 1)
    P.set(class4, class4, station1, sink, 1)
    model.link(P)
    return model


def gallery_mapmk(map=None, k=2, seed=23000):
    if map is None:
        map = MAP.rand(seed=seed)
    model = Network('MAP/M/k')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, map)
    queue.setService(oclass, Exp(2))
    queue.setNumberOfServers(k)
    model.link(Network.serial_routing([source, queue, sink]))
    return model


def gallery_mdk(k=2):
    model = Network('M/D/k')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp.fit_mean(1))
    queue.setService(oclass, Det(2 / k))
    queue.setNumberOfServers(k)
    model.link(Network.serial_routing([source, queue, sink]))
    return model


def gallery_merl1_linear(n=2, Umax=0.9):
    model = Network('M/Erl/1-Linear')
    line = [Source(model, 'mySource')]
    for i in range(n):
        line.append(Queue(model, f'Queue{i+1}', SchedStrategy.FCFS))
    line.append(Sink(model, 'mySink'))
    oclass = OpenClass(model, 'myClass')
    line[0].setArrival(oclass, Exp(1))
    means = [Umax] if n // 2 == 1 else list(np.linspace(0.1, Umax, n // 2))
    if n % 2 == 0:
        means = means + means[::-1]
    else:
        means = means + [Umax] + means[::-1]
    for i in range(n):
        line[i + 1].setService(oclass, Erlang.fit_mean_and_order(means[i], i + 1))
    model.link(Network.serial_routing(line))
    return model


def gallery_merl1_reentrant():
    model = Network('M/Erl/1-Reentrant')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass1 = OpenClass(model, 'Class1')
    oclass2 = OpenClass(model, 'Class2')
    source.setArrival(oclass1, Exp(1))
    source.setArrival(oclass2, Disabled())
    queue.setService(oclass1, Erlang.fit_mean_and_order(1 / 2, 5))
    queue.setService(oclass2, Exp(3))
    P = model.init_routing_matrix()
    P.set(oclass1, oclass1, source, queue, 1)
    P.set(oclass1, oclass2, queue, queue, 1)
    P.set(oclass2, oclass2, queue, sink, 1)
    model.link(P)
    return model


def gallery_merl1_tandem():
    return gallery_merl1_linear(2)


def gallery_merlk(k=2):
    model = Network(f'M/Erl/{k}')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp(1))
    queue.setService(oclass, Erlang.fit_mean_and_order(0.5, 2))
    queue.setNumberOfServers(k)
    model.link(Network.serial_routing([source, queue, sink]))
    return model


def gallery_mhyp1_linear(n=2, Umax=0.9):
    model = Network('M/Hyp/1-Linear')
    line = [Source(model, 'mySource')]
    for i in range(n):
        line.append(Queue(model, f'Queue{i+1}', SchedStrategy.FCFS))
    line.append(Sink(model, 'mySink'))
    oclass = OpenClass(model, 'myClass')
    line[0].setArrival(oclass, Exp(1))
    means = [Umax] if n // 2 == 1 else list(np.linspace(0.1, Umax, n // 2))
    if n % 2 == 0:
        means = means + means[::-1]
    else:
        means = means + [Umax] + means[::-1]
    for i in range(n):
        line[i + 1].setService(oclass, HyperExp.fit_mean_and_scv(means[i], n))
    model.link(Network.serial_routing(line))
    return model


def gallery_mhyp1_reentrant():
    model = Network('M/Hyper/1-Reentrant')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass1 = OpenClass(model, 'Class1')
    oclass2 = OpenClass(model, 'Class2')
    source.setArrival(oclass1, Exp(1))
    source.setArrival(oclass2, Disabled())
    queue.setService(oclass1, Coxian.fit_mean_and_scv(0.5, 4))
    queue.setService(oclass2, Exp(3))
    P = model.init_routing_matrix()
    P.set(oclass1, oclass1, source, queue, 1)
    P.set(oclass1, oclass2, queue, queue, 1)
    P.set(oclass2, oclass2, queue, sink, 1)
    model.link(P)
    return model


def gallery_mhyp1_tandem():
    return gallery_mhyp1_linear(2)


def gallery_mhypk(k=2):
    model = Network(f'M/Hyper/{k}')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp(1))
    queue.setService(oclass, Coxian.fit_mean_and_scv(0.5, 4))
    queue.setNumberOfServers(k)
    model.link(Network.serial_routing([source, queue, sink]))
    return model


def gallery_mm1_ps():
    model = Network('M/M/1-PS')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.PS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp(1))
    queue.setService(oclass, Exp(2))
    model.link(Network.serial_routing([source, queue, sink]))
    return model


def gallery_mm1_ps_feedback(p=1 / 3):
    model = Network('M/M/1-PS-Feedback')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.PS)
    sink = Sink(model, 'Sink')
    oclass1 = OpenClass(model, 'Class1')
    source.setArrival(oclass1, Exp.fit_mean(1))
    queue.setService(oclass1, Exp.fit_mean(0.5))
    P = model.init_routing_matrix()
    P.set(oclass1, oclass1, source, queue, 1)
    P.set(oclass1, oclass1, queue, queue, p)
    P.set(oclass1, oclass1, queue, sink, 1 - p)
    model.link(P)
    return model


def gallery_mm1_ps_multiclass():
    model = Network('M[2]/M[2]/1-PS')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.PS)
    sink = Sink(model, 'mySink')
    oclass1 = OpenClass(model, 'myClass1')
    source.setArrival(oclass1, Exp(1))
    queue.setService(oclass1, Exp(4))
    oclass2 = OpenClass(model, 'myClass2')
    source.setArrival(oclass2, Exp(0.5))
    queue.setService(oclass2, Exp(4))
    P = model.init_routing_matrix()
    P[oclass1] = Network.serial_routing([source, queue, sink])
    P[oclass2] = Network.serial_routing([source, queue, sink])
    model.link(P)
    return model


def gallery_mm1_ps_reentrant():
    model = Network('M/M/1')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.PS)
    sink = Sink(model, 'Sink')
    oclass1 = OpenClass(model, 'Class1')
    oclass2 = OpenClass(model, 'Class2')
    source.setArrival(oclass1, Exp(1))
    source.setArrival(oclass2, Disabled())
    queue.setService(oclass1, Exp(2))
    queue.setService(oclass2, Exp(3))
    P = model.init_routing_matrix()
    P.set(oclass1, oclass1, source, queue, 1)
    P.set(oclass1, oclass2, queue, queue, 1)
    P.set(oclass2, oclass2, queue, sink, 1)
    model.link(P)
    return model


def gallery_mm1_reentrant():
    model = Network('M/M/1-Reentrant')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass1 = OpenClass(model, 'Class1')
    oclass2 = OpenClass(model, 'Class2')
    source.setArrival(oclass1, Exp(1))
    source.setArrival(oclass2, Disabled())
    queue.setService(oclass1, Exp(2))
    queue.setService(oclass2, Exp(3))
    P = model.init_routing_matrix()
    P.set(oclass1, oclass1, source, queue, 1)
    P.set(oclass1, oclass2, queue, queue, 1)
    P.set(oclass2, oclass2, queue, sink, 1)
    model.link(P)
    return model


def gallery_mm1_tandem_multiclass():
    model = Network('M[2]/M[2]/1 -> -/M[2]/1')
    source = Source(model, 'Source')
    queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass1 = OpenClass(model, 'myClass1')
    source.setArrival(oclass1, Exp(1))
    queue1.setService(oclass1, Exp(4))
    queue2.setService(oclass1, Exp(6))
    oclass2 = OpenClass(model, 'myClass2')
    source.setArrival(oclass2, Exp(0.5))
    queue1.setService(oclass2, Exp(2))
    queue2.setService(oclass2, Exp(6))
    P = model.init_routing_matrix()
    P[oclass1] = Network.serial_routing([source, queue1, queue2, sink])
    P[oclass2] = Network.serial_routing([source, queue1, queue2, sink])
    model.link(P)
    return model


def gallery_cqn_multiclass(m=1, r=2, wantdelay=True, seed=23000):
    random.seed(seed)
    model = Network('Multi-class CQN')
    node = []
    for i in range(m):
        node.append(Queue(model, f'Queue {i+1}', SchedStrategy.PS))
    if wantdelay:
        node.append(Delay(model, 'Delay 1'))
    jobclass = []
    for s in range(r):
        jobclass.append(ClosedClass(model, f'Class{s+1}', 5, node[0], 0))
    for s in range(r):
        for i in range(m):
            node[i].setService(jobclass[s], Exp.fit_mean(round(50 * random.random())))
        if wantdelay:
            node[-1].setService(jobclass[s], Exp.fit_mean(round(100 * random.random())))
    P = model.init_routing_matrix()
    for s in range(r):
        P[jobclass[s], jobclass[s]] = Network.serial_routing(node)
    model.link(P)
    return model


def gallery_repairmen(nservers=1, seed=2300):
    model = Network('Finite repairmen CQN')
    M = 1
    random.seed(seed)
    node = []
    node.append(Queue(model, 'Queue1', SchedStrategy.PS))
    node[0].setNumberOfServers(nservers)
    node.append(Delay(model, 'Delay1'))
    jobclass = ClosedClass(model, 'Class1', round(random.random() * 10 * M + 3), node[0], 0)
    node[0].setService(jobclass, Exp.fit_mean(random.random() + 1))
    node[1].setService(jobclass, Exp.fit_mean(2.0))
    P = model.init_routing_matrix()
    P[jobclass, jobclass] = Network.serial_routing(node)
    model.link(P)
    return model


def gallery_mmap1(map=None, seed=23000):
    if map is None:
        map = MAP.rand(seed=seed)
        map = map.set_mean(0.5)
    model = Network('M/MAP/1')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp(1))
    queue.setService(oclass, map)
    model.link(Network.serial_routing([source, queue, sink]))
    return model


def gallery_mmap1_multiclass(map1=None, map2=None, seed=23000):
    if map1 is None:
        map1 = MAP.rand(n=2, seed=seed)
        map1 = map1.set_mean(0.5)
    if map2 is None:
        map2 = MAP.rand(n=3, seed=seed + 1)
        map2 = map2.set_mean(0.5)
    model = Network('M/MAP/1')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass1 = OpenClass(model, 'myClass1')
    source.setArrival(oclass1, Exp(0.35 / map1.get_mean()))
    queue.setService(oclass1, map1)
    oclass2 = OpenClass(model, 'myClass2')
    source.setArrival(oclass2, Exp(0.15 / map2.get_mean()))
    queue.setService(oclass2, map2)
    P = model.init_routing_matrix()
    P[oclass1] = Network.serial_routing([source, queue, sink])
    P[oclass2] = Network.serial_routing([source, queue, sink])
    model.link(P)
    return model


def gallery_mmapk(map=None, k=2, seed=23000):
    if map is None:
        map = MAP.rand(seed=seed)
    model = Network('M/MAP/k')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp(1))
    queue.setService(oclass, map)
    queue.setNumberOfServers(k)
    model.link(Network.serial_routing([source, queue, sink]))
    return model


def gallery_replayerm1(filename=None):
    if filename is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        possible_paths = [
            os.path.join(script_dir, '..', 'examples', 'gettingstarted', 'example_trace.txt'),
            os.path.join(script_dir, '..', 'examples', 'basic', 'openQN', 'example_trace.txt'),
        ]
        for path in possible_paths:
            if os.path.exists(path):
                filename = path
                break
        if filename is None:
            raise FileNotFoundError("example_trace.txt not found in expected locations")
    model = Network('Trace/M/1')
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')
    oclass = OpenClass(model, 'myClass')
    replayer = Replayer(filename)
    source.setArrival(oclass, replayer)
    queue.setService(oclass, Exp(3 / replayer.get_mean()))
    model.link(Network.serial_routing([source, queue, sink]))
    return model


def gallery_multitier():
    """
    Create a 4-tier J2EE Layered Queueing Network (client/app/database).

    Reference LQN with multiple entries per task and or-fork/or-join
    activity precedence. Ported from MATLAB gallery_multitier.

    Returns:
        LayeredNetwork: 3-layer (client, application, database) LQN model.
    """
    model = LayeredNetwork('testLQN3')

    # Layer 1: client
    P0 = Processor(model, 'P0', 1, SchedStrategy.PS)
    T0 = Task(model, 'T0', 1, SchedStrategy.REF).on(P0)
    E0 = Entry(model, 'E0').on(T0)

    # Layer 2: application server
    P1 = Processor(model, 'P1', 1, SchedStrategy.PS)
    T1 = Task(model, 'T1', 1, SchedStrategy.FCFS).on(P1)
    E10 = Entry(model, 'E10').on(T1)
    E11 = Entry(model, 'E11').on(T1)
    E12 = Entry(model, 'E12').on(T1)
    E13 = Entry(model, 'E13').on(T1)

    # Layer 3: database
    P2 = Processor(model, 'P2', 1, SchedStrategy.PS)
    T2 = Task(model, 'T2', 1, SchedStrategy.FCFS).on(P2)
    E20 = Entry(model, 'E20').on(T2)
    E21 = Entry(model, 'E21').on(T2)
    E22 = Entry(model, 'E22').on(T2)
    E23 = Entry(model, 'E23').on(T2)

    # Client activities
    A0 = Activity(model, 'A0', Exp(1.0)).on(T0).bound_to(E0).synch_call(E12, 1.0)
    A1 = Activity(model, 'A1', Exp(1.0)).on(T0).synch_call(E10, 1.0)
    A2 = Activity(model, 'A2', Exp(1.0)).on(T0).synch_call(E11, 1.0)
    A3 = Activity(model, 'A3', Exp(1.0)).on(T0).synch_call(E13, 1.0)

    # Application activities
    B0 = Activity(model, 'B0', Exp(1.0)).on(T1).bound_to(E10)
    B1 = Activity(model, 'B1', Exp(1.0)).on(T1).replies_to(E10)
    B2 = Activity(model, 'B2', Exp(1.0)).on(T1).bound_to(E11)
    B3 = Activity(model, 'B3', Exp(1.0)).on(T1).synch_call(E21, 1.0).replies_to(E11)
    B4 = Activity(model, 'B4', Exp(1.0)).on(T1).bound_to(E12).synch_call(E20, 1.0).replies_to(E12)
    B5 = Activity(model, 'B5', Exp(1.0)).on(T1).bound_to(E13)
    B6 = Activity(model, 'B6', Exp(1.0)).on(T1)
    B7 = Activity(model, 'B7', Exp(1.0)).on(T1).synch_call(E22, 1.0)
    B7a = Activity(model, 'B7a', Exp(1.0)).on(T1)
    B7b = Activity(model, 'B7b', Exp(1.0)).on(T1).synch_call(E23, 1.0)
    B8 = Activity(model, 'B8', Exp(1.0)).on(T1).replies_to(E13)

    # Database activities
    C0 = Activity(model, 'C0', Exp(1.0)).on(T2).bound_to(E20)
    C1 = Activity(model, 'C1', Exp(1.0)).on(T2).replies_to(E20)
    C2 = Activity(model, 'C2', Exp(1.0)).on(T2).bound_to(E21).replies_to(E21)
    C3 = Activity(model, 'C3', Exp(1.0)).on(T2).bound_to(E22)
    C4 = Activity(model, 'C4', Exp(1.0)).on(T2)
    C5 = Activity(model, 'C5', Exp(1.0)).on(T2).replies_to(E22)
    C6 = Activity(model, 'C6', Exp(1.0)).on(T2).bound_to(E23).replies_to(E23)

    # Precedences
    T0.add_precedence(ActivityPrecedence.Serial(A0, A1, A2, A3))
    T1.add_precedence(ActivityPrecedence.Serial(B0, B1))
    T1.add_precedence(ActivityPrecedence.Serial(B2, B3))
    T1.add_precedence(ActivityPrecedence.Serial(B5, B6, B7))
    T1.add_precedence(ActivityPrecedence.OrFork(B7, [B7a, B7b], [0.7, 0.3]))
    T1.add_precedence(ActivityPrecedence.OrJoin([B7a, B7b], B8))
    T2.add_precedence(ActivityPrecedence.Serial(C0, C1))
    T2.add_precedence(ActivityPrecedence.Serial(C3, C4, C5))

    return model


def gallery_multitier_storage():
    """
    Create a 4-tier J2EE LQN with a dedicated cache layer (LRU replacement).

    Extends gallery_multitier with a cache layer (CacheTask + ItemEntry) and
    cache-access hit/miss precedence. Ported from MATLAB gallery_multitier_storage.

    Returns:
        LayeredNetwork: 4-layer LQN model with cache layer.
    """
    model = LayeredNetwork('testLQN3_Cache')

    # Layer 1: client
    P0 = Processor(model, 'P0', 1, SchedStrategy.PS)
    T0 = Task(model, 'T0', 1, SchedStrategy.REF).on(P0)
    E0 = Entry(model, 'E0').on(T0)

    # Layer 2: application server
    P1 = Processor(model, 'P1', 1, SchedStrategy.PS)
    T1 = Task(model, 'T1', 1, SchedStrategy.FCFS).on(P1)
    E10 = Entry(model, 'E10').on(T1)
    E11 = Entry(model, 'E11').on(T1)
    E12 = Entry(model, 'E12').on(T1)
    E13 = Entry(model, 'E13').on(T1)

    # Layer 3: database
    P2 = Processor(model, 'P2', 1, SchedStrategy.PS)
    T2 = Task(model, 'T2', 1, SchedStrategy.FCFS).on(P2)
    E20 = Entry(model, 'E20').on(T2)
    E21 = Entry(model, 'E21').on(T2)
    E22 = Entry(model, 'E22').on(T2)
    E23 = Entry(model, 'E23').on(T2)

    # Layer 4: cache (10 items, capacity 2, LRU, uniform access)
    totalitems = 10
    cachecapacity = 2
    pAccess = DiscreteSampler([1.0 / totalitems] * totalitems)
    P3 = Processor(model, 'P3', 1, SchedStrategy.PS)
    T3 = CacheTask(model, 'T3', totalitems, cachecapacity, ReplacementStrategy.LRU, 1).on(P3)
    E3 = ItemEntry(model, 'E3', totalitems, pAccess).on(T3)

    # Client activities
    A0 = Activity(model, 'A0', Exp(1.0)).on(T0).bound_to(E0).synch_call(E12, 1.0)
    A1 = Activity(model, 'A1', Exp(1.0)).on(T0).synch_call(E10, 1.0)
    A2 = Activity(model, 'A2', Exp(1.0)).on(T0).synch_call(E11, 1.0)
    A3 = Activity(model, 'A3', Exp(1.0)).on(T0).synch_call(E13, 1.0)

    # Application activities
    B0 = Activity(model, 'B0', Exp(1.0)).on(T1).bound_to(E10)
    B1 = Activity(model, 'B1', Exp(1.0)).on(T1).replies_to(E10)
    B2 = Activity(model, 'B2', Exp(1.0)).on(T1).bound_to(E11)
    B3 = Activity(model, 'B3', Exp(1.0)).on(T1).synch_call(E21, 1.0).replies_to(E11)
    B4 = Activity(model, 'B4', Exp(1.0)).on(T1).bound_to(E12).synch_call(E20, 1.0).replies_to(E12)
    B5 = Activity(model, 'B5', Exp(1.0)).on(T1).bound_to(E13)
    B6 = Activity(model, 'B6', Exp(1.0)).on(T1)
    B7 = Activity(model, 'B7', Exp(1.0)).on(T1).synch_call(E22, 1.0)
    B7a = Activity(model, 'B7a', Exp(1.0)).on(T1).replies_to(E13)
    B7b = Activity(model, 'B7b', Exp(1.0)).on(T1).synch_call(E23, 1.0).replies_to(E13)

    # Database activities
    C0 = Activity(model, 'C0', Exp(1.0)).on(T2).bound_to(E20)
    C1 = Activity(model, 'C1', Exp(1.0)).on(T2).synch_call(E3, 1.0).replies_to(E20)
    C2 = Activity(model, 'C2', Exp(1.0)).on(T2).bound_to(E21).synch_call(E3, 1.0).replies_to(E21)
    C3 = Activity(model, 'C3', Exp(1.0)).on(T2).bound_to(E22)
    C4 = Activity(model, 'C4', Exp(1.0)).on(T2)
    C5 = Activity(model, 'C5', Exp(1.0)).on(T2).replies_to(E22)
    C6 = Activity(model, 'C6', Exp(1.0)).on(T2).bound_to(E23).replies_to(E23)

    # Cache activities
    D0 = Activity(model, 'D0', Immediate()).on(T3).bound_to(E3)
    D1a = Activity(model, 'D1a', Exp(1.0)).on(T3).replies_to(E3)
    D1b = Activity(model, 'D1b', Exp(0.5)).on(T3).replies_to(E3)

    # Precedences
    T0.add_precedence(ActivityPrecedence.Serial(A0, A1, A2, A3))
    T1.add_precedence(ActivityPrecedence.Serial(B0, B1))
    T1.add_precedence(ActivityPrecedence.Serial(B2, B3))
    T1.add_precedence(ActivityPrecedence.Serial(B5, B6, B7))
    T1.add_precedence(ActivityPrecedence.OrFork(B7, [B7a, B7b], [0.7, 0.3]))
    T2.add_precedence(ActivityPrecedence.Serial(C0, C1))
    T2.add_precedence(ActivityPrecedence.Serial(C3, C4, C5))
    T3.add_precedence(ActivityPrecedence.CacheAccess(D0, [D1a, D1b]))

    return model


def gallery_fj_open():
    """Open fork-join network (single class, two parallel tasks)."""
    model = Network('Fork-Join-Open')
    source = Source(model, 'Source')
    queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
    fork = Fork(model, 'Fork')
    join = Join(model, 'Join', fork)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'class1')
    source.setArrival(oclass, Exp(0.05))
    queue1.setService(oclass, Exp(1.0))
    queue2.setService(oclass, Exp(2.0))
    P = model.init_routing_matrix()
    P.set(oclass, oclass, source, fork, 1.0)
    P.set(oclass, oclass, fork, queue1, 1.0)
    P.set(oclass, oclass, fork, queue2, 1.0)
    P.set(oclass, oclass, queue1, join, 1.0)
    P.set(oclass, oclass, queue2, join, 1.0)
    P.set(oclass, oclass, join, sink, 1.0)
    model.link(P)
    return model


def gallery_fj_closed():
    """Closed fork-join network (single class, two parallel tasks)."""
    model = Network('Fork-Join-Closed')
    delay = Delay(model, 'Delay')
    queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.PS)
    fork = Fork(model, 'Fork')
    join = Join(model, 'Join', fork)
    oclass = ClosedClass(model, 'class1', 5, delay)
    delay.setService(oclass, Exp(1.0))
    queue1.setService(oclass, Exp(1.0))
    queue2.setService(oclass, Exp(1.0))
    P = model.init_routing_matrix()
    P.set(oclass, oclass, delay, fork, 1.0)
    P.set(oclass, oclass, fork, queue1, 1.0)
    P.set(oclass, oclass, fork, queue2, 1.0)
    P.set(oclass, oclass, queue1, join, 1.0)
    P.set(oclass, oclass, queue2, join, 1.0)
    P.set(oclass, oclass, join, delay, 1.0)
    model.link(P)
    return model


def gallery_cache_lru():
    """Closed cache model with LRU replacement (n=5 items, m=2 slots)."""
    model = Network('Cache-LRU')
    n = 5
    m = 2
    delay = Delay(model, 'Delay')
    cacheNode = Cache(model, 'Cache', n, m, ReplacementStrategy.LRU)
    jobClass = ClosedClass(model, 'JobClass', 1, delay, 0)
    hitClass = ClosedClass(model, 'HitClass', 0, delay, 0)
    missClass = ClosedClass(model, 'MissClass', 0, delay, 0)
    delay.setService(jobClass, Exp(1))
    pAccess = DiscreteSampler([1.0 / n] * n)
    cacheNode.setRead(jobClass, pAccess)
    cacheNode.setHitClass(jobClass, hitClass)
    cacheNode.setMissClass(jobClass, missClass)
    P = model.init_routing_matrix()
    P.set(jobClass, jobClass, delay, cacheNode, 1.0)
    P.set(hitClass, jobClass, cacheNode, delay, 1.0)
    P.set(missClass, jobClass, cacheNode, delay, 1.0)
    model.link(P)
    return model


def gallery_cache_routing():
    """Open cache with hit/miss routed to distinct queues."""
    model = Network('Cache-Routing')
    n = 4
    m = 2
    source = Source(model, 'Source')
    cacheNode = Cache(model, 'Cache', n, m, ReplacementStrategy.LRU)
    hitQueue = Queue(model, 'HitQueue', SchedStrategy.FCFS)
    missQueue = Queue(model, 'MissQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    jobClass = OpenClass(model, 'InitClass', 0)
    hitClass = OpenClass(model, 'HitClass', 0)
    missClass = OpenClass(model, 'MissClass', 0)
    source.setArrival(jobClass, Exp(1))
    hitQueue.setService(hitClass, Exp(2.0))
    missQueue.setService(missClass, Exp(1.0))
    pAccess = DiscreteSampler([1.0 / n] * n)
    cacheNode.setRead(jobClass, pAccess)
    cacheNode.setHitClass(jobClass, hitClass)
    cacheNode.setMissClass(jobClass, missClass)
    P = model.init_routing_matrix()
    P.set(jobClass, jobClass, source, cacheNode, 1.0)
    P.set(hitClass, hitClass, cacheNode, hitQueue, 1.0)
    P.set(hitClass, hitClass, hitQueue, sink, 1.0)
    P.set(missClass, missClass, cacheNode, missQueue, 1.0)
    P.set(missClass, missClass, missQueue, sink, 1.0)
    model.link(P)
    return model


def gallery_lqn_basic():
    """Basic 3-task layered queueing network (client/server chain)."""
    model = LayeredNetwork('LQN-Basic')
    P1 = Processor(model, 'P1', 2, SchedStrategy.PS)
    P2 = Processor(model, 'P2', 3, SchedStrategy.PS)
    T1 = Task(model, 'T1', 50, SchedStrategy.REF).on(P1).set_think_time(Exp(1 / 2))
    T2 = Task(model, 'T2', 50, SchedStrategy.FCFS).on(P1).set_think_time(Exp(1 / 3))
    T3 = Task(model, 'T3', 25, SchedStrategy.FCFS).on(P2).set_think_time(Exp(1 / 4))
    E1 = Entry(model, 'E1').on(T1)
    E2 = Entry(model, 'E2').on(T2)
    E3 = Entry(model, 'E3').on(T3)
    A1 = Activity(model, 'AS1', Exp(10)).on(T1).bound_to(E1).synch_call(E2, 1)
    A2 = Activity(model, 'AS2', Exp(20)).on(T2).bound_to(E2).synch_call(E3, 5).replies_to(E2)
    A3 = Activity(model, 'AS3', Exp(50)).on(T3).bound_to(E3).replies_to(E3)
    return model


def gallery_lqn_workflows():
    """Layered network with loop, and-fork/join and or-fork/join precedence."""
    model = LayeredNetwork('LQN-Workflows')
    P1 = Processor(model, 'P1', float('inf'), SchedStrategy.INF)
    T1 = Task(model, 'T1', 1, SchedStrategy.REF).on(P1)
    T1.set_think_time(Immediate())
    E1 = Entry(model, 'Entry').on(T1)

    P2 = Processor(model, 'P2', float('inf'), SchedStrategy.INF)
    T2 = Task(model, 'T2', float('inf'), SchedStrategy.INF).on(P2).set_think_time(Immediate())
    E2 = Entry(model, 'E2').on(T2)

    P3 = Processor(model, 'P3', 5, SchedStrategy.PS)
    T3 = Task(model, 'T3', float('inf'), SchedStrategy.INF).on(P3)
    T3.set_think_time(Exp.fit_mean(10))
    E3 = Entry(model, 'E3').on(T3)

    A1 = Activity(model, 'A1', Exp.fit_mean(1)).on(T1).bound_to(E1)
    A2 = Activity(model, 'A2', Exp.fit_mean(2)).on(T1)
    A3 = Activity(model, 'A3', Exp.fit_mean(3)).on(T1).synch_call(E2)

    B1 = Activity(model, 'B1', Exp.fit_mean(0.1)).on(T2).bound_to(E2)
    B2 = Activity(model, 'B2', Exp.fit_mean(0.2)).on(T2)
    B3 = Activity(model, 'B3', Exp.fit_mean(0.3)).on(T2)
    B4 = Activity(model, 'B4', Exp.fit_mean(0.4)).on(T2)
    B5 = Activity(model, 'B5', Exp.fit_mean(0.5)).on(T2)
    B6 = Activity(model, 'B6', Exp.fit_mean(0.6)).on(T2).synch_call(E3).replies_to(E2)

    C1 = Activity(model, 'C1', Exp.fit_mean(0.1)).on(T3).bound_to(E3)
    C2 = Activity(model, 'C2', Exp.fit_mean(0.2)).on(T3)
    C3 = Activity(model, 'C3', Exp.fit_mean(0.3)).on(T3)
    C4 = Activity(model, 'C4', Exp.fit_mean(0.4)).on(T3)
    C5 = Activity(model, 'C5', Exp.fit_mean(0.5)).on(T3).replies_to(E3)

    T1.add_precedence(ActivityPrecedence.Loop(A1, [A2, A3], 3))
    T2.add_precedence(ActivityPrecedence.Serial(B4, B5))
    T2.add_precedence(ActivityPrecedence.AndFork(B1, [B2, B3, B4]))
    T2.add_precedence(ActivityPrecedence.AndJoin([B2, B3, B5], B6))
    T3.add_precedence(ActivityPrecedence.OrFork(C1, [C2, C3, C4], [0.3, 0.3, 0.4]))
    T3.add_precedence(ActivityPrecedence.OrJoin([C2, C3, C4], C5))
    return model


def gallery_mm1k(K=3):
    """M/M/1/K queue with finite capacity K (blocking / loss)."""
    model = Network('M/M/1/K')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    queue.setNumberOfServers(1)
    queue.setCapacity(K)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'Class1', 0)
    source.setArrival(oclass, Exp(0.8))
    queue.setService(oclass, Exp(1.0))
    P = model.init_routing_matrix()
    P.set(oclass, oclass, source, queue, 1.0)
    P.set(oclass, oclass, queue, sink, 1.0)
    model.link(P)
    return model


def gallery_fcr(K=3):
    """Finite capacity region with dropping around a single queue (like M/M/1/K)."""
    model = Network('FCR-Dropping')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'Class1', 0)
    source.setArrival(oclass, Exp(0.8))
    queue.setService(oclass, Exp(1.0))
    P = model.init_routing_matrix()
    P.set(oclass, oclass, source, queue, 1.0)
    P.set(oclass, oclass, queue, sink, 1.0)
    model.link(P)
    fcr = model.add_region([queue])
    fcr.setGlobalMaxJobs(K)
    fcr.setDropRule(oclass, True)
    return model


def gallery_renv_breakdown():
    """Random environment: single server with breakdown/repair (UP/DOWN stages).

    Returns an Environment whose base model is an M/M/1 queue alternating
    between an UP stage (fast service) and a DOWN stage (degraded service).
    """
    model = Network('ServerWithFailures')
    source = Source(model, 'Arrivals')
    queue = Queue(model, 'Server', SchedStrategy.FCFS)
    sink = Sink(model, 'Departures')
    jobclass = OpenClass(model, 'Jobs')
    source.setArrival(jobclass, Exp(0.8))
    queue.setService(jobclass, Exp(2.0))
    queue.setNumberOfServers(1)
    P = model.init_routing_matrix()
    P.set(jobclass, jobclass, source, queue, 1.0)
    P.set(jobclass, jobclass, queue, sink, 1.0)
    model.link(P)
    env = Environment('ServerEnv')
    env.add_node_failure_repair(model, queue, Exp(0.1), Exp(1.0), Exp(0.5))
    env.init()
    return env


def gallery_qn_random(seed=23000):
    """Randomly generated mixed queueing network (reproducible).

    Uses NetworkGenerator with a deterministic (cyclic) topology and a fixed
    default seed so repeated calls yield the same model. Closed network:
    3 queues, 1 delay, 2 closed classes (always stable / solvable).
    """
    random.seed(seed)
    np.random.seed(seed)
    gen = NetworkGenerator(
        sched_strat='fcfs', routing_strat='Probabilities', distribution='Exp',
        cclass_job_load='medium', has_varying_service_rates=False,
        has_multi_server_queues=False, has_random_cs_nodes=False,
        has_multi_chain_cs=False, topology_fcn=cyclic_graph)
    return gen.generate(3, 1, 0, 2)


def gallery_lqn_random(seed=23000):
    """Randomly generated layered queueing network (reproducible).

    Uses LayeredNetworkGenerator with a fixed default seed so repeated calls
    yield the same model. 1 client, 2 levels, 4 tasks, 2 processors.
    """
    random.seed(seed)
    np.random.seed(seed)
    gen = LayeredNetworkGenerator()
    return gen.generate(1, 2, 4, 2)
