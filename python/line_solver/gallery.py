# Use relative imports to avoid circular import issues
from .lang import (
    Network, Source, Queue, Sink, Delay, OpenClass, ClosedClass,
    Fork, Join, ClassSwitch, Router
)
from .distributions import (
    APH, Cox2, Coxian, Det, Disabled, Erlang, Exp, Gamma, HyperExp,
    Immediate, MAP, Pareto, PH, Replayer, Uniform
)
from .constants import SchedStrategy, RoutingStrategy


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
    model = Network('M/M/1-Priority')
    source1 = Source(model, 'HighPrioSource')
    source2 = Source(model, 'LowPrioSource')
    queue = Queue(model, 'myQueue', SchedStrategy.HOL)
    sink = Sink(model, 'mySink')

    highPrioClass = OpenClass(model, 'HighPrio')
    lowPrioClass = OpenClass(model, 'LowPrio')

    queue.setPriorityClass(highPrioClass, 1)
    queue.setPriorityClass(lowPrioClass, 2)

    source1.setArrival(highPrioClass, Exp(0.3))
    source2.setArrival(lowPrioClass, Exp(0.5))
    queue.setService(highPrioClass, Exp(2))
    queue.setService(lowPrioClass, Exp(2))

    model.link(Network.serial_routing(source1, queue, sink))
    model.link(Network.serial_routing(source2, queue, sink))
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
