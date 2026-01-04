package jline.examples.java.advanced;

import jline.lang.Network;
import jline.lang.ClosedClass;
import jline.lang.OpenClass;
import jline.lang.Signal;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SignalType;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;

/**
 * Model factory for MAM (RCAT/INAP method) examples.
 *
 * These models demonstrate the capabilities of the RCAT algorithm for
 * analyzing queueing networks through agent decomposition.
 */
public class AgentModel {

    /**
     * Creates an open tandem queue model (M/M/1 -> M/M/1).
     *
     * Network structure: Source -> Queue1 -> Queue2 -> Sink
     *
     * Parameters:
     * - Arrival rate: 0.5
     * - Service rate at Queue1: 1.0
     * - Service rate at Queue2: 1.5
     *
     * @return the network model
     */
    public static Network tandemOpen() {
        Network model = new Network("Tandem-MM1");

        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(0.5));
        queue1.setService(oclass, new Exp(1.0));
        queue2.setService(oclass, new Exp(1.5));

        model.link(Network.serialRouting(source, queue1, queue2, sink));

        return model;
    }

    /**
     * Creates a closed network with two PS queues.
     *
     * Network structure: Queue1 <-> Queue2 (closed loop)
     *
     * Parameters:
     * - Number of jobs: 10
     * - Service rate at Queue1: 2.0
     * - Service rate at Queue2: 1.0
     *
     * @return the network model
     */
    public static Network closedNetwork() {
        Network model = new Network("Closed-2Q");

        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);

        ClosedClass cclass = new ClosedClass(model, "Class1", 10, queue1);
        queue1.setService(cclass, new Exp(2.0));
        queue2.setService(cclass, new Exp(1.0));

        model.link(Network.serialRouting(queue1, queue2));

        return model;
    }

    /**
     * Creates a multiclass closed network.
     *
     * Network structure: Queue1 <-> Queue2 (closed loop, 2 classes)
     *
     * Parameters:
     * - Class 1: 5 jobs
     * - Class 2: 3 jobs
     * - Different service rates per class at each queue
     *
     * @return the network model
     */
    public static Network multiclassClosed() {
        Network model = new Network("Multiclass-Closed");

        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);

        ClosedClass class1 = new ClosedClass(model, "Class1", 5, queue1);
        ClosedClass class2 = new ClosedClass(model, "Class2", 3, queue1);

        queue1.setService(class1, new Exp(2.0));
        queue1.setService(class2, new Exp(1.5));
        queue2.setService(class1, new Exp(1.0));
        queue2.setService(class2, new Exp(0.8));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, queue1, queue2, 1.0);
        P.set(class1, class1, queue2, queue1, 1.0);
        P.set(class2, class2, queue1, queue2, 1.0);
        P.set(class2, class2, queue2, queue1, 1.0);
        model.link(P);

        return model;
    }

    /**
     * Creates a Jackson network with probabilistic routing.
     *
     * Network structure: Source -> Queue1/2/3 (with feedback) -> Sink
     *
     * Parameters:
     * - Arrival rate: 1.0
     * - Service rates: [2.0, 3.0, 2.5]
     * - Routing with feedback between queues
     *
     * @return the network model
     */
    public static Network jacksonNetwork() {
        Network model = new Network("Jackson-3Q");

        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(1.0));
        queue1.setService(oclass, new Exp(2.0));
        queue2.setService(oclass, new Exp(3.0));
        queue3.setService(oclass, new Exp(2.5));

        // Routing probabilities
        double p12 = 0.4, p13 = 0.3, p1s = 0.3;
        double p21 = 0.2, p23 = 0.3, p2s = 0.5;
        double p3s = 1.0;

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(oclass, oclass, source, queue1, 1.0);
        P.set(oclass, oclass, queue1, queue2, p12);
        P.set(oclass, oclass, queue1, queue3, p13);
        P.set(oclass, oclass, queue1, sink, p1s);
        P.set(oclass, oclass, queue2, queue1, p21);
        P.set(oclass, oclass, queue2, queue3, p23);
        P.set(oclass, oclass, queue2, sink, p2s);
        P.set(oclass, oclass, queue3, sink, p3s);
        model.link(P);

        return model;
    }

    /**
     * Creates a G-network (Gelenbe network) with negative customers.
     *
     * G-networks extend standard queueing networks with "negative customers"
     * (signals) that remove jobs from queues upon arrival. This models
     * scenarios like job cancellations or service interrupts.
     *
     * Network structure:
     * - Source generates positive customers (jobs) and negative signals
     * - Positive customers: Source -> Queue1 -> Queue2 -> Sink
     * - Negative signals: Source -> Queue1 -> Queue2 (removes job) -> Sink
     *
     * Parameters:
     * - Positive arrival rate: 1.0
     * - Negative signal rate: 0.3
     * - Service rates: [2.0, 3.0]
     *
     * Reference: Gelenbe, E. (1991). "Product-form queueing networks with
     *            negative and positive customers", Journal of Applied Probability
     *
     * @return the network model with signal class configured
     */
    public static Network gNetwork() {
        double lambdaPos = 1.0;   // Positive customer arrival rate
        double lambdaNeg = 0.3;   // Negative signal arrival rate
        double mu1 = 2.0;         // Service rate at Queue1
        double mu2 = 3.0;         // Service rate at Queue2

        Network model = new Network("GNetwork-Example");

        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        // Positive customer class (normal jobs)
        OpenClass posClass = new OpenClass(model, "Positive");
        source.setArrival(posClass, new Exp(lambdaPos));
        queue1.setService(posClass, new Exp(mu1));
        queue2.setService(posClass, new Exp(mu2));

        // Negative signal class (removes jobs from target queue)
        // Using Signal class with SignalType.NEGATIVE for automatic G-network handling
        Signal negClass = new Signal(model, "Negative", SignalType.NEGATIVE);
        source.setArrival(negClass, new Exp(lambdaNeg));
        queue1.setService(negClass, new Exp(mu1));  // Signals also get "served"
        queue2.setService(negClass, new Exp(mu2));

        // Set routing matrix
        RoutingMatrix P = model.initRoutingMatrix();
        // Positive customers: Source -> Queue1 -> Queue2 -> Sink
        P.set(posClass, posClass, source, queue1, 1.0);
        P.set(posClass, posClass, queue1, queue2, 1.0);
        P.set(posClass, posClass, queue2, sink, 1.0);
        // Negative signals: Source -> Queue1 -> Queue2 (removes job) -> Sink
        P.set(negClass, negClass, source, queue1, 1.0);
        P.set(negClass, negClass, queue1, queue2, 1.0);
        P.set(negClass, negClass, queue2, sink, 1.0);
        model.link(P);

        // Signal class is automatically detected by refreshStruct

        return model;
    }
}
