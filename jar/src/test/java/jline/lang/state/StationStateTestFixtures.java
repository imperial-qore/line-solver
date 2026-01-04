package jline.lang.state;

import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Exp;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Source;
import jline.lang.nodes.Sink;
import jline.lang.ClosedClass;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;

/**
 * Test fixtures for NetworkStateTest.
 * Provides methods to create various network configurations for testing.
 */
public class StationStateTestFixtures {
    
    /**
     * Creates a closed network model with a single class and 3 jobs.
     * Network structure: Delay -> Queue -> Delay
     */
    public static Network createClosedNetwork(SchedStrategy schedStrategy) {
        // Create network model
        Network model = new Network("Closed Network with " + schedStrategy);
        
        // Create nodes
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", schedStrategy);
        queue.setNumberOfServers(1);
        
        // Create closed class with 3 jobs
        ClosedClass jobClass = new ClosedClass(model, "Class1", 3, delay, 0);
        
        // Set service times to match MATLAB exactly
        delay.setService(jobClass, new Exp(1.0));  // Rate = 1.0
        
        // For DPS/GPS, set weights; for HOL, set priority in class definition
        if (schedStrategy == SchedStrategy.DPS || schedStrategy == SchedStrategy.GPS) {
            queue.setService(jobClass, new Exp(0.5), 1.0);  // Rate = 0.5, weight = 1.0
        } else {
            queue.setService(jobClass, new Exp(0.5));  // Rate = 0.5
        }
        
        // Set routing - simple loop: Delay -> Queue -> Delay
        model.link(Network.serialRouting(delay, queue));
        
        return model;
    }
    
    /**
     * Creates a closed network model with two classes, each with 3 jobs.
     * Network structure: Delay -> Queue -> Delay
     */
    public static Network createTwoClassClosedNetwork(SchedStrategy schedStrategy) {
        // Create network model
        Network model = new Network("Two Class Closed Network with " + schedStrategy);
        
        // Create nodes
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", schedStrategy);
        queue.setNumberOfServers(1);
        
        // Create two closed classes with 3 jobs each
        // For HOL, assign different priorities (0 = highest priority)
        ClosedClass jobClass1;
        ClosedClass jobClass2;
        if (schedStrategy == SchedStrategy.HOL) {
            jobClass1 = new ClosedClass(model, "Class1", 3, delay, 0);  // Higher priority
            jobClass2 = new ClosedClass(model, "Class2", 3, delay, 1);  // Lower priority
        } else {
            jobClass1 = new ClosedClass(model, "Class1", 3, delay, 0);
            jobClass2 = new ClosedClass(model, "Class2", 3, delay, 0);
        }
        
        // Set service times for both classes
        delay.setService(jobClass1, Exp.fitMean(1.0));  // Mean service time = 1.0
        delay.setService(jobClass2, Exp.fitMean(1.5));  // Mean service time = 1.5
        
        // For DPS/GPS, set different weights for the classes
        if (schedStrategy == SchedStrategy.DPS || schedStrategy == SchedStrategy.GPS) {
            queue.setService(jobClass1, Exp.fitMean(2.0), 2.0);  // Mean service time = 2.0, weight = 2.0
            queue.setService(jobClass2, Exp.fitMean(2.5), 1.0);  // Mean service time = 2.5, weight = 1.0
        } else {
            queue.setService(jobClass1, Exp.fitMean(2.0));  // Mean service time = 2.0
            queue.setService(jobClass2, Exp.fitMean(2.5));  // Mean service time = 2.5
        }
        
        // Set routing - for multiclass models, we need to use RoutingMatrix
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass1, delay, queue, 1.0);
        P.set(jobClass1, queue, delay, 1.0);
        P.set(jobClass2, delay, queue, 1.0);
        P.set(jobClass2, queue, delay, 1.0);
        model.link(P);
        
        return model;
    }
    
    /**
     * Creates a mixed network model with one closed class (3 jobs) and one open class.
     * Open class has exponential arrivals with rate 0.1 (mean inter-arrival time = 10)
     * Network structure: Source -> Queue -> Sink (for open class)
     *                    Delay -> Queue -> Delay (for closed class)
     */
    public static Network createMixedNetwork(SchedStrategy schedStrategy) {
        // Create network model
        Network model = new Network("Mixed Network with " + schedStrategy);
        
        // Create nodes
        Source source = new Source(model, "Source");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", schedStrategy);
        queue.setNumberOfServers(1);
        Sink sink = new Sink(model, "Sink");
        
        // Create job classes
        // For HOL, assign different priorities (0 = highest priority)
        OpenClass openClass;
        ClosedClass closedClass;
        if (schedStrategy == SchedStrategy.HOL) {
            openClass = new OpenClass(model, "OpenClass", 0);  // Higher priority
            closedClass = new ClosedClass(model, "ClosedClass", 3, delay, 1);  // Lower priority
        } else {
            openClass = new OpenClass(model, "OpenClass", 0);
            closedClass = new ClosedClass(model, "ClosedClass", 3, delay, 0);
        }
        
        // Set arrival process for open class
        source.setArrival(openClass, Exp.fitMean(10.0)); // λ = 0.1
        
        // Set service times
        delay.setService(closedClass, Exp.fitMean(1.0));  // Mean service time = 1.0
        
        // For DPS/GPS, set different weights for the classes
        if (schedStrategy == SchedStrategy.DPS || schedStrategy == SchedStrategy.GPS) {
            queue.setService(openClass, Exp.fitMean(2.0), 1.5);    // Mean service time = 2.0, weight = 1.5
            queue.setService(closedClass, Exp.fitMean(2.5), 1.0);  // Mean service time = 2.5, weight = 1.0
        } else {
            queue.setService(openClass, Exp.fitMean(2.0));    // Mean service time = 2.0
            queue.setService(closedClass, Exp.fitMean(2.5));  // Mean service time = 2.5
        }
        
        // Set routing - for multiclass models, we need to use RoutingMatrix
        RoutingMatrix P = model.initRoutingMatrix();
        // Open class: Source -> Queue -> Sink
        P.set(openClass, source, queue, 1.0);
        P.set(openClass, queue, sink, 1.0);
        // Closed class: Delay -> Queue -> Delay
        P.set(closedClass, delay, queue, 1.0);
        P.set(closedClass, queue, delay, 1.0);
        model.link(P);
        
        return model;
    }
    
    /**
     * Creates an open network model with one open class.
     * Open class has exponential arrivals with rate 0.1 (mean inter-arrival time = 10)
     * Network structure: Source -> Queue -> Sink
     */
    public static Network createOpenNetwork(SchedStrategy schedStrategy) {
        // Create network model
        Network model = new Network("Open Network with " + schedStrategy);
        
        // Create nodes
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", schedStrategy);
        queue.setNumberOfServers(1);
        Sink sink = new Sink(model, "Sink");
        
        // Create open class
        OpenClass openClass = new OpenClass(model, "OpenClass", 0);
        
        // Set arrival process
        source.setArrival(openClass, Exp.fitMean(10.0)); // λ = 0.1
        
        // Set service times
        queue.setService(openClass, Exp.fitMean(2.0));  // Mean service time = 2.0
        
        // Set routing: Source -> Queue -> Sink
        model.link(Network.serialRouting(source, queue, sink));
        
        return model;
    }
}