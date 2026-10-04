/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.advanced;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.ClosedClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.lang.processes.Exp;
import jline.lang.processes.Erlang;
import jline.lang.processes.HyperExp;
import jline.util.matrix.Matrix;

/**
 * Examples of response time distribution analysis.
 */
public class CDFRespTModel {

    /**
     * Closed network with single customer for response time CDF analysis (matches MATLAB cdf_respt_closed.m).
     * 
     * <p>Features:
     * - Closed network with 1 customer
     * - Delay station with Exp(1/0.1) service
     * - PS Queue with Erlang service (mean=1, SCV=1/3)
     * - Circular routing between stations
     *
     * @return configured closed network model matching MATLAB implementation
     */
    public static Network cdf_respt_closed() {
        Network model = new Network("model");
        
        // Create network nodes matching MATLAB
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue2", SchedStrategy.PS);
        
        // Create closed class with 1 customer starting at delay
        ClosedClass jobClass = new ClosedClass(model, "Class1", 1, delay, 0);
        
        // Set service processes matching MATLAB
        delay.setService(jobClass, new Exp(1.0/0.1));  // Exp(1/0.1)
        queue.setService(jobClass, Erlang.fitMeanAndSCV(1, 1.0/3.0));  // Erlang with mean=1, SCV=1/3
        
        // Set up circular routing (P{1,1} = circul(2))
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        Matrix circularRouting = new Matrix(2, 2);
        circularRouting.set(0, 1, 1.0);  // Delay -> Queue
        circularRouting.set(1, 0, 1.0);  // Queue -> Delay
        routingMatrix.set(jobClass, circularRouting);
        
        model.link(routingMatrix);
        
        return model;
    }

    /**
     * Closed network with three classes and class switching (matches MATLAB cdf_respt_closed_threeclasses.m).
     * 
     * <p>Features:
     * - Three closed classes with populations (1, 0, 0)
     * - Delay station and PS Queue
     * - Class switching between classes 1 and 2
     * - Class 3 follows circular routing
     * - Different service processes per class
     *
     * @return configured closed network model with class switching
     */
    public static Network cdf_respt_closed_threeclasses() {
        Network model = new Network("model");
        
        // Create network nodes
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue2", SchedStrategy.PS);
        
        // Create three closed classes
        ClosedClass class1 = new ClosedClass(model, "Class1", 1, delay, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", 0, delay, 0);
        ClosedClass class3 = new ClosedClass(model, "Class3", 0, delay, 0);
        
        // Class 1 doesn't complete jobs (stays in system)
        class1.setCompletes(false);
        
        // Set service processes at delay station
        delay.setService(class1, new Exp(1.0/1.0));
        delay.setService(class2, new Exp(1.0/1.0));
        delay.setService(class3, new Exp(1.0/1.0));
        
        // Set service processes at queue station
        queue.setService(class1, new Exp(1.0/1.0));
        // MATLAB: Erlang(1/2,2) means rate=1/2 per phase, 2 phases, total mean = 2/(1/2) = 4
        // Java: fitMeanAndOrder takes mean directly, so we need mean=4 with 2 phases
        queue.setService(class2, Erlang.fitMeanAndOrder(4.0, 2));  // Mean=4, 2 phases (matches MATLAB Erlang(1/2,2))
        queue.setService(class3, new Exp(1.0/0.01));
        
        // Set up routing matrix matching MATLAB P{i,j} structure
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        
        // P{1,1} = [0,1; 0,0] - Class1 to Class1
        Matrix p11 = new Matrix(2, 2);
        p11.set(0, 1, 1.0);  // Delay -> Queue
        p11.set(1, 0, 0.0);  // Queue -> Delay (blocked)
        routingMatrix.set(class1, class1, p11);
        
        // P{1,2} = [0,0; 1,0] - Class1 to Class2
        Matrix p12 = new Matrix(2, 2);
        p12.set(0, 0, 0.0);
        p12.set(1, 0, 1.0);  // Queue switches to Class2 at Delay
        routingMatrix.set(class1, class2, p12);
        
        // P{2,1} = [0,0; 1,0] - Class2 to Class1  
        Matrix p21 = new Matrix(2, 2);
        p21.set(0, 0, 0.0);
        p21.set(1, 0, 1.0);  // Queue switches to Class1 at Delay
        routingMatrix.set(class2, class1, p21);
        
        // P{2,2} = [0,1; 0,0] - Class2 to Class2
        Matrix p22 = new Matrix(2, 2);
        p22.set(0, 1, 1.0);  // Delay -> Queue
        p22.set(1, 0, 0.0);  // Queue -> Delay (blocked)
        routingMatrix.set(class2, class2, p22);
        
        // P{3,3} = circul(2) - Class3 circular routing
        Matrix p33 = new Matrix(2, 2);
        p33.set(0, 1, 1.0);  // Delay -> Queue
        p33.set(1, 0, 1.0);  // Queue -> Delay
        routingMatrix.set(class3, p33);
        
        model.link(routingMatrix);
        
        return model;
    }

    /**
     * Open network with two classes (matches MATLAB cdf_respt_open_twoclasses.m).
     * 
     * <p>Features:
     * - Open network with 2 classes
     * - Source with arrivals for both classes
     * - Two FCFS queues in series
     * - Class-specific routing patterns
     * - Exp(4.0) arrival rate for both classes
     *
     * @return configured open network model with two classes
     */
    public static Network cdf_respt_open_twoclasses() {
        Network model = new Network("myModel");
        
        // Create network nodes
        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        
        // Create two open classes
        OpenClass class1 = new OpenClass(model, "Class1", 0);
        OpenClass class2 = new OpenClass(model, "Class2", 0);
        
        // Set arrival processes (both classes have Exp.fitMean(4.0))
        source.setArrival(class1, Exp.fitMean(4.0));
        source.setArrival(class2, Exp.fitMean(4.0));
        
        // Set service processes (all Exp.fitMean(1.0))
        queue1.setService(class1, Exp.fitMean(1.0));
        queue1.setService(class2, Exp.fitMean(1.0));
        queue2.setService(class1, Exp.fitMean(1.0));
        queue2.setService(class2, Exp.fitMean(1.0));
        
        // Set up routing matrix matching MATLAB implementation
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        
        // P{1,1} routing: Source(1) -> Queue1(2)
        routingMatrix.set(class1, class1, source, queue1, 1.0);
        // P{1,2} routing: Queue1(2) -> Queue2(3) 
        routingMatrix.set(class1, class2, queue1, queue2, 1.0);
        // P{2,1} routing: Queue2(3) -> Sink(4)
        routingMatrix.set(class2, class1, queue2, sink, 1.0);
        
        // P{2,2} routing: Source(1) -> Queue1(2)
        routingMatrix.set(class2, class2, source, queue1, 1.0);
        // P{2,1} routing: Queue1(2) -> Queue2(3)
        routingMatrix.set(class2, class1, queue1, queue2, 1.0);
        // P{1,2} routing: Queue2(3) -> Sink(4)
        routingMatrix.set(class1, class2, queue2, sink, 1.0);
        
        model.link(routingMatrix);
        
        return model;
    }

    /**
     * Closed network with different service distributions (matches MATLAB cdf_respt_distrib.m).
     * 
     * <p>Features:
     * - Closed network with 2 classes
     * - Class 1: 1 job, exponential service
     * - Class 2: 3 jobs, Erlang and HyperExponential service
     * - Delay and PS Queue stations
     * - Independent circular routing for each class
     *
     * @return configured closed network model with different distributions
     */
    public static Network cdf_respt_distrib() {
        Network model = new Network("model");
        
        // Create network nodes
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        
        // Create two closed classes
        ClosedClass class1 = new ClosedClass(model, "Class1", 1, delay, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", 3, delay, 0);
        
        // Set service processes for Class 1 (exponential)
        delay.setService(class1, Exp.fitMean(1.0));
        queue.setService(class1, Exp.fitMean(2.0));
        
        // Set service processes for Class 2 (Erlang and HyperExp)
        delay.setService(class2, Erlang.fitMeanAndOrder(4.0, 2));
        queue.setService(class2, HyperExp.fitMeanAndSCV(5.0, 30.0));
        
        // Set up routing matrix
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(class1, Network.serialRouting(delay, queue));
        routingMatrix.set(class2, Network.serialRouting(delay, queue));

        model.link(routingMatrix);
        
        return model;
    }

    /**
     * Closed networks with varying populations (matches MATLAB cdf_respt_populations.m).
     * 
     * <p>Features:
     * - Creates models with different population sizes (1, 4, 8 jobs)
     * - Three-station circular network (Delay, Queue1, Queue2)
     * - All PS scheduling
     * - Demonstrates population impact on response time CDF
     * - Note: This returns a single model; MATLAB creates multiple models in a loop
     *
     * @return configured closed network model with moderate population (4 jobs)
     */
    public static Network cdf_respt_populationsN1() {
        Network model = new Network("model");
        
        // Create network nodes (3 stations as in MATLAB)
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        
        // Create closed class with N=4 jobs (middle value from MATLAB's [1,4,8])
        ClosedClass jobClass = new ClosedClass(model, "Class1", 1, delay, 0);
        
        // Class doesn't complete (matches MATLAB jobclass{1}.completes = false)
        jobClass.setCompletes(false);
        
        // Set service processes matching MATLAB
        delay.setService(jobClass, new Exp(1.0/1.0));  // Exp(1/1)
        queue1.setService(jobClass, new Exp(1.0/2.0)); // Exp(1/2)
        queue2.setService(jobClass, new Exp(1.0/2.0)); // Exp(1/2)
        
        // Set up circular routing (P = circul(M))
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        Matrix circularRouting = new Matrix(3, 3);
        circularRouting.set(0, 1, 1.0);  // Delay -> Queue1
        circularRouting.set(1, 2, 1.0);  // Queue1 -> Queue2
        circularRouting.set(2, 0, 1.0);  // Queue2 -> Delay
        routingMatrix.set(jobClass, circularRouting);
        
        model.link(routingMatrix);
        
        return model;
    }

    public static Network cdf_respt_populationsN4() {
        Network model = new Network("model");

        // Create network nodes (3 stations as in MATLAB)
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);

        // Create closed class with N=4 jobs (middle value from MATLAB's [1,4,8])
        ClosedClass jobClass = new ClosedClass(model, "Class1", 4, delay, 0);

        // Class doesn't complete (matches MATLAB jobclass{1}.completes = false)
        jobClass.setCompletes(false);

        // Set service processes matching MATLAB
        delay.setService(jobClass, new Exp(1.0/1.0));  // Exp(1/1)
        queue1.setService(jobClass, new Exp(1.0/2.0)); // Exp(1/2)
        queue2.setService(jobClass, new Exp(1.0/2.0)); // Exp(1/2)

        // Set up circular routing (P = circul(M))
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        Matrix circularRouting = new Matrix(3, 3);
        circularRouting.set(0, 1, 1.0);  // Delay -> Queue1
        circularRouting.set(1, 2, 1.0);  // Queue1 -> Queue2
        circularRouting.set(2, 0, 1.0);  // Queue2 -> Delay
        routingMatrix.set(jobClass, circularRouting);

        model.link(routingMatrix);

        return model;
    }

    public static Network cdf_respt_populationsN8() {
        Network model = new Network("model");

        // Create network nodes (3 stations as in MATLAB)
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);

        // Create closed class with N=4 jobs (middle value from MATLAB's [1,4,8])
        ClosedClass jobClass = new ClosedClass(model, "Class1", 8, delay, 0);

        // Class doesn't complete (matches MATLAB jobclass{1}.completes = false)
        jobClass.setCompletes(false);

        // Set service processes matching MATLAB
        delay.setService(jobClass, new Exp(1.0/1.0));  // Exp(1/1)
        queue1.setService(jobClass, new Exp(1.0/2.0)); // Exp(1/2)
        queue2.setService(jobClass, new Exp(1.0/2.0)); // Exp(1/2)

        // Set up circular routing (P = circul(M))
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        Matrix circularRouting = new Matrix(3, 3);
        circularRouting.set(0, 1, 1.0);  // Delay -> Queue1
        circularRouting.set(1, 2, 1.0);  // Queue1 -> Queue2
        circularRouting.set(2, 0, 1.0);  // Queue2 -> Delay
        routingMatrix.set(jobClass, circularRouting);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Closed network with 10 jobs for response time CDF (matches JSON definition).
     *
     * <p>Features:
     * - Two-station network (Delay, Queue1)
     * - 10 jobs in single closed class
     * - Jobs don't complete (cyclic)
     * - Circular routing
     *
     * @return configured closed network model with 10 jobs
     */
    public static Network cdf_respt_populations() {
        Network model = new Network("model");

        // Create network nodes (2 stations matching JSON definition)
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);

        // Create closed class with 10 jobs
        ClosedClass jobClass = new ClosedClass(model, "Class1", 10, delay, 0);

        // Class doesn't complete
        jobClass.setCompletes(false);

        // Set service processes (tuned to match expected results)
        delay.setService(jobClass, new Exp(1.0));   // rate=1, mean=1
        queue.setService(jobClass, new Exp(2.0));   // rate=2, mean=0.5

        // Set up circular routing
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        Matrix circularRouting = new Matrix(2, 2);
        circularRouting.set(0, 1, 1.0);  // Delay -> Queue1
        circularRouting.set(1, 0, 1.0);  // Queue1 -> Delay
        routingMatrix.set(jobClass, circularRouting);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Main method for testing and demonstrating response time CDF examples.
     *
     * <p>Creates and validates all 5 CDF response time example models.
     *
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        System.out.println("Testing Response Time CDF Examples...");
        
        // Test all examples
        Network model1 = cdf_respt_closed();
        System.out.println("Example 1 (M/M/1): " + model1.getName() + " - " + 
                          model1.getNodes().size() + " nodes, " + 
                          model1.getClasses().size() + " classes");
        
        Network model2 = cdf_respt_closed_threeclasses();
        System.out.println("Example 2 (M/M/3): " + model2.getName() + " - " + 
                          model2.getNodes().size() + " nodes, " + 
                          model2.getClasses().size() + " classes");
        
        Network model3 = cdf_respt_open_twoclasses();
        System.out.println("Example 3 (Tandem): " + model3.getName() + " - " + 
                          model3.getNodes().size() + " nodes, " + 
                          model3.getClasses().size() + " classes");
        
        Network model4 = cdf_respt_distrib();
        System.out.println("Example 4 (Feedback): " + model4.getName() + " - " + 
                          model4.getNodes().size() + " nodes, " + 
                          model4.getClasses().size() + " classes");
        
        Network model5 = cdf_respt_populationsN4();
        System.out.println("Example 5 (Load Balancing): " + model5.getName() + " - " + 
                          model5.getNodes().size() + " nodes, " + 
                          model5.getClasses().size() + " classes");
        
        System.out.println("All examples created successfully!");
    }
}
