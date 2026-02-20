/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.advanced;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.PollingType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Det;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;

/**
 * Examples of models with polling
 */
public class CyclicPollingModel {

    /**
     * M[2]/M[2]/1-Gated polling system example.
     * 
     * This example demonstrates a polling system with:
     * - 2 job classes
     * - Exponential arrivals and service times
     * - Exhaustive polling strategy
     *
     * @return network model with polling queue
     */
    public static Network polling_exhaustive_exp() {
        Network model = new Network("M[2]/M[2]/1-Gated");
        
        // Block 1: nodes
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.POLLING);
        Sink sink = new Sink(model, "mySink");
        
        // Block 2: classes
        OpenClass oclass1 = new OpenClass(model, "myClass1");
        source.setArrival(oclass1, new Exp(0.1));
        queue.setService(oclass1, new Exp(1.0));
        
        OpenClass oclass2 = new OpenClass(model, "myClass2");
        source.setArrival(oclass2, new Exp(0.1));
        queue.setService(oclass2, new Exp(1.5));
        
        queue.setPollingType(PollingType.EXHAUSTIVE);
        
        // Block 3: topology
        model.addLink(source, queue);
        model.addLink(queue, sink);
        
        source.setProbRouting(oclass1, queue, 1.0);
        queue.setProbRouting(oclass1, sink, 1.0);
        
        source.setProbRouting(oclass2, queue, 1.0);
        queue.setProbRouting(oclass2, sink, 1.0);
        
        return model;
    }

    /**
     * M[2]/M[2]/1-Gated polling system with switchover times.
     * 
     * This example demonstrates a polling system with:
     * - 2 job classes
     * - Exponential arrivals and service times
     * - Gated polling strategy
     * - Exponential switchover times between classes
     *
     * @return network model with gated polling and switchover times
     */
    public static Network polling_gated() {
        Network model = new Network("M[2]/M[2]/1-Exhaustive");
        
        // Block 1: nodes
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.POLLING);
        Sink sink = new Sink(model, "mySink");
        
        // Block 2: classes
        OpenClass oclass1 = new OpenClass(model, "myClass1");
        source.setArrival(oclass1, new Exp(1.0));
        queue.setService(oclass1, new Exp(4.0));
        
        OpenClass oclass2 = new OpenClass(model, "myClass2");
        source.setArrival(oclass2, new Exp(0.8));
        queue.setService(oclass2, new Exp(1.5));
        
        queue.setPollingType(PollingType.GATED);
        queue.setSwitchover(oclass1, new Exp(1.0));
        queue.setSwitchover(oclass2, new Exp(0.5));
        
        // Block 3: topology
        model.addLink(source, queue);
        model.addLink(queue, sink);
        
        source.setProbRouting(oclass1, queue, 1.0);
        queue.setProbRouting(oclass1, sink, 1.0);
        
        source.setProbRouting(oclass2, queue, 1.0);
        queue.setProbRouting(oclass2, sink, 1.0);
        
        return model;
    }

    /**
     * M[2]/M[2]/1-K-Limited polling system with mixed switchover.
     * 
     * This example demonstrates a polling system with:
     * - 2 job classes
     * - Exponential arrivals and service times
     * - K-Limited polling strategy with K=1
     * - Mixed switchover times (Exponential for class 1, Immediate for class 2)
     *
     * The setPollingType method now accepts an optional K parameter for K-LIMITED polling.
     *
     * @return network model with K-limited polling and mixed switchover
     */
    public static Network polling_klimited() {
        Network model = new Network("M[2]/M[2]/1-K-Limited");
        
        // Block 1: nodes
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.POLLING);
        Sink sink = new Sink(model, "mySink");
        
        // Block 2: classes
        OpenClass oclass1 = new OpenClass(model, "myClass1");
        source.setArrival(oclass1, new Exp(0.2));
        queue.setService(oclass1, new Exp(1.0));
        
        OpenClass oclass2 = new OpenClass(model, "myClass2");
        source.setArrival(oclass2, new Exp(0.3));
        queue.setService(oclass2, new Exp(1.5));
        
        queue.setPollingType(PollingType.KLIMITED, 1); // Set K-LIMITED polling with K=1
        queue.setSwitchover(oclass1, new Exp(1));
        // Use Exp(1) for switchover (same as oclass1)
        queue.setSwitchover(oclass2, new Exp(1));
        
        // Block 3: topology
        model.addLink(source, queue);
        model.addLink(queue, sink);
        
        source.setProbRouting(oclass1, queue, 1.0);
        queue.setProbRouting(oclass1, sink, 1.0);
        
        source.setProbRouting(oclass2, queue, 1.0);
        queue.setProbRouting(oclass2, sink, 1.0);
        
        return model;
    }

    /**
     * M[2]/M[2]/1-Exhaustive polling system with deterministic arrivals and immediate switchover.
     * 
     * This example demonstrates a polling system with:
     * - 2 job classes
     * - Deterministic arrivals and service times
     * - Exhaustive polling strategy
     * - Immediate (zero-time) switchover between classes
     *
     * @return network model with exhaustive polling and immediate switchover
     */
    public static Network polling_exhaustive_det() {
        Network model = new Network("M[2]/M[2]/1-Gated");
        
        // Block 1: nodes
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.POLLING);
        Sink sink = new Sink(model, "mySink");
        
        // Block 2: classes
        OpenClass oclass1 = new OpenClass(model, "myClass1");
        source.setArrival(oclass1, Det.fitMean(1.0));
        queue.setService(oclass1, Det.fitMean(0.001));
        
        OpenClass oclass2 = new OpenClass(model, "myClass2");
        source.setArrival(oclass2, Det.fitMean(1.0));
        queue.setService(oclass2, Det.fitMean(0.001));
        
        queue.setPollingType(PollingType.EXHAUSTIVE);
        // Use Exp(1e6) for very fast switchover (approximating immediate)
        queue.setSwitchover(oclass2, new Exp(1e6));
        queue.setSwitchover(oclass1, new Exp(1e6));
        
        // Block 3: topology
        model.addLink(source, queue);
        model.addLink(queue, sink);
        
        source.setProbRouting(oclass1, queue, 1.0);
        queue.setProbRouting(oclass1, sink, 1.0);
        
        source.setProbRouting(oclass2, queue, 1.0);
        queue.setProbRouting(oclass2, sink, 1.0);
        
        return model;
    }

    /**
     * Main method for testing and demonstrating cyclic polling examples.
     *
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        System.out.println("Testing Cyclic Polling Examples:");
        
        // Test Example 1: EXHAUSTIVE polling
        System.out.println("\nExample 1: EXHAUSTIVE polling");
        Network model1 = polling_exhaustive_exp();
        System.out.println("Model created successfully: " + model1.getName());
        System.out.println("Number of nodes: " + model1.getNodes().size());
        System.out.println("Number of classes: " + model1.getClasses().size());
        
        // Test Example 2: GATED polling with switchover times
        System.out.println("\nExample 2: GATED polling with switchover times");
        Network model2 = polling_gated();
        System.out.println("Model created successfully: " + model2.getName());
        System.out.println("Number of nodes: " + model2.getNodes().size());
        System.out.println("Number of classes: " + model2.getClasses().size());
        
        // Test Example 3: K-LIMITED polling with mixed switchover
        System.out.println("\nExample 3: K-LIMITED polling with mixed switchover");
        Network model3 = polling_klimited();
        System.out.println("Model created successfully: " + model3.getName());
        System.out.println("Number of nodes: " + model3.getNodes().size());
        System.out.println("Number of classes: " + model3.getClasses().size());
        
        // Test Example 4: EXHAUSTIVE polling with deterministic arrivals
        System.out.println("\nExample 4: EXHAUSTIVE polling with deterministic arrivals");
        Network model4 = polling_exhaustive_det();
        System.out.println("Model created successfully: " + model4.getName());
        System.out.println("Number of nodes: " + model4.getNodes().size());
        System.out.println("Number of classes: " + model4.getClasses().size());
        
        System.out.println("\nAll cyclic polling examples completed successfully!");
    }
}
