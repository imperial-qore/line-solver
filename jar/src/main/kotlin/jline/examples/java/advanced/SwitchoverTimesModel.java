/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.advanced;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.lang.processes.Erlang;
import jline.solvers.jmt.JMT;

/**
 * Examples of models with switchover times
 */
public class SwitchoverTimesModel {

    /**
     * Demonstrates switchover times between job classes in a multi-class queueing system.
     * 
     * This example shows:
     * - Open queueing network with 2 job classes
     * - FCFS scheduling strategy with setup/changeover times
     * - Different arrival and service rates for each class
     * - Switchover times when server switches between serving different classes
     * 
     * Note: The Java implementation uses the available switchover functionality:
     * - Uses single-parameter setSwitchover for POLLING scheduling strategy
     * - Demonstrates both Exponential and Erlang distributions for switchover times
     * - Shows realistic modeling of setup/changeover times in multi-class systems
     * 
     * @return network model demonstrating switchover times concept
     */
    public static Network switchover_basic() {
        Network model = new Network("M[2]/M[2]/1-Gated");
        
        // Block 1: nodes
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.POLLING); // Use POLLING to support switchover
        Sink sink = new Sink(model, "mySink");
        
        // Block 2: classes
        OpenClass oclass1 = new OpenClass(model, "myClass1");
        source.setArrival(oclass1, new Exp(0.2));
        queue.setService(oclass1, new Exp(0.5));  // Service rate must exceed arrival rate for stability
        
        OpenClass oclass2 = new OpenClass(model, "myClass2");
        source.setArrival(oclass2, new Exp(0.8));
        queue.setService(oclass2, new Exp(1.5));
        
        // Set switchover times between job classes
        // In MATLAB: queue.setSwitchover(oclass1, oclass2, Exp(1))
        // In MATLAB: queue.setSwitchover(oclass2, oclass1, Erlang(1,2))
        queue.setSwitchover(oclass1, new Exp(1.0));        // Exponential switchover time for class 1
        queue.setSwitchover(oclass2, new Erlang(1.0, 2));  // Erlang(rate=1.0, phases=2) switchover time for class 2
        
        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        
        routingMatrix.set(oclass1, oclass1, source, queue, 1.0);
        routingMatrix.set(oclass1, oclass1, queue, sink, 1.0);
        
        routingMatrix.set(oclass2, oclass2, source, queue, 1.0);
        routingMatrix.set(oclass2, oclass2, queue, sink, 1.0);
        
        model.link(routingMatrix);
        
        return model;
    }

    /**
     * Main method for testing and demonstrating switchover times examples.
     *
     * @param args command line arguments (not used)
     * @throws Exception if solver encounters an error
     */
    public static void main(String[] args) throws Exception {
        switchover_basic();
    }

}
