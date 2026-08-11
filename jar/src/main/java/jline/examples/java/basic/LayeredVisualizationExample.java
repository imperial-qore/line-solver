/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.lang.layered.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Exp;

/**
 * Demonstrates JUNG visualization of LayeredNetwork models using the plot() method.
 *
 * <p>This example creates layered queueing network models and displays them
 * using the built-in plot() method. The visualization showcases the structure
 * of complex multi-tier software architectures.
 *
 * <p>The visualization uses shape coding:
 * <ul>
 * <li>Red pyramids: Processors/Hosts</li>
 * <li>Red parallelograms: Tasks</li>
 * <li>Red rectangles: Entries</li>
 * <li>Red circles: Activities</li>
 * </ul>
 *
 * <p>Edge styles indicate relationship types (all edges are black):
 * <ul>
 * <li>Solid lines: Parent-child containment, entry-to-activity binding, activity precedence</li>
 * <li>Dashed lines: Synchronous calls</li>
 * <li>Dot-dashed lines: Asynchronous calls</li>
 * <li>Dotted lines: Forwarding calls</li>
 * </ul>
 *
 * <p>The window includes menus for switching layouts and mouse modes.
 */
public class LayeredVisualizationExample {

    /**
     * Creates a simple model demonstrating forwarding calls.
     *
     * <p>In this model:
     * <ul>
     * <li>Client calls Server1</li>
     * <li>Server1 forwards the reply to Server2 (50% probability)</li>
     * <li>Server2 sends the reply back to Client</li>
     * </ul>
     *
     * @return LayeredNetwork with forwarding calls
     */
    public static LayeredNetwork createForwardingExample() {
        LayeredNetwork model = new LayeredNetwork("ForwardingExample");

        // Client processor and task
        Processor pClient = new Processor(model, "PClient", 1, SchedStrategy.INF);
        Task tClient = new Task(model, "TClient", 1, SchedStrategy.REF).on(pClient);
        Entry eClient = new Entry(model, "EClient").on(tClient);

        // Server1 processor and task - forwards to Server2
        Processor pServer1 = new Processor(model, "PServer1", 1, SchedStrategy.PS);
        Task tServer1 = new Task(model, "TServer1", 1, SchedStrategy.FCFS).on(pServer1);
        Entry eServer1 = new Entry(model, "EServer1").on(tServer1);

        // Server2 processor and task - receives forwarded requests
        Processor pServer2 = new Processor(model, "PServer2", 1, SchedStrategy.PS);
        Task tServer2 = new Task(model, "TServer2", 1, SchedStrategy.FCFS).on(pServer2);
        Entry eServer2 = new Entry(model, "EServer2").on(tServer2);

        // Client activity calls Server1
        Activity aClient = new Activity(model, "AClient", new Exp(1.0)).on(tClient);
        aClient.boundTo(eClient);
        aClient.synchCall(eServer1, 1);

        // Server1 activity
        Activity aServer1 = new Activity(model, "AServer1", new Exp(2.0)).on(tServer1);
        aServer1.boundTo(eServer1);
        aServer1.repliesTo(eServer1);

        // Server2 activity
        Activity aServer2 = new Activity(model, "AServer2", new Exp(3.0)).on(tServer2);
        aServer2.boundTo(eServer2);
        aServer2.repliesTo(eServer2);

        // Server1 forwards to Server2 with 50% probability
        eServer1.forward(eServer2, 0.5);

        return model;
    }

    /**
     * Main method demonstrating the plot() method on LayeredNetwork.
     *
     * @param args command line arguments (not used)
     * @throws Exception if model creation fails
     */
    public static void main(String[] args) throws Exception {
        System.out.println("=== LayeredNetwork Visualization Examples ===\n");

        // // Show serial model (contains synchronous calls - shown as dashed lines)
        // System.out.println("Visualizing lqn_serial model (synchronous calls - dashed lines)...");
        // LayeredNetwork serialModel = LayeredModel.lqn_serial();
        // serialModel.plot("LQN Serial - Synchronous Calls (dashed)");
        // Thread.sleep(1000);

        // // Show async cache model (contains asynchronous calls - shown as dot-dashed lines)
        // System.out.println("Visualizing lcq_async_prefetch model (asynchronous calls - dot-dashed lines)...");
        // LayeredNetwork asyncModel = CacheModel.lcq_async_prefetch();
        // asyncModel.plot("Async Cache - Asynchronous Calls (dot-dashed)");
        // Thread.sleep(1000);

        // // Show forwarding example (contains forwarding calls - shown as dotted lines)
        // System.out.println("Visualizing forwarding example (forwarding calls - dotted lines)...");
        // LayeredNetwork forwardingModel = createForwardingExample();
        // forwardingModel.plot("Forwarding Example - Forwarding Calls (dotted)");
        // Thread.sleep(1000);

        // // Show basic 3-tier model (small, ~8 elements)
        // System.out.println("Visualizing lqn_basic model...");
        // LayeredNetwork basicModel = LayeredModel.lqn_basic();
        // basicModel.plot("LQN Basic - 3-Tier Architecture");
        // Thread.sleep(1000);

        // Show OFBiz enterprise-scale model (largest)
        System.out.println("Visualizing OFBiz enterprise model (large-scale)...");
        LayeredNetwork ofbizModel = LayeredModel.lqn_ofbiz();
        if (ofbizModel != null) {
            ofbizModel.plot("OFBiz Enterprise Model - Large Scale", 1200, 800);
        } else {
            System.out.println("  (OFBiz model not found in resources)");
        }

        // // Show BPMN workflow model (complex fork-join patterns)
        // System.out.println("Visualizing BPMN workflow model (complex fork-join)...");
        // LayeredNetwork bpmnModel = LayeredModel.lqn_bpmn();
        // bpmnModel.plot("BPMN Workflow - Fork/Join Patterns", 1100, 700);
        // Thread.sleep(1000);

        System.out.println("\nVisualization windows opened. Use Layout menu to switch layouts.");
        System.out.println("Edge styles: Solid=containment/binding/precedence, Dashed=sync, Dot-dashed=async, Dotted=forwarding");
        System.out.println("Close windows to exit.");
    }
}
