package jline.examples;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;

/**
 * Test script to reproduce the environment example from the manual
 */
public class TestEnvPrintStageTable {

    public static void main(String[] args) {
        // Create network1 (Online stage)
        Network network1 = new Network("network1");
        Delay delay1 = new Delay(network1, "Delay");
        Queue queue1 = new Queue(network1, "Queue", SchedStrategy.FCFS);

        OpenClass class1 = new OpenClass(network1, "Class1");
        delay1.setService(class1, new Exp(1.0));
        queue1.setService(class1, new Exp(2.0));

        RoutingMatrix routing1 = network1.initRoutingMatrix();
        routing1.set(delay1, queue1, 1.0);
        routing1.set(queue1, delay1, 1.0);
        network1.link(routing1);

        // Create network2 (Offline stage)
        Network network2 = new Network("network2");
        Delay delay2 = new Delay(network2, "Delay");
        Queue queue2 = new Queue(network2, "Queue", SchedStrategy.FCFS);

        OpenClass class2 = new OpenClass(network2, "Class1");
        delay2.setService(class2, new Exp(0.5));
        queue2.setService(class2, new Exp(1.0));

        RoutingMatrix routing2 = network2.initRoutingMatrix();
        routing2.set(delay2, queue2, 1.0);
        routing2.set(queue2, delay2, 1.0);
        network2.link(routing2);

        // Create environment
        int E = 2;
        Environment envModel = new Environment("UnreliableEnv", E);

        // Add stages
        envModel.addStage(0, "Online", "UP", network1);
        envModel.addStage(1, "Offline", "DOWN", network2);

        // Add transitions
        envModel.addTransition(0, 1, new Exp(1));
        envModel.addTransition(1, 0, new Exp(2));

        // Initialize environment
        envModel.init();

        // Print stage table
        System.out.println("=== Java Output ===");
        envModel.printStageTable();
    }
}
