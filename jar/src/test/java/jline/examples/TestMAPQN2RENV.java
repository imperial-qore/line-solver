package jline.examples;

import jline.io.MAPQN2RENV;
import jline.lang.*;
import jline.lang.nodes.*;
import jline.lang.processes.*;
import jline.lang.constant.SchedStrategy;

/**
 * Test for MAPQN2RENV conversion.
 */
public class TestMAPQN2RENV {
    public static void main(String[] args) {
        // Create closed queueing network
        Network model = new Network("ClosedQN");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);

        // Create closed job class
        int N = 3;
        ClosedClass jobClass = new ClosedClass(model, "Class1", N, queue1, 0);

        // Set MMPP2 service (Markov-modulated service)
        // Parameters: slow_rate, fast_rate, slow_to_fast, fast_to_slow
        queue1.setService(jobClass, new MMPP2(2.0, 5.0, 0.5, 0.3));
        queue2.setService(jobClass, new Exp(3.0));

        // Create routing
        model.link(Network.serialRouting(queue1, queue2));

        System.out.println("Original network: " + model.getName());
        System.out.println("Queue1 service: MMPP2(2.0, 5.0, 0.5, 0.3)");
        System.out.println("Queue2 service: Exp(3.0)");
        System.out.println();

        // Convert MMPP2 network to random environment
        Environment env = MAPQN2RENV.mapqn2renv(model);
        env.init();
        env.printStageTable();
    }
}
