package jline.examples.java.advanced;

import jline.lang.Network;
import jline.lang.ClosedClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.DropStrategy;
import jline.lang.processes.Exp;
import jline.lang.nodes.Queue;
import jline.solvers.des.SolverDES;
import jline.solvers.des.DESOptions;
import jline.solvers.NetworkAvgTable;

/**
 * Debug test for BAS blocking in closed queueing networks.
 *
 * This test creates a closed network with 2 queues and 2 jobs, where:
 * - Queue1 has infinite capacity and BAS blocking rule
 * - Queue2 has finite capacity = 1
 *
 * The test runs a small number of samples to trace the BAS blocking behavior.
 */
public class BASBlockingDebugTest {

    public static void main(String[] args) {
        System.out.println("=== BAS Blocking Debug Test ===");
        System.out.println("Closed network: 2 queues, 2 jobs");
        System.out.println("Queue1: infinite cap, BAS blocking");
        System.out.println("Queue2: capacity = 1");
        System.out.println();

        // Create model
        Network model = new Network("des_bas_closed_debug");

        // Create two queues
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);

        // Create closed class with 2 jobs, reference station is Queue1
        ClosedClass jobclass = new ClosedClass(model, "Class1", 2, queue1, 0);

        // Service rates
        queue1.setService(jobclass, new Exp(1.0));  // Mean service time = 1
        queue2.setService(jobclass, new Exp(0.8));  // Mean service time = 1.25

        // Set number of servers
        queue1.setNumberOfServers(1);
        queue2.setNumberOfServers(1);

        // Set finite capacity for Queue2
        queue2.setCap(1);

        // Set BAS blocking rule at Queue1
        queue1.setDropRule(jobclass, DropStrategy.BlockingAfterService);

        // Serial routing: Queue1 -> Queue2 -> Queue1 (cyclic)
        model.link(Network.serialRouting(queue1, queue2));

        // Configure DES options
        DESOptions options = new DESOptions();
        options.seed = 23000;
        options.samples = 1000000;  // Large sample for accuracy

        System.out.println("Running DES with " + options.samples + " samples...\n");

        // Run solver
        SolverDES solver = new SolverDES(model, options);
        NetworkAvgTable avgTable = (NetworkAvgTable) solver.getAvgTable(true);

        System.out.println("\n=== DES Results ===");
        System.out.println("Queue1 QLen: " + avgTable.getQLen().get(0));
        System.out.println("Queue2 QLen: " + avgTable.getQLen().get(1));
        System.out.println("Queue1 Util: " + avgTable.getUtil().get(0));
        System.out.println("Queue2 Util: " + avgTable.getUtil().get(1));
        System.out.println("Queue1 RespT: " + avgTable.getRespT().get(0));
        System.out.println("Queue2 RespT: " + avgTable.getRespT().get(1));
        System.out.println("Queue1 Tput: " + avgTable.getTput().get(0));
        System.out.println("Queue2 Tput: " + avgTable.getTput().get(1));
        System.out.println("\nTotal QLen (should be 2.0): " + (avgTable.getQLen().get(0) + avgTable.getQLen().get(1)));

        System.out.println("\n=== Expected JMT Results (for comparison) ===");
        System.out.println("Queue1 QLen: ~0.85");
        System.out.println("Queue2 QLen: ~1.15");
    }
}
