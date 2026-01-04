/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.advanced;

import jline.VerboseLevel;
import jline.lang.Network;
import jline.solvers.NetworkAvgTable;
import jline.solvers.jmt.JMT;
import jline.solvers.nc.SolverNC;

/**
 * Examples demonstrating Finite Capacity Region (FCR) usage.
 */
public class FCRegionExamples {

    /**
     * Run the multiclass FCR blocking example.
     */
    public static void fcr_oqnwaitq() {
        System.out.println("=== FCR Multiclass Blocking Example ===\n");

        Network model = FCRegionModel.fcr_oqnwaitq();
        JMT solver = new JMT(model, "seed", 23000, "samples", 50000, "verbose", VerboseLevel.SILENT);
        NetworkAvgTable avgTable = solver.getAvgTable();

        System.out.println(avgTable);
    }

    /**
     * Run the multiclass FCR dropping example.
     */
    public static void fcr_oqndrop() {
        System.out.println("=== FCR Multiclass Dropping Example ===\n");

        Network model = FCRegionModel.fcr_oqndrop();
        JMT solver = new JMT(model, "seed", 23000, "samples", 50000, "verbose", VerboseLevel.SILENT);
        NetworkAvgTable avgTable = solver.getAvgTable();

        System.out.println(avgTable);
    }

    /**
     * Compare FCR blocking with standard M/M/1.
     * Demonstrates that FCR with blocking behaves identically to M/M/1.
     */
    public static void fcr_mm1waitq() {
        System.out.println("=== Comparison: FCR Blocking vs M/M/1 ===\n");

        // FCR blocking model
        Network model1 = FCRegionModel.fcr_mm1waitq();
        JMT solver1 = new JMT(model1, "seed", 23000, "samples", 100000, "verbose", VerboseLevel.SILENT);
        NetworkAvgTable avgTable1 = solver1.getAvgTable();

        // Standard M/M/1
        Network model2 = FCRegionModel.mm1();
        JMT solver2 = new JMT(model2, "seed", 23000, "samples", 100000, "verbose", VerboseLevel.SILENT);
        NetworkAvgTable avgTable2 = solver2.getAvgTable();

        // Compare results (queue is at index 1)
        int queueIdx = 1;
        System.out.println("Queue Metrics:");
        System.out.printf("                    FCR Blocking    M/M/1%n");
        System.out.printf("Queue Length:       %.4f          %.4f%n", avgTable1.getQLen().get(queueIdx), avgTable2.getQLen().get(queueIdx));
        System.out.printf("Utilization:        %.4f          %.4f%n", avgTable1.getUtil().get(queueIdx), avgTable2.getUtil().get(queueIdx));
        System.out.printf("Response Time:      %.4f          %.4f%n", avgTable1.getRespT().get(queueIdx), avgTable2.getRespT().get(queueIdx));
        System.out.printf("Throughput:         %.4f          %.4f%n", avgTable1.getTput().get(queueIdx), avgTable2.getTput().get(queueIdx));

        // Theoretical M/M/1 values
        double arrivalRate = 0.5;
        double serviceRate = 1.0;
        double rho = arrivalRate / serviceRate;
        double theoretical_Q = rho / (1 - rho);
        double theoretical_U = rho;
        double theoretical_R = 1 / (serviceRate - arrivalRate);
        double theoretical_X = arrivalRate;

        System.out.println("\nTheoretical M/M/1:");
        System.out.printf("Queue Length:       %.4f%n", theoretical_Q);
        System.out.printf("Utilization:        %.4f%n", theoretical_U);
        System.out.printf("Response Time:      %.4f%n", theoretical_R);
        System.out.printf("Throughput:         %.4f%n", theoretical_X);

        System.out.println("\n==> FCR with blocking matches M/M/1 (jobs wait, no loss)");
    }

    /**
     * Compare FCR dropping with M/M/1/K.
     * Demonstrates that FCR with dropping behaves like M/M/1/K.
     */
    public static void fcr_mm1kdrop() {
        int K = 2;
        System.out.printf("=== Comparison: FCR Dropping vs M/M/1/K (K=%d) ===%n%n", K);

        // FCR dropping model
        Network model1 = FCRegionModel.fcr_mm1kdrop();
        JMT solver1 = new JMT(model1, "seed", 23000, "samples", 100000, "verbose", VerboseLevel.SILENT);
        NetworkAvgTable avgTable1 = solver1.getAvgTable();

        // Standard M/M/1/K
        Network model2 = FCRegionModel.mm1k();
        JMT solver2 = new JMT(model2, "seed", 23000, "samples", 100000, "verbose", VerboseLevel.SILENT);
        NetworkAvgTable avgTable2 = solver2.getAvgTable();

        // Compare results (queue is at index 1)
        int queueIdx = 1;
        System.out.println("Queue Metrics:");
        System.out.printf("                    FCR Dropping    M/M/1/K%n");
        System.out.printf("Queue Length:       %.4f          %.4f%n", avgTable1.getQLen().get(queueIdx), avgTable2.getQLen().get(queueIdx));
        System.out.printf("Utilization:        %.4f          %.4f%n", avgTable1.getUtil().get(queueIdx), avgTable2.getUtil().get(queueIdx));
        System.out.printf("Response Time:      %.4f          %.4f%n", avgTable1.getRespT().get(queueIdx), avgTable2.getRespT().get(queueIdx));
        System.out.printf("Throughput:         %.4f          %.4f%n", avgTable1.getTput().get(queueIdx), avgTable2.getTput().get(queueIdx));

        // Theoretical M/M/1/K values
        double arrivalRate = 0.8;
        double serviceRate = 1.0;
        double rho = arrivalRate / serviceRate;
        double p0 = (1 - rho) / (1 - Math.pow(rho, K + 1));
        double pK = p0 * Math.pow(rho, K);
        double theoretical_Q = rho / (1 - rho) - (K + 1) * Math.pow(rho, K + 1) / (1 - Math.pow(rho, K + 1));
        double lambda_eff = arrivalRate * (1 - pK);
        double theoretical_U = lambda_eff / serviceRate;
        double theoretical_X = lambda_eff;
        double theoretical_R = theoretical_Q / theoretical_X;

        System.out.println("\nTheoretical M/M/1/K:");
        System.out.printf("Queue Length:       %.4f%n", theoretical_Q);
        System.out.printf("Utilization:        %.4f%n", theoretical_U);
        System.out.printf("Response Time:      %.4f%n", theoretical_R);
        System.out.printf("Throughput:         %.4f%n", theoretical_X);
        System.out.printf("Blocking Prob:      %.4f%n", pK);

        System.out.println("\n==> FCR with dropping matches M/M/1/K (jobs lost when full)");
    }

    /**
     * Demonstrate FCR with multiple constraint types.
     */
    public static void fcr_constraints() {
        System.out.println("=== FCR Constraint Types Demo ===\n");
        System.out.println("Region covers: Queue1, Queue2");
        System.out.println("Global max jobs: 6");
        System.out.println("HighPriority max: 4 (blocking)");
        System.out.println("LowPriority max: 3 (dropping)\n");

        Network model = FCRegionModel.fcr_constraints();
        JMT solver = new JMT(model, "seed", 23000, "samples", 100000, "verbose", VerboseLevel.SILENT);
        NetworkAvgTable avgTable = solver.getAvgTable();

        System.out.println(avgTable);

        System.out.println("\nNote: HighPriority jobs experience delays when region is full.");
        System.out.println("LowPriority jobs are dropped when their class limit or global limit is reached.");
    }

    /**
     * Demonstrate loss network analysis with NC solver.
     * NC solver uses Erlang fixed-point approximation for open
     * loss networks with a single Delay inside an FCR with DROP policy.
     */
    public static void fcr_lossn() {
        System.out.println("=== FCR Loss Network (NC Solver) ===\n");
        System.out.println("Model: Single Delay node in FCR with DROP policy");
        System.out.println("Global max: 5, Class1 max: 3, Class2 max: 3\n");

        Network model = FCRegionModel.fcr_lossn();

        // Run NC solver (uses lossn method)
        System.out.println("Running NC solver (lossn method)...");
        SolverNC solverNC = new SolverNC(model, "verbose", VerboseLevel.SILENT);
        NetworkAvgTable avgTableNC = solverNC.getAvgTable();
        System.out.println("\nNC Results:");
        System.out.println(avgTableNC);

        // Run JMT for comparison
        System.out.println("\nRunning JMT for comparison...");
        JMT solverJMT = new JMT(model, "seed", 23000, "samples", 500000, "verbose", VerboseLevel.SILENT);
        NetworkAvgTable avgTableJMT = solverJMT.getAvgTable();
        System.out.println("\nJMT Results:");
        System.out.println(avgTableJMT);

        // Compare throughputs at Delay node (index 1)
        int delayIdx = 1;
        System.out.println("\n=== Comparison (Delay node) ===");
        System.out.printf("%-15s %12s %12s %12s%n", "Class", "NC Tput", "JMT Tput", "Rel Err");
        for (int r = 0; r < 2; r++) {
            double tputNC = avgTableNC.getTput().get(delayIdx + r);
            double tputJMT = avgTableJMT.getTput().get(delayIdx + r);
            double relErr = Math.abs(tputNC - tputJMT) / tputJMT * 100;
            System.out.printf("%-15s %12.4f %12.4f %10.2f%%%n",
                    "Class" + (r + 1), tputNC, tputJMT, relErr);
        }

        System.out.println("\nNote: NC uses analytical Erlang fixed-point approximation.");
        System.out.println("JMT uses discrete-event simulation.");
    }

    public static void main(String[] args) {
        fcr_oqnwaitq();
        System.out.println();
        fcr_oqndrop();
        System.out.println();
        fcr_mm1waitq();
        System.out.println();
        fcr_mm1kdrop();
        System.out.println();
        fcr_constraints();
        System.out.println();
        fcr_lossn();
    }
}
