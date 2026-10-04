package jline.examples.java.advanced;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import jline.solvers.NetworkAvgTable;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.ldes.SolverLDES;
import jline.util.matrix.Matrix;

/**
 * Pass-and-swap (PAS) / order-independent queue examples.
 *
 * A PAS station is parameterized by the total service-rate function mu(c) of the
 * ordered state vector c (a 1xN row Matrix of 0-based class indices, c(0) the
 * oldest job) and by a swapping graph G; both are properties of the Queue.
 * See Dorsman and Gardner (2024), "New directions in pass-and-swap queues",
 * Queueing Systems 107:205-256. The stationary distribution is the
 * order-independent product form and is invariant to the swapping graph.
 *
 * Both the CTMC solver (exact) and the LDES discrete-event simulator support
 * PAS stations.
 */
public class PassAndSwapExample {

    /** M/M/K order-independent queue model: mu(c) = min(n, K). */
    public static Network mmkModel() {
        final int K = 2;
        double[] lam = {0.7, 0.5};
        int C = 4;
        Network model = new Network("PASmmk");
        Source src = new Source(model, "Source");
        Queue q = new Queue(model, "PASQueue", SchedStrategy.PAS);
        Sink snk = new Sink(model, "Sink");
        OpenClass c1 = new OpenClass(model, "Class1");
        OpenClass c2 = new OpenClass(model, "Class2");
        src.setArrival(c1, new Exp(lam[0]));
        src.setArrival(c2, new Exp(lam[1]));
        q.setService((Matrix c) -> (double) Math.min(c.getNumCols(), K));
        q.setSwapGraph(new Matrix(2, 2));   // empty graph -> plain OI queue
        q.setNumberOfServers(K);
        q.setCap(C);
        model.link(Network.serialRouting(src, q, snk));
        return model;
    }

    /**
     * Five-class, three-server compatibility model (paper Figs 1-2) with the
     * swapping graph E = {(0,2),(0,4),(1,3),(2,3),(3,4)} (0-based).
     * mu(c) = number of servers compatible with at least one class present in c.
     */
    public static Network compatibilityModel() {
        final int[][] comp = {
                {1, 0, 0, 1, 0},   // server 1: classes {0,3}
                {0, 1, 0, 1, 0},   // server 2: classes {1,3}
                {0, 0, 1, 0, 1}};  // server 3: classes {2,4}
        double[] lam = {0.5, 0.4, 0.3, 0.2, 0.1};
        int C = 3;
        Network model = new Network("PAScompatibility");
        Source src = new Source(model, "Source");
        Queue q = new Queue(model, "PASQueue", SchedStrategy.PAS);
        Sink snk = new Sink(model, "Sink");
        OpenClass[] cls = new OpenClass[5];
        for (int r = 0; r < 5; r++) {
            cls[r] = new OpenClass(model, "Class" + (r + 1));
            src.setArrival(cls[r], new Exp(lam[r]));
        }
        q.setService((Matrix c) -> {
            double total = 0;
            for (int s = 0; s < 3; s++) {
                boolean any = false;
                for (int i = 0; i < c.getNumCols(); i++) {
                    if (comp[s][(int) c.get(0, i)] != 0) { any = true; break; }
                }
                if (any) total += 1;
            }
            return total;
        });
        Matrix G = new Matrix(5, 5);
        int[][] edges = {{0, 2}, {0, 4}, {1, 3}, {2, 3}, {3, 4}};
        for (int[] e : edges) { G.set(e[0], e[1], 1); G.set(e[1], e[0], 1); }
        q.setSwapGraph(G);
        q.setNumberOfServers(3);
        q.setCap(C);
        model.link(Network.serialRouting(src, q, snk));
        return model;
    }

    /**
     * Three-class PAS queue whose swapping graph contains a self-loop (class 1)
     * and an edge (2,3); one dedicated server per class so mu(c) = number of
     * distinct classes present. Self-loops are permitted (Dorsman & Gardner
     * 2024, Sect. 2.3) and the product form remains invariant to the graph.
     */
    public static Network selfloopModel() {
        double[] lam = {0.6, 0.4, 0.3};
        int C = 3;
        Network model = new Network("PASselfloop");
        Source src = new Source(model, "Source");
        Queue q = new Queue(model, "PASQueue", SchedStrategy.PAS);
        Sink snk = new Sink(model, "Sink");
        OpenClass[] cls = new OpenClass[3];
        for (int r = 0; r < 3; r++) {
            cls[r] = new OpenClass(model, "Class" + (r + 1));
            src.setArrival(cls[r], new Exp(lam[r]));
        }
        // mu(c) = number of distinct classes present (one server per class)
        q.setService((Matrix c) -> {
            boolean[] seen = new boolean[3];
            for (int i = 0; i < c.getNumCols(); i++) seen[(int) c.get(0, i)] = true;
            int d = 0;
            for (boolean b : seen) if (b) d++;
            return (double) d;
        });
        Matrix G = new Matrix(3, 3);
        G.set(0, 0, 1);              // self-loop on class 1
        G.set(1, 2, 1); G.set(2, 1, 1);
        q.setSwapGraph(G);
        q.setNumberOfServers(3);
        q.setCap(C);
        model.link(Network.serialRouting(src, q, snk));
        return model;
    }

    /** M/M/K order-independent queue solved by CTMC. */
    public static NetworkAvgTable pas_mmk() {
        Network model = mmkModel();
        SolverOptions opt = SolverCTMC.defaultOptions();
        opt.method = "exact";   // pin the state-space path
        opt.cutoff = Matrix.singleton(4);
        return new SolverCTMC(model, opt).getAvgTable();
    }

    /** Five-class compatibility model solved by CTMC. */
    public static NetworkAvgTable pas_compatibility_5class() {
        Network model = compatibilityModel();
        SolverOptions opt = SolverCTMC.defaultOptions();
        opt.method = "exact";   // pin the state-space path
        opt.cutoff = Matrix.singleton(3);
        return new SolverCTMC(model, opt).getAvgTable();
    }

    /** Self-loop swapping-graph model solved by CTMC. */
    public static NetworkAvgTable pas_selfloop() {
        Network model = selfloopModel();
        SolverOptions opt = SolverCTMC.defaultOptions();
        opt.method = "exact";   // pin the state-space path
        opt.cutoff = Matrix.singleton(3);
        return new SolverCTMC(model, opt).getAvgTable();
    }

    /**
     * Pass-and-swap saturation handling for large / unbounded buffers
     * (pas_saturation.m).
     *
     * <p>An order-independent service rate mu(c) usually SATURATES: beyond some
     * per-class count threshold tau_r, adding more class-r jobs no longer
     * changes mu (a compatibility / rank function saturates once each present
     * class has one job; an M/M/K rate saturates at K jobs). A large finite
     * buffer therefore approximates an unbounded queue cheaply.</p>
     *
     * <p>Uses the Comte (thesis Sect. 4.1) compatibility queue: I = 2 classes,
     * S = 3 servers, S_1 = {1,2}, S_2 = {2,3}, so mu({1}) = mu({2}) = 3 and
     * mu({1,2}) = 5, constant once both classes are present. It checks that (i)
     * LDES matches CTMC on a finite buffer and (ii) an unset/infinite buffer
     * raises a clean, actionable error.</p>
     */
    public static void pas_saturation() {
        // compatibility comp(s,i) = 1 iff server s serves class i; capacities cap_s
        final int[][] comp = {
                {1, 0},   // server 1 -> class 1
                {1, 1},   // server 2 -> classes 1, 2 (shared)
                {0, 1}};  // server 3 -> class 2
        final double[] capS = {2, 1, 2};
        final double[] lambda = {1.5, 1.0};
        final int I = 2;
        final int S = 3;

        // Rank function: total capacity of servers compatible with a present class.
        ServiceRate muRank = c -> {
            boolean[] present = new boolean[I];
            for (int j = 0; j < c.getNumCols(); j++) {
                int r = (int) c.get(0, j);
                if (r >= 0 && r < I) {
                    present[r] = true;
                }
            }
            double total = 0;
            for (int s = 0; s < S; s++) {
                for (int i = 0; i < I; i++) {
                    if (present[i] && comp[s][i] == 1) {
                        total += capS[s];
                        break;
                    }
                }
            }
            return total;
        };

        // (1) a large buffer approximates the unbounded queue; LDES matches CTMC
        final int CAP = 8;
        SolverOptions copt = SolverCTMC.defaultOptions();
        copt.method = "exact";
        copt.cutoff = Matrix.singleton(CAP);
        double[] Qc = saturationQLen(
                new SolverCTMC(saturationModel(muRank, CAP, I, S, lambda), copt), I);
        SolverOptions lopt = SolverLDES.defaultOptions();
        lopt.samples = 300000;
        lopt.seed = 23000;
        double[] Ql = saturationQLen(
                new SolverLDES(saturationModel(muRank, CAP, I, S, lambda), lopt), I);
        System.out.printf("open compatibility PAS, buffer=%d:%n", CAP);
        System.out.println("  CTMC  QLen = " + java.util.Arrays.toString(Qc));
        System.out.println("  LDES  QLen = " + java.util.Arrays.toString(Ql));
        double rel = 0;
        for (int r = 0; r < I; r++) {
            rel = Math.max(rel, Math.abs(Qc[r] - Ql[r]) / Math.max(Qc[r], 1e-12));
        }
        System.out.printf("  LDES vs CTMC max rel|dQ| = %.2f%%%n", 100 * rel);
        if (rel > 0.03) {
            throw new RuntimeException("LDES does not match CTMC within simulation noise");
        }

        // (2) clean error for an unset / infinite buffer
        System.out.println("\ninfinite-buffer PAS raises a clean, actionable error:");
        try {
            SolverOptions bad = SolverLDES.defaultOptions();
            bad.samples = 1000;
            bad.seed = 1;
            new SolverLDES(saturationModel(muRank, -1, I, S, lambda), bad).getAvgTable();
            throw new IllegalStateException("expected a finite-buffer error");
        } catch (IllegalStateException e) {
            throw e;
        } catch (RuntimeException e) {
            System.out.println("  " + e.getMessage());
        }

        System.out.println("\nPASS: PAS saturation handling matches CTMC on a large buffer "
                + "and errors cleanly on infinite buffers.");
    }

    /** The total order-independent service rate of an ordered state. */
    private interface ServiceRate {
        double mu(Matrix c);
    }

    /** Open compatibility PAS queue; a negative cap leaves the buffer unbounded. */
    private static Network saturationModel(ServiceRate mu, int cap, int I, int S,
                                           double[] lambda) {
        Network model = new Network("PASsaturation");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "PASQueue", SchedStrategy.PAS);
        Sink sink = new Sink(model, "Sink");
        for (int r = 0; r < I; r++) {
            OpenClass jobclass = new OpenClass(model, "Class" + (r + 1));
            source.setArrival(jobclass, new Exp(lambda[r]));
        }
        queue.setService((Matrix c) -> mu.mu(c));
        queue.setSwapGraph(new Matrix(I, I));   // empty graph: plain OI queue
        queue.setNumberOfServers(S);
        if (cap > 0) {
            queue.setCap(cap);
        }
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    /** Mean queue length of the PAS station, one entry per class. */
    private static double[] saturationQLen(jline.solvers.NetworkSolver solver, int I) {
        Matrix QN = solver.getAvgQLen();
        jline.lang.NetworkStruct sn = solver.getModel().getStruct();
        int ind = sn.nodenames.indexOf("PASQueue");
        int sidx = (int) sn.nodeToStation.get(ind);
        double[] q = new double[I];
        for (int r = 0; r < I; r++) {
            q[r] = Math.round(QN.get(sidx, r) * 1e5) / 1e5;
        }
        return q;
    }

    public static void main(String[] args) {
        System.out.println("=== PAS M/M/K order-independent queue (CTMC) ===");
        pas_mmk().print();
        System.out.println("=== PAS M/M/K order-independent queue (LDES) ===");
        SolverOptions lopt = SolverLDES.defaultOptions();
        lopt.seed = 12345;
        lopt.samples = 2000000;
        new SolverLDES(mmkModel(), lopt).getAvgTable().print();

        System.out.println("=== PAS five-class compatibility (CTMC) ===");
        pas_compatibility_5class().print();
        System.out.println("=== PAS five-class compatibility (LDES) ===");
        SolverOptions lopt2 = SolverLDES.defaultOptions();
        lopt2.seed = 12345;
        lopt2.samples = 4000000;
        new SolverLDES(compatibilityModel(), lopt2).getAvgTable().print();

        System.out.println("=== PAS self-loop swapping graph (CTMC) ===");
        pas_selfloop().print();
        System.out.println("=== PAS self-loop swapping graph (LDES) ===");
        SolverOptions lopt3 = SolverLDES.defaultOptions();
        lopt3.seed = 12345;
        lopt3.samples = 4000000;
        new SolverLDES(selfloopModel(), lopt3).getAvgTable().print();
    }
}
