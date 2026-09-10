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
    public static NetworkAvgTable mmk() {
        Network model = mmkModel();
        SolverOptions opt = SolverCTMC.defaultOptions();
        opt.method = "exact";   // pin the state-space path
        opt.cutoff = Matrix.singleton(4);
        return new SolverCTMC(model, opt).getAvgTable();
    }

    /** Five-class compatibility model solved by CTMC. */
    public static NetworkAvgTable compatibility() {
        Network model = compatibilityModel();
        SolverOptions opt = SolverCTMC.defaultOptions();
        opt.method = "exact";   // pin the state-space path
        opt.cutoff = Matrix.singleton(3);
        return new SolverCTMC(model, opt).getAvgTable();
    }

    /** Self-loop swapping-graph model solved by CTMC. */
    public static NetworkAvgTable selfloop() {
        Network model = selfloopModel();
        SolverOptions opt = SolverCTMC.defaultOptions();
        opt.method = "exact";   // pin the state-space path
        opt.cutoff = Matrix.singleton(3);
        return new SolverCTMC(model, opt).getAvgTable();
    }

    public static void main(String[] args) {
        System.out.println("=== PAS M/M/K order-independent queue (CTMC) ===");
        mmk().print();
        System.out.println("=== PAS M/M/K order-independent queue (LDES) ===");
        SolverOptions lopt = SolverLDES.defaultOptions();
        lopt.seed = 12345;
        lopt.samples = 2000000;
        new SolverLDES(mmkModel(), lopt).getAvgTable().print();

        System.out.println("=== PAS five-class compatibility (CTMC) ===");
        compatibility().print();
        System.out.println("=== PAS five-class compatibility (LDES) ===");
        SolverOptions lopt2 = SolverLDES.defaultOptions();
        lopt2.seed = 12345;
        lopt2.samples = 4000000;
        new SolverLDES(compatibilityModel(), lopt2).getAvgTable().print();

        System.out.println("=== PAS self-loop swapping graph (CTMC) ===");
        selfloop().print();
        System.out.println("=== PAS self-loop swapping graph (LDES) ===");
        SolverOptions lopt3 = SolverLDES.defaultOptions();
        lopt3.seed = 12345;
        lopt3.samples = 4000000;
        new SolverLDES(selfloopModel(), lopt3).getAvgTable().print();
    }
}
