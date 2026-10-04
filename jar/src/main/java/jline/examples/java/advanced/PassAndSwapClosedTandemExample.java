package jline.examples.java.advanced;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.nc.SolverNC;
import jline.util.matrix.Matrix;
import jline.VerboseLevel;

import java.util.List;

/**
 * Closed tandem of two pass-and-swap (PAS) queues, reproducing Figures 5 and 6
 * of Comte and Dorsman, "Pass-and-Swap Queues" (2021, arXiv:2009.12299).
 *
 *   Topology (Fig. 6):   --> PASQueue1 --> PASQueue2 -->   (closed tandem)
 *
 * Six classes, one customer of each; both queues share the swapping graph of
 * Figure 5 (edges 1-3, 1-4, 2-4, 2-5, 3-6, 4-6, 5-6; 0-based at the JAR
 * boundary). Head-only single-server service. The departing customer is chosen
 * by the pass-and-swap scan and routed to the back of the other queue.
 *
 * A closed PAS network with a non-empty swapping graph is REDUCIBLE (the
 * pass-and-swap mechanism conserves the placement order, so each placement
 * order is a separate recurrent class). The initial job placement is therefore
 * a REQUIRED model input: it selects the recurrent component on which the
 * stationary distribution is supported. Here we use the placement of Fig. 6a:
 * queue 1 = (1,2,3,4,5,6) (class 1 oldest, at the head), queue 2 empty. The
 * PAS state list is 1-based (jobClass+1, 0 = empty), so setState uses 1..6.
 */
public class PassAndSwapClosedTandemExample {

    public static Network model() {
        final double mu1 = 1.0;
        final double mu2 = 1.3;          // head-only single-server rates
        int[][] edges = {{0, 2}, {0, 3}, {1, 3}, {1, 4}, {2, 5}, {3, 5}, {4, 5}};
        Matrix G = new Matrix(6, 6);
        for (int[] e : edges) { G.set(e[0], e[1], 1); G.set(e[1], e[0], 1); }

        Network model = new Network("PASclosedTandem");
        Queue q1 = new Queue(model, "PASQueue1", SchedStrategy.PAS);
        Queue q2 = new Queue(model, "PASQueue2", SchedStrategy.PAS);
        ClosedClass[] cls = new ClosedClass[6];
        for (int r = 0; r < 6; r++) {
            cls[r] = new ClosedClass(model, "Class" + (r + 1), 1, q1);
        }
        q1.setService((Matrix c) -> mu1);   // head-only: mu(prefix)=mu1
        q2.setService((Matrix c) -> mu2);
        q1.setSwapGraph(G); q1.setNumberOfServers(1); q1.setCap(6);
        q2.setSwapGraph(G); q2.setNumberOfServers(1); q2.setCap(6);
        RoutingMatrix P = model.initRoutingMatrix();
        for (int r = 0; r < 6; r++) {
            P.set(cls[r], cls[r], q1, q2, 1.0);
            P.set(cls[r], cls[r], q2, q1, 1.0);
        }
        model.link(P);

        // Required initial placement (Fig. 6a): queue 1 = (1,2,3,4,5,6) ordered,
        // class 1 oldest at the head. State list is 1-based.
        Matrix placement = new Matrix(1, 6);
        for (int i = 0; i < 6; i++) placement.set(0, i, i + 1);
        q1.setState(placement);
        return model;
    }

    public static NetworkAvgTable ctmc() {
        SolverOptions opt = SolverCTMC.defaultOptions();
        opt.method = "exact";   // pin the state-space path
        opt.cutoff = Matrix.singleton(6);
        return new SolverCTMC(model(), opt).getAvgTable();
    }

    /**
     * Closed tandem of two pass-and-swap queues, Fig. 6 (pas_closed_tandem_fig6.m).
     */
    public static void pas_closed_tandem_fig6() {
        System.out.println("=== Closed tandem of two PAS queues (Fig 5/6), CTMC ===");
        ctmc().print();
    }

    /**
     * Importance-sampling normalizing-constant analysis of the same closed
     * pass-and-swap tandem (pas_nc_sampling.m).
     *
     * <p>A non-empty swap graph makes the ordered-state chain reducible; the
     * recurrent communicating class carries the per-class product form
     * pi(c) = Phi_1 Phi_2 / G_C. SolverNC importance sampling estimates G_C and
     * the mean queue lengths by auto-normalized importance sampling
     * (Pfqn_pas_is), scaling to populations where the exact ordered-state CTMC
     * is intractable. Here it is validated against that exact CTMC.</p>
     *
     * <p>Reference: Casale, Comte and Dorsman (2026).</p>
     */
    public static void pas_nc_sampling() {
        System.out.println("=== CTMC (exact) ===");
        NetworkAvgTable Tc = ctmc();
        System.out.println(
                "=== SolverNC method='sampling' (importance sampling, 5e5 samples) ===");
        SolverOptions ncOpt = SolverNC.defaultOptions();
        ncOpt.method = "sampling";
        ncOpt.samples = 500000;
        ncOpt.seed = 777;
        ncOpt.verbose = VerboseLevel.SILENT;
        NetworkAvgTable Tn = new SolverNC(model(), ncOpt).getAvgTable();

        List<Double> qc = Tc.getQLen();
        List<Double> qn = Tn.getQLen();
        List<String> stations = Tc.getStationNames();
        List<String> classes = Tc.getClassNames();
        System.out.printf("%n  station     class    CTMC      NC-samp%n");
        double err = 0;
        for (int i = 0; i < qc.size(); i++) {
            System.out.printf("  %-9s  %-6s %9.5f %9.5f%n",
                    stations.get(i), classes.get(i), qc.get(i), qn.get(i));
            err = Math.max(err, Math.abs(qn.get(i) - qc.get(i)));
        }
        System.out.printf(
                "%nNC-samp vs CTMC:  max|dQ| = %.3e (importance-sampling noise)%n", err);
        if (err > 5e-2) {
            throw new RuntimeException(
                    "NC sampling does not match the exact CTMC within IS tolerance");
        }
        System.out.println("PASS: SolverNC 'sampling' matches the exact CTMC within "
                + "importance-sampling noise.");
    }

    public static void main(String[] args) {
        pas_closed_tandem_fig6();
        pas_nc_sampling();
    }
}
