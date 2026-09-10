package jline.solvers.ctmc;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * CTMC validation on CLOSED pass-and-swap (PAS) networks, i.e. the REDUCIBLE case.
 *
 * The existing PAS coverage ({@code PassAndSwapExamplesTest}) is all open
 * (Source -> PASQueue -> Sink) and hence irreducible, so it does not exercise the
 * analyzer's reducible branch at all. A closed PAS network with a non-empty swapping
 * graph conserves the placement order (Comte and Dorsman, "Pass-and-Swap Queues", 2021,
 * arXiv:2009.12299, Prop. 2), so its generator carries one bottom SCC per reachable
 * placement order and the stationary distribution is component-dependent. The declared
 * initial placement selects the component, and the analyzer must honour it by seeding
 * the absorption probabilities from that state rather than mixing the components.
 *
 * The oracle is INDEPENDENT rather than a recorded golden: the pass-and-swap Markov
 * chain is built here directly from the algorithm over ordered placements and solved on
 * its own. A recorded golden cannot detect a consistent error in the thing it records.
 */
public class SolverCTMCClosedPasTest {

    private static final double TOL = 1e-6;
    private static final double MU1 = 1.0;
    private static final double MU2 = 1.3;
    private static final int R = 6;

    /** Figure 5 swapping graph, 0-based: edges 1-3, 1-4, 2-4, 2-5, 3-6, 4-6, 5-6. */
    private static Matrix swapGraph() {
        int[][] edges = {{0, 2}, {0, 3}, {1, 3}, {1, 4}, {2, 5}, {3, 5}, {4, 5}};
        Matrix G = new Matrix(R, R);
        for (int[] e : edges) {
            G.set(e[0], e[1], 1);
            G.set(e[1], e[0], 1);
        }
        return G;
    }

    /**
     * Closed tandem of Fig. 6. When {@code placement} is null no initial placement is
     * declared, which leaves the model ill-posed.
     */
    private static Network buildClosedTandem(int[] placement) {
        Matrix G = swapGraph();
        Network model = new Network("PASclosedTandem");
        Queue q1 = new Queue(model, "PASQueue1", SchedStrategy.PAS);
        Queue q2 = new Queue(model, "PASQueue2", SchedStrategy.PAS);
        ClosedClass[] cls = new ClosedClass[R];
        for (int r = 0; r < R; r++) {
            cls[r] = new ClosedClass(model, "Class" + (r + 1), 1, q1);
        }
        q1.setService((Matrix c) -> MU1);   // head-only service: mu(prefix) = mu1
        q2.setService((Matrix c) -> MU2);
        q1.setSwapGraph(G); q1.setNumberOfServers(1); q1.setCap(R);
        q2.setSwapGraph(G); q2.setNumberOfServers(1); q2.setCap(R);
        RoutingMatrix P = model.initRoutingMatrix();
        for (int r = 0; r < R; r++) {
            P.set(cls[r], cls[r], q1, q2, 1.0);
            P.set(cls[r], cls[r], q2, q1, 1.0);
        }
        model.link(P);
        if (placement != null) {
            // The PAS state list is 1-based (jobClass+1, 0 = empty)
            Matrix m = new Matrix(1, placement.length);
            for (int i = 0; i < placement.length; i++) {
                m.set(0, i, placement[i]);
            }
            q1.setState(m);
        }
        return model;
    }

    private static Matrix ctmcQLen(Network model) {
        SolverOptions opt = SolverCTMC.defaultOptions();
        opt.cutoff = Matrix.singleton(R);
        return new SolverCTMC(model, opt).getAvgQLen();
    }

    // ---------------------------------------------------------------------
    // Independent reference: the pass-and-swap chain over ordered placements
    // ---------------------------------------------------------------------

    /**
     * One pass-and-swap completion at position p of the ordered list: classes shift one
     * step along the swap chain and the served slot is removed. Returns the new list
     * with the departing class appended as its last element.
     */
    private static int[] psStep(int[] lst, int p, Matrix G) {
        List<Integer> chain = new ArrayList<Integer>();
        chain.add(p);
        int cur = p;
        while (true) {
            int nxt = -1;
            for (int j = cur + 1; j < lst.length; j++) {
                if (G.get(lst[cur] - 1, lst[j] - 1) > 0) {
                    nxt = j;
                    break;
                }
            }
            if (nxt < 0) {
                break;
            }
            chain.add(nxt);
            cur = nxt;
        }
        int dep = lst[chain.get(chain.size() - 1)];
        int[] tmp = Arrays.copyOf(lst, lst.length);
        for (int i = 0; i < chain.size() - 1; i++) {
            tmp[chain.get(i + 1)] = lst[chain.get(i)];
        }
        int removed = chain.get(0);
        int[] out = new int[lst.length];        // last slot carries the departing class
        int w = 0;
        for (int i = 0; i < tmp.length; i++) {
            if (i != removed) {
                out[w++] = tmp[i];
            }
        }
        out[lst.length - 1] = dep;
        return out;
    }

    private static String key(int[] l1, int[] l2) {
        return Arrays.toString(l1) + "|" + Arrays.toString(l2);
    }

    private static int[] append(int[] l, int v) {
        int[] out = Arrays.copyOf(l, l.length + 1);
        out[l.length] = v;
        return out;
    }

    /**
     * Enumerate the pass-and-swap chain reachable from the Fig. 6a placement, solve it,
     * and return the mean per-class occupancy of the two queues as a flat station-major
     * vector of length 2*R.
     */
    private static double[] psReference(Matrix G) {
        List<int[]> q1 = new ArrayList<int[]>();
        List<int[]> q2 = new ArrayList<int[]>();
        Map<String, Integer> index = new HashMap<String, Integer>();
        int[] init1 = new int[R];
        for (int i = 0; i < R; i++) {
            init1[i] = i + 1;
        }
        int[] init2 = new int[0];
        q1.add(init1);
        q2.add(init2);
        index.put(key(init1, init2), 0);

        List<int[]> edges = new ArrayList<int[]>();   // {from, to}
        List<Double> rates = new ArrayList<Double>();
        int fr = 0;
        while (fr < q1.size()) {
            int[] l1 = q1.get(fr);
            int[] l2 = q2.get(fr);
            if (l1.length > 0) {
                int[] res = psStep(l1, 0, G);
                int[] nl = Arrays.copyOf(res, l1.length - 1);
                int dep = res[l1.length - 1];
                int[] t2 = append(l2, dep);
                String k = key(nl, t2);
                Integer j = index.get(k);
                if (j == null) {
                    j = q1.size();
                    index.put(k, j);
                    q1.add(nl);
                    q2.add(t2);
                }
                edges.add(new int[]{fr, j});
                rates.add(MU1);
            }
            if (l2.length > 0) {
                int[] res = psStep(l2, 0, G);
                int[] nl = Arrays.copyOf(res, l2.length - 1);
                int dep = res[l2.length - 1];
                int[] t1 = append(l1, dep);
                String k = key(t1, nl);
                Integer j = index.get(k);
                if (j == null) {
                    j = q1.size();
                    index.put(k, j);
                    q1.add(t1);
                    q2.add(nl);
                }
                edges.add(new int[]{fr, j});
                rates.add(MU2);
            }
            fr++;
        }

        int ns = q1.size();
        Matrix Q = new Matrix(ns, ns);
        for (int e = 0; e < edges.size(); e++) {
            int[] ij = edges.get(e);
            Q.set(ij[0], ij[1], Q.get(ij[0], ij[1]) + rates.get(e));
        }
        for (int i = 0; i < ns; i++) {
            double s = 0.0;
            for (int j = 0; j < ns; j++) {
                if (j != i) {
                    s += Q.get(i, j);
                }
            }
            Q.set(i, i, -s);
        }
        Matrix pi = jline.api.mc.Ctmc_solve.ctmc_solve(Q);

        double[] out = new double[2 * R];
        for (int s = 0; s < ns; s++) {
            double p = pi.get(0, s);
            for (int c : q1.get(s)) {
                out[c - 1] += p;
            }
            for (int c : q2.get(s)) {
                out[R + c - 1] += p;
            }
        }
        return out;
    }

    // ---------------------------------------------------------------------
    // Tests
    // ---------------------------------------------------------------------

    @Test
    public void testClosedTandemMatchesPassAndSwapChain() {
        Network model = buildClosedTandem(new int[]{1, 2, 3, 4, 5, 6});
        Matrix QN = ctmcQLen(model);
        double[] ref = psReference(swapGraph());
        for (int i = 0; i < 2; i++) {
            for (int r = 0; r < R; r++) {
                assertEquals(ref[i * R + r], QN.get(i, r), TOL,
                        "closed PAS QLen at station " + (i + 1) + " class " + (r + 1));
            }
        }
    }

    @Test
    public void testClosedTandemConservesPopulation() {
        Network model = buildClosedTandem(new int[]{1, 2, 3, 4, 5, 6});
        Matrix QN = ctmcQLen(model);
        for (int r = 0; r < R; r++) {
            assertEquals(1.0, QN.get(0, r) + QN.get(1, r), TOL,
                    "class " + (r + 1) + " population not conserved across the two queues");
        }
    }

    @Test
    public void testPlacementSelectsTheComponent() {
        // The reversed placement lies in a different recurrent component, so it must
        // give a different stationary vector. If the analyzer mixed the components (or
        // picked one by size) the two would agree, which is exactly the old defect.
        Matrix asc = ctmcQLen(buildClosedTandem(new int[]{1, 2, 3, 4, 5, 6}));
        Matrix desc = ctmcQLen(buildClosedTandem(new int[]{6, 5, 4, 3, 2, 1}));
        double maxDiff = 0.0;
        for (int i = 0; i < 2; i++) {
            for (int r = 0; r < R; r++) {
                maxDiff = Math.max(maxDiff, Math.abs(asc.get(i, r) - desc.get(i, r)));
            }
        }
        assertTrue(maxDiff > 1e-3,
                "the declared placement must select the recurrent component, but the "
                        + "ascending and descending placements agree to " + maxDiff);
    }

    @Test
    public void testClosedTandemWithoutPlacementIsRejected() {
        // Several recurrent components and no declared placement is ill-posed: refuse
        // rather than fabricate a component.
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            @Override
            public void execute() {
                ctmcQLen(buildClosedTandem(null));
            }
        });
    }
}
