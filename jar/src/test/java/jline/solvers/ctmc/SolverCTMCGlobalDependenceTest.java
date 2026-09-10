package jline.solvers.ctmc;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.nodes.Queue;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.lang.constant.SolverType;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Regression tests for {@link Network#setGlobalDependence}, the globally
 * state-dependent rate scaling phi(n) over the FULL (nstations x nclasses)
 * population matrix, and for the CTMC generator that consumes it.
 *
 * <p>WHY IT IS NOT setJointDependence WITH A WIDER ARGUMENT. Every existing hook
 * -- lldscaling, cdscaling, jdscaling -- is handed the population of ONE station,
 * so none can express a rate that reads the whole state. A Whittle network needs
 * exactly that, and so does bandwidth sharing, where one route holds several
 * links at once.
 *
 * <p>The oracles, in increasing strength: an identity phi must change nothing; a
 * phi reproducing a per-station alpha_i(n_i) must equal setLoadDependence, which
 * pins the fold point; two stations sharing one unit of capacity by
 * phi_s(n) = n_s/|n| conserve the population; and the same balanced model with
 * Erlang service of equal mean must give the same means, which is insensitivity
 * -- the defining property of a Whittle network and the only oracle here that
 * fails if PHASE transitions go unscaled.
 *
 * <p>The MATLAB twin is line-test.git/test_ctmc_global_dependence.m, the Python
 * twin is python/tests/test_ctmc_global_dependence.py and the C++ twin is
 * cpp/tests/test_ctmc_global_dependence.cpp.
 */
public class SolverCTMCGlobalDependenceTest {

    private static final double TOL = 1e-9;

    /** Closed cyclic PS pair, one class; erlang switches to an Erlang-3 of equal mean. */
    private static Network closedPair(int njobs, boolean erlang) {
        Network model = new Network("gdPair");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "C1", njobs, q1, 0);
        if (erlang) {
            q1.setService(c1, Erlang.fitMeanAndOrder(1.0, 3));
            q2.setService(c1, Erlang.fitMeanAndOrder(0.5, 3));
        } else {
            q1.setService(c1, new Exp(1.0));
            q2.setService(c1, new Exp(2.0));
        }
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c1, q1, q2, 1.0);
        P.set(c1, c1, q2, q1, 1.0);
        model.link(P);
        return model;
    }

    private static Matrix scalar(double v) {
        Matrix m = new Matrix(1, 1);
        m.set(0, 0, v);
        return m;
    }

    private static double[] qlen(Network model) {
        NetworkAvgTable t = new SolverCTMC(model).getAvgTable();
        java.util.List<Double> q = t.getQLen();
        double[] out = new double[q.size()];
        for (int i = 0; i < q.size(); i++) out[i] = q.get(i);
        return out;
    }

    @Test
    public void identityScalingChangesNothing() {
        double[] plain = qlen(closedPair(3, false));
        Network withGd = closedPair(3, false);
        withGd.setGlobalDependence(n -> scalar(1.0), scalar(1.0));
        double[] scaled = qlen(withGd);
        assertEquals(plain.length, scaled.length);
        for (int i = 0; i < plain.length; i++) {
            assertEquals(plain[i], scaled[i], 1e-12);
        }
    }

    @Test
    public void reproducesPerStationLoadDependence() {
        Network ld = closedPair(3, false);
        Matrix alpha = new Matrix(1, 4);
        alpha.set(0, 0, 1.0);
        alpha.set(0, 1, 2.0);
        alpha.set(0, 2, 2.0);
        alpha.set(0, 3, 2.0);
        ((Queue) ld.getNodeByName("Q2")).setLoadDependence(alpha);
        double[] a = qlen(ld);

        Network gd = closedPair(3, false);
        final int M = gd.getNumberOfStations();
        final int K = gd.getNumberOfClasses();
        gd.setGlobalDependence(n -> {
            // alpha(n) = min(n, 2) at station 2, expressed through the global handle
            Matrix v = Matrix.ones(M, K);
            double nj = 0;
            for (int r = 0; r < K; r++) nj += n.get(1, r);
            if (nj > 0) {
                for (int r = 0; r < K; r++) v.set(1, r, Math.min(nj, 2.0));
            }
            return v;
        }, scalar(2.0));
        double[] b = qlen(gd);

        for (int i = 0; i < a.length; i++) {
            assertEquals(a[i], b[i], TOL);
        }
    }

    @Test
    public void sharedLinkConservesPopulationAndIsInsensitive() {
        double[] exp = qlen(sharedLink(3, false));
        double[] erl = qlen(sharedLink(3, true));
        double tot = 0;
        for (int i = 0; i < exp.length; i++) tot += exp[i];
        assertEquals(3.0, tot, TOL);
        // Insensitivity: an Erlang-3 of the same mean must not move the means
        for (int i = 0; i < exp.length; i++) {
            assertEquals(exp[i], erl[i], 1e-6);
        }
    }

    /** phi_s(n) = n_s/|n|: the two stations share one unit of capacity (balanced). */
    private static Network sharedLink(int njobs, boolean erlang) {
        Network model = closedPair(njobs, erlang);
        final int M = model.getNumberOfStations();
        final int K = model.getNumberOfClasses();
        model.setGlobalDependence(n -> {
            Matrix v = Matrix.ones(M, K);
            double tot = 0;
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) tot += n.get(i, r);
            }
            if (tot > 0) {
                for (int i = 0; i < M; i++) {
                    double ni = 0;
                    for (int r = 0; r < K; r++) ni += n.get(i, r);
                    for (int r = 0; r < K; r++) v.set(i, r, ni / tot);
                }
            }
            return v;
        }, scalar(1.0));
        return model;
    }

    @Test
    public void malformedDeclarationIsRefused() {
        Network model = closedPair(3, false);
        assertThrows(IllegalArgumentException.class,
                () -> model.setGlobalDependence(null, scalar(1.0)));
        assertThrows(IllegalArgumentException.class,
                () -> model.setGlobalDependence(n -> scalar(1.0), null));
        assertThrows(IllegalArgumentException.class,
                () -> model.setGlobalDependence(n -> scalar(1.0), scalar(-1.0)));
        // wrong output shape, refused at declaration rather than mid-generation
        assertThrows(IllegalArgumentException.class,
                () -> model.setGlobalDependence(n -> Matrix.ones(7, 3), scalar(1.0)));
    }

    /**
     * SolverSSA carries the same factorization on the sample path: phi(n) is a
     * CONSTANT within a state, so it is evaluated once per state and multiplies
     * every station service rate there. An unscaled run would report [1.5, 1.5]
     * against the exact [2, 1], so the means ARE the oracle.
     */
    @Test
    public void ssaSerialMatchesCtmc() {
        Network model = sharedLink(3, false);
        double[] exact = qlen(model);
        SolverOptions opt = new SolverOptions(SolverType.SSA);
        opt.seed = 23000;
        opt.samples = 200000;
        NetworkAvgTable t = new jline.solvers.ssa.SolverSSA(model, opt).getAvgTable();
        java.util.List<Double> q = t.getQLen();
        for (int i = 0; i < exact.length; i++) {
            assertEquals(exact[i], q.get(i), 3e-2 * Math.max(exact[i], 1.0));
        }
        // Util normalizes by the declared gd peak, as it does for cd/jd. With
        // peak 1 the two stations share one unit of capacity, so the utilizations
        // are the CTMC's 2/3 and 1/3.
        NetworkAvgTable tc = new SolverCTMC(model).getAvgTable();
        java.util.List<Double> uc = tc.getUtil();
        java.util.List<Double> us = t.getUtil();
        for (int i = 0; i < uc.size(); i++) {
            assertEquals(uc.get(i), us.get(i), 3e-2);
        }
    }

    /**
     * The NRM builds its propensities from the per-station population slice and
     * never sees the whole population matrix phi reads, so an explicit
     * method='nrm' must divert to the serial engine rather than run unscaled.
     */
    @Test
    public void ssaNrmFallsBackToSerial() {
        Network model = sharedLink(3, false);
        double[] exact = qlen(model);
        SolverOptions opt = new SolverOptions(SolverType.SSA);
        opt.method = "nrm";
        opt.seed = 23000;
        opt.samples = 200000;
        NetworkAvgTable t = new jline.solvers.ssa.SolverSSA(model, opt).getAvgTable();
        java.util.List<Double> q = t.getQLen();
        for (int i = 0; i < exact.length; i++) {
            assertEquals(exact[i], q.get(i), 3e-2 * Math.max(exact[i], 1.0));
        }
    }

    /**
     * phi(n) is a handle, so it reaches the wire only as a TABLE: the writer
     * materializes it over the lattice of the whole network state, restricted to
     * the (station,class) slots a class can occupy. The oracle is the ANSWER, not
     * the document -- a table read back at the wrong coordinate still parses.
     */
    @Test
    public void jsonRoundTripPreservesTheAnswer() throws Exception {
        Network model = sharedLink(3, false);
        com.google.gson.JsonObject doc = jline.io.LineModelIO.toJsonObject(model);
        com.google.gson.JsonObject blk =
                doc.getAsJsonObject("model").getAsJsonObject("globalDependence");
        assertEquals("globalDependent", blk.get("type").getAsString());
        assertEquals(2, blk.getAsJsonArray("slots").size());
        assertEquals(16, blk.getAsJsonObject("scaling").size());  // the (3+1)^2 box

        java.io.File f = java.io.File.createTempFile("gdwire", ".json");
        f.deleteOnExit();
        jline.io.LineModelIO.save(model, f.getAbsolutePath());
        Network back = (Network) jline.io.LineModelIO.load(f.getAbsolutePath());
        assertNotNull(back.getGlobalDependence());
        Matrix a = new SolverCTMC(model, "exact").getAvgQLen();
        Matrix b = new SolverCTMC(back, "exact").getAvgQLen();
        for (int i = 0; i < a.getNumRows(); i++) {
            assertEquals(a.get(i, 0), b.get(i, 0), 1e-10);
        }
    }
}
