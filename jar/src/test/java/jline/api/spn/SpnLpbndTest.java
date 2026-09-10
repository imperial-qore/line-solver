package jline.api.spn;

import jline.api.spn.Spn_lpbnd.SpnLpBounds;
import jline.api.spn.Spn_lpbnd.SpnLpOptions;
import jline.lang.ClosedClass;
import jline.lang.Mode;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.nodes.Place;
import jline.lang.nodes.Transition;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import jline.solvers.ba.SolverBA;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Regression tests for {@link Spn_lpbnd} and the SolverBA {@code spnlp.*} family.
 *
 * <p>THE ACCEPTANCE TEST IS THE REFERENCE'S OWN TABLE 2. Liu (1998) publishes
 * four bound columns for the four-server production line of its Fig. 2b, on
 * five rate vectors, and this class asserts all four to three decimals. Those
 * numbers are what a transcription error in any constraint family would move:
 * the polytope has twelve of them and the published optimum is a function of
 * the whole set, so a bracket check alone would not catch a family that was
 * dropped or mis-signed.</p>
 *
 * <p>The bracket checks are the other half. A bound that agreed with the paper
 * and still failed to contain the exact answer would be wrong whatever it
 * agreed with. The exact values here are the ones {@code SpnPfTest} already
 * pins for the fork-join net, which shares no code path with this LP.</p>
 *
 * <p>Reference: Z. Liu, "Performance Analysis of Stochastic Timed Petri Nets
 * Using Linear Programming Approach", IEEE Trans. Software Engineering 24(11),
 * 1998, 1014-1030.</p>
 */
public class SpnLpbndTest {

    /** Liu Table 2, p. 1023: mu, l.b., u.b.2, o.l.b., o.u.b. */
    private static final double[][] TABLE2 = {
            {1.000, 1.25, 2.00, 0.50, 1.165, 2.000, 0.930, 2.000},
            {1.000, 1.25, 2.00, 2.50, 1.829, 3.529, 1.481, 4.000},
            {1.000, 1.25, 1.25, 2.50, 1.581, 3.333, 1.333, 4.000},
            {1.000, 1.25, 1.25, 1.00, 1.359, 3.333, 1.111, 4.000},
            {1.111, 1.111, 1.111, 1.111, 1.350, 2.963, 1.111, 4.444}};

    /** The exact fork-join token counts SpnPfTest pins, from the product form. */
    private static final double[] FJ_TOKENS =
            {0.710262429604, 1.864198441415, 1.864198441415, 0.425539128981};
    private static final double FJ_MODEX = 0.607360882508;

    private static double sum(double[] v) {
        double s = 0;
        for (int i = 0; i < v.length; i++) {
            s += v[i];
        }
        return s;
    }

    /**
     * Fig. 2b: four servers, blocking before service, buffers of 3, 2 and 4.
     *
     * <p>(p5, p2), (p4, p1) and (p3, p0) are the three buffer pairs, each
     * conserved at its capacity, so the net is a strongly connected marked graph
     * and every transition carries the same throughput.</p>
     */
    private static Network prodLine(double[] mu) {
        Network model = new Network("liu98");
        Place p5 = new Place(model, "p5");
        Place p4 = new Place(model, "p4");
        Place p3 = new Place(model, "p3");
        Place p2 = new Place(model, "p2");
        Place p1 = new Place(model, "p1");
        Place p0 = new Place(model, "p0");
        Transition t1 = new Transition(model, "t1");
        Transition t2 = new Transition(model, "t2");
        Transition t3 = new Transition(model, "t3");
        Transition t4 = new Transition(model, "t4");
        ClosedClass jc = new ClosedClass(model, "Class1", 9, p2);

        Mode m1 = t1.addMode("m1");
        t1.setDistribution(m1, new Exp(mu[0]));
        t1.setNumberOfServers(m1, 1);
        t1.setEnablingConditions(m1, jc, p2, 1);
        t1.setFiringOutcome(m1, jc, p5, 1);

        Mode m2 = t2.addMode("m2");
        t2.setDistribution(m2, new Exp(mu[1]));
        t2.setNumberOfServers(m2, 1);
        t2.setEnablingConditions(m2, jc, p5, 1);
        t2.setEnablingConditions(m2, jc, p1, 1);
        t2.setFiringOutcome(m2, jc, p4, 1);
        t2.setFiringOutcome(m2, jc, p2, 1);

        Mode m3 = t3.addMode("m3");
        t3.setDistribution(m3, new Exp(mu[2]));
        t3.setNumberOfServers(m3, 1);
        t3.setEnablingConditions(m3, jc, p4, 1);
        t3.setEnablingConditions(m3, jc, p0, 1);
        t3.setFiringOutcome(m3, jc, p3, 1);
        t3.setFiringOutcome(m3, jc, p1, 1);

        Mode m4 = t4.addMode("m4");
        t4.setDistribution(m4, new Exp(mu[3]));
        t4.setNumberOfServers(m4, 1);
        t4.setEnablingConditions(m4, jc, p3, 1);
        t4.setFiringOutcome(m4, jc, p0, 1);

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jc, jc, p2, t1, 1.0);
        P.set(jc, jc, t1, p5, 1.0);
        P.set(jc, jc, p5, t2, 1.0);
        P.set(jc, jc, p1, t2, 1.0);
        P.set(jc, jc, t2, p4, 1.0);
        P.set(jc, jc, t2, p2, 1.0);
        P.set(jc, jc, p4, t3, 1.0);
        P.set(jc, jc, p0, t3, 1.0);
        P.set(jc, jc, t3, p3, 1.0);
        P.set(jc, jc, t3, p1, 1.0);
        P.set(jc, jc, p3, t4, 1.0);
        P.set(jc, jc, t4, p0, 1.0);
        model.link(P);

        p5.setState(Matrix.singleton(0));
        p4.setState(Matrix.singleton(0));
        p3.setState(Matrix.singleton(0));
        p2.setState(Matrix.singleton(3));
        p1.setState(Matrix.singleton(2));
        p0.setState(Matrix.singleton(4));
        return model;
    }

    /**
     * P0 -(Tf)-&gt; P1 + P2 -(Tj)-&gt; P3 -(Tb)-&gt; P0.
     *
     * <p>Tf consumes one token and produces two, so the marking is not a
     * conserved population; the place invariant is (2, 1, 1, 2). This is what
     * exercises the WEIGHTED form of the invariant family -- the reference
     * writes that family for unweighted cycles, and an implementation that only
     * summed token counts would produce an unbounded polytope here.</p>
     */
    private static Network forkJoinSpn(int ntokens) {
        Network model = new Network("fj");
        Place[] pl = new Place[4];
        for (int i = 0; i < 4; i++) {
            pl[i] = new Place(model, "P" + i);
        }
        Transition tf = new Transition(model, "Tf");
        Transition tj = new Transition(model, "Tj");
        Transition tb = new Transition(model, "Tb");
        ClosedClass jc = new ClosedClass(model, "C", ntokens, pl[0]);

        Mode mf = tf.addMode("f");
        tf.setDistribution(mf, new Exp(1.3));
        tf.setNumberOfServers(mf, 1);
        tf.setEnablingConditions(mf, jc, pl[0], 1);
        tf.setFiringOutcome(mf, jc, pl[1], 1);
        tf.setFiringOutcome(mf, jc, pl[2], 1);

        Mode mj = tj.addMode("j");
        tj.setDistribution(mj, new Exp(0.7));
        tj.setNumberOfServers(mj, 1);
        tj.setEnablingConditions(mj, jc, pl[1], 1);
        tj.setEnablingConditions(mj, jc, pl[2], 1);
        tj.setFiringOutcome(mj, jc, pl[3], 1);

        Mode mb = tb.addMode("b");
        tb.setDistribution(mb, new Exp(1.9));
        tb.setNumberOfServers(mb, 1);
        tb.setEnablingConditions(mb, jc, pl[3], 1);
        tb.setFiringOutcome(mb, jc, pl[0], 1);

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jc, jc, pl[0], tf, 1.0);
        P.set(jc, jc, tf, pl[1], 1.0);
        P.set(jc, jc, tf, pl[2], 1.0);
        P.set(jc, jc, pl[1], tj, 1.0);
        P.set(jc, jc, pl[2], tj, 1.0);
        P.set(jc, jc, tj, pl[3], 1.0);
        P.set(jc, jc, pl[3], tb, 1.0);
        P.set(jc, jc, tb, pl[0], 1.0);
        model.link(P);
        for (int i = 0; i < 4; i++) {
            pl[i].setState(Matrix.singleton(i == 0 ? ntokens : 0));
        }
        return model;
    }

    /** One place, one mode, a self-loop: nothing can move. */
    private static Network selfLoop(int ntokens) {
        Network model = new Network("one");
        Place p = new Place(model, "P");
        Transition t = new Transition(model, "T");
        ClosedClass jc = new ClosedClass(model, "C", ntokens, p);
        Mode m = t.addMode("fire");
        t.setDistribution(m, new Exp(2.0));
        t.setNumberOfServers(m, 1);
        t.setEnablingConditions(m, jc, p, 1);
        t.setFiringOutcome(m, jc, p, 1);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jc, jc, p, t, 1.0);
        P.set(jc, jc, t, p, 1.0);
        model.link(P);
        p.setState(Matrix.singleton(ntokens));
        return model;
    }

    // ------------------------------------------------------------------
    @Test
    public void testLiu98Table2UpperMatchesPublished() {
        for (int c = 0; c < TABLE2.length; c++) {
            double[] mu = Arrays.copyOfRange(TABLE2[c], 0, 4);
            SpnLpOptions o = new SpnLpOptions();
            o.markovian = true;
            SpnLpBounds b = Spn_lpbnd.spn_lpbnd(prodLine(mu).getStruct(), o);
            assertEquals(TABLE2[c][5], sum(b.modeTput[1]), 5e-4,
                    "Liu Table 2 case " + (c + 1) + " u.b.2");
        }
    }

    /**
     * The published l.b. column is the lower side WITH the liveness rows.
     *
     * <p>They are opt-in because they hold only on a live net, and this one is:
     * a strongly connected marked graph with a token on every cycle. Without
     * them the lower side falls back to the operational value, which is the
     * second assertion here and is why the default is not a silent loss.</p>
     */
    @Test
    public void testLiu98Table2LowerNeedsLiveness() {
        for (int c = 0; c < TABLE2.length; c++) {
            double[] mu = Arrays.copyOfRange(TABLE2[c], 0, 4);
            SpnLpOptions live = new SpnLpOptions();
            live.markovian = true;
            live.assumelive = true;
            SpnLpBounds bl = Spn_lpbnd.spn_lpbnd(prodLine(mu).getStruct(), live);
            assertEquals(TABLE2[c][4], sum(bl.modeTput[0]), 5e-4,
                    "Liu Table 2 case " + (c + 1) + " l.b.");
            SpnLpBounds bp = Spn_lpbnd.spn_lpbnd(prodLine(mu).getStruct(), new SpnLpOptions());
            assertTrue(sum(bp.modeTput[0]) <= sum(bl.modeTput[0]) + 1e-9);
        }
    }

    @Test
    public void testLiu98Table2OperationalMatchesPublished() {
        for (int c = 0; c < TABLE2.length; c++) {
            double[] mu = Arrays.copyOfRange(TABLE2[c], 0, 4);
            SpnLpOptions o = new SpnLpOptions();
            o.markovian = false;
            o.assumelive = true;
            SpnLpBounds b = Spn_lpbnd.spn_lpbnd(prodLine(mu).getStruct(), o);
            assertEquals(TABLE2[c][6], sum(b.modeTput[0]), 5e-4,
                    "Liu Table 2 case " + (c + 1) + " o.l.b.");
            assertEquals(TABLE2[c][7], sum(b.modeTput[1]), 5e-4,
                    "Liu Table 2 case " + (c + 1) + " o.u.b.");
        }
    }

    // ------------------------------------------------------------------
    @Test
    public void testForkJoinBracketContainsTheProductForm() {
        for (int pass = 0; pass < 2; pass++) {
            SpnLpOptions o = new SpnLpOptions();
            o.markovian = pass == 0;
            SpnLpBounds b = Spn_lpbnd.spn_lpbnd(forkJoinSpn(3).getStruct(), o);
            assertEquals(4, b.nplacelevels);
            for (int l = 0; l < 4; l++) {
                assertTrue(b.tokens[0][l] <= FJ_TOKENS[l] + 1e-6,
                        "lower " + b.tokens[0][l] + " > exact " + FJ_TOKENS[l]);
                assertTrue(FJ_TOKENS[l] <= b.tokens[1][l] + 1e-6,
                        "exact " + FJ_TOKENS[l] + " > upper " + b.tokens[1][l]);
                assertTrue(b.placeTput[0][l] <= FJ_MODEX + 1e-6);
                assertTrue(FJ_MODEX <= b.placeTput[1][l] + 1e-6);
            }
        }
    }

    /**
     * The a priori caps come off the WEIGHTED invariants, and the fork-join net
     * is where an unweighted reading would be wrong: its two minimal-support
     * invariants are (1,1,0,1) and (1,0,1,1) at value 3, so every level is
     * capped at 3 rather than at the token count.
     */
    @Test
    public void testForkJoinLevelCapsComeFromTheInvariants() {
        SpnLpBounds b = Spn_lpbnd.spn_lpbnd(forkJoinSpn(3).getStruct(), new SpnLpOptions());
        for (int l = 0; l < 4; l++) {
            assertEquals(3.0, b.bound[l], 1e-12);
        }
    }

    /**
     * A [0, B] bracket is a missing constraint family, not a loose bound: it
     * means the polytope does not constrain the objective at all.
     */
    @Test
    public void testBracketIsNotVacuous() {
        SpnLpBounds b = Spn_lpbnd.spn_lpbnd(forkJoinSpn(3).getStruct(), new SpnLpOptions());
        for (int l = 0; l < b.nplacelevels; l++) {
            assertTrue(b.tokens[1][l] - b.tokens[0][l] < b.bound[l] - 1e-9,
                    "level " + l + " bracket is the whole box");
        }
    }

    @Test
    public void testMarkovianIsAtLeastAsTightAsOperational() {
        SpnLpOptions mk = new SpnLpOptions();
        SpnLpOptions op = new SpnLpOptions();
        op.markovian = false;
        SpnLpBounds bm = Spn_lpbnd.spn_lpbnd(forkJoinSpn(3).getStruct(), mk);
        SpnLpBounds bo = Spn_lpbnd.spn_lpbnd(forkJoinSpn(3).getStruct(), op);
        for (int l = 0; l < bm.nplacelevels; l++) {
            assertTrue(bo.tokens[0][l] <= bm.tokens[0][l] + 1e-9);
            assertTrue(bm.tokens[1][l] <= bo.tokens[1][l] + 1e-9);
        }
    }

    @Test
    public void testDegenerateSelfLoopIsExact() {
        SpnLpBounds b = Spn_lpbnd.spn_lpbnd(selfLoop(3).getStruct(), new SpnLpOptions());
        assertEquals(3.0, b.tokens[0][0], 1e-9);
        assertEquals(3.0, b.tokens[1][0], 1e-9);
    }

    // ------------------------------------------------------------------
    @Test
    public void testSolverBaOffersOnlySpnlpOnAPetriNet() {
        String[] valid = new SolverBA(forkJoinSpn(3)).listValidMethods();
        List<String> got = Arrays.asList(valid);
        assertEquals(4, got.size(), "got " + got);
        assertTrue(got.contains("spnlp.upper") && got.contains("spnlp.lower")
                && got.contains("spnlp.op.upper") && got.contains("spnlp.op.lower"),
                "got " + got);
    }

    @Test
    public void testSolverBaEndToEnd() throws Exception {
        String[] methods = {"spnlp.upper", "spnlp.lower", "spnlp.op.upper", "spnlp.op.lower"};
        for (int i = 0; i < methods.length; i++) {
            SolverOptions o = SolverBA.defaultOptions();
            o.method = methods[i];
            SolverBA s = new SolverBA(forkJoinSpn(3), o);
            s.runAnalyzer();
            Matrix qn = s.result.QN;
            // U = Q at an INF station, which is what a Place is
            for (int r = 0; r < qn.getNumRows(); r++) {
                assertEquals(qn.get(r, 0), s.result.UN.get(r, 0), 1e-12);
                if (methods[i].endsWith("upper")) {
                    assertTrue(qn.get(r, 0) >= FJ_TOKENS[r] - 1e-6);
                } else {
                    assertTrue(qn.get(r, 0) <= FJ_TOKENS[r] + 1e-6);
                }
            }
        }
    }

    // ------------------------------------------------------------------
    /** An Erlang mode has no marking-only state, but it does have a mean. */
    @Test
    public void testPhaseTypeRefusedByMarkovianAcceptedByOperational() {
        Network model = new Network("ph");
        Place p0 = new Place(model, "P0");
        Place p1 = new Place(model, "P1");
        Transition t0 = new Transition(model, "T0");
        Transition t1 = new Transition(model, "T1");
        ClosedClass jc = new ClosedClass(model, "C", 3, p0);
        Mode m0 = t0.addMode("a");
        t0.setDistribution(m0, Erlang.fitMeanAndOrder(1.0, 2));
        t0.setNumberOfServers(m0, 1);
        t0.setEnablingConditions(m0, jc, p0, 1);
        t0.setFiringOutcome(m0, jc, p1, 1);
        Mode m1 = t1.addMode("b");
        t1.setDistribution(m1, new Exp(1.5));
        t1.setNumberOfServers(m1, 1);
        t1.setEnablingConditions(m1, jc, p1, 1);
        t1.setFiringOutcome(m1, jc, p0, 1);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jc, jc, p0, t0, 1.0);
        P.set(jc, jc, t0, p1, 1.0);
        P.set(jc, jc, p1, t1, 1.0);
        P.set(jc, jc, t1, p0, 1.0);
        model.link(P);
        p0.setState(Matrix.singleton(3));
        p1.setState(Matrix.singleton(0));

        RuntimeException ex = assertThrows(RuntimeException.class,
                () -> Spn_lpbnd.spn_lpbnd(model.getStruct(), new SpnLpOptions()));
        assertTrue(ex.getMessage().contains("phase-type"), ex.getMessage());

        SpnLpOptions op = new SpnLpOptions();
        op.markovian = false;
        SpnLpBounds b = Spn_lpbnd.spn_lpbnd(model.getStruct(), op);
        for (int l = 0; l < b.nplacelevels; l++) {
            assertTrue(Double.isFinite(b.tokens[0][l]) && Double.isFinite(b.tokens[1][l]));
        }
    }

    @Test
    public void testInfiniteServerModeRefusedByName() {
        Network model = new Network("is");
        Place p = new Place(model, "P");
        Transition t = new Transition(model, "T");
        ClosedClass jc = new ClosedClass(model, "C", 2, p);
        Mode m = t.addMode("fire");
        t.setDistribution(m, new Exp(1.0));
        t.setNumberOfServers(m, Integer.valueOf(4));
        t.setEnablingConditions(m, jc, p, 1);
        t.setFiringOutcome(m, jc, p, 1);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jc, jc, p, t, 1.0);
        P.set(jc, jc, t, p, 1.0);
        model.link(P);
        p.setState(Matrix.singleton(2));
        RuntimeException ex = assertThrows(RuntimeException.class,
                () -> Spn_lpbnd.spn_lpbnd(model.getStruct(), new SpnLpOptions()));
        assertTrue(ex.getMessage().contains("servers"), ex.getMessage());
    }
}
