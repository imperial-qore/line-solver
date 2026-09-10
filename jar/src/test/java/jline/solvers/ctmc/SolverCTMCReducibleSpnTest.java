package jline.solvers.ctmc;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import jline.GlobalConstants;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.Mode;
import jline.lang.processes.Exp;
import jline.lang.nodes.Place;
import jline.lang.nodes.Transition;
import jline.solvers.SolverOptions;
import jline.util.graph.DirectedGraph;
import jline.util.matrix.Matrix;

/**
 * Reducible SPN with a transient SCC feeding two disjoint recurrent classes.
 *
 * <p>Three tokens sit in P1. T1 (rate 1) moves all three into cycle A (P2 &lt;-&gt; P4),
 * T2 (rate 2) moves all three into cycle B (P3 &lt;-&gt; P5). Nothing ever returns a token
 * to P1, so the choice is irreversible and the chain decomposes into one transient
 * SCC (the initial marking) and two BSCCs.
 *
 * <p>The reference is CLOSED FORM, not a recorded output. The branch is the T1/T2
 * race, P(A)=1/3 and P(B)=2/3, and each cycle is a closed two-place exponential
 * cycle whose marginal is geometric in the rate ratio:
 * E[P2]=108/175, E[P4]=67/175, E[P3]=772/671, E[P5]=570/671.
 * An independent external GSPN tool reproduces these to five decimals.
 *
 * <p>This model is the regression for three defects that were each silent: a
 * zero-visit station had its capacity zeroed, deleting the declared initial marking
 * from the state space; the SCC support graph was built with a sign test, which
 * shatters a matrix-exponential generator; and a chain reducible within a single
 * weak component was handed to the plain solve as a singular system.
 *
 * <p>See _kb/11-conventions-and-gotchas.md.
 */
public class SolverCTMCReducibleSpnTest {

    private static final double TOL = 1e-6;
    private static final double EXACT_P2 = 108.0 / 175.0;
    private static final double EXACT_P3 = 772.0 / 671.0;
    private static final double EXACT_P4 = 67.0 / 175.0;
    private static final double EXACT_P5 = 570.0 / 671.0;

    private static Network buildModel() {
        Network model = new Network("spn_twobscc");
        Place p1 = new Place(model, "P1");
        Place p2 = new Place(model, "P2");
        Place p3 = new Place(model, "P3");
        Place p4 = new Place(model, "P4");
        Place p5 = new Place(model, "P5");
        Transition t1 = new Transition(model, "T1");
        Transition t2 = new Transition(model, "T2");
        Transition ta = new Transition(model, "TA");
        Transition ta2 = new Transition(model, "TA2");
        Transition tb = new Transition(model, "TB");
        Transition tb2 = new Transition(model, "TB2");
        ClosedClass jc = new ClosedClass(model, "C", 3, p1, 0);

        addMode(t1, "mA", 1.0, jc, p1, 3, p2, 3);
        addMode(t2, "mB", 2.0, jc, p1, 3, p3, 3);
        addMode(ta, "a1", 3.0, jc, p2, 1, p4, 1);
        addMode(ta2, "a2", 4.0, jc, p4, 1, p2, 1);
        addMode(tb, "b1", 5.0, jc, p3, 1, p5, 1);
        addMode(tb2, "b2", 6.0, jc, p5, 1, p3, 1);

        RoutingMatrix r = model.initRoutingMatrix();
        r.set(jc, jc, p1, t1, 0.5);
        r.set(jc, jc, p1, t2, 0.5);
        r.set(jc, jc, t1, p2, 1.0);
        r.set(jc, jc, t2, p3, 1.0);
        r.set(jc, jc, p2, ta, 1.0);
        r.set(jc, jc, ta, p4, 1.0);
        r.set(jc, jc, p4, ta2, 1.0);
        r.set(jc, jc, ta2, p2, 1.0);
        r.set(jc, jc, p3, tb, 1.0);
        r.set(jc, jc, tb, p5, 1.0);
        r.set(jc, jc, p5, tb2, 1.0);
        r.set(jc, jc, tb2, p3, 1.0);
        model.link(r);
        return model;
    }

    private static void addMode(Transition t, String name, double rate, ClosedClass jc,
                                Place in, int nIn, Place out, int nOut) {
        Mode mode = t.addMode(name);
        t.setDistribution(mode, new Exp(rate));
        t.setEnablingConditions(mode, jc, in, nIn);
        t.setFiringOutcome(mode, jc, out, nOut);
    }

    private static SolverOptions options() {
        SolverOptions o = new SolverCTMC(buildModel()).defaultOptions();
        o.cutoff = new Matrix(1, 1);
        o.cutoff.set(0, 0, 3);
        o.keep = true;
        return o;
    }

    @Test
    public void testMeansMatchTheClosedFormMixture() {
        SolverCTMC solver = new SolverCTMC(buildModel(), options());
        Matrix qlen = solver.getAvgQLen();
        // station order follows the Place declaration order: P1..P5
        assertEquals(EXACT_P2, qlen.get(1, 0), TOL, "P2");
        assertEquals(EXACT_P3, qlen.get(2, 0), TOL, "P3");
        assertEquals(EXACT_P4, qlen.get(3, 0), TOL, "P4");
        assertEquals(EXACT_P5, qlen.get(4, 0), TOL, "P5");
    }

    @Test
    public void testMassSplitsByTheFiringRaceNotUniformly() {
        // P(A)=1/3 comes from the T1/T2 rates; a uniform 50/50 fallback must not pass
        SolverCTMC solver = new SolverCTMC(buildModel(), options());
        Matrix qlen = solver.getAvgQLen();
        double massA = qlen.get(1, 0) + qlen.get(3, 0);
        double massB = qlen.get(2, 0) + qlen.get(4, 0);
        assertEquals(1.0, massA, TOL, "cycle A mass");
        assertEquals(2.0, massB, TOL, "cycle B mass");
    }

    @Test
    public void testGeneratorKeepsItsTransientScc() {
        SolverCTMC solver = new SolverCTMC(buildModel(), options());
        solver.getAvgQLen();
        Matrix q = solver.getGenerator().infGen;
        int n = q.getNumRows();
        Matrix adj = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                // arc by MAGNITUDE, never by sign
                if (i != j && Math.abs(q.get(i, j)) > GlobalConstants.ArcTol) {
                    adj.set(i, j, 1.0);
                }
            }
        }
        DirectedGraph.SCCResult scc = new DirectedGraph(adj).stronglyconncomp();
        int nscc = scc.recurrent.length;
        int nbscc = 0;
        for (int c = 0; c < nscc; c++) {
            if (scc.recurrent[c]) {
                nbscc++;
            }
        }
        assertTrue(nscc > 1, "the chain must stay reducible, got a single SCC");
        assertTrue(nbscc >= 2, "expected at least two recurrent classes, got " + nbscc);
        assertEquals(1, nscc - nbscc,
                "expected exactly one transient SCC (the initial marking); zero means "
                        + "the declared marking was dropped from the state space");
    }
}
