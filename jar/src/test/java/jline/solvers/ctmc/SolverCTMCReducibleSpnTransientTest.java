package jline.solvers.ctmc;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.Mode;
import jline.lang.processes.Exp;
import jline.lang.nodes.Place;
import jline.lang.nodes.Transition;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

/**
 * Transient analysis of the reducible SPN with a transient SCC.
 *
 * <p>Same model as {@link SolverCTMCReducibleSpnTest}: three tokens in P1, T1 (rate 1)
 * drains them into cycle A, T2 (rate 2) into cycle B, and nothing returns to P1. The
 * transient has a closed form independent of the internal cycle dynamics:
 * E[P1](t) = 3 exp(-3t), E[P2]+E[P4] = (1-exp(-3t)), E[P3]+E[P5] = 2(1-exp(-3t)).
 *
 * <p>See _kb/11-conventions-and-gotchas.md.
 */
public class SolverCTMCReducibleSpnTransientTest {

    private static final double TSPAN = 4.0;
    private static final double RACE_RATE = 3.0;
    private static final double POP = 3.0;

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

    private static SolverCTMC solver() {
        Network model = buildModel();
        SolverOptions o = new SolverCTMC(model).defaultOptions();
        o.cutoff = new Matrix(1, 1);
        o.cutoff.set(0, 0, 3);
        o.timespan = new double[]{0.0, TSPAN};
        return new SolverCTMC(model, o);
    }

    @Test
    public void testDrainedMassSplitsOneToTwoAtEveryTime() {
        // the 1:2 split is set by the T1/T2 rates and holds at every t; a uniform
        // branch would give 1:1
        SolverCTMC s = solver();
        s.getTranAvg();
        jline.solvers.SolverResult res = s.getResults();
        Matrix t = res.t;
        double maxErrA = 0.0;
        double maxErrB = 0.0;
        for (int i = 0; i < t.getNumRows(); i++) {
            double drained = POP * (1.0 - Math.exp(-RACE_RATE * t.get(i, 0)));
            double massA = res.QNt[1][0].get(i, 0) + res.QNt[3][0].get(i, 0);
            double massB = res.QNt[2][0].get(i, 0) + res.QNt[4][0].get(i, 0);
            maxErrA = Math.max(maxErrA, Math.abs(massA - drained / 3.0));
            maxErrB = Math.max(maxErrB, Math.abs(massB - 2.0 * drained / 3.0));
        }
        assertTrue(maxErrA < 1e-4, "cycle A fill curve off by " + maxErrA);
        assertTrue(maxErrB < 1e-4, "cycle B fill curve off by " + maxErrB);
    }

    @Test
    public void testTheCurvesAreNotConstant() {
        // guard against the silent fallback: seeding with the stationary vector gives a
        // flat curve that satisfies conservation and the steady-state limit vacuously
        SolverCTMC s = solver();
        s.getTranAvg();
        jline.solvers.SolverResult res = s.getResults();
        double maxSpan = 0.0;
        for (int ist = 0; ist < 5; ist++) {
            Matrix c = res.QNt[ist][0];
            maxSpan = Math.max(maxSpan, c.elementMax() - c.elementMin());
        }
        assertTrue(maxSpan > 1e-3, "every curve is flat: no transient analysis ran");
    }

    @Test
    public void testTransientConvergesToTheSteadyStateMixture() {
        SolverCTMC s = solver();
        s.getTranAvg();
        jline.solvers.SolverResult res = s.getResults();
        int last = res.t.getNumRows() - 1;
        assertEquals(108.0 / 175.0, res.QNt[1][0].get(last, 0), 1e-3, "P2");
        assertEquals(772.0 / 671.0, res.QNt[2][0].get(last, 0), 1e-3, "P3");
        assertEquals(67.0 / 175.0, res.QNt[3][0].get(last, 0), 1e-3, "P4");
        assertEquals(570.0 / 671.0, res.QNt[4][0].get(last, 0), 1e-3, "P5");
    }

    @Test
    public void testP1DrainsAtTheRaceRate() {
        SolverCTMC s = solver();
        s.getTranAvg();
        jline.solvers.SolverResult res = s.getResults();
        Matrix t = res.t;
        Matrix p1 = res.QNt[0][0];
        double maxErr = 0.0;
        for (int i = 0; i < t.getNumRows(); i++) {
            double exact = POP * Math.exp(-RACE_RATE * t.get(i, 0));
            maxErr = Math.max(maxErr, Math.abs(p1.get(i, 0) - exact));
        }
        assertTrue(maxErr < 1e-4, "E[P1](t) deviates from 3exp(-3t) by " + maxErr);
    }

    @Test
    public void testTokensAreConservedAtEveryTime() {
        SolverCTMC s = solver();
        s.getTranAvg();
        jline.solvers.SolverResult res = s.getResults();
        Matrix t = res.t;
        double maxDev = 0.0;
        for (int i = 0; i < t.getNumRows(); i++) {
            double total = 0.0;
            for (int ist = 0; ist < 5; ist++) {
                total += res.QNt[ist][0].get(i, 0);
            }
            maxDev = Math.max(maxDev, Math.abs(total - POP));
        }
        assertTrue(maxDev < 1e-4, "token count drifts by " + maxDev);
    }
}
