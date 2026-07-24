package jline.solvers.ssa;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Coxian;
import jline.lang.processes.Distribution;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the NRM phase expansion at the INF/PS-family stations against an
 * exact analytic oracle.
 * <p>
 * The fixture is a closed two-job network whose station under test has mean
 * service 0.5 and whose Delay has mean service 1.0. Both PS and INF are
 * insensitive stations under BCMP, so the exact means depend on the service
 * distribution only through its mean:
 * </p>
 * <ul>
 *   <li>PS  station: Q = 0.8,     U = 0.6,     X = 1.2</li>
 *   <li>INF station: Q = 2/3,     U = 2/3,     X = 4/3</li>
 * </ul>
 * <p>
 * Insensitivity is what makes this a usable regression: every distribution
 * below has mean 0.5, so all of them must return the SAME exact value, and no
 * CTMC solve is needed to obtain it. The Exp row is the control -- with
 * nph == 1 the slot map degenerates to the pre-expansion flat class index, so
 * Exp must reproduce its pre-expansion value exactly and any Exp deviation
 * indicts the layout rather than the phase logic.
 * </p>
 * <p>
 * These tests exist because the expansion was previously covered only by
 * all-exponential models, in which the dependency graph's arithmetic decode of
 * a state index is still correct and its defects are therefore invisible: the
 * failure mode is a quiet bias that GROWS with distance from exponential (see
 * _kb/log.md [2026-07-17]). The Erlang/HyperExp/Coxian rows are the ones that
 * carry the signal; a green all-exponential suite proves nothing about them.
 * </p>
 * <p>
 * Tolerance: at SAMPLES the single-seed noise floor on Q is sd ~= 0.001
 * (0.1%), calibrated by an 8-seed spread on the Exp control. RTOL is set to 1%
 * -- several sd above that floor, so the fixed-seed assertion is not flaky,
 * yet below the signature of the decode defect this covers (Erlang-3 +1.3%,
 * HyperExp +4.0%), so a recurrence fails rather than passing inside a loose
 * band.
 * </p>
 */
public class SolverSSANrmPhaseTest {

    private static final int SAMPLES = 2000000;
    private static final int SEED = 23000;
    /** Relative tolerance on the simulated means at a fixed seed. */
    private static final double RTOL = 0.01;

    private static final double MEAN = 0.5;

    private static Distribution[] dists() {
        return new Distribution[]{
                new Exp(1.0 / MEAN),
                Erlang.fitMeanAndOrder(MEAN, 2),
                Erlang.fitMeanAndOrder(MEAN, 3),
                HyperExp.fitMeanAndSCV(MEAN, 4.0),
                Coxian.fitMeanAndSCV(MEAN, 3.0),
        };
    }

    private static String[] labels() {
        return new String[]{"Exp", "Erlang-2", "Erlang-3", "HyperExp-4", "Coxian-3"};
    }

    private static Network model(SchedStrategy sched, Distribution svc) {
        Network model = new Network("phase_" + sched);
        Queue queue = new Queue(model, "Q1", sched);
        Delay delay = new Delay(model, "Delay");
        ClosedClass c1 = new ClosedClass(model, "C1", 2, queue);
        queue.setService(c1, svc);
        delay.setService(c1, new Exp(1.0));
        model.link(model.serialRouting(queue, delay));
        return model;
    }

    private static double nrmQLenAtQueue(Network model) {
        SolverSSA solver = new SolverSSA(model);
        solver.options.method = "nrm";
        solver.options.samples = SAMPLES;
        solver.options.seed = SEED;
        Matrix q = solver.getAvgQLen();
        // Station 0 is the Queue: serialRouting(queue, delay) registers it first.
        return q.get(0, 0);
    }

    private static void assertInsensitive(SchedStrategy sched, double exact, String label) {
        Distribution[] ds = dists();
        String[] names = labels();
        for (int i = 0; i < ds.length; i++) {
            double got = nrmQLenAtQueue(model(sched, ds[i]));
            double err = Math.abs(got - exact) / exact;
            assertTrue(err < RTOL, label + "/" + names[i]
                    + ": NRM QLen " + got + " deviates from the exact insensitive value "
                    + exact + " by " + (100.0 * err) + "%");
        }
    }

    /** PS is insensitive: every mean-0.5 distribution must return Q = 0.8. */
    @Test
    public void testPhaseExpansionPsInsensitive() {
        assertInsensitive(SchedStrategy.PS, 0.8, "PS");
    }

    /** INF is insensitive: every mean-0.5 distribution must return Q = 2/3. */
    @Test
    public void testPhaseExpansionInfInsensitive() {
        assertInsensitive(SchedStrategy.INF, 2.0 / 3.0, "INF");
    }
}
