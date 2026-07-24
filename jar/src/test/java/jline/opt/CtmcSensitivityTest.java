package jline.opt;

import jline.api.mc.Ctmc_solve;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.opt.results.SensitivityData;
import jline.opt.sensitivity.SensitivityTable;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the generator-derivative CTMC sensitivity fallback (only reachable
 * for models the product-form branches decline) against central finite
 * differences of the exact CTMC mean queue length. Closed network: Think delay
 * Exp(1) and a 2-server PS queue Exp(mu), population 3. The multiserver queue
 * makes {@code SensitivityTable.compute(model)} return null, so the CTMC branch
 * must supply d(QLen)/d(rate).
 */
public class CtmcSensitivityTest {

    private static Network build(double muQ) {
        Network model = new Network("m");
        Delay think = new Delay(model, "Think");
        Queue q = new Queue(model, "Q", SchedStrategy.PS);
        q.setNumberOfServers(2);
        ClosedClass c = new ClosedClass(model, "C", 3, think);
        think.setService(c, new Exp(1.0));
        q.setService(c, new Exp(muQ));
        model.link(Network.serialRouting(think, q));
        return model;
    }

    /** Exact E[QLen] of class C at the queue, straight from the CTMC. */
    private static double meanQlenAtQueue(double muQ) {
        Network m = build(muQ);
        SolverCTMC solver = new SolverCTMC(m);
        Matrix Q = solver.getGenerator().infGen;
        Matrix pi = Ctmc_solve.ctmc_solve(Q);
        Matrix sa = solver.getStateSpaceAggr();
        // station-major, class-minor: Q is station index 1, class C is 0, K=1 -> col 1
        int K = m.getStruct().nclasses;
        int col = 1 * K + 0;
        double e = 0.0;
        for (int s = 0; s < pi.length(); s++) {
            e += pi.get(0, s) * sa.get(s, col);
        }
        return e;
    }

    @Test
    public void testCtmcQlenSensitivityMatchesFiniteDifference() {
        double mu = 2.0;
        Network model = build(mu);

        // product-form branches must decline (multiserver), forcing the CTMC path
        assertNull(SensitivityTable.compute(model, false),
                "product-form sensitivity should be null for a multiserver queue");

        SensitivityData sens = SensitivityTable.compute(model, true);
        assertNotNull(sens, "CTMC fallback should produce a sensitivity table");
        Double analytic = sens.forKind("QLen")
                .get(SensitivityData.metricKey("Q", "C"))
                .get(SensitivityData.paramKey("Q", "C"));
        assertNotNull(analytic, "d(QLen_Q,C)/d(rate_Q,C) should be present");

        double h = 1e-4;
        double fd = (meanQlenAtQueue(mu + h) - meanQlenAtQueue(mu - h)) / (2.0 * h);

        assertTrue(Math.abs(analytic - fd) <= 1e-6 * Math.max(Math.abs(fd), 1.0),
                "analytic=" + analytic + " fd=" + fd);
    }
}
