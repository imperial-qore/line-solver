package jline.solvers.mva;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.mva.SolverMVA;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * The joint-dependence evaluation point of solver_amvald: eta_k is read at
 * (Q_1, ..., 1 + delta_r Q_r, ..., Q_R), only the arriving class shifted.
 *
 * <p>Model: the IS+2xOI note model of Casale/Comte/Dorsman, closed under the
 * QD-AMVA closure -- an IS delay with think rates (0.7, 1.1) and two identical
 * PS stations carrying the support-rank eta of a three-server compatibility
 * station, mu = (2, 1, 2) with server 1 dedicated to class 1, server 2 shared
 * and server 3 dedicated to class 2. At N = (1, 1) the values are eta = 3 with
 * one class present and eta = 5 with both, bilinearly interpolated in between.
 *
 * <p>The expected queue lengths are python's, which has always evaluated eta at
 * this point. Incrementing every coordinate instead -- what the JAR and MATLAB
 * did before -- forces full support on a support-rank eta and moves these
 * numbers well outside the tolerance below.
 */
public class JdEvalPointNoteTest {

    private static final double[] MU = {2.0, 2.0};   // eta with exactly one class present
    private static final double MU_BOTH = 5.0;       // eta with both classes present

    /** svc(n) = 1/mu(supp n) on the unit box, bilinearly interpolated, returned as a rate. */
    private static final SerializableFunction<Matrix, Matrix> ETA = (Matrix n) -> {
        double n1 = Math.min(Math.max(n.get(0), 0.0), 1.0);
        double n2 = Math.min(Math.max(n.get(1), 0.0), 1.0);
        double s00 = 0.0, s10 = 1.0 / 3.0, s01 = 1.0 / 3.0, s11 = 1.0 / MU_BOTH;
        double svc = (1 - n1) * (1 - n2) * s00 + n1 * (1 - n2) * s10
                   + (1 - n1) * n2 * s01 + n1 * n2 * s11;
        Matrix out = new Matrix(1, 1);
        out.set(0, 0, svc > 0.0 ? 1.0 / svc : 1.0);
        return out;
    };

    private static Matrix peak(double v) {
        Matrix m = new Matrix(1, 2);
        m.set(0, 0, v);
        m.set(0, 1, v);
        return m;
    }

    @Test
    public void testNoteModelMatchesPython() {
        Network model = new Network("IS+2xOI");
        Delay is = new Delay(model, "IS");
        Queue cd1 = new Queue(model, "CD1", SchedStrategy.PS);
        Queue cd2 = new Queue(model, "CD2", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "C1", 1, is, 0);
        ClosedClass c2 = new ClosedClass(model, "C2", 1, is, 0);
        is.setService(c1, new Exp(0.7));
        is.setService(c2, new Exp(1.1));
        for (Queue q : new Queue[]{cd1, cd2}) {
            q.setService(c1, new Exp(1.0));
            q.setService(c2, new Exp(1.0));
            q.setNumberOfServers(1);
            q.setLimitedJointDependence(ETA, peak(MU_BOTH));
        }
        model.link(model.serialRouting(is, cd1, cd2));

        Matrix Q = new SolverMVA(model, "method", "qd").getAvgQLen();
        // python line_solver, SolverMVA(method='qd'), same model
        double[][] expected = {
                {0.6589337427, 0.5565433634},   // IS
                {0.1705331157, 0.2217281442},   // CD1
                {0.1705331157, 0.2217281442},   // CD2
        };
        for (int k = 0; k < 3; k++) {
            for (int r = 0; r < 2; r++) {
                assertEquals(expected[k][r], Q.get(k, r), 1e-5,
                        "QLen[" + k + "][" + r + "] disagrees with python");
            }
        }
        assertEquals(2.0, Q.elementSum(), 1e-5, "population is not conserved");
    }
}
