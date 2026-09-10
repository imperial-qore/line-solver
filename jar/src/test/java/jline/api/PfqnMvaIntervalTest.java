package jline.api;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import jline.api.pfqn.mva.Pfqn_mva;
import jline.api.pfqn.mva.Pfqn_mva_interval;
import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * Numerical validation of pfqn_mva_interval against the worked figures of Luthi and Haring,
 * "Mean value analysis for queueing network models with intervals as input parameters",
 * Performance Evaluation 32(3):185-215, 1998, and against a brute-force sweep of the box.
 */
public class PfqnMvaIntervalTest {

    private static Matrix box(double[][] rows) {
        Matrix m = new Matrix(rows.length, 2);
        for (int i = 0; i < rows.length; i++) {
            m.set(i, 0, rows[i][0]);
            m.set(i, 1, rows[i][1]);
        }
        return m;
    }

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
        }
        return m;
    }

    private static Ret.pfqnMVA mva(double dcpu, double n, double z) {
        Matrix L = new Matrix(2, 1);
        L.set(0, 0, dcpu);
        L.set(1, 0, 10.0);
        return Pfqn_mva.pfqn_mva(L, row(n), row(z), null);
    }

    /** Section 3.2: D_cpu = [12,16], D_disk = 10, Z = [15,20], n = 10 gives X in [0.0620, 0.0799]. */
    @Test
    public void throughputIntervalMatchesThePaper() {
        Pfqn_mva_interval.Result res = Pfqn_mva_interval.pfqn_mva_interval(
                box(new double[][]{{12, 16}, {10, 10}}), row(10), row(15, 20));
        assertEquals(0.0620, res.X.get(0, 0), 5e-4);
        assertEquals(0.0799, res.X.get(0, 1), 5e-4);
    }

    /** The hull must contain every interior point of the box and must be attained at its corners. */
    @Test
    public void hullContainsTheGridAndIsAttained() {
        Pfqn_mva_interval.Result res = Pfqn_mva_interval.pfqn_mva_interval(
                box(new double[][]{{12, 16}, {10, 10}}), row(10), row(15, 20));

        double xlo = Double.POSITIVE_INFINITY, xup = Double.NEGATIVE_INFINITY;
        double[] qlo = {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        double[] qup = {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
        double[] rlo = {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        double[] rup = {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
        double rtlo = Double.POSITIVE_INFINITY, rtup = Double.NEGATIVE_INFINITY;
        for (int a = 0; a <= 8; a++) {
            double dcpu = 12.0 + a * 0.5;
            for (int b = 0; b <= 5; b++) {
                double z = 15.0 + b;
                Ret.pfqnMVA point = mva(dcpu, 10, z);
                double x = point.X.get(0);
                xlo = Math.min(xlo, x);
                xup = Math.max(xup, x);
                double rt = 0.0;
                for (int i = 0; i < 2; i++) {
                    qlo[i] = Math.min(qlo[i], point.Q.get(i));
                    qup[i] = Math.max(qup[i], point.Q.get(i));
                    rlo[i] = Math.min(rlo[i], point.R.get(i));
                    rup[i] = Math.max(rup[i], point.R.get(i));
                    rt += point.R.get(i);
                }
                rtlo = Math.min(rtlo, rt);
                rtup = Math.max(rtup, rt);
            }
        }

        assertEquals(xlo, res.X.get(0, 0), 1e-9);
        assertEquals(xup, res.X.get(0, 1), 1e-9);
        for (int i = 0; i < 2; i++) {
            assertEquals(qlo[i], res.Q.get(i, 0), 1e-9);
            assertEquals(qup[i], res.Q.get(i, 1), 1e-9);
            assertEquals(rlo[i], res.R.get(i, 0), 1e-9);
            assertEquals(rup[i], res.R.get(i, 1), 1e-9);
        }
        assertTrue(rtlo >= res.Rtot.get(0, 0) - 1e-9 && rtup <= res.Rtot.get(0, 1) + 1e-9);
        assertEquals(rtlo, res.Rtot.get(0, 0), 1e-9);
        assertEquals(rtup, res.Rtot.get(0, 1), 1e-9);
    }

    /** A thin box must reproduce ordinary MVA on both endpoints. */
    @Test
    public void thinBoxDegeneratesToExactMva() {
        Pfqn_mva_interval.Result res = Pfqn_mva_interval.pfqn_mva_interval(
                box(new double[][]{{14, 14}, {10, 10}}), row(10), row(15, 15));
        Ret.pfqnMVA exact = mva(14.0, 10, 15.0);
        assertEquals(exact.X.get(0), res.X.get(0, 0), 1e-12);
        assertEquals(exact.X.get(0), res.X.get(0, 1), 1e-12);
        for (int i = 0; i < 2; i++) {
            assertEquals(exact.Q.get(i), res.Q.get(i, 0), 1e-12);
            assertEquals(exact.Q.get(i), res.Q.get(i, 1), 1e-12);
            assertEquals(exact.R.get(i), res.R.get(i, 0), 1e-12);
            assertEquals(exact.R.get(i), res.R.get(i, 1), 1e-12);
        }
    }

    /** The population interval widens the throughput on the lower side only. */
    @Test
    public void populationIntervalWidensTheCorrectSide() {
        Matrix L = box(new double[][]{{12, 16}, {10, 10}});
        Pfqn_mva_interval.Result wide = Pfqn_mva_interval.pfqn_mva_interval(L, row(5, 10), row(15, 20));
        Pfqn_mva_interval.Result thin = Pfqn_mva_interval.pfqn_mva_interval(L, row(10), row(15, 20));
        assertEquals(mva(16.0, 5, 20.0).X.get(0), wide.X.get(0, 0), 1e-12);
        assertEquals(thin.X.get(0, 1), wide.X.get(0, 1), 1e-12);
        assertTrue(wide.Rtot.get(0, 0) <= thin.Rtot.get(0, 0) + 1e-12);
    }
}
