package jline.api;

import jline.api.sum.Sum_closed;
import jline.api.sum.Sum_closing;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Tests for the summation method (SUM/ESUM) and the closing method,
 * validated against the worked examples in Bolch, Greiner, de Meer,
 * Trivedi, "Queueing Networks and Markov Chains", 2nd ed., Wiley, 2006.
 */
public class SumMethodTest {

    private static final double INF = Double.POSITIVE_INFINITY;

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
        }
        return m;
    }

    private static Matrix col(double... v) {
        Matrix m = new Matrix(v.length, 1);
        for (int i = 0; i < v.length; i++) {
            m.set(i, 0, v[i]);
        }
        return m;
    }

    @Test
    public void testExample95ProductFormSum() {
        // Bolch Example 9.5: closed PF network, N=4 nodes, K=3, node 1 is -/M/2
        Sum_closed.Result res = Sum_closed.sum_closed(
                col(0.5, 0.3, 0.4, 1.0), row(3), row(0),
                col(2, 1, 1, INF), col(1, 1, 1, 1), 1e-3, 10000);
        assertEquals(1.193, res.XN.get(0), 1e-3);
        assertEquals(0.637, res.QN.get(0, 0), 1e-3);
        assertEquals(0.470, res.QN.get(1, 0), 1e-3);
        assertEquals(0.700, res.QN.get(2, 0), 1e-3);
        assertEquals(1.193, res.QN.get(3, 0), 1e-3);
    }

    @Test
    public void testExample1010Esum() {
        // Bolch Example 10.10: closed NPF network, N=5 nodes, K=17 (Table 10.15)
        Sum_closed.Result res = Sum_closed.sum_closed(
                col(1 / 13.5, 0.2 / 1.15, 0.4 / 1.2, 0.3 / 1.2, 0.1 / 1.7),
                row(17), row(0), col(1, 3, INF, 4, 1),
                col(1.0, 0.8, 1.0, 3.0, 1.6), 1e-6, 10000);
        assertEquals(11.34, res.XN.get(0), 5e-3);
        assertEquals(1.50, 17 / res.XN.get(0), 5e-3);
    }

    @Test
    public void testExample1013ClosingMixed() {
        // Bolch Example 10.13: mixed NPF network, class 1 closed (K1=9),
        // class 2 open (lambda=5, ca2=0.7), closed with K2=500 (Table 10.21)
        double[] e1 = {1, 1.428571, 1.428571, 0.428571, 0.428571};
        double[] e2 = {1, 1.25, 1.428571, 0.25, 0.428571};
        double[] mu = {4, 4, 6, 4, 5};
        double[] c2 = {0.4, 0.3, 0.3, 0.4, 0.5};
        Matrix L = new Matrix(5, 2);
        Matrix scv = new Matrix(5, 2);
        for (int i = 0; i < 5; i++) {
            L.set(i, 0, e1[i] / mu[i]);
            L.set(i, 1, e2[i] / mu[i]);
            scv.set(i, 0, c2[i]);
            scv.set(i, 1, c2[i]);
        }
        Sum_closing.Result res = Sum_closing.sum_closing(
                row(0, 5), row(1, 0.7), L, col(3, 4, 3, 2, 2), scv,
                row(9, INF), row(0, 0), 500, 1e-6, 10000);
        assertEquals(5.0, res.XN.get(1), 1e-2);
        double[] rhoBook = {0.84, 0.85, 0.80, 0.43, 0.43};
        double[] kBook = {5.2, 5.8, 4.1, 1.0, 1.0};
        for (int i = 0; i < 5; i++) {
            assertEquals(rhoBook[i], res.UN.get(i, 0) + res.UN.get(i, 1), 5e-3);
            assertEquals(kBook[i], res.QN.get(i, 0) + res.QN.get(i, 1), 5e-2);
        }
    }

    @Test
    public void testClosingConvergesToArrivalRate() {
        // closing method: open-class throughput approaches lambda0 from below
        Matrix L = col(1 / 0.9 / 9, 3.0 / 7 / 0.9 / 10, 7.0 / 9 / 0.9 / 12, 1 / 0.9 / 4);
        Sum_closing.Result res = Sum_closing.sum_closing(
                row(3), row(1.5), L, col(1, 1, 1, 1), col(0.5, 0.8, 2.4, 4.0),
                row(INF), row(0), 5000, 1e-6, 10000);
        assertEquals(3.0, res.XN.get(0), 1e-3);
    }
}
