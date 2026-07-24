package jline.api.pfqn;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Random;

import org.junit.jupiter.api.Test;

import jline.api.pfqn.ld.Pfqn_dac;
import jline.api.pfqn.ld.Pfqn_mvald;
import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * Tests for the DAC (Distribution Analysis by Chain) algorithm.
 *
 * The reference values are those published in the availability-model example of
 * E. de Souza e Silva, "Distribution Analysis of Product Form Queueing Networks",
 * UCLA CSD-870023, 1987, Section 3 and Appendix.
 */
public class Pfqn_dacTest {

    private static final double TOL = 1e-4;

    /** Availability model of Figure 1: 1 CPU, 2 active memory modules + 1 spare, 1 repairman. */
    private static Matrix demands() {
        Matrix L = new Matrix(3, 2);
        L.set(0, 0, 5.0);
        L.set(0, 1, 0.0);
        L.set(1, 0, 0.0);
        L.set(1, 1, 10.0);
        L.set(2, 0, 2.0);
        L.set(2, 1, 1.0);
        return L;
    }

    private static Matrix pop() {
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 1.0);
        N.set(0, 1, 3.0);
        return N;
    }

    /** Locate the row of the aggregate state (a,b,c) and return its probability. */
    private static double at(Ret.pfqnDAC d, int a, int b, int c) {
        for (int i = 0; i < d.states.getNumRows(); i++) {
            if ((int) d.states.get(i, 0) == a && (int) d.states.get(i, 1) == b
                    && (int) d.states.get(i, 2) == c) {
                return d.Pjoint.get(i, 0);
            }
        }
        throw new IllegalStateException("state not enumerated");
    }

    @Test
    public void testColdSpareJointDistribution() {
        // Center 2 is queue dependent: only 2 modules can fail, the spare cannot.
        Matrix mu = new Matrix(3, 4);
        for (int n = 0; n < 4; n++) {
            mu.set(0, n, 1.0);
            mu.set(1, n, n == 0 ? 1.0 : 2.0);
            mu.set(2, n, 1.0);
        }
        Ret.pfqnDAC d = Pfqn_dac.pfqn_dac(demands(), pop(), null, mu);

        // Joint distribution, Appendix step k=4.
        assertEquals(0.5381, at(d, 1, 3, 0), TOL);
        assertEquals(0.1076, at(d, 1, 2, 1), TOL);
        assertEquals(0.02152, at(d, 1, 1, 2), TOL);
        assertEquals(0.002152, at(d, 1, 0, 3), TOL);
        assertEquals(0.2153, at(d, 0, 3, 1), TOL);
        assertEquals(0.08609, at(d, 0, 2, 2), TOL);
        assertEquals(0.02582, at(d, 0, 1, 3), TOL);
        assertEquals(0.003442, at(d, 0, 0, 4), TOL);

        // States unreachable because a chain does not visit the center.
        assertEquals(0.0, at(d, 4, 0, 0), TOL);
        assertEquals(0.0, at(d, 2, 2, 0), TOL);
        assertEquals(0.0, at(d, 0, 4, 0), TOL);

        assertEquals(1.0, d.Pjoint.elementSum(), 1e-10);

        // Probability that the CPU is working, Section 5 example.
        assertEquals(0.6694, d.pi.get(0, 1), TOL);

        // Availability: at least one CPU and one memory module operational.
        double av = at(d, 1, 3, 0) + at(d, 1, 2, 1) + at(d, 1, 1, 2);
        assertEquals(0.6672, av, TOL);
    }

    @Test
    public void testColdSpareMeanMeasures() {
        Matrix mu = new Matrix(3, 4);
        for (int n = 0; n < 4; n++) {
            mu.set(0, n, 1.0);
            mu.set(1, n, n == 0 ? 1.0 : 2.0);
            mu.set(2, n, 1.0);
        }
        Ret.pfqnDAC d = Pfqn_dac.pfqn_dac(demands(), pop(), null, mu);

        // Appendix step k=4 reports the per-customer queue lengths of the memory
        // chain as L(2)=0.8983 and L(3)=0.1017. The chain holds 3 exchangeable
        // customers, so compare per customer, at the precision the paper prints.
        assertEquals(0.8983, d.Q.get(1, 1) / 3.0, TOL);
        assertEquals(0.1017, d.Q.get(2, 1) / 3.0, TOL);
        assertEquals(0.0, d.Q.get(0, 1), TOL);

        // The memory chain never visits the CPU failure center, and the CPU chain
        // holds a single customer split between its own center and the repairman.
        assertEquals(1.0, d.Q.get(0, 0) + d.Q.get(2, 0), 1e-10);
    }

    @Test
    public void testHotStandbyJointDistribution() {
        // The spare is powered and can also fail, so mu_2(3) grows accordingly.
        Matrix mu = new Matrix(3, 4);
        for (int n = 0; n < 4; n++) {
            mu.set(0, n, 1.0);
            mu.set(2, n, 1.0);
        }
        mu.set(1, 0, 1.0);
        mu.set(1, 1, 2.0);
        mu.set(1, 2, 2.5);
        mu.set(1, 3, 2.5);
        Ret.pfqnDAC d = Pfqn_dac.pfqn_dac(demands(), pop(), null, mu);

        assertEquals(0.5068, at(d, 1, 3, 0), TOL);
        assertEquals(0.1267, at(d, 1, 2, 1), TOL);
        assertEquals(0.02535, at(d, 1, 1, 2), TOL);
        assertEquals(0.002536, at(d, 1, 0, 3), TOL);
        assertEquals(0.2027, at(d, 0, 3, 1), TOL);
        assertEquals(0.1014, at(d, 0, 2, 2), TOL);
        assertEquals(0.03041, at(d, 0, 1, 3), TOL);
        assertEquals(0.004054, at(d, 0, 0, 4), TOL);
        assertEquals(1.0, d.Pjoint.elementSum(), 1e-10);
    }

    /** DAC must reproduce the mean measures of load-dependent MVA exactly. */
    @Test
    public void testAgreesWithMvaldOnLoadDependentNetworks() {
        Random rng = new Random(7);
        for (int trial = 0; trial < 5; trial++) {
            int J = 3;
            int R = 2;
            Matrix L = new Matrix(J, R);
            for (int j = 0; j < J; j++) {
                for (int r = 0; r < R; r++) {
                    L.set(j, r, 0.5 + 3.5 * rng.nextDouble());
                }
            }
            Matrix N = new Matrix(1, R);
            N.set(0, 0, 2.0);
            N.set(0, 1, 3.0);
            int K = 5;
            Matrix mu = new Matrix(J, K);
            for (int n = 0; n < K; n++) {
                mu.set(0, n, Math.min(n + 1, 2)); // two-server station
                mu.set(1, n, 1.0);                // single-server fixed rate
                mu.set(2, n, n + 1);              // infinite server
            }
            Matrix Z = Matrix.zeros(1, R);

            Ret.pfqnDAC d = Pfqn_dac.pfqn_dac(L, N, Z, mu);
            Ret.pfqnMVALD m = Pfqn_mvald.pfqn_mvald(L, N, Z, mu);

            assertEquals(1.0, d.Pjoint.elementSum(), 1e-10);
            for (int r = 0; r < R; r++) {
                assertEquals(m.X.get(r), d.X.get(0, r), 1e-8);
                for (int j = 0; j < J; j++) {
                    assertEquals(m.Q.get(j, r), d.Q.get(j, r), 1e-8);
                }
            }
            for (int j = 0; j < J; j++) {
                for (int n = 0; n <= K; n++) {
                    assertEquals(m.pi.get(j, n), d.pi.get(j, n), 1e-8);
                }
            }
        }
    }

    /** A non-zero think time is modelled as an appended infinite-server center. */
    @Test
    public void testThinkTimeAppendsDelayStation() {
        Matrix L = new Matrix(2, 2);
        L.set(0, 0, 1.0);
        L.set(0, 1, 2.0);
        L.set(1, 0, 3.0);
        L.set(1, 1, 1.0);
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 2.0);
        N.set(0, 1, 2.0);
        Matrix Z = new Matrix(1, 2);
        Z.set(0, 0, 4.0);
        Z.set(0, 1, 5.0);

        Ret.pfqnDAC d = Pfqn_dac.pfqn_dac(L, N, Z, null);
        assertEquals(3, d.states.getNumCols());
        assertEquals(1.0, d.Pjoint.elementSum(), 1e-10);

        // Equivalent model with the delay written out as an explicit IS station.
        Matrix Lx = new Matrix(3, 2);
        Lx.set(0, 0, 1.0);
        Lx.set(0, 1, 2.0);
        Lx.set(1, 0, 3.0);
        Lx.set(1, 1, 1.0);
        Lx.set(2, 0, 4.0);
        Lx.set(2, 1, 5.0);
        Matrix mux = new Matrix(3, 4);
        for (int n = 0; n < 4; n++) {
            mux.set(0, n, 1.0);
            mux.set(1, n, 1.0);
            mux.set(2, n, n + 1);
        }
        Ret.pfqnMVALD m = Pfqn_mvald.pfqn_mvald(Lx, N, Matrix.zeros(1, 2), mux);
        for (int r = 0; r < 2; r++) {
            assertEquals(m.X.get(r), d.X.get(0, r), 1e-8);
            for (int j = 0; j < 2; j++) {
                assertEquals(m.Q.get(j, r), d.Q.get(j, r), 1e-8);
            }
        }
    }

    /** Every state must carry non-negative mass: the recursion is numerically stable. */
    @Test
    public void testProbabilitiesAreNonNegative() {
        Matrix L = new Matrix(2, 1);
        L.set(0, 0, 1.0);
        L.set(1, 0, 1e-3);
        Matrix N = new Matrix(1, 1);
        N.set(0, 0, 8.0);
        Ret.pfqnDAC d = Pfqn_dac.pfqn_dac(L, N);
        for (int i = 0; i < d.Pjoint.getNumRows(); i++) {
            assertTrue(d.Pjoint.get(i, 0) >= 0, "negative probability at state " + i);
        }
        assertEquals(1.0, d.Pjoint.elementSum(), 1e-10);
    }
}
