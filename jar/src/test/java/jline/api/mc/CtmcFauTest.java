package jline.api.mc;

import org.junit.jupiter.api.Test;

import jline.util.matrix.Matrix;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Fast adaptive uniformization against ANALYTICAL and cross-method oracles.
 *
 * The two-state chain has a closed-form transient, which is the primary oracle. The
 * queueing cases are checked against {@link Ctmc_foxglynn} at a tolerance three orders
 * tighter, which is a genuinely different algorithm on the same generator: full-space
 * uniformization at max_i |q_ii| against the adaptive rate sequence. Two structural
 * properties are asserted alongside the numbers, because they are what the method
 * promises and what a plausible-looking wrong answer would break: the result is a
 * componentwise lower bound on the exact distribution, and its missing mass IS its L1
 * error.
 *
 * Twin of python/tests/test_ctmc_fau.py and of cpp/tests/test_ctmc_fau.cpp.
 */
public class CtmcFauTest {

    private static Matrix mm1kGenerator(double lambda, double mu, int K) {
        int n = K + 1;
        Matrix Q = new Matrix(n, n);
        for (int i = 0; i + 1 < n; i++) {
            Q.set(i, i + 1, lambda);
            Q.set(i + 1, i, mu);
        }
        for (int i = 0; i < n; i++) {
            double off = 0.0;
            if (i + 1 < n) off += lambda;
            if (i > 0) off += mu;
            Q.set(i, i, -off);
        }
        return Q;
    }

    private static Matrix unitVector(int n, int i) {
        Matrix pi0 = new Matrix(1, n);
        pi0.set(0, i, 1.0);
        return pi0;
    }

    private static double l1(Matrix a, Matrix b) {
        double s = 0.0;
        for (int i = 0; i < a.getNumElements(); i++) {
            s += Math.abs(a.get(i) - b.get(i));
        }
        return s;
    }

    @Test
    public void twoStateChainMatchesTheClosedFormTransient() {
        double a = 3.0;
        double b = 1.5;
        double t = 0.7;
        Matrix Q = new Matrix(2, 2);
        Q.set(0, 0, -a);
        Q.set(0, 1, a);
        Q.set(1, 0, b);
        Q.set(1, 1, -b);
        Matrix pi0 = unitVector(2, 0);

        Ctmc_fau.CtmcFauResult r = Ctmc_fau.ctmc_fau(pi0, Q, t);

        double decay = Math.exp(-(a + b) * t);
        double p0 = b / (a + b) + a / (a + b) * decay;
        double p1 = a / (a + b) * (1.0 - decay);
        assertEquals(p0, r.pit.get(0), 1e-6);
        assertEquals(p1, r.pit.get(1), 1e-6);
        assertTrue(r.steps > 1);
    }

    @Test
    public void theErrorIsTheMissingMass() {
        Matrix Q = mm1kGenerator(1.0, 2.0, 30);
        Matrix pi0 = unitVector(31, 0);
        Ctmc_fau.CtmcFauResult r = Ctmc_fau.ctmc_fau(pi0, Q, 5.0);
        Matrix reference = Ctmc_foxglynn.ctmc_foxglynn(pi0, Q, 5.0, 1e-14, -1);

        // The three approximations only remove mass, so the defect IS the error.
        assertEquals(l1(r.pit, reference), r.errorBound, 1e-9);
        assertTrue(r.errorBound <= 1e-5, "errorBound " + r.errorBound);
    }

    @Test
    public void resultIsAComponentwiseLowerBound() {
        Matrix Q = mm1kGenerator(1.0, 2.0, 30);
        Matrix pi0 = unitVector(31, 0);
        Ctmc_fau.CtmcFauResult r = Ctmc_fau.ctmc_fau(pi0, Q, 5.0);
        Matrix reference = Ctmc_foxglynn.ctmc_foxglynn(pi0, Q, 5.0, 1e-14, -1);
        for (int i = 0; i < 31; i++) {
            assertTrue(r.pit.get(i) <= reference.get(i) + 1e-12, "state " + i);
        }
    }

    @Test
    public void tighteningTheToleranceTightensTheError() {
        Matrix Q = mm1kGenerator(1.0, 2.0, 20);
        Matrix pi0 = unitVector(21, 0);
        Ctmc_fau.CtmcFauResult loose = Ctmc_fau.ctmc_fau(pi0, Q, 4.0, 1e-4, 1e-12, -1);
        Ctmc_fau.CtmcFauResult tight = Ctmc_fau.ctmc_fau(pi0, Q, 4.0, 1e-10, 1e-12, -1);
        assertTrue(tight.errorBound < loose.errorBound);
        assertTrue(tight.steps > loose.steps);
    }

    @Test
    public void fastStatesCarryingNoMassDoNotSetTheCost() {
        // Six slow states, then four states of rate 1e6 that the initial distribution
        // cannot reach within the horizon. Ordinary uniformization pays for the fast
        // ones, adaptive uniformization does not.
        int ns = 6;
        int n = ns + 4;
        double t = 1.0;
        Matrix Q = new Matrix(n, n);
        for (int i = 0; i + 1 < ns; i++) {
            Q.set(i, i + 1, 0.5);
            Q.set(i + 1, i, 0.4);
        }
        for (int i = ns; i + 1 < n; i++) {
            Q.set(i, i + 1, 1e6);
            Q.set(i + 1, i, 1e6);
        }
        Q.set(n - 1, ns, 1e6);
        for (int i = 0; i < n; i++) {
            double off = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) off += Q.get(i, j);
            }
            Q.set(i, i, -off);
        }
        Matrix pi0 = unitVector(n, 0);

        Ctmc_fau.CtmcFauResult r = Ctmc_fau.ctmc_fau(pi0, Q, t);
        Matrix reference = Ctmc_foxglynn.ctmc_foxglynn(pi0, Q, t, 1e-14, -1);

        assertTrue(l1(r.pit, reference) < 1e-6, "l1 " + l1(r.pit, reference));
        assertTrue(r.lambdaMax <= 1.0, "lambdaMax " + r.lambdaMax);
        assertTrue(r.uniformRate >= 1e6, "uniformRate " + r.uniformRate);
        // The step count follows the visited rate, not the global one.
        assertTrue(r.steps < 50, "steps " + r.steps);
    }

    @Test
    public void occupancyThresholdBoundsTheSupport() {
        int K = 400;
        Matrix Q = mm1kGenerator(1.0, 3.0, K);
        Matrix pi0 = unitVector(K + 1, 0);
        Ctmc_fau.CtmcFauResult r = Ctmc_fau.ctmc_fau(pi0, Q, 4.0, 1e-8, 1e-10, -1);
        Matrix reference = Ctmc_foxglynn.ctmc_foxglynn(pi0, Q, 4.0, 1e-14, -1);

        assertTrue(r.supportMax < K + 1, "supportMax " + r.supportMax);
        assertTrue(l1(r.pit, reference) <= r.errorBound + 1e-12);
        assertTrue(l1(r.pit, reference) < 1e-7);
    }

    @Test
    public void absorbingChainTerminates() {
        Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, -2.0);
        Q.set(0, 1, 2.0);
        Q.set(1, 1, -1.0);
        Q.set(1, 2, 1.0);
        Matrix pi0 = unitVector(3, 0);
        Ctmc_fau.CtmcFauResult r = Ctmc_fau.ctmc_fau(pi0, Q, 3.0);
        Matrix reference = Ctmc_foxglynn.ctmc_foxglynn(pi0, Q, 3.0, 1e-14, -1);
        assertTrue(l1(r.pit, reference) < 1e-6);
        assertTrue(r.absorbed);
        assertEquals(3, r.steps);
    }

    @Test
    public void zeroHorizonReturnsTheInitialDistribution() {
        Matrix Q = mm1kGenerator(1.0, 2.0, 5);
        Matrix pi0 = unitVector(6, 2);
        Ctmc_fau.CtmcFauResult r = Ctmc_fau.ctmc_fau(pi0, Q, 0.0);
        for (int i = 0; i < 6; i++) {
            assertEquals(pi0.get(i), r.pit.get(i), 0.0);
        }
        assertEquals(0.0, r.errorBound, 0.0);
    }
}
