package jline.api.pfqn.nc;

import jline.io.Ret;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The Lam scaling of pfqn_ca must be read off a LOWER bound on log G that a
 * MIXED state can reach, not off one station holding every class at once.
 *
 * An LQN layer reaches SolverNC as a think-time class beside a zero-Z call
 * class, and the per-configuration estimate discarded the delay entirely as
 * soon as one class had no think time. On L=[1e-9,1], N=[99,1], Z=[1,0] -- the
 * layer lqn_twotasks hands the NC solver -- it returned the all-at-the-queue
 * -2051.6 against a true log G of -359.134, so kscale was -30 and the scaling
 * went the WRONG WAY: Z/2^-30 = 1.07e9 made the delay column Z^n/n! peak at
 * e^1699 and overflow at n=[40,0]. This class returned G = Inf, lG = Inf there
 * (python raised OverflowError on the same input, which is how it was found).
 */
public class Pfqn_caScalingTest {

    /** 120-digit mpmath convolution of the same recursion. */
    private static final double LG_EXACT = -359.13420517157538927;

    @Test
    public void testAZeroThinkTimeClassDoesNotDiscardTheDelay() {
        Matrix L = new Matrix(1, 2);
        L.set(0, 0, 1e-9);
        L.set(0, 1, 1.0);
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 99.0);
        N.set(0, 1, 1.0);
        Matrix Z = new Matrix(1, 2);
        Z.set(0, 0, 1.0);
        Z.set(0, 1, 0.0);

        Ret.pfqnNc r = Pfqn_ca.pfqn_ca(L, N, Z);
        assertTrue(Double.isFinite(r.lG), "lG must stay finite: the delay column overflowed");
        assertEquals(LG_EXACT, r.lG, 1e-9, "lG must match the exact convolution");
        assertTrue(r.G > 0 && Double.isFinite(r.G), "G must stay in range");
    }

    @Test
    public void testTheEstimateStaysALowerBoundOnLogG() {
        // A scaling read off an estimate ABOVE log G would drive the recursion
        // to underflow, the other half of the range problem Reiser reports.
        Matrix L = new Matrix(2, 2);
        L.set(0, 0, 1e-6);
        L.set(0, 1, 3.0);
        L.set(1, 0, 2.0);
        L.set(1, 1, 1e-4);
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 5.0);
        N.set(0, 1, 4.0);
        Matrix Z = new Matrix(1, 2);
        Z.set(0, 0, 7.0);
        Z.set(0, 1, 0.0);

        Ret.pfqnNc r = Pfqn_ca.pfqn_ca(L, N, Z);
        double best = 0.0;
        for (int cls = 0; cls < 2; cls++) {
            double perClass = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < 2; i++) {
                perClass = Math.max(perClass, N.get(cls) * Math.log(L.get(i, cls)));
            }
            if (Z.get(cls) > 0) {
                perClass = Math.max(perClass,
                    N.get(cls) * Math.log(Z.get(cls)) - jline.util.Maths.factln(N.get(cls)));
            }
            best += perClass;
        }
        assertTrue(best <= r.lG + 1e-9, "the per-class bound must not exceed log G");
    }
}
