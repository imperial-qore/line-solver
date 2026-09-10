/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mc;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.function.UnaryOperator;

import org.apache.commons.math3.complex.Complex;
import org.junit.jupiter.api.Test;

import jline.api.lti.Laplace_invert;
import jline.api.lti.WeeksParams;
import jline.api.pfqn.nc.Pfqn_cyclet_ofree;
import jline.api.pfqn.nc.PfqnCycletResult;
import jline.util.matrix.Matrix;

/**
 * First passage times, the Laguerre inverter, and the overtake-free cycle time.
 *
 * <p>THE ORACLES SHARE NOTHING WITH THE FORMULAS. The passage moments are
 * checked against the dense closed form n! alpha (-S)^-n 1, which the
 * implementation deliberately does NOT use; the semi-Markov recursion against
 * the Markov one on a chain that is both; the cycle time against the moments
 * Harrison and Knottenbelt 2002 publishes for four models. The MATLAB, native
 * Python and C++ twins report the same numbers.
 */
public class CtmcPassageTest {

    /** M/M/1/K generator, K+1 states. */
    private static Matrix mm1k(int K, double lambda, double mu) {
        int n = K + 1;
        Matrix Q = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            if (i + 1 < n) {
                Q.set(i, i + 1, lambda);
            }
            if (i > 0) {
                Q.set(i, i - 1, mu);
            }
        }
        for (int i = 0; i < n; i++) {
            double r = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    r += Q.get(i, j);
                }
            }
            Q.set(i, i, -r);
        }
        return Q;
    }

    private static Matrix unit(int n, int i) {
        Matrix p = new Matrix(1, n);
        p.set(0, i, 1.0);
        return p;
    }

    @Test
    public void passageLawMatchesTheAbsorptionLawBuiltIndependently() {
        int K = 6;
        int n = K + 1;
        Matrix Q = mm1k(K, 1.0, 1.5);
        int[] target = new int[]{n - 1};
        double[] tset = new double[81];
        for (int i = 0; i < tset.length; i++) {
            tset[i] = 0.25 * i;
        }
        PassageCurve c = Ctmc_passage_time.ctmc_passage_time(Q, unit(n, 0), target, tset);
        assertEquals(0.0, c.atom, 1e-14);

        // Oracle: assemble S by hand and exponentiate at each point.
        Matrix S = new Matrix(n - 1, n - 1);
        for (int i = 0; i + 1 < n; i++) {
            for (int j = 0; j + 1 < n; j++) {
                S.set(i, j, Q.get(i, j));
            }
        }
        double[] s0 = new double[n - 1];
        for (int i = 0; i + 1 < n; i++) {
            double r = 0.0;
            for (int j = 0; j + 1 < n; j++) {
                r += S.get(i, j);
            }
            s0[i] = -r;
        }
        for (int it = 0; it < tset.length; it++) {
            Matrix E = S.scale(tset[it]).expm();
            double sF = 0.0;
            double sf = 0.0;
            for (int j = 0; j + 1 < n; j++) {
                double a = E.get(0, j);
                sF += a;
                sf += a * s0[j];
            }
            assertEquals(1.0 - sF, c.F[it], 1e-10);
            assertEquals(sf, c.f[it], 1e-10);
        }
    }

    @Test
    public void eq3ReproducesTheDenseClosedFormItDoesNotUse() {
        int K = 6;
        int n = K + 1;
        Matrix Q = mm1k(K, 1.0, 1.5);
        PassageMomentsResult pm =
                Ctmc_passage_moments.ctmc_passage_moments(Q, unit(n, 0), new int[]{n - 1}, 4);

        int nA = n - 1;
        Matrix A = new Matrix(nA, nA);
        for (int i = 0; i < nA; i++) {
            for (int j = 0; j < nA; j++) {
                A.set(i, j, -Q.get(i, j));
            }
        }
        Matrix Ainv = new Matrix(nA, nA);
        for (int c = 0; c < nA; c++) {
            Matrix e = new Matrix(nA, 1);
            e.set(c, 0, 1.0);
            Matrix col = new Matrix(nA, 1);
            Matrix.solve(A, e, col);
            for (int r = 0; r < nA; r++) {
                Ainv.set(r, c, col.get(r, 0));
            }
        }
        Matrix P = new Matrix(nA, nA);
        for (int i = 0; i < nA; i++) {
            P.set(i, i, 1.0);
        }
        double fact = 1.0;
        for (int k = 1; k <= 4; k++) {
            P = P.mult(Ainv);
            fact *= k;
            double m = 0.0;
            for (int j = 0; j < nA; j++) {
                m += P.get(0, j);
            }
            assertEquals(fact * m, pm.m[k - 1], Math.abs(fact * m) * 1e-9);
        }
        // The value the MATLAB, Python and C++ twins report for the mean.
        assertEquals(50.34375, pm.m[0], 1e-9);
    }

    @Test
    public void theAtomAtZeroIsReportedNotDropped() {
        int K = 6;
        int n = K + 1;
        Matrix Q = mm1k(K, 1.0, 1.5);
        Matrix pi0 = new Matrix(1, n);
        pi0.set(0, 0, 0.7);
        pi0.set(0, n - 1, 0.3);
        PassageCurve c = Ctmc_passage_time.ctmc_passage_time(Q, pi0, new int[]{n - 1},
                new double[]{0.0, 1.0, 5.0});
        assertEquals(0.3, c.atom, 1e-14);
        // F(0) IS the atom: a passage started inside the target completed at once.
        assertEquals(0.3, c.F[0], 1e-12);
    }

    @Test
    public void hittingTimeIsTheFirstMomentAndSaysInfinityWhenItMust() {
        int K = 6;
        int n = K + 1;
        Matrix Q = mm1k(K, 1.0, 1.5);
        Matrix h = Ctmc_hitting_time.ctmc_hitting_time(Q, new int[]{n - 1});
        assertEquals(50.34375, h.get(0, 0), 1e-9);

        // An isolated absorbing state cannot reach state 0.
        Matrix Q2 = new Matrix(n + 1, n + 1);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Q2.set(i, j, Q.get(i, j));
            }
        }
        Matrix h2 = Ctmc_hitting_time.ctmc_hitting_time(Q2, new int[]{0});
        assertTrue(Double.isInfinite(h2.get(n, 0)));
    }

    @Test
    public void theTransformRouteAgreesWithTheExponentialOne() {
        int K = 6;
        int n = K + 1;
        Matrix Q = mm1k(K, 1.0, 1.5);
        int[] target = new int[]{n - 1};
        double[] tset = new double[20];
        for (int i = 0; i < tset.length; i++) {
            tset[i] = 0.5 * (i + 1);
        }
        PassageCurve a = Ctmc_passage_time.ctmc_passage_time(Q, unit(n, 0), target, tset,
                "expm", "euler");
        PassageCurve b = Ctmc_passage_time.ctmc_passage_time(Q, unit(n, 0), target, tset,
                "lt", "euler");
        for (int i = 0; i < tset.length; i++) {
            assertEquals(a.F[i], b.F[i], 1e-7);
            assertEquals(a.f[i], b.f[i], 1e-7);
        }
        Complex L0 = Ctmc_passage_lst.ctmc_passage_lst(Q, unit(n, 0), target,
                new Complex(0.0, 0.0));
        assertEquals(1.0, L0.getReal(), 1e-12);
    }

    @Test
    public void theSemiMarkovRecursionReproducesTheMarkovOne() {
        int K = 6;
        int n = K + 1;
        int nmax = 4;
        Matrix Q = mm1k(K, 1.0, 1.5);
        int[] target = new int[]{n - 1};

        Matrix P = new Matrix(n, n);
        double[] rate = new double[n];
        for (int i = 0; i < n; i++) {
            rate[i] = -Q.get(i, i);
            if (rate[i] > 0.0) {
                for (int j = 0; j < n; j++) {
                    P.set(i, j, (j == i) ? 0.0 : Q.get(i, j) / rate[i]);
                }
            } else {
                P.set(i, i, 1.0);
            }
        }
        Matrix hmom = new Matrix(n, nmax);
        for (int i = 0; i < n; i++) {
            double f = 1.0;
            double p = 1.0;
            for (int r = 1; r <= nmax; r++) {
                f *= r;
                p *= rate[i];
                hmom.set(i, r - 1, rate[i] > 0.0 ? f / p : 0.0);
            }
        }

        PassageMomentsResult mC =
                Ctmc_passage_moments.ctmc_passage_moments(Q, unit(n, 0), target, nmax);
        PassageMomentsResult mS =
                Smp_passage_moments.smp_passage_moments(P, hmom, unit(n, 0), target, nmax);
        for (int q = 0; q < nmax; q++) {
            assertEquals(mC.m[q], mS.m[q], Math.abs(mC.m[q]) * 1e-9);
        }

        List<UnaryOperator<Complex>> hlst = new ArrayList<UnaryOperator<Complex>>();
        for (int i = 0; i < n; i++) {
            final double r = rate[i];
            hlst.add(new UnaryOperator<Complex>() {
                @Override
                public Complex apply(Complex s) {
                    return new Complex(r, 0.0).divide(s.add(new Complex(r, 0.0)));
                }
            });
        }
        Complex Ls = Smp_passage_lst.smp_passage_lst(P, hlst, unit(n, 0), target,
                new Complex(0.05, 0.0));
        Complex Lc = Ctmc_passage_lst.ctmc_passage_lst(Q, unit(n, 0), target,
                new Complex(0.05, 0.0));
        assertEquals(Lc.getReal(), Ls.getReal(), 1e-11);
    }

    @Test
    public void nonMarkovHoldingTimesGiveTheirExactMoments() {
        // 1 -> 2 with a deterministic sojourn of 2: the passage IS 2.
        Matrix P = new Matrix(2, 2);
        P.set(0, 1, 1.0);
        P.set(1, 1, 1.0);
        Matrix hmom = new Matrix(2, 3);
        hmom.set(0, 0, 2.0);
        hmom.set(0, 1, 4.0);
        hmom.set(0, 2, 8.0);
        PassageMomentsResult m =
                Smp_passage_moments.smp_passage_moments(P, hmom, unit(2, 0), new int[]{1}, 3);
        assertEquals(2.0, m.m[0], 1e-12);
        assertEquals(4.0, m.m[1], 1e-12);
        assertEquals(8.0, m.m[2], 1e-12);

        // 1 -> 2 -> 3 with Exp(2) sojourns: the passage 1 -> 3 is Erlang(2,2).
        Matrix P3 = new Matrix(3, 3);
        P3.set(0, 1, 1.0);
        P3.set(1, 2, 1.0);
        P3.set(2, 2, 1.0);
        Matrix h3 = new Matrix(3, 3);
        for (int i = 0; i < 2; i++) {
            h3.set(i, 0, 0.5);
            h3.set(i, 1, 0.5);
            h3.set(i, 2, 0.75);
        }
        PassageMomentsResult m3 =
                Smp_passage_moments.smp_passage_moments(P3, h3, unit(3, 0), new int[]{2}, 3);
        assertEquals(1.0, m3.m[0], 1e-12);
        assertEquals(1.5, m3.m[1], 1e-12);
        assertEquals(3.0, m3.m[2], 1e-12);
    }

    @Test
    public void weeksIsTheSharpestInverterOnARationalTransform() {
        // Measured worst absolute error on 2/(s+2) over t in {0.1,...,5}:
        // euler 2.5e-10, talbot 2.2e-11, gaver 1.2e-04, cme 4.2e-04, weeks 3.3e-14.
        UnaryOperator<Complex> F = new UnaryOperator<Complex>() {
            @Override
            public Complex apply(Complex s) {
                return new Complex(2.0, 0.0).divide(s.add(new Complex(2.0, 0.0)));
            }
        };
        WeeksParams w = Laplace_invert.laplace_weeks_scaling(F);
        assertEquals(0.0, w.sigma, 1e-15);
        assertEquals(1.0, w.b, 1e-15);
        double[] tt = new double[]{0.1, 0.5, 1.0, 2.0, 5.0};
        for (double t : tt) {
            assertEquals(2.0 * Math.exp(-2.0 * t), Laplace_invert.laplace_invert_weeks(w, t), 1e-11);
        }
        // The true coefficients are q_n = 0.8 * 0.6^n, which is what the
        // 1/(1-z) factor produces; the factor (1-z) printed in Eq. 10 starts
        // 0.8, -1.12, 0.128 and is wrong at every t.
        double p = 0.8;
        for (int n = 0; n < 6; n++) {
            assertEquals(p, w.q[n], 1e-10);
            p *= 0.6;
        }
    }

    @Test
    public void theWeeksScalingSearchRefusesANonSmoothDensityByName() {
        // A deterministic delay has transform exp(-s); its density is a point
        // mass and has no Laguerre representation. Fig. 1 must SAY SO rather
        // than return the last iterate, which would be noise as an answer.
        UnaryOperator<Complex> det = new UnaryOperator<Complex>() {
            @Override
            public Complex apply(Complex s) {
                return s.negate().exp();
            }
        };
        RuntimeException e = assertThrows(RuntimeException.class,
                () -> Laplace_invert.laplace_weeks_scaling(det));
        assertTrue(e.getMessage().contains("no suitable scaling"));
    }

    @Test
    public void theCycleTimeReproducesThePapersPublishedMoments() {
        double[] mu = new double[]{3.0, 5.0, 4.0, 6.0, 2.0, 1.0};
        double p12 = 0.2;
        double p13 = 0.5;
        double p14 = 0.3;
        double[] v = new double[]{1.0, p12, p13, p14, p12, p14};
        List<int[]> paths = new ArrayList<int[]>();
        paths.add(new int[]{0, 2});
        paths.add(new int[]{0, 1, 4});
        paths.add(new int[]{0, 3, 5});
        double[] tset = new double[161];
        for (int i = 0; i < tset.length; i++) {
            tset[i] = 0.25 * i;
        }
        PfqnCycletResult r = Pfqn_cyclet_ofree.pfqn_cyclet_ofree(v, mu, 18, paths, tset, "auto", 3,
                new double[]{p13, p12, p14}, "euler", Pfqn_cyclet_ofree.DEFAULT_TOL);
        // Sec. 7.2: "6.12717, 53.3067 and 612.887".
        assertEquals(6.12717, r.mom[0], 1e-5);
        assertEquals(53.3067, r.mom[1], 1e-4);
        assertEquals(612.887, r.mom[2], 1e-3);
        for (String m : r.method) {
            assertEquals("exact", m);
        }
        assertEquals(1.0, r.F[r.F.length - 1], 1e-5);
    }

    @Test
    public void theErlangAndBranchingErlangExamplesOfSec6() {
        // Sec. 6.1: a 3-stage Erlang with lambda = 2. Equal rates, so 'auto'
        // MUST leave Theorem 2 alone: its partial fractions divide by the rate
        // differences.
        List<int[]> one = new ArrayList<int[]>();
        one.add(new int[]{0, 1, 2});
        PfqnCycletResult e3 = Pfqn_cyclet_ofree.pfqn_cyclet_ofree(
                new double[]{1.0, 1.0, 1.0}, new double[]{2.0, 2.0, 2.0}, 1, one,
                new double[]{1.0});
        assertEquals(1.5, e3.mom[0], 1e-10);
        assertEquals(3.0, e3.mom[1], 1e-10);
        assertEquals(7.5, e3.mom[2], 1e-10);
        assertEquals("lt", e3.method.get(0));

        // Sec. 6.2: a 3-stage Erlang(1) and a 12-stage Erlang(2) branch, each
        // with probability 1/2. Moments 4.5, 25.5 and 166.5.
        double[] vb = new double[15];
        double[] mub = new double[15];
        for (int i = 0; i < 15; i++) {
            vb[i] = 0.5;
            mub[i] = (i < 3) ? 1.0 : 2.0;
        }
        List<int[]> two = new ArrayList<int[]>();
        two.add(new int[]{0, 1, 2});
        int[] up = new int[12];
        for (int i = 0; i < 12; i++) {
            up[i] = 3 + i;
        }
        two.add(up);
        PfqnCycletResult br = Pfqn_cyclet_ofree.pfqn_cyclet_ofree(vb, mub, 1, two,
                new double[]{1.0}, "auto", 3, new double[]{0.5, 0.5}, "euler",
                Pfqn_cyclet_ofree.DEFAULT_TOL);
        assertEquals(4.5, br.mom[0], 1e-10);
        assertEquals(25.5, br.mom[1], 1e-10);
        assertEquals(166.5, br.mom[2], 1e-10);

        // Sec. 6.2, the first-to-last passage of the upper branch: an Erlang-11.
        List<int[]> u11 = new ArrayList<int[]>();
        int[] p11 = new int[11];
        for (int i = 0; i < 11; i++) {
            p11[i] = 3 + i;
        }
        u11.add(p11);
        PfqnCycletResult r11 =
                Pfqn_cyclet_ofree.pfqn_cyclet_ofree(vb, mub, 1, u11, new double[]{1.0});
        assertEquals(5.5, r11.mom[0], 1e-10);
        assertEquals(33.0, r11.mom[1], 1e-10);
        assertEquals(214.5, r11.mom[2], 1e-10);
    }

    @Test
    public void theRefusalsAreByName() {
        Matrix Q = mm1k(3, 1.0, 1.5);
        Matrix pi0 = unit(4, 0);
        assertThrows(RuntimeException.class,
                () -> Ctmc_passage_ph.ctmc_passage_ph(Q, pi0, new int[]{}));
        assertThrows(RuntimeException.class,
                () -> Ctmc_passage_ph.ctmc_passage_ph(Q, pi0, new int[]{99}));
        // A matrix whose rows do not sum to zero is not a generator, and is
        // refused rather than silently repaired.
        Matrix bad = mm1k(3, 1.0, 1.5);
        bad.set(0, 0, bad.get(0, 0) + 1.0);
        assertThrows(RuntimeException.class,
                () -> Ctmc_passage_ph.ctmc_passage_ph(bad, pi0, new int[]{3}));
        assertThrows(RuntimeException.class, () -> Ctmc_passage_time.ctmc_passage_time(
                Q, pi0, new int[]{3}, new double[]{1.0}, "nope", "euler"));
        // Equal rates on the path are refused by the exact route by name, not
        // silently divided by zero.
        List<int[]> pp = new ArrayList<int[]>();
        pp.add(new int[]{0, 1});
        assertThrows(RuntimeException.class,
                () -> Pfqn_cyclet_ofree.pfqn_cyclet_ofree(new double[]{1.0, 1.0},
                        new double[]{2.0, 2.0}, 2, pp, new double[]{1.0}, "exact", 3, null,
                        "euler", Pfqn_cyclet_ofree.DEFAULT_TOL));
    }
}
