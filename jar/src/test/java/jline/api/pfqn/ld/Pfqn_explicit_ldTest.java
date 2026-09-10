/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.ld;

import jline.api.pfqn.nc.Pfqn_explicit;
import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validation of {@link Pfqn_explicit_ld}, i.e. Theorem 1 of Casale, Harrison and Ong
 * (Perform. Eval. 2021) carried over the divided-difference form of Casale (SIGMETRICS
 * 2017), Corollary 3.2.
 *
 * <p>The closed form is EXACT, so it is checked against oracles rather than for
 * self-consistency:</p>
 * <ol>
 *   <li>{@link Pfqn_gld}, the load-dependent convolution, on multi-server, fixed-rate,
 *       arbitrary limited load-dependent and never-settling rate lattices.</li>
 *   <li>{@link Pfqn_explicit} itself on the fixed-rate degeneration mu = 1, which the
 *       load-dependent route must reproduce to the bit.</li>
 *   <li>The dispatch: the {@code divdiff} method name has to reach here through
 *       {@link Pfqn_ncld} and agree with its exact route.</li>
 *   <li>The near-tie contract: an induced tie that lands two ulps apart is MISSED at
 *       tol = eps, and the loss report has to say so; any looser tolerance merges the
 *       pair and Eq. (16) is exact.</li>
 * </ol>
 */
public class Pfqn_explicit_ldTest {

    private static final double TOL = 1e-8;

    private static Matrix mat(double[][] v) {
        Matrix m = new Matrix(v.length, v[0].length);
        for (int i = 0; i < v.length; i++) {
            for (int j = 0; j < v[0].length; j++) {
                m.set(i, j, v[i][j]);
            }
        }
        return m;
    }

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int j = 0; j < v.length; j++) {
            m.set(0, j, v[j]);
        }
        return m;
    }

    /** Rates of a c-server queue at each of M centers, mu(n) = min(n,c). */
    private static Matrix msRates(int M, int Nt, int c) {
        Matrix mu = new Matrix(M, Nt);
        for (int i = 0; i < M; i++) {
            for (int n = 1; n <= Nt; n++) {
                mu.set(i, n - 1, Math.min(n, c));
            }
        }
        return mu;
    }

    private static Matrix ones(int M, int Nt) {
        Matrix mu = new Matrix(M, Nt);
        mu.fill(1.0);
        return mu;
    }

    private static void agreesWithGld(String tag, Matrix L, Matrix N, Matrix mu) {
        Ret.pfqnNc ref = Pfqn_gld.pfqn_gld(L, N, mu, null);
        Pfqn_explicit_ld.Result got = Pfqn_explicit_ld.pfqn_explicit_ld(L, N, mu);
        assertTrue(Double.isFinite(got.lG), tag + ": the closed form refused");
        assertEquals(ref.lG, got.lG, TOL * Math.max(1.0, Math.abs(ref.lG)), tag);
    }

    @Test
    public void multiserverMulticlassMatchesGld() {
        Matrix L = mat(new double[][] {{1.2, 0.7}, {0.4, 1.9}});
        Matrix N = row(3, 2);
        Matrix mu = new Matrix(2, 5);
        for (int n = 1; n <= 5; n++) {
            mu.set(0, n - 1, Math.min(n, 2));
            mu.set(1, n - 1, Math.min(n, 3));
        }
        agreesWithGld("multiserver R=2", L, N, mu);
    }

    @Test
    public void fixedRateMulticlassMatchesGld() {
        Matrix L = mat(new double[][] {{1.0, 0.7, 0.3}, {0.5, 1.3, 0.9}, {0.9, 0.4, 1.7}});
        Matrix N = row(2, 1, 2);
        agreesWithGld("fixed rate R=3", L, N, ones(3, 5));
    }

    @Test
    public void arbitraryLimitedLoadDependenceMatchesGld() {
        // rates that neither increase nor follow a multi-server shape, settling at s=3
        Matrix L = mat(new double[][] {{1.1}, {0.6}, {2.3}});
        Matrix N = row(6);
        Matrix mu = mat(new double[][] {{0.5, 1.4, 2.2, 2.2, 2.2, 2.2},
                                        {1.0, 0.8, 1.7, 1.7, 1.7, 1.7},
                                        {2.0, 1.1, 0.9, 0.9, 0.9, 0.9}});
        agreesWithGld("arbitrary LLD R=1", L, N, mu);
    }

    @Test
    public void ratesThatNeverSettleAreStillExact() {
        // an infinite server never satisfies alpha(n)=alpha(s) for a finite s, but
        // s = |N| is admissible because no larger population occurs
        Matrix L = mat(new double[][] {{0.9, 1.4}, {1.1, 0.5}});
        Matrix N = row(2, 2);
        Matrix mu = new Matrix(2, 4);
        for (int n = 1; n <= 4; n++) {
            mu.set(0, n - 1, n);
            mu.set(1, n - 1, Math.min(n, 2));
        }
        agreesWithGld("never settles R=2", L, N, mu);
    }

    @Test
    public void aStationAClassNeverVisits() {
        // The oracle is the rational convolution over the state space, G = 0.755325,
        // computed offline: it is independent of Pfqn_gld, which used to be WRONG on
        // exactly this shape (see pfqnGldAgreesOnAZeroDemandLoadDependentModel).
        Matrix L = mat(new double[][] {{0.0, 1.3}, {0.9, 0.4}});
        Matrix N = row(2, 2);
        Pfqn_explicit_ld.Result got = Pfqn_explicit_ld.pfqn_explicit_ld(L, N, msRates(2, 4, 2));
        assertEquals(Math.log(0.755325), got.lG, 1e-9);
    }

    /**
     * Regression for the {@link Pfqn_gld} base case (fixed 2026-09-03).
     *
     * <p>A class with jobs and no demand at the ONLY station used to be dropped from the
     * single-station sum instead of zeroing the constant, and the recursion reaches that
     * base case with the full population every time it peels a station. So any
     * load-dependent model carrying a zero demand came back wrong. Fixed rate was
     * unaffected, which is why it hid for so long.</p>
     */
    @Test
    public void pfqnGldAgreesOnAZeroDemandLoadDependentModel() {
        Matrix N = row(2, 2);
        Matrix mu = msRates(2, 4, 2);
        Matrix L = mat(new double[][] {{0.0, 1.3}, {0.9, 0.4}});
        assertEquals(Math.log(0.755325), Pfqn_gld.pfqn_gld(L, N, mu, null).lG, 1e-9);
        // a whole station a class never visits, and the N=[2 3] fixture
        Matrix Lrow = mat(new double[][] {{0.0, 0.0}, {0.9, 0.4}});
        assertEquals(-2.3309845675200272, Pfqn_gld.pfqn_gld(Lrow, N, mu, null).lG, 1e-9);
        Matrix L2 = mat(new double[][] {{0.0, 1.3}, {0.9, 0.7}});
        assertEquals(Math.log(1.14240375),
                Pfqn_gld.pfqn_gld(L2, row(2, 3), msRates(2, 5, 2), null).lG, 1e-9);
        // and the fixed-rate model must be untouched by the fix
        Matrix ones = new Matrix(2, 4);
        ones.fill(1.0);
        assertEquals(1.2267416163786375, Pfqn_gld.pfqn_gld(L, N, ones, null).lG, 1e-9);
    }

    @Test
    public void repeatedScaledDemandsTakeEquationSixteen() {
        Matrix L = mat(new double[][] {{1.3, 0.8}, {1.3, 0.8}, {1.3, 0.8}});
        Matrix N = row(2, 2);
        Matrix mu = msRates(3, 4, 2);
        Pfqn_explicit_ld.Result got = Pfqn_explicit_ld.pfqn_explicit_ld(L, N, mu);
        assertEquals("repeated", got.method);
        agreesWithGld("repeated demands", L, N, mu);
    }

    @Test
    public void fixedRateDegenerationReproducesPfqnExplicit() {
        Matrix L = mat(new double[][] {{1.0, 0.7}, {0.5, 1.3}, {0.9, 0.4}});
        Matrix N = row(3, 2);
        Pfqn_explicit.Result a = Pfqn_explicit.pfqn_explicit(L, N);
        Pfqn_explicit_ld.Result b = Pfqn_explicit_ld.pfqn_explicit_ld(L, N, ones(3, 5));
        assertEquals(a.method, b.method);
        assertEquals(a.lG, b.lG, 1e-12);
    }

    @Test
    public void singleClassSkipsTheOuterSum() {
        // h_theta(N) is homogeneous of degree N in theta, so the divided difference is
        // the identity at R=1 and the constant must still match the convolution
        Matrix L = mat(new double[][] {{1.4}, {0.9}, {0.35}});
        Matrix N = row(7);
        agreesWithGld("single class", L, N, msRates(3, 7, 2));
    }

    @Test
    public void anUlpWideTieIsMissedAtEpsAndReported() {
        // both scaled demands are 1.95 at t=[2 3], but land two ulps apart in doubles
        Matrix L = mat(new double[][] {{0.0, 1.3}, {0.9, 0.7}});
        Matrix N = row(2, 3);
        Matrix mu = msRates(2, 5, 2);
        // the rational convolution over the state space, G = 1.14240375, pinned rather
        // than read off Pfqn_gld so this case does not depend on that recursion
        double lref = Math.log(1.14240375);

        Pfqn_explicit_ld.Result eps = Pfqn_explicit_ld.pfqn_explicit_ld(L, N, mu);
        assertEquals("distinct", eps.method);
        assertTrue(eps.lossDigits > 15,
                "an unmerged ulp-wide tie must be reported, got " + eps.lossDigits);
        assertTrue(Math.abs(lref - eps.lG) > 1.0,
                "the missed tie must actually move the answer");

        Pfqn_explicit_ld.Result merged = Pfqn_explicit_ld.pfqn_explicit_ld(
                L, N, mu, 1e-12, "auto", Double.POSITIVE_INFINITY);
        assertEquals("repeated", merged.method);
        assertEquals(lref, merged.lG, TOL * Math.max(1.0, Math.abs(lref)));

        // a caller holding a budget is refused rather than warned
        Pfqn_explicit_ld.Result budgeted = Pfqn_explicit_ld.pfqn_explicit_ld(
                L, N, mu, Math.ulp(1.0), "auto", 8.0);
        assertTrue(Double.isNaN(budgeted.lG));
    }

    @Test
    public void inadmissibleArgumentsAreRefused() {
        Matrix L = mat(new double[][] {{1.0, 0.5}, {0.5, 1.0}});
        Matrix N = row(2, 2);
        assertThrows(IllegalArgumentException.class, () -> Pfqn_explicit_ld.pfqn_explicit_ld(
                L, N, ones(2, 4), Math.ulp(1.0), "bogus", Double.POSITIVE_INFINITY));
        // a rate lattice shorter than the population cannot answer at |N|
        assertThrows(IllegalArgumentException.class,
                () -> Pfqn_explicit_ld.pfqn_explicit_ld(L, N, ones(2, 3)));
        // and a rate of zero would divide by zero in phi_k
        Matrix bad = ones(2, 4);
        bad.set(0, 2, 0.0);
        assertThrows(IllegalArgumentException.class,
                () -> Pfqn_explicit_ld.pfqn_explicit_ld(L, N, bad));
    }

    @Test
    public void theExplicitTokenRoutesThroughPfqnNcld() {
        Matrix L = mat(new double[][] {{1.2, 0.7}, {0.4, 1.9}});
        Matrix N = row(3, 2);
        Matrix mu = new Matrix(2, 5);
        for (int n = 1; n <= 5; n++) {
            mu.set(0, n - 1, Math.min(n, 2));
            mu.set(1, n - 1, Math.min(n, 3));
        }
        Matrix Z = new Matrix(1, 2);
        Z.fill(0.0);

        SolverOptions exact = SolverNC.defaultOptions();
        exact.method = "exact";
        double lref = Pfqn_ncld.pfqn_ncld(L, N, Z, mu, exact).lG;

        SolverOptions explicit = SolverNC.defaultOptions();
        explicit.method = "divdiff";
        Ret.pfqnNc got = Pfqn_ncld.pfqn_ncld(L, N, Z, mu, explicit);
        assertEquals("divdiff.ld/distinct", got.method);
        assertEquals(lref, got.lG, TOL * Math.max(1.0, Math.abs(lref)));

        // a think time would have to enter g_sigma, whose closed form covers queues only
        Matrix Zt = new Matrix(1, 2);
        Zt.fill(0.5);
        assertThrows(RuntimeException.class,
                () -> Pfqn_ncld.pfqn_ncld(L, N, Zt, mu, explicit));
    }
}
