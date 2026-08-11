package jline.api;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static jline.api.moment.Moment_binomial_from_factorial.moment_binomial_from_factorial;
import static jline.api.moment.Moment_binomial_from_negbinomial.moment_binomial_from_negbinomial;
import static jline.api.moment.Moment_binotrans.moment_binotrans;
import static jline.api.moment.Moment_binotransinv.moment_binotransinv;
import static jline.api.moment.Moment_central_from_raw.moment_central_from_raw;
import static jline.api.moment.Moment_factorial_from_binomial.moment_factorial_from_binomial;
import static jline.api.moment.Moment_factorial_from_raw.moment_factorial_from_raw;
import static jline.api.moment.Moment_factorial_from_upfactorial.moment_factorial_from_upfactorial;
import static jline.api.moment.Moment_lah.moment_lah;
import static jline.api.moment.Moment_cumulant_from_raw.moment_cumulant_from_raw;
import static jline.api.moment.Moment_factcumulant_from_factorial.moment_factcumulant_from_factorial;
import static jline.api.moment.Moment_factorial_from_factcumulant.moment_factorial_from_factcumulant;
import static jline.api.moment.Moment_housematrix.moment_housematrix;
import static jline.api.moment.Moment_joint_aggregate.moment_joint_aggregate;
import static jline.api.moment.Moment_joint_binomial_from_factorial.moment_joint_binomial_from_factorial;
import static jline.api.moment.Moment_joint_binomial_from_negbinomial.moment_joint_binomial_from_negbinomial;
import static jline.api.moment.Moment_joint_central_from_raw.moment_joint_central_from_raw;
import static jline.api.moment.Moment_joint_central_from_raw_mean.moment_joint_central_from_raw_mean;
import static jline.api.moment.Moment_joint_cumulant_from_raw.moment_joint_cumulant_from_raw;
import static jline.api.moment.Moment_joint_factcumulant_from_factorial.moment_joint_factcumulant_from_factorial;
import static jline.api.moment.Moment_joint_factorial_from_binomial.moment_joint_factorial_from_binomial;
import static jline.api.moment.Moment_joint_factorial_from_factcumulant.moment_joint_factorial_from_factcumulant;
import static jline.api.moment.Moment_joint_factorial_from_raw.moment_joint_factorial_from_raw;
import static jline.api.moment.Moment_joint_factorial_from_upfactorial.moment_joint_factorial_from_upfactorial;
import static jline.api.moment.Moment_joint_marking.moment_joint_marking;
import static jline.api.moment.Moment_joint_negbinomial_from_binomial.moment_joint_negbinomial_from_binomial;
import static jline.api.moment.Moment_joint_negbinomial_from_upfactorial.moment_joint_negbinomial_from_upfactorial;
import static jline.api.moment.Moment_joint_raw_from_central.moment_joint_raw_from_central;
import static jline.api.moment.Moment_joint_raw_from_cumulant.moment_joint_raw_from_cumulant;
import static jline.api.moment.Moment_joint_raw_from_factorial.moment_joint_raw_from_factorial;
import static jline.api.moment.Moment_joint_raw_from_upfactorial.moment_joint_raw_from_upfactorial;
import static jline.api.moment.Moment_joint_upfactorial_from_factorial.moment_joint_upfactorial_from_factorial;
import static jline.api.moment.Moment_joint_upfactorial_from_negbinomial.moment_joint_upfactorial_from_negbinomial;
import static jline.api.moment.Moment_joint_upfactorial_from_raw.moment_joint_upfactorial_from_raw;
import static jline.api.moment.Moment_raw_from_cumulant.moment_raw_from_cumulant;
import static jline.api.moment.Moment_tensortrans.moment_tensortrans;
import static jline.api.moment.Moment_negbinomial_from_binomial.moment_negbinomial_from_binomial;
import static jline.api.moment.Moment_negbinomial_from_upfactorial.moment_negbinomial_from_upfactorial;
import static jline.api.moment.Moment_raw_from_central.moment_raw_from_central;
import static jline.api.moment.Moment_raw_from_factorial.moment_raw_from_factorial;
import static jline.api.moment.Moment_raw_from_upfactorial.moment_raw_from_upfactorial;
import static jline.api.moment.Moment_stirling1.moment_stirling1;
import static jline.api.moment.Moment_stirling2.moment_stirling2;
import static jline.api.moment.Moment_stirlingcycle.moment_stirlingcycle;
import static jline.api.moment.Moment_upfactorial_from_factorial.moment_upfactorial_from_factorial;
import static jline.api.moment.Moment_upfactorial_from_negbinomial.moment_upfactorial_from_negbinomial;
import static jline.api.moment.Moment_upfactorial_from_raw.moment_upfactorial_from_raw;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Unit tests for the moment-conversion API, which implements the "house of
 * moments" of
 *
 * <p>A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * <p>The conversions link six families of moments of a discrete random variable
 * N: the power (raw) moments m_n = E[N^n], the central moments
 * m_n^c = E[(N-m_1)^n], the factorial moments f_n = E[N(N-1)...(N-n+1)], the
 * binomial moments b_n = E[C(N,n)], and the two families introduced by the
 * reference, the upward-factorial moments f_n^+ = E[N(N+1)...(N+n-1)] and the
 * negative-binomial moments b_n^- = E[C(N+n-1,n)].
 *
 * <p>Every vector uses the harmonized 0-based order convention of the
 * reference: a vector of length n+1 holds the moments of order 0..n, element 0
 * being the order-0 moment, which equals 1 for every family.
 *
 * <p>The conversions are refereed by a brute-force oracle that evaluates each
 * family by direct summation over the pmf, and by closed forms holding for
 * specific distributions. These tests verify that the JAR implementation
 * matches the MATLAB version; see test_moment.m in line-test.git and
 * python/tests/test_moment.py.
 */
public class MomentAPITest {

    private static final double TOL = 1e-10;

    // ---------- helpers --------------------------------------------------

    private static Matrix vec(double[] v) {
        Matrix m = new Matrix(v.length, 1);
        for (int i = 0; i < v.length; i++) {
            m.set(i, 0, v[i]);
        }
        return m;
    }

    private static double[] arr(Matrix m) {
        double[] v = new double[m.length()];
        for (int i = 0; i < v.length; i++) {
            v[i] = m.get(i);
        }
        return v;
    }

    /**
     * Normwise relative error max|a-b| / max(1,max|b|).
     *
     * <p>A componentwise relative error is not a usable criterion here. The
     * triangles carry alternating signs and large binomial coefficients, so an
     * entry that is exactly zero in exact arithmetic is reached by cancellation
     * between terms of the magnitude of the input: for Uniform(0..5) the
     * falling factorial annihilates every atom at order 6, so f_6 = 0 exactly,
     * yet it is computed by cancelling terms of size m_6 = 3419 and lands on
     * ~5e-13 of roundoff. Dividing that by |b_i| = 0 is meaningless. The
     * normwise residual is the standard conditioning-aware criterion for a
     * linear transform and still resolves a genuine defect to ~1e-14 here.
     */
    private static double relerr(double[] a, double[] b) {
        double num = 0.0;
        double den = 1.0;
        for (int i = 0; i < b.length; i++) {
            num = Math.max(num, Math.abs(a[i] - b[i]));
            den = Math.max(den, Math.abs(b[i]));
        }
        return num / den;
    }

    private static double relerr(Matrix a, double[] b) {
        return relerr(arr(a), b);
    }

    /** Exact binomial coefficient by the incremental product. */
    private static double binom(long n, long k) {
        if (k < 0 || k > n) {
            return 0.0;
        }
        double c = 1.0;
        for (long i = 0; i < k; i++) {
            c = c * (n - i) / (i + 1);
        }
        return c;
    }

    /** Every moment family by direct summation over the pmf, orders 0..N. */
    private static final class Oracle {
        double[] m;
        double[] f;
        double[] fp;
        double[] b;
        double[] bm;
        double[] mc;
        double m1;
    }

    private static Oracle pmfOracle(double[] kk, double[] pk, int N) {
        Oracle o = new Oracle();
        o.m = new double[N + 1];
        o.f = new double[N + 1];
        o.fp = new double[N + 1];
        o.b = new double[N + 1];
        o.bm = new double[N + 1];
        o.mc = new double[N + 1];
        for (int n = 0; n <= N; n++) {
            for (int ii = 0; ii < kk.length; ii++) {
                double k = kk[ii];
                o.m[n] += Math.pow(k, n) * pk[ii];
                double pf = 1.0;
                double pfp = 1.0;
                for (int j = 0; j < n; j++) {
                    pf *= (k - j);      // falling factorial k(k-1)...(k-n+1)
                    pfp *= (k + j);     // rising  factorial k(k+1)...(k+n-1)
                }
                o.f[n] += pf * pk[ii];
                o.fp[n] += pfp * pk[ii];
                long ki = Math.round(k);
                if (ki >= n) {
                    o.b[n] += binom(ki, n) * pk[ii];
                }
                if (n == 0) {
                    o.bm[n] += pk[ii];
                } else {
                    // the k=0 term contributes C(n-1,n) = 0, cf. eq. (7)
                    o.bm[n] += binom(ki + n - 1, n) * pk[ii];
                }
            }
        }
        o.m1 = o.m[1];
        for (int n = 0; n <= N; n++) {
            for (int ii = 0; ii < kk.length; ii++) {
                o.mc[n] += Math.pow(kk[ii] - o.m1, n) * pk[ii];
            }
        }
        return o;
    }

    /**
     * Four pmfs with distinct structure. All have mass at k=0, exercising the
     * k=0 term of the negative-binomial moments, and all have finite support so
     * the oracle is exact (the Poisson case is truncated far beyond its mean).
     */
    private static Oracle[] cases(int N) {
        Oracle[] out = new Oracle[4];

        double[] kk1 = new double[11];
        double[] pk1 = new double[11];
        for (int k = 0; k <= 10; k++) {
            kk1[k] = k;
            pk1[k] = binom(10, k) * Math.pow(0.3, k) * Math.pow(0.7, 10 - k);
        }
        out[0] = pmfOracle(kk1, pk1, N);

        double[] kk2 = new double[6];
        double[] pk2 = new double[6];
        for (int k = 0; k <= 5; k++) {
            kk2[k] = k;
            pk2[k] = 1.0 / 6.0;
        }
        out[1] = pmfOracle(kk2, pk2, N);

        double[] kk3 = {0, 1, 2, 3, 4, 5};
        double[] pk3 = {0.10, 0.20, 0.05, 0.30, 0.15, 0.20};
        out[2] = pmfOracle(kk3, pk3, N);

        double[] kk4 = new double[61];
        double[] pk4 = new double[61];
        double term = Math.exp(-2.0);
        for (int k = 0; k <= 60; k++) {
            kk4[k] = k;
            pk4[k] = term;
            term = term * 2.0 / (k + 1);
        }
        out[3] = pmfOracle(kk4, pk4, N);

        return out;
    }

    private static double[] row(Matrix m, int i) {
        double[] r = new double[m.getNumCols()];
        for (int j = 0; j < r.length; j++) {
            r[j] = m.get(i, j);
        }
        return r;
    }

    // ---------- triangles ------------------------------------------------

    @Test
    public void testStirling1Triangle() {
        Matrix s = moment_stirling1(4);
        assertArrayEquals(new double[]{1, 0, 0, 0, 0}, row(s, 0), 1e-12);
        assertArrayEquals(new double[]{0, 1, 0, 0, 0}, row(s, 1), 1e-12);
        assertArrayEquals(new double[]{0, -1, 1, 0, 0}, row(s, 2), 1e-12);
        assertArrayEquals(new double[]{0, 2, -3, 1, 0}, row(s, 3), 1e-12);
        assertArrayEquals(new double[]{0, -6, 11, -6, 1}, row(s, 4), 1e-12);
    }

    @Test
    public void testStirling2Triangle() {
        Matrix S = moment_stirling2(4);
        assertArrayEquals(new double[]{1, 0, 0, 0, 0}, row(S, 0), 1e-12);
        assertArrayEquals(new double[]{0, 1, 0, 0, 0}, row(S, 1), 1e-12);
        assertArrayEquals(new double[]{0, 1, 1, 0, 0}, row(S, 2), 1e-12);
        assertArrayEquals(new double[]{0, 1, 3, 1, 0}, row(S, 3), 1e-12);
        assertArrayEquals(new double[]{0, 1, 7, 6, 1}, row(S, 4), 1e-12);
    }

    @Test
    public void testStirlingCycleTriangle() {
        Matrix sigma = moment_stirlingcycle(4);
        assertArrayEquals(new double[]{1, 0, 0, 0, 0}, row(sigma, 0), 1e-12);
        assertArrayEquals(new double[]{0, 1, 0, 0, 0}, row(sigma, 1), 1e-12);
        assertArrayEquals(new double[]{0, 1, 1, 0, 0}, row(sigma, 2), 1e-12);
        assertArrayEquals(new double[]{0, 2, 3, 1, 0}, row(sigma, 3), 1e-12);
        assertArrayEquals(new double[]{0, 6, 11, 6, 1}, row(sigma, 4), 1e-12);
        // row sums equal n!, the permutations being partitioned by cycle count
        double fact = 1.0;
        for (int n = 0; n <= 4; n++) {
            if (n > 0) {
                fact *= n;
            }
            double sum = 0.0;
            for (double v : row(sigma, n)) {
                sum += v;
            }
            assertTrue(Math.abs(sum - fact) < 1e-12, "row sum of sigma at n=" + n);
        }
    }

    @Test
    public void testLahTriangle() {
        Matrix L = moment_lah(4);
        assertArrayEquals(new double[]{1, 0, 0, 0, 0}, row(L, 0), 1e-12);
        assertArrayEquals(new double[]{0, 1, 0, 0, 0}, row(L, 1), 1e-12);
        assertArrayEquals(new double[]{0, 2, 1, 0, 0}, row(L, 2), 1e-12);
        assertArrayEquals(new double[]{0, 6, 6, 1, 0}, row(L, 3), 1e-12);
        assertArrayEquals(new double[]{0, 24, 36, 12, 1}, row(L, 4), 1e-12);
    }

    @Test
    public void testLahMatchesExplicitFactorialForm() {
        // the recursion must agree with L(n,k) = (n!/k!)*C(n-1,k-1)
        int n = 8;
        Matrix L = moment_lah(n);
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= i; j++) {
                double fi = 1.0;
                for (int t = 2; t <= i; t++) {
                    fi *= t;
                }
                double fj = 1.0;
                for (int t = 2; t <= j; t++) {
                    fj *= t;
                }
                double expected = (fi / fj) * binom(i - 1, j - 1);
                assertTrue(Math.abs(L.get(i, j) - expected) <= 1e-12 * Math.max(1.0, expected),
                        "L(" + i + "," + j + ")");
            }
        }
    }

    @Test
    public void testStirlingSignRelation() {
        // s(n,k) = (-1)^(n-k) * sigma(n,k), eq. (12)
        int n = 7;
        Matrix s = moment_stirling1(n);
        Matrix sigma = moment_stirlingcycle(n);
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= i; j++) {
                double expected = (((i - j) % 2 == 0) ? 1.0 : -1.0) * sigma.get(i, j);
                assertTrue(Math.abs(s.get(i, j) - expected) <= 1e-12 * Math.max(1.0, Math.abs(expected)),
                        "s(" + i + "," + j + ")");
            }
        }
    }

    @Test
    public void testStirlingTrianglesAreInverse() {
        // eqs. (10) and (11) make the two Stirling triangles mutually inverse
        int n = 7;
        Matrix s = moment_stirling1(n);
        Matrix S = moment_stirling2(n);
        Matrix p = S.mult(s);
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= n; j++) {
                assertTrue(Math.abs(p.get(i, j) - (i == j ? 1.0 : 0.0)) < 1e-10,
                        "S*s at (" + i + "," + j + ")");
            }
        }
    }

    @Test
    public void testStirling1GeneratesFallingFactorial() {
        // defining identity (10): sum_k s(n,k) x^k = x(x-1)...(x-n+1)
        int n = 6;
        Matrix s = moment_stirling1(n);
        double[] xs = {-3, -0.5, 0, 1, 2.5, 4, 9};
        for (double x : xs) {
            for (int i = 0; i <= n; i++) {
                double lhs = 0.0;
                for (int j = 0; j <= i; j++) {
                    lhs += s.get(i, j) * Math.pow(x, j);
                }
                double rhs = 1.0;
                for (int j = 0; j < i; j++) {
                    rhs *= (x - j);
                }
                assertTrue(Math.abs(lhs - rhs) <= 1e-8 + 1e-10 * Math.abs(rhs),
                        "x=" + x + " n=" + i);
            }
        }
    }

    @Test
    public void testStirling2ExpandsPowerIntoFallingFactorials() {
        // defining identity (11): x^n = sum_k S(n,k) x(x-1)...(x-k+1)
        int n = 6;
        Matrix S = moment_stirling2(n);
        double[] xs = {-3, -0.5, 0, 1, 2.5, 4, 9};
        for (double x : xs) {
            for (int i = 0; i <= n; i++) {
                double rhs = 0.0;
                for (int j = 0; j <= i; j++) {
                    double ff = 1.0;
                    for (int l = 0; l < j; l++) {
                        ff *= (x - l);
                    }
                    rhs += S.get(i, j) * ff;
                }
                assertTrue(Math.abs(Math.pow(x, i) - rhs) <= 1e-8 + 1e-10 * Math.abs(Math.pow(x, i)),
                        "x=" + x + " n=" + i);
            }
        }
    }

    @Test
    public void testTriangleInputValidation() {
        assertThrows(IllegalArgumentException.class, () -> moment_stirling1(-1));
        assertThrows(IllegalArgumentException.class, () -> moment_stirling2(-1));
        assertThrows(IllegalArgumentException.class, () -> moment_stirlingcycle(-1));
        assertThrows(IllegalArgumentException.class, () -> moment_lah(-4));
    }

    // ---------- binomial transform ---------------------------------------

    @Test
    public void testBinotransKnownValues() {
        // eq. (8) evaluated by hand for x = [1,2,5,15]:
        //   y_0 =  1
        //   y_1 = -1 + 2            =  1
        //   y_2 =  1 - 4 + 5        =  2
        //   y_3 = -1 + 6 - 15 + 15  =  5
        Matrix y = moment_binotrans(vec(new double[]{1, 2, 5, 15}));
        assertArrayEquals(new double[]{1, 1, 2, 5}, arr(y), 1e-12);
    }

    @Test
    public void testBinotransIsInvolution() {
        // eq. (8) and its inverse (9) undo one another in both orders
        double[] x = {1, 2.5, 7, -3, 11, 0.5};
        assertTrue(relerr(moment_binotransinv(moment_binotrans(vec(x))), x) < TOL);
        assertTrue(relerr(moment_binotrans(moment_binotransinv(vec(x))), x) < TOL);
    }

    // ---------- conversions against the pmf oracle ------------------------

    @Test
    public void testAllConversionsVsPmfOracle() {
        // the acceptance gate: all 14 conversions refereed by direct summation
        Oracle[] cs = cases(6);
        for (int c = 0; c < cs.length; c++) {
            Oracle o = cs[c];
            String at = "case " + c;
            // power <-> factorial, eq. (13)
            assertTrue(relerr(moment_factorial_from_raw(vec(o.m)), o.f) < TOL, at);
            assertTrue(relerr(moment_raw_from_factorial(vec(o.f)), o.m) < TOL, at);
            // power <-> upward-factorial, via the Stirling cycle numbers
            assertTrue(relerr(moment_upfactorial_from_raw(vec(o.m)), o.fp) < TOL, at);
            assertTrue(relerr(moment_raw_from_upfactorial(vec(o.fp)), o.m) < TOL, at);
            // factorial <-> binomial
            assertTrue(relerr(moment_binomial_from_factorial(vec(o.f)), o.b) < TOL, at);
            assertTrue(relerr(moment_factorial_from_binomial(vec(o.b)), o.f) < TOL, at);
            // upward-factorial <-> negative-binomial, eq. (7)
            assertTrue(relerr(moment_negbinomial_from_upfactorial(vec(o.fp)), o.bm) < TOL, at);
            assertTrue(relerr(moment_upfactorial_from_negbinomial(vec(o.bm)), o.fp) < TOL, at);
            // binomial <-> negative-binomial, the shifted binomial transform, eq. (14)
            assertTrue(relerr(moment_binomial_from_negbinomial(vec(o.bm)), o.b) < TOL, at);
            assertTrue(relerr(moment_negbinomial_from_binomial(vec(o.b)), o.bm) < TOL, at);
            // factorial <-> upward-factorial, via the Lah numbers
            assertTrue(relerr(moment_factorial_from_upfactorial(vec(o.fp)), o.f) < TOL, at);
            assertTrue(relerr(moment_upfactorial_from_factorial(vec(o.f)), o.fp) < TOL, at);
            // power <-> central
            assertTrue(relerr(moment_central_from_raw(vec(o.m)), o.mc) < TOL, at);
            assertTrue(relerr(moment_raw_from_central(vec(o.mc), o.m1), o.m) < TOL, at);
        }
    }

    @Test
    public void testOrderZeroAndOneInvariants() {
        // every family agrees at orders 0 and 1; m_1^c = 0 (cf. Section 2)
        Oracle[] cs = cases(5);
        for (Oracle o : cs) {
            Matrix f = moment_factorial_from_raw(vec(o.m));
            Matrix fp = moment_upfactorial_from_raw(vec(o.m));
            Matrix b = moment_binomial_from_factorial(f);
            Matrix bm = moment_negbinomial_from_upfactorial(fp);
            Matrix mc = moment_central_from_raw(vec(o.m));
            assertTrue(Math.abs(o.m[0] - 1.0) < 1e-12);
            assertTrue(Math.abs(f.get(0) - 1.0) < 1e-12);
            assertTrue(Math.abs(fp.get(0) - 1.0) < 1e-12);
            assertTrue(Math.abs(b.get(0) - 1.0) < 1e-12);
            assertTrue(Math.abs(bm.get(0) - 1.0) < 1e-12);
            assertTrue(Math.abs(mc.get(0) - 1.0) < 1e-12);
            double s = Math.max(1.0, Math.abs(o.m1));
            assertTrue(Math.abs(f.get(1) - o.m1) <= 1e-10 * s);
            assertTrue(Math.abs(fp.get(1) - o.m1) <= 1e-10 * s);
            assertTrue(Math.abs(b.get(1) - o.m1) <= 1e-10 * s);
            assertTrue(Math.abs(bm.get(1) - o.m1) <= 1e-10 * s);
            assertTrue(Math.abs(mc.get(1)) < 1e-10);
        }
    }

    // ---------- conversions against closed forms --------------------------

    @Test
    public void testPoissonFactorialMomentsClosedForm() {
        // Poisson(lambda) has f_n = lambda^n exactly
        int N = 6;
        double[] lambdas = {0.5, 2.0, 7.0};
        for (double lambda : lambdas) {
            double[] kk = new double[121];
            double[] pk = new double[121];
            double term = Math.exp(-lambda);
            for (int k = 0; k <= 120; k++) {
                kk[k] = k;
                pk[k] = term;
                term = term * lambda / (k + 1);
            }
            Matrix f = moment_factorial_from_raw(vec(pmfOracle(kk, pk, N).m));
            double[] expected = new double[N + 1];
            for (int n = 0; n <= N; n++) {
                expected[n] = Math.pow(lambda, n);
            }
            assertTrue(relerr(f, expected) < 1e-8, "lambda=" + lambda);
        }
    }

    @Test
    public void testBinomialFactorialMomentsClosedForm() {
        // Binomial(K,p) has f_n = (K!/(K-n)!)*p^n, and f_n = 0 for n > K
        int N = 6;
        int K = 5;
        double p = 0.4;
        double[] kk = new double[K + 1];
        double[] pk = new double[K + 1];
        for (int k = 0; k <= K; k++) {
            kk[k] = k;
            pk[k] = binom(K, k) * Math.pow(p, k) * Math.pow(1 - p, K - k);
        }
        Matrix f = moment_factorial_from_raw(vec(pmfOracle(kk, pk, N).m));
        double[] expected = new double[N + 1];
        for (int n = 0; n <= N; n++) {
            if (n <= K) {
                double ff = 1.0;
                for (int j = 0; j < n; j++) {
                    ff *= (K - j);
                }
                expected[n] = ff * Math.pow(p, n);
            } else {
                expected[n] = 0.0;   // the falling factorial annihilates orders above K
            }
        }
        assertTrue(relerr(f, expected) < 1e-8);
        assertTrue(Math.abs(f.get(K + 1)) < 1e-8);
    }

    @Test
    public void testBernoulliAllFamiliesByHand() {
        // Bernoulli(p): m_n = p, f_n = 0 (n>=2), f_n^+ = n!*p, b_n^- = p
        int N = 5;
        double p = 0.3;
        Oracle o = pmfOracle(new double[]{0, 1}, new double[]{1 - p, p}, N);
        Matrix f = moment_factorial_from_raw(vec(o.m));
        Matrix fp = moment_upfactorial_from_raw(vec(o.m));
        Matrix b = moment_binomial_from_factorial(f);
        Matrix bm = moment_negbinomial_from_upfactorial(fp);
        double fact = 1.0;
        for (int n = 0; n <= N; n++) {
            if (n > 0) {
                fact *= n;
            }
            double expectedM = (n == 0) ? 1.0 : p;
            assertTrue(Math.abs(o.m[n] - expectedM) < 1e-12, "m_" + n);
            double expectedF = (n == 0) ? 1.0 : (n == 1 ? p : 0.0);
            assertTrue(Math.abs(f.get(n) - expectedF) < 1e-10, "f_" + n);
            assertTrue(Math.abs(b.get(n) - expectedF) < 1e-10, "b_" + n);
            double expectedFp = (n == 0) ? 1.0 : fact * p;
            assertTrue(Math.abs(fp.get(n) - expectedFp) <= 1e-10 * Math.max(1.0, expectedFp), "fp_" + n);
            assertTrue(Math.abs(bm.get(n) - expectedM) < 1e-10, "bm_" + n);
        }
    }

    @Test
    public void testDeterministicPointMass() {
        // a point mass at c reduces the conversions to the defining identities
        int N = 6;
        double c = 4;
        Oracle o = pmfOracle(new double[]{c}, new double[]{1.0}, N);
        Matrix f = moment_factorial_from_raw(vec(o.m));
        Matrix fp = moment_upfactorial_from_raw(vec(o.m));
        for (int n = 0; n <= N; n++) {
            double ff = 1.0;
            double rf = 1.0;
            for (int j = 0; j < n; j++) {
                ff *= (c - j);
                rf *= (c + j);
            }
            assertTrue(Math.abs(f.get(n) - ff) <= 1e-8 + 1e-10 * Math.abs(ff), "f_" + n);
            assertTrue(Math.abs(fp.get(n) - rf) <= 1e-8 + 1e-10 * Math.abs(rf), "fp_" + n);
        }
    }

    // ---------- structural properties -------------------------------------

    @Test
    public void testHouseOfMomentsCommutes() {
        // Figure 1 is a commuting diagram: two routes to a family must agree.
        // This catches a wrong triangle on one edge that a pure round-trip test
        // (traversing the same edge both ways) cannot see.
        Oracle[] cs = cases(6);
        for (int c = 0; c < cs.length; c++) {
            Oracle o = cs[c];
            String at = "case " + c;
            Matrix f = moment_factorial_from_raw(vec(o.m));
            Matrix fp = moment_upfactorial_from_raw(vec(o.m));
            // raw -> upfactorial directly, vs raw -> factorial -> upfactorial
            assertTrue(relerr(moment_upfactorial_from_factorial(f), arr(fp)) < TOL, at);
            // raw -> factorial directly, vs raw -> upfactorial -> factorial
            assertTrue(relerr(moment_factorial_from_upfactorial(fp), arr(f)) < TOL, at);
            // binomial by the n! route vs by the shifted-binomial-transform route
            Matrix bm = moment_negbinomial_from_upfactorial(fp);
            assertTrue(relerr(moment_binomial_from_negbinomial(bm),
                    arr(moment_binomial_from_factorial(f))) < TOL, at);
            // negative-binomial by the n! route vs by the shifted-transform route
            Matrix b = moment_binomial_from_factorial(f);
            assertTrue(relerr(moment_negbinomial_from_binomial(b), arr(bm)) < TOL, at);
        }
    }

    @Test
    public void testRoundtrips() {
        // each conversion pair composes to the identity
        Oracle[] cs = cases(6);
        for (int c = 0; c < cs.length; c++) {
            Oracle o = cs[c];
            String at = "case " + c;
            assertTrue(relerr(moment_raw_from_factorial(moment_factorial_from_raw(vec(o.m))), o.m) < TOL, at);
            assertTrue(relerr(moment_raw_from_upfactorial(moment_upfactorial_from_raw(vec(o.m))), o.m) < TOL, at);
            assertTrue(relerr(moment_raw_from_central(moment_central_from_raw(vec(o.m)), o.m1), o.m) < TOL, at);
            Matrix f = moment_factorial_from_raw(vec(o.m));
            assertTrue(relerr(moment_factorial_from_binomial(moment_binomial_from_factorial(f)), arr(f)) < TOL, at);
            assertTrue(relerr(moment_factorial_from_upfactorial(moment_upfactorial_from_factorial(f)), arr(f)) < TOL, at);
            Matrix fp = moment_upfactorial_from_raw(vec(o.m));
            assertTrue(relerr(moment_upfactorial_from_negbinomial(
                    moment_negbinomial_from_upfactorial(fp)), arr(fp)) < TOL, at);
            Matrix b = moment_binomial_from_factorial(f);
            assertTrue(relerr(moment_binomial_from_negbinomial(
                    moment_negbinomial_from_binomial(b)), arr(b)) < TOL, at);
        }
    }

    @Test
    public void testCentralConversionHoldsForContinuousRv() {
        // Section 5: the power <-> central rules also hold for continuous r.v.s.
        // Refereed on Gamma(k,theta), whose raw moments are
        // m_n = theta^n * gamma(k+n)/gamma(k) and whose variance is k*theta^2.
        int N = 5;
        int k = 3;
        double theta = 2;
        double[] m = new double[N + 1];
        for (int n = 0; n <= N; n++) {
            double ratio = 1.0;   // gamma(k+n)/gamma(k) = (k)(k+1)...(k+n-1)
            for (int j = 0; j < n; j++) {
                ratio *= (k + j);
            }
            m[n] = Math.pow(theta, n) * ratio;
        }
        Matrix mc = moment_central_from_raw(vec(m));
        assertTrue(Math.abs(mc.get(0) - 1.0) < 1e-12);
        assertTrue(Math.abs(mc.get(1)) < 1e-10);
        assertTrue(Math.abs(mc.get(2) - k * theta * theta) <= 1e-10 * k * theta * theta);
        assertTrue(Math.abs(mc.get(3) - 2 * k * Math.pow(theta, 3)) <= 1e-10 * 2 * k * Math.pow(theta, 3));
        assertTrue(relerr(moment_raw_from_central(mc, m[1]), m) < TOL);
    }

    // ---------- interface behaviour ---------------------------------------

    @Test
    public void testScalarInputIsOrderZeroOnly() {
        // a length-1 input carries only the order-0 moment and passes through
        Matrix one = vec(new double[]{1.0});
        assertTrue(Math.abs(moment_factorial_from_raw(one).get(0) - 1.0) < 1e-12);
        assertTrue(Math.abs(moment_upfactorial_from_raw(one).get(0) - 1.0) < 1e-12);
        assertTrue(Math.abs(moment_binomial_from_factorial(one).get(0) - 1.0) < 1e-12);
        assertTrue(Math.abs(moment_binomial_from_negbinomial(one).get(0) - 1.0) < 1e-12);
        assertTrue(Math.abs(moment_factorial_from_upfactorial(one).get(0) - 1.0) < 1e-12);
    }

    @Test
    public void testCentralFromRawRequiresTheMean() {
        // m_n^c is defined relative to m_1, so a length-1 input is rejected
        // rather than silently treated as a zero mean
        assertThrows(IllegalArgumentException.class,
                () -> moment_central_from_raw(vec(new double[]{1.0})));
    }

    // ---------- cumulants and the multivariate house ---------------------

    /**
     * Brute-force joint moment families of a discrete law, by direct summation
     * over the pmf. Row-major layout, matching the joint API.
     *
     * @param prob probabilities of the atoms
     * @param val atom values, val[r][j] being the j-th coordinate of atom r
     * @param dims extents of the requested array
     * @param kind 0 for power, 1 for factorial, 2 for upward-factorial moments
     * @return flattened joint moment array
     */
    private static double[] jointOracle(double[] prob, double[][] val, int[] dims, int kind) {
        int d = dims.length;
        int nel = 1;
        for (int j = 0; j < d; j++) {
            nel *= dims[j];
        }
        double[] out = new double[nel];
        int[] a = new int[d];
        for (int ia = 0; ia < nel; ia++) {
            double acc = 0.0;
            for (int r = 0; r < prob.length; r++) {
                double term = prob[r];
                for (int j = 0; j < d; j++) {
                    double x = val[r][j];
                    if (kind == 0) {
                        term *= Math.pow(x, a[j]);
                    } else {
                        for (int t = 0; t < a[j]; t++) {
                            term *= (kind == 1) ? (x - t) : (x + t);
                        }
                    }
                }
                acc += term;
            }
            out[ia] = acc;
            for (int l = d - 1; l >= 0; l--) {
                a[l]++;
                if (a[l] < dims[l]) {
                    break;
                }
                a[l] = 0;
            }
        }
        return out;
    }

    private static final double[] LAW_P = {0.4, 0.3, 0.2, 0.1};
    private static final double[][] LAW_V = {{1, 2}, {3, 0}, {2, 5}, {0, 1}};

    @Test
    public void testCumulantsOfThePoisson() {
        // every cumulant of a Poisson equals its rate, and the raw moments are
        // the Touchard polynomials
        double[] k = {0, 2, 2, 2, 2, 2};
        Matrix m = moment_raw_from_cumulant(vec(k));
        assertTrue(relerr(m, new double[]{1, 2, 6, 22, 94, 454}) < TOL);
        assertTrue(relerr(moment_cumulant_from_raw(m), k) < TOL);
    }

    @Test
    public void testLowOrderCumulantsAreTheTextbookOnes() {
        double[] m = {1, 2, 6, 22, 94};
        double[] k = arr(moment_cumulant_from_raw(vec(m)));
        double[] mc = arr(moment_central_from_raw(vec(m)));
        assertTrue(Math.abs(k[0]) < TOL);
        assertTrue(Math.abs(k[1] - m[1]) < TOL);
        assertTrue(Math.abs(k[2] - (m[2] - m[1] * m[1])) < TOL);
        assertTrue(Math.abs(k[3] - mc[3]) < TOL);
        assertTrue(Math.abs(k[4] - (mc[4] - 3 * mc[2] * mc[2])) < TOL);
    }

    @Test
    public void testFactorialCumulantsOfThePoissonVanish() {
        // f_n = lambda^n for a Poisson, so only the first factorial cumulant
        // survives: they measure the departure from Poisson behaviour
        double lam = 1.7;
        double[] f = new double[6];
        for (int n = 0; n < 6; n++) {
            f[n] = Math.pow(lam, n);
        }
        double[] kf = arr(moment_factcumulant_from_factorial(vec(f)));
        assertTrue(Math.abs(kf[1] - lam) < TOL);
        for (int n = 2; n < 6; n++) {
            assertTrue(Math.abs(kf[n]) < TOL, "order " + n);
        }
        assertTrue(relerr(moment_factorial_from_factcumulant(vec(kf)), f) < TOL);
    }

    @Test
    public void testJointConversionsAgainstTheOracle() {
        int[] dims = {5, 5};
        double[] m = jointOracle(LAW_P, LAW_V, dims, 0);
        double[] f = jointOracle(LAW_P, LAW_V, dims, 1);
        double[] fp = jointOracle(LAW_P, LAW_V, dims, 2);
        assertTrue(relerr(moment_joint_factorial_from_raw(m, dims), f) < TOL);
        assertTrue(relerr(moment_joint_raw_from_factorial(f, dims), m) < TOL);
        assertTrue(relerr(moment_joint_upfactorial_from_raw(m, dims), fp) < TOL);
        assertTrue(relerr(moment_joint_raw_from_upfactorial(fp, dims), m) < TOL);
        assertTrue(relerr(moment_joint_factorial_from_upfactorial(fp, dims), f) < TOL);
        assertTrue(relerr(moment_joint_upfactorial_from_factorial(f, dims), fp) < TOL);
        double[] b = moment_joint_binomial_from_factorial(f, dims);
        double[] bn = moment_joint_negbinomial_from_upfactorial(fp, dims);
        assertTrue(relerr(moment_joint_factorial_from_binomial(b, dims), f) < TOL);
        assertTrue(relerr(moment_joint_upfactorial_from_negbinomial(bn, dims), fp) < TOL);
        assertTrue(relerr(moment_joint_negbinomial_from_binomial(b, dims), bn) < TOL);
        assertTrue(relerr(moment_joint_binomial_from_negbinomial(bn, dims), b) < TOL);
    }

    @Test
    public void testJointConversionsReduceToTheUnivariateOnes() {
        // at d = 1 the two APIs must agree exactly
        double[] m = {1, 2, 6, 22, 94};
        int[] dims = {5};
        assertTrue(relerr(moment_joint_factorial_from_raw(m, dims),
                arr(moment_factorial_from_raw(vec(m)))) < TOL);
        assertTrue(relerr(moment_joint_upfactorial_from_raw(m, dims),
                arr(moment_upfactorial_from_raw(vec(m)))) < TOL);
        assertTrue(relerr(moment_joint_central_from_raw(m, dims),
                arr(moment_central_from_raw(vec(m)))) < TOL);
    }

    @Test
    public void testJointCentralMomentsAndTheCovariance() {
        int[] dims = {4, 4};
        double[] m = jointOracle(LAW_P, LAW_V, dims, 0);
        double[] mc = moment_joint_central_from_raw(m, dims);
        // row-major: index (i,j) is i*4+j
        assertTrue(Math.abs(mc[1 * 4 + 1] - (m[1 * 4 + 1] - m[1 * 4] * m[1])) < TOL);
        assertTrue(Math.abs(mc[1 * 4]) < TOL && Math.abs(mc[1]) < TOL);
        double[] mu = {m[1 * 4], m[1]};
        assertTrue(relerr(moment_joint_raw_from_central(mc, dims, mu), m) < TOL);
        assertTrue(relerr(moment_joint_central_from_raw_mean(m, dims, mu), mc) < TOL);
    }

    @Test
    public void testJointCumulantsVanishUnderIndependence() {
        // for independent components the joint moment array is an outer
        // product and every mixed cumulant is zero: no separable transform has
        // this property, which is why the cumulant edge is not a Kronecker one
        int[] dims = {4, 4};
        double[] m1 = arr(moment_raw_from_cumulant(vec(new double[]{0, 2, 2, 2})));
        double[] m2 = arr(moment_raw_from_cumulant(vec(new double[]{0, 1, 3, 5})));
        double[] m = new double[16];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                m[i * 4 + j] = m1[i] * m2[j];
            }
        }
        double[] k = moment_joint_cumulant_from_raw(m, dims);
        for (int i = 1; i < 4; i++) {
            for (int j = 1; j < 4; j++) {
                assertTrue(Math.abs(k[i * 4 + j]) < 1e-8, "mixed cumulant " + i + "," + j);
            }
        }
        assertTrue(relerr(moment_joint_raw_from_cumulant(k, dims), m) < TOL);
    }

    @Test
    public void testJointCumulantRoundtripAndCovariance() {
        int[] dims = {4, 4};
        double[] m = jointOracle(LAW_P, LAW_V, dims, 0);
        double[] f = jointOracle(LAW_P, LAW_V, dims, 1);
        double[] k = moment_joint_cumulant_from_raw(m, dims);
        assertTrue(Math.abs(k[1 * 4 + 1] - (m[1 * 4 + 1] - m[1 * 4] * m[1])) < TOL);
        assertTrue(relerr(moment_joint_raw_from_cumulant(k, dims), m) < TOL);
        double[] kf = moment_joint_factcumulant_from_factorial(f, dims);
        assertTrue(relerr(moment_joint_factorial_from_factcumulant(kf, dims), f) < TOL);
    }

    @Test
    public void testJointConversionsInThreeDimensions() {
        int[] dims = {3, 3, 3};
        double[] p = {0.5, 0.3, 0.2};
        double[][] v = {{1, 2, 0}, {3, 0, 1}, {2, 1, 4}};
        double[] m = jointOracle(p, v, dims, 0);
        double[] f = jointOracle(p, v, dims, 1);
        double[] fp = jointOracle(p, v, dims, 2);
        assertTrue(relerr(moment_joint_factorial_from_raw(m, dims), f) < TOL);
        assertTrue(relerr(moment_joint_upfactorial_from_raw(m, dims), fp) < TOL);
        double[] k = moment_joint_cumulant_from_raw(m, dims);
        assertTrue(relerr(moment_joint_raw_from_cumulant(k, dims), m) < TOL);
    }

    @Test
    public void testMarkingOfAPoissonGivesIndependentPoissons() {
        double lam = 3.0;
        double[] f = new double[7];
        for (int n = 0; n < 7; n++) {
            f[n] = Math.pow(lam, n);
        }
        double[] p = {0.25, 0.75};
        int[] ord = {3, 3};
        int[] dims = {4, 4};
        double[] F = moment_joint_marking(vec(f), p, ord);
        double[] kf = moment_joint_factcumulant_from_factorial(F, dims);
        assertTrue(Math.abs(kf[1 * 4] - lam * p[0]) < TOL);
        assertTrue(Math.abs(kf[1] - lam * p[1]) < TOL);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (i + j >= 2) {
                    assertTrue(Math.abs(kf[i * 4 + j]) < 1e-8, "cumulant " + i + "," + j);
                }
            }
        }
    }

    @Test
    public void testMarkingAndAggregationAreInverse() {
        double[] f = {1, 2, 4.5, 11, 30};
        double[] p = {0.4, 0.6};
        double[] F = moment_joint_marking(vec(f), p, new int[]{2, 2});
        assertTrue(relerr(moment_joint_aggregate(F, new int[]{3, 3}),
                new double[]{f[0], f[1], f[2]}) < TOL);
    }

    @Test
    public void testAggregationHoldsForADependentJointLaw() {
        // the aggregation identity is Vandermonde, so it needs no independence
        int[] dims = {5, 5};
        double[] f = jointOracle(LAW_P, LAW_V, dims, 1);
        double[][] tot = {{3}, {3}, {7}, {1}};
        double[] fs = jointOracle(LAW_P, tot, new int[]{5}, 1);
        assertTrue(relerr(moment_joint_aggregate(f, dims), fs) < TOL);
    }

    @Test
    public void testHouseMatrixAndTensorTransform() {
        Matrix T = moment_housematrix("factorial_from_raw", 4);
        assertTrue(relerr(T, arr(moment_stirling1(4))) < TOL);
        int[] dims = {4, 5};
        double[] m = new double[20];
        for (int i = 0; i < 20; i++) {
            m[i] = i + 1.0;
        }
        double[] got = moment_tensortrans(m, dims,
                moment_housematrix("raw_from_factorial", 4), 1);
        for (int i = 0; i < 4; i++) {
            double[] row = new double[5];
            System.arraycopy(m, i * 5, row, 0, 5);
            double[] exp = arr(moment_raw_from_factorial(vec(row)));
            for (int j = 0; j < 5; j++) {
                assertTrue(Math.abs(got[i * 5 + j] - exp[j]) < TOL, "row " + i);
            }
        }
        assertThrows(IllegalArgumentException.class,
                () -> moment_tensortrans(m, dims, T, 5));
        assertThrows(IllegalArgumentException.class,
                () -> moment_housematrix("no_such_edge", 3));
        assertThrows(IllegalArgumentException.class,
                () -> moment_joint_central_from_raw(new double[]{1, 2, 3}, new int[]{3, 1}));
        assertThrows(IllegalArgumentException.class,
                () -> moment_joint_marking(vec(new double[]{1, 2}), new double[]{0.5, 0.5},
                        new int[]{2, 2}));
    }
}
