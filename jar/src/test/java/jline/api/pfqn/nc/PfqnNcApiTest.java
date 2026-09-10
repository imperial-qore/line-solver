/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.nc;

import jline.VerboseLevel;
import jline.io.Ret;
import jline.lang.constant.SolverType;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validation of the normalizing-constant algorithms (jline.api.pfqn.nc) on
 * single-class product-form networks whose normalizing constant is computed
 * inline by exact convolution (Buzen). Exact algorithms (CA, RECAL, COMOM)
 * must reproduce it to numerical precision; asymptotic and Monte Carlo
 * methods must approach it within their documented accuracy regimes.
 */
public class PfqnNcApiTest {

    private static final double EXACT_TOL = 1e-8;

    // Test network: two queueing stations, one delay; single class
    private static final double L1 = 0.5, L2 = 0.4, Z0 = 0.3;

    /**
     * Exact log normalizing constant by direct convolution:
     * G(N) = sum_{k0+k1+k2=N} Z^k0/k0! * L1^k1 * L2^k2.
     */
    private static double exactLogG(double l1, double l2, double z, int n) {
        double g = 0.0;
        for (int k0 = 0; k0 <= n; k0++) {
            double zTerm = Math.pow(z, k0);
            for (int f = 2; f <= k0; f++) {
                zTerm /= f;
            }
            for (int k1 = 0; k1 + k0 <= n; k1++) {
                int k2 = n - k0 - k1;
                g += zTerm * Math.pow(l1, k1) * Math.pow(l2, k2);
            }
        }
        return Math.log(g);
    }

    private static Matrix rowVector(double... values) {
        Matrix m = new Matrix(1, values.length);
        for (int i = 0; i < values.length; i++) {
            m.set(0, i, values[i]);
        }
        return m;
    }

    private static Matrix colVector(double... values) {
        Matrix m = new Matrix(values.length, 1);
        for (int i = 0; i < values.length; i++) {
            m.set(i, 0, values[i]);
        }
        return m;
    }

    @Test
    public void convolutionAlgorithmIsExact() {
        int n = 4;
        Matrix L = colVector(L1, L2);
        Ret.pfqnNc r = Pfqn_ca.pfqn_ca(L, rowVector(n), rowVector(Z0));
        assertEquals(exactLogG(L1, L2, Z0, n), r.lG, EXACT_TOL,
                "CA lG vs inline convolution");
    }

    @Test
    public void recalMatchesConvolution() {
        int n = 4;
        Matrix L = colVector(L1, L2);
        Ret.pfqnNc r = Pfqn_recal.pfqn_recal(L, rowVector(n), rowVector(Z0));
        assertEquals(exactLogG(L1, L2, Z0, n), r.lG, EXACT_TOL,
                "RECAL lG vs inline convolution");
    }

    // ------------------------------------------------------------------
    // Multiclass: exact reference by direct state-space enumeration
    // ------------------------------------------------------------------

    // Two stations x two classes with think times
    private static final double[][] LM = {{0.5, 0.3}, {0.4, 0.6}};
    private static final double[] ZM = {0.3, 0.2};
    private static final int[] NM = {2, 1};

    private static double binom(int n, int k) {
        double b = 1.0;
        for (int i = 1; i <= k; i++) {
            b = b * (n - k + i) / i;
        }
        return b;
    }

    /**
     * Exact multiclass normalizing constant by convolution over the two
     * single-server stations and the delay:
     * G = sum over splits of (n1,n2) into (delay, st1, st2) of
     *     Z1^d1/d1! Z2^d2/d2! * prod_m (k1m+k2m)!/(k1m! k2m!) L(m,1)^k1m L(m,2)^k2m.
     */
    private static double exactMulticlassLogG() {
        double g = 0.0;
        for (int d1 = 0; d1 <= NM[0]; d1++) {
            for (int d2 = 0; d2 <= NM[1]; d2++) {
                double zTerm = Math.pow(ZM[0], d1) * Math.pow(ZM[1], d2);
                for (int f = 2; f <= d1; f++) zTerm /= f;
                for (int f = 2; f <= d2; f++) zTerm /= f;
                int r1 = NM[0] - d1, r2 = NM[1] - d2;
                for (int a1 = 0; a1 <= r1; a1++) {
                    for (int a2 = 0; a2 <= r2; a2++) {
                        int b1 = r1 - a1, b2 = r2 - a2;
                        double st1 = binom(a1 + a2, a1)
                                * Math.pow(LM[0][0], a1) * Math.pow(LM[0][1], a2);
                        double st2 = binom(b1 + b2, b1)
                                * Math.pow(LM[1][0], b1) * Math.pow(LM[1][1], b2);
                        g += zTerm * st1 * st2;
                    }
                }
            }
        }
        return Math.log(g);
    }

    private static Matrix multiclassL() {
        Matrix l = new Matrix(2, 2);
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                l.set(i, j, LM[i][j]);
            }
        }
        return l;
    }

    @ParameterizedTest(name = "ncMethod={0}")
    @ValueSource(strings = {"exact", "ca"})
    public void dispatcherExactMethodsMatchEnumeration(String method) {
        SolverOptions options = new SolverOptions(SolverType.NC);
        options.verbose = VerboseLevel.SILENT;
        options.method = method;
        Ret.pfqnNcXQ r = Pfqn_nc.pfqn_nc(Matrix.zeros(1, 2), multiclassL(),
                rowVector(NM[0], NM[1]), rowVector(ZM[0], ZM[1]), options);
        assertEquals(exactMulticlassLogG(), r.lG, 1e-6,
                "pfqn_nc method '" + method + "' vs exact enumeration");
    }

    @Test
    public void dispatcherDefaultMethodIsNearExact() {
        // The adaptive default picks a cubature whose error is far below 1e-3
        SolverOptions options = new SolverOptions(SolverType.NC);
        options.verbose = VerboseLevel.SILENT;
        options.method = "default";
        Ret.pfqnNcXQ r = Pfqn_nc.pfqn_nc(Matrix.zeros(1, 2), multiclassL(),
                rowVector(NM[0], NM[1]), rowVector(ZM[0], ZM[1]), options);
        assertEquals(exactMulticlassLogG(), r.lG, 1e-3,
                "pfqn_nc default method vs exact enumeration");
    }

    @Test
    public void dispatcherComomHandlesRepairmanModels() {
        // The dispatcher's comom branch covers the repairman case (one
        // queueing station plus delay) for multiclass models
        Matrix l = new Matrix(1, 2);
        l.set(0, 0, LM[0][0]);
        l.set(0, 1, LM[0][1]);
        // Exact: single station + delay by enumeration
        double g = 0.0;
        for (int d1 = 0; d1 <= NM[0]; d1++) {
            for (int d2 = 0; d2 <= NM[1]; d2++) {
                double zTerm = Math.pow(ZM[0], d1) * Math.pow(ZM[1], d2);
                for (int f = 2; f <= d1; f++) zTerm /= f;
                for (int f = 2; f <= d2; f++) zTerm /= f;
                int a1 = NM[0] - d1, a2 = NM[1] - d2;
                g += zTerm * binom(a1 + a2, a1)
                        * Math.pow(LM[0][0], a1) * Math.pow(LM[0][1], a2);
            }
        }
        SolverOptions options = new SolverOptions(SolverType.NC);
        options.verbose = VerboseLevel.SILENT;
        options.method = "comom";
        Ret.pfqnNcXQ r = Pfqn_nc.pfqn_nc(Matrix.zeros(1, 2), l,
                rowVector(NM[0], NM[1]), rowVector(ZM[0], ZM[1]), options);
        assertEquals(Math.log(g), r.lG, 1e-4,
                "pfqn_nc comom on repairman model vs exact enumeration");
    }

    @ParameterizedTest(name = "ncApproxMethod={0}")
    @ValueSource(strings = {"ls", "sampling"})
    public void dispatcherApproximateMethodsAreClose(String method) {
        SolverOptions options = new SolverOptions(SolverType.NC);
        options.verbose = VerboseLevel.SILENT;
        options.method = method;
        Ret.pfqnNcXQ r = Pfqn_nc.pfqn_nc(Matrix.zeros(1, 2), multiclassL(),
                rowVector(NM[0], NM[1]), rowVector(ZM[0], ZM[1]), options);
        assertEquals(exactMulticlassLogG(), r.lG, 0.35,
                "pfqn_nc approximate method '" + method + "' deviates grossly");
    }

    @Test
    public void logisticExpansionAccurateAtLargePopulation() {
        // The logistic expansion carries an O(1/N) asymptotic error; at N=64
        // the MATLAB reference itself deviates by |dlG| ~ 0.09 from exact
        int n = 64;
        Matrix L = colVector(L1, L2);
        double exact = exactLogG(L1, L2, 0.0, n);
        Ret.pfqnNc le = Pfqn_le.pfqn_le(L, rowVector(n));
        assertEquals(exact, le.lG, 0.2,
                "LE asymptotic lG beyond its O(1/N) error regime at N=64");
    }

    @Test
    public void logisticSamplingWithinStatisticalTolerance() {
        int n = 8;
        Matrix L = colVector(L1, L2);
        double exact = exactLogG(L1, L2, Z0, n);
        Ret.pfqnNc ls = Pfqn_ls.pfqn_ls(L, rowVector(n), rowVector(Z0), 200000L, 23000L);
        assertEquals(exact, ls.lG, 0.10,
                "LS Monte Carlo lG within 10% of exact");
    }

    @Test
    public void mcintegrationWithinStatisticalTolerance() {
        int n = 8;
        Matrix L = colVector(L1, L2);
        double exact = exactLogG(L1, L2, Z0, n);
        Ret.pfqnNc mci = Pfqn_mci.pfqn_mci(L, rowVector(n), rowVector(Z0), 200000, "imci");
        assertEquals(exact, mci.lG, 0.10,
                "MCI lG within 10% of exact");
    }

    @Test
    public void mmintegrationMatchesRepairmanClosedForm() {
        // McKenna-Mitra integral applies to the repairman model: ONE queueing
        // station with per-class demands plus per-class think times.
        // Single class: G = sum_{k} Z^{N-k}/(N-k)! L^k, computed inline.
        int n = 8;
        double g = 0.0;
        for (int k = 0; k <= n; k++) {
            double term = Math.pow(L1, k) * Math.pow(Z0, n - k);
            for (int f = 2; f <= n - k; f++) {
                term /= f;
            }
            g += term;
        }
        double exact = Math.log(g);
        Ret.pfqnNc mm = Pfqn_mmint2.pfqn_mmint2(rowVector(L1), rowVector(n), rowVector(Z0));
        assertEquals(exact, mm.lG, 1e-4,
                "McKenna-Mitra integration vs repairman closed form");
    }

    @Test
    public void propfairThroughputNearExact() {
        int n = 8;
        Matrix L = colVector(L1, L2);
        // Exact throughput from the normalizing constant ratio X = G(N-1)/G(N)
        double exactX = Math.exp(exactLogG(L1, L2, Z0, n - 1) - exactLogG(L1, L2, Z0, n));
        Ret.pfqnNcXQ pf = Pfqn_propfair.pfqn_propfair(L, rowVector(n), rowVector(Z0));
        assertTrue(pf.X != null && pf.X.getNumElements() > 0,
                "propfair must return a throughput estimate");
        assertEquals(exactX, pf.X.get(0), 0.15 * exactX,
                "proportionally-fair throughput within 15% of exact");
    }

}
