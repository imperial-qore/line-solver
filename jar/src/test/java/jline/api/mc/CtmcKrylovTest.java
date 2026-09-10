package jline.api.mc;

import org.junit.jupiter.api.Test;

import jline.util.matrix.Matrix;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The two Krylov kernels of the CTMC path, against an ANALYTICAL oracle.
 *
 * {@link Ctmc_gmres} and {@link Ctmc_bicgstab} share the same equilibration, reverse
 * Cuthill-McKee reordering and ILUT preconditioner, so checking one against the other
 * proves nothing about either: a fault in the shared preparation would cancel. The oracle
 * here is the stationary vector of M/M/1/K, which is the truncated geometric in closed
 * form.
 *
 * Twin of python/tests/test_ctmc_krylov.py and of the ctmc_gmres/ctmc_bicgstab cases in
 * cpp/tests/test_mc_aggregation.cpp.
 */
public class CtmcKrylovTest {

    private static final double LAMBDA = 0.7;
    private static final double MU = 1.0;

    /** The system Ctmc_solve builds for M/M/1/K: last column of Q replaced by ones, transposed. */
    private static Matrix mm1kSystem(int K) {
        int n = K + 1;
        Matrix Q = new Matrix(n, n);
        for (int i = 0; i + 1 < n; i++) {
            Q.set(i, i + 1, LAMBDA);
            Q.set(i + 1, i, MU);
        }
        for (int i = 0; i < n; i++) {
            double off = 0.0;
            if (i + 1 < n) off += LAMBDA;
            if (i > 0) off += MU;
            Q.set(i, i, -off);
        }
        for (int i = 0; i < n; i++) Q.set(i, n - 1, 1.0);
        return Q.transpose();
    }

    private static Matrix unitLastEntry(int n) {
        Matrix b = new Matrix(n, 1);
        b.set(n - 1, 0, 1.0);
        return b;
    }

    private static double[] truncatedGeometric(int n) {
        double[] exact = new double[n];
        double total = 0.0;
        double power = 1.0;
        for (int i = 0; i < n; i++) {
            exact[i] = power;
            total += power;
            power *= LAMBDA / MU;
        }
        for (int i = 0; i < n; i++) exact[i] /= total;
        return exact;
    }

    private static double maxDeviation(Matrix x, double[] exact) {
        double dev = 0.0;
        for (int i = 0; i < exact.length; i++) dev = Math.max(dev, Math.abs(x.get(i, 0) - exact[i]));
        return dev;
    }

    @Test
    public void testGmresReproducesTheTruncatedGeometric() {
        int K = 2000;
        Matrix A = mm1kSystem(K);
        Ctmc_gmres.GmresResult r = Ctmc_gmres.ctmc_gmres(A, unitLastEntry(K + 1));
        assertEquals(0, r.flag, "GMRES did not converge, relative residual " + r.relres);
        assertTrue(maxDeviation(r.x, truncatedGeometric(K + 1)) < 1e-11,
                "GMRES deviates from the closed form by " + maxDeviation(r.x, truncatedGeometric(K + 1)));
    }

    @Test
    public void testBicgstabReproducesTheTruncatedGeometric() {
        int K = 2000;
        Matrix A = mm1kSystem(K);
        Ctmc_bicgstab.BicgstabResult r = Ctmc_bicgstab.ctmc_bicgstab(A, unitLastEntry(K + 1));
        assertEquals(0, r.flag, "BiCGSTAB did not converge, relative residual " + r.relres);
        assertTrue(maxDeviation(r.x, truncatedGeometric(K + 1)) < 1e-11,
                "BiCGSTAB deviates from the closed form by " + maxDeviation(r.x, truncatedGeometric(K + 1)));
    }

    @Test
    public void testBothKernelsAgreeOnTheSameSystem() {
        int K = 500;
        Matrix A = mm1kSystem(K);
        Matrix b = unitLastEntry(K + 1);
        Ctmc_gmres.GmresResult g = Ctmc_gmres.ctmc_gmres(A, b);
        Ctmc_bicgstab.BicgstabResult s = Ctmc_bicgstab.ctmc_bicgstab(A, b);
        assertEquals(0, g.flag);
        assertEquals(0, s.flag);
        for (int i = 0; i <= K; i++) {
            assertEquals(g.x.get(i, 0), s.x.get(i, 0), 1e-10, "the two kernels disagree at state " + i);
        }
    }

    @Test
    public void testBicgstabReportsMatvecsNotIterations() {
        // iter counts matrix-vector products with A -- two per complete iteration, and
        // one when the iteration converges at its half step, so an ODD count is normal
        // and only the bound is assertable. Counting products rather than iterations is
        // what makes the number comparable with Ctmc_gmres and across the four codebases.
        int maxit = 10;
        Matrix A = mm1kSystem(200);
        Ctmc_bicgstab.BicgstabResult r = Ctmc_bicgstab.ctmc_bicgstab(A, unitLastEntry(201), 0.0, maxit, null);
        assertEquals(0, r.flag);
        assertTrue(r.iter >= 1, "a converged solve performs at least one product");
        assertTrue(r.iter <= 2 * maxit, "iter=" + r.iter + " exceeds two products per iteration");
    }

    @Test
    public void testMultiColumnBicgstabSolvesEveryColumn() {
        // The shape of the stochastic complement: one factorization, many columns.
        int K = 300;
        int n = K + 1;
        Matrix A = mm1kSystem(K);
        Matrix B = new Matrix(n, 3);
        for (int i = 0; i < n; i++) {
            B.set(i, 0, 1.0);
            B.set(i, 1, i % 2 == 0 ? 1.0 : -1.0);
            B.set(i, 2, (double) i / n);
        }
        Matrix X = Ctmc_bicgstab.ctmc_bicgstab(A, B, 1e-10, 0);
        assertNotNull(X, "the multi-column solve reported failure");
        Matrix R = A.mult(X);
        for (int i = 0; i < n; i++) {
            for (int c = 0; c < 3; c++) {
                assertEquals(B.get(i, c), R.get(i, c), 1e-7, "column " + c + " does not satisfy its equation");
            }
        }
    }
}
