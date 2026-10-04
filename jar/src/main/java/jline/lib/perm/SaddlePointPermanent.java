package jline.lib.perm;

import jline.util.matrix.Matrix;

/**
 * Saddle-point (SPM) approximation of the permanent of a positive matrix.
 *
 * Twin of MATLAB perm_spm.m and of the native Python SaddlePointPermanent.
 * This is the HOMOGENEOUS variant of cache_spm: both evaluate the same Cauchy
 * integral by Laplace's method and differ only in the generating function whose
 * coefficient they extract,
 *
 * <pre>
 *   cache_spm  E(m) = prod_l m_l! [prod_l z_l^m_l] prod_k (1 + sum_l g_kl z_l)
 *   perm_spm   P    = prod_l m_l! [prod_l z_l^m_l] prod_k (    sum_l A_kl z_l)
 * </pre>
 *
 * The cache factor carries a "+1" because an item may stay out of the cache, so
 * the coefficient it extracts is a rectangular permanent over n items and
 * sum(m) &lt; n slots. Dropping the "+1" forces every row to be matched, which is
 * exactly the permanent and requires sum(m) == n. That is the one case
 * cache_spm cannot serve: at n == sum(m) its multipliers diverge and it falls
 * back on cache_erec. Here the integrand is homogeneous, and the saddle point
 * is interior in the h-1 directions that survive.
 *
 * METHOD. With z_l = xi_l exp(i th_l) the saddle point in xi solves
 *
 * <pre>
 *   sum_k A_kl xi_l / (sum_j A_kj xi_j) = m_l,   l = 1..h,
 * </pre>
 *
 * that is, P_kl = A_kl xi_l / S_k with S = A*xi is the diagonal scaling of A to
 * row sums 1 and column sums m (Sinkhorn scaling; doubly stochastic when m is
 * all ones). There phi = sum_k log S_k - sum_l m_l log xi_l is the log of the
 * Gurvits capacity, an upper bound on the log permanent. The Gaussian
 * correction uses H = diag(m) - P'P, a weighted graph Laplacian on the columns:
 * H*ones = 0, which is the invariance of the integrand under th -&gt; th + c*ones
 * that homogeneity creates. That direction is a full period rather than a
 * Gaussian, so it contributes 2*pi and leaves an (h-1)-dimensional Laplace
 * integral. Any principal (h-1) submatrix serves, since all cofactors of a
 * Laplacian are equal, and
 *
 * <pre>
 *   log P = sum_l log(m_l!) - (h-1)/2 log(2 pi) + phi - 1/2 log det(H_red).
 * </pre>
 *
 * ACCURACY. Exact for h == 1, where the permanent is n! prod_k A(k,0). It is a
 * genuine asymptotic expansion as min(m) grows with h fixed, the ratio to the
 * exact permanent falling from 1.11 at m = (2,2,2) to 1.02 at m = (3,3). At
 * m = ones the dimension of the integral grows with the expansion parameter and
 * the leading term keeps a systematic bias: on the n x n matrix of ones it
 * returns (2 pi)^(-(n-1)/2) n^(n+1/2) against the exact n!, a ratio tending to
 * (e/sqrt(2 pi))^n = 1.084^n, and random positive matrices track that closely
 * (1.31 at n = 4, 1.87 at n = 8). So at m = ones it OVERESTIMATES, with a spread
 * across matrices far tighter than the bias itself, and it is not a bound in
 * either direction. BethePermanent is a genuine lower bound.
 */
public class SaddlePointPermanent extends PermSolver {

    private final int[] m;
    private final double tolerance;
    private final int maxIterations;
    private double logValue = 0.0;
    private double logCapacity = 0.0;
    private double[] xi;

    /**
     * Constructor with all parameters.
     *
     * @param matrix        n x h strictly positive matrix
     * @param m             column multiplicities, non-negative and summing to n;
     *                      null means all ones, which requires a square matrix
     * @param tolerance     margin on the column sums at which the scaling stops
     * @param maxIterations maximum number of scaling sweeps
     * @param solve         whether to run solve() after construction
     * @throws IllegalArgumentException if the matrix is negative or not strictly
     *                      positive, or the multiplicities do not sum to the row count
     */
    public SaddlePointPermanent(Matrix matrix, int[] m, double tolerance, int maxIterations,
                                boolean solve) {
        super(matrix);
        this.tolerance = tolerance;
        this.maxIterations = maxIterations;
        final int rows = matrix.getNumRows();
        final int cols = matrix.getNumCols();
        if (rows == 0 || cols == 0) {
            this.m = new int[0];
            this.xi = new double[0];
            this.value = 1.0;   // the permanent of the empty matrix is 1
            return;
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (matrix.get(i, j) < 0.0) {
                    throw new IllegalArgumentException("Matrix must be non-negative. Found negative"
                            + " element at (" + i + ", " + j + "): " + matrix.get(i, j));
                }
            }
        }
        int[] mult = m;
        if (mult == null) {
            if (cols != rows) {
                throw new IllegalArgumentException("Without column multiplicities the matrix must"
                        + " be square; it is " + rows + "x" + cols + ".");
            }
            mult = new int[rows];
            for (int j = 0; j < rows; j++) {
                mult[j] = 1;
            }
        }
        if (mult.length != cols) {
            throw new IllegalArgumentException("The multiplicity vector has " + mult.length
                    + " entries against " + cols + " columns.");
        }
        int total = 0;
        for (int j = 0; j < cols; j++) {
            if (mult[j] < 0) {
                throw new IllegalArgumentException("Column multiplicities must be non-negative"
                        + " integers; entry " + j + " is " + mult[j] + ".");
            }
            total += mult[j];
        }
        if (total != rows) {
            throw new IllegalArgumentException("The column multiplicities must sum to the number"
                    + " of rows (sum(m) = " + total + " against " + rows + " rows). The integrand"
                    + " is homogeneous of degree " + rows + ", so every other coefficient of it is"
                    + " exactly zero.");
        }
        PermSupport.requireFullSupport(matrix, "spm");
        this.m = mult.clone();
        this.xi = new double[cols];

        if (solve) {
            solve();
        }
    }

    /** Multiplicities and defaults for tolerance (1e-11) and maxIterations (10000). */
    public SaddlePointPermanent(Matrix matrix, int[] m, boolean solve) {
        this(matrix, m, 1e-11, 10000, solve);
    }

    /** Unit multiplicities, which requires a square matrix. */
    public SaddlePointPermanent(Matrix matrix, boolean solve) {
        this(matrix, null, 1e-11, 10000, solve);
    }

    /** Unit multiplicities, not solved on construction. */
    public SaddlePointPermanent(Matrix matrix) {
        this(matrix, null, 1e-11, 10000, false);
    }

    /** Logarithm of the estimate, correct even when the estimate overflows. */
    public double getLogValue() {
        return logValue;
    }

    /** Log Gurvits capacity at the saddle point, an upper bound on the log permanent. */
    public double getLogCapacity() {
        return logCapacity;
    }

    /** Saddle point, unit geometric mean, zero on a column of multiplicity zero. */
    public double[] getXi() {
        return xi;
    }

    @Override
    public void compute() {
        if (matrix.getNumRows() == 0 || matrix.getNumCols() == 0) {
            value = 1.0;
            logValue = 0.0;
            return;
        }
        expand();
        value = Math.exp(logValue);
    }

    /** Locate the saddle point and evaluate the Laplace expansion there. */
    private void expand() {
        final int rows = matrix.getNumRows();
        final int cols = matrix.getNumCols();

        // A column repeated zero times leaves the permanent unchanged, and its xi is a
        // boundary of the Laplace integral rather than a direction of it, so it must
        // leave the expansion. Dropping it is exact: setting z_l = 0 removes the column,
        // and prod_l m_l! is unchanged because 0! = 1.
        int hk = 0;
        for (int j = 0; j < cols; j++) {
            if (m[j] > 0) {
                hk++;
            }
        }
        int[] keep = new int[hk];
        double[] mk = new double[hk];
        int at = 0;
        for (int j = 0; j < cols; j++) {
            if (m[j] > 0) {
                keep[at] = j;
                mk[at] = m[j];
                at++;
            }
        }
        double[][] a = new double[rows][hk];
        for (int i = 0; i < rows; i++) {
            for (int l = 0; l < hk; l++) {
                a[i][l] = matrix.get(i, keep[l]);
            }
        }

        double[] xik = scale(a, mk);
        xi = new double[cols];
        for (int l = 0; l < hk; l++) {
            xi[keep[l]] = xik[l];
        }

        double[] s = rowScaled(a, xik);
        double[][] p = new double[rows][hk];
        for (int i = 0; i < rows; i++) {
            for (int l = 0; l < hk; l++) {
                p[i][l] = a[i][l] * xik[l] / s[i];
            }
        }
        logCapacity = 0.0;
        for (int i = 0; i < rows; i++) {
            logCapacity += Math.log(s[i]);
        }
        for (int l = 0; l < hk; l++) {
            logCapacity -= mk[l] * Math.log(xik[l]);
        }

        // H = diag(mk) - P'P is a Laplacian, so it is singular along ones and all of its
        // principal cofactors are equal; the last index is dropped only because one has
        // to be. Strict positivity of A makes the column graph complete, hence H_red
        // positive definite and Cholesky the right factor.
        double logDet = 0.0;
        if (hk > 1) {
            double[][] hred = new double[hk - 1][hk - 1];
            for (int l = 0; l < hk - 1; l++) {
                for (int j = 0; j < hk - 1; j++) {
                    double dot = 0.0;
                    for (int i = 0; i < rows; i++) {
                        dot += p[i][l] * p[i][j];
                    }
                    hred[l][j] = (l == j ? mk[l] : 0.0) - dot;
                }
            }
            logDet = logDetCholesky(hred);
        }
        // no direction survives the homogeneity when hk == 1, and det of the empty matrix is 1

        double logFact = 0.0;
        for (int l = 0; l < hk; l++) {
            logFact += logFactorial((int) mk[l]);
        }
        logValue = logFact - 0.5 * (hk - 1) * Math.log(2.0 * Math.PI) + logCapacity - 0.5 * logDet;
    }

    /** Scale a to row sums 1 and column sums mk, returning the multipliers. */
    private double[] scale(double[][] a, double[] mk) {
        final int rows = a.length;
        final int hk = mk.length;
        double[] xik = new double[hk];
        for (int l = 0; l < hk; l++) {
            xik[l] = 1.0;
        }
        double margin = Double.POSITIVE_INFINITY;
        for (int it = 0; it < maxIterations; it++) {
            double[] s = rowScaled(a, xik);
            double[] colsum = new double[hk];
            for (int l = 0; l < hk; l++) {
                double acc = 0.0;
                for (int i = 0; i < rows; i++) {
                    acc += a[i][l] / s[i];
                }
                colsum[l] = xik[l] * acc;
            }
            margin = 0.0;
            for (int l = 0; l < hk; l++) {
                margin = Math.max(margin, Math.abs(colsum[l] - mk[l]));
            }
            if (margin < tolerance) {
                return xik;
            }
            double logmean = 0.0;
            for (int l = 0; l < hk; l++) {
                xik[l] = xik[l] * mk[l] / colsum[l];
                logmean += Math.log(xik[l]);
            }
            logmean /= hk;
            for (int l = 0; l < hk; l++) {
                xik[l] /= Math.exp(logmean);    // the saddle is a ray; pin its scale
            }
        }
        throw new IllegalArgumentException("The scaling to row sums 1 and column sums m did not"
                + " converge in " + maxIterations + " sweeps (margin error " + margin
                + " against a tolerance of " + tolerance + "). The expansion assumes the saddle"
                + " point, so no value is returned. The usual cause is a matrix without total"
                + " support.");
    }

    /** Row sums of a scaled by xi, that is S(i) = sum_l a(i,l) xi(l). */
    private static double[] rowScaled(double[][] a, double[] xik) {
        double[] s = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            double acc = 0.0;
            for (int l = 0; l < xik.length; l++) {
                acc += a[i][l] * xik[l];
            }
            s[i] = acc;
        }
        return s;
    }

    /** Log determinant of a symmetric positive definite matrix, by Cholesky. */
    private static double logDetCholesky(double[][] h) {
        final int k = h.length;
        double[][] l = new double[k][k];
        double logdet = 0.0;
        for (int i = 0; i < k; i++) {
            for (int j = 0; j <= i; j++) {
                double acc = h[i][j];
                for (int t = 0; t < j; t++) {
                    acc -= l[i][t] * l[j][t];
                }
                if (i == j) {
                    if (!(acc > 0.0)) {
                        throw new IllegalArgumentException("The reduced Hessian is not positive"
                                + " definite (leading minor " + (i + 1) + "), so the saddle point"
                                + " is degenerate and the Gaussian factor does not exist.");
                    }
                    l[i][i] = Math.sqrt(acc);
                    logdet += 2.0 * Math.log(l[i][i]);
                } else {
                    l[i][j] = acc / l[j][j];
                }
            }
        }
        return logdet;
    }

    /**
     * log(n!). Summation of logs below 1000, which is more accurate than any series
     * there; past it Stirling with the 1/(12n) and 1/(360n^3) terms, good to 1e-14
     * relative and O(1) rather than O(n).
     */
    private static double logFactorial(int n) {
        if (n <= 1) {
            return 0.0;
        }
        if (n <= 1000) {
            double acc = 0.0;
            for (int i = 2; i <= n; i++) {
                acc += Math.log(i);
            }
            return acc;
        }
        final double d = n;
        return 0.5 * Math.log(2.0 * Math.PI * d) + d * Math.log(d) - d + 1.0 / (12.0 * d)
                - 1.0 / (360.0 * d * d * d);
    }
}
