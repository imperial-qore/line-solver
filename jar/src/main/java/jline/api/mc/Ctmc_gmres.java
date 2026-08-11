/**
 * Restarted GMRES with ILUT preconditioning for sparse linear systems
 *
 * Solves the sparse nonsymmetric system A*x = b by the restarted generalized minimal
 * residual method of Saad and Schultz (1986), right-preconditioned by the threshold
 * incomplete LU factorization ILUT of Saad (1994), with a Jacobi preconditioner as
 * fallback when the incomplete factorization breaks down. This is the iterative
 * counterpart of the direct sparse solve used by Ctmc_solve, intended for generators
 * whose LU fill-in exceeds available memory.
 *
 * Two preparation steps are not optional on a generator. Rows are equilibrated to unit
 * max norm, so the O(1) normalization row does not mix with rows carrying rates of a
 * different magnitude. The states are then reordered by reverse Cuthill-McKee: in the
 * natural ordering of a birth-death chain the unpivoted elimination has growth factor
 * (mu/lambda)^n, which overflows by a few thousand states, and a fill-reducing ordering
 * rather than pivoting is what removes it.
 *
 * A is the already-assembled coefficient matrix; no CTMC-specific processing is
 * performed here, so the same kernel serves the stochastic complementation and
 * aggregation kernels.
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import org.ejml.data.DMatrixSparseCSC;

import jline.util.matrix.Matrix;

public final class Ctmc_gmres {
    private Ctmc_gmres() {}

    /**
     * Default linear-solve residual. Much tighter than the fixed-point tolerance
     * options.iter_tol (1e-4): switching from the direct solve to GMRES must not move a
     * reported performance metric.
     */
    public static final double GMRES_DEFAULT_TOL = 1e-12;

    /** Default Krylov subspace dimension between restarts. */
    public static final int GMRES_DEFAULT_RESTART = 50;

    /** Relative threshold below which ILUT discards a fill-in entry. */
    public static final double ILUT_DROP_TOL = 1e-4;

    /** Factor bounding the ILUT factors, as a multiple of the nonzeros per row of A. */
    public static final double ILUT_FILL_FACTOR = 10.0;

    /** Arnoldi breakdown threshold, relative to the norm of the vector being expanded. */
    private static final double BREAKDOWN_TOL = 1e-14;

    /**
     * Outcome of a GMRES solve. The flag follows the MATLAB gmres convention: 0
     * converged, 1 iteration limit reached, 2 preconditioner ill-conditioned, 3
     * stagnation or breakdown. Callers must check it and fall back to the direct solve
     * when it is nonzero.
     */
    public static final class GmresResult {
        public final Matrix x;
        public final int flag;
        public final double relres;
        public final int iter;

        GmresResult(Matrix x, int flag, double relres, int iter) {
            this.x = x;
            this.flag = flag;
            this.relres = relres;
            this.iter = iter;
        }
    }

    /**
     * Solve A*x = b by restarted GMRES with the default tolerance, restart and iteration
     * cap, starting from a uniform initial guess.
     *
     * @param A Coefficient matrix
     * @param b Right-hand side, as a column vector
     * @return Solution, convergence flag, relative residual and iteration count
     */
    public static GmresResult ctmc_gmres(Matrix A, Matrix b) {
        return ctmc_gmres(A, b, GMRES_DEFAULT_TOL, 0, 0, null);
    }

    /**
     * Solve A*x = b by restarted GMRES.
     *
     * @param A       Coefficient matrix
     * @param b       Right-hand side, as a column vector
     * @param tol     Relative residual tolerance, nonpositive for the default 1e-12
     * @param restart Krylov subspace dimension, nonpositive for the default min(n,50)
     * @param maxit   Maximum number of restart cycles, nonpositive for ceil(n/restart)
     * @param x0      Initial guess, null for the uniform vector 1/n
     * @return Solution, convergence flag, relative residual and iteration count
     */
    public static GmresResult ctmc_gmres(Matrix A, Matrix b, double tol, int restart, int maxit, Matrix x0) {
        int n = A.getNumRows();
        if (b.length() != n) {
            throw new IllegalArgumentException("Matrix dimensions incompatible: A is " + n + "x" + n
                    + " but b has " + b.length() + " entries");
        }
        double[] rhs = new double[n];
        for (int i = 0; i < n; i++) rhs[i] = b.get(i);
        double[] guess = null;
        if (x0 != null) {
            guess = new double[n];
            for (int i = 0; i < n; i++) guess[i] = x0.get(i);
        }
        return new Prepared(A).solve(rhs, guess, tol, restart, maxit);
    }

    /**
     * Solve A*X = B for every column of B, reusing one incomplete factorization across
     * all of them and starting each column from the previous solution. This is the shape
     * of the stochastic complement, whose right-hand side is a whole block of the
     * generator: refactorizing per column would cost more than the direct solve it
     * replaces.
     *
     * @param A       Coefficient matrix
     * @param B       Right-hand sides, one per column
     * @param tol     Relative residual tolerance, nonpositive for the default 1e-12
     * @param restart Krylov subspace dimension, nonpositive for the default min(n,50)
     * @param maxit   Maximum number of restart cycles, nonpositive for ceil(n/restart)
     * @return The solution block, or null if any column failed to converge. Returning
     *         null rather than a partial block keeps the caller's fallback all-or-nothing.
     */
    public static Matrix ctmc_gmres(Matrix A, Matrix B, double tol, int restart, int maxit) {
        int n = A.getNumRows();
        if (B.getNumRows() != n) {
            throw new IllegalArgumentException("Matrix dimensions incompatible: A is " + n + "x" + n
                    + " but B has " + B.getNumRows() + " rows");
        }
        Prepared prepared = new Prepared(A);
        int nrhs = B.getNumCols();
        Matrix X = new Matrix(n, nrhs);
        double[] guess = null;
        for (int c = 0; c < nrhs; c++) {
            double[] rhs = new double[n];
            for (int i = 0; i < n; i++) rhs[i] = B.get(i, c);
            GmresResult r = prepared.solve(rhs, guess, tol, restart, maxit);
            if (r.flag != 0) return null;
            guess = new double[n];
            for (int i = 0; i < n; i++) {
                guess[i] = r.x.get(i, 0);
                if (guess[i] != 0.0) X.set(i, c, guess[i]);
            }
        }
        return X;
    }

    /**
     * The equilibrated, reordered and preconditioned form of a coefficient matrix. It is
     * built once and solved against any number of right-hand sides.
     */
    private static final class Prepared {
        final int n;
        final Csr csr;
        final int[] perm;
        final double[] rowScale;
        final Precond M;

        Prepared(Matrix A) {
            this.n = A.getNumRows();
            if (A.getNumCols() != n) {
                throw new IllegalArgumentException("Matrix A must be square for linear system solving");
            }
            // see _kb/03-api-layer.md for rationale
            Csr c = Csr.of(A.toDMatrixSparseCSC());
            this.rowScale = c.equilibrate();

            // Reverse Cuthill-McKee. perm maps a new index to the old one, iperm the reverse.
            this.perm = Rcm.order(c);
            int[] iperm = new int[n];
            for (int i = 0; i < n; i++) iperm[perm[i]] = i;
            this.csr = c.permuteSymmetric(perm, iperm);
            this.M = Precond.of(this.csr);
        }

        GmresResult solve(double[] rhsIn, double[] x0In, double tol, int restart, int maxit) {
            double[] scaled = new double[n];
            for (int i = 0; i < n; i++) scaled[i] = rhsIn[i] / rowScale[i];
            double[] x = new double[n];
            if (x0In == null) {
                java.util.Arrays.fill(x, 1.0 / n);
            } else {
                System.arraycopy(x0In, 0, x, 0, n);
            }
            return gmres(this, permute(scaled, perm), permute(x, perm), tol, restart, maxit);
        }
    }

    private static GmresResult gmres(Prepared prep, double[] rhs, double[] x,
                                     double tol, int restart, int maxit) {
        int n = prep.n;
        Csr csr = prep.csr;
        Precond M = prep.M;
        int[] perm = prep.perm;

        if (tol <= 0) tol = GMRES_DEFAULT_TOL;
        if (restart <= 0) restart = Math.min(n, GMRES_DEFAULT_RESTART);
        restart = Math.min(restart, n);
        if (maxit <= 0) maxit = (int) Math.ceil((double) n / restart);
        maxit = Math.max(1, Math.min(maxit, n));

        double bnorm = norm2(rhs);
        if (bnorm == 0.0) bnorm = 1.0;

        double[] r = new double[n];
        double[] w = new double[n];
        double[] z = new double[n];

        // r = b - A*x
        csr.mult(x, r);
        for (int i = 0; i < n; i++) r[i] = rhs[i] - r[i];
        double beta = norm2(r);
        if (beta / bnorm <= tol) {
            return new GmresResult(asColumn(unpermute(x, perm)), 0, beta / bnorm, 0);
        }

        double[][] V = new double[restart + 1][n];
        double[][] H = new double[restart + 1][restart];
        double[] cs = new double[restart];
        double[] sn = new double[restart];
        double[] g = new double[restart + 1];
        double[] y = new double[restart];

        int iter = 0;
        int flag = 1;
        double relres = beta / bnorm;
        double initres = relres;
        double prevrelres = Double.POSITIVE_INFINITY;

        for (int cycle = 0; cycle < maxit; cycle++) {
            beta = norm2(r);
            if (beta == 0.0) {
                flag = 0;
                relres = 0.0;
                break;
            }
            for (int i = 0; i < n; i++) V[0][i] = r[i] / beta;
            java.util.Arrays.fill(g, 0.0);
            g[0] = beta;

            int k = 0;
            for (int j = 0; j < restart; j++) {
                // see _kb/03-api-layer.md for rationale
                M.apply(V[j], z);
                csr.mult(z, w);
                iter++;

                double wnorm0 = norm2(w);
                // see _kb/03-api-layer.md for rationale
                for (int pass = 0; pass < 2; pass++) {
                    for (int i = 0; i <= j; i++) {
                        double hij = dot(V[i], w);
                        H[i][j] += hij;
                        for (int t = 0; t < n; t++) w[t] -= hij * V[i][t];
                    }
                }
                double hnext = norm2(w);
                H[j + 1][j] = hnext;

                k = j + 1;
                boolean breakdown = hnext <= BREAKDOWN_TOL * wnorm0;
                if (!breakdown) {
                    for (int t = 0; t < n; t++) V[j + 1][t] = w[t] / hnext;
                }

                // Apply the accumulated Givens rotations to the new Hessenberg column,
                // then annihilate its subdiagonal entry with a fresh rotation.
                for (int i = 0; i < j; i++) {
                    double t1 = cs[i] * H[i][j] + sn[i] * H[i + 1][j];
                    H[i + 1][j] = -sn[i] * H[i][j] + cs[i] * H[i + 1][j];
                    H[i][j] = t1;
                }
                double denom = Math.hypot(H[j][j], H[j + 1][j]);
                if (denom == 0.0) {
                    cs[j] = 1.0;
                    sn[j] = 0.0;
                } else {
                    cs[j] = H[j][j] / denom;
                    sn[j] = H[j + 1][j] / denom;
                }
                H[j][j] = cs[j] * H[j][j] + sn[j] * H[j + 1][j];
                H[j + 1][j] = 0.0;
                g[j + 1] = -sn[j] * g[j];
                g[j] = cs[j] * g[j];

                relres = Math.abs(g[j + 1]) / bnorm;
                if (relres <= tol || breakdown) break;
            }

            // Back-substitute the least-squares solution on the rotated Hessenberg
            // system, then map the correction back through the preconditioner.
            for (int i = k - 1; i >= 0; i--) {
                double s = g[i];
                for (int t = i + 1; t < k; t++) s -= H[i][t] * y[t];
                y[i] = H[i][i] == 0.0 ? 0.0 : s / H[i][i];
            }
            double[] corr = new double[n];
            for (int i = 0; i < k; i++) {
                for (int t = 0; t < n; t++) corr[t] += y[i] * V[i][t];
            }
            M.apply(corr, z);
            for (int t = 0; t < n; t++) x[t] += z[t];

            csr.mult(x, r);
            for (int t = 0; t < n; t++) r[t] = rhs[t] - r[t];
            relres = norm2(r) / bnorm;

            for (int i = 0; i <= restart; i++) java.util.Arrays.fill(H[i], 0.0);

            if (relres <= tol) {
                flag = 0;
                break;
            }
            // see _kb/03-api-layer.md for rationale
            if (!(relres < 1e2 * initres)) {
                flag = 3;
                break;
            }
            // see _kb/03-api-layer.md for rationale
            if (relres >= prevrelres * (1.0 - 1e-12)) {
                flag = 3;
                break;
            }
            prevrelres = relres;
        }

        for (int i = 0; i < n; i++) {
            if (Double.isNaN(x[i]) || Double.isInfinite(x[i])) {
                return new GmresResult(asColumn(unpermute(x, perm)), 3, Double.POSITIVE_INFINITY, iter);
            }
        }
        if (relres <= tol) flag = 0;
        return new GmresResult(asColumn(unpermute(x, perm)), flag, relres, iter);
    }

    /** Returns out with out[i] = v[perm[i]]. */
    private static double[] permute(double[] v, int[] perm) {
        double[] out = new double[v.length];
        for (int i = 0; i < v.length; i++) out[i] = v[perm[i]];
        return out;
    }

    /** Returns out with out[perm[i]] = v[i], the inverse of {@link #permute}. */
    private static double[] unpermute(double[] v, int[] perm) {
        double[] out = new double[v.length];
        for (int i = 0; i < v.length; i++) out[perm[i]] = v[i];
        return out;
    }

    private static Matrix asColumn(double[] v) {
        Matrix out = new Matrix(v.length, 1);
        for (int i = 0; i < v.length; i++) out.set(i, 0, v[i]);
        return out;
    }

    private static double norm2(double[] v) {
        double s = 0.0;
        for (double vi : v) s += vi * vi;
        return Math.sqrt(s);
    }

    private static double dot(double[] a, double[] b) {
        double s = 0.0;
        for (int i = 0; i < a.length; i++) s += a[i] * b[i];
        return s;
    }

    /**
     * Compressed sparse row view of a matrix, with the column indices of each row in
     * increasing order and the position of the diagonal entry recorded. The incomplete
     * factorization is a row-oriented elimination, so it needs this rather than the CSC
     * layout EJML uses.
     */
    private static final class Csr {
        final int n;
        final int[] rowPtr;
        final int[] colIdx;
        final double[] val;
        final int[] diagPtr;

        private Csr(int n, int[] rowPtr, int[] colIdx, double[] val) {
            this.n = n;
            this.rowPtr = rowPtr;
            this.colIdx = colIdx;
            this.val = val;
            this.diagPtr = new int[n];
            for (int i = 0; i < n; i++) {
                diagPtr[i] = -1;
                for (int p = rowPtr[i]; p < rowPtr[i + 1]; p++) {
                    if (colIdx[p] == i) {
                        diagPtr[i] = p;
                        break;
                    }
                }
            }
        }

        static Csr of(DMatrixSparseCSC csc) {
            int n = csc.getNumRows();
            int[] colStart = csc.col_idx;
            int[] rows = csc.nz_rows;
            double[] values = csc.nz_values;
            int nnz = csc.nz_length;

            int[] rowPtr = new int[n + 1];
            for (int p = 0; p < nnz; p++) rowPtr[rows[p] + 1]++;
            for (int i = 0; i < n; i++) rowPtr[i + 1] += rowPtr[i];

            int[] colIdx = new int[nnz];
            double[] val = new double[nnz];
            int[] next = new int[n];
            System.arraycopy(rowPtr, 0, next, 0, n);
            // Walking the columns in increasing order leaves each row's column indices
            // sorted, which the elimination below relies on.
            for (int c = 0; c < n; c++) {
                for (int p = colStart[c]; p < colStart[c + 1]; p++) {
                    int r = rows[p];
                    colIdx[next[r]] = c;
                    val[next[r]] = values[p];
                    next[r]++;
                }
            }
            return new Csr(n, rowPtr, colIdx, val);
        }

        /** Computes out = this * v. */
        void mult(double[] v, double[] out) {
            for (int i = 0; i < n; i++) {
                double s = 0.0;
                for (int p = rowPtr[i]; p < rowPtr[i + 1]; p++) s += val[p] * v[colIdx[p]];
                out[i] = s;
            }
        }

        /**
         * Scales every row to unit max norm and returns the scaling applied, which the
         * caller must divide into each right-hand side. Row scaling leaves the solution
         * unchanged.
         */
        double[] equilibrate() {
            double[] scale = new double[n];
            for (int i = 0; i < n; i++) {
                double m = 0.0;
                for (int p = rowPtr[i]; p < rowPtr[i + 1]; p++) m = Math.max(m, Math.abs(val[p]));
                scale[i] = (m == 0.0) ? 1.0 : m;
                if (m == 0.0 || m == 1.0) continue;
                for (int p = rowPtr[i]; p < rowPtr[i + 1]; p++) val[p] /= m;
            }
            return scale;
        }

        /** Returns the matrix reordered so that entry (i,j) becomes (iperm[i], iperm[j]). */
        Csr permuteSymmetric(int[] perm, int[] iperm) {
            int nnz = rowPtr[n];
            int[] newRowPtr = new int[n + 1];
            for (int i = 0; i < n; i++) {
                int oi = perm[i];
                newRowPtr[i + 1] = newRowPtr[i] + (rowPtr[oi + 1] - rowPtr[oi]);
            }
            int[] newColIdx = new int[nnz];
            double[] newVal = new double[nnz];
            for (int i = 0; i < n; i++) {
                int oi = perm[i];
                int len = rowPtr[oi + 1] - rowPtr[oi];
                int base = newRowPtr[i];
                for (int t = 0; t < len; t++) {
                    newColIdx[base + t] = iperm[colIdx[rowPtr[oi] + t]];
                    newVal[base + t] = val[rowPtr[oi] + t];
                }
                sortRow(newColIdx, newVal, base, base + len);
            }
            return new Csr(n, newRowPtr, newColIdx, newVal);
        }

        /** Insertion sort of one row by column index; rows are short, so this is enough. */
        private static void sortRow(int[] cols, double[] vals, int from, int to) {
            for (int i = from + 1; i < to; i++) {
                int c = cols[i];
                double v = vals[i];
                int j = i - 1;
                while (j >= from && cols[j] > c) {
                    cols[j + 1] = cols[j];
                    vals[j + 1] = vals[j];
                    j--;
                }
                cols[j + 1] = c;
                vals[j + 1] = v;
            }
        }
    }

    /**
     * Reverse Cuthill-McKee ordering of the symmetrized sparsity pattern. The level
     * structure is grown from the lowest-degree unvisited node of each component, with
     * the neighbours of each node appended in order of increasing degree, and the result
     * is reversed.
     */
    private static final class Rcm {
        static int[] order(Csr a) {
            int n = a.n;
            // Symmetric adjacency of the pattern, diagonal excluded.
            int[] deg = new int[n];
            int[] adjPtr = new int[n + 1];
            for (int i = 0; i < n; i++) {
                for (int p = a.rowPtr[i]; p < a.rowPtr[i + 1]; p++) {
                    int j = a.colIdx[p];
                    if (j == i) continue;
                    deg[i]++;
                    deg[j]++;
                }
            }
            for (int i = 0; i < n; i++) adjPtr[i + 1] = adjPtr[i] + deg[i];
            int[] adj = new int[adjPtr[n]];
            int[] fill = new int[n];
            System.arraycopy(adjPtr, 0, fill, 0, n);
            for (int i = 0; i < n; i++) {
                for (int p = a.rowPtr[i]; p < a.rowPtr[i + 1]; p++) {
                    int j = a.colIdx[p];
                    if (j == i) continue;
                    adj[fill[i]++] = j;
                    adj[fill[j]++] = i;
                }
            }
            // Duplicate neighbours, produced when both (i,j) and (j,i) are nonzero, are
            // harmless here: they only bias the degree used for tie-breaking.

            boolean[] seen = new boolean[n];
            int[] result = new int[n];
            int count = 0;
            int[] queue = new int[n];

            while (count < n) {
                int start = -1;
                for (int i = 0; i < n; i++) {
                    if (!seen[i] && (start < 0 || deg[i] < deg[start])) start = i;
                }
                int head = count;
                int tail = count;
                queue[tail++] = start;
                seen[start] = true;
                while (head < tail) {
                    int v = queue[head++];
                    result[count++] = v;
                    int from = tail;
                    for (int p = adjPtr[v]; p < adjPtr[v + 1]; p++) {
                        int u = adj[p];
                        if (!seen[u]) {
                            seen[u] = true;
                            queue[tail++] = u;
                        }
                    }
                    // Order the newly discovered neighbours by increasing degree.
                    for (int i = from + 1; i < tail; i++) {
                        int c = queue[i];
                        int j = i - 1;
                        while (j >= from && deg[queue[j]] > deg[c]) {
                            queue[j + 1] = queue[j];
                            j--;
                        }
                        queue[j + 1] = c;
                    }
                }
            }

            // Reverse: Cuthill-McKee reversed is what reduces the profile.
            int[] perm = new int[n];
            for (int i = 0; i < n; i++) perm[i] = result[n - 1 - i];
            return perm;
        }
    }

    /**
     * Right preconditioner: either the ILUT factors of A, applied by a forward and a
     * backward triangular sweep, or a Jacobi diagonal when the incomplete factorization
     * breaks down.
     */
    private static final class Precond {
        private final Ilut lu;
        private final double[] dinv;

        private Precond(Ilut lu, double[] dinv) {
            this.lu = lu;
            this.dinv = dinv;
        }

        static Precond of(Csr a) {
            Ilut lu = Ilut.factorize(a, ILUT_DROP_TOL, ILUT_FILL_FACTOR);
            if (lu != null) return new Precond(lu, null);
            // Jacobi fallback. A zero diagonal entry would make the preconditioner
            // singular, so those rows are left unscaled rather than inverted.
            int n = a.n;
            double[] dinv = new double[n];
            for (int i = 0; i < n; i++) {
                double d = a.diagPtr[i] < 0 ? 0.0 : a.val[a.diagPtr[i]];
                dinv[i] = d == 0.0 ? 1.0 : 1.0 / d;
            }
            return new Precond(null, dinv);
        }

        /** Computes out = M^{-1} * v. */
        void apply(double[] v, double[] out) {
            if (dinv != null) {
                for (int i = 0; i < dinv.length; i++) out[i] = dinv[i] * v[i];
                return;
            }
            lu.solve(v, out);
        }
    }

    /**
     * Threshold incomplete LU factorization, ILUT(p, tau) of Saad (1994). Row i of the
     * IKJ elimination is expanded into a dense workspace, entries below tau times the
     * average magnitude of the original row are dropped as they are produced, and the p
     * largest survivors are retained in each of the L and U parts of the row. The
     * diagonal is never allowed to vanish: a pivot that drops below the row threshold is
     * replaced, which is Saad's remedy and keeps the factorization defined without
     * pivoting.
     */
    private static final class Ilut {
        private final int n;
        private final int[] lPtr;
        private final int[] lCol;
        private final double[] lVal;
        private final int[] uPtr;
        private final int[] uCol;
        private final double[] uVal;
        private final double[] dInv;

        private Ilut(int n, int[] lPtr, int[] lCol, double[] lVal,
                     int[] uPtr, int[] uCol, double[] uVal, double[] dInv) {
            this.n = n;
            this.lPtr = lPtr;
            this.lCol = lCol;
            this.lVal = lVal;
            this.uPtr = uPtr;
            this.uCol = uCol;
            this.uVal = uVal;
            this.dInv = dInv;
        }

        /** Computes out = U^{-1} L^{-1} v, with L unit lower triangular. */
        void solve(double[] v, double[] out) {
            for (int i = 0; i < n; i++) {
                double s = v[i];
                for (int p = lPtr[i]; p < lPtr[i + 1]; p++) s -= lVal[p] * out[lCol[p]];
                out[i] = s;
            }
            for (int i = n - 1; i >= 0; i--) {
                double s = out[i];
                for (int p = uPtr[i]; p < uPtr[i + 1]; p++) s -= uVal[p] * out[uCol[p]];
                out[i] = s * dInv[i];
            }
        }

        static Ilut factorize(Csr a, double dropTol, double fillFactor) {
            int n = a.n;
            int nnz = a.rowPtr[n];
            int lfil = Math.max(1, (int) Math.ceil(fillFactor * nnz / Math.max(1, n)));

            // Initial factor capacity; the arrays grow on demand, so this is only a hint.
            int cap = (int) Math.max(16L, Math.min((long) nnz * 2L, 1L << 24));
            int[] lPtr = new int[n + 1];
            int[] uPtr = new int[n + 1];
            int[] lCol = new int[cap];
            double[] lVal = new double[cap];
            int[] uCol = new int[cap];
            double[] uVal = new double[cap];
            double[] dInv = new double[n];
            int lLen = 0;
            int uLen = 0;

            // Dense accumulator for the current row, with an index list so that clearing
            // it costs only the number of entries touched.
            double[] w = new double[n];
            int[] wPos = new int[n];
            java.util.Arrays.fill(wPos, -1);
            int[] wIdx = new int[n];
            int wCount = 0;

            // see _kb/03-api-layer.md for rationale
            java.util.PriorityQueue<Integer> pending = new java.util.PriorityQueue<Integer>();

            int[] rowsL = new int[n];
            int[] rowsU = new int[n];

            for (int i = 0; i < n; i++) {
                double tnorm = 0.0;
                int rowLen = a.rowPtr[i + 1] - a.rowPtr[i];
                for (int p = a.rowPtr[i]; p < a.rowPtr[i + 1]; p++) {
                    int j = a.colIdx[p];
                    w[j] = a.val[p];
                    wPos[j] = wCount;
                    wIdx[wCount++] = j;
                    tnorm += Math.abs(a.val[p]);
                    if (j < i) pending.add(j);
                }
                if (rowLen == 0 || tnorm == 0.0) {
                    // An empty row leaves the factorization singular; report breakdown
                    // rather than invent a pivot for a row that carries no information.
                    return null;
                }
                double tau = dropTol * tnorm / rowLen;

                while (!pending.isEmpty()) {
                    int k = pending.poll();
                    double mult = w[k] * dInv[k];
                    if (Math.abs(mult) <= tau) {
                        w[k] = 0.0;
                        continue;
                    }
                    w[k] = mult;
                    for (int p = uPtr[k]; p < uPtr[k + 1]; p++) {
                        int j = uCol[p];
                        double upd = mult * uVal[p];
                        if (wPos[j] >= 0) {
                            w[j] -= upd;
                        } else {
                            if (Math.abs(upd) <= tau) continue;
                            w[j] = -upd;
                            wPos[j] = wCount;
                            wIdx[wCount++] = j;
                            if (j < i) pending.add(j);
                        }
                    }
                }

                // Split the row, keeping the lfil largest entries on each side.
                int nl = 0;
                int nu = 0;
                double diag = w[i];
                for (int t = 0; t < wCount; t++) {
                    int j = wIdx[t];
                    if (j == i) continue;
                    if (Math.abs(w[j]) <= tau) continue;
                    if (j < i) rowsL[nl++] = j;
                    else rowsU[nu++] = j;
                }
                nl = keepLargest(rowsL, w, nl, lfil);
                nu = keepLargest(rowsU, w, nu, lfil);
                java.util.Arrays.sort(rowsL, 0, nl);
                java.util.Arrays.sort(rowsU, 0, nu);

                if (lLen + nl > lCol.length) {
                    int grow = Math.max(lCol.length * 2, lLen + nl);
                    lCol = java.util.Arrays.copyOf(lCol, grow);
                    lVal = java.util.Arrays.copyOf(lVal, grow);
                }
                if (uLen + nu > uCol.length) {
                    int grow = Math.max(uCol.length * 2, uLen + nu);
                    uCol = java.util.Arrays.copyOf(uCol, grow);
                    uVal = java.util.Arrays.copyOf(uVal, grow);
                }
                for (int t = 0; t < nl; t++) {
                    lCol[lLen] = rowsL[t];
                    lVal[lLen++] = w[rowsL[t]];
                }
                for (int t = 0; t < nu; t++) {
                    uCol[uLen] = rowsU[t];
                    uVal[uLen++] = w[rowsU[t]];
                }
                lPtr[i + 1] = lLen;
                uPtr[i + 1] = uLen;

                // see _kb/03-api-layer.md for rationale
                if (!(Math.abs(diag) > tau) || Double.isNaN(diag)) {
                    double substitute = tau > 0 ? tau : 1e-8;
                    diag = diag < 0 ? -substitute : substitute;
                }
                dInv[i] = 1.0 / diag;
                if (Double.isInfinite(dInv[i]) || Double.isNaN(dInv[i])) return null;

                for (int t = 0; t < wCount; t++) {
                    w[wIdx[t]] = 0.0;
                    wPos[wIdx[t]] = -1;
                }
                wCount = 0;
                pending.clear();
            }

            return new Ilut(n, lPtr, lCol, lVal, uPtr, uCol, uVal, dInv);
        }

        /**
         * Reduces cols[0..len) to the keep entries of largest magnitude in w, by partial
         * selection. Order within the retained set is irrelevant; the caller sorts.
         */
        private static int keepLargest(int[] cols, double[] w, int len, int keep) {
            if (len <= keep) return len;
            for (int i = 0; i < keep; i++) {
                int best = i;
                for (int j = i + 1; j < len; j++) {
                    if (Math.abs(w[cols[j]]) > Math.abs(w[cols[best]])) best = j;
                }
                int tmp = cols[i];
                cols[i] = cols[best];
                cols[best] = tmp;
            }
            return keep;
        }
    }
}
