/**
 * @file MMAP(3,K) closed-form marking fit: a marked MAP of third order
 *
 * The MMAP(2,K) argument does not depend on the order. Two facts carry over
 * (verified symbolically, see io/sage/proofs/mmap3k_marking_inverse.py):
 *
 * 1. every per-class characteristic in which the class matrix appears exactly
 *    once is LINEAR in the marking fractions, so z = nnz(D1) fractions per
 *    class are determined by z characteristics through a square system;
 * 2. that system is BLOCK DIAGONAL in the classes, so one z x z block is built
 *    once and reused for every class: the cost does not grow with K.
 *
 * What changes with the order is WHICH characteristics are needed. At order two
 * (p_c, F_c, B_c) suffice; at order three the independent set of lowest total
 * order is (a,b) = (1,0), (1,1), (2,0), (3,0), that is p_c, F_c, B_c and the
 * second-order backward moment, with a the backward and b the forward order of
 * pie A^a D1c A^b 1. Alternating forward and backward orders does NOT stay
 * independent at higher orders.
 *
 * The block is assembled exactly by evaluating the linear map on unit markings,
 * so there are no finite differences and no computer algebra at run time. The
 * underlying MAP(3) is an input: unlike order two there is no canonical inverse
 * here that turns moments and autocorrelation into an order-3 MAP.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap3k_fit {
    private Mmap3k_fit() {}

    private static final double FEASTOL = 1e-8;

    /** Result of a fit: the marked MAP and whether the marking is feasible. */
    public static final class Result {
        public final MatrixCell mmap;
        public final boolean exact;

        public Result(MatrixCell mmap, boolean exact) {
            this.mmap = mmap;
            this.exact = exact;
        }
    }

    /** (backward, forward) orders of an independent characteristic set. */
    public static int[][] markingOrders(int n, int z) {
        List<int[]> out = new ArrayList<int[]>();
        out.add(new int[]{1, 0});
        out.add(new int[]{1, 1});
        int a = 2;
        while (out.size() < n + 1) {
            out.add(new int[]{a, 0});
            a++;
        }
        int[][] arr = new int[z][2];
        for (int i = 0; i < z; i++) {
            arr[i] = out.get(i);
        }
        return arr;
    }

    /**
     * Marks a given MAP so that the per-class characteristics are matched.
     *
     * @param D0 hidden transition matrix of the underlying MAP
     * @param D1 visible transition matrix of the underlying MAP
     * @param P class probabilities, summing to one
     * @param F first-order forward moments
     * @param B first-order backward moments
     * @param B2 second-order backward moments, required from order three
     */
    public static Result mmap3k_fit(Matrix D0, Matrix D1, double[] P, double[] F,
                                    double[] B, double[] B2) {
        int n = D0.getNumRows();
        int K = P.length;

        List<int[]> nz = new ArrayList<int[]>();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (D1.get(i, j) != 0.0) {
                    nz.add(new int[]{i, j});
                }
            }
        }
        int z = nz.size();
        int[][] orders = markingOrders(n, z);
        if (z > 3 && B2 == null) {
            throw new IllegalArgumentException(
                "mmap3k_fit: order " + n + " needs the second-order backward moments");
        }

        Matrix A = D0.scale(-1.0).inv();
        Matrix P_emb = A.mult(D1);
        Matrix T = P_emb.transpose().sub(1.0, Matrix.eye(n));
        for (int j = 0; j < n; j++) {
            T.set(n - 1, j, 1.0);
        }
        Matrix rhs = new Matrix(n, 1);
        rhs.set(n - 1, 0, 1.0);
        Matrix pieCol = new Matrix(n, 1);
        if (!Matrix.solve(T, rhs, pieCol)) {
            throw new RuntimeException("mmap3k_fit: the embedded chain is singular");
        }
        Matrix pie = pieCol.transpose();

        double[][] M = new double[z][z];
        for (int jj = 0; jj < z; jj++) {
            Matrix Dc = Matrix.zeros(n, n);
            int[] pos = nz.get(jj);
            Dc.set(pos[0], pos[1], D1.get(pos[0], pos[1]));
            for (int ii = 0; ii < z; ii++) {
                M[ii][jj] = functional(pie, A, Dc, orders[ii][0], orders[ii][1], n);
            }
        }

        double[][] q = new double[z][K];
        for (int c = 0; c < K; c++) {
            double[] y = new double[z];
            for (int ii = 0; ii < z; ii++) {
                int a = orders[ii][0];
                int b = orders[ii][1];
                if (a == 1 && b == 0) {
                    y[ii] = P[c];
                } else if (a == 1 && b == 1) {
                    y[ii] = P[c] * F[c];
                } else if (a == 2 && b == 0) {
                    y[ii] = P[c] * B[c];
                } else if (a == 3 && b == 0) {
                    y[ii] = P[c] * B2[c];
                } else {
                    throw new IllegalArgumentException(
                        "mmap3k_fit: no target for the characteristic (a=" + a + ", b=" + b + ")");
                }
            }
            double[] sol = solve(M, y);
            for (int ii = 0; ii < z; ii++) {
                q[ii][c] = sol[ii];
            }
        }

        double viol = 0.0;
        for (int j = 0; j < z; j++) {
            double sum = 0.0;
            for (int c = 0; c < K; c++) {
                viol = Math.max(viol, -q[j][c]);
                viol = Math.max(viol, q[j][c] - 1.0);
                sum += q[j][c];
            }
            viol = Math.max(viol, Math.abs(sum - 1.0));
        }

        MatrixCell out = new MatrixCell(2 + K);
        out.set(0, D0.copy());
        out.set(1, D1.copy());
        for (int c = 0; c < K; c++) {
            Matrix Dc = Matrix.zeros(n, n);
            for (int jj = 0; jj < z; jj++) {
                int[] pos = nz.get(jj);
                double v = Math.min(Math.max(q[jj][c], 0.0), 1.0);
                Dc.set(pos[0], pos[1], D1.get(pos[0], pos[1]) * v);
            }
            out.set(2 + c, Dc);
        }
        return new Result(out, viol <= FEASTOL);
    }

    private static double functional(Matrix pie, Matrix A, Matrix Dc, int a, int b, int n) {
        Matrix left = pie.copy();
        for (int k = 0; k < a; k++) {
            left = left.mult(A);
        }
        Matrix mid = left.mult(Dc);
        for (int k = 0; k < b; k++) {
            mid = mid.mult(A);
        }
        double s = 0.0;
        for (int j = 0; j < n; j++) {
            s += mid.get(0, j);
        }
        return s;
    }

    /** Gaussian elimination with partial pivoting on a small dense system. */
    private static double[] solve(double[][] Ain, double[] bin) {
        int m = bin.length;
        double[][] M = new double[m][m + 1];
        for (int i = 0; i < m; i++) {
            System.arraycopy(Ain[i], 0, M[i], 0, m);
            M[i][m] = bin[i];
        }
        for (int col = 0; col < m; col++) {
            int piv = col;
            for (int i = col + 1; i < m; i++) {
                if (Math.abs(M[i][col]) > Math.abs(M[piv][col])) {
                    piv = i;
                }
            }
            if (Math.abs(M[piv][col]) < 1e-14) {
                throw new RuntimeException("mmap3k_fit: the marking system is singular");
            }
            double[] tmp = M[col];
            M[col] = M[piv];
            M[piv] = tmp;
            for (int i = col + 1; i < m; i++) {
                double f = M[i][col] / M[col][col];
                for (int j = col; j <= m; j++) {
                    M[i][j] -= f * M[col][j];
                }
            }
        }
        double[] x = new double[m];
        for (int i = m - 1; i >= 0; i--) {
            double s = M[i][m];
            for (int j = i + 1; j < m; j++) {
                s -= M[i][j] * x[j];
            }
            x[i] = s / M[i][i];
        }
        return x;
    }
}
