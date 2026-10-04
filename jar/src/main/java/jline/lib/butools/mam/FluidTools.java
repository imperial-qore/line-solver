/*
 * Helpers for the BUTools-family level-dependent / priority fluid queue ports.
 */
package jline.lib.butools.mam;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import jline.util.matrix.Matrix;

/**
 * Small matrix helpers shared by the fluid-queue ports (diagonal extraction,
 * index-based submatrices, block assembly). Kept package-local to mirror the
 * self-contained BUTools thirdparty style.
 */
public final class FluidTools {
    private FluidTools() {}

    /** Diagonal of M as a double array. */
    public static double[] diagArray(Matrix M) {
        int n = Math.min(M.getNumRows(), M.getNumCols());
        double[] d = new double[n];
        for (int i = 0; i < n; i++) d[i] = M.get(i, i);
        return d;
    }

    /** Diagonal matrix from a double array. */
    public static Matrix diagFrom(double[] v) {
        Matrix D = Matrix.zeros(v.length, v.length);
        for (int i = 0; i < v.length; i++) D.set(i, i, v[i]);
        return D;
    }

    /** Diagonal matrix from a row/column vector. */
    public static Matrix diagFrom(Matrix v) {
        int n = v.getNumRows() * v.getNumCols();
        double[] a = new double[n];
        int k = 0;
        for (int i = 0; i < v.getNumRows(); i++)
            for (int j = 0; j < v.getNumCols(); j++)
                a[k++] = v.get(i, j);
        return diagFrom(a);
    }

    /** Submatrix M(ri, ci) with index arrays (0-based). */
    public static Matrix sub(Matrix M, int[] ri, int[] ci) {
        Matrix R = Matrix.zeros(ri.length, ci.length);
        for (int i = 0; i < ri.length; i++)
            for (int j = 0; j < ci.length; j++)
                R.set(i, j, M.get(ri[i], ci[j]));
        return R;
    }

    /** Assigns src into dst at rows ri, cols ci. */
    public static void setSub(Matrix dst, int[] ri, int[] ci, Matrix src) {
        for (int i = 0; i < ri.length; i++)
            for (int j = 0; j < ci.length; j++)
                dst.set(ri[i], ci[j], src.get(i, j));
    }

    /** Assigns src (rows x |cols|) into the given columns of dst (all rows). */
    public static void setCols(Matrix dst, int[] cols, Matrix src) {
        for (int i = 0; i < dst.getNumRows(); i++)
            for (int j = 0; j < cols.length; j++)
                dst.set(i, cols[j], src.get(i, j));
    }

    /** [A B] side by side (handles zero-dimension blocks). */
    public static Matrix hstack(Matrix A, Matrix B) {
        int rows = Math.max(A.getNumRows(), B.getNumRows());
        Matrix R = Matrix.zeros(rows, A.getNumCols() + B.getNumCols());
        for (int i = 0; i < A.getNumRows(); i++)
            for (int j = 0; j < A.getNumCols(); j++) R.set(i, j, A.get(i, j));
        for (int i = 0; i < B.getNumRows(); i++)
            for (int j = 0; j < B.getNumCols(); j++) R.set(i, A.getNumCols() + j, B.get(i, j));
        return R;
    }

    /** [A; B] stacked (handles zero-dimension blocks). */
    public static Matrix vstack(Matrix A, Matrix B) {
        int cols = Math.max(A.getNumCols(), B.getNumCols());
        Matrix R = Matrix.zeros(A.getNumRows() + B.getNumRows(), cols);
        for (int i = 0; i < A.getNumRows(); i++)
            for (int j = 0; j < A.getNumCols(); j++) R.set(i, j, A.get(i, j));
        for (int i = 0; i < B.getNumRows(); i++)
            for (int j = 0; j < B.getNumCols(); j++) R.set(A.getNumRows() + i, j, B.get(i, j));
        return R;
    }

    /** 2x2 block matrix [[A B];[C D]]. */
    public static Matrix block(Matrix A, Matrix B, Matrix C, Matrix D) {
        return vstack(hstack(A, B), hstack(C, D));
    }

    /** Block-diagonal [[A 0];[0 B]]. */
    public static Matrix blkdiag(Matrix A, Matrix B) {
        return block(A, Matrix.zeros(A.getNumRows(), B.getNumCols()),
                Matrix.zeros(B.getNumRows(), A.getNumCols()), B);
    }

    /** Concatenate two index arrays. */
    public static int[] concat(int[] a, int[] b) {
        int[] r = new int[a.length + b.length];
        System.arraycopy(a, 0, r, 0, a.length);
        System.arraycopy(b, 0, r, a.length, b.length);
        return r;
    }

    /** 0..n-1. */
    public static int[] range(int n) {
        int[] r = new int[n];
        for (int i = 0; i < n; i++) r[i] = i;
        return r;
    }

    /** Column-major reshape of M into a single column vector. */
    public static Matrix vec(Matrix M) {
        return M.columnMajorOrder();
    }

    /** Row vector of row-sums (each row summed across columns). */
    public static Matrix rowSums(Matrix M) {
        Matrix out = Matrix.zeros(M.getNumRows(), 1);
        for (int i = 0; i < M.getNumRows(); i++) {
            double s = 0;
            for (int j = 0; j < M.getNumCols(); j++) s += M.get(i, j);
            out.set(i, 0, s);
        }
        return out;
    }

    /** Left null vector of KA via the CRPSolve algorithm (column, no generator check). */
    public static Matrix nullvec(Matrix KA) {
        int n = KA.getNumRows();
        Matrix M = KA.copy();
        for (int i = 0; i < n; i++) M.set(i, 0, 1.0);
        Matrix e = Matrix.zeros(n, 1);
        e.set(0, 0, 1.0);
        return M.transpose().leftMatrixDivide(e); // (n,1)
    }

    /** Minimum modulus of the eigenvalues of M. */
    public static double minAbsEig(Matrix M) {
        double m = Double.POSITIVE_INFINITY;
        for (Complex c : M.eig()) m = Math.min(m, c.abs());
        return m;
    }

    /** int_0^L expm(KA u) du, robust when KA has a zero eigenvalue (deflation). */
    public static Matrix integExp(Matrix KA, double L) {
        int n = KA.getNumRows();
        Matrix lrow = nullvec(KA).transpose();       // (1,n)
        Matrix rcol = nullvec(KA.transpose());       // (n,1)
        double lr = lrow.mult(rcol).get(0, 0);
        lrow = lrow.scale(1.0 / lr);
        Matrix rl = rcol.mult(lrow);                 // (n,n)
        Matrix Kd = KA.sub(rl);
        return Kd.scale(-1.0).pinv().mult(Matrix.eye(n).sub(Kd.scale(L).expm()))
                .add(rl.scale(L + Math.exp(-L) - 1.0));
    }

    /** {int_0^L expm(KA u)du, int_0^L expm(KB u)du}, applying deflation to the near-singular one. */
    public static Matrix[] integExp2(Matrix KA, Matrix KB, double L) {
        double eKA = KA.getNumRows() > 0 ? minAbsEig(KA) : Double.POSITIVE_INFINITY;
        double eKB = KB.getNumRows() > 0 ? minAbsEig(KB) : Double.POSITIVE_INFINITY;
        Matrix KAi, KBi;
        if (eKA > eKB) {
            KAi = KA.scale(-1.0).pinv().mult(Matrix.eye(KA.getNumRows()).sub(KA.scale(L).expm()));
            KBi = integExp(KB, L);
        } else {
            KAi = integExp(KA, L);
            KBi = KB.scale(-1.0).pinv().mult(Matrix.eye(KB.getNumRows()).sub(KB.scale(L).expm()));
        }
        return new Matrix[]{KAi, KBi};
    }

    /** {J0, J1} with J0=int_0^L expm(M u)du, J1=int_0^L u expm(M u)du via nilpotent augmentation. */
    public static Matrix[] expIntMoments(Matrix M, double L) {
        int n = M.getNumRows();
        Matrix A = Matrix.zeros(3 * n, n * 3);
        setSub(A, range(n), range(n), M);
        setSub(A, range(n), shiftRange(n, n), Matrix.eye(n));
        setSub(A, shiftRange(n, n), shiftRange(2 * n, n), Matrix.eye(n));
        Matrix W = A.scale(L).expm();
        Matrix J0 = sub(W, range(n), shiftRange(n, n));
        Matrix W13 = sub(W, range(n), shiftRange(2 * n, n));
        Matrix J1 = J0.scale(L).sub(W13);
        return new Matrix[]{J0, J1};
    }

    /** start..start+n-1. */
    public static int[] shiftRange(int start, int n) {
        int[] r = new int[n];
        for (int i = 0; i < n; i++) r[i] = start + i;
        return r;
    }

    // ---------------- ordered real Schur decomposition ----------------

    /**
     * Real Schur decomposition of A reordered so that the eigenvalues appear
     * grouped as [zero real-part, negative real-part, positive real-part], the
     * ordering used by MATLAB's ordschur in the multi-regime fluid solver.
     * Returns {Z, T, int[]{zeroeig, negeig, poseig}} with Z orthogonal,
     * T quasi-upper-triangular and Z'*A*Z = T.
     */
    public static Object[] orderedSchur(Matrix A, double tol) {
        Map<String, Matrix> sd = A.schur();
        Matrix Z = sd.get("U").copy();
        Matrix T = sd.get("T").copy();
        List<int[]> blocks = blockList(T);
        int nb = blocks.size();
        int[] size = new int[nb], rank = new int[nb];
        for (int i = 0; i < nb; i++) {
            int s = blocks.get(i)[0], sz = blocks.get(i)[1];
            double re = sz == 1 ? T.get(s, s) : 0.5 * (T.get(s, s) + T.get(s + 1, s + 1));
            size[i] = sz;
            rank[i] = rankOf(re, tol);
        }
        // bubble sort blocks by descending rank (zero>neg>pos), fixed ranks -> terminates
        boolean swapped = true;
        while (swapped) {
            swapped = false;
            int pos = 0;
            for (int bi = 0; bi < nb - 1; bi++) {
                if (rank[bi + 1] > rank[bi]) {
                    swapAdjacentBlocks(T, Z, pos, size[bi], size[bi + 1]);
                    int ts = size[bi]; size[bi] = size[bi + 1]; size[bi + 1] = ts;
                    int tr = rank[bi]; rank[bi] = rank[bi + 1]; rank[bi + 1] = tr;
                    swapped = true;
                }
                pos += size[bi];
            }
        }
        int zeroeig = 0, negeig = 0, poseig = 0;
        for (int i = 0; i < nb; i++) {
            if (rank[i] == 2) zeroeig += size[i];
            else if (rank[i] == 1) negeig += size[i];
            else poseig += size[i];
        }
        return new Object[]{Z, T, new int[]{zeroeig, negeig, poseig}};
    }

    /**
     * Real Schur decomposition of A reordered so that eigenvalues with
     * positive real part appear in the leading block, matching MATLAB's
     * ordschur(U,T,'rhp'). Returns {Z, T} with Z orthogonal,
     * T quasi-upper-triangular and Z'*A*Z = T.
     */
    public static Matrix[] orderedSchurRhp(Matrix A) {
        Map<String, Matrix> sd = A.schur();
        Matrix Z = sd.get("U").copy();
        Matrix T = sd.get("T").copy();
        List<int[]> blocks = blockList(T);
        int nb = blocks.size();
        int[] size = new int[nb], rank = new int[nb];
        for (int i = 0; i < nb; i++) {
            int s = blocks.get(i)[0], sz = blocks.get(i)[1];
            double re = sz == 1 ? T.get(s, s) : 0.5 * (T.get(s, s) + T.get(s + 1, s + 1));
            size[i] = sz;
            rank[i] = re > 0 ? 1 : 0;
        }
        boolean swapped = true;
        while (swapped) {
            swapped = false;
            int pos = 0;
            for (int bi = 0; bi < nb - 1; bi++) {
                if (rank[bi + 1] > rank[bi]) {
                    swapAdjacentBlocks(T, Z, pos, size[bi], size[bi + 1]);
                    int ts = size[bi]; size[bi] = size[bi + 1]; size[bi + 1] = ts;
                    int tr = rank[bi]; rank[bi] = rank[bi + 1]; rank[bi + 1] = tr;
                    swapped = true;
                }
                pos += size[bi];
            }
        }
        return new Matrix[]{Z, T};
    }

    /** Diagonal 1x1/2x2 block partition of a real Schur form: each entry = {start, size}. */
    private static List<int[]> blockList(Matrix T) {
        int n = T.getNumRows();
        List<int[]> bl = new ArrayList<int[]>();
        int i = 0;
        while (i < n) {
            if (i < n - 1 && Math.abs(T.get(i + 1, i)) > 1e-13 * (1 + Math.abs(T.get(i, i)) + Math.abs(T.get(i + 1, i + 1)))) {
                bl.add(new int[]{i, 2});
                i += 2;
            } else {
                bl.add(new int[]{i, 1});
                i += 1;
            }
        }
        return bl;
    }

    /** Cluster rank of a real part: zero=2 (top), negative=1, positive=0. */
    private static int rankOf(double re, double tol) {
        if (Math.abs(re) < tol) return 2;
        return re < 0 ? 1 : 0;
    }

    /**
     * Swap the adjacent diagonal blocks of a real Schur form: the p-sized block
     * at [j0,j0+p) and the q-sized block at [j0+p,j0+p+q), updating T=Q'TQ and
     * Z=ZQ via an orthogonal similarity (Bai-Demmel direct swap).
     */
    private static void swapAdjacentBlocks(Matrix T, Matrix Z, int j0, int p, int q) {
        int n = T.getNumRows();
        int w = p + q;
        int[] rp = shiftRange(j0, p), rq = shiftRange(j0 + p, q), win = shiftRange(j0, w);
        Matrix A11 = sub(T, rp, rp);
        Matrix B22 = sub(T, rq, rq);
        Matrix C = sub(T, rp, rq);
        // solve A11 X - X B22 = -C  ->  X spans the B22-invariant subspace of [[A11,C],[0,B22]]
        Matrix X = triSylvester(A11, B22, C.scale(-1.0));
        Matrix stack = vstack(X, Matrix.eye(q)); // (p+q) x q
        Matrix Q = fullQ(stack);                  // full orthogonal (p+q)x(p+q)
        int[] allRows = range(n);
        // T <- Q' T Q applied to the window rows then columns
        setSub(T, win, allRows, Q.transpose().mult(sub(T, win, allRows)));
        setSub(T, allRows, win, sub(T, allRows, win).mult(Q));
        // Z <- Z Q on the window columns
        setSub(Z, allRows, win, sub(Z, allRows, win).mult(Q));
    }

    /** Full orthogonal factor Q (m x m) of a tall matrix (m >= n) via QR. */
    private static Matrix fullQ(Matrix tall) {
        RealMatrix rm = MatrixUtils.createRealMatrix(tall.toArray2D());
        RealMatrix q = new QRDecomposition(rm).getQ();
        return new Matrix(q.getData());
    }

    /** Solves A*X = B via partial-pivot LU (LAPACK-equivalent accuracy). */
    public static Matrix linsolve(Matrix A, Matrix B) {
        RealMatrix ra = MatrixUtils.createRealMatrix(A.toArray2D());
        RealMatrix rb = MatrixUtils.createRealMatrix(B.toArray2D());
        RealMatrix x = new LUDecomposition(ra).getSolver().solve(rb);
        return new Matrix(x.getData());
    }

    /**
     * Solves the Sylvester equation A*X - X*B = C where A (n x n) and B (m x m)
     * are quasi-upper-triangular (real Schur form), via column-wise
     * back-substitution (Bartels-Stewart), handling 1x1 and 2x2 diagonal blocks
     * of B. Accurate even when A/B are large.
     */
    public static Matrix triSylvester(Matrix A, Matrix B, Matrix C) {
        int n = A.getNumRows(), m = B.getNumRows();
        Matrix X = Matrix.zeros(n, m);
        Matrix In = Matrix.eye(n);
        int j = 0;
        while (j < m) {
            boolean two = (j < m - 1)
                    && Math.abs(B.get(j + 1, j)) > 1e-13 * (1 + Math.abs(B.get(j, j)) + Math.abs(B.get(j + 1, j + 1)));
            if (!two) {
                Matrix rhs = colVec(C, j);
                for (int k = 0; k < j; k++) rhs = rhs.add(colVec(X, k).scale(B.get(k, j)));
                Matrix xj = linsolve(A.sub(In.scale(B.get(j, j))), rhs);
                putCol(X, j, xj);
                j += 1;
            } else {
                Matrix rj = colVec(C, j), rj1 = colVec(C, j + 1);
                for (int k = 0; k < j; k++) {
                    rj = rj.add(colVec(X, k).scale(B.get(k, j)));
                    rj1 = rj1.add(colVec(X, k).scale(B.get(k, j + 1)));
                }
                Matrix LHS = Matrix.zeros(2 * n, 2 * n);
                setSub(LHS, range(n), range(n), A.sub(In.scale(B.get(j, j))));
                setSub(LHS, range(n), shiftRange(n, n), In.scale(-B.get(j + 1, j)));
                setSub(LHS, shiftRange(n, n), range(n), In.scale(-B.get(j, j + 1)));
                setSub(LHS, shiftRange(n, n), shiftRange(n, n), A.sub(In.scale(B.get(j + 1, j + 1))));
                Matrix xsol = linsolve(LHS, vstack(rj, rj1));
                putCol(X, j, sub(xsol, range(n), new int[]{0}));
                putCol(X, j + 1, sub(xsol, shiftRange(n, n), new int[]{0}));
                j += 2;
            }
        }
        return X;
    }

    private static Matrix colVec(Matrix M, int j) {
        Matrix c = Matrix.zeros(M.getNumRows(), 1);
        for (int i = 0; i < M.getNumRows(); i++) c.set(i, 0, M.get(i, j));
        return c;
    }

    private static void putCol(Matrix M, int j, Matrix col) {
        for (int i = 0; i < M.getNumRows(); i++) M.set(i, j, col.get(i, 0));
    }
}
