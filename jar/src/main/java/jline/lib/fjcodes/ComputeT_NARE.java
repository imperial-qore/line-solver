package jline.lib.fjcodes;

import jline.util.matrix.Matrix;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * Compute T-matrix using NARE (Nonsymmetric Algebraic Riccati Equation) method.
 */
public final class ComputeT_NARE {
    private ComputeT_NARE() {}

    public static Matrix computeT_NARE(Matrix D0, Matrix D1, Matrix S, Matrix A_jump) {
        int m = S.getNumRows();
        int ma = D0.getNumRows();
        int ms = m / ma;

        Matrix Ims = Matrix.eye(ms);
        Matrix Ima = Matrix.eye(ma);

        Matrix A = Ims.kron(D0);
        Matrix B = Ims.kron(D1);
        Matrix C = A_jump.kron(Ima);
        Matrix D = S;

        Matrix H = new Matrix(2 * m, 2 * m);
        setSubMatrix(H, 0, 0, A);
        setSubMatrix(H, 0, m, B);
        setSubMatrix(H, m, 0, C.scale(-1.0));
        setSubMatrix(H, m, m, D.scale(-1.0));

        try {
            // Reorder the real Schur form of H so the m stable (negative real-part)
            // eigenvalues occupy the leading columns, i.e. the leading m columns of Z
            // span the stable invariant subspace (MATLAB: ordschur selecting the m
            // smallest eigenvalues). FluidTools.orderedSchur groups eigenvalues as
            // [zero, negative, positive] with a block-aware swap that correctly
            // handles 2x2 (complex-conjugate) blocks, unlike the previous 1x1-only
            // reordering here, which corrupted T for Hamiltonians with complex
            // eigenvalues. For a stable FJ queue there are no zero eigenvalues and
            // exactly m negative ones, so columns [0, m) are the stable subspace.
            Object[] os = jline.lib.butools.mam.FluidTools.orderedSchur(H, 1e-8);
            Matrix Z = (Matrix) os[0];
            int[] counts = (int[]) os[2];
            int lead = counts[0];  // number of zero-real-part eigenvalues (normally 0)

            // Stable subspace = columns [lead, lead+m); rows split at the physical
            // midpoint m (upper/lower halves of the 2m x 2m orthogonal factor).
            Matrix Q1_11 = Matrix.getSubMatrix(Z, 0, m, lead, lead + m);
            Matrix Q1_21 = Matrix.getSubMatrix(Z, m, 2 * m, lead, lead + m);

            Matrix X = Q1_21.rightMatrixDivide(Q1_11);
            Matrix T2 = S.add(1.0, X.mult(Ims.kron(D1)));
            return T2;
        } catch (Exception e) {
            throw new RuntimeException("computeT_NARE failed", e);
        }
    }

    private static void setSubMatrix(Matrix dest, int rowOff, int colOff, Matrix src) {
        for (int i = 0; i < src.getNumRows(); i++) {
            for (int j = 0; j < src.getNumCols(); j++) {
                dest.set(rowOff + i, colOff + j, src.get(i, j));
            }
        }
    }

    private static List<Complex> extractSchurEigenvalues(Matrix Q) {
        int n = Q.getNumRows();
        List<Complex> eigenvalues = new ArrayList<Complex>();
        int i = 0;
        while (i < n) {
            if (i < n - 1 && Math.abs(Q.get(i + 1, i)) > 1e-10) {
                double a = Q.get(i, i);
                double b = Q.get(i, i + 1);
                double c = Q.get(i + 1, i);
                double d = Q.get(i + 1, i + 1);
                double trace = a + d;
                double det = a * d - b * c;
                double discriminant = trace * trace / 4.0 - det;
                if (discriminant >= 0) {
                    double sqrtDisc = Math.sqrt(discriminant);
                    eigenvalues.add(new Complex(trace / 2.0 + sqrtDisc, 0.0));
                    eigenvalues.add(new Complex(trace / 2.0 - sqrtDisc, 0.0));
                } else {
                    double sqrtDisc = Math.sqrt(-discriminant);
                    eigenvalues.add(new Complex(trace / 2.0, sqrtDisc));
                    eigenvalues.add(new Complex(trace / 2.0, -sqrtDisc));
                }
                i += 2;
            } else {
                eigenvalues.add(new Complex(Q.get(i, i), 0.0));
                i++;
            }
        }
        return eigenvalues;
    }

    private static Matrix orderedSchur(Matrix U, Matrix T, boolean[] sel) {
        int n = U.getNumRows();
        Matrix U1 = U.copy();
        Matrix T1 = T.copy();
        boolean changed = true;
        int iter = 0;
        int maxIter = n * n;
        while (changed && iter < maxIter) {
            changed = false;
            iter++;
            for (int i = 0; i < n - 1; i++) {
                if (!sel[i] && sel[i + 1]) {
                    swapAdjacentEigenvalues(U1, T1, i);
                    boolean tmp = sel[i];
                    sel[i] = sel[i + 1];
                    sel[i + 1] = tmp;
                    changed = true;
                }
            }
        }
        return U1;
    }

    private static void swapAdjacentEigenvalues(Matrix U, Matrix T, int k) {
        int n = T.getNumRows();
        double t11 = T.get(k, k);
        double t22 = T.get(k + 1, k + 1);
        if (Math.abs(T.get(k + 1, k)) < 1e-14) {
            double t12 = T.get(k, k + 1);
            if (Math.abs(t12) < 1e-14) {
                T.set(k, k, t22);
                T.set(k + 1, k + 1, t11);
                return;
            }
        }
        double t12 = T.get(k, k + 1);
        double cs;
        double sn;
        if (Math.abs(t11 - t22) < 1e-14) {
            cs = 1.0; sn = 0.0;
        } else {
            double sign = (t22 - t11 > 0) ? 1.0 : -1.0;
            double temp = (t22 - t11) / (2.0 * t12);
            double tau = sign / (Math.abs(temp) + Math.sqrt(1.0 + temp * temp));
            cs = 1.0 / Math.sqrt(1.0 + tau * tau);
            sn = tau * cs;
        }
        for (int i = 0; i < n; i++) {
            double temp1 = T.get(i, k);
            double temp2 = T.get(i, k + 1);
            T.set(i, k, cs * temp1 + sn * temp2);
            T.set(i, k + 1, -sn * temp1 + cs * temp2);
        }
        for (int j = 0; j < n; j++) {
            double temp1 = T.get(k, j);
            double temp2 = T.get(k + 1, j);
            T.set(k, j, cs * temp1 + sn * temp2);
            T.set(k + 1, j, -sn * temp1 + cs * temp2);
        }
        for (int i = 0; i < n; i++) {
            double temp1 = U.get(i, k);
            double temp2 = U.get(i, k + 1);
            U.set(i, k, cs * temp1 + sn * temp2);
            U.set(i, k + 1, -sn * temp1 + cs * temp2);
        }
    }

    private static double normInf(Matrix A) {
        double maxSum = 0.0;
        for (int i = 0; i < A.getNumRows(); i++) {
            double rowSum = 0.0;
            for (int j = 0; j < A.getNumCols(); j++) {
                rowSum += Math.abs(A.get(i, j));
            }
            maxSum = Math.max(maxSum, rowSum);
        }
        return maxSum;
    }
}
