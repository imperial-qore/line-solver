/**
 * @file CoMoM algorithm for computing the normalizing constant.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;

import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

public final class Pfqn_comom {
    private Pfqn_comom() {}

    /**
     * Triple of (A, B, DA) coefficient matrices.
     */
    public static final class GenMatrixResult {
        public final double[][] A;
        public final double[][] B;
        public final double[][] DA;

        public GenMatrixResult(double[][] A, double[][] B, double[][] DA) {
            this.A = A;
            this.B = B;
            this.DA = DA;
        }
    }

    public static double pfqn_comom(Matrix L, Matrix N, Matrix Z, double atol) {
        Matrix Zmat = (Z != null) ? Z : new Matrix(1, N.getNumCols());
        int M = L.getNumRows();
        int R = L.getNumCols();
        if (M > 1) {
            // Same contract as MATLAB pfqn_comom: the COMOM recursion here is
            // the repairman-model variant; for M>1 it silently returns wrong
            // values, so reject such inputs (use pfqn_ca / pfqn_recal instead)
            throw new IllegalArgumentException(
                    "pfqn_comom supports at most one queueing station (repairman models with a delay); "
                            + "use pfqn_ca or pfqn_recal for M>1.");
        }

        Matrix Lmax = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            double maxVal = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < M; i++) {
                if (L.get(i, r) > maxVal) maxVal = L.get(i, r);
            }
            if (maxVal < atol) maxVal = Zmat.get(r);
            Lmax.set(r, maxVal);
        }

        Matrix Ls = L.copy();
        Matrix Zs = Zmat.copy();
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                if (Lmax.get(r) > 0) {
                    Ls.set(i, r, Ls.get(i, r) / Lmax.get(r));
                }
            }
        }
        for (int r = 0; r < R; r++) {
            if (Lmax.get(r) > 0) {
                Zs.set(r, Zs.get(r) / Lmax.get(r));
            }
        }

        List<int[]> DnRaw = multichoose(R, M);
        for (int i = 0; i < DnRaw.size(); i++) {
            DnRaw.get(i)[R - 1] = 0;
        }
        List<int[]> Dn = sortByNnzPos(DnRaw, R);
        int numDn = Dn.size();

        int[] nvec = new int[R];
        int totalCols = numDn * (M + 1);
        double[] h = new double[totalCols];
        for (int i = 0; i <= M; i++) {
            int col = hash(Dn, nvec, nvec, i + 1, M, R);
            h[col] = 1.0;
        }

        int sumNvec = sumArray(nvec);
        double logScale = Maths.factln(sumNvec + M - 1) - sumFactln(nvec);
        for (int i = 0; i < totalCols; i++) {
            if (h[i] > 0) {
                h[i] = FastMath.exp(FastMath.log(h[i]) + logScale);
            }
        }

        double[] scale = new double[sumArrayInt(N)];

        for (int r = 0; r < R; r++) {
            double[][] A = null;
            double[][] B = null;
            double[][] DA = null;

            for (int Nr = 1; Nr <= (int) N.get(r); Nr++) {
                nvec[r]++;

                if (Nr == 1) {
                    GenMatrixResult matrices = genMatrix(Ls, nvec, Zs, r, Dn, M, R, N);
                    A = matrices.A;
                    B = matrices.B;
                    DA = matrices.DA;
                } else {
                    for (int i = 0; i < A.length; i++) {
                        for (int j = 0; j < A[0].length; j++) {
                            A[i][j] += DA[i][j];
                        }
                    }
                }

                int sz = A.length;
                double coeff = (double) nvec[r] / (sumArray(nvec) + M - 1);
                double[] b = new double[sz];
                for (int i = 0; i < sz; i++) {
                    double sum = 0.0;
                    for (int j = 0; j < sz; j++) {
                        sum += B[i][j] * h[j];
                    }
                    b[i] = sum * coeff;
                }

                double[] hNew = solveLinearSystem(A, b);
                for (int i = 0; i < sz; i++) {
                    h[i] = hNew[i];
                }

                int nt = sumArray(nvec);
                double[] sortedH = h.clone();
                Arrays.sort(sortedH);
                double absSum = 0.0;
                for (double v : sortedH) absSum += v;
                scale[nt - 1] = FastMath.abs(absSum);

                for (int i = 0; i < h.length; i++) {
                    h[i] = FastMath.abs(h[i]) / scale[nt - 1];
                }
            }
        }

        int sumN = sumArrayInt(N);
        double lG = FastMath.log(h[totalCols - R]);
        lG += Maths.factln(sumN + M - 1);
        for (int r = 0; r < R; r++) {
            lG -= Maths.factln((int) N.get(r));
        }
        for (int r = 0; r < R; r++) {
            lG += N.get(r) * FastMath.log(Lmax.get(r));
        }
        for (int i = 0; i < sumN; i++) {
            lG += FastMath.log(scale[i]);
        }

        return lG;
    }

    private static GenMatrixResult genMatrix(
            Matrix L, int[] nvec, Matrix Z, int r,
            List<int[]> Dn, int M, int R, Matrix N) {
        int numDn = Dn.size();
        int sz = numDn * (M + 1);
        double[][] A = new double[sz][sz];
        double[][] DA = new double[sz][sz];
        double[][] B = new double[sz][sz];

        int row = 0;
        for (int d = 0; d < numDn; d++) {
            int[] dn = Dn.get(d);

            int sumDnRToEnd = 0;
            for (int t = r; t < R - 1; t++) {
                sumDnRToEnd += dn[t];
            }

            if (sumDnRToEnd > 0) {
                for (int k = 0; k <= M; k++) {
                    int[] nMinusDn = new int[R];
                    for (int t = 0; t < R; t++) nMinusDn[t] = nvec[t] - dn[t];
                    int col = hash(Dn, nvec, nMinusDn, k + 1, M, R);
                    A[row][col] = 1.0;

                    int sumDnR1ToEnd = 0;
                    for (int t = r + 1; t < R - 1; t++) sumDnR1ToEnd += dn[t];

                    if (sumDnR1ToEnd > 0) {
                        B[row][col] = 1.0;
                    } else {
                        int[] er = new int[R];
                        er[r] = 1;
                        int[] nMinusDnPlusEr = new int[R];
                        for (int t = 0; t < R; t++) nMinusDnPlusEr[t] = nvec[t] - dn[t] + er[t];
                        int colB = hash(Dn, nvec, nMinusDnPlusEr, k + 1, M, R);
                        B[row][colB] = 1.0;
                    }
                    row++;
                }
            } else {
                // see _kb/03-api-layer.md for rationale
                int sumDn1ToR = 0;
                for (int t = 0; t <= r && t < R; t++) sumDn1ToR += dn[t];

                if (sumDn1ToR < M) {
                    for (int k = 1; k <= M; k++) {
                        int[] nMinusDn = new int[R];
                        for (int t = 0; t < R; t++) nMinusDn[t] = nvec[t] - dn[t];

                        int colKP1 = hash(Dn, nvec, nMinusDn, k + 1, M, R);
                        A[row][colKP1] = 1.0;

                        int col0P1 = hash(Dn, nvec, nMinusDn, 0 + 1, M, R);
                        A[row][col0P1] = -1.0;

                        for (int s = 0; s < r; s++) {
                            int[] nOner = new int[R];
                            for (int t = 0; t < R; t++) nOner[t] = nvec[t] - dn[t];
                            nOner[s]--;
                            int colOner = hash(Dn, nvec, nOner, k + 1, M, R);
                            A[row][colOner] = A[row][colOner] - L.get(k - 1, s);
                        }

                        B[row][colKP1] = L.get(k - 1, r);
                        row++;
                    }

                    for (int s = 0; s < r; s++) {
                        int[] nMinusDn = new int[R];
                        for (int t = 0; t < R; t++) nMinusDn[t] = nvec[t] - dn[t];

                        int col0P1 = hash(Dn, nvec, nMinusDn, 0 + 1, M, R);
                        A[row][col0P1] = (double) nMinusDn[s];

                        int[] nOner = new int[R];
                        for (int t = 0; t < R; t++) nOner[t] = nMinusDn[t];
                        nOner[s]--;

                        int colOner0P1 = hash(Dn, nvec, nOner, 0 + 1, M, R);
                        A[row][colOner0P1] = A[row][colOner0P1] - Z.get(s);
                        for (int k = 1; k <= M; k++) {
                            int colOnerKP1 = hash(Dn, nvec, nOner, k + 1, M, R);
                            A[row][colOnerKP1] = A[row][colOnerKP1] - L.get(k - 1, s);
                        }
                        row++;
                    }
                }

                int[] nMinusDn = new int[R];
                for (int t = 0; t < R; t++) nMinusDn[t] = nvec[t] - dn[t];
                int col0P1 = hash(Dn, nvec, nMinusDn, 0 + 1, M, R);
                A[row][col0P1] = (double) nMinusDn[r];
                DA[row][col0P1] = 1.0;
                B[row][col0P1] = Z.get(r);
                for (int k = 1; k <= M; k++) {
                    int colKP1 = hash(Dn, nvec, nMinusDn, k + 1, M, R);
                    B[row][colKP1] = L.get(k - 1, r);
                }
                row++;
            }
        }

        return new GenMatrixResult(A, B, DA);
    }

    /**
     * Column index of the normalizing constant G(ref - n) with i-1 extra jobs at
     * the queueing stations. The reference vector is the RUNNING population nvec,
     * not the target population N: MATLAB's genmatrix shadows the outer N with its
     * own parameter, which is called as genmatrix(L,nvec,Z,r), so every hash there
     * resolves ref - n to a row of Dn by construction. Passing the full N instead
     * makes matchRow miss on every call except the last iteration, and the -1 was
     * then swallowed by bounds guards, leaving h all zeros and lG = -Inf.
     */
    private static int hash(List<int[]> Dn, int[] ref, int[] n, int i, int M, int R) {
        int[] diff = new int[R];
        for (int t = 0; t < R; t++) {
            diff[t] = ref[t] - n[t];
        }
        int pos = matchRow(Dn, diff);
        if (pos < 0)
            throw new RuntimeException("pfqn_comom: population vector outside the CoMoM basis");
        if (i == 1) {
            return Dn.size() * M + pos;
        } else {
            return pos * M + i - 2;
        }
    }

    private static int matchRow(List<int[]> Dn, int[] row) {
        for (int i = 0; i < Dn.size(); i++) {
            int[] r = Dn.get(i);
            boolean match = true;
            for (int j = 0; j < row.length; j++) {
                if (r[j] != row[j]) {
                    match = false;
                    break;
                }
            }
            if (match) return i;
        }
        return -1;
    }

    private static List<int[]> multichoose(int R, int M) {
        List<int[]> result = new ArrayList<int[]>();
        if (R == 1) {
            result.add(new int[]{M});
            return result;
        }
        if (M == 0) {
            result.add(new int[R]);
            return result;
        }
        for (int i = 0; i <= M; i++) {
            List<int[]> sub = multichoose(R - 1, M - i);
            for (int[] s : sub) {
                int[] row = new int[R];
                row[0] = i;
                System.arraycopy(s, 0, row, 1, R - 1);
                result.add(row);
            }
        }
        return result;
    }

    private static List<int[]> sortByNnzPos(List<int[]> list, int R) {
        List<int[]> result = new ArrayList<int[]>(list);
        for (int i = 0; i < result.size() - 1; i++) {
            for (int j = i + 1; j < result.size(); j++) {
                if (nnzCmp(result.get(i), result.get(j)) == 1) {
                    int[] tmp = result.get(i);
                    result.set(i, result.get(j));
                    result.set(j, tmp);
                }
            }
        }
        return result;
    }

    private static int nnzCmp(int[] i1, int[] i2) {
        int nnz1 = 0;
        int nnz2 = 0;
        for (int v : i1) if (v != 0) nnz1++;
        for (int v : i2) if (v != 0) nnz2++;
        if (nnz1 > nnz2) return 1;
        if (nnz1 < nnz2) return 0;
        for (int j = 0; j < i1.length; j++) {
            if (i1[j] == 0 && i2[j] > 0) return 1;
            if (i1[j] > 0 && i2[j] == 0) return 0;
        }
        return 0;
    }

    private static int sumArray(int[] arr) {
        int s = 0;
        for (int v : arr) s += v;
        return s;
    }

    private static int sumArrayInt(Matrix N) {
        return (int) N.elementSum();
    }

    private static double sumFactln(int[] arr) {
        double s = 0.0;
        for (int v : arr) s += Maths.factln(v);
        return s;
    }

    private static double[] solveLinearSystem(double[][] A, double[] b) {
        int n = b.length;
        double[][] aug = new double[n][n + 1];
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, aug[i], 0, n);
            aug[i][n] = b[i];
        }

        for (int col = 0; col < n; col++) {
            double maxVal = FastMath.abs(aug[col][col]);
            int maxRow = col;
            for (int row = col + 1; row < n; row++) {
                if (FastMath.abs(aug[row][col]) > maxVal) {
                    maxVal = FastMath.abs(aug[row][col]);
                    maxRow = row;
                }
            }
            if (maxRow != col) {
                double[] tmp = aug[col];
                aug[col] = aug[maxRow];
                aug[maxRow] = tmp;
            }

            if (FastMath.abs(aug[col][col]) < 1e-30)
                throw new RuntimeException("pfqn_comom: singular CoMoM system");

            for (int row = col + 1; row < n; row++) {
                double factor = aug[row][col] / aug[col][col];
                for (int j = col; j <= n; j++) {
                    aug[row][j] -= factor * aug[col][j];
                }
            }
        }

        double[] x = new double[n];
        for (int row = n - 1; row >= 0; row--) {
            if (FastMath.abs(aug[row][row]) < 1e-30)
                throw new RuntimeException("pfqn_comom: singular CoMoM system");
            double sum = aug[row][n];
            for (int col = row + 1; col < n; col++) {
                sum -= aug[row][col] * x[col];
            }
            x[row] = sum / aug[row][row];
        }
        return x;
    }
}
