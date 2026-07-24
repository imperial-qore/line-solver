/**
 * ProCoMoM algorithm for computing marginal queue-length probabilities.
 */
package jline.api.pfqn.nc;

import jline.io.Ret;
import jline.util.Maths;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public final class Pfqn_procomom {

    private Pfqn_procomom() {
    }

    public static Ret.pfqnProcomom pfqn_procomom(Matrix L, Matrix N) {
        return pfqn_procomom(L, N, null, 1e-14);
    }

    public static Ret.pfqnProcomom pfqn_procomom(Matrix L, Matrix N, Matrix Z) {
        return pfqn_procomom(L, N, Z, 1e-14);
    }

    public static Ret.pfqnProcomom pfqn_procomom(Matrix L, Matrix N, Matrix Z, double atol) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        int sumN = (int) N.elementSum();

        Matrix Zmat;
        if (Z == null || Z.isEmpty()) {
            Zmat = new Matrix(1, R);
            Zmat.fill(0.0);
        } else {
            Zmat = Z.copy();
        }

        Matrix Lmax = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            double maxVal = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < M; i++) {
                if (L.get(i, r) > maxVal) maxVal = L.get(i, r);
            }
            if (maxVal < atol) maxVal = 1.0;
            Lmax.set(r, maxVal);
        }
        Matrix Ls = L.copy();
        Matrix Zs = Zmat.copy();
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                Ls.set(i, r, Ls.get(i, r) / Lmax.get(r));
            }
        }
        for (int r = 0; r < R; r++) {
            Zs.set(r, Zs.get(r) / Lmax.get(r));
        }

        List<int[]> DnRaw = multichooseP(R, M);
        for (int i = 0; i < DnRaw.size(); i++) {
            DnRaw.get(i)[R - 1] = 0;
        }
        List<int[]> Dn = sortByNnzPosP(DnRaw);
        int numDn = Dn.size();
        int basisSize = numDn * M;

        Pair<Matrix, Boolean> Pr = solveAll(Ls, Zs, Dn, N, M, R, sumN, basisSize, numDn);
        boolean rankdef = Pr.getRight();

        if (rankdef) {
            double Lscale = 0.0;
            for (int i = 0; i < Ls.getNumRows(); i++) {
                for (int j = 0; j < Ls.getNumCols(); j++) {
                    Lscale = FastMath.max(Lscale, FastMath.abs(Ls.get(i, j)));
                }
            }
            if (Lscale < atol) Lscale = 1.0;

            double[] deltaExps = new double[]{-10.0, -8.0, -6.0, -4.0};
            for (double deltaExp : deltaExps) {
                double delta = Lscale * FastMath.pow(10.0, deltaExp);
                Random rng2 = new Random(23000);
                Matrix Lp = Ls.copy();
                Matrix Zp = Zs.copy();
                for (int i = 0; i < M; i++) {
                    for (int r = 0; r < R; r++) {
                        Lp.set(i, r, Lp.get(i, r) + delta * (1 + rng2.nextDouble()));
                    }
                }
                for (int r = 0; r < R; r++) {
                    Zp.set(r, Zp.get(r) + delta * (1 + rng2.nextDouble()));
                }

                Pair<Matrix, Boolean> PrTry = solveAll(Lp, Zp, Dn, N, M, R, sumN, basisSize, numDn);
                if (!PrTry.getRight()) {
                    boolean allValid = true;
                    for (int i = 0; i < M; i++) {
                        for (int j = 0; j <= sumN; j++) {
                            if (PrTry.getLeft().get(i, j) < -1e-6) allValid = false;
                        }
                        double rowSum = 0.0;
                        for (int j = 0; j <= sumN; j++) rowSum += PrTry.getLeft().get(i, j);
                        if (FastMath.abs(rowSum - 1.0) > 0.01) allValid = false;
                    }
                    if (allValid) {
                        Pr = PrTry;
                        break;
                    }
                }
                Pr = PrTry;
            }
        }

        Matrix PrMat = Pr.getLeft();
        Matrix Q = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            double mean = 0.0;
            for (int j = 0; j <= sumN; j++) {
                mean += (double) j * PrMat.get(i, j);
            }
            Q.set(i, 0, mean);
        }

        return new Ret.pfqnProcomom(PrMat, Q);
    }

    private static Pair<Matrix, Boolean> solveAll(Matrix Ls, Matrix Zs, List<int[]> Dn,
                                                    Matrix N, int M, int R, int sumN, int basisSize, int numDn) {
        Matrix Pr = new Matrix(M, sumN + 1);
        boolean rankdef = false;

        for (int station = 0; station < M; station++) {
            Matrix Lrot = Ls.copy();
            for (int r = 0; r < R; r++) {
                double tmp = Lrot.get(station, r);
                Lrot.set(station, r, Lrot.get(M - 1, r));
                Lrot.set(M - 1, r, tmp);
            }

            Pair<double[], Boolean> result = solveStation(Lrot, Zs, Dn, N, M, R, sumN, basisSize, numDn);
            if (result.getRight()) rankdef = true;

            double[] dist = result.getLeft();
            double total = 0.0;
            for (double v : dist) total += v;
            if (FastMath.abs(total) > 0) {
                for (int j = 0; j <= sumN; j++) {
                    Pr.set(station, j, dist[j] / total);
                }
            }
        }

        return new Pair<Matrix, Boolean>(Pr, rankdef);
    }

    private static Pair<double[], Boolean> solveStation(Matrix Ls, Matrix Zs, List<int[]> Dn,
                                                          Matrix N, int M, int R, int sumN, int basisSize, int numDn) {
        boolean rankdef = false;

        double[][] pk = new double[basisSize][sumN + 1];

        int[] zeroDn = new int[R];
        for (int kk = 1; kk <= M; kk++) {
            int idx = phash(Dn, zeroDn, kk, M);
            if (idx >= 0 && idx < basisSize) {
                pk[idx][0] = 1.0;
            }
        }

        int[] Ncur = new int[R];
        for (int r = 0; r < R; r++) {
            for (int Nr = 1; Nr <= (int) N.get(r); Nr++) {
                Ncur[r] = Nr;
                double[][] pklast = new double[basisSize][];
                for (int i = 0; i < basisSize; i++) pklast[i] = pk[i].clone();
                pk = new double[basisSize][sumN + 1];

                double[][][] matrices = genPMatrix(Ls, Zs, Ncur, r, Dn, M, R, numDn, basisSize);
                double[][] Ag = matrices[0];
                double[][] Bg = matrices[1];
                double[][] DCg = matrices[2];
                double[][] DDg = matrices[3];
                int numRows = Ag.length;

                Triple<double[][], double[], double[][]> svdResult = computeSVD(Ag, numRows, basisSize);
                double[][] U = svdResult.first;
                double[] sv = svdResult.second;
                double[][] V = svdResult.third;

                double tol = (double) FastMath.max(numRows, basisSize) * FastMath.ulp(1.0) * sv[0];
                int rk = 0;
                for (double s : sv) {
                    if (s > tol) rk++;
                }

                int sumNcur = sumArrayP(Ncur);

                if (rk < basisSize) {
                    rankdef = true;

                    double[][] pB = pseudoInvMult(U, sv, V, rk, Bg, numRows, basisSize);
                    double[][] pDC = pseudoInvMult(U, sv, V, rk, DCg, numRows, basisSize);
                    double[][] pDD = pseudoInvMult(U, sv, V, rk, DDg, numRows, basisSize);

                    for (int i = 0; i < basisSize; i++) {
                        double sum = 0.0;
                        for (int j = 0; j < basisSize; j++) {
                            sum += pB[i][j] * pklast[j][0];
                        }
                        pk[i][0] = sum;
                    }

                    for (int n = 1; n <= sumNcur; n++) {
                        for (int i = 0; i < basisSize; i++) {
                            double sum = 0.0;
                            for (int j = 0; j < basisSize; j++) {
                                sum += pB[i][j] * pklast[j][n];
                                sum += (double) n * pDC[i][j] * pk[j][n - 1];
                            }
                            pk[i][n] = sum;
                        }
                        for (int i = 0; i < basisSize; i++) {
                            double sum = 0.0;
                            for (int j = 0; j < basisSize; j++) {
                                sum += (double) n * pDD[i][j] * pklast[j][n - 1];
                            }
                            pk[i][n] += sum;
                        }
                    }
                } else {
                    Pair<double[][], double[][]> qr = computeQR(Ag, numRows, basisSize);
                    double[][] Q = qr.getLeft();
                    double[][] Rmat = qr.getRight();

                    double[][] QtB = transposeMultiply(Q, Bg, numRows, basisSize);
                    double[][] QtDC = transposeMultiply(Q, DCg, numRows, basisSize);
                    double[][] QtDD = transposeMultiply(Q, DDg, numRows, basisSize);

                    double[] rhsInit = matVecMult(QtB, pklast, 0, basisSize);
                    double[] solInit = backSubstitute(Rmat, rhsInit, basisSize);
                    for (int i = 0; i < basisSize; i++) pk[i][0] = solInit[i];

                    for (int n = 1; n <= sumNcur; n++) {
                        double[] rhs = new double[basisSize];
                        double[] pklastVec = matVecMult(QtB, pklast, n, basisSize);
                        double[] pkPrevVec = matVecMult(QtDC, pk, n - 1, basisSize);
                        double[] pklastPrevVec = matVecMult(QtDD, pklast, n - 1, basisSize);
                        for (int i = 0; i < basisSize; i++) {
                            rhs[i] = pklastVec[i] + (double) n * pkPrevVec[i] + (double) n * pklastPrevVec[i];
                        }
                        double[] sol = backSubstitute(Rmat, rhs, basisSize);
                        for (int i = 0; i < basisSize; i++) pk[i][n] = sol[i];
                    }
                }

                double smax = 0.0;
                for (int i = 0; i < basisSize; i++) {
                    for (int j = 0; j <= sumN; j++) {
                        double v = FastMath.abs(pk[i][j]);
                        if (v > smax) smax = v;
                    }
                }
                if (smax > 0 && Double.isFinite(smax)) {
                    for (int i = 0; i < basisSize; i++) {
                        for (int j = 0; j <= sumN; j++) {
                            pk[i][j] /= smax;
                        }
                    }
                }
            }
        }

        int[] zeroDnFinal = new int[R];
        int idx = phash(Dn, zeroDnFinal, 1, M);
        double[] dist = new double[sumN + 1];
        if (idx >= 0 && idx < basisSize) {
            for (int j = 0; j <= sumN; j++) {
                dist[j] = pk[idx][j];
            }
        }

        return new Pair<double[], Boolean>(dist, rankdef);
    }

    private static class Triple<A, B, C> {
        final A first;
        final B second;
        final C third;
        Triple(A a, B b, C c) {
            this.first = a;
            this.second = b;
            this.third = c;
        }
    }

    private static double[][][] genPMatrix(Matrix Ls, Matrix Zs, int[] Ncur, int r,
                                             List<int[]> Dn, int M, int R, int numDn, int basisSize) {
        int numRows = countRows(Dn, Ncur, r, M, R, numDn);
        double[][] A = new double[numRows][basisSize];
        double[][] B = new double[numRows][basisSize];
        double[][] DC = new double[numRows][basisSize];
        double[][] DD = new double[numRows][basisSize];

        int row = 0;
        for (int d = 0; d < numDn; d++) {
            int[] dn = Dn.get(d);

            int sumRange = 0;
            if (r <= R - 2) {
                for (int t = r; t < R - 1; t++) sumRange += dn[t];
            }

            if (r <= R - 2 && sumRange > 0) {
                for (int k = 1; k <= M; k++) {
                    if (row >= numRows) break;
                    int colA = phash(Dn, dn, k, M);
                    if (colA >= 0 && colA < basisSize) A[row][colA] = 1.0;

                    int sumRange2 = 0;
                    if (r + 1 <= R - 2) {
                        for (int t = r + 1; t < R - 1; t++) sumRange2 += dn[t];
                    }

                    if (r + 1 <= R - 2 && sumRange2 > 0) {
                        if (colA >= 0 && colA < basisSize) B[row][colA] = 1.0;
                    } else {
                        int[] shifted = dn.clone();
                        shifted[r]--;
                        int colB = phash(Dn, shifted, k, M);
                        if (colB >= 0 && colB < basisSize) B[row][colB] = 1.0;
                    }
                    row++;
                }
            } else {
                int sumDn1ToR = 0;
                for (int t = 0; t <= r; t++) sumDn1ToR += dn[t];

                if (sumDn1ToR < M) {
                    for (int k = 1; k < M; k++) {
                        if (row >= numRows) break;
                        int colKP1 = phash(Dn, dn, k + 1, M);
                        if (colKP1 >= 0 && colKP1 < basisSize) A[row][colKP1] = 1.0;
                        int col1 = phash(Dn, dn, 1, M);
                        if (col1 >= 0 && col1 < basisSize) A[row][col1] = -1.0;

                        for (int s = 0; s < r; s++) {
                            int[] shifted = dn.clone();
                            shifted[s]++;
                            int col = phash(Dn, shifted, k + 1, M);
                            if (col >= 0 && col < basisSize) {
                                A[row][col] -= Ls.get(k - 1, s);
                            }
                        }
                        int colB = phash(Dn, dn, k + 1, M);
                        if (colB >= 0 && colB < basisSize) B[row][colB] = Ls.get(k - 1, r);
                        row++;
                    }

                    for (int s = 0; s < r; s++) {
                        if (row >= numRows) break;
                        int ndS = Ncur[s] - dn[s];
                        int col1 = phash(Dn, dn, 1, M);
                        if (col1 >= 0 && col1 < basisSize) A[row][col1] = (double) ndS;

                        int[] shifted = dn.clone();
                        shifted[s]++;
                        int colBase = phash(Dn, shifted, 1, M);
                        if (colBase >= 0 && colBase < basisSize) {
                            A[row][colBase] -= Zs.get(s);
                            DC[row][colBase] = Ls.get(M - 1, s);
                        }
                        for (int k = 1; k < M; k++) {
                            int colK = phash(Dn, shifted, k + 1, M);
                            if (colK >= 0 && colK < basisSize) {
                                A[row][colK] -= Ls.get(k - 1, s);
                            }
                        }
                        row++;
                    }
                }

                if (row < numRows) {
                    int ndR = Ncur[r] - dn[r];
                    int col1 = phash(Dn, dn, 1, M);
                    if (col1 >= 0 && col1 < basisSize) {
                        A[row][col1] = (double) ndR;
                        B[row][col1] = Zs.get(r);
                    }
                    for (int k = 1; k < M; k++) {
                        int colKP1 = phash(Dn, dn, k + 1, M);
                        if (colKP1 >= 0 && colKP1 < basisSize) {
                            B[row][colKP1] = Ls.get(k - 1, r);
                        }
                    }
                    // see _kb/03-api-layer.md for rationale
                    if (col1 < 0 || col1 >= basisSize)
                        throw new RuntimeException(
                                "pfqn_procomom: population vector outside the CoMoM basis");
                    DD[row][col1] = Ls.get(M - 1, r);
                    row++;
                }
            }
        }

        return new double[][][]{A, B, DC, DD};
    }

    private static int countRows(List<int[]> Dn, int[] Ncur, int r, int M, int R, int numDn) {
        int numRows = 0;
        for (int d = 0; d < numDn; d++) {
            int[] dn = Dn.get(d);
            int sumRange = 0;
            if (r <= R - 2) {
                for (int t = r; t < R - 1; t++) sumRange += dn[t];
            }
            if (r <= R - 2 && sumRange > 0) {
                numRows += M;
            } else {
                int sumDn1ToR = 0;
                for (int t = 0; t <= r; t++) sumDn1ToR += dn[t];
                if (sumDn1ToR < M) {
                    numRows += M + r - 1;
                }
                numRows += 1;
            }
        }
        return numRows;
    }

    private static int phash(List<int[]> Dn, int[] dn, int k, int M) {
        int pos = matchRowP(Dn, dn);
        if (pos < 0) return -1;
        return pos * M + k - 1;
    }

    private static int matchRowP(List<int[]> Dn, int[] row) {
        for (int i = 0; i < Dn.size(); i++) {
            boolean match = true;
            for (int j = 0; j < row.length; j++) {
                if (Dn.get(i)[j] != row[j]) {
                    match = false;
                    break;
                }
            }
            if (match) return i;
        }
        return -1;
    }

    private static List<int[]> multichooseP(int R, int M) {
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
            List<int[]> sub = multichooseP(R - 1, M - i);
            for (int[] s : sub) {
                int[] row = new int[R];
                row[0] = i;
                System.arraycopy(s, 0, row, 1, R - 1);
                result.add(row);
            }
        }
        return result;
    }

    private static List<int[]> sortByNnzPosP(List<int[]> list) {
        List<int[]> result = new ArrayList<int[]>(list);
        for (int i = 0; i < result.size() - 1; i++) {
            for (int j = i + 1; j < result.size(); j++) {
                if (nnzCmpP(result.get(i), result.get(j)) == 1) {
                    int[] tmp = result.get(i);
                    result.set(i, result.get(j));
                    result.set(j, tmp);
                }
            }
        }
        return result;
    }

    private static int nnzCmpP(int[] i1, int[] i2) {
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

    private static int sumArrayP(int[] arr) {
        int s = 0;
        for (int v : arr) s += v;
        return s;
    }

    private static Triple<double[][], double[], double[][]> computeSVD(double[][] A, int numRows, int numCols) {
        org.ejml.interfaces.decomposition.SingularValueDecomposition_F64<org.ejml.data.DMatrixRMaj> svd =
                org.ejml.dense.row.factory.DecompositionFactory_DDRM.svd(numRows, numCols, true, true, false);
        org.ejml.data.DMatrixRMaj ddrm = new org.ejml.data.DMatrixRMaj(numRows, numCols);
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                ddrm.set(i, j, A[i][j]);
            }
        }
        svd.decompose(ddrm);
        int minDim = FastMath.min(numRows, numCols);
        double[] sv = new double[minDim];
        org.ejml.data.DMatrixRMaj Umat = svd.getU(null, false);
        org.ejml.data.DMatrixRMaj Vmat = svd.getV(null, false);
        double[] svals = svd.getSingularValues();
        System.arraycopy(svals, 0, sv, 0, minDim);

        double[][] U = new double[numRows][minDim];
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < minDim; j++) U[i][j] = Umat.get(i, j);
        }
        double[][] V = new double[numCols][minDim];
        for (int i = 0; i < numCols; i++) {
            for (int j = 0; j < minDim; j++) V[i][j] = Vmat.get(i, j);
        }

        return new Triple<double[][], double[], double[][]>(U, sv, V);
    }

    private static Pair<double[][], double[][]> computeQR(double[][] A, int numRows, int numCols) {
        org.ejml.data.DMatrixRMaj ddrm = new org.ejml.data.DMatrixRMaj(numRows, numCols);
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                ddrm.set(i, j, A[i][j]);
            }
        }
        org.ejml.interfaces.decomposition.QRDecomposition<org.ejml.data.DMatrixRMaj> qr =
                org.ejml.dense.row.factory.DecompositionFactory_DDRM.qr(numRows, numCols);
        qr.decompose(ddrm);
        org.ejml.data.DMatrixRMaj Qmat = qr.getQ(null, true);
        org.ejml.data.DMatrixRMaj Rmat = qr.getR(null, true);

        int minDim = FastMath.min(numRows, numCols);
        double[][] Q = new double[numRows][minDim];
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < minDim; j++) Q[i][j] = Qmat.get(i, j);
        }
        double[][] R = new double[minDim][numCols];
        for (int i = 0; i < minDim; i++) {
            for (int j = 0; j < numCols; j++) R[i][j] = Rmat.get(i, j);
        }

        return new Pair<double[][], double[][]>(Q, R);
    }

    private static double[][] pseudoInvMult(double[][] U, double[] sv, double[][] V,
                                              int rk, double[][] B, int numRows, int numCols) {
        double[][] result = new double[numCols][numCols];

        double[][] UrTB = new double[rk][numCols];
        for (int i = 0; i < rk; i++) {
            for (int j = 0; j < numCols; j++) {
                double sum = 0.0;
                for (int k = 0; k < numRows; k++) {
                    sum += U[k][i] * B[k][j];
                }
                UrTB[i][j] = sum;
            }
        }

        for (int i = 0; i < rk; i++) {
            for (int j = 0; j < numCols; j++) {
                UrTB[i][j] /= sv[i];
            }
        }

        for (int i = 0; i < numCols; i++) {
            for (int j = 0; j < numCols; j++) {
                double sum = 0.0;
                for (int k = 0; k < rk; k++) {
                    sum += V[i][k] * UrTB[k][j];
                }
                result[i][j] = sum;
            }
        }

        return result;
    }

    private static double[][] transposeMultiply(double[][] Q, double[][] B, int numRows, int numCols) {
        double[][] result = new double[numCols][numCols];
        for (int i = 0; i < numCols; i++) {
            for (int j = 0; j < numCols; j++) {
                double sum = 0.0;
                for (int k = 0; k < numRows; k++) {
                    sum += Q[k][i] * B[k][j];
                }
                result[i][j] = sum;
            }
        }
        return result;
    }

    private static double[] matVecMult(double[][] mat, double[][] pk, int col, int n) {
        double[] result = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += mat[i][j] * pk[j][col];
            }
            result[i] = sum;
        }
        return result;
    }

    private static double[] backSubstitute(double[][] R, double[] b, int n) {
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = b[i];
            for (int j = i + 1; j < n; j++) {
                sum -= R[i][j] * x[j];
            }
            // see _kb/03-api-layer.md for rationale
            if (FastMath.abs(R[i][i]) < 1e-30)
                throw new RuntimeException("pfqn_procomom: singular triangular factor");
            x[i] = sum / R[i][i];
        }
        return x;
    }
}
