/**
 * Method of Moments (MOM) for exact normalizing constant computation
 *
 * Implements the Method of Moments using exact arithmetic with BigFraction.
 */
package jline.api.pfqn.nc;

import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.fraction.BigFraction;
import org.apache.commons.math3.fraction.BigFractionField;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;
import org.apache.commons.math3.linear.ArrayFieldVector;
import org.apache.commons.math3.linear.FieldLUDecomposition;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.FieldVector;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;

/**
 * Method of Moments (MOM) solver for product-form queueing networks.
 */
public final class Pfqn_mom {

    private Pfqn_mom() {
    }

    public static Ret.pfqnMom pfqn_mom(Matrix L, Matrix N, Matrix Z) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        // see _kb/03-api-layer.md for rationale
        boolean hasEmpty = false;
        for (int r = 0; r < R; r++) {
            if ((int) N.get(0, r) <= 0) { hasEmpty = true; break; }
        }
        if (hasEmpty) {
            return solveWithEmptyClasses(L, N, Z);
        }

        // see _kb/03-api-layer.md for rationale
        if (R == 1) {
            return singleClass(L, N, Z);
        }

        // see _kb/03-api-layer.md for rationale
        for (int r = 0; r < R; r++) {
            double zc = 0.0;
            for (int i = 0; i < Z.getNumRows(); i++) zc += Z.get(i, r);
            if (zc <= 0.0) {
                throw new IllegalArgumentException(
                        "pfqn_mom: multiclass (R=" + R + ") requires strictly positive think time Z"
                                + " for every class (class " + (r + 1) + " has Z=0); use pfqn_ca instead.");
            }
        }

        try {

            BigFraction[][] Lf = new BigFraction[M][R];
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < R; j++) {
                    Lf[i][j] = new BigFraction(L.get(i, j));
                }
            }

            BigFraction[] Nf = new BigFraction[R];
            for (int r = 0; r < R; r++) {
                Nf[r] = new BigFraction((int) N.get(0, r));
            }

            BigFraction[] Zf = new BigFraction[R];
            if (Z.getNumRows() == 1) {
                for (int r = 0; r < R; r++) {
                    Zf[r] = new BigFraction(Z.get(0, r));
                }
            } else {
                for (int r = 0; r < R; r++) {
                    BigFraction sum = BigFraction.ZERO;
                    for (int i = 0; i < Z.getNumRows(); i++) {
                        sum = sum.add(new BigFraction(Z.get(i, r)));
                    }
                    Zf[r] = sum;
                }
            }

            int[] n = new int[R];
            BigFraction[] g = null;
            BigFraction[] gr = null;

            FieldMatrix<BigFraction> C = null;
            FieldMatrix<BigFraction> Cg = null;
            FieldMatrix<BigFraction> D = null;
            FieldMatrix<BigFraction> Dr = null;
            FieldMatrix<BigFraction> CgDr = null;

            for (int r = 1; r <= R; r++) {
                Quadruple q = setupls(Lf, Nf, Zf, r);
                C = q.first;
                Cg = q.second;
                D = q.third;
                Dr = q.fourth;

                if (r == 1) {
                    g = new BigFraction[M + 1];
                    for (int i = 0; i < M + 1; i++) g[i] = BigFraction.ONE;
                } else {
                    int actualSize = Maths.nchoosek(M + r - 2, r - 1) * r;
                    BigFraction[] G = new BigFraction[actualSize];
                    for (int i = 0; i < actualSize; i++) G[i] = BigFraction.ZERO;

                    for (int i = 0; i < Maths.nchoosek(M + r - 2, r - 1); i++) {
                        for (int s = 0; s < r - 1; s++) {
                            int srcIdx = i * (r - 1) + s;
                            int dstIdx = i * r + s;
                            G[dstIdx] = g[srcIdx];
                        }
                        int srcIdx = i * (r - 1);
                        int dstIdx = i * r + r - 1;
                        G[dstIdx] = gr[srcIdx];
                    }

                    BigFraction[] Gk = new BigFraction[Maths.nchoosek(M + r - 1, r) * r];
                    for (int i = 0; i < Gk.length; i++) Gk[i] = BigFraction.ZERO;
                    BigFraction[] newGk = blocksolve(M, r, C,
                            matrixVectorMultiply(Cg.scalarMultiply(BigFraction.MINUS_ONE), G), Gk);
                    g = concatenate(newGk, G);
                }

                CgDr = C.createMatrix(Cg.getRowDimension(), Dr.getColumnDimension());
                for (int i = 0; i < Cg.getRowDimension(); i++) {
                    for (int j = 0; j < Dr.getColumnDimension(); j++) {
                        BigFraction sum = BigFraction.ZERO;
                        for (int k = 0; k < Cg.getColumnDimension(); k++) {
                            sum = sum.add(Cg.getEntry(i, k).multiply(Dr.getEntry(k, j)));
                        }
                        CgDr.setEntry(i, j, sum);
                    }
                }

                BigFraction[] Gk = new BigFraction[Maths.nchoosek(M + r - 1, r) * r];
                for (int i = 0; i < Gk.length; i++) Gk[i] = BigFraction.ZERO;
                int numIterations = (int) N.get(0, r - 1) - 1;
                for (int nr = 0; nr < numIterations; nr++) {
                    n[r - 1] = nr + 1;

                    BigFraction nrFrac = new BigFraction(n[r - 1]);
                    BigFraction[] gLocal = g;
                    BigFraction[] G = new BigFraction[Dr.getRowDimension()];
                    for (int i = 0; i < Dr.getRowDimension(); i++) {
                        BigFraction sum = BigFraction.ZERO;
                        for (int j = 0; j < Dr.getColumnDimension(); j++) {
                            if (j < gLocal.length) {
                                sum = sum.add(Dr.getEntry(i, j).multiply(gLocal[j]));
                            }
                        }
                        G[i] = sum.divide(nrFrac);
                    }

                    BigFraction[] b = matrixVectorMultiply(D.subtract(CgDr.scalarMultiply(nrFrac.reciprocal())), g);
                    Gk = blocksolve(M, r, C, b, Gk);
                    g = concatenate(Gk, G);
                }

                gr = g;
                n[r - 1] = (int) N.get(0, r - 1);

                BigFraction nrFrac = new BigFraction(n[r - 1]);
                BigFraction[] gLocal2 = g;
                BigFraction[] G = new BigFraction[Dr.getRowDimension()];
                for (int i = 0; i < Dr.getRowDimension(); i++) {
                    BigFraction sum = BigFraction.ZERO;
                    for (int j = 0; j < Dr.getColumnDimension(); j++) {
                        if (j < gLocal2.length) {
                            sum = sum.add(Dr.getEntry(i, j).multiply(gLocal2[j]));
                        }
                    }
                    G[i] = sum.divide(nrFrac);
                }

                BigFraction[] b = matrixVectorMultiply(D.subtract(CgDr.scalarMultiply(nrFrac.reciprocal())), g);
                Gk = blocksolve(M, r, C, b, Gk);

                FieldMatrix<BigFraction> Cinv = invertMatrix(C);
                FieldMatrix<BigFraction> F1 = concatenateVertical(Cinv.multiply(D), zeros(Dr.getRowDimension(), D.getColumnDimension()));
                FieldMatrix<BigFraction> F2 = concatenateVertical(Cinv.scalarMultiply(BigFraction.MINUS_ONE).multiply(CgDr), Dr);
                FieldMatrix<BigFraction> F = F1.add(F2.scalarMultiply(nrFrac.reciprocal()));

                g = concatenate(Gk, G);
            }

            double[][] Xdata = new double[1][R];
            double[][] Qdata = new double[M][R];
            BigFraction Gconst;

            if (R == 1) {
                Gconst = g[1];
                Xdata[0][0] = gr[1].divide(Gconst).doubleValue();
                Qdata[0][0] = g[0].divide(Gconst).subtract(BigFraction.ONE).doubleValue();
            } else {
                BigFraction[] Gk = new BigFraction[Maths.nchoosek(M + R - 2, R - 1) * (R + 1)];
                for (int i = 0; i < Gk.length; i++) Gk[i] = BigFraction.ZERO;

                for (int i = 0; i < Maths.nchoosek(M + R - 2, R - 1); i++) {
                    for (int s = 0; s < R; s++) {
                        int idx_from = Maths.nchoosek(M + R - 1, R) * R + i * R + s;
                        int idx_to = i * (R + 1) + s;
                        Gk[idx_to] = g[idx_from];
                    }
                    int idx_from_gr = Maths.nchoosek(M + R - 1, R) * R + i * R;
                    int idx_to_gr = i * (R + 1) + R;
                    Gk[idx_to_gr] = gr[idx_from_gr];
                }

                BigFraction[] finalG = null;

                for (int l = R - 1; l >= 1; l--) {
                    BigFraction[] G = new BigFraction[Maths.nchoosek(M + l - 2, l - 1) * (R + 1)];
                    for (int i = 0; i < G.length; i++) G[i] = BigFraction.ZERO;
                    List<int[]> Ik = Maths.sortByNnzPos(Maths.multichooseList(M, l));
                    List<int[]> I = Maths.sortByNnzPos(Maths.multichooseList(M, l - 1));

                    for (int i = 0; i < I.size(); i++) {
                        G[i * (R + 1)] = BigFraction.ZERO;
                        int[] Ii = I.get(i).clone();
                        Ii[0]++;
                        int t = Maths.matchRow(Ik, Ii);
                        G[i * (R + 1)] = G[i * (R + 1)].add(Gk[t * (R + 1)]);

                        for (int s = 0; s < R; s++) {
                            G[i * (R + 1)] = G[i * (R + 1)].subtract(Lf[0][s].multiply(Gk[t * (R + 1) + s + 1]));
                        }

                        for (int s = 0; s < R; s++) {
                            if (Z.get(0, s) == 0.0) {
                                G[i * (R + 1) + s + 1] = BigFraction.ZERO;
                            } else {
                                G[i * (R + 1) + s + 1] = Nf[s].divide(Zf[s]).multiply(G[i * (R + 1)]);

                                for (int j = 0; j < M; j++) {
                                    int[] Ij = I.get(i).clone();
                                    Ij[j]++;
                                    int tj = Maths.matchRow(Ik, Ij);
                                    G[i * (R + 1) + s + 1] = G[i * (R + 1) + s + 1].subtract(
                                            new BigFraction(1 + I.get(i)[j]).multiply(Lf[j][s]).divide(Zf[s])
                                                    .multiply(Gk[tj * (R + 1) + s + 1])
                                    );
                                }
                            }
                        }
                    }

                    if (l > 1) {
                        Gk = G;
                    } else {
                        finalG = G;
                    }
                }

                BigFraction[] G = finalG;
                Gconst = G[0];

                for (int s = 0; s < R; s++) {
                    Xdata[0][s] = G[s + 1].divide(Gconst).doubleValue();
                    for (int m = 0; m < M; m++) {
                        Qdata[m][s] = Lf[m][s].multiply(Gk[m * (R + 1) + s + 1]).divide(Gconst).doubleValue();
                    }
                }
            }

            Matrix X = new Matrix(Xdata);
            Matrix Q = new Matrix(Qdata);

            double lG;
            if (Gconst.getDenominator().equals(BigInteger.ONE)) {
                lG = logBigInteger(Gconst.getNumerator());
            } else {
                lG = logBigInteger(Gconst.getNumerator()) - logBigInteger(Gconst.getDenominator());
            }

            return new Ret.pfqnMom(X, Q, Gconst, lG, g, gr);
        } catch (org.apache.commons.math3.linear.SingularMatrixException e) {
            // see _kb/03-api-layer.md for rationale
            boolean zeroThink = false;
            for (int r = 0; r < Z.getNumCols(); r++) {
                double zc = 0.0;
                for (int i = 0; i < Z.getNumRows(); i++) zc += Z.get(i, r);
                if (zc == 0.0) { zeroThink = true; break; }
            }
            String hint = zeroThink
                    ? " (a class has zero think time Z; pfqn_mom requires Z > 0 for R >= 3 classes)"
                    : "";
            throw new IllegalArgumentException(
                    "pfqn_mom: the method-of-moments linear system is singular" + hint
                            + "; use pfqn_ca for an exact normalizing constant.", e);
        } catch (Exception e) {
            // The multichain moment recursion has data-dependent index/dimension
            // limitations for higher R and certain demand patterns. It never
            // returns a silently-wrong value (verified: whenever it returns it
            // matches pfqn_ca), so surface a clear diagnostic instead of an
            // opaque internal error.
            throw new IllegalArgumentException(
                    "pfqn_mom: the method-of-moments recursion failed for M=" + M + ", R=" + R
                            + " with these demands; use pfqn_ca for an exact normalizing constant.", e);
        }
    }

    private static class Quadruple {
        final FieldMatrix<BigFraction> first;
        final FieldMatrix<BigFraction> second;
        final FieldMatrix<BigFraction> third;
        final FieldMatrix<BigFraction> fourth;
        Quadruple(FieldMatrix<BigFraction> a, FieldMatrix<BigFraction> b,
                  FieldMatrix<BigFraction> c, FieldMatrix<BigFraction> d) {
            this.first = a;
            this.second = b;
            this.third = c;
            this.fourth = d;
        }
    }

    private static Quadruple setupls(BigFraction[][] L, BigFraction[] N, BigFraction[] Z, int R) {
        int M = L.length;
        int[] m = new int[M];
        for (int i = 0; i < M; i++) m[i] = 1;

        List<int[]> Ik = Maths.sortByNnzPos(Maths.multichooseList(M, R));
        List<int[]> I = Maths.sortByNnzPos(Maths.multichooseList(M, R - 1));

        int rows = Ik.size() * R;
        int colsC = Ik.size() * R;
        int colsCg = I.size() * R;
        int colsD = (Ik.size() + I.size()) * R;
        int rowsDr = I.size() * R;

        BigFractionField field = BigFractionField.getInstance();
        Array2DRowFieldMatrix<BigFraction> C = new Array2DRowFieldMatrix<BigFraction>(field, rows, colsC);
        Array2DRowFieldMatrix<BigFraction> Cg = new Array2DRowFieldMatrix<BigFraction>(field, rows, colsCg);
        Array2DRowFieldMatrix<BigFraction> D = new Array2DRowFieldMatrix<BigFraction>(field, rows, colsD);
        Array2DRowFieldMatrix<BigFraction> Dr = new Array2DRowFieldMatrix<BigFraction>(field, rowsDr, colsD);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < colsC; j++) C.setEntry(i, j, BigFraction.ZERO);
            for (int j = 0; j < colsCg; j++) Cg.setEntry(i, j, BigFraction.ZERO);
            for (int j = 0; j < colsD; j++) D.setEntry(i, j, BigFraction.ZERO);
        }
        for (int i = 0; i < rowsDr; i++) {
            for (int j = 0; j < colsD; j++) Dr.setEntry(i, j, BigFraction.ZERO);
        }

        List<Integer> pcpos = new ArrayList<Integer>();
        int currentRow = 0;

        for (int i = 0; i < Ik.size(); i++) {
            int h = 0;
            for (int v : Ik.get(i)) if (v > 0) h++;

            for (int j = 0; j < M; j++) {
                if (Ik.get(i)[j] > 0) {
                    C.setEntry(currentRow, i * R, BigFraction.ONE);
                    for (int s = 0; s < R - 1; s++) {
                        C.setEntry(currentRow, i * R + s + 1, L[j][s].negate());
                    }

                    int[] IkCopy = Ik.get(i).clone();
                    IkCopy[j]--;
                    int idx = Maths.matchRow(I, IkCopy);
                    Cg.setEntry(currentRow, idx * R, BigFraction.MINUS_ONE);

                    D.setEntry(currentRow, i * R, L[j][R - 1]);
                    currentRow++;
                }
            }

            for (int j = 0; j < R - h; j++) {
                pcpos.add(currentRow);
                currentRow++;
            }
        }

        int last = 0;
        for (int i = 0; i < I.size(); i++) {
            for (int s = 0; s < R; s++) {
                Dr.setEntry(i * R + s, Ik.size() * R + i * R + s, Z[R - 1]);
            }

            for (int j = 0; j < M; j++) {
                int[] ICopy = I.get(i).clone();
                ICopy[j]++;
                int t = Maths.matchRow(Ik, ICopy);
                Dr.setEntry(i * R, t * R, new BigFraction(m[j] + I.get(i)[j]).multiply(L[j][R - 1]));

                for (int s = 0; s < R - 1; s++) {
                    C.setEntry(pcpos.get(last + s), t * R + s + 1,
                            new BigFraction(m[j] + I.get(i)[j]).multiply(L[j][s]).negate());
                    Cg.setEntry(pcpos.get(last + s), i * R, N[s]);
                    Cg.setEntry(pcpos.get(last + s), i * R + s + 1, Z[s].negate());
                    Dr.setEntry(i * R + s + 1, t * R + s + 1, new BigFraction(m[j] + I.get(i)[j]).multiply(L[j][R - 1]));
                }
            }
            last += R - 1;
        }

        return new Quadruple(C, Cg, D, Dr);
    }

    private static BigFraction[] blocksolve(int M, int R, FieldMatrix<BigFraction> C,
                                              BigFraction[] b, BigFraction[] Gr) {
        int blockend = C.getColumnDimension() - 1;
        int H = Math.min(M, R);
        BigFraction[] x = new BigFraction[b.length];
        for (int i = 0; i < x.length; i++) x[i] = BigFraction.ZERO;
        BigFraction[] bAdjusted = b.clone();

        for (int i = 0; i < bAdjusted.length; i++) {
            for (int j = 0; j < Gr.length; j++) {
                bAdjusted[i] = bAdjusted[i].subtract(C.getEntry(i, j).multiply(Gr[j]));
            }
        }

        for (int h = H; h >= 1; h--) {
            for (int t = 0; t < Maths.nchoosek(M, h); t++) {
                int blockstart = blockend - Maths.nchoosek(R - 1, R - h) * R + 1;
                int blockSize = blockend - blockstart + 1;

                FieldMatrix<BigFraction> blockC = C.getSubMatrix(blockstart, blockend, blockstart, blockend);
                BigFraction[] sliceArr = new BigFraction[blockend - blockstart + 1];
                for (int k = 0; k < sliceArr.length; k++) sliceArr[k] = bAdjusted[blockstart + k];
                ArrayFieldVector<BigFraction> blockB = new ArrayFieldVector<BigFraction>(BigFractionField.getInstance(), sliceArr);

                if (blockend < C.getColumnDimension() - 1) {
                    for (int i = 0; i < blockSize; i++) {
                        for (int j = blockend + 1; j < C.getColumnDimension(); j++) {
                            BigFraction term = C.getEntry(blockstart + i, j).multiply(x[j]);
                            blockB.setEntry(i, blockB.getEntry(i).subtract(term));
                        }
                    }
                }

                FieldVector<BigFraction> blockX = new FieldLUDecomposition<BigFraction>(blockC).getSolver().solve(blockB);

                for (int i = 0; i < blockSize; i++) {
                    x[blockstart + i] = blockX.getEntry(i);
                }

                blockend = blockstart - 1;
            }
        }

        for (int i = 0; i < x.length; i++) {
            if (i < Gr.length) {
                x[i] = x[i].add(Gr[i]);
            }
        }

        return x;
    }

    private static BigFraction[] matrixVectorMultiply(FieldMatrix<BigFraction> A, BigFraction[] v) {
        BigFraction[] result = new BigFraction[A.getRowDimension()];
        for (int i = 0; i < result.length; i++) result[i] = BigFraction.ZERO;
        for (int i = 0; i < A.getRowDimension(); i++) {
            for (int j = 0; j < A.getColumnDimension(); j++) {
                result[i] = result[i].add(A.getEntry(i, j).multiply(v[j]));
            }
        }
        return result;
    }

    private static FieldMatrix<BigFraction> invertMatrix(FieldMatrix<BigFraction> A) {
        return new FieldLUDecomposition<BigFraction>(A).getSolver().getInverse();
    }

    private static BigFraction[] concatenate(BigFraction[] a, BigFraction[] b) {
        BigFraction[] r = new BigFraction[a.length + b.length];
        System.arraycopy(a, 0, r, 0, a.length);
        System.arraycopy(b, 0, r, a.length, b.length);
        return r;
    }

    private static FieldMatrix<BigFraction> concatenateVertical(FieldMatrix<BigFraction> A, FieldMatrix<BigFraction> B) {
        int rows = A.getRowDimension() + B.getRowDimension();
        int cols = A.getColumnDimension();
        BigFractionField field = BigFractionField.getInstance();
        Array2DRowFieldMatrix<BigFraction> result = new Array2DRowFieldMatrix<BigFraction>(field, rows, cols);

        for (int i = 0; i < A.getRowDimension(); i++) {
            for (int j = 0; j < A.getColumnDimension(); j++) {
                result.setEntry(i, j, A.getEntry(i, j));
            }
        }

        for (int i = 0; i < B.getRowDimension(); i++) {
            for (int j = 0; j < B.getColumnDimension(); j++) {
                result.setEntry(A.getRowDimension() + i, j, B.getEntry(i, j));
            }
        }

        return result;
    }

    private static FieldMatrix<BigFraction> zeros(int rows, int cols) {
        BigFractionField field = BigFractionField.getInstance();
        Array2DRowFieldMatrix<BigFraction> result = new Array2DRowFieldMatrix<BigFraction>(field, rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.setEntry(i, j, BigFraction.ZERO);
            }
        }
        return result;
    }

    private static double logBigInteger(BigInteger bigInt) {
        if (bigInt.signum() <= 0) {
            throw new IllegalArgumentException("Cannot compute log of non-positive number");
        }
        if (bigInt.equals(BigInteger.ONE)) {
            return 0.0;
        }
        int bitLength = bigInt.bitLength();
        if (bitLength <= 53) {
            return Math.log(bigInt.doubleValue());
        }
        BigInteger powerOfTwo = BigInteger.ONE.shiftLeft(bitLength - 1);
        double ratio = bigInt.doubleValue() / powerOfTwo.doubleValue();
        return (bitLength - 1) * Math.log(2.0) + Math.log(ratio);
    }

    /**
     * Exact single-class solution (any number of load-independent stations plus
     * an aggregate delay Z) via a BigFraction Mean Value Analysis recursion.
     * Returns the exact normalizing constant G(N) = prod_{n=1}^{N} 1/X(n),
     * throughput X(N), per-station queue lengths Q(N), and the g/g_1 sequences
     * g[k] = G(k) and g_1 = {G(N-1)}.
     */
    private static Ret.pfqnMom singleClass(Matrix L, Matrix N, Matrix Z) {
        int M = L.getNumRows();
        int Npop = (int) N.get(0, 0);

        BigFraction[] Lf = new BigFraction[M];
        for (int i = 0; i < M; i++) Lf[i] = new BigFraction(L.get(i, 0));

        BigFraction Zf = BigFraction.ZERO;
        for (int i = 0; i < Z.getNumRows(); i++) Zf = Zf.add(new BigFraction(Z.get(i, 0)));

        BigFraction[] Q = new BigFraction[M];
        for (int i = 0; i < M; i++) Q[i] = BigFraction.ZERO;

        BigFraction[] Gseq = new BigFraction[Npop + 1];
        Gseq[0] = BigFraction.ONE;
        BigFraction X = BigFraction.ZERO;

        for (int nn = 1; nn <= Npop; nn++) {
            BigFraction[] Rres = new BigFraction[M];
            BigFraction denom = Zf;
            for (int i = 0; i < M; i++) {
                Rres[i] = Lf[i].multiply(BigFraction.ONE.add(Q[i]));
                denom = denom.add(Rres[i]);
            }
            if (denom.equals(BigFraction.ZERO)) {
                throw new IllegalArgumentException(
                        "pfqn_mom: single-class network has zero total demand and think time.");
            }
            X = new BigFraction(nn).divide(denom);
            for (int i = 0; i < M; i++) Q[i] = X.multiply(Rres[i]);
            Gseq[nn] = Gseq[nn - 1].divide(X); // X(n) = G(n-1)/G(n)
        }

        BigFraction G = Gseq[Npop];
        BigFraction Gm1 = Gseq[Npop - 1];

        Matrix Xout = new Matrix(new double[][]{{X.doubleValue()}});
        double[][] Qdata = new double[M][1];
        for (int i = 0; i < M; i++) Qdata[i][0] = Q[i].doubleValue();
        Matrix Qout = new Matrix(Qdata);

        double lG;
        if (G.getDenominator().equals(BigInteger.ONE)) {
            lG = logBigInteger(G.getNumerator());
        } else {
            lG = logBigInteger(G.getNumerator()) - logBigInteger(G.getDenominator());
        }

        return new Ret.pfqnMom(Xout, Qout, G, lG, Gseq, new BigFraction[]{Gm1});
    }

    /**
     * Solves a model with one or more empty-population classes (N_r == 0) by
     * removing them, solving the reduced problem, and re-expanding the outputs
     * (X_r = Q_r = 0 for each removed class). The normalizing constant is
     * unaffected by classes with zero jobs.
     */
    private static Ret.pfqnMom solveWithEmptyClasses(Matrix L, Matrix N, Matrix Z) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        List<Integer> active = new ArrayList<Integer>();
        for (int r = 0; r < R; r++) {
            if ((int) N.get(0, r) > 0) active.add(r);
        }

        // All classes empty: G = 1, everything zero.
        if (active.isEmpty()) {
            Matrix X = new Matrix(1, R);
            X.zero();
            Matrix Q = new Matrix(M, R);
            Q.zero();
            BigFraction[] g = new BigFraction[]{BigFraction.ONE};
            return new Ret.pfqnMom(X, Q, BigFraction.ONE, 0.0, g, new BigFraction[]{BigFraction.ONE});
        }

        int Ra = active.size();
        Matrix La = new Matrix(M, Ra);
        Matrix Na = new Matrix(1, Ra);
        Matrix Za = new Matrix(Z.getNumRows(), Ra);
        for (int c = 0; c < Ra; c++) {
            int r = active.get(c);
            for (int i = 0; i < M; i++) La.set(i, c, L.get(i, r));
            Na.set(0, c, N.get(0, r));
            for (int i = 0; i < Z.getNumRows(); i++) Za.set(i, c, Z.get(i, r));
        }

        Ret.pfqnMom sub = pfqn_mom(La, Na, Za);

        Matrix X = new Matrix(1, R);
        X.zero();
        Matrix Q = new Matrix(M, R);
        Q.zero();
        for (int c = 0; c < Ra; c++) {
            int r = active.get(c);
            X.set(0, r, sub.X.get(0, c));
            for (int i = 0; i < M; i++) Q.set(i, r, sub.Q.get(i, c));
        }

        return new Ret.pfqnMom(X, Q, sub.G, sub.lG, sub.g, sub.g_1);
    }
}
