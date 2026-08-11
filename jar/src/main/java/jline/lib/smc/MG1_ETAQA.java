/**
 * ETAQA (Efficient Truncation and Aggregation for Queueing Analysis) algorithms
 * for M/G/1-type Markov chains.
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc;

import jline.util.matrix.Matrix;

public final class MG1_ETAQA {
    private MG1_ETAQA() {}

    /** Result of ETAQA computation containing aggregated probabilities. */
    public static final class ETAQAResult {
        public final Matrix pi;
        public final Matrix G;

        public ETAQAResult(Matrix pi, Matrix G) {
            this.pi = pi;
            this.G = G;
        }

        public Matrix getPi() { return pi; }
        public Matrix getG() { return G; }
    }

    /**
     * Computes the aggregated stationary probability vector for an M/G/1-type
     * Markov chain using the ETAQA method.
     */
    public static Matrix mg1_pi_etaqa(Matrix B, Matrix A, Matrix G, Matrix C0) {
        int m = A.getNumRows();
        int dega = A.getNumCols() / m - 1;

        Matrix boundary;
        int mb;
        int degb;

        if (B == null || B.isEmpty()) {
            mb = m;
            degb = dega;
            boundary = A.copy();
        } else {
            mb = B.getNumRows();
            degb = (B.getNumCols() - mb) / m;
            boundary = B.copy();
        }

        Matrix boundaryC0 = (C0 != null) ? C0 : A.extractCols(0, m);

        Matrix rowSums = boundary.sumRows();
        boolean isContinuous = true;
        for (int i = 0; i < mb; i++) {
            if (Math.abs(rowSums.get(0, i)) > 1e-12) {
                isContinuous = false;
                break;
            }
        }

        Matrix Bwork = boundary.copy();
        Matrix Awork = A.copy();
        if (!isContinuous) {
            Matrix discreteCheck = boundary.sumRows();
            boolean isDiscrete = true;
            for (int i = 0; i < mb; i++) {
                if (Math.abs(discreteCheck.get(0, i) - 1.0) > 1e-12) {
                    isDiscrete = false;
                    break;
                }
            }
            if (isDiscrete) {
                for (int i = 0; i < mb; i++) {
                    Bwork.set(i, i, Bwork.get(i, i) - 1.0);
                }
                for (int i = 0; i < m; i++) {
                    Awork.set(i, m + i, Awork.get(i, m + i) - 1.0);
                }
            }
        }

        Matrix Gmat = (G != null) ? G : mg1_g_etaqa(Awork);

        Matrix sumA = Awork.extractCols(dega * m, (dega + 1) * m).copy();
        Matrix alpha = sumA.sumCols();
        for (int i = dega - 1; i >= 1; i--) {
            sumA = sumA.add(Awork.extractCols(i * m, (i + 1) * m));
            alpha = alpha.add(sumA.sumCols());
        }
        sumA = sumA.add(Awork.extractCols(0, m));
        Matrix a = Stat.stat(sumA);
        double drift = a.mult(alpha.transpose()).get(0, 0);

        if (drift >= 1.0) {
            throw new IllegalStateException(
                    "The Markov chain characterized by A is not positive recurrent (drift = " + drift + ")");
        }

        Matrix Shat = Bwork.extractCols(mb + (degb - 1) * m, mb + degb * m).copy();
        if (degb > 1) {
            for (int i = degb - 1; i >= 1; i--) {
                Matrix temp = Bwork.extractCols(mb + (i - 1) * m, mb + i * m).add(
                        Shat.extractCols(0, m).mult(Gmat));
                Shat = temp.concatCols(Shat);
            }
        }

        Matrix S = Awork.extractCols(dega * m, (dega + 1) * m).copy();
        if (dega > 1) {
            for (int i = dega - 1; i >= 1; i--) {
                Matrix temp = Awork.extractCols(i * m, (i + 1) * m).add(
                        S.extractCols(0, m).mult(Gmat));
                S = temp.concatCols(S);
            }
        }

        Matrix firstc = Matrix.ones(mb + 2 * m, 1);

        Matrix secondc = new Matrix(mb + 2 * m, mb);
        for (int i = 0; i < mb; i++) {
            for (int j = 0; j < mb; j++) {
                secondc.set(i, j, Bwork.get(i, j));
            }
        }
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < mb; j++) {
                secondc.set(mb + i, j, boundaryC0.get(i, j));
            }
        }

        if (Shat.getNumCols() < 2 * m) {
            Shat = Shat.concatCols(new Matrix(mb, m));
        }
        if (S.getNumCols() < 2 * m) {
            S = S.concatCols(new Matrix(m, m));
        }

        Matrix thirdc = new Matrix(mb + 2 * m, m);
        Matrix B1 = Bwork.extractCols(mb, mb + m);
        Matrix Shat2G = Shat.extractCols(m, 2 * m).mult(Gmat);
        Matrix B1plusShat2G = B1.add(Shat2G);
        for (int i = 0; i < mb; i++) {
            for (int j = 0; j < m; j++) {
                thirdc.set(i, j, B1plusShat2G.get(i, j));
            }
        }
        Matrix A1 = Awork.extractCols(m, 2 * m);
        Matrix S1G = S.extractCols(m, 2 * m).mult(Gmat);
        Matrix A1plusS1G = A1.add(S1G);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                thirdc.set(mb + i, j, A1plusS1G.get(i, j));
            }
        }

        Matrix Bsum = new Matrix(mb, m);
        Matrix ShatSum = new Matrix(mb, m);
        if (degb > 2) {
            for (int i = 2; i < degb; i++) {
                Bsum = Bsum.add(Bwork.extractCols(mb + (i - 1) * m, mb + i * m));
                if (i < Shat.getNumCols() / m) {
                    ShatSum = ShatSum.add(Shat.extractCols(i * m, (i + 1) * m));
                }
            }
            Bsum = Bsum.add(Bwork.extractCols(mb + (degb - 1) * m, mb + degb * m));
        } else if (degb == 2) {
            Bsum = Bsum.add(Bwork.extractCols(mb + m, mb + 2 * m));
        }

        Matrix Asum = new Matrix(m, m);
        Matrix Ssum = new Matrix(m, m);
        if (dega >= 3) {
            for (int i = 2; i < dega; i++) {
                Ssum = Ssum.add(S.extractCols(i * m, (i + 1) * m));
                Asum = Asum.add(Awork.extractCols(i * m, (i + 1) * m));
            }
            Asum = Asum.add(Awork.extractCols(dega * m, (dega + 1) * m));
        } else if (dega == 2) {
            Asum = Asum.add(Awork.extractCols(dega * m, (dega + 1) * m));
        }

        Matrix fourthc = new Matrix(mb + 2 * m, m);
        Matrix BsumPlusShatSumG = Bsum.add(ShatSum.mult(Gmat));
        for (int i = 0; i < mb; i++) {
            for (int j = 0; j < m; j++) {
                fourthc.set(i, j, BsumPlusShatSumG.get(i, j));
            }
        }
        Matrix AsumPlusSsumG = Asum.add(Ssum.mult(Gmat));
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                fourthc.set(mb + i, j, AsumPlusSsumG.get(i, j));
            }
        }
        Matrix S1 = S.extractCols(m, 2 * m);
        Matrix AsumPlusA1PlusSsumPlusS1G = Asum.add(A1).add(Ssum.add(S1).mult(Gmat));
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                fourthc.set(mb + m + i, j, AsumPlusA1PlusSsumPlusS1G.get(i, j));
            }
        }

        Matrix Xtemp = secondc.concatCols(thirdc).concatCols(fourthc);
        int rankXtemp = Xtemp.rank();

        int redundantCol = 0;
        for (int i = 0; i < Xtemp.getNumCols(); i++) {
            Matrix withoutCol = removeColumn(Xtemp, i);
            if (withoutCol.rank() == rankXtemp) {
                redundantCol = i;
                break;
            }
        }
        Xtemp = removeColumn(Xtemp, redundantCol);

        Matrix Xnew = firstc.concatCols(Xtemp);
        Matrix rside = new Matrix(1, mb + 2 * m);
        rside.set(0, 0, 1.0);

        Matrix pi = rside.mult(Xnew.pinv());
        return pi;
    }

    public static Matrix mg1_pi_etaqa(Matrix B, Matrix A) {
        return mg1_pi_etaqa(B, A, null, null);
    }

    public static Matrix mg1_pi_etaqa(Matrix B, Matrix A, Matrix G) {
        return mg1_pi_etaqa(B, A, G, null);
    }

    /**
     * Computes the G matrix for M/G/1-type Markov chains.
     */
    public static Matrix mg1_g_etaqa(Matrix A) {
        int r = A.getNumRows();

        Matrix rowSums = A.sumRows();
        boolean isContinuous = true;
        for (int i = 0; i < r; i++) {
            if (Math.abs(rowSums.get(0, i)) > 1e-16) {
                isContinuous = false;
                break;
            }
        }

        Matrix Awork;
        if (isContinuous) {
            Matrix A1 = A.extractCols(r, 2 * r);
            double minDiag = Double.MAX_VALUE;
            for (int i = 0; i < r; i++) {
                if (A1.get(i, i) < minDiag) {
                    minDiag = A1.get(i, i);
                }
            }
            if (minDiag >= 0) {
                throw new IllegalArgumentException(
                        "Invalid generator: A1 diagonal must be negative for continuous time");
            }
            double lamb = -minDiag;
            Awork = A.scale(1.0 / lamb);
            for (int i = 0; i < r; i++) {
                Awork.set(i, r + i, Awork.get(i, r + i) + 1.0);
            }
        } else {
            Awork = A.copy();
        }

        return mg1_cr_internal(Awork);
    }

    private static Matrix mg1_cr_internal(Matrix A) {
        int m = A.getNumRows();
        int n = A.getNumCols() / m;

        if (n < 2) {
            throw new IllegalArgumentException("Need at least 2 blocks");
        }

        int maxIter = 100;
        double tol = 1e-14;

        Matrix G = new Matrix(m, m);
        for (int iter = 0; iter < maxIter; iter++) {
            Matrix Gnew = A.extractCols(0, m).copy();
            Matrix Gpow = G.copy();
            for (int i = 1; i < n; i++) {
                Gnew = Gnew.add(A.extractCols(i * m, (i + 1) * m).mult(Gpow));
                Gpow = Gpow.mult(G);
            }
            double diff = Gnew.sub(G).normFrobenius();
            if (diff < tol) {
                return Gnew;
            }
            G = Gnew;
        }
        return G;
    }

    /**
     * Computes the n-th moment of queue length using ETAQA.
     */
    public static double mg1_qlen_etaqa(Matrix B, Matrix A, Matrix pi, int n) {
        int m = A.getNumRows();
        int dega = A.getNumCols() / m - 1;

        int mb;
        int degb;
        if (B == null || B.isEmpty()) {
            mb = m;
            degb = dega;
        } else {
            mb = B.getNumRows();
            degb = (B.getNumCols() - mb) / m;
        }

        Matrix pi1 = pi.extractCols(mb, mb + m);
        Matrix piStar = pi.extractCols(mb + m, mb + 2 * m);

        if (n == 1) {
            Matrix e = Matrix.ones(m, 1);
            return pi1.mult(e).get(0, 0) + 2.0 * piStar.mult(e).get(0, 0);
        }

        double moment = 0.0;
        Matrix e = Matrix.ones(m, 1);
        moment += pi1.mult(e).get(0, 0);

        double piStarSum = piStar.mult(e).get(0, 0);
        for (int k = 2; k < 100; k++) {
            double levelProb = piStarSum * Math.pow(1.0 - piStarSum, k - 2);
            moment += Math.pow((double) k, n) * levelProb;
            if (levelProb < 1e-15) break;
        }
        return moment;
    }

    private static Matrix removeColumn(Matrix M, int col) {
        if (M.getNumCols() <= 1) return M;
        Matrix result = new Matrix(M.getNumRows(), M.getNumCols() - 1);
        int destCol = 0;
        for (int j = 0; j < M.getNumCols(); j++) {
            if (j != col) {
                for (int i = 0; i < M.getNumRows(); i++) {
                    result.set(i, destCol, M.get(i, j));
                }
                destCol++;
            }
        }
        return result;
    }
}
