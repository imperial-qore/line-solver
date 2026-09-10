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

        // sumRows() returns a COLUMN of row sums, so the index is (i,0). Reading
        // it as (0,i) threw "Outside of matrix bounds" for every boundary with
        // more than one row, i.e. for every model this solver actually builds.
        Matrix rowSums = boundary.sumRows();
        double total = 0.0;
        for (int i = 0; i < mb; i++) {
            total += rowSums.get(i, 0);
        }
        boolean isContinuous = !(total > 1e-12);

        Matrix Bwork = boundary.copy();
        Matrix Awork = A.copy();
        if (!isContinuous) {
            // The reference tests the TOTAL of the row sums against mb, not each
            // row against 1: MG1_pi_ETAQA.m's (sum(B,2) - ones(mb,1))'*ones(mb,1).
            boolean isDiscrete = (total - (double) mb) < 1e-12;
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

        // alpha is sum(sumA,2), the ROW sums: sumCols() returns the column sums,
        // which is a different vector and made the drift test answer about a
        // quantity the reference never forms.
        Matrix sumA = Awork.extractCols(dega * m, (dega + 1) * m).copy();
        Matrix alpha = sumA.sumRows();
        for (int i = dega - 1; i >= 1; i--) {
            sumA = sumA.add(Awork.extractCols(i * m, (i + 1) * m));
            alpha = alpha.add(sumA.sumRows());
        }
        sumA = sumA.add(Awork.extractCols(0, m));
        Matrix a = Stat.stat(sumA);
        double drift = a.mult(alpha).get(0, 0);

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
     * Computes the G matrix for M/G/1-type Markov chains. Port of
     * MG1_G_ETAQA.m: the generator is uniformized by dividing through by
     * -min(diag(A1)) and adding the identity back onto A1, then handed to
     * cyclic reduction.
     *
     * The uniformization is UNCONDITIONAL because the reference's own
     * discrete-time test is write-only: MG1_G_ETAQA.m assigns a flag named
     * isdicrete and reads one named isdiscrete, so the branch never changes.
     * For the generators LINE passes it that is the correct branch anyway.
     */
    public static Matrix mg1_g_etaqa(Matrix A) {
        int r = A.getNumRows();
        if (A.getNumCols() % r != 0) {
            throw new IllegalArgumentException(
                    "MG1_G_ETAQA: A is not a block sequence of A's width");
        }
        Matrix A1 = A.extractCols(r, 2 * r);
        double t = Double.MAX_VALUE;
        for (int i = 0; i < r; i++) {
            if (A1.get(i, i) < t) {
                t = A1.get(i, i);
            }
        }
        if (t > 0) {
            throw new IllegalArgumentException(
                    "MG1_G_ETAQA: this is not a stochastic matrix, neither continuous nor "
                            + "discrete; every row must sum to 0 or 1");
        }
        Matrix Awork = A.scale(1.0 / (-t));
        for (int i = 0; i < r; i++) {
            Awork.set(i, r + i, Awork.get(i, r + i) + 1.0);
        }
        return MG1_pi.mg1_cr(Awork);
    }

    /**
     * Computes the n-th moment of the level (the queue length) of an
     * M/G/1-type chain from the ETAQA aggregates. Port of MG1_qlen_ETAQA.m.
     *
     * The moment is NOT read off the aggregates directly: the vectors
     * r^(k) = sum_{j>=2} j^k pi_j are propagated through a recurrence whose
     * left-hand side is one fixed m x m system, so the n-th moment costs n
     * solves of that size however heavy the tail is.
     *
     * REPRODUCED REFERENCE DEFECT: F0j is built from block dega down to block
     * 2 but indexed from 1, so F0j(j) is the tail sum starting at j+1 rather
     * than at j.
     */
    public static double mg1_qlen_etaqa(Matrix B, Matrix A, Matrix pi, int n) {
        int m = A.getNumRows();
        int dega = A.getNumCols() / m - 1;

        int mb;
        int degb;
        Matrix Bw;
        if (B == null || B.isEmpty()) {
            mb = m;
            degb = dega;
            Bw = A.copy();
        } else {
            Bw = B.copy();
            mb = Bw.getNumRows();
            if ((Bw.getNumCols() - mb) % m != 0) {
                throw new IllegalArgumentException(
                        "MG1_qlen_ETAQA: matrix B has an incorrect number of columns");
            }
            degb = (Bw.getNumCols() - mb) / m;
        }

        double mass = 0.0;
        for (int j = 0; j < pi.getNumCols(); j++) {
            mass += pi.get(0, j);
        }
        if (Math.abs(mass - 1.0) > 1e-10) {
            throw new IllegalArgumentException(
                    "MG1_qlen_ETAQA: the input probability vector does not sum up to 1");
        }

        Matrix pi0 = pi.extractCols(0, mb);
        Matrix pi1 = pi.extractCols(mb, mb + m);
        Matrix pistar = pi.extractCols(mb + m, mb + 2 * m);

        // lsleft = [(A0+...+Amax) without its last column, (F11 - A0) e]
        Matrix Asum = A.extractCols(0, m).copy();
        for (int i = 1; i <= dega; i++) {
            Asum = Asum.add(A.extractCols(i * m, (i + 1) * m));
        }
        Matrix F11 = A.extractCols(2 * m, 3 * m).copy();
        for (int i = 3; i <= dega; i++) {
            F11 = F11.add(A.extractCols(i * m, (i + 1) * m).scale((double) (i - 1)));
        }
        Matrix A0 = A.extractCols(0, m);
        Matrix lsleft = new Matrix(m, m);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j + 1 < m; j++) {
                lsleft.set(i, j, Asum.get(i, j));
            }
            double s = 0.0;
            for (int j = 0; j < m; j++) {
                s += F11.get(i, j) - A0.get(i, j);
            }
            lsleft.set(i, m - 1, s);
        }

        // Fhat0(j) = sum_{l>=j} B(l), j = 1..degb
        Matrix[] fhat0j = new Matrix[degb];
        fhat0j[degb - 1] = Bw.extractCols(mb + (degb - 1) * m, Bw.getNumCols()).copy();
        for (int j = degb - 1; j >= 1; j--) {
            fhat0j[j - 1] = Bw.extractCols(mb + (j - 1) * m, mb + j * m).add(fhat0j[j]);
        }

        // F0(j) = sum_{l>=j+1} A(l), the reference's off-by-one indexing
        int nf0 = dega - 1;
        Matrix[] f0j = new Matrix[nf0 > 0 ? nf0 : 0];
        if (nf0 > 0) {
            f0j[nf0 - 1] = A.extractCols(dega * m, (dega + 1) * m).copy();
            for (int j = dega - 1; j >= 2; j--) {
                f0j[j - 2] = A.extractCols(j * m, (j + 1) * m).add(f0j[j - 1]);
            }
        }

        Matrix[] r = new Matrix[n + 1];
        r[0] = pistar.copy();

        Matrix[] frestsaver = new Matrix[n];
        Matrix[] fcrestsaver = new Matrix[n];
        for (int l = 1; l <= n; l++) {
            Matrix t1 = new Matrix(m, m);
            for (int j = 2; j <= dega; j++) {
                t1 = t1.add(A.extractCols(j * m, (j + 1) * m).scale(Math.pow((double) j, (double) l)));
            }
            frestsaver[l - 1] = t1;
            Matrix t2 = new Matrix(m, m);
            for (int j = 1; j <= dega; j++) {
                if (j <= f0j.length) {
                    t2 = t2.add(f0j[j - 1].scale(Math.pow((double) j, (double) l)));
                }
            }
            fcrestsaver[l - 1] = t2;
        }

        Matrix A1blk = A.extractCols(m, 2 * m);
        Matrix ones = Matrix.ones(m, 1);
        for (int k = 1; k <= n; k++) {
            double dk = (double) k;
            Matrix fhatkM = new Matrix(mb, m);
            for (int j = 1; j <= degb; j++) {
                fhatkM = fhatkM.add(Bw.extractCols(mb + (j - 1) * m, mb + j * m)
                        .scale(Math.pow((double) (j + 1), dk)));
            }
            Matrix fhatk = pi0.mult(fhatkM);

            Matrix fkM = new Matrix(m, m);
            for (int j = 2; j <= dega; j++) {
                fkM = fkM.add(A.extractCols(j * m, (j + 1) * m)
                        .scale(Math.pow((double) (j + 1), dk)));
            }
            fkM = A1blk.scale(Math.pow(2.0, dk)).add(fkM);
            Matrix fk = pi1.mult(fkM);

            Matrix frest = new Matrix(1, m);
            for (int l = 1; l <= k; l++) {
                Matrix M = A1blk.add(frestsaver[l - 1]);
                Matrix tvec = r[k - l].mult(M);
                frest = frest.add(tvec.scale(bino(k, l)));
            }

            Matrix bk = new Matrix(1, m);
            for (int j = 0; j < m; j++) {
                bk.set(0, j, -fhatk.get(0, j) - fk.get(0, j) - frest.get(0, j));
            }

            Matrix fchatkM = new Matrix(mb, m);
            for (int j = 2; j <= degb; j++) {
                fchatkM = fchatkM.add(fhat0j[j - 1].scale(Math.pow((double) j, dk)));
            }
            double fchatk = pi0.mult(fchatkM).mult(ones).get(0, 0);

            Matrix fckM = new Matrix(m, m);
            for (int j = 1; j + 1 <= dega; j++) {
                if (j <= f0j.length) {
                    fckM = fckM.add(f0j[j - 1].scale(Math.pow((double) (j + 1), dk)));
                }
            }
            double fck = pi1.mult(fckM).mult(ones).get(0, 0);

            double fcrest = 0.0;
            for (int l = 1; l <= k; l++) {
                fcrest += bino(k, l) * r[k - l].mult(fcrestsaver[l - 1]).mult(ones).get(0, 0);
            }

            double ck = -fchatk - fck - fcrest;

            Matrix rside = new Matrix(1, m);
            for (int j = 0; j + 1 < m; j++) {
                rside.set(0, j, bk.get(0, j));
            }
            rside.set(0, m - 1, ck);
            r[k] = rside.mult(lsleft.inv());
        }

        double qlen = 0.0;
        for (int j = 0; j < m; j++) {
            qlen += r[n].get(0, j);
        }
        for (int j = 0; j < m; j++) {
            qlen += pi1.get(0, j);
        }
        return qlen;
    }

    /** The binomial coefficient the reference computes as a ratio of factorials. */
    static double bino(int n, int k) {
        double num = 1.0;
        for (int i = 1; i <= n; i++) {
            num *= i;
        }
        double den = 1.0;
        for (int i = 1; i <= k; i++) {
            den *= i;
        }
        for (int i = 1; i <= n - k; i++) {
            den *= i;
        }
        return num / den;
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
