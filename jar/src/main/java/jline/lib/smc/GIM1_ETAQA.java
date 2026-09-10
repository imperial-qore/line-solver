/**
 * @file GI/M/1-type ETAQA
 *
 * The aggregated stationary vector and the queue-length moments of a
 * GI/M/1-type Markov chain. Port of MAMSolver's GIM1_R_ETAQA.m,
 * GIM1_pi_ETAQA.m and GIM1_qlen_ETAQA.m (matlab/lib/thirdparty/MAMSolver).
 *
 * ETAQA solves a FINITE system for exactly three aggregates -- pi_0, pi_1 and
 * pi* = sum_{j>=2} pi_j -- by replacing the infinitely many balance equations
 * for levels 2 and above by their sum. It is exact: no level is truncated and
 * no tail is fitted, which is what separates it from a level-truncated solve
 * or from reading the moments off a matrix-geometric tail.
 *
 * REPRODUCED REFERENCE DEFECT: GIM1_qlen_ETAQA.m initializes its accumulator
 * with A(3), a SCALAR read at column-major linear index 3 of the stacked A
 * where the third BLOCK is meant, and MATLAB broadcasts it over the m x m
 * accumulator. It corrupts the last column of the moment system whenever
 * m > 1 and can make the reported mean queue length NEGATIVE. Value and
 * broadcast are reproduced here, including the case where the loop that would
 * turn the scalar into a matrix never runs, so that this port agrees with
 * MATLAB rather than silently answering something else.
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc;

import jline.util.matrix.Matrix;

public final class GIM1_ETAQA {
    private GIM1_ETAQA() {}

    /**
     * R of a GI/M/1-type chain, uniformized first. Port of GIM1_R_ETAQA.m.
     *
     * `A` is the VERTICAL stack [A0; A1; ...; Amax], which is how a GI/M/1-type
     * sequence is written; the reference transposes that stack into the
     * horizontal one GIM1_R wants and asks for the automatic dual with
     * functional iterations.
     */
    public static Matrix gim1_r_etaqa(Matrix A) {
        int s = A.getNumCols();
        if (A.getNumRows() % s != 0) {
            throw new IllegalArgumentException(
                    "GIM1_R_ETAQA: A is not a vertical stack of square blocks");
        }
        int b = A.getNumRows() / s;
        Matrix An = A.copy();
        double t = Double.MAX_VALUE;
        for (int i = 0; i < s; i++) {
            if (An.get(s + i, i) < t) {
                t = An.get(s + i, i);
            }
        }
        if (!(t > 0)) {
            An = An.scale(1.0 / (-t));
            for (int i = 0; i < s; i++) {
                An.set(s + i, i, An.get(s + i, i) + 1.0);
            }
        }
        Matrix wide = new Matrix(s, s * b);
        for (int k = 0; k < b; k++) {
            for (int i = 0; i < s; i++) {
                for (int j = 0; j < s; j++) {
                    wide.set(i, k * s + j, An.get(k * s + i, j));
                }
            }
        }
        return GIM1_R.gim1_R_dual(wide, "A", "FI");
    }

    /**
     * Aggregated stationary vector [pi0, pi1, pi2+pi3+...] of a GI/M/1-type
     * chain. Port of GIM1_pi_ETAQA.m. B and A are vertical stacks; B0 is the
     * reference's 'Boundary' option (null for the default A0).
     */
    public static Matrix gim1_pi_etaqa(Matrix Bin, Matrix Ain, Matrix R, Matrix B0in) {
        Matrix B = Bin.copy();
        Matrix A = Ain.copy();
        int m = R.getNumRows();
        int mb = B.getNumCols();
        if ((B.getNumRows() - mb) % m != 0) {
            throw new IllegalArgumentException(
                    "GIM1_pi_ETAQA: input matrix B has an incorrect number of rows");
        }
        int degb = (B.getNumRows() - mb) / m;
        if (A.getNumRows() % m != 0) {
            throw new IllegalArgumentException(
                    "GIM1_pi_ETAQA: input matrix A has an incorrect number of rows");
        }
        int dega = A.getNumRows() / m - 1;

        checkSpectralRadius(R, "GIM1_pi_ETAQA");

        Matrix B0 = (B0in != null) ? B0in : A.extractRows(0, m);
        Matrix Btop = B.extractRows(0, mb);

        // A transition matrix is turned into a generator, as the reference tests it.
        Matrix test = Btop.add(B0);
        double tot = 0.0;
        for (int i = 0; i < mb; i++) {
            double rs = 0.0;
            for (int j = 0; j < test.getNumCols(); j++) {
                rs += test.get(i, j);
            }
            tot += rs - 1.0;
        }
        if (Math.abs(tot) < 1e-10) {
            for (int i = 0; i < mb; i++) {
                B.set(i, i, B.get(i, i) - 1.0);
            }
            for (int i = 0; i < m; i++) {
                A.set(m + i, i, A.get(m + i, i) - 1.0);
            }
            Btop = B.extractRows(0, mb);
        }

        Matrix firstc = Matrix.ones(mb + 2 * m, 1);

        // sum_{i>=2} R^(i-2) (I - R) B(i)
        Matrix temp = Matrix.eye(m).sub(R);
        Matrix tempsum = new Matrix(m, mb);
        for (int i = 2; i <= degb; i++) {
            tempsum = tempsum.add(temp.mult(vblock(B, mb, m, i)));
            temp = R.mult(temp);
        }
        Matrix secondc = stack(Btop, vblock(B, mb, m, 1), tempsum, mb + 2 * m, mb);

        // sum_{i>=2} R^(i-2) (I - R) A(i)
        temp = Matrix.eye(m).sub(R);
        tempsum = new Matrix(m, m);
        for (int i = 2; i <= dega; i++) {
            tempsum = tempsum.add(temp.mult(A.extractRows(i * m, (i + 1) * m)));
            temp = R.mult(temp);
        }
        Matrix thirdc = stack(B0, A.extractRows(m, 2 * m), tempsum, mb + 2 * m, m);

        // sum_{i>=2} R^(i-1) A(i)
        temp = R.copy();
        tempsum = new Matrix(m, m);
        for (int i = 2; i <= dega; i++) {
            tempsum = tempsum.add(temp.mult(A.extractRows(i * m, (i + 1) * m)));
            temp = R.mult(temp);
        }
        Matrix A0 = A.extractRows(0, m);
        Matrix A1 = A.extractRows(m, 2 * m);
        Matrix fourthFull = stack(new Matrix(mb, m), A0, A0.add(A1).add(tempsum), mb + 2 * m, m);
        Matrix fourthc = fourthFull.extractCols(0, m - 1);

        Matrix X = firstc.concatCols(secondc).concatCols(thirdc).concatCols(fourthc);
        Matrix rside = new Matrix(1, mb + 2 * m);
        rside.set(0, 0, 1.0);
        return rside.mult(X.inv());
    }

    /**
     * n-th moment of the level of a GI/M/1-type chain from the ETAQA
     * aggregates. Port of GIM1_qlen_ETAQA.m, including the scalar A(3) defect
     * described in the class comment.
     */
    public static double gim1_qlen_etaqa(Matrix Bin, Matrix Ain, Matrix R, Matrix pi, int n,
                                         Matrix B0in) {
        Matrix B = Bin.copy();
        Matrix A = Ain.copy();
        int m = R.getNumRows();
        int mb = B.getNumCols();
        if ((B.getNumRows() - mb) % m != 0) {
            throw new IllegalArgumentException(
                    "GIM1_qlen_ETAQA: input matrix B has an incorrect number of rows");
        }
        int degb = (B.getNumRows() - mb) / m;
        if (A.getNumRows() % m != 0) {
            throw new IllegalArgumentException(
                    "GIM1_qlen_ETAQA: input matrix A has an incorrect number of rows");
        }
        int dega = A.getNumRows() / m - 1;

        checkSpectralRadius(R, "GIM1_qlen_ETAQA");

        Matrix B0 = (B0in != null) ? B0in : A.extractRows(0, m);

        Matrix pi0 = pi.extractCols(0, mb);
        Matrix pi1 = pi.extractCols(mb, mb + m);
        Matrix pistar = pi.extractCols(mb + m, mb + 2 * m);

        if (n == 0) {
            return 1.0;
        }

        // The scalar of the reproduced defect, at column-major linear index 3.
        double a3 = A.get(2 % A.getNumRows(), 2 / A.getNumRows());

        Matrix Btop = B.extractRows(0, mb);
        Matrix test = Btop.add(B0);
        double tot = 0.0;
        for (int i = 0; i < mb; i++) {
            double rs = 0.0;
            for (int j = 0; j < test.getNumCols(); j++) {
                rs += test.get(i, j);
            }
            tot += rs - 1.0;
        }
        if (Math.abs(tot) < 1e-10) {
            for (int i = 0; i < mb; i++) {
                B.set(i, i, B.get(i, i) - 1.0);
            }
            for (int i = 0; i < m; i++) {
                A.set(m + i, i, A.get(m + i, i) - 1.0);
            }
        }

        Matrix A0 = A.extractRows(0, m);
        Matrix A1 = A.extractRows(m, 2 * m);

        Matrix rpower = Matrix.eye(m);
        Matrix lsum = A0.add(A1);
        for (int i = 2; i <= dega; i++) {
            lsum = lsum.add(rpower.mult(A.extractRows(i * m, (i + 1) * m)));
            rpower = R.mult(rpower);
        }

        double[] leftr = new double[m];
        if (degb >= 2 && dega >= 2) {
            boolean loopRuns = dega >= 3;
            Matrix part1 = new Matrix(m, m);
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    part1.set(i, j, a3);
                }
            }
            Matrix part2 = new Matrix(m, m);
            rpower = R.copy();
            for (int i = 1; i + 2 <= dega; i++) {
                Matrix blk = rpower.mult(A.extractRows((i + 2) * m, (i + 3) * m));
                part1 = part1.add(blk);
                part2 = part2.add(blk.scale((double) i));
                rpower = R.mult(rpower);
            }
            for (int i = 0; i < m; i++) {
                double s = 0.0;
                if (loopRuns) {
                    for (int j = 0; j < m; j++) {
                        s += part1.get(i, j) + part2.get(i, j);
                    }
                } else {
                    s = a3;  // scalar * ones(m,1), never expanded to a matrix
                }
                for (int j = 0; j < m; j++) {
                    s -= A0.get(i, j);
                }
                leftr[i] = s;
            }
        } else if (degb == 1 && dega != 1) {
            Matrix acc = new Matrix(m, m);
            rpower = Matrix.eye(m);
            for (int i = 2; i <= dega; i++) {
                acc = acc.add(rpower.mult(A.extractRows(i * m, (i + 1) * m)).scale((double) (i - 1)));
                rpower = R.mult(rpower);
            }
            for (int i = 0; i < m; i++) {
                double s = 0.0;
                for (int j = 0; j < m; j++) {
                    s += acc.get(i, j) - A0.get(i, j);
                }
                leftr[i] = s;
            }
        } else {
            throw new IllegalArgumentException(
                    "GIM1_qlen_ETAQA: the number of A blocks is not enough, this is a reducible "
                            + "Markov chain");
        }

        Matrix lsleft = new Matrix(m, m);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j + 1 < m; j++) {
                lsleft.set(i, j, lsum.get(i, j));
            }
            lsleft.set(i, m - 1, leftr[i]);
        }

        Matrix[] r = new Matrix[n + 1];
        r[0] = pistar.copy();
        Matrix[] reuse = new Matrix[n];
        Matrix onesM = Matrix.ones(m, 1);
        Matrix onesMb = Matrix.ones(mb, 1);

        for (int k = 1; k <= n; k++) {
            double dk = (double) k;
            Matrix bk = pi0.mult(B0).scale(Math.pow(2.0, dk))
                    .add(pi1.mult(A1.scale(Math.pow(2.0, dk)).add(A0.scale(Math.pow(3.0, dk)))))
                    .scale(-1.0);
            for (int l = 1; l <= k; l++) {
                Matrix Ml = A1.add(A0.scale(Math.pow(2.0, (double) l)));
                bk = bk.sub(r[k - l].mult(Ml).scale(MG1_ETAQA.bino(k, l)));
            }

            Matrix tempsum = new Matrix(m, m);
            rpower = R.copy();
            for (int i = 1; i + 2 <= dega; i++) {
                double t = 0.0;
                for (int z = 1; z <= i; z++) {
                    t += Math.pow((double) z, dk);
                }
                tempsum = tempsum.add(rpower.mult(A.extractRows((i + 2) * m, (i + 3) * m)).scale(t));
                rpower = rpower.mult(R);
            }
            reuse[k - 1] = A0.sub(tempsum).mult(onesM);

            double ck = Math.pow(2.0, dk) * pi1.mult(A0).mult(onesM).get(0, 0);
            Matrix tempsum2 = new Matrix(m, mb);
            rpower = R.copy();
            for (int i = 2; i <= degb; i++) {
                double t = 0.0;
                for (int z = 2; z <= i; z++) {
                    t += Math.pow((double) z, dk);
                }
                tempsum2 = tempsum2.add(rpower.mult(vblock(B, mb, m, i)).scale(t));
                rpower = R.mult(rpower);
            }
            ck -= pi1.mult(tempsum2).mult(onesMb).get(0, 0);
            for (int l = 1; l <= k; l++) {
                ck += MG1_ETAQA.bino(k, l) * r[k - l].mult(reuse[l - 1]).get(0, 0);
            }

            Matrix rside = new Matrix(1, m);
            for (int j = 0; j + 1 < m; j++) {
                rside.set(0, j, bk.get(0, j));
            }
            rside.set(0, m - 1, ck);
            r[k] = rside.mult(lsleft.inv());
        }

        double qlen = 0.0;
        for (int j = 0; j < m; j++) {
            qlen += r[n].get(0, j) + pi1.get(0, j);
        }
        return qlen;
    }

    /**
     * Block i (i >= 1) of a vertical stack whose boundary occupies the first mb
     * rows: M(mb+(i-1)*m+1 : mb+i*m, :). The blocks are m x mb, which is not
     * m x m when the boundary has a different size.
     */
    private static Matrix vblock(Matrix M, int mb, int m, int i) {
        return M.extractRows(mb + (i - 1) * m, mb + i * m);
    }

    /** Vertical concatenation of three blocks into a (rows x cols) column. */
    private static Matrix stack(Matrix top, Matrix mid, Matrix bot, int rows, int cols) {
        Matrix out = new Matrix(rows, cols);
        int off = 0;
        Matrix[] parts = new Matrix[]{top, mid, bot};
        for (int p = 0; p < parts.length; p++) {
            for (int i = 0; i < parts[p].getNumRows(); i++) {
                for (int j = 0; j < cols && j < parts[p].getNumCols(); j++) {
                    out.set(off + i, j, parts[p].get(i, j));
                }
            }
            off += parts[p].getNumRows();
        }
        return out;
    }

    /** The reference's own test that sp(R) < 1, column by column. */
    private static void checkSpectralRadius(Matrix R, String who) {
        int m = R.getNumRows();
        Matrix temp = Matrix.eye(m).sub(R).inv();
        boolean allCols = true;
        for (int j = 0; j < m && allCols; j++) {
            boolean any = false;
            for (int i = 0; i < m; i++) {
                if (temp.get(i, j) < -100.0 * 2.220446049250313e-16) {
                    any = true;
                }
            }
            allCols = any;
        }
        if (allCols) {
            throw new IllegalStateException(
                    who + ": the spectral radius of R is not below 1, GIM1 is not positive "
                            + "recurrent");
        }
    }
}
