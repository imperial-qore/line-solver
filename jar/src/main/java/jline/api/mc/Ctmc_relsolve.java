/**
 * @file Equilibrium distribution of a CTMC re-normalized w.r.t. reference state probability
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import java.util.Map;
import java.util.Set;

import jline.util.matrix.Matrix;

public final class Ctmc_relsolve {
    private Ctmc_relsolve() {}

    /**
     * Equilibrium distribution of a continuous-time Markov chain re-normalized with respect to
     * the probability of a reference state.
     *
     * @return Object[] containing [p, Q, nConnComp, connComp]
     */
    public static Object[] ctmc_relsolve(Matrix Q, int refstate, Map<String, Object> options) {
        Matrix Qmat = Q.copy();

        if (Qmat.length() > 6000) {
            System.out.println("ctmc_relsolve: the order of Q is large (" + Qmat.length() + ").");
        }

        if (Qmat.length() == 1) {
            Matrix p = new Matrix(1, 1);
            p.set(0, 0, 1.0);
            return new Object[] { p, Qmat, 1, Matrix.ones(1, 1) };
        }

        Qmat = Ctmc_makeinfgen.ctmc_makeinfgen(Qmat);
        int n = Qmat.length();

        // Check connectivity: B = abs(Q+Q')>0
        Matrix B = Qmat.add(1.0, Qmat.transpose());
        B.absEq();
        for (int colIdx = 0; colIdx < B.getNumCols(); colIdx++) {
            int col1 = B.getColIndexes()[colIdx];
            int col2 = B.getColIndexes()[colIdx + 1];

            for (int i = col1; i < col2; i++) {
                if (B.getNonZeroValues()[i] > 0) B.getNonZeroValues()[i] = 1.0;
            }
        }

        Set<Set<Integer>> sets = Matrix.weaklyConnect(B, null);
        int nConnComp = sets.size();

        Matrix connComp = new Matrix(1, n);
        int compIdx = 0;
        for (Set<Integer> set : sets) {
            for (int nodeIdx : set) {
                connComp.set(0, nodeIdx, (double) compIdx);
            }
            compIdx++;
        }

        if (nConnComp > 1) {
            Matrix p = new Matrix(1, n);

            int c = 0;
            for (Set<Integer> set_c : sets) {
                Matrix Qc = new Matrix(set_c.size(), set_c.size());
                int Qc_row = 0;
                for (int q_row : set_c) {
                    int Qc_col = 0;
                    for (int q_col : set_c) {
                        Qc.set(Qc_row, Qc_col, Qmat.get(q_row, q_col));
                        Qc_col++;
                    }
                    Qc_row++;
                }

                Matrix QcProper = Ctmc_makeinfgen.ctmc_makeinfgen(Qc);
                Matrix pc = Ctmc_solve.ctmc_solve(QcProper);

                int idx = 0;
                for (int i : set_c) {
                    p.set(0, i, pc.get(0, idx));
                    idx++;
                }
                c++;
            }

            double pSum = p.sumRows(0);
            p.divide(pSum, p, true);

            return new Object[] { p, Qmat, nConnComp, connComp };
        }

        if (Qmat.getNonZeroLength() == 0) {
            Matrix p = new Matrix(1, n);
            p.fill(1.0 / n);
            return new Object[] { p, Qmat, nConnComp, connComp };
        }

        Matrix p = new Matrix(1, n);
        Matrix b = new Matrix(n, 1);

        Matrix nnzel = new Matrix(1, n);
        for (int i = 0; i < n; i++) nnzel.set(0, i, (double) i);
        Matrix Qnnz = Qmat.copy();
        Matrix bnnz = b.copy();
        Matrix Qnnz_1 = Qmat.copy();
        Matrix bnnz_1 = bnnz.copy();

        boolean isReducible = false;
        boolean goon = true;

        while (goon) {
            Matrix Qnnz_abs = Qnnz.copy();
            Qnnz_abs.absEq();
            Matrix Qnnz_abs_sum_col = Qnnz_abs.sumCols();
            Matrix Qnnz_abs_sum_rows = Qnnz_abs.sumRows();
            Matrix find_res = new Matrix(1, Qnnz_abs_sum_col.getNumCols());

            for (int i = 0; i < Qnnz_abs_sum_col.getNumCols(); i++) {
                if (Qnnz_abs_sum_col.get(i) != 0.0 && Qnnz_abs_sum_rows.get(i) != 0.0) {
                    find_res.set(0, i, 1.0);
                }
            }
            nnzel = find_res.find().transpose();

            if (nnzel.length() < n && !isReducible) {
                isReducible = true;
            }

            Matrix new_Qnnz = new Matrix(nnzel.getNumCols(), nnzel.getNumCols());
            for (int i = 0; i < nnzel.getNumCols(); i++) {
                for (int j = 0; j < nnzel.getNumCols(); j++) {
                    double matrixValue = Qnnz.get((int) nnzel.get(0, i), (int) nnzel.get(0, j));
                    if (matrixValue != 0.0) new_Qnnz.set(i, j, matrixValue);
                }
            }
            Qnnz = new_Qnnz;

            Matrix new_bnnz = new Matrix(nnzel.getNumCols(), 1);
            for (int i = 0; i < nnzel.getNumCols(); i++) {
                new_bnnz.set(i, 0, bnnz.get((int) nnzel.get(0, i), 0));
            }
            bnnz = new_bnnz;

            Qnnz = Ctmc_makeinfgen.ctmc_makeinfgen(Qnnz);

            if ((Qnnz.getNumCols() * Qnnz.getNumRows() == Qnnz_1.getNumCols() * Qnnz_1.getNumRows())
                    && (bnnz.getNumCols() * bnnz.getNumRows() == bnnz_1.getNumCols() * bnnz_1.getNumRows())) {
                goon = false;
            } else {
                Qnnz_1 = Qnnz.copy();
                bnnz_1 = bnnz.copy();
                nnzel = new Matrix(1, Qnnz.length());
                for (int i = 0; i < Qnnz.length(); i++) nnzel.set(0, i, (double) i);
            }
        }

        if (Qnnz.isEmpty()) {
            p.fill(1.0 / n);
            return new Object[] { p, Qmat, nConnComp, connComp };
        }

        int lastCol = Qnnz.getNumCols() - 1;
        int adjustedRefstate = (refstate < Qnnz.getNumRows()) ? refstate : 0;

        for (int i = 0; i < Qnnz.getNumRows(); i++) {
            Qnnz.set(i, lastCol, 0.0);
        }
        Qnnz.set(adjustedRefstate, lastCol, 1.0);
        bnnz.set(lastCol, 0, 1.0);

        Matrix solve_res = new Matrix(0, 0);
        Matrix.solveSafe(Qnnz.transpose(), bnnz, solve_res);

        for (int i = 0; i < nnzel.getNumCols(); i++) {
            p.set(0, (int) nnzel.get(0, i), solve_res.get(i, 0));
        }

        return new Object[] { p, Qmat, nConnComp, connComp };
    }

    public static Object[] ctmc_relsolve(Matrix Q, int refstate) {
        return ctmc_relsolve(Q, refstate, null);
    }

    public static Object[] ctmc_relsolve(Matrix Q) {
        return ctmc_relsolve(Q, 0, null);
    }
}
