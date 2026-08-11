/**
 * Joint Queue-Length Probability for Product-Form Networks
 *
 * @since LINE 3.0
 */
package jline.api.pfqn;

import jline.api.pfqn.nc.Pfqn_ca;
import jline.io.Ret;
import jline.lib.perm.Permanent;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_joint {
    private Pfqn_joint() {}

    /**
     * Compute the joint queue-length probability for vector n
     */
    public static double pfqn_joint(Matrix n, Matrix L, Matrix N, Matrix Z, Double lGN) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        Matrix zThink = (Z != null) ? Z : Matrix.zeros(1, R);

        double lGn;
        if (lGN != null) {
            lGn = lGN;
        } else {
            Ret.pfqnNc result = Pfqn_ca.pfqn_ca(L, N, zThink);
            lGn = result.lG;
        }

        if (n.getNumCols() == 1) {
            return computeTotalQueueLengthProb(n, L, N, zThink, lGn);
        } else if (n.getNumCols() == R) {
            return computePerClassQueueLengthProb(n, L, N, zThink, lGn, M, R);
        } else {
            throw new IllegalArgumentException("Invalid argument to pfqn_joint: n must have 1 or " + R + " columns");
        }
    }

    public static double pfqn_joint(Matrix n, Matrix L, Matrix N, Matrix Z) {
        return pfqn_joint(n, L, N, Z, null);
    }

    public static double pfqn_joint(Matrix n, Matrix L, Matrix N) {
        return pfqn_joint(n, L, N, null, null);
    }

    /**
     * Compute joint probability for total queue lengths (n is M x 1)
     */
    private static double computeTotalQueueLengthProb(Matrix n, Matrix L, Matrix N, Matrix Z, double lGn) {
        boolean hasThinkTime = Z.elementSum() > 0.0;

        if (hasThinkTime) {
            double n0 = N.elementSum() - n.elementSum();
            double Fjoint = fper(Matrix.concatRows(L, Z, null), N,
                    Matrix.concatRows(n, Matrix.singleton(n0), null));
            double logFjoint = Math.log(Fjoint);
            double logFactorial = Maths.factln((int) n0);
            return Math.exp(logFjoint - lGn - logFactorial);
        } else {
            double Fjoint = fper(L, N, n);
            double logFjoint = Math.log(Fjoint);
            return Math.exp(logFjoint - lGn);
        }
    }

    /**
     * Compute joint probability for per-class queue lengths (n is M x R)
     */
    private static double computePerClassQueueLengthProb(Matrix n, Matrix L, Matrix N, Matrix Z,
                                                         double lGn, int M, int R) {
        Matrix n0 = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            double colSum = 0.0;
            for (int i = 0; i < M; i++) {
                colSum += n.get(i, r);
            }
            n0.set(0, r, N.get(0, r) - colSum);
        }

        double Fjoint = 0.0;

        boolean hasThinkTime = Z.elementSum() > 0.0;
        if (hasThinkTime) {
            for (int r = 0; r < R; r++) {
                if (n0.get(0, r) > 0) {
                    Fjoint += n0.get(0, r) * Math.log(Z.get(0, r));
                    Fjoint -= Maths.factln((int) n0.get(0, r));
                }
            }
        }

        for (int i = 0; i < M; i++) {
            Matrix nRow = new Matrix(1, R);
            for (int r = 0; r < R; r++) {
                nRow.set(0, r, n.get(i, r));
            }
            Fjoint += Maths.multinomialln(nRow);

            for (int r = 0; r < R; r++) {
                if (n.get(i, r) > 0 && L.get(i, r) > 0) {
                    Fjoint += n.get(i, r) * Math.log(L.get(i, r));
                }
            }
        }

        return Math.exp(Fjoint - lGn);
    }

    /**
     * Helper function F_per: computes permanent-based probability term
     */
    private static double fper(Matrix L, Matrix N, Matrix m) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        Matrix Ak = null;
        for (int r = 0; r < R; r++) {
            int nRep = (int) N.get(0, r);
            Matrix col = Matrix.extractColumn(L, r, null);
            Matrix replicatedCols = col.repmat(1, nRep);
            if (Ak == null) {
                Ak = replicatedCols;
            } else {
                Ak = Matrix.concatColumns(Ak, replicatedCols, null);
            }
        }

        Matrix A = null;
        for (int i = 0; i < M; i++) {
            int mi = (int) m.get(i, 0);
            if (mi > 0) {
                Matrix rowToReplicate = Matrix.extractRows(Ak, i, i + 1, null);
                Matrix replicatedRows = rowToReplicate.repmat(mi, 1);
                if (A == null) {
                    A = replicatedRows;
                } else {
                    A = Matrix.concatRows(A, replicatedRows, null);
                }
            }
        }

        if (A == null || A.getNumRows() == 0) {
            return 1.0;
        }

        if (A.getNumRows() != A.getNumCols()) {
            return 0.0;
        }

        Permanent permanent = new Permanent(A, true);
        double permValue = permanent.value;

        double logProdFactorial = 0.0;
        for (int r = 0; r < R; r++) {
            logProdFactorial += Maths.factln((int) N.get(0, r));
        }

        return permValue / Math.exp(logProdFactorial);
    }
}
