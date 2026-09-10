/**
 * Joint Queue-Length Probability for Product-Form Networks
 *
 * @since LINE 3.0
 */
package jline.api.pfqn;

import jline.api.pfqn.nc.Pfqn_ca;
import jline.io.Ret;
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
     * Compute joint probability for total queue lengths (n is M x 1).
     *
     * Delegated to Pfqn_jointmarg so that the permanent identity lives in one
     * place; here the think time is a single aggregated delay row, which is the
     * (M+1)-th station.
     */
    private static double computeTotalQueueLengthProb(Matrix n, Matrix L, Matrix N, Matrix Z, double lGn) {
        boolean hasThinkTime = Z.elementSum() > 0.0;

        if (hasThinkTime) {
            double n0 = N.elementSum() - n.elementSum();
            if (n0 < 0) {
                return 0.0;
            }
            Matrix Zrow = Matrix.zeros(1, L.getNumCols());
            for (int r = 0; r < L.getNumCols(); r++) {
                double zr = 0.0;
                for (int zi = 0; zi < Z.getNumRows(); zi++) zr += Z.get(zi, r);
                Zrow.set(0, r, zr);
            }
            Matrix Lext = Matrix.concatRows(L, Zrow, null);
            Matrix next = Matrix.concatRows(n, Matrix.singleton(n0), null);
            return Pfqn_jointmarg.pfqn_jointmarg(next, Lext, N,
                    new int[]{L.getNumRows()}, Double.valueOf(lGn)).pjoint;
        } else {
            return Pfqn_jointmarg.pfqn_jointmarg(n, L, N, null, Double.valueOf(lGn)).pjoint;
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
                    // Column sum, not row 0: Z may carry one row per delay node,
                    // and the reference sums it (pfqn_joint.m:66,77 write
                    // sum(Z)). Reproducing that rather than assuming the caller
                    // summed; no current caller passes more than one row.
                    double zr = 0.0;
                    for (int zi = 0; zi < Z.getNumRows(); zi++) zr += Z.get(zi, r);
                    Fjoint += n0.get(0, r) * Math.log(zr);
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
}
