/**
 * Exact MVA recursion carrying the interlocked-flow correction.
 *
 * The correction of Franks (1999), Ch. 4, Eq. (4.7) replaces the arrival theorem term
 * Q(n-1_s,i) by a per-class weighted sum, so the recursion has to carry per-class queue
 * lengths that {@link Pfqn_mva} does not need. Closed single-server models only.
 *
 * The discounted arrival-instant queue is floored at the in-service component, as in lqns
 * MVA::queueOnly_adjusted, so the correction damps itself out as a station saturates. That
 * is a self-limiting guard, NOT a hard capacity test: sum_s XN(s)*L(i,s) &lt;= mi(i) is
 * still asserted nowhere. See git show 449847e7b:_kb/log.md.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import org.apache.commons.math3.util.FastMath;

import jline.api.pfqn.Pfqn_replicas;
import jline.io.Ret;
import jline.util.matrix.Matrix;
import jline.util.Triple;

public final class Pfqn_mva_ilock {
    private Pfqn_mva_ilock() {}

    /**
     * Mean Value Analysis with the interlocked-flow correction of Franks (1999), Eq. (4.7).
     * IL(r,s) is the share of the class-s queue that a class-r arrival cannot see, because that
     * work was itself caused by the class-r request. IL is required; the model is not
     * product-form under it, so lGN comes back as NaN.
     */
    public static Ret.pfqnMVA pfqn_mva_ilock(Matrix L, Matrix N, Matrix Z, Matrix mi, Matrix IL) {
        N = N.copy(); // Create local copy to avoid modifying the original
        Z = Z.isEmpty() ? new Matrix(1, L.getNumCols()) : Z;
        Matrix XN; // throughputs
        Matrix QN;  // queue lengths
        Matrix UN;  // utilizations
        Matrix CN;  // residence times
        double lGN = 0.0; // log of the normalizing constant
        double InfServ = 1.0;
        if (Z.isEmpty() && mi == null) {
            InfServ = 0.0;
        }
        N = N.ceil();
        int M_original = L.getNumRows(); // Original number of stations
        int R = L.getNumCols(); // R Classes
        N = N.columnMajorOrder().transpose();

        // Normalize mi to row vector
        if (mi == null) {
            mi = Matrix.ones(1, M_original);
        } else {
            if (mi.getNumRows() > 1) {
                mi = mi.transpose();
            }
        }

        // see _kb/03-api-layer.md for rationale
        Matrix L_reduced = L;
        int[] mapping = new int[M_original];
        for (int i = 0; i < M_original; i++) mapping[i] = i;
        int M = M_original;
        if (!N.any()) {
            // see _kb/03-api-layer.md for rationale
            return new Ret.pfqnMVA(new Matrix(1, R), new Matrix(M, R), new Matrix(M, R), new Matrix(M, R), 0.0);
        }
        int NR = N.length();
        if (R != NR) {
            throw new RuntimeException("pfqn_mva_ilock: Demand matrix and population vector have different number of classes");
        }

        boolean hasIL = IL != null && !IL.isEmpty();
        if (!hasIL) {
            throw new RuntimeException("pfqn_mva_ilock: an interlock matrix is required; use pfqn_mva for the standard arrival theorem");
        }
        Matrix ILw = null;
        {
            if (IL.getNumRows() != R || IL.getNumCols() != R) {
                throw new RuntimeException("pfqn_mva_ilock: the interlock matrix must be nclasses x nclasses");
            }
            ILw = new Matrix(R, R);
            for (int r = 0; r < R; r++) {
                for (int s = 0; s < R; s++) {
                    double w = (r == s) ? 1.0 : 1.0 - IL.get(r, s);
                    ILw.set(r, s, Math.max(0.0, Math.min(1.0, w)));
                }
            }
        }

        XN = new Matrix(1, R);
        QN = new Matrix(M, R);
        CN = new Matrix(M, R);
        if (InfServ == 1.0) {
            Z = Z.columnMajorOrder().transpose();
        } else {
            Z = new Matrix(1, R);
        }
        Matrix prods = new Matrix(1, R - 1);
        for (int w = 0; w < R - 1; w++) {
            // ones(1, R-(w+1)+1)
            // (w + 2) instead of (w + 1) because Java indexing starts at 0
            Matrix o = Matrix.ones(1, R - (w + 2) + 1);
            // Addition: ones(1,R-(w+1)+1) + N(1, w+1:R)
            for (int i = 0; i < R - (w + 2) + 1; i++) {
                o.set(0, i, o.get(0, i) + N.get(0, w + 1 + i));
            }
            // Now take prod(o)
            prods.set(0, w, o.elementMult());
        }

        int firstNonEmpty = R - 1;
        while (firstNonEmpty >= 0 && N.get(0, firstNonEmpty) == 0.0) {
            firstNonEmpty--;
        }
        double totpop = Matrix.ones(1, N.getNumCols()).add(1.0, N).elementMult();
        double ctr = totpop;
        Matrix Q = new Matrix((int) totpop, M);
        Matrix[] Qc = new Matrix[R]; // per-class queue lengths, needed by the interlock
        Matrix[] Uc = new Matrix[R]; // per-class in-service component, the interlock's floor
        for (int r = 0; r < R; r++) {
            Qc[r] = new Matrix((int) totpop, M);
            Uc[r] = new Matrix((int) totpop, M);
        }
        int currentpop = 1;

        Matrix n = new Matrix(1, R);
        n.set(0, firstNonEmpty, 1);
        while (ctr > 0) {
            int s = 0;
            while (s < R) {
                int pos_n_1s = 0;
                if (n.get(0, s) > 0) {
                    n.set(0, s, n.get(0, s) - 1);
                    pos_n_1s = (int) n.get(0, R - 1);
                    // for w=1:R-1
                    int w = 0;
                    while (w < R - 1) {
                        pos_n_1s = (int) (pos_n_1s + n.get(0, w) * prods.get(0, w));
                        w++;
                    }
                    n.set(0, s, n.get(0, s) + 1);
                }
                double CNtot = 0.0;
                int i = 0;
                // Compute the residence times. Compute the total residence time as well to avoid another iteration
                // through all the stations.
                while (i < M) {
                    double Lis = L_reduced.get(i, s);
                    double qarv = 0.0;
                    for (int r = 0; r < R; r++) {
                        // In-service protection, as in lqns MVA::queueOnly_adjusted: the
                        // discount bites on the WAITING part only, never on the job already
                        // in service, so it damps itself out as the station saturates.
                        qarv += Math.max(ILw.get(s, r) * Qc[r].get(pos_n_1s, i), Uc[r].get(pos_n_1s, i));
                    }
                    CN.set(i, s, Lis * (mi.get(0, i) + qarv));
                    CNtot += CN.get(i, s);
                    i++;
                }
                // Compute the throughput for class s
                XN.set(0, s, n.get(0, s) / (Z.get(0, s) + CNtot));
                i = 0;
                // Compute the queue lengths
                while (i < M) {
                    QN.set(i, s, XN.get(0, s) * CN.get(i, s));
                    Q.set(currentpop, i, Q.get(currentpop, i) + QN.get(i, s));
                    Qc[s].set(currentpop, i, QN.get(i, s));
                    Uc[s].set(currentpop, i, XN.get(0, s) * L_reduced.get(i, s));
                    i++;
                }
                s++;
            }
            s = R - 1;
            while ((s >= 0 && (n.get(0, s) == N.get(0, s))) || s > firstNonEmpty) {
                s--;
            }
            Matrix nonZero = n.find();
            if (!nonZero.isEmpty()) {
                int nonZeroIdx = nonZero.getNumRows() - 1;
                while (nonZeroIdx >= 0 && n.get((int) nonZero.get(nonZeroIdx, 0)) <= 0) {
                    nonZeroIdx--;
                }
                int last_nnz = (int) nonZero.get(nonZeroIdx, 0);
                double sumn = 0.0;
                double sumN = 0.0;
                double sumnprime = 0.0;
                for (int i = 0; i < last_nnz; i++) {
                    sumn += n.get(0, i);
                    sumN += N.get(0, i);
                }
                for (int i = last_nnz + 1; i < R; i++) {
                    sumnprime += n.get(0, i);
                }
                if (sumn == sumN && sumnprime == 0.0) {
                    double logX = FastMath.log(XN.get(0, last_nnz));
                    lGN -= logX;
                }
            }
            if (s == -1) {
                break;
            }
            n.set(0, s, n.get(0, s) + 1);
            s++;
            while (s < R) {
                n.set(0, s, 0);
                s++;
            }
            ctr--;
            currentpop++;
        }

        lGN = Double.NaN; // the interlock leaves the model outside product form

        UN = new Matrix(M, R); // Utilizations
        for (int mIdx = 0; mIdx < M; mIdx++) {
            for (int r = 0; r < R; r++) {
                UN.set(mIdx, r, XN.get(0, r) * L_reduced.get(mIdx, r));
            }
        }

        // see _kb/03-api-layer.md for rationale

        // Expand results back to original dimensions if stations were consolidated
        if (M < M_original) {
            Triple<Matrix, Matrix, Matrix> expanded = Pfqn_replicas.pfqn_expand(QN, UN, CN, mapping);
            QN = expanded.getFirst();
            UN = expanded.getSecond();
            CN = expanded.getThird();
        }

        return new Ret.pfqnMVA(XN, QN, UN, CN, lGN);
    }
}
