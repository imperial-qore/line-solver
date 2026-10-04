/**
 * @file Exact MAP/PH/c FCFS queue
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

import java.util.ArrayList;
import java.util.List;

import jline.api.mam.LdqbdMphc;
import jline.api.mam.Qbd_R_logred;
import jline.util.matrix.Matrix;

/**
 * The MAP/PH/c FCFS queue, solved exactly.
 *
 * THE STATE SPACE. With c identical servers the server identities carry no
 * information, so the service phases are held as a MULTISET: a configuration is
 * n = (n_1..n_ms) with sum(n) = k servers busy in phase i. There are
 * binomial(ms+k-1,k) of them, the count of Asmussen and Moller (2001), against
 * ms^k for the ordered space. Levels 0..c-1 are the boundary (level = servers
 * busy), levels &gt;= c repeat and carry the queue.
 *
 * THE WAITING TIME. An arrival that finds j customers waiting ahead of it waits
 * for j+1 service completions, so Wq is the (j+1)-st event time of the
 * configuration MAP (Lc, Cdep) started at the arrival-epoch configuration.
 * Folding the matrix-geometric level distribution over j gives the LINEAR
 * matrix ODE G'(t) = G Lj + R G Cj with G(0) = (I-R)^-1 kron(D1,I)/lambda and
 * P(Wq &gt; t) = pi_c G(t) e, so Wq is matrix-exponential. Its transform obeys
 * the generalized Sylvester equation g(sI-Lj) - R g Cj = G(0), and every moment
 * reuses that one operator with a different right-hand side.
 *
 * References:
 * S. Asmussen and J.R. Moller, "Calculation of the steady state waiting time
 * distribution in GI/PH/c and MAP/PH/c queues", Queueing Systems 37(1):9-29,
 * 2001.
 * D.P. Gaver, P.A. Jacobs, G. Latouche, "Finite birth-and-death models in
 * randomly changing environments", Adv. Appl. Probab. 16:715-731, 1984.
 */
public final class Qsys_mapphc {
    private Qsys_mapphc() {}

    public static QsysMapPhcResult qsys_mapphc(Matrix D0, Matrix D1, Matrix alpha, Matrix S, int c) {
        return qsys_mapphc(D0, D1, alpha, S, c, 500, 3, null);
    }

    /**
     * @param D0         arrival MAP hidden block, order ma
     * @param D1         arrival MAP arrival block, order ma
     * @param alpha      PH service initial vector, order ms
     * @param S          PH service sub-generator, order ms
     * @param c          number of servers
     * @param maxNumComp cap on the queue length probabilities returned
     * @param numWMoms   how many waiting-time moments to return
     * @param wPoints    times at which to evaluate P(Wq &gt; t), or null
     * @return the exact solution
     */
    public static QsysMapPhcResult qsys_mapphc(Matrix D0, Matrix D1, Matrix alpha, Matrix S, int c,
                                               int maxNumComp, int numWMoms, Matrix wPoints) {
        final int ma = D0.getNumRows();
        final int ms = S.getNumRows();
        if (D0.getNumCols() != ma || D1.getNumRows() != ma || D1.getNumCols() != ma) {
            throw new IllegalArgumentException("D0 and D1 must be square and of equal order");
        }
        if (S.getNumCols() != ms || alpha.length() != ms) {
            throw new IllegalArgumentException("alpha and S must have matching order");
        }
        if (c < 1) {
            throw new IllegalArgumentException("The number of servers c must be at least one");
        }

        final Matrix onesMs = Matrix.ones(ms, 1);
        final Matrix s0 = S.mult(onesMs).scale(-1.0);

        final Matrix theta = statLeftNull(D0.add(D1));
        final double lambda = theta.mult(D1).mult(Matrix.ones(ma, 1)).get(0, 0);
        final double meanService = -alpha.mult(S.inv()).mult(onesMs).get(0, 0);
        final double rho = lambda * meanService / c;
        if (rho >= 1.0) {
            throw new IllegalArgumentException("The load " + rho + " of the system is not below one");
        }

        // Multiset configurations of every occupancy 0..c
        final List<int[][]> cfg = new ArrayList<int[][]>();
        for (int k = 0; k <= c; k++) {
            cfg.add(multisets(ms, k));
        }

        final List<Matrix> Lcfg = new ArrayList<Matrix>();
        final List<Matrix> Up = new ArrayList<Matrix>();
        final List<Matrix> Dn = new ArrayList<Matrix>();
        for (int k = 0; k <= c; k++) {
            final int[][] Ck = cfg.get(k);
            final int nk = Ck.length;
            final Matrix Lk = new Matrix(nk, nk);
            for (int row = 0; row < nk; row++) {
                final int[] n = Ck[row];
                for (int i = 0; i < ms; i++) {
                    if (n[i] == 0) continue;
                    for (int j = 0; j < ms; j++) {
                        if (j == i) continue;
                        final int col = find(Ck, moveOne(n, i, j));
                        Lk.set(row, col, Lk.get(row, col) + n[i] * S.get(i, j));
                    }
                    Lk.set(row, row, Lk.get(row, row) + n[i] * S.get(i, i));
                }
            }
            Lcfg.add(Lk);

            if (k < c) {
                final int[][] Ck1 = multisets(ms, k + 1);
                final Matrix Uk = new Matrix(nk, Ck1.length);
                for (int row = 0; row < nk; row++) {
                    for (int j = 0; j < ms; j++) {
                        final int col = find(Ck1, addOne(Ck[row], j));
                        Uk.set(row, col, Uk.get(row, col) + alpha.get(j));
                    }
                }
                Up.add(Uk);
            } else {
                Up.add(null);
            }

            if (k > 0) {
                final int[][] Ckm = cfg.get(k - 1);
                final Matrix Dk = new Matrix(nk, Ckm.length);
                for (int row = 0; row < nk; row++) {
                    final int[] n = Ck[row];
                    for (int i = 0; i < ms; i++) {
                        if (n[i] == 0) continue;
                        final int col = find(Ckm, removeOne(n, i));
                        Dk.set(row, col, Dk.get(row, col) + n[i] * s0.get(i));
                    }
                }
                Dn.add(Dk);
            } else {
                Dn.add(null);
            }
        }

        // Completion WITH an immediate restart: the repeating down block
        final int[][] Cc = cfg.get(c);
        final int nc = Cc.length;
        final Matrix Cdep = new Matrix(nc, nc);
        for (int row = 0; row < nc; row++) {
            final int[] n = Cc[row];
            for (int i = 0; i < ms; i++) {
                if (n[i] == 0) continue;
                for (int j = 0; j < ms; j++) {
                    final int col = find(Cc, moveOne(n, i, j));
                    Cdep.set(row, col, Cdep.get(row, col) + n[i] * s0.get(i) * alpha.get(j));
                }
            }
        }

        final Matrix Ima = Matrix.eye(ma);
        final Matrix Inc = Matrix.eye(nc);
        final Matrix A_up = D1.kron(Inc);
        final Matrix A_loc = D0.kron(Inc).add(Ima.kron(Lcfg.get(c)));
        final Matrix A_dn = Ima.kron(Cdep);
        final Matrix R = Qbd_R_logred.qbd_R_logred(A_dn, A_loc, A_up);

        // Boundary levels 0..c, with the tail folded into level c through R
        final int[] sz = new int[c + 1];
        final int[] off = new int[c + 2];
        for (int k = 0; k <= c; k++) {
            sz[k] = ma * cfg.get(k).length;
            off[k + 1] = off[k] + sz[k];
        }
        final int tot = off[c + 1];
        final Matrix Q = new Matrix(tot, tot);
        for (int k = 0; k <= c; k++) {
            final Matrix Ick = Matrix.eye(cfg.get(k).length);
            final Matrix diag = (k < c)
                    ? D0.kron(Ick).add(Ima.kron(Lcfg.get(k)))
                    : A_loc.add(R.mult(A_dn));
            placeBlock(Q, diag, off[k], off[k]);
            if (k < c) {
                placeBlock(Q, D1.kron(Up.get(k)), off[k], off[k + 1]);
            }
            if (k > 0) {
                placeBlock(Q, Ima.kron(Dn.get(k)), off[k], off[k - 1]);
            }
        }
        Matrix piVec = statLeftNull(Q);

        final int nOp = nc * ma;
        final Matrix ImR = Matrix.eye(nOp).sub(R);
        final Matrix ImRinv = ImR.inv();
        Matrix piC = rowSlice(piVec, off[c], sz[c]);
        double massBoundary = 0.0;
        for (int i = 0; i < off[c]; i++) massBoundary += piVec.get(i);
        final double massTail = piC.mult(ImRinv).mult(Matrix.ones(nOp, 1)).get(0, 0);
        piVec = piVec.scale(1.0 / (massBoundary + massTail));
        piC = rowSlice(piVec, off[c], sz[c]);

        // Queue length distribution
        final List<Double> ql = new ArrayList<Double>();
        for (int k = 0; k < c; k++) {
            double s = 0.0;
            for (int i = 0; i < sz[k]; i++) s += piVec.get(off[k] + i);
            ql.add(s);
        }
        Matrix tail = piC.copy();
        double acc = 0.0;
        for (int k = 0; k < c; k++) acc += ql.get(k);
        double head = tail.elementSum();
        ql.add(head);
        acc += head;
        while (acc < 1 - 1e-12 && ql.size() < maxNumComp) {
            tail = tail.mult(R);
            final double v = tail.elementSum();
            ql.add(v);
            acc += v;
        }
        // E[N] in CLOSED FORM. maxNumComp caps the probabilities RETURNED, not the
        // mean: summing the truncated list loses the matrix-geometric tail, which at
        // rho -> 1 carries a first-order share of the mass. With pi_{c+j} = pi_c R^j,
        // sum_j (c+j) pi_c R^j e = pi_c [c (I-R)^-1 + R (I-R)^-2] e.
        double meanQL = 0.0;
        for (int k = 0; k < c; k++) meanQL += k * ql.get(k);
        final Matrix uTail = ImRinv.mult(Matrix.ones(nOp, 1));
        meanQL += c * piC.mult(uTail).get(0, 0)
                + piC.mult(R.mult(ImRinv.mult(uTail))).get(0, 0);

        // Waiting time
        final Matrix Lj = Ima.kron(Lcfg.get(c));
        final Matrix Cj = Ima.kron(Cdep);
        final Matrix G0 = ImRinv.mult(D1.kron(Inc)).scale(1.0 / lambda);
        final double probWait = piC.mult(G0).mult(Matrix.ones(nOp, 1)).get(0, 0);

        // X(-Lj) - R X Cj = rhs, vectorized column-major
        final Matrix Kop = Lj.scale(-1.0).transpose().kron(Matrix.eye(nOp))
                .sub(Cj.transpose().kron(R));
        final Matrix wMoms = new Matrix(1, Math.max(numWMoms, 1));
        Matrix gPrev = null;
        for (int k = 1; k <= numWMoms; k++) {
            final Matrix rhs = (k == 1) ? G0 : gPrev.scale(-(k - 1));
            final Matrix g = unvec(Kop.inv().mult(vec(rhs)), nOp);
            final double mk = k * ((k % 2 == 1) ? 1.0 : -1.0)
                    * piC.mult(g).mult(Matrix.ones(nOp, 1)).get(0, 0);
            wMoms.set(0, k - 1, mk);
            gPrev = g;
        }
        final double meanWT = numWMoms >= 1 ? wMoms.get(0, 0) : Double.NaN;

        Matrix wCCDF = null;
        if (wPoints != null && wPoints.length() > 0) {
            final Matrix Kt = Lj.transpose().kron(Matrix.eye(nOp)).add(Cj.transpose().kron(R));
            final Matrix v0 = vec(G0);
            wCCDF = new Matrix(1, wPoints.length());
            for (int it = 0; it < wPoints.length(); it++) {
                final Matrix Gt = unvec(Kt.scale(wPoints.get(it)).expm().mult(v0), nOp);
                wCCDF.set(0, it, piC.mult(Gt).mult(Matrix.ones(nOp, 1)).get(0, 0));
            }
        }

        final Matrix qlM = new Matrix(1, ql.size());
        for (int i = 0; i < ql.size(); i++) qlM.set(0, i, ql.get(i));

        return new QsysMapPhcResult(meanQL, meanWT, meanWT + meanService, rho, qlM, wMoms,
                wCCDF, wPoints, probWait, nc, "LINE:MAP/PH/" + c);
    }

    /**
     * Left null vector of G normalized to sum one. The R-corrected level-c block
     * has nonzero row sums, so the augmented system is used rather than a
     * generator solve.
     */
    private static Matrix statLeftNull(Matrix G) {
        final int n = G.getNumRows();
        // G is singular of rank n-1, so one of the n equations x G = 0 is implied
        // by the others: replace that column by the normalization x e = 1.
        final Matrix A = G.copy();
        for (int i = 0; i < n; i++) A.set(i, n - 1, 1.0);
        final Matrix rhs = new Matrix(n, 1);
        rhs.set(n - 1, 0, 1.0);
        return A.transpose().inv().mult(rhs).transpose();
    }

    /**
     * Rows are the compositions of k into ms nonnegative parts, in a fixed order.
     * Shared with the exact M/PH/c LD-QBD blocks, so a configuration index means
     * the same thing in both.
     */
    static int[][] multisets(int ms, int k) {
        return LdqbdMphc.ph_multisets(ms, k);
    }

    private static int find(int[][] rows, int[] key) {
        for (int i = 0; i < rows.length; i++) {
            boolean eq = true;
            for (int j = 0; j < key.length; j++) {
                if (rows[i][j] != key[j]) { eq = false; break; }
            }
            if (eq) return i;
        }
        throw new IllegalStateException("configuration not found");
    }

    private static int[] moveOne(int[] n, int from, int to) {
        final int[] m = n.clone();
        m[from]--;
        m[to]++;
        return m;
    }

    private static int[] addOne(int[] n, int at) {
        final int[] m = n.clone();
        m[at]++;
        return m;
    }

    private static int[] removeOne(int[] n, int at) {
        final int[] m = n.clone();
        m[at]--;
        return m;
    }

    private static void placeBlock(Matrix target, Matrix block, int r0, int c0) {
        for (int i = 0; i < block.getNumRows(); i++) {
            for (int j = 0; j < block.getNumCols(); j++) {
                target.set(r0 + i, c0 + j, block.get(i, j));
            }
        }
    }

    private static Matrix rowSlice(Matrix v, int off, int len) {
        final Matrix out = new Matrix(1, len);
        for (int i = 0; i < len; i++) out.set(0, i, v.get(off + i));
        return out;
    }

    /** Column-major vectorization, matching the MATLAB reference. */
    private static Matrix vec(Matrix X) {
        final int n = X.getNumRows(), m = X.getNumCols();
        final Matrix out = new Matrix(n * m, 1);
        for (int j = 0; j < m; j++) {
            for (int i = 0; i < n; i++) out.set(j * n + i, 0, X.get(i, j));
        }
        return out;
    }

    private static Matrix unvec(Matrix v, int n) {
        final Matrix out = new Matrix(n, n);
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) out.set(i, j, v.get(j * n + i, 0));
        }
        return out;
    }
}
