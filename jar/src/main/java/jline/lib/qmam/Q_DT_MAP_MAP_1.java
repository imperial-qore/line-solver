/**
 * @file Q_DT_MAP_MAP_1 - Discrete-Time DMAP/DMAP/1 Queue Analyzer
 *
 * Computes the queue length distribution of a discrete-time DMAP/DMAP/1/FCFS
 * queue. Port of Q_DT_MAP_MAP_1.m of the Q-MAM library by J. F. Perez,
 * J. Van Velthoven and B. Van Houdt (VALUETOOLS 2008).
 *
 * @since LINE 3.1.0
 */
package jline.lib.qmam;

import java.util.Map;

import jline.lib.smc.QBD_CR;
import jline.lib.smc.QBD_pi;
import jline.lib.smc.Stat;
import jline.util.matrix.Matrix;

public final class Q_DT_MAP_MAP_1 {
    private Q_DT_MAP_MAP_1() {}

    /**
     * Queue length distribution of a discrete-time DMAP/DMAP/1/FCFS queue.
     *
     * <p>Convention: late arrival system with delayed access. The QBD blocks
     * encode it directly, {@code A_1 = kron(C1,D0)} being the statement that a
     * job arriving at the end of a slot cannot be served within it, and
     * {@code A_0 = kron(C0,D0) + kron(C1,D1)} that a simultaneous arrival and
     * completion leave the level unchanged. The boundary row
     * {@code B_k = kron(C_k, I)} freezes the service phase while the system is
     * empty. This matches the LDES slotted engine, whose intra-slot order is
     * completion, internal movement, arrival, bookkeeping.
     *
     * <p>Only the queue length is returned. The waiting and sojourn pmfs that
     * Q_DT_MAP_MAP_1.m also computes are deliberately not ported: no LINE
     * caller consumes them, and the MATLAB discrete-time solver path reads the
     * queue length alone. They should be added here, not approximated by a
     * caller, when a discrete-time getCdfRespT is wired.
     *
     * @param C0 ma x ma matrix of the arrival DMAP, transitions without arrival
     * @param C1 ma x ma matrix of the arrival DMAP, transitions with arrival
     * @param D0 ms x ms matrix of the service DMAP, transitions without completion
     * @param D1 ms x ms matrix of the service DMAP, transitions with completion
     * @param options mode, maximum number of components and verbosity
     * @return the queue length distribution, entry i being Prob[i in system]
     */
    public static DTQueueResult qDtMapMap1(Matrix C0, Matrix C1, Matrix D0, Matrix D1,
                                           MAPMAP1Options options) {
        if (options == null) {
            options = new MAPMAP1Options();
        }
        int ma = C0.getNumRows();
        int ms = D0.getNumRows();
        int mtot = ma * ms;

        if (C0.getNumCols() != ma || C1.getNumRows() != ma || C1.getNumCols() != ma) {
            throw new IllegalArgumentException("Arrival process matrices must be ma x ma");
        }
        if (D0.getNumCols() != ms || D1.getNumRows() != ms || D1.getNumCols() != ms) {
            throw new IllegalArgumentException("Service process matrices must be ms x ms");
        }

        // Per-slot event rates of the two DMAPs, from the phase chain D0+D1
        Matrix piA = Stat.stat(C0.add(1.0, C1));
        double avga = piA.mult(C1).mult(Matrix.ones(ma, 1)).get(0);
        Matrix piS = Stat.stat(D0.add(1.0, D1));
        double avgs = piS.mult(D1).mult(Matrix.ones(ms, 1)).get(0);

        double rho = avga / avgs;
        if (rho >= 1) {
            throw new RuntimeException("The load " + rho + " of the system exceeds one");
        }

        Matrix eyeMs = Matrix.eye(ms);
        Matrix Am1 = C0.kron(D1);
        Matrix A0 = C0.kron(D0).add(1.0, C1.kron(D1));
        Matrix A1 = C1.kron(D0);
        Matrix Bm1 = C0.kron(D1);
        Matrix B0 = C0.kron(eyeMs);
        Matrix B1 = C1.kron(eyeMs);

        Map<String, Matrix> qbd = QBD_CR.QBD_CR(Am1, A0, A1, null,
                options.verbose > 0 ? Integer.valueOf(1) : null, null, null);
        Matrix R = qbd.get("R");

        // Boundary [B1; A0 + R*Am1]: the empty system has its own local block,
        // so the general-boundary branch of QBD_pi is required here
        Matrix boundary = Matrix.concatRows(B1, A0.add(1.0, R.mult(Am1)), null);
        Matrix stv = QBD_pi.QBD_pi(Bm1, B0, R, options.maxNumComp, options.verbose, boundary, 0);

        int nlev = stv.getNumCols() / mtot;
        Matrix ql = new Matrix(1, nlev);
        double total = 0;
        for (int i = 0; i < nlev; i++) {
            double mass = 0;
            for (int j = 0; j < mtot; j++) {
                mass += stv.get(0, i * mtot + j);
            }
            ql.set(0, i, mass);
            total += mass;
        }
        if (total > 0) {
            ql.scaleEq(1.0 / total);
        }

        return new DTQueueResult(ql);
    }
}
