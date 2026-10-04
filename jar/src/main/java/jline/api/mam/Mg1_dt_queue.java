package jline.api.mam;

import jline.lib.smc.MG1PiOptions;
import jline.lib.smc.MG1_pi;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Discrete-time single-server queue with batch DMAP arrivals, DBMAP/DMAP/1.
 *
 * <p>Solved on a slotted time scale under the late arrival system with delayed
 * access (LAS-DA): within a slot the service completion resolves first,
 * arrivals are appended at the end of the slot and cannot enter service before
 * the next one, and the level is read after both. This is the convention of
 * Q_DT_MAP_MAP_1, whose QBD blocks this reproduces for a single arrival per
 * slot, and of the LDES slotted engine.
 *
 * <p>The chain is M/G/1-type because a slot may deliver a batch: with arrival
 * matrices A_k and service pair (S0,S1),
 * A^(-1) = kron(A_0,S1), A^(k) = kron(A_k,S0) + kron(A_{k+1},S1),
 * B^(k) = kron(A_k,I), the boundary row holding the empty system where no
 * service runs.
 *
 * <p>MATLAB twin: mg1_dt_queue.m
 */
public final class Mg1_dt_queue {
    private Mg1_dt_queue() {}

    /** Queue length, utilization, throughput, pmf and departure process. */
    public static final class Result {
        public final double QN;
        public final double UN;
        public final double TN;
        public final Matrix ql;
        public final MatrixCell dep;

        public Result(double QN, double UN, double TN, Matrix ql, MatrixCell dep) {
            this.QN = QN;
            this.UN = UN;
            this.TN = TN;
            this.ql = ql;
            this.dep = dep;
        }
    }

    /** Level mass left above the departure-process cut. */
    private static final double TAIL_TOL = 1e-10;

    public static Result mg1_dt_queue(MatrixCell ARV, MatrixCell SVC, int maxNumComp,
                                      boolean wantDeparture) {
        Matrix S0 = SVC.get(0);
        Matrix S1 = SVC.get(1);
        int ms = S0.getNumRows();
        int ma = ARV.get(0).getNumRows();
        int K = ARV.size() - 1;
        Matrix Ims = Matrix.eye(ms);

        double lambda = Dmap_batch.dmap_lambda(ARV);
        double mu = Dmap_batch.dmap_lambda(SVC);
        if (lambda >= mu) {
            throw new RuntimeException("The discrete-time load " + (lambda / mu)
                    + " of the station is not below one (" + lambda
                    + " arrivals per slot against " + mu + " completions per busy slot).");
        }

        int m = ma * ms;
        // A = [A^(-1) A^(0) ... A^(K)], the layout MG1_pi consumes
        Matrix Acat = new Matrix(m, m * (K + 2));
        setBlock(Acat, 0, ARV.get(0).kron(S1));
        for (int k = 0; k <= K; k++) {
            Matrix blk = ARV.get(k).kron(S0);
            if (k + 1 <= K) {
                blk = blk.add(1.0, ARV.get(k + 1).kron(S1));
            }
            setBlock(Acat, k + 1, blk);
        }
        Matrix Bcat = new Matrix(m, m * (K + 1));
        for (int k = 0; k <= K; k++) {
            setBlock(Bcat, k, ARV.get(k).kron(Ims));
        }

        // mg1_pi computes G itself from the solver named in the options; the
        // solver and mode strings are dereferenced there, so neither may be null.
        // FI rather than the reference's CR: the two agree to 1.6e-15 on these
        // blocks, but this port's mg1_cr costs 19.8 s against 25 ms at order 17.
        MG1PiOptions opts = new MG1PiOptions(null, maxNumComp, 0, "FI", false, "ShiftPWCR");
        Matrix pivec = MG1_pi.mg1_pi(Bcat, Acat, opts);
        Matrix ql = aggregateLevels(pivec, pivec.getNumCols() / m - 1, m);

        int nlev = ql.getNumCols();
        double QN = 0;
        for (int i = 0; i < nlev; i++) {
            QN += i * ql.get(0, i);
        }
        double UN = 1.0 - ql.get(0, 0);

        MatrixCell dep = null;
        if (wantDeparture) {
            double cum = 0;
            int L = nlev - 1;
            for (int i = 0; i < nlev; i++) {
                cum += ql.get(0, i);
                if (cum > 1 - TAIL_TOL) {
                    L = i;
                    break;
                }
            }
            dep = levelChain(ARV, S0, S1, Ims, m, K, Math.max(1, L));
        }
        return new Result(QN, UN, lambda, ql, dep);
    }

    /**
     * Transition matrix of the level-truncated chain, split by whether the slot
     * carries a departure: levels 0..L with arrivals that would cross L held at
     * L. (D0,D1) is therefore both the chain and the departure process seen by
     * the next station, so the two can never disagree.
     */
    private static MatrixCell levelChain(MatrixCell ARV, Matrix S0, Matrix S1, Matrix Ims,
                                         int m, int K, int L) {
        int nstates = (L + 1) * m;
        Matrix D0 = new Matrix(nstates, nstates);
        Matrix D1 = new Matrix(nstates, nstates);
        for (int i = 0; i <= L; i++) {
            for (int k = 0; k <= K; k++) {
                if (i == 0) {
                    // empty system: the slot carries no completion
                    addBlock(D0, i * m, Math.min(L, k) * m, ARV.get(k).kron(Ims));
                } else {
                    addBlock(D0, i * m, Math.min(L, i + k) * m, ARV.get(k).kron(S0));
                    addBlock(D1, i * m, Math.min(L, i - 1 + k) * m, ARV.get(k).kron(S1));
                }
            }
        }
        return new MatrixCell(D0, D1);
    }

    /** Level marginal of a phase-level stationary vector, renormalized. */
    private static Matrix aggregateLevels(Matrix pivec, int L, int m) {
        Matrix ql = new Matrix(1, L + 1);
        double total = 0;
        for (int i = 0; i <= L; i++) {
            double mass = 0;
            for (int j = 0; j < m; j++) {
                mass += pivec.get(i * m + j);
            }
            ql.set(0, i, mass);
            total += mass;
        }
        if (total > 0) {
            ql.scaleEq(1.0 / total);
        }
        return ql;
    }

    private static void setBlock(Matrix target, int blockIndex, Matrix block) {
        int m = block.getNumRows();
        for (int r = 0; r < m; r++) {
            for (int c = 0; c < block.getNumCols(); c++) {
                target.set(r, blockIndex * m + c, block.get(r, c));
            }
        }
    }

    private static void addBlock(Matrix target, int rowOffset, int colOffset, Matrix block) {
        for (int r = 0; r < block.getNumRows(); r++) {
            for (int c = 0; c < block.getNumCols(); c++) {
                target.set(rowOffset + r, colOffset + c,
                        target.get(rowOffset + r, colOffset + c) + block.get(r, c));
            }
        }
    }
}
