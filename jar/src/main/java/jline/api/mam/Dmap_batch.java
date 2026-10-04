package jline.api.mam;

import jline.lib.butools.dph.CanonicalFromDPH2.DPH2Representation;
import jline.lib.butools.dph.DPH2From3Moments;
import jline.lib.smc.Stat;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Discrete-time batch arrival streams: event rate, superposition, Bernoulli
 * splitting and order reduction.
 *
 * <p>A batch stream is a MatrixCell {A_0, A_1, ...} where A_k carries the slots
 * delivering k events. A plain DMAP is the two-entry case. These four
 * operations are what a slotted decomposition needs between stations.
 *
 * <p>MATLAB twins: dmap_lambda.m, dmap_super.m, dmap_thin.m, dmap_compress.m,
 * dmap_compress_batch.m
 */
public final class Dmap_batch {
    private Dmap_batch() {}

    /**
     * Mean number of EVENTS per slot, pi * sum_k k*A_k * e. A slot carrying a
     * batch of two counts twice, which is what Little's law consumes
     * downstream.
     */
    public static double dmap_lambda(MatrixCell A) {
        int m = A.get(0).getNumRows();
        Matrix P = new Matrix(m, m);
        for (int k = 0; k < A.size(); k++) {
            P = P.add(1.0, A.get(k));
        }
        Matrix pi = Stat.stat(P);
        Matrix W = new Matrix(m, m);
        for (int k = 1; k < A.size(); k++) {
            W = W.add((double) k, A.get(k));
        }
        return pi.mult(W).mult(Matrix.ones(m, 1)).get(0);
    }

    /**
     * Superposition, E_k = sum_{i+j=k} kron(A_i, B_j).
     *
     * <p>NOT closed on DMAPs: two slotted streams fire in the same slot with
     * positive probability, so the merged stream carries batches. Folding E_2
     * into E_1 would conserve neither the arrival rate nor the slot in which
     * the work appears, so the batch dimension is kept and the downstream
     * station is solved as an M/G/1-type chain instead of a QBD.
     */
    public static MatrixCell dmap_super(MatrixCell A, MatrixCell B) {
        int p = A.size() - 1;
        int q = B.size() - 1;
        MatrixCell E = new MatrixCell(p + q + 1);
        for (int k = 0; k <= p + q; k++) {
            Matrix Ek = null;
            for (int i = Math.max(0, k - q); i <= Math.min(p, k); i++) {
                Matrix term = A.get(i).kron(B.get(k - i));
                Ek = Ek == null ? term : Ek.add(1.0, term);
            }
            E.set(k, Ek);
        }
        return E;
    }

    /**
     * Bernoulli thinning, B_k = sum_{n>=k} C(n,k) p^k (1-p)^(n-k) A_n. The
     * phase process is untouched, so this is exact for PROB/RAND routing.
     */
    public static MatrixCell dmap_thin(MatrixCell A, double p) {
        if (p < 0 || p > 1) {
            throw new RuntimeException("The routing probability must lie in [0,1], got " + p);
        }
        int n = A.size() - 1;
        int rows = A.get(0).getNumRows();
        int cols = A.get(0).getNumCols();
        MatrixCell B = new MatrixCell(n + 1);
        for (int k = 0; k <= n; k++) {
            Matrix Bk = new Matrix(rows, cols);
            for (int j = k; j <= n; j++) {
                double w = binomial(j, k) * Math.pow(p, k) * Math.pow(1 - p, j - k);
                if (w > 0) {
                    Bk = Bk.add(w, A.get(j));
                }
            }
            B.set(k, Bk);
        }
        // trailing zero batch levels carry no mass and only inflate the blocks
        while (B.size() > 2 && B.get(B.size() - 1).elementMaxAbs() < 1e-14) {
            B.remove(B.size() - 1);
        }
        return B;
    }

    /**
     * Order reduction of a DMAP by matching three interevent moments through
     * BuTools DPH2From3Moments, leaving it untouched while it stays within
     * maxOrder. Correlation is NOT preserved, mirroring the continuous-time
     * 'mixture.order1' compression, and that is where the multi-station
     * discrete-time path becomes approximate. Outside the DPH(2) region the
     * fallback keeps the exact mean with a Geometric, so the arrival rate is
     * conserved in every branch.
     */
    public static MatrixCell dmap_compress(MatrixCell DMAP, int maxOrder) {
        if (DMAP.get(0).getNumRows() <= maxOrder) {
            return DMAP;
        }
        double m1 = Dmap_moment.dmap_moment(DMAP, 1);
        double[] moms = new double[]{m1,
                Dmap_moment.dmap_moment(DMAP, 2),
                Dmap_moment.dmap_moment(DMAP, 3)};
        try {
            DPH2Representation rep = DPH2From3Moments.dph2From3Moments(moms);
            if (rep != null) {
                MatrixCell cand = Dph_to_dmap.dph_to_dmap(rep.beta, rep.B);
                if (isFeasible(cand)) {
                    return cand;
                }
            }
        } catch (RuntimeException e) {
            // fall through to the rate-preserving surrogate below
        }
        double p = Math.min(1.0, Math.max(1e-14, 1.0 / m1));
        Matrix alpha = new Matrix(1, 1);
        alpha.set(0, 0, 1.0);
        Matrix A = new Matrix(1, 1);
        A.set(0, 0, 1.0 - p);
        return Dph_to_dmap.dph_to_dmap(alpha, A);
    }

    /**
     * Order reduction of a batch stream. Keeps the two features the downstream
     * M/G/1-type solve consumes: the law of the time between NONEMPTY slots,
     * matched to three moments, and the stationary batch-size distribution
     * conditional on a nonempty slot, kept exactly. The event rate of the
     * reduced stream equals the original one by construction.
     */
    public static MatrixCell dmap_compress_batch(MatrixCell B, int maxOrder) {
        if (B.get(0).getNumRows() <= maxOrder) {
            return B;
        }
        int nb = B.size() - 1;
        int m = B.get(0).getNumRows();

        Matrix Ptot = new Matrix(m, m);
        for (int k = 0; k < B.size(); k++) {
            Ptot = Ptot.add(1.0, B.get(k));
        }
        Matrix piPhase = Stat.stat(Ptot);
        Matrix e = Matrix.ones(m, 1);
        double[] qraw = new double[nb];
        double massNonEmpty = 0;
        for (int k = 1; k <= nb; k++) {
            qraw[k - 1] = piPhase.mult(B.get(k)).mult(e).get(0);
            massNonEmpty += qraw[k - 1];
        }
        if (massNonEmpty <= 1e-14) {
            throw new RuntimeException("The batch stream carries no events, so it cannot be compressed.");
        }

        MatrixCell marked = new MatrixCell(B.get(0).copy(), Ptot.add(-1.0, B.get(0)));
        MatrixCell markedC = dmap_compress(marked, maxOrder);

        MatrixCell Bc = new MatrixCell(nb + 1);
        Bc.set(0, markedC.get(0));
        for (int k = 1; k <= nb; k++) {
            Bc.set(k, markedC.get(1).scale(qraw[k - 1] / massNonEmpty));
        }
        return Bc;
    }

    private static boolean isFeasible(MatrixCell DMAP) {
        Matrix D0 = DMAP.get(0);
        Matrix D1 = DMAP.get(1);
        int m = D0.getNumRows();
        if (D0.elementMin() < -1e-10 || D1.elementMin() < -1e-10) {
            return false;
        }
        Matrix rs = D0.add(1.0, D1).mult(Matrix.ones(m, 1));
        for (int i = 0; i < m; i++) {
            if (Math.abs(rs.get(i, 0) - 1.0) > 1e-6) {
                return false;
            }
        }
        return true;
    }

    private static double binomial(int n, int k) {
        double c = 1;
        for (int i = 0; i < k; i++) {
            c = c * (n - i) / (i + 1);
        }
        return c;
    }
}
