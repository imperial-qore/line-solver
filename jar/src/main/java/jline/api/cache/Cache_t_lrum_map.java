/**
 * @file Characteristic times for LRU(m) caches under MAP request streams
 *
 * TTL approximation of LRU(m) with per-item Markovian arrival processes
 * (Gast and Van Houdt, Performance Evaluation 2017, Section 3.1.2). Each
 * item is modeled by an embedded Markov chain over (list, phase) states;
 * the level probability vectors obey pi_l = pi_0 prod_s R_s with the
 * R-recursions of eqs. (6)-(7), and pi_0 is the left Perron vector of
 * R_1 e^{D0 T_1} (level-0 balance). The characteristic times T_1..T_h are
 * fixed so that the expected occupancy of each list equals its capacity.
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import de.xypron.jcobyla.Calcfc;
import de.xypron.jcobyla.Cobyla;

import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Cache_t_lrum_map {
    private Cache_t_lrum_map() {}

    public static Matrix cache_t_lrum_map(final MatrixCell[] D0, final MatrixCell[] D1, final Matrix m) {
        final int n = D0.length;
        final int h = D0[0].size();

        double[] x = new double[h];
        for (int i = 0; i < h; i++) {
            x[i] = 1.0;
        }

        double rhobeg = 0.5;
        double rhoend = 1.0e-6;
        int maxFunEvals = (int) (500 * m.length()
                * Math.pow(10.0, Math.floor(Math.log10(n) / 4.0))); // maximum number of function evaluations

        Calcfc objectiveFunction = new Calcfc() {
            @Override
            public double compute(int comn, int comm, double[] comx, double[] comcon) {
                Matrix result = lrummapTime(comx, D0, D1, m, n, h);
                return result.norm();
            }
        };

        Cobyla.findMinimum(objectiveFunction, h, h, x, rhobeg, rhoend, 0, maxFunEvals); // 1: print result; 0: no print

        // Extract the optimized characteristic times
        Matrix t = new Matrix(1, h);
        for (int i = 0; i < h; i++) {
            t.set(0, i, x[i]);
        }

        return t;
    }

    /** Capacity residuals m_l - E[occupancy of list l] for times x. */
    public static Matrix lrummapTime(double[] x, MatrixCell[] D0Matrix, MatrixCell[] D1Matrix,
                                     Matrix m, int n, int h) {
        Matrix F = new Matrix(1, h);
        double[] occ = new double[h];
        for (int k = 0; k < n; k++) {
            LevelStats st = levelStats(x, D0Matrix[k], D1Matrix[k], h);
            for (int l = 0; l < h; l++) {
                occ[l] += st.occupancy[l];
            }
        }
        for (int l = 0; l < h; l++) {
            F.set(0, l, m.get(l) - occ[l]);
        }
        return F;
    }

    /** Per-item level statistics of the embedded (list, phase) chain. */
    static final class LevelStats {
        /** Time-stationary probability of residing in list l (index 0 = out of cache). */
        final double[] prob;
        /** Occupancy contribution of lists 1..h (prob columns 1..h). */
        final double[] occupancy;
        /** Fraction of the item's requests that hit in list l=1..h (index l-1). */
        final double[] hitfrac;

        LevelStats(double[] prob, double[] occupancy, double[] hitfrac) {
            this.prob = prob;
            this.occupancy = occupancy;
            this.hitfrac = hitfrac;
        }
    }

    /**
     * Level statistics for one item with per-list MAP cells (D0c, D1c) and
     * characteristic times x, following eqs. (5)-(9) of the reference.
     */
    static LevelStats levelStats(double[] x, MatrixCell D0c, MatrixCell D1c, int h) {
        int d = D0c.get(0).getNumRows();

        Matrix[] expD0 = new Matrix[h];
        Matrix[] trans = new Matrix[h]; // A_l of the paper, l=1..h at index l-1
        Matrix[] Nh = new Matrix[h];
        for (int l = 0; l < h; l++) {
            Matrix expD0l = Maths.matrixExp(Matrix.scaleMult(D0c.get(l), x[l]));
            expD0[l] = expD0l;
            Matrix iD0 = (Matrix.negative(D0c.get(l))).inv();
            Nh[l] = Matrix.oneMinusMatrix(expD0l).mult(iD0);
            trans[l] = Nh[l].mult(D1c.get(l));
        }
        Matrix A0 = (Matrix.negative(D0c.get(0))).inv().mult(D1c.get(0));
        Matrix N0 = (Matrix.negative(D0c.get(0))).inv();

        // R recursion, eqs. (6)-(7); R[l] holds the paper's R_{l+1}
        Matrix[] R = new Matrix[h];
        for (int l = h - 1; l >= 0; l--) {
            if (l == h - 1) {
                Matrix Aprev = (h == 1) ? A0 : trans[l - 1];
                R[l] = Aprev.mult(Matrix.oneMinusMatrix(trans[l]).inv());
            } else if (l == 0) {
                R[l] = A0.mult(Matrix.oneMinusMatrix(R[l + 1].mult(expD0[l + 1])).inv());
            } else {
                R[l] = trans[l - 1].mult(Matrix.oneMinusMatrix(R[l + 1].mult(expD0[l + 1])).inv());
            }
        }

        // pi_0: left Perron vector of R_1 e^{D0 T_1} (level-0 balance)
        Matrix pi0 = perronLeft(R[0].mult(expD0[0]), d);

        Matrix[] pih = new Matrix[h];
        pih[0] = pi0.mult(R[0]);
        for (int l = 1; l < h; l++) {
            pih[l] = pih[l - 1].mult(R[l]);
        }

        Matrix e = Matrix.ones(d, 1);
        double[] holding = new double[h + 1];
        holding[0] = pi0.mult(N0).mult(e).get(0);
        for (int l = 1; l <= h; l++) {
            holding[l] = pih[l - 1].mult(Nh[l - 1]).mult(e).get(0);
        }
        double denom = 0.0;
        for (int l = 0; l <= h; l++) {
            denom += holding[l];
        }

        double[] prob = new double[h + 1];
        double[] occupancy = new double[h];
        double[] hitfrac = new double[h];
        double lambdaK = mapRate(D0c.get(0), D1c.get(0), d);
        for (int l = 0; l <= h; l++) {
            prob[l] = holding[l] / denom;
        }
        for (int l = 1; l <= h; l++) {
            occupancy[l - 1] = prob[l];
            double hitTput = pih[l - 1].mult(Nh[l - 1]).mult(D1c.get(l - 1)).mult(e).get(0) / denom;
            hitfrac[l - 1] = lambdaK > 0 ? hitTput / lambdaK : 0.0;
        }
        return new LevelStats(prob, occupancy, hitfrac);
    }

    /** Stationary request rate of the item's MAP. */
    static double mapRate(Matrix D0, Matrix D1, int d) {
        // stationary phase vector of the generator D0+D1 via power iteration
        // on the uniformized transition matrix
        Matrix Q = D0.add(1.0, D1);
        double qmax = 0.0;
        for (int i = 0; i < d; i++) {
            qmax = Math.max(qmax, -Q.get(i, i));
        }
        Matrix P = Matrix.eye(d);
        if (qmax > 0) {
            P = P.add(1.0 / (1.1 * qmax), Q);
        }
        Matrix pi = perronLeft(P, d);
        return pi.mult(D1).mult(Matrix.ones(d, 1)).get(0);
    }

    /** Left Perron vector (normalized to sum 1) of a nonnegative matrix. */
    static Matrix perronLeft(Matrix M, int d) {
        Matrix pi = new Matrix(1, d);
        pi.fill(1.0 / d);
        for (int it = 0; it < 5000; it++) {
            Matrix next = pi.mult(M);
            double s = next.elementSum();
            if (s <= 0) {
                return pi;
            }
            next = Matrix.scaleMult(next, 1.0 / s);
            double diff = 0.0;
            for (int i = 0; i < d; i++) {
                diff = Math.max(diff, Math.abs(next.get(i) - pi.get(i)));
            }
            pi = next;
            if (diff < 1e-15) {
                break;
            }
        }
        return pi;
    }
}
