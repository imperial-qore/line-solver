/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.pfqn;

import java.util.ArrayList;
import java.util.List;

import jline.api.pfqn.nc.Pfqn_ca;
import jline.api.pfqn.nc.Pfqn_perm;
import jline.io.Ret;
import jline.lib.perm.AdaPartSampler;
import jline.lib.perm.BethePermanent;
import jline.lib.perm.HeuristicPermanent;
import jline.lib.perm.HuberLawSampler;
import jline.lib.perm.SaddlePointPermanent;
import jline.util.Maths;
import jline.util.matrix.Matrix;

/**
 * Joint probability of the per-station TOTAL queue lengths.
 *
 * <p>Joint probability that station i holds n(i) jobs in total, all classes
 * summed out, in a closed multiclass product-form network:</p>
 *
 * <pre>
 *   P(n_1,...,n_M) = perm(A) / ( prod_r N_r! * prod_{j in infset} n_j! * G(N) )
 * </pre>
 *
 * <p>with A the demand matrix whose column r is repeated N_r times and whose
 * row i is repeated n_i times, so A is square of order sum(N). Unlike
 * {@link Pfqn_joint}, which takes the delay as a single aggregated row, every
 * infinite-server station keeps its own row here and contributes its own
 * 1/n_j!: the queueing stations contribute the n_i! the permanent identity
 * supplies, the infinite servers do not.</p>
 *
 * <p>The identity holds for load-independent single-server queues plus infinite
 * servers. Multiserver and load-dependent stations break the n_i! factor and
 * are the caller's responsibility to exclude.</p>
 *
 * <p>ZERO ELEMENTS are safe under the exact engine and only under it: a station
 * holding no jobs contributes no row, a class with no jobs contributes no
 * column, a zero demand is an ordinary zero entry of A, and the permanent of
 * the empty matrix is 1. The approximate engines are REFUSED on a matrix with a
 * structural zero rather than having it floored at eps: Sinkhorn scaling needs
 * full support, and the Bethe gap is a state-dependent lower bound that does
 * not cancel when the estimates are normalized against each other.</p>
 *
 * <p>Reference: H. J. Ryser, "Combinatorial Mathematics", Carus Mathematical
 * Monographs 14, Mathematical Association of America, 1963.</p>
 */
public final class Pfqn_jointmarg {
    private Pfqn_jointmarg() {}

    /** Result of a joint total-queue-length evaluation. */
    public static class Ret_jointmarg {
        /** Joint probability. */
        public final double pjoint;
        /** Its logarithm, which survives populations pjoint underflows at. */
        public final double lpjoint;
        /** Log normalizing constant used, whether supplied or computed here. */
        public final double lG;

        public Ret_jointmarg(double pjoint, double lpjoint, double lG) {
            this.pjoint = pjoint;
            this.lpjoint = lpjoint;
            this.lG = lG;
        }
    }

    /**
     * Joint probability of the per-station total queue lengths.
     *
     * @param n      (1 x M) or (M x 1) per-station total job counts, infinite
     *               servers included
     * @param L      (M x R) demand matrix, infinite-server rows included
     * @param N      (1 x R) per-class populations
     * @param infset row indices of L that are infinite-server stations, may be
     *               null or empty
     * @param lGN    log normalizing constant, or null to compute it here
     * @param engine "exact" (default), "spm", "bethe", "heur", "huberlaw" or
     *               "adapart". "spm" is the only engine that does not expand the
     *               matrix to order sum(N): it takes the row-replicated matrix with
     *               the class populations as column multiplicities, which is the
     *               regime its saddle-point expansion is asymptotically exact in, so
     *               its cost does not grow with the population and its relative error
     *               is O((R-1)/min(N)). Measured on a 3-station 2-class model, 12.8%
     *               at N = (1,1), 4.2% at (3,3), 2.1% at (6,6); it degrades the other
     *               way round, when the class count grows at fixed population (2.7%
     *               at R = 2, 21% at R = 7, both at N_r = 3), because R-1 is the
     *               dimension being expanded in. The bias is nearly constant across
     *               the lattice, so a caller that renormalizes a full sweep keeps far
     *               less of it: total variation distance 5.0e-3 at N = (1,1), 8.4e-4
     *               at (3,3), 4.3e-4 at (5,5), better than "bethe" and "heur" at
     *               every population measured.
     * @return the probability, its logarithm and the constant used
     */
    public static Ret_jointmarg pfqn_jointmarg(Matrix n, Matrix L, Matrix N, int[] infset,
                                               Double lGN, String engine) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        String eng = (engine == null || engine.isEmpty()) ? "exact" : engine.toLowerCase();

        int[] nvec = toIntVector(n);
        int[] npop = toIntVector(N);
        if (nvec.length != M) {
            throw new IllegalArgumentException("pfqn_jointmarg: the occupancy vector has " + nvec.length
                    + " entries but L has " + M + " rows");
        }
        if (npop.length != R) {
            throw new IllegalArgumentException("pfqn_jointmarg: the population vector has " + npop.length
                    + " entries but L has " + R + " columns");
        }
        boolean[] isInf = new boolean[M];
        if (infset != null) {
            for (int k = 0; k < infset.length; k++) {
                if (infset[k] < 0 || infset[k] >= M) {
                    throw new IllegalArgumentException("pfqn_jointmarg: infset indexes a station outside 0.."
                            + (M - 1));
                }
                isInf[infset[k]] = true;
            }
        }

        int ntot = 0;
        for (int i = 0; i < M; i++) {
            if (nvec[i] < 0) {
                throw new IllegalArgumentException("pfqn_jointmarg: the occupancy vector has a negative entry");
            }
            ntot += nvec[i];
        }
        int npoptot = 0;
        for (int r = 0; r < R; r++) {
            npoptot += npop[r];
        }

        double lG;
        if (lGN != null && !Double.isNaN(lGN) && !Double.isInfinite(lGN)) {
            lG = lGN;
        } else {
            lG = logNormConst(L, N, isInf);
        }

        // Infeasible occupancies are not an error: the caller sweeps a lattice.
        if (ntot != npoptot) {
            return new Ret_jointmarg(0.0, Double.NEGATIVE_INFINITY, lG);
        }
        if (npoptot == 0) {
            return new Ret_jointmarg(Math.exp(-lG), -lG, lG);
        }

        // The expanded matrix is square of order sum(N). "spm" works on the
        // unexpanded form, and building this would throw away the very property
        // that makes it independent of the population.
        Matrix A = "spm".equals(eng) ? null : replicate(L, npop, nvec);

        if (!"exact".equals(eng)) {
            int[] zero = firstZero(L, npop, nvec);
            if (zero != null) {
                throw new IllegalArgumentException("pfqn_jointmarg: the '" + eng
                        + "' permanent engine cannot be applied: the demand of class " + (zero[1] + 1)
                        + " at station " + (zero[0] + 1)
                        + " is zero, so the replicated matrix has no full support. Use engine 'exact'.");
            }
        }

        double F;
        if ("exact".equals(eng)) {
            F = Pfqn_perm.pfqn_perm(A);
        } else if ("spm".equals(eng)) {
            // Never the expanded A: the saddle point is asymptotic in the column
            // multiplicities, which are the class populations themselves.
            int kept = 0;
            for (int r = 0; r < R; r++) {
                if (npop[r] > 0) {
                    kept++;
                }
            }
            int[] mult = new int[kept];
            int at = 0;
            for (int r = 0; r < R; r++) {
                if (npop[r] > 0) {
                    mult[at++] = npop[r];
                }
            }
            F = new SaddlePointPermanent(replicateRows(L, npop, nvec), mult, true).getValue();
        } else if ("bethe".equals(eng)) {
            F = new BethePermanent(A, 1e-3, 200000, true).value;
        } else if ("heur".equals(eng)) {
            F = new HeuristicPermanent(A, true).value;
        } else if ("huberlaw".equals(eng)) {
            HuberLawSampler huber = new HuberLawSampler(A);
            huber.solve();
            F = huber.value;
        } else if ("adapart".equals(eng)) {
            AdaPartSampler sampler = new AdaPartSampler(A);
            sampler.solve();
            F = sampler.value;
        } else {
            throw new IllegalArgumentException("pfqn_jointmarg: unrecognized permanent engine '" + engine
                    + "'. Use exact, spm, bethe, heur, huberlaw or adapart.");
        }

        if (F <= 0.0) {
            return new Ret_jointmarg(0.0, Double.NEGATIVE_INFINITY, lG);
        }

        double lp = Math.log(F) - lG;
        for (int r = 0; r < R; r++) {
            lp -= Maths.factln(npop[r]);
        }
        for (int i = 0; i < M; i++) {
            if (isInf[i]) {
                lp -= Maths.factln(nvec[i]);
            }
        }
        return new Ret_jointmarg(Math.exp(lp), lp, lG);
    }

    public static Ret_jointmarg pfqn_jointmarg(Matrix n, Matrix L, Matrix N, int[] infset, Double lGN) {
        return pfqn_jointmarg(n, L, N, infset, lGN, "exact");
    }

    public static Ret_jointmarg pfqn_jointmarg(Matrix n, Matrix L, Matrix N, int[] infset) {
        return pfqn_jointmarg(n, L, N, infset, null, "exact");
    }

    /**
     * Log normalizing constant with the infinite-server rows aggregated into
     * the think time. That aggregation is exact: the delay stations combine by
     * the multinomial theorem, so G does not depend on how they are split.
     */
    private static double logNormConst(Matrix L, Matrix N, boolean[] isInf) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        List<Integer> queues = new ArrayList<Integer>();
        Matrix Z = Matrix.zeros(1, R);
        for (int i = 0; i < M; i++) {
            if (isInf[i]) {
                for (int r = 0; r < R; r++) {
                    Z.set(0, r, Z.get(0, r) + L.get(i, r));
                }
            } else {
                queues.add(Integer.valueOf(i));
            }
        }
        Matrix Lq = new Matrix(queues.size(), R);
        for (int k = 0; k < queues.size(); k++) {
            int i = queues.get(k).intValue();
            for (int r = 0; r < R; r++) {
                Lq.set(k, r, L.get(i, r));
            }
        }
        Ret.pfqnNc res = Pfqn_ca.pfqn_ca(Lq, N, Z);
        return res.lG;
    }

    /**
     * Column r of L repeated N(r) times, then row i of that repeated n(i)
     * times. A station holding no jobs and a class holding no jobs each drop
     * out here, which is what makes a zero entry of the occupancy vector free
     * of any special case: the result stays square of order sum(N).
     */
    private static Matrix replicate(Matrix L, int[] npop, int[] nvec) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        Matrix Ak = null;
        for (int r = 0; r < R; r++) {
            if (npop[r] <= 0) {
                continue;
            }
            Matrix col = Matrix.extractColumn(L, r, null);
            Matrix rep = col.repmat(1, npop[r]);
            Ak = (Ak == null) ? rep : Matrix.concatColumns(Ak, rep, null);
        }
        if (Ak == null) {
            return null;
        }
        Matrix A = null;
        for (int i = 0; i < M; i++) {
            if (nvec[i] <= 0) {
                continue;
            }
            Matrix row = Matrix.extractRows(Ak, i, i + 1, null);
            Matrix rep = row.repmat(nvec[i], 1);
            A = (A == null) ? rep : Matrix.concatRows(A, rep, null);
        }
        return A;
    }

    /**
     * Row i of L repeated n(i) times, keeping only the classes that hold jobs.
     *
     * The same matrix replicate() expands, one step earlier: perm(Ar, m) equals
     * perm(A), and the saddle point wants the unexpanded form because its
     * expansion is asymptotic in m. A class with no jobs is dropped rather than
     * passed with multiplicity zero, so a zero demand in such a column cannot
     * trip the full-support check.
     */
    private static Matrix replicateRows(Matrix L, int[] npop, int[] nvec) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        List<Integer> keepc = new ArrayList<Integer>();
        for (int r = 0; r < R; r++) {
            if (npop[r] > 0) {
                keepc.add(r);
            }
        }
        int rows = 0;
        for (int i = 0; i < M; i++) {
            if (nvec[i] > 0) {
                rows += nvec[i];
            }
        }
        Matrix Ar = new Matrix(rows, keepc.size());
        int at = 0;
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < nvec[i]; c++) {
                for (int l = 0; l < keepc.size(); l++) {
                    Ar.set(at, l, L.get(i, keepc.get(l)));
                }
                at++;
            }
        }
        return Ar;
    }

    /**
     * First (station, class) whose zero demand actually reaches the replicated
     * matrix. A class with no jobs or a station with no jobs contributes
     * nothing, so its zeros are irrelevant.
     */
    private static int[] firstZero(Matrix L, int[] npop, int[] nvec) {
        for (int i = 0; i < L.getNumRows(); i++) {
            if (nvec[i] == 0) {
                continue;
            }
            for (int r = 0; r < L.getNumCols(); r++) {
                if (npop[r] > 0 && L.get(i, r) <= 0.0) {
                    return new int[]{i, r};
                }
            }
        }
        return null;
    }

    private static int[] toIntVector(Matrix v) {
        int len = v.getNumRows() * v.getNumCols();
        int[] out = new int[len];
        int k = 0;
        for (int i = 0; i < v.getNumRows(); i++) {
            for (int j = 0; j < v.getNumCols(); j++) {
                out[k++] = (int) Math.round(v.get(i, j));
            }
        }
        return out;
    }
}
