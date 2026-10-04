/**
 * QRF BAS-blocking NLP bounds: qrf.bas.mmi, qrf.bas.mem and qrf.bas.bethe.
 *
 * @since LINE 3.0
 */
package jline.api.mapqn;

import java.util.HashMap;
import java.util.Map;

/**
 * Minimum-mutual-information, maximum-entropy and tree-reweighted (Bethe)
 * bounds on the BAS-blocking polytope, the JAR twin of MATLAB
 * {@code qrf_bas_mmi.m} / {@code qrf_bas_mem.m} / {@code qrf_bas_bethe.m} and of
 * native Python {@code api/mapqn/qrf_bas_nlp.py}.
 *
 * <p>The polytope is NOT re-derived here. It is the one
 * {@link Mapqn_qr_bounds_bas#buildSystem} assembles and that the LP token
 * {@code qrf.bas} is validated on against the AMPL model
 * {@code qrboundsbas_skel.mod}; only the objective differs. Re-transcribing the
 * fifteen families would duplicate the very index conventions that twice went
 * wrong in the MATLAB twin: the THM30/THM3 population-to-index shift, and a
 * MARGINALS sum whose upper limit was taken on the 1-based index rather than
 * the population.
 */
public final class Mapqn_qrf_bas_nlp {
    private Mapqn_qrf_bas_nlp() {}

    /** Shift keeping the logarithms finite, as in every QRF objective. */
    private static final double LOGTOL = 1.0e-6;

    private static final int KIND_MMI = 0;
    private static final int KIND_MEM = 1;
    private static final int KIND_BETHE = 2;

    /**
     * Minimise mutual information over the BAS polytope.
     *
     * @return a solution carrying {@code UN_i} and {@code QN_i}, 1-based
     */
    public static Mapqn_solution solveMmi(Mapqn_qr_bounds_bas_parameters params) {
        return solve(params, KIND_MMI);
    }

    /** Minimise negative entropy (maximum entropy) over the BAS polytope. */
    public static Mapqn_solution solveMem(Mapqn_qr_bounds_bas_parameters params) {
        return solve(params, KIND_MEM);
    }

    /**
     * Minimise the tree-reweighted (Bethe) free entropy over the BAS polytope.
     *
     * <p>With {@code lambda = 1/M} the objective is
     * {@code lambda*sum_{i!=j} I(n_i;n_j) - sum_i H(n_i)}, i.e. lambda times the
     * MMI body plus the MEM body, both taken over populations from ZERO. lambda
     * is half the uniform point {@code 2/M} of the spanning-tree polytope of
     * K_M, the largest uniform edge weight at which the tree-reweighted entropy
     * is concave and the program therefore convex. This is the objective of
     * {@code qrf.bethe} evaluated over the BAS decision vector, so the blocking
     * configurations and the per-station capacities enter through the ranges
     * alone.
     */
    public static Mapqn_solution solveBethe(Mapqn_qr_bounds_bas_parameters params) {
        return solve(params, KIND_BETHE);
    }

    private static Mapqn_solution solve(Mapqn_qr_bounds_bas_parameters params, int kind) {
        final Mapqn_qr_bounds_bas.BasSystem sys = Mapqn_qr_bounds_bas.buildSystem(params);
        final int n = sys.nCols;
        double[][] Aeq = sys.equalityMatrix();
        double[] beq = sys.equalityRhs();
        double[][] Aub = sys.inequalityMatrix();
        double[] bub = sys.inequalityRhs();
        if (Aub.length == 0) {
            Aub = null;
            bub = null;
        }

        double[] lb = new double[n];
        double[] ub = new double[n];
        for (int i = 0; i < n; i++) ub[i] = 1.0;

        double[] x0 = Mapqn_nlp_solver.feasibleStart(Aeq, beq, Aub, bub, n);
        if (x0 == null) {
            throw new IllegalStateException(
                    "The BAS polytope is infeasible: no phase-1 point was found. Check the "
                            + "blocking parameters (F, BB, MM, ZZ, ZM).");
        }

        // The MI block spans 0..F, the AMPL source's own range; MEM's spans
        // 1..F, as ITS AMPL statement does, except under BETHE, which combines
        // the two and needs one range across both.
        final int[][] mi = (kind == KIND_MEM) ? null : mmiTerms(sys, params, 0);
        final int[][] diag =
                (kind == KIND_MMI) ? null : memTerms(sys, params, (kind == KIND_BETHE) ? 0 : 1);
        if ((mi != null && mi[0].length == 0) || (diag != null && diag[0].length == 0)) {
            throw new IllegalStateException(
                    "The objective is empty: the model registered no joint variables at "
                            + "population one or above.");
        }

        Mapqn_nlp_solver.ObjectiveFn objective;
        Mapqn_nlp_solver.GradientFn gradient;
        if (kind == KIND_MMI) {
            objective = mmiObjective(mi);
            gradient = mmiGradient(mi);
        } else if (kind == KIND_MEM) {
            objective = memObjective(diag);
            gradient = memGradient(diag);
        } else {
            final double lam = 1.0 / params.M;
            objective = betheObjective(mi, diag, lam);
            gradient = betheGradient(mi, diag, lam);
        }
        double[] x = Mapqn_nlp_solver.solve(objective, gradient, n, Aeq, beq, Aub, bub, lb, ub, x0);
        return metrics(x, sys, params);
    }

    /**
     * Column triples (ij, ii, jj) of the MI objective, i != j and ni,nj &gt;= nFrom.
     *
     * <p>nFrom is 0 for both callers: that is the AMPL source's range
     * ({@code sum {nj in 0..F[j]} sum {ni in 0..F[i]}}) and MATLAB's. This port
     * used 1 until 2026-09-03, which dropped the idle cells and made
     * {@code qrf.bas.mmi} disagree with MATLAB.
     */
    private static int[][] mmiTerms(Mapqn_qr_bounds_bas.BasSystem sys,
                                    Mapqn_qr_bounds_bas_parameters params, int nFrom) {
        int M = params.M;
        int[] K = params.K;
        int[] F = params.F;
        int MR = params.MR;
        java.util.List<int[]> triples = new java.util.ArrayList<int[]>();
        for (int m = 0; m < MR; m++) {
            for (int i = 0; i < M; i++) {
                for (int ki = 0; ki < K[i]; ki++) {
                    for (int j = 0; j < M; j++) {
                        if (i == j) continue;
                        for (int kj = 0; kj < K[j]; kj++) {
                            for (int ni = nFrom; ni <= F[i]; ni++) {
                                int cii = sys.p2Column(i, ni, ki, i, ni, ki, m);
                                if (cii < 0) continue;
                                for (int nj = nFrom; nj <= F[j]; nj++) {
                                    int cij = sys.p2Column(i, ni, ki, j, nj, kj, m);
                                    int cjj = sys.p2Column(j, nj, kj, j, nj, kj, m);
                                    if (cij < 0 || cjj < 0) continue;
                                    triples.add(new int[]{cij, cii, cjj});
                                }
                            }
                        }
                    }
                }
            }
        }
        int[][] out = new int[3][triples.size()];
        for (int t = 0; t < triples.size(); t++) {
            int[] tr = triples.get(t);
            out[0][t] = tr[0];
            out[1][t] = tr[1];
            out[2][t] = tr[2];
        }
        return out;
    }

    /** Diagonal columns of the MEM objective, ni &gt;= nFrom (1 for MEM, 0 for BETHE). */
    private static int[][] memTerms(Mapqn_qr_bounds_bas.BasSystem sys,
                                    Mapqn_qr_bounds_bas_parameters params, int nFrom) {
        int M = params.M;
        int[] K = params.K;
        int[] F = params.F;
        int MR = params.MR;
        java.util.List<Integer> cols = new java.util.ArrayList<Integer>();
        for (int m = 0; m < MR; m++) {
            for (int i = 0; i < M; i++) {
                for (int k = 0; k < K[i]; k++) {
                    for (int ni = nFrom; ni <= F[i]; ni++) {
                        int c = sys.p2Column(i, ni, k, i, ni, k, m);
                        if (c >= 0) cols.add(Integer.valueOf(c));
                    }
                }
            }
        }
        int[][] out = new int[1][cols.size()];
        for (int t = 0; t < cols.size(); t++) out[0][t] = cols.get(t).intValue();
        return out;
    }

    private static Mapqn_nlp_solver.ObjectiveFn mmiObjective(final int[][] terms) {
        return new Mapqn_nlp_solver.ObjectiveFn() {
            @Override
            public double apply(double[] x) {
                double f = 0.0;
                for (int t = 0; t < terms[0].length; t++) {
                    double pij = x[terms[0][t]];
                    double pii = x[terms[1][t]];
                    double pjj = x[terms[2][t]];
                    f += pij * (Math.log(LOGTOL + pij) - Math.log(LOGTOL + pii)
                            - Math.log(LOGTOL + pjj));
                }
                return f;
            }
        };
    }

    private static Mapqn_nlp_solver.GradientFn mmiGradient(final int[][] terms) {
        return new Mapqn_nlp_solver.GradientFn() {
            @Override
            public void apply(double[] x, double[] gradOut) {
                for (int t = 0; t < terms[0].length; t++) {
                    int aij = terms[0][t];
                    int aii = terms[1][t];
                    int ajj = terms[2][t];
                    double pij = x[aij];
                    double pii = x[aii];
                    double pjj = x[ajj];
                    gradOut[aij] += Math.log(LOGTOL + pij) - Math.log(LOGTOL + pii)
                            - Math.log(LOGTOL + pjj) + pij / (LOGTOL + pij);
                    gradOut[aii] -= pij / (LOGTOL + pii);
                    gradOut[ajj] -= pij / (LOGTOL + pjj);
                }
            }
        };
    }

    private static Mapqn_nlp_solver.ObjectiveFn memObjective(final int[][] terms) {
        return new Mapqn_nlp_solver.ObjectiveFn() {
            @Override
            public double apply(double[] x) {
                double f = 0.0;
                for (int t = 0; t < terms[0].length; t++) {
                    double p = x[terms[0][t]];
                    f += p * Math.log(LOGTOL + p);
                }
                return f;
            }
        };
    }

    private static Mapqn_nlp_solver.GradientFn memGradient(final int[][] terms) {
        return new Mapqn_nlp_solver.GradientFn() {
            @Override
            public void apply(double[] x, double[] gradOut) {
                for (int t = 0; t < terms[0].length; t++) {
                    int c = terms[0][t];
                    double p = x[c];
                    gradOut[c] += Math.log(LOGTOL + p) + p / (LOGTOL + p);
                }
            }
        };
    }

    private static Mapqn_nlp_solver.ObjectiveFn betheObjective(final int[][] mi,
                                                              final int[][] diag,
                                                              final double lam) {
        return new Mapqn_nlp_solver.ObjectiveFn() {
            @Override
            public double apply(double[] x) {
                double f = 0.0;
                for (int t = 0; t < mi[0].length; t++) {
                    double pij = x[mi[0][t]];
                    double pii = x[mi[1][t]];
                    double pjj = x[mi[2][t]];
                    f += lam * pij * (Math.log(LOGTOL + pij) - Math.log(LOGTOL + pii)
                            - Math.log(LOGTOL + pjj));
                }
                for (int t = 0; t < diag[0].length; t++) {
                    double p = x[diag[0][t]];
                    f += p * Math.log(LOGTOL + p);
                }
                return f;
            }
        };
    }

    private static Mapqn_nlp_solver.GradientFn betheGradient(final int[][] mi,
                                                             final int[][] diag,
                                                             final double lam) {
        return new Mapqn_nlp_solver.GradientFn() {
            @Override
            public void apply(double[] x, double[] gradOut) {
                for (int t = 0; t < mi[0].length; t++) {
                    int aij = mi[0][t];
                    int aii = mi[1][t];
                    int ajj = mi[2][t];
                    double pij = x[aij];
                    double pii = x[aii];
                    double pjj = x[ajj];
                    gradOut[aij] += lam * (Math.log(LOGTOL + pij) - Math.log(LOGTOL + pii)
                            - Math.log(LOGTOL + pjj) + pij / (LOGTOL + pij));
                    gradOut[aii] -= lam * pij / (LOGTOL + pii);
                    gradOut[ajj] -= lam * pij / (LOGTOL + pjj);
                }
                for (int t = 0; t < diag[0].length; t++) {
                    int c = diag[0][t];
                    double p = x[c];
                    gradOut[c] += Math.log(LOGTOL + p) + p / (LOGTOL + p);
                }
            }
        };
    }

    /**
     * Utilization and queue length at the optimal point.
     *
     * UTILIZATION, not occupancy. UN used to sum the diagonal p2 over ALL
     * blocking configurations, i.e. P(n_i &gt;= 1) with the BLOCKED ones
     * included. A blocked BAS server holds a job it has already finished and
     * does no work, so that is occupancy: on the M=2, N=3, F=[2 3] cyclic model
     * it reported U2 = 1 where the exact utilization is 7/15, which the LP over
     * the SAME polytope already returns. The e variables carry the right
     * quantity -- UEFF pins e(i,ki) to the mass with n_i &gt;= 1 in the
     * configurations where i is NOT blocked.
     *
     * NO 1/M, matching the LP builder this shares: {@code Mapqn_qr_bounds_bas}
     * emits UEFF as one row per (j, i, ki), which leaves e unscaled. QN stays on
     * the diagonal p2 over every configuration, blocked included, because a
     * blocked job is still held at the station and counts towards its
     * population.
     */
    private static Mapqn_solution metrics(double[] x, Mapqn_qr_bounds_bas.BasSystem sys,
                                          Mapqn_qr_bounds_bas_parameters params) {
        int M = params.M;
        int[] K = params.K;
        int[] F = params.F;
        int MR = params.MR;
        Map<String, Double> vars = new HashMap<String, Double>();
        for (int i = 0; i < M; i++) {
            double un = 0.0;
            double qn = 0.0;
            for (int ki = 0; ki < K[i]; ki++) {
                int ec = sys.eColumn(i, ki);
                if (ec >= 0) un += x[ec];
            }
            for (int m = 0; m < MR; m++) {
                for (int ni = 1; ni <= F[i]; ni++) {
                    for (int ki = 0; ki < K[i]; ki++) {
                        int c = sys.p2Column(i, ni, ki, i, ni, ki, m);
                        if (c < 0) continue;
                        qn += ni * x[c];
                    }
                }
            }
            vars.put("UN_" + (i + 1), Double.valueOf(un));
            vars.put("QN_" + (i + 1), Double.valueOf(qn));
        }
        return new Mapqn_solution(0.0, vars);
    }
}
