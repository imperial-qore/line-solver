/**
 * QRF No-Blocking BETHE approximation: the tree-reweighted (Bethe) free entropy
 * minimized over the same no-blocking polytope {@code qrf.mmi} uses.
 * Port of MATLAB qrf_noblo_bethe.m.
 * @since LINE 3.0
 */
package jline.api.mapqn;

public final class Mapqn_qrf_noblo_bethe {
    private Mapqn_qrf_noblo_bethe() {}

    /**
     * Solve the no-blocking quadratic reduction under the tree-reweighted free
     * entropy.
     *
     * <p>The polytope, the phase-1 feasible start and the NLP call are exactly
     * those of {@link Mapqn_qrf_noblo_mmi#solve}; the objective is the only
     * difference. See {@link Mapqn_qrf_noblo_mmi#betheObjective} for what it is
     * and why the weight is {@code lambda = 1/M}.
     *
     * <p>ONE SOLVE, NO RESTARTS. The objective is convex on this polytope, so
     * there is no second local minimum for a restart to find; the single solve
     * from the phase-1 point returns the global optimum.
     *
     * @param M  number of queues
     * @param MR_input ignored -- no blocking means one configuration, by definition
     * @param K  phases per queue
     * @param N  total population
     * @param mu completion rates, [M][Kmax][Kmax]
     * @param v  background rates, [M][Kmax][Kmax]
     * @param rt routing matrix, [M][M]
     */
    public static Mapqn_solution solve(int M, int MR_input, int[] K, int N,
                                       double[][][] mu, double[][][] v, double[][] rt) {
        final int MR = 1;
        final double[][] BB = new double[1][M];
        final int[] F = new int[M];
        for (int i = 0; i < M; i++) F[i] = N;

        int Kmax = 0;
        for (int kv : K) if (kv > Kmax) Kmax = kv;
        final int KmaxF = Kmax;

        final double[][][][] q = new double[M][M][Kmax][Kmax];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                for (int k = 0; k < K[i]; k++) {
                    for (int h = 0; h < K[i]; h++) {
                        q[i][j][k][h] = (j != i)
                                ? rt[i][j] * mu[i][k][h]
                                : v[i][k][h] + rt[i][i] * mu[i][k][h];
                    }
                }
            }
        }

        final int numVars = Mapqn_qrf_noblo_mmi.computeNumVars(M, N, Kmax, MR);
        double[] x0 = new double[numVars];
        double[] lb = new double[numVars];
        double[] ub = new double[numVars];
        for (int i = 0; i < numVars; i++) {
            lb[i] = 0.0;
            ub[i] = 1.0;
        }

        Mapqn_qrf_noblo_mmi.ConstraintSet cs =
                Mapqn_qrf_noblo_mmi.buildConstraints(q, M, MR, BB, F, N, K, Kmax, numVars);

        double[][] Aeq = cs.aeq.toArray(new double[0][]);
        double[] beq = new double[cs.beq.size()];
        for (int i = 0; i < cs.beq.size(); i++) beq[i] = cs.beq.get(i);

        double[][] Aub = cs.aub.isEmpty() ? null : cs.aub.toArray(new double[0][]);
        double[] bub = null;
        if (!cs.bub.isEmpty()) {
            bub = new double[cs.bub.size()];
            for (int i = 0; i < cs.bub.size(); i++) bub[i] = cs.bub.get(i);
        }

        final int Mf = M;
        final int Nf = N;
        final int MRf = MR;
        final int[] Kf = K;
        final int[] Ff = F;

        Mapqn_nlp_solver.ObjectiveFn objective = new Mapqn_nlp_solver.ObjectiveFn() {
            public double apply(double[] x) {
                return Mapqn_qrf_noblo_mmi.betheObjective(x, Mf, Nf, Kf, KmaxF, Ff, MRf);
            }
        };
        Mapqn_nlp_solver.GradientFn gradient = new Mapqn_nlp_solver.GradientFn() {
            public void apply(double[] x, double[] gradOut) {
                Mapqn_qrf_noblo_mmi.betheGradient(x, gradOut, Mf, Nf, Kf, KmaxF, Ff, MRf);
            }
        };

        double[] start = Mapqn_nlp_solver.feasibleStart(Aeq, beq, Aub, bub, numVars);
        if (start != null) x0 = start;
        double[] xOpt = Mapqn_nlp_solver.solve(objective, gradient, numVars,
                Aeq, beq, Aub, bub, lb, ub, x0);

        return Mapqn_qrf_noblo_mmi.extractResults(xOpt, M, N, K, Kmax, F, MR);
    }
}
