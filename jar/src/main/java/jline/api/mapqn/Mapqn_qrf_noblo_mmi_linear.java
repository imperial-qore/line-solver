/**
 * QRF No-Blocking NLP with Linear Constraint Matrices.
 * Port of MATLAB qrf_noblo_mmi_linear.m.
 * The "linear" is about HOW the constraints are built (explicit sparse Aeq/beq
 * rather than a residual callback), not about the objective, which is MMI. Until
 * 2026-08-29 the MATLAB reference called its own mem() here and this port
 * mirrored that.
 * @since LINE 3.0
 */
package jline.api.mapqn;

import java.util.List;

public final class Mapqn_qrf_noblo_mmi_linear {
    private Mapqn_qrf_noblo_mmi_linear() {}


    public static Mapqn_solution solve(double[][][][] MAPs, int N, double[][] rt, double[][] alpha) {
        int M = MAPs.length;
        int[] K = new int[M];
        for (int i = 0; i < M; i++) {
            K[i] = MAPs[i][0].length;
        }
        int Kmax = 0;
        for (int i = 0; i < M; i++) {
            if (K[i] > Kmax) Kmax = K[i];
        }

        double[][] alphaEff;
        if (alpha != null) {
            alphaEff = alpha;
        } else {
            alphaEff = new double[M][N];
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < N; j++) alphaEff[i][j] = 1.0;
            }
        }

        double[][][] mu = new double[M][Kmax][Kmax];
        double[][][] v = new double[M][Kmax][Kmax];
        for (int i = 0; i < M; i++) {
            double[][] D0 = MAPs[i][0];
            double[][] D1 = MAPs[i][1];
            for (int h = 0; h < K[i]; h++) {
                for (int k = 0; k < K[i]; k++) {
                    mu[i][h][k] = D1[h][k];
                    v[i][h][k] = (h == k) ? 0.0 : D0[h][k]; // (from, to), as mu is
                }
            }
        }

        int MR = 1;
        double[][] BB = new double[1][M];
        int[] F = new int[M];
        for (int i = 0; i < M; i++) F[i] = N;

        // see _kb/03-api-layer.md for rationale
        double[][][][][] q = new double[M][M][Kmax][Kmax][N + 1];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                for (int k = 0; k < K[i]; k++) {
                    for (int h = 0; h < K[i]; h++) {
                        for (int n = 1; n <= N; n++) {
                            double a = alphaEff[i][n - 1];
                            q[i][j][k][h][n] = (j != i)
                                    ? rt[i][j] * mu[i][k][h] * a
                                    : v[i][k][h] * a + rt[i][i] * mu[i][k][h] * a;
                        }
                    }
                }
            }
        }

        final int numVars = Mapqn_qrf_noblo_mmi.computeNumVars(M, N, Kmax, MR);
        double[] x0 = new double[numVars];
        for (int i = 0; i < numVars; i++) x0[i] = 1.0;
        double totalSum = 0.0;
        for (int i = 0; i < x0.length; i++) totalSum += x0[i];
        for (int i = 0; i < x0.length; i++) x0[i] /= totalSum;

        double[] lb = new double[numVars];
        double[] ub = new double[numVars];
        for (int i = 0; i < numVars; i++) {
            lb[i] = 0.0;
            ub[i] = 1.0;
        }

        Mapqn_qrf_noblo_mmi.ConstraintSet cs =
                Mapqn_qrf_noblo_mmi.buildConstraintsLd(q, M, MR, BB, F, N, K, Kmax, numVars);

        List<double[]> aeqList = cs.aeq;
        List<Double> beqList = cs.beq;
        List<double[]> aubList = cs.aub;
        List<Double> bubList = cs.bub;

        double[][] Aeq = aeqList.toArray(new double[0][]);
        double[] beq = new double[beqList.size()];
        for (int i = 0; i < beq.length; i++) beq[i] = beqList.get(i);

        double[][] Aub = aubList.isEmpty() ? null : aubList.toArray(new double[0][]);
        double[] bub;
        if (bubList.isEmpty()) {
            bub = null;
        } else {
            bub = new double[bubList.size()];
            for (int i = 0; i < bub.length; i++) bub[i] = bubList.get(i);
        }

        final int Mf = M;
        final int Nf = N;
        final int KmaxF = Kmax;
        final int MRf = MR;
        final int[] Kf = K;
        final int[] Ff = F;

        // MMI objective. The 'linear' in the name is about HOW the constraints are
        // built, not about the objective; until 2026-08-29 the MATLAB reference called
        // its own mem() here and this port mirrored that.
        Mapqn_nlp_solver.ObjectiveFn objective = new Mapqn_nlp_solver.ObjectiveFn() {
            @Override
            public double apply(double[] x) {
                return Mapqn_qrf_noblo_mmi.mmiObjective(x, Mf, Nf, Kf, KmaxF, Ff, MRf);
            }
        };
        Mapqn_nlp_solver.GradientFn gradient = new Mapqn_nlp_solver.GradientFn() {
            @Override
            public void apply(double[] x, double[] gradOut) {
                Mapqn_qrf_noblo_mmi.mmiGradient(x, gradOut, Mf, Nf, Kf, KmaxF, Ff, MRf);
            }
        };

        double[] start = Mapqn_nlp_solver.feasibleStart(Aeq, beq, Aub, bub, numVars);
        if (start != null) x0 = start;
        double[] xOpt = Mapqn_nlp_solver.solve(objective, gradient, numVars,
                Aeq, beq, Aub, bub, lb, ub, x0);
        return Mapqn_qrf_noblo_mmi.extractResults(xOpt, M, N, K, Kmax, F, MR, alphaEff);
    }

    public static Mapqn_solution solve(double[][][][] MAPs, int N, double[][] rt) {
        return solve(MAPs, N, rt, null);
    }
}
