/**
 * QRF No-Blocking MEM (Maximum Entropy Method) Approximation.
 * Port of MATLAB qrf_noblo_mem.m.
 * @since LINE 3.0
 */
package jline.api.mapqn;

import java.util.List;

public final class Mapqn_qrf_noblo_mem {
    private Mapqn_qrf_noblo_mem() {}


    public static Mapqn_solution solve(double[][][][] MAPs, int N, double[][] rt) {
        int M = MAPs.length;
        int[] K = new int[M];
        for (int i = 0; i < M; i++) {
            K[i] = MAPs[i][0].length;
        }
        int Kmax = 0;
        for (int i = 0; i < M; i++) {
            if (K[i] > Kmax) Kmax = K[i];
        }

        // Extract mu, v from MAPs
        final double[][][] mu = new double[M][Kmax][Kmax];
        final double[][][] v = new double[M][Kmax][Kmax];
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

        // Build q
        double[][][][] q = new double[M][M][Kmax][Kmax];
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
            x0[i] = 0.0;
            lb[i] = 0.0;
            ub[i] = 1.0;
        }

        Mapqn_qrf_noblo_mmi.ConstraintSet cs =
                Mapqn_qrf_noblo_mmi.buildConstraints(q, M, MR, BB, F, N, K, Kmax, numVars);

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

        // MEM objective: memObjective returns the NEGATIVE entropy, so minimizing it
        // maximizes H, which is what the AMPL model asks for.
        Mapqn_nlp_solver.ObjectiveFn objective = new Mapqn_nlp_solver.ObjectiveFn() {
            @Override
            public double apply(double[] x) {
                return Mapqn_qrf_noblo_mmi.memObjective(x, Mf, Nf, Kf, KmaxF, Ff, MRf);
            }
        };
        Mapqn_nlp_solver.GradientFn gradient = new Mapqn_nlp_solver.GradientFn() {
            @Override
            public void apply(double[] x, double[] gradOut) {
                Mapqn_qrf_noblo_mmi.memGradient(x, gradOut, Mf, Nf, Kf, KmaxF, Ff, MRf);
            }
        };

        double[] start = Mapqn_nlp_solver.feasibleStart(Aeq, beq, Aub, bub, numVars);
        if (start != null) x0 = start;
        double[] xOpt = Mapqn_nlp_solver.solve(objective, gradient, numVars,
                Aeq, beq, Aub, bub, lb, ub, x0);
        return Mapqn_qrf_noblo_mmi.extractResults(xOpt, M, N, K, Kmax, F, MR);
    }
}
