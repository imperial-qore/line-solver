/**
 * QRF No-Blocking MMI (Minimize Mutual Information) Approximation.
 * Port of MATLAB qrf_noblo_mmi.m.
 * @since LINE 3.0
 */
package jline.api.mapqn;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public final class Mapqn_qrf_noblo_mmi {
    private Mapqn_qrf_noblo_mmi() {}

    private static final double LOGTOL = 1e-6;

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

        final int numVars = computeNumVars(M, N, Kmax, MR);
        double[] x0 = new double[numVars];
        double[] lb = new double[numVars];
        double[] ub = new double[numVars];
        for (int i = 0; i < numVars; i++) {
            lb[i] = 0.0;
            ub[i] = 1.0;
        }

        ConstraintSet cs = buildConstraints(q, M, MR, BB, F, N, K, Kmax, numVars);

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
                Double[][][][][][][] p2 = unflattenP2(x, Mf, Nf, KmaxF, MRf);
                double fobj = 0.0;
                for (int m = 0; m < MRf; m++) {
                    for (int i = 0; i < Mf; i++) {
                        for (int ki = 0; ki < Kf[i]; ki++) {
                            for (int j = 0; j < Mf; j++) {
                                if (i != j) {
                                    for (int kj = 0; kj < Kf[j]; kj++) {
                                        for (int ni = 1; ni <= Ff[i]; ni++) {
                                            for (int nj = 1; nj <= Ff[j]; nj++) {
                                                double pij = p2[i][ni][ki][j][nj][kj][m];
                                                double pii = p2[i][ni][ki][i][ni][ki][m];
                                                double pjj = p2[j][nj][kj][j][nj][kj][m];
                                                fobj += pij * (Math.log(LOGTOL + pij)
                                                        - Math.log(LOGTOL + pii) - Math.log(LOGTOL + pjj));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                return fobj;
            }
        };

        double[] xOpt = Mapqn_nlp_solver.solve(objective, numVars, Aeq, beq, Aub, bub, lb, ub, x0);

        return extractResults(xOpt, M, N, K, Kmax, F, MR);
    }

    public static int computeNumVars(int M, int N, int Kmax, int MR) {
        return M * (N + 1) * Kmax * M * (N + 1) * Kmax * MR + M * Kmax;
    }

    public static int p2Index(int j, int nj, int k, int i, int ni, int h, int m,
                              int M, int N, int Kmax, int MR) {
        int idx = (N + 1) * Kmax * M * (N + 1) * Kmax * MR * j;
        idx += Kmax * M * (N + 1) * Kmax * MR * nj;
        idx += M * (N + 1) * Kmax * MR * k;
        idx += (N + 1) * Kmax * MR * i;
        idx += Kmax * MR * ni;
        idx += MR * h;
        idx += m;
        return idx;
    }

    public static int eIndex(int i, int k, int M, int N, int Kmax, int MR) {
        int p2Size = M * (N + 1) * Kmax * M * (N + 1) * Kmax * MR;
        return p2Size + Kmax * i + k;
    }

    public static Double[][][][][][][] unflattenP2(double[] x, int M, int N, int Kmax, int MR) {
        Double[][][][][][][] p2 = new Double[M][N + 1][Kmax][M][N + 1][Kmax][MR];
        int ctr = 0;
        for (int j = 0; j < M; j++) {
            for (int nj = 0; nj <= N; nj++) {
                for (int k = 0; k < Kmax; k++) {
                    for (int i = 0; i < M; i++) {
                        for (int ni = 0; ni <= N; ni++) {
                            for (int h = 0; h < Kmax; h++) {
                                for (int m = 0; m < MR; m++) {
                                    p2[j][nj][k][i][ni][h][m] = Double.valueOf(x[ctr++]);
                                }
                            }
                        }
                    }
                }
            }
        }
        return p2;
    }

    public static Mapqn_solution extractResults(double[] xOpt, int M, int N, int[] K, int Kmax,
                                                int[] F, int MR) {
        Double[][][][][][][] p2 = unflattenP2(xOpt, M, N, Kmax, MR);
        double[] UN = new double[M];
        double[] QN = new double[M];
        for (int ti = 0; ti < M; ti++) {
            for (int m = 0; m < MR; m++) {
                for (int ni = 1; ni <= F[ti]; ni++) {
                    for (int ki = 0; ki < K[ti]; ki++) {
                        UN[ti] += p2[ti][ni][ki][ti][ni][ki][m];
                        QN[ti] += ni * p2[ti][ni][ki][ti][ni][ki][m];
                    }
                }
            }
        }
        Map<String, Double> vars = new HashMap<String, Double>();
        for (int i = 0; i < M; i++) {
            vars.put("UN_" + (i + 1), Double.valueOf(UN[i]));
            vars.put("QN_" + (i + 1), Double.valueOf(QN[i]));
        }
        return new Mapqn_solution(0.0, vars);
    }

    public static final class ConstraintSet {
        public final List<double[]> aeq;
        public final List<Double> beq;
        public final List<double[]> aub;
        public final List<Double> bub;

        public ConstraintSet(List<double[]> aeq, List<Double> beq,
                             List<double[]> aub, List<Double> bub) {
            this.aeq = aeq;
            this.beq = beq;
            this.aub = aub;
            this.bub = bub;
        }
    }

    /**
     * Population-free form, matching the q of qrboundsbas_skel.mod. Used by the
     * qrf_noblo_mmi / qrf_noblo_mem variants, which have no load dependence.
     */
    public static ConstraintSet buildConstraints(double[][][][] q, int M, int MR, double[][] BB,
                                                 int[] F, int N, int[] K, int Kmax, int numVars) {
        return buildConstraintsCore(q, null, M, MR, BB, F, N, K, Kmax, numVars);
    }

    /**
     * Load-dependent form, matching the five-index q of qrboundsrsrd_skel.mod:11
     * whose last index is the population of the EMITTING station. A 4D q cannot
     * represent load dependence: the balance families would be unable to tell
     * the rate at which station i empties at population n from the rate at n'.
     */
    public static ConstraintSet buildConstraintsLd(double[][][][][] q, int M, int MR, double[][] BB,
                                                   int[] F, int N, int[] K, int Kmax, int numVars) {
        return buildConstraintsCore(null, q, M, MR, BB, F, N, K, Kmax, numVars);
    }

    private static double qAt(double[][][][] q4, double[][][][][] q5,
                              int i, int j, int k, int h, int n) {
        return (q5 != null) ? q5[i][j][k][h][n] : q4[i][j][k][h];
    }

    private static ConstraintSet buildConstraintsCore(double[][][][] q4, double[][][][][] q5,
                                                      int M, int MR, double[][] BB,
                                                      int[] F, int N, int[] K, int Kmax,
                                                      int numVars) {
        final boolean loadDependent = (q5 != null);
        List<double[]> aeq = new ArrayList<double[]>();
        List<Double> beq = new ArrayList<Double>();
        List<double[]> aub = new ArrayList<double[]>();
        List<Double> bub = new ArrayList<Double>();

        // ONE: normalization
        for (int j = 0; j < M; j++) {
            double[] row = new double[numVars];
            for (int nj = 0; nj <= N; nj++) {
                for (int k = 0; k < K[j]; k++) {
                    for (int m = 0; m < MR; m++) {
                        row[p2Index(j, nj, k, j, nj, k, m, M, N, Kmax, MR)] += 1.0;
                    }
                }
            }
            aeq.add(row);
            beq.add(Double.valueOf(1.0));
        }

        // ZERO1: i==j and ni==nj and h!=k
        for (int j = 0; j < M; j++)
            for (int k = 0; k < K[j]; k++)
                for (int nj = 0; nj <= N; nj++)
                    for (int i = 0; i < M; i++)
                        for (int h = 0; h < K[i]; h++)
                            for (int ni = 0; ni <= N; ni++)
                                for (int m = 0; m < MR; m++) {
                                    if (i == j && nj == ni && h != k) {
                                        double[] row = new double[numVars];
                                        row[p2Index(j, nj, k, i, ni, h, m, M, N, Kmax, MR)] = 1.0;
                                        aeq.add(row);
                                        beq.add(Double.valueOf(0.0));
                                    }
                                }

        // ZERO2: i==j and nj!=ni
        for (int j = 0; j < M; j++)
            for (int k = 0; k < K[j]; k++)
                for (int nj = 0; nj <= N; nj++)
                    for (int i = 0; i < M; i++)
                        for (int h = 0; h < K[i]; h++)
                            for (int ni = 0; ni <= N; ni++)
                                for (int m = 0; m < MR; m++) {
                                    if (i == j && nj != ni) {
                                        double[] row = new double[numVars];
                                        row[p2Index(j, nj, k, i, ni, h, m, M, N, Kmax, MR)] = 1.0;
                                        aeq.add(row);
                                        beq.add(Double.valueOf(0.0));
                                    }
                                }

        // ZERO3: i!=j and nj+ni>N
        for (int j = 0; j < M; j++)
            for (int k = 0; k < K[j]; k++)
                for (int nj = 0; nj <= N; nj++)
                    for (int i = 0; i < M; i++)
                        for (int h = 0; h < K[i]; h++)
                            for (int ni = 0; ni <= N; ni++)
                                for (int m = 0; m < MR; m++) {
                                    if (i != j && nj + ni > N) {
                                        double[] row = new double[numVars];
                                        row[p2Index(j, nj, k, i, ni, h, m, M, N, Kmax, MR)] = 1.0;
                                        aeq.add(row);
                                        beq.add(Double.valueOf(0.0));
                                    }
                                }

        // ZERO5: BB[m,j]==1 for m >= 1 (0-based) => p2[j,0,...]=0.
        // Vacuous at MR==1, live for the blocking configurations.
        for (int j = 0; j < M; j++)
            for (int k = 0; k < K[j]; k++)
                for (int i = 0; i < M; i++)
                    for (int h = 0; h < K[i]; h++)
                        for (int ni = 0; ni <= F[i]; ni++)
                            for (int m = 1; m < MR; m++) {
                                if (BB[m][j] == 1.0) {
                                    double[] row = new double[numVars];
                                    row[p2Index(j, 0, k, i, ni, h, m, M, N, Kmax, MR)] = 1.0;
                                    aeq.add(row);
                                    beq.add(Double.valueOf(0.0));
                                }
                            }

        // ZERO6: nj > F[j]
        for (int j = 0; j < M; j++)
            for (int k = 0; k < K[j]; k++)
                for (int nj = F[j] + 1; nj <= N; nj++)
                    for (int i = 0; i < M; i++)
                        for (int h = 0; h < K[i]; h++)
                            for (int ni = 0; ni <= N; ni++)
                                for (int m = 0; m < MR; m++) {
                                    double[] row = new double[numVars];
                                    row[p2Index(j, nj, k, i, ni, h, m, M, N, Kmax, MR)] = 1.0;
                                    aeq.add(row);
                                    beq.add(Double.valueOf(0.0));
                                }

        // SYMMETRY
        for (int j = 0; j < M; j++)
            for (int nj = 0; nj <= N; nj++)
                for (int k = 0; k < K[j]; k++)
                    for (int i = 0; i < M; i++)
                        for (int ni = 0; ni <= N; ni++)
                            for (int h = 0; h < K[i]; h++)
                                for (int m = 0; m < MR; m++) {
                                    double[] row = new double[numVars];
                                    row[p2Index(i, ni, h, j, nj, k, m, M, N, Kmax, MR)] += 1.0;
                                    row[p2Index(j, nj, k, i, ni, h, m, M, N, Kmax, MR)] -= 1.0;
                                    aeq.add(row);
                                    beq.add(Double.valueOf(0.0));
                                }

        // MARGINALS
        for (int j = 0; j < M; j++)
            for (int k = 0; k < K[j]; k++)
                for (int nj = 0; nj <= N; nj++)
                    for (int i = 0; i < M; i++)
                        for (int m = 0; m < MR; m++) {
                            if (i != j) {
                                double[] row = new double[numVars];
                                row[p2Index(j, nj, k, j, nj, k, m, M, N, Kmax, MR)] += 1.0;
                                // see _kb/03-api-layer.md for rationale
                                for (int ni = 0; ni <= N; ni++) {
                                    for (int h = 0; h < K[i]; h++) {
                                        row[p2Index(j, nj, k, i, ni, h, m, M, N, Kmax, MR)] -= 1.0;
                                    }
                                }
                                aeq.add(row);
                                beq.add(Double.valueOf(0.0));
                            }
                        }

        // UEFF
        for (int j = 0; j < M; j++)
            for (int i = 0; i < M; i++)
                for (int ki = 0; ki < K[i]; ki++) {
                    double[] row = new double[numVars];
                    row[eIndex(i, ki, M, N, Kmax, MR)] += 1.0;
                    for (int nj = 0; nj <= N; nj++)
                        for (int kj = 0; kj < K[j]; kj++)
                            for (int m = 0; m < MR; m++)
                                for (int ni = 1; ni <= N; ni++) {
                                    if (BB[m][i] == 0.0) {
                                        row[p2Index(j, nj, kj, i, ni, ki, m, M, N, Kmax, MR)] -= 1.0;
                                    }
                                }
                    aeq.add(row);
                    beq.add(Double.valueOf(0.0));
                }

        // see _kb/03-api-layer.md for rationale
        for (int i = 0; i < M; i++)
            for (int k = 0; k < K[i]; k++) {
                double[] row = new double[numVars];
                if (loadDependent) {
                    for (int ni = 1; ni <= F[i]; ni++)
                        for (int m = 0; m < MR; m++)
                            for (int j = 0; j < M; j++)
                                for (int h = 0; h < K[i]; h++) {
                                    row[p2Index(i, ni, k, i, ni, k, m, M, N, Kmax, MR)]
                                            += q5[i][j][k][h][ni];
                                    row[p2Index(i, ni, h, i, ni, h, m, M, N, Kmax, MR)]
                                            -= q5[i][j][h][k][ni];
                                }
                } else {
                    for (int j = 0; j < M; j++)
                        for (int h = 0; h < K[i]; h++) {
                            row[eIndex(i, k, M, N, Kmax, MR)] += q4[i][j][k][h];
                        }
                    for (int j = 0; j < M; j++)
                        for (int h = 0; h < K[i]; h++) {
                            row[eIndex(i, h, M, N, Kmax, MR)] -= q4[i][j][h][k];
                        }
                }
                aeq.add(row);
                beq.add(Double.valueOf(0.0));
            }

        // THM2
        for (int j = 0; j < M; j++)
            for (int k = 0; k < K[j]; k++)
                for (int nj = 0; nj <= F[j]; nj++)
                    for (int m = 0; m < MR; m++) {
                        double[] row = new double[numVars];
                        for (int i = 0; i < M; i++)
                            for (int ni = 1; ni <= F[i]; ni++)
                                for (int ki = 0; ki < K[i]; ki++) {
                                    row[p2Index(j, nj, k, i, ni, ki, m, M, N, Kmax, MR)] += (double) ni;
                                }
                        row[p2Index(j, nj, k, j, nj, k, m, M, N, Kmax, MR)] -= (double) N;
                        aeq.add(row);
                        beq.add(Double.valueOf(0.0));
                    }

        // COR1
        double[] rowCor = new double[numVars];
        for (int m = 0; m < MR; m++)
            for (int i = 0; i < M; i++)
                for (int j = 0; j < M; j++)
                    for (int nj = 1; nj <= F[j]; nj++)
                        for (int ni = 1; ni <= F[i]; ni++)
                            for (int ki = 0; ki < K[i]; ki++)
                                for (int kj = 0; kj < K[j]; kj++) {
                                    rowCor[p2Index(j, nj, kj, i, ni, ki, m, M, N, Kmax, MR)] += (double) (ni * nj);
                                }
        aeq.add(rowCor);
        beq.add(Double.valueOf((double) (N * N)));

        // see _kb/03-api-layer.md for rationale

        // THM30 {i, u}: balance across the ni = 0 boundary of station i.
        for (int i = 0; i < M; i++)
            for (int u = 0; u < K[i]; u++) {
                double[] row = new double[numVars];
                for (int j = 0; j < M; j++) {
                    if (j == i) continue;
                    for (int nj = 1; nj <= F[j]; nj++)
                        for (int k = 0; k < K[j]; k++) {
                            double coef = 0.0;
                            for (int h = 0; h < K[j]; h++) coef += qAt(q4, q5, j, i, k, h, nj);
                            if (coef == 0.0) continue;
                            for (int m = 0; m < MR; m++) {
                                if (BB[m][j] != 0.0) continue;
                                row[p2Index(j, nj, k, i, 0, u, m, M, N, Kmax, MR)] += coef;
                            }
                        }
                }
                for (int j = 0; j < M; j++) {
                    if (j == i) continue;
                    for (int nj = 0; nj <= F[j]; nj++)
                        for (int k = 0; k < K[i]; k++) {
                            // Population index 1, not 0: the AMPL right-hand side
                            // is p2[j,nj,h, i,0+1,k, m].
                            double coef = qAt(q4, q5, i, j, k, u, 1);
                            if (coef == 0.0) continue;
                            for (int h = 0; h < K[j]; h++)
                                for (int m = 0; m < MR; m++) {
                                    if (BB[m][i] != 0.0) continue;
                                    row[p2Index(j, nj, h, i, 1, k, m, M, N, Kmax, MR)] -= coef;
                                }
                        }
                }
                aeq.add(row);
                beq.add(Double.valueOf(0.0));
            }

        // THM3 {i, ni in 0..F[i]-1}: balance across the ni -> ni+1 boundary.
        for (int i = 0; i < M; i++)
            for (int ni = 0; ni < F[i]; ni++) {
                double[] row = new double[numVars];
                for (int j = 0; j < M; j++) {
                    if (j == i) continue;
                    for (int nj = 1; nj <= F[j]; nj++)
                        for (int k = 0; k < K[j]; k++) {
                            double coef = 0.0;
                            for (int h = 0; h < K[j]; h++) coef += qAt(q4, q5, j, i, k, h, nj);
                            if (coef == 0.0) continue;
                            for (int u = 0; u < K[i]; u++)
                                for (int m = 0; m < MR; m++) {
                                    if (BB[m][j] != 0.0) continue;
                                    row[p2Index(j, nj, k, i, ni, u, m, M, N, Kmax, MR)] += coef;
                                }
                        }
                }
                for (int j = 0; j < M; j++) {
                    if (j == i) continue;
                    for (int nj = 0; nj <= F[j]; nj++)
                        for (int k = 0; k < K[i]; k++) {
                            double coef = 0.0;
                            for (int h = 0; h < K[i]; h++) coef += qAt(q4, q5, i, j, k, h, ni + 1);
                            if (coef == 0.0) continue;
                            for (int u = 0; u < K[j]; u++)
                                for (int m = 0; m < MR; m++) {
                                    if (BB[m][i] != 0.0) continue;
                                    row[p2Index(j, nj, u, i, ni + 1, k, m, M, N, Kmax, MR)] -= coef;
                                }
                        }
                }
                aeq.add(row);
                beq.add(Double.valueOf(0.0));
            }

        // THM4
        for (int j = 0; j < M; j++)
            for (int k = 0; k < K[j]; k++)
                for (int i = 0; i < M; i++)
                    for (int m = 0; m < MR; m++) {
                        double[] row = new double[numVars];
                        for (int t = 0; t < M; t++)
                            for (int h = 0; h < K[t]; h++)
                                for (int njIdx = 0; njIdx <= N; njIdx++)
                                    for (int nt = 0; nt <= N; nt++) {
                                        row[p2Index(j, njIdx, k, t, nt, h, m, M, N, Kmax, MR)] -= (double) nt;
                                    }
                        for (int h = 0; h < K[i]; h++)
                            for (int njIdx = 0; njIdx <= N; njIdx++)
                                for (int ni = 1; ni <= N; ni++) {
                                    row[p2Index(j, njIdx, k, i, ni, h, m, M, N, Kmax, MR)] += (double) N;
                                }
                        aub.add(row);
                        bub.add(Double.valueOf(0.0));
                    }

        return new ConstraintSet(aeq, beq, aub, bub);
    }
}
