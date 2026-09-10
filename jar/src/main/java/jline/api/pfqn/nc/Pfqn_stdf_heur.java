/**
 * @file Heuristic Sojourn Time Distribution Function for multi-server FCFS stations
 *
 * Implements a heuristic variant of McKenna's 1987 method for computing sojourn time
 * distributions at multi-server FCFS stations.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import jline.GlobalConstants;
import jline.api.mam.Map_cdf;
import jline.api.mam.Map_exponential;
import jline.api.mam.Map_sumind;
import jline.api.pfqn.ld.Pfqn_comomrm_ld;
import jline.api.pfqn.ld.Pfqn_mu_ms;
import jline.api.pfqn.ld.Pfqn_mushift;
import jline.api.pfqn.ld.Pfqn_mvald;
import jline.io.InputOutput;
import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Pfqn_stdf_heur {
    private Pfqn_stdf_heur() {}

    /**
     * Heuristic sojourn time distribution analysis at multiserver FCFS nodes.
     */
    public static Matrix[][] pfqn_stdf_heur(Matrix L, Matrix N, Matrix Z, Matrix S,
                                            Matrix fcfsNodes, Matrix rates, Matrix tset) {
        boolean stabilityWarnIssued = false;
        int M = L.getNumRows();
        int R = L.getNumCols();
        int T = tset.length();

        // Initialize service rate matrix
        Matrix mu = new Matrix(M, (int) N.elementSum());
        for (int k = 0; k < M; k++) {
            for (int n = 0; n < (int) N.elementSum(); n++) {
                mu.set(k, n, Math.min(S.get(k), (double) (n + 1)));
            }
        }

        // Ensure t > 0 for stability of hkc function
        Matrix tsetMod = tset.copy();
        for (int i = 0; i < tsetMod.length(); i++) {
            if (tsetMod.get(i) == 0.0) {
                tsetMod.set(i, GlobalConstants.FineTol);
            }
        }

        Matrix[][] RD = new Matrix[M][R];
        Matrix[] hkc = new Matrix[R];

        for (int k = 0; k < fcfsNodes.length(); k++) {
            int kIdx = (int) fcfsNodes.get(k);
            for (int r = 0; r < R; r++) {
                if (L.get(kIdx, r) > GlobalConstants.FineTol) {
                    Matrix Nr = Matrix.oner(N, new ArrayList<Integer>(Collections.singletonList(r)));

                    double lGr;
                    Matrix Q1;
                    boolean isNumStable;

                    if (M == 1) {
                        SolverOptions options = new SolverOptions();
                        Ret.pfqnComomrmLd result1 = Pfqn_comomrm_ld.pfqn_comomrm_ld(L, Nr, Z, mu, options);
                        int m = 2;
                        Matrix muMs = Pfqn_mu_ms.pfqn_mu_ms((int) N.elementSum(), m, (int) S.get(kIdx));
                        Ret.pfqnComomrmLd result2 = Pfqn_comomrm_ld.pfqn_comomrm_ld(L, Nr, Z, muMs, options);

                        Matrix Q1Local = new Matrix(M, R);
                        for (int s = 0; s < R; s++) {
                            Q1Local.set(kIdx, s, L.get(kIdx, s) * Math.exp(result2.lG - result1.lG));
                        }
                        lGr = result1.lG;
                        Q1 = Q1Local;
                        isNumStable = true;
                    } else {
                        Ret.pfqnMVALD result = Pfqn_mvald.pfqn_mvald(L, Nr, Z, mu);
                        lGr = result.lG.get(result.lG.size() - 1);
                        Q1 = result.Q;
                        isNumStable = result.isNumStable;
                    }

                    // Build heuristic MAPs for each population level
                    MatrixCell[] hMAPr = new MatrixCell[(int) N.elementSum() + 1];
                    for (int n = 0; n <= (int) N.elementSum(); n++) {
                        if (n < S.get(kIdx)) {
                            // Use average rates for n < S(k)
                            hMAPr[n] = Map_exponential.map_exponential(1.0 / rates.get(kIdx, r));
                        } else {
                            // Use rates per class for n >= S(k)
                            List<MatrixCell> mapList = new ArrayList<MatrixCell>();
                            mapList.add(Map_exponential.map_exponential(1.0 / rates.get(kIdx, r)));

                            for (int s = 0; s < R; s++) {
                                double qSum = Q1.sumRows(kIdx);
                                if (qSum > GlobalConstants.FineTol) {
                                    double rate = Q1.get(kIdx, s) * (n - S.get(kIdx) + 1) / qSum / rates.get(kIdx, s);
                                    mapList.add(Map_exponential.map_exponential(rate));
                                }
                            }
                            hMAPr[n] = Map_sumind.map_sumind(mapList.toArray(new MatrixCell[0]));
                        }
                    }

                    // Compute CDFs
                    hkc[r] = new Matrix(T, (int) N.elementSum() + 1);
                    for (int n = 0; n <= (int) N.elementSum(); n++) {
                        Matrix cdfValues = Map_cdf.map_cdf(hMAPr[n], tsetMod);
                        for (int t = 0; t < T; t++) {
                            hkc[r].set(t, n, cdfValues.get(t));
                        }
                    }

                    if (!isNumStable && !stabilityWarnIssued) {
                        stabilityWarnIssued = true;
                        InputOutput.line_warning("pfqn_stdf_heur",
                                "The computation of the sojourn time distribution is numerically unstable.");
                    }

                    RD[kIdx][r] = new Matrix(tset.length(), 2);
                    for (int t = 0; t < tset.length(); t++) {
                        RD[kIdx][r].set(t, 1, tset.get(t));
                    }

                    // Compute response time distribution using recursive form
                    double[] Hkrt = new double[T];
                    Matrix LReduced = new Matrix(L.getNumRows() - 1, L.getNumCols());
                    Matrix muReduced = new Matrix(mu.getNumRows() - 1, mu.getNumCols());
                    int rowIdx = 0;
                    for (int i = 0; i < M; i++) {
                        if (i != kIdx) {
                            for (int j = 0; j < R; j++) {
                                LReduced.set(rowIdx, j, L.get(i, j));
                            }
                            for (int j = 0; j < mu.getNumCols(); j++) {
                                muReduced.set(rowIdx, j, mu.get(i, j));
                            }
                            rowIdx++;
                        }
                    }

                    double lGk;
                    if (M == 1) {
                        // The network WITHOUT station k is EMPTY here, and its
                        // normalizing constant is the think-only one, NOT 0:
                        // G = prod_r Z_r^n_r / n_r!. Short-circuiting to 0.0
                        // asserts G = 1 and makes the returned CDF identically
                        // 1.0. Same defect and same fix as Pfqn_stdf.
                        // Z HAS ONE ROW PER DELAY NODE (snGetProductFormParams
                        // builds it (max(1,Mz) x R)), so class r's think time is
                        // the COLUMN SUM. Z.get(0,r) used the first delay only.
                        lGk = 0.0;
                        for (int rr = 0; rr < R; rr++) {
                            double nr = Nr.get(0, rr);
                            if (nr > 0) {
                                double Ztot = 0.0;
                                for (int zi = 0; zi < Z.getNumRows(); zi++) Ztot += Z.get(zi, rr);
                                lGk += nr * Math.log(Ztot) - Maths.factln((int) nr);
                            }
                        }
                    } else {
                        Ret.pfqnMVALD result = Pfqn_mvald.pfqn_mvald(LReduced, Nr, Z, muReduced);
                        lGk = result.lG.get(result.lG.size() - 1);
                    }

                    for (int t = 0; t < T; t++) {
                        Matrix gammat = mu.copy();
                        for (int m = 0; m < (int) Nr.elementSum(); m++) {
                            if (m + 1 < hkc[r].getNumCols()) {
                                double denominator = hkc[r].get(t, m + 1);
                                // > 0, NOT > FineTol: pfqn_stdf_heur.m:109-111
                                // divides unconditionally, so the wider guard is
                                // a port invention that discards legitimately
                                // small hkc values at small t.
                                if (denominator > 0.0) {
                                    gammat.set(kIdx, m, mu.get(kIdx, m) * hkc[r].get(t, m) / denominator);
                                }
                            }
                        }

                        Matrix gammak = Pfqn_mushift.pfqn_mushift(gammat, kIdx);
                        Hkrt[t] = hkc[r].get(t, 0) * Math.exp(lGk);

                        for (int s = 0; s < R; s++) {
                            if (Nr.get(s) > 0) {
                                Ret.pfqnRd lYksRd = Pfqn_rd.pfqn_rd(L,
                                        Matrix.oner(Nr, new ArrayList<Integer>(Collections.singletonList(s))),
                                        Z, gammak, null);
                                double lYks_t = lYksRd.lG;
                                if (gammat.get(kIdx, 0) > GlobalConstants.FineTol) {
                                    Hkrt[t] += (L.get(kIdx, s) * hkc[r].get(t, 0) / gammat.get(kIdx, 0))
                                            * Math.exp(lYks_t);
                                }
                            }
                        }
                    }

                    // Handle NaN values
                    for (int t = 0; t < T; t++) {
                        if (Double.isNaN(Hkrt[t])) {
                            Hkrt[t] = GlobalConstants.FineTol;
                        }
                    }

                    double[] lHkrt = new double[T];
                    for (int t = 0; t < T; t++) lHkrt[t] = Math.log(Hkrt[t]);

                    for (int t = 0; t < T; t++) {
                        RD[kIdx][r].set(t, 0, Math.exp(lHkrt[t] - lGr));
                    }
                }
            }
        }

        return RD;
    }

    /**
     * Helper function Fm(m,x) - not used in current implementation but kept for reference.
     */
    @SuppressWarnings("unused")
    private static double Fm(int m, double x) {
        if (m == 1) {
            return 1.0 - Math.exp(-x);
        } else {
            double A = 0.0;
            for (int j = 0; j < m; j++) {
                A += Math.pow(x, j) / Maths.fact(j);
            }
            return 1.0 - Math.exp(-x) * A;
        }
    }
}
