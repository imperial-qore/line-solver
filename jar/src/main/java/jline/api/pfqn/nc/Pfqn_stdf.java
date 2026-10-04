/**
 * @file Sojourn Time Distribution Function (STDF) for multi-server FCFS stations
 *
 * Implements McKenna's 1987 method for computing sojourn time distributions at
 * multi-server FCFS stations in closed queueing networks.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import java.util.ArrayList;
import java.util.Collections;

import jline.GlobalConstants;
import jline.api.mam.Map_cdf;
import jline.api.mam.Map_erlang;
import jline.api.mam.Map_exponential;
import jline.api.mam.Map_sumind;
import jline.api.pfqn.ld.Pfqn_comomrm_ld;
import jline.api.pfqn.ld.Pfqn_mushift;
import jline.api.pfqn.ld.Pfqn_mvald;
import jline.io.InputOutput;
import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Pfqn_stdf {
    private Pfqn_stdf() {}

    /**
     * Sojourn time distribution function at multiserver FCFS nodes (McKenna 1987 JACM).
     */
    public static Matrix[][] pfqn_stdf(Matrix L, Matrix N, Matrix Z, Matrix S,
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

        for (int k = 0; k < fcfsNodes.length(); k++) {
            int kIdx = (int) fcfsNodes.get(k);
            // Check that FCFS station has uniform service rates
            Matrix rateRow = new Matrix(1, R);
            Matrix.extract(rates, kIdx, kIdx + 1, 0, R, rateRow, 0, 0);
            double rateRange = rateRow.elementMax() - rateRow.elementMin();
            if (rateRange > GlobalConstants.FineTol) {
                throw new IllegalArgumentException("The FCFS stations has distinct service rates, the model is invalid.");
            }

            // Build MAPs for each population level
            MatrixCell[] hMAP = new MatrixCell[(int) N.elementSum() + 1];
            for (int n = 0; n <= (int) N.elementSum(); n++) {
                if (n < S.get(kIdx)) {
                    hMAP[n] = Map_exponential.map_exponential(1.0 / rates.get(kIdx, 0));
                } else {
                    MatrixCell[] mapList = new MatrixCell[]{
                            Map_exponential.map_exponential(1.0 / rates.get(kIdx, 0)),
                            Map_erlang.map_erlang((n - S.get(kIdx) + 1) / (S.get(kIdx) * rates.get(kIdx, 0)),
                                    (int) (n - S.get(kIdx) + 1))
                    };
                    hMAP[n] = Map_sumind.map_sumind(mapList);
                }
            }

            // Compute CDFs
            Matrix hkc = new Matrix(T, (int) N.elementSum() + 1);
            for (int n = 0; n <= (int) N.elementSum(); n++) {
                Matrix cdfValues = Map_cdf.map_cdf(hMAP[n], tsetMod);
                for (int t = 0; t < T; t++) {
                    hkc.set(t, n, cdfValues.get(t));
                }
            }

            for (int r = 0; r < R; r++) {
                if (L.get(kIdx, r) > GlobalConstants.FineTol) {
                    Matrix Nr = Matrix.oner(N, new ArrayList<Integer>(Collections.singletonList(r)));

                    double lGr;
                    boolean isNumStable;

                    if (L.getNumRows() == 1) {
                        SolverOptions options = new SolverOptions();
                        Ret.pfqnComomrmLd result = Pfqn_comomrm_ld.pfqn_comomrm_ld(L, Nr, Z, mu, options);
                        lGr = result.lG;
                        isNumStable = true;
                    } else {
                        Ret.pfqnMVALD result = Pfqn_mvald.pfqn_mvald(L, Nr, Z, mu);
                        lGr = result.lG.get(result.lG.size() - 1);
                        isNumStable = result.isNumStable;
                    }

                    if (!isNumStable && !stabilityWarnIssued) {
                        stabilityWarnIssued = true;
                        InputOutput.line_warning("pfqn_stdf",
                                "The computation of the sojourn time distribution is numerically unstable.");
                    }

                    RD[kIdx][r] = new Matrix(tset.length(), 2);
                    for (int t = 0; t < tset.length(); t++) {
                        RD[kIdx][r].set(t, 1, tset.get(t));
                    }

                    // Compute response time distribution using recursive form
                    double[] Hkrt = new double[T];

                    // Create reduced matrices excluding station k
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
                    if (L.getNumRows() == 1) {
                        // The network WITHOUT station k is EMPTY here, and its
                        // normalizing constant is the think-only one, NOT 0:
                        // G = prod_r Z_r^n_r / n_r!. Short-circuiting to 0.0
                        // asserts G = 1 and makes the returned CDF identically
                        // 1.0. The reference calls pfqn_comomrm_ld on the empty
                        // slice (pfqn_stdf.m:73-77), which returns exactly the
                        // closed form used below (verified to 0.0 over
                        // n = 1..4, Z = 0.5..2).
                        // Z HAS ONE ROW PER DELAY NODE, not one row overall:
                        // snGetProductFormParams builds it as (max(1,Mz) x R).
                        // The think time of class r is therefore the COLUMN SUM,
                        // and reading Z.get(0,r) silently used the first delay
                        // only -- correct on the single-delay models this
                        // shortcut was validated against, wrong on any other.
                        // The caller passes Z through untouched, so the summing
                        // has to happen here; pfqn_comomrm_ld and pfqn_mvald,
                        // which the other branch delegates to, accept the column
                        // and are why only this inline shortcut was affected.
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
                            if (m + 1 < hkc.getNumCols()) {
                                double denominator = hkc.get(t, m + 1);
                                // > 0, NOT > FineTol: the reference divides
                                // unconditionally, and hkc is a CDF value that is
                                // legitimately far below FineTol at small t, so
                                // the wider guard discarded real terms and left
                                // gamma at mu. Zero is still excluded because hkc
                                // underflows to exactly 0 at extreme t, where the
                                // ratio is undefined rather than merely small.
                                if (denominator > 0.0) {
                                    gammat.set(kIdx, m, mu.get(kIdx, m) * hkc.get(t, m) / denominator);
                                }
                            }
                        }

                        Matrix gammak = Pfqn_mushift.pfqn_mushift(gammat, kIdx);
                        // Extract reduced gammak matrix (excluding one column)
                        // Clamped at 0: with a single job the reduced population
                        // is empty, sum(Nr)-1 is -1, and the Matrix constructor
                        // below throws on a negative column count. MATLAB's
                        // `gammak(:,1:(sum(Nr)-1))` yields an EMPTY slice there
                        // (pfqn_stdf.m:86), which is what 0 columns reproduces.
                        // The existing `if (gammakCols > 0)` guard sits AFTER
                        // the constructor and so never prevented the throw.
                        int gammakCols = Math.max(0,
                                Math.min(gammak.getNumCols(), (int) Nr.elementSum() - 1));
                        Matrix gammakReduced = new Matrix(gammak.getNumRows(), gammakCols);
                        if (gammakCols > 0) {
                            Matrix.extract(gammak, 0, gammak.getNumRows(), 0, gammakCols, gammakReduced, 0, 0);
                        }

                        Hkrt[t] = hkc.get(t, 0) * Math.exp(lGk);

                        for (int s = 0; s < R; s++) {
                            if (Nr.get(s) > 0) {
                                double lYks_t;
                                if (L.getNumRows() == 1) {
                                    SolverOptions options = new SolverOptions();
                                    Ret.pfqnComomrmLd result =
                                            Pfqn_comomrm_ld.pfqn_comomrm_ld(L,
                                                    Matrix.oner(Nr, new ArrayList<Integer>(Collections.singletonList(s))),
                                                    Z, gammakReduced, options);
                                    lYks_t = result.lG;
                                } else {
                                    Ret.pfqnMVALD result = Pfqn_mvald.pfqn_mvald(L,
                                            Matrix.oner(Nr, new ArrayList<Integer>(Collections.singletonList(s))),
                                            Z, gammakReduced);
                                    lYks_t = result.lG.get(result.lG.size() - 1);
                                }

                                if (gammat.get(kIdx, 0) > GlobalConstants.FineTol) {
                                    Hkrt[t] += (L.get(kIdx, s) * hkc.get(t, 0) / gammat.get(kIdx, 0)) * Math.exp(lYks_t);
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
                        RD[kIdx][r].set(t, 0, Math.min(1.0, Math.exp(lHkrt[t] - lGr)));
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
