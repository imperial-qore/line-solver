package jline.lib.butools;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;

import jline.api.mc.Ctmc_solve;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class MMAPPH1FCFS {
    private MMAPPH1FCFS() {}

    /**
     * Convenience alias for callers that expect a {@code solve(...)} entry point
     * on the MMAPPH1FCFS class. Delegates to {@link #MMAPPH1FCFS}.
     */
    public static Map<String, Map<Integer, Matrix>> solve(MatrixCell D,
                                                          Map<Integer, Matrix> sigma,
                                                          Map<Integer, Matrix> S,
                                                          Integer numOfQLMoms,
                                                          Integer numOfQLProbs,
                                                          Integer numOfSTMoms,
                                                          Matrix stDistr,
                                                          boolean stDistrME,
                                                          boolean stDistrPH,
                                                          Double prec,
                                                          Matrix classes_) {
        return MMAPPH1FCFS(D, sigma, S, numOfQLMoms, numOfQLProbs, numOfSTMoms,
                           stDistr, stDistrME, stDistrPH, prec, classes_);
    }

    public static Map<String, Map<Integer, Matrix>> MMAPPH1FCFS(MatrixCell D,
                                                                  Map<Integer, Matrix> sigma,
                                                                  Map<Integer, Matrix> S,
                                                                  Integer numOfQLMoms,
                                                                  Integer numOfQLProbs,
                                                                  Integer numOfSTMoms,
                                                                  Matrix stDistr,
                                                                  boolean stDistrME,
                                                                  boolean stDistrPH,
                                                                  Double prec,
                                                                  Matrix classes_) {
        int K = D.size() - 1;
        double precision = 1e-14;
        Matrix classes = new Matrix(1, K, K);
        for (int i = 0; i < K; i++) {
            classes.set(i, (double) i);
        }
        if (prec != null) {
            precision = prec.doubleValue();
        }

        if (classes_ != null) {
            classes = classes_;
        }

        Matrix D0 = D.get(0);
        int N = D0.getNumRows();
        Matrix Ia = Matrix.eye(N);
        Matrix Da = new Matrix(N, N, N * N);
        for (int q = 0; q < K; q++) {
            Da = Da.add(1.0, D.get(q + 1));
        }
        Map<Integer, Matrix> beta = new HashMap<Integer, Matrix>();
        Matrix theta = Ctmc_solve.ctmc_solve(D0.add(1.0, Da));
        Matrix lambda = new Matrix(K, 1, K);
        Matrix mu = new Matrix(K, 1, K);
        Matrix Nsk = new Matrix(1, K, K);
        double ro = 0.0;
        for (int k = 0; k < K; k++) {
            lambda.set(k, theta.mult(D.get(k + 1)).elementSum());
            beta.put(Integer.valueOf(k),
                    Ctmc_solve.ctmc_solve(S.get(Integer.valueOf(k))
                            .add(-1.0, S.get(Integer.valueOf(k)).sumRows().mult(sigma.get(Integer.valueOf(k))))));
            Matrix neg_sk = S.get(Integer.valueOf(k)).copy();
            neg_sk.scaleEq(-1.0);
            mu.set(k, beta.get(Integer.valueOf(k)).mult(neg_sk).elementSum());
            Nsk.set(k, (double) S.get(Integer.valueOf(k)).getNumRows());
            ro = ro + lambda.get(k) / mu.get(k);
        }

        Matrix alpha = theta.mult(Da);
        alpha.scaleEq(1 / lambda.elementSum());
        // MATLAB BUTools computes D0i = inv(-D0) with plain inv(), which merely
        // warns (not errors) when D0 is ill-conditioned to working precision and
        // still returns a usable LU inverse. Matrix.inv() instead hard-throws at
        // det < 1e-14, breaking mixed multiserver-FCFS dec.source models whose
        // aggregate D0 is near-singular. robustLeftDivide uses the same LU solve
        // (matching MATLAB) and only falls back to a pseudo-inverse if genuinely
        // singular, so it reproduces MATLAB's behaviour without spurious errors.
        Matrix D0i = Matrix.robustLeftDivide(D0, Matrix.eye(D0.getNumRows())).copy();
        D0i.scaleEq(-1.0);

        Matrix Sa = S.get(Integer.valueOf(0));
        Map<Integer, Matrix> sa = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> ba = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> sv = new HashMap<Integer, Matrix>();
        sa.put(Integer.valueOf(0), sigma.get(Integer.valueOf(0)));
        ba.put(Integer.valueOf(0), beta.get(Integer.valueOf(0)));
        Matrix sv0 = S.get(Integer.valueOf(0)).sumRows();
        sv0.scaleEq(-1.0);
        sv.put(Integer.valueOf(0), sv0);

        Map<Integer, Matrix> Pk = new HashMap<Integer, Matrix>();
        Pk.put(Integer.valueOf(0), D0i.mult(D.get(1)));

        for (int q = 1; q < K; q++) {
            int sigLen = sigma.get(Integer.valueOf(0)).length();
            int betaLen = beta.get(Integer.valueOf(0)).length();
            sa.put(Integer.valueOf(q), new Matrix(1, sigLen, sigLen));
            ba.put(Integer.valueOf(q), new Matrix(1, betaLen, betaLen));
            sv.put(Integer.valueOf(q), new Matrix(sigLen, 1, sigLen));
            Pk.put(Integer.valueOf(q), D0i.mult(D.get(q + 1)));
        }

        for (int k = 1; k < K; k++) {
            Sa = Sa.createBlockDiagonal(S.get(Integer.valueOf(k)));
            for (int q = 0; q < K; q++) {
                if (q == k) {
                    sa.put(Integer.valueOf(q),
                            Matrix.concatColumns(sa.get(Integer.valueOf(q)), sigma.get(Integer.valueOf(k)), null));
                    ba.put(Integer.valueOf(q),
                            Matrix.concatColumns(ba.get(Integer.valueOf(q)), beta.get(Integer.valueOf(k)), null));
                    Matrix sk_neg_sum = S.get(Integer.valueOf(k)).sumRows();
                    sk_neg_sum.scaleEq(-1.0);
                    sv.put(Integer.valueOf(q),
                            Matrix.concatRows(sv.get(Integer.valueOf(q)), sk_neg_sum, null));
                } else {
                    Matrix sigK = sigma.get(Integer.valueOf(k));
                    Matrix betK = beta.get(Integer.valueOf(k));
                    sa.put(Integer.valueOf(q),
                            Matrix.concatColumns(sa.get(Integer.valueOf(q)),
                                    new Matrix(sigK.getNumRows(), sigK.getNumCols(), 0), null));
                    ba.put(Integer.valueOf(q),
                            Matrix.concatColumns(ba.get(Integer.valueOf(q)),
                                    new Matrix(betK.getNumRows(), betK.getNumCols(), 0), null));
                    sv.put(Integer.valueOf(q),
                            Matrix.concatRows(sv.get(Integer.valueOf(q)),
                                    new Matrix(sigK.length(), 1, 0), null));
                }
            }
        }
        D0i.mult(Da);
        Matrix iVec = D.get(1).kron(sa.get(Integer.valueOf(0)));

        for (int k = 1; k < K; k++) {
            iVec = iVec.add(1.0, D.get(k + 1).kron(sa.get(Integer.valueOf(k))));
        }

        int Ns = Sa.getNumRows();
        Matrix Is = Matrix.eye(Ns);
        Matrix neg_Sa_row_sum = Sa.sumRows();
        neg_Sa_row_sum.scaleEq(-1.0);
        Matrix Y0 = FluidFundamentalMatrices.FluidFundamentalMatrices(Ia.kron(Sa),
                Ia.kron(neg_Sa_row_sum),
                iVec,
                D0,
                precision,
                null,
                null).get("P");
        Matrix T = Ia.kron(Sa).add(1.0, Y0.mult(iVec));
        Matrix pi0 = new Matrix(1, T.getNumRows(), T.getNumRows());
        for (int k = 0; k < K; k++) {
            Matrix ba_mu = ba.get(Integer.valueOf(k)).copy();
            ba_mu.scaleEq(1 / mu.get(k));
            pi0 = pi0.add(1.0, theta.mult(D.get(k + 1)).kron(ba_mu));
        }
        pi0 = pi0.mult(T);
        pi0.scaleEq(-1.0);

        Matrix iT = T.inv();
        iT.scaleEq(-1.0);
        Matrix oa = Matrix.ones(N, 1);

        Map<String, Map<Integer, Matrix>> result = new HashMap<String, Map<Integer, Matrix>>();
        if (numOfSTMoms != null) {
            result.put("stNoms", new HashMap<Integer, Matrix>());
        }
        if (stDistr != null) {
            result.put("stDistr", new HashMap<Integer, Matrix>());
        }
        if (stDistrME) {
            result.put("stDistrME_alpha", new HashMap<Integer, Matrix>());
            result.put("stDistrME_A", new HashMap<Integer, Matrix>());
        }
        if (stDistrPH) {
            result.put("stDistrPH_alpha", new HashMap<Integer, Matrix>());
            result.put("stDistrPH_A", new HashMap<Integer, Matrix>());
        }
        if (numOfQLProbs != null) {
            result.put("ncDistr", new HashMap<Integer, Matrix>());
        }
        if (numOfQLMoms != null) {
            result.put("ncMoms", new HashMap<Integer, Matrix>());
        }

        for (int i = 0; i < classes.length(); i++) {
            int k = (int) classes.get(i);
            Matrix clo = iT.mult(oa.kron(sv.get(Integer.valueOf(k))));
            if (numOfSTMoms != null) {
                Matrix rtMoms = new Matrix(1, numOfSTMoms.intValue());
                for (int m = 0; m < numOfSTMoms.intValue(); m++) {
                    // (m+1)-th sojourn moment: factorial(m+1)*pi0*iT^(m+1)*clo / (pi0*clo)
                    // (BUTools MMAPPH1FCFS "stMoms"; note clo already contains one iT factor)
                    double rtMoms_m = pi0.mult(Matrix.pow(iT, m + 1)).mult(clo).get(0) / pi0.mult(clo).get(0);
                    rtMoms.set(m, CombinatoricsUtils.factorial(m + 1) * rtMoms_m);
                }
                result.get("stNoms").put(Integer.valueOf(k), rtMoms);
            }

            if (stDistr != null) {
                Matrix cdf = new Matrix(0, 0, 0);
                for (int p = 0; p < stDistr.length(); p++) {
                    Matrix Tt = T.copy();
                    Tt.scaleEq(stDistr.get(p));
                    double pr = 1 - pi0.mult(Maths.matrixExp(Tt)).mult(clo).get(0) / pi0.mult(clo).get(0);
                    Matrix cdf_ = new Matrix(1, 1, 1);
                    cdf_.set(0, 0, pr);
                    if (p == 0) {
                        cdf = cdf_;
                    } else {
                        cdf = Matrix.concatColumns(cdf, cdf_, null);
                    }
                }
                result.get("stDistr").put(Integer.valueOf(k), cdf);
            }

            if (stDistrME) {
                Matrix Bm = SimilarityMatrixForVectors.SimilarityMatrixForVectors(
                        clo.mult(pi0.mult(clo).inv()), Matrix.ones(N * Ns, 1));
                Matrix Bmi = Bm.inv();
                Matrix A = Bm.mult(T).mult(Bmi);
                pi0.mult(Bmi);
                result.get("stDistrME_alpha").put(Integer.valueOf(k), alpha);
                result.get("stDistrME_A").put(Integer.valueOf(k), A);
            }

            if (stDistrPH) {
                Matrix vv = pi0.mult(iT);

                List<Double> nz = new ArrayList<Double>();
                List<Integer> nz_index = new ArrayList<Integer>();
                for (int n = 0; n < vv.length(); n++) {
                    if (vv.get(n) > precision) {
                        nz.add(Double.valueOf(vv.get(n)));
                        nz_index.add(Integer.valueOf(n));
                    }
                }
                double[] nzArr = new double[nz.size()];
                for (int n = 0; n < nz.size(); n++) nzArr[n] = nz.get(n).doubleValue();
                Matrix delta = Matrix.diag(nzArr);
                Matrix neg_T = T.copy();
                neg_T.scaleEq(-1.0);
                Matrix cl = neg_T.mult(clo);
                cl.scaleEq(pi0.mult(clo).value());
                Matrix alpha_stDistrPH = new Matrix(0, 0, 0);
                for (int n = 0; n < nz_index.size(); n++) {
                    if (alpha_stDistrPH.length() == 0) {
                        alpha_stDistrPH = Matrix.extractRows(cl, nz_index.get(n).intValue(),
                                nz_index.get(n).intValue() + 1, null);
                    } else {
                        Matrix.concatRows(alpha_stDistrPH,
                                Matrix.extractRows(cl, nz_index.get(n).intValue(),
                                        nz_index.get(n).intValue() + 1, null),
                                null);
                    }
                }
                Matrix A = new Matrix(nz_index.size(), nz_index.size(),
                        (int) Math.pow((double) nz_index.size(), 2.0));
                for (int n = 0; n < nz_index.size(); n++) {
                    for (int m = 0; m < nz_index.size(); m++) {
                        A.set(nz_index.get(n).intValue(), nz_index.get(m).intValue(),
                                T.get(nz_index.get(n).intValue(), nz_index.get(m).intValue()));
                    }
                }
                A = delta.inv().mult(A.transpose()).mult(delta);
                result.get("stDistrPH_alpha").put(Integer.valueOf(k), alpha);
                result.get("stDistrPH_A").put(Integer.valueOf(k), A);
            }

            if (numOfQLProbs != null) {
                Matrix value = new Matrix(1, numOfQLProbs.intValue(), numOfQLProbs.intValue());
                Matrix jm = new Matrix(Ns, 1, Ns);
                {
                    int n = (int) Nsk.sumRows(0, k).elementSum();
                    while (n < Nsk.sumRows(0, k + 1).elementSum()) {
                        jm.set(n, 0, 1.0);
                        n++;
                    }
                }
                Matrix jmc = Matrix.ones(Ns, 1);
                jmc = jmc.add(-1.0, jm);
                Matrix LmCurr = Matrix.lyap(T, D0.add(1.0, Da).add(-1.0, D.get(k + 1)).kron(Is),
                        Matrix.eye(N * Ns), null);
                value.set(0, 1 - ro + pi0.mult(LmCurr).mult(oa.kron(jmc)).get(0));
                for (int n = 0; n < numOfQLProbs.intValue() - 1; n++) {
                    Matrix LmPrev = LmCurr.copy();
                    LmCurr = Matrix.lyap(T, D0.add(1.0, Da).add(-1.0, D.get(k + 1)).kron(Is),
                            LmPrev.mult(D.get(k + 1).kron(Is)), null);
                    value.set(n + 1, pi0.mult(LmCurr).mult(oa.kron(jmc)).get(0)
                            + pi0.mult(LmPrev).mult(oa.kron(jm)).get(0));
                }
                result.get("ncDistr").put(Integer.valueOf(k), value);
            }

            if (numOfQLMoms != null) {
                Matrix jm = new Matrix(Ns, 1, Ns);
                {
                    int n = (int) Nsk.sumRows(0, k).elementSum();
                    while (n < Nsk.sumRows(0, k + 1).elementSum()) {
                        jm.set(n, 0, 1.0);
                        n++;
                    }
                }
                Map<Integer, Matrix> ELn = new HashMap<Integer, Matrix>();
                ELn.put(Integer.valueOf(0), Matrix.lyap(T, D0.add(1.0, Da).kron(Is), Matrix.eye(N * Ns), null));
                Matrix qlMoms = new Matrix(1, numOfQLMoms.intValue(), numOfQLMoms.intValue());
                for (int n = 0; n < numOfQLMoms.intValue(); n++) {
                    int bino = 1;
                    Matrix Btag = new Matrix(N * Ns, N * Ns, (int) FastMath.pow((double) (N * Ns), 2));
                    for (int m = -1; m < n; m++) {
                        Btag = Btag.add((double) bino, ELn.get(Integer.valueOf(m + 1)));
                        bino = bino * (n - m) / (m + 2);
                    }
                    ELn.put(Integer.valueOf(n + 1),
                            Matrix.lyap(T, D0.add(1.0, Da).kron(Is), Btag.mult(D.get(k + 1).kron(Is)), null));
                    qlMoms.set(n, pi0.mult(ELn.get(Integer.valueOf(n + 1))).elementSum()
                            + pi0.mult(Btag).mult(oa.kron(jm)).get(0));
                }
                result.get("ncMoms").put(Integer.valueOf(k), qlMoms);
            }
        }

        return result;
    }
}
