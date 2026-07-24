/**
 * Analyzes multi-class MMAP[K]/PH[K]/1 queue with preemptive-resume (PRPR) priority.
 */
package jline.lib.butools;

import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

import jline.api.mc.Ctmc_solve;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.CombinatoricsUtils;

public final class MMAPPH1PRPR {
    private MMAPPH1PRPR() {}

    public static Map<String, Map<Integer, Matrix>> MMAPPH1PRPR(
            MatrixCell D,
            MatrixCell sigma,
            MatrixCell S,
            Integer numOfQLMoms,
            Integer numOfQLProbs,
            Integer numOfSTMoms,
            Matrix stCdfPoints,
            Double prec,
            Integer erlMaxOrder_,
            Matrix classes_) {
        int K = D.size() - 1;

        int erlMaxOrder = (erlMaxOrder_ != null) ? erlMaxOrder_ : 200;
        double precision = (prec != null) ? prec : 1e-14;

        Matrix classes = new Matrix(1, K, K);
        for (int i = 0; i < K; i++) {
            classes.set(i, (double) i);
        }
        if (classes_ != null) {
            classes = classes_;
        }

        Matrix D0 = D.get(0);
        int N = D0.getNumRows();
        Matrix I = Matrix.eye(N);

        Matrix sD = new Matrix(N, N, N * N);
        for (int i = 0; i <= K; i++) {
            sD = sD.add(1.0, D.get(i));
        }

        Map<Integer, Matrix> s = new HashMap<Integer, Matrix>();
        Matrix M = new Matrix(1, K, K);
        for (int i = 0; i < K; i++) {
            Matrix neg_Si = S.get(i).copy();
            neg_Si.scaleEq(-1.0);
            s.put(i, neg_Si.sumCols());
            M.set(i, (double) sigma.get(i).length());
        }

        Map<String, Map<Integer, Matrix>> result = new HashMap<String, Map<Integer, Matrix>>();
        if (numOfSTMoms != null) result.put("stMoms", new HashMap<Integer, Matrix>());
        if (stCdfPoints != null) result.put("stDistr", new HashMap<Integer, Matrix>());
        if (numOfQLMoms != null) result.put("ncMoms", new HashMap<Integer, Matrix>());
        if (numOfQLProbs != null) result.put("ncDistr", new HashMap<Integer, Matrix>());

        for (int g = 0; g < classes.length(); g++) {
            int k = (int) classes.get(g);

            int sM = 0;
            for (int i = k; i < K; i++) sM += (int) M.get(i);

            Matrix Qwmm = D0.copy();
            for (int i = 0; i < k; i++) {
                Qwmm = Qwmm.add(1.0, D.get(i + 1));
            }

            Matrix Qwpm = new Matrix(N * sM, N, N * sM * N);
            Matrix Qwmp = new Matrix(N, N * sM, N * N * sM);
            Matrix Qwpp = new Matrix(N * sM, N * sM, N * sM * N * sM);

            int kix = 0;
            for (int i = k; i < K; i++) {
                int bs = N * (int) M.get(i);
                Qwmp.insertSubMatrix(0, kix, N, kix + bs, D.get(i + 1).kron(sigma.get(i)));
                Qwpm.insertSubMatrix(kix, 0, kix + bs, N, I.kron(s.get(i)));
                Qwpp.insertSubMatrix(kix, kix, kix + bs, kix + bs, I.kron(S.get(i)));
                kix += bs;
            }

            Map<String, Matrix> fluidResult = FluidFundamentalMatrices.FluidFundamentalMatrices(
                    Qwpp, Qwpm, Qwmp, Qwmm, precision, null, null);
            Matrix Psiw = fluidResult.get("P");
            Matrix Kw = fluidResult.get("K");
            Matrix Uw = fluidResult.get("U");

            Matrix neg_Kw = Kw.copy();
            neg_Kw.scaleEq(-1.0);
            Matrix Ua = Matrix.ones(N, 1).add(2.0, Qwmp.mult(neg_Kw.pinv()).sumRows());

            Matrix A_ls = Matrix.concatColumns(Uw, Ua, null).transpose();
            Matrix b_ls = Matrix.concatColumns(new Matrix(1, N, 0), Matrix.singleton(1.0), null).transpose();
            Matrix pm = A_ls.leftMatrixDivide(b_ls).transpose();

            Matrix Bw = new Matrix(N * sM, N, N * sM * N);
            Bw.insertSubMatrix(0, 0, N * (int) M.get(k), N, I.kron(s.get(k)));

            double denom = pm.mult(Qwmp).mult(neg_Kw.pinv()).mult(Bw).elementSum();
            Matrix kappa = pm.mult(Qwmp);
            kappa.scaleEq(1.0 / denom);

            if (k < K - 1) {
                Matrix Qsmm = D0.copy();
                for (int i = 0; i <= k; i++) {
                    Qsmm = Qsmm.add(1.0, D.get(i + 1));
                }

                int Np = Kw.getNumRows();
                int sizeHigh = 0;
                for (int i = k + 1; i < K; i++) {
                    sizeHigh += N * (int) M.get(i);
                }

                Matrix Qspm = new Matrix(Np + sizeHigh, N, (Np + sizeHigh) * N);
                Matrix Qsmp = new Matrix(N, Np + sizeHigh, N * (Np + sizeHigh));
                Matrix Qspp = new Matrix(Np + sizeHigh, Np + sizeHigh, (Np + sizeHigh) * (Np + sizeHigh));

                Qspp.insertSubMatrix(0, 0, Np, Np, Kw);
                Qspm.insertSubMatrix(0, 0, Np, N, Bw);

                kix = Np;
                for (int i = k + 1; i < K; i++) {
                    int bs = N * (int) M.get(i);
                    Qsmp.insertSubMatrix(0, kix, N, kix + bs, D.get(i + 1).kron(sigma.get(i)));
                    Qspm.insertSubMatrix(kix, 0, kix + bs, N, I.kron(s.get(i)));
                    Qspp.insertSubMatrix(kix, kix, kix + bs, kix + bs, I.kron(S.get(i)));
                    kix += bs;
                }

                Matrix inis = Matrix.concatColumns(kappa, new Matrix(1, sizeHigh, 0), null);

                Matrix Psis = FluidFundamentalMatrices.FluidFundamentalMatrices(
                        Qspp, Qspm, Qsmp, Qsmm, precision, null, null).get("P");

                if (numOfSTMoms != null) {
                    Map<Integer, Matrix> Pn = new HashMap<Integer, Matrix>();
                    Pn.put(0, Psis);
                    Matrix rtMoms = new Matrix(1, numOfSTMoms, numOfSTMoms);

                    for (int n = 1; n <= numOfSTMoms; n++) {
                        Matrix A = Qspp.add(1.0, Psis.mult(Qsmp));
                        Matrix B = Qsmm.add(1.0, Qsmp.mult(Psis));
                        Matrix C = Pn.get(n - 1).copy();
                        C.scaleEq(-2.0 * n);

                        int bino = 1;
                        for (int i = 1; i < n; i++) {
                            bino = bino * (n - i + 1) / i;
                            C = C.add((double) bino, Pn.get(i).mult(Qsmp).mult(Pn.get(n - i)));
                        }

                        Matrix P = Matrix.lyap(A, B, C, null);
                        Pn.put(n, P);
                        rtMoms.set(n - 1, inis.mult(P).elementSum() * Math.pow(-1.0, n) / Math.pow(2.0, n));
                    }
                    result.get("stMoms").put(k, rtMoms);
                }

                if (stCdfPoints != null) {
                    Matrix res = new Matrix(1, stCdfPoints.length(), stCdfPoints.length());
                    for (int o = 0; o < stCdfPoints.length(); o++) {
                        double t = stCdfPoints.get(o);
                        int L = erlMaxOrder;
                        double lambda = L / t / 2.0;

                        Matrix Psie = FluidFundamentalMatrices.FluidFundamentalMatrices(
                                Qspp.add(-lambda, Matrix.eye(Qspp.getNumRows())),
                                Qspm,
                                Qsmp,
                                Qsmm.add(-lambda, Matrix.eye(Qsmm.getNumRows())),
                                precision, null, null).get("P");

                        Map<Integer, Matrix> Pn = new HashMap<Integer, Matrix>();
                        Pn.put(0, Psie);
                        double pr = inis.mult(Psie).elementSum();

                        for (int n = 1; n < L; n++) {
                            Matrix A = Qspp.add(1.0, Psie.mult(Qsmp)).add(-lambda, Matrix.eye(Qspp.getNumRows()));
                            Matrix B = Qsmm.add(1.0, Qsmp.mult(Psie)).add(-lambda, Matrix.eye(Qsmm.getNumRows()));
                            Matrix C = Pn.get(n - 1).copy();
                            C.scaleEq(2.0 * lambda);
                            for (int i = 1; i < n; i++) {
                                C = C.add(1.0, Pn.get(i).mult(Qsmp).mult(Pn.get(n - i)));
                            }
                            Matrix P = Matrix.lyap(A, B, C, null);
                            Pn.put(n, P);
                            pr += inis.mult(P).elementSum();
                        }
                        res.set(o, pr);
                    }
                    result.get("stDistr").put(k, res);
                }

                if (numOfQLMoms != null) {
                    Map<Integer, Matrix> QLDPn = new HashMap<Integer, Matrix>();
                    QLDPn.put(0, Psis);
                    Matrix dqlMoms = new Matrix(1, numOfQLMoms, numOfQLMoms);

                    for (int n = 1; n <= numOfQLMoms; n++) {
                        Matrix A = Qspp.add(1.0, Psis.mult(Qsmp));
                        Matrix B = Qsmm.add(1.0, Qsmp.mult(Psis));
                        Matrix C = QLDPn.get(n - 1).mult(D.get(k + 1));
                        C.scaleEq((double) n);

                        int bino = 1;
                        for (int i = 1; i < n; i++) {
                            bino = bino * (n - i + 1) / i;
                            C = C.add((double) bino, QLDPn.get(i).mult(Qsmp).mult(QLDPn.get(n - i)));
                        }
                        Matrix P = Matrix.lyap(A, B, C, null);
                        QLDPn.put(n, P);
                        dqlMoms.set(n - 1, inis.mult(P).elementSum());
                    }
                    Matrix dqlMomsConverted = MomsFromFactorialMoms.MomsFromFactorialMoms(dqlMoms);

                    Matrix pi = Ctmc_solve.ctmc_solve(sD);
                    double lambdak = pi.mult(D.get(k + 1)).elementSum();
                    Map<Integer, Matrix> QLPn = new HashMap<Integer, Matrix>();
                    QLPn.put(0, pi);
                    Matrix qlMoms = new Matrix(1, numOfQLMoms, numOfQLMoms);

                    Matrix iTerm = Matrix.ones(N, 1).mult(pi).add(-1.0, sD).pinv();

                    for (int n = 1; n <= numOfQLMoms; n++) {
                        double sumP = inis.mult(QLDPn.get(n)).elementSum()
                                + n * (inis.mult(QLDPn.get(n - 1)).add(-1.0 / lambdak,
                                        QLPn.get(n - 1).mult(D.get(k + 1))))
                                .mult(iTerm).mult(D.get(k + 1).sumRows()).get(0);
                        Matrix P = Matrix.scaleMult(pi, sumP).add(
                                (double) n,
                                QLPn.get(n - 1).mult(D.get(k + 1))
                                        .add(-lambdak, inis.mult(QLDPn.get(n - 1))).mult(iTerm));
                        QLPn.put(n, P);
                        qlMoms.set(n - 1, P.elementSum());
                    }
                    result.get("ncMoms").put(k, MomsFromFactorialMoms.MomsFromFactorialMoms(qlMoms));
                }

                if (numOfQLProbs != null) {
                    Matrix sDk = D0.copy();
                    for (int i = 0; i < k; i++) {
                        sDk = sDk.add(1.0, D.get(i + 1));
                    }

                    Matrix Psid = FluidFundamentalMatrices.FluidFundamentalMatrices(
                            Qspp, Qspm, Qsmp, sDk, precision, null, null).get("P");
                    Map<Integer, Matrix> Pn = new HashMap<Integer, Matrix>();
                    Pn.put(0, Psid);
                    Matrix dqlProbs = inis.mult(Psid);

                    for (int n = 1; n < numOfQLProbs; n++) {
                        Matrix A = Qspp.add(1.0, Psid.mult(Qsmp));
                        Matrix B = sDk.add(1.0, Qsmp.mult(Psid));
                        Matrix C = Pn.get(n - 1).mult(D.get(k + 1));
                        for (int i = 1; i < n; i++) {
                            C = C.add(1.0, Pn.get(i).mult(Qsmp).mult(Pn.get(n - i)));
                        }
                        Matrix P = Matrix.lyap(A, B, C, null);
                        Pn.put(n, P);
                        dqlProbs = Matrix.concatRows(dqlProbs, inis.mult(P), null);
                    }

                    Matrix pi = Ctmc_solve.ctmc_solve(sD);
                    double lambdak = pi.mult(D.get(k + 1)).elementSum();
                    Matrix iTerm = Matrix.negative(sD.add(-1.0, D.get(k + 1))).pinv();

                    Matrix qlProbs = Matrix.scaleMult(
                            Matrix.extractRows(dqlProbs, 0, 1, null).mult(iTerm), lambdak);
                    for (int n = 1; n < numOfQLProbs; n++) {
                        Matrix P = Matrix.extractRows(qlProbs, n - 1, n, null).mult(D.get(k + 1))
                                .add(lambdak, Matrix.extractRows(dqlProbs, n, n + 1, null)
                                        .add(-1.0, Matrix.extractRows(dqlProbs, n - 1, n, null)))
                                .mult(iTerm);
                        qlProbs = Matrix.concatRows(qlProbs, P, null);
                    }
                    result.get("ncDistr").put(k, qlProbs.sumRows().transpose());
                }
            } else if (k == K - 1) {
                if (numOfSTMoms != null) {
                    Matrix neg_Kw_inv = neg_Kw.pinv();
                    Matrix rtMoms = new Matrix(1, numOfSTMoms, numOfSTMoms);
                    for (int i = 1; i <= numOfSTMoms; i++) {
                        rtMoms.set(i - 1, CombinatoricsUtils.factorial(i)
                                * kappa.mult(Matrix.pow(neg_Kw_inv, i + 1)).mult(Bw.sumRows()).get(0));
                    }
                    result.get("stMoms").put(k, rtMoms);
                }

                if (stCdfPoints != null) {
                    Matrix neg_Kw_inv = neg_Kw.pinv();
                    Matrix rtDistr = new Matrix(1, stCdfPoints.length(), stCdfPoints.length());
                    for (int o = 0; o < stCdfPoints.length(); o++) {
                        double t = stCdfPoints.get(o);
                        Matrix Kwt = Kw.copy();
                        Kwt.scaleEq(t);
                        rtDistr.set(o, kappa.mult(neg_Kw_inv)
                                .mult(Matrix.eye(Kw.getNumRows()).add(-1.0, Maths.matrixExp(Kwt)))
                                .mult(Bw.sumRows()).get(0));
                    }
                    result.get("stDistr").put(k, rtDistr);
                }

                if (numOfQLMoms != null || numOfQLProbs != null) {
                    int Mk = (int) M.get(k);
                    int qbdSize = N * Mk;

                    Matrix L = sD.add(-1.0, D.get(k + 1)).kron(Matrix.eye(Mk)).add(1.0, I.kron(S.get(k)));
                    Matrix B = I.kron(s.get(k).mult(sigma.get(k)));
                    Matrix F = D.get(k + 1).kron(Matrix.eye(Mk));
                    Matrix L0 = sD.add(-1.0, D.get(k + 1)).kron(Matrix.eye(Mk));

                    Matrix R = QBDFundamentalMatrices.QBDFundamentalMatrices(B, L, F, precision, null, null, null).get("R");

                    Matrix p0 = Ctmc_solve.ctmc_solve(L0.add(1.0, R.mult(B)));
                    Matrix IminusR_inv = Matrix.eye(R.getNumRows()).add(-1.0, R).pinv();
                    p0 = Matrix.scaleMult(p0, 1.0 / p0.mult(IminusR_inv).elementSum());

                    if (numOfQLMoms != null) {
                        Matrix qlMoms = new Matrix(1, numOfQLMoms, numOfQLMoms);
                        for (int i = 1; i <= numOfQLMoms; i++) {
                            qlMoms.set(i - 1, CombinatoricsUtils.factorial(i)
                                    * p0.mult(Matrix.pow(R, i)).mult(Matrix.pow(IminusR_inv, i + 1)).elementSum());
                        }
                        result.get("ncMoms").put(k, MomsFromFactorialMoms.MomsFromFactorialMoms(qlMoms));
                    }

                    if (numOfQLProbs != null) {
                        Matrix qlProbs = p0.copy();
                        for (int i = 1; i < numOfQLProbs; i++) {
                            qlProbs = Matrix.concatRows(qlProbs, p0.mult(Matrix.pow(R, i)), null);
                        }
                        result.get("ncDistr").put(k, qlProbs.sumRows().transpose());
                    }
                }
            }
        }

        return result;
    }
}
