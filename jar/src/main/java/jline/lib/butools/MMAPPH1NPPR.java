package jline.lib.butools;

import jline.api.mc.Ctmc_solve;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;

import java.util.HashMap;
import java.util.Map;

/**
 * Top-level functions for MMAPPH1NPPR analysis (ported from Kotlin).
 */
public final class MMAPPH1NPPR {

    private MMAPPH1NPPR() {
    }

    public static Map<String, Map<Integer, Matrix>> MMAPPH1NPPR(MatrixCell D,
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

        for (int i = 0; i < K + 1; i++) {
            sD = sD.add(1.0, D.get(i));
        }

        Matrix M = new Matrix(1, K, K);
        Map<Integer, Matrix> s = new HashMap<Integer, Matrix>();
        for (int i = 0; i < K; i++) {
            Matrix neg_si = S.get(i).copy();
            neg_si.scaleEq(-1.0);
            s.put(i, neg_si.sumCols());
            M.set(i, M.get(i) + sigma.get(i).length());
        }
        Matrix QWMM = D0.copy();
        Matrix QWPP = new Matrix(N * (int) M.elementSum(),
                N * (int) M.elementSum(),
                N * (int) M.elementSum() * N * (int) M.elementSum());
        Matrix QWMP = new Matrix(N, N * (int) M.elementSum(), N * N * (int) M.elementSum());
        Matrix QWPM = new Matrix(N * (int) M.elementSum(), N, N * N * (int) M.elementSum());
        int kix = 0;
        for (int i = 0; i < K; i++) {
            int bs = N * (int) M.get(i);
            QWPP.insertSubMatrix(kix, kix, kix + bs, kix + bs, Matrix.eye(N).kron(S.get(i)));
            QWMP.insertSubMatrix(0, kix, QWMP.getNumRows(), kix + bs, D.get(i + 1).kron(sigma.get(i)));
            QWPM.insertSubMatrix(kix, 0, kix + bs, QWPM.getNumCols(), Matrix.eye(N).kron(s.get(i)));
            kix = kix + bs;
        }

        Map<String, Matrix> ffmRoot = FluidFundamentalMatrices.FluidFundamentalMatrices(QWPP, QWPM, QWMP, QWMM, precision, null, null);

        Matrix Kw = ffmRoot.get("K");
        Matrix Uw = ffmRoot.get("U");
        Matrix neg_Kw = Kw.copy();
        neg_Kw.scaleEq(-1.0);

        Matrix Ua = Matrix.ones(N, 1).add(2.0, QWMP.mult(neg_Kw.pinv()).sumRows());
        Matrix A_ls = Matrix.concatColumns(Uw, Ua, null).transpose();
        Matrix b_ls = Matrix.concatColumns(new Matrix(1, N, 0), Matrix.singleton(1.0), null).transpose();
        Matrix pm = A_ls.leftMatrixDivide(b_ls).transpose();

        double ro = ((1 - pm.elementSum()) / 2) / (pm.elementSum() + (1 - pm.elementSum()) / 2);
        Matrix kappa = pm.copy();
        kappa.scaleEq(1 / pm.elementSum());

        Matrix pi = Ctmc_solve.ctmc_solve(sD);
        Matrix lambda = new Matrix(K, 1, 1);
        for (int i = 0; i < K; i++) {
            lambda.set(i, pi.mult(D.get(i + 1)).elementSum());
        }

        Map<Integer, Matrix> Psiw = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> Qwmp = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> Qwzp = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> Qwpp = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> Qwmz = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> Qwpz = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> Qwzz = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> Qwmm = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> Qwpm = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> Qwzm = new HashMap<Integer, Matrix>();

        for (int k = 0; k < K; k++) {
            double Mlo = (k == 0) ? 0.0 : M.sumSubMatrix(0, 1, 0, k);
            double Mhi = M.elementSum() - Mlo;
            int a = (int) (N * Mlo * Mhi + N * Mhi);
            Matrix Qkwpp = new Matrix(a, a, a * a);
            Matrix Qkwpz = new Matrix(a, N * (int) Mlo, a * N * (int) Mlo);
            Matrix Qkwpm = new Matrix(a, N, a * N);
            Matrix Qkwmz = new Matrix(N, N * (int) Mlo);
            Matrix Qkwmp = new Matrix(N, a, a * N);
            Matrix Dlo = D0.copy();
            for (int i = 0; i < k; i++) {
                Dlo = Dlo.add(1.0, D.get(i + 1));
            }
            Matrix Qkwmm = Dlo;
            Matrix Qkwzp = new Matrix(N * (int) Mlo, a, N * (int) Mlo * a);
            Matrix Qkwzm = new Matrix(N * (int) Mlo, N, N * (int) Mlo * N);
            Matrix Qkwzz = new Matrix(N * (int) Mlo, N * (int) Mlo, N * (int) Mlo * N * (int) Mlo);
            kix = 0;
            for (int i = k; i < K; i++) {
                int kix2 = 0;
                for (int j = 0; j < k; j++) {
                    int bs = (int) (N * M.get(j) * M.get(i));
                    int bs2 = N * (int) M.get(j);
                    Qkwpp.insertSubMatrix(kix, kix, kix + bs, kix + bs,
                            Matrix.eye(N).kron(Matrix.eye((int) M.get(j)).kron(S.get(i))));
                    Qkwpz.insertSubMatrix(kix, kix2, kix + bs, kix2 + bs2,
                            Matrix.eye(N).kron(Matrix.eye((int) M.get(j)).kron(s.get(i))));
                    Qkwzp.insertSubMatrix(kix2, kix, kix2 + bs2, kix + bs,
                            D.get(i + 1).kron(Matrix.eye((int) M.get(j)).kron(sigma.get(i))));
                    kix = kix + bs;
                    kix2 = kix2 + bs2;
                }
            }
            for (int i = k; i < K; i++) {
                int bs = N * (int) M.get(i);
                Qkwpp.insertSubMatrix(kix, kix, kix + bs, kix + bs, Matrix.eye(N).kron(S.get(i)));
                Qkwpm.insertSubMatrix(kix, 0, kix + bs, Qkwpm.getNumCols(), Matrix.eye(N).kron(s.get(i)));
                Qkwmp.insertSubMatrix(0, kix, Qkwmp.getNumRows(), kix + bs, D.get(i + 1).kron(sigma.get(i)));
                kix = kix + bs;
            }
            kix = 0;
            for (int j = 0; j < k; j++) {
                int bs = N * (int) M.get(j);
                Qkwzz.insertSubMatrix(kix, kix, kix + bs, kix + bs,
                        Dlo.kron(Matrix.eye((int) M.get(j))).add(1.0, Matrix.eye(N).kron(S.get(j))));
                Qkwzm.insertSubMatrix(kix, 0, kix + bs, Qkwzm.getNumCols(), Matrix.eye(N).kron(s.get(j)));
                kix = kix + bs;
            }
            Matrix Psikw;
            if (Mlo > 0) {
                Matrix neg_Qkwzz = Qkwzz.copy();
                neg_Qkwzz.scaleEq(-1.0);
                Psikw = FluidFundamentalMatrices.FluidFundamentalMatrices(
                        Qkwpp.add(1.0, Qkwpz.mult(neg_Qkwzz.pinv()).mult(Qkwzp)),
                        Qkwpm.add(1.0, Qkwpz.mult(neg_Qkwzz.pinv()).mult(Qkwzm)),
                        Qkwmp, Qkwmm, precision, null, null).get("P");
            } else {
                Psikw = FluidFundamentalMatrices.FluidFundamentalMatrices(
                        Qkwpp, Qkwpm, Qkwmp, Qkwmm, precision, null, null).get("P");
            }
            Psiw.put(k, Psikw);
            Qwzp.put(k, Qkwzp);
            Qwmp.put(k, Qkwmp);
            Qwpp.put(k, Qkwpp);
            Qwmz.put(k, Qkwmz);
            Qwpz.put(k, Qkwpz);
            Qwzz.put(k, Qkwzz);
            Qwmm.put(k, Qkwmm);
            Qwpm.put(k, Qkwpm);
            Qwzm.put(k, Qkwzm);
        }

        double lambdaS = lambda.elementSum();
        Map<Integer, Matrix> phi = new HashMap<Integer, Matrix>();
        Matrix neg_D0 = D0.copy();
        neg_D0.scaleEq(-1.0);
        Matrix phi0 = kappa.mult(neg_D0);
        phi0.scaleEq((1 - ro) / lambdaS);
        phi.put(0, phi0);

        Map<Integer, Matrix> q0 = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> qL = new HashMap<Integer, Matrix>();
        q0.put(0, new Matrix(0, 0, 0));
        qL.put(0, new Matrix(0, 0, 0));

        for (int k = 0; k < K - 1; k++) {
            Matrix sDk = D.get(0);
            for (int j = 0; j <= k; j++) {
                sDk = sDk.add(1.0, D.get(j + 1));
            }
            double pk = 0.0;
            for (int j = 0; j <= k; j++) {
                pk = pk + lambda.get(j);
            }
            pk = pk / lambdaS - (1 - ro) * kappa.mult(sDk.sumRows()).get(0) / lambdaS;
            Matrix Qwzpk = Qwzp.get(k + 1);
            int vix = 0;
            Map<Integer, Matrix> Ak = new HashMap<Integer, Matrix>();
            for (int ii = 0; ii <= k; ii++) {
                int bs = (int) (N * M.get(ii));
                Matrix V1 = Matrix.extractRows(Qwzpk, vix, vix + bs, null);
                Matrix a = sDk.kron(Matrix.eye((int) M.get(ii)));
                a.scaleEq(-1.0);
                a.add(-1.0, I.kron(S.get(ii)));
                Ak.put(ii, I.kron(sigma.get(ii)).mult(a.pinv()).mult(I.kron(s.get(ii))).add(1.0, V1.mult(Psiw.get(k + 1))));
                vix = vix + bs;
            }
            Matrix Qwmpk = Qwmp.get(k + 1);
            Matrix Bk = Qwmpk.mult(Psiw.get(k + 1));
            Matrix ztag = phi.get(0).mult(neg_D0.pinv().mult(D.get(k + 1)).mult(Ak.get(k))
                    .add(-1.0, Ak.get(0)).add(1.0, neg_D0.pinv().mult(Bk)));
            for (int i = 0; i < k; i++) {
                ztag = ztag.add(1.0, phi.get(i + 1).mult(Ak.get(i).add(-1.0, Ak.get(i + 1))))
                        .add(1.0, phi.get(0).mult(neg_D0.pinv()).mult(D.get(i + 1)).mult(Ak.get(i)));
            }
            Matrix Mx = Matrix.eye(Ak.get(k).getNumCols()).add(-1.0, Ak.get(k));
            Mx.insertSubMatrix(0, 0, Mx.getNumRows(), 1, Matrix.ones(N, 1));
            phi.put(k + 1, Matrix.concatColumns(Matrix.singleton(pk),
                    Matrix.extractRows(ztag.transpose(), 1, ztag.length(), null).transpose(),
                    null).mult(Mx.pinv()));
            q0.put(k + 1, phi.get(0).mult(neg_D0.pinv()));
            qL.put(k + 1, new Matrix(0, 0, 0));
            for (int ii = 0; ii <= k; ii++) {
                Matrix a = sDk.kron(Matrix.eye((int) M.get(ii)));
                a.scaleEq(-1.0);
                Matrix qLii = phi.get(ii + 1).add(-1.0, phi.get(ii)).add(1.0, phi.get(0).mult(neg_D0.pinv()).mult(D.get(ii + 1)))
                        .mult(I.kron(sigma.get(ii))).mult(a.add(-1.0, I.kron(S.get(ii))).pinv());
                qL.put(k + 1, Matrix.concatColumns(qL.get(k + 1), qLii, null));
            }
        }

        Map<String, Map<Integer, Matrix>> result = new HashMap<String, Map<Integer, Matrix>>();
        if (numOfSTMoms != null) result.put("stMoms", new HashMap<Integer, Matrix>());
        if (stCdfPoints != null) result.put("stDistr", new HashMap<Integer, Matrix>());
        if (numOfQLMoms != null) result.put("ncMoms", new HashMap<Integer, Matrix>());
        if (numOfQLProbs != null) result.put("ncDistr", new HashMap<Integer, Matrix>());

        for (int g = 0; g < classes.length(); g++) {
            int k = (int) classes.get(g);
            Matrix sD0k = D0;
            for (int i = 0; i < k; i++) {
                sD0k = sD0k.add(1.0, D.get(i + 1));
            }
            if (k < K - 1) {
                double Mlo_k = (k == 0) ? 0.0 : M.sumSubMatrix(0, 1, 0, k);
                if (Mlo_k > 0) {
                    Matrix neg_Qwzz_k = Qwzz.get(k).copy();
                    neg_Qwzz_k.scaleEq(-1.0);
                    Kw = Qwpp.get(k).add(1.0, Qwpz.get(k).mult(neg_Qwzz_k.pinv()).mult(Qwzp.get(k)))
                            .add(1.0, Psiw.get(k).mult(Qwmp.get(k)));
                } else {
                    Kw = Qwpp.get(k).add(1.0, Psiw.get(k).mult(Qwmp.get(k)));
                }
                Matrix BM = new Matrix(0, 0, 0);
                Matrix CM = new Matrix(0, N, 0);
                Matrix DM = new Matrix(0, 0, 0);
                for (int i = 0; i < k; i++) {
                    BM = BM.createBlockDiagonal(I.kron(S.get(i)));
                    DM = DM.createBlockDiagonal(D.get(k + 1).kron(Matrix.eye((int) M.get(i))));
                    CM = Matrix.concatRows(CM, I.kron(s.get(i)), null);
                }
                Matrix Kwu;
                if (Mlo_k > 0) {
                    Kwu = Matrix.concatRows(Matrix.concatColumns(Kw,
                                    Qwpz.get(k).add(1.0, Psiw.get(k).mult(Qwmz.get(k))).mult(Matrix.negative(Qwzz.get(k)).pinv()).mult(DM),
                                    null),
                            Matrix.concatColumns(new Matrix(BM.getNumRows(), Kw.getNumRows(), 0), BM, null), null);
                } else {
                    Kwu = Kw;
                }
                Matrix Bwu;
                if (Mlo_k > 0) {
                    Bwu = Matrix.concatRows(Psiw.get(k).mult(D.get(k + 1)), CM, null);
                } else {
                    Bwu = Psiw.get(k).mult(D.get(k + 1));
                }
                Matrix iniw;
                Matrix pwu;
                if (k > 0) {
                    iniw = Matrix.concatColumns(q0.get(k).mult(Qwmp.get(k)).add(1.0, qL.get(k).mult(Qwzp.get(k))),
                            qL.get(k).mult(DM), null);
                    pwu = q0.get(k).mult(D.get(k + 1));
                } else {
                    iniw = pm.mult(Qwmp.get(k));
                    pwu = pm.mult(D.get(k + 1));
                }
                double norm = pwu.elementSum() + iniw.mult(Matrix.negative(Kwu).pinv()).mult(Bwu).elementSum();
                pwu.scaleEq(1 / norm);
                iniw.scaleEq(1 / norm);
                int KN = Kwu.getNumRows();
                double sizeHigh = M.sumSubMatrix(0, 1, k + 1, K);
                Matrix Qspp = new Matrix((int) (KN + N * sizeHigh),
                        (int) (KN + N * sizeHigh),
                        (int) Math.pow(KN + N * sizeHigh, 2.0));
                Matrix Qspm = new Matrix((int) (KN + N * sizeHigh), N, (int) (N * (KN + N * sizeHigh)));
                Matrix Qsmp = new Matrix(N, (int) (KN + N * sizeHigh), (int) (N * (KN + N * sizeHigh)));
                Matrix Qsmm = sD0k.add(1.0, D.get(k + 1));
                kix = 0;
                for (int i = k + 1; i < K; i++) {
                    int bs = N * (int) M.get(i);
                    Qspp.insertSubMatrix(KN + kix, KN + kix, KN + kix + bs, KN + kix + bs, I.kron(S.get(i)));
                    Qspm.insertSubMatrix(KN + kix, 0, KN + kix + bs, Qspm.getNumCols(), I.kron(s.get(i)));
                    Qsmp.insertSubMatrix(0, KN + kix, Qsmp.getNumRows(), KN + kix + bs, D.get(i + 1).kron(sigma.get(i)));
                    kix = kix + bs;
                }
                Qspp.insertSubMatrix(0, 0, KN, KN, Kwu);
                Qspm.insertSubMatrix(0, 0, KN, Qspm.getNumCols(), Bwu);
                Matrix inis = Matrix.concatColumns(iniw, new Matrix(1, (int) (N * sizeHigh)), null);

                Matrix Psis = FluidFundamentalMatrices.FluidFundamentalMatrices(Qspp, Qspm, Qsmp, Qsmm, precision, null, null).get("P");

                if (numOfSTMoms != null) {
                    Map<Integer, Matrix> Pn = new HashMap<Integer, Matrix>();
                    Pn.put(0, Psis);
                    Matrix wtMoms = new Matrix(1, numOfSTMoms, numOfSTMoms);
                    for (int n = 1; n <= numOfSTMoms; n++) {
                        Matrix A = Qspp.add(1.0, Psis.mult(Qsmp));
                        Matrix B = Qsmm.add(1.0, Qsmp.mult(Psis));
                        Matrix C = Pn.get(n - 1).copy();
                        C.scaleEq((double) (-2 * n));
                        long bino = 1;
                        for (int i = 1; i <= n - 1; i++) {
                            bino = bino * (n - i + 1) / i;
                            C = C.add((double) bino, Pn.get(i).mult(Qsmp).mult(Pn.get(n - i)));
                        }
                        Matrix P = Matrix.lyap(A, B, C, null);
                        Pn.put(n, P);
                        wtMoms.set(n - 1, inis.mult(P).elementSum() * Math.pow(-1.0, (double) n) / Math.pow(2.0, (double) n));
                    }
                    Map<Integer, Matrix> Pnr = new HashMap<Integer, Matrix>();
                    Pnr.put(0, Matrix.scaleMult(sigma.get(k), inis.mult(Pn.get(0)).elementSum()));
                    Matrix rtMoms = new Matrix(1, numOfSTMoms, numOfSTMoms);
                    Matrix negSk = S.get(k).copy();
                    negSk.scaleEq(-1.0);
                    for (int n = 1; n <= numOfSTMoms; n++) {
                        Matrix term1 = Pnr.get(n - 1).mult(negSk.pinv());
                        term1.scaleEq((double) n);
                        double coeff = Math.pow(-1.0, (double) n) / Math.pow(2.0, (double) n) * inis.mult(Pn.get(n)).elementSum();
                        Matrix P = term1.add(coeff, sigma.get(k));
                        Pnr.put(n, P);
                        rtMoms.set(n - 1, P.elementSum() + pwu.elementSum() * CombinatoricsUtils.factorial(n) *
                                sigma.get(k).mult(Matrix.pow(negSk.pinv(), n)).elementSum());
                    }
                    result.get("stMoms").put(k, rtMoms);
                }

                if (stCdfPoints != null) {
                    Matrix res = new Matrix(1, 0, 0);
                    for (int o = 0; o < stCdfPoints.length(); o++) {
                        int t = (int) stCdfPoints.get(o);
                        int L = erlMaxOrder;
                        double lambdae = (double) L / (double) t / 2.0;
                        Matrix Psie = FluidFundamentalMatrices.FluidFundamentalMatrices(
                                Qspp.add(-lambdae, Matrix.eye(Qspp.getNumRows())),
                                Qspm, Qsmp,
                                Qsmm.add(-lambdae, Matrix.eye(Qsmm.getNumRows())),
                                precision, null, null).get("P");
                        Map<Integer, Matrix> Pn = new HashMap<Integer, Matrix>();
                        Pn.put(0, Psis);
                        double pr = (pwu.elementSum() + inis.mult(Psie).elementSum()) *
                                (1 - sigma.get(k).mult(Matrix.pow(Matrix.eye(S.get(k).getNumRows()).pinv()
                                        .add(-1 / 2.0 / lambdae, S.get(k)), L)).elementSum());
                        for (int n = 0; n < L; n++) {
                            Matrix A = Qspp.add(1.0, Psie.mult(Qsmp)).add(-lambdae, Matrix.eye(Qspp.length()));
                            Matrix B = Qsmm.add(1.0, Qsmp.mult(Psie)).add(-lambdae, Matrix.eye(Qsmm.length()));
                            Matrix C = Pn.get(n);
                            C.scaleEq(2 * lambdae);
                            for (int i = 0; i < n - 1; i++) {
                                C = C.add(1.0, Pn.get(i + 1).mult(Qsmp).mult(Pn.get(n - i + 1)));
                            }
                            Matrix P = Matrix.lyap(A, B, C, null);
                            Pn.put(n + 1, P);
                            pr = pr + inis.mult(P).elementSum() * (1 - sigma.get(k).mult(Matrix.pow(Matrix.eye(S.get(k).length())
                                    .pinv().add(-1 / 2.0 / lambdae, S.get(k)), L - n)).elementSum());
                        }
                        res = Matrix.concatColumns(res, Matrix.singleton(pr), null);
                    }
                    result.get("stDistr").put(k, res);
                }
                if (numOfQLMoms != null || numOfQLProbs != null) {
                    Matrix W = Matrix.negative(sD.add(-1.0, D.get(k + 1)).kron(Matrix.eye((int) M.get(k))))
                            .add(-1.0, I.kron(S.get(k))).pinv()
                            .mult(D.get(k + 1).kron(Matrix.eye((int) M.get(k))));
                    Matrix iW = Matrix.eye(W.getNumCols()).add(-1.0, W).pinv();
                    Matrix w = Matrix.eye(N).kron(sigma.get(k));
                    Matrix omega = Matrix.negative(sD.add(-1.0, D.get(k + 1)).kron(Matrix.eye((int) M.get(k))))
                            .add(-1.0, I.kron(S.get(k))).pinv()
                            .mult(I.kron(s.get(k)));
                    if (numOfQLMoms != null) {
                        Map<Integer, Matrix> Psii = new HashMap<Integer, Matrix>();
                        Psii.put(0, Psis);
                        Map<Integer, Matrix> QLDPn = new HashMap<Integer, Matrix>();
                        QLDPn.put(0, inis.mult(Psii.get(0)).mult(w).mult(iW));
                        for (int n = 1; n <= numOfQLMoms; n++) {
                            Matrix A = Qspp.add(1.0, Psis.mult(Qsmp));
                            Matrix B = Qsmm.add(1.0, Qsmp.mult(Psis));
                            Matrix C = Psii.get(n - 1).mult(D.get(k + 1));
                            C.scaleEq((double) n);
                            long bino = 1;
                            for (int i = 1; i <= n - 1; i++) {
                                bino = bino * (n - i + 1) / i;
                                C = C.add((double) bino, Psii.get(i).mult(Qsmp).mult(Psii.get(n - i)));
                            }
                            Matrix P = Matrix.lyap(A, B, C, null);
                            Psii.put(n, P);
                            Matrix qlDTerm = QLDPn.get(n - 1).mult(iW).mult(W);
                            qlDTerm.scaleEq((double) n);
                            QLDPn.put(n, qlDTerm.add(1.0, inis.mult(P).mult(w).mult(iW)));
                        }
                        for (int n = 0; n <= numOfQLMoms; n++) {
                            QLDPn.put(n, QLDPn.get(n).add(1.0,
                                    pwu.mult(w).mult(Matrix.pow(iW, n + 1)).mult(Matrix.pow(W, n))).mult(omega));
                        }
                        Map<Integer, Matrix> QLPn = new HashMap<Integer, Matrix>();
                        QLPn.put(0, pi);
                        Matrix qlMOms = new Matrix(1, numOfQLMoms, numOfQLMoms);
                        Matrix iTerm = Matrix.ones(N, 1).mult(pi).add(-1.0, sD).pinv();
                        for (int n = 1; n <= numOfQLMoms; n++) {
                            double sumP = QLDPn.get(n).elementSum() + n * (QLDPn.get(n - 1)
                                    .add(-1.0 / lambda.get(k), QLPn.get(n - 1).mult(D.get(k + 1)))
                                    .mult(iTerm).mult(D.get(k + 1).sumRows())).get(0);
                            Matrix P = Matrix.scaleMult(pi, sumP)
                                    .add((double) n, QLPn.get(n - 1).mult(D.get(k + 1)).add(-lambda.get(k), QLDPn.get(n - 1)).mult(iTerm));
                            QLPn.put(n, P);
                            qlMOms.set(n - 1, P.elementSum());
                        }
                        qlMOms = MomsFromFactorialMoms.MomsFromFactorialMoms(qlMOms);
                        result.get("ncMoms").put(k, qlMOms);
                    }
                    if (numOfQLProbs != null) {
                        Matrix Psid = FluidFundamentalMatrices.FluidFundamentalMatrices(Qspp, Qspm, Qsmp, sD0k, precision, null, null).get("P");
                        Map<Integer, Matrix> Pn = new HashMap<Integer, Matrix>();
                        Pn.put(0, Psid);
                        Matrix XDn = inis.mult(Psid).mult(w);
                        Matrix dqlProbs = XDn.add(1.0, pwu.mult(w)).mult(omega);
                        for (int n = 1; n < numOfQLProbs; n++) {
                            Matrix A = Qspp.add(1.0, Psid.mult(Qsmp));
                            Matrix B = sD0k.add(1.0, Qsmp.mult(Psid));
                            Matrix C = Pn.get(n - 1).mult(D.get(k + 1));
                            for (int i = 1; i <= n - 1; i++) {
                                C = C.add(1.0, Pn.get(i).mult(Qsmp).mult(Pn.get(n - i)));
                            }
                            Matrix P = Matrix.lyap(A, B, C, null);
                            Pn.put(n, P);
                            XDn = XDn.mult(W).add(1.0, inis.mult(P).mult(w));
                            dqlProbs = Matrix.concatRows(dqlProbs,
                                    XDn.add(1.0, pwu.mult(w).mult(Matrix.pow(W, n))).mult(omega), null);
                        }
                        Matrix iTerm = Matrix.negative(sD.add(-1.0, D.get(k + 1))).pinv();
                        Matrix qlProbs = Matrix.scaleMult(Matrix.extractRows(dqlProbs, 0, 1, null).mult(iTerm), lambda.get(k));
                        for (int n = 1; n < numOfQLProbs; n++) {
                            Matrix P = Matrix.extractRows(qlProbs, n - 1, n, null).mult(D.get(k + 1)).add(1.0,
                                    Matrix.scaleMult(Matrix.extractRows(dqlProbs, n, n + 1, null)
                                            .add(-1.0, Matrix.extractRows(dqlProbs, n - 1, n, null)), lambda.get(k))).mult(iTerm);
                            qlProbs = Matrix.concatRows(qlProbs, P, null);
                        }
                        result.get("ncDistr").put(k, qlProbs.sumRows().transpose());
                    }
                }
            } else if (k == K - 1) {
                if (numOfSTMoms != null || stCdfPoints != null) {
                    Kw = Qwpp.get(k).add(1.0, Qwpz.get(k).mult(Matrix.negative(Qwzz.get(k)).pinv()).mult(Qwzp.get(k)))
                            .add(1.0, Psiw.get(k).mult(Qwmp.get(k)));
                    Matrix AM = new Matrix(0, 0, 0);
                    Matrix BM = new Matrix(0, 0, 0);
                    Matrix CM = new Matrix(0, s.get(0).getNumCols(), 0);
                    Matrix DM = new Matrix(0, 0, 0);
                    for (int i = 0; i < k; i++) {
                        AM = AM.createBlockDiagonal(new Matrix(N, 1).kron(Matrix.eye((int) M.get(i)).kron(s.get(k))));
                        BM = BM.createBlockDiagonal(S.get(i));
                        CM = Matrix.concatRows(CM, s.get(i), null);
                        DM = DM.createBlockDiagonal(D.get(k + 1).kron(Matrix.eye((int) M.get(i))));
                    }
                    Matrix Z = Matrix.concatRows(Matrix.concatColumns(Kw,
                                    Matrix.concatRows(AM, new Matrix(N * (int) M.get(k), AM.getNumCols()), null), null),
                            Matrix.concatColumns(new Matrix(BM.getNumRows(), Kw.getNumCols()), BM, null), null);
                    Matrix z = Matrix.concatRows(
                            Matrix.concatRows(new Matrix(AM.getNumRows(), 1, 0), Matrix.ones(N, 1).kron(s.get(k)), null),
                            CM, null);
                    Matrix iniw = Matrix.concatColumns(q0.get(k).mult(Qwmp.get(k)).add(1.0, qL.get(k).mult(Qwzp.get(k))),
                            new Matrix(1, BM.getNumRows()), null);
                    Matrix zeta = Matrix.scaleMult(iniw, 1 / iniw.mult(Matrix.negative(Z).pinv()).mult(z).elementSum());
                    if (numOfSTMoms != null) {
                        Matrix rtMomsH = new Matrix(1, numOfSTMoms, numOfSTMoms);
                        for (int i = 0; i < numOfSTMoms; i++) {
                            rtMomsH.set(i, CombinatoricsUtils.factorial(i + 1) *
                                    zeta.mult(Matrix.pow(Matrix.negative(Z).pinv(), i + 2))
                                            .mult(z).elementSum());
                        }
                        result.get("stMoms").put(k, rtMomsH);
                    }
                    if (stCdfPoints != null) {
                        Matrix rtDistr = zeta.mult(Matrix.negative(Z).pinv())
                                .mult(Matrix.eye(Z.getNumCols()).add(-1.0, Z.mult(Maths.matrixExp(Matrix.scaleMult(Z, stCdfPoints.get(0))))))
                                .mult(z);
                        for (int i = 1; i < stCdfPoints.length(); i++) {
                            rtDistr = Matrix.concatColumns(rtDistr,
                                    zeta.mult(Matrix.negative(Z).pinv())
                                            .mult(Matrix.eye(Z.getNumCols()).add(-1.0, Z.mult(Maths.matrixExp(Matrix.scaleMult(Z, stCdfPoints.get(i))))))
                                            .mult(z), null);
                        }
                        result.get("stDistr").put(k, rtDistr);
                    }
                }

                if (numOfQLMoms != null || numOfQLProbs != null) {
                    Matrix L = new Matrix(N * (int) M.elementSum(),
                            N * (int) M.elementSum(),
                            (int) FastMath.pow(N * M.elementSum(), 2));
                    Matrix B = new Matrix(N * (int) M.elementSum(),
                            N * (int) M.elementSum(),
                            (int) FastMath.pow(N * M.elementSum(), 2));
                    Matrix F = new Matrix(N * (int) M.elementSum(),
                            N * (int) M.elementSum(),
                            (int) FastMath.pow(N * M.elementSum(), 2));
                    kix = 0;
                    for (int i = 0; i < K; i++) {
                        int bs = N * (int) M.get(i);
                        F.insertSubMatrix(kix, kix, kix + bs, kix + bs, D.get(k + 1).kron(Matrix.eye((int) M.get(i))));
                        L.insertSubMatrix(kix, kix, kix + bs, kix + bs,
                                sD0k.kron(Matrix.eye((int) M.get(i))).add(1.0, I.kron(S.get(i))));
                        if (i < K - 1) {
                            L.insertSubMatrix(kix, N * (int) M.sumSubMatrix(0, 1, 0, k),
                                    kix + bs, L.getNumCols(),
                                    I.kron(s.get(i).mult(sigma.get(k))));
                        } else {
                            B.insertSubMatrix(kix, N * (int) M.sumSubMatrix(0, 1, 0, k),
                                    kix + bs, B.getNumCols(),
                                    I.kron(s.get(i).mult(sigma.get(k))));
                        }
                        kix = kix + bs;
                    }
                    Matrix R = QBDFundamentalMatrices.QBDFundamentalMatrices(B, L, F, precision, null, null, null).get("R");
                    Matrix P0 = Matrix.concatColumns(qL.get(k), q0.get(k).mult(I.kron(sigma.get(k))), null);
                    P0.scaleEq(1 / P0.mult(Matrix.eye(R.getNumRows()).add(-1.0, R).pinv()).elementSum());

                    if (numOfQLMoms != null) {
                        Matrix qlMoms = new Matrix(1, numOfQLMoms, numOfQLMoms);
                        for (int i = 0; i < numOfQLMoms; i++) {
                            qlMoms.set(i, Matrix.scaleMult(P0.mult(Matrix.pow(R, i + 1))
                                            .mult(Matrix.pow(Matrix.eye(R.getNumRows()).add(-1.0, R), i + 2).pinv()),
                                    (double) CombinatoricsUtils.factorial(i + 1)).elementSum());
                        }
                        result.get("ncMoms").put(k, MomsFromFactorialMoms.MomsFromFactorialMoms(qlMoms));
                    }

                    if (numOfQLProbs != null) {
                        Matrix qlProbs = P0.copy();
                        for (int i = 1; i < numOfQLProbs; i++) {
                            qlProbs = Matrix.concatRows(qlProbs, P0.mult(Matrix.pow(R, i)), null);
                        }
                        result.get("ncDistr").put(k, qlProbs.sumRows());
                    }
                }
            }
        }
        return result;
    }
}
