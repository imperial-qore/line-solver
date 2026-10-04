/*
 * Ported from BUTools-family fluid tools (G. Horvath).
 *
 * Performance measures of a continuous-time fluid priority queue.
 *
 * Reference: G. Horvath, "Efficient analysis of the MMAP[K]/PH[K]/1 priority
 * queue", European Journal of Operational Research, 246(1):128-139, 2015.
 */
package jline.lib.butools.queues;

import java.util.ArrayList;
import java.util.List;

import jline.lib.butools.FluidFundamentalMatrices;
import jline.lib.butools.mam.FluidStationaryDistr;
import jline.lib.butools.mam.GeneralFluidSolution;
import jline.lib.butools.mam.GeneralFluidSolve;
import jline.lib.butools.mc.CTMCSolve;
import jline.util.matrix.Matrix;

import static jline.lib.butools.mam.FluidTools.*;

public final class FluidPrioQueue {
    private FluidPrioQueue() {}

    /**
     * Returns performance measures of a continuous-time fluid priority queue.
     * Q (N,N) is the modulating generator, R (K,N) the per-class input rates, d
     * the constant service rate. measures alternates a name ("flMoms","stMoms",
     * "flDistr","stDistr") and a value (Integer number of moments, or double[]
     * points). classes (1-based) selects the priority classes to analyze.
     *
     * @return one double[] per (class, measure) in class-major order.
     */
    public static List<double[]> fluidPrioQueue(Matrix Q, Matrix R, double d, int[] classes,
                                                double prec, int erlMaxOrder, Object... measures) {
        int K = R.getNumRows();
        int N = Q.getNumRows();
        if (classes == null) {
            classes = new int[K];
            for (int i = 0; i < K; i++) classes[i] = i + 1;
        }
        Matrix pi = CTMCSolve.ctmcSolve(Q);
        Matrix lambda = pi.mult(R.transpose()); // 1 x K

        List<double[]> Ret = new ArrayList<double[]>();
        for (int kOne : classes) {
            int k = kOne - 1;
            Matrix RsumHi = sumRowsRange(R, k, K);         // 1 x N
            GeneralFluidSolution gs = GeneralFluidSolve.generalFluidSolve(
                    Q, diagFrom(RsumHi).scale(1.0 / d).sub(Matrix.eye(N)), null, prec);
            Matrix mass0 = gs.getMass0(), ini = gs.getIni(), Km = gs.getK(), clo = gs.getClo();
            int KN = Km.getNumRows();
            double lamk = lambda.get(0, k);
            Matrix Rk = diagFrom(rowVec(R, k));
            Matrix clok = clo.mult(Rk).scale(1.0 / lamk);

            Matrix Delta = diagFrom(Km.scale(-1.0).pinv().mult(rowSums(clok)));
            Km = Delta.pinv().mult(Km).mult(Delta);
            clok = Delta.pinv().mult(clok);
            ini = ini.mult(Delta);

            if (k < K - 1) {
                Matrix RsumLo = sumRowsRange(R, k + 1, K);
                Matrix F = block(Km, clok, Matrix.zeros(N, KN), Q);
                Matrix Clo = block(Matrix.eye(KN), Matrix.zeros(KN, N),
                        Matrix.zeros(N, KN), diagFrom(RsumLo).scale(1.0 / d).sub(Matrix.eye(N)));
                Matrix inis = hstack(ini, Matrix.zeros(1, N));
                for (int mi = 0; mi < measures.length; mi += 2) {
                    String opt = (String) measures[mi];
                    Object val = measures[mi + 1];
                    if (opt.equals("stMoms")) {
                        int nm = (Integer) val;
                        Matrix D = block(Matrix.zeros(KN, KN), Matrix.zeros(KN, N), Matrix.zeros(N, KN), Matrix.eye(N));
                        List<Matrix> Tmp = busyPeriodRewardMoms(F, Clo, D, nm, prec);
                        double[] st = new double[nm];
                        for (int i = 1; i < Tmp.size(); i++)
                            st[i - 1] = Math.pow(-1, i) * inis.mult(Tmp.get(i)).elementSum();
                        Ret.add(st);
                    } else if (opt.equals("flMoms")) {
                        int nm = (Integer) val;
                        Matrix D = block(Matrix.zeros(KN, KN), Matrix.zeros(KN, N), Matrix.zeros(N, KN), Rk);
                        List<Matrix> FLDPn = busyPeriodRewardMoms(F, Clo, D, nm, prec);
                        List<Matrix> FLDv = new ArrayList<Matrix>();
                        for (int i = 0; i < FLDPn.size(); i++) {
                            Matrix Xi = inis.mult(sub(FLDPn.get(i), range(FLDPn.get(i).getNumRows()), shiftRange(KN, N)));
                            if (i == 0) Xi = Xi.add(mass0.mult(Rk).scale(1.0 / lamk));
                            FLDv.add(Xi);
                        }
                        double[] fld = new double[nm];
                        List<Matrix> FLPn = new ArrayList<Matrix>();
                        FLPn.add(pi);
                        Matrix iTerm = Matrix.ones(N, 1).mult(pi).sub(Q).pinv();
                        Matrix RkCol = rowVec(R, k).transpose(); // (N,1)
                        for (int n = 1; n <= nm; n++) {
                            double sumP = FLDv.get(n).elementSum()
                                    + n * FLDv.get(n - 1).scale(-1.0).add(FLPn.get(n - 1).mult(Rk).scale(1.0 / lamk))
                                    .mult(iTerm).mult(RkCol).get(0, 0);
                            Matrix P = pi.scale(sumP).add(
                                    FLPn.get(n - 1).mult(Rk).scale(-1.0).add(FLDv.get(n - 1).scale(lamk)).mult(iTerm).scale(n));
                            FLPn.add(P);
                            fld[n - 1] = Math.pow(-1, n) * P.elementSum();
                        }
                        Ret.add(fld);
                    } else if (opt.equals("stDistr")) {
                        double[] pts = (double[]) val;
                        Matrix D = block(Matrix.zeros(KN, KN), Matrix.zeros(KN, N), Matrix.zeros(N, KN), Matrix.eye(N));
                        double[] res = new double[pts.length];
                        for (int xi = 0; xi < pts.length; xi++) {
                            Matrix[] br = busyPeriodRewardDistr(F, Clo, D, pts[xi], prec, erlMaxOrder);
                            res[xi] = mass0.mult(Rk).scale(1.0 / lamk).elementSum() + inis.mult(br[0]).elementSum();
                        }
                        Ret.add(res);
                    } else if (opt.equals("flDistr")) {
                        double[] pts = (double[]) val;
                        Matrix D = block(Matrix.zeros(KN, KN), Matrix.zeros(KN, N), Matrix.zeros(N, KN), Rk);
                        double[] res = new double[pts.length];
                        for (int xi = 0; xi < pts.length; xi++) {
                            double nu = erlMaxOrder / pts[xi];
                            Matrix[] br = busyPeriodRewardDistrFull(F, Clo, D, pts[xi], prec, erlMaxOrder);
                            List<Matrix> Psix = splitList(br);
                            Matrix Psiy = mass0.mult(Rk).scale(1.0 / lamk)
                                    .add(inis.mult(sub(Psix.get(0), range(Psix.get(0).getNumRows()), shiftRange(KN, N))))
                                    .scale(lamk * nu).mult(Rk.scale(nu).sub(Q).pinv());
                            for (int i = 1; i < Psix.size(); i++) {
                                Psiy = inis.mult(sub(Psix.get(i), range(Psix.get(i).getNumRows()), shiftRange(KN, N))).scale(lamk)
                                        .add(Psiy.mult(Rk)).scale(nu).mult(Rk.scale(nu).sub(Q).pinv());
                            }
                            res[xi] = Psiy.elementSum();
                        }
                        Ret.add(res);
                    } else {
                        throw new IllegalArgumentException("FluidPrioQueue: Unknown parameter " + opt);
                    }
                }
            } else { // k == K (lowest priority)
                for (int mi = 0; mi < measures.length; mi += 2) {
                    String opt = (String) measures[mi];
                    Object val = measures[mi + 1];
                    if (opt.equals("stMoms")) {
                        int nm = (Integer) val;
                        double[] st = new double[nm];
                        for (int i = 1; i <= nm; i++)
                            st[i - 1] = factorial(i) * ini.mult(Matrix.pow(Km.scale(-1.0).pinv(), i + 1)).mult(clok).elementSum();
                        Ret.add(st);
                    } else if (opt.equals("stDistr")) {
                        double[] pts = (double[]) val;
                        double[] res = new double[pts.length];
                        for (int xi = 0; xi < pts.length; xi++)
                            res[xi] = mass0.mult(Rk).scale(1.0 / lamk).elementSum()
                                    + ini.mult(Km.scale(-1.0).pinv()).mult(Matrix.eye(KN).sub(Km.scale(pts[xi]).expm())).mult(clok).elementSum();
                        Ret.add(res);
                    } else if (opt.equals("flMoms")) {
                        int nm = (Integer) val;
                        jline.lib.butools.queues.FluFluResult r = jline.lib.butools.queues.FluFluQueue.fluFluQueue(
                                Q, Rk, Matrix.zeros(1, 1), new Matrix(new double[][]{{d}}), false, nm, 0, prec);
                        Ret.add(r.getFluidMoments());
                    } else if (opt.equals("flDistr")) {
                        double[] pts = (double[]) val;
                        jline.lib.butools.queues.FluFluResult r = jline.lib.butools.queues.FluFluQueue.fluFluQueue(
                                Q, Rk, Matrix.zeros(1, 1), new Matrix(new double[][]{{d}}), false, 1, 0, prec);
                        GeneralFluidSolution fs = r.getFluidSolution();
                        Matrix y = FluidStationaryDistr.fluidStationaryDistr(fs.getMass0(), fs.getIni(), fs.getK(), fs.getClo(), pts);
                        double[] res = new double[pts.length];
                        for (int xi = 0; xi < pts.length; xi++) {
                            double s = 0;
                            for (int j = 0; j < y.getNumCols(); j++) s += y.get(xi, j);
                            res[xi] = s;
                        }
                        Ret.add(res);
                    } else {
                        throw new IllegalArgumentException("FluidPrioQueue: Unknown parameter " + opt);
                    }
                }
            }
        }
        return Ret;
    }

    // ---------------- internal algorithms ----------------

    private static Matrix dReward(Matrix Qm, Matrix Rm, int n, double prec) {
        int NQ = Qm.getNumRows();
        if (NQ == 0) return Matrix.zeros(0, 0);
        double[] dR = diagArray(Rm);
        int[] ixz = idxAbsLe(dR, prec);
        int[] ixp = concat(idxGt(dR, prec), idxLt(dR, -prec));
        int Nz = ixz.length, Np = ixp.length;
        Matrix Per = Matrix.zeros(NQ, NQ);
        for (int w = 0; w < Nz; w++) Per.set(w, ixz[w], 1.0);
        for (int w = 0; w < Np; w++) Per.set(Nz + w, ixp[w], 1.0);
        Matrix iPer = Per.transpose(); // permutation inverse = transpose
        Matrix Rp = sub(Rm, ixp, ixp), Qpp = sub(Qm, ixp, ixp), Qpz = sub(Qm, ixp, ixz),
                Qzp = sub(Qm, ixz, ixp), Qzz = sub(Qm, ixz, ixz);
        Matrix iRp = Rp.pinv();
        Matrix inner = Nz > 0
                ? iRp.mult(Qpp.scale(-1.0).sub(Qpz.mult(Qzz.scale(-1.0).pinv()).mult(Qzp)))
                : iRp.mult(Qpp.scale(-1.0));
        Matrix dXvn = Matrix.pow(inner.pinv(), n + 1).mult(iRp).scale(Math.pow(-1, n) * factorial(n));
        Matrix drpar;
        if (Nz > 0) {
            Matrix iQzz = Qzz.scale(-1.0).pinv();
            drpar = block(iQzz.mult(Qzp).mult(dXvn).mult(Qpz).mult(iQzz), iQzz.mult(Qzp).mult(dXvn),
                    dXvn.mult(Qpz).mult(iQzz), dXvn);
        } else {
            drpar = dXvn;
        }
        return iPer.mult(drpar).mult(Per);
    }

    private static Matrix perm(int NF, int[] ixz, int[] ixp, int[] ixn) {
        Matrix Per = Matrix.zeros(NF, NF);
        for (int i = 0; i < ixz.length; i++) Per.set(i, ixz[i], 1.0);
        for (int i = 0; i < ixp.length; i++) Per.set(ixz.length + i, ixp[i], 1.0);
        for (int i = 0; i < ixn.length; i++) Per.set(ixz.length + ixp.length + i, ixn[i], 1.0);
        return Per;
    }

    private static List<Matrix> busyPeriodRewardMoms(Matrix F, Matrix C, Matrix D, int numOfMoms, double prec) {
        int NF = F.getNumRows();
        double[] dC = diagArray(C);
        int[] ixz = idxAbsLe(dC, prec), ixp = idxGt(dC, prec), ixn = idxLt(dC, -prec);
        int Nz = ixz.length, Np = ixp.length, Nn = ixn.length;
        Matrix Per = perm(NF, ixz, ixp, ixn);
        Matrix iPer = Per.transpose();
        Matrix PF = Per.mult(F).mult(iPer);
        int[] rz = range(Nz), rp = shiftRange(Nz, Np), rn = shiftRange(Nz + Np, Nn);
        Matrix Fzz = sub(PF, rz, rz), Fpz = sub(PF, rp, rz), Fmz = sub(PF, rn, rz);
        Matrix Fzp = sub(PF, rz, rp), Fpp = sub(PF, rp, rp), Fmp = sub(PF, rn, rp);
        Matrix Fzm = sub(PF, rz, rn), Fpm = sub(PF, rp, rn), Fmm = sub(PF, rn, rn);
        Matrix Cm = sub(C, ixn, ixn), Cp = sub(C, ixp, ixp);
        Matrix Dm = sub(D, ixn, ixn), Dp = sub(D, ixp, ixp), Dz = sub(D, ixz, ixz);
        Matrix iFzz = Nz > 0 ? Fzz.scale(-1.0).pinv() : Matrix.zeros(0, 0);
        Matrix iCp = Np > 0 ? Cp.pinv() : Matrix.zeros(0, 0);
        Matrix iCm = Nn > 0 ? Cm.scale(-1.0).pinv() : Matrix.zeros(0, 0);

        Matrix[] Fppd = new Matrix[numOfMoms + 1], Fpmd = new Matrix[numOfMoms + 1];
        Matrix[] Fmpd = new Matrix[numOfMoms + 1], Fmmd = new Matrix[numOfMoms + 1];
        Fppd[0] = iCp.mult(Fpp.add(Fpz.mult(iFzz).mult(Fzp)));
        Fpmd[0] = iCp.mult(Fpm.add(Fpz.mult(iFzz).mult(Fzm)));
        Fmpd[0] = iCm.mult(Fmp.add(Fmz.mult(iFzz).mult(Fzp)));
        Fmmd[0] = iCm.mult(Fmm.add(Fmz.mult(iFzz).mult(Fzm)));
        for (int i = 1; i <= numOfMoms; i++) {
            Matrix dr = dReward(Fzz, Dz, i, prec);
            Fppd[i] = iCp.mult(Fpz).mult(dr).mult(Fzp);
            Fpmd[i] = iCp.mult(Fpz).mult(dr).mult(Fzm);
            Fmpd[i] = iCm.mult(Fmz).mult(dr).mult(Fzp);
            Fmmd[i] = iCm.mult(Fmz).mult(dr).mult(Fzm);
            if (i == 1) {
                Fppd[i] = Fppd[i].sub(iCp.mult(Dp));
                Fmmd[i] = Fmmd[i].sub(iCm.mult(Dm));
            }
        }
        Matrix Psi = FluidFundamentalMatrices.FluidFundamentalMatrices(Fppd[0], Fpmd[0], Fmpd[0], Fmmd[0], prec, null, null).get("P");
        Matrix[] BPM = new Matrix[numOfMoms + 1];
        BPM[0] = Psi;
        for (int i = 1; i <= numOfMoms; i++) {
            Matrix X = Psi.scale(-1.0).mult(Fmpd[i]).mult(Psi).add(Fpmd[i]);
            for (int m = 0; m < i; m++) {
                double c = comb(i, m);
                X = X.add(Fppd[i - m].add(Psi.mult(Fmpd[i - m])).mult(BPM[m])
                        .add(BPM[m].mult(Fmmd[i - m].add(Fmpd[i - m].mult(Psi)))).scale(c));
            }
            for (int l = 1; l < i; l++)
                for (int m = 1; m <= i - l; m++) {
                    double c = comb(i, l) * comb(i - l, m);
                    X = X.add(BPM[l].mult(Fmpd[i - l - m]).mult(BPM[m]).scale(c));
                }
            BPM[i] = Matrix.lyap(Fppd[0].add(Psi.mult(Fmpd[0])), Fmmd[0].add(Fmpd[0].mult(Psi)), X, null);
        }
        List<Matrix> out = new ArrayList<Matrix>();
        for (int i = 0; i < BPM.length; i++) {
            Matrix emb = Matrix.zeros(NF, NF);
            setSub(emb, rp, rn, BPM[i]);
            out.add(iPer.mult(emb).mult(Per));
        }
        return out;
    }

    // returns {pr}
    private static Matrix[] busyPeriodRewardDistr(Matrix F, Matrix C, Matrix D, double t, double prec, int erlMaxOrder) {
        return new Matrix[]{busyPeriodRewardDistrFull(F, C, D, t, prec, erlMaxOrder)[0]};
    }

    // returns [pr, Pn_0, Pn_1, ...] (full reordered)
    private static Matrix[] busyPeriodRewardDistrFull(Matrix F, Matrix C, Matrix D, double t, double prec, int erlMaxOrder) {
        int NF = F.getNumRows();
        double[] dC = diagArray(C);
        int[] ixz = idxAbsLe(dC, prec), ixp = idxGt(dC, prec), ixn = idxLt(dC, -prec);
        int Nz = ixz.length, Np = ixp.length, Nn = ixn.length;
        Matrix Per = perm(NF, ixz, ixp, ixn);
        Matrix iPer = Per.transpose();
        Matrix PF = Per.mult(F).mult(iPer);
        int[] rz = range(Nz), rp = shiftRange(Nz, Np), rn = shiftRange(Nz + Np, Nn);
        Matrix Fzz = sub(PF, rz, rz), Fpz = sub(PF, rp, rz), Fmz = sub(PF, rn, rz);
        Matrix Fzp = sub(PF, rz, rp), Fpp = sub(PF, rp, rp), Fmp = sub(PF, rn, rp);
        Matrix Fzm = sub(PF, rz, rn), Fpm = sub(PF, rp, rn), Fmm = sub(PF, rn, rn);
        Matrix Cm = sub(C, ixn, ixn), Cp = sub(C, ixp, ixp);
        Matrix Dm = sub(D, ixn, ixn), Dp = sub(D, ixp, ixp), Dz = sub(D, ixz, ixz);
        Matrix iCp = Np > 0 ? Cp.pinv() : Matrix.zeros(0, 0);
        Matrix iCm = Nn > 0 ? Cm.scale(-1.0).pinv() : Matrix.zeros(0, 0);

        int L = erlMaxOrder;
        double nu = L / t;
        Matrix Z = Nz > 0 ? Dz.scale(nu).sub(Fzz).pinv() : Matrix.zeros(0, 0);
        Matrix AFpp = iCp.mult(Fpp.sub(Dp.scale(nu)).add(Fpz.mult(Z).mult(Fzp)));
        Matrix AFpm = iCp.mult(Fpm.add(Fpz.mult(Z).mult(Fzm)));
        Matrix AFmp = iCm.mult(Fmp.add(Fmz.mult(Z).mult(Fzp)));
        Matrix AFmm = iCm.mult(Fmm.sub(Dm.scale(nu)).add(Fmz.mult(Z).mult(Fzm)));
        Matrix Psie = FluidFundamentalMatrices.FluidFundamentalMatrices(AFpp, AFpm, AFmp, AFmm, prec, null, null).get("P");
        List<Matrix> Pn = new ArrayList<Matrix>();
        Pn.add(Psie);
        Matrix pr = Psie.copy();
        Matrix AM = AFpp.add(Psie.mult(AFmp));
        Matrix BM = AFmm.add(AFmp.mult(Psie));
        Matrix nuDzZ = Nz > 0 ? Dz.scale(nu).mult(Z) : Matrix.zeros(0, 0);
        for (int n = 1; n < L; n++) {
            Matrix CM = iCp.mult(Dp.scale(nu)).mult(Pn.get(n - 1)).add(Pn.get(n - 1).mult(iCm.mult(Dm.scale(nu))));
            for (int i = 1; i < n; i++)
                CM = CM.add(Pn.get(i).mult(iCm.mult(Fmp)).mult(Pn.get(n - i)));
            if (Nz > 0) {
                Matrix nuDzZn = Matrix.pow(nuDzZ, n);
                CM = CM.add(iCp.mult(Fpz).mult(Z).mult(nuDzZn).mult(Fzm))
                        .sub(Psie.mult(iCm).mult(Fmz).mult(Z).mult(nuDzZn).mult(Fzp).mult(Psie));
                for (int i = 0; i < n; i++) {
                    Matrix nuDzZni = Matrix.pow(nuDzZ, n - i);
                    CM = CM.add(Pn.get(i).mult(iCm).mult(Fmz).mult(Z).mult(nuDzZni).mult(Fzm.add(Fzp.mult(Psie))));
                    CM = CM.add(iCp.mult(Fpz).add(Psie.mult(iCm).mult(Fmz)).mult(Z).mult(nuDzZni).mult(Fzp).mult(Pn.get(i)));
                }
                for (int i = 1; i < n; i++)
                    for (int j = 1; j <= n - i; j++) {
                        Matrix nuDzZnij = Matrix.pow(nuDzZ, n - i - j);
                        CM = CM.add(Pn.get(i).mult(iCm).mult(Fmz).mult(Z).mult(nuDzZnij).mult(Fzp).mult(Pn.get(j)));
                    }
            }
            Matrix PM = Matrix.lyap(AM, BM, CM, null);
            Pn.add(PM);
            pr = pr.add(PM);
        }
        // re-order
        Matrix[] out = new Matrix[1 + Pn.size()];
        out[0] = reorder(pr, iPer, Per, NF, rp, rn);
        for (int i = 0; i < Pn.size(); i++) out[1 + i] = reorder(Pn.get(i), iPer, Per, NF, rp, rn);
        return out;
    }

    private static Matrix reorder(Matrix M, Matrix iPer, Matrix Per, int NF, int[] rp, int[] rn) {
        Matrix emb = Matrix.zeros(NF, NF);
        setSub(emb, rp, rn, M);
        return iPer.mult(emb).mult(Per);
    }

    private static List<Matrix> splitList(Matrix[] full) {
        List<Matrix> l = new ArrayList<Matrix>();
        for (int i = 1; i < full.length; i++) l.add(full[i]); // Pn entries
        return l;
    }

    // ---------- small helpers ----------
    private static double factorial(int n) { double r = 1; for (int i = 2; i <= n; i++) r *= i; return r; }
    private static double comb(int n, int k) { if (k < 0 || k > n) return 0; double r = 1; for (int i = 0; i < k; i++) r = r * (n - i) / (i + 1); return r; }
    private static Matrix rowVec(Matrix M, int r) {
        Matrix v = Matrix.zeros(1, M.getNumCols());
        for (int j = 0; j < M.getNumCols(); j++) v.set(0, j, M.get(r, j));
        return v;
    }
    private static Matrix sumRowsRange(Matrix M, int r0, int r1) {
        Matrix v = Matrix.zeros(1, M.getNumCols());
        for (int r = r0; r < r1; r++) for (int j = 0; j < M.getNumCols(); j++) v.set(0, j, v.get(0, j) + M.get(r, j));
        return v;
    }
    private static int[] idxAbsLe(double[] a, double prec) {
        List<Integer> l = new ArrayList<Integer>();
        for (int i = 0; i < a.length; i++) if (Math.abs(a[i]) <= prec) l.add(i);
        return toArr(l);
    }
    private static int[] idxGt(double[] a, double prec) {
        List<Integer> l = new ArrayList<Integer>();
        for (int i = 0; i < a.length; i++) if (a[i] > prec) l.add(i);
        return toArr(l);
    }
    private static int[] idxLt(double[] a, double prec) {
        List<Integer> l = new ArrayList<Integer>();
        for (int i = 0; i < a.length; i++) if (a[i] < prec) l.add(i);
        return toArr(l);
    }
    private static int[] toArr(List<Integer> l) { int[] r = new int[l.size()]; for (int i = 0; i < r.length; i++) r[i] = l.get(i); return r; }
}
