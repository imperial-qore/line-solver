/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 *
 * Transient analysis of a single-class open queue via the Laplace-domain
 * transient QBD method plus numerical inverse Laplace transform. Supports
 * single-server MAP/MAP/1 (infinite buffer) and MAP/MAP/1/N (finite buffer);
 * arrival and service are read uniformly as (D0, D1) MAPs. This subsumes M/M/1
 * and M/PH/1, and covers correlated-MAP arrival/service that the libQBD/expm
 * fast path cannot represent. Port of the MATLAB solver_mam_transient_qbd.
 */
package jline.solvers.mam.handlers;

import java.util.Map;
import java.util.function.Function;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.FieldMatrix;

import jline.api.mam.Map_pie;
import jline.api.mam.Map_prob;
import jline.api.mam.TransientQbd;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.lib.lti.iltcme;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Solver_mam_transient_qbd {

    private Solver_mam_transient_qbd() {}

    // ---------------------------------------------------------------------
    // Applicability: single-server open queue with non-Poisson arrival or
    // correlated-MAP service selects the Laplace transient QBD solver.
    // ---------------------------------------------------------------------
    public static boolean applicable(NetworkStruct sn) {
        if (sn.nclasses != 1) {
            return false;
        }
        if (!Double.isInfinite(sn.njobs.get(0))) {
            return false;
        }
        int[] sq = findSourceQueue(sn);
        if (sq[0] < 0 || sq[1] < 0) {
            return false;
        }
        if ((int) sn.nservers.get(sq[1], 0) != 1) {
            return false;
        }
        Matrix[] arr = readMap(sn, sq[0]);
        Matrix[] svc = readMap(sn, sq[1]);
        if (arr == null || svc == null) {
            return false;
        }
        boolean arrivalIsPoisson = (arr[0].getNumRows() == 1);
        boolean serviceIsRenewal = isRenewalMap(svc[0], svc[1]);
        return !arrivalIsPoisson || !serviceIsRenewal;
    }

    private static boolean isRenewalMap(Matrix D0, Matrix D1) {
        int ns = D0.getNumRows();
        if (ns == 1) {
            return true;
        }
        Matrix tExit = D1.mult(ones(ns, 1));
        Matrix alpha = Map_pie.map_pie(D0, D1);       // 1 x ns
        double diff = 0.0, ref = 0.0;
        for (int i = 0; i < ns; i++) {
            for (int j = 0; j < ns; j++) {
                double d = D1.get(i, j) - tExit.get(i, 0) * alpha.get(0, j);
                diff += d * d;
                ref += D1.get(i, j) * D1.get(i, j);
            }
        }
        return Math.sqrt(diff) < 1e-9 * Math.max(1.0, Math.sqrt(ref));
    }

    // ---------------------------------------------------------------------
    // Main solver
    // ---------------------------------------------------------------------
    public static TransientResult solver_mam_transient_qbd(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;
        int K = sn.nclasses;
        int[] sq = findSourceQueue(sn);
        int sourceIdx = sq[0], queueIdx = sq[1];

        Matrix[] arr = readMap(sn, sourceIdx);
        Matrix[] svc = readMap(sn, queueIdx);
        Matrix Da0 = arr[0], Da1 = arr[1];
        Matrix Ds0 = svc[0], Ds1 = svc[1];
        int na = Da0.getNumRows();
        int ns = Ds0.getNumRows();
        Matrix Ina = eye(na);
        Matrix Ins = eye(ns);

        Matrix Lrep = Da0.krons(Ds0);
        Matrix Frep = Da1.kron(Ins);
        Matrix Brep = Ina.kron(Ds1);
        Matrix Lv0 = Da0.kron(Ins);
        Matrix F0 = Da1.kron(Ins);
        Matrix B0 = Ina.kron(Ds1);

        Matrix piArr = Map_prob.map_prob(Da0, Da1);   // 1 x na
        Matrix piSvc = Map_prob.map_prob(Ds0, Ds1);   // 1 x ns
        Matrix pi0 = piArr.kron(piSvc);               // 1 x na*ns

        Matrix sExit = Ds1.mult(ones(ns, 1));
        Matrix wDep = ones(na, 1).kron(sExit);        // na*ns x 1
        final FieldMatrix<Complex> wDepC = TransientQbd.cm(wDep);
        final FieldMatrix<Complex> wOneC = TransientQbd.cm(ones(na * ns, 1));

        double bufCap = sn.cap.get(queueIdx, 0);
        boolean isFinite = !Double.isInfinite(bufCap);

        double tStart = options.timespan[0];
        double tEnd = options.timespan[1];
        double dur = tEnd - tStart;
        int nTime = (int) Math.min(101, Math.max(11, Math.round(dur * 10)));
        double[] times = linspace(tStart, tEnd, nTime);
        int maxFnEvals = 100;

        // Positive-time mask (ILT is singular at t=0; system starts empty).
        int nPos = 0;
        for (int i = 0; i < nTime; i++) {
            if (times[i] > 0) {
                nPos++;
            }
        }
        double[] tpos = new double[nPos];
        int[] posIdx = new int[nPos];
        int c = 0;
        for (int i = 0; i < nTime; i++) {
            if (times[i] > 0) {
                tpos[c] = times[i];
                posIdx[c] = i;
                c++;
            }
        }

        double[] EN = new double[nTime];
        double[] DEP = new double[nTime];
        double[] Uval = new double[nTime];

        // Single matrix-valued inverse Laplace pass returns [E[N]; Tput; P0] per
        // time point, so the transient transform is evaluated once per node
        // (rather than three times for three separate scalar inversions).
        final Function<Complex, Complex[][]> metricsFun;
        if (isFinite) {
            final int Ncap = (int) bufCap;
            Matrix LvTop = Da0.add(1.0, Da1).kron(Ins).add(1.0, Ina.kron(Ds0));
            final FieldMatrix<Complex>[] B = fmArr(null, TransientQbd.cm(Brep));
            final FieldMatrix<Complex>[] L = fmArr(null, TransientQbd.cm(Lrep));
            final FieldMatrix<Complex>[] F = fmArr(null, TransientQbd.cm(Frep));
            final FieldMatrix<Complex>[] Lv = fmArr(null, TransientQbd.cm(Lv0), TransientQbd.cm(LvTop));
            final int[] T = {0, 0, Ncap};
            metricsFun = s -> {
                Complex en = Complex.ZERO;
                Complex dep = Complex.ZERO;
                for (int mm = 1; mm <= Ncap; mm++) {
                    FieldMatrix<Complex> V = TransientQbd.transient2(B, L, F, Lv, T, 0, mm, s);
                    en = en.add(piMV(pi0, V, wOneC).multiply(mm));
                    dep = dep.add(piMV(pi0, V, wDepC));
                }
                Complex p0 = piMV(pi0, TransientQbd.transient2(B, L, F, Lv, T, 0, 0, s), wOneC);
                return new Complex[][]{{en}, {dep}, {p0}};
            };
        } else {
            final FieldMatrix<Complex>[] B = fmArr(null, TransientQbd.cm(B0), TransientQbd.cm(Brep));
            final FieldMatrix<Complex>[] L = fmArr(null, null, TransientQbd.cm(Lrep));
            final FieldMatrix<Complex>[] F = fmArr(null, TransientQbd.cm(F0), TransientQbd.cm(Frep));
            final FieldMatrix<Complex>[] Lv = fmArr(null, TransientQbd.cm(Lv0), TransientQbd.cm(Lrep));
            final int[] T = {0, 0, 1};
            final int Kreg = 2;
            final FieldMatrix<Complex> wOne = wOneC;
            final FieldMatrix<Complex> wDepF = wDepC;
            // see _kb/06-solver-catalog.md for rationale
            metricsFun = s -> {
                FieldMatrix<Complex> Lk = L[Kreg].subtract(TransientQbd.sEye(L[Kreg].getRowDimension(), s));
                FieldMatrix<Complex> R = TransientQbd.qbdFundMat(B[Kreg], Lk, F[Kreg])[1];
                int nk = R.getRowDimension();
                FieldMatrix<Complex> ImR = TransientQbd.eye(nk).subtract(R);
                FieldMatrix<Complex> V1 = TransientQbd.transient2Open(B, L, F, Lv, T, 0, 1, s);
                FieldMatrix<Complex> V0 = TransientQbd.transient2Open(B, L, F, Lv, T, 0, 0, s);
                FieldMatrix<Complex> onesCol = TransientQbd.cm(ones(nk, 1));
                Complex en = piMV(pi0, V1, TransientQbd.ldiv(ImR, TransientQbd.ldiv(ImR, onesCol)));
                Complex dep = piMV(pi0, V1, TransientQbd.ldiv(ImR, wDepF));
                Complex p0 = piMV(pi0, V0, wOne);
                return new Complex[][]{{en}, {dep}, {p0}};
            };
        }

        double[][][] out = iltcme.iltMatrix(metricsFun, tpos, maxFnEvals);
        for (int i = 0; i < nPos; i++) {
            EN[posIdx[i]] = out[i][0][0];
            DEP[posIdx[i]] = out[i][1][0];
            Uval[posIdx[i]] = 1.0 - out[i][2][0];
        }

        Matrix[][] Qt = new Matrix[M][K];
        Matrix[][] Ut = new Matrix[M][K];
        Matrix[][] Tt = new Matrix[M][K];
        Matrix qResult = new Matrix(nTime, 2);
        Matrix uResult = new Matrix(nTime, 2);
        Matrix tResult = new Matrix(nTime, 2);
        for (int t = 0; t < nTime; t++) {
            qResult.set(t, 0, EN[t]);   qResult.set(t, 1, times[t]);
            uResult.set(t, 0, Uval[t]); uResult.set(t, 1, times[t]);
            tResult.set(t, 0, DEP[t]);  tResult.set(t, 1, times[t]);
        }
        Qt[queueIdx][0] = qResult;
        Ut[queueIdx][0] = uResult;
        Tt[queueIdx][0] = tResult;
        return new TransientResult(Qt, Ut, Tt);
    }

    // ---------------------------------------------------------------------
    // Helpers
    // ---------------------------------------------------------------------
    private static Complex piMV(Matrix pi, FieldMatrix<Complex> V, FieldMatrix<Complex> w) {
        FieldMatrix<Complex> Vw = V.multiply(w);   // N x 1
        Complex acc = Complex.ZERO;
        for (int i = 0; i < Vw.getRowDimension(); i++) {
            acc = acc.add(Vw.getEntry(i, 0).multiply(pi.get(0, i)));
        }
        return acc;
    }

    @SafeVarargs
    private static FieldMatrix<Complex>[] fmArr(FieldMatrix<Complex>... a) {
        return a;
    }

    private static int[] findSourceQueue(NetworkStruct sn) {
        int sourceIdx = -1, queueIdx = -1;
        for (int i = 0; i < sn.nstations; i++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(i));
            if (sched == SchedStrategy.EXT) {
                sourceIdx = i;
            } else if (sched == SchedStrategy.FCFS) {
                queueIdx = i;
            }
        }
        return new int[]{sourceIdx, queueIdx};
    }

    private static Matrix[] readMap(NetworkStruct sn, int idx) {
        Station station = sn.stations.get(idx);
        JobClass jobClass = sn.jobclasses.get(0);
        Map<JobClass, MatrixCell> procMap = sn.proc.get(station);
        if (procMap == null) {
            return null;
        }
        MatrixCell cell = procMap.get(jobClass);
        if (cell == null) {
            return null;
        }
        Matrix D0 = cell.get(0);
        Matrix D1 = cell.get(1);
        if (D0 == null || D1 == null) {
            return null;
        }
        return new Matrix[]{D0, D1};
    }

    private static Matrix eye(int n) {
        Matrix I = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            I.set(i, i, 1.0);
        }
        return I;
    }

    private static Matrix ones(int r, int col) {
        Matrix o = new Matrix(r, col);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < col; j++) {
                o.set(i, j, 1.0);
            }
        }
        return o;
    }

    private static double[] linspace(double a, double b, int n) {
        double[] out = new double[n];
        if (n == 1) {
            out[0] = a;
            return out;
        }
        double step = (b - a) / (n - 1);
        for (int i = 0; i < n; i++) {
            out[i] = a + i * step;
        }
        return out;
    }

}
