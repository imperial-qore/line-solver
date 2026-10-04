/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.special.Gamma;

import jline.api.mam.Map_lambda;
import jline.api.mc.Dtmc_solve;
import jline.util.matrix.Matrix;

import static jline.io.InputOutput.line_error;

/**
 * The MAP/G/1/K queue with tail drop: Markovian arrivals, an arbitrary service
 * law F, and a buffer of K packets counting the one in transmission.
 *
 * <p>Port of {@code matlab/src/api/qsys/qsys_mapg1k.m}, the twin of
 * {@code cpp/include/line/api/qsys/qsys_mapg1k.h} and of the native Python
 * {@code api.qsys.mapg1k}. Unlike {@link Qsys_mapg1} the service law is NOT
 * fitted to a phase-type distribution: F enters exactly through the functionals
 * A_m and Q_m, evaluated by uniformizing the arrival MAP at
 * theta = max_i (-D0(i,i)). Unlike {@link Qsys_mg1k_loss}, which embeds the same
 * way but assumes Poisson input, arrivals may be a general MAP, so flows of
 * equal rate but different interarrival variability or autocorrelation are told
 * apart.
 *
 * <p>METHOD. The chain embedded at departure epochs has state (n, j) with
 * n = 0..K-1 the packets left behind by a departure and j the MAP phase. With
 * A_m the matrix of "m arrivals during one service, phase i to phase j",
 * <pre>
 *   n &gt;= 1: n' = n - 1 + min(m, K - n), the overflow being sum_{m &gt;= K-n} A_m
 *   n == 0: the phase first jumps through Psi = (-D0)^-1 D1, the idle period
 *           ending at an arrival, and the service proceeds as from n = 1.
 * </pre>
 * Its stationary law sigma gives, by Markov renewal reward over one
 * inter-departure cycle, E[cycle] = S + sigma_0 (-D0)^-1 e, T = 1/E[cycle] and
 * p0 = 1 - T S; the level-holding times come from Q_m, the expected time within
 * a service during which exactly m arrivals have occurred. No PASTA argument is
 * used anywhere -- the MAP phase resolution does that work instead, which is
 * what lets {@link Qsys_mmapg1k} read exact PER-CLASS loss ratios off pKvec.
 *
 * <p>This is not the transform solution of Theorem 1 of [1], stated through a
 * sequence R_m obeying R(z) = z (A(z) - z I)^-1. That sequence grows
 * geometrically while the quantity extracted from it stays O(s) as s -&gt; 0+, so
 * the conditioning degrades like that ratio^K and crosses the double-precision
 * ceiling near K = 20 under gamma service with CV = 2. Reference [1] evaluates
 * its formulae in arbitrary precision, so the restriction is invisible there.
 * The embedded chain used here has every entry a probability or a time.
 *
 * <p>The uniformization coefficients c_n are built by RECURSION rather than from
 * log-gamma, as in the C++ twin: for the gamma law they are the negative
 * binomial pmf, for the deterministic law the Poisson pmf, and for a PH law
 * c_n = theta^n alpha M^{n+1} t with M = (theta I - T)^-1 accumulated as a
 * running row vector. The recursions are monotone in n and cannot lose leading
 * digits to cancellation.
 *
 * <p>References:
 * <br>[1] Chydzinski, A. Per-Flow Throughput of a FIFO Buffer. Applied System
 * Innovation 2026, 9, 112.
 * <br>[2] Niu, Z.; Cooper, R.B. Transform-Free Analysis of M/G/1/K and Related
 * Queues. Mathematics of Operations Research 1993, 18, 486-510.
 */
public final class Qsys_mapg1k {

    private Qsys_mapg1k() {
    }

    /** Default uniformization truncation tolerance. */
    public static final double DEFAULT_TOL = 1e-12;
    /** Default cap on the uniformization order. */
    public static final int DEFAULT_NMAX = 200000;

    /** Coefficients c_n of the uniformized service law, with the mean and residual. */
    static final class ServiceCoefficients {
        double[] cn;
        double mean;
        double residual;
    }

    /**
     * E[w(S)] for a density-specified service law, under the substitution
     * x = exp(u): the integrand becomes w(e^u) f(e^u) e^u, which tends to zero at
     * both ends because w is bounded and integrability of f forces x f(x) -&gt; 0.
     * The window widens geometrically until a whole new panel adds less than the
     * tolerance, so the truncation is measured rather than assumed.
     */
    private static double densityMoment(final QsysServiceLaw svc, final UnivariateFunction w,
                                        double reltol) {
        final UnivariateFunction g = new UnivariateFunction() {
            @Override
            public double value(double u) {
                double x = Math.exp(u);
                double v = w.value(x) * svc.pdf.value(x) * x;
                // exp(u) underflows while pdf(exp(u)) overflows: the NaN is an
                // artifact of the substitution, not of the integrand.
                return (Double.isNaN(v) || Double.isInfinite(v)) ? 0.0 : v;
            }
        };
        IterativeLegendreGaussIntegrator integ =
                new IterativeLegendreGaussIntegrator(15, reltol, 1e-300, 2, 64);
        double hi0 = svc.tmaxFinite ? Math.log(svc.tmax) : 1.0;
        double lo = hi0 - 2.0;
        double hi = hi0;
        double total = panel(integ, g, lo, hi);
        for (int k = 0; k < 60; k++) {
            double lonew = lo - 4.0;
            double addLo = panel(integ, g, lonew, lo);
            lo = lonew;
            double addHi = 0.0;
            if (!svc.tmaxFinite) {
                double hinew = hi + 4.0;
                addHi = panel(integ, g, hi, hinew);
                hi = hinew;
            }
            total += addLo + addHi;
            if (Math.abs(addLo) + Math.abs(addHi) <= reltol * Math.abs(total)) {
                break;
            }
        }
        return total;
    }

    private static double panel(IterativeLegendreGaussIntegrator integ, UnivariateFunction g,
                                double a, double b) {
        double v = integ.integrate(1 << 20, g, a, b);
        return Double.isNaN(v) ? 0.0 : v;
    }

    /** Mean of the service law; the density case is handled by the caller. */
    private static double serviceMean(QsysServiceLaw svc) {
        switch (svc.kind) {
            case GAMMA:
                return svc.shape * svc.scale;
            case DETERMINISTIC:
                return svc.det;
            case PHASE_TYPE: {
                int p = svc.phT.getNumRows();
                Matrix e = new Matrix(p, 1);
                for (int i = 0; i < p; i++) {
                    e.set(i, 0, 1.0);
                }
                Matrix x = svc.phT.inv().mult(e);   // T^-1 e
                double s = 0.0;
                for (int i = 0; i < p; i++) {
                    s -= svc.phAlpha.get(i) * x.get(i, 0);
                }
                return s;
            }
            default:
                return Double.NaN;
        }
    }

    /**
     * c_n = E[e^{-theta S} (theta S)^n / n!] for n = 0..N, plus the mean.
     * sum_n c_n = E[e^{-theta S} e^{theta S}] = 1 exactly, which both sets the
     * truncation order and certifies it.
     */
    static ServiceCoefficients serviceCoefficients(QsysServiceLaw svc, double theta, double tol,
                                                   int nmaxCap) {
        final double qreltol = 1e-13;
        ServiceCoefficients out = new ServiceCoefficients();
        out.mean = (svc.kind == QsysServiceLaw.Kind.DENSITY)
                ? densityMoment(svc, new UnivariateFunction() {
                    @Override
                    public double value(double x) {
                        return x;
                    }
                }, qreltol)
                : serviceMean(svc);
        if (!(out.mean > 0.0)) {
            line_error("qsys_mapg1k", "the service law has a non-positive mean.");
        }

        double md = theta * out.mean;
        double guess = md + 10.0 * Math.sqrt(Math.max(md, 1.0)) + 32.0;
        int n0 = (int) Math.max(32.0, Math.ceil(guess));
        if (n0 + 1 > nmaxCap) {
            n0 = Math.max(nmaxCap - 1, 0);
        }

        double[] cn = new double[n0 + 1];
        double[][] phState = null;   // carried row alpha*M^(n+1) and the exit vector
        switch (svc.kind) {
            case GAMMA: {
                if (!(svc.shape > 0.0) || !(svc.scale > 0.0)) {
                    line_error("qsys_mapg1k", "the gamma shape and scale must be positive.");
                }
                double q = 1.0 + svc.scale * theta;
                double p = svc.scale * theta / q;
                cn[0] = Math.exp(-svc.shape * Math.log(q));
                for (int n = 1; n <= n0; n++) {
                    cn[n] = cn[n - 1] * p * (svc.shape + n - 1.0) / n;
                }
                break;
            }
            case DETERMINISTIC: {
                if (!(svc.det > 0.0)) {
                    line_error("qsys_mapg1k", "the deterministic service time must be positive.");
                }
                double m = theta * svc.det;
                cn[0] = Math.exp(-m);
                for (int n = 1; n <= n0; n++) {
                    cn[n] = cn[n - 1] * m / n;
                }
                break;
            }
            case PHASE_TYPE: {
                phState = phInit(svc, theta);
                for (int n = 0; n <= n0; n++) {
                    cn[n] = phNext(phState, theta);
                }
                break;
            }
            case DENSITY: {
                if (svc.pdf == null) {
                    line_error("qsys_mapg1k", "the density service law has no pdf.");
                }
                for (int n = 0; n <= n0; n++) {
                    cn[n] = densityMoment(svc, poissonWeight(theta, n), qreltol);
                }
                break;
            }
            default:
                line_error("qsys_mapg1k", "unsupported service law kind.");
        }

        double total = 0.0;
        for (int i = 0; i < cn.length; i++) {
            total += cn[i];
        }
        // Grow in blocks of 64 until the series sums to 1 within tol, or until a
        // whole block adds nothing in floating point (the representable series is
        // exhausted, which is the only stop available to a quadrature path).
        while (cn.length < nmaxCap) {
            if (Math.abs(1.0 - total) <= tol) {
                break;
            }
            int first = cn.length;
            int last = Math.min(first + 63, nmaxCap - 1);
            if (last < first) {
                break;
            }
            double[] grown = new double[last + 1];
            System.arraycopy(cn, 0, grown, 0, cn.length);
            double added = 0.0;
            switch (svc.kind) {
                case GAMMA: {
                    double q = 1.0 + svc.scale * theta;
                    double p = svc.scale * theta / q;
                    for (int n = first; n <= last; n++) {
                        grown[n] = grown[n - 1] * p * (svc.shape + n - 1.0) / n;
                        added += grown[n];
                    }
                    break;
                }
                case DETERMINISTIC: {
                    double m = theta * svc.det;
                    for (int n = first; n <= last; n++) {
                        grown[n] = grown[n - 1] * m / n;
                        added += grown[n];
                    }
                    break;
                }
                case PHASE_TYPE: {
                    for (int n = first; n <= last; n++) {
                        grown[n] = phNext(phState, theta);
                        added += grown[n];
                    }
                    break;
                }
                default: {
                    for (int n = first; n <= last; n++) {
                        grown[n] = densityMoment(svc, poissonWeight(theta, n), qreltol);
                        added += grown[n];
                    }
                    break;
                }
            }
            cn = grown;
            total += added;
            if (added <= Math.ulp(1.0) * total) {
                break;
            }
        }
        out.cn = cn;
        out.residual = Math.abs(1.0 - total);
        return out;
    }

    /** State of the PH recursion: {row alpha*M^(n+1), exit vector t}. */
    private static double[][] phInit(QsysServiceLaw svc, double theta) {
        int p = svc.phT.getNumRows();
        if (svc.phAlpha.length() != p || svc.phT.getNumCols() != p) {
            line_error("qsys_mapg1k", "the PH alpha and T are inconsistent.");
        }
        Matrix ThI = new Matrix(p, p);
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < p; j++) {
                ThI.set(i, j, (i == j ? theta : 0.0) - svc.phT.get(i, j));
            }
        }
        Matrix Minv = ThI.inv();
        double[] t = new double[p];
        for (int i = 0; i < p; i++) {
            double s = 0.0;
            for (int j = 0; j < p; j++) {
                s += svc.phT.get(i, j);
            }
            t[i] = -s;
        }
        double[] row = new double[p];
        for (int j = 0; j < p; j++) {
            double s = 0.0;
            for (int i = 0; i < p; i++) {
                s += svc.phAlpha.get(i) * Minv.get(i, j);
            }
            row[j] = s;
        }
        double[] flat = new double[p * p];
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < p; j++) {
                flat[i * p + j] = Minv.get(i, j);
            }
        }
        return new double[][]{row, t, flat};
    }

    /** c_n = row*t, then advance the row to theta*row*M for the next n. */
    private static double phNext(double[][] state, double theta) {
        double[] row = state[0];
        double[] t = state[1];
        double[] Minv = state[2];
        int p = row.length;
        double v = 0.0;
        for (int i = 0; i < p; i++) {
            v += row[i] * t[i];
        }
        double[] next = new double[p];
        for (int j = 0; j < p; j++) {
            double s = 0.0;
            for (int i = 0; i < p; i++) {
                s += row[i] * Minv[i * p + j];
            }
            next[j] = theta * s;
        }
        state[0] = next;
        return v;
    }

    /** x -&gt; e^{-theta x} (theta x)^n / n!, the Poisson weight of order n. */
    private static UnivariateFunction poissonWeight(final double theta, final int n) {
        final double lfact = Gamma.logGamma(n + 1.0);
        return new UnivariateFunction() {
            @Override
            public double value(double x) {
                return Math.exp(-theta * x + n * Math.log(theta * x) - lfact);
            }
        };
    }

    /** qsys_mapg1k with the reference defaults tol = 1e-12, nmax = 200000. */
    public static QsysMapG1kResult qsys_mapg1k(Matrix D0, Matrix D1, QsysServiceLaw svc, int K) {
        return qsys_mapg1k(D0, D1, svc, K, DEFAULT_TOL, DEFAULT_NMAX);
    }

    /**
     * MAP/G/1/K with tail drop.
     *
     * @param D0      M x M hidden transition matrix of the arrival MAP
     * @param D1      M x M arrival matrix; D0 + D1 must be an irreducible generator
     * @param svc     service law
     * @param K       buffer size in packets, K &gt;= 1, the one in service included
     * @param tol     uniformization truncation tolerance
     * @param nmaxCap cap on the uniformization order
     */
    public static QsysMapG1kResult qsys_mapg1k(Matrix D0, Matrix D1, QsysServiceLaw svc, int K,
                                               double tol, int nmaxCap) {
        int M = D0.getNumRows();
        if (D0.getNumCols() != M || D1.getNumRows() != M || D1.getNumCols() != M) {
            line_error("qsys_mapg1k", "D0 and D1 must be square matrices of equal size.");
        }
        if (K < 1) {
            line_error("qsys_mapg1k", "the buffer size K must be a positive integer.");
        }
        double theta = 0.0;
        for (int i = 0; i < M; i++) {
            double b = -D0.get(i, i);
            if (!(b > 0.0)) {
                line_error("qsys_mapg1k", "D0 must have strictly negative diagonal entries.");
            }
            theta = Math.max(theta, b);
        }

        ServiceCoefficients sc = serviceCoefficients(svc, theta, tol, nmaxCap);
        double[] cn = sc.cn;
        double Smean = sc.mean;
        int nmax = cn.length - 1;
        if (sc.residual > 1e-6) {
            line_error("qsys_mapg1k", "the uniformization series for c_n was truncated with "
                    + "residual " + sc.residual + "; raise nmax.");
        }

        // d_n = (sum_{k > n} c_k)/theta, the expected time within a service with
        // exactly n arrivals so far; accumulated as a telescoping tail.
        double[] dn = new double[nmax + 1];
        double tail = 0.0;
        for (int n = nmax; n >= 0; n--) {
            dn[n] = tail / theta;
            tail += cn[n];
        }

        int mmax = K - 1;
        Matrix[] A = new Matrix[mmax + 1];
        Matrix[] Q = new Matrix[mmax + 1];
        Matrix[] Sn = new Matrix[mmax + 1];
        for (int m = 0; m <= mmax; m++) {
            A[m] = new Matrix(M, M);
            Q[m] = new Matrix(M, M);
            Sn[m] = new Matrix(M, M);
        }
        for (int i = 0; i < M; i++) {
            Sn[0].set(i, i, 1.0);
        }
        Matrix B0 = new Matrix(M, M);
        Matrix Qtot = new Matrix(M, M);
        Matrix Pn = Matrix.eye(M);
        Matrix Pt0 = new Matrix(M, M);
        Matrix Pt1 = new Matrix(M, M);
        Matrix PD = new Matrix(M, M);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                Pt0.set(i, j, (i == j ? 1.0 : 0.0) + D0.get(i, j) / theta);
                Pt1.set(i, j, D1.get(i, j) / theta);
                PD.set(i, j, (i == j ? 1.0 : 0.0) + (D0.get(i, j) + D1.get(i, j)) / theta);
            }
        }
        for (int n = 0; n <= nmax; n++) {
            int mtop = Math.min(n, mmax);
            for (int m = 0; m <= mtop; m++) {
                for (int i = 0; i < M; i++) {
                    for (int j = 0; j < M; j++) {
                        A[m].set(i, j, A[m].get(i, j) + Sn[m].get(i, j) * cn[n]);
                        Q[m].set(i, j, Q[m].get(i, j) + Sn[m].get(i, j) * dn[n]);
                    }
                }
            }
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    B0.set(i, j, B0.get(i, j) + Pn.get(i, j) * cn[n]);
                    Qtot.set(i, j, Qtot.get(i, j) + Pn.get(i, j) * dn[n]);
                }
            }
            if (n < nmax) {
                Matrix[] Snew = new Matrix[mmax + 1];
                int mt = Math.min(n + 1, mmax);
                for (int m = 0; m <= mmax; m++) {
                    Snew[m] = new Matrix(M, M);
                }
                for (int m = 0; m <= mt; m++) {
                    Matrix acc = new Matrix(M, M);
                    if (m <= n) {
                        acc = Sn[m].mult(Pt0);
                    }
                    if (m >= 1) {
                        Matrix add = Sn[m - 1].mult(Pt1);
                        for (int i = 0; i < M; i++) {
                            for (int j = 0; j < M; j++) {
                                acc.set(i, j, acc.get(i, j) + add.get(i, j));
                            }
                        }
                    }
                    Snew[m] = acc;
                }
                Sn = Snew;
                Pn = Pn.mult(PD);
            }
        }

        Matrix negD0 = new Matrix(M, M);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                negD0.set(i, j, -D0.get(i, j));
            }
        }
        Matrix negD0inv = negD0.inv();
        Matrix Psi = negD0inv.mult(D1);        // phase at the arrival ending an idle period
        double[] idle = new double[M];
        for (int i = 0; i < M; i++) {
            double s = 0.0;
            for (int j = 0; j < M; j++) {
                s += negD0inv.get(i, j);
            }
            idle[i] = s;
        }

        // Embedded chain at departure epochs, state (n,j) -> index n*M + j.
        Matrix P = new Matrix(K * M, K * M);
        int lastblk = (K - 1) * M;
        for (int n = 1; n <= K - 1; n++) {
            Matrix Bacc = B0.copy();
            for (int m = 0; m <= K - n - 1; m++) {
                int col = (n - 1 + m) * M;
                for (int i = 0; i < M; i++) {
                    for (int j = 0; j < M; j++) {
                        P.set(n * M + i, col + j, P.get(n * M + i, col + j) + A[m].get(i, j));
                        Bacc.set(i, j, Bacc.get(i, j) - A[m].get(i, j));
                    }
                }
            }
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    P.set(n * M + i, lastblk + j, P.get(n * M + i, lastblk + j) + Bacc.get(i, j));
                }
            }
        }
        Matrix Bacc0 = B0.copy();
        for (int m = 0; m <= K - 2; m++) {
            Matrix blk = Psi.mult(A[m]);
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    P.set(i, m * M + j, P.get(i, m * M + j) + blk.get(i, j));
                    Bacc0.set(i, j, Bacc0.get(i, j) - A[m].get(i, j));
                }
            }
        }
        Matrix tailblk = Psi.mult(Bacc0);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                P.set(i, lastblk + j, P.get(i, lastblk + j) + tailblk.get(i, j));
            }
        }

        double rowdev = 0.0;
        for (int i = 0; i < K * M; i++) {
            double s = 0.0;
            for (int j = 0; j < K * M; j++) {
                s += P.get(i, j);
            }
            rowdev = Math.max(rowdev, Math.abs(s - 1.0));
        }
        if (rowdev > 1e-8) {
            line_error("qsys_mapg1k", "the embedded chain rows deviate from 1 by " + rowdev
                    + "; the uniformization series for A_m has not converged, raise nmax.");
        }

        Matrix sigmaM = Dtmc_solve.dtmc_solve(P);
        double[] sigma = new double[K * M];
        for (int i = 0; i < K * M; i++) {
            sigma[i] = sigmaM.get(i);
        }

        double idleTime = 0.0;
        for (int j = 0; j < M; j++) {
            idleTime += sigma[j] * idle[j];
        }
        double Ecyc = Smean + idleTime;
        double Tput = 1.0 / Ecyc;
        double p0 = idleTime / Ecyc;

        Matrix[] Qcum = new Matrix[mmax + 1];
        Matrix acc = new Matrix(M, M);
        for (int m = 0; m <= mmax; m++) {
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    acc.set(i, j, acc.get(i, j) + Q[m].get(i, j));
                }
            }
            Qcum[m] = acc.copy();
        }
        double[] timeKvec = new double[M];
        for (int n = 1; n <= K - 1; n++) {
            int rr = K - n - 1;
            for (int j = 0; j < M; j++) {
                double s = 0.0;
                for (int i = 0; i < M; i++) {
                    s += sigma[n * M + i] * (Qtot.get(i, j) - Qcum[rr].get(i, j));
                }
                timeKvec[j] += s;
            }
        }
        double[] s0Psi = new double[M];
        for (int j = 0; j < M; j++) {
            double s = 0.0;
            for (int i = 0; i < M; i++) {
                s += sigma[i] * Psi.get(i, j);
            }
            s0Psi[j] = s;
        }
        for (int j = 0; j < M; j++) {
            double s = 0.0;
            for (int i = 0; i < M; i++) {
                double q = Qtot.get(i, j) - (K >= 2 ? Qcum[K - 2].get(i, j) : 0.0);
                s += s0Psi[i] * q;
            }
            timeKvec[j] += s;
        }
        Matrix pKvec = new Matrix(1, M);
        double pK = 0.0;
        for (int j = 0; j < M; j++) {
            pKvec.set(0, j, timeKvec[j] / Ecyc);
            pK += pKvec.get(0, j);
        }

        double[] timeL = new double[K + 1];
        timeL[0] = idleTime;
        for (int n = 1; n <= K - 1; n++) {
            for (int l = n; l <= K - 1; l++) {
                double s = 0.0;
                for (int i = 0; i < M; i++) {
                    for (int j = 0; j < M; j++) {
                        s += sigma[n * M + i] * Q[l - n].get(i, j);
                    }
                }
                timeL[l] += s;
            }
        }
        for (int l = 1; l <= K - 1; l++) {
            double s = 0.0;
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    s += s0Psi[i] * Q[l - 1].get(i, j);
                }
            }
            timeL[l] += s;
        }
        double sK = 0.0;
        for (int j = 0; j < M; j++) {
            sK += timeKvec[j];
        }
        timeL[K] = sK;

        Matrix plevel = new Matrix(1, K + 1);
        double mass = 0.0;
        for (int l = 0; l <= K; l++) {
            plevel.set(0, l, timeL[l] / Ecyc);
            mass += plevel.get(0, l);
        }
        if (Math.abs(mass - 1.0) > 1e-8) {
            line_error("qsys_mapg1k", "the level distribution has mass " + mass
                    + "; the Q_m series has not converged, raise nmax.");
        }
        double meanQ = 0.0;
        for (int l = 0; l <= K; l++) {
            meanQ += l * plevel.get(0, l);
        }

        Matrix p0vec = new Matrix(1, M);
        for (int j = 0; j < M; j++) {
            double s = 0.0;
            for (int i = 0; i < M; i++) {
                s += sigma[i] * negD0inv.get(i, j);
            }
            p0vec.set(0, j, s / Ecyc);
        }

        double lambda = Map_lambda.map_lambda(D0, D1);

        QsysMapG1kResult r = new QsysMapG1kResult();
        r.p0 = p0;
        r.pK = pK;
        r.throughput = Tput;
        r.lossProbability = 1.0 - Tput / lambda;
        r.lambda = lambda;
        r.meanServiceTime = Smean;
        r.utilization = 1.0 - p0;
        r.rho = lambda * Smean;
        r.meanQueueLength = meanQ;
        r.nmax = nmax;
        r.countingResidual = sc.residual;
        Matrix sig = new Matrix(1, K * M);
        for (int i = 0; i < K * M; i++) {
            sig.set(0, i, sigma[i]);
        }
        r.sigma = sig;
        r.pKvec = pKvec;
        r.p0vec = p0vec;
        r.plevel = plevel;
        return r;
    }
}
