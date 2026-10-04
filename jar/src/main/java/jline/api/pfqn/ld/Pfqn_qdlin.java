package jline.api.pfqn.ld;

import jline.util.matrix.Matrix;

/**
 * QD-LIN: the Linearizer arm of AMVA-LD, on a plain demand matrix.
 *
 * <p>Array-level twin of what SolverMVA computes for {@code method='qdlin'}: the Linearizer of
 * Chandy and Neuse, Commun. ACM 25(2), 1982, run inside the queue-dependent AMVA framework of
 * Casale, Perez and Wang (IFIP PERFORMANCE 2015), so the load-dependent term g_k is evaluated at
 * the CORRECTED arrival-instant queue rather than at the plain one.
 *
 * <p>THIS IS A TRANSCRIPTION of {@code solver_amvald.m} together with
 * {@code solver_amvald_forward.m}, restricted to the domain a demand matrix describes: closed
 * classes only, one chain per class, unit visits, PS queueing stations and one optional delay
 * carrying Z. It is NOT an independent re-derivation. Ported from
 * {@code python/line_solver/api/pfqn/qdlin.py}, which reproduces the native-Python
 * {@code SolverMVA(model,'qdlin')} to machine precision over 640 random closed models.
 *
 * <p>FOUR PROPERTIES OF THE REFERENCE ARE REPRODUCED DELIBERATELY:
 *
 * <ol>
 *   <li>THE GAMMA CORRECTION IS CLASS-AGGREGATE, IN SLICE 0. {@code solver_amvald.m} allocates the
 *       (K,M,K) per-class Linearizer array for qdlin but writes
 *       {@code gamma(s,k) = sum_r Q_s(k,r)/(Nt-1) - sum_r Q(k,r)/Nt} into it with two subscripts,
 *       which MATLAB linear-indexes to (s,k,1); the other class slices stay zero while every
 *       reader indexes gamma per class. The correction that reaches the residence time is
 *       {@code N_0*gamma(r,k,0) - [r==0]*gamma(r,k,0)}, which coincides with the queue-dependent
 *       AMVA form {@code (Nt-1)*gamma_agg} iff K == 1. {@code method='lin'} takes the per-class
 *       form instead.</li>
 *   <li>A SINGLE-SERVER STATION STILL CARRIES A SOFTMIN TERM. The multiserver factor is
 *       {@link Pfqn_lldfun} at the arrival-instant total with the server counts, and its softmin at
 *       c = 1 is not exactly 1, so qdlin does not reduce to a textbook single-server AMVA even when
 *       every station has one server.</li>
 *   <li>THE WAIT FACTOR IS FLOORED AT {@code wtol}, LINE's {@code options.tol}. That is a DIFFERENT
 *       knob from the convergence tolerance and SolverMVA never sets it, so it stays at the
 *       lineDefaults 1e-4 while the fixed point converges to {@code iter_tol} 1e-6. MATLAB and the
 *       JAR solver do NOT carry this floor; native Python does, and removing it there was tried and
 *       reverted on 2026-09-04 because the unfloored Python recursion diverges where MATLAB and C++
 *       do not. See {@code _kb/06-solver-catalog.md}.</li>
 *   <li>WHICH UTILIZATION IS REPORTED DEPENDS ON THE MODEL. The analyzer forwards the iterated
 *       Uchain to the deaggregation ONLY under lld, cd or jd scaling; with none of those the
 *       deaggregation recomputes T*S/c from the NOMINAL demand, and the two differ by the iteration
 *       residual.</li>
 * </ol>
 *
 * <p>MU AND NSERVERS ARE DIFFERENT MECHANISMS, unlike in {@link Pfqn_qdamva}, which folds the
 * multiserver curve into mu. Here mu is {@code sn.lldscaling}, an interpolated rate multiplier per
 * station, and nservers is the server count feeding the softmin. A c-server station is
 * {@code nservers(k)=c}, NOT a mu row of {@code min(1..smax,c)}; the latter reproduces
 * {@code Queue.setLoadDependence}, a different station.
 */
public final class Pfqn_qdlin {
    private Pfqn_qdlin() {}

    /** Under-relaxation parameter of solver_amvald. */
    private static final double OMICRON = 0.5;

    /** The fixed point {@link #pfqn_qdlin} reaches, and how many sweeps it took. */
    public static class Result {
        /** (M x R) mean queue lengths at the queueing stations. */
        public final Matrix Q;
        /** (M x R) per-class utilizations. */
        public final Matrix U;
        /** (M x R) per-class residence times. */
        public final Matrix R;
        /** (1 x R) per-class throughputs. */
        public final Matrix X;
        /** (1 x R) per-class cycle times, think time included. */
        public final Matrix C;
        /** Number of forward evaluations performed. */
        public final int iter;

        public Result(Matrix Q, Matrix U, Matrix R, Matrix X, Matrix C, int iter) {
            this.Q = Q;
            this.U = U;
            this.R = R;
            this.X = X;
            this.C = C;
            this.iter = iter;
        }
    }

    /**
     * QD-LIN with the reference's default tolerances and iteration cap.
     *
     * @param L        (M x R) service demand matrix, queueing stations only
     * @param N        (1 x R) population vector, finite
     * @param Z        (1 x R) think time vector, or null for none
     * @param mu       (M x smax) load-dependent rate multipliers, or null for none
     * @param nservers (M x 1) server counts, or null for one server everywhere
     * @return the fixed point
     */
    public static Result pfqn_qdlin(Matrix L, Matrix N, Matrix Z, Matrix mu, Matrix nservers) {
        return pfqn_qdlin(L, N, Z, mu, nservers, 1e-6, 1000, 1e-4);
    }

    /**
     * QD-LIN.
     *
     * @param L        (M x R) service demand matrix, queueing stations only
     * @param N        (1 x R) population vector, finite
     * @param Z        (1 x R) think time vector, or null for none; a delay station carrying it is
     *                 prepended to the station list when any entry is positive, exactly as the
     *                 equivalent Network would hold one
     * @param mu       (M x smax) load-dependent rate multipliers, or null for none
     * @param nservers (M x 1) server counts, or null for one server everywhere
     * @param tol      convergence tolerance on the queue lengths, LINE's iter_tol
     * @param maxiter  iteration budget, LINE's iter_max; the outer sweep and each inner sweep are
     *                 capped at sqrt(maxiter) and the forward evaluations at min(maxiter,10000)
     * @param wtol     floor on the wait factor, LINE's options.tol
     * @return the fixed point
     */
    public static Result pfqn_qdlin(Matrix L, Matrix N, Matrix Z, Matrix mu, Matrix nservers,
                                    double tol, int maxiter, double wtol) {
        final int M = L.getNumRows();
        final int K = L.getNumCols();
        if (N == null || N.length() != K) {
            throw new IllegalArgumentException(
                    "pfqn_qdlin: the population vector must have one entry per class");
        }
        for (int r = 0; r < K; r++) {
            if (Double.isInfinite(N.get(r))) {
                throw new IllegalArgumentException(
                        "pfqn_qdlin: an infinite population is not supported, closed classes only");
            }
        }
        if (Z != null && !Z.isEmpty() && Z.length() != K) {
            throw new IllegalArgumentException(
                    "pfqn_qdlin: the think-time vector must have one entry per class");
        }
        if (nservers != null && !nservers.isEmpty() && nservers.length() != M) {
            throw new IllegalArgumentException(
                    "pfqn_qdlin: the server-count vector must have one entry per station");
        }

        double Nt = 0.0;
        for (int r = 0; r < K; r++) Nt += N.get(r);
        if (!(Nt > 0.0)) {
            return new Result(new Matrix(M, K), new Matrix(M, K), new Matrix(M, K),
                    new Matrix(1, K), new Matrix(1, K), 0);
        }

        // station list: the delay, when there is one, then the queueing stations
        boolean hasDelay = false;
        if (Z != null && !Z.isEmpty()) {
            for (int r = 0; r < K; r++) {
                if (Z.get(r) > 0.0) hasDelay = true;
            }
        }
        final int off = hasDelay ? 1 : 0;
        final int Ms = M + off;

        Matrix ST = new Matrix(Ms, K);
        double[] srv = new double[Ms];
        boolean[] isdelay = new boolean[Ms];
        if (hasDelay) {
            for (int r = 0; r < K; r++) ST.set(0, r, Z.get(r));
            srv[0] = Double.POSITIVE_INFINITY;
            isdelay[0] = true;
        }
        for (int k = 0; k < M; k++) {
            for (int r = 0; r < K; r++) ST.set(k + off, r, L.get(k, r));
            srv[k + off] = (nservers == null || nservers.isEmpty()) ? 1.0 : nservers.get(k);
        }
        Matrix srvMat = new Matrix(Ms, 1);
        for (int k = 0; k < Ms; k++) srvMat.set(k, 0, srv[k]);

        final boolean hasMu = mu != null && !mu.isEmpty();
        Matrix muFull = new Matrix(0, 0);
        if (hasMu) {
            final int smax = mu.getNumCols();
            muFull = new Matrix(Ms, smax);
            for (int j = 0; j < smax; j++) {
                if (hasDelay) muFull.set(0, j, 1.0);
                for (int k = 0; k < M; k++) muFull.set(k + off, j, mu.get(k, j));
            }
        }

        // the classes with jobs; every loop below runs over these
        int nnzCount = 0;
        for (int r = 0; r < K; r++) {
            if (N.get(r) > 0.0) nnzCount++;
        }
        int[] nnz = new int[nnzCount];
        int at = 0;
        for (int r = 0; r < K; r++) {
            if (N.get(r) > 0.0) nnz[at++] = r;
        }

        // balanced initialization, as in solver_amvald
        Matrix Q = new Matrix(Ms, K);
        final double share = 1.0 / Ms;
        for (int i = 0; i < nnz.length; i++) {
            final int r = nnz[i];
            for (int k = 0; k < Ms; k++) Q.set(k, r, share * N.get(r));
        }
        double[] X = new double[K];
        for (int r = 0; r < K; r++) {
            double col = 0.0;
            for (int k = 0; k < Ms; k++) col += ST.get(k, r);
            X[r] = col > 0.0 ? 1.0 / col : 0.0;
        }
        Matrix U = new Matrix(Ms, K);
        for (int k = 0; k < Ms; k++) {
            for (int i = 0; i < nnz.length; i++) {
                final int r = nnz[i];
                U.set(k, r, Double.isInfinite(srv[k]) ? ST.get(k, r) * X[r]
                                                      : ST.get(k, r) * X[r] / srv[k]);
            }
        }

        // gamma(s,i,r) as one (Ms x K) matrix per reduced class s
        Matrix[] gamma = new Matrix[K];
        for (int s = 0; s < K; s++) gamma[s] = new Matrix(Ms, K);
        Matrix Tput = new Matrix(Ms, K);
        Matrix Cout = new Matrix(1, K);
        Matrix STeff = new Matrix(Ms, K);

        final double maxSweep = Math.sqrt(maxiter);
        final int maxTotiter = Math.min(maxiter, 10000);
        int totiter = 0;
        int outerIter = 0;
        Matrix QouterPrev = null;

        while (true) {
            if (outerIter >= 2) {
                double gap = 0.0;
                for (int k = 0; k < Ms; k++) {
                    for (int r = 0; r < K; r++) {
                        gap = Math.max(gap, Math.abs(Q.get(k, r) - QouterPrev.get(k, r)));
                    }
                }
                if (!(gap > tol)) break;
            }
            if (!(outerIter < maxSweep) || totiter > maxTotiter) break;
            outerIter++;
            QouterPrev = Q.copy();

            // Linearizer recursion: one sweep at each reduced population N - 1_s
            boolean exhausted = false;
            for (int s = 0; s < K && !exhausted; s++) {
                if (!(N.get(s) > 0.0)) continue;
                double[] Ns = new double[K];
                for (int r = 0; r < K; r++) Ns[r] = N.get(r);
                Ns[s] -= 1.0;
                final double shrink = (Nt - 1.0) / Nt;
                Matrix Qs = Q.copy();
                for (int k = 0; k < Ms; k++) {
                    for (int r = 0; r < K; r++) Qs.set(k, r, Qs.get(k, r) * shrink);
                }
                double[] Xs = new double[K];
                for (int r = 0; r < K; r++) Xs[r] = X[r] * shrink;

                int iterS = 0;
                Matrix QsPrev = null;
                while (true) {
                    if (iterS >= 2) {
                        double gap = 0.0;
                        for (int k = 0; k < Ms; k++) {
                            for (int r = 0; r < K; r++) {
                                gap = Math.max(gap, Math.abs(Qs.get(k, r) - QsPrev.get(k, r)));
                            }
                        }
                        if (!(gap > tol)) break;
                    }
                    if (!(iterS <= maxSweep)) break;
                    iterS++;
                    QsPrev = Qs.copy();
                    final double[] XsPrev = new double[K];
                    System.arraycopy(Xs, 0, XsPrev, 0, K);

                    Matrix[] fw = forward(ST, srvMat, srv, isdelay, muFull, gamma, QsPrev, Ns, nnz,
                            K, wtol);
                    final Matrix Ws = fw[0];
                    totiter++;
                    if (totiter >= maxTotiter) {
                        exhausted = true;
                        break;
                    }

                    for (int i = 0; i < nnz.length; i++) {
                        final int r = nnz[i];
                        double wsum = 0.0;
                        for (int k = 0; k < Ms; k++) wsum += Ws.get(k, r);
                        if (wsum == 0.0 || Ns[r] == 0.0) {
                            Xs[r] = 0.0;
                        } else if (wsum > 1e-14) {
                            Xs[r] = OMICRON * Ns[r] / wsum + (1 - OMICRON) * XsPrev[r];
                        } else {
                            Xs[r] = XsPrev[r];
                        }
                        for (int k = 0; k < Ms; k++) {
                            Qs.set(k, r, OMICRON * Xs[r] * Ws.get(k, r)
                                    + (1 - OMICRON) * QsPrev.get(k, r));
                        }
                    }
                }

                // class-aggregate correction into column 0, see the class comment
                if (Nt > 1.0) {
                    for (int k = 0; k < Ms; k++) {
                        double qs = 0.0;
                        double qo = 0.0;
                        for (int r = 0; r < K; r++) {
                            qs += QsPrev.get(k, r);
                            qo += QouterPrev.get(k, r);
                        }
                        gamma[s].set(k, 0, qs / (Nt - 1.0) - qo / Nt);
                    }
                } else {
                    for (int k = 0; k < Ms; k++) gamma[s].set(k, 0, 0.0);
                }
            }
            if (exhausted) break;

            // sweep at the full population N
            double[] Nfull = new double[K];
            for (int r = 0; r < K; r++) Nfull[r] = N.get(r);
            int innerIter = 0;
            Matrix Qprev = null;
            while (true) {
                if (innerIter >= 2) {
                    double gap = 0.0;
                    for (int k = 0; k < Ms; k++) {
                        for (int r = 0; r < K; r++) {
                            gap = Math.max(gap, Math.abs(Q.get(k, r) - Qprev.get(k, r)));
                        }
                    }
                    if (!(gap > tol)) break;
                }
                if (!(innerIter <= maxSweep)) break;
                innerIter++;
                Qprev = Q.copy();
                final double[] Xprev = new double[K];
                System.arraycopy(X, 0, Xprev, 0, K);
                final Matrix Uprev = U.copy();

                Matrix[] fw = forward(ST, srvMat, srv, isdelay, muFull, gamma, Qprev, Nfull, nnz,
                        K, wtol);
                final Matrix W = fw[0];
                STeff = fw[1];
                totiter++;
                if (totiter >= maxTotiter) {
                    exhausted = true;
                    break;
                }

                for (int i = 0; i < nnz.length; i++) {
                    final int r = nnz[i];
                    double wsum = 0.0;
                    for (int k = 0; k < Ms; k++) wsum += W.get(k, r);
                    if (wsum == 0.0) {
                        X[r] = 0.0;
                    } else {
                        Cout.set(0, r, wsum);
                        if (wsum > 1e-14) {
                            X[r] = OMICRON * Nfull[r] / wsum + (1 - OMICRON) * Xprev[r];
                        } else {
                            X[r] = Xprev[r];
                        }
                    }
                    for (int k = 0; k < Ms; k++) {
                        Q.set(k, r, OMICRON * X[r] * W.get(k, r) + (1 - OMICRON) * Qprev.get(k, r));
                        Tput.set(k, r, X[r]);
                        U.set(k, r, OMICRON * STeff.get(k, r) * X[r]
                                + (1 - OMICRON) * Uprev.get(k, r));
                    }
                }
            }
            if (exhausted) break;
        }

        // utilization capping, as in solver_amvald: a queueing station whose class utilizations
        // sum above one has them renormalized in proportion to STeff. Delay stations are exempt.
        for (int k = 0; k < Ms; k++) {
            if (isdelay[k]) continue;
            double usum = 0.0;
            for (int r = 0; r < K; r++) usum += U.get(k, r);
            if (!(usum > 1.0)) continue;
            double denom = 0.0;
            for (int r = 0; r < K; r++) denom += STeff.get(k, r) * X[r];
            if (!(denom > 0.0)) continue;
            for (int r = 0; r < K; r++) {
                if (STeff.get(k, r) > 0.0) {
                    U.set(k, r, Math.min(1.0, usum) * STeff.get(k, r) * X[r] / denom);
                }
            }
        }

        // the analyzer forwards the iterated Uchain to the deaggregation ONLY under lld, cd or jd
        // scaling; with none of those it recomputes T*S/c from the NOMINAL demand. mu is the only
        // one of the three a demand matrix can carry.
        if (!hasMu) {
            for (int k = 0; k < Ms; k++) {
                for (int i = 0; i < nnz.length; i++) {
                    final int r = nnz[i];
                    U.set(k, r, Double.isInfinite(srv[k]) ? ST.get(k, r) * X[r]
                                                          : ST.get(k, r) * X[r] / srv[k]);
                }
            }
        }

        Matrix Qout = new Matrix(M, K);
        Matrix Uout = new Matrix(M, K);
        Matrix Rout = new Matrix(M, K);
        // A class with no jobs keeps its 1/sum(ST) SEED in X unless it is left out
        // here: the sweeps only ever write the classes in nnz, so the initial value
        // would otherwise be reported as that class's throughput.
        Matrix Xout = new Matrix(1, K);
        for (int i = 0; i < nnz.length; i++) Xout.set(0, nnz[i], X[nnz[i]]);
        for (int k = 0; k < M; k++) {
            final int ks = k + off;
            for (int r = 0; r < K; r++) {
                Qout.set(k, r, Q.get(ks, r));
                Uout.set(k, r, U.get(ks, r));
                Rout.set(k, r, Tput.get(ks, r) > 0.0 ? Q.get(ks, r) / Tput.get(ks, r) : 0.0);
            }
        }
        return new Result(Qout, Uout, Rout, Xout, Cout, totiter);
    }

    /**
     * One forward evaluation, solver_amvald_forward restricted to PS and INF.
     *
     * @return {W, STeff}, both (Ms x K)
     */
    private static Matrix[] forward(Matrix ST, Matrix srvMat, double[] srv, boolean[] isdelay,
                                    Matrix mu, Matrix[] gamma, Matrix Qin, double[] Nin, int[] nnz,
                                    int K, double wtol) {
        final int Ms = ST.getNumRows();
        double Ntin = 0.0;
        for (int r = 0; r < K; r++) Ntin += Nin[r];
        final double delta = Ntin > 0.0 ? (Ntin - 1.0) / Ntin : 1.0;

        double[] dcl = new double[K];
        for (int r = 0; r < K; r++) dcl[r] = 1.0;
        for (int i = 0; i < nnz.length; i++) {
            final int r = nnz[i];
            dcl[r] = (Nin[r] - 1.0) / Nin[r];
        }

        // arrival-instant queue lengths, class-aggregate and per class
        double[] interp = new double[Ms];
        Matrix totArvl = new Matrix(Ms, K);
        for (int k = 0; k < Ms; k++) {
            double sumQk = 0.0;
            for (int i = 0; i < nnz.length; i++) sumQk += Qin.get(k, nnz[i]);
            interp[k] = delta * sumQk;
            for (int i = 0; i < nnz.length; i++) {
                final int r = nnz[i];
                totArvl.set(k, r, dcl[r] * Qin.get(k, r) + sumQk - Qin.get(k, r));
            }
        }

        // the load-dependent term, at the gamma-corrected arrival-instant queue
        Matrix lldterm = new Matrix(Ms, K);
        lldterm.fill(1.0);
        final Matrix noServers = new Matrix(0, 0);
        if (nnz.length > 0) {
            for (int i = 0; i < nnz.length; i++) {
                final int r = nnz[i];
                Matrix arg = new Matrix(Ms, 1);
                for (int k = 0; k < Ms; k++) {
                    double corr = 0.0;
                    for (int j = 0; j < nnz.length; j++) {
                        corr += Nin[nnz[j]] * gamma[r].get(k, nnz[j]);
                    }
                    corr -= gamma[r].get(k, r);
                    arg.set(k, 0, 1.0 + interp[k] + corr);
                }
                Matrix v = Pfqn_lldfun.pfqn_lldfun(arg, mu, noServers);
                for (int k = 0; k < Ms; k++) lldterm.set(k, r, v.get(k, 0));
            }
        } else {
            Matrix arg = new Matrix(Ms, 1);
            for (int k = 0; k < Ms; k++) arg.set(k, 0, 1.0 + interp[k]);
            Matrix v = Pfqn_lldfun.pfqn_lldfun(arg, mu, noServers);
            for (int k = 0; k < Ms; k++) {
                for (int r = 0; r < K; r++) lldterm.set(k, r, v.get(k, 0));
            }
        }

        // the multiserver term; the 'default' rule leaves PS on the softmin arm
        Matrix msarg = new Matrix(Ms, 1);
        if (nnz.length > 0 && Ntin > 0.0) {
            final double shrink = (Ntin - 1.0) / Ntin;
            for (int k = 0; k < Ms; k++) {
                double acc = 0.0;
                for (int a = 0; a < nnz.length; a++) {
                    final int s = nnz[a];
                    double gs = 0.0;
                    for (int b = 0; b < nnz.length; b++) {
                        final int r = nnz[b];
                        gs += shrink * Nin[r] * gamma[s].get(k, r);
                    }
                    acc += gs;
                }
                msarg.set(k, 0, 1.0 + interp[k] + acc / nnz.length);
            }
        } else {
            for (int k = 0; k < Ms; k++) msarg.set(k, 0, 1.0 + interp[k]);
        }
        final Matrix noLld = new Matrix(0, 0);
        Matrix msterm = Pfqn_lldfun.pfqn_lldfun(msarg, noLld, srvMat);

        Matrix STeff = new Matrix(Ms, K);
        for (int i = 0; i < nnz.length; i++) {
            final int r = nnz[i];
            for (int k = 0; k < Ms; k++) {
                STeff.set(k, r, ST.get(k, r) * lldterm.get(k, r) * msterm.get(k, 0));
            }
        }

        Matrix W = new Matrix(Ms, K);
        for (int i = 0; i < nnz.length; i++) {
            final int r = nnz[i];
            for (int k = 0; k < Ms; k++) {
                if (isdelay[k]) {
                    W.set(k, r, STeff.get(k, r));
                } else {
                    double corr = 0.0;
                    for (int j = 0; j < nnz.length; j++) {
                        corr += Nin[nnz[j]] * gamma[r].get(k, nnz[j]);
                    }
                    corr -= gamma[r].get(k, r);
                    W.set(k, r, STeff.get(k, r)
                            * Math.max(wtol, 1.0 + totArvl.get(k, r) + corr));
                }
            }
        }
        return new Matrix[] {W, STeff};
    }
}
