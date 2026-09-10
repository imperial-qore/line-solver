/**
 * @file Mean value analysis of closed networks with shortest-job-next stations
 *
 * Approximate MVA for closed queueing networks in which a subset of the single-server stations
 * schedules non-preemptively by shortest job next (SJN/SJF), the job size being known on arrival.
 * Ported at parity from MATLAB pfqn_mvasjn.m.
 *
 * Reference: K. Kant, "MVA approximations for SJN scheduling", Performance Evaluation 15(1):41-61,
 * 1992.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.util.matrix.Matrix;

import java.util.Arrays;

public final class Pfqn_mvasjn {
    private Pfqn_mvasjn() {}

    /** Conditional waiting time profile of one SJN station at the target population. */
    public static final class Profile {
        /** Station index within the rows of L. */
        public final int station;
        /** Job size grid. */
        public final double[] x;
        /** Conditional waiting times (ngrid x R). */
        public final double[][] W;
        /** Tail parameters a, b, c per class (R x 3). */
        public final double[][] tail;

        public Profile(int station, double[] x, double[][] W, double[][] tail) {
            this.station = station;
            this.x = x;
            this.W = W;
            this.tail = tail;
        }
    }

    /** Mean performance measures of the network. */
    public static final class Result {
        /** Throughputs (1 x R). */
        public final Matrix X;
        /** Mean queue lengths (M x R). */
        public final Matrix Q;
        /** Utilizations (M x R). */
        public final Matrix U;
        /** Residence times (M x R). */
        public final Matrix C;
        /** Conditional waiting time profiles, one per SJN station. */
        public final Profile[] profiles;
        /** Iterations performed, one for the lattice recursion. */
        public final int iter;

        public Result(Matrix X, Matrix Q, Matrix U, Matrix C, Profile[] profiles, int iter) {
            this.X = X;
            this.Q = Q;
            this.U = U;
            this.C = C;
            this.profiles = profiles;
            this.iter = iter;
        }
    }

    /**
     * Mean value analysis with shortest-job-next stations, over the whole population lattice.
     *
     * <p>The SJN station is modelled by the conditional waiting time W(x,n) of a tagged customer
     * whose service requirement is x, obtained from the arrival theorem as the sum of the residual
     * life of the job in service, the work of the queued jobs that will be served before the
     * tagged one, and the work of the jobs that overtake it while it waits:</p>
     *
     * <pre>
     *   W(x,n) = [ (1+CV^2) s U(n-1)/2 + X(n-1) phi(x,n-1) ] / [ 1 - X(n-1) theta(x) ]
     *   theta(x) = int_0^x t f(t) dt,   phi(x,n) = int_0^x W(t,n) t f(t) dt
     *   R(n) = s + int_0^inf W(x,n) f(x) dx
     * </pre>
     *
     * <p>The recursion is explicit: W(.,n) needs only phi(.,n-1), so it is carried alongside the
     * population recursion of exact MVA. This costs prod(N+1) steps;
     * {@link Pfqn_amvasjn} is the fixed-point counterpart that trades the lattice for a
     * Schweitzer closure on the same profile.</p>
     *
     * <p>The service time density is not an input: only its mean and squared coefficient of
     * variation are, and the density is reconstructed by the two-moment Erlang-mixture fit the
     * reference prescribes. The x-integrals run on a fixed grid by composite Simpson, W(.,n) being
     * needed at the next population so that quadrature rules sampling at arbitrary abscissae
     * cannot be used; beyond the grid the profile is closed by the analytic tail
     * W(x,n) = a - b exp(-c (x - Lx)).</p>
     *
     * @param L service demand matrix (M x R) of the queueing stations
     * @param N population vector (1 x R)
     * @param Z think time vector (1 x R), may be null
     * @param scv squared coefficients of variation of the service times (M x R), may be null
     * @param sjnset zero-based indices of the stations scheduling by SJN, may be null
     * @param V visit ratios (M x R), so that the per-visit service time is L./V, may be null
     * @param options grid, priority and cap options, may be null
     * @return mean performance measures and the conditional waiting time profiles
     */
    public static Result pfqn_mvasjn(Matrix L, Matrix N, Matrix Z, Matrix scv, int[] sjnset,
            Matrix V, SjnOptions options) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        if (options == null) {
            options = new SjnOptions();
        }
        options.validate(R);
        double[] Nv = SjnArgs.vector(N, R, "population vector");
        for (int r = 0; r < R; r++) {
            Nv[r] = Math.round(Nv[r]);
            if (Nv[r] < 0) {
                throw new IllegalArgumentException("negative class populations");
            }
        }
        double[] Zv = Z == null || Z.isEmpty() ? new double[R] : SjnArgs.vector(Z, R, "think times");
        double[][] Lm = SjnArgs.matrix(L, M, R);
        double[][] Vm = V == null || V.isEmpty() ? SjnArgs.ones(M, R) : SjnArgs.matrix(V, M, R);
        double[][] scvm = scv == null || scv.isEmpty() ? SjnArgs.ones(M, R)
                : SjnArgs.matrix(scv, M, R);
        double[][] S = SjnArgs.perVisit(Lm, Vm);
        int[] sjn = SjnArgs.stationSet(sjnset, M);
        boolean useprio = options.prio != null;

        int nsjn = sjn.length;
        int ngrid = options.ns + 1;
        SjnSupport.Grid[] G = new SjnSupport.Grid[nsjn];
        for (int q = 0; q < nsjn; q++) {
            G[q] = SjnSupport.setup(S[sjn[q]], scvm[sjn[q]], options.ns, options.Lfactor);
        }

        int[] stride = new int[R];
        int npop = 1;
        for (int r = 0; r < R; r++) {
            stride[r] = npop;
            npop *= (int) Nv[r] + 1;
        }
        double[][] Xp = new double[npop][R];
        double[][][] Qp = new double[npop][M][R];
        double[][][] Up = new double[npop][M][R];
        double[][][] Cp = new double[npop][M][R];
        double[][][][] Wp = new double[nsjn][npop][ngrid][R];
        double[][][][] Pp = new double[nsjn][npop][ngrid][R];
        double[][][] Ip = new double[nsjn][npop][R];
        double[][][] Tp = new double[nsjn][npop][R * 3];

        double[] beta = SjnArgs.onesVector(R);
        for (int idx = 1; idx < npop; idx++) {
            double[] n = SjnArgs.decode(idx, stride, Nv);
            double[][] Call = new double[M][R];
            for (int r = 0; r < R; r++) {
                if (n[r] == 0) {
                    continue;
                }
                int iprev = idx - stride[r];
                for (int m = 0; m < M; m++) {
                    int q = SjnArgs.indexOf(sjn, m);
                    if (q < 0) {
                        double qsum = 0;
                        for (int k = 0; k < R; k++) {
                            qsum += Qp[iprev][m][k];
                        }
                        Call[m][r] = Lm[m][r] * (1 + qsum);
                        continue;
                    }
                    // the population step already supplies the neighbouring profile, no deflation
                    double[] lam = new double[R];
                    for (int k = 0; k < R; k++) {
                        lam[k] = Xp[iprev][k] * Vm[m][k];
                    }
                    SjnSupport.State st = new SjnSupport.State(lam, Up[iprev][m], Qp[iprev][m],
                            Wp[q][iprev], Pp[q][iprev], Ip[q][iprev]);
                    SjnSupport.StationResult sr = SjnSupport.station(m, r, G[q], S[m], scvm[m],
                            Vm[m], st, beta, useprio, options.prio);
                    Call[m][r] = sr.C;
                    for (int i = 0; i < ngrid; i++) {
                        Wp[q][idx][i][r] = sr.W[i];
                        Pp[q][idx][i][r] = sr.phi[i];
                    }
                    Ip[q][idx][r] = sr.phiinf;
                    Tp[q][idx][3 * r] = sr.tail[0];
                    Tp[q][idx][3 * r + 1] = sr.tail[1];
                    Tp[q][idx][3 * r + 2] = sr.tail[2];
                }
            }
            SjnSupport.CapResult cr = SjnSupport.cap(Call, Lm, n, Zv, sjn, options.umax);
            if (cr.bound) {
                // the cap has invalidated the profile the next population step reads back
                throw new SjnStarvationException(String.format(
                        "the utilization cap of %g was binding at an SJN station at population %s: "
                        + "the station is in the starvation regime, where the conditional waiting "
                        + "time equation has no solution and the population lattice no valid "
                        + "continuation. Use the Schweitzer fixed point (pfqn_amvasjn, method "
                        + "'amva'), SolverCTMC or SolverLDES.",
                        options.umax, Arrays.toString(n)));
            }
            for (int r = 0; r < R; r++) {
                Xp[idx][r] = cr.X[r];
                for (int m = 0; m < M; m++) {
                    Cp[idx][m][r] = cr.C[m][r];
                    Qp[idx][m][r] = cr.X[r] * cr.C[m][r];
                    Up[idx][m][r] = cr.X[r] * Lm[m][r];
                }
            }
        }
        int last = npop - 1;
        Profile[] profiles = new Profile[nsjn];
        for (int q = 0; q < nsjn; q++) {
            double[][] tail = new double[R][3];
            for (int r = 0; r < R; r++) {
                tail[r][0] = Tp[q][last][3 * r];
                tail[r][1] = Tp[q][last][3 * r + 1];
                tail[r][2] = Tp[q][last][3 * r + 2];
            }
            profiles[q] = new Profile(sjn[q], G[q].x, Wp[q][last], tail);
        }
        return new Result(SjnArgs.rowMatrix(Xp[last]), SjnArgs.toMatrix(Qp[last]),
                SjnArgs.toMatrix(Up[last]), SjnArgs.toMatrix(Cp[last]), profiles, 1);
    }
}
