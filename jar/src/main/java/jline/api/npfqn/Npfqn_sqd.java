/**
 * Smith Queue Decomposition (SQD): approximate MVA for closed Blocking-After-Service (BAS) networks.
 *
 * Solves a finite-buffer closed queueing network under Blocking-After-Service
 * (manufacturing/transfer) blocking directly from its {@link NetworkStruct}. The
 * method is an AMVA-style population recursion in which each finite-capacity station
 * is described by a load-dependent effective service rate calibrated from an M/M/1/K
 * blocking probability; downstream blocking is propagated through the effective
 * routing between service stations. Delay (INF/EXT) stations are treated as
 * infinite-capacity pure-delay nodes.
 *
 * Single-chain (chain-aggregated) demands only. Returns per-station throughput,
 * queue length, utilization, and residence time.
 *
 * Originally contributed as {@code SolverDBT} by Avinash Bommareddy (Imperial College
 * London FYP, 2026); refactored here into an {@code sn}-based API function.
 *
 * @since LINE 3.0
 */
package jline.api.npfqn;

import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.util.matrix.Matrix;

import static jline.api.sn.SnGetDemandsChain.snGetDemandsChain;

public final class Npfqn_sqd {

    private Npfqn_sqd() {}

    /** Selects how downstream blocking is aggregated for a station. */
    public enum NeighborMode { DOWNSTREAM, OWNSERVER }

    /** Selects how the load-dependent rate scale V1 is updated across populations. */
    public enum V1Policy { COMPOUND, FRESH }

    private static final int     DEFAULT_CALIBRATION_MODE   = 0;
    private static final boolean DEFAULT_SERVER_BLOCKING    = true;

    private static final double INITIAL_V1          = 692.192;
    private static final double CALIBRATION_EPSILON = 0.05;

    /**
     * Solve a BAS closed network from its NetworkStruct using the default population
     * {@code sn.nclosedjobs} and default tuning (calibration mode 0, server blocking
     * time on, downstream neighbor mode, compound V1 policy).
     */
    public static Ret.pfqnMVA npfqn_sqd(NetworkStruct sn) {
        return npfqn_sqd(sn, sn.nclosedjobs);
    }

    /** Solve a BAS closed network from its NetworkStruct for a given total population N. */
    public static Ret.pfqnMVA npfqn_sqd(NetworkStruct sn, int N) {
        return npfqn_sqd(sn, N, DEFAULT_CALIBRATION_MODE, DEFAULT_SERVER_BLOCKING,
                NeighborMode.DOWNSTREAM, V1Policy.COMPOUND, null);
    }

    /**
     * Solve a BAS closed network from its NetworkStruct.
     *
     * @param sn                 the chain-aggregated network structure
     * @param N                  total closed-class population
     * @param calibrationMode    0 = base (beta=K, gamma=1), 1 = fixed heuristic,
     *                           2 = blocking-aware (beta,gamma) calibration
     * @param serverBlockingTime if true, add a manufacturing-blocking term to server time
     * @param neighborMode       downstream-routed vs own-server blocking aggregation
     * @param v1Policy           compound vs fresh load-dependent rate-scale update
     * @param initialV1          optional per-station initial V1 (null = INITIAL_V1)
     * @return per-station throughput X, queue length Q, utilization U, residence time R
     */
    public static Ret.pfqnMVA npfqn_sqd(NetworkStruct sn, int N,
                                        int calibrationMode, boolean serverBlockingTime,
                                        NeighborMode neighborMode, V1Policy v1Policy,
                                        double[] initialV1) {
        int M = sn.nstations;
        int R = sn.nclasses;

        Ret.snGetDemands demands = snGetDemandsChain(sn);
        Matrix Vchain  = demands.Vchain;
        Matrix STchain = demands.STchain;

        double[] V  = new double[M];
        double[] ST = new double[M];
        for (int i = 0; i < M; i++) {
            V[i]  = Vchain.get(i, 0);
            ST[i] = STchain.get(i, 0);
        }

        boolean[] isDelay = new boolean[M];
        int[]     cap     = new int[M];
        for (int i = 0; i < M; i++) {
            Station st = sn.stations.get(i);
            SchedStrategy sched = sn.sched.get(st);
            isDelay[i] = (sched == SchedStrategy.INF || sched == SchedStrategy.EXT);
            if (isDelay[i]) {
                cap[i] = Integer.MAX_VALUE;
            } else {
                double c = sn.cap.get(i, 0);
                cap[i] = (Double.isInfinite(c) || c > 1e14) ? Integer.MAX_VALUE : (int) c;
            }
        }

        double[][] pEff = computeEffectiveRouting(sn, M, R, isDelay);

        // see _kb/03-api-layer.md for rationale
        if (M == 2 && !isDelay[0] && !isDelay[1]
                && cap[0] < Integer.MAX_VALUE && cap[1] < Integer.MAX_VALUE
                && pEff[0][1] > 1.0 - 1e-9 && pEff[1][0] > 1.0 - 1e-9) {
            int C = cap[0] + cap[1] + 2;
            int Ndual = C - N;
            if (2 * N > C && Ndual >= 1 && Ndual < N) {
                Ret.pfqnMVA dual = npfqn_sqd(sn, Ndual, calibrationMode, serverBlockingTime,
                        neighborMode, v1Policy, initialV1);
                Matrix XN = new Matrix(M, 1);
                Matrix QN = new Matrix(M, 1);
                Matrix UN = new Matrix(M, 1);
                Matrix WN = new Matrix(M, 1);
                for (int i = 0; i < M; i++) {
                    double T_i = dual.X.get(i, 0);
                    double Q_i = (cap[i] + 1) - dual.Q.get(1 - i, 0);
                    double U_i = Math.min(1.0, T_i * ST[i]);
                    double W_i = (T_i > 1e-15) ? Q_i / T_i : 0.0;
                    XN.set(i, 0, T_i);
                    QN.set(i, 0, Q_i);
                    UN.set(i, 0, U_i);
                    WN.set(i, 0, W_i);
                }
                return new Ret.pfqnMVA(XN, QN, UN, WN, Double.NaN);
            }
        }

        double[] V1     = new double[M];
        double[] v1init = new double[M];
        double[] L_buf  = new double[M];
        double[] L_svr  = new double[M];
        for (int i = 0; i < M; i++) {
            v1init[i] = (initialV1 == null) ? INITIAL_V1 : initialV1[i];
            V1[i] = v1init[i];
        }

        double   X     = 0.0;
        double[] W_buf = new double[M];
        double[] W_svr = new double[M];

        for (int pop = 1; pop <= N; pop++) {

            // wait times
            for (int i = 0; i < M; i++) {
                if (isDelay[i]) {
                    W_buf[i] = 0.0;
                    W_svr[i] = ST[i];
                    continue;
                }

                double pBlockDown;
                if (neighborMode == NeighborMode.OWNSERVER) {
                    pBlockDown = (cap[i] < Integer.MAX_VALUE)
                            ? mm1kBlocking(cap[i], X * V[i] * ST[i])
                            : 0.0;
                } else {
                    pBlockDown = 0.0;
                    for (int j = 0; j < M; j++) {
                        if (!isDelay[j] && pEff[i][j] > 0 && cap[j] < Integer.MAX_VALUE) {
                            double rho_j = X * V[j] * ST[j];
                            pBlockDown += pEff[i][j] * mm1kBlocking(cap[j], rho_j);
                        }
                    }
                }

                double n;
                if (neighborMode == NeighborMode.OWNSERVER) {
                    n = L_svr[i];
                } else {
                    n = 0.0;
                    for (int j = 0; j < M; j++) {
                        n += pEff[i][j] * L_svr[j];
                    }
                }

                double mu_n;
                if (n < 1e-10) {
                    mu_n = V1[i];
                } else {
                    double[] bg     = computeBetaGamma(cap[i], pBlockDown, calibrationMode);
                    double   beta   = bg[0];
                    double   gamma  = bg[1];
                    double   base   = Math.max(0.0, (n - 1.0) / beta);
                    double   expArg = Math.pow(base, gamma);
                    mu_n = n * V1[i] * Math.exp(-expArg);   // Eq.13
                }
                if (mu_n < 1e-10) mu_n = 1e-10;

                W_buf[i] = (1.0 / mu_n) * (1.0 + n);        // Eq.18
                W_svr[i] = ST[i] * (1.0 + L_svr[i]);        // Eq.17

                if (serverBlockingTime) {
                    double bt = 0.0;
                    for (int j = 0; j < M; j++) {
                        if (!isDelay[j] && pEff[i][j] > 0 && cap[j] < Integer.MAX_VALUE) {
                            double rho_j = X * V[j] * ST[j];
                            double pBj   = mm1kBlocking(cap[j], rho_j);
                            double denom = ST[i] + ST[j];
                            double theta = (denom > 1e-15) ? ST[j] / denom : 0.0;
                            bt += pEff[i][j] * pBj * ST[j] * theta;
                        }
                    }
                    W_svr[i] += bt;
                }
            }

            // throughput
            double sumVW = 0.0;
            for (int i = 0; i < M; i++) {
                sumVW += V[i] * (W_buf[i] + W_svr[i]);
            }
            X = (sumVW > 1e-15) ? (double) pop / sumVW : 0.0;

            // queue lengths
            for (int i = 0; i < M; i++) {
                L_buf[i] = X * V[i] * W_buf[i];
                L_svr[i] = X * V[i] * W_svr[i];
            }

            // adjust V1
            if (pop < N) {
                for (int i = 0; i < M; i++) {
                    if (isDelay[i]) continue;
                    double pBlock;
                    if (neighborMode == NeighborMode.OWNSERVER) {
                        pBlock = (cap[i] < Integer.MAX_VALUE)
                                ? mm1kBlocking(cap[i], X * V[i] * ST[i])
                                : 0.0;
                    } else {
                        pBlock = 0.0;
                        for (int j = 0; j < M; j++) {
                            if (!isDelay[j] && pEff[i][j] > 0 && cap[j] < Integer.MAX_VALUE) {
                                double rho_j = X * V[j] * ST[j];
                                pBlock += pEff[i][j] * mm1kBlocking(cap[j], rho_j);
                            }
                        }
                    }
                    if (v1Policy == V1Policy.FRESH) {
                        V1[i] = v1init[i] * (1.0 - pBlock);
                    } else {
                        V1[i] = V1[i] * (1.0 - pBlock);
                    }
                    if (V1[i] < 1e-10) V1[i] = 1e-10;
                }
            }
        }

        Matrix XN = new Matrix(M, 1);
        Matrix QN = new Matrix(M, 1);
        Matrix UN = new Matrix(M, 1);
        Matrix WN = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            double T_i   = X * V[i];
            double Q_i   = L_buf[i] + L_svr[i];
            double W_tot = W_buf[i] + W_svr[i];
            double R_i   = (T_i > 1e-15) ? Q_i / T_i : W_tot;
            double U_i   = Math.min(1.0, T_i * ST[i]);
            XN.set(i, 0, T_i);
            QN.set(i, 0, Q_i);
            UN.set(i, 0, U_i);
            WN.set(i, 0, R_i);
        }

        return new Ret.pfqnMVA(XN, QN, UN, WN, Double.NaN);
    }

    /**
     * Calibrate the (beta, gamma) shape parameters of the load-dependent effective
     * service rate from the station capacity K and a downstream blocking probability.
     */
    static double[] computeBetaGamma(int K, double pBlockDown, int mode) {
        if (K <= 2 || K == Integer.MAX_VALUE) {
            return new double[]{K, 1.0};
        }

        final double a = 2.0;
        final double b = (double) K;

        double Va, Vb;
        switch (mode) {
            case 2: {   // blocking-aware
                double pK = Math.max(1e-6, Math.min(pBlockDown, 1.0 - 1e-6));
                Va = 1.0 - pK * (a - 1.0) / (b - 1.0);
                Vb = 1.0 - pK;
                break;
            }
            case 1: {   // fixed heuristic
                Va = (b - a) / b;
                Vb = CALIBRATION_EPSILON;
                break;
            }
            default:    // mode 0: base
                return new double[]{K, 1.0};
        }

        Va = Math.min(Va, 0.999);
        Va = Math.max(Va, 0.01);
        Vb = Math.max(Vb, 1e-6);
        if (Vb >= Va) {
            return new double[]{K, 1.0};
        }

        double lnVa = Math.log(Va);
        double lnVb = Math.log(Vb);
        double gamma = Math.log(lnVa / lnVb) / Math.log((a - 1.0) / (b - 1.0));
        gamma = Math.max(0.5, Math.min(gamma, 10.0));
        double beta = (a - 1.0) / Math.pow(-lnVa, 1.0 / gamma);

        if (Double.isNaN(gamma) || Double.isInfinite(gamma) ||
                Double.isNaN(beta)  || Double.isInfinite(beta)  || beta <= 0) {
            return new double[]{K, 1.0};
        }

        return new double[]{beta, gamma};
    }

    /**
     * Build station-to-station effective routing among service stations, collapsing
     * pass-through delay (INF/EXT) stations into a single hop.
     */
    private static double[][] computeEffectiveRouting(NetworkStruct sn, int M, int R,
                                                      boolean[] isDelay) {
        double[][] p = new double[M][M];
        for (int i = 0; i < M; i++) {
            if (isDelay[i]) continue;
            int sf_i = (int) sn.stationToStateful.get(i);
            for (int j = 0; j < M; j++) {
                int    sf_j = (int) sn.stationToStateful.get(j);
                double p_ij = sn.rt.get(sf_i * R, sf_j * R);
                if (p_ij <= 0) continue;
                if (!isDelay[j]) {
                    p[i][j] += p_ij;
                } else {
                    for (int k = 0; k < M; k++) {
                        if (!isDelay[k]) {
                            int    sf_k = (int) sn.stationToStateful.get(k);
                            double p_jk = sn.rt.get(sf_j * R, sf_k * R);
                            if (p_jk > 0) p[i][k] += p_ij * p_jk;
                        }
                    }
                }
            }
        }
        return p;
    }

    /** Steady-state blocking probability of an M/M/1/K queue at load rho. */
    private static double mm1kBlocking(int K, double rho) {
        if (rho <= 1e-15) return 0.0;
        if (Math.abs(rho - 1.0) < 1e-9) return 1.0 / (K + 1.0);
        return (1.0 - rho) * Math.pow(rho, K) / (1.0 - Math.pow(rho, K + 1));
    }
}
