/**
 * Maximum Entropy algorithm for single-class open queueing networks with
 * finite buffers, loss and transfer blocking.
 *
 * Extends {@link Me_oqn} to open networks in which a station has a finite
 * buffer. Two per-station policies are supported. Under loss a job that
 * finds the destination full is discarded, so each station is a censored
 * GE/GE/c/0;N queue and the network is the ME decomposition of Kouvatsos
 * (1994), Section 4. Under transfer blocking a job that completes service
 * at station i and finds the destination j full is held in i's server,
 * which cannot serve anyone else until j has room (blocking after
 * service, BAS).
 *
 * Transfer blocking is not work conserving, so a product-form
 * approximation cannot be applied to the network as it stands. Following
 * Tahilramani, Manjunath and Bose (1999) the network is first made work
 * conserving by inserting a GE/GE/inf holding node on every routing pair
 * with a finite-buffer destination. The holding node absorbs the blocked
 * job, releasing station i's server; the delay it introduces is the
 * residual life of the minimum of the c_j service times in progress at j,
 * inflated geometrically because the released job may find j full again.
 * Station i's own service time is inflated by the same blocking
 * probability so that the jobs queued behind the blocked one still see
 * the server as busy. The expanded network is work conserving and is
 * solved node by node with the censored ME queue of {@link Me_gegecn},
 * iterating over the blocking probabilities and the first two moments of
 * the flows until they converge.
 *
 * @since LINE 3.0
 */
package jline.api.nc;

import jline.io.InputOutput;
import jline.util.matrix.Matrix;

public final class Me_oqn_blk {
    private Me_oqn_blk() {}

    /** Drop rule marker: a blocked job is discarded. */
    public static final int RULE_LOSS = 0;
    /** Drop rule marker: a blocked job is held in the upstream server. */
    public static final int RULE_BAS = 1;

    public static MeOqnBlkResult me_oqn_blk(int M, Matrix lambda0, Matrix Ca0, Matrix mu,
                                            Matrix Cs, Matrix P, Matrix c, Matrix N,
                                            int[] blockrule) {
        return me_oqn_blk(M, lambda0, Ca0, mu, Cs, P, c, N, blockrule, new MeOqnOptions(), 0.5);
    }

    /**
     * Maximum Entropy algorithm for a single-class open network with finite
     * buffers.
     *
     * @param M         number of stations
     * @param lambda0   external arrival rates [M x 1]
     * @param Ca0       external interarrival scvs [M x 1]
     * @param mu        service rates [M x 1]
     * @param Cs        service scvs [M x 1]
     * @param P         routing probabilities [M x M]; a row sum below one
     *                  sends the residual flow out of the network
     * @param c         servers per station [M x 1], infinite marking an IS
     * @param N         buffer capacity per station [M x 1] in jobs, in
     *                  service included, infinite marking an unbounded buffer
     * @param blockrule policy at each finite buffer, {@link #RULE_LOSS} or
     *                  {@link #RULE_BAS}
     * @param options   tolerance, iteration cap and verbosity
     * @param damping   relaxation weight applied to the blocking
     *                  probabilities
     * @return station metrics and the flow moments at the fixed point
     */
    public static MeOqnBlkResult me_oqn_blk(int M, Matrix lambda0, Matrix Ca0, Matrix mu,
                                            Matrix Cs, Matrix P, Matrix c, Matrix N,
                                            int[] blockrule, MeOqnOptions options, double damping) {
        boolean[] finiteBuf = new boolean[M];
        boolean[] bas = new boolean[M];
        for (int i = 0; i < M; i++) {
            finiteBuf[i] = !Double.isInfinite(N.get(i, 0)) && !Double.isInfinite(c.get(i, 0));
            bas[i] = finiteBuf[i] && blockrule != null && blockrule[i] == RULE_BAS;
            if (finiteBuf[i] && Cs.get(i, 0) < 1 - 1e-12) {
                InputOutput.line_error("me_oqn_blk", "MEM with finite buffers requires a service scv of at least 1 at station "
                        + (i + 1) + ": the GE distribution is not defined for scv < 1.");
            }
            if (lambda0.get(i, 0) > 0 && Ca0.get(i, 0) < 1 - 1e-12) {
                InputOutput.line_error("me_oqn_blk", "MEM with finite buffers requires an external interarrival scv of at least 1 at station "
                        + (i + 1) + ".");
            }
        }

        // see _kb/03-api-layer.md for rationale
        Matrix Pf = P.copy();
        double[] muf = new double[M];
        double[] Csf = new double[M];
        for (int i = 0; i < M; i++) {
            muf[i] = mu.get(i, 0);
            Csf[i] = Cs.get(i, 0);
            double pii = P.get(i, i);
            if (pii > 0) {
                muf[i] = mu.get(i, 0) * (1 - pii);
                Csf[i] = pii + (1 - pii) * Cs.get(i, 0);
                for (int j = 0; j < M; j++) {
                    Pf.set(i, j, P.get(i, j) / (1 - pii));
                }
                Pf.set(i, i, 0.0);
            }
        }

        // see _kb/03-api-layer.md for rationale
        double[] muRes = new double[M];
        for (int j = 0; j < M; j++) {
            double sigma = 2.0 / (Csf[j] + 1.0);
            muRes[j] = c.get(j, 0) * muf[j] * sigma;
        }

        double[] Ca = new double[M];
        double[] Cd = new double[M];
        double[] PBe = new double[M];
        double[] PBa = new double[M];
        double[] Qv = new double[M];
        double[] Uv = new double[M];
        double[] Tv = new double[M];
        double[] lam = new double[M];
        double[][] PBs = new double[M][M];
        double[][] PBh = new double[M][M];
        double[][] Lhold = new double[M][M];
        double[][] CaStreamInt = new double[M][M];
        for (int i = 0; i < M; i++) {
            Ca[i] = 1.0;
            Cd[i] = Csf[i];
        }

        int iter = 0;
        double delta = Double.POSITIVE_INFINITY;
        double[] attExt = new double[M];
        double[][] attInt = new double[M][M];
        for (iter = 1; iter <= options.getMaxIter(); iter++) {
            double[] Ca_old = Ca.clone();
            double[] PBe_old = PBe.clone();
            double[][] PBs_old = new double[M][M];
            for (int i = 0; i < M; i++) {
                PBs_old[i] = PBs[i].clone();
            }

            // see _kb/03-api-layer.md for rationale
            double[] PBf = new double[M];
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    if (Pf.get(i, j) > 0 && bas[j]) {
                        PBf[i] += Pf.get(i, j) * PBs[i][j];
                    }
                }
                if (PBf[i] >= 1 - 1e-9) {
                    InputOutput.line_error("me_oqn_blk",
                            "MEM transfer-blocking fixed point saturates: a station is blocked with probability one. The network has no stable operating point under BAS.");
                }
            }
            double[] muEff = new double[M];
            double[] CsEff = new double[M];
            for (int i = 0; i < M; i++) {
                muEff[i] = muf[i] * (1 - PBf[i]);
                CsEff[i] = PBf[i] + Csf[i] * (1 - PBf[i]);
            }

            // see _kb/03-api-layer.md for rationale
            Matrix A = new Matrix(M, M);
            Matrix b = new Matrix(M, 1);
            for (int j = 0; j < M; j++) {
                b.set(j, 0, lambda0.get(j, 0) * (1 - PBe[j]));
                for (int i = 0; i < M; i++) {
                    if (Pf.get(i, j) > 0) {
                        if (finiteBuf[j] && !bas[j]) {
                            A.set(i, j, Pf.get(i, j) * (1 - PBs[i][j]));
                        } else {
                            A.set(i, j, Pf.get(i, j));
                        }
                    }
                }
            }
            Matrix sys = new Matrix(M, M);
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    sys.set(i, j, (i == j ? 1.0 : 0.0) - A.get(j, i));
                }
            }
            Matrix Tm = new Matrix(M, 1);
            if (!Matrix.solveSafe(sys, b, Tm)) {
                InputOutput.line_error("me_oqn_blk", "MEM flow balance is singular: the routing matrix is reducible under blocking.");
            }
            for (int i = 0; i < M; i++) {
                Tv[i] = Math.max(0.0, Tm.get(i, 0));
            }

            // see _kb/03-api-layer.md for rationale
            for (int j = 0; j < M; j++) {
                attExt[j] = lambda0.get(j, 0);
                for (int i = 0; i < M; i++) {
                    attInt[i][j] = 0.0;
                    if (Pf.get(i, j) > 0) {
                        if (finiteBuf[j] && bas[j]) {
                            attInt[i][j] = Tv[i] * Pf.get(i, j) / Math.max(1 - PBs[i][j], 1e-12);
                        } else {
                            attInt[i][j] = Tv[i] * Pf.get(i, j);
                        }
                    }
                }
            }
            for (int j = 0; j < M; j++) {
                double tot = attExt[j];
                double wsum = attExt[j] * PBe[j];
                for (int i = 0; i < M; i++) {
                    tot += attInt[i][j];
                    wsum += attInt[i][j] * PBs[i][j];
                }
                lam[j] = tot;
                PBa[j] = tot > 0 ? wsum / tot : 0.0;
            }

            // see _kb/03-api-layer.md for rationale
            for (int j = 0; j < M; j++) {
                if (lam[j] <= 0) {
                    continue;
                }
                double sumInv = 0.0;
                if (attExt[j] > 0) {
                    sumInv += (attExt[j] / lam[j]) / (Ca0.get(j, 0) + 1);
                }
                for (int i = 0; i < M; i++) {
                    CaStreamInt[i][j] = 1.0;
                    if (attInt[i][j] > 0) {
                        CaStreamInt[i][j] = 1 - Pf.get(i, j) + Pf.get(i, j) * Cd[i];
                        sumInv += (attInt[i][j] / lam[j]) / (CaStreamInt[i][j] + 1);
                    }
                }
                if (sumInv > 0) {
                    Ca[j] = -1 + 1 / sumInv;
                }
            }

            // see _kb/03-api-layer.md for rationale
            double[] PBeNew = new double[M];
            double[][] PBsNew = new double[M][M];
            double[][] PBhNew = new double[M][M];
            for (int j = 0; j < M; j++) {
                if (Double.isInfinite(c.get(j, 0))) {
                    Qv[j] = muEff[j] > 0 ? lam[j] / muEff[j] : 0.0;
                    Uv[j] = Qv[j];
                    Cd[j] = Ca[j];
                    continue;
                }
                if (!finiteBuf[j]) {
                    double rho = muEff[j] > 0 ? lam[j] / (c.get(j, 0) * muEff[j]) : 0.0;
                    if (rho >= 1) {
                        Qv[j] = Double.POSITIVE_INFINITY;
                        Uv[j] = 1.0;
                        Cd[j] = CsEff[j];
                    } else if (c.get(j, 0) == 1) {
                        Qv[j] = rho * (Ca[j] + 1) / 2 + rho * rho * (CsEff[j] + Ca[j]) / (2 * (1 - rho));
                        Uv[j] = rho;
                        Cd[j] = rho * rho * CsEff[j] + (1 - rho) * Ca[j] + rho * (1 - rho);
                    } else {
                        Qv[j] = Me_oqn.geGecMql(lam[j], Ca[j], muEff[j], CsEff[j], (int) c.get(j, 0));
                        Uv[j] = rho;
                        Cd[j] = rho * rho * CsEff[j] + (1 - rho) * Ca[j] + rho * (1 - rho);
                    }
                    continue;
                }
                int cj = (int) c.get(j, 0);
                int Nj = (int) N.get(j, 0);
                MeGegecnResult nodeRes = Me_gegecn.me_gegecn(lam[j], Ca[j], muEff[j], CsEff[j], cj, 0, Nj);
                Qv[j] = nodeRes.getL();
                Uv[j] = nodeRes.getU();
                // see _kb/03-api-layer.md for rationale
                double uj = nodeRes.getU();
                Cd[j] = uj * uj * CsEff[j] + (1 - uj) * Ca[j] + uj * (1 - uj);
                if (attExt[j] > 0) {
                    PBeNew[j] = Me_gegecn.me_gegecn_pb(nodeRes.getP(), 0, Nj, cj, CsEff[j], Ca0.get(j, 0));
                }
                for (int i = 0; i < M; i++) {
                    if (attInt[i][j] > 0) {
                        PBsNew[i][j] = Me_gegecn.me_gegecn_pb(nodeRes.getP(), 0, Nj, cj, CsEff[j], CaStreamInt[i][j]);
                        if (bas[j]) {
                            // see _kb/03-api-layer.md for rationale
                            double q = Pf.get(i, j) * PBs[i][j];
                            double CaH = 1 - q + q * Cd[i];
                            PBhNew[i][j] = Me_gegecn.me_gegecn_pb(nodeRes.getP(), 0, Nj, cj, CsEff[j], CaH);
                        }
                    }
                }
            }

            // Relaxation on the blocking probabilities
            delta = 0.0;
            for (int j = 0; j < M; j++) {
                PBe[j] = (1 - damping) * PBe[j] + damping * PBeNew[j];
                for (int i = 0; i < M; i++) {
                    PBs[i][j] = (1 - damping) * PBs[i][j] + damping * PBsNew[i][j];
                    PBh[i][j] = (1 - damping) * PBh[i][j] + damping * PBhNew[i][j];
                    delta = Math.max(delta, Math.abs(PBs[i][j] - PBs_old[i][j]));
                }
                delta = Math.max(delta, Math.abs(Ca[j] - Ca_old[j]));
                delta = Math.max(delta, Math.abs(PBe[j] - PBe_old[j]));
            }
            if (options.getVerbose()) {
                System.out.printf("Iteration %d: max delta = %e%n", iter, delta);
            }
            if (delta < options.getTol()) {
                break;
            }
        }
        if (iter > options.getMaxIter() && delta >= options.getTol()) {
            InputOutput.line_warning("me_oqn_blk", "Did not converge within %d iterations (delta=%e)\n",
                    options.getMaxIter(), delta);
        }

        // see _kb/03-api-layer.md for rationale
        for (int i = 0; i < M; i++) {
            double add = 0.0;
            for (int j = 0; j < M; j++) {
                Lhold[i][j] = 0.0;
                if (bas[j] && Pf.get(i, j) > 0 && PBs[i][j] > 0) {
                    double rateH = Tv[i] * Pf.get(i, j) * PBs[i][j];
                    double muH = muRes[j] * (1 - PBh[i][j]);
                    if (muH > 0) {
                        Lhold[i][j] = rateH / muH;
                    }
                }
                add += Lhold[i][j];
            }
            Qv[i] += add;
        }

        // see _kb/03-api-layer.md for rationale
        Matrix Q = new Matrix(M, 1);
        Matrix W = new Matrix(M, 1);
        Matrix T = new Matrix(M, 1);
        Matrix U = new Matrix(M, 1);
        Matrix CaOut = new Matrix(M, 1);
        Matrix CdOut = new Matrix(M, 1);
        Matrix PBaOut = new Matrix(M, 1);
        Matrix lamOut = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            double u = 0.0;
            if (muf[i] > 0) {
                u = Double.isInfinite(c.get(i, 0)) ? Tv[i] / muf[i] : Tv[i] / (c.get(i, 0) * muf[i]);
            }
            Uv[i] = u;
            Q.set(i, 0, Qv[i]);
            T.set(i, 0, Tv[i]);
            U.set(i, 0, u);
            W.set(i, 0, Tv[i] > 0 ? Qv[i] / Tv[i] : 0.0);
            CaOut.set(i, 0, Ca[i]);
            CdOut.set(i, 0, Cd[i]);
            PBaOut.set(i, 0, PBa[i]);
            lamOut.set(i, 0, lam[i]);
        }
        return new MeOqnBlkResult(Q, W, T, U, CaOut, CdOut, PBaOut, lamOut, Math.min(iter, options.getMaxIter()));
    }
}
