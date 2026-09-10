/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import java.util.ArrayList;
import java.util.List;

import jline.GlobalConstants;
import jline.api.mam.Ldqbd;
import jline.api.mam.LdqbdOptions;
import jline.api.mam.LdqbdResult;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

/**
 * Open classes of one station, solved as a level-dependent QBD MODULATED by the
 * background chain: level = number of open jobs held by the station, phase =
 * (arrival MAP phase, environment state, service phase).
 *
 * <p>The station is a MAP/PH/c queue whose server is shared with the closed
 * jobs. With k open and e closed jobs present the open aggregate completes at
 * rate phi(k,e) = min(k+e,c) * k/(k+e) times the phase-type completion rate of
 * one busy server: the open class receives the share k/(k+e) of the min(k+e,c)
 * busy servers. The dependence on k is what makes the QBD level-dependent, the
 * dependence on e is what makes it modulated. The c parallel servers are
 * collapsed into a single phase-type process scaled by phi, which is exact for
 * exponential service at any c and for phase-type service at c = 1, and
 * approximates the multiset of in-service phases otherwise -- the same collapse
 * {@link Solver_mam_ldqbd} documents.</p>
 *
 * <p>The environment is level-dependent too, and for the same reason. A lumped
 * transition that LOWERS the closed occupancy of this station is a closed
 * completion here, so it carries the closed share min(e+k,c)*e/(e+k) of the
 * server; the background chain built it at the averaged share gref, and level k
 * rescales it by the ratio of the two. A transition that RAISES the occupancy is
 * an arrival from elsewhere and is left alone. Without this the closed jobs
 * would drain at their mean-field rate however long the open queue is, and the
 * positive correlation between the two occupancies -- the very thing a congested
 * station produces -- would be lost.</p>
 *
 * <p>The level space is truncated at Kmax. An arrival at the top level is lost
 * but still advances the arrival phase, so the arrival process keeps its exact
 * marginal and autocorrelation and only the queue tail is cut.</p>
 *
 * @see Solver_mam_bgchain
 * @see Mam_bgchain_env
 */
public final class Mam_bgchain_station {
    private Mam_bgchain_station() {}

    /** Solved open station. */
    public static final class StationResult {
        /** Mean number of open jobs at the station. */
        public double QLen;
        /** Mean fraction of the servers held by open jobs. */
        public double Util;
        /** Open departure rate. */
        public double Tput;
        /** Stationary probability of the truncation level. */
        public double ploss;
        /** Stationary probability of each environment state. */
        public double[] penv;
        /** E[min(e+k,c)*e/(e+k) | e], the share the background chain reads back. */
        public double[] cshare;
        /** Environment support, ascending. */
        public int[] esup;
    }

    /**
     * @param Da0      D0 of the aggregate open arrival MAP
     * @param Da1      D1 of the aggregate open arrival MAP
     * @param alphaS   initial vector of the open service phase-type law
     * @param T        subgenerator of the open service phase-type law
     * @param A        environment generator from {@link Mam_bgchain_env}
     * @param esup     closed jobs each environment state stands for, ascending
     * @param nservers number of servers
     * @param gref     closed capacity share A was built at, one per environment state
     * @param Kmax     truncation level of the open queue
     * @param options  solver options
     * @return the solved station
     */
    public static StationResult mam_bgchain_station(Matrix Da0, Matrix Da1, Matrix alphaS, Matrix T,
                                                    Matrix A, int[] esup, double nservers, double[] gref,
                                                    int Kmax, SolverOptions options) {
        int ma = Da0.getNumRows();
        int me = esup.length;
        int ms = alphaS.getNumCols();

        Matrix t = T.mult(Matrix.ones(ms, 1)).scale(-1.0);
        Matrix Ima = Matrix.eye(ma);
        Matrix Ime = Matrix.eye(me);
        Matrix Ims = Matrix.eye(ms);

        if (Kmax < 1) Kmax = 1;

        // Split the environment into the closed departures from this station
        // (which the open level throttles) and the arrivals to it (which it does not).
        Matrix Adown = new Matrix(me, me);
        Matrix Aup = new Matrix(me, me);
        for (int e = 0; e < me; e++) {
            for (int ep = 0; ep < me; ep++) {
                if (ep < e) {
                    Adown.set(e, ep, A.get(e, ep));
                } else if (ep > e) {
                    Aup.set(e, ep, A.get(e, ep));
                }
            }
        }

        List<Matrix> Q0 = new ArrayList<Matrix>();
        List<Matrix> Q1 = new ArrayList<Matrix>();
        List<Matrix> Q2 = new ArrayList<Matrix>();
        List<double[]> phiae = new ArrayList<double[]>();

        Matrix A0 = envAtLevel(Aup, Adown, esup, gref, 0, nservers);
        Q1.add(Da0.kron(Ime).add(1.0, Ima.kron(A0)));
        Q0.add(Da1.kron(Ime).kron(alphaS));

        Matrix Da0kron = Da0.kron(Ime).kron(Ims);
        Matrix Da1kron = Da1.kron(Ime).kron(Ims);
        for (int k = 1; k <= Kmax; k++) {
            double[] share = openShare(k, esup, nservers);
            double[] rep = new double[ma * me];
            for (int a = 0; a < ma; a++) {
                System.arraycopy(share, 0, rep, a * me, me);
            }
            phiae.add(rep);
            Matrix Ak = envAtLevel(Aup, Adown, esup, gref, k, nservers);
            Matrix local = Da0kron.add(1.0, Ima.kron(Ak).kron(Ims)).add(1.0, diag(rep).kron(T));
            Q1.add(local);
            if (k < Kmax) {
                Q0.add(Da1kron);
            }
            if (k == 1) {
                Q2.add(diag(rep).kron(t));
            } else {
                Q2.add(diag(rep).kron(t.mult(alphaS)));
            }
        }
        // truncation: an arrival at the top level is lost, its phase transition is kept
        Q1.set(Kmax, Q1.get(Kmax).add(1.0, Da1kron));

        LdqbdOptions ldopts = new LdqbdOptions(options.tol, options.iter_max, false);
        LdqbdResult ld = Ldqbd.ldqbd(Q0, Q1, Q2, ldopts);
        Matrix plev = ld.getPi();
        List<Matrix> pcell = ld.getPiCells();

        double tot = 0;
        for (int k = 0; k <= Kmax; k++) {
            double v = Math.max(plev.get(0, k), 0.0);
            plev.set(0, k, v);
            tot += v;
        }
        if (tot > 0) plev.scaleEq(1.0 / tot);

        StationResult res = new StationResult();
        res.QLen = 0;
        for (int k = 0; k <= Kmax; k++) res.QLen += k * plev.get(0, k);
        res.ploss = plev.get(0, Kmax);
        res.esup = esup;

        double[] penv = new double[me];
        double[] gacc = new double[me];
        double util = 0;
        double tput = 0;
        for (int k = 0; k <= Kmax; k++) {
            Matrix pk = pcell.get(k);
            double[] marg = new double[me];
            if (k == 0) {
                for (int a = 0; a < ma; a++) {
                    for (int e = 0; e < me; e++) {
                        marg[e] += Math.max(pk.get(0, a * me + e), 0.0);
                    }
                }
            } else {
                double[] rep = phiae.get(k - 1);
                for (int a = 0; a < ma; a++) {
                    for (int e = 0; e < me; e++) {
                        double blockSum = 0;
                        for (int s = 0; s < ms; s++) {
                            double v = Math.max(pk.get(0, (a * me + e) * ms + s), 0.0);
                            blockSum += v;
                            tput += v * rep[a * me + e] * t.get(s, 0);
                        }
                        marg[e] += blockSum;
                        util += blockSum * rep[a * me + e];
                    }
                }
            }
            double[] closedShare = closedShare(esup, k, nservers);
            for (int e = 0; e < me; e++) {
                penv[e] += marg[e];
                gacc[e] += marg[e] * closedShare[e];
            }
        }
        double psum = 0;
        for (int e = 0; e < me; e++) psum += penv[e];
        if (psum > 0) {
            for (int e = 0; e < me; e++) {
                penv[e] /= psum;
                gacc[e] /= psum;
            }
        }
        res.cshare = new double[me];
        for (int e = 0; e < me; e++) {
            res.cshare[e] = penv[e] > GlobalConstants.Zero ? gacc[e] / penv[e] : 0.0;
        }
        res.penv = penv;
        res.Util = util / nservers;
        res.Tput = tput;
        return res;
    }

    /**
     * Share of the server capacity that k open jobs hold when e closed jobs are
     * also present: min(k+e,c) busy servers times the open fraction k/(k+e).
     */
    static double[] openShare(int k, int[] esup, double nservers) {
        double[] v = new double[esup.length];
        for (int e = 0; e < esup.length; e++) {
            double tot = k + esup[e];
            if (tot > 0) {
                v[e] = Math.min(tot, nservers) * (k / tot);
            }
        }
        return v;
    }

    /** The mirror image of {@link #openShare}, for the closed jobs. */
    static double[] closedShare(int[] esup, int k, double nservers) {
        double[] v = new double[esup.length];
        for (int e = 0; e < esup.length; e++) {
            double tot = esup[e] + k;
            if (tot > 0) {
                v[e] = Math.min(tot, nservers) * (esup[e] / tot);
            }
        }
        return v;
    }

    /**
     * Environment generator seen at open level k: the closed departures are
     * rescaled from the mean-field share gref to the share they hold against k
     * open jobs, the closed arrivals are unchanged, the diagonal is rebuilt.
     */
    static Matrix envAtLevel(Matrix Aup, Matrix Adown, int[] esup, double[] gref, int k, double nservers) {
        int me = esup.length;
        double[] g = closedShare(esup, k, nservers);
        Matrix Ak = new Matrix(me, me);
        for (int e = 0; e < me; e++) {
            double ratio = (gref[e] > 0) ? g[e] / gref[e] : 1.0;
            double diag = 0;
            for (int ep = 0; ep < me; ep++) {
                if (ep == e) continue;
                double v = Aup.get(e, ep) + ratio * Adown.get(e, ep);
                if (v != 0) {
                    Ak.set(e, ep, v);
                    diag += v;
                }
            }
            Ak.set(e, e, -diag);
        }
        return Ak;
    }

    private static Matrix diag(double[] v) {
        Matrix D = new Matrix(v.length, v.length);
        for (int i = 0; i < v.length; i++) {
            if (v[i] != 0) D.set(i, i, v[i]);
        }
        return D;
    }
}
