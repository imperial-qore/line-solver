/**
 * @file Method of Layers on the SRVN decomposition of an entry-only LQN.
 *
 * Port of {@code matlab/src/api/lqn/lqn_mol.m}.
 *
 * @since LINE 3.0
 */
package jline.api.lqn;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import jline.api.pfqn.ld.Pfqn_qdamva;
import jline.lang.constant.CallType;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.LayeredNetworkStruct;
import jline.util.matrix.Matrix;

/**
 * Method of Layers on the SRVN decomposition of a layered queueing network whose entries carry no
 * activity graph.
 *
 * <p>A compact, self-contained reimplementation of the layered fixed point {@code SolverLN} runs,
 * restricted to LQNs in which every entry binds exactly one activity and there are no activity
 * precedences. It decomposes the model the way {@code lqns --srvn-layering} does -- one submodel
 * per processor and one per called task -- and sweeps them in the two phases of Rolia-Sevcik's
 * Method of Layers: all software (task) submodels, then all hardware (processor) ones.
 *
 * <p>Every submodel is a closed multiclass queueing network with ONE station and one class per
 * client task, so it is solved by {@link Pfqn_qdamva} rather than by building a Network. The
 * surrogate client delay of {@code SolverLN} collapses into the think-time vector Z of that call.
 *
 * <p>WHERE THIS DIFFERS FROM {@code SolverLN}'s {@code srvn.cs}: a submodel here carries one class
 * per client TASK with visit-weighted demands, where {@code srvn.cs} carries one class per activity
 * and encodes the call multiplicities as routing. On an entry-only model the two agree on the
 * structure and differ only in the aggregation, so the throughputs and processor utilizations track
 * closely while entry response times spread more.
 *
 * <p>SCOPE. Entry-only models. Activity graphs (fork/join, OR-branches, loops, second phases,
 * forwarding), asynchronous calls, caches, setup tasks, admission constraints, replication and open
 * arrivals are REFUSED, not approximated, and named when they are.
 *
 * <p>INDEX BASE. {@link LayeredNetworkStruct} is 0-BASED, unlike MATLAB's {@code lsn}, which
 * numbers from 1: hosts run {@code 0..nhosts-1}, tasks {@code tshift..tshift+ntasks-1}, and so on,
 * with {@code nidx} the element COUNT. Every loop here is written in that convention.
 */
public final class Lqn_mol {
    private Lqn_mol() {}

    /** GlobalConstants.FineTol, the reference's own "effectively zero". */
    private static final double FINE_TOL = 1e-8;

    /** Tuning of the outer fixed point. */
    public static class Options {
        public int iter_max = 200;
        public double iter_tol = 1e-6;
        public double relax_factor = 0.5;
    }

    /**
     * The four (nidx) vectors in the column convention {@code SolverLN} and LQNS report, so they
     * line up with {@code LN(model).getAvgTable()} cell for cell.
     *
     * <table>
     *   <caption>Reported measures by element kind</caption>
     *   <tr><th>index<th>QN (QLen)<th>UN (Util)<th>RN (RespT)<th>TN (Tput)
     *   <tr><td>host<td>NaN<td>processor utilization<td>NaN<td>NaN
     *   <tr><td>task<td>sum of entry T*S<td>sum of entry proc util<td>NaN<td>cycle rate
     *   <tr><td>entry<td>T*S<td>processor utilization<td>response time<td>throughput
     *   <tr><td>activity<td>as its entry<td>as its entry<td>as its entry<td>as its entry
     * </table>
     */
    public static class Result {
        public final double[] QN;
        public final double[] UN;
        public final double[] RN;
        public final double[] TN;
        /** Iterations performed. */
        public final int iter;
        /** The infinity norm of the last relative change, over servt and thinkt alike. */
        public final double resid;
        /** Entry response time seen by a caller. */
        public final double[] servt;
        /** Entry processor residence. */
        public final double[] residt;
        /** Blocking time per call. */
        public final double[] callservt;
        /** Task surrogate idle time. */
        public final double[] thinkt;
        /** Entry share of its task's invocations. */
        public final double[] share;

        Result(double[] QN, double[] UN, double[] RN, double[] TN, int iter, double resid,
               double[] servt, double[] residt, double[] callservt, double[] thinkt,
               double[] share) {
            this.QN = QN;
            this.UN = UN;
            this.RN = RN;
            this.TN = TN;
            this.iter = iter;
            this.resid = resid;
            this.servt = servt;
            this.residt = residt;
            this.callservt = callservt;
            this.thinkt = thinkt;
            this.share = share;
        }
    }

    /** Method of Layers with the reference's default tuning. */
    public static Result lqn_mol(LayeredNetworkStruct lsn) {
        return lqn_mol(lsn, new Options());
    }

    /**
     * Method of Layers on an entry-only, closed, synchronous LQN.
     *
     * @param lsn     the layered struct, from {@code LayeredNetwork.getStruct()}
     * @param options iteration cap, tolerance and under-relaxation factor
     * @return the four measure vectors and the fixed-point state behind them
     */
    public static Result lqn_mol(LayeredNetworkStruct lsn, Options options) {
        if (options == null) {
            options = new Options();
        }
        assertSupported(lsn);

        final int nidx = lsn.nidx;
        final int ncalls = lsn.ncalls;
        final int e0 = lsn.eshift, e1 = lsn.eshift + lsn.nentries;
        final int t0 = lsn.tshift, t1 = lsn.tshift + lsn.ntasks;
        final double om = options.relax_factor;

        // ---- static per-entry data: bound activity, host demand, owning task ----
        int[] actof = new int[nidx];
        int[] taskof = new int[nidx];
        int[] hostof = new int[nidx];
        double[] dem = new double[nidx];
        for (int tidx = t0; tidx < t1; tidx++) {
            hostof[tidx] = (int) lsn.parent.get(tidx);
        }
        for (int eidx = e0; eidx < e1; eidx++) {
            actof[eidx] = lsn.actsof.get(eidx).get(0);
            Double d = lsn.hostdem_mean.get(actof[eidx]);
            dem[eidx] = (d == null || Double.isNaN(d)) ? 0.0 : d;
            taskof[eidx] = (int) lsn.parent.get(eidx);
        }

        // ---- static per-call data ----------------------------------------------
        int[] callsrc = new int[ncalls];
        int[] calldst = new int[ncalls];
        double[] cally = new double[ncalls];
        int[] entryOfAct = new int[nidx];
        for (int eidx = e0; eidx < e1; eidx++) {
            entryOfAct[actof[eidx]] = eidx;
        }
        List<List<Integer>> callsFrom = new ArrayList<>();
        List<List<Integer>> callsTo = new ArrayList<>();
        for (int i = 0; i < nidx; i++) {
            callsFrom.add(new ArrayList<Integer>());
            callsTo.add(new ArrayList<Integer>());
        }
        for (int c = 0; c < ncalls; c++) {
            callsrc[c] = entryOfAct[(int) lsn.callpair.get(c, 0)];
            calldst[c] = (int) lsn.callpair.get(c, 1);
            Double y = lsn.callproc_mean.get(c);
            cally[c] = (y == null || Double.isNaN(y)) ? 0.0 : y;
            callsFrom.get(callsrc[c]).add(c);
            callsTo.get(calldst[c]).add(c);
        }

        // ---- populations, from maxmult (mult is wrong for INF tasks) ------------
        double[] npop = new double[nidx];
        Arrays.fill(npop, 1.0);
        for (int idx = 0; idx < t1; idx++) {
            double m = lsn.maxmult.get(idx);
            if (!Double.isFinite(m) || m < 1.0) {
                m = 1.0;
            }
            npop[idx] = m;
        }

        // ---- layer sets --------------------------------------------------------
        // One hardware layer per populated host, one software layer per called
        // non-reference task, as buildLayers.m draws them.
        List<Integer> hostLayers = new ArrayList<>();
        for (int h = 0; h < lsn.nhosts; h++) {
            List<Integer> ts = lsn.tasksof.get(h);
            if (ts != null && !ts.isEmpty()) {
                hostLayers.add(h);
            }
        }
        boolean[] isCalled = new boolean[nidx];
        for (int c = 0; c < ncalls; c++) {
            isCalled[taskof[calldst[c]]] = true;
        }
        List<Integer> taskLayers = new ArrayList<>();
        for (int tidx = t0; tidx < t1; tidx++) {
            if (lsn.isref.get(tidx) == 0 && isCalled[tidx]) {
                taskLayers.add(tidx);
            }
        }

        // ---- fixed-point state -------------------------------------------------
        double[] residt = dem.clone();
        double[] servt = new double[nidx];
        double[] callservt = new double[ncalls];
        double[] thinkt = new double[nidx];
        double[] share = new double[nidx];
        double[] Xtask = new double[nidx];
        double[] Xentry = new double[nidx];
        double[] busyth = new double[nidx];
        double[] zref = new double[nidx];
        for (int tidx = t0; tidx < t1; tidx++) {
            // lqn_ref_thinktime: a reference task's declared think time, and zero
            // for every other task and for a negative or non-finite one.
            double z = 0.0;
            if (lsn.isref.get(tidx) != 0) {
                Double zz = lsn.think_mean.get(tidx);
                z = (zz == null) ? 0.0 : zz;
                if (!Double.isFinite(z) || z < 0.0) {
                    z = 0.0;
                }
            }
            zref[tidx] = z;
            thinkt[tidx] = z;
            List<Integer> es = entriesOf(lsn, tidx);
            for (int eidx : es) {
                share[eidx] = 1.0 / es.size();
            }
        }
        // Seed servt bottom-up over the call graph so a callee is priced before
        // its caller; a cycle just leaves the residual demand seeded at 0.
        for (int eidx = e0; eidx < e1; eidx++) {
            servt[eidx] = dem[eidx];
        }
        for (int pass = 0; pass < Math.max(1, lsn.nentries); pass++) {
            for (int eidx = e0; eidx < e1; eidx++) {
                double s = residt[eidx];
                for (int c : callsFrom.get(eidx)) {
                    s += cally[c] * servt[calldst[c]];
                }
                servt[eidx] = s;
            }
        }
        for (int c = 0; c < ncalls; c++) {
            callservt[c] = servt[calldst[c]];
        }

        int iter = 0;
        double resid = Double.POSITIVE_INFINITY;
        while (iter < options.iter_max) {
            iter++;
            double[] servt_prev = servt.clone();
            double[] thinkt_prev = thinkt.clone();

            // ---- phase 1: software layers (thread contention at each called task)
            for (int tidx : taskLayers) {
                // The task is the station, its caller tasks the classes.
                Set<Integer> callerset = new TreeSet<>();
                for (int eidx : entriesOf(lsn, tidx)) {
                    for (int c : callsTo.get(eidx)) {
                        callerset.add(taskof[callsrc[c]]);
                    }
                }
                List<Integer> callers = new ArrayList<>(callerset);
                final int K = callers.size();
                if (K == 0) {
                    continue;
                }
                double[] gcl = new double[K];
                Arrays.fill(gcl, 1.0);
                if (lsn.sched.get(tidx) != SchedStrategy.INF) {
                    // An infinite-thread task never queues for a thread.
                    Matrix L = new Matrix(1, K);
                    Matrix N = new Matrix(1, K);
                    Matrix Z = new Matrix(1, K);
                    for (int k = 0; k < K; k++) {
                        int ctask = callers.get(k);
                        N.set(0, k, npop[ctask]);
                        double d = 0.0;
                        for (int eidx : entriesOf(lsn, ctask)) {
                            for (int c : callsFrom.get(eidx)) {
                                if (taskof[calldst[c]] == tidx) {
                                    d += share[eidx] * cally[c] * servt[calldst[c]];
                                }
                            }
                        }
                        L.set(0, k, d);
                        Z.set(0, k, cycleOutside(lsn, ctask, tidx, thinkt, share, residt, callservt,
                                cally, taskof, calldst, callsFrom));
                    }
                    // The AMVA U output is X*L*g, where g is the RECIPROCAL RATE
                    // MULTIPLIER at the current congestion, not 1/c -- it is not
                    // a busy-server count, so nothing here reads it. Occupancy
                    // comes from Little's law in throughputs().
                    Pfqn_qdamva.Result r =
                            Pfqn_qdamva.pfqn_qdamva(L, N, Z, molMu(N, npop[tidx]));
                    for (int k = 0; k < K; k++) {
                        if (L.get(0, k) > FINE_TOL) {
                            gcl[k] = r.R.get(0, k) / L.get(0, k);
                        }
                    }
                }
                for (int k = 0; k < K; k++) {
                    int ctask = callers.get(k);
                    for (int eidx : entriesOf(lsn, ctask)) {
                        for (int c : callsFrom.get(eidx)) {
                            if (taskof[calldst[c]] != tidx) {
                                continue;
                            }
                            double newv = gcl[k] * servt[calldst[c]];
                            callservt[c] = om * newv + (1.0 - om) * callservt[c];
                        }
                    }
                }
            }

            // ---- phase 2: hardware layers (processor contention at each host) ---
            for (int hidx : hostLayers) {
                // The processor is the station, its tasks the classes.
                List<Integer> tsks = lsn.tasksof.get(hidx);
                final int K = tsks.size();
                double[] f = new double[K];
                Arrays.fill(f, 1.0);
                if (lsn.sched.get(hidx) != SchedStrategy.INF) {
                    // A delay processor never queues.
                    Matrix L = new Matrix(1, K);
                    Matrix N = new Matrix(1, K);
                    Matrix Z = new Matrix(1, K);
                    for (int k = 0; k < K; k++) {
                        int tidx = tsks.get(k);
                        N.set(0, k, npop[tidx]);
                        double d = 0.0;
                        double z = thinkt[tidx];
                        for (int eidx : entriesOf(lsn, tidx)) {
                            d += share[eidx] * dem[eidx];
                            for (int c : callsFrom.get(eidx)) {
                                z += share[eidx] * cally[c] * callservt[c];
                            }
                        }
                        L.set(0, k, d);
                        Z.set(0, k, z);
                    }
                    Pfqn_qdamva.Result r =
                            Pfqn_qdamva.pfqn_qdamva(L, N, Z, molMu(N, npop[hidx]));
                    for (int k = 0; k < K; k++) {
                        if (L.get(0, k) > FINE_TOL) {
                            f[k] = r.R.get(0, k) / L.get(0, k);
                        }
                    }
                }
                for (int k = 0; k < K; k++) {
                    for (int eidx : entriesOf(lsn, tsks.get(k))) {
                        residt[eidx] = f[k] * dem[eidx];
                    }
                }
            }

            // ---- recompose entry service times ---------------------------------
            for (int eidx = e0; eidx < e1; eidx++) {
                double s = residt[eidx];
                for (int c : callsFrom.get(eidx)) {
                    s += cally[c] * callservt[c];
                }
                servt[eidx] = om * s + (1.0 - om) * servt[eidx];
            }

            // ---- throughputs, entry shares, think-time closure ------------------
            throughputs(lsn, e0, e1, t0, t1, servt, share, thinkt, npop, Xtask, Xentry, busyth,
                    taskof, callsrc, callsTo, calldst, cally);
            for (int tidx = t0; tidx < t1; tidx++) {
                if (lsn.isref.get(tidx) != 0) {
                    thinkt[tidx] = zref[tidx];
                    continue;
                }
                if (Xtask[tidx] <= FINE_TOL) {
                    continue;
                }
                // Idle time of a thread per cycle. updateThinkTimes splits this
                // into an INF arm (njobs - util) and a finite arm
                // (njobs*abs(1-util)) only because LINE reports busy SERVERS at
                // an infinite server and a busy FRACTION at a finite one;
                // carrying the count in both cases makes the two arms the same
                // expression.
                double newz = Math.max(0.0,
                        Math.abs(npop[tidx] - busyth[tidx]) / Xtask[tidx] - zref[tidx]);
                thinkt[tidx] = om * newz + (1.0 - om) * thinkt[tidx];
            }

            // Both halves of the state must settle: servt alone can sit still for
            // an iteration while the think times are still moving.
            resid = 0.0;
            for (int eidx = e0; eidx < e1; eidx++) {
                resid = Math.max(resid, Math.abs(servt[eidx] - servt_prev[eidx])
                        / Math.max(1.0, Math.abs(servt[eidx])));
            }
            for (int tidx = t0; tidx < t1; tidx++) {
                resid = Math.max(resid, Math.abs(thinkt[tidx] - thinkt_prev[tidx])
                        / Math.max(1.0, Math.abs(thinkt[tidx])));
            }
            if (resid < options.iter_tol) {
                break;
            }
        }
        throughputs(lsn, e0, e1, t0, t1, servt, share, thinkt, npop, Xtask, Xentry, busyth, taskof,
                callsrc, callsTo, calldst, cally);

        // ---- assemble the reported vectors -------------------------------------
        double[] QN = new double[nidx];
        double[] UN = new double[nidx];
        double[] RN = new double[nidx];
        double[] TN = new double[nidx];
        Arrays.fill(QN, Double.NaN);
        Arrays.fill(UN, Double.NaN);
        Arrays.fill(RN, Double.NaN);
        Arrays.fill(TN, Double.NaN);
        for (int eidx = e0; eidx < e1; eidx++) {
            int hidx = hostof[taskof[eidx]];
            double procutil = Xentry[eidx] * dem[eidx] / hostServers(lsn, hidx, npop);
            QN[eidx] = Xentry[eidx] * servt[eidx];
            UN[eidx] = procutil;
            RN[eidx] = servt[eidx];
            TN[eidx] = Xentry[eidx];
            int aidx = actof[eidx];
            QN[aidx] = QN[eidx];
            UN[aidx] = UN[eidx];
            RN[aidx] = RN[eidx];
            TN[aidx] = TN[eidx];
        }
        for (int tidx = t0; tidx < t1; tidx++) {
            double q = 0.0;
            double u = 0.0;
            for (int eidx : entriesOf(lsn, tidx)) {
                q += QN[eidx];
                u += UN[eidx];
            }
            QN[tidx] = q;
            UN[tidx] = u;
            RN[tidx] = Double.NaN;
            TN[tidx] = Xtask[tidx];
        }
        for (int hidx = 0; hidx < lsn.nhosts; hidx++) {
            double u = 0.0;
            List<Integer> ts = lsn.tasksof.get(hidx);
            if (ts != null) {
                for (int tidx : ts) {
                    u += UN[tidx];
                }
            }
            QN[hidx] = Double.NaN;
            UN[hidx] = u;
            RN[hidx] = Double.NaN;
            // No throughput is defined at a processor, as in LQNS.
            TN[hidx] = Double.NaN;
        }
        return new Result(QN, UN, RN, TN, iter, resid, servt, residt, callservt, thinkt, share);
    }

    /** The entries of a task, never null. */
    private static List<Integer> entriesOf(LayeredNetworkStruct lsn, int tidx) {
        List<Integer> es = lsn.entriesof.get(tidx);
        return (es == null) ? new ArrayList<Integer>() : es;
    }

    /**
     * LINE scales a station utilization into [0,1] whatever its multiplicity, and reports busy
     * SERVERS at an infinite server.
     */
    private static double hostServers(LayeredNetworkStruct lsn, int hidx, double[] npop) {
        return (lsn.sched.get(hidx) == SchedStrategy.INF) ? 1.0 : npop[hidx];
    }

    /**
     * Time a thread of {@code tidx} spends away from {@code excl} in one cycle: its think time, its
     * own processor residence, and its blocking at every callee other than {@code excl}.
     */
    private static double cycleOutside(LayeredNetworkStruct lsn, int tidx, int excl, double[] thinkt,
                                       double[] share, double[] residt, double[] callservt,
                                       double[] cally, int[] taskof, int[] calldst,
                                       List<List<Integer>> callsFrom) {
        double z = thinkt[tidx];
        for (int eidx : entriesOf(lsn, tidx)) {
            double w = share[eidx] * residt[eidx];
            for (int c : callsFrom.get(eidx)) {
                if (taskof[calldst[c]] != excl) {
                    w += share[eidx] * cally[c] * callservt[c];
                }
            }
            z += w;
        }
        return z;
    }

    /**
     * Reference tasks set the pace; every other rate follows from the call rates, so the entries
     * are visited in call-graph order until stable.
     */
    private static void throughputs(LayeredNetworkStruct lsn, int e0, int e1, int t0, int t1,
                                    double[] servt, double[] share, double[] thinkt, double[] npop,
                                    double[] Xtask, double[] Xentry, double[] busyth, int[] taskof,
                                    int[] callsrc, List<List<Integer>> callsTo, int[] calldst,
                                    double[] cally) {
        for (int t = t0; t < t1; t++) {
            if (lsn.isref.get(t) != 0) {
                double cyc = thinkt[t];
                for (int eidx : entriesOf(lsn, t)) {
                    cyc += share[eidx] * servt[eidx];
                }
                Xtask[t] = (cyc > FINE_TOL) ? npop[t] / cyc : 0.0;
                for (int eidx : entriesOf(lsn, t)) {
                    Xentry[eidx] = Xtask[t] * share[eidx];
                }
            } else {
                Xtask[t] = 0.0;
                for (int eidx : entriesOf(lsn, t)) {
                    Xentry[eidx] = 0.0;
                }
            }
        }
        for (int pass = 0; pass < Math.max(1, lsn.ntasks); pass++) {
            for (int eidx = e0; eidx < e1; eidx++) {
                if (lsn.isref.get(taskof[eidx]) != 0) {
                    continue;
                }
                double x = 0.0;
                for (int c : callsTo.get(eidx)) {
                    x += cally[c] * Xentry[callsrc[c]];
                }
                Xentry[eidx] = x;
            }
            for (int t = t0; t < t1; t++) {
                if (lsn.isref.get(t) != 0) {
                    continue;
                }
                List<Integer> es = entriesOf(lsn, t);
                double sum = 0.0;
                for (int eidx : es) {
                    sum += Xentry[eidx];
                }
                Xtask[t] = sum;
                if (sum > FINE_TOL) {
                    for (int eidx : es) {
                        share[eidx] = Xentry[eidx] / sum;
                    }
                }
            }
        }
        // Mean busy threads, by Little's law over the entries the task serves. This is the
        // occupancy the think-time closure needs, and it is exact given the throughputs -- unlike
        // the layer AMVA's own U.
        for (int t = t0; t < t1; t++) {
            double u = 0.0;
            for (int eidx : entriesOf(lsn, t)) {
                u += Xentry[eidx] * servt[eidx];
            }
            busyth[t] = u;
        }
    }

    /**
     * The queue-dependent rate multiplier row of a c-server station over a population of sum(N).
     *
     * <p>{@code Pfqn_lldfun} SKIPS a constant row, so a single server must come back as a row of
     * ones and not as a scalar 1, or the multiserver term is never applied.
     */
    private static Matrix molMu(Matrix N, double c) {
        double tot = 0.0;
        for (int k = 0; k < N.length(); k++) {
            tot += N.get(k);
        }
        int smax = Math.max(2, (int) Math.ceil(Math.max(0.0, tot)));
        Matrix mu = new Matrix(1, smax);
        for (int n = 1; n <= smax; n++) {
            mu.set(0, n - 1, (!Double.isFinite(c) || c <= 1.0) ? 1.0 : Math.min(n, c));
        }
        return mu;
    }

    /**
     * Refuse every feature this decomposition does not represent, NAMING the element, rather than
     * returning a number that quietly ignores it.
     */
    private static void assertSupported(LayeredNetworkStruct lsn) {
        for (int eidx = lsn.eshift; eidx < lsn.eshift + lsn.nentries; eidx++) {
            List<Integer> acts = lsn.actsof.get(eidx);
            int n = (acts == null) ? 0 : acts.size();
            if (n != 1) {
                throw new UnsupportedOperationException(String.format(
                        "lqn_mol: entry %s binds %d activities. lqn_mol solves entry-only models; "
                                + "use SolverLN for an activity graph.",
                        lsn.hashnames.get(eidx), n));
            }
        }
        for (int a = 0; a < lsn.nacts; a++) {
            int aidx = lsn.ashift + a;
            // A successor INSIDE the activity band is a precedence; an edge to an entry or a task
            // is the ordinary binding every entry-only model has.
            for (int b = 0; b < lsn.nacts; b++) {
                if (lsn.graph.get(aidx, lsn.ashift + b) != 0) {
                    throw new UnsupportedOperationException(String.format(
                            "lqn_mol: activity %s has an activity precedence. lqn_mol solves "
                                    + "entry-only models; use SolverLN for an activity graph.",
                            lsn.hashnames.get(aidx)));
                }
            }
            if (lsn.actphase != null && a < lsn.actphase.length()
                    && (int) lsn.actphase.get(a) != 1) {
                throw new UnsupportedOperationException(String.format(
                        "lqn_mol: activity %s is in phase %d. lqn_mol supports phase 1 only.",
                        lsn.hashnames.get(aidx), (int) lsn.actphase.get(a)));
            }
        }
        for (int c = 0; c < lsn.ncalls; c++) {
            CallType ct = lsn.calltype.get(c);
            if (ct != CallType.SYNC) {
                throw new UnsupportedOperationException(String.format(
                        "lqn_mol: call %s is %s. lqn_mol supports synchronous calls only.",
                        lsn.callhashnames.get(c), String.valueOf(ct)));
            }
        }
        for (int idx = 0; idx < lsn.tshift + lsn.ntasks; idx++) {
            if (lsn.iscache != null && idx < lsn.iscache.length() && lsn.iscache.get(idx) != 0) {
                throw new UnsupportedOperationException(String.format(
                        "lqn_mol: %s is a cache task, which lqn_mol does not model.",
                        lsn.hashnames.get(idx)));
            }
            if (lsn.hassetup != null && idx < lsn.hassetup.length()
                    && lsn.hassetup.get(idx) != 0) {
                throw new UnsupportedOperationException(String.format(
                        "lqn_mol: %s has a setup time, which lqn_mol does not model.",
                        lsn.hashnames.get(idx)));
            }
            if (lsn.repl != null && idx < lsn.repl.length() && lsn.repl.get(idx) != 1.0) {
                throw new UnsupportedOperationException(String.format(
                        "lqn_mol: %s is replicated %d times, which lqn_mol does not model.",
                        lsn.hashnames.get(idx), (int) lsn.repl.get(idx)));
            }
            if (lsn.lincon != null && lsn.lincon.get(idx) != null) {
                throw new UnsupportedOperationException(String.format(
                        "lqn_mol: %s carries an admission constraint, which lqn_mol does not "
                                + "model.", lsn.hashnames.get(idx)));
            }
            SchedStrategy s = lsn.sched.get(idx);
            boolean ok = (s == SchedStrategy.PS || s == SchedStrategy.FCFS
                    || s == SchedStrategy.INF
                    || (idx >= lsn.nhosts && s == SchedStrategy.REF));
            if (!ok) {
                throw new UnsupportedOperationException(String.format(
                        "lqn_mol: %s is scheduled %s, which lqn_mol does not model.",
                        lsn.hashnames.get(idx), String.valueOf(s)));
            }
        }
        if (lsn.callgroups != null && !lsn.callgroups.isEmpty()) {
            throw new UnsupportedOperationException(
                    "lqn_mol: this model uses routed call groups, which lqn_mol does not model.");
        }
        if (lsn.arrival != null) {
            for (int eidx = lsn.eshift; eidx < lsn.eshift + lsn.nentries; eidx++) {
                if (lsn.arrival.get(eidx) != null) {
                    throw new UnsupportedOperationException(String.format(
                            "lqn_mol: entry %s has an open arrival. lqn_mol solves closed models "
                                    + "only.", lsn.hashnames.get(eidx)));
                }
            }
        }
    }
}
