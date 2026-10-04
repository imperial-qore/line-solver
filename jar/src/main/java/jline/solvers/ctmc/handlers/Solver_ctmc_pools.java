package jline.solvers.ctmc.handlers;

import jline.io.Ret;
import jline.lang.Event;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.Sync;
import jline.lang.constant.EventType;
import jline.lang.constant.HeteroSchedPolicy;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodeparam.ServiceNodeParam;
import jline.lang.nodes.StatefulNode;
import jline.lang.nodes.Station;
import jline.lang.state.AfterEventStationPool;
import jline.lang.state.CtmcPool;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.GlobalConstants;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Heterogeneous server pools (Queue.addServerType) for the CTMC state space; port of MATLAB
 * solver_ctmc_pools.m, following the LDES semantics.
 *
 * <p>A class-r job in service at a pooled station occupies one server of ONE compatible pool
 * t and is served by that pool's law (setHeteroService), or by the station's own law for r
 * when the pool declares none. The per-class service process is replaced by a block-diagonal
 * phase-type law with one block per compatible pool, in ascending pool order, and the pool
 * bookkeeping read by {@link AfterEventStationPool} is stored in {@code sn.ctmcpool}. The
 * pools are the server bank, so nservers becomes their total. PHASE synchronizations are
 * added where a pooled law has more than one phase, and sn.state is rebuilt in the pooled
 * layout.</p>
 *
 * <p>The rewrite works on a shallow copy whose mutated containers are replaced by copies, so
 * the caller's struct (typically the model's cached one) is left untouched.</p>
 */
public final class Solver_ctmc_pools {

    private Solver_ctmc_pools() {
    }

    /** The struct with every heterogeneous-server station rewritten, or sn itself when there is none. */
    public static NetworkStruct apply(NetworkStruct sn) {
        if (sn.nodeparam == null) {
            return sn;
        }
        int R = sn.nclasses;
        List<Integer> todo = new ArrayList<Integer>();
        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.isstation.get(ind) != 1) continue;
            int ist = (int) sn.nodeToStation.get(ind);
            ServiceNodeParam p = sn.getServiceParam(sn.stations.get(ist));
            if (p == null || p.nservertypes <= 0) continue;
            if (sn.ctmcpool != null && sn.ctmcpool.containsKey(ind)) continue;
            SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
            // PAS/OI stations model compatible servers through the OI rank rate, not through pools
            if (sched == SchedStrategy.PAS || sched == SchedStrategy.OI) continue;
            todo.add(ind);
        }
        if (todo.isEmpty()) {
            return sn;
        }
        NetworkStruct old = sn;
        sn = old.shallowCopy();
        sn.proc = copy2(old.proc);
        sn.pie = copy2(old.pie);
        sn.mu = copy2(old.mu);
        sn.phi = copy2(old.phi);
        sn.procid = copy2(old.procid);
        sn.isph = old.isph == null ? null : copy2(old.isph);
        sn.phases = old.phases.copy();
        sn.phasessz = old.phasessz.copy();
        sn.phaseshift = old.phaseshift.copy();
        sn.nservers = old.nservers.copy();
        sn.nvars = old.nvars.copy();
        sn.sync = new HashMap<Integer, Sync>(old.sync);
        sn.state = old.state == null ? null : new HashMap<StatefulNode, Matrix>(old.state);
        sn.ctmcpool = old.ctmcpool == null ? new HashMap<Integer, CtmcPool>()
                : new HashMap<Integer, CtmcPool>(old.ctmcpool);
        int local = sn.nnodes;
        for (int ind : todo) {
            int ist = (int) sn.nodeToStation.get(ind);
            poolStation(sn, ind, ist, R);
            for (int r = 0; r < R; r++) {
                if (sn.phases.get(ist, r) > 1 && old.phases.get(ist, r) <= 1 && !hasPhaseSync(sn, ind, r)) {
                    // refreshSync only emits a PHASE action where the station's own law has phases
                    Sync s = new Sync();
                    s.active.put(0, new Event(EventType.PHASE, ind, r, Double.NaN, new Matrix(0, 0), Double.NaN, Double.NaN));
                    s.passive.put(0, new Event(EventType.LOCAL, local, r, 1.0, new Matrix(0, 0), Double.NaN, Double.NaN));
                    sn.sync.put(sn.sync.size(), s);
                }
            }
            int isf = (int) sn.nodeToStateful.get(ind);
            StatefulNode node = sn.stateful.get(isf);
            if (sn.state != null && sn.state.get(node) != null && sn.state.get(node).getNumRows() > 0) {
                sn.state.put(node, rebuildState(old, sn, ind, ist, Matrix.extractRows(old.state.get(node), 0, 1, null)));
            }
        }
        return sn;
    }

    private static <A, B, C> Map<A, Map<B, C>> copy2(Map<A, Map<B, C>> m) {
        if (m == null) return null;
        Map<A, Map<B, C>> out = new HashMap<A, Map<B, C>>();
        for (Map.Entry<A, Map<B, C>> e : m.entrySet()) {
            out.put(e.getKey(), e.getValue() == null ? null : new HashMap<B, C>(e.getValue()));
        }
        return out;
    }

    private static String cname(NetworkStruct sn, int r) {
        return sn.classnames != null && r < sn.classnames.size() ? sn.classnames.get(r) : String.valueOf(r);
    }

    /** (D0, D1) of a served law, or null when the class is disabled or absent. */
    private static MatrixCell law(MatrixCell pc) {
        if (pc == null || pc.size() < 2 || pc.get(0) == null || pc.get(0).isEmpty() || pc.get(0).hasNaN()) {
            return null;
        }
        return pc;
    }

    private static double[] entry(Matrix D1) {
        int n = D1.getNumRows();
        double[] a = new double[n];
        double s = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                a[j] += D1.get(i, j);
            }
        }
        for (double v : a) s += v;
        for (int j = 0; j < n; j++) a[j] /= Math.max(s, GlobalConstants.FineTol);
        return a;
    }

    private static double[] exitRates(Matrix D0) {
        int n = D0.getNumRows();
        double[] ex = new double[n];
        for (int i = 0; i < n; i++) {
            double s = 0;
            for (int j = 0; j < n; j++) s += D0.get(i, j);
            ex[i] = -s;
        }
        return ex;
    }

    /** Mean alpha (-D0)^-1 1 of a phase-type law, by Gaussian elimination. */
    private static double phMean(Matrix D0, double[] a) {
        int n = D0.getNumRows();
        double[][] A = new double[n][n + 1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) A[i][j] = -D0.get(i, j);
            A[i][n] = 1.0;
        }
        for (int c = 0; c < n; c++) {
            int piv = c;
            for (int r = c + 1; r < n; r++) if (Math.abs(A[r][c]) > Math.abs(A[piv][c])) piv = r;
            double[] tmp = A[c];
            A[c] = A[piv];
            A[piv] = tmp;
            for (int r = 0; r < n; r++) {
                if (r == c) continue;
                double f = A[r][c] / A[c][c];
                for (int k = c; k <= n; k++) A[r][k] -= f * A[c][k];
            }
        }
        double m = 0;
        for (int i = 0; i < n; i++) m += a[i] * A[i][n] / A[i][i];
        return m;
    }

    private static void requirePh(MatrixCell law, ProcessType procid, NetworkStruct sn, int ind, int r, String what) {
        if (procid == ProcessType.MAP || procid == ProcessType.MMPP2 || procid == ProcessType.MMAP
                || procid == ProcessType.ME || procid == ProcessType.RAP) {
            throw new RuntimeException("SolverCTMC serves a heterogeneous server pool with a renewal phase-type law; class '"
                    + cname(sn, r) + "' at station '" + sn.nodenames.get(ind) + "' has a " + procid + " law in " + what + ".");
        }
        Matrix D0 = law.get(0);
        Matrix D1 = law.get(1);
        int n = D0.getNumRows();
        double[] ex = exitRates(D0);
        double[] a = entry(D1);
        boolean bad = false;
        double diff = 0;
        double norm = 0;
        for (int i = 0; i < n; i++) {
            if (ex[i] < -GlobalConstants.FineTol) bad = true;
            for (int j = 0; j < n; j++) {
                if (i != j && D0.get(i, j) < -GlobalConstants.FineTol) bad = true;
                diff += Math.abs(D1.get(i, j) - ex[i] * a[j]);
                norm += Math.abs(D1.get(i, j));
            }
        }
        if (bad || diff > GlobalConstants.CoarseTol * Math.max(1.0, norm)) {
            throw new RuntimeException("SolverCTMC serves a heterogeneous server pool with a renewal phase-type law; class '"
                    + cname(sn, r) + "' at station '" + sn.nodenames.get(ind) + "' has a correlated or matrix-exponential law in "
                    + what + ".");
        }
    }

    private static void refuseFeatures(NetworkStruct sn, int ind, int ist, int R, ServiceNodeParam p) {
        Station st = sn.stations.get(ist);
        String why = "";
        if (sn.lldscaling != null && !sn.lldscaling.isEmpty() && sn.lldscaling.getNumRows() > ist
                && anyNotOne(sn.lldscaling, ist)) {
            why = "load-dependent service";
        } else if (sn.cdscaling != null && sn.cdscaling.get(st) != null) {
            why = "class-dependent service";
        } else if (sn.jdscaling != null && sn.jdscaling.get(st) != null) {
            why = "joint-dependent service";
        } else if (sn.gdscaling != null) {
            why = "global dependence";
        } else if (sn.hasbreakdown != null && sn.hasbreakdown.length() > ind && sn.hasbreakdown.get(ind) == 1) {
            why = "server breakdowns";
        } else if (anyValue(sn.retrialProc, st)) {
            why = "retrial";
        } else if (anyValue(sn.balkingStrategy, st)) {
            why = "balking";
        } else if (anyValue(sn.impatienceClass, st)) {
            why = "reneging";
        } else if (sn.isbasblocking != null && sn.isbasblocking.length() > ind && sn.isbasblocking.get(ind) == 1) {
            why = "BAS blocking";
        } else if (sn.immfeed != null && !sn.immfeed.isEmpty() && sn.immfeed.getNumRows() > ist && rowPositive(sn.immfeed, ist)) {
            why = "immediate feedback";
        } else if (sn.replyblock != null && !sn.replyblock.isEmpty() && sn.replyblock.getNumRows() > ind
                && rowPositive(sn.replyblock, ind)) {
            why = "synchronous calls";
        } else if (p.serverparallelism != null && p.serverparallelism.elementMax() > 1) {
            why = "server parallelism";
        } else if (sn.nvars != null && Matrix.extractRows(sn.nvars, ind, ind + 1, null).elementSum() > 0) {
            why = "a feature that keeps local variables at the station (round-robin routing, MAP service, polling)";
        } else if (sn.issignal != null && sn.rtnodes != null && !sn.rtnodes.isEmpty()) {
            for (int s = 0; s < R && why.isEmpty(); s++) {
                if (sn.issignal.get(s) <= 0) continue;
                for (int row = 0; row < sn.rtnodes.getNumRows(); row++) {
                    if (sn.rtnodes.get(row, ind * R + s) > 0) {
                        why = "signals";
                        break;
                    }
                }
            }
        }
        if (!why.isEmpty()) {
            throw new RuntimeException("SolverCTMC does not combine heterogeneous server pools with " + why + " (station '"
                    + sn.nodenames.get(ind) + "').");
        }
    }

    private static boolean anyNotOne(Matrix m, int row) {
        for (int c = 0; c < m.getNumCols(); c++) if (m.get(row, c) != 1) return true;
        return false;
    }

    private static boolean rowPositive(Matrix m, int row) {
        for (int c = 0; c < m.getNumCols(); c++) if (m.get(row, c) > 0) return true;
        return false;
    }

    private static <V> boolean anyValue(Map<Station, Map<JobClass, V>> m, Station st) {
        if (m == null || m.get(st) == null) return false;
        for (V v : m.get(st).values()) {
            if (v instanceof MatrixCell ? ((MatrixCell) v).size() > 0 : v != null) return true;
        }
        return false;
    }

    private static void poolStation(NetworkStruct sn, int ind, int ist, int R) {
        Station station = sn.stations.get(ist);
        ServiceNodeParam p = sn.getServiceParam(station);
        String name = sn.nodenames.get(ind);
        int T = p.nservertypes;
        SchedStrategy sched = sn.sched.get(station);
        if (sched != SchedStrategy.FCFS && sched != SchedStrategy.HOL && sched != SchedStrategy.FCFSPRIO
                && sched != SchedStrategy.LCFS && sched != SchedStrategy.LCFSPRIO && sched != SchedStrategy.SIRO) {
            throw new RuntimeException("SolverCTMC supports heterogeneous server pools under FCFS, HOL, FCFSPRIO, LCFS, LCFSPRIO and "
                    + "SIRO; station '" + name + "' uses " + sched + ".");
        }
        refuseFeatures(sn, ind, ist, R, p);
        CtmcPool pool = new CtmcPool();
        pool.ntypes = T;
        pool.count = new int[T];
        pool.compat = new boolean[T][R];
        for (int t = 0; t < T; t++) {
            pool.count[t] = (int) p.serverspertype.get(t);
            for (int r = 0; r < R; r++) pool.compat[t][r] = p.servercompat.get(t, r) > 0;
        }
        pool.policy = p.heteroschedpolicy == null ? HeteroSchedPolicy.ORDER : p.heteroschedpolicy;
        pool.pools = new int[R][];
        pool.off = new int[R][];
        pool.len = new int[R][];
        pool.alpha = new double[R][][];
        pool.exit = new double[R][][];
        pool.D0 = new Matrix[R][];
        pool.fsfrate = new double[T][R];
        for (int r = 0; r < R; r++) {
            JobClass jc = sn.jobclasses.get(r);
            MatrixCell base = law(sn.proc.get(station) == null ? null : sn.proc.get(station).get(jc));
            boolean hasPoolLaw = false;
            for (int t = 0; t < T; t++) {
                hasPoolLaw = hasPoolLaw || (pool.compat[t][r] && poolLaw(p, t, r) != null);
            }
            if (base == null && !hasPoolLaw) {
                pool.pools[r] = new int[0];
                pool.off[r] = new int[0];
                pool.len[r] = new int[0];
                pool.alpha[r] = new double[0][];
                pool.exit[r] = new double[0][];
                pool.D0[r] = new Matrix[0];
                continue; // the class is not served here
            }
            List<Integer> pl = new ArrayList<Integer>();
            for (int t = 0; t < T; t++) if (pool.compat[t][r]) pl.add(t);
            if (pl.isEmpty()) {
                throw new RuntimeException("Station '" + name + "' declares no server pool compatible with class '" + cname(sn, r)
                        + "', so a job of that class would wait forever.");
            }
            if (base != null) {
                ProcessType pid = sn.procid.get(station) == null ? null : sn.procid.get(station).get(jc);
                requirePh(base, pid, sn, ind, r, "its default service");
            }
            int nk = pl.size();
            pool.pools[r] = new int[nk];
            pool.off[r] = new int[nk];
            pool.len[r] = new int[nk];
            pool.alpha[r] = new double[nk][];
            pool.exit[r] = new double[nk][];
            pool.D0[r] = new Matrix[nk];
            int shift = 0;
            for (int k = 0; k < nk; k++) {
                int t = pl.get(k);
                MatrixCell lw = poolLaw(p, t, r);
                if (lw == null) {
                    if (base == null) {
                        throw new RuntimeException("Server pool '" + p.servertypenames.get(t) + "' of station '" + name
                                + "' accepts class '" + cname(sn, r) + "' but declares no law for it, and the station's own "
                                + "service for the class is disabled.");
                    }
                    lw = base;
                } else {
                    requirePh(lw, null, sn, ind, r, "pool '" + p.servertypenames.get(t) + "'");
                }
                Matrix D0 = lw.get(0);
                double[] a = entry(lw.get(1));
                pool.pools[r][k] = t;
                pool.off[r][k] = shift;
                pool.len[r][k] = D0.getNumRows();
                pool.alpha[r][k] = a;
                pool.exit[r][k] = exitRates(D0);
                pool.D0[r][k] = D0.copy();
                pool.fsfrate[t][r] = 1.0 / phMean(D0, a);
                shift += D0.getNumRows();
            }
            int n = shift;
            Matrix D0x = new Matrix(n, n);
            for (int k = 0; k < nk; k++) {
                Matrix B = pool.D0[r][k];
                int o = pool.off[r][k];
                for (int i = 0; i < B.getNumRows(); i++)
                    for (int j = 0; j < B.getNumCols(); j++)
                        if (B.get(i, j) != 0) D0x.set(o + i, o + j, B.get(i, j));
            }
            double[] exx = exitRates(D0x);
            Matrix piex = new Matrix(1, n);
            for (int j = 0; j < pool.len[r][0]; j++) piex.set(0, j, pool.alpha[r][0][j]);
            Matrix D1x = new Matrix(n, n);
            Matrix mu = new Matrix(n, 1);
            Matrix phi = new Matrix(n, 1);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) if (exx[i] * piex.get(0, j) != 0) D1x.set(i, j, exx[i] * piex.get(0, j));
                mu.set(i, 0, -D0x.get(i, i));
                phi.set(i, 0, exx[i] / (-D0x.get(i, i)));
            }
            sn.proc.get(station).put(jc, new MatrixCell(D0x, D1x));
            sn.pie.get(station).put(jc, piex);
            sn.mu.get(station).put(jc, mu);
            sn.phi.get(station).put(jc, phi);
            sn.procid.get(station).put(jc, ProcessType.PH);
            if (sn.isph != null && sn.isph.get(station) != null) sn.isph.get(station).put(jc, true);
            sn.phases.set(ist, r, n);
            sn.phasessz.set(ist, r, n);
        }
        double[] ncls = new double[T];
        for (int t = 0; t < T; t++) for (int r = 0; r < R; r++) if (pool.compat[t][r]) ncls[t] += 1;
        Integer[] ord = new Integer[T];
        for (int t = 0; t < T; t++) ord[t] = t;
        final double[] nc = ncls;
        Arrays.sort(ord, (x, y) -> Double.compare(nc[x], nc[y])); // stable
        pool.alfsorder = new int[T];
        for (int t = 0; t < T; t++) pool.alfsorder[t] = ord[t];
        boolean multi = false;
        for (int r = 0; r < R; r++) multi = multi || pool.pools[r].length >= 2;
        pool.rotate = (pool.policy == HeteroSchedPolicy.ALIS || pool.policy == HeteroSchedPolicy.FAIRNESS) && multi;
        pool.perms = pool.rotate ? permutations(T) : new int[0][];
        pool.varpos = 0;
        if (pool.rotate) {
            sn.nvars.set(ind, 2 * R, 1);
        }
        sn.ctmcpool.put(ind, pool);
        int tot = 0;
        for (int c : pool.count) tot += c;
        sn.nservers.set(ist, 0, tot); // the pools are the server bank, as in LDES
        int shift = 0;
        for (int r = 0; r < R; r++) {
            sn.phaseshift.set(ist, r, shift);
            shift += (int) sn.phasessz.get(ist, r);
        }
        if (sn.phaseshift.getNumCols() > R) sn.phaseshift.set(ist, R, shift);
    }

    private static MatrixCell poolLaw(ServiceNodeParam p, int t, int r) {
        if (p.heteroproc == null || p.heteroproc.get(t) == null) return null;
        MatrixCell c = p.heteroproc.get(t).get(r);
        return (c == null || c.size() < 2) ? null : c;
    }

    /** All permutations of 0..T-1 in lexicographic order (identity first). */
    private static int[][] permutations(int T) {
        List<int[]> out = new ArrayList<int[]>();
        permRec(new int[T], new boolean[T], 0, out);
        return out.toArray(new int[0][]);
    }

    private static void permRec(int[] cur, boolean[] used, int pos, List<int[]> out) {
        if (pos == cur.length) {
            out.add(cur.clone());
            return;
        }
        for (int v = 0; v < cur.length; v++) {
            if (used[v]) continue;
            used[v] = true;
            cur[pos] = v;
            permRec(cur, used, pos + 1, out);
            used[v] = false;
        }
    }

    private static boolean hasPhaseSync(NetworkStruct sn, int ind, int r) {
        for (Sync s : sn.sync.values()) {
            Event ev = s.active.get(0);
            if (ev != null && ev.getEvent() == EventType.PHASE && ev.getNode() == ind && ev.getJobClass() == r) return true;
        }
        return false;
    }

    /**
     * Replay the declared initial jobs as arrivals into an empty pooled station: in-service
     * jobs first, by class, then the waiting jobs from the oldest to the newest, each taking
     * the most likely outcome (RAIS: the first candidate pool).
     */
    private static Matrix rebuildState(NetworkStruct old, NetworkStruct sn, int ind, int ist, Matrix oldrow) {
        int R = sn.nclasses;
        int nsrvOld = 0;
        int[] Kold = new int[R];
        int[] Ksold = new int[R];
        for (int r = 0; r < R; r++) {
            Kold[r] = (int) old.phasessz.get(ist, r);
            Ksold[r] = (int) old.phaseshift.get(ist, r);
            nsrvOld += Kold[r];
        }
        int V = (int) Matrix.extractRows(old.nvars, ind, ind + 1, null).elementSum();
        int W = oldrow.getNumCols() - nsrvOld - V;
        List<Integer> seq = new ArrayList<Integer>();
        for (int r = 0; r < R; r++) {
            int cnt = 0;
            for (int j = 0; j < Kold[r]; j++) cnt += (int) oldrow.get(0, W + Ksold[r] + j);
            for (int c = 0; c < cnt; c++) seq.add(r);
        }
        boolean isSiro = sn.sched.get(sn.stations.get(ist)) == SchedStrategy.SIRO;
        if (isSiro) {
            for (int r = 0; r < R; r++) for (int c = 0; c < (int) oldrow.get(0, r); c++) seq.add(r);
        } else {
            for (int c = W - 1; c >= 0; c--) if (oldrow.get(0, c) > 0) seq.add((int) oldrow.get(0, c) - 1); // rightmost is the oldest
        }
        CtmcPool pool = sn.ctmcpool.get(ind);
        int Wb = isSiro ? R : Math.max(1, seq.size());
        int nsrv = 0;
        for (int r = 0; r < R; r++) nsrv += (int) sn.phasessz.get(ist, r);
        int Vn = pool.rotate ? 1 : 0;
        Matrix st = new Matrix(1, Wb + nsrv + Vn);
        if (pool.rotate) st.set(0, Wb + nsrv, 1);
        Matrix savedCap = sn.cap;
        Matrix savedClassCap = sn.classcap;
        sn.cap = savedCap.copy();
        sn.classcap = savedClassCap.copy();
        sn.cap.set(ist, 0, Double.POSITIVE_INFINITY); // the declared state is placed without the cutoff
        for (int r = 0; r < R; r++) sn.classcap.set(ist, r, Double.POSITIVE_INFINITY);
        try {
            for (int c : seq) {
                Ret.EventResult res = AfterEventStationPool.afterEventStationPool(sn, ind, ist, st, EventType.ARV, c, R, Vn, false);
                if (res.outspace.getNumRows() == 0) {
                    throw new RuntimeException("The initial state of station '" + sn.nodenames.get(ind)
                            + "' cannot be placed on its server pools.");
                }
                int best = 0;
                for (int i = 1; i < res.outprob.getNumRows(); i++) {
                    if (res.outprob.get(i, 0) > res.outprob.get(best, 0)) best = i;
                }
                st = Matrix.extractRows(res.outspace, best, best + 1, null);
            }
        } finally {
            sn.cap = savedCap;
            sn.classcap = savedClassCap;
        }
        return st;
    }

    /**
     * Pooled utilization over the whole server bank, as in LDES: UN(i,k) = E[sir_k]/nservers.
     * Overwrites row i of UN; the columns of station i in STATESPACE are c0..c1-1.
     */
    public static void utilization(NetworkStruct sn, int i, Matrix UN, Matrix stateSpace, int c0, int c1,
                                   Matrix wset, Matrix probSysState, int K) {
        int ind = (int) sn.stationToNode.get(i);
        for (int k = 0; k < K; k++) UN.set(i, k, 0);
        double S = sn.nservers.get(i);
        for (int index = 0; index < wset.length(); index++) {
            int st = (int) wset.get(index);
            State.StateMarginalStatistics tm = ToMarginal.toMarginal(sn, ind,
                    Matrix.extract(stateSpace, st, st + 1, c0, c1), null, null, null, null, null);
            for (int k = 0; k < K; k++) {
                UN.set(i, k, UN.get(i, k) + probSysState.get(st) * tm.sir.get(k) / S);
            }
        }
    }
}
