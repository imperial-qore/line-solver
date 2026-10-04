package jline.lang.state;

import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.EventType;
import jline.lang.constant.HeteroSchedPolicy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.util.Maths;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.TreeSet;

/**
 * Event handler and state enumerator of a station with heterogeneous server pools, as
 * rewritten for the CTMC by {@code jline.solvers.ctmc.handlers.Solver_ctmc_pools}. Port of
 * MATLAB State.afterEventStationPool and State.fromMarginalPool.
 *
 * <p>The local state is [buffer | servers | local variables]. The buffer holds the WAITING
 * jobs only (1-based class ids right-aligned, the rightmost is the oldest; per-class counts
 * under SIRO), and the class-r server block holds one sub-block of phases per compatible
 * pool. The semantics are those of LDES:</p>
 * <ul>
 * <li>ARV: the candidates are the compatible pools with a free server. None: the job waits.
 * Otherwise ORDER takes the lowest index, ALFS the least flexible pool, FSF the pool with the
 * highest rate for the class (ties by index), RAIS each candidate with equal probability, and
 * ALIS/FAIRNESS the first candidate of the rotating pool order, which then moves to the back
 * when there was more than one candidate.</li>
 * <li>DEP: the freed server of pool t takes the first waiting job, in the discipline's
 * service order, that pool t can serve.</li>
 * <li>PHASE: a transition inside the block of the job's pool.</li>
 * </ul>
 * <p>Invariant: no waiting job has a free compatible server.</p>
 */
public final class AfterEventStationPool {

    private AfterEventStationPool() {
    }

    /** True when node ind is a pooled station of this (CTMC) struct. */
    public static boolean isPooled(NetworkStruct sn, int ind) {
        return sn.ctmcpool != null && sn.ctmcpool.containsKey(ind);
    }

    private static final class Out {
        final List<double[]> rows = new ArrayList<double[]>();
        final List<Double> rates = new ArrayList<Double>();
        final List<Double> probs = new ArrayList<Double>();
        final List<Integer> starts = new ArrayList<Integer>();

        void add(double[] row, double rate, double prob, int start) {
            rows.add(row);
            rates.add(rate);
            probs.add(prob);
            starts.add(start);
        }
    }

    private static double[] concat(double[] a, double[] b, double[] c) {
        double[] out = new double[a.length + b.length + c.length];
        System.arraycopy(a, 0, out, 0, a.length);
        System.arraycopy(b, 0, out, a.length, b.length);
        System.arraycopy(c, 0, out, a.length + b.length, c.length);
        return out;
    }

    private static int[] busy(CtmcPool pool, double[] srv, int[] Ks) {
        int[] b = new int[pool.ntypes];
        for (int r = 0; r < pool.pools.length; r++) {
            for (int k = 0; k < pool.pools[r].length; k++) {
                int c0 = Ks[r] + pool.off[r][k];
                for (int j = 0; j < pool.len[r][k]; j++) {
                    b[pool.pools[r][k]] += (int) srv[c0 + j];
                }
            }
        }
        return b;
    }

    /**
     * Successors of EVENT for class JOBCLASS at pooled station IND (0-based).
     */
    public static Ret.EventResult afterEventStationPool(NetworkStruct sn, int ind, int ist, Matrix inspace,
                                                        EventType event, int jobClass, int R, int V,
                                                        boolean isSimulation) {
        CtmcPool pool = sn.ctmcpool.get(ind);
        int[] K = new int[R];
        int[] Ks = new int[R];
        for (int r = 0; r < R; r++) {
            K[r] = (int) sn.phasessz.get(ist, r);
            Ks[r] = (int) sn.phaseshift.get(ist, r);
        }
        Station station = sn.stations.get(ist);
        SchedStrategy sched = sn.sched.get(station);
        boolean isSiro = sched == SchedStrategy.SIRO;
        int nsrv = 0;
        for (int r = 0; r < R; r++) {
            nsrv += K[r];
        }
        Out out = new Out();
        for (int row = 0; row < inspace.getNumRows(); row++) {
            int ncols = inspace.getNumCols();
            int W = ncols - nsrv - V;
            double[] buf = new double[W];
            double[] srv = new double[nsrv];
            double[] var = new double[V];
            for (int c = 0; c < W; c++) buf[c] = inspace.get(row, c);
            for (int c = 0; c < nsrv; c++) srv[c] = inspace.get(row, W + c);
            for (int c = 0; c < V; c++) var[c] = inspace.get(row, W + nsrv + c);
            if (event == EventType.ARV) {
                arrival(sn, ist, pool, buf, srv, var, jobClass, K, Ks, isSiro, out);
            } else if (event == EventType.DEP) {
                for (int k = 0; k < pool.pools[jobClass].length; k++) {
                    int t = pool.pools[jobClass][k];
                    for (int j = 0; j < pool.len[jobClass][k]; j++) {
                        int col = Ks[jobClass] + pool.off[jobClass][k] + j;
                        double nj = srv[col];
                        double ex = pool.exit[jobClass][k][j];
                        if (nj <= 0 || ex <= 0) continue;
                        double[] sd = srv.clone();
                        sd[col] -= 1;
                        List<Object[]> nxt = serveNext(sn, pool, buf, sd, t, Ks, isSiro, sched);
                        for (Object[] o : nxt) {
                            out.add(concat((double[]) o[0], (double[]) o[1], var), ex * nj * (Double) o[2], 1.0,
                                    (Integer) o[3]);
                        }
                    }
                }
            } else if (event == EventType.PHASE) {
                for (int k = 0; k < pool.pools[jobClass].length; k++) {
                    Matrix D0 = pool.D0[jobClass][k];
                    int c0 = Ks[jobClass] + pool.off[jobClass][k];
                    for (int j = 0; j < pool.len[jobClass][k]; j++) {
                        double nj = srv[c0 + j];
                        if (nj <= 0) continue;
                        for (int jd = 0; jd < pool.len[jobClass][k]; jd++) {
                            if (jd == j || D0.get(j, jd) <= 0) continue;
                            double[] sp = srv.clone();
                            sp[c0 + j] -= 1;
                            sp[c0 + jd] += 1;
                            out.add(concat(buf, sp, var), D0.get(j, jd) * nj, 1.0, -1);
                        }
                    }
                }
            }
        }
        int n = out.rows.size();
        if (n == 0) {
            return new Ret.EventResult(new Matrix(0, 0), new Matrix(0, 0), new Matrix(0, 0),
                    new Matrix(0, R), new Matrix(0, R));
        }
        int width = 0;
        for (double[] r : out.rows) width = Math.max(width, r.length);
        Matrix outspace = new Matrix(n, width);
        Matrix outrate = new Matrix(n, 1);
        Matrix outprob = new Matrix(n, 1);
        Matrix outstart = new Matrix(n, R);
        Matrix outpreempt = new Matrix(n, R);
        for (int i = 0; i < n; i++) {
            double[] r = out.rows.get(i);
            int pad = width - r.length; // left-pad narrower buffers
            for (int c = 0; c < r.length; c++) {
                if (r[c] != 0) outspace.set(i, pad + c, r[c]);
            }
            outrate.set(i, 0, event == EventType.ARV ? -1.0 : out.rates.get(i));
            outprob.set(i, 0, out.probs.get(i));
            if (out.starts.get(i) >= 0) outstart.set(i, out.starts.get(i), 1.0);
        }
        if (isSimulation && n > 1) {
            double tot = 0;
            double[] w = new double[n];
            for (int i = 0; i < n; i++) {
                w[i] = event == EventType.ARV ? outprob.get(i, 0) : outrate.get(i, 0);
                tot += w[i];
            }
            double u = Maths.rand() * tot;
            int fc = n - 1;
            double acc = 0;
            for (int i = 0; i < n; i++) {
                acc += w[i];
                if (u < acc) {
                    fc = i;
                    break;
                }
            }
            double rsum = event == EventType.ARV ? -1.0 : tot;
            Matrix os = Matrix.extractRows(outspace, fc, fc + 1, null);
            Matrix orate = new Matrix(1, 1);
            orate.set(0, 0, rsum);
            Matrix oprob = new Matrix(1, 1);
            oprob.set(0, 0, 1.0);
            return new Ret.EventResult(os, orate, oprob, Matrix.extractRows(outstart, fc, fc + 1, null),
                    Matrix.extractRows(outpreempt, fc, fc + 1, null));
        }
        return new Ret.EventResult(outspace, outrate, outprob, outstart, outpreempt);
    }

    private static void arrival(NetworkStruct sn, int ist, CtmcPool pool, double[] buf, double[] srv, double[] var,
                                int cls, int[] K, int[] Ks, boolean isSiro, Out out) {
        int R = sn.nclasses;
        if (pool.pools[cls].length == 0) {
            return; // the class is not served here
        }
        double[] nir = new double[R];
        double ni = 0;
        for (int r = 0; r < R; r++) {
            for (int j = 0; j < K[r]; j++) nir[r] += srv[Ks[r] + j];
        }
        if (isSiro) {
            for (int r = 0; r < R; r++) nir[r] += buf[r];
        } else {
            for (double b : buf) if (b > 0) nir[(int) b - 1] += 1;
        }
        for (int r = 0; r < R; r++) ni += nir[r];
        if (ni >= sn.cap.get(ist) || nir[cls] >= sn.classcap.get(ist, cls)) {
            if (State.isPhysicalCapacity(sn, ist, cls) && State.arrivalIsLost(sn, ist, cls)) {
                out.add(concat(buf, srv, var), -1, 1.0, -1); // lost: the state is unchanged
            }
            return; // blocked, or beyond the state-space cutoff
        }
        int[] busy = busy(pool, srv, Ks);
        List<Integer> cand = new ArrayList<Integer>();
        for (int t : pool.pools[cls]) {
            if (busy[t] < pool.count[t]) cand.add(t);
        }
        if (cand.isEmpty()) {
            double[] b;
            if (isSiro) {
                b = buf.clone();
                b[cls] += 1;
            } else {
                int slot = -1;
                for (int c = buf.length - 1; c >= 0; c--) {
                    if (buf[c] == 0) {
                        slot = c;
                        break;
                    }
                }
                if (slot < 0) {
                    b = new double[buf.length + 1];
                    System.arraycopy(buf, 0, b, 1, buf.length);
                    slot = 0;
                } else {
                    b = buf.clone();
                }
                b[slot] = cls + 1;
            }
            out.add(concat(b, srv, var), -1, 1.0, -1);
            return;
        }
        List<Integer> choice = new ArrayList<Integer>();
        List<Double> pch = new ArrayList<Double>();
        choice.add(cand.get(0));
        pch.add(1.0);
        double[] varc = var.clone();
        if (cand.size() > 1) {
            HeteroSchedPolicy pol = pool.policy;
            if (pol == HeteroSchedPolicy.ALFS) {
                for (int t : pool.alfsorder) {
                    if (cand.contains(t)) {
                        choice.set(0, t);
                        break;
                    }
                }
            } else if (pol == HeteroSchedPolicy.FSF) {
                int best = cand.get(0);
                for (int t : cand) {
                    if (pool.fsfrate[t][cls] > pool.fsfrate[best][cls]) best = t; // first maximum wins ties
                }
                choice.set(0, best);
            } else if (pol == HeteroSchedPolicy.RAIS) {
                choice.clear();
                pch.clear();
                for (int t : cand) {
                    choice.add(t);
                    pch.add(1.0 / cand.size());
                }
            } else if ((pol == HeteroSchedPolicy.ALIS || pol == HeteroSchedPolicy.FAIRNESS) && pool.rotate) {
                int[] ord = pool.perms[(int) var[pool.varpos] - 1];
                for (int k = 0; k < ord.length; k++) {
                    if (cand.contains(ord[k])) {
                        choice.set(0, ord[k]);
                        int[] no = new int[ord.length];
                        int p = 0;
                        for (int q = 0; q < ord.length; q++) {
                            if (q != k) no[p++] = ord[q];
                        }
                        no[p] = ord[k];
                        varc[pool.varpos] = pool.permIndexOf(no);
                        break;
                    }
                }
            }
        }
        for (int c = 0; c < choice.size(); c++) {
            int t = choice.get(c);
            int k = pool.blockOf(cls, t);
            double[] a = pool.alpha[cls][k];
            for (int j = 0; j < a.length; j++) {
                if (a[j] <= 0) continue;
                double[] sa = srv.clone();
                sa[Ks[cls] + pool.off[cls][k] + j] += 1;
                out.add(concat(buf, sa, varc), -1, pch.get(c) * a[j], cls);
            }
        }
    }

    /** Entries {buf, srv, prob, startedClass} of the freed pool-t server taking its next job. */
    private static List<Object[]> serveNext(NetworkStruct sn, CtmcPool pool, double[] buf, double[] srv, int t,
                                            int[] Ks, boolean isSiro, SchedStrategy sched) {
        List<Object[]> res = new ArrayList<Object[]>();
        if (isSiro) {
            double tot = 0;
            List<Integer> elig = new ArrayList<Integer>();
            for (int s = 0; s < buf.length; s++) {
                if (buf[s] > 0 && pool.compat[t][s]) {
                    elig.add(s);
                    tot += buf[s];
                }
            }
            if (elig.isEmpty()) {
                res.add(new Object[]{buf, srv, 1.0, -1});
                return res;
            }
            for (int s : elig) {
                double[] b = buf.clone();
                b[s] -= 1;
                for (Object[] o : start(pool, b, srv, s, t, Ks)) {
                    o[2] = (Double) o[2] * buf[s] / tot;
                    res.add(o);
                }
            }
            return res;
        }
        List<Integer> pos = new ArrayList<Integer>();
        for (int p = 0; p < buf.length; p++) {
            if (buf[p] > 0 && pool.compat[t][(int) buf[p] - 1]) pos.add(p);
        }
        if (pos.isEmpty()) {
            res.add(new Object[]{buf, srv, 1.0, -1});
            return res;
        }
        if (sched == SchedStrategy.HOL || sched == SchedStrategy.FCFSPRIO || sched == SchedStrategy.LCFSPRIO) {
            double best = Double.POSITIVE_INFINITY;
            for (int p : pos) best = Math.min(best, sn.classprio.get((int) buf[p] - 1));
            List<Integer> keep = new ArrayList<Integer>();
            for (int p : pos) {
                if (sn.classprio.get((int) buf[p] - 1) == best) keep.add(p); // a lower value is a higher priority
            }
            pos = keep;
        }
        int p = (sched == SchedStrategy.LCFS || sched == SchedStrategy.LCFSPRIO)
                ? pos.get(0) // leftmost is the newest
                : pos.get(pos.size() - 1); // rightmost is the oldest
        int s = (int) buf[p] - 1;
        double[] b = new double[buf.length];
        System.arraycopy(buf, 0, b, 1, p);
        System.arraycopy(buf, p + 1, b, p + 1, buf.length - p - 1);
        return start(pool, b, srv, s, t, Ks);
    }

    private static List<Object[]> start(CtmcPool pool, double[] buf, double[] srv, int s, int t, int[] Ks) {
        List<Object[]> res = new ArrayList<Object[]>();
        int k = pool.blockOf(s, t);
        double[] a = pool.alpha[s][k];
        for (int j = 0; j < a.length; j++) {
            if (a[j] <= 0) continue;
            double[] ss = srv.clone();
            ss[Ks[s] + pool.off[s][k] + j] += 1;
            res.add(new Object[]{buf, ss, a[j], s});
        }
        return res;
    }

    /**
     * Local states [buffer | servers | pool order] with marginal n at pooled station ind:
     * every split of the jobs into in-service (per class and pool) and waiting ones such that
     * no waiting job has a free compatible server, every phase assignment of the jobs in
     * service, every order of the waiting jobs (per-class counts under SIRO), and every pool
     * order index when the order rotates. The ordered buffer is max(1, min(sum(n), cap)) wide.
     */
    public static Matrix fromMarginalPool(NetworkStruct sn, int ind, Matrix n) {
        int ist = (int) sn.nodeToStation.get(ind);
        int R = sn.nclasses;
        CtmcPool pool = sn.ctmcpool.get(ind);
        int V = (int) Matrix.extractRows(sn.nvars, ind, ind + 1, null).elementSum();
        if (V != (pool.rotate ? 1 : 0)) {
            throw new RuntimeException("Station '" + sn.nodenames.get(ind)
                    + "' carries local variables that a heterogeneous server pool cannot hold.");
        }
        int[] K = new int[R];
        int nsrv = 0;
        int ntot = 0;
        int[] nn = new int[R];
        for (int r = 0; r < R; r++) {
            K[r] = (int) sn.phasessz.get(ist, r);
            nsrv += K[r];
            nn[r] = (int) n.get(r);
            ntot += nn[r];
        }
        boolean isSiro = sn.sched.get(sn.stations.get(ist)) == SchedStrategy.SIRO;
        double cap = sn.cap.get(ist);
        int W = isSiro ? R : (int) Math.max(1, Math.min(ntot, cap));
        boolean unservedBusy = false;
        for (int r = 0; r < R; r++) {
            if (pool.pools[r].length == 0 && nn[r] > 0) unservedBusy = true;
        }
        if (ntot > cap || unservedBusy) {
            return new Matrix(0, W + nsrv + V);
        }
        List<int[][]> allocs = new ArrayList<int[][]>();
        int[][] x = new int[R][];
        for (int r = 0; r < R; r++) x[r] = new int[pool.pools[r].length];
        rec(pool, nn, R, 0, 0, x, new int[pool.ntypes], allocs);
        Comparator<double[]> lex = new Comparator<double[]>() {
            @Override
            public int compare(double[] a, double[] b) {
                for (int c = 0; c < a.length; c++) {
                    int v = Double.compare(a[c], b[c]);
                    if (v != 0) return v;
                }
                return 0;
            }
        };
        TreeSet<double[]> rows = new TreeSet<double[]>(lex);
        int nperm = pool.rotate ? pool.perms.length : 0;
        for (int[][] xa : allocs) {
            int[] w = nn.clone();
            for (int r = 0; r < R; r++) {
                for (int c : xa[r]) w[r] -= c;
            }
            List<double[]> srvset = new ArrayList<double[]>();
            srvset.add(new double[0]);
            for (int r = 0; r < R; r++) {
                List<double[]> blk = new ArrayList<double[]>();
                blk.add(new double[0]);
                for (int k = 0; k < pool.pools[r].length; k++) {
                    List<double[]> comps = compositions(xa[r][k], pool.len[r][k]);
                    blk = product(blk, comps);
                }
                List<double[]> padded = new ArrayList<double[]>();
                for (double[] b : blk) padded.add(Arrays.copyOf(b, K[r]));
                srvset = product(srvset, padded);
            }
            List<double[]> bufset = new ArrayList<double[]>();
            int wsum = 0;
            for (int r = 0; r < R; r++) wsum += w[r];
            if (isSiro) {
                double[] b = new double[R];
                for (int r = 0; r < R; r++) b[r] = w[r];
                bufset.add(b);
            } else if (wsum == 0) {
                bufset.add(new double[W]);
            } else {
                Matrix vi = new Matrix(1, wsum);
                int p = 0;
                for (int r = 0; r < R; r++) {
                    for (int c = 0; c < w[r]; c++) vi.set(0, p++, r + 1);
                }
                Matrix mi = Maths.multisetPerms(vi);
                for (int i = 0; i < mi.getNumRows(); i++) {
                    double[] b = new double[W];
                    for (int c = 0; c < mi.getNumCols(); c++) b[W - mi.getNumCols() + c] = mi.get(i, c);
                    bufset.add(b);
                }
            }
            for (double[] b : bufset) {
                for (double[] s : srvset) {
                    double[] row = concat(b, s, new double[0]);
                    if (pool.rotate) {
                        for (int q = 1; q <= nperm; q++) {
                            double[] rv = Arrays.copyOf(row, row.length + 1);
                            rv[row.length] = q;
                            rows.add(rv);
                        }
                    } else {
                        rows.add(row);
                    }
                }
            }
        }
        Matrix space = new Matrix(rows.size(), W + nsrv + V);
        int i = 0;
        for (double[] row : rows.descendingSet()) {
            for (int c = 0; c < row.length; c++) {
                if (row[c] != 0) space.set(i, c, row[c]);
            }
            i++;
        }
        return space;
    }

    private static void rec(CtmcPool pool, int[] n, int R, int r, int k, int[][] x, int[] load, List<int[][]> out) {
        if (r >= R) {
            for (int s = 0; s < R; s++) {
                if (pool.pools[s].length == 0) continue;
                int used = 0;
                for (int c : x[s]) used += c;
                if (n[s] - used > 0) {
                    for (int t : pool.pools[s]) {
                        if (load[t] < pool.count[t]) return; // a waiting job would have a free server
                    }
                }
            }
            int[][] cp = new int[R][];
            for (int s = 0; s < R; s++) cp[s] = x[s].clone();
            out.add(cp);
            return;
        }
        if (k >= pool.pools[r].length) {
            rec(pool, n, R, r + 1, 0, x, load, out);
            return;
        }
        int t = pool.pools[r][k];
        int used = 0;
        for (int q = 0; q < k; q++) used += x[r][q];
        int cmax = Math.min(n[r] - used, pool.count[t] - load[t]);
        for (int c = 0; c <= cmax; c++) {
            x[r][k] = c;
            load[t] += c;
            rec(pool, n, R, r, k + 1, x, load, out);
            load[t] -= c;
        }
        x[r][k] = 0;
    }

    private static List<double[]> product(List<double[]> a, List<double[]> b) {
        List<double[]> res = new ArrayList<double[]>();
        for (double[] u : a) {
            for (double[] v : b) {
                res.add(concat(u, v, new double[0]));
            }
        }
        return res;
    }

    /** All vectors of PARTS non-negative integers summing to TOTAL. */
    private static List<double[]> compositions(int total, int parts) {
        List<double[]> res = new ArrayList<double[]>();
        if (parts == 1) {
            res.add(new double[]{total});
            return res;
        }
        for (int first = total; first >= 0; first--) {
            for (double[] rest : compositions(total - first, parts - 1)) {
                double[] v = new double[parts];
                v[0] = first;
                System.arraycopy(rest, 0, v, 1, parts - 1);
                res.add(v);
            }
        }
        return res;
    }
}
