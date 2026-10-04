package jline.solvers.ctmc.handlers;

import jline.api.mam.Map_pie;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.StatefulNode;
import jline.lang.nodes.Station;
import jline.lang.state.CtmcMapc;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Pair form of a multiserver MAP (and MMPP2) service for the CTMC state space; port of MATLAB
 * solver_ctmc_mapc.m.
 *
 * <p>JMT and both LDES engines sample a MAP service through ONE sampler per (station, class):
 * every draw starts in the phase the previous draw ENDED in, draws chained in service-start
 * order, idle periods included. With c &gt; 1 servers the next start can happen while earlier
 * draws are still in progress, so the landing phase of a draw must be known when it starts.
 * With V = (-D0)^-1 D1 a draw started in h ends in j w.p. V(h,j); a busy server is a PAIR
 * (i,j) that moves i-&gt;k at D0(i,k)V(k,j)/V(i,j) and completes at D1(i,j)/V(i,j). The class
 * memory variable is the landing of the most recently STARTED draw: a start from h enters
 * (h,j) w.p. V(h,j) and sets h := j, while moves and completions keep it. For a renewal MAP
 * the chain reduces in law to PH/c; single-server stations are left untouched.</p>
 *
 * <p>The pair law is stored as the lifted MAP D0p (conditioned moves) and
 * D1p((i,j),(j,j')) = D1(i,j)/V(i,j)*V(j,j'), equivalent in law to the original; the
 * bookkeeping read by the state handlers is {@code sn.ctmcmapc}. The rewrite works on a
 * shallow copy whose mutated containers are replaced by copies.</p>
 */
public final class Solver_ctmc_mapc {

    private Solver_ctmc_mapc() {
    }

    /** The struct with every multiserver MAP service in pair form, or sn itself when there is none. */
    public static NetworkStruct apply(NetworkStruct sn) {
        int R = sn.nclasses;
        List<int[]> todo = new ArrayList<int[]>();
        for (int ist = 0; ist < sn.nstations; ist++) {
            int ind = (int) sn.stationToNode.get(ist);
            Station station = sn.stations.get(ist);
            double c = sn.nservers.get(ist);
            if (sn.nodetype.get(ind) == NodeType.Source || sn.sched.get(station) != SchedStrategy.FCFS
                    || Double.isInfinite(c) || c <= 1) {
                continue;
            }
            for (int r = 0; r < R; r++) {
                ProcessType pt = sn.procid.get(station).get(sn.jobclasses.get(r));
                if ((pt == ProcessType.MAP || pt == ProcessType.MMPP2) && CtmcMapc.get(sn, ist, r) == null) {
                    todo.add(new int[]{ist, r});
                }
            }
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
        sn.phases = old.phases.copy();
        sn.phasessz = old.phasessz.copy();
        sn.phaseshift = old.phaseshift.copy();
        sn.state = old.state == null ? null : new HashMap<StatefulNode, Matrix>(old.state);
        sn.ctmcmapc = new HashMap<Integer, CtmcMapc[]>();
        if (old.ctmcmapc != null) {
            for (Map.Entry<Integer, CtmcMapc[]> e : old.ctmcmapc.entrySet()) {
                sn.ctmcmapc.put(e.getKey(), e.getValue().clone());
            }
        }
        List<Integer> changed = new ArrayList<Integer>();
        for (int[] ir : todo) {
            int ist = ir[0];
            int r = ir[1];
            Station station = sn.stations.get(ist);
            JobClass jc = sn.jobclasses.get(r);
            MatrixCell law = sn.proc.get(station).get(jc);
            Matrix D0 = law.get(0);
            Matrix D1 = law.get(1);
            int p = D0.getNumRows();
            Matrix V = D0.scale(-1.0).inv().mult(D1);
            List<int[]> pl = new ArrayList<int[]>();
            for (int i = 0; i < p; i++) {
                for (int j = 0; j < p; j++) {
                    if (Math.abs(V.get(i, j)) < 1e-14) V.set(i, j, 0);
                    if (V.get(i, j) > 0) pl.add(new int[]{i, j});
                }
            }
            int T = pl.size();
            int[][] pairs = pl.toArray(new int[T][]);
            int[][] tp = new int[p][p];
            for (int[] row : tp) java.util.Arrays.fill(row, -1);
            for (int t = 0; t < T; t++) tp[pairs[t][0]][pairs[t][1]] = t;
            double[] done = new double[T];
            Matrix D0p = new Matrix(T, T);
            Matrix D1p = new Matrix(T, T);
            Matrix mu = new Matrix(T, 1);
            Matrix phi = new Matrix(T, 1);
            for (int t = 0; t < T; t++) {
                int i = pairs[t][0];
                int j = pairs[t][1];
                done[t] = D1.get(i, j) / V.get(i, j);
                D0p.set(t, t, D0.get(i, i));
                for (int k = 0; k < p; k++) {
                    if (k != i && D0.get(i, k) != 0 && tp[k][j] >= 0) {
                        D0p.set(t, tp[k][j], D0.get(i, k) * V.get(k, j) / V.get(i, j));
                    }
                }
                for (int jn = 0; jn < p; jn++) {
                    if (tp[j][jn] >= 0) D1p.set(t, tp[j][jn], done[t] * V.get(j, jn));
                }
                mu.set(t, 0, -D0.get(i, i));
                phi.set(t, 0, done[t] / (-D0.get(i, i)));
            }
            sn.proc.get(station).put(jc, new MatrixCell(D0p, D1p));
            sn.pie.get(station).put(jc, Map_pie.map_pie(D0p, D1p));
            sn.mu.get(station).put(jc, mu);
            sn.phi.get(station).put(jc, phi);
            sn.phases.set(ist, r, T);
            sn.phasessz.set(ist, r, T);
            CtmcMapc[] row = sn.ctmcmapc.get(ist);
            if (row == null) {
                row = new CtmcMapc[R];
                sn.ctmcmapc.put(ist, row);
            }
            row[r] = new CtmcMapc(p, pairs, V, done);
            if (!changed.contains(ist)) changed.add(ist);
        }
        for (int ist : changed) {
            int shift = 0;
            for (int r = 0; r < R; r++) {
                sn.phaseshift.set(ist, r, shift);
                shift += (int) sn.phasessz.get(ist, r);
            }
            if (sn.phaseshift.getNumCols() > R) sn.phaseshift.set(ist, R, shift);
            int ind = (int) sn.stationToNode.get(ist);
            int isf = (int) sn.nodeToStateful.get(ind);
            StatefulNode node = sn.stateful.get(isf);
            if (sn.state != null && sn.state.get(node) != null && sn.state.get(node).getNumRows() > 0) {
                sn.state.put(node, rebuildState(old, sn, ind, ist, sn.state.get(node)));
            }
        }
        return sn;
    }

    /**
     * A job in service in MAP phase i is placed in the first pair (i,j); the buffer and the
     * local variables, the carried phase included, are unchanged.
     */
    private static Matrix rebuildState(NetworkStruct old, NetworkStruct sn, int ind, int ist, Matrix oldst) {
        int R = sn.nclasses;
        int nsrvOld = 0;
        int nsrv = 0;
        for (int r = 0; r < R; r++) {
            nsrvOld += (int) old.phasessz.get(ist, r);
            nsrv += (int) sn.phasessz.get(ist, r);
        }
        int V = (int) Matrix.extractRows(old.nvars, ind, ind + 1, null).elementSum();
        int W = oldst.getNumCols() - nsrvOld - V;
        Matrix st = new Matrix(oldst.getNumRows(), W + nsrv + V);
        for (int row = 0; row < oldst.getNumRows(); row++) {
            for (int c = 0; c < W; c++) st.set(row, c, oldst.get(row, c));
            for (int c = 0; c < V; c++) st.set(row, W + nsrv + c, oldst.get(row, W + nsrvOld + c));
            for (int r = 0; r < R; r++) {
                int oo = W + (int) old.phaseshift.get(ist, r);
                int no = W + (int) sn.phaseshift.get(ist, r);
                CtmcMapc mc = CtmcMapc.get(sn, ist, r);
                for (int i = 0; i < (int) old.phasessz.get(ist, r); i++) {
                    double cnt = oldst.get(row, oo + i);
                    if (cnt == 0) continue;
                    int t = i;
                    if (mc != null) {
                        for (t = 0; t < mc.pairs.length; t++) {
                            if (mc.pairs[t][0] == i) break;
                        }
                    }
                    st.set(row, no + t, st.get(row, no + t) + cnt);
                }
            }
        }
        return st;
    }

    private static <A, B, C> Map<A, Map<B, C>> copy2(Map<A, Map<B, C>> m) {
        if (m == null) return null;
        Map<A, Map<B, C>> out = new HashMap<A, Map<B, C>>();
        for (Map.Entry<A, Map<B, C>> e : m.entrySet()) {
            out.put(e.getKey(), e.getValue() == null ? null : new HashMap<B, C>(e.getValue()));
        }
        return out;
    }
}
