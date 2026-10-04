package jline.solvers.ssa.handlers;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import java.util.function.BiFunction;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;

public final class Solver_ssa_nrm_space {
    private Solver_ssa_nrm_space() {}

    public static class SolverSSAResultNRMSpace {
        public final Matrix pi;
        public final Matrix outspace;
        public final Matrix depRates;
        public final NetworkStruct sn;
        public SolverSSAResultNRMSpace(Matrix pi, Matrix outspace, Matrix depRates, NetworkStruct sn) {
            this.pi = pi; this.outspace = outspace; this.depRates = depRates; this.sn = sn;
        }
    }

    /** Holder for the three outputs of {@link #nrm_space}. */
    public static class NrmSpaceResult {
        public final List<Double> times;
        public final Matrix states;
        public final List<ArrayDeque<Integer>[]> bufferStates;
        public NrmSpaceResult(List<Double> times, Matrix states, List<ArrayDeque<Integer>[]> bufferStates) {
            this.times = times; this.states = states; this.bufferStates = bufferStates;
        }
    }

    public static SolverSSAResultNRMSpace solver_ssa_nrm_space(final NetworkStruct sn, final SolverOptions options) {
        RandomManager.setMasterSeed(options.seed);
        int samples = options.samples;
        final int R = sn.nclasses;
        final int I = sn.nnodes;
        Map<jline.lang.nodes.StatefulNode, Matrix> state = sn.state;

        final java.util.Set<SchedStrategy> bufferedSched = new java.util.HashSet<SchedStrategy>(
                java.util.Arrays.asList(SchedStrategy.FCFS, SchedStrategy.LCFS));

        Matrix S = new Matrix(0, I * R);
        final List<Integer> fromIdx = new ArrayList<Integer>();
        final List<List<Integer>> toIdx = new ArrayList<List<Integer>>();
        final List<int[]> fromIR = new ArrayList<int[]>();
        final List<List<Double>> probIR = new ArrayList<List<Double>>();

        int k = 0;
        for (int ind = 0; ind < I; ind++) {
            for (int r = 0; r < R; r++) {
                k++;
                fromIR.add(new int[]{ind, r});
                fromIdx.add(ind * R + r);
                probIR.add(new ArrayList<Double>());
                toIdx.add(new ArrayList<Integer>());
                double[] Srow = new double[I * R];
                if (sn.isslc.get(r) != 0.0) {
                    Srow[fromIdx.get(k - 1)] = GlobalConstants.NegInf;
                } else {
                    Srow[fromIdx.get(k - 1)] = -1.0;
                    for (int jnd = 0; jnd < I; jnd++) {
                        for (int s = 0; s < R; s++) {
                            if (sn.rtnodes.get(ind * R + r, jnd * R + s) > 0) {
                                toIdx.get(k - 1).add(jnd * R + s);
                                double p = sn.rtnodes.get(ind * R + r, jnd * R + s);
                                probIR.get(k - 1).add(p);
                                Srow[jnd * R + s] = Srow[jnd * R + s] + p;
                            }
                        }
                    }
                }
                if (k > S.getNumRows()) {
                    Matrix newS = new Matrix(k, I * R);
                    for (int i = 0; i < S.getNumRows(); i++) {
                        for (int j = 0; j < S.getNumCols(); j++) newS.set(i, j, S.get(i, j));
                    }
                    S = newS;
                }
                for (int j = 0; j < Srow.length; j++) S.set(k - 1, j, Srow[j]);
            }
        }
        S = S.transpose();

        Matrix nvec0 = new Matrix(I * R, 1);
        @SuppressWarnings("unchecked")
        final ArrayDeque<Integer>[] buffers0 = new ArrayDeque[I];
        for (int ind = 0; ind < I; ind++) buffers0[ind] = new ArrayDeque<Integer>();
        for (int ind = 0; ind < I; ind++) {
            if (sn.isstateful.get(ind) != 0.0) {
                Matrix state_i = state.get(sn.stateful.get((int) sn.nodeToStateful.get(ind)));
                if (state_i == null) {
                    throw new RuntimeException("State matrix for stateful node " + ind + " is null");
                }
                State.StateMarginalStatistics aggr = ToMarginal.toMarginalAggr(sn,
                        ind, state_i, null, null, null, null, null);
                for (int r = 0; r < R; r++) {
                    double nir = aggr.nir.get(r);
                    if (Double.isInfinite(nir)) {
                        if (sn.nodetype.get(ind) == NodeType.Source) nir = 1.0;
                        else throw new RuntimeException("Infinite population error.");
                    }
                    nvec0.set(ind * R + r, 0, nir);
                }

                // see _kb/06-solver-catalog.md for rationale
                int ist = (int) sn.nodeToStation.get(ind);
                if (ist >= 0 && bufferedSched.contains(sn.sched.get(sn.stations.get(ist)))) {
                    Matrix Kmat = new Matrix(1, sn.phasessz.getNumCols());
                    Matrix.extract(sn.phasessz, ist, ist + 1, 0, sn.phasessz.getNumCols(), Kmat, 0, 0);
                    int sumK = (int) Kmat.elementSum();
                    int sumNvars = (int) sn.nvars.sumRows(ind);
                    int bufCols = state_i.getNumCols() - sumK - sumNvars;
                    for (int pos = 0; pos < bufCols; pos++) {
                        int classId = (int) state_i.get(0, pos);
                        if (classId >= 1 && classId <= R) {
                            buffers0[ind].addLast(classId);
                        }
                        // classId == 0 means empty position, skip
                    }
                }
            }
        }

        final double[] mi = new double[I];
        final double[][] rates = new double[I][R];
        for (int ind = 0; ind < I; ind++) {
            if (sn.isstation.get(ind) != 0.0) {
                int ist = (int) sn.nodeToStation.get(ind);
                for (int r = 0; r < R; r++) {
                    double muir = sn.rates.get(ist, r);
                    if (!Double.isNaN(muir)) rates[ind][r] = muir;
                }
                mi[ind] = sn.nservers.get(ist);
            } else {
                for (int r = 0; r < R; r++) {
                    rates[ind][r] = GlobalConstants.Immediate;
                    mi[ind] = (double) GlobalConstants.MaxInt;
                }
            }
            if (Double.isInfinite(mi[ind])) mi[ind] = (double) GlobalConstants.MaxInt;
        }

        final double epstol = GlobalConstants.Zero;
        @SuppressWarnings("unchecked")
        final BiFunction<Matrix, ArrayDeque<Integer>[], Double>[] a = new BiFunction[fromIdx.size()];
        for (int j = 0; j < fromIdx.size(); j++) {
            final int idxJ = j;
            a[j] = new BiFunction<Matrix, ArrayDeque<Integer>[], Double>() {
                @Override
                public Double apply(Matrix X, ArrayDeque<Integer>[] bufs) {
                    int ind = fromIR.get(idxJ)[0];
                    int r = fromIR.get(idxJ)[1];
                    double base = rates[ind][r];
                    if (sn.isstation.get(ind) == 0.0) {
                        return base * Math.min(1.0, X.get(fromIdx.get(idxJ), 0));
                    }
                    int ist = (int) sn.nodeToStation.get(ind);
                    SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
                    if (sched == SchedStrategy.EXT) return base;
                    if (sched == SchedStrategy.INF)
                        return base * X.get(fromIdx.get(idxJ), 0)
                                * Solver_ssa_nrm.cdfac(sn, ist, X, ind, R, r);
                    if (sched == SchedStrategy.PS || sched == SchedStrategy.DPS || sched == SchedStrategy.GPS) {
                        if (R == 1) return base * Math.min(mi[ind], X.get(fromIdx.get(idxJ), 0))
                                * Solver_ssa_nrm.lldfac(sn, ist, X.get(fromIdx.get(idxJ), 0))
                                * Solver_ssa_nrm.cdfac(sn, ist, X, ind, R, r);
                        double total = 0.0;
                        for (int rr = 0; rr < R; rr++) total += X.get(ind * R + rr, 0);
                        total += epstol;
                        return base * (X.get(fromIdx.get(idxJ), 0) / total) * Math.min(mi[ind], total)
                                * Solver_ssa_nrm.lldfac(sn, ist, total)
                                * Solver_ssa_nrm.cdfac(sn, ist, X, ind, R, r);
                    }
                    if (sched == SchedStrategy.FCFS || sched == SchedStrategy.LCFS) {
                        // Invariant: buffers[ind].size == max(0, total - mi[ind])
                        // State update must ensure that buffer only fills when all servers busy.
                        int waiting = 0;
                        for (Integer cl : bufs[ind]) {
                            if (cl == r + 1) waiting++;
                        }
                        double inService = X.get(fromIdx.get(idxJ), 0) - waiting;
                        if (inService <= 0) return 0.0;
                        double total = 0.0;
                        for (int rr = 0; rr < R; rr++) total += X.get(ind * R + rr, 0);
                        // rate proportional to jobs actually being served, scaled by
                        // the load-dependent factor at the total station population
                        return base * inService * Solver_ssa_nrm.lldfac(sn, ist, total)
                                * Solver_ssa_nrm.cdfac(sn, ist, X, ind, R, r);
                    }
                    if (sched == SchedStrategy.PSPRIO || sched == SchedStrategy.DPSPRIO || sched == SchedStrategy.GPSPRIO) {
                        if (R == 1) return base * Math.min(mi[ind], X.get(fromIdx.get(idxJ), 0));
                        double total = 0.0;
                        for (int rr = 0; rr < R; rr++) total += X.get(ind * R + rr, 0);
                        total += epstol;
                        int minPrio = Integer.MAX_VALUE;
                        for (int rr = 0; rr < R; rr++) {
                            if (X.get(ind * R + rr, 0) > 0) {
                                int classPrio = (int) sn.classprio.get(rr);
                                if (classPrio < minPrio) minPrio = classPrio;
                            }
                        }
                        int classPrio = (int) sn.classprio.get(r);
                        if (total <= mi[ind] || classPrio == minPrio) {
                            return base * (X.get(fromIdx.get(idxJ), 0) / total) * Math.min(mi[ind], total);
                        }
                        return 0.0;
                    }
                    return base * X.get(fromIdx.get(idxJ), 0);
                }
            };
        }

        List<List<Integer>> D = new ArrayList<List<Integer>>(S.getNumCols());
        for (int kk = 0; kk < S.getNumCols(); kk++) D.add(new ArrayList<Integer>());
        for (int kk = 0; kk < D.size(); kk++) {
            List<Integer> J = new ArrayList<Integer>();
            for (int i = 0; i < S.getNumRows(); i++) {
                if (S.get(i, kk) != 0.0) J.add(i);
            }
            List<Integer> vecd = new ArrayList<Integer>();
            for (int j = 0; j < J.size(); j++) {
                int pos = J.get(j);
                int r = pos % R;
                int ind = (pos - r) / R;
                for (int rr = 0; rr < R; rr++) vecd.add(ind * R + rr);
            }
            if (!vecd.isEmpty()) {
                List<Integer> vecdUnique = new ArrayList<Integer>(new java.util.LinkedHashSet<Integer>(vecd));
                List<Integer> vecs = new ArrayList<Integer>();
                for (int j = 0; j < vecdUnique.size(); j++) {
                    for (int kk2 = 0; kk2 < S.getNumCols(); kk2++) {
                        if (S.get(vecdUnique.get(j), kk2) < 0.0) vecs.add(kk2);
                    }
                }
                List<Integer> vecsUnique = new ArrayList<Integer>(new java.util.LinkedHashSet<Integer>(vecs));
                D.get(kk).addAll(vecsUnique);
            }
        }

        for (int i = 0; i < S.getNumRows(); i++) {
            for (int j = 0; j < S.getNumCols(); j++) {
                if (Double.isInfinite(S.get(i, j))) S.set(i, j, 0.0);
            }
        }

        Map<String, double[]> reactCache = new HashMap<String, double[]>();

        NrmSpaceResult result = nrm_space(S, D, a, nvec0, buffers0, samples, options, reactCache, fromIR, mi, R, sn);
        List<Double> t = result.times;
        Matrix nvecsim = result.states;
        List<ArrayDeque<Integer>[]> bufferStates = result.bufferStates;

        double[] dt = new double[t.size() - 1];
        for (int i = 0; i < dt.length; i++) dt[i] = t.get(i + 1) - t.get(i);

        List<double[]> stateCols = new ArrayList<double[]>();
        for (int c = 0; c < nvecsim.getNumCols() - 1; c++) {
            double[] arr = new double[nvecsim.getNumRows()];
            for (int rIdx = 0; rIdx < nvecsim.getNumRows(); rIdx++) arr[rIdx] = nvecsim.get(rIdx, c);
            stateCols.add(arr);
        }

        // Find unique states keyed on both nvec and buffer contents
        LinkedHashMap<String, Integer> uniq = new LinkedHashMap<String, Integer>();
        List<double[]> outspaceRows = new ArrayList<double[]>();
        List<ArrayDeque<Integer>[]> outspaceBuffers = new ArrayList<ArrayDeque<Integer>[]>();
        int[] ic = new int[stateCols.size()];
        for (int i = 0; i < stateCols.size(); i++) {
            double[] arr = stateCols.get(i);
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < arr.length; j++) {
                if (j > 0) sb.append(',');
                sb.append(arr[j]);
            }
            String bufKey = bufferHashAll(bufferStates.get(i));
            String key = sb.toString() + "|" + bufKey;
            Integer existing = uniq.get(key);
            if (existing == null) {
                outspaceRows.add(arr);
                outspaceBuffers.add(bufferStates.get(i));
                int idx = outspaceRows.size() - 1;
                uniq.put(key, idx);
                ic[i] = idx;
            } else {
                ic[i] = existing;
            }
        }

        double[] timeAccum = new double[outspaceRows.size()];
        for (int i = 0; i < ic.length; i++) timeAccum[ic[i]] += dt[i];
        double total = 0.0;
        for (double v : timeAccum) total += v;
        Matrix pi = new Matrix(1, timeAccum.length);
        for (int i = 0; i < timeAccum.length; i++) pi.set(0, i, timeAccum[i] / total);

        Matrix outspace = new Matrix(outspaceRows.size(), I * R);
        for (int i = 0; i < outspaceRows.size(); i++) {
            double[] row = outspaceRows.get(i);
            for (int j = 0; j < row.length; j++) outspace.set(i, j, row[j]);
        }

        int numStates = outspace.getNumRows();
        Matrix depRates = new Matrix(numStates, I * R);
        for (int st = 0; st < numStates; st++) {
            Matrix stateVec = new Matrix(I * R, 1);
            for (int i = 0; i < I * R; i++) stateVec.set(i, 0, outspace.get(st, i));
            double[] a_state = reactCache.get(hashState(stateVec, outspaceBuffers.get(st), sn));
            if (a_state == null) continue;
            for (int j = 0; j < fromIdx.size(); j++) {
                depRates.set(st, fromIdx.get(j), depRates.get(st, fromIdx.get(j)) + a_state[j]);
            }
        }
        return new SolverSSAResultNRMSpace(pi, outspace, depRates, sn);
    }

    // ======================================================================
    // Next-Reaction Method core
    // ======================================================================
    public static NrmSpaceResult nrm_space(Matrix S,
                                           List<List<Integer>> D,
                                           BiFunction<Matrix, ArrayDeque<Integer>[], Double>[] a,
                                           Matrix nvec0,
                                           ArrayDeque<Integer>[] buffers0,
                                           int samples,
                                           SolverOptions options,
                                           Map<String, double[]> reactcache,
                                           List<int[]> fromIR,
                                           double[] mi,
                                           int R,
                                           NetworkStruct sn) {
        int numReactions = S.getNumCols();
        @SuppressWarnings("unchecked")
        ArrayDeque<Integer>[] buffers = new ArrayDeque[buffers0.length];
        for (int i = 0; i < buffers0.length; i++) buffers[i] = new ArrayDeque<Integer>(buffers0[i]);

        double[] Ak = new double[numReactions];
        for (int i = 0; i < numReactions; i++) Ak[i] = a[i].apply(nvec0, buffers0);
        double[] Pk = new double[numReactions];
        for (int i = 0; i < numReactions; i++) Pk[i] = -Math.log(Maths.rand());
        double[] Tk = new double[numReactions];
        List<Double> times = new ArrayList<Double>();
        times.add(0.0);
        List<double[]> states = new ArrayList<double[]>();
        List<ArrayDeque<Integer>[]> bufferStates = new ArrayList<ArrayDeque<Integer>[]>();
        bufferStates.add(copyBuffers(buffers));
        Matrix nvec = nvec0.copy();
        reactcache.clear();
        reactcache.put(hashState(nvec, buffers, sn), Ak.clone());

        final Routing rt = buildRouting(S);

        int n = 0;
        while (n < samples) {
            double[] tau = new double[numReactions];
            for (int i = 0; i < numReactions; i++) {
                tau[i] = (Ak[i] > 0) ? (Pk[i] - Tk[i]) / Ak[i] : Double.POSITIVE_INFINITY;
            }
            int kfire = 0;
            double minVal = tau[0];
            for (int i = 1; i < numReactions; i++) {
                if (tau[i] < minVal) { minVal = tau[i]; kfire = i; }
            }
            double dt = tau[kfire];
            times.add(times.get(times.size() - 1) + dt);

            int srcRow = fromIR.get(kfire)[0] * R + fromIR.get(kfire)[1];
            int destPos = fireReaction(kfire, nvec, S, rt, srcRow);

            updateBuffers(kfire, nvec, buffers, fromIR, destPos, mi, R, sn);

            for (int i = 0; i < numReactions; i++) Tk[i] += Ak[i] * dt;
            for (int i : D.get(kfire)) Ak[i] = a[i].apply(nvec, buffers);
            Pk[kfire] -= Math.log(Maths.rand());

            reactcache.put(hashState(nvec, buffers, sn), Ak.clone());
            double[] colState = new double[nvec.getNumRows()];
            for (int rIdx = 0; rIdx < nvec.getNumRows(); rIdx++) colState[rIdx] = nvec.get(rIdx, 0);
            states.add(colState);
            bufferStates.add(copyBuffers(buffers));
            n++;
            printProgress(options, n);
        }

        Matrix result = new Matrix(S.getNumRows(), states.size() + 1);
        for (int rIdx = 0; rIdx < nvec0.getNumRows(); rIdx++) result.set(rIdx, 0, nvec0.get(rIdx, 0));
        for (int c = 0; c < states.size(); c++) {
            double[] col = states.get(c);
            for (int rIdx = 0; rIdx < col.length; rIdx++) result.set(rIdx, c + 1, col[rIdx]);
        }
        return new NrmSpaceResult(times, result, bufferStates);
    }

    private static ArrayDeque<Integer>[] copyBuffers(ArrayDeque<Integer>[] buffers) {
        @SuppressWarnings("unchecked")
        ArrayDeque<Integer>[] copy = new ArrayDeque[buffers.length];
        for (int i = 0; i < buffers.length; i++) copy[i] = new ArrayDeque<Integer>(buffers[i]);
        return copy;
    }

    /** Hash key combining the state vector and the buffer contents of FCFS/LCFS nodes. */
    private static String hashState(Matrix v, ArrayDeque<Integer>[] bufs, NetworkStruct sn) {
        StringBuilder sb = new StringBuilder();
        boolean first = true;
        for (int ind = 0; ind < bufs.length; ind++) {
            if (isFCFS(ind, sn) || isLCFS(ind, sn)) {
                if (!first) sb.append('|');
                first = false;
                sb.append(ind).append(":[").append(joinDeque(bufs[ind])).append(']');
            }
        }
        return v.toString() + "|" + sb.toString();
    }

    /** Buffer hash over all node indices (used to distinguish unique states). */
    private static String bufferHashAll(ArrayDeque<Integer>[] bufs) {
        StringBuilder sb = new StringBuilder();
        for (int ind = 0; ind < bufs.length; ind++) {
            if (ind > 0) sb.append('|');
            sb.append(ind).append(":[").append(joinDeque(bufs[ind])).append(']');
        }
        return sb.toString();
    }

    private static String joinDeque(ArrayDeque<Integer> deque) {
        StringBuilder sb = new StringBuilder();
        boolean first = true;
        for (Integer x : deque) {
            if (!first) sb.append(',');
            first = false;
            sb.append(x);
        }
        return sb.toString();
    }

    private static void updateBuffers(
            int kfire,
            Matrix nvec,
            ArrayDeque<Integer>[] buffers,
            List<int[]> fromIR,
            int destPos,
            double[] mi,
            int R,
            NetworkStruct sn) {
        int ind = fromIR.get(kfire)[0];

        // Handle departure from FCFS/LCFS source node
        if (isFCFS(ind, sn)) {
            buffers[ind].pollLast();
        } else if (isLCFS(ind, sn)) {
            buffers[ind].pollFirst();
        }

        // Destination is the (state-row) position selected by the routing draw
        if (destPos < 0) {
            return;
        }
        int jnd = destPos / R;
        int r = destPos % R;

        // Handle arrival at buffered destination node
        if (isFCFS(jnd, sn) || isLCFS(jnd, sn)) {
            double totalAtDest = 0.0;
            for (int i = 0; i < R; i++) totalAtDest += nvec.get(jnd * R + i, 0);
            if (totalAtDest > mi[jnd]) {
                // All servers busy - arriving job joins back of buffer
                buffers[jnd].addFirst(r + 1);
            }
            // Otherwise job went straight into service, buffer unchanged
        }
    }

    // Per-reaction routing for inverse-CDF destination sampling. Mirrors MATLAB
    // next_reaction_method_direct: each firing moves exactly one job to a single
    // stochastically chosen destination (keeping the marginal state integer), so
    // immediate pass-through nodes (e.g. ClassSwitch) are not drained below zero.
    private static final class Routing {
        final int[] nnzP;
        final int[][] destRow;
        final double[][] cdf;
        Routing(int n) { nnzP = new int[n]; destRow = new int[n][]; cdf = new double[n][]; }
    }

    private static Routing buildRouting(Matrix S) {
        int numReactions = S.getNumCols();
        int rows = S.getNumRows();
        Routing rt = new Routing(numReactions);
        for (int k = 0; k < numReactions; k++) {
            List<Integer> dest = new ArrayList<Integer>();
            List<Double> pr = new ArrayList<Double>();
            for (int row = 0; row < rows; row++) {
                double v = S.get(row, k);
                double p = (v < 0.0) ? v + 1.0 : v; // P = S; P(P<0) += 1
                if (p > 0.0) { dest.add(row); pr.add(p); }
            }
            int nd = dest.size();
            rt.nnzP[k] = nd;
            rt.destRow[k] = new int[nd];
            rt.cdf[k] = new double[nd];
            double cum = 0.0;
            for (int x = 0; x < nd; x++) {
                rt.destRow[k][x] = dest.get(x);
                cum += pr.get(x);
                rt.cdf[k][x] = cum;
            }
        }
        return rt;
    }

    private static int fireReaction(int kfire, Matrix nvec, Matrix S, Routing rt, int srcRow) {
        if (rt.nnzP[kfire] > 1) {
            double u = Maths.rand();
            int sel = rt.cdf[kfire].length - 1;
            for (int x = 0; x < rt.cdf[kfire].length; x++) {
                if (rt.cdf[kfire][x] > u) { sel = x; break; }
            }
            nvec.set(srcRow, 0, nvec.get(srcRow, 0) - 1.0);
            int destPos = rt.destRow[kfire][sel];
            nvec.set(destPos, 0, nvec.get(destPos, 0) + 1.0);
            return destPos;
        }
        nvec.addEq(S.getColumn(kfire));
        return (rt.nnzP[kfire] == 1) ? rt.destRow[kfire][0] : -1;
    }

    private static boolean isFCFS(int ind, NetworkStruct sn) {
        if (sn.isstation.get(ind) == 0.0) return false;
        int ist = (int) sn.nodeToStation.get(ind);
        return sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFS;
    }

    private static boolean isLCFS(int ind, NetworkStruct sn) {
        if (sn.isstation.get(ind) == 0.0) return false;
        int ist = (int) sn.nodeToStation.get(ind);
        return sn.sched.get(sn.stations.get(ist)) == SchedStrategy.LCFS;
    }

    private static void printProgress(SolverOptions options, int samples_collected) {
        // The solver console owns the line while it narrates: the in-place
        // backspace counter below cannot be rewritten in a paged log, so the
        // progress is reported as decimated rows instead.
        if (jline.io.LineConsole.ownsLog()) {
            final long every = Math.max(1L, options.samples / 20L);
            if (samples_collected % every == 0) {
                jline.io.LineConsole.iter(samples_collected / every,
                        "simulated %d of %d samples (%.0f%%)", samples_collected,
                        options.samples, 100.0 * samples_collected / options.samples);
            }
            return;
        }
        if (System.console() != null && !"parallel".equals(options.method)
                && (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG)) {
            if (samples_collected == 2) {
                System.out.printf("\nSSA samples: %9d ", samples_collected);
                System.out.flush();
            } else if (samples_collected % 1000 == 0) {
                System.out.printf("\b\b\b\b\b\b\b\b\b\b %9d", samples_collected);
                System.out.flush();
            }
            if (samples_collected == options.samples) System.out.println();
        }
    }
}
