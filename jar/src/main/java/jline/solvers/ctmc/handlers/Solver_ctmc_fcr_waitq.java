package jline.solvers.ctmc.handlers;

import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.Sync;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.EventType;
import jline.lang.constant.NodeType;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Reachability-based state space and per-action rate filters for models with
 * a finite capacity region (FCR) whose drop rule is WAITQ (waiting queue).
 *
 * Port of MATLAB solver_ctmc_fcr_waitq.m. JMT semantics (mirrored by LDES): a
 * job refused entry to a full region leaves the upstream station and waits in
 * a per-region FIFO of (class, destination) tokens outside the region; after
 * every transition that frees region capacity, tokens are released strictly
 * in FIFO order (head-of-line) as long as the admission constraints (global
 * cap, per-class caps, memory budget, linear constraints A*x&lt;=b) permit; a
 * fresh cap-admissible arrival overtakes a stuck head. Blocked jobs are
 * counted neither in the region occupancy nor in any station state.
 * Class-switching hops between two members of the same region are an exit
 * (with FIFO release) followed by a gated re-entry of the new class. Classes
 * whose region rule is DROP keep the transition-censoring behavior of the
 * default generator.
 *
 * The CTMC state is augmented as [h(0:nstateful-1), buf_0, ..., buf_{F-1}]
 * where h are the per-node hashed states and buf_f is the token FIFO of
 * region f, padded with -1 to its maximum length.
 */
public final class Solver_ctmc_fcr_waitq {

    private Solver_ctmc_fcr_waitq() {
    }

    /**
     * Result bundle: augmented spaces, the per-action rate filters, and the true-BAS
     * become-blocked arcs. The latter are part of the generator but are NOT departures
     * of any action, so they are kept out of Dfilt and folded straight into Q.
     */
    public static final class Result {
        public final Matrix stateSpace;
        public final Matrix stateSpaceAggr;
        public final Matrix stateSpaceHashed;
        public final MatrixCell Dfilt;
        public final Matrix basBlockQ;

        Result(Matrix stateSpace, Matrix stateSpaceAggr, Matrix stateSpaceHashed, MatrixCell Dfilt,
               Matrix basBlockQ) {
            this.stateSpace = stateSpace;
            this.stateSpaceAggr = stateSpaceAggr;
            this.stateSpaceHashed = stateSpaceHashed;
            this.Dfilt = Dfilt;
            this.basBlockQ = basBlockQ;
        }
    }

    public static Result build(NetworkStruct sn, SolverOptions options) {
        final int nstateful = sn.nstateful;
        final int K = sn.nclasses;
        final int nnodes = sn.nnodes;
        final int F = sn.nregions;
        final Map<Integer, Sync> sync = sn.sync;
        final int A = sync.size();
        final int local = sn.nnodes + 1;

        // feature gates
        if (sn.gsync != null && !sn.gsync.isEmpty()) {
            throw new RuntimeException("WAITQ finite capacity regions are not supported together with stochastic Petri net transitions in SolverCTMC.");
        }
        if (sn.fjsync != null && !sn.fjsync.isEmpty()) {
            throw new RuntimeException("WAITQ finite capacity regions are not supported together with fork-join in SolverCTMC.");
        }
        if (sn.isstatedep != null) {
            for (int ind = 0; ind < nnodes; ind++) {
                if (sn.isstatedep.get(ind, 2) != 0.0) {
                    throw new RuntimeException("WAITQ finite capacity regions are not supported together with state-dependent routing in SolverCTMC.");
                }
            }
        }

        // region data
        final int M = sn.nstations;
        final boolean[][] memberMask = new boolean[F][M];
        final double[][] ccap = new double[F][K];
        final double[] gcap = new double[F];
        final double[] memcap = new double[F];
        final double[][] szrow = new double[F][K];
        final Matrix[] linA = new Matrix[F];
        final Matrix[] linb = new Matrix[F];
        final boolean[][] iswaitq = new boolean[F][K];
        final int dropId = DropStrategy.Drop.getID();
        for (int f = 0; f < F; f++) {
            Matrix Rmat = sn.region.get(f); // M x (K+1)
            Matrix memMat = (sn.regionmaxmem != null && sn.regionmaxmem.size() > f) ? sn.regionmaxmem.get(f) : null;
            List<Integer> members = new ArrayList<Integer>();
            for (int i = 0; i < M; i++) {
                boolean isMember = false;
                for (int col = 0; col <= K; col++) {
                    if (Rmat.get(i, col) != -1) {
                        isMember = true;
                        break;
                    }
                }
                // membership: any job-count cap OR the region memory budget set
                // on the station row (a memory-only region has all caps at -1)
                if (!isMember && memMat != null && memMat.get(i, 0) != -1) {
                    isMember = true;
                }
                if (isMember) {
                    members.add(i);
                    memberMask[f][i] = true;
                }
            }
            gcap[f] = Double.POSITIVE_INFINITY;
            memcap[f] = Double.POSITIVE_INFINITY;
            for (int r = 0; r < K; r++) {
                ccap[f][r] = Double.POSITIVE_INFINITY;
                for (int ii = 0; ii < members.size(); ii++) {
                    double v = Rmat.get(members.get(ii), r);
                    if (v != -1) {
                        ccap[f][r] = Math.min(ccap[f][r], v);
                    }
                }
                iswaitq[f][r] = (sn.regionrule == null || sn.regionrule.isEmpty())
                        || sn.regionrule.get(f, r) != dropId;
                szrow[f][r] = (sn.regionsz != null && !sn.regionsz.isEmpty()) ? sn.regionsz.get(f, r) : 1.0;
            }
            for (int ii = 0; ii < members.size(); ii++) {
                double v = Rmat.get(members.get(ii), K);
                if (v != -1) {
                    gcap[f] = Math.min(gcap[f], v);
                }
                if (memMat != null) {
                    double mv = memMat.get(members.get(ii), 0);
                    if (mv != -1) {
                        memcap[f] = Math.min(memcap[f], mv);
                    }
                }
            }
            if (sn.regionlincon != null && sn.regionlincon.containsKey(f)) {
                MatrixCell ab = sn.regionlincon.get(f);
                if (ab != null && ab.size() >= 2 && ab.get(0) != null && ab.get(1) != null) {
                    linA[f] = ab.get(0);
                    linb[f] = ab.get(1);
                }
            }
        }

        // see _kb/06-solver-catalog.md for rationale
        Matrix cutoffMat = options.getCutoffMatrix(M, K);
        double[] tokbound = new double[K];
        for (int r = 0; r < K; r++) {
            boolean anyWaitq = false;
            for (int f = 0; f < F; f++) {
                if (iswaitq[f][r]) {
                    anyWaitq = true;
                    break;
                }
            }
            if (!anyWaitq) {
                continue;
            }
            int chainOf = -1;
            if (sn.chains != null && !sn.chains.isEmpty()) {
                for (int c = 0; c < sn.chains.getNumRows(); c++) {
                    if (sn.chains.get(c, r) > 0) {
                        chainOf = c;
                        break;
                    }
                }
            }
            double chainpop;
            if (chainOf >= 0) {
                chainpop = 0;
                for (int r2 = 0; r2 < K; r2++) {
                    if (sn.chains.get(chainOf, r2) > 0) {
                        chainpop += sn.njobs.get(r2);
                    }
                }
            } else {
                chainpop = sn.njobs.get(r);
            }
            if (Double.isFinite(chainpop)) {
                tokbound[r] = chainpop;
            } else {
                double mx = 0;
                for (int i = 0; i < cutoffMat.getNumRows(); i++) {
                    double cv = cutoffMat.get(i, r);
                    if (Double.isFinite(cv)) {
                        mx = Math.max(mx, cv);
                    }
                }
                tokbound[r] = mx;
            }
        }
        final int[] Lmax = new int[F];
        final int[] bufoff = new int[F];
        int widthAcc = nstateful;
        for (int f = 0; f < F; f++) {
            double lsum = 0;
            for (int r = 0; r < K; r++) {
                if (iswaitq[f][r]) {
                    lsum += tokbound[r];
                }
            }
            Lmax[f] = (int) lsum;
            bufoff[f] = widthAcc;
            widthAcc += Lmax[f];
        }
        final int width = widthAcc;

        // initial augmented state (buffers empty, pad -1)
        int[] h0 = new int[nstateful];
        for (int ind = 0; ind < nnodes; ind++) {
            if (sn.isstateful.get(ind, 0) != 1.0) {
                continue;
            }
            int isf = (int) sn.nodeToStateful.get(ind);
            Matrix spc = sn.space.get(sn.stateful.get(isf));
            if (spc.getNumRows() == 1) {
                h0[isf] = 0;
                continue;
            }
            Matrix st = (sn.state != null) ? sn.state.get(sn.stateful.get(isf)) : null;
            if (st == null || st.isEmpty()) {
                throw new RuntimeException("WAITQ finite capacity regions need the initial state (sn.state) to be set.");
            }
            Matrix row = st.getRow(0);
            int w = spc.getNumCols();
            Matrix padded = new Matrix(1, w);
            int shift = w - row.getNumCols();
            for (int j = 0; j < row.getNumCols(); j++) {
                padded.set(0, Math.max(0, shift) + j, row.get(0, j));
            }
            int idx = Matrix.matchrow(spc, padded);
            if (idx < 0) {
                throw new RuntimeException("Initial state of a stateful node not found in its local state space (WAITQ FCR).");
            }
            h0[isf] = idx;
        }
        int[] row0 = new int[width];
        Arrays.fill(row0, -1);
        for (int isf = 0; isf < nstateful; isf++) {
            row0[isf] = h0[isf];
        }

        final BuilderCtx ctx = new BuilderCtx(sn, K, F, nstateful, width, bufoff, Lmax,
                memberMask, ccap, gcap, memcap, szrow, linA, linb);

        for (int f = 0; f < F; f++) {
            if (ctx.violates(f, ctx.regionAggr(h0, f))) {
                throw new RuntimeException("The initial state violates the finite capacity region constraints.");
            }
        }

        ctx.register(row0);

        // BFS over the augmented reachable space
        while (!ctx.frontier.isEmpty()) {
            int s = ctx.frontier.poll();
            int[] row = ctx.SSH.get(s);
            int[] h = Arrays.copyOfRange(row, 0, nstateful);
            List<List<Integer>> bufs = ctx.decodeBufs(row);
            double[][] xf = new double[F][];
            for (int f = 0; f < F; f++) {
                xf[f] = ctx.regionAggr(h, f);
            }
            for (int a = 0; a < A; a++) {
                Sync syncA = sync.get(a);
                int node_a = syncA.active.get(0).getNode();
                if (sn.isstateful.get(node_a, 0) != 1.0) {
                    continue;
                }
                int isf_a = (int) sn.nodeToStateful.get(node_a);
                int class_a = syncA.active.get(0).getJobClass();
                EventType event_a = syncA.active.get(0).getEvent();
                Ret.EventResult resA = State.afterEventHashed(sn, node_a, (double) h[isf_a], event_a, class_a);
                Matrix new_state_a = resA.outspace;
                Matrix rate_a = resA.outrate;
                boolean allInvalid = true;
                for (int c = 0; c < new_state_a.length(); c++) {
                    if (new_state_a.get(c) != -1.0) {
                        allInvalid = false;
                        break;
                    }
                }
                if (allInvalid) {
                    continue;
                }
                int node_p = syncA.passive.get(0).getNode();
                boolean isLocal = (node_p + 1 == local);
                for (int ia = 0; ia < new_state_a.length(); ia++) {
                    double ra = rate_a.get(ia);
                    if (Double.isNaN(ra) || ra <= 0 || new_state_a.get(ia) == -1.0) {
                        continue;
                    }
                    if (isLocal) {
                        int[] newh = h.clone();
                        newh[isf_a] = (int) new_state_a.get(ia);
                        ctx.emit(a, s, newh, bufs, ra, null);
                        continue;
                    }
                    int class_p = syncA.passive.get(0).getJobClass();
                    EventType event_p = syncA.passive.get(0).getEvent();
                    if (sn.isstateful.get(node_p, 0) != 1.0) {
                        continue;
                    }
                    int isf_p = (int) sn.nodeToStateful.get(node_p);
                    int stat_a = (sn.isstation.get(node_a, 0) == 1.0) ? (int) sn.nodeToStation.get(node_a) : -1;
                    int stat_p = (sn.isstation.get(node_p, 0) == 1.0) ? (int) sn.nodeToStation.get(node_p) : -1;
                    int blockedf = -1;
                    int switchf = -1;
                    int droppedf = -1;
                    if (event_p == EventType.ARV && stat_p >= 0) {
                        for (int f = 0; f < F; f++) {
                            if (memberMask[f][stat_p] && (stat_a < 0 || !memberMask[f][stat_a])) {
                                double[] xn = xf[f].clone();
                                xn[class_p] += 1;
                                if (ctx.violates(f, xn)) {
                                    if (!iswaitq[f][class_p]) {
                                        droppedf = f; // DROP: the job is destroyed
                                    } else {
                                        blockedf = f;
                                    }
                                    break;
                                }
                            } else if (memberMask[f][stat_p] && stat_a >= 0 && memberMask[f][stat_a]
                                    && class_p != class_a) {
                                switchf = f;
                                break;
                            }
                        }
                    }
                    double prob_p = syncA.passive.get(0).getProb();
                    if (droppedf >= 0) {
                        // DROP rule (JMT): only the active part applies, job vanishes
                        int[] newh = h.clone();
                        newh[isf_a] = (int) new_state_a.get(ia);
                        ctx.emit(a, s, newh, bufs, ra * prob_p, null);
                        continue;
                    }
                    if (switchf >= 0) {
                        int[] newh = h.clone();
                        newh[isf_a] = (int) new_state_a.get(ia);
                        ctx.emit(a, s, newh, bufs, ra * prob_p,
                                new int[]{switchf, class_p, node_p, iswaitq[switchf][class_p] ? 1 : 0});
                        continue;
                    }
                    if (blockedf >= 0) {
                        if (bufs.get(blockedf).size() >= Lmax[blockedf]) {
                            continue; // FIFO truncation boundary (open-class cutoff)
                        }
                        int[] newh = h.clone();
                        newh[isf_a] = (int) new_state_a.get(ia);
                        List<List<Integer>> newbufs = ctx.copyBufs(bufs);
                        newbufs.get(blockedf).add(node_p * K + class_p);
                        ctx.emit(a, s, newh, newbufs, ra * prob_p, null);
                        continue;
                    }
                    // normal passive application
                    Ret.EventResult resP;
                    if (node_p == node_a) {
                        resP = State.afterEventHashed(sn, node_p, new_state_a.get(ia), event_p, class_p);
                    } else {
                        resP = State.afterEventHashed(sn, node_p, (double) h[isf_p], event_p, class_p);
                    }
                    if (resP == null || resP.outspace == null || resP.outspace.isEmpty()
                            || allInvalidRows(resP.outspace)) {
                        // see _kb/06-solver-catalog.md for rationale
                        if (event_a == EventType.DEP && sn.isbasblocking != null
                                && node_a < sn.isbasblocking.length()
                                && sn.isbasblocking.get(node_a) == 1) {
                            Matrix spaceA = sn.space.get(sn.stateful.get(isf_a));
                            Matrix curVecA = spaceA.getRow(h[isf_a]);
                            int bcolA = curVecA.getNumCols() - 1;
                            if (bcolA >= 0 && curVecA.get(0, bcolA) == 0.0) {
                                Matrix blockedVec = curVecA.copy();
                                blockedVec.set(0, bcolA, 1.0);
                                int blockedIdx = Matrix.matchrow(spaceA, blockedVec);
                                if (blockedIdx >= 0) {
                                    int[] newh = h.clone();
                                    newh[isf_a] = blockedIdx;
                                    ctx.emit(-1, s, newh, bufs, ra * prob_p, null);
                                }
                            }
                        }
                        continue;
                    }
                    Matrix new_state_p = resP.outspace;
                    Matrix outprob_p = resP.outprob;
                    for (int ip = 0; ip < new_state_p.getNumRows(); ip++) {
                        if (new_state_p.get(ip) == -1.0 || ip >= outprob_p.length()) {
                            continue;
                        }
                        double psync = prob_p * outprob_p.get(ip);
                        if (psync <= 0) {
                            continue;
                        }
                        int[] newh = h.clone();
                        newh[isf_a] = (int) new_state_a.get(ia);
                        newh[isf_p] = (int) new_state_p.get(ip);
                        ctx.emit(a, s, newh, bufs, ra * psync, null);
                    }
                }
            }
        }

        // assemble outputs
        int n = ctx.SSH.size();
        MatrixCell Dfilt = new MatrixCell();
        for (int a = 0; a < A; a++) {
            Dfilt.set(a, new Matrix(n, n));
        }
        // Sentinel action -1 collects the true-BAS become-blocked arcs: part of the
        // generator, but not a departure of any action, so kept out of Dfilt.
        Matrix basBlockQ = new Matrix(n, n);
        for (int t = 0; t < ctx.tripA.size(); t++) {
            int a = ctx.tripA.get(t);
            int i = ctx.tripI.get(t);
            int j = ctx.tripJ.get(t);
            double v = ctx.tripV.get(t);
            if (a < 0) {
                basBlockQ.set(i, j, basBlockQ.get(i, j) + v);
            } else {
                Dfilt.get(a).set(i, j, Dfilt.get(a).get(i, j) + v);
            }
        }
        Matrix stateSpaceHashed = new Matrix(n, width);
        for (int s = 0; s < n; s++) {
            int[] r = ctx.SSH.get(s);
            for (int j = 0; j < width; j++) {
                stateSpaceHashed.set(s, j, r[j]);
            }
        }
        int cols = 0;
        int[] nodeW = new int[nstateful];
        for (int isf = 0; isf < nstateful; isf++) {
            nodeW[isf] = sn.space.get(sn.stateful.get(isf)).getNumCols();
            cols += nodeW[isf];
        }
        Matrix stateSpace = new Matrix(n, cols + (width - nstateful));
        Matrix stateSpaceAggr = new Matrix(n, M * K);
        for (int s = 0; s < n; s++) {
            int[] r = ctx.SSH.get(s);
            int pos = 0;
            for (int ind = 0; ind < nnodes; ind++) {
                if (sn.isstateful.get(ind, 0) != 1.0) {
                    continue;
                }
                int isf = (int) sn.nodeToStateful.get(ind);
                Matrix srow = sn.space.get(sn.stateful.get(isf)).getRow(r[isf]);
                for (int j = 0; j < nodeW[isf]; j++) {
                    stateSpace.set(s, pos + j, srow.get(0, j));
                }
                pos += nodeW[isf];
                if (sn.isstation.get(ind, 0) == 1.0) {
                    int ist = (int) sn.nodeToStation.get(ind);
                    State.StateMarginalStatistics ms = ToMarginal.toMarginal(sn, ind, srow, null, null, null, null, null);
                    if (ms != null && ms.nir != null) {
                        for (int rr = 0; rr < K; rr++) {
                            stateSpaceAggr.set(s, ist * K + rr, ms.nir.get(rr));
                        }
                    }
                }
            }
            for (int j = nstateful; j < width; j++) {
                stateSpace.set(s, cols + (j - nstateful), r[j]);
            }
        }
        return new Result(stateSpace, stateSpaceAggr, stateSpaceHashed, Dfilt, basBlockQ);
    }

    /**
     * True when every row index of an outspace is the -1 sentinel, i.e. the event was
     * refused at every branch. Matches the MATLAB isequal(new_state_p,-1) test.
     *
     * @param outspace the outspace returned by State.afterEventHashed
     * @return true when no valid successor row is present
     */
    private static boolean allInvalidRows(Matrix outspace) {
        for (int i = 0; i < outspace.length(); i++) {
            if (outspace.get(i) != -1.0) {
                return false;
            }
        }
        return true;
    }

    /** Mutable BFS state shared by build() and emit(). */
    private static final class BuilderCtx {
        final NetworkStruct sn;
        final int K;
        final int F;
        final int nstateful;
        final int width;
        final int[] bufoff;
        final int[] Lmax;
        final boolean[][] memberMask;
        final double[][] ccap;
        final double[] gcap;
        final double[] memcap;
        final double[][] szrow;
        final Matrix[] linA;
        final Matrix[] linb;

        final Map<String, Integer> keymap = new HashMap<String, Integer>();
        final List<int[]> SSH = new ArrayList<int[]>();
        final Deque<Integer> frontier = new ArrayDeque<Integer>();
        final List<Integer> tripA = new ArrayList<Integer>();
        final List<Integer> tripI = new ArrayList<Integer>();
        final List<Integer> tripJ = new ArrayList<Integer>();
        final List<Double> tripV = new ArrayList<Double>();

        BuilderCtx(NetworkStruct sn, int K, int F, int nstateful, int width, int[] bufoff, int[] Lmax,
                   boolean[][] memberMask, double[][] ccap, double[] gcap, double[] memcap,
                   double[][] szrow, Matrix[] linA, Matrix[] linb) {
            this.sn = sn;
            this.K = K;
            this.F = F;
            this.nstateful = nstateful;
            this.width = width;
            this.bufoff = bufoff;
            this.Lmax = Lmax;
            this.memberMask = memberMask;
            this.ccap = ccap;
            this.gcap = gcap;
            this.memcap = memcap;
            this.szrow = szrow;
            this.linA = linA;
            this.linb = linb;
        }

        boolean violates(int f, double[] x) {
            double tot = 0;
            double mem = 0;
            for (int r = 0; r < K; r++) {
                if (x[r] > ccap[f][r]) {
                    return true;
                }
                tot += x[r];
                mem += x[r] * szrow[f][r];
            }
            if (tot > gcap[f] || mem > memcap[f]) {
                return true;
            }
            if (linA[f] != null && linb[f] != null) {
                int C = linA[f].getNumRows();
                for (int c = 0; c < C; c++) {
                    double lhs = 0;
                    for (int r = 0; r < K; r++) {
                        lhs += linA[f].get(c, r) * x[r];
                    }
                    if (lhs > linb[f].get(c, 0)) {
                        return true;
                    }
                }
            }
            return false;
        }

        double[] regionAggr(int[] h, int f) {
            double[] x = new double[K];
            for (int ist = 0; ist < memberMask[f].length; ist++) {
                if (!memberMask[f][ist]) {
                    continue;
                }
                int ind = (int) sn.stationToNode.get(ist);
                int isf = (int) sn.nodeToStateful.get(ind);
                Matrix srow = sn.space.get(sn.stateful.get(isf)).getRow(h[isf]);
                State.StateMarginalStatistics ms = ToMarginal.toMarginal(sn, ind, srow, null, null, null, null, null);
                if (ms != null && ms.nir != null) {
                    for (int r = 0; r < K; r++) {
                        x[r] += ms.nir.get(r);
                    }
                }
            }
            return x;
        }

        List<List<Integer>> decodeBufs(int[] row) {
            List<List<Integer>> bufs = new ArrayList<List<Integer>>();
            for (int f = 0; f < F; f++) {
                List<Integer> bf = new ArrayList<Integer>();
                for (int j = 0; j < Lmax[f]; j++) {
                    int t = row[bufoff[f] + j];
                    if (t >= 0) {
                        bf.add(t);
                    }
                }
                bufs.add(bf);
            }
            return bufs;
        }

        List<List<Integer>> copyBufs(List<List<Integer>> bufs) {
            List<List<Integer>> out = new ArrayList<List<Integer>>();
            for (int f = 0; f < bufs.size(); f++) {
                out.add(new ArrayList<Integer>(bufs.get(f)));
            }
            return out;
        }

        int register(int[] row) {
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < row.length; j++) {
                sb.append(row[j]).append(',');
            }
            String key = sb.toString();
            Integer idx = keymap.get(key);
            if (idx != null) {
                return idx;
            }
            int newIdx = SSH.size();
            SSH.add(row);
            keymap.put(key, newIdx);
            frontier.add(newIdx);
            return newIdx;
        }

        /**
         * Applies the FIFO release cascade (and an optional pending gated
         * re-entry from a class-switching hop) to the tentative augmented
         * state, then records the transitions of action a.
         */
        void emit(int a, int src, int[] newh, List<List<Integer>> newbufs, double w, int[] pend) {
            if (w <= 0) {
                return;
            }
            Deque<Object[]> work = new ArrayDeque<Object[]>();
            work.add(new Object[]{newh.clone(), copyBufs(newbufs), 1.0, pend});
            while (!work.isEmpty()) {
                Object[] it = work.poll();
                int[] hh = (int[]) it[0];
                @SuppressWarnings("unchecked")
                List<List<Integer>> bb = (List<List<Integer>>) it[1];
                double pw = (Double) it[2];
                int[] pd = (int[]) it[3];
                boolean progressed = false;
                for (int f = 0; f < F; f++) {
                    if (bb.get(f).isEmpty()) {
                        continue;
                    }
                    double[] x = regionAggr(hh, f);
                    int tok = bb.get(f).get(0);
                    int dest = tok / K;
                    int r = tok % K;
                    double[] xn = x.clone();
                    xn[r] += 1;
                    if (violates(f, xn)) {
                        continue; // head-of-line: this region's FIFO stays blocked
                    }
                    int isf_d = (int) sn.nodeToStateful.get(dest);
                    Ret.EventResult resD = State.afterEventHashed(sn, dest, (double) hh[isf_d], EventType.ARV, r);
                    if (resD == null || resD.outspace == null || resD.outspace.isEmpty()) {
                        continue;
                    }
                    Matrix hd = resD.outspace;
                    Matrix opd = resD.outprob;
                    boolean anyBranch = false;
                    for (int id = 0; id < hd.length(); id++) {
                        double op = (id < opd.length()) ? opd.get(id) : 1.0;
                        if (hd.get(id) == -1.0 || op <= 0) {
                            continue;
                        }
                        int[] hh2 = hh.clone();
                        hh2[isf_d] = (int) hd.get(id);
                        List<List<Integer>> bb2 = copyBufs(bb);
                        bb2.get(f).remove(0);
                        work.add(new Object[]{hh2, bb2, pw * op, pd});
                        anyBranch = true;
                    }
                    if (!anyBranch) {
                        continue;
                    }
                    progressed = true;
                    break;
                }
                if (progressed) {
                    continue;
                }
                if (pd != null) {
                    // pending gated re-entry after the cascade settled
                    int f = pd[0];
                    int cls = pd[1];
                    int dest = pd[2];
                    double[] x = regionAggr(hh, f);
                    double[] xn = x.clone();
                    xn[cls] += 1;
                    boolean admitted = false;
                    if (!violates(f, xn)) {
                        int isf_d = (int) sn.nodeToStateful.get(dest);
                        Ret.EventResult resD = State.afterEventHashed(sn, dest, (double) hh[isf_d], EventType.ARV, cls);
                        if (resD != null && resD.outspace != null && !resD.outspace.isEmpty()) {
                            Matrix hd = resD.outspace;
                            Matrix opd = resD.outprob;
                            for (int id = 0; id < hd.length(); id++) {
                                double op = (id < opd.length()) ? opd.get(id) : 1.0;
                                if (hd.get(id) == -1.0 || op <= 0) {
                                    continue;
                                }
                                int[] hh2 = hh.clone();
                                hh2[isf_d] = (int) hd.get(id);
                                work.add(new Object[]{hh2, copyBufs(bb), pw * op, null});
                                admitted = true;
                            }
                        }
                    }
                    if (!admitted) {
                        if (pd.length < 4 || pd[3] != 0) {
                            List<List<Integer>> bb2 = copyBufs(bb);
                            bb2.get(f).add(dest * K + cls);
                            work.add(new Object[]{hh.clone(), bb2, pw, null});
                        } else {
                            // DROP: the switching job is destroyed
                            work.add(new Object[]{hh.clone(), copyBufs(bb), pw, null});
                        }
                    }
                    continue;
                }
                // settled: register the augmented state and the transition
                int[] rowNew = new int[width];
                Arrays.fill(rowNew, -1);
                for (int isf = 0; isf < nstateful; isf++) {
                    rowNew[isf] = hh[isf];
                }
                for (int f = 0; f < F; f++) {
                    for (int j = 0; j < bb.get(f).size(); j++) {
                        rowNew[bufoff[f] + j] = bb.get(f).get(j);
                    }
                }
                int dst = register(rowNew);
                tripA.add(a);
                tripI.add(src);
                tripJ.add(dst);
                tripV.add(w * pw);
            }
        }
    }
}
