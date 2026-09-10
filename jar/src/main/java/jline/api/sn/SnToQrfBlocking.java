/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.sn;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.lang.NetworkStruct;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.solvers.QrfParams;
import jline.util.matrix.Matrix;

/**
 * Derives the QRF BAS blocking tables (f, MR, BB, MM, ZZ, MM1) from an sn.
 *
 * <p>Port of {@code matlab/src/api/sn/sn_to_qrf_blocking.m}.
 *
 * <p>{@code qrf_bas} describes a Blocking-After-Service network by a
 * finite-capacity queue f and an enumeration of the BLOCKING CONFIGURATIONS
 * reachable behind it. Everything in that enumeration is implied by the model,
 * so it is derived here rather than demanded from the caller;
 * {@code options.qrfParams} remains an explicit override.
 *
 * <p>The tables, and the constraint that reads each one in {@code qrf_bas}:
 * <ul>
 *   <li>{@code f} - the ONE finite-capacity queue. The formulation carries a
 *       scalar f (ZERO4/ZERO7/ZERO8, THM30, THM3I, THM3L all index it), so a
 *       model with two binding buffers is refused here.</li>
 *   <li>{@code F(i)} - min(buffer size, N) for every queue; N where the buffer
 *       is unbounded, since no queue can hold more than the population.</li>
 *   <li>{@code BB(m,i)} - 1 iff queue i is blocked in configuration m.</li>
 *   <li>{@code ZZ(m)} - the blocking depth of configuration m.</li>
 *   <li>{@code MM(m,0)} - head of the FIFO blocking order: the queue that takes
 *       the slot when f completes. {@code qrf_bas} reads ONLY column 0.</li>
 *   <li>{@code MM1(m,j)} - index of the configuration reached from m when j
 *       becomes blocked. Read by THM3L alone, at depth ZM-1.</li>
 * </ul>
 *
 * <p>THREE INVARIANTS, each a correctness condition rather than a convention:
 * <ol>
 *   <li>Configuration 1 MUST be the empty one. ZERO4 iterates {@code m = 2:MR}
 *       and ZERO5/ZERO7/ZERO8 test {@code m >= 2} to mean "some queue is
 *       blocked".</li>
 *   <li>{@code ZM = max(ZZ)} MUST be the reachable maximum. {@code qrf_bas}
 *       recomputes ZM from ZZ and closes the depth ladder there, so a truncated
 *       enumeration excises states the real chain visits and the polytope stops
 *       containing the true distribution -- the bound stops bounding. The size
 *       guard therefore REFUSES; it never truncates.</li>
 *   <li>Blocking APPENDS at the tail: a queue that becomes blocked joins behind
 *       those already waiting, so MM1's successor is the configuration with j
 *       appended, and the head MM(m,0) names never moves.</li>
 * </ol>
 *
 * <p>THE ENUMERATION IS THE FULL ORDERED ONE, and it has to be. A (set, head)
 * collapse looks sound -- the LP reads configurations only through BB, ZZ,
 * MM(:,0) and MM1, and both objective and readout sum over m -- and it would
 * shrink MR from {@code sum_z P(B,z)} to {@code 1 + sum_z C(B,z)*z}. It was
 * tried and it is WRONG. Merging the depth-ZM configurations that share a set
 * and a head makes several THM3L rows, one per depth-(ZM-1) predecessor,
 * reference the SAME merged successor block. That is extra coupling the fine
 * system does not have, so the collapsed polytope is strictly SMALLER, not a
 * projection of the fine one, and it can cut off the true distribution.
 * Measured on a 4-station model with three feeders (B=3, ZM=3, MR 13 collapsed
 * vs 16 full), the collapse reported upper bounds of 0.681/0.979/0.768 where
 * the full enumeration gives 0.709/0.982/0.800: tighter, from a coarser state
 * space, which is the signature of a cut that is not valid.
 *
 * <p>So MR is factorial in the number of feeders B, and the size guard is what
 * keeps that honest: it REFUSES an oversized instance rather than trimming the
 * enumeration, because trimming is the same unsound cut by another name.
 *
 * <p>WHO CAN BE BLOCKED is read from {@code sn.isbasblocking}, not from
 * {@code sn.droprule}: LINE accepts the BAS declaration on the upstream station
 * or on the full destination, and reading droprule at the capped station sees
 * only the second (BUG-83).
 *
 * @since LINE 3.0
 */
public final class SnToQrfBlocking {

    /**
     * Variable-count ceiling of the derived LP. A guard, not a tuning knob: the
     * enumeration cannot be truncated (invariant 2), so an oversized model is
     * refused rather than approximated.
     */
    public static final double DEFAULT_MAXVARS = 5e5;

    private SnToQrfBlocking() {
    }

    /** The derived tables, or the reason they cannot be built. */
    public static final class Result {
        /** Blocking tables, null when msg is non-empty. */
        public final QrfParams params;
        /** Maximum reachable blocking depth. */
        public final int ZM;
        /** 1-based station indices that can be blocked behind f. */
        public final int[] blockers;
        /** Empty on success, otherwise why the tables cannot be derived. */
        public final String msg;

        Result(QrfParams params, int ZM, int[] blockers, String msg) {
            this.params = params;
            this.ZM = ZM;
            this.blockers = blockers;
            this.msg = msg;
        }
    }

    /**
     * @param sn      network structure
     * @param maxVars variable-count ceiling; pass {@link #DEFAULT_MAXVARS} for the default
     * @return the derived tables, or a Result carrying a non-empty msg
     */
    public static Result snToQrfBlocking(NetworkStruct sn, double maxVars) {
        int M = sn.nstations;
        int N = (int) Math.round(sn.njobs == null ? 0.0 : sn.njobs.elementSum());

        SnToQrfCapacity.Result cap = SnToQrfCapacity.snToQrfCapacity(sn);
        if (!cap.msg.isEmpty()) {
            return new Result(null, 0, new int[0], cap.msg);
        }
        int[] F = cap.F;

        List<Integer> fcand = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            if (cap.binding[i]) {
                fcand.add(i);
            }
        }
        if (fcand.isEmpty()) {
            // No binding buffer: callers gate on snHasBlocking first, so this is
            // a defensive branch rather than a normal path.
            return new Result(emptyBlocking(F, 1, M), 0, new int[0], "");
        }
        if (fcand.size() > 1) {
            StringBuilder names = new StringBuilder();
            for (int c = 0; c < fcand.size(); c++) {
                if (c > 0) {
                    names.append(", ");
                }
                names.append(stationName(sn, fcand.get(c)));
            }
            return new Result(null, 0, new int[0], String.format(
                    "'qrf.bas' models a single finite-capacity queue (its f is a scalar), but %d "
                            + "stations have a binding buffer: %s. Use 'qrf.rsrd', whose PBB "
                            + "constraint sums over every full queue and therefore admits several, "
                            + "or cap only one station.", fcand.size(), names.toString()));
        }
        int f = fcand.get(0);

        int[] blockers = qrfBlockers(sn, f, M);

        // Blocking needs f at capacity plus one held job per blocked queue, so
        // the population caps the depth as tightly as the feeder count does.
        int ZM = Math.min(blockers.length, N - F[f]);
        if (ZM < 0) {
            ZM = 0;
        }
        if (ZM == 0) {
            return new Result(emptyBlocking(F, f + 1, M), 0, oneBased(blockers), "");
        }

        List<int[]> cfg = enumeratePermutations(blockers, ZM);
        int MR = cfg.size();

        // Size guard: refuse, never truncate (invariant 2).
        int Ktot = totalPhases(sn, M);
        double nVars = (double) MR * (N + 1) * (N + 1) * Ktot * Ktot + Ktot;
        if (nVars > maxVars) {
            return new Result(null, ZM, oneBased(blockers), String.format(
                    "the QRF BAS linear program for this model would carry %.3g variables (MR=%d "
                            + "blocking configurations, N=%d, %d service phases in total), above "
                            + "the qrfMaxVars limit of %.3g. The enumeration cannot be truncated "
                            + "-- a depth below the reachable maximum ZM=%d excises states the "
                            + "chain visits, and the result would no longer bound. Reduce the "
                            + "population, the number of stations feeding %s, or the phase counts; "
                            + "or raise the limit deliberately.",
                    nVars, MR, N, Ktot, maxVars, ZM, stationName(sn, f)));
        }

        Matrix BB = new Matrix(MR, M);
        BB.zero();
        Matrix MM = new Matrix(MR, Math.max(2, blockers.length));
        MM.zero();
        Matrix MM1 = new Matrix(MR, M);
        MM1.zero();
        int[] ZZ = new int[MR];

        Map<String, Integer> index = new HashMap<String, Integer>();
        for (int m = 0; m < MR; m++) {
            int[] seq = cfg.get(m);
            ZZ[m] = seq.length;
            for (int z = 0; z < seq.length; z++) {
                BB.set(m, seq[z], 1.0);
                // only column 0 is read; the rest records the full order
                MM.set(m, z, seq[z] + 1);
            }
            index.put(cfgKey(seq), m);
        }
        for (int m = 0; m < MR; m++) {
            if (ZZ[m] >= ZM) {
                continue; // THM3L reads MM1 only below ZM
            }
            int[] seq = cfg.get(m);
            for (int b = 0; b < blockers.length; b++) {
                int j = blockers[b];
                if (BB.get(m, j) == 1.0) {
                    continue;
                }
                int[] succ = new int[seq.length + 1];
                System.arraycopy(seq, 0, succ, 0, seq.length);
                succ[seq.length] = j;
                Integer mp = index.get(cfgKey(succ));
                if (mp != null) {
                    MM1.set(m, j, mp + 1);
                }
            }
        }

        QrfParams qp = new QrfParams();
        qp.f = f + 1;
        qp.F = F;
        qp.MR = MR;
        qp.BB = BB;
        qp.MM = MM;
        qp.MM1 = MM1;
        qp.ZZ = ZZ;
        qp.ZM = ZM;
        return new Result(qp, ZM, oneBased(blockers), "");
    }

    /**
     * The one-configuration table for a model in which no blocking state is
     * reachable. {@code qrf_bas} reads it as a plain finite-buffer network.
     */
    private static QrfParams emptyBlocking(int[] F, int fOneBased, int M) {
        QrfParams qp = new QrfParams();
        qp.f = fOneBased;
        qp.F = F;
        qp.MR = 1;
        qp.BB = new Matrix(1, M);
        qp.BB.zero();
        qp.MM = new Matrix(1, 2);
        qp.MM.zero();
        qp.MM1 = new Matrix(1, M);
        qp.MM1.zero();
        qp.ZZ = new int[]{0};
        qp.ZM = 0;
        return qp;
    }

    /**
     * 0-based station indices that hold a completed job when f is full. A
     * blocker must route into f, must not be f, and must not be an infinite
     * server (which has a server per job and cannot be held). BAS itself is read
     * from the BUG-83 field, with a structural fallback for an sn built without
     * the local-variable refresh.
     */
    private static int[] qrfBlockers(NetworkStruct sn, int f, int M) {
        boolean[] declared = new boolean[M];
        boolean any = false;
        if (sn.isbasblocking != null && sn.stationToNode != null) {
            for (int i = 0; i < M && i < sn.stationToNode.getNumRows(); i++) {
                int ind = (int) sn.stationToNode.get(i, 0);
                if (ind >= 0 && ind < sn.isbasblocking.length()
                        && sn.isbasblocking.get(ind) == 1.0) {
                    declared[i] = true;
                    any = true;
                }
            }
        }
        if (!any && sn.droprule != null && sn.stations != null) {
            // Fallback: BAS declared on the upstream station or on the full
            // destination, the same two forms declaresBlockedMarker resolves.
            boolean destBAS = declaresBas(sn, f);
            for (int i = 0; i < M; i++) {
                if (i == f) {
                    continue;
                }
                if (destBAS || declaresBas(sn, i)) {
                    declared[i] = true;
                }
            }
        }

        int R = sn.nclasses;
        List<Integer> out = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            if (i == f || !declared[i]) {
                continue;
            }
            if (sn.sched != null && sn.stations != null && i < sn.stations.size()
                    && sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                continue;
            }
            // Summed over class pairs so the test survives a multiclass sn, even
            // though the QRF gate upstream admits one class only.
            boolean routes = false;
            for (int r = 0; r < R && !routes; r++) {
                for (int s = 0; s < R; s++) {
                    if (sn.rt.get(i * R + r, f * R + s) > 0) {
                        routes = true;
                        break;
                    }
                }
            }
            if (routes) {
                out.add(i);
            }
        }
        Collections.sort(out);
        int[] arr = new int[out.size()];
        for (int i = 0; i < arr.length; i++) {
            arr[i] = out.get(i);
        }
        return arr;
    }

    private static boolean declaresBas(NetworkStruct sn, int ist) {
        if (sn.stations == null || ist >= sn.stations.size()) {
            return false;
        }
        Station st = sn.stations.get(ist);
        Map<jline.lang.JobClass, DropStrategy> row = sn.droprule.get(st);
        if (row == null) {
            return false;
        }
        for (DropStrategy d : row.values()) {
            if (d == DropStrategy.BlockingAfterService) {
                return true;
            }
        }
        return false;
    }

    /**
     * Every ordered sequence of distinct blockers up to length ZM, first entry
     * the head. Deterministic order -- depth ascending, then subsets
     * lexicographic by ascending station index, then the orders of each subset
     * sorted -- so every codebase emits identical tables. The empty
     * configuration is first (invariant 1).
     */
    private static List<int[]> enumeratePermutations(int[] blockers, int ZM) {
        List<int[]> cfg = new ArrayList<int[]>();
        cfg.add(new int[0]);
        int nb = blockers.length;
        for (int z = 1; z <= ZM; z++) {
            List<int[]> subsets = new ArrayList<int[]>();
            combinations(blockers, z, 0, new int[z], 0, subsets);
            for (int s = 0; s < subsets.size(); s++) {
                List<int[]> orders = new ArrayList<int[]>();
                permute(subsets.get(s), 0, orders);
                sortLexicographic(orders);
                cfg.addAll(orders);
            }
        }
        return cfg;
    }

    private static void combinations(int[] src, int k, int start, int[] buf, int depth,
                                     List<int[]> out) {
        if (depth == k) {
            out.add(buf.clone());
            return;
        }
        for (int i = start; i < src.length; i++) {
            buf[depth] = src[i];
            combinations(src, k, i + 1, buf, depth + 1, out);
        }
    }

    private static void permute(int[] arr, int depth, List<int[]> out) {
        if (depth == arr.length) {
            out.add(arr.clone());
            return;
        }
        for (int i = depth; i < arr.length; i++) {
            int t = arr[depth];
            arr[depth] = arr[i];
            arr[i] = t;
            permute(arr, depth + 1, out);
            t = arr[depth];
            arr[depth] = arr[i];
            arr[i] = t;
        }
    }

    private static void sortLexicographic(List<int[]> orders) {
        Collections.sort(orders, new java.util.Comparator<int[]>() {
            @Override
            public int compare(int[] a, int[] b) {
                for (int i = 0; i < Math.min(a.length, b.length); i++) {
                    if (a[i] != b[i]) {
                        return a[i] < b[i] ? -1 : 1;
                    }
                }
                return a.length - b.length;
            }
        });
    }

    /**
     * Identity of a configuration: the whole blocking order, since that is what
     * distinguishes configurations in the enumeration {@code qrf_bas} is
     * entitled to.
     */
    private static String cfgKey(int[] seq) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < seq.length; i++) {
            sb.append(seq[i]).append(',');
        }
        return sb.toString();
    }

    /**
     * Total service phases across stations, which sizes the QRF variable space
     * together with MR and the population.
     */
    private static int totalPhases(NetworkStruct sn, int M) {
        int total = 0;
        for (int i = 0; i < M; i++) {
            int ki = 1;
            if (sn.proc != null && sn.stations != null && i < sn.stations.size()) {
                Map<jline.lang.JobClass, jline.util.matrix.MatrixCell> row =
                        sn.proc.get(sn.stations.get(i));
                if (row != null && !row.isEmpty()) {
                    jline.util.matrix.MatrixCell mc = row.values().iterator().next();
                    if (mc != null && mc.size() > 0 && mc.get(0) != null) {
                        ki = mc.get(0).getNumRows();
                    }
                }
            }
            total += Math.max(1, ki);
        }
        return total;
    }

    private static String stationName(NetworkStruct sn, int ist) {
        if (sn.nodenames != null && sn.stationToNode != null
                && ist < sn.stationToNode.getNumRows()) {
            int ind = (int) sn.stationToNode.get(ist, 0);
            if (ind >= 0 && ind < sn.nodenames.size()) {
                return sn.nodenames.get(ind);
            }
        }
        return "station " + (ist + 1);
    }

    private static int[] oneBased(int[] zeroBased) {
        int[] out = new int[zeroBased.length];
        for (int i = 0; i < out.length; i++) {
            out[i] = zeroBased[i] + 1;
        }
        return out;
    }
}
