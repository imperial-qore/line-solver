package jline.solvers.ssa.handlers;

import java.util.HashMap;
import java.util.Map;

import jline.GlobalConstants;
import jline.lang.NetworkStruct;
import jline.lang.Sync;
import jline.lang.constant.EventType;
import jline.lang.constant.NodeType;
import jline.lang.state.AfterEventContext;
import jline.lang.state.EventCache;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.solvers.ssa.SolverSSA;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Solver_ssa_findenabled {
    private Solver_ssa_findenabled() {}

    public static void solver_ssa_findenabled(NetworkStruct sn,
                                              EventCache eventCache,
                                              int A,
                                              Map<Integer, Integer> node_a,
                                              Map<Integer, Map<Integer, Matrix>> next_state,
                                              Map<Integer, Matrix> stateCell,
                                              Map<Integer, EventType> event_a,
                                              Map<Integer, Integer> class_a,
                                              boolean isSimulation,
                                              Map<Integer, Double> outprob_a,
                                              Map<Integer, Integer> node_p,
                                              int local,
                                              Map<Integer, EventType> event_p,
                                              Map<Integer, Integer> class_p,
                                              Map<Integer, Double> outprob_p,
                                              Map<Integer, Double> prob_sync_p,
                                              Map<Integer, Sync> sync,
                                              Map<Integer, Integer> node_a_sf,
                                              Map<Integer, Integer> node_p_sf,
                                              Map<Integer, Matrix> depRatesSamples,
                                              int samples_collected,
                                              Map<Integer, Matrix> arvRatesSamples,
                                              Matrix csmask,
                                              Map<Integer, Double> enabled_rates,
                                              Map<Integer, Integer> enabled_sync,
                                              Map<Integer, int[]> enabled_fcr,
                                              SolverSSA solverSSA,
                                              AfterEventContext aectx) {
        // see _kb/06-solver-catalog.md for rationale
        int nreg = sn.nregions;
        boolean fcrOn = nreg > 0;
        double[][] fcrClassCap = null, fcrXcur = null;
        double[] fcrGlobalCap = null, fcrMemCap = null;
        boolean[][] fcrMemberMask = null;
        Matrix[] fcrA = null, fcrB = null;
        if (fcrOn) {
            int Kf = sn.nclasses;
            fcrClassCap = new double[nreg][Kf];
            fcrXcur = new double[nreg][Kf];
            fcrGlobalCap = new double[nreg];
            fcrMemCap = new double[nreg];
            fcrMemberMask = new boolean[nreg][];
            fcrA = new Matrix[nreg];
            fcrB = new Matrix[nreg];
            for (int f = 0; f < nreg; f++) {
                Matrix Rmat = sn.region.get(f);
                int M = Rmat.getNumRows();
                boolean[] mask = new boolean[M];
                java.util.List<Integer> members = new java.util.ArrayList<Integer>();
                Matrix memMatMember = (sn.regionmaxmem != null && sn.regionmaxmem.size() > f) ? sn.regionmaxmem.get(f) : null;
                for (int i = 0; i < M; i++) {
                    boolean m = false;
                    for (int c = 0; c <= Kf; c++) {
                        if (Rmat.get(i, c) != -1) { m = true; break; }
                    }
                    // membership: any job-count cap OR the region memory budget set
                    // on the station row (a memory-only region has all caps at -1)
                    if (!m && memMatMember != null && memMatMember.get(i, 0) != -1) { m = true; }
                    mask[i] = m;
                    if (m) { members.add(i); }
                }
                fcrMemberMask[f] = mask;
                for (int r = 0; r < Kf; r++) {
                    double v = Double.POSITIVE_INFINITY;
                    for (int i : members) { double x = Rmat.get(i, r); if (x != -1) { v = Math.min(v, x); } }
                    fcrClassCap[f][r] = v;
                }
                double g = Double.POSITIVE_INFINITY;
                for (int i : members) { double x = Rmat.get(i, Kf); if (x != -1) { g = Math.min(g, x); } }
                fcrGlobalCap[f] = g;
                double mc = Double.POSITIVE_INFINITY;
                Matrix mm = (sn.regionmaxmem != null && sn.regionmaxmem.size() > f) ? sn.regionmaxmem.get(f) : null;
                if (mm != null) {
                    for (int i : members) { double x = mm.get(i, 0); if (x != -1) { mc = Math.min(mc, x); } }
                }
                fcrMemCap[f] = mc;
                if (sn.regionlincon != null && sn.regionlincon.containsKey(f)) {
                    MatrixCell ab = sn.regionlincon.get(f);
                    if (ab != null && ab.size() >= 2 && ab.get(0) != null && ab.get(1) != null) {
                        fcrA[f] = ab.get(0);
                        fcrB[f] = ab.get(1);
                    }
                }
                for (int i : members) {
                    int ind_i = (int) sn.stationToNode.get(i);
                    int isf_i = (int) sn.stationToStateful.get(i);
                    Matrix nirM = ToMarginal.toMarginal(sn, ind_i, stateCell.get(isf_i), null, null, null, null, null).nir;
                    for (int r = 0; r < Kf; r++) { fcrXcur[f][r] += nirM.get(0, r); }
                }
            }
        }
        for (int act = 0; act < A; act++) {
            Map<Integer, Matrix> rate_a = new HashMap<Integer, Matrix>();
            int isf_a = (int) sn.nodeToStateful.get(node_a.get(act));
            int isf_p;
            // next_state[act] = stateCell.mapValues { it.copy() }.toMutableMap()
            Map<Integer, Matrix> initial = new HashMap<Integer, Matrix>();
            for (Map.Entry<Integer, Matrix> entry : stateCell.entrySet()) {
                Matrix v = entry.getValue();
                initial.put(entry.getKey(), v != null ? v.copy() : null);
            }
            next_state.put(act, initial);

            // see _kb/06-solver-catalog.md for rationale
            boolean noPromote = false;
            if (event_a.get(act) == EventType.DEP
                    && node_p.get(act).equals(node_a.get(act))
                    && sn.isstation.get(node_a.get(act)) == 1.0
                    && sn.immfeed != null && !sn.immfeed.isEmpty()) {
                int istIf = (int) sn.nodeToStation.get(node_a.get(act));
                int cpIf = class_p.get(act);
                if (istIf >= 0 && istIf < sn.immfeed.getNumRows()
                        && cpIf >= 0 && cpIf < sn.immfeed.getNumCols()
                        && sn.immfeed.get(istIf, cpIf) > 0.0) {
                    noPromote = true;
                }
            }

            // solverSSA.run { ... } — inline the block (Kotlin scope receiver)
            {
                jline.io.Ret.EventResult eventResult = State.afterEvent(sn,
                        node_a.get(act),
                        stateCell.get(isf_a),
                        event_a.get(act),
                        class_a.get(act),
                        isSimulation,
                        eventCache,
                        aectx,
                        noPromote);
                if (!eventResult.outspace.isEmpty()) {
                    next_state.get(act).put((int) sn.nodeToStateful.get(node_a.get(act)), eventResult.outspace);
                } else {
                    next_state.get(act).remove((int) sn.nodeToStateful.get(node_a.get(act)));
                }
                if (!eventResult.outrate.isEmpty()) {
                    rate_a.put(act, eventResult.outrate);
                } else {
                    rate_a.remove(act);
                }
                if (!eventResult.outprob.isEmpty()) {
                    outprob_a.put(act, eventResult.outprob.toDouble());
                } else {
                    outprob_a.remove(act);
                }
            }

            if (!next_state.get(act).containsKey(isf_a) || !rate_a.containsKey(act)) {
                continue;
            }

            // see _kb/06-solver-catalog.md for rationale
            int numRowsA = next_state.get(act).get(isf_a).getNumRows();
            for (int ia = 0; ia < numRowsA; ia++) {
                if (Double.isNaN(rate_a.get(act).get(ia)) || rate_a.get(act).get(ia) == 0.0) {
                    // handling degenerate rate values
                    rate_a.get(act).set(ia, GlobalConstants.Zero);
                }

                Matrix hash_check = next_state.get(act).get(isf_a);
                boolean hash_found = false;
                for (int col = 0; col < hash_check.getNumCols(); col++) {
                    if (hash_check.get(ia, col) != -1.0) {
                        hash_found = true;
                        break;
                    }
                }

                if (!hash_found) {
                    continue;
                }

                // boolean update_cond = true;
                boolean becomeBlocked = false;   // true BAS: this DEP holds a completed job (not a departure)
                if (rate_a.get(act).get(ia) > 0) {
                    if (!node_p.get(act).equals(local)) {
                        isf_p = (int) sn.nodeToStateful.get(node_p.get(act));
                        if (node_p.get(act).equals(node_a.get(act))) {
                            // self-loop

                            jline.io.Ret.EventResult eventResult = State.afterEvent(sn,
                                    node_p.get(act),
                                    next_state.get(act).get(isf_a),
                                    event_p.get(act),
                                    class_p.get(act),
                                    isSimulation,
                                    eventCache,
                                    aectx);
                            if (!eventResult.outspace.isEmpty()) {
                                next_state.get(act).put(isf_p, eventResult.outspace);
                            } else {
                                next_state.get(act).remove(isf_p);
                            }
                            if (!eventResult.outprob.isEmpty()) {
                                outprob_p.put(act, eventResult.outprob.toDouble());
                            }
                        } else {
                            // departure
                            jline.io.Ret.EventResult eventResult = State.afterEvent(sn,
                                    node_p.get(act),
                                    next_state.get(act).get(isf_p),
                                    event_p.get(act),
                                    class_p.get(act),
                                    isSimulation,
                                    eventCache,
                                    aectx);

                            if (!eventResult.outspace.isEmpty()) {
                                next_state.get(act).put(isf_p, eventResult.outspace);
                            } else {
                                next_state.get(act).remove(isf_p);
                            }
                            if (!eventResult.outprob.isEmpty()) {
                                outprob_p.put(act, eventResult.outprob.toDouble());
                            }
                        }

                        if (next_state.get(act).containsKey(isf_p)) {
                            // Check if the source node (node_a) has state-dependent routing
                            if (node_a.get(act) < sn.nnodes && sn.isstatedep.get(node_a.get(act), 2) == 1.0) {
                                // see _kb/06-solver-catalog.md for rationale
                                Map<jline.lang.nodes.Node, Matrix> stateCell_node = new HashMap<jline.lang.nodes.Node, Matrix>();
                                for (Map.Entry<Integer, Matrix> entry : stateCell.entrySet()) {
                                    Integer stateful_index = entry.getKey();
                                    Matrix matrix = entry.getValue();
                                    if (stateful_index != null && stateful_index < sn.stateful.size()) {
                                        jline.lang.nodes.Node node = sn.stateful.get(stateful_index);
                                        if (node != null) {
                                            stateCell_node.put(node, matrix);
                                        }
                                    }
                                }
                                Map<jline.lang.nodes.Node, Matrix> nextState_node = new HashMap<jline.lang.nodes.Node, Matrix>();
                                for (Map.Entry<Integer, Matrix> entry : next_state.get(act).entrySet()) {
                                    Integer stateful_index = entry.getKey();
                                    Matrix matrix = entry.getValue();
                                    if (stateful_index != null && stateful_index < sn.stateful.size()) {
                                        jline.lang.nodes.Node node = sn.stateful.get(stateful_index);
                                        if (node != null) {
                                            nextState_node.put(node, matrix);
                                        }
                                    }
                                }
                                jline.util.Pair<Map<jline.lang.nodes.Node, Matrix>, Map<jline.lang.nodes.Node, Matrix>> nodePairs =
                                        new jline.util.Pair<Map<jline.lang.nodes.Node, Matrix>, Map<jline.lang.nodes.Node, Matrix>>(
                                                stateCell_node, nextState_node);
                                prob_sync_p.put(act, sync.get(act).passive.get(0).getProb(nodePairs));
                            } else {
                                prob_sync_p.put(act, sync.get(act).passive.get(0).getProb());
                            }
                        } else {
                            prob_sync_p.put(act, 0.0);
                            // see _kb/06-solver-catalog.md for rationale
                            int Rc = sn.nclasses;
                            if (event_a.get(act) == EventType.DEP
                                    && sn.nvars != null && sn.nvars.getNumCols() > 2 * Rc
                                    && sn.nvars.get(node_a.get(act), 2 * Rc) == 1) {
                                Matrix curA = stateCell.get(isf_a);
                                int bcol = curA.getNumCols() - 1;
                                if (curA.get(0, bcol) == 0.0) {
                                    Matrix blockedA = curA.copy();
                                    blockedA.set(0, bcol, 1.0);
                                    next_state.get(act).put(isf_a, blockedA);
                                    next_state.get(act).put(isf_p, stateCell.get(isf_p).copy());
                                    prob_sync_p.put(act, 1.0);
                                    becomeBlocked = true;
                                }
                            }
                        }
                    }
                    if (next_state.get(act).containsKey(isf_a)) {
                        if (node_p.get(act).equals(local)) {
                            prob_sync_p.put(act, 1.0);
                        }
                        if (!Double.isNaN(rate_a.get(act).toDouble())) {
                            if (next_state.get(act).size() == stateCell.size()) {
                                // see _kb/06-solver-catalog.md for rationale
                                boolean blockFCR = false;
                                int[] fcrMark = null;
                                if (fcrOn && !node_p.get(act).equals(local) && node_p.get(act) < sn.nnodes) {
                                    int jp = (int) sn.nodeToStation.get(node_p.get(act));
                                    if (jp >= 0) {
                                        int ja = (int) sn.nodeToStation.get(node_a.get(act));
                                        int cc = class_p.get(act);
                                        int dropIdFcr = jline.lang.constant.DropStrategy.Drop.getID();
                                        for (int f = 0; f < nreg; f++) {
                                            boolean[] mask = fcrMemberMask[f];
                                            boolean waitqCc = (sn.regionrule == null || sn.regionrule.isEmpty())
                                                    || sn.regionrule.get(f, cc) != dropIdFcr;
                                            if (jp < mask.length && mask[jp] && (ja < 0 || ja >= mask.length || !mask[ja])) {
                                                double[] xn = fcrXcur[f].clone();
                                                xn[cc] += 1;
                                                boolean bad = false;
                                                double tot = 0;
                                                for (int r = 0; r < sn.nclasses; r++) {
                                                    tot += xn[r];
                                                    if (xn[r] > fcrClassCap[f][r]) { bad = true; break; }
                                                }
                                                if (!bad && tot > fcrGlobalCap[f]) { bad = true; }
                                                if (!bad && !Double.isInfinite(fcrMemCap[f])) {
                                                    double mem = 0;
                                                    for (int r = 0; r < sn.nclasses; r++) { mem += sn.regionsz.get(f, r) * xn[r]; }
                                                    if (mem > fcrMemCap[f]) { bad = true; }
                                                }
                                                if (!bad && fcrA[f] != null && fcrB[f] != null) {
                                                    int C = fcrA[f].getNumRows();
                                                    for (int q = 0; q < C; q++) {
                                                        double lhs = 0;
                                                        for (int r = 0; r < sn.nclasses; r++) { lhs += fcrA[f].get(q, r) * xn[r]; }
                                                        if (lhs > fcrB[f].get(q, 0)) { bad = true; break; }
                                                    }
                                                }
                                                if (bad) {
                                                    if (waitqCc) {
                                                        fcrMark = new int[]{f, cc, node_p.get(act), 0}; // park
                                                    } else {
                                                        fcrMark = new int[]{f, cc, node_p.get(act), 2}; // DROP: destroyed
                                                    }
                                                    break;
                                                }
                                            } else if (jp < mask.length && mask[jp]
                                                    && ja >= 0 && ja < mask.length && mask[ja]
                                                    && cc != class_a.get(act)) {
                                                // exit + gated re-entry; DROP destroys on refusal
                                                fcrMark = new int[]{f, cc, node_p.get(act), waitqCc ? 1 : 3};
                                                break;
                                            }
                                        }
                                    }
                                }
                                if (fcrMark != null) {
                                    if (node_a.get(act) < sn.nnodes && sn.isstatedep.get(node_a.get(act), 2) == 1.0) {
                                        throw new RuntimeException("WAITQ finite capacity regions are not supported together with state-dependent routing in SolverSSA.");
                                    }
                                    // suppress the passive application: the job leaves the
                                    // upstream node; the entry is resolved at application time
                                    int isf_pk = (int) sn.nodeToStateful.get(node_p.get(act));
                                    next_state.get(act).put(isf_pk, stateCell.get(isf_pk).copy());
                                    prob_sync_p.put(act, sync.get(act).passive.get(0).getProb());
                                }
                                // A true-BAS become-blocked outcome holds the job at the
                                // active station, so it is not a departure either.
                                if (event_a.get(act) == EventType.DEP && !blockFCR && !becomeBlocked) {
                                    isf_p = (int) sn.nodeToStateful.get(node_p.get(act));
                                    node_a_sf.put(act, isf_a);
                                    node_p_sf.put(act, isf_p);

                                    // Matrix original_departure = depRatesSamples.get(class_a.get(act));
                                    // Matrix original_arrival = arvRatesSamples.get(class_p.get(act));
                                    double added_value = (outprob_a.get(act) * outprob_p.get(act)
                                            * rate_a.get(act).get(ia) * prob_sync_p.get(act));

                                    int a_sf_act = node_a_sf.get(act);
                                    int p_sf_act = node_p_sf.get(act);

                                    double dep_value = (depRatesSamples.get(samples_collected - 1)
                                            .get(class_a.get(act), a_sf_act) + added_value);
                                    double arv_val = (arvRatesSamples.get(samples_collected - 1)
                                            .get(class_p.get(act), p_sf_act) + added_value);

                                    depRatesSamples.get(samples_collected - 1)
                                            .set(class_a.get(act), a_sf_act, dep_value);
                                    arvRatesSamples.get(samples_collected - 1)
                                            .set(class_p.get(act), p_sf_act, arv_val);
                                }
                                if (node_p.get(act) < local
                                        && csmask.get(class_a.get(act), class_p.get(act)) != 1.0
                                        && sn.nodetype.get(node_p.get(act)) != NodeType.Source
                                        && (rate_a.get(act).get(ia) * prob_sync_p.get(act) > 0)) {
                                    // Error: state-dependent routing violates the class switching mask
                                    throw new RuntimeException("Error: state-dependent routing at node "
                                            + node_a.get(act) + " violates the class switching mask (node "
                                            + node_a.get(act) + " -> node " + node_p.get(act) + ", class "
                                            + class_a.get(act) + " -> class " + class_p.get(act) + ").");
                                }

                                if (!blockFCR) {
                                    int ctr = enabled_rates.size();
                                    enabled_rates.put(ctr, rate_a.get(act).get(ia) * prob_sync_p.get(act));
                                    enabled_sync.put(ctr, act);
                                    if (fcrMark != null && enabled_fcr != null) {
                                        enabled_fcr.put(ctr, fcrMark);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
