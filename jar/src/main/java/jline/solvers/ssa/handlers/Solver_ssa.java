package jline.solvers.ssa.handlers;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.util.FastMath;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.DMatrixSparseTriplet;
import org.ejml.ops.DConvertMatrixStruct;

import jline.VerboseLevel;
import jline.lang.Event;
import jline.lang.GlobalSync;
import jline.lang.ModeEvent;
import jline.lang.NetworkStruct;
import jline.lang.Sync;
import jline.lang.constant.EventType;
import jline.lang.constant.NodeType;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.NodeParam;
import jline.lang.nodes.StatefulNode;
import jline.lang.state.AfterEventContext;
import jline.lang.state.AfterGlobalEvent;
import jline.lang.state.EventCache;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.solvers.SolverOptions;
import jline.solvers.ssa.SSAValues;
import jline.solvers.ssa.SolverSSA;
import jline.streaming.Collector;
import jline.util.Maths;
import jline.util.RandomManager;
import jline.util.Utils;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Solver_ssa {

    private Solver_ssa() {}

    public static SSAValues solver_ssa(NetworkStruct sn_in,
                                       EventCache eventCache,
                                       Map<StatefulNode, Matrix> init_state,
                                       SolverOptions optionsIn,
                                       SolverSSA solverSSA) {
        SolverOptions options = solverSSA.getOptions();
        NetworkStruct sn = sn_in.copy();

        // Set master seed for reproducible SSA simulation
        RandomManager.setMasterSeed(options.seed);

        int nstateful = sn.nstateful;
        int R = sn.nclasses;
        Matrix N = sn.njobs.transpose();
        Map<Integer, Sync> sync = sn.sync;
        Matrix csmask = sn.csmask;

        Matrix cutoff = options.getCutoffMatrix(sn.nstations, sn.nclasses);

        Matrix Np = N.transpose();
        Matrix capacityc = new Matrix(sn.nnodes, sn.nclasses);
        capacityc.zero();

        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.isstation.get(ind) == 1.0) {
                int ist = (int) sn.nodeToStation.get(ind);
                for (int r = 0; r < sn.nclasses; r++) {
                    int c = 0;
                    for (int i = 0; i < sn.chains.getNumRows(); i++) {
                        if (sn.chains.get(i, r) == 1.0) {
                            c = i;
                        }
                    }

                    Matrix proc_m = sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0);
                    boolean disabled = false;
                    for (int row = 0; row < proc_m.getNumRows(); row++) {
                        for (int col = 0; col < proc_m.getNumCols(); col++) {
                            if (Double.isNaN(proc_m.get(row, col)) && sn.nodetype.get(ind) != NodeType.Place) {
                                disabled = true;
                            }
                        }
                    }

                    if (sn.fjclassmap != null && !sn.fjclassmap.isEmpty()
                            && r < sn.fjclassmap.getNumCols() && sn.fjclassmap.get(0, r) >= 0) {
                        // see _kb/06-solver-catalog.md for rationale
                        capacityc.set(ind, r, sn.classcap.get(ist, r));
                    } else if (!sn.visits.get(c).isEmpty() && sn.visits.get(c).get(ist, r) == 0.0) {
                        capacityc.set(ind, r, 0);
                    } else if (!sn.proc.isEmpty() && !sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r))
                            .isEmpty() && disabled) {
                        capacityc.set(ind, r, 0);
                    } else {
                        if (Utils.isInf(N.get(r))) {
                            capacityc.set(ind, r, Maths.min(cutoff.get(ist, r), sn.classcap.get(ist, r)));
                        } else {
                            int njobs_sum = 0;
                            for (int i = 0; i < sn.njobs.getNumCols(); i++) {
                                if (sn.chains.get(c, i) == 1.0) {
                                    njobs_sum += (int) sn.njobs.get(i);
                                }
                            }
                            // closed classes: enumerate up to the chain population, but never
                            // beyond the class capacity at this station (finite-buffer stations)
                            capacityc.set(ind, r, Maths.min(njobs_sum, sn.classcap.get(ist, r)));
                        }
                    }
                }
                // never raise the station capacity above its configured total capacity
                int capacity_sum = (int) Maths.min(capacityc.sumRows().get(ind), sn.cap.get(ist));
                if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.PAS || sn.sched.get(sn.stations.get(ist)) == SchedStrategy.OI) {
                    // see _kb/06-solver-catalog.md for rationale
                    capacity_sum = (int) sn.cap.get(ist);
                }
                if (Utils.isInf(sn.nservers.get(ist))) {
                    sn.nservers.set(ist, (double) capacity_sum);
                }
                for (int col = 0; col < sn.cap.getNumCols(); col++) {
                    sn.cap.set(ist, col, capacity_sum);
                }
                for (int col = 0; col < sn.classcap.getNumCols(); col++) {
                    sn.classcap.set(ist, col, capacityc.get(ind, col));
                }
            }
        }

        // Signal (G-network) classes never occupy a station; cap their per-station
        // capacity at 0 (except at EXT/Source). Mirrors MATLAB solver_ssa.m.
        if (sn.issignal != null && sn.classcap != null) {
            for (int ii = 0; ii < sn.nstations; ii++) {
                if (sn.sched.get(sn.stations.get(ii)) != jline.lang.constant.SchedStrategy.EXT) {
                    for (int r = 0; r < sn.nclasses; r++) {
                        if (sn.issignal.get(r) > 0) sn.classcap.set(ii, r, 0);
                    }
                }
            }
        }
        // Heterogeneous servers (single-class, ORDER policy) -> load-dependent rate
        // mu(n) = sum of the first min(n,c) server rates. Mirrors MATLAB solver_ssa.m.
        for (int ii = 0; ii < sn.nstations; ii++) {
            jline.lang.nodes.Station stt = sn.stations.get(ii);
            jline.lang.nodeparam.ServiceNodeParam snp = sn.getServiceParam(stt);
            if (snp == null || snp.nservertypes <= 0) continue;
            java.util.List<Integer> served = new java.util.ArrayList<Integer>();
            for (int r = 0; r < sn.nclasses; r++) {
                jline.util.matrix.MatrixCell pc = (sn.proc != null && sn.proc.get(stt) != null) ? sn.proc.get(stt).get(sn.jobclasses.get(r)) : null;
                boolean disabled = (pc == null || pc.isEmpty() || pc.get(0).hasNaN());
                if (!disabled && sn.rates.get(ii, r) > 0) served.add(r);
            }
            if (served.size() > 1) {
                throw new RuntimeException("SolverSSA supports heterogeneous servers only for single-class stations. Use SolverJMT or SolverLDES for multi-class heterogeneous servers.");
            }
            if (served.size() == 1) {
                int r = served.get(0);
                java.util.List<Double> srvrates = new java.util.ArrayList<Double>();
                for (int t = 0; t < snp.nservertypes; t++) {
                    boolean compat = snp.servercompat.get(t, r) > 0;
                    Double rate = (snp.heterorates != null && snp.heterorates.get(t) != null) ? snp.heterorates.get(t).get(r) : null;
                    if (compat && rate != null && rate > 0) {
                        int spt = (int) snp.serverspertype.get(t);
                        for (int s = 0; s < spt; s++) srvrates.add(rate);
                    }
                }
                int c = srvrates.size();
                double mu_base = sn.rates.get(ii, r);
                if (c > 0 && mu_base > 0) {
                    int njobsSum = 0;
                    for (int r2 = 0; r2 < sn.nclasses; r2++) { double nj = sn.njobs.get(r2); if (!Double.isInfinite(nj)) njobsSum += (int) nj; }
                    int Lh = Math.max(c, Math.max(njobsSum, 1));
                    if (sn.lldscaling == null || sn.lldscaling.isEmpty()) {
                        sn.lldscaling = new Matrix(sn.nstations, Lh);
                        sn.lldscaling.ones();
                    } else if (sn.lldscaling.getNumCols() < c) {
                        Matrix ext = new Matrix(sn.nstations, c);
                        for (int a = 0; a < sn.nstations; a++)
                            for (int b = 0; b < c; b++)
                                ext.set(a, b, (b < sn.lldscaling.getNumCols()) ? sn.lldscaling.get(a, b) : sn.lldscaling.get(a, sn.lldscaling.getNumCols() - 1));
                        sn.lldscaling = ext;
                    }
                    for (int n = 1; n <= sn.lldscaling.getNumCols(); n++) {
                        int lim = Math.min(n, c);
                        double mun = 0;
                        for (int s = 0; s < lim; s++) mun += srvrates.get(s);
                        sn.lldscaling.set(ii, n - 1, mun / (mu_base * lim));
                    }
                }
            }
        }

        if (Np.hasInfinite()) {
            for (int col = 0; col < Np.getNumCols(); col++) {
                if (Utils.isInf(Np.get(0, col))) {
                    Np.set(0, col, 0);
                }
            }
        }

        // see _kb/06-solver-catalog.md for rationale
        AfterEventContext aectx = State.afterEventInit(sn);

        Matrix init_state_hashed = new Matrix(1, nstateful);
        init_state_hashed.zero();

        Map<Integer, Matrix> arvRatesSamples = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> depRatesSamples = new HashMap<Integer, Matrix>();
        for (int r = 0; r < options.samples; r++) {
            Matrix m = new Matrix(R, nstateful);
            m.zero();
            arvRatesSamples.put(r, m);
            depRatesSamples.put(r, m.copy());
        }
        int A = sync.size();
        int samples_collected = 1;
        Matrix state = init_state_hashed.copy();
        Map<Integer, Matrix> cur_state = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> nir = new HashMap<Integer, Matrix>();

        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.isstateful.get(ind) == 1.0) {
                int isf = (int) sn.nodeToStateful.get(ind);
                // see _kb/06-solver-catalog.md for rationale
                StatefulNode sfKey = sn.isfjaugmented ? sn.stateful.get(isf) : solverSSA.model.getStatefulNodes().get(isf);
                Matrix initRow = init_state.get(sfKey);
                if (initRow == null && sn.isfjaugmented) {
                    for (StatefulNode kk : init_state.keySet()) {
                        if (kk.getName().equals(sfKey.getName())) {
                            initRow = init_state.get(kk);
                            break;
                        }
                    }
                }
                if (initRow == null) {
                    throw new RuntimeException("solver_ssa: no initial state for stateful node '" + sfKey.getName() + "'.");
                }
                Matrix state_space = Matrix.extractRows(initRow,
                        (int) state.get(isf),
                        (int) state.get(isf) + 1,
                        null);
                cur_state.put(isf, state_space);

                if (sn.isstation.get(ind) == 1.0) {
                    int ist = (int) sn.nodeToStation.get(ind);
                    nir.put(ist, ToMarginal.toMarginal(sn, ind, state_space, null, null, null, null, null).nir.transpose());
                }
            }
        }

        state = new Matrix(0, 0);
        Matrix statelen = new Matrix(cur_state.size(), 1);

        for (int ind = 0; ind < cur_state.size(); ind++) {
            if (cur_state.containsKey(ind)) {
                Matrix row = cur_state.get(ind);
                if (state.isEmpty()) {
                    state = row;
                } else {
                    state = Matrix.concatColumns(state, row, null);
                }
                statelen.set(ind, (double) row.getNumElements());
            } else {
                statelen.set(ind, 0.0);
            }
        }

        int nSamples = options.samples - 1;

        double[] tranSyncData = new double[nSamples];

        Matrix z = new Matrix(1, 1);
        z.zero();
        Matrix tranStateInit = Matrix.concatColumns(z, state, null).transpose();
        int tranStateRows = tranStateInit.getNumRows();

        DMatrixRMaj tranStateD = new DMatrixRMaj(tranStateRows, nSamples);
        for (int row = 0; row < tranStateRows; row++) {
            tranStateD.set(row, 0, tranStateInit.get(row, 0));
        }

        samples_collected = 1;

        int SSqRows = 0;
        for (int ist = 0; ist < sn.nstations; ist++) {
            if (nir.containsKey(ist)) {
                SSqRows += nir.get(ist).getNumRows();
            }
        }

        DMatrixRMaj SSqD = new DMatrixRMaj(SSqRows, nSamples);
        int rowOff = 0;
        for (int ist = 0; ist < sn.nstations; ist++) {
            if (nir.containsKey(ist)) {
                Matrix col = nir.get(ist);
                for (int i = 0; i < col.getNumRows(); i++) {
                    SSqD.set(rowOff + i, 0, col.get(i));
                }
                rowOff += col.getNumRows();
            }
        }

        int local = sn.nnodes;

        Map<Integer, Integer> node_a = new HashMap<Integer, Integer>();
        Map<Integer, Integer> node_p = new HashMap<Integer, Integer>();
        Map<Integer, Integer> class_a = new HashMap<Integer, Integer>();
        Map<Integer, Integer> class_p = new HashMap<Integer, Integer>();
        Map<Integer, EventType> event_a = new HashMap<Integer, EventType>();
        Map<Integer, EventType> event_p = new HashMap<Integer, EventType>();
        Map<Integer, Double> outprob_a = new HashMap<Integer, Double>();
        Map<Integer, Double> outprob_p = new HashMap<Integer, Double>();
        for (int act = 0; act < A; act++) {
            Sync s = sync.get(act);
            Event active = s.active.get(0);
            Event passive = s.passive.get(0);
            node_a.put(act, active.getNode());
            node_p.put(act, passive.getNode());
            class_a.put(act, active.getJobClass());
            class_p.put(act, passive.getJobClass());
            event_a.put(act, active.getEvent());
            event_p.put(act, passive.getEvent());
        }

        Map<Integer, Map<Integer, Matrix>> next_state = new HashMap<Integer, Map<Integer, Matrix>>();
        boolean isSimulation = true;
        double cur_time = 0.0;
        Map<Integer, Double> enabled_rates = new LinkedHashMap<Integer, Double>();
        Map<Integer, Integer> enabled_sync = new LinkedHashMap<Integer, Integer>();
        Map<Integer, int[]> enabled_fcr = new LinkedHashMap<Integer, int[]>();
        // FCR WAITQ: per-region FIFO of parked (class, destination) tokens plus
        // the caps used by the release cascade (JMT waiting-queue semantics)
        FcrData[] fcrData = (sn.nregions > 0) ? fcrPrep(sn) : null;
        java.util.List<java.util.List<Integer>> fcrBuf = new java.util.ArrayList<java.util.List<Integer>>();
        for (int f = 0; f < sn.nregions; f++) {
            fcrBuf.add(new java.util.ArrayList<Integer>());
        }
        Map<Integer, Matrix> cur_state_1 = new HashMap<Integer, Matrix>();
        while (samples_collected < options.samples && cur_time <= options.timespan[1]) {
            if (samples_collected == 1) {
                if ("parallel".equals(options.method)) {
                    // pass
                }
            }

            Map<Integer, Integer> node_a_sf = new HashMap<Integer, Integer>();
            Map<Integer, Integer> node_p_sf = new HashMap<Integer, Integer>();
            Map<Integer, Double> prob_sync_p = new HashMap<Integer, Double>();
            enabled_rates.clear();
            enabled_sync.clear();
            enabled_fcr.clear();

            Solver_ssa_findenabled.solver_ssa_findenabled(sn,
                    eventCache,
                    A,
                    node_a,
                    next_state,
                    cur_state,
                    event_a,
                    class_a,
                    isSimulation,
                    outprob_a,
                    node_p,
                    local,
                    event_p,
                    class_p,
                    outprob_p,
                    prob_sync_p,
                    sync,
                    node_a_sf,
                    node_p_sf,
                    depRatesSamples,
                    samples_collected,
                    arvRatesSamples,
                    csmask,
                    enabled_rates,
                    enabled_sync,
                    enabled_fcr,
                    solverSSA,
                    aectx);

            // Handle global synchronizations
            Map<Integer, GlobalSync> gsync = sn.gsync;
            if (gsync != null && !gsync.isEmpty()) {
                int G = gsync.size();
                int ctr = enabled_rates.size();

                for (int gact = 0; gact < G; gact++) {
                    GlobalSync gSync = gsync.get(gact);
                    if (gSync == null) continue;
                    if (gSync.active.isEmpty()) continue;

                    int gind = gSync.active.get(0).getNode();

                    List<Matrix> glspace = new ArrayList<Matrix>();
                    for (int isf = 0; isf < sn.nstateful; isf++) {
                        if (cur_state.get(isf) != null) {
                            glspace.add(cur_state.get(isf).copy());
                        } else {
                            glspace.add(new Matrix(1, 1));
                        }
                    }

                    AfterGlobalEvent.AfterGlobalEventResult result = AfterGlobalEvent.afterGlobalEvent(sn, gind, glspace, gSync, isSimulation);
                    Matrix outrate = result.outrate;
                    Matrix outprob = result.outprob;
                    List<Matrix> outglspace = result.outglspace;

                    for (int ia = 0; ia < outrate.getNumRows(); ia++) {
                        double rate = outrate.get(ia, 0);
                        double prob = (outprob.getNumRows() > ia) ? outprob.get(ia, 0) : 1.0;

                        if (!Double.isNaN(rate) && rate > 0 && !Double.isNaN(prob) && prob > 0) {
                            double combinedRate = rate * prob;
                            enabled_rates.put(ctr, combinedRate);
                            enabled_sync.put(ctr, A + gact);

                            Map<Integer, Matrix> nextStateForGsync = new HashMap<Integer, Matrix>();
                            for (int isf = 0; isf < sn.nstateful; isf++) {
                                if (outglspace.size() > isf && outglspace.get(isf) != null) {
                                    nextStateForGsync.put(isf, outglspace.get(isf).copy());
                                } else if (cur_state.containsKey(isf)) {
                                    nextStateForGsync.put(isf, cur_state.get(isf).copy());
                                }
                            }
                            next_state.put(A + gact, nextStateForGsync);

                            if (gSync.active.get(0).getEvent() == EventType.FIRE) {
                                for (ModeEvent pev : gSync.passive) {
                                    int pevNode = pev.getNode();
                                    if (pevNode >= 0 && pevNode < sn.nnodes &&
                                            !Double.isNaN(sn.nodeToStateful.get(pevNode)) &&
                                            sn.nodeToStateful.get(pevNode) >= 0) {
                                        int pevIsf = (int) sn.nodeToStateful.get(pevNode);
                                        int pevClass = pev.mode % R;

                                        // see _kb/06-solver-catalog.md for rationale
                                        int gsyncRateIdx = samples_collected - 1;
                                        if (pev.getEvent() == EventType.PRE) {
                                            Matrix depMatrix = depRatesSamples.get(gsyncRateIdx);
                                            if (depMatrix != null) {
                                                double currentVal = depMatrix.get(pevClass, pevIsf);
                                                depMatrix.set(pevClass, pevIsf, currentVal + combinedRate);
                                            }
                                        } else if (pev.getEvent() == EventType.POST) {
                                            Matrix arvMatrix = arvRatesSamples.get(gsyncRateIdx);
                                            if (arvMatrix != null) {
                                                double currentVal = arvMatrix.get(pevClass, pevIsf);
                                                arvMatrix.set(pevClass, pevIsf, currentVal + combinedRate);
                                            }
                                        }
                                    }
                                }
                            }

                            ctr++;
                        }
                    }
                }
            }

            // see _kb/06-solver-catalog.md for rationale
            if (sn.fjsync != null && !sn.fjsync.isEmpty()) {
                int G_fj = (sn.gsync != null) ? sn.gsync.size() : 0;
                int FJ = sn.fjsync.size();
                int ctrFj = enabled_rates.size();
                for (int kfj = 0; kfj < FJ; kfj++) {
                    jline.lang.FJSync entry = sn.fjsync.get(kfj);
                    List<Matrix> glspaceFj = new ArrayList<Matrix>();
                    for (int isf = 0; isf < sn.nstateful; isf++) {
                        glspaceFj.add(cur_state.get(isf) != null ? cur_state.get(isf).copy() : new Matrix(1, 1));
                    }
                    jline.lang.state.AfterFJEvent.AfterFJEventResult fjRes =
                            jline.lang.state.AfterFJEvent.afterFJEvent(sn, entry, glspaceFj, true, eventCache);
                    if (fjRes.outGlobalStates.isEmpty()) {
                        continue;
                    }
                    double effFj = fjRes.outrate.get(0, 0) * fjRes.outprob.get(0, 0);
                    if (Double.isNaN(effFj) || effFj <= 0) {
                        continue;
                    }
                    enabled_rates.put(ctrFj, effFj);
                    enabled_sync.put(ctrFj, A + G_fj + kfj);
                    Map<Integer, Matrix> nextStateForFj = new HashMap<Integer, Matrix>();
                    for (int isf = 0; isf < sn.nstateful; isf++) {
                        nextStateForFj.put(isf, fjRes.outGlobalStates.get(0).get(isf).copy());
                    }
                    next_state.put(A + G_fj + kfj, nextStateForFj);
                    // rate accounting: parent departure at the fork, one sibling
                    // arrival per branch head in the tag's auxiliary classes
                    Matrix depMatrixFj = depRatesSamples.get(samples_collected);
                    Matrix arvMatrixFj = arvRatesSamples.get(samples_collected);
                    int isfForkFj = (int) sn.nodeToStateful.get(entry.fork);
                    if (depMatrixFj != null) {
                        depMatrixFj.set(entry.jobclass, isfForkFj, depMatrixFj.get(entry.jobclass, isfForkFj) + effFj);
                    }
                    if (arvMatrixFj != null) {
                        for (int b = 0; b < entry.branchheads.length; b++) {
                            int isfBhFj = (int) sn.nodeToStateful.get(entry.branchheads[b]);
                            arvMatrixFj.set(entry.auxclasses[b], isfBhFj, arvMatrixFj.get(entry.auxclasses[b], isfBhFj) + effFj);
                        }
                    }
                    ctrFj++;
                }
            }

            Matrix enabled_rates_m = new Matrix(1, enabled_rates.size());
            for (int i = 0; i < enabled_rates.size(); i++) {
                enabled_rates_m.set(0, i, enabled_rates.get(i));
            }
            double tot_rate = enabled_rates_m.elementSum();
            Matrix cum_sum = enabled_rates_m.cumsumViaRow();
            Matrix cum_rate = Matrix.scaleMult(cum_sum, 1.0 / tot_rate);

            double rand = Maths.rand();
            int firing_ctr = -1;
            for (int i = 0; i < cum_rate.getNumElements(); i++) {
                if (rand > cum_rate.get(i)) {
                    firing_ctr = i;
                } else {
                    break;
                }
            }
            firing_ctr++;
            if (enabled_sync.isEmpty()) {
                throw new RuntimeException("SSA simulation entered a deadlock before collecting all samples, no synchronization is enabled.");
            }

            tranStateD = update_paddings_dense(sn, cur_state, statelen, tranStateD);
            tranStateRows = tranStateD.getNumRows();

            double dt = -(FastMath.log(Maths.rand()) / tot_rate);
            cur_time += dt;

            save_log_dense(dt,
                    cur_state,
                    tranStateD,
                    samples_collected,
                    tranSyncData,
                    enabled_sync,
                    firing_ctr,
                    sn,
                    nir,
                    cur_state,
                    SSqD,
                    solverSSA.getStreamingCollector(),
                    cur_time);

            cur_state = next_state.get(enabled_sync.get(firing_ctr));

            // see _kb/06-solver-catalog.md for rationale
            if (fcrData != null) {
                int[] mk = enabled_fcr.get(firing_ctr);
                if (mk != null && mk[3] == 0) {
                    fcrBuf.get(mk[0]).add(mk[2] * sn.nclasses + mk[1]);
                }
                // mk[3] == 2: DROP, the refused job was destroyed (active only)
                fcrRelease(sn, fcrData, fcrBuf, cur_state, eventCache, aectx);
                if (mk != null && (mk[3] == 1 || mk[3] == 3)) {
                    // class-switching hop: gate the re-entry after the cascade
                    double[] x = fcrRegionPop(sn, fcrData[mk[0]], cur_state);
                    double[] xn = x.clone();
                    xn[mk[1]] += 1;
                    boolean admitted = false;
                    if (!fcrViolates(sn, fcrData[mk[0]], mk[0], xn)) {
                        int isf_d = (int) sn.nodeToStateful.get(mk[2]);
                        jline.io.Ret.EventResult res = State.afterEvent(sn, mk[2], cur_state.get(isf_d),
                                EventType.ARV, mk[1], true, eventCache, aectx);
                        if (res != null && res.outspace != null && !res.outspace.isEmpty()) {
                            cur_state.put(isf_d, res.outspace);
                            admitted = true;
                        }
                    }
                    if (!admitted && mk[3] == 1) {
                        fcrBuf.get(mk[0]).add(mk[2] * sn.nclasses + mk[1]);
                    }
                    // mk[3] == 3 refused: DROP, the switching job is destroyed
                }
            }

            samples_collected++;
            print_progress(options, samples_collected);
        }

        // Copy final state
        cur_state_1 = new HashMap<Integer, Matrix>();
        for (Map.Entry<Integer, Matrix> e : cur_state.entrySet()) {
            cur_state_1.put(e.getKey(), e.getValue().copy());
        }

        // see _kb/06-solver-catalog.md for rationale
        int firstCol = 0;
        if (options.config != null && options.config.warmupfrac != null && options.config.warmupfrac > 0) {
            double wf = Math.max(0.0, Math.min(0.99, options.config.warmupfrac));
            int nDrop = (int) Math.floor(wf * nSamples);
            if (nDrop > 0 && nDrop < nSamples) {
                firstCol = nDrop;
            }
        }
        int nKeep = nSamples - firstCol;

        // Cumulative sum of dwell times (over the kept samples)
        DMatrixRMaj timesCumSumD = new DMatrixRMaj(nKeep, 1);
        double cumSum = 0.0;
        for (int i = firstCol; i < nSamples; i++) {
            cumSum += tranStateD.get(0, i);
            timesCumSumD.set(i - firstCol, 0, cumSum);
        }
        Matrix timesCumSum = denseToSparseMatrix(timesCumSumD);

        // Find unique state rows
        Map<String, List<Integer>> rowKeyMap = new HashMap<String, List<Integer>>();
        List<String> uniqueKeysOrdered = new ArrayList<String>();
        int estimatedCapacity = (tranStateRows - 1) * 20;

        for (int col = firstCol; col < nSamples; col++) {
            StringBuilder sb = new StringBuilder(estimatedCapacity);
            for (int row = 1; row < tranStateRows; row++) {
                if (row > 1) sb.append(',');
                sb.append(tranStateD.get(row, col));
            }
            String key = sb.toString();
            List<Integer> indices = rowKeyMap.get(key);
            if (indices == null) {
                List<Integer> newList = new ArrayList<Integer>();
                newList.add(col);
                rowKeyMap.put(key, newList);
                uniqueKeysOrdered.add(key);
            } else {
                indices.add(col);
            }
        }

        java.util.Collections.sort(uniqueKeysOrdered);

        int numUniqueStates = uniqueKeysOrdered.size();
        int[] ui = new int[numUniqueStates];
        @SuppressWarnings("unchecked")
        List<Integer>[] uj = new List[numUniqueStates];

        for (int i = 0; i < numUniqueStates; i++) {
            String key = uniqueKeysOrdered.get(i);
            List<Integer> indices = rowKeyMap.get(key);
            ui[i] = indices.get(0);
            uj[i] = indices;
        }

        int stateSizeCount = cur_state_1.size();
        int[] statesz = new int[stateSizeCount];
        for (int i = 0; i < stateSizeCount; i++) {
            statesz[i] = cur_state_1.get(i).getNumElements();
        }

        Map<Integer, Matrix> tranSysState = new HashMap<Integer, Matrix>();
        tranSysState.put(0, timesCumSum);

        int start_index = 1;
        for (int j = 0; j < stateSizeCount; j++) {
            int size = statesz[j];
            int end_index = start_index + size;
            DMatrixRMaj tmpD = new DMatrixRMaj(nKeep, size);
            for (int samp = firstCol; samp < nSamples; samp++) {
                for (int k = start_index; k < end_index; k++) {
                    tmpD.set(samp - firstCol, k - start_index, tranStateD.get(k, samp));
                }
            }
            tranSysState.put(j + 1, denseToSparseMatrix(tmpD));
            start_index = end_index;
        }

        // arvRates and depRates
        DMatrixRMaj[] arvRatesD = new DMatrixRMaj[R];
        DMatrixRMaj[] depRatesD = new DMatrixRMaj[R];
        for (int r = 0; r < R; r++) {
            arvRatesD[r] = new DMatrixRMaj(numUniqueStates, sn.nstateful);
            depRatesD[r] = new DMatrixRMaj(numUniqueStates, sn.nstateful);
        }

        DMatrixRMaj piD = new DMatrixRMaj(1, numUniqueStates);
        for (int s = 0; s < numUniqueStates; s++) {
            List<Integer> stateIndexes = uj[s];
            double dwellSum = 0.0;
            for (Integer idx : stateIndexes) {
                dwellSum += tranStateD.get(0, idx);
            }
            piD.set(0, s, dwellSum);
        }

        DMatrixRMaj SSqNewD = new DMatrixRMaj(numUniqueStates, SSqRows);
        for (int s = 0; s < numUniqueStates; s++) {
            int srcCol = ui[s];
            for (int row = 0; row < SSqRows; row++) {
                SSqNewD.set(s, row, SSqD.get(row, srcCol));
            }
        }
        Matrix SSq = denseToSparseMatrix(SSqNewD);

        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.isstateful.get(ind) == 1.0) {
                int isf = (int) sn.nodeToStateful.get(ind);
                for (int s = 0; s < numUniqueStates; s++) {
                    int uis = ui[s];
                    for (int r = 0; r < R; r++) {
                        arvRatesD[r].set(s, isf, arvRatesSamples.get(uis).get(r, isf));
                        depRatesD[r].set(s, isf, depRatesSamples.get(uis).get(r, isf));
                    }
                }
            }
        }

        Map<Integer, Matrix> arvRates = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> depRates = new HashMap<Integer, Matrix>();
        for (int i = 0; i < R; i++) {
            arvRates.put(i, denseToSparseMatrix(arvRatesD[i]));
            depRates.put(i, denseToSparseMatrix(depRatesD[i]));
        }

        // Normalize pi
        double piSum = 0.0;
        for (int i = 0; i < piD.data.length; i++) piSum += piD.data[i];
        for (int i = 0; i < numUniqueStates; i++) {
            piD.set(0, i, piD.get(0, i) / piSum);
        }
        Matrix pi = denseToSparseMatrix(piD);

        DMatrixRMaj tranSyncRMaj = new DMatrixRMaj(nKeep, 1);
        for (int i = firstCol; i < nSamples; i++) {
            tranSyncRMaj.set(i - firstCol, 0, tranSyncData[i]);
        }
        Matrix tranSync = denseToSparseMatrix(tranSyncRMaj);

        return new SSAValues(pi, SSq, arvRates, depRates, tranSysState, tranSync, sn);
    }

    public static void save_log_dense(double dt,
                                      Map<Integer, Matrix> cur_state,
                                      DMatrixRMaj tranStateD,
                                      int samples_collected,
                                      double[] tranSyncData,
                                      Map<Integer, Integer> enabled_sync,
                                      int firing_ctr,
                                      NetworkStruct sn,
                                      Map<Integer, Matrix> nir,
                                      Map<Integer, Matrix> stateCell,
                                      DMatrixRMaj SSqD,
                                      Collector streamingCollector,
                                      double curTime) {
        int colIdx = samples_collected - 1;
        tranStateD.set(0, colIdx, dt);
        int offset = 1;
        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (cur_state.containsKey(ind)) {
                Matrix row = cur_state.get(ind);
                for (int i = 0; i < row.getNumElements(); i++) {
                    tranStateD.set(offset + i, colIdx, row.get(i));
                }
                offset += row.getNumElements();
            }
        }

        tranSyncData[samples_collected - 1] = (double) (enabled_sync.get(firing_ctr) + 1);

        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.isstation.get(ind) == 1.0) {
                int isf = (int) sn.nodeToStateful.get(ind);
                int ist = (int) sn.nodeToStation.get(ind);
                nir.put(ist, ToMarginal.toMarginal(sn, ind, stateCell.get(isf), null, null, null, null, null).nir.transpose());
            }
        }

        int rowOffset = 0;
        for (int ist = 0; ist < sn.nstations; ist++) {
            if (nir.containsKey(ist)) {
                Matrix col = nir.get(ist);
                for (int i = 0; i < col.getNumRows(); i++) {
                    SSqD.set(rowOffset + i, colIdx, col.get(i));
                }
                rowOffset += col.getNumRows();
            }
        }

        if (streamingCollector != null) {
            Matrix nir_col = new Matrix(0, 0);
            for (int ist = 0; ist < sn.nstations; ist++) {
                if (nir.containsKey(ist)) {
                    Matrix col = nir.get(ist);
                    if (nir_col.isEmpty()) {
                        nir_col = col;
                    } else {
                        nir_col = Matrix.concatRows(nir_col, col, null);
                    }
                }
            }
            streamingCollector.recordState(curTime, dt, nir_col, null, null);
        }
    }

    public static void save_log_dense(double dt,
                                      Map<Integer, Matrix> cur_state,
                                      DMatrixRMaj tranStateD,
                                      int samples_collected,
                                      double[] tranSyncData,
                                      Map<Integer, Integer> enabled_sync,
                                      int firing_ctr,
                                      NetworkStruct sn,
                                      Map<Integer, Matrix> nir,
                                      Map<Integer, Matrix> stateCell,
                                      DMatrixRMaj SSqD) {
        save_log_dense(dt, cur_state, tranStateD, samples_collected, tranSyncData,
                enabled_sync, firing_ctr, sn, nir, stateCell, SSqD, null, 0.0);
    }

    /** Per-region FCR data for the WAITQ release cascade. */
    private static final class FcrData {
        boolean[] mask;
        double[] ccap;
        double gcap;
        double mcap;
        Matrix A;
        Matrix b;
    }

    /** Builds the per-region member masks and caps (memory-only regions included). */
    private static FcrData[] fcrPrep(NetworkStruct sn) {
        int F = sn.nregions;
        int K = sn.nclasses;
        FcrData[] out = new FcrData[F];
        for (int f = 0; f < F; f++) {
            FcrData d = new FcrData();
            Matrix Rmat = sn.region.get(f);
            int M = Rmat.getNumRows();
            Matrix memMat = (sn.regionmaxmem != null && sn.regionmaxmem.size() > f) ? sn.regionmaxmem.get(f) : null;
            d.mask = new boolean[M];
            java.util.List<Integer> members = new java.util.ArrayList<Integer>();
            for (int i = 0; i < M; i++) {
                boolean m = false;
                for (int c = 0; c <= K; c++) {
                    if (Rmat.get(i, c) != -1) { m = true; break; }
                }
                if (!m && memMat != null && memMat.get(i, 0) != -1) { m = true; }
                d.mask[i] = m;
                if (m) { members.add(i); }
            }
            d.ccap = new double[K];
            for (int r = 0; r < K; r++) {
                double v = Double.POSITIVE_INFINITY;
                for (int ii = 0; ii < members.size(); ii++) {
                    double x = Rmat.get(members.get(ii), r);
                    if (x != -1) { v = Math.min(v, x); }
                }
                d.ccap[r] = v;
            }
            d.gcap = Double.POSITIVE_INFINITY;
            d.mcap = Double.POSITIVE_INFINITY;
            for (int ii = 0; ii < members.size(); ii++) {
                double x = Rmat.get(members.get(ii), K);
                if (x != -1) { d.gcap = Math.min(d.gcap, x); }
                if (memMat != null) {
                    double mv = memMat.get(members.get(ii), 0);
                    if (mv != -1) { d.mcap = Math.min(d.mcap, mv); }
                }
            }
            if (sn.regionlincon != null && sn.regionlincon.containsKey(f)) {
                MatrixCell ab = sn.regionlincon.get(f);
                if (ab != null && ab.size() >= 2 && ab.get(0) != null && ab.get(1) != null) {
                    d.A = ab.get(0);
                    d.b = ab.get(1);
                }
            }
            out[f] = d;
        }
        return out;
    }

    /** True if population vector xn breaks any admission constraint of region f. */
    private static boolean fcrViolates(NetworkStruct sn, FcrData d, int f, double[] xn) {
        double tot = 0;
        double mem = 0;
        for (int r = 0; r < sn.nclasses; r++) {
            if (xn[r] > d.ccap[r]) { return true; }
            tot += xn[r];
            mem += ((sn.regionsz != null && !sn.regionsz.isEmpty()) ? sn.regionsz.get(f, r) : 1.0) * xn[r];
        }
        if (tot > d.gcap || mem > d.mcap) { return true; }
        if (d.A != null && d.b != null) {
            int C = d.A.getNumRows();
            for (int c = 0; c < C; c++) {
                double lhs = 0;
                for (int r = 0; r < sn.nclasses; r++) { lhs += d.A.get(c, r) * xn[r]; }
                if (lhs > d.b.get(c, 0)) { return true; }
            }
        }
        return false;
    }

    /** Per-class population of the region under the current state cells. */
    private static double[] fcrRegionPop(NetworkStruct sn, FcrData d, Map<Integer, Matrix> curState) {
        double[] x = new double[sn.nclasses];
        for (int ist = 0; ist < d.mask.length; ist++) {
            if (!d.mask[ist]) { continue; }
            int ind = (int) sn.stationToNode.get(ist);
            int isf = (int) sn.stationToStateful.get(ist);
            Matrix nirM = ToMarginal.toMarginal(sn, ind, curState.get(isf), null, null, null, null, null).nir;
            for (int r = 0; r < sn.nclasses; r++) { x[r] += nirM.get(0, r); }
        }
        return x;
    }

    /**
     * Strict-FIFO head-of-line release of parked region tokens: admits FIFO
     * heads while the admission constraints permit, applying the arrival to
     * the destination station state.
     */
    private static void fcrRelease(NetworkStruct sn, FcrData[] fcr, java.util.List<java.util.List<Integer>> bufs,
                                   Map<Integer, Matrix> curState, EventCache eventCache, AfterEventContext aectx) {
        int K = sn.nclasses;
        boolean progress = true;
        while (progress) {
            progress = false;
            for (int f = 0; f < fcr.length; f++) {
                if (bufs.get(f).isEmpty()) { continue; }
                double[] x = fcrRegionPop(sn, fcr[f], curState);
                int tok = bufs.get(f).get(0);
                int dest = tok / K;
                int r = tok % K;
                double[] xn = x.clone();
                xn[r] += 1;
                if (fcrViolates(sn, fcr[f], f, xn)) {
                    continue; // head-of-line: this region's FIFO stays blocked
                }
                int isf_d = (int) sn.nodeToStateful.get(dest);
                jline.io.Ret.EventResult res = State.afterEvent(sn, dest, curState.get(isf_d),
                        EventType.ARV, r, true, eventCache, aectx);
                if (res == null || res.outspace == null || res.outspace.isEmpty()) {
                    continue; // destination cannot accept (e.g. station capacity)
                }
                curState.put(isf_d, res.outspace);
                bufs.get(f).remove(0);
                progress = true;
            }
        }
    }

    public static void print_progress(SolverOptions options, int samples_collected) {
        if (System.console() != null && !"parallel".equals(options.method) && (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG)) {
            if (samples_collected == 2) {
                System.out.printf("\nSSA samples: %6d ", samples_collected);
                System.out.flush();
            } else if (samples_collected % 100 == 0) {
                System.out.printf("\b\b\b\b\b\b\b %6d", samples_collected);
                System.out.flush();
            }
            if (samples_collected == options.samples) {
                System.out.println();
            }
        }
    }

    private static Matrix denseToSparseMatrix(DMatrixRMaj d) {
        int rows = d.getNumRows();
        int cols = d.getNumCols();
        double[] data = d.data;

        int nnz = 0;
        for (int i = 0; i < data.length; i++) {
            if (data[i] != 0.0) nnz++;
        }

        DMatrixSparseTriplet triplet = new DMatrixSparseTriplet(rows, cols, nnz);
        for (int i = 0; i < rows; i++) {
            int rowOff = i * cols;
            for (int j = 0; j < cols; j++) {
                double v = data[rowOff + j];
                if (v != 0.0) {
                    triplet.addItem(i, j, v);
                }
            }
        }

        DMatrixSparseCSC csc = DConvertMatrixStruct.convert(triplet, (DMatrixSparseCSC) null);
        return new Matrix((org.ejml.data.DMatrix) csc);
    }

    public static DMatrixRMaj update_paddings_dense(NetworkStruct sn,
                                                    Map<Integer, Matrix> stateCell,
                                                    Matrix statelen,
                                                    DMatrixRMaj tranStateD) {
        DMatrixRMaj result = tranStateD;
        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.isstation.get(ind) == 1.0) {
                int isf = (int) sn.nodeToStateful.get(ind);
                boolean deltalen = (stateCell.get(isf).getNumElements() > statelen.get(isf));
                if (deltalen) {
                    statelen.set(isf, (double) stateCell.get(isf).getNumElements());
                    int shift = 0;
                    if (ind > 0) {
                        for (int col = 0; col < isf; col++) {
                            shift += (int) statelen.get(col);
                        }
                    }
                    int oldRows = result.getNumRows();
                    int nCols = result.getNumCols();
                    DMatrixRMaj newResult = new DMatrixRMaj(oldRows + 1, nCols);
                    for (int r = 0; r <= shift; r++) {
                        for (int c = 0; c < nCols; c++) {
                            newResult.set(r, c, result.get(r, c));
                        }
                    }
                    for (int r = shift + 1; r < oldRows; r++) {
                        for (int c = 0; c < nCols; c++) {
                            newResult.set(r + 1, c, result.get(r, c));
                        }
                    }
                    result = newResult;
                }
            }
        }
        // see _kb/06-solver-catalog.md for rationale
        int needed = 1;
        for (int isf = 0; isf < sn.nstateful; isf++) {
            if (stateCell.containsKey(isf)) {
                needed += stateCell.get(isf).getNumElements();
            }
        }
        if (result.getNumRows() < needed) {
            int oldRows = result.getNumRows();
            int nCols = result.getNumCols();
            DMatrixRMaj grown = new DMatrixRMaj(needed, nCols);
            for (int r = 0; r < oldRows; r++) {
                for (int c = 0; c < nCols; c++) {
                    grown.set(r, c, result.get(r, c));
                }
            }
            result = grown;
        }
        return result;
    }
}
