package jline.solvers.ctmc.handlers;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.interfaces.linsol.LinearSolverDense;

import jline.api.mc.Ctmc_makeinfgen;
import jline.api.mc.Ctmc_ssg;
import jline.api.mc.Ctmc_ssg_reachability;
import jline.api.mc.Ctmc_stochcomp;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.EventType;
import jline.lang.constant.NodeType;
import jline.lang.constant.TimingStrategy;
import jline.lang.nodeparam.TransitionNodeParam;
import jline.lang.nodes.Node;
import jline.lang.state.AfterGlobalEvent;
import jline.lang.state.State;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.ResultCTMC;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import jline.util.matrix.MatrixEntry;

public class Solver_ctmc {

    private final SolverCTMC solverCTMC;

    public Solver_ctmc(SolverCTMC solverCTMC) {
        this.solverCTMC = solverCTMC;
    }

    public static ResultCTMC solver_ctmc(NetworkStruct snIn, SolverOptions options) {
        NetworkStruct sn = snIn;
        int nstateful = sn.nstateful;
        int nclasses = sn.nclasses;
        java.util.Map<Integer, jline.lang.Sync> sync = sn.sync;
        int A = sync.size();
        // True when at least one service or arrival process is a matrix exponential, so
        // that the generator legitimately carries negative off-diagonal entries.
        boolean hasMEproc = false;
        if (sn.isph != null) {
            for (java.util.Map<jline.lang.JobClass, Boolean> row : sn.isph.values()) {
                for (Boolean v : row.values()) {
                    if (v != null && !v) {
                        hasMEproc = true;
                    }
                }
            }
        }
        Matrix csmask = sn.csmask;
        // see _kb/06-solver-catalog.md for rationale (incl. the REPLY-signal exception)
        if (sn.issignal != null && sn.classcap != null) {
            for (int ii = 0; ii < sn.nstations; ii++) {
                if (sn.sched.get(sn.stations.get(ii)) != jline.lang.constant.SchedStrategy.EXT) {
                    for (int r = 0; r < nclasses; r++) {
                        boolean isReply = sn.signaltype != null && r < sn.signaltype.size()
                                && sn.signaltype.get(r) == jline.lang.constant.SignalType.REPLY;
                        if (sn.issignal.get(r) > 0 && !isReply) sn.classcap.set(ii, r, 0);
                    }
                }
            }
        }
        // see _kb/06-solver-catalog.md for rationale
        for (int ii = 0; ii < sn.nstations; ii++) {
            jline.lang.nodes.Station stt = sn.stations.get(ii);
            jline.lang.nodeparam.ServiceNodeParam snp = sn.getServiceParam(stt);
            if (snp == null || snp.nservertypes <= 0) continue;
            java.util.List<Integer> served = new java.util.ArrayList<Integer>();
            for (int r = 0; r < nclasses; r++) {
                MatrixCell pc = (sn.proc != null && sn.proc.get(stt) != null) ? sn.proc.get(stt).get(sn.jobclasses.get(r)) : null;
                boolean disabled = (pc == null || pc.isEmpty() || pc.get(0).hasNaN());
                if (!disabled && sn.rates.get(ii, r) > 0) served.add(r);
            }
            if (served.size() > 1) {
                throw new RuntimeException("SolverCTMC supports heterogeneous servers only for single-class stations. Use SolverJMT or SolverLDES for multi-class heterogeneous servers.");
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
                    for (int r2 = 0; r2 < nclasses; r2++) { double nj = sn.njobs.get(r2); if (!Double.isInfinite(nj)) njobsSum += (int) nj; }
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
        if (options.config.state_space_gen == null) {
            options.config.state_space_gen = "default";
        }

        Matrix stateSpace;
        Matrix stateSpaceAggr;
        Matrix stateSpaceHashed;
        if (sn.fjsync != null && !sn.fjsync.isEmpty()) {
            // fork firings break per-chain population conservation, so the
            // population-lattice enumeration cannot generate FJ state spaces
            jline.api.mc.CtmcSsgReachabilityResult ssgResult = jline.api.mc.Ctmc_ssg_fj.ctmc_ssg_fj(sn, options);
            stateSpace = ssgResult.getStateSpace();
            stateSpaceAggr = ssgResult.getStateSpaceAggr();
            stateSpaceHashed = ssgResult.getStateSpaceHashed();
            sn = ssgResult.getSn();
        } else if ("reachable".equals(options.config.state_space_gen)) {
            jline.api.mc.CtmcSsgReachabilityResult ssgResult = Ctmc_ssg_reachability.ctmc_ssg_reachability(sn, options);
            stateSpace = ssgResult.getStateSpace();
            stateSpaceAggr = ssgResult.getStateSpaceAggr();
            stateSpaceHashed = ssgResult.getStateSpaceHashed();
            sn = ssgResult.getSn();
        } else if ("default".equals(options.config.state_space_gen) || "full".equals(options.config.state_space_gen)) {
            jline.solvers.ctmc.SolverCTMC.CtmcSsgResult ssgResult = Ctmc_ssg.ctmc_ssg(sn, options);
            stateSpace = ssgResult.getStateSpace();
            stateSpaceAggr = ssgResult.getStateSpaceAggr();
            stateSpaceHashed = ssgResult.getStateSpaceHashed();
            sn = ssgResult.getSn();
        } else {
            throw new IllegalArgumentException("Unknown state space generation method: " + options.config.state_space_gen);
        }

        // see _kb/06-solver-catalog.md for rationale
        boolean fcrWaitq = sn.nregions > 0;

        MatrixCell Dfilt = new MatrixCell();
        // see _kb/06-solver-catalog.md for rationale
        Matrix basBlockQ = null;
        int local = sn.nnodes + 1;
        if (fcrWaitq) {
            // see _kb/06-solver-catalog.md for rationale
            Solver_ctmc_fcr_waitq.Result wq = Solver_ctmc_fcr_waitq.build(sn, options);
            stateSpace = wq.stateSpace;
            stateSpaceAggr = wq.stateSpaceAggr;
            stateSpaceHashed = wq.stateSpaceHashed;
            Dfilt = wq.Dfilt;
            basBlockQ = wq.basBlockQ;
        } else {
        int sizeGen = stateSpaceHashed.getNumRows();
        basBlockQ = new Matrix(sizeGen, sizeGen);
        for (int a = 0; a < A; a++) {
            Dfilt.set(a, new Matrix(sizeGen, sizeGen));
        }
        for (int a = 0; a < A; a++) {
            MatrixCell stateCell = new MatrixCell();

            for (int s = 0; s < stateSpaceHashed.getNumRows(); s++) {
                Matrix state = stateSpaceHashed.getRow(s);

                for (int ind = 0; ind < sn.nnodes; ind++) {
                    if (sn.isstateful.get(ind, 0) == 1.0) {
                        int isf = (int) sn.nodeToStateful.get(ind);
                        int stateIndex = (int) state.get(isf);

                        Matrix spaceMatrix = sn.space.get(sn.stateful.get(isf));
                        Matrix stateRow = spaceMatrix.getRow(stateIndex);
                        stateCell.set(isf, stateRow);
                    }
                }

                jline.lang.Sync syncA = (jline.lang.Sync) sync.get(a);
                int node_a = syncA.active.get(0).getNode();
                double state_a = state.get((int) sn.nodeToStateful.get(node_a));
                int class_a = syncA.active.get(0).getJobClass();
                EventType event_a = syncA.active.get(0).getEvent();

                Ret.EventResult eventResult = State.afterEventHashed(sn, node_a, state_a, event_a, class_a);

                Matrix new_state_a = eventResult.outspace;
                Matrix rate_a = eventResult.outrate;

                boolean allInvalid = true;
                for (int checkIdx = 0; checkIdx < new_state_a.length(); checkIdx++) {
                    if (new_state_a.get(checkIdx) != -1.0) {
                        allInvalid = false;
                        break;
                    }
                }
                if (allInvalid) {
                    continue;
                }

                for (int ia = 0; ia < new_state_a.length(); ia++) {
                    // A matrix-exponential process embeds in the generator exactly as a
                    // phase-type does, except that the off-diagonal entries of D0 and the
                    // completion vector -A*e may be negative. Those transitions are part of
                    // the balance equations: dropping them leaves the diagonal to absorb
                    // their mass and silently answers a different model (an M/CME/1 lost 3%
                    // of its mean queue length). The stationary vector is then a signed
                    // measure whose aggregates over each phase block are still the exact
                    // probabilities. See sn.isph and _kb/04-networkstruct.md.
                    if (rate_a.get(ia) != 0 && (rate_a.get(ia) > 0 || hasMEproc)) {
                        int node_p = syncA.passive.get(0).getNode();
                        if (node_p + 1 != local) {
                            int state_p = (int) state.get((int) sn.nodeToStateful.get(node_p));
                            int class_p = syncA.passive.get(0).getJobClass();
                            EventType event_p = syncA.passive.get(0).getEvent();

                            Ret.EventResult afterEventResult = null;
                            if (node_p == node_a) {
                                if (new_state_a.get(ia) != -1.0) {
                                    afterEventResult = State.afterEventHashed(sn, node_p, new_state_a.get(ia), event_p, class_p);
                                }
                            } else {
                                if (new_state_a.get(ia) != -1.0) {
                                    afterEventResult = State.afterEventHashed(sn, node_p, (double) state_p, event_p, class_p);
                                }
                            }

                            if (afterEventResult == null) {
                                continue;
                            }

                            Matrix new_state_p = afterEventResult.outspace;
                            Matrix outprob_p = afterEventResult.outprob;

                            for (int ip = 0; ip < new_state_p.getNumRows(); ip++) {
                                double prob_sync_p = 0.0;

                                if (new_state_p.get(ip) == -1.0) {
                                    continue;
                                }

                                if (ip >= outprob_p.length()) {
                                    continue;
                                }

                                double outprob_ip;
                                try {
                                    outprob_ip = outprob_p.get(ip);
                                } catch (Exception e) {
                                    continue;
                                }

                                if (node_p + 1 != local) {
                                    if (new_state_p.get(ip) != -1.0) {
                                        if (sn.isstatedep.get(node_a, 2) != 0.0) {
                                            MatrixCell newStateCell = new MatrixCell(stateCell);
                                            int statefulNodeA = (int) sn.nodeToStateful.get(node_a);
                                            int statefulNodeP = (int) sn.nodeToStateful.get(node_p);

                                            Matrix spaceMatrixA = sn.space.get(sn.stateful.get(statefulNodeA))
                                                    .getRow((int) new_state_a.get(ia));
                                            Matrix spaceMatrixP = sn.space.get(sn.stateful.get(statefulNodeP))
                                                    .getRow((int) new_state_p.get(ip));

                                            newStateCell.set(statefulNodeA, spaceMatrixA);
                                            newStateCell.set(statefulNodeP, spaceMatrixP);

                                            Map<Node, Matrix> stateCell_node = new HashMap<Node, Matrix>();
                                            for (Map.Entry<Integer, Matrix> entry : stateCell.toMap().entrySet()) {
                                                int isf_index = entry.getKey();
                                                Matrix matrix = entry.getValue();
                                                Node node = sn.stateful.get(isf_index);
                                                if (node != null) {
                                                    if (!stateCell_node.containsKey(node)) {
                                                        stateCell_node.put(node, matrix);
                                                    }
                                                }
                                            }

                                            Map<Node, Matrix> newStateCell_node = new HashMap<Node, Matrix>();
                                            for (Map.Entry<Integer, Matrix> entry : newStateCell.toMap().entrySet()) {
                                                int isf_index = entry.getKey();
                                                Matrix matrix = entry.getValue();
                                                Node node = sn.stateful.get(isf_index);
                                                if (node != null) {
                                                    if (!newStateCell_node.containsKey(node)) {
                                                        newStateCell_node.put(node, matrix);
                                                    }
                                                }
                                            }

                                            Pair<Map<Node, Matrix>, Map<Node, Matrix>> nodePairs =
                                                    new Pair<Map<Node, Matrix>, Map<Node, Matrix>>(stateCell_node, newStateCell_node);

                                            prob_sync_p = syncA.passive.get(0).getProb(nodePairs) * outprob_ip;
                                        } else {
                                            prob_sync_p = syncA.passive.get(0).getProb() * outprob_ip;
                                        }
                                    } else {
                                        prob_sync_p = 0.0;
                                    }
                                }

                                Matrix new_state = null;
                                if (!Double.isNaN(new_state_a.get(ia))) {
                                    if (node_p + 1 == local) {
                                        new_state = state.copy();
                                        new_state.set((int) sn.nodeToStateful.get(node_a), new_state_a.get(ia));
                                        prob_sync_p = outprob_p.get(ip);
                                    } else if (!new_state_p.isEmpty()) {
                                        new_state = state.copy();
                                        new_state.set((int) sn.nodeToStateful.get(node_a), new_state_a.get(ia));
                                        new_state.set((int) sn.nodeToStateful.get(node_p), new_state_p.get(ip));
                                    }
                                    java.util.Objects.requireNonNull(new_state);
                                    int ns = Matrix.matchrow(stateSpaceHashed, new_state);

                                    if (ns >= 0) {
                                        if (!rate_a.isEmpty()) {
                                            if (node_p + 1 < local && csmask.get(class_a, class_p) == 0.0
                                                    && rate_a.get(ia) * prob_sync_p > 0
                                                    && sn.nodetype.get(node_p) != NodeType.Source) {
                                                System.err.printf(
                                                        "Error: state-dependent routing at node %d (%s) violates the class switching mask (node %s -> node %s, class %s -> class %s).",
                                                        node_a, sn.nodenames.get(node_a),
                                                        sn.nodenames.get(node_a), sn.nodenames.get(node_p),
                                                        sn.classnames.get(class_a), sn.classnames.get(class_p));
                                            }

                                            double finalRate = rate_a.get(ia) * prob_sync_p;

                                            if (Dfilt.get(a).getNumRows() >= s + 1 && Dfilt.get(a).getNumCols() >= ns + 1) {
                                                Dfilt.get(a).set(s, ns, Dfilt.get(a).get(s, ns) + finalRate);
                                            } else {
                                                Dfilt.get(a).set(s, ns, finalRate);
                                            }
                                        }
                                    }
                                }
                            }
                            // see _kb/06-solver-catalog.md for rationale
                            if (event_a == EventType.DEP && rate_a.get(ia) > 0
                                    && sn.nvars != null && sn.nvars.getNumCols() > 2 * sn.nclasses
                                    && sn.isbasblocking != null && node_a < sn.isbasblocking.length()
                                    && sn.isbasblocking.get(node_a) == 1) {
                                boolean destFull = new_state_p.isEmpty();
                                if (!destFull) {
                                    destFull = true;
                                    for (int ipc = 0; ipc < new_state_p.length(); ipc++) {
                                        if (new_state_p.get(ipc) != -1.0) { destFull = false; break; }
                                    }
                                }
                                int isfA = (int) sn.nodeToStateful.get(node_a);
                                Matrix curVecA = sn.space.get(sn.stateful.get(isfA)).getRow((int) state.get(isfA));
                                int bcolA = curVecA.getNumCols() - 1;
                                if (destFull && curVecA.get(0, bcolA) == 0.0) {
                                    Matrix blockedVec = curVecA.copy();
                                    blockedVec.set(0, bcolA, 1.0);
                                    int blockedIdx = Matrix.matchrow(sn.space.get(sn.stateful.get(isfA)), blockedVec);
                                    if (blockedIdx >= 0) {
                                        Matrix new_state_b = state.copy();
                                        new_state_b.set(isfA, (double) blockedIdx);
                                        int nsb = Matrix.matchrow(stateSpaceHashed, new_state_b);
                                        if (nsb >= 0) {
                                            basBlockQ.set(s, nsb, basBlockQ.get(s, nsb) + rate_a.get(ia));
                                        }
                                    }
                                }
                            }
                        } else {
                            if (!Double.isNaN(new_state_a.get(ia))) {
                                Matrix new_state = state.copy();
                                new_state.set((int) sn.nodeToStateful.get(node_a), new_state_a.get(ia));
                                double prob_sync_p = 1.0;
                                int ns = Matrix.matchrow(stateSpaceHashed, new_state);
                                if (ns > -1) {
                                    if (!rate_a.hasNaN()) {
                                        if (Dfilt.get(a).getNumRows() >= s + 1 && Dfilt.get(a).getNumCols() >= ns + 1) {
                                            Dfilt.get(a).set(s, ns, Dfilt.get(a).get(s, ns) + rate_a.get(ia) * prob_sync_p);
                                        } else {
                                            Dfilt.get(a).set(s, ns, rate_a.get(ia) * prob_sync_p);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        } // if fcrWaitq

        int size = stateSpaceHashed.getNumRows();
        Matrix Q = Matrix.eye(size);
        // see _kb/06-solver-catalog.md for rationale
        Matrix Qimm = new Matrix(size, size);
        // see _kb/06-solver-catalog.md for rationale
        boolean isfjaug = sn.fjsync != null && !sn.fjsync.isEmpty();
        boolean[] immAction = new boolean[A];
        for (int a = 0; a < A; a++) {
            jline.lang.Sync syncImm = (jline.lang.Sync) sync.get(a);
            NodeType nt_a = sn.nodetype.get(syncImm.active.get(0).getNode());
            immAction[a] = (nt_a == NodeType.Router) || (nt_a == NodeType.Fork)
                    || (isfjaug && nt_a == NodeType.Join);
        }
        for (int a = 0; a < A; a++) {
            Q = Q.add(1.0, Dfilt.get(a));
            if (immAction[a]) {
                Qimm = Qimm.add(1.0, Dfilt.get(a));
            }
        }
        // Fold in true-BAS become-blocked transitions (not counted as departures).
        if (basBlockQ != null) {
            Q = Q.add(1.0, basBlockQ);
        }

        // Process global synchronization events
        Map<Integer, jline.lang.GlobalSync> gsyncEvents = sn.gsync;
        int G = (gsyncEvents != null) ? gsyncEvents.size() : 0;
        Matrix[] DfiltGsyncComp = new Matrix[G];
        for (int g = 0; g < G; g++) DfiltGsyncComp[g] = new Matrix(size, size);
        // ENABLE phase moves and firings of a TimingStrategy.IMMEDIATE mode are the
        // two gsync sources emitted at the GlobalConstants.Immediate scale.
        boolean[] immGsync = new boolean[G];

        if (gsyncEvents != null && G > 0) {
            for (int g = 0; g < G; g++) {
                Object globj = gsyncEvents.get(g);
                if (globj == null) continue;
                jline.lang.GlobalSync glevent = (jline.lang.GlobalSync) globj;
                if (glevent.active.isEmpty()) continue;
                int gind = glevent.active.get(0).getNode();
                int isf_transition = (int) sn.nodeToStateful.get(gind);
                Object npObj = sn.nodeparam.get(sn.nodes.get(gind));
                TransitionNodeParam transParam = (npObj instanceof TransitionNodeParam) ? (TransitionNodeParam) npObj : null;
                if (transParam == null) continue;
                int nmodes_g = transParam.nmodes;

                if (glevent.active.get(0).getEvent() == EventType.ENABLE) {
                    immGsync[g] = true;
                } else if (glevent.active.get(0).getEvent() == EventType.FIRE) {
                    int mode_g = glevent.active.get(0).getMode();
                    immGsync[g] = transParam.timing != null && mode_g >= 0
                            && mode_g < transParam.timing.size()
                            && transParam.timing.get(mode_g) == TimingStrategy.IMMEDIATE;
                }

                for (int s = 0; s < stateSpaceHashed.getNumRows(); s++) {
                    Matrix state = stateSpaceHashed.getRow(s);

                    List<Matrix> glspace = new ArrayList<Matrix>(nstateful);
                    for (int isf = 0; isf < nstateful; isf++) {
                        int stateIndex = (int) state.get(isf);
                        Matrix spaceMatrix = sn.space.get(sn.stateful.get(isf));
                        glspace.add(spaceMatrix.getRow(stateIndex));
                    }

                    AfterGlobalEvent.AfterGlobalEventResult result = AfterGlobalEvent.afterGlobalEvent(sn, gind, glspace, glevent, false);
                    Matrix outrate = result.outrate;
                    Matrix outprob = result.outprob;
                    Matrix outcomp = result.outcomp;
                    List<Matrix> outglspace = result.outglspace;

                    if (outrate.isEmpty() || outrate.length() == 0) continue;
                    boolean allZero = true;
                    for (int io = 0; io < outrate.length(); io++) {
                        if (outrate.get(io) != 0.0) { allZero = false; break; }
                    }
                    if (allZero) continue;

                    for (int io = 0; io < outrate.length(); io++) {
                        if (outrate.get(io) == 0.0) continue;

                        Matrix new_state = state.copy();

                        Matrix transOutSpace = outglspace.get(isf_transition);
                        if (transOutSpace.getNumRows() <= io) continue;
                        Matrix trans_state = transOutSpace.getRow(io);
                        int hash_t = Matrix.matchrow(sn.space.get(sn.stateful.get(isf_transition)), trans_state);
                        if (hash_t < 0) continue;
                        new_state.set(isf_transition, (double) hash_t);

                        // see _kb/06-solver-catalog.md for rationale
                        boolean is_comp = false;
                        if (glevent.active.get(0).getEvent() == EventType.FIRE) {
                            if (outcomp.getNumRows() > io) {
                                is_comp = outcomp.get(io, 0) != 0.0;
                            }
                            for (int isf = 0; isf < nstateful; isf++) {
                                if (isf != isf_transition) {
                                    Matrix origPlace = glspace.get(isf);
                                    Matrix newPlace = outglspace.get(isf);
                                    if (!origPlace.isEqualTo(newPlace)) {
                                        int hash_p = Matrix.matchrow(sn.space.get(sn.stateful.get(isf)), newPlace);
                                        if (hash_p < 0) continue;
                                        new_state.set(isf, (double) hash_p);
                                    }
                                }
                            }
                        }

                        int ns = Matrix.matchrow(stateSpaceHashed, new_state);
                        if (ns >= 0) {
                            double prob_val = 1.0;
                            if (!outprob.isEmpty() && io < outprob.length()) {
                                prob_val = outprob.get(io);
                            }
                            double rate_val = outrate.get(io) * prob_val;
                            Q.set(s, ns, Q.get(s, ns) + rate_val);
                            if (immGsync[g]) {
                                Qimm.set(s, ns, Qimm.get(s, ns) + rate_val);
                            }
                            if (is_comp) {
                                DfiltGsyncComp[g].set(s, ns, DfiltGsyncComp[g].get(s, ns) + rate_val);
                            }
                        }
                    }
                }
            }
        }

        // see _kb/06-solver-catalog.md for rationale
        int FJ = (sn.fjsync != null) ? sn.fjsync.size() : 0;
        Matrix[] DfiltFjsync = new Matrix[FJ];
        for (int k = 0; k < FJ; k++) DfiltFjsync[k] = new Matrix(size, size);
        if (FJ > 0) {
            jline.lang.state.EventCache fjEventCache = new jline.lang.state.EventCache(false, false);
            for (int k = 0; k < FJ; k++) {
                jline.lang.FJSync entry = sn.fjsync.get(k);
                for (int s = 0; s < stateSpaceHashed.getNumRows(); s++) {
                    Matrix state = stateSpaceHashed.getRow(s);
                    List<Matrix> glspace = new ArrayList<Matrix>(nstateful);
                    for (int isf = 0; isf < nstateful; isf++) {
                        glspace.add(sn.space.get(sn.stateful.get(isf)).getRow((int) state.get(isf)));
                    }
                    jline.lang.state.AfterFJEvent.AfterFJEventResult fjRes =
                            jline.lang.state.AfterFJEvent.afterFJEvent(sn, entry, glspace, false, fjEventCache);
                    for (int io = 0; io < fjRes.outGlobalStates.size(); io++) {
                        if (fjRes.outprob.get(io, 0) <= 0) continue;
                        Matrix new_state = state.copy();
                        boolean ok = true;
                        List<Matrix> gl_io = fjRes.outGlobalStates.get(io);
                        for (int isf = 0; isf < nstateful; isf++) {
                            Matrix newRow = gl_io.get(isf);
                            if (newRow.isEqualTo(glspace.get(isf))) continue;
                            Matrix spaceIsf = sn.space.get(sn.stateful.get(isf));
                            Matrix prow = newRow;
                            if (prow.getNumCols() < spaceIsf.getNumCols()) {
                                Matrix padded = new Matrix(1, spaceIsf.getNumCols());
                                padded.zero();
                                int shift = spaceIsf.getNumCols() - prow.getNumCols();
                                for (int c = 0; c < prow.getNumCols(); c++) {
                                    padded.set(0, shift + c, prow.get(0, c));
                                }
                                prow = padded;
                            }
                            int hash_p = Matrix.matchrow(spaceIsf, prow);
                            if (hash_p < 0) { ok = false; break; }
                            new_state.set(isf, (double) hash_p);
                        }
                        if (!ok) continue;
                        int ns = Matrix.matchrow(stateSpaceHashed, new_state);
                        if (ns >= 0) {
                            double rate_val = fjRes.outrate.get(io, 0) * fjRes.outprob.get(io, 0);
                            Q.set(s, ns, Q.get(s, ns) + rate_val);
                            Qimm.set(s, ns, Qimm.get(s, ns) + rate_val);
                            DfiltFjsync[k].set(s, ns, DfiltFjsync[k].get(s, ns) + rate_val);
                        }
                    }
                }
            }
        }

        // see _kb/06-solver-catalog.md for rationale
        List<Double> immPurged = null;
        if (options.config.hide_immediate) {
            immPurged = ctmcFindVanishingStates(sn, stateSpaceHashed, nclasses, nstateful, FJ, gsyncEvents, G);
            // see _kb/06-solver-catalog.md for rationale
            List<Double> immRows = new ArrayList<Double>();
            int immGap = 0;
            Matrix qimmRowSums = Qimm.sumRows();
            for (int i = 0; i < immPurged.size(); i++) {
                int r = immPurged.get(i).intValue();
                if (qimmRowSums.get(r, 0) <= 0) {
                    immGap++;
                } else {
                    immRows.add(immPurged.get(i));
                }
            }
            if (immGap > 0) {
                jline.io.InputOutput.line_warning_always("solver_ctmc",
                        "CTMC: %d vanishing state(s) have no immediate outgoing arc; the vanishing predicate and the immediate-arc tagging disagree, so those rows keep their timed arcs.",
                        immGap);
            }
            // see _kb/06-solver-catalog.md for rationale
            Set<Integer> immRowSet = new HashSet<Integer>();
            for (int i = 0; i < immRows.size(); i++) {
                immRowSet.add(Integer.valueOf(immRows.get(i).intValue()));
            }
            zeroRows(Q, immRowSet);
            Iterator<MatrixEntry> qimmIt = Qimm.nonZeroIterator();
            while (qimmIt.hasNext()) {
                MatrixEntry e = qimmIt.next();
                if (immRowSet.contains(Integer.valueOf(e.row))) {
                    Q.set(e.row, e.col, e.value);
                }
            }
            for (int a = 0; a < A; a++) {
                if (immAction[a]) continue;
                zeroRows(Dfilt.get(a), immRowSet);
            }
            for (int g = 0; g < G; g++) {
                if (immGsync[g]) continue;
                zeroRows(DfiltGsyncComp[g], immRowSet);
            }
        }

        Matrix diag_Q = new Matrix(Q);
        Matrix.extractDiag(Q, diag_Q);
        Matrix colMatrix = diag_Q.getColumn(0);
        diag_Q = Matrix.diag(colMatrix.toArray1D());
        Q = Q.sub(diag_Q);

        double[][][] arvRates = new double[stateSpaceHashed.getNumRows()][nstateful][nclasses];
        double[][][] depRates = new double[stateSpaceHashed.getNumRows()][nstateful][nclasses];

        for (int a = 0; a < A; a++) {
            jline.lang.Sync syncA = (jline.lang.Sync) sync.get(a);
            int node_a = syncA.active.get(0).getNode();
            int class_a = syncA.active.get(0).getJobClass();
            EventType event_a = syncA.active.get(0).getEvent();

            int node_p = syncA.passive.get(0).getNode();
            int class_p = syncA.passive.get(0).getJobClass();
            if (event_a == EventType.DEP) {
                int node_a_sf = (int) sn.nodeToStateful.get(node_a);
                int node_p_sf = (int) sn.nodeToStateful.get(node_p);
                Matrix rowSums = Dfilt.get(a).sumRows();
                for (int s = 0; s < stateSpaceHashed.getNumRows(); s++) {
                    double rate = rowSums.get(s, 0);
                    depRates[s][node_a_sf][class_a] += rate;
                    arvRates[s][node_p_sf][class_p] += rate;
                }
            }
        }

        if (gsyncEvents != null && G > 0) {
            for (int g = 0; g < G; g++) {
                Object globj = gsyncEvents.get(g);
                if (globj == null) continue;
                jline.lang.GlobalSync glevent = (jline.lang.GlobalSync) globj;
                if (glevent.active.isEmpty()) continue;
                if (glevent.active.get(0).getEvent() == EventType.FIRE) {
                    int gind_dep = glevent.active.get(0).getNode();
                    for (int j = 0; j < glevent.passive.size(); j++) {
                        jline.lang.ModeEvent pev = glevent.passive.get(j);
                        // see _kb/06-solver-catalog.md for rationale
                        int pev_node = pev.getNode();
                        if (pev_node < 0 || pev_node >= sn.nnodes) continue;
                        if (sn.isstateful.get(pev_node, 0) == 0.0) continue;
                        int pev_isf = (int) sn.nodeToStateful.get(pev_node);
                        int pev_mode = pev.mode;
                        jline.lang.nodeparam.TransitionNodeParam tnp =
                                (jline.lang.nodeparam.TransitionNodeParam) sn.nodeparam.get(sn.nodes.get(gind_dep));
                        Matrix arcs = null;
                        if (pev.getEvent() == EventType.PRE) {
                            arcs = tnp.enabling.get(pev_mode);
                        } else if (pev.getEvent() == EventType.POST) {
                            arcs = tnp.firing.get(pev_mode);
                        }
                        if (arcs == null) continue;
                        Matrix rowSums = DfiltGsyncComp[g].sumRows();
                        for (int pev_class = 0; pev_class < nclasses && pev_class < arcs.getNumCols(); pev_class++) {
                            if (arcs.get(pev_node, pev_class) <= 0) continue;
                            for (int s = 0; s < stateSpaceHashed.getNumRows(); s++) {
                                if (pev.getEvent() == EventType.PRE) {
                                    depRates[s][pev_isf][pev_class] += rowSums.get(s, 0);
                                } else {
                                    arvRates[s][pev_isf][pev_class] += rowSums.get(s, 0);
                                }
                            }
                        }
                    }
                }
            }
        }

        // fork firing rate accounting: parent departure at the fork, one
        // sibling arrival per branch head in the tag's auxiliary classes
        if (FJ > 0) {
            for (int k = 0; k < FJ; k++) {
                jline.lang.FJSync entry = sn.fjsync.get(k);
                int isf_fork = (int) sn.nodeToStateful.get(entry.fork);
                Matrix rowSums = DfiltFjsync[k].sumRows();
                for (int s = 0; s < stateSpaceHashed.getNumRows(); s++) {
                    double rate = rowSums.get(s, 0);
                    if (rate > 0) {
                        depRates[s][isf_fork][entry.jobclass] += rate;
                        for (int b = 0; b < entry.branchheads.length; b++) {
                            int isf_bh = (int) sn.nodeToStateful.get(entry.branchheads[b]);
                            arvRates[s][isf_bh][entry.auxclasses[b]] += rate;
                        }
                    }
                }
            }
        }

        List<Integer> zero_row = Matrix.findIndexWithZeroSum(Q, true);
        List<Integer> zero_col = Matrix.findIndexWithZeroSum(Q, false);

        Q.expandMatrixToSquare();
        Matrix negativeIdentity_row = Matrix.eye(zero_row.size());
        negativeIdentity_row.mulByMinusOne();

        for (int rowIdx = 0; rowIdx < zero_row.size(); rowIdx++) {
            int row = zero_row.get(rowIdx);
            for (int colIdx = 0; colIdx < zero_row.size(); colIdx++) {
                int col = zero_row.get(colIdx);
                Q.set(row, col, negativeIdentity_row.get(rowIdx, colIdx));
            }
        }

        Matrix negativeIdentity_col = Matrix.eye(zero_col.size());
        negativeIdentity_col.mulByMinusOne();

        for (int rowIdx = 0; rowIdx < zero_col.size(); rowIdx++) {
            int row = zero_col.get(rowIdx);
            for (int colIdx = 0; colIdx < zero_col.size(); colIdx++) {
                int col = zero_col.get(colIdx);
                Q.set(row, col, negativeIdentity_col.get(rowIdx, colIdx));
            }
        }

        for (int a = 0; a < A; a++) {
            int colBound = Dfilt.get(a).getNumRows() - 1;
            for (int col = Dfilt.get(a).getNumCols(); col < colBound; col++) {
                for (int row = 0; row < Dfilt.get(a).getNumRows(); row++) {
                    Dfilt.get(a).set(row, col, 0);
                }
            }
        }

        Q = Ctmc_makeinfgen.ctmc_makeinfgen(Q);

        // Drop states unreachable from the initial state.
        // The default state space generator enumerates the whole population lattice, so
        // for a model whose reachable set is constrained -- a Petri net with P-invariants
        // and a synchronous call, whose held-server counters are enumerated independently
        // of the marginals, are the clearest cases -- it also produces states that cannot
        // be reached. Some of those enable nothing at all, or form a closed class of their
        // own, which leaves the generator with several recurrent classes and no unique
        // stationary distribution. Such states carry zero probability by definition, and
        // no reachable state has an arc into them, so restricting every quantity to the
        // reachable set is exact rather than an approximation.
        // This runs before the immediate-state removal below, since the initial state may
        // itself be vanishing and is then absent from the complemented chain.
        // Mirrors MATLAB solver_ctmc.m.
        if (sn.state != null && !sn.state.isEmpty()) {
            boolean allStatesSet = true;
            for (int isf = 0; isf < nstateful; isf++) {
                Matrix row_isf = sn.state.get(sn.stateful.get(isf));
                if (row_isf == null || row_isf.isEmpty()) {
                    allStatesSet = false;
                    break;
                }
            }
            if (allStatesSet) {
                // The per-station initial rows carry only as many buffer slots as the
                // initial population needs, while the enumerated local space is sized for
                // the full capacity. Left-pad each row to its space width (empty buffer
                // slots pad the left, so the server-phase and local-variable tail stays
                // aligned) before matching; without this the lookup fails and the pruning
                // is silently skipped, leaving any enumerated-but-unreachable state to
                // break the stationary solve.
                Matrix initRow = new Matrix(0, 0);
                for (int isf = 0; isf < nstateful; isf++) {
                    Matrix row_isf = Matrix.extractRows(sn.state.get(sn.stateful.get(isf)), 0, 1, null);
                    Matrix space_isf = sn.space.get(sn.stateful.get(isf));
                    int w_isf = space_isf == null ? row_isf.getNumCols() : space_isf.getNumCols();
                    if (row_isf.getNumCols() < w_isf) {
                        Matrix pad = new Matrix(1, w_isf - row_isf.getNumCols());
                        pad.zero();
                        row_isf = pad.concatCols(row_isf);
                    }
                    initRow = initRow.isEmpty() ? row_isf : initRow.concatCols(row_isf);
                }
                int initState = -1;
                if (initRow.getNumCols() == stateSpace.getNumCols()) {
                    initState = Matrix.matchrow(stateSpace, initRow);
                }
                if (initState >= 0) {
                    int nQ = Q.getNumRows();
                    boolean[] reach = new boolean[nQ];
                    reach[initState] = true;
                    List<Integer> frontier = new ArrayList<Integer>();
                    frontier.add(initState);
                    while (!frontier.isEmpty()) {
                        List<Integer> next = new ArrayList<Integer>();
                        for (int fi = 0; fi < frontier.size(); fi++) {
                            int s = frontier.get(fi).intValue();
                            for (int ns = 0; ns < nQ; ns++) {
                                if (ns == s || reach[ns]) continue;
                                // any nonzero off-diagonal is an arc: an ME embeds with negative ones
                                if (Math.abs(Q.get(s, ns)) > 1e-12) {
                                    reach[ns] = true;
                                    next.add(ns);
                                }
                            }
                        }
                        frontier = next;
                    }
                    List<Integer> keep = new ArrayList<Integer>();
                    for (int s = 0; s < nQ; s++) {
                        if (reach[s]) keep.add(s);
                    }
                    if (keep.size() < nQ) {
                        Q = Ctmc_makeinfgen.ctmc_makeinfgen(submatrix(Q, keep, keep));
                        stateSpace = subrows(stateSpace, keep);
                        stateSpaceAggr = subrows(stateSpaceAggr, keep);
                        stateSpaceHashed = subrows(stateSpaceHashed, keep);
                        arvRates = subrates(arvRates, keep, nstateful, nclasses);
                        depRates = subrates(depRates, keep, nstateful, nclasses);
                        for (int a = 0; a < A; a++) {
                            Dfilt.set(a, submatrix(Dfilt.get(a), keep, keep));
                        }
                        for (int k = 0; k < FJ; k++) {
                            if (DfiltFjsync[k] != null) {
                                DfiltFjsync[k] = submatrix(DfiltFjsync[k], keep, keep);
                            }
                        }
                        for (int g = 0; g < G; g++) {
                            if (DfiltGsyncComp[g] != null) {
                                DfiltGsyncComp[g] = submatrix(DfiltGsyncComp[g], keep, keep);
                            }
                        }
                    }
                }
            }
        }

        if (options.config.hide_immediate) {
            List<Double> imm_unique = (immPurged != null) ? immPurged
                    : ctmcFindVanishingStates(sn, stateSpaceHashed, nclasses, nstateful, FJ, gsyncEvents, G);
            Matrix imm = new Matrix(imm_unique.isEmpty() ? 0 : imm_unique.size(), 1);
            for (int i = 0; i < imm_unique.size(); i++) {
                imm.set(i, 0, imm_unique.get(i));
            }

            List<Double> nonimm = new ArrayList<Double>();
            List<Double> allStates = new ArrayList<Double>();
            for (int i = 0; i < Q.getNumRows(); i++) allStates.add((double) i);

            if (!imm_unique.isEmpty()) {
                for (Double state2 : allStates) {
                    if (!imm_unique.contains(state2)) {
                        nonimm.add(state2);
                    }
                }
            } else {
                for (Double state2 : allStates) {
                    nonimm.add(state2);
                }
            }

            if (!imm_unique.isEmpty()) {
                Matrix newStateSpace = new Matrix(nonimm.size(), stateSpace.getNumCols());
                Matrix newStateSpaceAggr = new Matrix(nonimm.size(), stateSpaceAggr.getNumCols());
                Matrix newStateSpaceHashed = new Matrix(nonimm.size(), stateSpaceHashed.getNumCols());

                for (int i = 0; i < nonimm.size(); i++) {
                    int origRow = nonimm.get(i).intValue();
                    for (int col = 0; col < stateSpace.getNumCols(); col++) {
                        newStateSpace.set(i, col, stateSpace.get(origRow, col));
                    }
                    for (int col = 0; col < stateSpaceAggr.getNumCols(); col++) {
                        newStateSpaceAggr.set(i, col, stateSpaceAggr.get(origRow, col));
                    }
                    for (int col = 0; col < stateSpaceHashed.getNumCols(); col++) {
                        newStateSpaceHashed.set(i, col, stateSpaceHashed.get(origRow, col));
                    }
                }

                stateSpace = newStateSpace;
                stateSpaceAggr = newStateSpaceAggr;
                stateSpaceHashed = newStateSpaceHashed;

                jline.solvers.ctmc.SolverCTMC.StochCompResult stochcompResult = Ctmc_stochcomp.ctmc_stochcomp(Q, nonimm);
                Q = stochcompResult.S;
                Matrix Q12 = stochcompResult.Q12;

                Map<Integer, Integer> immSetMap = new HashMap<Integer, Integer>(imm_unique.size() * 2);
                for (int i = 0; i < imm_unique.size(); i++) immSetMap.put(imm_unique.get(i).intValue(), i);
                Map<Integer, Integer> nonimmSet = new HashMap<Integer, Integer>(nonimm.size() * 2);
                for (int i = 0; i < nonimm.size(); i++) nonimmSet.put(nonimm.get(i).intValue(), i);

                @SuppressWarnings("unchecked")
                LinearSolverDense<DMatrixRMaj> denseLU =
                        (stochcompResult.denseLUSolver instanceof LinearSolverDense)
                                ? (LinearSolverDense<DMatrixRMaj>) stochcompResult.denseLUSolver : null;
                int nImm = imm_unique.size();
                int nNonimm = nonimm.size();

                if (FJ > 0 || nImm > 0) {
                    // see _kb/06-solver-catalog.md for rationale
                    double[][][] fjDepRates = new double[nNonimm][nstateful][nclasses];
                    double[][][] fjArvRates = new double[nNonimm][nstateful][nclasses];
                    java.util.function.Function<Matrix, double[]> actionRate = (Matrix D) -> {
                        Matrix rs = D.sumRows();
                        Matrix bImm = new Matrix(nImm, 1);
                        for (int i = 0; i < nImm; i++) {
                            bImm.set(i, 0, rs.get(imm_unique.get(i).intValue(), 0));
                        }
                        Matrix xImm = new Matrix(nImm, 1);
                        if (bImm.elementSum() != 0.0) {
                            if (denseLU != null) {
                                DMatrixRMaj denseB = new DMatrixRMaj(nImm, 1);
                                for (int i = 0; i < nImm; i++) denseB.set(i, 0, bImm.get(i, 0));
                                DMatrixRMaj denseX = new DMatrixRMaj(nImm, 1);
                                denseLU.solve(denseB, denseX);
                                for (int i = 0; i < nImm; i++) xImm.set(i, 0, denseX.get(i, 0));
                            } else {
                                Matrix.solve(stochcompResult.Q22.neg(), bImm, xImm);
                            }
                        }
                        Matrix chain = Q12.mult(xImm);
                        double[] out = new double[nNonimm];
                        for (int i = 0; i < nNonimm; i++) {
                            out[i] = rs.get(nonimm.get(i).intValue(), 0) + chain.get(i, 0);
                        }
                        return out;
                    };
                    for (int a = 0; a < A; a++) {
                        jline.lang.Sync syncA2 = (jline.lang.Sync) sync.get(a);
                        if (syncA2.active.get(0).getEvent() != EventType.DEP) continue;
                        int nodeA2 = syncA2.active.get(0).getNode();
                        if (sn.isstateful.get(nodeA2, 0) == 0.0) continue;
                        // see _kb/06-solver-catalog.md for rationale
                        double[] rA;
                        if (immAction[a] && isfjaug && sn.nodetype.get(nodeA2) == NodeType.Join) {
                            rA = actionRate.apply(Dfilt.get(a));
                        } else {
                            Matrix rowSumsA2 = Dfilt.get(a).sumRows();
                            rA = new double[nNonimm];
                            for (int i = 0; i < nNonimm; i++) {
                                rA[i] = rowSumsA2.get(nonimm.get(i).intValue(), 0);
                            }
                        }
                        int isfA2 = (int) sn.nodeToStateful.get(nodeA2);
                        int classA2 = syncA2.active.get(0).getJobClass();
                        int nodeP2 = syncA2.passive.get(0).getNode();
                        for (int i = 0; i < nNonimm; i++) {
                            fjDepRates[i][isfA2][classA2] += rA[i];
                        }
                        if (nodeP2 < sn.nnodes && sn.isstateful.get(nodeP2, 0) == 1.0) {
                            int isfP2 = (int) sn.nodeToStateful.get(nodeP2);
                            int classP2 = syncA2.passive.get(0).getJobClass();
                            for (int i = 0; i < nNonimm; i++) {
                                fjArvRates[i][isfP2][classP2] += rA[i];
                            }
                        }
                    }
                    for (int k = 0; k < FJ; k++) {
                        jline.lang.FJSync entry = sn.fjsync.get(k);
                        double[] rK = actionRate.apply(DfiltFjsync[k]);
                        int isfFork2 = (int) sn.nodeToStateful.get(entry.fork);
                        for (int i = 0; i < nNonimm; i++) {
                            fjDepRates[i][isfFork2][entry.jobclass] += rK[i];
                        }
                        for (int b = 0; b < entry.branchheads.length; b++) {
                            int isfBh2 = (int) sn.nodeToStateful.get(entry.branchheads[b]);
                            for (int i = 0; i < nNonimm; i++) {
                                fjArvRates[i][isfBh2][entry.auxclasses[b]] += rK[i];
                            }
                        }
                    }
                    // see _kb/06-solver-catalog.md for rationale
                    if (gsyncEvents != null && G > 0) {
                        for (int g_rc = 0; g_rc < G; g_rc++) {
                            Object gObj = gsyncEvents.get(g_rc);
                            if (gObj == null) continue;
                            jline.lang.GlobalSync glevent_rc = (jline.lang.GlobalSync) gObj;
                            if (glevent_rc.active.isEmpty()) continue;
                            if (glevent_rc.active.get(0).getEvent() != EventType.FIRE) continue;
                            int gind_rc = glevent_rc.active.get(0).getNode();
                            double[] rG = actionRate.apply(DfiltGsyncComp[g_rc]);
                            for (int j = 0; j < glevent_rc.passive.size(); j++) {
                                jline.lang.ModeEvent pev = (jline.lang.ModeEvent) glevent_rc.passive.get(j);
                                int pevNode = pev.getNode();
                                if (pevNode < 0 || pevNode >= sn.nnodes) continue;
                                if (sn.isstateful.get(pevNode, 0) == 0.0) continue;
                                int pevIsf = (int) sn.nodeToStateful.get(pevNode);
                                Object tnpObj = sn.nodeparam.get(sn.nodes.get(gind_rc));
                                if (!(tnpObj instanceof jline.lang.nodeparam.TransitionNodeParam)) continue;
                                jline.lang.nodeparam.TransitionNodeParam tnp2 =
                                        (jline.lang.nodeparam.TransitionNodeParam) tnpObj;
                                // see _kb/06-solver-catalog.md for rationale
                                Matrix arcs2 = null;
                                if (pev.getEvent() == EventType.PRE) {
                                    arcs2 = tnp2.enabling.get(pev.mode);
                                } else if (pev.getEvent() == EventType.POST) {
                                    arcs2 = tnp2.firing.get(pev.mode);
                                }
                                if (arcs2 == null) continue;
                                for (int pevClass = 0; pevClass < nclasses && pevClass < arcs2.getNumCols(); pevClass++) {
                                    if (arcs2.get(pevNode, pevClass) <= 0) continue;
                                    if (pev.getEvent() == EventType.PRE) {
                                        for (int i = 0; i < nNonimm; i++) fjDepRates[i][pevIsf][pevClass] += rG[i];
                                    } else {
                                        for (int i = 0; i < nNonimm; i++) fjArvRates[i][pevIsf][pevClass] += rG[i];
                                    }
                                }
                            }
                        }
                    }
                    depRates = fjDepRates;
                    arvRates = fjArvRates;
                }

                for (int a = 0; a < A; a++) {
                    Matrix dfiltA = Dfilt.get(a);
                    DMatrixSparseCSC sparseD = dfiltA.toDMatrixSparseCSC();
                    Matrix Q21a = new Matrix(nImm, nNonimm);
                    for (int c = 0; c < sparseD.getNumCols(); c++) {
                        Integer newCol = nonimmSet.get(c);
                        if (newCol == null) continue;
                        int idx0 = sparseD.col_idx[c];
                        int idx1 = sparseD.col_idx[c + 1];
                        for (int idx = idx0; idx < idx1; idx++) {
                            int r = sparseD.nz_rows[idx];
                            Integer newRow = immSetMap.get(r);
                            if (newRow == null) continue;
                            Q21a.set(newRow, newCol, sparseD.nz_values[idx]);
                        }
                    }

                    Matrix T_intermediate;
                    if (denseLU != null) {
                        DMatrixSparseCSC sparseQ21a = Q21a.toDMatrixSparseCSC();
                        List<Integer> nzCols = new ArrayList<Integer>();
                        for (int c = 0; c < sparseQ21a.getNumCols(); c++) {
                            if (sparseQ21a.col_idx[c + 1] > sparseQ21a.col_idx[c]) {
                                nzCols.add(c);
                            }
                        }

                        if (nzCols.isEmpty()) {
                            T_intermediate = new Matrix(nImm, nNonimm);
                        } else {
                            int k = nzCols.size();
                            DMatrixRMaj denseB = new DMatrixRMaj(nImm, k);
                            for (int ci = 0; ci < k; ci++) {
                                int c = nzCols.get(ci);
                                int ci0 = sparseQ21a.col_idx[c];
                                int ci1 = sparseQ21a.col_idx[c + 1];
                                for (int idx = ci0; idx < ci1; idx++) {
                                    denseB.set(sparseQ21a.nz_rows[idx], ci, sparseQ21a.nz_values[idx]);
                                }
                            }
                            DMatrixRMaj denseX = new DMatrixRMaj(nImm, k);
                            denseLU.solve(denseB, denseX);

                            T_intermediate = new Matrix(nImm, nNonimm);
                            for (int ci = 0; ci < k; ci++) {
                                int c = nzCols.get(ci);
                                for (int r = 0; r < nImm; r++) {
                                    double v = denseX.get(r, ci);
                                    if (Math.abs(v) > 1e-15) {
                                        T_intermediate.set(r, c, v);
                                    }
                                }
                            }
                        }
                    } else {
                        T_intermediate = new Matrix(nImm, nNonimm);
                        Matrix.solve(stochcompResult.Q22.neg(), Q21a, T_intermediate);
                    }

                    Matrix Ta = Q12.mult(T_intermediate);

                    Matrix dfilt_value = new Matrix(nNonimm, nNonimm);
                    for (int c = 0; c < sparseD.getNumCols(); c++) {
                        Integer newCol = nonimmSet.get(c);
                        if (newCol == null) continue;
                        int idx0 = sparseD.col_idx[c];
                        int idx1 = sparseD.col_idx[c + 1];
                        for (int idx = idx0; idx < idx1; idx++) {
                            int r = sparseD.nz_rows[idx];
                            Integer newRow = nonimmSet.get(r);
                            if (newRow == null) continue;
                            dfilt_value.set(newRow, newCol, sparseD.nz_values[idx]);
                        }
                    }

                    if (Ta.getNumNonZeros() > 0) {
                        dfilt_value.add(Ta);
                    }

                    Dfilt.set(a, dfilt_value);
                }

                if (FJ == 0 && nImm == 0) {
                    // (when there are vanishing states the rates were already
                    // recomputed above on the nonimm-restricted index set)
                    double[][][] newDepRates = new double[nonimm.size()][nstateful][nclasses];
                    double[][][] newArvRates = new double[nonimm.size()][nstateful][nclasses];

                    for (int i = 0; i < nonimm.size(); i++) {
                        int origRow = nonimm.get(i).intValue();
                        for (int j = 0; j < nstateful; j++) {
                            for (int k = 0; k < nclasses; k++) {
                                newDepRates[i][j][k] = depRates[origRow][j][k];
                                newArvRates[i][j][k] = arvRates[origRow][j][k];
                            }
                        }
                    }

                    depRates = newDepRates;
                    arvRates = newArvRates;
                }
            }
        }

        return new ResultCTMC(Q, stateSpace, stateSpaceAggr, Dfilt, arvRates, depRates, sn);
    }

    /**
     * Clears every entry of A that lies in one of the given rows, touching only the
     * entries that actually exist. Scalar scatter into CSC is O(nnz) per write, so a
     * dense column scan over the vanishing rows would be quadratic; one pass over the
     * nonzeros is not. Entries are collected before being written because the iterator
     * must not be invalidated mid-traversal.
     */
    private static void zeroRows(Matrix A, Set<Integer> rows) {
        if (A == null || rows.isEmpty()) {
            return;
        }
        List<int[]> victims = new ArrayList<int[]>();
        Iterator<MatrixEntry> it = A.nonZeroIterator();
        while (it.hasNext()) {
            MatrixEntry e = it.next();
            if (rows.contains(Integer.valueOf(e.row))) {
                victims.add(new int[]{e.row, e.col});
            }
        }
        for (int i = 0; i < victims.size(); i++) {
            A.set(victims.get(i)[0], victims.get(i)[1], 0);
        }
    }

    /**
     * Indices of the vanishing (zero-sojourn) global states: Router/Fork pass-through
     * occupancy, firable Join sibling sets, SPN markings from which an ENABLE event
     * moves the Transition row, and markings enabling a TimingStrategy.IMMEDIATE mode.
     * Extracted so the same predicate drives both the vanishing-row purge and the
     * stochastic complementation. Mirrors the MATLAB local function
     * ctmc_find_vanishing_states in matlab/src/solvers/CTMC/solver_ctmc.m.
     */
    private static List<Double> ctmcFindVanishingStates(NetworkStruct sn, Matrix stateSpaceHashed,
            int nclasses, int nstateful, int FJ, Map<Integer, jline.lang.GlobalSync> gsyncEvents, int G) {
        List<Double> imm_list = new ArrayList<Double>();

        for (int ind = 0; ind < sn.nnodes; ind++) {
            NodeType nt = sn.nodetype.get(ind);
            // Fork qualifies only on FJ tag-augmented structs, where it is
            // stateful and holds the parent for one vanishing state
            boolean isImmediatePassThrough = (nt == NodeType.Router) || (nt == NodeType.Fork);
            if (sn.isstateful.get(ind) != 0.0 && sn.isstation.get(ind) == 0.0 && isImmediatePassThrough) {
                int isf = (int) sn.nodeToStateful.get(ind);

                Matrix space_slice = Matrix.extract(
                        sn.space.get(sn.stateful.get(isf)),
                        0,
                        sn.space.get(sn.stateful.get(isf)).getNumRows(),
                        0,
                        nclasses);
                Matrix rowSum = space_slice.sumRows();
                List<Integer> imm_st = new ArrayList<Integer>();
                for (int row = 0; row < rowSum.getNumRows(); row++) {
                    if (rowSum.get(row, 0) > 0) {
                        imm_st.add(row);
                    }
                }

                for (int s = 0; s < stateSpaceHashed.getNumRows(); s++) {
                    double hashValue = stateSpaceHashed.get(s, isf);
                    boolean anyMatch = false;
                    for (Integer immStIndex : imm_st) {
                        if (hashValue == (double) immStIndex.intValue()) {
                            anyMatch = true;
                            break;
                        }
                    }
                    if (anyMatch) {
                        imm_list.add((double) s);
                    }
                }
            }
        }

        // Transition immediate states
        if (gsyncEvents != null && G > 0) {
            Set<Double> immSet = new HashSet<Double>();
            for (Double d : imm_list) immSet.add(d);
            for (int s = 0; s < stateSpaceHashed.getNumRows(); s++) {
                if (immSet.contains((double) s)) continue;
                Matrix state_sc = stateSpaceHashed.getRow(s);
                List<Matrix> glspace_sc = new ArrayList<Matrix>(nstateful);
                for (int isf = 0; isf < nstateful; isf++) {
                    int stateIndex = (int) state_sc.get(isf);
                    Matrix spaceMatrix = sn.space.get(sn.stateful.get(isf));
                    glspace_sc.add(spaceMatrix.getRow(stateIndex));
                }
                for (int g_sc = 0; g_sc < G; g_sc++) {
                    Object glscObj = gsyncEvents.get(g_sc);
                    if (glscObj == null) continue;
                    jline.lang.GlobalSync glevent_sc = (jline.lang.GlobalSync) glscObj;
                    if (glevent_sc.active.isEmpty()) continue;
                    if (glevent_sc.active.get(0).getEvent() == EventType.ENABLE) {
                        int gind_sc = glevent_sc.active.get(0).getNode();
                        int isf_sc = (int) sn.nodeToStateful.get(gind_sc);
                        Matrix origTransState = glspace_sc.get(isf_sc);
                        AfterGlobalEvent.AfterGlobalEventResult result_sc = AfterGlobalEvent.afterGlobalEvent(sn, gind_sc, glspace_sc, glevent_sc, false);
                        Matrix outrate_sc = result_sc.outrate;
                        if (!outrate_sc.isEmpty() && outrate_sc.length() > 0) {
                            // see _kb/06-solver-catalog.md for rationale
                            Matrix newTransState = result_sc.outglspace.get(isf_sc);
                            boolean rowChanged = false;
                            int nOutcomes = Math.min(outrate_sc.length(), newTransState.getNumRows());
                            for (int io = 0; io < nOutcomes; io++) {
                                if (outrate_sc.get(io) > 0
                                        && !origTransState.isEqualTo(newTransState.getRow(io))) {
                                    rowChanged = true;
                                    break;
                                }
                            }
                            if (rowChanged) {
                                imm_list.add((double) s);
                                break;
                            }
                        }
                    } else if (glevent_sc.active.get(0).getEvent() == EventType.FIRE) {
                        // see _kb/06-solver-catalog.md for rationale
                        int gind_sc = glevent_sc.active.get(0).getNode();
                        int mode_sc = glevent_sc.active.get(0).getMode();
                        Object tpObj = sn.nodeparam.get(sn.nodes.get(gind_sc));
                        if (!(tpObj instanceof jline.lang.nodeparam.TransitionNodeParam)) continue;
                        jline.lang.nodeparam.TransitionNodeParam tp =
                                (jline.lang.nodeparam.TransitionNodeParam) tpObj;
                        if (tp.timing == null || mode_sc >= tp.timing.size()
                                || tp.timing.get(mode_sc) != TimingStrategy.IMMEDIATE) {
                            continue;
                        }
                        AfterGlobalEvent.AfterGlobalEventResult result_im =
                                AfterGlobalEvent.afterGlobalEvent(sn, gind_sc, glspace_sc, glevent_sc, false);
                        Matrix outrate_im = result_im.outrate;
                        if (outrate_im != null && !outrate_im.isEmpty()) {
                            boolean fires = false;
                            for (int io = 0; io < outrate_im.length(); io++) {
                                if (outrate_im.get(io) > 0.0) { fires = true; break; }
                            }
                            if (fires) {
                                imm_list.add((double) s);
                                break;
                            }
                        }
                    }
                }
            }
        }

        // see _kb/06-solver-catalog.md for rationale
        if (FJ > 0) {
            for (int ind = 0; ind < sn.nnodes; ind++) {
                if (sn.nodetype.get(ind) != NodeType.Join) continue;
                Object jnpObj = sn.nodeparam.get(sn.nodes.get(ind));
                if (!(jnpObj instanceof jline.lang.nodeparam.JoinNodeParam)) continue;
                jline.lang.nodeparam.JoinNodeParam jnp = (jline.lang.nodeparam.JoinNodeParam) jnpObj;
                if (jnp.fjOrigclasses == null) continue;
                int isf = (int) sn.nodeToStateful.get(ind);
                Matrix spaceJ = sn.space.get(sn.stateful.get(isf));
                Set<Integer> firableRows = new HashSet<Integer>();
                for (int row = 0; row < spaceJ.getNumRows(); row++) {
                    for (int rr : jnp.fjOrigclasses) {
                        Ret.EventResult jr = State.afterEvent(sn, ind, spaceJ.getRow(row), EventType.DEP, rr, false, new jline.lang.state.EventCache(false, false));
                        if (jr.outspace != null && !jr.outspace.isEmpty()) {
                            firableRows.add(row);
                            break;
                        }
                    }
                }
                if (!firableRows.isEmpty()) {
                    for (int s = 0; s < stateSpaceHashed.getNumRows(); s++) {
                        if (firableRows.contains((int) stateSpaceHashed.get(s, isf))) {
                            imm_list.add((double) s);
                        }
                    }
                }
            }
        }

        // Convert imm_list to unique sorted values
        TreeSet<Double> uniqueSet = new TreeSet<Double>(imm_list);
        List<Double> imm_unique = new ArrayList<Double>(uniqueSet);
        return imm_unique;
    }

    /** Rows {@code keep} of {@code m}, in the given order. */
    private static Matrix subrows(Matrix m, List<Integer> keep) {
        Matrix out = new Matrix(keep.size(), m.getNumCols());
        for (int i = 0; i < keep.size(); i++) {
            int r = keep.get(i).intValue();
            for (int c = 0; c < m.getNumCols(); c++) {
                out.set(i, c, m.get(r, c));
            }
        }
        return out;
    }

    /** The {@code rows} x {@code cols} submatrix of {@code m}. */
    private static Matrix submatrix(Matrix m, List<Integer> rows, List<Integer> cols) {
        Matrix out = new Matrix(rows.size(), cols.size());
        for (int i = 0; i < rows.size(); i++) {
            int r = rows.get(i).intValue();
            for (int j = 0; j < cols.size(); j++) {
                double v = m.get(r, cols.get(j).intValue());
                if (v != 0) {
                    out.set(i, j, v);
                }
            }
        }
        return out;
    }

    /** The state slice {@code keep} of a (states x stateful x classes) rate array. */
    private static double[][][] subrates(double[][][] rates, List<Integer> keep, int nstateful, int nclasses) {
        double[][][] out = new double[keep.size()][nstateful][nclasses];
        for (int i = 0; i < keep.size(); i++) {
            int r = keep.get(i).intValue();
            for (int j = 0; j < nstateful; j++) {
                for (int k = 0; k < nclasses; k++) {
                    out[i][j][k] = rates[r][j][k];
                }
            }
        }
        return out;
    }
}
