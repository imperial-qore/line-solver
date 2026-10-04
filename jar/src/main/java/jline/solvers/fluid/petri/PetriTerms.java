/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.petri;

import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.ProcessType;
import jline.lang.constant.TimingStrategy;
import jline.lang.nodeparam.TransitionNodeParam;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

/**
 * Event-based representation of the fluid marking process of a stochastic Petri
 * net. Java twin of the MATLAB {@code fluid_petri_terms}.
 *
 * <p>A GSPN is already a density-dependent Markov population process, which is
 * the object the moment-closure family of SolverFluid is built on: the marking
 * is the population, a transition mode is a reaction, its incidence column is
 * the jump, and the rate law {@code lambda*min(enabling degree, servers)} is the
 * same min() non-linearity the min-normal closure exists to smooth. Nothing
 * about the closure changes here; only where the drift comes from.
 *
 * <pre>    dx/dt = D * r(x, Sigma, phi, mu)</pre>
 *
 * <p>THE STATE, x = [ m ; y ].
 * <ul>
 *   <li>{@code m(p,k)} token mass of class k at place p. One coordinate per
 *   (place, class) pair some arc touches, the initial marking loads, or a Source
 *   feeds; a pair nothing reaches is dropped rather than carried as a null
 *   direction of the Newton system.</li>
 *   <li>{@code y(j,h)} the number of mode-j servers running in phase h, for a
 *   mode whose firing time has more than one phase. Their SUM is not free: the
 *   ENABLE synchronization latches it instantaneously to min(enabling degree,
 *   servers), so the latch is an ALGEBRAIC row with one free-sign unknown mu_j
 *   and the phase split evolves differentially. CARRYING THE DISTRIBUTION
 *   INSTEAD OF THE COUNT LOOKS TIDIER AND IS WRONG: the resulting equation is
 *   missing a term and agrees with the count form only AT a fixed point.</li>
 * </ul>
 *
 * <p>There is no Source coordinate: an exogenous arrival is a CONSTANT-propensity
 * event depositing one token, as it is in the NRM SPN runner. There is no Sink
 * coordinate either: a firing arc into a sink is mass leaving the net.
 *
 * <p>THE EVENTS, one column of D each: kind 1 a firing of mode j out of phase h
 * into phase h'; kind 2 an internal phase change; kind 3 an exogenous arrival;
 * kind 4 a firing of an IMMEDIATE mode, at the algebraic flow phi_j; kind 5 the
 * server latch of a multi-phase mode, at the free-sign unknown mu_j.
 */
public class PetriTerms {

    public int M;
    public int K;
    public int I;
    public List<Integer> places = new ArrayList<Integer>();
    public List<Integer> transitions = new ArrayList<Integer>();
    public List<String> namesNode;
    public int nstate;
    public int nm;
    /** (node, class) -> state coordinate, -1 where the pair carries none. */
    public int[][] pidx;
    public int[] coordNode;
    public int[] coordClass;
    public int[] coordStation;
    public List<PetriMode> modes = new ArrayList<PetriMode>();
    public List<Integer> timedIdx = new ArrayList<Integer>();
    public List<Integer> immIdx = new ArrayList<Integer>();
    public Matrix D;
    public double[] rateBase;
    public int nev;
    public int[] evKind;
    public int[] evMode;
    public int[] evPhase;
    public int[] evTo;
    public int[] evStation;
    public int[] evClass;
    public int[] immCol;
    public int[] latchCol;
    /**
     * The DIFFUSION counts the stochastic events only: an immediate flow and a
     * server latch are both the limit of an infinitely fast mechanism whose
     * fluctuation is slaved, not a Poisson stream with an intensity.
     */
    public int[] stochCol;
    public List<Integer> latchMode = new ArrayList<Integer>();
    /**
     * THE COVARIANCE COVERS EVERY COORDINATE, phases included: the rate of a
     * multi-phase mode is linear in y and reads no marking, so dropping the
     * phases would sever that mode's whole restoring force.
     */
    public int[] covIdx;
    public int[][] covPairs;
    public int npair;
    /** (a,b) -> closure unknown index, symmetric, -1 where the drift never reads it. */
    public int[][] pairIndex;
    public int[] immColOf;
    public Map<Integer, List<Integer>> consumers = new HashMap<Integer, List<Integer>>();
    public Map<Integer, List<Double>> consumerW = new HashMap<Integer, List<Double>>();
    public Map<Integer, List<Integer>> producers = new HashMap<Integer, List<Integer>>();
    public double[] x0;
    public double[][] m0full;
    public SolverOptions options;

    /** Key for the per-(station,class) accumulators. */
    private int key(int ist, int k) {
        return ist * K + k;
    }

    public List<Integer> consumersOf(int ist, int k) {
        return consumers.get(key(ist, k));
    }

    public List<Double> consumerWOf(int ist, int k) {
        return consumerW.get(key(ist, k));
    }

    public List<Integer> producersOf(int ist, int k) {
        return producers.get(key(ist, k));
    }

    /** Pad an arc matrix to (nnodes x nclasses); addMode sizes them at creation. */
    private static double[][] pad(Matrix A, int I, int K, double fill) {
        double[][] out = new double[I][K];
        for (int i = 0; i < I; i++) {
            for (int j = 0; j < K; j++) {
                out[i][j] = fill;
            }
        }
        if (A == null) {
            return out;
        }
        int rows = Math.min(I, A.getNumRows());
        int cols = Math.min(K, A.getNumCols());
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                out[i][j] = A.get(i, j);
            }
        }
        return out;
    }

    public static PetriTerms build(NetworkStruct sn, SolverOptions options) {
        PetriTerms t = new PetriTerms();
        t.options = options;
        t.I = sn.nnodes;
        t.K = sn.nclasses;
        t.M = sn.nstations;
        t.namesNode = new ArrayList<String>();
        for (int i = 0; i < t.I; i++) {
            t.namesNode.add(sn.nodenames.get(i));
        }
        for (int i = 0; i < t.I; i++) {
            if (sn.nodetype.get(i) == NodeType.Place) {
                t.places.add(i);
            } else if (sn.nodetype.get(i) == NodeType.Transition) {
                t.transitions.add(i);
            }
        }

        // ---- the initial marking, read off the model state exactly as the NRM does
        t.m0full = new double[t.I][t.K];
        for (int idx = 0; idx < t.places.size(); idx++) {
            int ind = t.places.get(idx);
            Matrix stateI = sn.state.get(sn.stateful.get((int) sn.nodeToStateful.get(ind)));
            State.StateMarginalStatistics aggr =
                    ToMarginal.toMarginalAggr(sn, ind, stateI, null, null, null, null, null);
            for (int k = 0; k < t.K; k++) {
                double nir = aggr.nir.get(k);
                if (Double.isInfinite(nir)) {
                    line_error(mfilename(new Object() {
                    }), String.format("Place %s holds an infinite initial marking of class %d.",
                            t.namesNode.get(ind), k + 1));
                }
                t.m0full[ind][k] = nir;
            }
        }

        // ---- which (place, class) pairs carry a coordinate
        boolean[][] touched = new boolean[t.I][t.K];
        for (int ti = 0; ti < t.transitions.size(); ti++) {
            int ind = t.transitions.get(ti);
            TransitionNodeParam tp = (TransitionNodeParam) sn.nodeparam.get(sn.nodes.get(ind));
            for (int m = 0; m < tp.nmodes; m++) {
                double[][] en = pad(tp.enabling.get(m), t.I, t.K, 0.0);
                double[][] fir = pad(tp.firing.get(m), t.I, t.K, 0.0);
                double[][] inh = pad(tp.inhibiting.get(m), t.I, t.K, Double.POSITIVE_INFINITY);
                for (int i = 0; i < t.I; i++) {
                    for (int k = 0; k < t.K; k++) {
                        if (en[i][k] > 0 || fir[i][k] != 0
                                || (!Double.isInfinite(inh[i][k]) && inh[i][k] > 0)) {
                            touched[i][k] = true;
                        }
                    }
                }
            }
        }

        // an arrival makes its target a coordinate even when no arc mentions it
        List<int[]> srcArr = new ArrayList<int[]>();
        for (int ind = 0; ind < t.I; ind++) {
            if (sn.nodetype.get(ind) != NodeType.Source) {
                continue;
            }
            int ist = (int) sn.nodeToStation.get(ind);
            for (int r = 0; r < t.K; r++) {
                double lambda = sn.rates.get(ist, r);
                if (Double.isNaN(lambda) || lambda <= 0) {
                    continue;
                }
                for (int pi = 0; pi < t.places.size(); pi++) {
                    int jnd = t.places.get(pi);
                    for (int s = 0; s < t.K; s++) {
                        double p = sn.rtnodes.get(ind * t.K + r, jnd * t.K + s);
                        if (p > 0) {
                            touched[jnd][s] = true;
                            srcArr.add(new int[]{ind, r, jnd, s});
                        }
                    }
                }
            }
        }

        t.pidx = new int[t.I][t.K];
        for (int i = 0; i < t.I; i++) {
            for (int k = 0; k < t.K; k++) {
                t.pidx[i][k] = -1;
            }
        }
        List<Integer> cn = new ArrayList<Integer>();
        List<Integer> cc = new ArrayList<Integer>();
        List<Integer> cs = new ArrayList<Integer>();
        int nm = 0;
        for (int pi = 0; pi < t.places.size(); pi++) {
            int ind = t.places.get(pi);
            for (int k = 0; k < t.K; k++) {
                if (!(touched[ind][k] || t.m0full[ind][k] > 0)) {
                    continue;
                }
                t.pidx[ind][k] = nm;
                cn.add(ind);
                cc.add(k);
                cs.add((int) sn.nodeToStation.get(ind));
                nm++;
            }
        }
        t.nm = nm;
        t.coordNode = toIntArray(cn);
        t.coordClass = toIntArray(cc);
        t.coordStation = toIntArray(cs);

        // ---- the modes, and the phase coordinates of the multi-phase ones
        int nstate = nm;
        for (int ti = 0; ti < t.transitions.size(); ti++) {
            int ind = t.transitions.get(ti);
            TransitionNodeParam tp = (TransitionNodeParam) sn.nodeparam.get(sn.nodes.get(ind));
            for (int m = 0; m < tp.nmodes; m++) {
                PetriMode rec = buildMode(sn, tp, ind, m, t, nm);
                if (rec.nph > 1) {
                    rec.zblk = new int[rec.nph];
                    for (int h = 0; h < rec.nph; h++) {
                        rec.zblk[h] = nstate + h;
                    }
                    nstate += rec.nph;
                }
                t.modes.add(rec);
            }
        }
        t.nstate = nstate;
        // The phase coordinates are appended after every marking coordinate, so a
        // jump column built at nm width has to be GROWN once the total is known.
        for (int j = 0; j < t.modes.size(); j++) {
            PetriMode md = t.modes.get(j);
            if (md.cvec.length < nstate) {
                double[] grown = new double[nstate];
                System.arraycopy(md.cvec, 0, grown, 0, md.cvec.length);
                md.cvec = grown;
            }
        }
        for (int j = 0; j < t.modes.size(); j++) {
            if (t.modes.get(j).timing == TimingStrategy.TIMED) {
                t.timedIdx.add(j);
            } else {
                t.immIdx.add(j);
            }
        }

        // ---- the event columns
        List<double[]> cols = new ArrayList<double[]>();
        List<Double> rateBase = new ArrayList<Double>();
        List<Integer> evKind = new ArrayList<Integer>();
        List<Integer> evMode = new ArrayList<Integer>();
        List<Integer> evPhase = new ArrayList<Integer>();
        List<Integer> evTo = new ArrayList<Integer>();
        List<Integer> evStation = new ArrayList<Integer>();
        List<Integer> evClass = new ArrayList<Integer>();

        for (int q = 0; q < t.timedIdx.size(); q++) {
            int j = t.timedIdx.get(q);
            PetriMode md = t.modes.get(j);
            if (md.nph == 1) {
                cols.add(md.cvec.clone());
                rateBase.add(md.d1[0]);
                evKind.add(1); evMode.add(j); evPhase.add(0); evTo.add(0);
                evStation.add(-1); evClass.add(-1);
            } else {
                for (int h = 0; h < md.nph; h++) {
                    for (int hp = 0; hp < md.nph; hp++) {
                        double w = md.D1.get(h, hp);
                        if (w <= 0) {
                            continue;
                        }
                        double[] col = md.cvec.clone();
                        col[md.zblk[hp]] += 1.0;
                        col[md.zblk[h]] -= 1.0;
                        cols.add(col);
                        rateBase.add(w);
                        evKind.add(1); evMode.add(j); evPhase.add(h); evTo.add(hp);
                        evStation.add(-1); evClass.add(-1);
                    }
                }
                for (int h = 0; h < md.nph; h++) {
                    for (int hp = 0; hp < md.nph; hp++) {
                        if (hp == h) {
                            continue;
                        }
                        double w = md.D0.get(h, hp);
                        if (w <= 0) {
                            continue;
                        }
                        double[] col = new double[nstate];
                        col[md.zblk[hp]] = 1.0;
                        col[md.zblk[h]] = -1.0;
                        cols.add(col);
                        rateBase.add(w);
                        evKind.add(2); evMode.add(j); evPhase.add(h); evTo.add(hp);
                        evStation.add(-1); evClass.add(-1);
                    }
                }
            }
        }
        // One latch column per multi-phase mode: mu_j servers per unit time enter
        // at the firing process's own entry distribution. The rate is FREE IN
        // SIGN -- a mode whose enabling degree drops stops servers rather than
        // starting them -- and it is zero at any fixed point.
        for (int q = 0; q < t.timedIdx.size(); q++) {
            int j = t.timedIdx.get(q);
            PetriMode md = t.modes.get(j);
            if (md.nph <= 1) {
                continue;
            }
            double[] col = new double[nstate];
            for (int h = 0; h < md.nph; h++) {
                col[md.zblk[h]] = md.pie[h];
            }
            cols.add(col);
            rateBase.add(1.0);
            evKind.add(5); evMode.add(j); evPhase.add(0); evTo.add(0);
            evStation.add(-1); evClass.add(-1);
        }
        for (int q = 0; q < t.immIdx.size(); q++) {
            int j = t.immIdx.get(q);
            cols.add(t.modes.get(j).cvec.clone());
            rateBase.add(1.0);
            evKind.add(4); evMode.add(j); evPhase.add(0); evTo.add(0);
            evStation.add(-1); evClass.add(-1);
        }
        for (int a = 0; a < srcArr.size(); a++) {
            int[] rec = srcArr.get(a);
            int snd = rec[0], r = rec[1], qnd = rec[2], l = rec[3];
            int ist = (int) sn.nodeToStation.get(snd);
            if (sn.procid != null && sn.procid.containsKey(sn.stations.get(ist))
                    && sn.procid.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)) != ProcessType.EXP) {
                line_error(mfilename(new Object() {
                }), String.format(
                        "Source %s has a non-exponential arrival for class %d. The fluid Petri route models "
                                + "an arrival as a constant-propensity event, which a renewal stream with memory "
                                + "is not; use SolverCTMC, SolverJMT or SolverSSA.", t.namesNode.get(snd), r + 1));
            }
            double[] col = new double[nstate];
            col[t.pidx[qnd][l]] = 1.0;
            cols.add(col);
            rateBase.add(sn.rates.get(ist, r) * sn.rtnodes.get(snd * t.K + r, qnd * t.K + l));
            evKind.add(3); evMode.add(-1); evPhase.add(0); evTo.add(0);
            evStation.add(ist); evClass.add(r);
        }

        t.nev = cols.size();
        t.D = new Matrix(nstate, Math.max(t.nev, 1));
        t.D.fill(0.0);
        for (int e = 0; e < t.nev; e++) {
            double[] col = cols.get(e);
            for (int s = 0; s < nstate; s++) {
                if (col[s] != 0) {
                    t.D.set(s, e, col[s]);
                }
            }
        }
        t.rateBase = toDoubleArray(rateBase);
        t.evKind = toIntArray(evKind);
        t.evMode = toIntArray(evMode);
        t.evPhase = toIntArray(evPhase);
        t.evTo = toIntArray(evTo);
        t.evStation = toIntArray(evStation);
        t.evClass = toIntArray(evClass);

        // ---- which Sigma entries the closure reads
        Map<Long, Integer> pairKey = new HashMap<Long, Integer>();
        List<int[]> pairs = new ArrayList<int[]>();
        for (int j = 0; j < t.modes.size(); j++) {
            PetriMode md = t.modes.get(j);
            if (md.closable) {
                for (int a = 0; a < md.arcSlot.size(); a++) {
                    for (int b = a; b < md.arcSlot.size(); b++) {
                        addPair(pairKey, pairs, md.arcSlot.get(a), md.arcSlot.get(b));
                    }
                }
            }
            if (md.timing == TimingStrategy.TIMED) {
                for (int b = 0; b < md.inhSlot.size(); b++) {
                    addPair(pairKey, pairs, md.inhSlot.get(b), md.inhSlot.get(b));
                }
            }
        }
        t.npair = pairs.size();
        t.covPairs = new int[t.npair][2];
        for (int i = 0; i < t.npair; i++) {
            t.covPairs[i] = pairs.get(i);
        }
        int pn = Math.max(nm, 1);
        t.pairIndex = new int[pn][pn];
        for (int i = 0; i < pn; i++) {
            for (int j = 0; j < pn; j++) {
                t.pairIndex[i][j] = -1;
            }
        }
        for (int i = 0; i < t.npair; i++) {
            t.pairIndex[t.covPairs[i][0]][t.covPairs[i][1]] = i;
            t.pairIndex[t.covPairs[i][1]][t.covPairs[i][0]] = i;
        }

        List<Integer> imm = new ArrayList<Integer>();
        List<Integer> latch = new ArrayList<Integer>();
        List<Integer> stoch = new ArrayList<Integer>();
        for (int e = 0; e < t.nev; e++) {
            if (t.evKind[e] == 4) {
                imm.add(e);
            } else if (t.evKind[e] == 5) {
                latch.add(e);
            }
            if (t.evKind[e] != 4 && t.evKind[e] != 5) {
                stoch.add(e);
            }
        }
        t.immCol = toIntArray(imm);
        t.latchCol = toIntArray(latch);
        t.stochCol = toIntArray(stoch);
        for (int q = 0; q < t.timedIdx.size(); q++) {
            int j = t.timedIdx.get(q);
            if (t.modes.get(j).nph > 1) {
                t.latchMode.add(j);
            }
        }
        t.immColOf = new int[t.modes.size()];
        for (int j = 0; j < t.modes.size(); j++) {
            t.immColOf[j] = -1;
        }
        for (int q = 0; q < t.immCol.length; q++) {
            t.immColOf[t.evMode[t.immCol[q]]] = t.immCol[q];
        }
        t.covIdx = new int[nstate];
        for (int s = 0; s < nstate; s++) {
            t.covIdx[s] = s;
        }

        // ---- per-place consumption and production, for the metric reader
        // A Place's throughput is the rate at which TOKENS leave it, so each
        // consuming mode contributes its firing rate times the multiplicity of
        // the arc it takes them through. That is SolverCTMC's convention and the
        // one Little's law needs; SolverSSA's NRM sums the UNWEIGHTED propensity,
        // so the two disagree wherever an input arc has multiplicity above one.
        for (int e = 0; e < t.nev; e++) {
            if (t.evKind[e] == 3) {
                int kk = t.evStation[e] * t.K + t.evClass[e];
                if (!t.producers.containsKey(kk)) {
                    t.producers.put(kk, new ArrayList<Integer>());
                }
                t.producers.get(kk).add(e);
                continue;
            }
            if (t.evKind[e] != 1 && t.evKind[e] != 4) {
                continue;
            }
            PetriMode md = t.modes.get(t.evMode[e]);
            for (int a = 0; a < md.arcSlot.size(); a++) {
                int s = md.arcSlot.get(a);
                int kk = t.coordStation[s] * t.K + t.coordClass[s];
                if (!t.consumers.containsKey(kk)) {
                    t.consumers.put(kk, new ArrayList<Integer>());
                    t.consumerW.put(kk, new ArrayList<Double>());
                }
                t.consumers.get(kk).add(e);
                t.consumerW.get(kk).add(md.arcW.get(a));
            }
        }

        // ---- the initial state
        t.x0 = new double[nstate];
        for (int s = 0; s < nm; s++) {
            t.x0[s] = t.m0full[t.coordNode[s]][t.coordClass[s]];
        }
        for (int j = 0; j < t.modes.size(); j++) {
            PetriMode md = t.modes.get(j);
            if (md.nph > 1) {
                double e = Double.POSITIVE_INFINITY;
                for (int a = 0; a < md.arcSlot.size(); a++) {
                    e = Math.min(e, t.x0[md.arcSlot.get(a)] / md.arcW.get(a));
                }
                if (md.arcSlot.isEmpty()) {
                    e = 1.0;
                }
                for (int h = 0; h < md.nph; h++) {
                    t.x0[md.zblk[h]] = Math.min(e, md.c) * md.pie[h];
                }
            }
        }
        return t;
    }

    private static void addPair(Map<Long, Integer> key, List<int[]> pairs, int a, int b) {
        int lo = Math.min(a, b);
        int hi = Math.max(a, b);
        Long k = Long.valueOf(((long) lo << 32) | (hi & 0xffffffffL));
        if (!key.containsKey(k)) {
            key.put(k, pairs.size());
            pairs.add(new int[]{lo, hi});
        }
    }

    private static PetriMode buildMode(NetworkStruct sn, TransitionNodeParam tp, int ind, int m,
                                       PetriTerms t, int nm) {
        String modeName = (tp.modenames != null && tp.modenames.size() > m && tp.modenames.get(m) != null)
                ? tp.modenames.get(m) : ("Mode" + (m + 1));
        TimingStrategy timing = (tp.timing != null && tp.timing.size() > m)
                ? tp.timing.get(m) : TimingStrategy.TIMED;
        PetriMode rec = new PetriMode(ind, m, timing,
                t.namesNode.get(ind) + "." + modeName, nm);

        double[][] en = pad(tp.enabling.get(m), t.I, t.K, 0.0);
        double[][] fir = pad(tp.firing.get(m), t.I, t.K, 0.0);
        double[][] inh = pad(tp.inhibiting.get(m), t.I, t.K, Double.POSITIVE_INFINITY);

        double[] cvec = new double[nm];
        for (int p = 0; p < t.I; p++) {
            for (int k = 0; k < t.K; k++) {
                double w = en[p][k];
                if (!(w > 0)) {
                    continue;
                }
                if (Double.isInfinite(w)) {
                    line_error(mfilename(new Object() {
                    }), String.format("Mode %s has a non-finite enabling arc weight at %s. An arc that no "
                            + "marking can satisfy disables the mode; declare a finite multiplicity.",
                            rec.label, t.namesNode.get(p)));
                }
                if (t.pidx[p][k] < 0) {
                    line_error(mfilename(new Object() {
                    }), String.format("Mode %s takes an enabling arc from %s, which is not a Place.",
                            rec.label, t.namesNode.get(p)));
                }
                rec.arcSlot.add(t.pidx[p][k]);
                rec.arcW.add(w);
                cvec[t.pidx[p][k]] -= w;
            }
        }
        for (int p = 0; p < t.I; p++) {
            for (int k = 0; k < t.K; k++) {
                if (fir[p][k] <= 0) {
                    continue; // a negative entry marks an input place, whose token the PRE already removed
                }
                if (t.pidx[p][k] < 0) {
                    continue; // a firing arc into a Sink is mass leaving the net
                }
                cvec[t.pidx[p][k]] += fir[p][k];
            }
        }
        for (int p = 0; p < t.I; p++) {
            for (int k = 0; k < t.K; k++) {
                double th = inh[p][k];
                if (Double.isInfinite(th) || th <= 0 || t.pidx[p][k] < 0) {
                    continue; // JMT writes a missing inhibitor arc as 0 or -1, never as a threshold
                }
                rec.inhSlot.add(t.pidx[p][k]);
                rec.inhThr.add(th);
            }
        }
        rec.cvec = cvec;

        double c = (tp.nmodeservers != null) ? tp.nmodeservers.get(0, m) : 1.0;
        rec.c = Double.isNaN(c) ? 1.0 : c;
        rec.prio = (tp.firingprio != null) ? (int) tp.firingprio.get(0, m) : 1;
        rec.weight = (tp.fireweight != null) ? tp.fireweight.get(0, m) : 1.0;
        if (tp.firingdep != null && tp.firingdep.size() > m) {
            rec.dep = tp.firingdep.get(m);
        }

        if (timing == TimingStrategy.IMMEDIATE) {
            // An immediate mode carries no firing process: its flow is an
            // algebraic unknown of the DAE, not a rate.
            rec.nph = 0;
            rec.d1 = new double[]{0.0};
            rec.closable = false;
            return rec;
        }

        // firingproc / firingpie are keyed by the Mode OBJECT, not the index
        jline.lang.Mode modeObj = ((jline.lang.nodes.Transition) sn.nodes.get(ind)).getModes().get(m);
        MatrixCell proc = (tp.firingproc != null) ? tp.firingproc.get(modeObj) : null;
        if (proc == null || proc.get(0) == null) {
            line_error(mfilename(new Object() {
            }), String.format("Mode %s has no Markovian firing process. The renewal families are converted "
                    + "to phase type before the solver runs, so this is a distribution the fluid Petri "
                    + "route cannot time; use SolverCTMC or SolverLDES.", rec.label));
        }
        Matrix D0 = proc.get(0);
        Matrix D1 = proc.get(1);
        rec.nph = D0.getNumRows();
        rec.D0 = D0;
        rec.D1 = D1;
        rec.d1 = new double[rec.nph];
        for (int h = 0; h < rec.nph; h++) {
            double acc = 0;
            for (int hp = 0; hp < D1.getNumCols(); hp++) {
                acc += D1.get(h, hp);
            }
            rec.d1[h] = acc;
        }
        if (rec.nph > 1) {
            Matrix pieM = (tp.firingpie != null) ? tp.firingpie.get(modeObj) : null;
            rec.pie = new double[rec.nph];
            double tot = 0;
            if (pieM == null || pieM.getNumElements() == 0) {
                rec.pie[0] = 1.0;
                tot = 1.0;
            } else {
                for (int h = 0; h < rec.nph; h++) {
                    rec.pie[h] = pieM.get(h);
                    tot += rec.pie[h];
                }
            }
            if (tot > 0) {
                for (int h = 0; h < rec.nph; h++) {
                    rec.pie[h] /= tot;
                }
            }
        }
        rec.closable = !(rec.arcSlot.size() <= 1 && Double.isInfinite(rec.c));
        return rec;
    }

    private static int[] toIntArray(List<Integer> l) {
        int[] out = new int[l.size()];
        for (int i = 0; i < l.size(); i++) {
            out[i] = l.get(i);
        }
        return out;
    }

    private static double[] toDoubleArray(List<Double> l) {
        double[] out = new double[l.size()];
        for (int i = 0; i < l.size(); i++) {
            out[i] = l.get(i);
        }
        return out;
    }
}
