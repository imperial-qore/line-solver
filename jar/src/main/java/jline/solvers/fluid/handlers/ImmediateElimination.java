/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.handlers;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;

import static jline.io.InputOutput.line_warning;

/**
 * Stochastic complementation of the IMMEDIATE coordinates of the fluid ODE.
 *
 * <p>A coordinate whose exit rate is {@link GlobalConstants#Immediate} (= 1/FineTol = 1e8) is not a
 * fast coordinate, it is an INSTANTANEOUS one: the rate is LINE's stand-in for infinity, written by
 * SolverLN for the branch of an activity that takes no time (an entry called with probability
 * y &lt; 1 carries a second PH phase at InfRate entered with probability 1-y). Integrating it
 * numerically is meaningless work -- the mode relaxes 1e8 times faster than anything else in the
 * model -- and it is what made a two-station LN layer take 145 s here and 70 s in MATLAB for an
 * answer identical, to every digit, to the one the reduced system returns in 0.35 s.
 *
 * <p>THE REDUCTION IS EXACT, not an approximation. It is the ODE twin of ctmc_stochcomp: the
 * instantaneous coordinates F are absorbed into the timed ones S by the absorption probabilities of
 * the embedded jump chain restricted to F, so the flow that would enter F is routed straight to
 * where F would have sent it.
 *
 * <p>WHY THIS IS A STRUCTURAL COMPOSITION AND NOT A GENERATOR ROUND TRIP. Every event of the fluid
 * event set is a single -1 at its source coordinate and a single +1 at its destination, so a path
 * through F composes to one event, -1 at the original source and +1 at the absorbing coordinate,
 * that keeps the original source's GATING. Rebuilding the events from a reduced generator instead
 * (the shape this class had before) loses that identity, and with it the event ORDER that
 * FluidMomentTerms reads throughputs off -- which is why the moment-closure methods used to refuse
 * the reduction outright. {@code Emap} carries the identity across instead: {@code Emap(e,o)} is the
 * expected number of times the ORIGINAL event o fires per firing of the reduced event e, so a caller
 * maps any per-event quantity with {@code newAttr = Emap * oldAttr}. It is the identity when nothing
 * is eliminated.
 *
 * <p>A COMPOSED EVENT CAN BE A DEPARTURE AT TWO STATIONS AT ONCE, which is why a single
 * evIsDeparture flag cannot survive the composition: a job that leaves the delay, passes through the
 * queue's immediate phase and returns has completed at BOTH, and both throughputs must count it.
 *
 * @see jline.api.mc.Ctmc_stochcomp
 */
public class ImmediateElimination {

    /**
     * Result of immediate elimination containing the reduced system, the surviving coordinates and
     * the two maps that carry per-event and per-coordinate quantities across the reduction.
     */
    public static class EliminationResult {
        public Matrix allJumpsReduced;
        public Matrix rateBaseReduced;
        public Matrix eventIdxReduced;
        public int[] stateMap;
        /** [nEventsReduced x nEventsOriginal] expected firings of each original event. */
        public Matrix Emap;
        /**
         * [nStates x nStates] projector for the initial condition: identity on the timed rows, the
         * absorption distribution on the immediate ones. Mass parked on an eliminated coordinate
         * would otherwise be frozen there for the whole integration, because nothing moves it.
         */
        public Matrix absorb;

        public EliminationResult(Matrix allJumpsReduced, Matrix rateBaseReduced,
                                 Matrix eventIdxReduced, int[] stateMap,
                                 Matrix Emap, Matrix absorb) {
            this.allJumpsReduced = allJumpsReduced;
            this.rateBaseReduced = rateBaseReduced;
            this.eventIdxReduced = eventIdxReduced;
            this.stateMap = stateMap;
            this.Emap = Emap;
            this.absorb = absorb;
        }
    }

    private static EliminationResult identity(Matrix allJumps, Matrix rateBase, Matrix eventIdx) {
        int nStates = (int) allJumps.getNumRows();
        int nEvents = (int) rateBase.getNumRows();
        int[] map = new int[nStates];
        for (int i = 0; i < nStates; i++) {
            map[i] = i;
        }
        return new EliminationResult(allJumps, rateBase, eventIdx, map,
                Matrix.eye(nEvents), Matrix.eye(nStates));
    }

    /**
     * Eliminate the instantaneous coordinates from the fluid event set.
     *
     * @param allJumps [nStates x nEvents] jump matrix
     * @param rateBase [nEvents x 1] fixed part of each event rate
     * @param eventIdx [nEvents x 1] source coordinate of each event, which also gates it
     * @param sn       model structure, unused, kept for the caller's signature
     * @param options  solver options; {@code options.config immediate_tol} overrides the threshold
     * @return the reduced system together with Emap and absorb
     */
    public static EliminationResult eliminateImmediate(
            Matrix allJumps,
            Matrix rateBase,
            Matrix eventIdx,
            NetworkStruct sn,
            SolverOptions options) {
        try {
            return eliminateStructural(allJumps, rateBase, eventIdx, options);
        } catch (Exception e) {
            line_warning("ImmediateElimination",
                    String.format("Immediate coordinate elimination failed: %s. Using original system.",
                            e.getMessage()));
            return identity(allJumps, rateBase, eventIdx);
        }
    }

    private static EliminationResult eliminateStructural(
            Matrix allJumps, Matrix rateBase, Matrix eventIdx, SolverOptions options) {

        int nStates = (int) allJumps.getNumRows();
        int nEvents = (int) rateBase.getNumRows();

        // A rate at or above the threshold is the InfRate sentinel, not a fast rate the user wrote:
        // the default sits just under GlobalConstants.Immediate so that only the sentinel qualifies.
        double immTol = GlobalConstants.Immediate * (1 - 1e-2);
        if (options.config != null && options.config.containsKey("immediate_tol")) {
            immTol = (Double) options.config.get("immediate_tol");
        }

        boolean[] isImm = new boolean[nStates];
        boolean anyImm = false;
        for (int e = 0; e < nEvents; e++) {
            if (rateBase.get(e, 0) >= immTol) {
                isImm[(int) eventIdx.get(e, 0)] = true;
                anyImm = true;
            }
        }
        if (!anyImm) {
            return identity(allJumps, rateBase, eventIdx);
        }

        // Destination of each event. Every event is -1 at its source and +1 at one destination; a
        // departure that re-enters its own coordinate cancels to an all-zero column, whose
        // destination is that same coordinate.
        int[] src = new int[nEvents];
        int[] dst = new int[nEvents];
        for (int e = 0; e < nEvents; e++) {
            src[e] = (int) eventIdx.get(e, 0);
            dst[e] = src[e];
            for (int s = 0; s < nStates; s++) {
                if (allJumps.get(s, e) > 0) {
                    dst[e] = s;
                    break;
                }
            }
        }

        // A coordinate with no outflow at all cannot be complemented away, and one whose outflow is
        // entirely a self-loop would make the fundamental matrix singular. Both are dropped from F
        // rather than guessed at.
        for (int f = 0; f < nStates; f++) {
            if (!isImm[f]) {
                continue;
            }
            double tot = 0;
            boolean leaves = false;
            for (int e = 0; e < nEvents; e++) {
                if (src[e] == f) {
                    tot += rateBase.get(e, 0);
                    if (dst[e] != f) {
                        leaves = true;
                    }
                }
            }
            if (tot <= 0 || !leaves) {
                isImm[f] = false;
            }
        }

        List<Integer> Flist = new ArrayList<Integer>();
        List<Integer> Slist = new ArrayList<Integer>();
        for (int i = 0; i < nStates; i++) {
            if (isImm[i]) {
                Flist.add(i);
            } else {
                Slist.add(i);
            }
        }
        if (Flist.isEmpty() || Slist.isEmpty()) {
            return identity(allJumps, rateBase, eventIdx);
        }
        int nF = Flist.size();
        int nS = Slist.size();
        int[] posF = new int[nStates];
        int[] posS = new int[nStates];
        for (int a = 0; a < nF; a++) {
            posF[Flist.get(a)] = a;
        }
        for (int b = 0; b < nS; b++) {
            posS[Slist.get(b)] = b;
        }

        // Branching of the embedded jump chain out of each immediate coordinate. The probabilities
        // are the rate shares, so a coordinate carrying both an immediate and an ordinary exit gives
        // the ordinary one its (vanishing) share rather than being special-cased.
        Matrix PFF = new Matrix(nF, nF);
        Matrix PFS = new Matrix(nF, nS);
        Matrix cnt = new Matrix(nF, nEvents);
        for (int a = 0; a < nF; a++) {
            int f = Flist.get(a);
            double tot = 0;
            for (int e = 0; e < nEvents; e++) {
                if (src[e] == f) {
                    tot += rateBase.get(e, 0);
                }
            }
            for (int e = 0; e < nEvents; e++) {
                if (src[e] != f) {
                    continue;
                }
                double p = rateBase.get(e, 0) / tot;
                cnt.set(a, e, cnt.get(a, e) + p);
                if (isImm[dst[e]]) {
                    int c = posF[dst[e]];
                    PFF.set(a, c, PFF.get(a, c) + p);
                } else {
                    int c = posS[dst[e]];
                    PFS.set(a, c, PFS.get(a, c) + p);
                }
            }
        }

        // Fundamental matrix of the instantaneous chain. (I-PFF) is invertible whenever every
        // immediate coordinate reaches a timed one; a closed cycle of immediate coordinates is a
        // modelling error and is left to the original system rather than solved.
        Matrix ImP = Matrix.eye(nF);
        for (int a = 0; a < nF; a++) {
            for (int b = 0; b < nF; b++) {
                ImP.set(a, b, ImP.get(a, b) - PFF.get(a, b));
            }
        }
        Matrix Nfm;
        try {
            Nfm = ImP.inv();
        } catch (Exception ex) {
            line_warning("ImmediateElimination", "the immediate coordinates form a closed cycle, "
                    + "so they have no absorption distribution; integrating the unreduced system instead");
            return identity(allJumps, rateBase, eventIdx);
        }
        Matrix Aabs = Nfm.mult(PFS);
        Matrix expCnt = Nfm.mult(cnt);
        for (int a = 0; a < nF; a++) {
            double rowSum = 0;
            for (int b = 0; b < nS; b++) {
                rowSum += Aabs.get(a, b);
            }
            if (!(rowSum > 0.5) || Double.isNaN(rowSum)) {
                line_warning("ImmediateElimination", "the immediate coordinates form a closed cycle, "
                        + "so they have no absorption distribution; integrating the unreduced system instead");
                return identity(allJumps, rateBase, eventIdx);
            }
        }

        // Compose the event list. An event sourced in F is dropped: its flow is already carried by
        // whichever event feeds F.
        List<double[]> jumpsNew = new ArrayList<double[]>();
        List<Double> rateNew = new ArrayList<Double>();
        List<Integer> evidxNew = new ArrayList<Integer>();
        List<int[]> emapIdx = new ArrayList<int[]>();
        List<Double> emapVal = new ArrayList<Double>();
        for (int e = 0; e < nEvents; e++) {
            if (isImm[src[e]]) {
                continue;
            }
            if (!isImm[dst[e]]) {
                double[] col = new double[nStates];
                for (int s = 0; s < nStates; s++) {
                    col[s] = allJumps.get(s, e);
                }
                int row = jumpsNew.size();
                jumpsNew.add(col);
                rateNew.add(rateBase.get(e, 0));
                evidxNew.add(src[e]);
                emapIdx.add(new int[]{row, e});
                emapVal.add(1.0);
                continue;
            }
            // The event feeds an immediate coordinate: replace it by one event per absorbing
            // destination, keeping the original source and so the original gating, since the rate of
            // the composed flow IS the rate of the inflow.
            int a = posF[dst[e]];
            for (int b = 0; b < nS; b++) {
                double pa = Aabs.get(a, b);
                if (pa <= 0) {
                    continue;
                }
                int s = Slist.get(b);
                double[] col = new double[nStates];
                col[src[e]] -= 1;
                col[s] += 1;
                int row = jumpsNew.size();
                jumpsNew.add(col);
                rateNew.add(rateBase.get(e, 0) * pa);
                evidxNew.add(src[e]);
                // Weighting every absorbing branch by the SAME unconditional expected counts is what
                // makes the rate accounting exact: the branch rates sum back to rateBase(e), so the
                // mapped total is rateBase(e) times the counts.
                emapIdx.add(new int[]{row, e});
                emapVal.add(1.0);
                for (int o = 0; o < nEvents; o++) {
                    double v = expCnt.get(a, o);
                    if (v != 0) {
                        emapIdx.add(new int[]{row, o});
                        emapVal.add(v);
                    }
                }
            }
        }

        int nNew = jumpsNew.size();
        Matrix allJumpsReduced = new Matrix(nStates, nNew);
        Matrix rateBaseReduced = new Matrix(nNew, 1);
        Matrix eventIdxReduced = new Matrix(nNew, 1);
        for (int e = 0; e < nNew; e++) {
            double[] col = jumpsNew.get(e);
            for (int s = 0; s < nStates; s++) {
                if (col[s] != 0) {
                    allJumpsReduced.set(s, e, col[s]);
                }
            }
            rateBaseReduced.set(e, 0, rateNew.get(e));
            eventIdxReduced.set(e, 0, evidxNew.get(e));
        }
        Matrix Emap = new Matrix(nNew, nEvents);
        for (int i = 0; i < emapIdx.size(); i++) {
            int[] ij = emapIdx.get(i);
            Emap.set(ij[0], ij[1], Emap.get(ij[0], ij[1]) + emapVal.get(i));
        }

        Matrix absorb = Matrix.eye(nStates);
        for (int a = 0; a < nF; a++) {
            int f = Flist.get(a);
            for (int j = 0; j < nStates; j++) {
                absorb.set(f, j, 0);
            }
            for (int b = 0; b < nS; b++) {
                if (Aabs.get(a, b) > 0) {
                    absorb.set(f, Slist.get(b), Aabs.get(a, b));
                }
            }
        }

        int[] stateMap = new int[nS];
        for (int b = 0; b < nS; b++) {
            stateMap[b] = Slist.get(b);
        }

        if (options.verbose == VerboseLevel.DEBUG) {
            System.out.println(String.format(
                    "Eliminated %d immediate coordinates of %d, %d events of %d",
                    nF, nStates, nEvents - nNew, nEvents));
        }
        return new EliminationResult(allJumpsReduced, rateBaseReduced, eventIdxReduced,
                stateMap, Emap, absorb);
    }
}
