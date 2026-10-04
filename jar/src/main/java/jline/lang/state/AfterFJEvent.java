/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.state;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.lang.FJSync;
import jline.lang.NetworkStruct;
import jline.lang.constant.EventType;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import static jline.util.Maths.rand;

/**
 * Fork firing synchronization handler (FJ tag-augmented structs, see
 * ModelAdapter.fjtag). The firing atomically consumes one parent job of
 * the entry's class held at the (stateful) Fork node and emits weight
 * (tasksPerLink) siblings per branch, in the auxiliary classes of the
 * entry's tag, at the branch head nodes.
 *
 * Enabling condition (evaluated on the global state): the Fork holds at
 * least one parent job, and the entry's tag is the LOWEST free tag for
 * this (fork, class): a tag is free iff its auxiliary classes have zero
 * occupancy network-wide. Canonical lowest-free-tag allocation ensures
 * exactly one fjsync entry per (fork, class) is enabled in any state.
 */
public class AfterFJEvent implements Serializable {

    public static class AfterFJEventResult {
        public final List<List<Matrix>> outGlobalStates;
        public final Matrix outrate;
        public final Matrix outprob;

        public AfterFJEventResult(List<List<Matrix>> outGlobalStates, Matrix outrate, Matrix outprob) {
            this.outGlobalStates = outGlobalStates;
            this.outrate = outrate;
            this.outprob = outprob;
        }
    }

    public static AfterFJEventResult afterFJEvent(NetworkStruct sn, FJSync fjentry, List<Matrix> glspace, boolean isSimulation, EventCache eventCache) {
        int R = sn.nclasses;
        int f = fjentry.fork;
        int r = fjentry.jobclass;
        int isfF = (int) sn.nodeToStateful.get(f);
        List<List<Matrix>> empty = new ArrayList<List<Matrix>>();

        // 1. parent job held at the fork
        Matrix forkstate = glspace.get(isfF);
        if (forkstate.get(0, forkstate.getNumCols() - R + r) < 1) {
            return new AfterFJEventResult(empty, new Matrix(0, 0), new Matrix(0, 0));
        }

        // 2. lowest-free-tag test: occupancy of each tag's auxiliary classes,
        // scanned across all stateful nodes
        Matrix auxall = fjentry.auxall; // B x T
        int T = auxall.getNumCols();
        int t = fjentry.tag;
        double[] nglobal = new double[R];
        for (int isf = 0; isf < sn.nstateful; isf++) {
            int ind = (int) sn.statefulToNode.get(isf);
            State.StateMarginalStatistics stats = ToMarginal.toMarginalAggr(sn, ind, glspace.get(isf), null, null, null, null, null);
            for (int rr = 0; rr < Math.min(R, stats.nir.getNumCols()); rr++) {
                nglobal[rr] += stats.nir.get(0, rr);
            }
        }
        double[] occ = new double[T];
        for (int tt = 0; tt < T; tt++) {
            for (int b = 0; b < auxall.getNumRows(); b++) {
                occ[tt] += nglobal[(int) auxall.get(b, tt)];
            }
        }
        if (occ[t] > 0) {
            return new AfterFJEventResult(empty, new Matrix(0, 0), new Matrix(0, 0)); // tag in use
        }
        for (int tt = 0; tt < t; tt++) {
            if (occ[tt] == 0) {
                return new AfterFJEventResult(empty, new Matrix(0, 0), new Matrix(0, 0)); // a lower tag is free
            }
        }

        // consume the parent job at the fork
        List<Matrix> newgl = new ArrayList<Matrix>();
        for (Matrix m : glspace) {
            newgl.add(m.copy());
        }
        Matrix fk = newgl.get(isfF);
        fk.set(0, fk.getNumCols() - R + r, fk.get(0, fk.getNumCols() - R + r) - 1);

        // emit weight (tasksPerLink) siblings per branch: sequential
        // application over the partial outcome list handles branches sharing
        // the same head node, repeated emissions, and phase-entry mixtures of
        // non-exponential sibling services
        int B = fjentry.branchheads.length;
        int weight = (int) fjentry.weight;
        List<List<Matrix>> partials = new ArrayList<List<Matrix>>();
        partials.add(newgl);
        List<Double> partprob = new ArrayList<Double>();
        partprob.add(1.0);

        for (int w = 0; w < weight; w++) {
            for (int b = 0; b < B; b++) {
                int bh = fjentry.branchheads[b];
                int isfB = (int) sn.nodeToStateful.get(bh);
                int a = fjentry.auxclasses[b];
                List<List<Matrix>> newpartials = new ArrayList<List<Matrix>>();
                List<Double> newpartprob = new ArrayList<Double>();
                for (int pp = 0; pp < partials.size(); pp++) {
                    List<Matrix> curgl = partials.get(pp);
                    Ret.EventResult arv = State.afterEvent(sn, bh, curgl.get(isfB), EventType.ARV, a, isSimulation, eventCache);
                    if (arv.outspace.isEmpty()) {
                        // sibling arrival blocked (cannot occur under the per-tag
                        // auxiliary class capacity invariant); disable the firing
                        return new AfterFJEventResult(empty, new Matrix(0, 0), new Matrix(0, 0));
                    }
                    for (int io = 0; io < arv.outspace.getNumRows(); io++) {
                        List<Matrix> nextgl = new ArrayList<Matrix>();
                        for (Matrix m : curgl) {
                            nextgl.add(m.copy());
                        }
                        nextgl.set(isfB, Matrix.extractRows(arv.outspace, io, io + 1, null));
                        newpartials.add(nextgl);
                        double pio = (io < arv.outprob.getNumRows()) ? arv.outprob.get(io, 0) : 1.0;
                        newpartprob.add(partprob.get(pp) * pio);
                    }
                }
                partials = newpartials;
                partprob = newpartprob;
            }
        }

        Matrix outprob = new Matrix(partials.size(), 1);
        for (int i = 0; i < partials.size(); i++) {
            outprob.set(i, 0, partprob.get(i) * fjentry.prob);
        }
        Matrix outrate = new Matrix(partials.size(), 1);
        outrate.fill(GlobalConstants.Immediate);

        if (isSimulation) {
            if (partials.size() > 1) {
                double tot = outprob.elementSum();
                double u = rand() * (tot > 0 ? tot : 1.0);
                double cum = 0.0;
                int sel = partials.size() - 1;
                for (int i = 0; i < partials.size(); i++) {
                    cum += outprob.get(i, 0);
                    if (u <= cum) {
                        sel = i;
                        break;
                    }
                }
                List<List<Matrix>> selected = new ArrayList<List<Matrix>>();
                selected.add(partials.get(sel));
                partials = selected;
                Matrix rr1 = new Matrix(1, 1);
                rr1.fill(GlobalConstants.Immediate);
                outrate = rr1;
            }
            // the phase-entry choice has already been sampled inside
            // afterEvent, so the firing competes at the full immediate rate
            outprob = new Matrix(partials.size(), 1);
            outprob.ones();
        }

        return new AfterFJEventResult(partials, outrate, outprob);
    }
}
