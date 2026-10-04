/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.lqn;

import java.util.Iterator;
import jline.GlobalConstants;
import jline.lang.constant.CallType;
import jline.lang.layered.ActivityPrecedence;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.LayeredNetworkStruct;
import jline.lang.layered.Task;
import jline.lang.processes.Immediate;
import jline.lang.workflow.Workflow;
import jline.lang.workflow.WorkflowActivity;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixEntry;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Phase-type composition of an LQN activity graph, the machinery behind
 * SolverLN method 'srvn.ph'.
 * <p>
 * An entry becomes a {@link Workflow} whose leaves are its activities and, when
 * requested, its synchronous calls; the series-parallel reduction of that
 * workflow is then the exact law of the entry service time. Mirrors the MATLAB
 * lqn_entry_workflow, lqn_ph_serial_law and lqn_ph_moments.
 * </p>
 */
public class LqnPh {

    private LqnPh() {
    }

    /** Activity graph of one entry, as a workflow plus its execution counts. */
    public static class EntryWorkflow {
        /** Workflow whose leaves are the activities (and calls) of the entry. */
        public Workflow wf;
        /** Workflow activity index of each LQN activity, -1 when absent. */
        public int[] actIdxOf;
        /** Workflow activity index of each call, -1 when absent. */
        public int[] callIdxOf;
        /** Expected executions of each LQN activity per entry invocation. */
        public double[] execs;
        /** Expected executions of each call per entry invocation. */
        public double[] callexecs;
    }

    /**
     * Activity graph of LQN entry EIDX as a workflow.
     * <p>
     * The activity precedences are read from the task object rather than
     * reconstructed from lqn.graph, whose loop back-edges carry probabilities
     * and not counts.
     * </p>
     *
     * @param model     the layered network the struct was obtained from
     * @param lqn       the layered network struct
     * @param eidx      absolute index of the entry
     * @param withCalls true to expand every synchronous call into a leaf of its
     *                  own, placed in series after its activity; false to keep
     *                  only the host demands, which is the processor-demand law
     * @return the workflow and its execution counts
     */
    public static EntryWorkflow entryWorkflow(LayeredNetwork model, LayeredNetworkStruct lqn,
                                              int eidx, boolean withCalls) {
        int tidx = (int) lqn.parent.get(0, eidx);
        Task task = taskNamed(model, lqn.names.get(tidx));
        if (task == null) {
            throw new IllegalStateException("Entry " + lqn.hashnames.get(eidx)
                    + " has no task in the layered model.");
        }
        List<Integer> acts = lqn.actsof.get(eidx);
        if (acts == null || acts.isEmpty()) {
            throw new IllegalStateException("Entry " + lqn.hashnames.get(eidx) + " binds no activity.");
        }

        EntryWorkflow out = new EntryWorkflow();
        out.wf = new Workflow(lqn.names.get(eidx) + ".Workflow");
        out.actIdxOf = new int[lqn.nidx];
        out.callIdxOf = new int[lqn.ncalls];
        java.util.Arrays.fill(out.actIdxOf, -1);
        java.util.Arrays.fill(out.callIdxOf, -1);

        Map<Integer, String> headName = new HashMap<Integer, String>();
        Map<Integer, String> tailName = new HashMap<Integer, String>();
        Map<String, Integer> actOfName = new HashMap<String, Integer>();

        for (int aidx : acts) {
            String nm = lqn.names.get(aidx);
            WorkflowActivity a = out.wf.addActivity(nm, lqn.hostdem.get(aidx));
            out.actIdxOf[aidx] = a.getIndex();
            headName.put(aidx, nm);
            tailName.put(aidx, nm);
            actOfName.put(nm, aidx);
        }

        if (withCalls) {
            for (int aidx : acts) {
                List<String> chain = new ArrayList<String>();
                chain.add(headName.get(aidx));
                List<Integer> calls = lqn.callsof.get(aidx);
                if (calls != null) {
                    for (int cidx : calls) {
                        if (lqn.calltype.get(cidx) != CallType.SYNC) {
                            continue; // an asynchronous call blocks the caller for no time
                        }
                        String cnm = lqn.callhashnames.get(cidx);
                        WorkflowActivity c = out.wf.addActivity(cnm, Immediate.getInstance());
                        out.callIdxOf[cidx] = c.getIndex();
                        chain.add(cnm);
                    }
                }
                if (chain.size() > 1) {
                    for (int k = 1; k < chain.size(); k++) {
                        out.wf.addPrecedence(ActivityPrecedence.Serial(chain.get(k - 1), chain.get(k)));
                    }
                    tailName.put(aidx, chain.get(chain.size() - 1));
                }
            }
        }

        // Precedences of the task, restricted to the activities of this entry and
        // rewritten so that a predecessor is entered at its head and left at its tail
        for (ActivityPrecedence prec : task.getPrecedences()) {
            List<Integer> preIdx = resolveNames(prec.getPreActs(), actOfName);
            List<Integer> postIdx = resolveNames(prec.getPostActs(), actOfName);
            if (preIdx == null || postIdx == null) {
                continue; // the precedence belongs to another entry of the same task
            }
            List<String> preNames = new ArrayList<String>();
            for (int k : preIdx) {
                preNames.add(tailName.get(k));
            }
            List<String> postNames = new ArrayList<String>();
            for (int k : postIdx) {
                postNames.add(headName.get(k));
            }
            out.wf.addPrecedence(new ActivityPrecedence(preNames, postNames,
                    prec.getPreType(), prec.getPostType(), prec.getPreParams(), prec.getPostParams()));
        }

        Pair<Boolean, String> valid = out.wf.validate();
        if (!valid.getLeft()) {
            throw new IllegalStateException("Entry " + lqn.hashnames.get(eidx)
                    + " cannot be composed into a phase-type law: " + valid.getRight());
        }
        Workflow.SPTree tree = out.wf.getSPTree();
        if (tree == null) {
            throw new IllegalStateException("Entry " + lqn.hashnames.get(eidx)
                    + " has a precedence graph that is not series-parallel, so its activity graph has no "
                    + "exact phase-type reduction. Use method='default'.");
        }

        out.execs = new double[lqn.nidx];
        out.callexecs = new double[lqn.ncalls];
        for (int aidx : acts) {
            out.execs[aidx] = tree.execs[tree.leafOf[out.actIdxOf[aidx]]];
        }
        for (int cidx = 0; cidx < lqn.ncalls; cidx++) {
            if (out.callIdxOf[cidx] >= 0) {
                out.callexecs[cidx] = tree.execs[tree.leafOf[out.callIdxOf[cidx]]];
            }
        }
        return out;
    }

    /**
     * Composed law of a workflow in which the branches of an AND fork are
     * SERIAL rather than concurrent, that is, the total work the branches
     * request rather than the elapsed time until the last of them finishes.
     * <p>
     * This is the law of the PROCESSOR demand of an LQN entry. Two branches of
     * an AND fork are two activity threads of the same task instance: they
     * overlap in time, so the entry response time is the maximum of the
     * branches, but they run on ONE processor, so the demand they place on it
     * is the sum. Composing the host law with Workflow.toPH would charge the
     * processor the maximum and let the layer report a utilization below the
     * true one, which no amount of iterating recovers.
     * </p>
     *
     * @param wf the workflow
     * @return the initial vector and subgenerator of the composed law
     */
    public static Pair<Matrix, Matrix> serialLaw(Workflow wf) {
        Workflow.SPTree tree = wf.getSPTree();
        if (tree == null) {
            throw new IllegalStateException("Workflow " + wf.getName()
                    + " is not series-parallel, so it has no exact phase-type reduction.");
        }
        return composeSerialized(wf, tree, tree.root);
    }

    private static Pair<Matrix, Matrix> composeSerialized(Workflow wf, Workflow.SPTree tree, int k) {
        Workflow.SPNode node = tree.nodes.get(k);
        int[] kids = node.kids;
        if ("leaf".equals(node.type)) {
            return wf.getActivities().get(node.act).getPHRepresentation();
        } else if ("serial".equals(node.type) || "par".equals(node.type)) {
            Pair<Matrix, Matrix> acc = composeSerialized(wf, tree, kids[0]);
            for (int i = 1; i < kids.length; i++) {
                Pair<Matrix, Matrix> nx = composeSerialized(wf, tree, kids[i]);
                acc = Workflow.composeSerial(acc.getLeft(), acc.getRight(), nx.getLeft(), nx.getRight());
            }
            return acc;
        } else if ("or".equals(node.type)) {
            List<Matrix> alphas = new ArrayList<Matrix>();
            List<Matrix> Ts = new ArrayList<Matrix>();
            for (int i = 0; i < kids.length; i++) {
                Pair<Matrix, Matrix> br = composeSerialized(wf, tree, kids[i]);
                alphas.add(br.getLeft());
                Ts.add(br.getRight());
            }
            return Workflow.composeMixture(alphas, Ts, node.probs);
        } else if ("loop".equals(node.type)) {
            Pair<Matrix, Matrix> body = composeSerialized(wf, tree, kids[0]);
            return Workflow.composeLoopGeometric(body.getLeft(), body.getRight(), node.count);
        }
        throw new IllegalStateException("Unknown series-parallel node type \"" + node.type + "\".");
    }

    /**
     * First two moments of a phase-type law without building a Distribution
     * object, which is what the layered fixed point needs at every iteration
     * for every composed entry law. A defective ALPHA carries an atom at zero
     * and contributes nothing to either moment.
     *
     * @param alpha initial vector
     * @param T     subgenerator
     * @return {mean, squared coefficient of variation}
     */
    public static double[] moments(Matrix alpha, Matrix T) {
        int n = T.getNumRows();
        Matrix Tinv = T.inv();
        // Tinv is read n^2 times below, twice over, and every read of a sparse-backed
        // matrix is a binary search in a CSC column. Materialise it row-major once.
        // The accumulation order below is UNCHANGED -- reassociating these sums would
        // move the last bit -- and an absent CSC entry reads as the same 0.0 that
        // Tinv.get returns, so every partial sum is the same double as before.
        double[][] tinv = new double[n][n];
        Iterator<MatrixEntry> tit = Tinv.nonZeroIterator();
        while (tit.hasNext()) {
            MatrixEntry e = tit.next();
            tinv[e.row][e.col] = e.value;
        }
        double[] x1 = new double[n];
        for (int i = 0; i < n; i++) {
            double s = 0;
            double[] ti = tinv[i];
            for (int j = 0; j < n; j++) {
                s += ti[j];
            }
            x1[i] = -s;
        }
        double m1 = 0;
        for (int i = 0; i < n; i++) {
            m1 += alphaAt(alpha, i) * x1[i];
        }
        double[] x2 = new double[n];
        for (int i = 0; i < n; i++) {
            double s = 0;
            double[] ti = tinv[i];
            for (int j = 0; j < n; j++) {
                s += ti[j] * x1[j];
            }
            x2[i] = -s;
        }
        double m2 = 0;
        for (int i = 0; i < n; i++) {
            m2 += alphaAt(alpha, i) * x2[i];
        }
        m2 = 2 * m2;
        if (Double.isNaN(m1) || Double.isInfinite(m1) || m1 <= GlobalConstants.FineTol) {
            return new double[]{GlobalConstants.FineTol, 1.0};
        }
        double scv = m2 / (m1 * m1) - 1;
        if (Double.isNaN(scv) || Double.isInfinite(scv) || scv <= GlobalConstants.FineTol) {
            scv = GlobalConstants.FineTol;
        }
        return new double[]{m1, scv};
    }

    private static double alphaAt(Matrix alpha, int i) {
        if (alpha.getNumRows() == 1) {
            return alpha.get(0, i);
        }
        return alpha.get(i, 0);
    }

    private static Task taskNamed(LayeredNetwork model, String name) {
        for (Task t : model.getTasks().values()) {
            if (t.getName().equals(name)) {
                return t;
            }
        }
        return null;
    }

    private static List<Integer> resolveNames(List<String> names, Map<String, Integer> actOfName) {
        List<Integer> out = new ArrayList<Integer>();
        for (String nm : names) {
            Integer v = actOfName.get(nm);
            if (v == null) {
                return null;
            }
            out.add(v);
        }
        return out;
    }
}
