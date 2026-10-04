/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.workflow;

import jline.GlobalConstants;
import jline.io.WfCommonsLoader;
import jline.io.WfCommonsOptions;
import jline.lang.Model;
import jline.lang.constant.ActivityPrecedenceType;
import jline.lang.layered.ActivityPrecedence;
import jline.lang.processes.APH;
import jline.lang.processes.Distribution;
import jline.lang.processes.Markovian;
import jline.lang.processes.PH;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixEntry;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

/**
 * A computational workflow that can be converted to a phase-type distribution.
 * <p>
 * A workflow whose precedence graph is series-parallel is reduced exactly, by
 * recursive composition of the series-parallel tree, which handles arbitrary
 * nesting (a fork inside a loop, a branch that is itself a fork-join). Graphs
 * that are not series-parallel fall back to the block-based composition, which
 * is a heuristic.
 * </p>
 * <p>
 * A loop repeats its body a geometric number of times of mean COUNT, the
 * semantics of the POST_LOOP precedence of an activity graph, rather than
 * COUNT times deterministically.
 * </p>
 */
public class Workflow extends Model {

    private List<WorkflowActivity> activities;
    private List<ActivityPrecedence> precedences;
    private Map<String, Integer> activityMap;
    private Markovian cachedPH;
    private SPTree spTree;
    private boolean spFailed;

    public Workflow(String name) {
        super(name);
        this.activities = new ArrayList<WorkflowActivity>();
        this.precedences = new ArrayList<ActivityPrecedence>();
        this.activityMap = new HashMap<String, Integer>();
        this.cachedPH = null;
        this.spTree = null;
        this.spFailed = false;
    }

    public WorkflowActivity addActivity(String name, double meanServiceTime) {
        WorkflowActivity act = new WorkflowActivity(this, name, meanServiceTime);
        activities.add(act);
        act.setIndex(activities.size() - 1);
        activityMap.put(name, act.getIndex());
        invalidateTopology();
        return act;
    }

    public WorkflowActivity addActivity(String name, Distribution hostDemand) {
        WorkflowActivity act = new WorkflowActivity(this, name, hostDemand);
        activities.add(act);
        act.setIndex(activities.size() - 1);
        activityMap.put(name, act.getIndex());
        invalidateTopology();
        return act;
    }

    public void addPrecedence(ActivityPrecedence prec) {
        precedences.add(prec);
        invalidateTopology();
    }

    public void addPrecedence(ActivityPrecedence[] precs) {
        for (ActivityPrecedence prec : precs) {
            precedences.add(prec);
        }
        invalidateTopology();
    }

    public WorkflowActivity getActivity(String name) {
        Integer idx = activityMap.get(name);
        if (idx == null) {
            return null;
        }
        return activities.get(idx);
    }

    public List<WorkflowActivity> getActivities() {
        return activities;
    }

    public List<ActivityPrecedence> getPrecedences() {
        return precedences;
    }

    public Pair<Boolean, String> validate() {
        if (activities.isEmpty()) {
            return new Pair<Boolean, String>(false, "Workflow must have at least one activity.");
        }

        for (ActivityPrecedence prec : precedences) {
            for (String actName : prec.getPreActs()) {
                if (!activityMap.containsKey(actName)) {
                    return new Pair<Boolean, String>(false,
                            "Activity '" + actName + "' referenced in precedence not found.");
                }
            }
            for (String actName : prec.getPostActs()) {
                if (!activityMap.containsKey(actName)) {
                    return new Pair<Boolean, String>(false,
                            "Activity '" + actName + "' referenced in precedence not found.");
                }
            }
        }

        for (ActivityPrecedence prec : precedences) {
            if (ActivityPrecedenceType.POST_OR.equals(prec.getPostType())) {
                if (prec.getPostParams() == null) {
                    return new Pair<Boolean, String>(false, "OR-fork must have probabilities.");
                }
                double sum = 0;
                for (int i = 0; i < prec.getPostParams().getNumElements(); i++) {
                    sum += prec.getPostParams().get(i);
                }
                if (Math.abs(sum - 1.0) > GlobalConstants.FineTol) {
                    return new Pair<Boolean, String>(false, "OR-fork probabilities must sum to 1.");
                }
            }
            if (ActivityPrecedenceType.POST_LOOP.equals(prec.getPostType())) {
                Matrix counts = prec.getPostParams();
                if (counts == null || counts.getNumElements() != 1) {
                    return new Pair<Boolean, String>(false, "Loop count must be a single positive number.");
                }
                if (counts.get(0) <= 0) {
                    return new Pair<Boolean, String>(false, "Loop count must be a positive number.");
                }
            }
            // A partial (quorum) AND-join is refused rather than silently
            // served as a full join, which would be a different law
            if (ActivityPrecedenceType.PRE_AND.equals(prec.getPreType()) && prec.getPreParams() != null
                    && prec.getPreParams().getNumElements() > 0) {
                double quorum = prec.getPreParams().get(0);
                if (quorum > 0 && quorum < prec.getPreActs().size()) {
                    return new Pair<Boolean, String>(false,
                            "AND-join with quorum " + (int) quorum + " of " + prec.getPreActs().size()
                                    + " is not supported by Workflow: a partial join is not the maximum of the "
                                    + "branches. Use a full join, or SolverLN with method='default', which routes "
                                    + "the join explicitly.");
                }
            }
        }

        return new Pair<Boolean, String>(true, "");
    }

    /**
     * Composed law of the workflow execution time.
     *
     * @return an APH, or a PH when the composed generator is cyclic, which a
     *         geometric loop over a multi-phase body makes it
     */
    public Markovian toPH() {
        if (cachedPH != null) {
            return cachedPH;
        }

        Pair<Boolean, String> validation = validate();
        if (!validation.getLeft()) {
            throw new IllegalStateException(validation.getRight());
        }

        Pair<Matrix, Matrix> law = composeSeriesParallel();
        if (law == null) {
            law = buildCTMC();
        }

        if (isAcyclicGenerator(law.getRight())) {
            cachedPH = new APH(law.getLeft(), law.getRight());
        } else {
            cachedPH = new PH(law.getLeft(), law.getRight());
        }
        return cachedPH;
    }

    /**
     * Recompose the workflow law after a demand change, reusing every cached
     * series-parallel node whose subtree is unchanged.
     *
     * @return the composed law
     */
    public Markovian refreshPH() {
        return toPH();
    }

    /**
     * Change the host demand of one activity, marking only that leaf dirty.
     *
     * @param name       activity name
     * @param hostDemand new host demand law
     */
    public void setActivityDemand(String name, Distribution hostDemand) {
        WorkflowActivity act = getActivity(name);
        if (act == null) {
            throw new IllegalArgumentException("Activity '" + name + "' not found in workflow.");
        }
        act.setHostDemand(hostDemand);
    }

    /**
     * Change only the mean of one activity, preserving its shape and order.
     *
     * @param name      activity name
     * @param meanValue new mean
     */
    public void setActivityDemandMean(String name, double meanValue) {
        WorkflowActivity act = getActivity(name);
        if (act == null) {
            throw new IllegalArgumentException("Activity '" + name + "' not found in workflow.");
        }
        act.setHostDemandMean(meanValue);
    }

    /**
     * Discard the cached law and decomposition, after a change that can alter
     * the shape of the series-parallel tree.
     */
    public void invalidateTopology() {
        cachedPH = null;
        spTree = null;
        spFailed = false;
    }

    /**
     * Mark the law of one activity dirty, keeping the topology and every other
     * cached block.
     *
     * @param actIdx activity index
     */
    public void invalidateActivity(int actIdx) {
        cachedPH = null;
        if (spTree == null) {
            return;
        }
        if (actIdx < 0 || actIdx >= spTree.leafOf.length || spTree.leafOf[actIdx] < 0) {
            spTree = null;
            spFailed = false;
            return;
        }
        invalidateBranch(spTree, spTree.leafOf[actIdx]);
    }

    /**
     * Time-scale a cached leaf in place, as T -&gt; T*factor with alpha fixed.
     *
     * @param actIdx activity index
     * @param factor rate scaling factor
     */
    public void rescaleActivityLeaf(int actIdx, double factor) {
        cachedPH = null;
        if (spTree == null) {
            return;
        }
        if (actIdx < 0 || actIdx >= spTree.leafOf.length || spTree.leafOf[actIdx] < 0) {
            spTree = null;
            spFailed = false;
            return;
        }
        int k = spTree.leafOf[actIdx];
        invalidateBranch(spTree, k);
        SPNode leaf = spTree.nodes.get(k);
        if (leaf.T != null) {
            leaf.T = leaf.T.scale(factor);
            leaf.valid = true;
        }
    }

    /**
     * Cached series-parallel decomposition, or null when the precedence graph
     * is not series-parallel. Field execs carries the expected number of
     * executions of each node per workflow execution.
     *
     * @return the decomposition
     */
    public SPTree getSPTree() {
        if (spTree == null && !spFailed) {
            buildSPTree();
        }
        return spTree;
    }

    private Pair<Matrix, Matrix> buildCTMC() {
        int n = activities.size();

        if (n == 1) {
            return activities.get(0).getPHRepresentation();
        }

        WorkflowStructure structure = analyzeStructure();

        if (structure.forkInfo.isEmpty() && structure.joinInfo.isEmpty() && structure.loopInfo.isEmpty()) {
            return composeSerialWorkflow(structure.adjList);
        }

        return composeComplexWorkflow(structure);
    }

    /** A node of the series-parallel decomposition. */
    public static class SPNode {
        /** One of leaf, serial, par, or, loop. */
        public String type;
        /** Activity index of a leaf, -1 otherwise. */
        public int act = -1;
        /** Child node indices. */
        public int[] kids = new int[0];
        /** Parent node index, -1 at the root. */
        public int parent = -1;
        /** Branch probabilities of an or node. */
        public double[] probs;
        /** Mean number of executions of a loop body. */
        public double count;
        /** Cached initial vector. */
        public Matrix alpha;
        /** Cached subgenerator. */
        public Matrix T;
        /** True when the cached law is up to date. */
        public boolean valid;
    }

    /** The series-parallel decomposition of the precedence graph. */
    public static class SPTree {
        /** Nodes, in creation order. */
        public List<SPNode> nodes = new ArrayList<SPNode>();
        /** Root node index. */
        public int root = -1;
        /** Node index of each activity leaf, -1 when absent. */
        public int[] leafOf;
        /** Expected executions of each node per workflow execution. */
        public double[] execs;
    }

    private static class SPParse {
        int[] outP;
        int[] inP;
        boolean[] consumed;
        List<SPNode> nodes = new ArrayList<SPNode>();
        boolean ok = true;
        String status = "end";
        int stopAt = -1;
    }

    private Pair<Matrix, Matrix> composeSeriesParallel() {
        if (spTree == null) {
            if (spFailed) {
                return null;
            }
            buildSPTree();
            if (spTree == null) {
                return null;
            }
        }
        return composeNode(spTree.root);
    }

    private boolean buildSPTree() {
        spTree = null;
        spFailed = true;

        int n = activities.size();
        if (n == 0) {
            return false;
        }

        SPParse S = new SPParse();
        S.outP = new int[n];
        S.inP = new int[n];
        S.consumed = new boolean[n];
        Arrays.fill(S.outP, -1);
        Arrays.fill(S.inP, -1);

        // An activity may head at most one precedence and be reached by at
        // most one precedence; otherwise the graph is not series-parallel
        for (int p = 0; p < precedences.size(); p++) {
            ActivityPrecedence prec = precedences.get(p);
            for (String name : prec.getPreActs()) {
                Integer i = activityMap.get(name);
                if (i == null || S.outP[i] >= 0) {
                    return false;
                }
                S.outP[i] = p;
            }
            for (String name : prec.getPostActs()) {
                Integer j = activityMap.get(name);
                if (j == null || S.inP[j] >= 0) {
                    return false;
                }
                S.inP[j] = p;
            }
        }

        int start = -1;
        for (int i = 0; i < n; i++) {
            if (S.inP[i] < 0) {
                if (start >= 0) {
                    return false;
                }
                start = i;
            }
        }
        if (start < 0) {
            return false;
        }

        List<Integer> kids = spParseSeq(S, start, new HashSet<Integer>());
        if (!S.ok || !"end".equals(S.status)) {
            return false;
        }
        for (int i = 0; i < n; i++) {
            if (!S.consumed[i]) {
                return false;
            }
        }

        int root = spSerialNode(S, kids);
        if (root < 0) {
            return false;
        }

        SPTree tree = new SPTree();
        tree.nodes = S.nodes;
        tree.root = root;
        tree.leafOf = new int[n];
        Arrays.fill(tree.leafOf, -1);
        for (int k = 0; k < tree.nodes.size(); k++) {
            SPNode node = tree.nodes.get(k);
            if ("leaf".equals(node.type)) {
                tree.leafOf[node.act] = k;
            }
        }
        tree.execs = spExecutionCounts(tree);

        spTree = tree;
        spFailed = false;
        return true;
    }

    private List<Integer> spParseSeq(SPParse S, int cur, Set<Integer> stopSet) {
        List<Integer> kids = new ArrayList<Integer>();
        S.status = "end";
        S.stopAt = -1;

        while (true) {
            if (cur < 0) {
                S.status = "end";
                return kids;
            }
            if (stopSet.contains(cur)) {
                S.status = "stop";
                S.stopAt = cur;
                return kids;
            }
            if (S.consumed[cur]) {
                S.ok = false;
                return kids;
            }
            S.consumed[cur] = true;
            kids.add(spAddNode(S, "leaf", cur, new int[0], null, 0));

            int p = S.outP[cur];
            if (p < 0) {
                S.status = "end";
                return kids;
            }
            ActivityPrecedence prec = precedences.get(p);
            if (prec.getPreActs().size() > 1) {
                // CUR is the tail of a branch: the caller composes the join
                S.status = "join";
                S.stopAt = p;
                return kids;
            }

            int[] postInds = spIndicesOf(prec.getPostActs());
            if (postInds == null) {
                S.ok = false;
                return kids;
            }

            String postType = prec.getPostType();
            if (ActivityPrecedenceType.POST_AND.equals(postType)) {
                int knode = spParseFork(S, postInds, null, stopSet, true);
                if (!S.ok) {
                    return kids;
                }
                kids.add(knode);
                cur = S.stopAt;
            } else if (ActivityPrecedenceType.POST_OR.equals(postType)) {
                Matrix params = prec.getPostParams();
                if (params == null || params.getNumElements() != postInds.length) {
                    S.ok = false;
                    return kids;
                }
                double[] probs = new double[postInds.length];
                for (int i = 0; i < postInds.length; i++) {
                    probs[i] = params.get(i);
                }
                int knode = spParseFork(S, postInds, probs, stopSet, false);
                if (!S.ok) {
                    return kids;
                }
                kids.add(knode);
                cur = S.stopAt;
            } else if (ActivityPrecedenceType.POST_LOOP.equals(postType)) {
                Matrix params = prec.getPostParams();
                if (params == null || params.getNumElements() != 1) {
                    S.ok = false;
                    return kids;
                }
                int knode = spParseLoop(S, postInds, params.get(0), stopSet);
                if (!S.ok) {
                    return kids;
                }
                kids.add(knode);
                cur = S.stopAt;
            } else if (ActivityPrecedenceType.POST_SEQ.equals(postType)) {
                if (postInds.length != 1) {
                    S.ok = false;
                    return kids;
                }
                cur = postInds[0];
            } else {
                // POST_CACHE and any other pattern is not a composition rule
                S.ok = false;
                return kids;
            }
        }
    }

    private int spParseFork(SPParse S, int[] branchHeads, double[] probs, Set<Integer> stopSet, boolean isAnd) {
        int nb = branchHeads.length;
        int[] branchNodes = new int[nb];
        String[] bstatus = new String[nb];
        int[] bstop = new int[nb];

        for (int b = 0; b < nb; b++) {
            List<Integer> bkids = spParseSeq(S, branchHeads[b], stopSet);
            if (!S.ok) {
                return -1;
            }
            bstatus[b] = S.status;
            bstop[b] = S.stopAt;
            int bn = spSerialNode(S, bkids);
            if (bn < 0) {
                S.ok = false;
                return -1;
            }
            branchNodes[b] = bn;
        }

        boolean allJoin = true;
        boolean allEnd = true;
        boolean allStop = true;
        for (int b = 0; b < nb; b++) {
            allJoin = allJoin && "join".equals(bstatus[b]) && bstop[b] == bstop[0];
            allEnd = allEnd && "end".equals(bstatus[b]);
            allStop = allStop && "stop".equals(bstatus[b]) && bstop[b] == bstop[0];
        }

        int nextAct;
        if (allJoin) {
            ActivityPrecedence joinPrec = precedences.get(bstop[0]);
            if (joinPrec.getPreActs().size() != nb) {
                S.ok = false;
                return -1;
            }
            String needed = isAnd ? ActivityPrecedenceType.PRE_AND : ActivityPrecedenceType.PRE_OR;
            if (!needed.equals(joinPrec.getPreType())) {
                S.ok = false;
                return -1;
            }
            int[] postInds = spIndicesOf(joinPrec.getPostActs());
            if (postInds == null || postInds.length != 1) {
                S.ok = false;
                return -1;
            }
            nextAct = postInds[0];
        } else if (allEnd) {
            nextAct = -1;
        } else if (!isAnd && allStop) {
            nextAct = bstop[0];
        } else {
            S.ok = false;
            return -1;
        }

        S.stopAt = nextAct;
        if (isAnd) {
            return spAddNode(S, "par", -1, branchNodes, null, 0);
        }
        return spAddNode(S, "or", -1, branchNodes, probs, 0);
    }

    private int spParseLoop(SPParse S, int[] postInds, double count, Set<Integer> stopSet) {
        int[] bodyActs;
        int endAct;
        if (postInds.length >= 2) {
            bodyActs = Arrays.copyOf(postInds, postInds.length - 1);
            endAct = postInds[postInds.length - 1];
        } else {
            bodyActs = postInds;
            endAct = -1;
        }

        Set<Integer> loopStop = new HashSet<Integer>(stopSet);
        for (int a : bodyActs) {
            loopStop.add(a);
        }
        if (endAct >= 0) {
            loopStop.add(endAct);
        }

        List<Integer> bodyKids = new ArrayList<Integer>();
        int j = 0;
        while (j < bodyActs.length) {
            int a = bodyActs[j];
            if (S.consumed[a]) {
                j++;
                continue;
            }
            Set<Integer> thisStop = new HashSet<Integer>(loopStop);
            thisStop.remove(a);
            List<Integer> kk = spParseSeq(S, a, thisStop);
            if (!S.ok) {
                return -1;
            }
            String st = S.status;
            int sa = S.stopAt;
            bodyKids.addAll(kk);
            if ("end".equals(st)) {
                j++;
            } else if ("stop".equals(st)) {
                int idx = -1;
                for (int i = 0; i < bodyActs.length; i++) {
                    if (bodyActs[i] == sa) {
                        idx = i;
                        break;
                    }
                }
                if (idx < 0) {
                    if (endAct >= 0 && sa == endAct) {
                        j = bodyActs.length;
                    } else {
                        S.ok = false;
                        return -1;
                    }
                } else {
                    j = idx;
                }
            } else {
                // A join reached from inside the body crosses the loop boundary
                S.ok = false;
                return -1;
            }
        }

        int bodyNode = spSerialNode(S, bodyKids);
        if (bodyNode < 0) {
            S.ok = false;
            return -1;
        }

        S.stopAt = endAct;
        return spAddNode(S, "loop", -1, new int[]{bodyNode}, null, count);
    }

    private int spSerialNode(SPParse S, List<Integer> kids) {
        if (kids.isEmpty()) {
            return -1;
        }
        if (kids.size() == 1) {
            return kids.get(0);
        }
        int[] arr = new int[kids.size()];
        for (int i = 0; i < kids.size(); i++) {
            arr[i] = kids.get(i);
        }
        return spAddNode(S, "serial", -1, arr, null, 0);
    }

    private int spAddNode(SPParse S, String type, int act, int[] kids, double[] probs, double count) {
        SPNode node = new SPNode();
        node.type = type;
        node.act = act;
        node.kids = kids;
        node.probs = probs;
        node.count = count;
        node.valid = false;
        S.nodes.add(node);
        int k = S.nodes.size() - 1;
        for (int c : kids) {
            S.nodes.get(c).parent = k;
        }
        return k;
    }

    private int[] spIndicesOf(List<String> names) {
        int[] idx = new int[names.size()];
        for (int a = 0; a < names.size(); a++) {
            Integer i = activityMap.get(names.get(a));
            if (i == null) {
                return null;
            }
            idx[a] = i;
        }
        return idx;
    }

    private Pair<Matrix, Matrix> composeNode(int k) {
        SPNode node = spTree.nodes.get(k);
        if (node.valid) {
            return new Pair<Matrix, Matrix>(node.alpha, node.T);
        }

        Matrix alpha;
        Matrix T;
        if ("leaf".equals(node.type)) {
            Pair<Matrix, Matrix> ph = activities.get(node.act).getPHRepresentation();
            alpha = ph.getLeft();
            T = ph.getRight();
        } else if ("serial".equals(node.type)) {
            Pair<Matrix, Matrix> acc = composeNode(node.kids[0]);
            for (int i = 1; i < node.kids.length; i++) {
                Pair<Matrix, Matrix> next = composeNode(node.kids[i]);
                acc = composeSerial(acc.getLeft(), acc.getRight(), next.getLeft(), next.getRight());
            }
            alpha = acc.getLeft();
            T = acc.getRight();
        } else if ("par".equals(node.type)) {
            Pair<Matrix, Matrix> acc = composeNode(node.kids[0]);
            for (int i = 1; i < node.kids.length; i++) {
                Pair<Matrix, Matrix> next = composeNode(node.kids[i]);
                acc = composeParallel(acc.getLeft(), acc.getRight(), next.getLeft(), next.getRight());
            }
            alpha = acc.getLeft();
            T = acc.getRight();
        } else if ("or".equals(node.type)) {
            List<Matrix> alphas = new ArrayList<Matrix>();
            List<Matrix> Ts = new ArrayList<Matrix>();
            for (int i = 0; i < node.kids.length; i++) {
                Pair<Matrix, Matrix> kid = composeNode(node.kids[i]);
                alphas.add(kid.getLeft());
                Ts.add(kid.getRight());
            }
            Pair<Matrix, Matrix> mix = composeMixture(alphas, Ts, node.probs);
            alpha = mix.getLeft();
            T = mix.getRight();
        } else if ("loop".equals(node.type)) {
            Pair<Matrix, Matrix> body = composeNode(node.kids[0]);
            Pair<Matrix, Matrix> rep = composeLoopGeometric(body.getLeft(), body.getRight(), node.count);
            alpha = rep.getLeft();
            T = rep.getRight();
        } else {
            throw new IllegalStateException("Unknown series-parallel node type '" + node.type + "'.");
        }

        node.alpha = alpha;
        node.T = T;
        node.valid = true;
        return new Pair<Matrix, Matrix>(alpha, T);
    }

    private static void invalidateBranch(SPTree tree, int k) {
        while (k >= 0) {
            tree.nodes.get(k).valid = false;
            k = tree.nodes.get(k).parent;
        }
    }

    private static double[] spExecutionCounts(SPTree tree) {
        double[] execs = new double[tree.nodes.size()];
        execs[tree.root] = 1.0;
        LinkedList<Integer> stack = new LinkedList<Integer>();
        stack.add(tree.root);
        while (!stack.isEmpty()) {
            int k = stack.removeLast();
            SPNode node = tree.nodes.get(k);
            for (int i = 0; i < node.kids.length; i++) {
                if ("or".equals(node.type)) {
                    execs[node.kids[i]] = execs[k] * node.probs[i];
                } else if ("loop".equals(node.type)) {
                    execs[node.kids[i]] = execs[k] * node.count;
                } else {
                    execs[node.kids[i]] = execs[k];
                }
                stack.add(node.kids[i]);
            }
        }
        return execs;
    }

    private static class WorkflowStructure {
        List<List<Integer>> adjList;
        int[] inDegree;
        int[] outDegree;
        List<ForkInfo> forkInfo;
        List<JoinInfo> joinInfo;
        List<LoopInfo> loopInfo;
    }

    private static class ForkInfo {
        String type;
        int preAct;
        int[] postActs;
        double[] probs;
    }

    private static class JoinInfo {
        String type;
        int[] preActs;
        int postAct;
    }

    private static class LoopInfo {
        int preAct;
        int[] loopActs;
        int endAct;
        double count;
    }

    private WorkflowStructure analyzeStructure() {
        int n = activities.size();
        WorkflowStructure structure = new WorkflowStructure();

        structure.adjList = new ArrayList<List<Integer>>();
        for (int i = 0; i < n; i++) {
            structure.adjList.add(new ArrayList<Integer>());
        }
        structure.inDegree = new int[n];
        structure.outDegree = new int[n];
        structure.forkInfo = new ArrayList<ForkInfo>();
        structure.joinInfo = new ArrayList<JoinInfo>();
        structure.loopInfo = new ArrayList<LoopInfo>();

        for (ActivityPrecedence prec : precedences) {
            List<String> preActNames = prec.getPreActs();
            List<String> postActNames = prec.getPostActs();

            int[] preInds = new int[preActNames.size()];
            for (int i = 0; i < preActNames.size(); i++) {
                preInds[i] = activityMap.get(preActNames.get(i));
            }

            int[] postInds = new int[postActNames.size()];
            for (int i = 0; i < postActNames.size(); i++) {
                postInds[i] = activityMap.get(postActNames.get(i));
            }

            for (int preIdx : preInds) {
                for (int postIdx : postInds) {
                    structure.adjList.get(preIdx).add(postIdx);
                    structure.outDegree[preIdx]++;
                    structure.inDegree[postIdx]++;
                }
            }

            String postType = prec.getPostType();
            String preType = prec.getPreType();
            Matrix postParams = prec.getPostParams();

            if (ActivityPrecedenceType.POST_AND.equals(postType)) {
                ForkInfo fork = new ForkInfo();
                fork.type = "and";
                fork.preAct = preInds[0];
                fork.postActs = postInds;
                structure.forkInfo.add(fork);
            } else if (ActivityPrecedenceType.POST_OR.equals(postType)) {
                ForkInfo fork = new ForkInfo();
                fork.type = "or";
                fork.preAct = preInds[0];
                fork.postActs = postInds;
                fork.probs = new double[postParams.getNumElements()];
                for (int i = 0; i < postParams.getNumElements(); i++) {
                    fork.probs[i] = postParams.get(i);
                }
                structure.forkInfo.add(fork);
            } else if (ActivityPrecedenceType.POST_LOOP.equals(postType)) {
                LoopInfo loop = new LoopInfo();
                loop.preAct = preInds[0];
                if (postInds.length >= 2) {
                    loop.loopActs = Arrays.copyOf(postInds, postInds.length - 1);
                    loop.endAct = postInds[postInds.length - 1];
                } else {
                    loop.loopActs = postInds;
                    loop.endAct = -1;
                }
                loop.count = postParams.get(0);
                structure.loopInfo.add(loop);
            }

            if (ActivityPrecedenceType.PRE_AND.equals(preType)) {
                JoinInfo join = new JoinInfo();
                join.type = "and";
                join.preActs = preInds;
                join.postAct = postInds[0];
                structure.joinInfo.add(join);
            } else if (ActivityPrecedenceType.PRE_OR.equals(preType)) {
                JoinInfo join = new JoinInfo();
                join.type = "or";
                join.preActs = preInds;
                join.postAct = postInds[0];
                structure.joinInfo.add(join);
            }
        }

        return structure;
    }

    private Pair<Matrix, Matrix> composeSerialWorkflow(List<List<Integer>> adjList) {
        int[] order = topologicalSort(adjList);

        Pair<Matrix, Matrix> result = activities.get(order[0]).getPHRepresentation();
        Matrix alpha = result.getLeft();
        Matrix T = result.getRight();

        for (int i = 1; i < order.length; i++) {
            Pair<Matrix, Matrix> next = activities.get(order[i]).getPHRepresentation();
            Pair<Matrix, Matrix> composed = composeSerial(alpha, T, next.getLeft(), next.getRight());
            alpha = composed.getLeft();
            T = composed.getRight();
        }

        return new Pair<Matrix, Matrix>(alpha, T);
    }

    private int[] topologicalSort(List<List<Integer>> adjList) {
        int n = activities.size();
        int[] inDeg = new int[n];

        for (int i = 0; i < n; i++) {
            for (int j : adjList.get(i)) {
                inDeg[j]++;
            }
        }

        Queue<Integer> queue = new LinkedList<Integer>();
        for (int i = 0; i < n; i++) {
            if (inDeg[i] == 0) {
                queue.add(i);
            }
        }

        List<Integer> order = new ArrayList<Integer>();
        while (!queue.isEmpty()) {
            int curr = queue.poll();
            order.add(curr);

            for (int next : adjList.get(curr)) {
                inDeg[next]--;
                if (inDeg[next] == 0) {
                    queue.add(next);
                }
            }
        }

        Set<Integer> orderSet = new HashSet<Integer>(order);
        for (int i = 0; i < n; i++) {
            if (!orderSet.contains(i)) {
                order.add(i);
            }
        }

        int[] result = new int[order.size()];
        for (int i = 0; i < order.size(); i++) {
            result[i] = order.get(i);
        }
        return result;
    }

    private Pair<Matrix, Matrix> composeComplexWorkflow(WorkflowStructure structure) {
        int n = activities.size();

        Matrix[] blockAlpha = new Matrix[n];
        Matrix[] blockT = new Matrix[n];
        // isConsumed tracks activities that have been absorbed into another block
        // (should be skipped in final composition)
        boolean[] isConsumed = new boolean[n];
        // isProcessed tracks activities whose blocks have been updated with composite results
        boolean[] isProcessed = new boolean[n];

        for (int i = 0; i < n; i++) {
            Pair<Matrix, Matrix> ph = activities.get(i).getPHRepresentation();
            blockAlpha[i] = ph.getLeft();
            blockT[i] = ph.getRight();
        }

        for (LoopInfo loop : structure.loopInfo) {
            int preIdx = loop.preAct;

            Matrix alphaLoop, TLoop;
            if (loop.loopActs.length == 1) {
                Pair<Matrix, Matrix> ph = activities.get(loop.loopActs[0]).getPHRepresentation();
                alphaLoop = ph.getLeft();
                TLoop = ph.getRight();
            } else {
                Pair<Matrix, Matrix> ph = activities.get(loop.loopActs[0]).getPHRepresentation();
                alphaLoop = ph.getLeft();
                TLoop = ph.getRight();
                for (int j = 1; j < loop.loopActs.length; j++) {
                    Pair<Matrix, Matrix> next = activities.get(loop.loopActs[j]).getPHRepresentation();
                    Pair<Matrix, Matrix> composed = composeSerial(alphaLoop, TLoop, next.getLeft(), next.getRight());
                    alphaLoop = composed.getLeft();
                    TLoop = composed.getRight();
                }
            }

            // Repeat the body a geometric number of times of mean COUNT, which
            // is the POST_LOOP semantics of the activity graph
            Pair<Matrix, Matrix> conv = composeLoopGeometric(alphaLoop, TLoop, loop.count);
            Pair<Matrix, Matrix> result = composeSerial(blockAlpha[preIdx], blockT[preIdx], conv.getLeft(), conv.getRight());

            if (loop.endAct >= 0) {
                Pair<Matrix, Matrix> endPh = activities.get(loop.endAct).getPHRepresentation();
                result = composeSerial(result.getLeft(), result.getRight(), endPh.getLeft(), endPh.getRight());
                isConsumed[loop.endAct] = true;
            }

            blockAlpha[preIdx] = result.getLeft();
            blockT[preIdx] = result.getRight();
            isProcessed[preIdx] = true;
            for (int idx : loop.loopActs) {
                isConsumed[idx] = true;
            }
        }

        for (ForkInfo fork : structure.forkInfo) {
            if ("and".equals(fork.type)) {
                JoinInfo matchingJoin = findMatchingJoin(fork.postActs, structure.joinInfo, "and");

                if (matchingJoin != null) {
                    int preIdx = fork.preAct;
                    int postIdx = matchingJoin.postAct;

                    Pair<Matrix, Matrix> parResult = composeAndForkBlock(fork.postActs, blockAlpha, blockT);

                    Pair<Matrix, Matrix> result;
                    if (!isProcessed[preIdx]) {
                        result = composeSerial(blockAlpha[preIdx], blockT[preIdx], parResult.getLeft(), parResult.getRight());
                    } else {
                        result = parResult;
                    }

                    if (!isConsumed[postIdx] && !isProcessed[postIdx]) {
                        result = composeSerial(result.getLeft(), result.getRight(), blockAlpha[postIdx], blockT[postIdx]);
                    }

                    blockAlpha[preIdx] = result.getLeft();
                    blockT[preIdx] = result.getRight();
                    isProcessed[preIdx] = true;
                    for (int idx : fork.postActs) {
                        isConsumed[idx] = true;
                    }
                    isConsumed[postIdx] = true;
                }
            }
        }

        for (ForkInfo fork : structure.forkInfo) {
            if ("or".equals(fork.type)) {
                JoinInfo matchingJoin = findMatchingJoin(fork.postActs, structure.joinInfo, "or");

                int preIdx = fork.preAct;

                Pair<Matrix, Matrix> orResult = composeOrForkBlock(fork.postActs, fork.probs, blockAlpha, blockT);

                Pair<Matrix, Matrix> result;
                if (!isProcessed[preIdx]) {
                    result = composeSerial(blockAlpha[preIdx], blockT[preIdx], orResult.getLeft(), orResult.getRight());
                } else {
                    result = orResult;
                }

                if (matchingJoin != null) {
                    int postIdx = matchingJoin.postAct;
                    if (!isConsumed[postIdx] && !isProcessed[postIdx]) {
                        result = composeSerial(result.getLeft(), result.getRight(), blockAlpha[postIdx], blockT[postIdx]);
                        isConsumed[postIdx] = true;
                    }
                }

                blockAlpha[preIdx] = result.getLeft();
                blockT[preIdx] = result.getRight();
                isProcessed[preIdx] = true;
                for (int idx : fork.postActs) {
                    isConsumed[idx] = true;
                }
            }
        }

        int[] order = topologicalSort(structure.adjList);
        Matrix alpha = null;
        Matrix T = null;

        // Compose all non-consumed activities in topological order
        // Block roots (isProcessed but not isConsumed) contain composite results
        for (int idx : order) {
            if (!isConsumed[idx]) {
                if (alpha == null) {
                    alpha = blockAlpha[idx];
                    T = blockT[idx];
                } else {
                    Pair<Matrix, Matrix> composed = composeSerial(alpha, T, blockAlpha[idx], blockT[idx]);
                    alpha = composed.getLeft();
                    T = composed.getRight();
                }
            }
        }

        if (alpha == null) {
            Pair<Matrix, Matrix> ph = activities.get(0).getPHRepresentation();
            alpha = ph.getLeft();
            T = ph.getRight();
        }

        return new Pair<Matrix, Matrix>(alpha, T);
    }

    private JoinInfo findMatchingJoin(int[] postActs, List<JoinInfo> joinInfo, String joinType) {
        Set<Integer> postSet = new HashSet<Integer>();
        for (int idx : postActs) {
            postSet.add(idx);
        }

        for (JoinInfo join : joinInfo) {
            if (joinType.equals(join.type)) {
                Set<Integer> preSet = new HashSet<Integer>();
                for (int idx : join.preActs) {
                    preSet.add(idx);
                }
                if (postSet.equals(preSet)) {
                    return join;
                }
            }
        }
        return null;
    }

    private Pair<Matrix, Matrix> composeAndForkBlock(int[] parallelInds, Matrix[] blockAlpha, Matrix[] blockT) {
        Matrix alpha = blockAlpha[parallelInds[0]];
        Matrix T = blockT[parallelInds[0]];

        for (int i = 1; i < parallelInds.length; i++) {
            Pair<Matrix, Matrix> composed = composeParallel(alpha, T, blockAlpha[parallelInds[i]], blockT[parallelInds[i]]);
            alpha = composed.getLeft();
            T = composed.getRight();
        }

        return new Pair<Matrix, Matrix>(alpha, T);
    }

    private Pair<Matrix, Matrix> composeOrForkBlock(int[] branchInds, double[] probs, Matrix[] blockAlpha, Matrix[] blockT) {
        int totalPhases = 0;
        for (int idx : branchInds) {
            totalPhases += blockT[idx].getNumRows();
        }

        Matrix T = new Matrix(totalPhases, totalPhases, totalPhases * totalPhases);
        Matrix alpha = new Matrix(1, totalPhases, totalPhases);

        int offset = 0;
        for (int i = 0; i < branchInds.length; i++) {
            int idx = branchInds[i];
            Matrix alphaI = blockAlpha[idx];
            Matrix TI = blockT[idx];
            int nI = TI.getNumRows();

            for (int r = 0; r < nI; r++) {
                for (int c = 0; c < nI; c++) {
                    double val = TI.get(r, c);
                    if (Math.abs(val) > GlobalConstants.Zero) {
                        T.set(offset + r, offset + c, val);
                    }
                }
            }

            for (int j = 0; j < alphaI.getNumElements(); j++) {
                alpha.set(0, offset + j, probs[i] * alphaI.get(j));
            }

            offset += nI;
        }

        return new Pair<Matrix, Matrix>(alpha, T);
    }

    public static Pair<Matrix, Matrix> composeSerial(Matrix alpha1, Matrix T1, Matrix alpha2, Matrix T2) {
        int n1 = T1.getNumRows();
        int n2 = T2.getNumRows();

        Matrix e1 = Matrix.ones(n1, 1);
        Matrix absRate1 = T1.mult(e1).scale(-1);

        Matrix TOut = new Matrix(n1 + n2, n1 + n2, (n1 + n2) * (n1 + n2));

        // The two generator blocks are copied by walking their NONZEROS. Probing all
        // n1^2 (resp. n2^2) cells cost a sparse get -- a binary search in a CSC column
        // -- per cell, and every cell it found absent held 0.0, which the
        // |val| > Zero guard rejected anyway. Same entries, same guard, so the result
        // is bit-identical; only the wasted probes are gone.
        Iterator<MatrixEntry> it1 = T1.nonZeroIterator();
        while (it1.hasNext()) {
            MatrixEntry e = it1.next();
            if (Math.abs(e.value) > GlobalConstants.Zero) {
                TOut.set(e.row, e.col, e.value);
            }
        }

        // absRate1(r) and alpha2(c) are read n1*n2 times between them; hoist each out
        // of the loop that does not index it.
        double[] a2 = new double[n2];
        for (int c = 0; c < n2; c++) a2[c] = alpha2.get(c);
        for (int r = 0; r < n1; r++) {
            double ar = absRate1.get(r, 0);
            for (int c = 0; c < n2; c++) {
                double val = ar * a2[c];
                if (Math.abs(val) > GlobalConstants.Zero) {
                    TOut.set(r, n1 + c, val);
                }
            }
        }

        Iterator<MatrixEntry> it2 = T2.nonZeroIterator();
        while (it2.hasNext()) {
            MatrixEntry e = it2.next();
            if (Math.abs(e.value) > GlobalConstants.Zero) {
                TOut.set(n1 + e.row, n1 + e.col, e.value);
            }
        }

        Matrix alphaOut = new Matrix(1, n1 + n2, n1 + n2);
        double mass1 = 0;
        for (int i = 0; i < n1; i++) {
            alphaOut.set(0, i, alpha1.get(i));
            mass1 += alpha1.get(i);
        }
        // A defective alpha1 carries an atom at zero, which starts the second
        // law immediately; this is aph_simplify pattern 1
        double defect1 = 1.0 - mass1;
        if (Math.abs(defect1) > GlobalConstants.Zero) {
            for (int c = 0; c < n2; c++) {
                alphaOut.set(0, n1 + c, defect1 * alpha2.get(c));
            }
        }

        return new Pair<Matrix, Matrix>(alphaOut, TOut);
    }

    /**
     * Probabilistic mixture of several PH laws, which is aph_simplify pattern
     * 3 generalised to any number of branches.
     *
     * @param alphas branch initial vectors
     * @param Ts     branch subgenerators
     * @param probs  branch probabilities
     * @return the mixture law
     */
    public static Pair<Matrix, Matrix> composeMixture(List<Matrix> alphas, List<Matrix> Ts, double[] probs) {
        int total = 0;
        for (Matrix T : Ts) {
            total += T.getNumRows();
        }

        Matrix TOut = new Matrix(total, total, total * total);
        Matrix alphaOut = new Matrix(1, total, total);

        int offset = 0;
        for (int i = 0; i < Ts.size(); i++) {
            Matrix TI = Ts.get(i);
            Matrix alphaI = alphas.get(i);
            int nI = TI.getNumRows();
            for (int r = 0; r < nI; r++) {
                for (int c = 0; c < nI; c++) {
                    double val = TI.get(r, c);
                    if (Math.abs(val) > GlobalConstants.Zero) {
                        TOut.set(offset + r, offset + c, val);
                    }
                }
            }
            for (int j = 0; j < nI; j++) {
                alphaOut.set(0, offset + j, probs[i] * alphaI.get(j));
            }
            offset += nI;
        }

        return new Pair<Matrix, Matrix>(alphaOut, TOut);
    }

    /**
     * Geometric repetition of a PH law, the POST_LOOP semantics of an activity
     * graph: the number of executions of the body is geometric of mean COUNT.
     * <p>
     * For COUNT&gt;=1 the body runs at least once and repeats on absorption
     * with probability P = 1-1/COUNT, so T_OUT = T + P/D*(-T*e)*ALPHA and
     * ALPHA_OUT = ALPHA/D with D = 1-P*(1-ALPHA*e) the correction for an atom
     * at zero in ALPHA. The order is that of the body, unlike the COUNT-fold
     * convolution of composeRepeat, and the mean is COUNT times the mean of
     * the body in both cases. For COUNT&lt;1 the body is executed at most once,
     * with probability COUNT.
     * </p>
     *
     * @param alpha body initial vector
     * @param T     body subgenerator
     * @param count mean number of executions
     * @return the repeated law
     */
    public static Pair<Matrix, Matrix> composeLoopGeometric(Matrix alpha, Matrix T, double count) {
        int n = T.getNumRows();

        if (count <= 0) {
            return new Pair<Matrix, Matrix>(Matrix.singleton(1.0),
                    Matrix.singleton(-GlobalConstants.Immediate));
        }

        if (Math.abs(count - 1.0) <= GlobalConstants.FineTol) {
            return new Pair<Matrix, Matrix>(alpha, T);
        }

        if (count < 1) {
            // Executed with probability COUNT, skipped otherwise
            Matrix TOut = new Matrix(n + 1, n + 1, (n + 1) * (n + 1));
            Matrix alphaOut = new Matrix(1, n + 1, n + 1);
            for (int r = 0; r < n; r++) {
                for (int c = 0; c < n; c++) {
                    double val = T.get(r, c);
                    if (Math.abs(val) > GlobalConstants.Zero) {
                        TOut.set(r, c, val);
                    }
                }
                alphaOut.set(0, r, count * alpha.get(r));
            }
            TOut.set(n, n, -GlobalConstants.Immediate);
            alphaOut.set(0, n, 1.0 - count);
            return new Pair<Matrix, Matrix>(alphaOut, TOut);
        }

        double p = 1.0 - 1.0 / count;
        double mass = 0;
        for (int i = 0; i < n; i++) {
            mass += alpha.get(i);
        }
        double denom = 1.0 - p * (1.0 - mass);

        Matrix e = Matrix.ones(n, 1);
        Matrix absRate = T.mult(e).scale(-1);

        Matrix TOut = new Matrix(n, n, n * n);
        for (int r = 0; r < n; r++) {
            for (int c = 0; c < n; c++) {
                double val = T.get(r, c) + (p / denom) * absRate.get(r, 0) * alpha.get(c);
                if (Math.abs(val) > GlobalConstants.Zero) {
                    TOut.set(r, c, val);
                }
            }
        }

        Matrix alphaOut = new Matrix(1, n, n);
        for (int i = 0; i < n; i++) {
            alphaOut.set(0, i, alpha.get(i) / denom);
        }

        return new Pair<Matrix, Matrix>(alphaOut, TOut);
    }

    /**
     * True when the phase graph of T has no cycle. A geometric loop over a
     * body of two or more phases closes a cycle, so the composed law is a PH
     * and not an APH.
     *
     * @param T subgenerator
     * @return true when acyclic
     */
    public static boolean isAcyclicGenerator(Matrix T) {
        int n = T.getNumRows();
        int[] inDeg = new int[n];
        for (int r = 0; r < n; r++) {
            for (int c = 0; c < n; c++) {
                if (r != c && Math.abs(T.get(r, c)) > GlobalConstants.ArcTol) {
                    inDeg[c]++;
                }
            }
        }
        Queue<Integer> queue = new LinkedList<Integer>();
        for (int i = 0; i < n; i++) {
            if (inDeg[i] == 0) {
                queue.add(i);
            }
        }
        int visited = 0;
        while (!queue.isEmpty()) {
            int curr = queue.poll();
            visited++;
            for (int c = 0; c < n; c++) {
                if (c != curr && Math.abs(T.get(curr, c)) > GlobalConstants.ArcTol) {
                    inDeg[c]--;
                    if (inDeg[c] == 0) {
                        queue.add(c);
                    }
                }
            }
        }
        return visited == n;
    }

    public static Pair<Matrix, Matrix> composeParallel(Matrix alpha1, Matrix T1, Matrix alpha2, Matrix T2) {
        int n1 = T1.getNumRows();
        int n2 = T2.getNumRows();

        Matrix e1 = Matrix.ones(n1, 1);
        Matrix e2 = Matrix.ones(n2, 1);
        Matrix absRate1 = T1.mult(e1).scale(-1);
        Matrix absRate2 = T2.mult(e2).scale(-1);

        int nBoth = n1 * n2;
        int nOnly1 = n1;
        int nOnly2 = n2;
        int nTotal = nBoth + nOnly1 + nOnly2;

        Matrix TOut = new Matrix(nTotal, nTotal, nTotal * nTotal);

        Matrix TBoth = T1.krons(T2);
        for (int r = 0; r < nBoth; r++) {
            for (int c = 0; c < nBoth; c++) {
                double val = TBoth.get(r, c);
                if (Math.abs(val) > GlobalConstants.Zero) {
                    TOut.set(r, c, val);
                }
            }
        }

        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n2; j++) {
                int bothIdx = i * n2 + j;
                int only1Idx = nBoth + i;
                double val = TOut.get(bothIdx, only1Idx) + absRate2.get(j, 0);
                TOut.set(bothIdx, only1Idx, val);
            }
        }

        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n2; j++) {
                int bothIdx = i * n2 + j;
                int only2Idx = nBoth + nOnly1 + j;
                double val = TOut.get(bothIdx, only2Idx) + absRate1.get(i, 0);
                TOut.set(bothIdx, only2Idx, val);
            }
        }

        for (int r = 0; r < n1; r++) {
            for (int c = 0; c < n1; c++) {
                double val = T1.get(r, c);
                if (Math.abs(val) > GlobalConstants.Zero) {
                    TOut.set(nBoth + r, nBoth + c, val);
                }
            }
        }

        for (int r = 0; r < n2; r++) {
            for (int c = 0; c < n2; c++) {
                double val = T2.get(r, c);
                if (Math.abs(val) > GlobalConstants.Zero) {
                    TOut.set(nBoth + nOnly1 + r, nBoth + nOnly1 + c, val);
                }
            }
        }

        Matrix alphaOut = new Matrix(1, nTotal, nTotal);
        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n2; j++) {
                int bothIdx = i * n2 + j;
                alphaOut.set(0, bothIdx, alpha1.get(i) * alpha2.get(j));
            }
        }

        return new Pair<Matrix, Matrix>(alphaOut, TOut);
    }

    /**
     * Deterministic COUNT-fold convolution of a PH law. The POST_LOOP
     * precedence of an activity graph is instead geometric; use
     * composeLoopGeometric for it.
     *
     * @param alpha initial vector
     * @param T     subgenerator
     * @param count number of executions
     * @return the convolved law
     */
    public static Pair<Matrix, Matrix> composeRepeat(Matrix alpha, Matrix T, int count) {
        if (count <= 0) {
            return new Pair<Matrix, Matrix>(Matrix.singleton(1.0),
                    Matrix.singleton(-GlobalConstants.Immediate));
        }

        if (count == 1) {
            return new Pair<Matrix, Matrix>(alpha, T);
        }

        Matrix alphaOut = alpha;
        Matrix TOut = T;
        for (int i = 1; i < count; i++) {
            Pair<Matrix, Matrix> composed = composeSerial(alphaOut, TOut, alpha, T);
            alphaOut = composed.getLeft();
            TOut = composed.getRight();
        }

        return new Pair<Matrix, Matrix>(alphaOut, TOut);
    }

    public static ActivityPrecedence[] Serial(WorkflowActivity... activities) {
        List<String> names = new ArrayList<String>();
        for (WorkflowActivity act : activities) {
            names.add(act.getName());
        }
        return ActivityPrecedence.Serial(names);
    }

    public static ActivityPrecedence AndFork(WorkflowActivity preAct, List<WorkflowActivity> postActs) {
        List<String> postNames = new ArrayList<String>();
        for (WorkflowActivity act : postActs) {
            postNames.add(act.getName());
        }
        return ActivityPrecedence.AndFork(preAct.getName(), postNames);
    }

    public static ActivityPrecedence AndJoin(List<WorkflowActivity> preActs, WorkflowActivity postAct) {
        List<String> preNames = new ArrayList<String>();
        for (WorkflowActivity act : preActs) {
            preNames.add(act.getName());
        }
        return ActivityPrecedence.AndJoin(preNames, postAct.getName());
    }

    /**
     * AND-join waiting for QUORUM of the branches. Workflow refuses a partial
     * join, since it is not the maximum of the branches; the overload exists so
     * that such a graph is rejected by name rather than mis-composed.
     *
     * @param preActs branch tail activities
     * @param postAct activity after the join
     * @param quorum  number of branches to wait for
     * @return the precedence
     */
    public static ActivityPrecedence AndJoin(List<WorkflowActivity> preActs, WorkflowActivity postAct, int quorum) {
        List<String> preNames = new ArrayList<String>();
        for (WorkflowActivity act : preActs) {
            preNames.add(act.getName());
        }
        return ActivityPrecedence.AndJoin(preNames, postAct.getName(), Matrix.singleton(quorum));
    }

    public static ActivityPrecedence OrFork(WorkflowActivity preAct, List<WorkflowActivity> postActs, double[] probs) {
        List<String> postNames = new ArrayList<String>();
        for (WorkflowActivity act : postActs) {
            postNames.add(act.getName());
        }
        return ActivityPrecedence.OrFork(preAct.getName(), postNames, new Matrix(probs));
    }

    public static ActivityPrecedence OrJoin(List<WorkflowActivity> preActs, WorkflowActivity postAct) {
        List<String> preNames = new ArrayList<String>();
        for (WorkflowActivity act : preActs) {
            preNames.add(act.getName());
        }
        return ActivityPrecedence.OrJoin(preNames, postAct.getName());
    }

    public static ActivityPrecedence Loop(WorkflowActivity preAct, List<WorkflowActivity> postActs, double count) {
        List<String> postNames = new ArrayList<String>();
        for (WorkflowActivity act : postActs) {
            postNames.add(act.getName());
        }
        return ActivityPrecedence.Loop(preAct.getName(), postNames, count);
    }

    /**
     * Load a workflow from a WfCommons JSON file.
     * <p>
     * WfCommons (<a href="https://github.com/wfcommons/workflow-schema">https://github.com/wfcommons/workflow-schema</a>)
     * is a standard format for representing scientific workflow traces.
     * </p>
     *
     * @param jsonFile Path to the WfCommons JSON file
     * @return Workflow object
     * @throws IOException If the file cannot be read
     */
    public static Workflow fromWfCommons(String jsonFile) throws IOException {
        return WfCommonsLoader.load(jsonFile);
    }

    /**
     * Load a workflow from a WfCommons JSON file with options.
     *
     * @param jsonFile Path to the WfCommons JSON file
     * @param options  Loader options
     * @return Workflow object
     * @throws IOException If the file cannot be read
     */
    public static Workflow fromWfCommons(String jsonFile, WfCommonsOptions options) throws IOException {
        return WfCommonsLoader.load(jsonFile, options);
    }
}
