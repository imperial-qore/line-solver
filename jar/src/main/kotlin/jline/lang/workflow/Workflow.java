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
import jline.util.Pair;
import jline.util.matrix.Matrix;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

/**
 * A computational workflow that can be converted to a phase-type distribution.
 */
public class Workflow extends Model {

    private List<WorkflowActivity> activities;
    private List<ActivityPrecedence> precedences;
    private Map<String, Integer> activityMap;
    private APH cachedPH;

    public Workflow(String name) {
        super(name);
        this.activities = new ArrayList<WorkflowActivity>();
        this.precedences = new ArrayList<ActivityPrecedence>();
        this.activityMap = new HashMap<String, Integer>();
        this.cachedPH = null;
    }

    public WorkflowActivity addActivity(String name, double meanServiceTime) {
        WorkflowActivity act = new WorkflowActivity(this, name, meanServiceTime);
        activities.add(act);
        act.setIndex(activities.size() - 1);
        activityMap.put(name, act.getIndex());
        cachedPH = null;
        return act;
    }

    public WorkflowActivity addActivity(String name, Distribution hostDemand) {
        WorkflowActivity act = new WorkflowActivity(this, name, hostDemand);
        activities.add(act);
        act.setIndex(activities.size() - 1);
        activityMap.put(name, act.getIndex());
        cachedPH = null;
        return act;
    }

    public void addPrecedence(ActivityPrecedence prec) {
        precedences.add(prec);
        cachedPH = null;
    }

    public void addPrecedence(ActivityPrecedence[] precs) {
        for (ActivityPrecedence prec : precs) {
            precedences.add(prec);
        }
        cachedPH = null;
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
        }

        return new Pair<Boolean, String>(true, "");
    }

    public APH toPH() {
        if (cachedPH != null) {
            return cachedPH;
        }

        Pair<Boolean, String> validation = validate();
        if (!validation.getLeft()) {
            throw new IllegalStateException(validation.getRight());
        }

        Pair<Matrix, Matrix> ctmc = buildCTMC();
        cachedPH = new APH(ctmc.getLeft(), ctmc.getRight());
        return cachedPH;
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

            Pair<Matrix, Matrix> conv = composeRepeat(alphaLoop, TLoop, (int) loop.count);
            Pair<Matrix, Matrix> result = composeSerial(blockAlpha[preIdx], blockT[preIdx], conv.getLeft(), conv.getRight());

            if (loop.endAct >= 0) {
                Pair<Matrix, Matrix> endPh = activities.get(loop.endAct).getPHRepresentation();
                result = composeSerial(result.getLeft(), result.getRight(), endPh.getLeft(), endPh.getRight());
                isProcessed[loop.endAct] = true;
            }

            blockAlpha[preIdx] = result.getLeft();
            blockT[preIdx] = result.getRight();
            isProcessed[preIdx] = true;
            for (int idx : loop.loopActs) {
                isProcessed[idx] = true;
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

                    if (!isProcessed[postIdx]) {
                        result = composeSerial(result.getLeft(), result.getRight(), blockAlpha[postIdx], blockT[postIdx]);
                    }

                    blockAlpha[preIdx] = result.getLeft();
                    blockT[preIdx] = result.getRight();
                    isProcessed[preIdx] = true;
                    for (int idx : fork.postActs) {
                        isProcessed[idx] = true;
                    }
                    isProcessed[postIdx] = true;
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
                    if (!isProcessed[postIdx]) {
                        result = composeSerial(result.getLeft(), result.getRight(), blockAlpha[postIdx], blockT[postIdx]);
                        isProcessed[postIdx] = true;
                    }
                }

                blockAlpha[preIdx] = result.getLeft();
                blockT[preIdx] = result.getRight();
                isProcessed[preIdx] = true;
                for (int idx : fork.postActs) {
                    isProcessed[idx] = true;
                }
            }
        }

        int[] order = topologicalSort(structure.adjList);
        Matrix alpha = null;
        Matrix T = null;

        for (int idx : order) {
            if (!isProcessed[idx] || alpha == null) {
                if (alpha == null) {
                    alpha = blockAlpha[idx];
                    T = blockT[idx];
                } else if (!isProcessed[idx]) {
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

        for (int r = 0; r < n1; r++) {
            for (int c = 0; c < n1; c++) {
                double val = T1.get(r, c);
                if (Math.abs(val) > GlobalConstants.Zero) {
                    TOut.set(r, c, val);
                }
            }
        }

        for (int r = 0; r < n1; r++) {
            for (int c = 0; c < n2; c++) {
                double val = absRate1.get(r, 0) * alpha2.get(c);
                if (Math.abs(val) > GlobalConstants.Zero) {
                    TOut.set(r, n1 + c, val);
                }
            }
        }

        for (int r = 0; r < n2; r++) {
            for (int c = 0; c < n2; c++) {
                double val = T2.get(r, c);
                if (Math.abs(val) > GlobalConstants.Zero) {
                    TOut.set(n1 + r, n1 + c, val);
                }
            }
        }

        Matrix alphaOut = new Matrix(1, n1 + n2, n1 + n2);
        for (int i = 0; i < n1; i++) {
            alphaOut.set(0, i, alpha1.get(i));
        }

        return new Pair<Matrix, Matrix>(alphaOut, TOut);
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

    public static Pair<Matrix, Matrix> composeRepeat(Matrix alpha, Matrix T, int count) {
        if (count <= 0) {
            return new Pair<Matrix, Matrix>(Matrix.singleton(1.0), Matrix.singleton(-1e10));
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
