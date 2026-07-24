/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io.tikz;

import jline.lang.layered.*;
import jline.lang.constant.ActivityPrecedenceType;
import jline.util.matrix.Matrix;

import java.util.*;

/**
 * Traverses a LayeredNetwork model to extract interactions and fragments
 * for UML sequence diagram generation.
 */
public class SequenceDiagramTraverser {

    /**
     * Represents a single interaction (message or activation) in the sequence diagram.
     */
    public static class Interaction {
        public enum Type {
            SYNCH_CALL,    // Synchronous call (blocking)
            ASYNCH_CALL,   // Asynchronous call (fire-and-forget)
            REPLY          // Reply message
        }

        public Type type;
        public Task sourceTask;
        public Task targetTask;
        public Entry targetEntry;
        public Activity sourceActivity;
        public double callMean;
        public int level;

        public Interaction(Type type, Task sourceTask, Task targetTask, Entry targetEntry,
                           Activity sourceActivity, double callMean, int level) {
            this.type = type;
            this.sourceTask = sourceTask;
            this.targetTask = targetTask;
            this.targetEntry = targetEntry;
            this.sourceActivity = sourceActivity;
            this.callMean = callMean;
            this.level = level;
        }
    }

    /**
     * Represents a combined fragment (loop, par, alt) in the sequence diagram.
     */
    public static class Fragment {
        public enum FragmentType {
            LOOP,   // Loop construct
            PAR,    // Parallel execution (fork-join)
            ALT,    // Alternative paths
            OPT     // Optional execution
        }

        public FragmentType type;
        public int startLevel;
        public int endLevel;
        public List<Task> involvedTasks;
        public String label;
        public List<Fragment> operands;  // For ALT fragments, each alternative
        public Matrix params;  // Loop count or probabilities

        public Fragment(FragmentType type, int startLevel) {
            this.type = type;
            this.startLevel = startLevel;
            this.endLevel = startLevel;
            this.involvedTasks = new ArrayList<>();
            this.operands = new ArrayList<>();
        }
    }

    private final LayeredNetwork model;
    private final SequenceDiagramLayoutEngine layoutEngine;
    private final List<Interaction> interactions;
    private final List<Fragment> fragments;

    // Entry name -> Task mapping
    private Map<String, Task> entryToTask;

    // Entry name -> Entry mapping
    private Map<String, Entry> entryByName;

    public SequenceDiagramTraverser(LayeredNetwork model, SequenceDiagramLayoutEngine layoutEngine) {
        this.model = model;
        this.layoutEngine = layoutEngine;
        this.interactions = new ArrayList<>();
        this.fragments = new ArrayList<>();
        this.entryToTask = new HashMap<>();
        this.entryByName = new HashMap<>();

        buildMaps();
    }

    /**
     * Builds internal maps for efficient lookups.
     */
    private void buildMaps() {
        for (Map.Entry<Integer, Task> taskEntry : model.getTasks().entrySet()) {
            Task task = taskEntry.getValue();
            for (Entry entry : task.getEntries()) {
                entryToTask.put(entry.getName(), task);
                entryByName.put(entry.getName(), entry);
            }
        }
    }

    /**
     * Extracts all interactions from the model.
     */
    public List<Interaction> extractInteractions() {
        interactions.clear();

        for (Activity act : layoutEngine.getSortedActivities()) {
            int level = layoutEngine.getActivityLevel(act);
            Task sourceTask = act.getParent();
            if (sourceTask == null) continue;

            // Process synchronous calls
            Map<Integer, String> syncDests = act.getSyncCallDests();
            Matrix syncMeans = act.getSyncCallMeans();
            for (Map.Entry<Integer, String> dest : syncDests.entrySet()) {
                int idx = dest.getKey();
                String destEntryName = dest.getValue();
                Entry destEntry = entryByName.get(destEntryName);
                Task destTask = entryToTask.get(destEntryName);

                double callMean = 1.0;
                if (syncMeans != null && idx < syncMeans.getNumCols()) {
                    callMean = syncMeans.get(0, idx);
                }

                if (destTask != null && destEntry != null) {
                    interactions.add(new Interaction(
                            Interaction.Type.SYNCH_CALL,
                            sourceTask, destTask, destEntry, act, callMean, level
                    ));
                }
            }

            // Process asynchronous calls
            Map<Integer, String> asyncDests = act.getAsyncCallDests();
            Matrix asyncMeans = act.getAsyncCallMeans();
            for (Map.Entry<Integer, String> dest : asyncDests.entrySet()) {
                int idx = dest.getKey();
                String destEntryName = dest.getValue();
                Entry destEntry = entryByName.get(destEntryName);
                Task destTask = entryToTask.get(destEntryName);

                double callMean = 1.0;
                if (asyncMeans != null && idx < asyncMeans.getNumCols()) {
                    callMean = asyncMeans.get(0, idx);
                }

                if (destTask != null && destEntry != null) {
                    interactions.add(new Interaction(
                            Interaction.Type.ASYNCH_CALL,
                            sourceTask, destTask, destEntry, act, callMean, level
                    ));
                }
            }
        }

        return interactions;
    }

    /**
     * Extracts combined fragments from activity precedences.
     */
    public List<Fragment> extractFragments() {
        fragments.clear();

        // Track fork-join pairs for PAR fragments
        Map<String, ActivityPrecedence> forksByPostAct = new HashMap<>();
        Map<String, ActivityPrecedence> joinsByPreAct = new HashMap<>();

        for (Map.Entry<Integer, Task> taskEntry : model.getTasks().entrySet()) {
            Task task = taskEntry.getValue();

            for (ActivityPrecedence prec : task.getPrecedences()) {
                String postType = prec.getPostType();
                String preType = prec.getPreType();

                // Handle loops
                if (ActivityPrecedenceType.POST_LOOP.equals(postType)) {
                    Fragment loopFrag = createLoopFragment(task, prec);
                    if (loopFrag != null) {
                        fragments.add(loopFrag);
                    }
                }

                // Handle AND-fork (parallel start)
                if (ActivityPrecedenceType.POST_AND.equals(postType)) {
                    for (String postAct : prec.getPostActs()) {
                        forksByPostAct.put(postAct, prec);
                    }
                }

                // Handle AND-join (parallel end)
                if (ActivityPrecedenceType.PRE_AND.equals(preType)) {
                    for (String preAct : prec.getPreActs()) {
                        joinsByPreAct.put(preAct, prec);
                    }
                }

                // Handle OR-fork (alternatives)
                if (ActivityPrecedenceType.POST_OR.equals(postType)) {
                    Fragment altFrag = createAltFragment(task, prec);
                    if (altFrag != null) {
                        fragments.add(altFrag);
                    }
                }
            }
        }

        // Match fork-join pairs to create PAR fragments
        Set<ActivityPrecedence> processedForks = new HashSet<>();
        for (Map.Entry<Integer, Task> taskEntry : model.getTasks().entrySet()) {
            Task task = taskEntry.getValue();

            for (ActivityPrecedence prec : task.getPrecedences()) {
                if (ActivityPrecedenceType.POST_AND.equals(prec.getPostType())
                        && !processedForks.contains(prec)) {
                    Fragment parFrag = createParFragment(task, prec, joinsByPreAct);
                    if (parFrag != null) {
                        fragments.add(parFrag);
                        processedForks.add(prec);
                    }
                }
            }
        }

        return fragments;
    }

    /**
     * Creates a LOOP fragment from a POST_LOOP precedence.
     */
    private Fragment createLoopFragment(Task task, ActivityPrecedence prec) {
        List<String> postActs = prec.getPostActs();
        if (postActs.isEmpty()) return null;

        String loopActName = postActs.get(0);
        Activity loopAct = findActivityByName(task, loopActName);
        if (loopAct == null) return null;

        int startLevel = layoutEngine.getActivityLevel(loopAct);
        Fragment frag = new Fragment(Fragment.FragmentType.LOOP, startLevel);

        // Get loop count from postParams
        Matrix postParams = prec.getPostParams();
        double loopCount = 1.0;
        if (postParams != null && !postParams.isEmpty()) {
            loopCount = postParams.get(0, 0);
        }
        frag.label = String.format("loop [%.0f]", loopCount);
        frag.params = postParams;

        // Find end level based on last activity in loop
        int endLevel = startLevel;
        for (int i = 0; i < postActs.size() - 1; i++) {  // Last act is end, not in loop
            Activity act = findActivityByName(task, postActs.get(i));
            if (act != null) {
                int level = layoutEngine.getActivityLevel(act);
                if (level > endLevel) endLevel = level;
                if (!frag.involvedTasks.contains(task)) {
                    frag.involvedTasks.add(task);
                }
            }
        }
        frag.endLevel = endLevel;

        return frag;
    }

    /**
     * Creates a PAR fragment from a POST_AND fork precedence.
     */
    private Fragment createParFragment(Task task, ActivityPrecedence forkPrec,
                                        Map<String, ActivityPrecedence> joinsByPreAct) {
        List<String> parallelActs = forkPrec.getPostActs();
        if (parallelActs.isEmpty()) return null;

        // Find start level from fork activity
        String forkActName = forkPrec.getPreActs().get(0);
        Activity forkAct = findActivityByName(task, forkActName);
        int startLevel = forkAct != null ? layoutEngine.getActivityLevel(forkAct) : 0;

        Fragment frag = new Fragment(Fragment.FragmentType.PAR, startLevel);
        frag.label = "par";

        int endLevel = startLevel;
        for (String actName : parallelActs) {
            Activity act = findActivityByName(task, actName);
            if (act != null) {
                int level = layoutEngine.getActivityLevel(act);
                if (level > endLevel) endLevel = level;
                if (!frag.involvedTasks.contains(task)) {
                    frag.involvedTasks.add(task);
                }
            }
        }

        // Look for matching join
        for (String actName : parallelActs) {
            ActivityPrecedence joinPrec = joinsByPreAct.get(actName);
            if (joinPrec != null && !joinPrec.getPostActs().isEmpty()) {
                String joinActName = joinPrec.getPostActs().get(0);
                Activity joinAct = findActivityByName(task, joinActName);
                if (joinAct != null) {
                    int joinLevel = layoutEngine.getActivityLevel(joinAct);
                    if (joinLevel > endLevel) endLevel = joinLevel;
                }
            }
        }

        frag.endLevel = endLevel;
        return frag;
    }

    /**
     * Creates an ALT fragment from a POST_OR fork precedence.
     */
    private Fragment createAltFragment(Task task, ActivityPrecedence prec) {
        List<String> altActs = prec.getPostActs();
        if (altActs.isEmpty()) return null;

        String preActName = prec.getPreActs().get(0);
        Activity preAct = findActivityByName(task, preActName);
        int startLevel = preAct != null ? layoutEngine.getActivityLevel(preAct) : 0;

        Fragment frag = new Fragment(Fragment.FragmentType.ALT, startLevel);
        frag.label = "alt";
        frag.params = prec.getPostParams();

        int endLevel = startLevel;
        Matrix probs = prec.getPostParams();

        for (int i = 0; i < altActs.size(); i++) {
            String actName = altActs.get(i);
            Activity act = findActivityByName(task, actName);
            if (act != null) {
                int level = layoutEngine.getActivityLevel(act);
                if (level > endLevel) endLevel = level;

                // Create operand fragment for this alternative
                Fragment operand = new Fragment(Fragment.FragmentType.OPT, level);
                if (probs != null && i < probs.getNumCols()) {
                    operand.label = String.format("[%.2f]", probs.get(0, i));
                } else {
                    operand.label = "[else]";
                }
                operand.involvedTasks.add(task);
                frag.operands.add(operand);
            }
        }

        frag.endLevel = endLevel;
        if (!frag.involvedTasks.contains(task)) {
            frag.involvedTasks.add(task);
        }

        return frag;
    }

    /**
     * Finds an activity by name within a task.
     */
    private Activity findActivityByName(Task task, String name) {
        for (Activity act : task.getActivities()) {
            if (act.getName().equals(name)) {
                return act;
            }
        }
        // Also search model-wide
        for (Map.Entry<Integer, Activity> entry : model.getActivities().entrySet()) {
            if (entry.getValue().getName().equals(name)) {
                return entry.getValue();
            }
        }
        return null;
    }

    /**
     * Gets all interactions sorted by level.
     */
    public List<Interaction> getInteractions() {
        return interactions;
    }

    /**
     * Gets all fragments sorted by start level.
     */
    public List<Fragment> getFragments() {
        return fragments;
    }
}
