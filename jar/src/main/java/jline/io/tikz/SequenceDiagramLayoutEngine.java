/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io.tikz;

import jline.lang.layered.*;
import jline.lang.constant.ActivityPrecedenceType;

import java.util.*;

/**
 * Computes layout positions for UML sequence diagram elements.
 * Groups tasks by processor and determines vertical ordering of activities.
 */
public class SequenceDiagramLayoutEngine {

    private final LayeredNetwork model;
    private final SequenceDiagramOptions options;

    // Processor -> List of Tasks mapping for swim lane organization
    private Map<Host, List<Task>> processorTaskGroups;

    // Task -> X position mapping
    private Map<Task, Double> lifelinePositions;

    // Processor -> [minX, maxX] bounds
    private Map<Host, double[]> processorBounds;

    // Activity -> vertical level (time order)
    private Map<Activity, Integer> activityLevels;

    // Entry -> Task mapping for resolving call destinations
    private Map<String, Task> entryToTask;

    // Activity -> Entry mapping for determining bound entry
    private Map<Activity, Entry> activityToEntry;

    // Sorted list of activities in execution order
    private List<Activity> sortedActivities;

    public SequenceDiagramLayoutEngine(LayeredNetwork model, SequenceDiagramOptions options) {
        this.model = model;
        this.options = options;
        this.processorTaskGroups = new LinkedHashMap<>();
        this.lifelinePositions = new HashMap<>();
        this.processorBounds = new HashMap<>();
        this.activityLevels = new HashMap<>();
        this.entryToTask = new HashMap<>();
        this.activityToEntry = new HashMap<>();
        this.sortedActivities = new ArrayList<>();
    }

    /**
     * Computes the complete layout for the sequence diagram.
     */
    public void computeLayout() {
        buildStructureMaps();
        computeLifelinePositions();
        computeActivityLevels();
    }

    /**
     * Builds internal maps for efficient lookups.
     */
    private void buildStructureMaps() {
        // Group tasks by processor
        for (Map.Entry<Integer, Host> entry : model.getHosts().entrySet()) {
            Host host = entry.getValue();
            List<Task> tasks = new ArrayList<>(host.getTasks());
            processorTaskGroups.put(host, tasks);
        }

        // Build entry-to-task mapping
        for (Map.Entry<Integer, Task> entry : model.getTasks().entrySet()) {
            Task task = entry.getValue();
            for (Entry e : task.getEntries()) {
                entryToTask.put(e.getName(), task);
            }
        }

        // Build activity-to-entry mapping
        for (Map.Entry<Integer, Task> entry : model.getTasks().entrySet()) {
            Task task = entry.getValue();
            for (Entry e : task.getEntries()) {
                for (Map.Entry<Integer, String> bound : e.getBoundToActivity().entrySet()) {
                    String actName = bound.getValue();
                    Activity act = findActivityByName(task, actName);
                    if (act != null) {
                        activityToEntry.put(act, e);
                    }
                }
            }
        }
    }

    /**
     * Computes horizontal positions for each lifeline (task) with processor grouping.
     */
    private void computeLifelinePositions() {
        double currentX = 0;
        double lifelineSpacing = options.getLifelineSpacing();

        for (Map.Entry<Host, List<Task>> entry : processorTaskGroups.entrySet()) {
            Host processor = entry.getKey();
            List<Task> tasks = entry.getValue();
            if (tasks.isEmpty()) continue;

            double processorStartX = currentX;

            for (Task task : tasks) {
                lifelinePositions.put(task, currentX);
                currentX += lifelineSpacing;
            }

            double processorEndX = currentX - lifelineSpacing;
            processorBounds.put(processor, new double[]{processorStartX, processorEndX});

            // Add gap between processors
            currentX += lifelineSpacing * 0.5;
        }
    }

    /**
     * Computes vertical levels for activities based on execution order.
     * Uses topological sort considering precedences and call relationships.
     */
    private void computeActivityLevels() {
        // Build adjacency list for topological sort
        Map<Activity, List<Activity>> graph = new HashMap<>();
        Map<Activity, Integer> inDegree = new HashMap<>();

        // Initialize all activities
        for (Map.Entry<Integer, Activity> entry : model.getActivities().entrySet()) {
            Activity act = entry.getValue();
            graph.put(act, new ArrayList<>());
            inDegree.put(act, 0);
        }

        // Add edges from precedences
        for (Map.Entry<Integer, Task> taskEntry : model.getTasks().entrySet()) {
            Task task = taskEntry.getValue();
            for (ActivityPrecedence prec : task.getPrecedences()) {
                for (String preActName : prec.getPreActs()) {
                    Activity preAct = findActivityByName(task, preActName);
                    for (String postActName : prec.getPostActs()) {
                        Activity postAct = findActivityByName(task, postActName);
                        if (preAct != null && postAct != null) {
                            graph.get(preAct).add(postAct);
                            inDegree.put(postAct, inDegree.get(postAct) + 1);
                        }
                    }
                }
            }
        }

        // Add edges for synch calls (caller activity -> callee's bound activity)
        for (Map.Entry<Integer, Activity> entry : model.getActivities().entrySet()) {
            Activity act = entry.getValue();
            for (Map.Entry<Integer, String> dest : act.getSyncCallDests().entrySet()) {
                String destEntryName = dest.getValue();
                Entry destEntry = findEntryByName(destEntryName);
                if (destEntry != null) {
                    Activity destAct = findBoundActivity(destEntry);
                    if (destAct != null && graph.containsKey(destAct)) {
                        graph.get(act).add(destAct);
                        inDegree.put(destAct, inDegree.get(destAct) + 1);
                    }
                }
            }
        }

        // Kahn's algorithm for topological sort
        Queue<Activity> queue = new LinkedList<>();
        for (Map.Entry<Activity, Integer> entry : inDegree.entrySet()) {
            if (entry.getValue() == 0) {
                queue.add(entry.getKey());
            }
        }

        int level = 0;
        while (!queue.isEmpty()) {
            Activity current = queue.poll();
            activityLevels.put(current, level++);
            sortedActivities.add(current);

            for (Activity neighbor : graph.get(current)) {
                inDegree.put(neighbor, inDegree.get(neighbor) - 1);
                if (inDegree.get(neighbor) == 0) {
                    queue.add(neighbor);
                }
            }
        }

        // Handle activities not in the graph (disconnected)
        for (Map.Entry<Integer, Activity> entry : model.getActivities().entrySet()) {
            Activity act = entry.getValue();
            if (!activityLevels.containsKey(act)) {
                activityLevels.put(act, level++);
                sortedActivities.add(act);
            }
        }
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
        return null;
    }

    /**
     * Finds an entry by name in the model.
     */
    private Entry findEntryByName(String name) {
        for (Map.Entry<Integer, Entry> entry : model.getEntries().entrySet()) {
            if (entry.getValue().getName().equals(name)) {
                return entry.getValue();
            }
        }
        return null;
    }

    /**
     * Finds the activity bound to an entry.
     */
    private Activity findBoundActivity(Entry entry) {
        if (entry.getBoundToActivity().isEmpty()) {
            return null;
        }
        String actName = entry.getBoundToActivity().values().iterator().next();
        if (entry.getParent() != null) {
            return findActivityByName(entry.getParent(), actName);
        }
        return null;
    }

    // Getters for layout information

    public Map<Host, List<Task>> getProcessorTaskGroups() {
        return processorTaskGroups;
    }

    public double getLifelineX(Task task) {
        return lifelinePositions.getOrDefault(task, 0.0);
    }

    public int getActivityLevel(Activity activity) {
        return activityLevels.getOrDefault(activity, 0);
    }

    public double[] getProcessorBounds(Host processor) {
        return processorBounds.get(processor);
    }

    public List<Activity> getSortedActivities() {
        return sortedActivities;
    }

    public Task getTaskForEntry(String entryName) {
        return entryToTask.get(entryName);
    }

    public Entry getEntryForActivity(Activity activity) {
        return activityToEntry.get(activity);
    }

    public int getMaxLevel() {
        int max = 0;
        for (int level : activityLevels.values()) {
            if (level > max) max = level;
        }
        return max;
    }
}
