/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import static jline.io.InputOutputKt.line_warning;

import com.google.gson.Gson;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import jline.GlobalConstants;
import jline.lang.constant.ActivityPrecedenceType;
import jline.lang.layered.ActivityPrecedence;
import jline.lang.processes.APH;
import jline.lang.processes.Det;
import jline.lang.processes.Distribution;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.lang.processes.Immediate;
import jline.lang.workflow.Workflow;
import jline.lang.workflow.WorkflowActivity;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
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
 * Loader for WfCommons workflow JSON files.
 * <p>
 * WfCommons (<a href="https://github.com/wfcommons/workflow-schema">https://github.com/wfcommons/workflow-schema</a>)
 * is a standard format for representing scientific workflow traces.
 * This loader converts WfCommons JSON files into LINE Workflow objects.
 * </p>
 * <p>
 * Supported schema versions: 1.4, 1.5
 * </p>
 * <p>
 * Example usage:
 * <pre>
 * Workflow wf = WfCommonsLoader.load("workflow.json");
 * APH ph = wf.toPH();
 * System.out.println("Mean execution time: " + ph.getMean());
 * </pre>
 * </p>
 */
public class WfCommonsLoader {

    private static final List<String> SUPPORTED_SCHEMA_VERSIONS = Arrays.asList("1.3", "1.4", "1.5");

    private WfCommonsLoader() {
        // Static utility class
    }

    /**
     * Load a WfCommons JSON file into a Workflow object.
     *
     * @param jsonFile Path to the WfCommons JSON file
     * @return Workflow object
     * @throws IOException If the file cannot be read
     * @throws IllegalArgumentException If the JSON is invalid
     */
    public static Workflow load(String jsonFile) throws IOException {
        return load(jsonFile, new WfCommonsOptions());
    }

    /**
     * Load a WfCommons JSON file with options.
     *
     * @param jsonFile Path to the WfCommons JSON file
     * @param options  Loader options
     * @return Workflow object
     * @throws IOException If the file cannot be read
     * @throws IllegalArgumentException If the JSON is invalid
     */
    public static Workflow load(String jsonFile, WfCommonsOptions options) throws IOException {
        FileReader reader = new FileReader(jsonFile);
        JsonObject data = JsonParser.parseReader(reader).getAsJsonObject();
        reader.close();

        validateSchema(data);
        String workflowName = extractName(data, jsonFile);
        return buildWorkflow(data, workflowName, options);
    }

    /**
     * Load a WfCommons JSON file from a URL.
     * <p>
     * Useful for loading workflows directly from repositories like wfcommons/pegasus-instances.
     * </p>
     *
     * @param urlString URL pointing to a WfCommons JSON file
     * @return Workflow object
     * @throws IOException If the URL cannot be read
     * @throws IllegalArgumentException If the JSON is invalid
     */
    public static Workflow loadFromUrl(String urlString) throws IOException {
        return loadFromUrl(urlString, new WfCommonsOptions());
    }

    /**
     * Load a WfCommons JSON file from a URL with options.
     *
     * @param urlString URL pointing to a WfCommons JSON file
     * @param options   Loader options
     * @return Workflow object
     * @throws IOException If the URL cannot be read
     * @throws IllegalArgumentException If the JSON is invalid
     */
    public static Workflow loadFromUrl(String urlString, WfCommonsOptions options) throws IOException {
        URL url = new URL(urlString);
        BufferedReader reader = new BufferedReader(new InputStreamReader(url.openStream()));
        StringBuilder jsonBuilder = new StringBuilder();
        String line;
        while ((line = reader.readLine()) != null) {
            jsonBuilder.append(line);
        }
        reader.close();

        JsonObject data = JsonParser.parseString(jsonBuilder.toString()).getAsJsonObject();
        validateSchema(data);

        // Extract name from URL if not present in JSON
        String defaultName = urlString;
        int lastSlash = urlString.lastIndexOf('/');
        if (lastSlash >= 0) {
            defaultName = urlString.substring(lastSlash + 1);
        }

        String workflowName = extractName(data, defaultName);
        return buildWorkflow(data, workflowName, options);
    }

    /**
     * Load a workflow from a JSON string.
     *
     * @param json JSON string in WfCommons format
     * @return Workflow object
     * @throws IllegalArgumentException If the JSON is invalid
     */
    public static Workflow loadFromString(String json) {
        return loadFromString(json, new WfCommonsOptions());
    }

    /**
     * Load a workflow from a JSON string with options.
     *
     * @param json    JSON string in WfCommons format
     * @param options Loader options
     * @return Workflow object
     * @throws IllegalArgumentException If the JSON is invalid
     */
    public static Workflow loadFromString(String json, WfCommonsOptions options) {
        JsonObject data = JsonParser.parseString(json).getAsJsonObject();
        validateSchema(data);
        String workflowName = extractName(data, "string_input");
        return buildWorkflow(data, workflowName, options);
    }

    /**
     * Validate if a file is a valid WfCommons JSON.
     *
     * @param jsonFile Path to the file
     * @return true if valid
     */
    public static boolean validateFile(String jsonFile) {
        try {
            FileReader reader = new FileReader(jsonFile);
            JsonObject data = JsonParser.parseReader(reader).getAsJsonObject();
            reader.close();
            validateSchema(data);
            return true;
        } catch (Exception e) {
            return false;
        }
    }

    private static void validateSchema(JsonObject data) {
        if (!data.has("schemaVersion")) {
            throw new IllegalArgumentException("Missing schemaVersion field in WfCommons JSON.");
        }

        String version = data.get("schemaVersion").getAsString();
        if (!SUPPORTED_SCHEMA_VERSIONS.contains(version)) {
            line_warning("WfCommonsLoader", "Schema version %s may not be fully supported.", version);
        }

        if (!data.has("workflow")) {
            throw new IllegalArgumentException("Missing workflow field in WfCommons JSON.");
        }

        JsonObject workflow = data.getAsJsonObject("workflow");

        // Schema 1.4 and earlier use workflow.tasks directly
        // Schema 1.5+ uses workflow.specification.tasks
        JsonArray tasks = null;
        if (workflow.has("specification")) {
            JsonObject specification = workflow.getAsJsonObject("specification");
            if (specification.has("tasks")) {
                tasks = specification.getAsJsonArray("tasks");
            }
        } else if (workflow.has("tasks")) {
            // Schema 1.4 format: tasks directly under workflow
            tasks = workflow.getAsJsonArray("tasks");
        }

        if (tasks == null || tasks.size() == 0) {
            throw new IllegalArgumentException("Workflow must have at least one task.");
        }
    }

    private static String extractName(JsonObject data, String defaultName) {
        if (data.has("name") && !data.get("name").isJsonNull()) {
            String name = data.get("name").getAsString();
            if (!name.isEmpty()) {
                return sanitizeName(name);
            }
        }
        // Extract filename without extension
        int lastSlash = Math.max(defaultName.lastIndexOf('/'), defaultName.lastIndexOf('\\'));
        String filename = lastSlash >= 0 ? defaultName.substring(lastSlash + 1) : defaultName;
        int lastDot = filename.lastIndexOf('.');
        String baseName = lastDot >= 0 ? filename.substring(0, lastDot) : filename;
        return sanitizeName(baseName);
    }

    private static String sanitizeName(String name) {
        String sanitized = name.replaceAll("[^a-zA-Z0-9_]", "_");
        if (sanitized.isEmpty()) {
            sanitized = "Workflow";
        }
        return sanitized;
    }

    private static Workflow buildWorkflow(JsonObject data, String workflowName, WfCommonsOptions options) {
        Workflow wf = new Workflow(workflowName);

        JsonObject workflow = data.getAsJsonObject("workflow");

        // Detect schema format: 1.5+ uses specification.tasks, 1.4 uses tasks directly
        JsonArray tasks;
        boolean isLegacySchema = !workflow.has("specification");
        if (isLegacySchema) {
            // Schema 1.4 format: tasks directly under workflow with embedded runtime
            tasks = workflow.getAsJsonArray("tasks");
        } else {
            // Schema 1.5+ format: tasks under specification
            JsonObject specification = workflow.getAsJsonObject("specification");
            tasks = specification.getAsJsonArray("tasks");
        }

        // Build execution data map if available
        Map<String, JsonObject> execMap = new HashMap<String, JsonObject>();
        if (options.isUseExecutionData()) {
            if (isLegacySchema) {
                // Schema 1.4: runtime data is embedded in task objects
                for (JsonElement elem : tasks) {
                    JsonObject task = elem.getAsJsonObject();
                    String id = task.has("id") ? task.get("id").getAsString() :
                               (task.has("name") ? task.get("name").getAsString() : null);
                    if (id != null) {
                        execMap.put(id, task);
                    }
                }
            } else if (workflow.has("execution")) {
                // Schema 1.5+: separate execution section
                JsonObject execution = workflow.getAsJsonObject("execution");
                if (execution.has("tasks")) {
                    JsonArray execTasks = execution.getAsJsonArray("tasks");
                    for (JsonElement elem : execTasks) {
                        JsonObject execTask = elem.getAsJsonObject();
                        String id = execTask.get("id").getAsString();
                        execMap.put(id, execTask);
                    }
                }
            }
        }

        // Phase 1: Create all activities
        Map<String, WorkflowActivity> taskMap = new HashMap<String, WorkflowActivity>();
        Map<String, Integer> taskIdxMap = new HashMap<String, Integer>();
        int n = tasks.size();

        for (int i = 0; i < n; i++) {
            JsonObject task = tasks.get(i).getAsJsonObject();
            // Schema 1.4 may use "name" as identifier, Schema 1.5+ uses "id"
            String taskId = task.has("id") ? task.get("id").getAsString() :
                           (task.has("name") ? task.get("name").getAsString() : "task_" + i);

            // Get runtime from execution data
            double runtime = options.getDefaultRuntime();
            if (execMap.containsKey(taskId)) {
                JsonObject execData = execMap.get(taskId);
                if (execData.has("runtimeInSeconds")) {
                    runtime = execData.get("runtimeInSeconds").getAsDouble();
                }
            }

            // Fit distribution
            Distribution dist = fitDistribution(runtime, options);

            // Create activity
            WorkflowActivity act = wf.addActivity(taskId, dist);
            taskMap.put(taskId, act);
            taskIdxMap.put(taskId, i);

            // Store metadata if requested
            if (options.isStoreMetadata()) {
                Map<String, Object> metadata = extractMetadata(task, execMap.get(taskId));
                act.setMetadata(metadata);
            }
        }

        // Phase 2: Build adjacency structures
        List<List<Integer>> adjList = new ArrayList<List<Integer>>();
        List<List<Integer>> predList = new ArrayList<List<Integer>>();
        int[] inDeg = new int[n];
        int[] outDeg = new int[n];

        for (int i = 0; i < n; i++) {
            adjList.add(new ArrayList<Integer>());
            predList.add(new ArrayList<Integer>());
        }

        for (int i = 0; i < n; i++) {
            JsonObject task = tasks.get(i).getAsJsonObject();

            // Process children
            if (task.has("children") && !task.get("children").isJsonNull()) {
                JsonArray children = task.getAsJsonArray("children");
                for (JsonElement childElem : children) {
                    String childId = childElem.getAsString();
                    if (taskIdxMap.containsKey(childId)) {
                        int childIdx = taskIdxMap.get(childId);
                        adjList.get(i).add(childIdx);
                        predList.get(childIdx).add(i);
                        outDeg[i]++;
                        inDeg[childIdx]++;
                    }
                }
            }
        }

        // Phase 3: Detect and add precedences
        addPrecedences(wf, tasks, taskMap, adjList, inDeg, outDeg, predList);

        return wf;
    }

    private static Distribution fitDistribution(double runtime, WfCommonsOptions options) {
        if (runtime <= GlobalConstants.FineTol) {
            return new Immediate();
        }

        switch (options.getDistributionType()) {
            case DET:
                return new Det(runtime);
            case APH:
                return APH.fitMeanAndSCV(runtime, options.getDefaultSCV());
            case HYPEREXP:
                if (options.getDefaultSCV() > 1.0) {
                    return HyperExp.fitMeanAndSCV(runtime, options.getDefaultSCV());
                }
                return Exp.fitMean(runtime);
            case EXP:
            default:
                return Exp.fitMean(runtime);
        }
    }

    private static Map<String, Object> extractMetadata(JsonObject task, JsonObject execData) {
        Map<String, Object> metadata = new HashMap<String, Object>();

        metadata.put("taskId", task.get("id").getAsString());

        if (task.has("name") && !task.get("name").isJsonNull()) {
            metadata.put("name", task.get("name").getAsString());
        }
        if (task.has("inputFiles") && !task.get("inputFiles").isJsonNull()) {
            List<String> files = new ArrayList<String>();
            for (JsonElement elem : task.getAsJsonArray("inputFiles")) {
                files.add(elem.getAsString());
            }
            metadata.put("inputFiles", files);
        }
        if (task.has("outputFiles") && !task.get("outputFiles").isJsonNull()) {
            List<String> files = new ArrayList<String>();
            for (JsonElement elem : task.getAsJsonArray("outputFiles")) {
                files.add(elem.getAsString());
            }
            metadata.put("outputFiles", files);
        }

        if (execData != null) {
            String[] execFields = {"executedAt", "coreCount", "avgCPU", "readBytes",
                                   "writtenBytes", "memoryInBytes", "energyInKWh",
                                   "avgPowerInW", "priority"};
            for (String field : execFields) {
                if (execData.has(field) && !execData.get(field).isJsonNull()) {
                    JsonElement value = execData.get(field);
                    if (value.isJsonPrimitive()) {
                        if (value.getAsJsonPrimitive().isNumber()) {
                            metadata.put(field, value.getAsDouble());
                        } else {
                            metadata.put(field, value.getAsString());
                        }
                    }
                }
            }
            if (execData.has("machines") && !execData.get("machines").isJsonNull()) {
                List<String> machines = new ArrayList<String>();
                for (JsonElement elem : execData.getAsJsonArray("machines")) {
                    machines.add(elem.getAsString());
                }
                metadata.put("machines", machines);
            }
            if (execData.has("command") && !execData.get("command").isJsonNull()) {
                metadata.put("command", execData.getAsJsonObject("command").toString());
            }
        }

        return metadata;
    }

    private static void addPrecedences(Workflow wf, JsonArray tasks,
                                       Map<String, WorkflowActivity> taskMap,
                                       List<List<Integer>> adjList,
                                       int[] inDeg, int[] outDeg,
                                       List<List<Integer>> predList) {
        int n = tasks.size();
        boolean[][] processed = new boolean[n][n];

        // Step 1: Identify and process fork-join pairs
        for (int i = 0; i < n; i++) {
            if (outDeg[i] > 1) {
                // Potential AND-fork
                List<Integer> children = adjList.get(i);

                // Find common join point
                Integer joinPoint = findCommonJoin(children, adjList, inDeg, n);

                if (joinPoint != null) {
                    // Valid fork-join structure
                    String preTaskId = tasks.get(i).getAsJsonObject().get("id").getAsString();
                    String postTaskId = tasks.get(joinPoint).getAsJsonObject().get("id").getAsString();

                    WorkflowActivity preAct = taskMap.get(preTaskId);
                    List<WorkflowActivity> postActs = new ArrayList<WorkflowActivity>();
                    for (int childIdx : children) {
                        String childId = tasks.get(childIdx).getAsJsonObject().get("id").getAsString();
                        postActs.add(taskMap.get(childId));
                    }
                    WorkflowActivity postAct = taskMap.get(postTaskId);

                    // Add fork and join
                    wf.addPrecedence(Workflow.AndFork(preAct, postActs));
                    wf.addPrecedence(Workflow.AndJoin(postActs, postAct));

                    // Mark edges as processed
                    for (int childIdx : children) {
                        processed[i][childIdx] = true;
                        processed[childIdx][joinPoint] = true;
                    }
                }
            }
        }

        // Step 2: Add remaining edges as serial connections
        for (int i = 0; i < n; i++) {
            for (int j : adjList.get(i)) {
                if (!processed[i][j]) {
                    String preTaskId = tasks.get(i).getAsJsonObject().get("id").getAsString();
                    String postTaskId = tasks.get(j).getAsJsonObject().get("id").getAsString();

                    WorkflowActivity preAct = taskMap.get(preTaskId);
                    WorkflowActivity postAct = taskMap.get(postTaskId);

                    wf.addPrecedence(Workflow.Serial(preAct, postAct));
                }
            }
        }
    }

    private static Integer findCommonJoin(List<Integer> children,
                                          List<List<Integer>> adjList,
                                          int[] inDeg, int n) {
        if (children.isEmpty()) {
            return null;
        }

        // Use BFS from each child to find reachable nodes
        List<Set<Integer>> reachableSets = new ArrayList<Set<Integer>>();
        for (int child : children) {
            reachableSets.add(getReachableNodes(child, adjList, n));
        }

        // Find intersection
        Set<Integer> common = new HashSet<Integer>(reachableSets.get(0));
        for (int i = 1; i < reachableSets.size(); i++) {
            common.retainAll(reachableSets.get(i));
        }

        if (common.isEmpty()) {
            return null;
        }

        // Find the first common node that has all children as direct predecessors
        for (int node : common) {
            if (inDeg[node] >= children.size()) {
                boolean allDirectChild = true;
                for (int child : children) {
                    if (!adjList.get(child).contains(node)) {
                        allDirectChild = false;
                        break;
                    }
                }
                if (allDirectChild) {
                    return node;
                }
            }
        }

        return null;
    }

    private static Set<Integer> getReachableNodes(int start,
                                                   List<List<Integer>> adjList,
                                                   int n) {
        Set<Integer> visited = new HashSet<Integer>();
        Queue<Integer> queue = new LinkedList<Integer>();
        visited.add(start);
        queue.add(start);
        Set<Integer> reachable = new HashSet<Integer>();

        while (!queue.isEmpty()) {
            int current = queue.poll();
            for (int next : adjList.get(current)) {
                if (!visited.contains(next)) {
                    visited.add(next);
                    queue.add(next);
                    reachable.add(next);
                }
            }
        }

        return reachable;
    }
}
