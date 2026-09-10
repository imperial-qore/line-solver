/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import com.google.gson.*;

import jline.GlobalConstants;
import jline.lang.*;
import jline.lang.constant.*;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.HeteroSchedPolicy;
import jline.lang.constant.JoinStrategy;
import jline.lang.constant.ServerType;
import jline.lang.layered.*;
import jline.lang.nodes.*;
import jline.lang.processes.*;
import jline.lang.reward.Reward;
import jline.lang.reward.RewardDescriptor;
import jline.lang.reward.RewardFunction;
import jline.lang.sections.ClassSwitcher;
import jline.lang.sections.Forker;
import jline.lang.sections.Joiner;
import jline.lang.sections.PollingServer;
import jline.lang.workflow.Workflow;
import jline.lang.workflow.WorkflowActivity;
import jline.util.matrix.Matrix;
import jline.util.SerializableFunction;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static jline.io.InputOutput.line_warning;
import static jline.io.InputOutput.mfilename;

/**
 * Provides save/load functionality for LINE queueing network models to/from JSON format.
 *
 * <p>Supports both {@link Network} (queueing network) and {@link LayeredNetwork} (layered
 * queueing network) models. The JSON format follows the line-model.schema.json specification.</p>
 *
 * <p>Usage:
 * <pre>
 *   // Save a Network model
 *   LineModelIO.save(model, "mymodel.json");
 *
 *   // Save a LayeredNetwork model
 *   LineModelIO.save(lqnModel, "mylqn.json");
 *
 *   // Load any model
 *   Object model = LineModelIO.load("mymodel.json");
 *   if (model instanceof Network) { ... }
 *   if (model instanceof LayeredNetwork) { ... }
 * </pre>
 * </p>
 */
public class LineModelIO {

    private static final String FORMAT_NAME = "line-model";
    private static final String FORMAT_VERSION = "1.0";

    // ========================================================================
    // SAVE: Network
    // ========================================================================

    /**
     * Builds a {@link JsonObject} representation of a {@link Network} model.
     *
     * @param model the queueing network model to serialize
     * @return the JSON document representing the model
     */
    public static JsonObject toJsonObject(Network model) {
        JsonObject doc = new JsonObject();
        doc.addProperty("format", FORMAT_NAME);
        doc.addProperty("version", FORMAT_VERSION);

        JsonObject modelObj = new JsonObject();
        modelObj.addProperty("type", "Network");
        modelObj.addProperty("name", model.getName());

        // Nodes
        modelObj.add("nodes", serializeNetworkNodes(model));

        // Classes
        modelObj.add("classes", serializeNetworkClasses(model));

        // Routing
        JsonObject routingResult = serializeNetworkRouting(model);
        // Move routingStrategies/routingWeights to model level for cross-codebase compatibility
        JsonObject routingStrategies = null;
        JsonObject routingWeights = null;
        JsonObject routingParams = null;
        if (routingResult.has("routingStrategies")) {
            routingStrategies = routingResult.getAsJsonObject("routingStrategies");
            routingResult.remove("routingStrategies");
        }
        if (routingResult.has("routingWeights")) {
            routingWeights = routingResult.getAsJsonObject("routingWeights");
            routingResult.remove("routingWeights");
        }
        if (routingResult.has("routingParams")) {
            routingParams = routingResult.getAsJsonObject("routingParams");
            routingResult.remove("routingParams");
        }
        modelObj.add("routing", routingResult);
        if (routingStrategies != null) {
            modelObj.add("routingStrategies", routingStrategies);
        }
        if (routingWeights != null) {
            modelObj.add("routingWeights", routingWeights);
        }
        if (routingParams != null) {
            modelObj.add("routingParams", routingParams);
        }

        // Krzesinski state-dependent routing. Every center travels by NODE NAME,
        // so the block is language independent and a node reordering on either
        // side cannot shift a center. See _kb/16-state-dependent-routing.md
        JsonObject sdrObj = serializeStateDepRouting(model);
        if (sdrObj != null) {
            modelObj.add("stateDepRouting", sdrObj);
        }

        // Global (Whittle) dependence phi(n): materialized over the lattice of the
        // WHOLE network state, unlike the per-station classDependence /
        // jointDependence blocks. Only the (station,class) slots a class can
        // actually occupy carry a coordinate, which keeps the lattice finite.
        if (model.getGlobalDependence() != null) {
            modelObj.add("globalDependence", serializeGlobalDependence(model));
        }

        // Finite capacity regions
        List<Region> regions = model.getRegions();
        if (regions != null && !regions.isEmpty()) {
            JsonArray fcrArr = new JsonArray();
            List<JobClass> classes = model.getClasses();
            for (Region region : regions) {
                JsonObject rj = new JsonObject();
                rj.addProperty("name", region.getName());
                // Stations
                JsonArray stationsArr = new JsonArray();
                for (Node n : region.getNodes()) {
                    JsonObject sj = new JsonObject();
                    sj.addProperty("node", n.getName());
                    // Per-class classWeight
                    JsonObject cwObj = new JsonObject();
                    for (JobClass jc : classes) {
                        double w = region.getClassWeight(jc);
                        if (w != 1.0) {
                            cwObj.addProperty(jc.getName(), w);
                        }
                    }
                    if (cwObj.size() > 0) sj.add("classWeight", cwObj);
                    // Per-class classSize
                    JsonObject csObj = new JsonObject();
                    for (JobClass jc : classes) {
                        double sz = region.getClassSize(jc);
                        if (sz != 1.0) {
                            csObj.addProperty(jc.getName(), sz);
                        }
                    }
                    if (csObj.size() > 0) sj.add("classSize", csObj);
                    stationsArr.add(sj);
                }
                rj.add("stations", stationsArr);
                // globalMaxJobs
                int gmj = region.getGlobalMaxJobs();
                if (gmj >= 0) {
                    rj.addProperty("globalMaxJobs", gmj);
                }
                // globalMaxMemory
                double gmm = region.getGlobalMaxMemory();
                if (gmm >= 0) {
                    rj.addProperty("globalMaxMemory", gmm);
                }
                // classMaxJobs
                JsonObject cmjObj = new JsonObject();
                for (JobClass jc : classes) {
                    int cmj = region.getClassMaxJobs(jc);
                    if (cmj >= 0 && cmj < Integer.MAX_VALUE) {
                        cmjObj.addProperty(jc.getName(), cmj);
                    }
                }
                if (cmjObj.size() > 0) rj.add("classMaxJobs", cmjObj);
                // classMaxMemory: see _kb/09-ldes-and-cache.md ("More model.json writer/reader notes")
                JsonObject cmmObj = new JsonObject();
                for (JobClass jc : classes) {
                    int cmm = region.getClassMaxMemory(jc);
                    if (cmm >= 0 && cmm < Integer.MAX_VALUE) {
                        cmmObj.addProperty(jc.getName(), cmm);
                    }
                }
                if (cmmObj.size() > 0) rj.add("classMaxMemory", cmmObj);
                // dropRule (only non-default entries; default is WaitingQueue)
                JsonObject drObj = new JsonObject();
                for (JobClass jc : classes) {
                    DropStrategy ds = region.getDropStrategy(jc);
                    if (ds != DropStrategy.WaitingQueue) {
                        drObj.addProperty(jc.getName(), dropStrategyToStr(ds));
                    }
                }
                if (drObj.size() > 0) rj.add("dropRule", drObj);
                // linear admission constraints A*n<=b, row-major A / flat b: see _kb/09-ldes-and-cache.md
                if (region.hasLinearConstraints()) {
                    Matrix[] lincon = region.getLinearConstraints();
                    rj.add("constraintA", matrixToJson2D(lincon[0]));
                    rj.add("constraintB", matrixToJsonRowVector(lincon[1]));
                }
                fcrArr.add(rj);
            }
            modelObj.add("finiteCapacityRegions", fcrArr);
        }

        JsonArray rewardsArr = serializeNetworkRewards(model);
        if (rewardsArr.size() > 0) {
            modelObj.add("rewards", rewardsArr);
        }

        doc.add("model", modelObj);
        return doc;
    }

    /**
     * Serializes the model's rewards in the declarative {name, type, node, class} form.
     *
     * Only rewards created through a {@link jline.lang.reward.Reward} template carry the
     * structural metadata needed to reproduce them. A reward defined from a bare lambda
     * (or via Reward.custom) is not reproducible from JSON: warn and omit it rather than
     * emit a reward that would be wrong on reload.
     */
    private static JsonArray serializeNetworkRewards(Network model) {
        JsonArray rewardsArr = new JsonArray();
        Map<String, RewardFunction> rewards = model.getRewards();
        if (rewards == null || rewards.isEmpty()) {
            return rewardsArr;
        }
        // sorted by name for determinism (HashMap iteration order is unstable): see _kb/09-ldes-and-cache.md
        List<String> rewardNames = new ArrayList<String>(rewards.keySet());
        Collections.sort(rewardNames);
        for (String name : rewardNames) {
            RewardFunction fn = rewards.get(name);
            if (!(fn instanceof RewardDescriptor)) {
                System.err.println("Warning: reward \"" + name + "\" is defined by a bare reward function and "
                        + "cannot be serialized to JSON; it is omitted from the saved model. Use a "
                        + "jline.lang.reward.Reward template (queueLength/utilization/blocking) for a "
                        + "serializable reward.");
                continue;
            }
            RewardDescriptor descriptor = (RewardDescriptor) fn;
            if (descriptor.getKind() == RewardDescriptor.Kind.Custom) {
                System.err.println("Warning: reward \"" + name + "\" is a custom reward wrapping an arbitrary "
                        + "function and cannot be serialized to JSON; it is omitted from the saved model.");
                continue;
            }
            if (descriptor.getNode() == null) {
                System.err.println("Warning: reward \"" + name + "\" of type " + descriptor.getKind()
                        + " has no associated node and cannot be serialized to JSON; it is omitted from the "
                        + "saved model.");
                continue;
            }
            JsonObject rj = new JsonObject();
            rj.addProperty("name", name);
            rj.addProperty("type", descriptor.getKind().name());
            rj.addProperty("node", descriptor.getNode().getName());
            if (descriptor.getJobClass() != null) {
                rj.addProperty("class", descriptor.getJobClass().getName());
            }
            rewardsArr.add(rj);
        }
        return rewardsArr;
    }

    /**
     * Restores declarative rewards onto a freshly loaded Network.
     */
    private static void loadNetworkRewards(Network model, JsonObject modelObj,
                                           Map<String, Node> nodeMap, Map<String, JobClass> classMap) {
        if (!modelObj.has("rewards")) {
            return;
        }
        JsonArray rewardsArr = modelObj.getAsJsonArray("rewards");
        for (int i = 0; i < rewardsArr.size(); i++) {
            JsonObject rw = rewardsArr.get(i).getAsJsonObject();
            if (!rw.has("name") || !rw.has("type")) {
                System.err.println("Warning: ignoring a reward entry without a \"name\" or \"type\" field.");
                continue;
            }
            String name = rw.get("name").getAsString();
            String type = rw.get("type").getAsString();
            String nodeName = rw.has("node") ? rw.get("node").getAsString() : null;
            Node node = (nodeName == null) ? null : nodeMap.get(nodeName);
            if (node == null) {
                System.err.println("Warning: reward \"" + name + "\" refers to node \"" + nodeName
                        + "\", which is not defined in this model; the reward is ignored.");
                continue;
            }
            JobClass jobClass = null;
            if (rw.has("class")) {
                String className = rw.get("class").getAsString();
                jobClass = classMap.get(className);
                if (jobClass == null) {
                    System.err.println("Warning: reward \"" + name + "\" refers to class \"" + className
                            + "\", which is not defined in this model; the reward is ignored.");
                    continue;
                }
            }
            if ("QLen".equals(type)) {
                model.setReward(name, Reward.queueLength(node, jobClass));
            } else if ("Util".equals(type)) {
                model.setReward(name, Reward.utilization(node, jobClass));
            } else if ("Blocking".equals(type)) {
                model.setReward(name, Reward.blocking(node));
            } else {
                System.err.println("Warning: reward \"" + name + "\" has type \"" + type
                        + "\", for which no reward template is implemented; the reward is ignored.");
            }
        }
    }

    /**
     * Rewrites every non-finite number in a model document into the spelling
     * the wire format uses, then writes it.
     * <p>
     * JSON HAS NO INFINITY LITERAL. Gson's serializeSpecialFloatingPointValues
     * writes a bare {@code Infinity} / {@code NaN} anyway, which no strict
     * parser takes: {@code linemodel_save} (the reference) writes an infinite
     * scalar as the STRING "Infinity" / "-Infinity" and a NaN as {@code null},
     * for every numeric field rather than a named few, and the C++ reader
     * {@code num_from_json} decodes exactly that pair. The spelling is
     * reachable whenever a model carries its declared state, since an open
     * class holds an infinite population and the Source's {@code initialState}
     * row then carries the sentinel. Left as a Gson literal, such a document is
     * refused by nlohmann at PARSE time, so the native {@code common/ldes}
     * engine cannot read it and SolverLDES silently falls back to the JVM
     * engine -- a different sample path, not a different spelling.
     * <p>
     * The builder no longer permits the special values, so a non-finite that
     * escapes this rewrite is reported here rather than written.
     *
     * @param doc      the model document
     * @param filename the output file path
     * @throws IOException if the file cannot be written
     */
    private static void writeJson(JsonObject doc, String filename) throws IOException {
        Gson gson = new GsonBuilder().setPrettyPrinting().create();
        Writer writer = new BufferedWriter(new FileWriter(filename));
        try {
            gson.toJson(wireNonFinite(doc), writer);
        } finally {
            writer.close();
        }
    }

    /**
     * The document with every non-finite number replaced by its wire spelling.
     * Containers are rewritten in place; see {@link #writeJson}.
     *
     * @param el the element to rewrite
     * @return the rewritten element
     */
    private static JsonElement wireNonFinite(JsonElement el) {
        if (el == null || el.isJsonNull()) {
            return el;
        }
        if (el.isJsonObject()) {
            JsonObject obj = el.getAsJsonObject();
            List<String> keys = new ArrayList<String>(obj.keySet());
            for (String key : keys) {
                obj.add(key, wireNonFinite(obj.get(key)));
            }
            return obj;
        }
        if (el.isJsonArray()) {
            JsonArray arr = el.getAsJsonArray();
            for (int i = 0; i < arr.size(); i++) {
                arr.set(i, wireNonFinite(arr.get(i)));
            }
            return arr;
        }
        if (el.isJsonPrimitive() && el.getAsJsonPrimitive().isNumber()) {
            double v = el.getAsDouble();
            if (Double.isNaN(v)) {
                return JsonNull.INSTANCE;
            }
            if (Double.isInfinite(v)) {
                return new JsonPrimitive(v > 0 ? "Infinity" : "-Infinity");
            }
        }
        return el;
    }

    /**
     * Saves a {@link Network} model to a JSON file.
     *
     * @param model    the queueing network model to save
     * @param filename the output file path (should end with .json)
     * @throws IOException if the file cannot be written
     */
    public static void save(Network model, String filename) throws IOException {
        JsonObject doc = toJsonObject(model);
        writeJson(doc, filename);
    }

    // ========================================================================
    // SAVE: LayeredNetwork
    // ========================================================================

    /**
     * Builds a {@link JsonObject} representation of a {@link LayeredNetwork} model.
     *
     * @param model the layered queueing network model to serialize
     * @return the JSON document representing the model
     */
    public static JsonObject toJsonObject(LayeredNetwork model) {
        JsonObject doc = new JsonObject();
        doc.addProperty("format", FORMAT_NAME);
        doc.addProperty("version", FORMAT_VERSION);

        JsonObject modelObj = new JsonObject();
        modelObj.addProperty("type", "LayeredNetwork");
        modelObj.addProperty("name", model.getName());

        // Hosts
        modelObj.add("hosts", serializeHosts(model));
        // Tasks
        modelObj.add("tasks", serializeTasks(model));
        // Entries
        modelObj.add("entries", serializeEntries(model));
        // Activities
        modelObj.add("activities", serializeActivities(model));
        // Precedences
        modelObj.add("precedences", serializePrecedences(model));

        doc.add("model", modelObj);
        return doc;
    }

    /**
     * Saves a {@link LayeredNetwork} model to a JSON file.
     *
     * @param model    the layered queueing network model to save
     * @param filename the output file path (should end with .json)
     * @throws IOException if the file cannot be written
     */
    public static void save(LayeredNetwork model, String filename) throws IOException {
        JsonObject doc = toJsonObject(model);
        writeJson(doc, filename);
    }

    // ========================================================================
    // SAVE: Workflow
    // ========================================================================

    /**
     * Saves a {@link Workflow} model to a JSON file.
     *
     * @param model    the workflow model to save
     * @param filename the output file path (should end with .json)
     * @throws IOException if the file cannot be written
     */
    public static void save(Workflow model, String filename) throws IOException {
        JsonObject doc = new JsonObject();
        doc.addProperty("format", FORMAT_NAME);
        doc.addProperty("version", FORMAT_VERSION);

        JsonObject modelObj = new JsonObject();
        modelObj.addProperty("type", "Workflow");
        modelObj.addProperty("name", model.getName());

        // Activities
        JsonArray actsArr = new JsonArray();
        for (WorkflowActivity act : model.getActivities()) {
            JsonObject actObj = new JsonObject();
            actObj.addProperty("name", act.getName());
            if (act.getHostDemand() != null) {
                actObj.add("hostDemand", serializeDistribution(act.getHostDemand()));
            }
            actsArr.add(actObj);
        }
        modelObj.add("activities", actsArr);

        // Precedences
        JsonArray precsArr = new JsonArray();
        for (ActivityPrecedence prec : model.getPrecedences()) {
            JsonObject precObj = new JsonObject();
            JsonArray preActs = new JsonArray();
            for (String name : prec.getPreActs()) {
                preActs.add(name);
            }
            precObj.add("preActs", preActs);
            JsonArray postActs = new JsonArray();
            for (String name : prec.getPostActs()) {
                postActs.add(name);
            }
            precObj.add("postActs", postActs);
            precObj.addProperty("preType", prec.getPreType());
            precObj.addProperty("postType", prec.getPostType());
            if (prec.getPreParams() != null && !prec.getPreParams().isEmpty()) {
                precObj.add("preParams", matrixToJsonArray(prec.getPreParams()));
            }
            if (prec.getPostParams() != null && !prec.getPostParams().isEmpty()) {
                precObj.add("postParams", matrixToJsonArray(prec.getPostParams()));
            }
            precsArr.add(precObj);
        }
        modelObj.add("precedences", precsArr);

        doc.add("model", modelObj);

        writeJson(doc, filename);
    }

    // ========================================================================
    // SAVE: Environment
    // ========================================================================

    /**
     * Saves an {@link Environment} model to a JSON file.
     *
     * @param env      the environment model to save
     * @param filename the output file path (should end with .json)
     * @throws IOException if the file cannot be written
     */
    public static void save(Environment env, String filename) throws IOException {
        JsonObject doc = new JsonObject();
        doc.addProperty("format", FORMAT_NAME);
        doc.addProperty("version", FORMAT_VERSION);

        JsonObject modelObj = new JsonObject();
        modelObj.addProperty("type", "Environment");
        modelObj.addProperty("name", env.getName());

        int nStages = env.getNumberOfStages();
        modelObj.addProperty("numStages", nStages);

        // Stages
        JsonArray stagesArr = new JsonArray();
        for (int i = 0; i < nStages; i++) {
            String stageName = env.getStageName(i);
            if (stageName == null) {
                continue;
            }
            JsonObject stageObj = new JsonObject();
            stageObj.addProperty("name", stageName);
            // The stage TYPE, which every reader already looks for and no writer
            // emitted: an Environment round-tripped through JSON came back with its
            // stage types blanked, so getStageTable and any consumer keying off
            // UP/DOWN read a different environment from the one that was saved.
            String stageType = env.getStageType(i);
            if (stageType != null && stageType.length() > 0) {
                stageObj.addProperty("type", stageType);
            }
            // Serialize the stage's Network model using existing Network serialization
            Network stageModel = env.getEnsemble().size() > i ? env.getModel(i) : null;
            if (stageModel != null) {
                JsonObject netJson = new JsonObject();
                netJson.addProperty("type", "Network");
                netJson.addProperty("name", stageModel.getName());
                netJson.add("nodes", serializeNetworkNodes(stageModel));
                netJson.add("classes", serializeNetworkClasses(stageModel));
                JsonObject routingResult = serializeNetworkRouting(stageModel);
                JsonObject routingStrategies = null;
                JsonObject routingWeights = null;
                if (routingResult.has("routingStrategies")) {
                    routingStrategies = routingResult.getAsJsonObject("routingStrategies");
                    routingResult.remove("routingStrategies");
                }
                if (routingResult.has("routingWeights")) {
                    routingWeights = routingResult.getAsJsonObject("routingWeights");
                    routingResult.remove("routingWeights");
                }
                netJson.add("routing", routingResult);
                if (routingStrategies != null) {
                    netJson.add("routingStrategies", routingStrategies);
                }
                if (routingWeights != null) {
                    netJson.add("routingWeights", routingWeights);
                }
                stageObj.add("model", netJson);
            }
            stagesArr.add(stageObj);
        }
        modelObj.add("stages", stagesArr);

        // Transitions
        JsonArray transArr = new JsonArray();
        for (int i = 0; i < nStages; i++) {
            for (int j = 0; j < nStages; j++) {
                if (env.env[i][j] != null) {
                    JsonObject transObj = new JsonObject();
                    transObj.addProperty("from", i);
                    transObj.addProperty("to", j);
                    transObj.add("distribution", serializeDistribution(env.env[i][j]));
                    transArr.add(transObj);
                }
            }
        }
        modelObj.add("transitions", transArr);

        // nodeFailures: see _kb/09-ldes-and-cache.md / _kb/12-interfaces-and-docs.md ("nodeFailures")
        JsonArray nodeFailuresArr = serializeNodeFailures(env);
        if (nodeFailuresArr.size() > 0) {
            modelObj.add("nodeFailures", nodeFailuresArr);
        }

        doc.add("model", modelObj);

        writeJson(doc, filename);
    }

    /**
     * Serializes the environment's node breakdown/repair descriptors.
     */
    private static JsonArray serializeNodeFailures(Environment env) {
        JsonArray nfArr = new JsonArray();
        List<Environment.NodeFailure> nodeFailures = env.getNodeFailures();
        if (nodeFailures == null || nodeFailures.isEmpty()) {
            return nfArr;
        }
        for (Environment.NodeFailure nf : nodeFailures) {
            JsonObject breakdown = serializeDistribution(nf.breakdown);
            if (breakdown == null) {
                System.err.println("Warning: node failure on \"" + nf.node + "\" has a breakdown distribution "
                        + "that cannot be serialized; the nodeFailures entry is omitted.");
                continue;
            }
            JsonObject downService = serializeDistribution(nf.downService);
            if (downService == null) {
                System.err.println("Warning: node failure on \"" + nf.node + "\" has a down-service "
                        + "distribution that cannot be serialized; the nodeFailures entry is omitted.");
                continue;
            }
            JsonObject nj = new JsonObject();
            nj.addProperty("node", nf.node);
            nj.add("breakdownRate", breakdown);
            if (nf.repair != null) {
                JsonObject repair = serializeDistribution(nf.repair);
                if (repair == null) {
                    System.err.println("Warning: node failure on \"" + nf.node + "\" has a repair distribution "
                            + "that cannot be serialized; the nodeFailures entry is omitted.");
                    continue;
                }
                nj.add("repairRate", repair);
            }
            nj.add("downService", downService);
            if (Environment.RESET_POLICY_CUSTOM.equals(nf.breakdownResetPolicy)) {
                System.err.println("Warning: node failure on \"" + nf.node + "\" uses a custom breakdown reset "
                        + "function, which cannot be serialized to JSON; the saved model falls back to the "
                        + "\"keep\" policy on reload.");
            } else {
                nj.addProperty("breakdownResetPolicy", nf.breakdownResetPolicy);
            }
            if (nf.repairResetPolicy != null && !nf.repairResetPolicy.isEmpty()) {
                if (Environment.RESET_POLICY_CUSTOM.equals(nf.repairResetPolicy)) {
                    System.err.println("Warning: node failure on \"" + nf.node + "\" uses a custom repair reset "
                            + "function, which cannot be serialized to JSON; the saved model falls back to the "
                            + "\"keep\" policy on reload.");
                } else {
                    nj.addProperty("repairResetPolicy", nf.repairResetPolicy);
                }
            }
            nfArr.add(nj);
        }
        return nfArr;
    }

    // ========================================================================
    // LOAD
    // ========================================================================

    /**
     * Loads a model from a JSON file.
     *
     * @param filename the input file path
     * @return a {@link Network} or {@link LayeredNetwork} depending on the model type
     * @throws IOException if the file cannot be read or the format is invalid
     */
    public static Object load(String filename) throws IOException {
        Reader reader = new BufferedReader(new FileReader(filename));
        JsonObject doc;
        try {
            doc = JsonParser.parseReader(reader).getAsJsonObject();
        } finally {
            reader.close();
        }

        String format = doc.has("format") ? doc.get("format").getAsString() : "";
        if (!FORMAT_NAME.equals(format)) {
            throw new IOException("Unknown format: " + format + " (expected " + FORMAT_NAME + ")");
        }

        JsonObject modelObj = doc.getAsJsonObject("model");
        String type = modelObj.get("type").getAsString();

        if ("Network".equals(type)) {
            return loadNetwork(modelObj);
        } else if ("LayeredNetwork".equals(type)) {
            return loadLayeredNetwork(modelObj);
        } else if ("Workflow".equals(type)) {
            return loadWorkflow(modelObj);
        } else if ("Environment".equals(type)) {
            return loadEnvironment(modelObj);
        } else {
            throw new IOException("Unknown model type: " + type);
        }
    }

    // ========================================================================
    // NETWORK SERIALIZATION HELPERS
    // ========================================================================

    private static JsonArray serializeNetworkNodes(Network model) {
        JsonArray nodesArr = new JsonArray();
        List<Node> nodes = model.getNodes();
        List<JobClass> classes = model.getClasses();

        // A STATE ON A STRICT SUBSET OF THE STATEFUL NODES DOES NOT TRAVEL: it is
        // not an initialization, and getState() runs initDefault() for exactly
        // that reason, so this side answers for the DEFAULT marking. A document
        // naming only the node the caller moved made the reader combine that row
        // with default markings for the rest -- on Think -> Q1 with 2 jobs and Q1
        // alone set to 2, the C++ read (Think=2, Q1=2) and answered
        // getProbSysAggr 0 for a joint state holding 4 of 2 jobs, against 0.4
        // here. A PAS placement is the exception initDefault() itself makes,
        // its ordering being a required input rather than a default.
        boolean fullyInitialized = model.hasInitState();

        for (Node node : nodes) {
            boolean nodeStateTravels = fullyInitialized
                    || (node instanceof Queue && ((Queue) node).getSchedStrategy() == SchedStrategy.PAS);
            // Skip auto-added ClassSwitch nodes (recreated by link())
            if (node instanceof ClassSwitch && ((ClassSwitch) node).autoAdded) {
                continue;
            }

            JsonObject nodeObj = new JsonObject();
            nodeObj.addProperty("name", node.getName());

            if (node instanceof Source) {
                nodeObj.addProperty("type", "Source");
                Source src = (Source) node;
                JsonObject arrivals = new JsonObject();
                for (JobClass jc : classes) {
                    Distribution dist = src.getArrivalDistribution(jc);
                    if (dist != null && !(dist instanceof Disabled)) {
                        arrivals.add(jc.getName(), serializeDistribution(dist));
                    }
                }
                nodeObj.add("service", arrivals);
                // arrivalBatch wire key: see _kb/09-ldes-and-cache.md ("Batch arrivals")
                JsonObject arrivalBatch = new JsonObject();
                for (JobClass jc : classes) {
                    jline.lang.processes.DiscreteDistribution batch = src.getArrivalBatch(jc);
                    if (batch != null) {
                        arrivalBatch.add(jc.getName(), serializeDistribution(batch));
                    }
                }
                if (arrivalBatch.size() > 0) {
                    nodeObj.add("arrivalBatch", arrivalBatch);
                }
                // Marked (MMAP) arrival binding: class names ordered by mark
                if (src.getMarkedClasses() != null) {
                    JsonArray markedArr = new JsonArray();
                    for (JobClass mc : src.getMarkedClasses()) {
                        markedArr.add(mc.getName());
                    }
                    nodeObj.add("markedClasses", markedArr);
                }

            } else if (node instanceof Sink) {
                nodeObj.addProperty("type", "Sink");

            } else if (node instanceof Queue) {
                // Delay extends Queue, shares this branch: see _kb/09-ldes-and-cache.md ("More model.json writer/reader notes")
                Queue queue = (Queue) node;
                boolean isDelay = node instanceof Delay;
                nodeObj.addProperty("type", isDelay ? "Delay" : "Queue");
                if (!isDelay) {
                    nodeObj.addProperty("scheduling", queue.getSchedStrategy().toString());
                    int nServers = queue.getNumberOfServers();
                    if (nServers != Integer.MAX_VALUE && nServers > 1) {
                        nodeObj.addProperty("servers", nServers);
                    }
                }
                JsonObject services = new JsonObject();
                for (JobClass jc : classes) {
                    Distribution dist = queue.getService(jc);
                    if (dist != null && !(dist instanceof Disabled)) {
                        services.add(jc.getName(), serializeDistribution(dist));
                    }
                }
                nodeObj.add("service", services);

                // Scheduling parameters (for DPS, GPS, etc.)
                boolean hasSchedPar = false;
                JsonObject schedParObj = new JsonObject();
                for (JobClass jc : classes) {
                    double par = queue.getSchedStrategyPar(jc);
                    if (par != 0.0 && !Double.isNaN(par)) {
                        schedParObj.addProperty(jc.getName(), par);
                        hasSchedPar = true;
                    }
                }
                if (hasSchedPar && !isDelay) {
                    nodeObj.add("schedParams", schedParObj);
                }

                // Buffer capacity
                int cap = (int) queue.getCap();
                if (cap < Integer.MAX_VALUE && cap > 0) {
                    nodeObj.addProperty("buffer", cap);
                }

                // Per-class buffer capacity
                JsonObject classCapObj = new JsonObject();
                for (JobClass jc : classes) {
                    int ccap = (int) queue.getClassCap(jc);
                    if (ccap < Integer.MAX_VALUE && ccap > 0 && ccap != cap) {
                        classCapObj.addProperty(jc.getName(), ccap);
                    }
                }
                if (classCapObj.size() > 0) {
                    nodeObj.add("classCap", classCapObj);
                }

                // Drop rules
                JsonObject dropRuleObj = new JsonObject();
                for (JobClass jc : classes) {
                    DropStrategy dr = queue.getDropRule(jc);
                    if (dr != null && dr != DropStrategy.WaitingQueue) {
                        dropRuleObj.addProperty(jc.getName(), dropStrategyToStr(dr));
                    }
                }
                if (dropRuleObj.size() > 0) {
                    nodeObj.add("dropRule", dropRuleObj);
                }

                // Load-dependent scaling
                Matrix lld = queue.getLimitedLoadDependence();
                if (lld != null && !lld.isEmpty()) {
                    JsonObject ldObj = new JsonObject();
                    ldObj.addProperty("type", "loadDependent");
                    JsonArray scalingArr = new JsonArray();
                    for (int k = 0; k < lld.getNumElements(); k++) {
                        scalingArr.add(lld.get(k));
                    }
                    ldObj.add("scaling", scalingArr);
                    nodeObj.add("loadDependence", ldObj);
                }
                // classDependence lattice materialization: see _kb/12-interfaces-and-docs.md
                SerializableFunction<Matrix, Matrix> lcd = queue.getLimitedClassDependence();
                if (lcd != null) {
                    int Kcd = model.getNumberOfClasses();
                    int[] maxc = classDependenceCutoffs(model, Kcd);
                    JsonObject cdObj = new JsonObject();
                    cdObj.addProperty("type", "classDependent");
                    JsonArray cutArr = new JsonArray();
                    for (int r = 0; r < Kcd; r++) {
                        cutArr.add(maxc[r]);
                    }
                    cdObj.add("cutoffs", cutArr);
                    JsonObject tblObj = new JsonObject();
                    int[] n = new int[Kcd];
                    int total = 1;
                    for (int r = 0; r < Kcd; r++) {
                        total *= (maxc[r] + 1);
                    }
                    for (int li = 0; li < total; li++) {
                        int rem = li;
                        for (int r = 0; r < Kcd; r++) {
                            n[r] = rem % (maxc[r] + 1);
                            rem /= (maxc[r] + 1);
                        }
                        Matrix nv = new Matrix(1, Kcd);
                        for (int r = 0; r < Kcd; r++) {
                            nv.set(0, r, n[r]);
                        }
                        Matrix bv = lcd.apply(nv);
                        JsonArray va = new JsonArray();
                        for (int r = 0; r < Kcd; r++) {
                            // A scalar-valued handle is broadcast to K entries so
                            // the reader is uniform.
                            double x = (bv.getNumElements() == 1) ? bv.get(0) : bv.get(r);
                            va.add(Double.isFinite(x) ? x : 0.0);
                        }
                        StringBuilder sb = new StringBuilder();
                        for (int r = 0; r < Kcd; r++) {
                            if (r > 0) sb.append(",");
                            sb.append(n[r]);
                        }
                        tblObj.add(sb.toString(), va);
                    }
                    cdObj.add("scaling", tblObj);
                    // Declared peak rate scaling per class (Util = T*S/peak),
                    // broadcast a 1x1 declaration to Kcd entries.
                    Matrix lcdPeak = queue.getLimitedClassDependencePeak();
                    JsonArray peakArr = new JsonArray();
                    for (int r = 0; r < Kcd; r++) {
                        double pv = 1.0;
                        if (lcdPeak != null && !lcdPeak.isEmpty()) {
                            pv = (lcdPeak.length() == 1) ? lcdPeak.get(0) : lcdPeak.get(r);
                        }
                        peakArr.add(pv);
                    }
                    cdObj.add("peak", peakArr);
                    nodeObj.add("classDependence", cdObj);
                }
                // jointDependence lattice materialization: non-product-form eta_i,
                // twin of classDependence above (see _kb/12-interfaces-and-docs.md).
                SerializableFunction<Matrix, Matrix> ljd = queue.getLimitedJointDependence();
                if (ljd != null) {
                    int Kjd = model.getNumberOfClasses();
                    int[] maxc = classDependenceCutoffs(model, Kjd);
                    JsonObject jdObj = new JsonObject();
                    jdObj.addProperty("type", "jointDependent");
                    JsonArray cutArr = new JsonArray();
                    for (int r = 0; r < Kjd; r++) {
                        cutArr.add(maxc[r]);
                    }
                    jdObj.add("cutoffs", cutArr);
                    JsonObject tblObj = new JsonObject();
                    int[] n = new int[Kjd];
                    int total = 1;
                    for (int r = 0; r < Kjd; r++) {
                        total *= (maxc[r] + 1);
                    }
                    for (int li = 0; li < total; li++) {
                        int rem = li;
                        for (int r = 0; r < Kjd; r++) {
                            n[r] = rem % (maxc[r] + 1);
                            rem /= (maxc[r] + 1);
                        }
                        Matrix nv = new Matrix(1, Kjd);
                        for (int r = 0; r < Kjd; r++) {
                            nv.set(0, r, n[r]);
                        }
                        Matrix bv = ljd.apply(nv);
                        JsonArray va = new JsonArray();
                        for (int r = 0; r < Kjd; r++) {
                            double x = (bv.getNumElements() == 1) ? bv.get(0) : bv.get(r);
                            va.add(Double.isFinite(x) ? x : 0.0);
                        }
                        StringBuilder sb = new StringBuilder();
                        for (int r = 0; r < Kjd; r++) {
                            if (r > 0) sb.append(",");
                            sb.append(n[r]);
                        }
                        tblObj.add(sb.toString(), va);
                    }
                    jdObj.add("scaling", tblObj);
                    Matrix ljdPeak = queue.getLimitedJointDependencePeak();
                    JsonArray peakArr = new JsonArray();
                    for (int r = 0; r < Kjd; r++) {
                        double pv = 1.0;
                        if (ljdPeak != null && !ljdPeak.isEmpty()) {
                            pv = (ljdPeak.length() == 1) ? ljdPeak.get(0) : ljdPeak.get(r);
                        }
                        peakArr.add(pv);
                    }
                    jdObj.add("peak", peakArr);
                    nodeObj.add("jointDependence", jdObj);
                }

                // Heterogeneous server types
                if (queue.isHeterogeneous()) {
                    JsonArray stArr = new JsonArray();
                    for (ServerType st : queue.getServerTypes()) {
                        JsonObject stj = new JsonObject();
                        stj.addProperty("name", st.getName());
                        stj.addProperty("count", st.getNumOfServers());
                        // Compatible classes
                        JsonArray ccArr = new JsonArray();
                        for (JobClass cc : st.getCompatibleClasses()) {
                            ccArr.add(cc.getName());
                        }
                        if (ccArr.size() > 0) stj.add("compatibleClasses", ccArr);
                        // Per-class service distributions
                        JsonObject svcObj = new JsonObject();
                        for (JobClass jc : classes) {
                            Distribution svcDist = queue.getService(jc, st);
                            if (svcDist != null) {
                                svcObj.add(jc.getName(), serializeDistribution(svcDist));
                            }
                        }
                        if (svcObj.size() > 0) stj.add("service", svcObj);
                        stArr.add(stj);
                    }
                    nodeObj.add("serverTypes", stArr);
                    HeteroSchedPolicy hsp = queue.getHeteroSchedPolicy();
                    if (hsp != null && hsp != HeteroSchedPolicy.ORDER) {
                        nodeObj.addProperty("heteroSchedPolicy", hsp.toText());
                    }
                }

                // Balking
                JsonObject balkObj = new JsonObject();
                for (JobClass jc : classes) {
                    if (queue.hasBalking(jc)) {
                        BalkingStrategy bs = queue.getBalkingStrategy(jc);
                        List<BalkingThreshold> thresholds = queue.getBalkingThresholds(jc);
                        JsonObject bjc = new JsonObject();
                        bjc.addProperty("strategy", bs.name());
                        JsonArray thArr = new JsonArray();
                        for (BalkingThreshold th : thresholds) {
                            JsonObject tjson = new JsonObject();
                            tjson.addProperty("minJobs", th.getMinJobs());
                            tjson.addProperty("maxJobs", th.getMaxJobs() == Integer.MAX_VALUE ? -1 : th.getMaxJobs());
                            tjson.addProperty("probability", th.getProbability());
                            thArr.add(tjson);
                        }
                        bjc.add("thresholds", thArr);
                        balkObj.add(jc.getName(), bjc);
                    }
                }
                if (balkObj.size() > 0) {
                    nodeObj.add("balking", balkObj);
                }

                // Retrial
                JsonObject retObj = new JsonObject();
                for (JobClass jc : classes) {
                    if (queue.hasRetrial(jc)) {
                        Distribution delayDist = queue.getRetrialDelayDistribution(jc);
                        int maxAttempts = queue.getMaxRetrialAttempts(jc);
                        JsonObject rjc = new JsonObject();
                        rjc.add("delay", serializeDistribution(delayDist));
                        rjc.addProperty("maxAttempts", maxAttempts);
                        retObj.add(jc.getName(), rjc);
                    }
                }
                if (retObj.size() > 0) {
                    nodeObj.add("retrial", retObj);
                }

                // Patience
                JsonObject patObj = new JsonObject();
                for (JobClass jc : classes) {
                    if (queue.hasPatience(jc)) {
                        Distribution patDist = queue.getPatience(jc);
                        JsonObject pjc = new JsonObject();
                        pjc.add("distribution", serializeDistribution(patDist));
                        ImpatienceType impType = queue.getImpatienceType(jc);
                        if (impType != null) {
                            pjc.addProperty("impatienceType", ImpatienceType.toText(impType));
                        }
                        patObj.add(jc.getName(), pjc);
                    }
                }
                if (patObj.size() > 0) {
                    nodeObj.add("patience", patObj);
                }

                // Orbit impatience (abandonment from the retrial orbit). Declared
                // on Station, so it applies to a Delay as well.
                JsonObject orbitObj = new JsonObject();
                for (JobClass jc : classes) {
                    if (queue.hasOrbitImpatience(jc)) {
                        orbitObj.add(jc.getName(),
                                serializeDistribution(queue.getOrbitImpatience(jc)));
                    }
                }
                if (orbitObj.size() > 0) {
                    nodeObj.add("orbitImpatience", orbitObj);
                }

                // Batch rejection probability (retrial queues). Declared on
                // Station, so it applies to a Delay as well.
                JsonObject batchRejectObj = new JsonObject();
                for (JobClass jc : classes) {
                    double brp = queue.getBatchRejectProbability(jc);
                    if (brp > 0) {
                        batchRejectObj.addProperty(jc.getName(), brp);
                    }
                }
                if (batchRejectObj.size() > 0) {
                    nodeObj.add("batchRejectProb", batchRejectObj);
                }

                // Job parallelism: servers seized at once by a job, per class
                JsonObject parallelismObj = new JsonObject();
                for (JobClass jc : classes) {
                    int npar = queue.getServerParallelism(jc);
                    if (npar > 1) {
                        parallelismObj.addProperty(jc.getName(), npar);
                    }
                }
                if (parallelismObj.size() > 0) {
                    nodeObj.add("serverParallelism", parallelismObj);
                }

                // immediateFeedback: see _kb/09-ldes-and-cache.md ("More model.json writer/reader notes")
                JsonObject immFeedObj = new JsonObject();
                if (queue.isImmediateFeedbackAll()) {
                    for (JobClass jc : classes) {
                        immFeedObj.addProperty(jc.getName(), true);
                    }
                } else {
                    Set<Integer> immFeedIdx = queue.getImmediateFeedbackClasses();
                    if (immFeedIdx != null) {
                        for (JobClass jc : classes) {
                            if (immFeedIdx.contains(jc.getIndex())) {
                                immFeedObj.addProperty(jc.getName(), true);
                            }
                        }
                    }
                }
                if (immFeedObj.size() > 0) {
                    nodeObj.add("immediateFeedback", immFeedObj);
                }

                // Setup / delay-off (server vacation). Written per class so the
                // LDES engine can rebuild isDelayOffEnabled() after a round trip.
                JsonObject setupObj = new JsonObject();
                JsonObject delayOffObj = new JsonObject();
                for (JobClass jc : classes) {
                    Distribution su = queue.getSetupTime(jc);
                    if (su != null && !(su instanceof Disabled)) {
                        setupObj.add(jc.getName(), serializeDistribution(su));
                    }
                    Distribution doff = queue.getDelayOffTime(jc);
                    if (doff != null && !(doff instanceof Disabled)) {
                        delayOffObj.add(jc.getName(), serializeDistribution(doff));
                    }
                }
                if (setupObj.size() > 0) {
                    nodeObj.add("setupTime", setupObj);
                }
                if (delayOffObj.size() > 0) {
                    nodeObj.add("delayOffTime", delayOffObj);
                }

                // Server breakdown/repair. The degraded down-server service is
                // written per class, which flattens the class-independent form
                // setBreakdown also accepts: both rebuild the same
                // sn.downServiceRates row.
                if (queue.hasBreakdown()) {
                    JsonObject bdObj = new JsonObject();
                    bdObj.add("failure", serializeDistribution(queue.getBreakdownFailure()));
                    bdObj.add("repair", serializeDistribution(queue.getBreakdownRepair()));
                    JsonObject downSvcObj = new JsonObject();
                    for (JobClass jc : classes) {
                        Distribution ds = queue.getDownService(jc);
                        if (ds != null && !(ds instanceof Disabled)) {
                            downSvcObj.add(jc.getName(), serializeDistribution(ds));
                        }
                    }
                    if (downSvcObj.size() > 0) {
                        bdObj.add("downService", downSvcObj);
                    }
                    nodeObj.add("breakdown", bdObj);
                }

                // pollingType written by name (Python auto() differs numerically): see _kb/09-ldes-and-cache.md
                if (queue.getSchedStrategy() == SchedStrategy.POLLING
                        && queue.getServer() instanceof PollingServer) {
                    PollingServer ps = (PollingServer) queue.getServer();
                    PollingType pt = ps.getPollingType();
                    if (pt != null) {
                        nodeObj.addProperty("pollingType", pt.name());
                        if (pt == PollingType.KLIMITED) {
                            nodeObj.addProperty("pollingPar", ps.getPollingK());
                        }
                    }
                }
                // POLLING switchover indexed by departing class alone: see _kb/09-ldes-and-cache.md
                JsonArray soArr = new JsonArray();
                if (queue.getSchedStrategy() == SchedStrategy.POLLING) {
                    for (JobClass fromJc : classes) {
                        Distribution so = queue.getSwitchover(fromJc);
                        if (so != null && !(so instanceof Disabled)) {
                            JsonObject soj = new JsonObject();
                            soj.addProperty("from", fromJc.getName());
                            soj.add("distribution", serializeDistribution(so));
                            soArr.add(soj);
                        }
                    }
                } else {
                    for (JobClass fromJc : classes) {
                        for (JobClass toJc : classes) {
                            Distribution so = queue.getSwitchoverTime(fromJc, toJc);
                            if (so != null && !(so instanceof Disabled)) {
                                JsonObject soj = new JsonObject();
                                soj.addProperty("from", fromJc.getName());
                                soj.addProperty("to", toJc.getName());
                                soj.add("distribution", serializeDistribution(so));
                                soArr.add(soj);
                            }
                        }
                    }
                }
                if (soArr.size() > 0) {
                    nodeObj.add("switchoverTimes", soArr);
                }

            } else if (node instanceof Fork) {
                nodeObj.addProperty("type", "Fork");
                Fork f = (Fork) node;
                Forker fkr = (Forker) f.getOutput();
                int tpl = (int) fkr.tasksPerLink;
                if (tpl > 1) {
                    nodeObj.addProperty("tasksPerLink", tpl);
                }

            } else if (node instanceof Join) {
                nodeObj.addProperty("type", "Join");
                Join join = (Join) node;
                if (join.joinOf != null) {
                    nodeObj.addProperty("forkNode", join.joinOf.getName());
                }
                // Join strategy and quorum
                Joiner joiner = (Joiner) join.getInput();
                if (joiner != null && joiner.joinStrategy != null) {
                    for (Map.Entry<JobClass, JoinStrategy> jsEntry : joiner.joinStrategy.entrySet()) {
                        if (jsEntry.getValue() != JoinStrategy.STD) {
                            // Save as PARTIAL for MATLAB compatibility (Quorum→PARTIAL)
                            String jsName = jsEntry.getValue() == JoinStrategy.Quorum ? "PARTIAL" : jsEntry.getValue().name();
                            nodeObj.addProperty("joinStrategy", jsName);
                            break;
                        }
                    }
                }
                if (joiner != null && joiner.joinRequired != null) {
                    for (Map.Entry<JobClass, Double> jqEntry : joiner.joinRequired.entrySet()) {
                        if (jqEntry.getValue() > 0) {
                            nodeObj.addProperty("joinQuorum", jqEntry.getValue().intValue());
                            break;
                        }
                    }
                }

            } else if (node instanceof Router) {
                nodeObj.addProperty("type", "Router");

            } else if (node instanceof ClassSwitch) {
                ClassSwitch cs = (ClassSwitch) node;
                nodeObj.addProperty("type", "ClassSwitch");
                // Serialize the class-switch matrix as nested dict (compatible with Python/MATLAB)
                int K = classes.size();
                ClassSwitcher switcher = (ClassSwitcher) cs.getServer();
                JsonObject csDict = new JsonObject();
                for (int r = 0; r < K; r++) {
                    JsonObject rowObj = new JsonObject();
                    for (int s = 0; s < K; s++) {
                        double val = switcher.applyCsFun(r, s);
                        if (val != 0.0) {
                            rowObj.addProperty(classes.get(s).getName(), val);
                        }
                    }
                    if (rowObj.size() > 0) {
                        csDict.add(classes.get(r).getName(), rowObj);
                    }
                }
                if (csDict.size() > 0) {
                    nodeObj.add("classSwitchMatrix", csDict);
                }

            } else if (node instanceof Cache) {
                Cache cache = (Cache) node;
                nodeObj.addProperty("type", "Cache");
                nodeObj.addProperty("numItems", cache.getNumberOfItems());
                // see _kb/09-ldes-and-cache.md (CLIMB is rewritten at solve time, never on the wire)
                nodeObj.add("itemLevelCap", matrixToJsonArray(cache.getItemLevelCap()));
                nodeObj.addProperty("replacementStrategy", cache.getReplacementStrategy().toString());
                nodeObj.addProperty("admissionProb", cache.getAdmissionProb());
                // per-item storage costs and per-list cost caps (ton21cache Sec. IX)
                if (cache.getItemSizes() != null && !cache.getItemSizes().isEmpty()) {
                    nodeObj.add("itemSizes", matrixToJsonArray(cache.getItemSizes()));
                }
                if (cache.getCostCaps() != null && !cache.getCostCaps().isEmpty()) {
                    if (cache.isCostCapGlobal()) {
                        nodeObj.addProperty("costCaps", cache.getCostCaps().get(0));
                    } else {
                        nodeObj.add("costCaps", matrixToJsonArray(cache.getCostCaps()));
                    }
                }

                // Hit/miss class mapping
                Matrix hitClass = cache.getHitClass();
                Matrix missClass = cache.getMissClass();
                if (hitClass != null && hitClass.getNumCols() > 0) {
                    JsonObject hitMap = new JsonObject();
                    for (int r = 0; r < hitClass.getNumCols(); r++) {
                        int hc = (int) hitClass.get(r);
                        if (hc >= 0 && hc < classes.size()) {
                            hitMap.addProperty(classes.get(r).getName(), classes.get(hc).getName());
                        }
                    }
                    nodeObj.add("hitClass", hitMap);
                }
                if (missClass != null && missClass.getNumCols() > 0) {
                    JsonObject missMap = new JsonObject();
                    for (int r = 0; r < missClass.getNumCols(); r++) {
                        int mc = (int) missClass.get(r);
                        if (mc >= 0 && mc < classes.size()) {
                            missMap.addProperty(classes.get(r).getName(), classes.get(mc).getName());
                        }
                    }
                    nodeObj.add("missClass", missMap);
                }

                // Popularity distributions
                JsonObject popObj = new JsonObject();
                for (JobClass jc : classes) {
                    int jcIdx = jc.getIndex() - 1;
                    Distribution pop = cache.popularityGet(0, jcIdx);
                    if (pop != null) {
                        popObj.add(jc.getName(), serializeDistribution(pop));
                    }
                }
                if (popObj.size() > 0) {
                    nodeObj.add("popularity", popObj);
                }

                // see _kb/09-ldes-and-cache.md for the access-cost (accessGraph/accessProb) rationale
                if (cache.getGraph() != null && cache.getGraph().length > 0) {
                    JsonArray gArr = new JsonArray();
                    for (Matrix g : cache.getGraph()) {
                        gArr.add(g != null ? matrixToJson2D(g) : new JsonArray());
                    }
                    nodeObj.add("accessGraph", gArr);
                } else if (cache.accessProb != null && cache.accessProb.length > 0) {
                    JsonArray apArr = new JsonArray();
                    for (Matrix[] classRow : cache.accessProb) {
                        JsonArray rowArr = new JsonArray();
                        if (classRow != null) {
                            for (Matrix m : classRow) {
                                rowArr.add(m != null ? matrixToJson2D(m) : new JsonArray());
                            }
                        }
                        apArr.add(rowArr);
                    }
                    nodeObj.add("accessProb", apArr);
                }

                // Initial cache state [class counts | contents | retrieval bitmap]
                Matrix cacheState = cache.getState();
                if (cacheState != null && !cacheState.isEmpty()) {
                    JsonArray stArr = new JsonArray();
                    for (int k = 0; k < cacheState.getNumElements(); k++) {
                        stArr.add(cacheState.get(k));
                    }
                    nodeObj.add("initialState", stArr);
                }

                // Retrieval system (delayed hits): capacity, per-arrival-class
                // retrieval queues and per-item retrieval classes.
                if (cache.getRetrievalSystemCapacity() > 0) {
                    JsonObject rsObj = new JsonObject();
                    rsObj.addProperty("capacity", cache.getRetrievalSystemCapacity());
                    JsonObject byClass = new JsonObject();
                    Matrix rc = cache.getRetrievalClasses();   // [items x classes]
                    for (Map.Entry<Integer, List<Integer>> qe : cache.getRetrievalSystemQueueIndices().entrySet()) {
                        int jobinIdx0 = qe.getKey();
                        if (jobinIdx0 < 0 || jobinIdx0 >= classes.size()) continue;
                        JsonObject entry = new JsonObject();
                        JsonArray qArr = new JsonArray();
                        for (Integer qi : qe.getValue()) {
                            if (qi >= 0 && qi < nodes.size()) qArr.add(nodes.get(qi).getName());
                        }
                        entry.add("queues", qArr);
                        JsonObject itemsObj = new JsonObject();
                        if (rc != null) {
                            for (int it = 0; it < rc.getNumRows(); it++) {
                                int rci = (int) rc.get(it, jobinIdx0);
                                if (rci >= 0 && rci < classes.size()) {
                                    itemsObj.addProperty(String.valueOf(it), classes.get(rci).getName());
                                }
                            }
                        }
                        entry.add("items", itemsObj);
                        byClass.add(classes.get(jobinIdx0).getName(), entry);
                    }
                    rsObj.add("byClass", byClass);
                    nodeObj.add("retrievalSystem", rsObj);
                }

            } else if (node instanceof Place) {
                nodeObj.addProperty("type", "Place");
                Place place = (Place) node;
                // see _kb/09-ldes-and-cache.md for the queueing-place (QPN) wire-format rationale
                if (place.isQueueing()) {
                    nodeObj.addProperty("scheduling", place.getSchedStrategy().toString());
                    int pns = place.getNumberOfServers();
                    if (pns != 1 && pns != Integer.MAX_VALUE) {
                        nodeObj.addProperty("servers", pns);
                    }
                    JsonObject services = new JsonObject();
                    JsonObject depDisc = new JsonObject();
                    for (JobClass jc : model.getClasses()) {
                        Distribution svc = place.getService(jc);
                        if (svc != null) {
                            services.add(jc.getName(), serializeDistribution(svc));
                            depDisc.addProperty(jc.getName(),
                                    place.getDepartureDiscipline(jc).toString());
                        }
                    }
                    if (services.size() > 0) {
                        nodeObj.add("service", services);
                        nodeObj.add("departureDiscipline", depDisc);
                    }
                }
                Matrix st = place.getState();
                if (nodeStateTravels && st != null && !st.isEmpty()) {
                    JsonArray stArr = new JsonArray();
                    for (int k = 0; k < st.getNumElements(); k++) {
                        stArr.add(st.get(k));
                    }
                    nodeObj.add("initialState", stArr);
                }

            } else if (node instanceof Transition) {
                Transition trans = (Transition) node;
                nodeObj.addProperty("type", "Transition");
                List<Mode> modesList = trans.getModes();
                if (!modesList.isEmpty()) {
                    List<Node> allNodes = model.getNodes();
                    JsonArray modesArr = new JsonArray();
                    for (int mi = 0; mi < modesList.size(); mi++) {
                        Mode mode = modesList.get(mi);
                        JsonObject modeObj = new JsonObject();
                        modeObj.addProperty("name", mode.getName());
                        // Distribution
                        Distribution dist = trans.getFiringDistribution(mode);
                        if (dist != null && !(dist instanceof Disabled)) {
                            modeObj.add("distribution", serializeDistribution(dist));
                        }
                        // Timing strategy
                        TimingStrategy ts = trans.timingStrategies.get(mode);
                        if (ts != null) {
                            modeObj.addProperty("timingStrategy",
                                ts == TimingStrategy.IMMEDIATE ? "IMMEDIATE" : "TIMED");
                        }
                        // Number of servers
                        int numSrv = trans.getNumberOfModeServers(mode);
                        if (numSrv != 1) {
                            if (numSrv == Integer.MAX_VALUE) {
                                modeObj.addProperty("numServers", "Infinity");
                            } else {
                                modeObj.addProperty("numServers", numSrv);
                            }
                        }
                        // Firing priority
                        if (trans.firingPriorities.getNumElements() > mi) {
                            double fp = trans.firingPriorities.get(mi);
                            // Omit only when it equals the builder default of 1: an explicit 0 is a
                            // legal JMT firing priority and has to survive the round trip (BUG-90).
                            if (fp != 1) {
                                modeObj.addProperty("firingPriority", fp);
                            }
                        }
                        // Firing weight
                        if (trans.firingWeights.getNumElements() > mi) {
                            double fw = trans.firingWeights.get(mi);
                            if (fw != 1.0) {
                                modeObj.addProperty("firingWeight", fw);
                            }
                        }
                        // Marking-dependent firing-rate multiplier g_mode(marking),
                        // materialized over the enabling (place,class) box lattice.
                        // Timed modes only; a null handle is omitted (unit multiplier).
                        SerializableFunction<Matrix, Double> gmod = trans.getFiringRateDependence(mode);
                        Matrix ecMatMod = trans.enablingConditions.get(mode);
                        if (gmod != null && ts != TimingStrategy.IMMEDIATE && ecMatMod != null) {
                            List<int[]> slotIdx = new ArrayList<>();
                            List<int[]> slotCap = new ArrayList<>();
                            JsonArray slotArr = new JsonArray();
                            JsonArray cutArr = new JsonArray();
                            for (int ni = 0; ni < ecMatMod.getNumRows(); ni++) {
                                for (int ci = 0; ci < ecMatMod.getNumCols(); ci++) {
                                    if (ecMatMod.get(ni, ci) > 0) {
                                        double pcap = allNodes.get(ni).getCap();
                                        int cap = (Double.isFinite(pcap) && pcap < Integer.MAX_VALUE) ? (int) Math.round(pcap) : 10;
                                        slotIdx.add(new int[]{ni, ci});
                                        slotCap.add(new int[]{cap});
                                        JsonObject sm = new JsonObject();
                                        sm.addProperty("node", allNodes.get(ni).getName());
                                        sm.addProperty("class", classes.get(ci).getName());
                                        slotArr.add(sm);
                                        cutArr.add(cap);
                                    }
                                }
                            }
                            if (!slotIdx.isEmpty()) {
                                int P = slotIdx.size();
                                int total = 1;
                                for (int s = 0; s < P; s++) total *= (slotCap.get(s)[0] + 1);
                                int nnodesAll = allNodes.size();
                                int nclassesAll = classes.size();
                                JsonObject tblObj = new JsonObject();
                                for (int li = 0; li < total; li++) {
                                    int rem = li;
                                    int[] c = new int[P];
                                    for (int s = 0; s < P; s++) {
                                        c[s] = rem % (slotCap.get(s)[0] + 1);
                                        rem /= (slotCap.get(s)[0] + 1);
                                    }
                                    Matrix mm = new Matrix(nnodesAll, nclassesAll);
                                    mm.zero();
                                    for (int s = 0; s < P; s++) mm.set(slotIdx.get(s)[0], slotIdx.get(s)[1], c[s]);
                                    double v = gmod.apply(mm);
                                    StringBuilder sb = new StringBuilder();
                                    for (int s = 0; s < P; s++) {
                                        if (s > 0) sb.append(",");
                                        sb.append(c[s]);
                                    }
                                    tblObj.addProperty(sb.toString(), Double.isFinite(v) ? v : 0.0);
                                }
                                JsonObject frmObj = new JsonObject();
                                frmObj.add("slots", slotArr);
                                frmObj.add("cutoffs", cutArr);
                                frmObj.add("scaling", tblObj);
                                modeObj.add("firingRateDependence", frmObj);
                            }
                        }
                        // Enabling conditions
                        Matrix ecMat = trans.enablingConditions.get(mode);
                        if (ecMat != null) {
                            JsonArray ecArr = new JsonArray();
                            for (int ni = 0; ni < ecMat.getNumRows(); ni++) {
                                for (int ci = 0; ci < ecMat.getNumCols(); ci++) {
                                    if (ecMat.get(ni, ci) > 0) {
                                        JsonObject ec = new JsonObject();
                                        ec.addProperty("node", allNodes.get(ni).getName());
                                        ec.addProperty("class", classes.get(ci).getName());
                                        ec.addProperty("count", ecMat.get(ni, ci));
                                        ecArr.add(ec);
                                    }
                                }
                            }
                            if (ecArr.size() > 0) {
                                modeObj.add("enablingConditions", ecArr);
                            }
                        }
                        // Inhibiting conditions
                        Matrix icMat = trans.inhibitingConditions.get(mode);
                        if (icMat != null) {
                            JsonArray icArr = new JsonArray();
                            for (int ni = 0; ni < icMat.getNumRows(); ni++) {
                                for (int ci = 0; ci < icMat.getNumCols(); ci++) {
                                    double val = icMat.get(ni, ci);
                                    if (Double.isFinite(val) && val > 0) {
                                        JsonObject ic = new JsonObject();
                                        ic.addProperty("node", allNodes.get(ni).getName());
                                        ic.addProperty("class", classes.get(ci).getName());
                                        ic.addProperty("count", val);
                                        icArr.add(ic);
                                    }
                                }
                            }
                            if (icArr.size() > 0) {
                                modeObj.add("inhibitingConditions", icArr);
                            }
                        }
                        // Firing outcomes
                        Matrix foMat = trans.firingOutcomes.get(mode);
                        if (foMat != null) {
                            JsonArray foArr = new JsonArray();
                            for (int ni = 0; ni < foMat.getNumRows(); ni++) {
                                for (int ci = 0; ci < foMat.getNumCols(); ci++) {
                                    if (foMat.get(ni, ci) != 0) {
                                        JsonObject fo = new JsonObject();
                                        fo.addProperty("node", allNodes.get(ni).getName());
                                        fo.addProperty("class", classes.get(ci).getName());
                                        fo.addProperty("count", foMat.get(ni, ci));
                                        foArr.add(fo);
                                    }
                                }
                            }
                            if (foArr.size() > 0) {
                                modeObj.add("firingOutcomes", foArr);
                            }
                        }
                        modesArr.add(modeObj);
                    }
                    nodeObj.add("modes", modesArr);
                }

            } else {
                nodeObj.addProperty("type", node.getClass().getSimpleName());
            }

            // see _kb/09-ldes-and-cache.md (statePrior is emitted only as a pair with stateSpace)
            if (node instanceof StatefulNode) {
                StatefulNode sfNode = (StatefulNode) node;
                // The initial state itself, for every stateful node that has one:
                // the reader decides the model is initialized only when EVERY
                // stateful node carries a state, so naming just the ones the
                // caller moved makes the document read as uninitialized and the
                // state space and prior below are then rebuilt from scratch. The
                // Place and Cache branches above write the same key for the nodes
                // whose state is their whole model, so this fills it in only where
                // it is still absent.
                Matrix nodeState = sfNode.getState();
                if (nodeStateTravels && !nodeObj.has("initialState") && nodeState != null && !nodeState.isEmpty()) {
                    nodeObj.add("initialState", matrixToJsonRowVector(nodeState));
                }
                Matrix prior = sfNode.getStatePrior();
                boolean trivialPrior = prior == null || prior.isEmpty()
                        || (prior.getNumRows() == 1 && Math.abs(prior.get(0) - 1.0) < 1e-12);
                if (!trivialPrior && nodeStateTravels) {
                    Matrix space = sfNode.getStateSpace();
                    if (space == null || space.getNumRows() != prior.getNumRows()) {
                        line_warning(mfilename(new Object() {
                        }), "Node %s carries a state prior over %d states but a state space of %d "
                                + "rows; the prior is not saved.", node.getName(),
                                prior.getNumRows(), space == null ? 0 : space.getNumRows());
                    } else {
                        nodeObj.add("stateSpace", matrixToJson2D(space));
                        nodeObj.add("statePrior", matrixToJsonRowVector(prior));
                    }
                }
            }

            nodesArr.add(nodeObj);
        }
        return nodesArr;
    }

    private static JsonArray serializeNetworkClasses(Network model) {
        JsonArray classesArr = new JsonArray();
        List<JobClass> classes = model.getClasses();

        for (JobClass jc : classes) {
            JsonObject classObj = new JsonObject();
            classObj.addProperty("name", jc.getName());

            if (jc instanceof OpenSignal) {
                OpenSignal os = (OpenSignal) jc;
                classObj.addProperty("type", "Signal");
                classObj.addProperty("openOrClosed", "Open");
                classObj.addProperty("signalType", SignalType.toText(os.getSignalType()));
                if (os.getTargetJobClass() != null) {
                    classObj.addProperty("targetClass", os.getTargetJobClass().getName());
                }
                if (os.getRemovalDistribution() != null) {
                    classObj.add("removalDistribution", serializeDistribution(os.getRemovalDistribution()));
                }
                if (os.getRemovalPolicy() != null && os.getRemovalPolicy() != RemovalPolicy.RANDOM) {
                    classObj.addProperty("removalPolicy", RemovalPolicy.toText(os.getRemovalPolicy()));
                }
            } else if (jc instanceof ClosedSignal) {
                ClosedSignal cs = (ClosedSignal) jc;
                classObj.addProperty("type", "Signal");
                classObj.addProperty("openOrClosed", "Closed");
                classObj.addProperty("signalType", SignalType.toText(cs.getSignalType()));
                if (cs.getReferenceStation() != null) {
                    classObj.addProperty("refNode", cs.getReferenceStation().getName());
                }
                if (cs.getTargetJobClass() != null) {
                    classObj.addProperty("targetClass", cs.getTargetJobClass().getName());
                }
                if (cs.getRemovalDistribution() != null) {
                    classObj.add("removalDistribution", serializeDistribution(cs.getRemovalDistribution()));
                }
                if (cs.getRemovalPolicy() != null && cs.getRemovalPolicy() != RemovalPolicy.RANDOM) {
                    classObj.addProperty("removalPolicy", RemovalPolicy.toText(cs.getRemovalPolicy()));
                }
            } else if (jc instanceof OpenClass) {
                classObj.addProperty("type", "Open");
                // see _kb/09-ldes-and-cache.md for the OpenClass refNode default-omission rationale
                Station openRef = jc.getReferenceStation();
                int srcIdx = model.getIndexSourceNode();
                Node defaultSrc = (srcIdx >= 0) ? model.getNodes().get(srcIdx) : null;
                if (openRef != null && openRef != defaultSrc) {
                    classObj.addProperty("refNode", openRef.getName());
                }
            } else if (jc instanceof SelfLoopingClass) {
                // see _kb/09-ldes-and-cache.md (SelfLoopingClass must be tested before ClosedClass)
                SelfLoopingClass slc = (SelfLoopingClass) jc;
                classObj.addProperty("type", "SelfLooping");
                classObj.addProperty("population", slc.getPopulation());
                if (slc.getReferenceStation() != null) {
                    classObj.addProperty("refNode", slc.getReferenceStation().getName());
                }
            } else if (jc instanceof ClosedClass) {
                ClosedClass cc = (ClosedClass) jc;
                classObj.addProperty("type", "Closed");
                classObj.addProperty("population", cc.getPopulation());
                if (cc.getReferenceStation() != null) {
                    classObj.addProperty("refNode", cc.getReferenceStation().getName());
                }
            }

            int priority = jc.getPriority();
            if (priority != 0) {
                classObj.addProperty("priority", priority);
            }

            double deadline = jc.getDeadline();
            if (Double.isFinite(deadline)) {
                classObj.addProperty("deadline", deadline);
            }

            if (jc.isReferenceClass()) {
                classObj.addProperty("isReferenceClass", true);
            }

            // Class-scoped patience, distinct from the node-scoped "patience"
            // that overrides it.
            if (jc.hasPatience()) {
                classObj.add("patience", serializeDistribution(jc.getPatience()));
                if (jc.getImpatienceType() != null) {
                    classObj.addProperty("impatienceType",
                            ImpatienceType.toText(jc.getImpatienceType()));
                }
            }

            // Reply signal class (syncreply). Stored 1-based, as the setter
            // documents; without it a REPLY signal is inert after a round trip.
            int replyIdx = jc.getReplySignalClassIndex();
            if (replyIdx >= 1 && replyIdx <= classes.size()) {
                classObj.addProperty("replySignalClass", classes.get(replyIdx - 1).getName());
            }

            // Spawn-on-completion binding (classspawn). Stored 1-based, as the
            // setter documents.
            int spawnIdx = jc.getSpawnClassIndex();
            if (spawnIdx >= 1 && spawnIdx <= classes.size()) {
                classObj.addProperty("spawnClass", classes.get(spawnIdx - 1).getName());
            }

            classesArr.add(classObj);
        }
        return classesArr;
    }

    private static JsonObject serializeNetworkRouting(Network model) {
        JsonObject routingObj = new JsonObject();
        routingObj.addProperty("type", "matrix");
        JsonObject matrixObj = new JsonObject();

        Map<JobClass, Map<JobClass, Matrix>> rtMap = model.getLinkedRoutingMatrix();
        // Use only original nodes (not auto-added ClassSwitch nodes) to match rtorig dimensions
        List<Node> allNodes = model.getNodes();
        List<Node> nodeList = new ArrayList<Node>();
        for (Node n : allNodes) {
            if (n instanceof ClassSwitch && ((ClassSwitch) n).autoAdded) {
                continue;
            }
            nodeList.add(n);
        }
        int M = nodeList.size();

        // see _kb/09-ldes-and-cache.md for the class-switch-node routing serialization rationale
        List<JobClass> rtClasses = model.getClasses();
        int Kcls = rtClasses.size();
        Set<Integer> csMatrixSwitchIndices = new HashSet<Integer>();
        for (int i = 0; i < M; i++) {
            Node n = nodeList.get(i);
            if (n instanceof ClassSwitch && !((ClassSwitch) n).autoAdded
                    && n.getServer() instanceof ClassSwitcher) {
                ClassSwitcher sw = (ClassSwitcher) n.getServer();
                boolean identity = true;
                for (int r = 0; r < Kcls && identity; r++) {
                    for (int s = 0; s < Kcls; s++) {
                        double expected = (r == s) ? 1.0 : 0.0;
                        if (Math.abs(sw.applyCsFun(r, s) - expected) > 1e-12) {
                            identity = false;
                            break;
                        }
                    }
                }
                if (!identity) {
                    csMatrixSwitchIndices.add(i);
                }
            }
        }
        // csSame[r][i][j] = sum over target class s of rt[r][s].get(i,j), for the
        // ClassSwitch rows whose switch is carried by the matrix.
        double[][][] csSame = null;
        if (rtMap != null && !csMatrixSwitchIndices.isEmpty()) {
            csSame = new double[Kcls][M][M];
            for (int r = 0; r < Kcls; r++) {
                Map<JobClass, Matrix> destMap = rtMap.get(rtClasses.get(r));
                if (destMap == null) {
                    continue;
                }
                for (int s = 0; s < Kcls; s++) {
                    Matrix rt = destMap.get(rtClasses.get(s));
                    if (rt == null) {
                        continue;
                    }
                    for (int ii : csMatrixSwitchIndices) {
                        for (int jj = 0; jj < M; jj++) {
                            csSame[r][ii][jj] += rt.get(ii, jj);
                        }
                    }
                }
            }
        }
        if (rtMap != null) {
            // Iterate classes in model order (not map entry order, which is a
            // HashMap and would make the JSON key order nondeterministic)
            for (int rIdx = 0; rIdx < Kcls; rIdx++) {
                JobClass cs = rtClasses.get(rIdx);
                Map<JobClass, Matrix> destMap = rtMap.get(cs);
                if (destMap == null) {
                    continue;
                }
                for (int sIdx = 0; sIdx < Kcls; sIdx++) {
                    JobClass cd = rtClasses.get(sIdx);
                    Matrix rt = destMap.get(cd);
                    if (rt == null) {
                        continue;
                    }
                    String key = cs.getName() + "," + cd.getName();
                    JsonObject fromTo = new JsonObject();
                    for (int i = 0; i < M; i++) {
                        if (csMatrixSwitchIndices.contains(i)) {
                            // Same-class topological routing only; the switch is
                            // reapplied on load from the classSwitchMatrix.
                            if (csSame == null || rIdx != sIdx || rIdx < 0) {
                                continue;
                            }
                            JsonObject dests = null;
                            for (int j = 0; j < M; j++) {
                                double val = csSame[rIdx][i][j];
                                if (val > 1e-14) {
                                    if (dests == null) {
                                        dests = new JsonObject();
                                    }
                                    dests.addProperty(nodeList.get(j).getName(), val);
                                }
                            }
                            if (dests != null) {
                                fromTo.add(nodeList.get(i).getName(), dests);
                            }
                            continue;
                        }
                        JsonObject dests = null;
                        for (int j = 0; j < M; j++) {
                            double val = rt.get(i, j);
                            if (val > 1e-14) {
                                if (dests == null) {
                                    dests = new JsonObject();
                                }
                                dests.addProperty(nodeList.get(j).getName(), val);
                            }
                        }
                        if (dests != null) {
                            fromTo.add(nodeList.get(i).getName(), dests);
                        }
                    }
                    if (fromTo.size() > 0) {
                        matrixObj.add(key, fromTo);
                    }
                }
            }
        }
        routingObj.add("matrix", matrixObj);

        // see _kb/09-ldes-and-cache.md for the non-default routing-strategy save rationale
        JsonObject routingStrategies = new JsonObject();
        JsonObject routingWeights = new JsonObject();
        JsonObject routingParams = new JsonObject();
        for (Node n : allNodes) {
            JsonObject nodeStrats = null;
            JsonObject nodeWeights = null;
            JsonObject nodeParams = null;
            for (JobClass jc : model.getClasses()) {
                RoutingStrategy rs = n.getRoutingStrategy(jc);
                // RAND is DECLARED not derived: dropping it wrote JMT Empirical where model wants Random. Implicit ClassSwitch not in `nodes`, so naming it dangles.
                boolean implicitCs = n instanceof ClassSwitch && ((ClassSwitch) n).autoAdded;
                if (rs != null && rs != RoutingStrategy.PROB && !implicitCs) {
                    if (nodeStrats == null) {
                        nodeStrats = new JsonObject();
                    }
                    nodeStrats.addProperty(jc.getName(), rs.toString());
                }
                // Power-of-k-choices parameters. Without them the reader
                // rebuilds the OutputStrategy default k=2, no memory.
                if (rs == RoutingStrategy.SQ) {
                    List<OutputStrategy> kosList = n.getOutput().getOutputStrategyByClass(jc);
                    if (kosList != null) {
                        for (OutputStrategy os : kosList) {
                            if (os.getRoutingStrategy() != RoutingStrategy.SQ) {
                                continue;
                            }
                            JsonObject kp = new JsonObject();
                            kp.addProperty("d", os.getSqD());
                            if (nodeParams == null) {
                                nodeParams = new JsonObject();
                            }
                            nodeParams.add(jc.getName(), kp);
                            break;
                        }
                    }
                }
                // Save WRROBIN weights
                if (rs == RoutingStrategy.WRROBIN) {
                    List<OutputStrategy> osList = n.getOutput().getOutputStrategyByClass(jc);
                    if (osList != null) {
                        JsonObject destWeights = null;
                        for (OutputStrategy os : osList) {
                            Node dest = os.getDestination();
                            if (dest != null) {
                                if (destWeights == null) {
                                    destWeights = new JsonObject();
                                }
                                destWeights.addProperty(dest.getName(), os.getProbability());
                            }
                        }
                        if (destWeights != null) {
                            if (nodeWeights == null) {
                                nodeWeights = new JsonObject();
                            }
                            nodeWeights.add(jc.getName(), destWeights);
                        }
                    }
                }
            }
            if (nodeStrats != null) {
                routingStrategies.add(n.getName(), nodeStrats);
            }
            if (nodeWeights != null) {
                routingWeights.add(n.getName(), nodeWeights);
            }
            if (nodeParams != null) {
                routingParams.add(n.getName(), nodeParams);
            }
        }
        if (routingStrategies.size() > 0) {
            routingObj.add("routingStrategies", routingStrategies);
        }
        if (routingWeights.size() > 0) {
            routingObj.add("routingWeights", routingWeights);
        }
        if (routingParams.size() > 0) {
            routingObj.add("routingParams", routingParams);
        }

        return routingObj;
    }

    // ========================================================================
    // LAYERED NETWORK SERIALIZATION HELPERS
    // ========================================================================

    private static JsonArray serializeHosts(LayeredNetwork model) {
        JsonArray hostsArr = new JsonArray();
        Map<Integer, Host> hosts = model.getHosts();
        for (Map.Entry<Integer, Host> entry : hosts.entrySet()) {
            Host host = entry.getValue();
            JsonObject hostObj = new JsonObject();
            hostObj.addProperty("name", host.getName());
            hostObj.addProperty("multiplicity", host.getMultiplicity());
            hostObj.addProperty("scheduling", host.getScheduling().toString());
            hostObj.addProperty("quantum", host.getQuantum());
            hostObj.addProperty("speedFactor", host.getSpeedFactor());
            int repl = host.getReplication();
            if (repl > 1) {
                hostObj.addProperty("replication", repl);
            }
            // Admission constraints: columns of a host constraint are its tasks
            List<String> hostCols = new ArrayList<String>();
            for (Task t : host.getTasks()) {
                hostCols.add(t.getName());
            }
            JsonArray hostRows = serializeLincon(host, hostCols);
            if (hostRows.size() > 0) {
                hostObj.add("admissionConstraints", hostRows);
            }
            hostsArr.add(hostObj);
        }
        return hostsArr;
    }

    private static JsonArray serializeTasks(LayeredNetwork model) {
        JsonArray tasksArr = new JsonArray();
        Map<Integer, Task> tasks = model.getTasks();
        for (Map.Entry<Integer, Task> entry : tasks.entrySet()) {
            Task task = entry.getValue();
            JsonObject taskObj = new JsonObject();
            taskObj.addProperty("name", task.getName());
            if (task.getParent() != null) {
                taskObj.addProperty("host", task.getParent().getName());
            }
            taskObj.addProperty("multiplicity", task.getMultiplicity());
            taskObj.addProperty("scheduling", task.getScheduling().toString());

            int repl = task.getReplication();
            if (repl > 1) {
                taskObj.addProperty("replication", repl);
            }

            // Both spellings are written: the object form carries the whole
            // distribution, the mean/SCV pair is kept so a reader that only knows
            // the older keys still loads the document. See the reader, which
            // accepts both and prefers the object.
            double thinkMean = task.getThinkTimeMean();
            if (thinkMean > 1e-8) {
                Distribution thinkDist = task.getThinkTime();
                if (thinkDist != null) {
                    taskObj.add("thinkTime", serializeDistribution(thinkDist));
                }
                taskObj.addProperty("thinkTimeMean", thinkMean);
                taskObj.addProperty("thinkTimeSCV", task.getThinkTimeSCV());
            }

            int priority = task.getPriority();
            if (priority != 0) {
                taskObj.addProperty("priority", priority);
            }

            // Fan-in, keyed by SOURCE task exactly as fanOut is keyed by dest.
            // MATLAB linemodel_save and the Python writer both emit the map, so
            // the pair form this used to write was unreadable by either.
            String fanInSrc = task.getFanInSource();
            if (fanInSrc != null && !fanInSrc.isEmpty()) {
                JsonObject fanInObj = new JsonObject();
                fanInObj.addProperty(fanInSrc, task.getFanInValue());
                taskObj.add("fanIn", fanInObj);
            }

            // Fan-out
            Map<String, Integer> fanOutMap = task.getFanOutMap();
            if (fanOutMap != null && !fanOutMap.isEmpty()) {
                JsonObject fanOutObj = new JsonObject();
                for (Map.Entry<String, Integer> fo : fanOutMap.entrySet()) {
                    fanOutObj.addProperty(fo.getKey(), fo.getValue());
                }
                taskObj.add("fanOut", fanOutObj);
            }

            // CacheTask / SetupTask properties
            if (task instanceof CacheTask) {
                CacheTask ct = (CacheTask) task;
                taskObj.addProperty("taskType", "CacheTask");
                taskObj.addProperty("totalItems", ct.getItems());
                JsonArray capArr = new JsonArray();
                for (int c : ct.getItemLevelCap()) {
                    capArr.add(c);
                }
                taskObj.add("cacheCapacity", capArr);
                taskObj.addProperty("replacementStrategy", ct.getReplacestrategy().toString());
            } else if (task instanceof SetupTask) {
                taskObj.addProperty("taskType", "SetupTask");
            }
            Distribution setupDist = task.getSetupTime();
            if (setupDist != null && !(setupDist instanceof Immediate)
                    && task.getSetupTimeMean() > GlobalConstants.FineTol) {
                taskObj.add("setupTime", serializeDistribution(setupDist));
                taskObj.addProperty("setupTimeMean", task.getSetupTimeMean());
                taskObj.addProperty("setupTimeSCV", task.getSetupTimeSCV());
            }
            Distribution delayOffDist = task.getDelayOffTime();
            if (delayOffDist != null && !(delayOffDist instanceof Immediate)
                    && task.getDelayOffTimeMean() > GlobalConstants.FineTol) {
                taskObj.add("delayOffTime", serializeDistribution(delayOffDist));
                taskObj.addProperty("delayOffTimeMean", task.getDelayOffTimeMean());
                taskObj.addProperty("delayOffTimeSCV", task.getDelayOffTimeSCV());
            }

            // Admission constraints: columns of a task constraint are its entries
            List<String> taskCols = new ArrayList<String>();
            for (Entry e : task.getEntries()) {
                taskCols.add(e.getName());
            }
            JsonArray taskRows = serializeLincon(task, taskCols);
            if (taskRows.size() > 0) {
                taskObj.add("admissionConstraints", taskRows);
            }

            tasksArr.add(taskObj);
        }
        return tasksArr;
    }

    /**
     * <p>Admission constraint rows of a Task or Host on the wire.</p>
     *
     * <p>Both declaration forms are normalised to the named form, so the wire is
     * order-independent: a positional setConstraint(A,b) matrix is resolved
     * against colNames (the element's entries, or its tasks) at write time.
     * colNames must be in the same declaration order the positional columns
     * assume.</p>
     *
     * @param elem     the Task or Host declaring the constraint
     * @param colNames operand names in positional-column order
     * @return the rows, empty when the element declares no constraint
     */
    private static JsonArray serializeLincon(LayeredNetworkElement elem, List<String> colNames) {
        JsonArray rows = new JsonArray();
        if (!elem.hasLinearConstraints()) {
            return rows;
        }
        Matrix A = elem.linConA;
        Matrix b = elem.linConB;
        if (A != null && b != null) {
            for (int r = 0; r < A.getNumRows(); r++) {
                JsonArray operands = new JsonArray();
                JsonArray coeffs = new JsonArray();
                for (int j = 0; j < A.getNumCols(); j++) {
                    if (A.get(r, j) == 0) {
                        continue;
                    }
                    if (j >= colNames.size()) {
                        throw new IllegalArgumentException("Admission constraint on " + elem.getName()
                                + " references column " + (j + 1) + " but the element has only "
                                + colNames.size() + " operands.");
                    }
                    operands.add(colNames.get(j));
                    coeffs.add(A.get(r, j));
                }
                if (operands.size() == 0) {
                    continue;
                }
                JsonObject row = new JsonObject();
                row.add("operands", operands);
                row.add("coeffs", coeffs);
                row.addProperty("cap", b.get(r, 0));
                rows.add(row);
            }
        }
        for (LayeredNetworkElement.LinConRow named : elem.linConRows) {
            JsonArray operands = new JsonArray();
            JsonArray coeffs = new JsonArray();
            for (int k = 0; k < named.names.size(); k++) {
                operands.add(named.names.get(k));
                coeffs.add(named.coeffs[k]);
            }
            JsonObject row = new JsonObject();
            row.add("operands", operands);
            row.add("coeffs", coeffs);
            row.addProperty("cap", named.cap);
            rows.add(row);
        }
        return rows;
    }

    /**
     * Replays admission constraint rows from the wire onto a Task or Host. Rows
     * name their operands, so no column order is assumed and the referenced
     * entries or tasks need not exist yet.
     *
     * @param elem the Task or Host to configure
     * @param rows the wire rows, may be null
     */
    private static void applyLincon(LayeredNetworkElement elem, JsonArray rows) {
        if (rows == null) {
            return;
        }
        for (int r = 0; r < rows.size(); r++) {
            JsonObject row = rows.get(r).getAsJsonObject();
            JsonArray operands = row.getAsJsonArray("operands");
            List<String> names = new ArrayList<String>();
            for (int k = 0; k < operands.size(); k++) {
                names.add(operands.get(k).getAsString());
            }
            double[] coeffs = null;
            if (row.has("coeffs")) {
                JsonArray cf = row.getAsJsonArray("coeffs");
                coeffs = new double[cf.size()];
                for (int k = 0; k < cf.size(); k++) {
                    coeffs[k] = cf.get(k).getAsDouble();
                }
            }
            elem.addConstraintByName(names, coeffs, row.get("cap").getAsDouble());
        }
    }

    private static JsonArray serializeEntries(LayeredNetwork model) {
        JsonArray entriesArr = new JsonArray();
        Map<Integer, Entry> entries = model.getEntries();
        for (Map.Entry<Integer, Entry> entry : entries.entrySet()) {
            Entry e = entry.getValue();
            JsonObject entryObj = new JsonObject();
            entryObj.addProperty("name", e.getName());
            if (e.getParent() != null) {
                entryObj.addProperty("task", e.getParent().getName());
            }
            if (e.getArrival() != null) {
                entryObj.add("arrival", serializeDistribution(e.getArrival()));
            }

            // Forwarding
            Map<Integer, String> fwDests = e.getForwardingDests();
            Matrix fwProbs = e.getForwardingProbs();
            if (fwDests != null && !fwDests.isEmpty()) {
                JsonArray fwArr = new JsonArray();
                for (Map.Entry<Integer, String> fwEntry : fwDests.entrySet()) {
                    JsonObject fwObj = new JsonObject();
                    fwObj.addProperty("dest", fwEntry.getValue());
                    int idx = fwEntry.getKey();
                    if (fwProbs != null && !fwProbs.isEmpty() && idx < fwProbs.getNumCols()) {
                        fwObj.addProperty("prob", fwProbs.get(0, idx));
                    } else {
                        fwObj.addProperty("prob", 1.0);
                    }
                    fwArr.add(fwObj);
                }
                entryObj.add("forwarding", fwArr);
            }

            // ItemEntry properties
            if (e instanceof ItemEntry) {
                ItemEntry ie = (ItemEntry) e;
                entryObj.addProperty("entryType", "ItemEntry");
                entryObj.addProperty("totalItems", ie.getCardinality());
                Distribution pop = ie.getPopularity();
                if (pop != null) {
                    entryObj.add("accessProb", serializeDistribution(pop));
                }
            }

            entriesArr.add(entryObj);
        }
        return entriesArr;
    }

    private static JsonArray serializeActivities(LayeredNetwork model) {
        JsonArray actsArr = new JsonArray();
        Map<Integer, Activity> activities = model.getActivities();
        Map<Integer, Entry> entries = model.getEntries();

        for (Map.Entry<Integer, Activity> entry : activities.entrySet()) {
            Activity act = entry.getValue();
            int actIdx = entry.getKey();
            JsonObject actObj = new JsonObject();
            actObj.addProperty("name", act.getName());
            if (act.getParent() != null) {
                actObj.addProperty("task", act.getParent().getName());
            }

            // Host demand
            Distribution hostDem = act.getHostDemand();
            if (hostDem != null) {
                actObj.add("hostDemand", serializeDistribution(hostDem));
            }

            // Bound-to entry
            String boundTo = act.getBoundToEntry();
            if (boundTo != null && !boundTo.isEmpty()) {
                actObj.addProperty("boundToEntry", boundTo);
            }

            // Replies-to entries
            for (Map.Entry<Integer, Entry> entryEntry : entries.entrySet()) {
                Entry e = entryEntry.getValue();
                Map<Integer, String> replyMap = e.getReplyActivity();
                if (replyMap != null) {
                    for (Map.Entry<Integer, String> replyEntry : replyMap.entrySet()) {
                        if (replyEntry.getValue().equals(act.getName())) {
                            actObj.addProperty("repliesTo", e.getName());
                        }
                    }
                }
            }

            // Think time
            double thinkMean = act.getThinkTimeMean();
            if (thinkMean > 1e-8) {
                Distribution thinkDist = act.getThinkTime();
                if (thinkDist != null) {
                    actObj.add("thinkTime", serializeDistribution(thinkDist));
                }
            }

            // Synchronous calls
            Map<Integer, String> syncDests = act.getSyncCallDests();
            Matrix syncMeans = act.getSyncCallMeans();
            if (syncDests != null && !syncDests.isEmpty()) {
                JsonArray syncArr = new JsonArray();
                for (Map.Entry<Integer, String> syncEntry : syncDests.entrySet()) {
                    JsonObject callObj = new JsonObject();
                    callObj.addProperty("dest", syncEntry.getValue());
                    int idx = syncEntry.getKey();
                    double mean = 1.0;
                    if (syncMeans != null && idx < syncMeans.length()) {
                        mean = syncMeans.get(idx);
                    }
                    callObj.addProperty("mean", mean);
                    syncArr.add(callObj);
                }
                actObj.add("synchCalls", syncArr);
            }

            // Asynchronous calls
            Map<Integer, String> asyncDests = act.getAsyncCallDests();
            Matrix asyncMeans = act.getAsyncCallMeans();
            if (asyncDests != null && !asyncDests.isEmpty()) {
                JsonArray asyncArr = new JsonArray();
                for (Map.Entry<Integer, String> asyncEntry : asyncDests.entrySet()) {
                    JsonObject callObj = new JsonObject();
                    callObj.addProperty("dest", asyncEntry.getValue());
                    int idx = asyncEntry.getKey();
                    double mean = 1.0;
                    if (asyncMeans != null && idx < asyncMeans.length()) {
                        mean = asyncMeans.get(idx);
                    }
                    callObj.addProperty("mean", mean);
                    asyncArr.add(callObj);
                }
                actObj.add("asynchCalls", asyncArr);
            }

            // Call order
            String callOrder = act.getCallOrder();
            if (callOrder != null && !"STOCHASTIC".equals(callOrder)) {
                actObj.addProperty("callOrder", callOrder);
            }

            actsArr.add(actObj);
        }
        return actsArr;
    }

    private static JsonArray serializePrecedences(LayeredNetwork model) {
        JsonArray precsArr = new JsonArray();
        Map<Integer, Task> tasks = model.getTasks();

        for (Map.Entry<Integer, Task> taskEntry : tasks.entrySet()) {
            Task task = taskEntry.getValue();
            List<ActivityPrecedence> precs = task.getPrecedences();
            if (precs == null) {
                continue;
            }
            for (ActivityPrecedence prec : precs) {
                JsonObject precObj = new JsonObject();
                precObj.addProperty("task", task.getName());

                JsonArray preActs = new JsonArray();
                for (String name : prec.getPreActs()) {
                    preActs.add(name);
                }
                precObj.add("preActs", preActs);

                JsonArray postActs = new JsonArray();
                for (String name : prec.getPostActs()) {
                    postActs.add(name);
                }
                precObj.add("postActs", postActs);

                precObj.addProperty("preType", prec.getPreType());
                precObj.addProperty("postType", prec.getPostType());

                if (prec.getPreParams() != null && !prec.getPreParams().isEmpty()) {
                    precObj.add("preParams", matrixToJsonArray(prec.getPreParams()));
                }
                if (prec.getPostParams() != null && !prec.getPostParams().isEmpty()) {
                    precObj.add("postParams", matrixToJsonArray(prec.getPostParams()));
                }

                precsArr.add(precObj);
            }
        }
        return precsArr;
    }

    // ========================================================================
    // DISTRIBUTION SERIALIZATION
    // ========================================================================

    private static JsonObject serializeDistribution(Distribution dist) {
        JsonObject obj = new JsonObject();
        String name = dist.getName();
        obj.addProperty("type", name);

        if (dist instanceof Immediate) {
            // No params
        } else if (dist instanceof Disabled) {
            obj.addProperty("type", "Disabled");
        } else if ("Exp".equals(name)) {
            JsonObject params = new JsonObject();
            params.addProperty("lambda", ((Number) dist.getParam(1).getValue()).doubleValue());
            obj.add("params", params);
        } else if ("Det".equals(name)) {
            JsonObject params = new JsonObject();
            params.addProperty("value", ((Number) dist.getParam(1).getValue()).doubleValue());
            obj.add("params", params);
        } else if ("Erlang".equals(name)) {
            JsonObject params = new JsonObject();
            params.addProperty("lambda", ((Number) dist.getParam(1).getValue()).doubleValue());
            params.addProperty("k", ((Number) dist.getParam(2).getValue()).intValue());
            obj.add("params", params);
        } else if ("HyperExp".equals(name)) {
            // see _kb/09-ldes-and-cache.md (HyperExp wire form: p/lambda both length n)
            JsonObject params = new JsonObject();
            double[] p = ((HyperExp) dist).getP();
            double[] lambda = ((HyperExp) dist).getLambda();
            JsonArray pArr = new JsonArray();
            for (int i = 0; i < p.length; i++) {
                pArr.add(p[i]);
            }
            params.add("p", pArr);
            JsonArray lambdaArr = new JsonArray();
            for (int i = 0; i < lambda.length; i++) {
                lambdaArr.add(lambda[i]);
            }
            params.add("lambda", lambdaArr);
            obj.add("params", params);
        } else if ("Gamma".equals(name)) {
            JsonObject params = new JsonObject();
            params.addProperty("alpha", ((Number) dist.getParam(1).getValue()).doubleValue());
            params.addProperty("beta", ((Number) dist.getParam(2).getValue()).doubleValue());
            obj.add("params", params);
        } else if ("Lognormal".equals(name)) {
            JsonObject params = new JsonObject();
            params.addProperty("mu", ((Number) dist.getParam(1).getValue()).doubleValue());
            params.addProperty("sigma", ((Number) dist.getParam(2).getValue()).doubleValue());
            obj.add("params", params);
        } else if ("Uniform".equals(name)) {
            JsonObject params = new JsonObject();
            params.addProperty("a", ((Number) dist.getParam(1).getValue()).doubleValue());
            params.addProperty("b", ((Number) dist.getParam(2).getValue()).doubleValue());
            obj.add("params", params);
        } else if ("Weibull".equals(name)) {
            JsonObject params = new JsonObject();
            params.addProperty("alpha", ((Number) dist.getParam(1).getValue()).doubleValue());
            params.addProperty("beta", ((Number) dist.getParam(2).getValue()).doubleValue());
            obj.add("params", params);
        } else if ("Pareto".equals(name)) {
            JsonObject params = new JsonObject();
            params.addProperty("alpha", ((Number) dist.getParam(1).getValue()).doubleValue());
            params.addProperty("scale", ((Number) dist.getParam(2).getValue()).doubleValue());
            obj.add("params", params);
        } else if ("Normal".equals(name)) {
            JsonObject params = new JsonObject();
            params.addProperty("mu", ((Number) dist.getParam(1).getValue()).doubleValue());
            params.addProperty("sigma", ((Number) dist.getParam(2).getValue()).doubleValue());
            obj.add("params", params);
        } else if ("Geometric".equals(name)) {
            JsonObject params = new JsonObject();
            params.addProperty("p", ((Number) dist.getParam(1).getValue()).doubleValue());
            obj.add("params", params);
        } else if ("Binomial".equals(name)) {
            JsonObject params = new JsonObject();
            params.addProperty("n", ((Number) dist.getParam(1).getValue()).intValue());
            params.addProperty("p", ((Number) dist.getParam(2).getValue()).doubleValue());
            obj.add("params", params);
        } else if ("Poisson".equals(name)) {
            JsonObject params = new JsonObject();
            params.addProperty("lambda", ((Number) dist.getParam(1).getValue()).doubleValue());
            obj.add("params", params);
        } else if ("Bernoulli".equals(name)) {
            JsonObject params = new JsonObject();
            params.addProperty("p", ((Number) dist.getParam(1).getValue()).doubleValue());
            obj.add("params", params);
        } else if ("DiscreteUniform".equals(name)) {
            JsonObject params = new JsonObject();
            params.addProperty("min", ((Number) dist.getParam(1).getValue()).doubleValue());
            params.addProperty("max", ((Number) dist.getParam(2).getValue()).doubleValue());
            obj.add("params", params);
        } else if ("Zipf".equals(name)) {
            JsonObject params = new JsonObject();
            // Zipf params: 1=p (Matrix), 2=x (Matrix), 3=s (double), 4=n (int)
            params.addProperty("s", ((Number) dist.getParam(3).getValue()).doubleValue());
            params.addProperty("n", ((Number) dist.getParam(4).getValue()).intValue());
            obj.add("params", params);
        } else if (dist instanceof DiscreteSampler) {
            obj.addProperty("type", "DiscreteSampler");
            JsonObject params = new JsonObject();
            Matrix pMat = (Matrix) dist.getParam(1).getValue();
            Matrix xMat = (Matrix) dist.getParam(2).getValue();
            JsonArray pArr = new JsonArray();
            JsonArray xArr = new JsonArray();
            for (int k = 0; k < pMat.length(); k++) {
                pArr.add(pMat.get(k));
            }
            for (int k = 0; k < xMat.length(); k++) {
                xArr.add(xMat.get(k));
            }
            params.add("p", pArr);
            params.add("x", xArr);
            obj.add("params", params);
        } else if (dist instanceof Replayer) {
            // see _kb/09-ldes-and-cache.md for the Replayer/Trace wire-format rationale
            obj.addProperty("type", "Replayer");
            Replayer rep = (Replayer) dist;
            JsonObject params = new JsonObject();
            String fileName = rep.getFileName();
            if (fileName == null || fileName.isEmpty()) {
                line_warning(mfilename(new Object() {
                }), "%s distribution has no backing trace file; only its mean and an "
                        + "APH fit are saved, and the sample path will not survive the round trip.",
                        dist.getName());
            } else {
                params.addProperty("fileName", fileName);
            }
            params.addProperty("mean", dist.getMean());
            obj.add("params", params);
            // APH fallback, mirroring linemodel_save.m and the Python writer: it
            // is what a reader uses when the trace file is not reachable.
            APH aphFit = APH.fitMeanAndSCV(dist.getMean(), dist.getSCV());
            JsonObject phObj = new JsonObject();
            phObj.add("alpha", matrixToJsonRowVector(aphFit.getInitProb()));
            phObj.add("T", matrixToJson2D(aphFit.getSubgenerator()));
            obj.add("ph", phObj);
        } else if (dist instanceof EmpiricalCDF) {
            // Stored as an n x 2 table [F, x]; the wire form splits the columns.
            obj.addProperty("type", "EmpiricalCDF");
            Matrix data = ((EmpiricalCDF) dist).getData();
            JsonObject params = new JsonObject();
            JsonArray fArr = new JsonArray();
            JsonArray xArr = new JsonArray();
            for (int k = 0; k < data.getNumRows(); k++) {
                fArr.add(data.get(k, 0));
                xArr.add(data.getNumCols() > 1 ? data.get(k, 1) : data.get(k, 0));
            }
            params.add("F", fArr);
            params.add("x", xArr);
            obj.add("params", params);
        } else if (dist instanceof MMDP2) {
            // MMDP2 extends MMDP; the scalar parameterization is stored verbatim.
            obj.addProperty("type", "MMDP2");
            JsonObject params = new JsonObject();
            params.addProperty("r0", ((Number) dist.getParam(1).getValue()).doubleValue());
            params.addProperty("r1", ((Number) dist.getParam(2).getValue()).doubleValue());
            params.addProperty("sigma0", ((Number) dist.getParam(3).getValue()).doubleValue());
            params.addProperty("sigma1", ((Number) dist.getParam(4).getValue()).doubleValue());
            obj.add("params", params);
        } else if (dist instanceof MMPP2) {
            obj.addProperty("type", "MMPP2");
            JsonObject params = new JsonObject();
            params.addProperty("lambda0", ((Number) dist.getParam(1).getValue()).doubleValue());
            params.addProperty("lambda1", ((Number) dist.getParam(2).getValue()).doubleValue());
            params.addProperty("sigma0", ((Number) dist.getParam(3).getValue()).doubleValue());
            params.addProperty("sigma1", ((Number) dist.getParam(4).getValue()).doubleValue());
            obj.add("params", params);
        } else if (dist instanceof jline.lang.processes.BMAP) {
            // see _kb/09-ldes-and-cache.md (BMAP must precede the MarkedMAP branch it subclasses)
            obj.addProperty("type", "BMAP");
            jline.lang.processes.BMAP bm = (jline.lang.processes.BMAP) dist;
            JsonArray dArr = new JsonArray();
            dArr.add(matrixToJson2D(bm.D(0)));
            for (int b = 1; b <= bm.getMaxBatchSize(); b++) {
                dArr.add(matrixToJson2D(bm.getProcess().get(1 + b)));
            }
            JsonObject bmapParams = new JsonObject();
            bmapParams.add("D", dArr);
            obj.add("params", bmapParams);
        } else if (dist instanceof jline.lang.processes.MarkedMMPP) {
            // see _kb/09-ldes-and-cache.md for the MarkedMMPP (M3PP) wire-format rationale
            obj.addProperty("type", "MarkedMMPP");
            jline.lang.processes.MarkedMMPP mm = (jline.lang.processes.MarkedMMPP) dist;
            JsonArray dArr = new JsonArray();
            for (int k = 0; k < mm.getProcess().size(); k++) {
                dArr.add(matrixToJson2D(mm.getProcess().get(k)));
            }
            JsonObject mmppParams = new JsonObject();
            mmppParams.add("D", dArr);
            mmppParams.addProperty("K", Math.max(0, mm.getProcess().size() - 2));
            obj.add("params", mmppParams);
        } else if (dist instanceof jline.lang.processes.MarkedMAP
                && !(dist instanceof jline.lang.processes.BMAP)) {
            // see _kb/09-ldes-and-cache.md for the MMAP (Marked MAP) wire-format rationale
            obj.addProperty("type", "MMAP");
            jline.lang.processes.MarkedMAP mk = (jline.lang.processes.MarkedMAP) dist;
            JsonObject mmapObj = new JsonObject();
            mmapObj.add("D0", matrixToJson2D(mk.D(0)));
            JsonArray d1kArr = new JsonArray();
            for (int k = 1; k <= mk.getNumberOfTypes(); k++) {
                d1kArr.add(matrixToJson2D(mk.getProcess().get(1 + k)));
            }
            mmapObj.add("D1k", d1kArr);
            obj.add("mmap", mmapObj);
        } else if (dist instanceof DMAP) {
            // see _kb/09-ldes-and-cache.md (DMAP extends MarkovModulated, not MAP)
            obj.addProperty("type", "DMAP");
            DMAP dm = (DMAP) dist;
            JsonObject params = new JsonObject();
            params.add("D0", matrixToJson2D(dm.D(0)));
            params.add("D1", matrixToJson2D(dm.D(1)));
            obj.add("params", params);
        } else if (dist instanceof MAP) {
            obj.addProperty("type", "MAP");
            MAP mapDist = (MAP) dist;
            JsonObject mapObj = new JsonObject();
            mapObj.add("D0", matrixToJson2D(mapDist.D(0)));
            mapObj.add("D1", matrixToJson2D(mapDist.D(1)));
            obj.add("map", mapObj);
        } else if (dist instanceof RAP) {
            obj.addProperty("type", "RAP");
            RAP rap = (RAP) dist;
            JsonObject params = new JsonObject();
            params.add("H0", matrixToJson2D(rap.getH0()));
            params.add("H1", matrixToJson2D(rap.getH1()));
            obj.add("params", params);
        } else if (dist instanceof ME) {
            obj.addProperty("type", "ME");
            ME me = (ME) dist;
            JsonObject params = new JsonObject();
            params.add("alpha", matrixToJsonRowVector(me.getAlpha()));
            params.add("A", matrixToJson2D(me.getA()));
            obj.add("params", params);
        } else if (dist instanceof APH) {
            obj.addProperty("type", "APH");
            APH aphDist = (APH) dist;
            JsonObject phObj = new JsonObject();
            phObj.add("alpha", matrixToJsonRowVector(aphDist.getInitProb()));
            phObj.add("T", matrixToJson2D(aphDist.getSubgenerator()));
            obj.add("ph", phObj);
        } else if (dist instanceof Coxian) {
            obj.addProperty("type", "Coxian");
            Coxian coxDist = (Coxian) dist;
            JsonObject params = new JsonObject();
            Matrix mu = coxDist.getMu();
            Matrix phi = coxDist.getPhi();
            JsonArray muArr = new JsonArray();
            JsonArray phiArr = new JsonArray();
            for (int k = 0; k < mu.getNumElements(); k++) {
                muArr.add(mu.get(k));
            }
            for (int k = 0; k < phi.getNumElements(); k++) {
                phiArr.add(phi.get(k));
            }
            params.add("mu", muArr);
            params.add("phi", phiArr);
            obj.add("params", params);
        } else if (dist instanceof PH) {
            obj.addProperty("type", "PH");
            PH phDist = (PH) dist;
            JsonObject phObj = new JsonObject();
            phObj.add("alpha", matrixToJsonRowVector(phDist.getInitProb()));
            phObj.add("T", matrixToJson2D(phDist.getSubgenerator()));
            obj.add("ph", phObj);
        } else if (dist instanceof Expolynomial) {
            // see _kb/09-ldes-and-cache.md for the Expolynomial wire-format rationale
            obj.addProperty("type", "Expolynomial");
            Expolynomial expDist = (Expolynomial) dist;
            JsonObject ep = new JsonObject();
            ep.addProperty("density", expDist.getDensity());
            ep.addProperty("eft", expDist.getEft());
            double lft = expDist.getLft();
            if (Double.isInfinite(lft) || lft >= GlobalConstants.Inf) {
                ep.addProperty("lft", "Inf");
            } else {
                ep.addProperty("lft", lft);
            }
            obj.add("expolynomial", ep);
        } else if (dist instanceof NHPP) {
            obj.addProperty("type", "NHPP");
            NHPP nhppDist = (NHPP) dist;
            JsonObject nhppParams = new JsonObject();
            JsonArray bpArr = new JsonArray();
            JsonArray nhppRatesArr = new JsonArray();
            for (double b : nhppDist.getBreakpoints()) {
                bpArr.add(b);
            }
            for (double r : nhppDist.getRates()) {
                nhppRatesArr.add(r);
            }
            nhppParams.add("breakpoints", bpArr);
            nhppParams.add("rates", nhppRatesArr);
            nhppParams.addProperty("cyclic", nhppDist.isCyclic());
            obj.add("params", nhppParams);
        } else if (dist instanceof MAPt || dist instanceof PHt) {
            // Segment matrices go out as a JSON array of matrices, one per segment, so a
            // single-segment schedule keeps its nesting instead of collapsing to a vector.
            boolean isMapt = dist instanceof MAPt;
            obj.addProperty("type", isMapt ? "MAPt" : "PHt");
            JsonObject schedParams = new JsonObject();
            JsonArray schedBp = new JsonArray();
            double[] bps = isMapt ? ((MAPt) dist).getBreakpoints() : ((PHt) dist).getBreakpoints();
            for (double b : bps) {
                schedBp.add(b);
            }
            schedParams.add("breakpoints", schedBp);
            java.util.List<Matrix> firstSeg = isMapt
                    ? ((MAPt) dist).getD0Segments() : ((PHt) dist).getAlphaSegments();
            java.util.List<Matrix> secondSeg = isMapt
                    ? ((MAPt) dist).getD1Segments() : ((PHt) dist).getSSegments();
            JsonArray firstArr = new JsonArray();
            JsonArray secondArr = new JsonArray();
            for (int k = 0; k < firstSeg.size(); k++) {
                firstArr.add(matrixToJsonArray(firstSeg.get(k)));
                secondArr.add(matrixToJsonArray(secondSeg.get(k)));
            }
            schedParams.add(isMapt ? "D0" : "alpha", firstArr);
            schedParams.add(isMapt ? "D1" : "S", secondArr);
            schedParams.addProperty("cyclic",
                    isMapt ? ((MAPt) dist).isCyclic() : ((PHt) dist).isCyclic());
            obj.add("params", schedParams);
        } else if (dist instanceof Prior) {
            obj.addProperty("type", "Prior");
            obj.addProperty("kind", "discrete");
            Prior priorDist = (Prior) dist;
            JsonArray distsArr = new JsonArray();
            for (int i = 0; i < priorDist.getNumAlternatives(); i++) {
                distsArr.add(serializeDistribution(priorDist.getAlternative(i)));
            }
            obj.add("distributions", distsArr);
            JsonArray probsArr = new JsonArray();
            double[] probs = priorDist.getProbabilities();
            for (double p : probs) {
                probsArr.add(p);
            }
            obj.add("probabilities", probsArr);
        } else {
            // see _kb/09-ldes-and-cache.md (unhandled distribution type: type name + mean/SCV, never silently Exp)
            line_warning(mfilename(new Object() {
            }), "No JSON encoding is implemented for distribution type %s; saving its mean "
                    + "and SCV only. The reloaded model will use an APH matching those two moments.",
                    name);
            JsonObject params = new JsonObject();
            params.addProperty("mean", dist.getMean());
            params.addProperty("scv", dist.getSCV());
            obj.add("params", params);
        }

        return obj;
    }

    // ========================================================================
    // DISTRIBUTION DESERIALIZATION
    // ========================================================================

    /**
     * Rebuilds a distribution from a fit specification, i.e. from its moments
     * rather than from explicit parameters. The declared type selects the family
     * to fit; when that family cannot honour the requested moments, the fit falls
     * back to an exponential of the given mean, as the reference Python reader
     * does.
     *
     * @param type the declared distribution type
     * @param fit  the fit block, holding the method and the target moments
     * @return the fitted distribution
     */
    private static Distribution deserializeFit(String type, JsonObject fit) {
        String method = fit.get("method").getAsString();

        // see _kb/09-ldes-and-cache.md for the fitCentral/fitRawMoments deserialization rationale
        if ("fitCentral".equals(method) || "fitRawMoments".equals(method)) {
            double[] m = fitMoments(fit, method);
            if ("fitRawMoments".equals(method)) {
                if ("Coxian".equals(type) || "Cox2".equals(type)) {
                    // Cox2 is parameterized centrally; convert m1,m2,m3 -> mean,var,skew.
                    double mean1 = m[0];
                    double var = m[1] - mean1 * mean1;
                    double sd = Math.sqrt(Math.max(var, 0.0));
                    double skew = (sd > 0)
                            ? (m[2] - 3 * mean1 * m[1] + 2 * mean1 * mean1 * mean1) / (sd * sd * sd)
                            : 0.0;
                    return Cox2.fitCentral(mean1, var, skew);
                }
                return APH.fitRawMoments(m[0], m[1], m[2]);
            }
            if ("Coxian".equals(type) || "Cox2".equals(type)) {
                return Cox2.fitCentral(m[0], m[1], m[2]);
            }
            return APH.fitCentral(m[0], m[1], m[2]);
        }

        double mean = fit.get("mean").getAsDouble();

        if ("fitMeanAndSCV".equals(method)) {
            double scv = fit.get("scv").getAsDouble();
            if ("Erlang".equals(type)) {
                return Erlang.fitMeanAndSCV(mean, scv);
            } else if ("HyperExp".equals(type)) {
                return HyperExp.fitMeanAndSCV(mean, scv);
            } else if ("Gamma".equals(type)) {
                return Gamma.fitMeanAndSCV(mean, scv);
            } else if ("Lognormal".equals(type)) {
                return Lognormal.fitMeanAndSCV(mean, scv);
            }
            return Exp.fitMean(mean);
        } else if ("fitMeanAndOrder".equals(method)) {
            long order = fit.get("order").getAsLong();
            if ("Erlang".equals(type)) {
                return Erlang.fitMeanAndOrder(mean, order);
            }
            return Exp.fitMean(mean);
        } else if (!"fitMean".equals(method)) {
            line_warning(mfilename(new Object() {
            }), "Unrecognized fit method '%s' for distribution type %s; fitting the mean alone.",
                    method, type);
        }

        // fitMean, and any unrecognised method, are honoured on the mean alone.
        if ("Erlang".equals(type)) {
            long order = fit.has("order") ? fit.get("order").getAsLong() : 1L;
            return Erlang.fitMeanAndOrder(mean, order);
        } else if ("Det".equals(type)) {
            return new Det(mean);
        }
        return Exp.fitMean(mean);
    }

    /**
     * Extracts the three target moments of a {@code fitCentral} /
     * {@code fitRawMoments} fit block. They may be given as a {@code moments}
     * array, or as named scalars: {@code mean}/{@code var}/{@code skew} for a
     * central fit, {@code m1}/{@code m2}/{@code m3} for a raw-moment fit.
     *
     * @param fit    the fit block
     * @param method the fit method name, used in the error message
     * @return the three moments in the order the corresponding APH factory expects
     */
    private static double[] fitMoments(JsonObject fit, String method) {
        if (fit.has("moments")) {
            JsonArray ma = fit.getAsJsonArray("moments");
            if (ma.size() < 3) {
                throw new IllegalArgumentException(
                        method + " requires 3 moments, the file carries " + ma.size());
            }
            return new double[]{ma.get(0).getAsDouble(), ma.get(1).getAsDouble(),
                    ma.get(2).getAsDouble()};
        }
        if ("fitRawMoments".equals(method)) {
            if (fit.has("m1") && fit.has("m2") && fit.has("m3")) {
                return new double[]{fit.get("m1").getAsDouble(), fit.get("m2").getAsDouble(),
                        fit.get("m3").getAsDouble()};
            }
        } else if (fit.has("mean") && fit.has("var") && fit.has("skew")) {
            return new double[]{fit.get("mean").getAsDouble(), fit.get("var").getAsDouble(),
                    fit.get("skew").getAsDouble()};
        }
        throw new IllegalArgumentException(
                method + " requires either a 'moments' array of length 3 or the corresponding "
                        + "named moment fields");
    }

    /**
     * <p>The {@code params} object of a distribution record, or a named failure.</p>
     *
     * <p>Nearly every branch of {@link #deserializeDistribution} reads {@code params}
     * and then dereferences a key inside it. A record written in a foreign dialect
     * carries no {@code params} at all, and {@code getAsJsonObject} answers null
     * rather than throwing, so the first key read surfaced as a bare
     * NullPointerException naming neither the field nor the distribution. A
     * wire-format mismatch has to arrive as a diagnosis, not a stack trace.</p>
     *
     * @param obj  the distribution record
     * @param type the record's declared type, so the message names which one failed
     * @return the params object, never null
     */
    private static JsonObject requireParams(JsonObject obj, String type) {
        JsonElement el = obj.get("params");
        if (el == null || !el.isJsonObject()) {
            throw new IllegalArgumentException(
                    "distribution of type '" + type + "' carries no 'params' object");
        }
        return el.getAsJsonObject();
    }

    private static Distribution deserializeDistribution(JsonObject obj) {
        String type = obj.get("type").getAsString();

        if ("Immediate".equals(type)) {
            return Immediate.getInstance();
        } else if ("Disabled".equals(type)) {
            return new Disabled();
        } else if ("Expolynomial".equals(type)) {
            // see _kb/09-ldes-and-cache.md for the Expolynomial deserialization rationale
            if (!obj.has("expolynomial")) {
                throw new IllegalArgumentException(
                        "Expolynomial distribution record has no 'expolynomial' object: the "
                                + "density expression, eft and lft cannot be recovered");
            }
            JsonObject ep = obj.getAsJsonObject("expolynomial");
            if (!ep.has("density") || !ep.has("eft") || !ep.has("lft")) {
                throw new IllegalArgumentException(
                        "Expolynomial distribution record must carry 'density', 'eft' and 'lft'");
            }
            String density = ep.get("density").getAsString();
            double eft = ep.get("eft").getAsDouble();
            JsonElement lftElem = ep.get("lft");
            double lft;
            if (lftElem.isJsonPrimitive() && lftElem.getAsJsonPrimitive().isString()) {
                lft = Double.POSITIVE_INFINITY;
            } else {
                lft = lftElem.getAsDouble();
            }
            return new Expolynomial(density, eft, lft);
        }

        // see _kb/09-ldes-and-cache.md (fit-specification vs explicit-parameters precedence)
        if (!obj.has("params") && obj.has("fit")) {
            return deserializeFit(type, obj.getAsJsonObject("fit"));
        }

        if ("Exp".equals(type)) {
            JsonObject params = requireParams(obj, type);
            // "lambda" is canonical; "rate" is accepted as a read-side alias
            // because the published manual example and the Python reader use it.
            JsonElement lambdaEl = params.has("lambda") ? params.get("lambda") : params.get("rate");
            if (lambdaEl == null) {
                throw new IllegalArgumentException("Exp distribution declares neither lambda nor rate");
            }
            return new Exp(lambdaEl.getAsDouble());
        } else if ("Det".equals(type)) {
            JsonObject params = requireParams(obj, type);
            double value = params.get("value").getAsDouble();
            return new Det(value);
        } else if ("Erlang".equals(type)) {
            JsonObject params = requireParams(obj, type);
            double lambda = params.get("lambda").getAsDouble();
            int k = params.get("k").getAsInt();
            return new Erlang(lambda, k);
        } else if ("HyperExp".equals(type)) {
            JsonObject params = requireParams(obj, type);
            JsonArray pArr = params.getAsJsonArray("p");
            JsonArray lambdaArr = params.getAsJsonArray("lambda");
            int n = lambdaArr.size();
            if (n == 2) {
                return new HyperExp(pArr.get(0).getAsDouble(),
                        lambdaArr.get(0).getAsDouble(), lambdaArr.get(1).getAsDouble());
            }
            // see _kb/09-ldes-and-cache.md (HyperExp wire form: p/lambda both length n)
            if (pArr.size() != n) {
                throw new IllegalArgumentException(
                        "HyperExp declares " + pArr.size() + " probabilities for " + n + " rates");
            }
            double[] pVec = new double[n];
            double[] lambdaVec = new double[n];
            for (int k = 0; k < n; k++) {
                pVec[k] = pArr.get(k).getAsDouble();
                lambdaVec[k] = lambdaArr.get(k).getAsDouble();
            }
            return new HyperExp(pVec, lambdaVec);
        } else if ("Gamma".equals(type)) {
            JsonObject params = requireParams(obj, type);
            double alpha = params.get("alpha").getAsDouble();
            double beta = params.get("beta").getAsDouble();
            return new Gamma(alpha, beta);
        } else if ("Lognormal".equals(type)) {
            JsonObject params = requireParams(obj, type);
            double mu = params.get("mu").getAsDouble();
            double sigma = params.get("sigma").getAsDouble();
            return new Lognormal(mu, sigma);
        } else if ("Uniform".equals(type)) {
            JsonObject params = requireParams(obj, type);
            double a = params.get("a").getAsDouble();
            double b = params.get("b").getAsDouble();
            return new Uniform(a, b);
        } else if ("Weibull".equals(type)) {
            JsonObject params = requireParams(obj, type);
            double alpha = params.get("alpha").getAsDouble();  // scale
            double beta = params.get("beta").getAsDouble();    // shape
            return new Weibull(beta, alpha);  // constructor: Weibull(shape, scale)
        } else if ("Pareto".equals(type)) {
            JsonObject params = requireParams(obj, type);
            double alpha = params.get("alpha").getAsDouble();
            double scale = params.has("scale") ? params.get("scale").getAsDouble()
                    : params.get("beta").getAsDouble();
            return new Pareto(alpha, scale);
        } else if ("Normal".equals(type)) {
            JsonObject params = requireParams(obj, type);
            double mu = params.get("mu").getAsDouble();
            double sigma = params.get("sigma").getAsDouble();
            return new Normal(mu, sigma);
        } else if ("Geometric".equals(type)) {
            JsonObject params = requireParams(obj, type);
            double p = params.get("p").getAsDouble();
            return new Geometric(p);
        } else if ("Binomial".equals(type)) {
            JsonObject params = requireParams(obj, type);
            int n = params.get("n").getAsInt();
            double p = params.get("p").getAsDouble();
            return new Binomial(n, p);
        } else if ("Poisson".equals(type)) {
            JsonObject params = requireParams(obj, type);
            double lambda = params.get("lambda").getAsDouble();
            return new Poisson(lambda);
        } else if ("Bernoulli".equals(type)) {
            JsonObject params = requireParams(obj, type);
            double p = params.get("p").getAsDouble();
            return new Bernoulli(p);
        } else if ("DiscreteUniform".equals(type)) {
            JsonObject params = requireParams(obj, type);
            double min = params.get("min").getAsDouble();
            double max = params.get("max").getAsDouble();
            return new DiscreteUniform(min, max);
        } else if ("Zipf".equals(type)) {
            JsonObject params = requireParams(obj, type);
            double s = params.get("s").getAsDouble();
            int n = params.get("n").getAsInt();
            return new Zipf(s, n);
        } else if ("DiscreteSampler".equals(type)) {
            JsonObject params = requireParams(obj, type);
            JsonArray pArr = params.getAsJsonArray("p");
            JsonArray xArr = params.getAsJsonArray("x");
            int n = pArr.size();
            Matrix pMat = new Matrix(1, n);
            Matrix xMat = new Matrix(1, n);
            for (int k = 0; k < n; k++) {
                pMat.set(0, k, pArr.get(k).getAsDouble());
                xMat.set(0, k, xArr.get(k).getAsDouble());
            }
            return new DiscreteSampler(pMat, xMat);
        } else if ("MMPP2".equals(type)) {
            JsonObject p = requireParams(obj, type);
            return new MMPP2(
                p.get("lambda0").getAsDouble(),
                p.get("lambda1").getAsDouble(),
                p.get("sigma0").getAsDouble(),
                p.get("sigma1").getAsDouble()
            );
        } else if ("NHPP".equals(type)) {
            JsonObject params = requireParams(obj, type);
            JsonArray bpArr = params.getAsJsonArray("breakpoints");
            JsonArray rateArr = params.getAsJsonArray("rates");
            int n = rateArr.size();
            double[] breakpoints = new double[n + 1];
            double[] rates = new double[n];
            for (int k = 0; k < n; k++) {
                breakpoints[k] = bpArr.get(k).getAsDouble();
                rates[k] = rateArr.get(k).getAsDouble();
            }
            breakpoints[n] = bpArr.get(n).getAsDouble();
            // Absent 'cyclic' means cyclic, matching the constructor default.
            boolean cyclic = !params.has("cyclic") || params.get("cyclic").getAsBoolean();
            return new NHPP(breakpoints, rates, cyclic);
        } else if ("MAPt".equals(type) || "PHt".equals(type)) {
            JsonObject params = requireParams(obj, type);
            JsonArray bpArr = params.getAsJsonArray("breakpoints");
            boolean isMapt = "MAPt".equals(type);
            JsonArray firstArr = params.getAsJsonArray(isMapt ? "D0" : "alpha");
            JsonArray secondArr = params.getAsJsonArray(isMapt ? "D1" : "S");
            int n = secondArr.size();
            double[] breakpoints = new double[n + 1];
            for (int k = 0; k <= n; k++) {
                breakpoints[k] = bpArr.get(k).getAsDouble();
            }
            java.util.List<Matrix> firstSeg = new java.util.ArrayList<Matrix>();
            java.util.List<Matrix> secondSeg = new java.util.ArrayList<Matrix>();
            for (int k = 0; k < n; k++) {
                firstSeg.add(jsonToMatrix2DLenient(firstArr.get(k)));
                secondSeg.add(jsonToMatrix2DLenient(secondArr.get(k)));
            }
            // Absent 'cyclic' means cyclic, matching the constructor default.
            boolean cyclic = !params.has("cyclic") || params.get("cyclic").getAsBoolean();
            if (isMapt) {
                return new MAPt(breakpoints, firstSeg, secondSeg, cyclic);
            }
            return new PHt(breakpoints, firstSeg, secondSeg, cyclic);
        } else if ("MMDP2".equals(type)) {
            JsonObject p = requireParams(obj, type);
            return new MMDP2(
                p.get("r0").getAsDouble(),
                p.get("r1").getAsDouble(),
                p.get("sigma0").getAsDouble(),
                p.get("sigma1").getAsDouble()
            );
        } else if ("MAP".equals(type)) {
            JsonObject mapObj = obj.getAsJsonObject("map");
            Matrix D0 = jsonToMatrix2DLenient(mapObj.get("D0"));
            Matrix D1 = jsonToMatrix2DLenient(mapObj.get("D1"));
            return new MAP(D0, D1);
        } else if ("DMAP".equals(type)) {
            JsonObject params = requireParams(obj, type);
            Matrix D0 = jsonToMatrix2DLenient(params.get("D0"));
            Matrix D1 = jsonToMatrix2DLenient(params.get("D1"));
            return new DMAP(D0, D1);
        } else if ("RAP".equals(type)) {
            JsonObject params = requireParams(obj, type);
            Matrix H0 = jsonToMatrix2DLenient(params.get("H0"));
            Matrix H1 = jsonToMatrix2DLenient(params.get("H1"));
            return new RAP(H0, H1);
        } else if ("ME".equals(type)) {
            JsonObject params = requireParams(obj, type);
            Matrix alpha = jsonToRowVector(params.getAsJsonArray("alpha"));
            Matrix A = jsonToMatrix2DLenient(params.get("A"));
            return new ME(alpha, A);
        } else if ("BMAP".equals(type)) {
            // D[0]=D0, D[b]=D_b for batch size b; the aggregate D1 is recomputed
            // by the constructor. The array lives under "params", which is what
            // linemodel_save.m and linemodel_io.py write and what
            // linemodel_io.py reads; the top-level form is what this class used
            // to write and is accepted so old model.json files still load.
            JsonArray dArr = blockArray(obj, "D");
            if (dArr == null || dArr.size() < 2) {
                throw new IllegalArgumentException("BMAP requires D0 and at least one batch matrix");
            }
            Matrix D0 = jsonToMatrix2DLenient(dArr.get(0));
            Matrix[] Db = new Matrix[dArr.size() - 1];
            for (int b = 1; b < dArr.size(); b++) {
                Db[b - 1] = jsonToMatrix2DLenient(dArr.get(b));
            }
            return new jline.lang.processes.BMAP(D0, Db);
        } else if ("MarkedMMPP".equals(type)) {
            // Process cell in M3A layout {D0, D1_agg, D11..D1K}, taken verbatim.
            // Under "params", as above.
            JsonArray dArr = blockArray(obj, "D");
            if (dArr == null || dArr.size() < 2) {
                throw new IllegalArgumentException("MarkedMMPP requires at least D0 and D1");
            }
            jline.util.matrix.MatrixCell cell = new jline.util.matrix.MatrixCell();
            for (int k = 0; k < dArr.size(); k++) {
                cell.set(k, jsonToMatrix2DLenient(dArr.get(k)));
            }
            return new jline.lang.processes.MarkedMMPP(cell);
        } else if ("EmpiricalCDF".equals(type) || "EmpiricalCdf".equals(type)) {
            JsonObject params = requireParams(obj, type);
            JsonArray fArr = params.getAsJsonArray("F");
            JsonArray xArr = params.getAsJsonArray("x");
            int n = Math.min(fArr.size(), xArr.size());
            Matrix cdfData = new Matrix(n, 1);
            Matrix xData = new Matrix(n, 1);
            for (int k = 0; k < n; k++) {
                cdfData.set(k, 0, fArr.get(k).getAsDouble());
                xData.set(k, 0, xArr.get(k).getAsDouble());
            }
            return new EmpiricalCDF(cdfData, xData);
        } else if ("Replayer".equals(type) || "Trace".equals(type)) {
            // see _kb/09-ldes-and-cache.md for the Replayer/Trace fallback rationale
            JsonObject params = obj.has("params") ? obj.getAsJsonObject("params") : null;
            String fileName = (params != null && params.has("fileName"))
                    ? params.get("fileName").getAsString() : null;
            if (fileName != null && new File(fileName).isFile()) {
                return new Replayer(fileName);
            }
            if (fileName != null) {
                line_warning(mfilename(new Object() {
                }), "Trace file '%s' is not readable; falling back to the phase-type fit "
                        + "saved alongside it.", fileName);
            } else {
                line_warning(mfilename(new Object() {
                }), "Replayer carries no trace file name; falling back to the phase-type fit "
                        + "saved alongside it.");
            }
            if (obj.has("ph")) {
                JsonObject phObj = obj.getAsJsonObject("ph");
                return new PH(jsonToRowVector(phObj.getAsJsonArray("alpha")),
                        jsonToMatrix2D(phObj.getAsJsonArray("T")));
            }
            if (params != null && params.has("mean")) {
                double mean = params.get("mean").getAsDouble();
                double scv = params.has("scv") ? params.get("scv").getAsDouble() : 1.0;
                return APH.fitMeanAndSCV(mean, scv);
            }
            throw new IllegalArgumentException(
                    "Replayer carries neither a readable trace file, a phase-type fit nor a mean");
        } else if ("MMAP".equals(type)) {
            // see _kb/09-ldes-and-cache.md for the MMAP (Marked MAP) wire-format rationale
            JsonObject mmapObj = obj.getAsJsonObject("mmap");
            Matrix D0 = jsonToMatrix2DLenient(mmapObj.get("D0"));
            JsonArray d1kArr = mmapObj.getAsJsonArray("D1k");
            Matrix[] Dk = new Matrix[d1kArr.size()];
            for (int k = 0; k < d1kArr.size(); k++) {
                Dk[k] = jsonToMatrix2DLenient(d1kArr.get(k));
            }
            return new jline.lang.processes.MMAP(D0, Dk);
        } else if ("APH".equals(type)) {
            JsonObject phObj = obj.getAsJsonObject("ph");
            Matrix alpha = jsonToRowVector(phObj.getAsJsonArray("alpha"));
            Matrix T = jsonToMatrix2D(phObj.getAsJsonArray("T"));
            return new APH(alpha, T);
        } else if ("Coxian".equals(type) || "Cox2".equals(type)) {
            // see _kb/09-ldes-and-cache.md for the Cox2/Coxian alias rationale
            if ("Cox2".equals(type) && obj.has("params")) {
                JsonObject cx = requireParams(obj, type);
                if (cx.has("mu1") && cx.has("mu2") && cx.has("phi1")) {
                    return new Cox2(cx.get("mu1").getAsDouble(),
                            cx.get("mu2").getAsDouble(), cx.get("phi1").getAsDouble());
                }
            }
            // Try mu/phi params format first
            if (obj.has("params")) {
                JsonObject params = requireParams(obj, type);
                if (params.has("mu") && params.has("phi")) {
                    JsonArray muArr = params.getAsJsonArray("mu");
                    JsonArray phiArr = params.getAsJsonArray("phi");
                    java.util.List<Double> muList = new java.util.ArrayList<Double>();
                    java.util.List<Double> phiList = new java.util.ArrayList<Double>();
                    for (int k = 0; k < muArr.size(); k++) {
                        muList.add(muArr.get(k).getAsDouble());
                    }
                    for (int k = 0; k < phiArr.size(); k++) {
                        phiList.add(phiArr.get(k).getAsDouble());
                    }
                    return new Coxian(muList, phiList);
                }
            }
            // Legacy: D0/D1 or alpha/T format
            JsonObject phObj = obj.getAsJsonObject("ph");
            if (phObj != null) {
                if (phObj.has("alpha") && phObj.has("T")) {
                    Matrix alpha = jsonToRowVector(phObj.getAsJsonArray("alpha"));
                    Matrix T = jsonToMatrix2D(phObj.getAsJsonArray("T"));
                    return new PH(alpha, T);
                }
                if (phObj.has("D0") && phObj.has("D1")) {
                    Matrix D0 = jsonToMatrix2D(phObj.getAsJsonArray("D0"));
                    int n = D0.getNumRows();
                    Matrix alpha = new Matrix(1, n);
                    alpha.set(0, 0, 1.0);
                    return new PH(alpha, D0);
                }
            }
            // Final fallback: mean/scv
            if (obj.has("params")) {
                JsonObject params = requireParams(obj, type);
                double mean = params.get("mean").getAsDouble();
                double scv = params.get("scv").getAsDouble();
                return Coxian.fitMeanAndSCV(mean, scv);
            }
            throw new IllegalArgumentException(
                    "Coxian distribution declares neither mu/phi, a phase-type nor mean/scv");
        } else if ("PH".equals(type)) {
            JsonObject phObj = obj.getAsJsonObject("ph");
            Matrix alpha = jsonToRowVector(phObj.getAsJsonArray("alpha"));
            Matrix T = jsonToMatrix2D(phObj.getAsJsonArray("T"));
            return new PH(alpha, T);
        } else if ("Prior".equals(type)) {
            // The continuous form (parameter density plus factory template) is
            // MATLAB/C++ only; the JAR Prior is discrete, and its keys read as a
            // discrete set would be a null alternative array
            String priorKind = obj.has("kind") ? obj.get("kind").getAsString() : "discrete";
            if ("continuous".equals(priorKind)) {
                throw new RuntimeException(
                        "a continuous Prior (paramDist plus a factory template) is carried by the "
                        + "MATLAB and C++ codebases only; the JAR Prior is discrete. Re-save the "
                        + "model with the Prior expanded by Prior.discretize, or solve it in "
                        + "MATLAB or C++");
            }
            JsonArray distsArr = obj.getAsJsonArray("distributions");
            JsonArray probsArr = obj.getAsJsonArray("probabilities");
            java.util.List<Distribution> alternatives = new java.util.ArrayList<Distribution>();
            for (int i = 0; i < distsArr.size(); i++) {
                alternatives.add(deserializeDistribution(distsArr.get(i).getAsJsonObject()));
            }
            double[] probs = new double[probsArr.size()];
            for (int i = 0; i < probsArr.size(); i++) {
                probs[i] = probsArr.get(i).getAsDouble();
            }
            return new Prior(alternatives, probs);
        }

        // see _kb/09-ldes-and-cache.md (unrecognized type on read: APH fit from mean/SCV, never Immediate/Exp)
        if (obj.has("params")) {
            JsonObject params = requireParams(obj, type);
            if (params.has("mean")) {
                double mean = params.get("mean").getAsDouble();
                double scv = params.has("scv") ? params.get("scv").getAsDouble() : 1.0;
                line_warning(mfilename(new Object() {
                }), "No JSON decoding is implemented for distribution type %s; rebuilding an APH "
                        + "matching its saved mean (%g) and SCV (%g).", type, mean, scv);
                if (mean > 0 && scv > 0) {
                    return APH.fitMeanAndSCV(mean, scv);
                }
                if (mean > 0) {
                    return new Det(mean);
                }
                return Immediate.getInstance();
            }
        }
        throw new IllegalArgumentException("Distribution type '" + type
                + "' is not recognized and carries no mean to fit; the model cannot be loaded.");
    }

    // ========================================================================
    // NETWORK LOAD
    // ========================================================================

    private static Network loadNetwork(JsonObject modelObj) throws IOException {
        String name = modelObj.has("name") ? modelObj.get("name").getAsString() : "model";
        Network model = new Network(name);

        JsonArray nodesArr = modelObj.getAsJsonArray("nodes");
        JsonArray classesArr = modelObj.getAsJsonArray("classes");
        JsonObject routingRaw = modelObj.has("routing") ? modelObj.getAsJsonObject("routing") : null;
        JsonObject routingObj = null;
        if (routingRaw != null) {
            if (routingRaw.has("matrix")) {
                routingObj = routingRaw.getAsJsonObject("matrix");
            } else {
                // Flat format (no wrapper)
                routingObj = routingRaw;
            }
        }

        // Two-pass node creation: first pass creates nodes, second sets services
        // Phase 1: Create nodes
        Map<String, Node> nodeMap = new LinkedHashMap<String, Node>();
        Map<String, JsonObject> nodeJsonMap = new LinkedHashMap<String, JsonObject>();

        for (JsonElement nodeEl : nodesArr) {
            JsonObject nodeObj = nodeEl.getAsJsonObject();
            String nodeName = nodeObj.get("name").getAsString();
            String nodeType = nodeObj.get("type").getAsString();
            nodeJsonMap.put(nodeName, nodeObj);

            Node node;
            if ("Source".equals(nodeType)) {
                node = new Source(model, nodeName);
            } else if ("Sink".equals(nodeType)) {
                node = new Sink(model, nodeName);
            } else if ("Delay".equals(nodeType)) {
                node = new Delay(model, nodeName);
                // see _kb/09-ldes-and-cache.md (LineModelIO.serializeNetworkNodes handles Delay/Queue jointly)
                if (nodeObj.has("buffer")) {
                    ((Delay) node).setCapacity(nodeObj.get("buffer").getAsInt());
                }
            } else if ("Queue".equals(nodeType)) {
                String schedStr = nodeObj.has("scheduling") ? nodeObj.get("scheduling").getAsString() : "PS";
                SchedStrategy sched = parseSchedStrategy(schedStr);
                node = new Queue(model, nodeName, sched);
                if (nodeObj.has("servers")) {
                    int servers = nodeObj.get("servers").getAsInt();
                    if (servers > 0 && servers != Integer.MAX_VALUE) {
                        ((Queue) node).setNumberOfServers(servers);
                    }
                }
                if (nodeObj.has("buffer")) {
                    ((Queue) node).setCapacity(nodeObj.get("buffer").getAsInt());
                }
            } else if ("Fork".equals(nodeType)) {
                node = new Fork(model, nodeName);
                if (nodeObj.has("tasksPerLink")) {
                    ((Fork) node).setTasksPerLink(nodeObj.get("tasksPerLink").getAsInt());
                }
            } else if ("Join".equals(nodeType)) {
                node = new Join(model, nodeName);
            } else if ("Router".equals(nodeType)) {
                node = new Router(model, nodeName);
            } else if ("ClassSwitch".equals(nodeType)) {
                node = new ClassSwitch(model, nodeName);
            } else if ("Cache".equals(nodeType)) {
                // see _kb/12-interfaces-and-docs.md (cache config emitted with flat JAR-compatible keys)
                JsonElement itemsEl = cacheField(nodeObj, "items", "numItems");
                if (itemsEl == null) {
                    throw new IOException("Cache node '" + nodeName
                            + "' declares neither cache.items nor numItems");
                }
                int numItems = itemsEl.getAsInt();
                JsonElement capEl = cacheField(nodeObj, "capacity", "itemLevelCap");
                if (capEl == null) {
                    throw new IOException("Cache node '" + nodeName
                            + "' declares neither cache.capacity nor itemLevelCap");
                }
                // A single-level cache may arrive as a bare scalar.
                Matrix itemLevelCap = capEl.isJsonArray()
                        ? jsonToRowVector(capEl.getAsJsonArray())
                        : Matrix.singleton(capEl.getAsDouble());
                JsonElement replEl = cacheField(nodeObj, "replacement", "replacementStrategy");
                if (replEl == null) {
                    replEl = cacheField(nodeObj, "replacementStrategy", "replacement");
                }
                ReplacementStrategy replPolicy = parseReplacementStrategy(
                        replEl != null ? replEl.getAsString() : "LRU");
                Cache cacheNode = new Cache(model, nodeName, numItems, itemLevelCap, replPolicy);
                JsonElement apEl = cacheField(nodeObj, "admissionProb", "admissionProb");
                if (apEl != null) {
                    cacheNode.setAdmissionProb(apEl.getAsDouble());
                }
                JsonElement szEl = cacheField(nodeObj, "itemSizes", "itemSizes");
                if (szEl != null) {
                    cacheNode.setItemSizes(szEl.isJsonArray()
                            ? jsonToRowVector(szEl.getAsJsonArray())
                            : Matrix.singleton(szEl.getAsDouble()));
                }
                JsonElement ccEl = cacheField(nodeObj, "costCaps", "costCaps");
                if (ccEl != null) {
                    cacheNode.setCostCaps(ccEl.isJsonArray()
                            ? jsonToRowVector(ccEl.getAsJsonArray())
                            : Matrix.singleton(ccEl.getAsDouble()));
                }
                node = cacheNode;
            } else if ("Place".equals(nodeType)) {
                // see _kb/09-ldes-and-cache.md for the queueing-place (QPN) wire-format rationale
                if (nodeObj.has("service")) {
                    String pSchedStr = nodeObj.has("scheduling")
                            ? nodeObj.get("scheduling").getAsString() : "FCFS";
                    node = new Place(model, nodeName, parseSchedStrategy(pSchedStr));
                } else {
                    node = new Place(model, nodeName);
                }
            } else if ("Transition".equals(nodeType)) {
                node = new Transition(model, nodeName);
            } else {
                // Unknown node type; skip
                continue;
            }
            nodeMap.put(nodeName, node);
        }

        // Phase 2: Create classes (after all nodes exist)
        Map<String, JobClass> classMap = new LinkedHashMap<String, JobClass>();
        for (JsonElement classEl : classesArr) {
            JsonObject classObj = classEl.getAsJsonObject();
            String className = classObj.get("name").getAsString();
            String classType = classObj.get("type").getAsString();

            JobClass jc;
            if ("Open".equals(classType)) {
                int priority = classObj.has("priority") ? classObj.get("priority").getAsInt() : 0;
                jc = new OpenClass(model, className, priority);
                // An explicit reference-station override (another Source)
                if (classObj.has("refNode")) {
                    Node refNode = nodeMap.get(classObj.get("refNode").getAsString());
                    if (refNode instanceof Source) {
                        try {
                            jc.setReferenceStation((Source) refNode);
                        } catch (Exception e) {
                            throw new IOException("Cannot set reference station '"
                                    + refNode.getName() + "' for open class '" + className + "'", e);
                        }
                    } else {
                        throw new IOException("Reference station '"
                                + classObj.get("refNode").getAsString()
                                + "' of open class '" + className + "' is not a Source");
                    }
                }
            } else if ("Closed".equals(classType) || "SelfLooping".equals(classType)) {
                double population = classObj.get("population").getAsDouble();
                String refNodeName = classObj.has("refNode") ? classObj.get("refNode").getAsString() : null;
                Station refStat = null;
                if (refNodeName != null) {
                    Node refNode = nodeMap.get(refNodeName);
                    if (refNode instanceof Station) {
                        refStat = (Station) refNode;
                    }
                }
                int priority = classObj.has("priority") ? classObj.get("priority").getAsInt() : 0;
                if (refStat == null) {
                    throw new IOException("Reference station '" + refNodeName + "' not found for closed class '" + className + "'");
                }
                if ("SelfLooping".equals(classType)) {
                    jc = new SelfLoopingClass(model, className, Math.round(population), refStat, priority);
                } else {
                    jc = new ClosedClass(model, className, population, refStat, priority);
                }
            } else if ("Signal".equals(classType)) {
                int priority = classObj.has("priority") ? classObj.get("priority").getAsInt() : 0;
                SignalType sigType = SignalType.NEGATIVE;
                if (classObj.has("signalType")) {
                    sigType = SignalType.fromText(classObj.get("signalType").getAsString());
                }
                // see _kb/09-ldes-and-cache.md (a bare Signal class is resolved open-vs-closed from the network)
                String openOrClosed = classObj.has("openOrClosed") ? classObj.get("openOrClosed").getAsString() : null;
                if (openOrClosed == null) {
                    boolean hasSource = false;
                    for (Node n : nodeMap.values()) {
                        if (n instanceof Source) {
                            hasSource = true;
                            break;
                        }
                    }
                    openOrClosed = hasSource ? "Open" : "Closed";
                }
                if ("Closed".equals(openOrClosed)) {
                    String refNodeName = classObj.has("refNode") ? classObj.get("refNode").getAsString() : null;
                    Station refStat = null;
                    if (refNodeName != null) {
                        Node refNode = nodeMap.get(refNodeName);
                        if (refNode instanceof Station) {
                            refStat = (Station) refNode;
                        }
                    }
                    // see _kb/09-ldes-and-cache.md (a bare Signal class is resolved open-vs-closed from the network)
                    if (refStat == null) {
                        for (Node n : nodeMap.values()) {
                            if (n instanceof Station && !(n instanceof Source)) {
                                refStat = (Station) n;
                                break;
                            }
                        }
                    }
                    if (refStat == null) {
                        throw new IOException("Reference station not found for closed signal '" + className + "'");
                    }
                    jc = new ClosedSignal(model, className, sigType, refStat, priority);
                } else {
                    jc = new OpenSignal(model, className, sigType, priority);
                }
                // Removal distribution
                if (classObj.has("removalDistribution")) {
                    Distribution remDist = deserializeDistribution(classObj.getAsJsonObject("removalDistribution"));
                    if (remDist instanceof DiscreteDistribution) {
                        DiscreteDistribution dremDist = (DiscreteDistribution) remDist;
                        if (jc instanceof OpenSignal) ((OpenSignal) jc).setRemovalDistribution(dremDist);
                        else if (jc instanceof ClosedSignal) ((ClosedSignal) jc).setRemovalDistribution(dremDist);
                    }
                }
                // Removal policy
                if (classObj.has("removalPolicy")) {
                    RemovalPolicy rp = RemovalPolicy.fromText(classObj.get("removalPolicy").getAsString());
                    if (jc instanceof OpenSignal) ((OpenSignal) jc).setRemovalPolicy(rp);
                    else if (jc instanceof ClosedSignal) ((ClosedSignal) jc).setRemovalPolicy(rp);
                }
            } else {
                throw new IOException("Unknown class type: " + classType);
            }
            // Class deadline
            if (classObj.has("deadline")) {
                double deadline = classObj.get("deadline").getAsDouble();
                if (Double.isFinite(deadline)) {
                    jc.setDeadline(deadline);
                }
            }

            // Reference class within its chain
            if (classObj.has("isReferenceClass")) {
                jc.setReferenceClass(classObj.get("isReferenceClass").getAsBoolean());
            }

            // Class-scoped patience (a node-scoped "patience" overrides it)
            if (classObj.has("patience")) {
                Distribution patDist = deserializeDistribution(classObj.getAsJsonObject("patience"));
                ImpatienceType impType = ImpatienceType.RENEGING;
                if (classObj.has("impatienceType")) {
                    String itStr = classObj.get("impatienceType").getAsString();
                    if ("balking".equals(itStr)) impType = ImpatienceType.BALKING;
                    else if ("retrial".equals(itStr)) impType = ImpatienceType.RETRIAL;
                }
                jc.setPatience(impType, patDist);
            }

            classMap.put(className, jc);
        }

        // Resolve replySignalClass associations (after all classes are created).
        // The index is 1-based, as JobClass.setReplySignalClassIndex documents.
        List<JobClass> loadedClasses = model.getClasses();
        for (JsonElement classEl : classesArr) {
            JsonObject classObj = classEl.getAsJsonObject();
            if (!classObj.has("replySignalClass")) continue;
            JobClass jc = classMap.get(classObj.get("name").getAsString());
            JobClass replyCls = classMap.get(classObj.get("replySignalClass").getAsString());
            if (jc == null || replyCls == null) continue;
            int replyIdx = loadedClasses.indexOf(replyCls);
            if (replyIdx >= 0) {
                jc.setReplySignalClassIndex(replyIdx + 1);
            }
        }

        // Resolve spawnClass associations (after all classes are created).
        // The index is 1-based, as JobClass.setSpawnClassIndex documents.
        for (JsonElement classEl : classesArr) {
            JsonObject classObj = classEl.getAsJsonObject();
            if (!classObj.has("spawnClass")) continue;
            JobClass jc = classMap.get(classObj.get("name").getAsString());
            JobClass spawnCls = classMap.get(classObj.get("spawnClass").getAsString());
            if (jc == null || spawnCls == null) continue;
            int spawnIdx = loadedClasses.indexOf(spawnCls);
            if (spawnIdx >= 0) {
                jc.setSpawnClassIndex(spawnIdx + 1);
            }
        }

        // Resolve signal targetClass associations (after all classes are created)
        for (JsonElement classEl : classesArr) {
            JsonObject classObj = classEl.getAsJsonObject();
            if (!"Signal".equals(classObj.get("type").getAsString())) continue;
            if (!classObj.has("targetClass")) continue;
            String sigName = classObj.get("name").getAsString();
            String targetName = classObj.get("targetClass").getAsString();
            JobClass sigCls = classMap.get(sigName);
            JobClass targetCls = classMap.get(targetName);
            if (sigCls != null && targetCls != null) {
                if (sigCls instanceof OpenSignal) ((OpenSignal) sigCls).forJobClass(targetCls);
                else if (sigCls instanceof ClosedSignal) {
                    ((ClosedSignal) sigCls).forJobClass(targetCls);
                    // see _kb/09-ldes-and-cache.md (a bare Signal class is resolved open-vs-closed from the network)
                    if (!classObj.has("refNode") && targetCls.getReferenceStation() != null) {
                        try {
                            sigCls.setReferenceStation(targetCls.getReferenceStation());
                        } catch (Exception e) {
                            // keep the fallback reference station
                        }
                    }
                }
            }
        }

        // Phase 3: Set service distributions and node-specific parameters
        List<JobClass> jobClasses = model.getClasses();
        for (Map.Entry<String, JsonObject> njEntry : nodeJsonMap.entrySet()) {
            String nodeName = njEntry.getKey();
            JsonObject nodeObj = njEntry.getValue();
            String nodeType = nodeObj.get("type").getAsString();
            Node node = nodeMap.get(nodeName);
            if (node == null) {
                continue;
            }

            if ("Source".equals(nodeType)) {
                Source src = (Source) node;
                JsonObject arrivals = null;
                if (nodeObj.has("service")) {
                    arrivals = nodeObj.getAsJsonObject("service");
                } else if (nodeObj.has("arrivals")) {
                    arrivals = nodeObj.getAsJsonObject("arrivals");
                }
                if (arrivals != null) {
                    for (Map.Entry<String, JsonElement> ae : arrivals.entrySet()) {
                        JobClass jc = classMap.get(ae.getKey());
                        if (jc != null) {
                            Distribution dist = deserializeDistribution(ae.getValue().getAsJsonObject());
                            src.setArrival(jc, dist);
                        }
                    }
                }
                // see _kb/09-ldes-and-cache.md (Batch arrivals: BatchArrival, Geo^X)
                if (nodeObj.has("arrivalBatch")) {
                    JsonObject batchObj = nodeObj.getAsJsonObject("arrivalBatch");
                    for (Map.Entry<String, JsonElement> be : batchObj.entrySet()) {
                        JobClass jc = classMap.get(be.getKey());
                        if (jc != null) {
                            Distribution batch = deserializeDistribution(be.getValue().getAsJsonObject());
                            // A DETERMINISTIC BATCH IS A LEGITIMATE BATCH SIZE. MATLAB's Det
                            // is declared `ContinuousDistribution & DiscreteDistribution`, so
                            // setArrivalBatch accepts it and the writer emits it; this class
                            // hierarchy has Det outside DiscreteDistribution, which refused
                            // the very models the reference writes. A Det on a positive
                            // integer is the degenerate DiscreteUniform on that value, so it
                            // is converted rather than rejected.
                            if (batch instanceof jline.lang.processes.Det) {
                                double v = batch.getMean();
                                if (v < 1 || Math.abs(v - Math.rint(v)) > 1e-9) {
                                    throw new RuntimeException("arrivalBatch for class '" + be.getKey()
                                            + "' is a Det on " + v + "; a batch must carry a positive "
                                            + "integer number of jobs");
                                }
                                batch = new jline.lang.processes.DiscreteUniform(Math.rint(v), Math.rint(v));
                            }
                            if (!(batch instanceof jline.lang.processes.DiscreteDistribution)) {
                                throw new RuntimeException("arrivalBatch for class '" + be.getKey()
                                        + "' deserialized to " + batch.getClass().getSimpleName()
                                        + ", which is not a DiscreteDistribution; a batch size must "
                                        + "be integer valued");
                            }
                            src.setArrivalBatch(jc, (jline.lang.processes.DiscreteDistribution) batch);
                        }
                    }
                }
                // see _kb/09-ldes-and-cache.md for the marked (MMAP) arrival-binding rationale
                if (nodeObj.has("markedClasses")) {
                    JsonArray markedArr = nodeObj.getAsJsonArray("markedClasses");
                    List<JobClass> markedClasses = new ArrayList<JobClass>();
                    for (int mi = 0; mi < markedArr.size(); mi++) {
                        JobClass mc = classMap.get(markedArr.get(mi).getAsString());
                        if (mc != null) {
                            markedClasses.add(mc);
                        }
                    }
                    if (!markedClasses.isEmpty()) {
                        Distribution first = src.getArrivalDistribution(markedClasses.get(0));
                        if (first instanceof jline.lang.processes.MarkedMAP) {
                            src.setMarkedArrival((jline.lang.processes.MarkedMAP) first, markedClasses);
                        }
                    }
                }
            } else if ("Place".equals(nodeType)) {
                Place place = (Place) node;
                // Queueing place: reconstruct the embedded-queue service
                // processes, server count and depository departure discipline.
                if (nodeObj.has("service")) {
                    JsonObject services = nodeObj.getAsJsonObject("service");
                    for (Map.Entry<String, JsonElement> se : services.entrySet()) {
                        JobClass jc = classMap.get(se.getKey());
                        if (jc != null) {
                            Distribution dist = deserializeDistribution(se.getValue().getAsJsonObject());
                            place.setService(jc, dist);
                        }
                    }
                }
                if (nodeObj.has("servers")) {
                    int servers = nodeObj.get("servers").getAsInt();
                    if (servers > 0 && servers != Integer.MAX_VALUE) {
                        place.setNumberOfServers(servers);
                    }
                }
                if (nodeObj.has("departureDiscipline")) {
                    JsonObject ddObj = nodeObj.getAsJsonObject("departureDiscipline");
                    for (Map.Entry<String, JsonElement> de : ddObj.entrySet()) {
                        JobClass jc = classMap.get(de.getKey());
                        if (jc != null) {
                            String ddStr = de.getValue().getAsString();
                            jline.lang.constant.DepartureDiscipline disc =
                                    "FIFO".equalsIgnoreCase(ddStr)
                                            ? jline.lang.constant.DepartureDiscipline.FIFO
                                            : jline.lang.constant.DepartureDiscipline.Normal;
                            place.setDepartureDiscipline(jc, disc);
                        }
                    }
                }
                // see _kb/12-interfaces-and-docs.md (per-class buffer capacity must be serialized for Place nodes too)
                if (nodeObj.has("classCap")) {
                    JsonObject classCapObj = nodeObj.getAsJsonObject("classCap");
                    for (Map.Entry<String, JsonElement> ccEntry : classCapObj.entrySet()) {
                        JobClass jc = classMap.get(ccEntry.getKey());
                        if (jc != null) {
                            place.setClassCapacity(jc, ccEntry.getValue().getAsInt());
                        }
                    }
                }
            } else if ("Delay".equals(nodeType) || "Queue".equals(nodeType)) {
                Queue queue = (Queue) node;
                if (nodeObj.has("service")) {
                    JsonObject services = nodeObj.getAsJsonObject("service");
                    for (Map.Entry<String, JsonElement> se : services.entrySet()) {
                        JobClass jc = classMap.get(se.getKey());
                        if (jc != null) {
                            Distribution dist = deserializeDistribution(se.getValue().getAsJsonObject());
                            queue.setService(jc, dist);
                        }
                    }
                }
                // see _kb/12-interfaces-and-docs.md for the PAS/OI service-rate reconstruction rationale
                if (nodeObj.has("oiServiceRate") &&
                        (queue.getSchedStrategy() == SchedStrategy.PAS
                                || queue.getSchedStrategy() == SchedStrategy.OI)) {
                    final int K = model.getNumberOfClasses();
                    // OI queues keep a fixed zero swap graph (setSwapGraph rejects them).
                    if (nodeObj.has("swapGraph") && queue.getSchedStrategy() == SchedStrategy.PAS) {
                        JsonArray sg = nodeObj.getAsJsonArray("swapGraph");
                        Matrix G = new Matrix(K, K);
                        for (int i = 0; i < sg.size(); i++) {
                            JsonArray row = sg.get(i).getAsJsonArray();
                            for (int j = 0; j < row.size(); j++) {
                                G.set(i, j, row.get(j).getAsDouble());
                            }
                        }
                        queue.setSwapGraph(G);
                    }
                    JsonObject tblJson = nodeObj.getAsJsonObject("oiServiceRate");
                    final HashMap<String, Double> rateTbl = new HashMap<String, Double>();
                    double maxRateTmp = 0.0;
                    for (Map.Entry<String, JsonElement> e : tblJson.entrySet()) {
                        double v = e.getValue().getAsDouble();
                        rateTbl.put(e.getKey(), v);
                        if (v > maxRateTmp) maxRateTmp = v;
                    }
                    final double maxRate = maxRateTmp;
                    // see _kb/12-interfaces-and-docs.md for the PAS/OI service-rate reconstruction rationale
                    final int[] cutoffs;
                    if (nodeObj.has("oiCutoffs")) {
                        JsonArray ca = nodeObj.getAsJsonArray("oiCutoffs");
                        cutoffs = new int[ca.size()];
                        for (int i = 0; i < ca.size(); i++) {
                            cutoffs[i] = ca.get(i).getAsInt();
                        }
                    } else {
                        cutoffs = null;
                    }
                    SerializableFunction<Matrix, Double> muFun = (Matrix prefix) -> {
                        int n = prefix.getNumCols();
                        if (n <= 0) {
                            return 0.0;
                        }
                        int[] cnt = new int[K];
                        for (int i = 0; i < n; i++) {
                            int cid = (int) prefix.get(0, i);
                            if (cid >= 0 && cid < K) cnt[cid]++;
                        }
                        if (cutoffs != null) {
                            for (int r = 0; r < K; r++) {
                                if (r < cutoffs.length && cnt[r] > cutoffs[r]) {
                                    cnt[r] = cutoffs[r];
                                }
                            }
                        }
                        StringBuilder sb = new StringBuilder();
                        for (int i = 0; i < K; i++) {
                            if (i > 0) sb.append(",");
                            sb.append(cnt[i]);
                        }
                        Double rate = rateTbl.get(sb.toString());
                        return rate != null ? rate : maxRate;
                    };
                    queue.setService(muFun);
                }
                if (nodeObj.has("schedParams")) {
                    JsonObject schedParams = nodeObj.getAsJsonObject("schedParams");
                    for (Map.Entry<String, JsonElement> sp : schedParams.entrySet()) {
                        JobClass jc = classMap.get(sp.getKey());
                        if (jc != null) {
                            queue.setSchedStrategyPar(jc, sp.getValue().getAsDouble());
                        }
                    }
                }
                // Per-class buffer capacity
                if (nodeObj.has("classCap")) {
                    JsonObject classCapObj = nodeObj.getAsJsonObject("classCap");
                    for (Map.Entry<String, JsonElement> ccEntry : classCapObj.entrySet()) {
                        JobClass jc = classMap.get(ccEntry.getKey());
                        if (jc != null) {
                            queue.setClassCap(jc, ccEntry.getValue().getAsInt());
                        }
                    }
                }
                // Drop rules
                if (nodeObj.has("dropRule")) {
                    JsonObject dropRuleObj = nodeObj.getAsJsonObject("dropRule");
                    for (Map.Entry<String, JsonElement> drEntry : dropRuleObj.entrySet()) {
                        JobClass jc = classMap.get(drEntry.getKey());
                        if (jc != null) {
                            DropStrategy ds = parseDropStrategy(drEntry.getValue().getAsString());
                            queue.setDropRule(jc, ds);
                        }
                    }
                }
                // Load-dependent scaling
                if (nodeObj.has("loadDependence")) {
                    JsonObject ldObj = nodeObj.getAsJsonObject("loadDependence");
                    String ldType = ldObj.has("type") ? ldObj.get("type").getAsString() : "";
                    if ("loadDependent".equals(ldType) && ldObj.has("scaling")) {
                        JsonArray scalingArr = ldObj.getAsJsonArray("scaling");
                        Matrix alpha = new Matrix(1, scalingArr.size());
                        for (int k = 0; k < scalingArr.size(); k++) {
                            alpha.set(0, k, scalingArr.get(k).getAsDouble());
                        }
                        queue.setLoadDependence(alpha);
                    }
                }
                // see _kb/12-interfaces-and-docs.md (a class-dependent rate-scaling callable is materialized over the macrostate lattice)
                if (nodeObj.has("classDependence")) {
                    JsonObject cdObj = nodeObj.getAsJsonObject("classDependence");
                    String cdType = cdObj.has("type") ? cdObj.get("type").getAsString() : "";
                    if ("classDependent".equals(cdType) && cdObj.has("scaling")) {
                        final int Kcd = model.getNumberOfClasses();
                        final HashMap<String, double[]> betaTbl = new HashMap<String, double[]>();
                        for (Map.Entry<String, JsonElement> e : cdObj.getAsJsonObject("scaling").entrySet()) {
                            JsonArray va = e.getValue().getAsJsonArray();
                            double[] v = new double[va.size()];
                            for (int r = 0; r < va.size(); r++) {
                                v[r] = va.get(r).getAsDouble();
                            }
                            betaTbl.put(e.getKey(), v);
                        }
                        final int[] cdCutoffs;
                        if (cdObj.has("cutoffs")) {
                            JsonArray ca = cdObj.getAsJsonArray("cutoffs");
                            cdCutoffs = new int[ca.size()];
                            for (int i = 0; i < ca.size(); i++) {
                                cdCutoffs[i] = ca.get(i).getAsInt();
                            }
                        } else {
                            cdCutoffs = null;
                        }
                        SerializableFunction<Matrix, Matrix> betaFun = (Matrix nvec) -> {
                            int[] n = new int[Kcd];
                            for (int r = 0; r < Kcd && r < nvec.getNumElements(); r++) {
                                n[r] = (int) Math.round(nvec.get(r));
                                if (n[r] < 0) n[r] = 0;
                                if (cdCutoffs != null && r < cdCutoffs.length && n[r] > cdCutoffs[r]) {
                                    n[r] = cdCutoffs[r];
                                }
                            }
                            StringBuilder sb = new StringBuilder();
                            for (int r = 0; r < Kcd; r++) {
                                if (r > 0) sb.append(",");
                                sb.append(n[r]);
                            }
                            double[] v = betaTbl.get(sb.toString());
                            Matrix out = new Matrix(1, Kcd);
                            for (int r = 0; r < Kcd; r++) {
                                // A state absent from the table is neutral (beta=1),
                                // which leaves the nominal service rate unscaled.
                                out.set(0, r, (v != null && r < v.length) ? v[r] : 1.0);
                            }
                            return out;
                        };
                        // see _kb/12-interfaces-and-docs.md (a class-dependent rate-scaling callable is materialized over the macrostate lattice)
                        if (cdObj.has("peak")) {
                            JsonArray pa = cdObj.getAsJsonArray("peak");
                            Matrix peakVec = new Matrix(1, pa.size());
                            for (int i = 0; i < pa.size(); i++) {
                                peakVec.set(0, i, pa.get(i).getAsDouble());
                            }
                            queue.setLimitedClassDependence(betaFun, peakVec);
                        } else {
                            // LEGACY JSON, written before the peak became a wire
                            // key. The table IS the whole lattice, bounded by
                            // `cutoffs`, so max_n beta(n) is recoverable here in
                            // a way it never is from a user's handle -- the same
                            // recovery linemodel_load.m and linemodel_io.py make.
                            // Dropping the peak instead would defer the failure
                            // to the solver, or worse to a zeroed utilization.
                            final int[] peakNK = new int[Kcd];
                            for (int i = 0; i < Kcd; i++) {
                                peakNK[i] = (cdCutoffs != null && i < cdCutoffs.length) ? cdCutoffs[i] : 1;
                            }
                            double bmax = jline.api.pfqn.ld.CdPeakScaling.cd_peak_scaling(betaFun, peakNK, Kcd);
                            if (!(bmax > 0)) {
                                throw new RuntimeException("classDependence at node '" + queue.getName()
                                        + "' carries no peak and none can be recovered from its lattice table;"
                                        + " write the \"peak\" key.");
                            }
                            Matrix derived = new Matrix(1, Kcd);
                            for (int i = 0; i < Kcd; i++) derived.set(0, i, bmax);
                            queue.setLimitedClassDependence(betaFun, derived);
                        }
                    }
                }
                // jointDependence: non-product-form eta_i, twin of classDependence.
                if (nodeObj.has("jointDependence")) {
                    JsonObject jdObj = nodeObj.getAsJsonObject("jointDependence");
                    String jdType = jdObj.has("type") ? jdObj.get("type").getAsString() : "";
                    if ("jointDependent".equals(jdType) && jdObj.has("scaling")) {
                        final int Kjd = model.getNumberOfClasses();
                        final HashMap<String, double[]> etaTbl = new HashMap<String, double[]>();
                        for (Map.Entry<String, JsonElement> e : jdObj.getAsJsonObject("scaling").entrySet()) {
                            JsonArray va = e.getValue().getAsJsonArray();
                            double[] v = new double[va.size()];
                            for (int r = 0; r < va.size(); r++) {
                                v[r] = va.get(r).getAsDouble();
                            }
                            etaTbl.put(e.getKey(), v);
                        }
                        final int[] jdCutoffs;
                        if (jdObj.has("cutoffs")) {
                            JsonArray ca = jdObj.getAsJsonArray("cutoffs");
                            jdCutoffs = new int[ca.size()];
                            for (int i = 0; i < ca.size(); i++) {
                                jdCutoffs[i] = ca.get(i).getAsInt();
                            }
                        } else {
                            jdCutoffs = null;
                        }
                        SerializableFunction<Matrix, Matrix> etaFun = (Matrix nvec) -> {
                            int[] n = new int[Kjd];
                            for (int r = 0; r < Kjd && r < nvec.getNumElements(); r++) {
                                n[r] = (int) Math.round(nvec.get(r));
                                if (n[r] < 0) n[r] = 0;
                                if (jdCutoffs != null && r < jdCutoffs.length && n[r] > jdCutoffs[r]) {
                                    n[r] = jdCutoffs[r];
                                }
                            }
                            StringBuilder sb = new StringBuilder();
                            for (int r = 0; r < Kjd; r++) {
                                if (r > 0) sb.append(",");
                                sb.append(n[r]);
                            }
                            double[] v = etaTbl.get(sb.toString());
                            Matrix out = new Matrix(1, Kjd);
                            for (int r = 0; r < Kjd; r++) {
                                out.set(0, r, (v != null && r < v.length) ? v[r] : 1.0);
                            }
                            return out;
                        };
                        if (jdObj.has("peak")) {
                            JsonArray pa = jdObj.getAsJsonArray("peak");
                            Matrix peakVec = new Matrix(1, pa.size());
                            for (int i = 0; i < pa.size(); i++) {
                                peakVec.set(0, i, pa.get(i).getAsDouble());
                            }
                            queue.setLimitedJointDependence(etaFun, peakVec);
                        } else {
                            // Legacy JSON: recover the peak from the lattice
                            // table, as the classDependence branch above does.
                            final int[] peakNK = new int[Kjd];
                            for (int i = 0; i < Kjd; i++) {
                                peakNK[i] = (jdCutoffs != null && i < jdCutoffs.length) ? jdCutoffs[i] : 1;
                            }
                            double emax = jline.api.pfqn.ld.CdPeakScaling.cd_peak_scaling(etaFun, peakNK, Kjd);
                            if (!(emax > 0)) {
                                throw new RuntimeException("jointDependence at node '" + queue.getName()
                                        + "' carries no peak and none can be recovered from its lattice table;"
                                        + " write the \"peak\" key.");
                            }
                            Matrix derived = new Matrix(1, Kjd);
                            for (int i = 0; i < Kjd; i++) derived.set(0, i, emax);
                            queue.setLimitedJointDependence(etaFun, derived);
                        }
                    }
                }
                // Heterogeneous server types
                if (nodeObj.has("serverTypes")) {
                    for (JsonElement stEl : nodeObj.getAsJsonArray("serverTypes")) {
                        JsonObject stObj = stEl.getAsJsonObject();
                        String stName = stObj.get("name").getAsString();
                        int stCount = stObj.get("count").getAsInt();
                        ServerType st = new ServerType(stName, stCount);
                        // Compatible classes
                        if (stObj.has("compatibleClasses")) {
                            for (JsonElement ccEl : stObj.getAsJsonArray("compatibleClasses")) {
                                JobClass jc = classMap.get(ccEl.getAsString());
                                if (jc != null) {
                                    st.addCompatibleClass(jc);
                                }
                            }
                        }
                        queue.addServerType(st);
                        // Per-class service distributions
                        if (stObj.has("service")) {
                            JsonObject svcObj = stObj.getAsJsonObject("service");
                            for (Map.Entry<String, JsonElement> svcEntry : svcObj.entrySet()) {
                                JobClass jc = classMap.get(svcEntry.getKey());
                                if (jc != null) {
                                    Distribution dist = deserializeDistribution(svcEntry.getValue().getAsJsonObject());
                                    queue.setService(jc, st, dist);
                                }
                            }
                        }
                    }
                    // Scheduling policy
                    if (nodeObj.has("heteroSchedPolicy")) {
                        HeteroSchedPolicy hsp = HeteroSchedPolicy.fromText(nodeObj.get("heteroSchedPolicy").getAsString());
                        queue.setHeteroSchedPolicy(hsp);
                    }
                }

                // Balking
                if (nodeObj.has("balking")) {
                    JsonObject balkObj = nodeObj.getAsJsonObject("balking");
                    for (Map.Entry<String, JsonElement> balkEntry : balkObj.entrySet()) {
                        JobClass jc = classMap.get(balkEntry.getKey());
                        if (jc == null) continue;
                        JsonObject bjc = balkEntry.getValue().getAsJsonObject();
                        String stratStr = bjc.get("strategy").getAsString();
                        BalkingStrategy bs = BalkingStrategy.valueOf(stratStr);
                        List<BalkingThreshold> thresholds = new ArrayList<BalkingThreshold>();
                        JsonArray thArr = bjc.getAsJsonArray("thresholds");
                        for (JsonElement thEl : thArr) {
                            JsonObject tObj = thEl.getAsJsonObject();
                            int minJobs = tObj.get("minJobs").getAsInt();
                            int maxJobs = tObj.get("maxJobs").getAsInt();
                            if (maxJobs < 0) maxJobs = Integer.MAX_VALUE;
                            double prob = tObj.get("probability").getAsDouble();
                            thresholds.add(new BalkingThreshold(minJobs, maxJobs, prob));
                        }
                        queue.setBalking(jc, bs, thresholds);
                    }
                }

                // Retrial
                if (nodeObj.has("retrial")) {
                    JsonObject retObj = nodeObj.getAsJsonObject("retrial");
                    for (Map.Entry<String, JsonElement> retEntry : retObj.entrySet()) {
                        JobClass jc = classMap.get(retEntry.getKey());
                        if (jc == null) continue;
                        JsonObject rjc = retEntry.getValue().getAsJsonObject();
                        Distribution delayDist = deserializeDistribution(rjc.getAsJsonObject("delay"));
                        int maxAttempts = -1;
                        if (rjc.has("maxAttempts")) {
                            maxAttempts = rjc.get("maxAttempts").getAsInt();
                        }
                        queue.setRetrial(jc, delayDist, maxAttempts);
                    }
                }

                // Patience
                if (nodeObj.has("patience")) {
                    JsonObject patObj = nodeObj.getAsJsonObject("patience");
                    for (Map.Entry<String, JsonElement> patEntry : patObj.entrySet()) {
                        JobClass jc = classMap.get(patEntry.getKey());
                        if (jc == null) continue;
                        JsonObject pjc = patEntry.getValue().getAsJsonObject();
                        Distribution patDist = deserializeDistribution(pjc.getAsJsonObject("distribution"));
                        ImpatienceType impType = ImpatienceType.RENEGING;
                        if (pjc.has("impatienceType")) {
                            String itStr = pjc.get("impatienceType").getAsString();
                            if ("balking".equals(itStr)) impType = ImpatienceType.BALKING;
                            else if ("retrial".equals(itStr)) impType = ImpatienceType.RETRIAL;
                        }
                        queue.setPatience(jc, impType, patDist);
                    }
                }

                // Orbit impatience (abandonment from the retrial orbit)
                if (nodeObj.has("orbitImpatience")) {
                    JsonObject orbitObj = nodeObj.getAsJsonObject("orbitImpatience");
                    for (Map.Entry<String, JsonElement> oe : orbitObj.entrySet()) {
                        JobClass jc = classMap.get(oe.getKey());
                        if (jc == null) continue;
                        queue.setOrbitImpatience(jc,
                                deserializeDistribution(oe.getValue().getAsJsonObject()));
                    }
                }

                // Batch rejection probability (retrial queues)
                if (nodeObj.has("batchRejectProb")) {
                    JsonObject batchRejectObj = nodeObj.getAsJsonObject("batchRejectProb");
                    for (Map.Entry<String, JsonElement> be : batchRejectObj.entrySet()) {
                        JobClass jc = classMap.get(be.getKey());
                        if (jc == null) continue;
                        queue.setBatchRejectProbability(jc, be.getValue().getAsDouble());
                    }
                }

                // Job parallelism: servers seized at once by a job, per class
                if (nodeObj.has("serverParallelism")) {
                    JsonObject parallelismObj = nodeObj.getAsJsonObject("serverParallelism");
                    for (Map.Entry<String, JsonElement> pe : parallelismObj.entrySet()) {
                        JobClass jc = classMap.get(pe.getKey());
                        if (jc == null) continue;
                        queue.setServerParallelism(jc, pe.getValue().getAsInt());
                    }
                }

                // Immediate feedback
                if (nodeObj.has("immediateFeedback")) {
                    JsonObject immFeedObj = nodeObj.getAsJsonObject("immediateFeedback");
                    for (Map.Entry<String, JsonElement> fe : immFeedObj.entrySet()) {
                        JobClass jc = classMap.get(fe.getKey());
                        if (jc == null || !fe.getValue().getAsBoolean()) continue;
                        queue.setImmediateFeedback(jc);
                    }
                }

                // see _kb/09-ldes-and-cache.md (setupTime/delayOffTime are always emitted as a pair)
                if (nodeObj.has("setupTime") || nodeObj.has("delayOffTime")) {
                    JsonObject setupObj = nodeObj.has("setupTime")
                            ? nodeObj.getAsJsonObject("setupTime") : new JsonObject();
                    JsonObject delayOffObj = nodeObj.has("delayOffTime")
                            ? nodeObj.getAsJsonObject("delayOffTime") : new JsonObject();
                    for (Map.Entry<String, JsonElement> suEntry : setupObj.entrySet()) {
                        JobClass jc = classMap.get(suEntry.getKey());
                        if (jc == null) continue;
                        if (!delayOffObj.has(suEntry.getKey())) continue;
                        Distribution suDist = deserializeDistribution(suEntry.getValue().getAsJsonObject());
                        Distribution doffDist = deserializeDistribution(
                                delayOffObj.getAsJsonObject(suEntry.getKey()));
                        queue.setDelayOff(jc, suDist, doffDist);
                    }
                }

                // Server breakdown/repair, with the optional per-class degraded
                // down-server service.
                if (nodeObj.has("breakdown")) {
                    JsonObject bdObj = nodeObj.getAsJsonObject("breakdown");
                    if (!bdObj.has("failure") || !bdObj.has("repair")) {
                        throw new RuntimeException("Node \"" + queue.getName() + "\": \"breakdown\" "
                                + "requires both a \"failure\" and a \"repair\" distribution.");
                    }
                    Distribution failDist = deserializeDistribution(bdObj.getAsJsonObject("failure"));
                    Distribution repairDist = deserializeDistribution(bdObj.getAsJsonObject("repair"));
                    queue.setBreakdown(failDist, repairDist);
                    if (bdObj.has("downService")) {
                        JsonObject downSvcObj = bdObj.getAsJsonObject("downService");
                        for (Map.Entry<String, JsonElement> de : downSvcObj.entrySet()) {
                            JobClass jc = classMap.get(de.getKey());
                            if (jc == null) continue;
                            queue.setDownService(jc,
                                    deserializeDistribution(de.getValue().getAsJsonObject()));
                        }
                    }
                }

                // see _kb/09-ldes-and-cache.md (readers must restore pollingType before switchoverTimes)
                if (nodeObj.has("pollingType")) {
                    PollingType pt = PollingType.valueOf(nodeObj.get("pollingType").getAsString());
                    if (pt == PollingType.KLIMITED) {
                        int k = nodeObj.has("pollingPar") ? nodeObj.get("pollingPar").getAsInt() : 1;
                        queue.setPollingType(pt, k);
                    } else {
                        queue.setPollingType(pt);
                    }
                }

                // Switchover times: entries without "to" carry the per-class
                // polling form, entries with "to" the (from,to) pair form.
                if (nodeObj.has("switchoverTimes")) {
                    JsonArray soArr = nodeObj.getAsJsonArray("switchoverTimes");
                    for (int si = 0; si < soArr.size(); si++) {
                        JsonObject soj = soArr.get(si).getAsJsonObject();
                        JobClass fromJc = classMap.get(soj.get("from").getAsString());
                        if (fromJc == null || !soj.has("distribution")) continue;
                        Distribution soDist = deserializeDistribution(soj.getAsJsonObject("distribution"));
                        if (soj.has("to")) {
                            JobClass toJc = classMap.get(soj.get("to").getAsString());
                            if (toJc == null) continue;
                            queue.setSwitchover(fromJc, toJc, soDist);
                        } else {
                            queue.setSwitchover(fromJc, soDist);
                        }
                    }
                }
            } else if ("Join".equals(nodeType)) {
                Join join = (Join) node;
                // Support both "forkNode" (MATLAB/Python format) and "joinOf" (JAR format)
                String forkNodeKey = nodeObj.has("forkNode") ? "forkNode" : "joinOf";
                if (nodeObj.has(forkNodeKey)) {
                    String forkName = nodeObj.get(forkNodeKey).getAsString();
                    Node forkNode = nodeMap.get(forkName);
                    if (forkNode != null) {
                        join.joinOf = forkNode;
                    }
                }
                // Join strategy
                if (nodeObj.has("joinStrategy")) {
                    String jsStr = nodeObj.get("joinStrategy").getAsString();
                    JoinStrategy js = JoinStrategy.STD;
                    try {
                        // Map aliases from other codebases
                        if ("PARTIAL".equals(jsStr) || "QUORUM".equals(jsStr)) {
                            js = JoinStrategy.Quorum;
                        } else {
                            js = JoinStrategy.valueOf(jsStr);
                        }
                    } catch (IllegalArgumentException e) {
                        // keep STD as default
                    }
                    for (JobClass jc : jobClasses) {
                        join.setStrategy(jc, js);
                    }
                }
                // Join quorum
                if (nodeObj.has("joinQuorum")) {
                    int jq = nodeObj.get("joinQuorum").getAsInt();
                    Joiner joiner = (Joiner) join.getInput();
                    if (joiner != null) {
                        for (JobClass jc : jobClasses) {
                            joiner.setRequired(jc, jq);
                        }
                    }
                }
            } else if ("ClassSwitch".equals(nodeType)) {
                ClassSwitch cs = (ClassSwitch) node;
                if (nodeObj.has("classSwitchMatrix")) {
                    // Dict format: {"Class1": {"Class1": 0.3, "Class2": 0.7}, ...}
                    JsonObject csDict = nodeObj.getAsJsonObject("classSwitchMatrix");
                    ClassSwitchMatrix csm = cs.initClassSwitchMatrix();
                    for (Map.Entry<String, JsonElement> fromEntry : csDict.entrySet()) {
                        JobClass fromClass = classMap.get(fromEntry.getKey());
                        if (fromClass == null) continue;
                        int r = jobClasses.indexOf(fromClass);
                        if (r < 0) continue;
                        JsonObject toObj = fromEntry.getValue().getAsJsonObject();
                        for (Map.Entry<String, JsonElement> toEntry : toObj.entrySet()) {
                            JobClass toClass = classMap.get(toEntry.getKey());
                            if (toClass == null) continue;
                            int s = jobClasses.indexOf(toClass);
                            if (s < 0) continue;
                            csm.set(r, s, toEntry.getValue().getAsDouble());
                        }
                    }
                    cs.setClassSwitchingMatrix(csm);
                } else if (nodeObj.has("csMatrix")) {
                    // Legacy 2D array format: [[0.3, 0.7], [1.0, 0.0]]
                    JsonArray csMatrixArr = nodeObj.getAsJsonArray("csMatrix");
                    int K = csMatrixArr.size();
                    ClassSwitchMatrix csm = cs.initClassSwitchMatrix();
                    for (int r = 0; r < K; r++) {
                        JsonArray row = csMatrixArr.get(r).getAsJsonArray();
                        for (int s = 0; s < row.size(); s++) {
                            csm.set(r, s, row.get(s).getAsDouble());
                        }
                    }
                    cs.setClassSwitchingMatrix(csm);
                }
            } else if ("Cache".equals(nodeType)) {
                Cache cache = (Cache) node;
                // Hit/miss class mapping
                JsonElement hitEl = cacheField(nodeObj, "hitClass", "hitClass");
                if (hitEl != null) {
                    JsonObject hitMap = hitEl.getAsJsonObject();
                    for (Map.Entry<String, JsonElement> he : hitMap.entrySet()) {
                        JobClass inClass = classMap.get(he.getKey());
                        JobClass outClass = classMap.get(he.getValue().getAsString());
                        if (inClass != null && outClass != null) {
                            cache.setHitClass(inClass, outClass);
                        }
                    }
                }
                JsonElement missEl = cacheField(nodeObj, "missClass", "missClass");
                if (missEl != null) {
                    JsonObject missMap = missEl.getAsJsonObject();
                    for (Map.Entry<String, JsonElement> me : missMap.entrySet()) {
                        JobClass inClass = classMap.get(me.getKey());
                        JobClass outClass = classMap.get(me.getValue().getAsString());
                        if (inClass != null && outClass != null) {
                            cache.setMissClass(inClass, outClass);
                        }
                    }
                }
                // Popularity distributions
                JsonElement popEl = cacheField(nodeObj, "popularity", "popularity");
                if (popEl != null) {
                    JsonObject popObj = popEl.getAsJsonObject();
                    for (Map.Entry<String, JsonElement> pe : popObj.entrySet()) {
                        JobClass jc = classMap.get(pe.getKey());
                        if (jc != null) {
                            Distribution popDist = deserializeDistribution(pe.getValue().getAsJsonObject());
                            // Skip Disabled distributions (hit/miss classes have no popularity)
                            if (popDist.isDiscrete()) {
                                cache.setRead(jc, popDist);
                            }
                        }
                    }
                }
                // `itemClass`: for a cache network, the item each per-item class reads
                // (Cache.setItemReadClasses). Recorded directly: popularity and the
                // hit/miss switches are read from their own keys above, so only the
                // item mapping itself is missing, and inferring it from a one-hot
                // popularity would be ambiguous against a single-item popularity.
                JsonElement icEl = cacheField(nodeObj, "itemClass", "itemClass");
                if (icEl != null) {
                    JsonObject icObj = icEl.getAsJsonObject();
                    for (Map.Entry<String, JsonElement> ie : icObj.entrySet()) {
                        JobClass jc = classMap.get(ie.getKey());
                        if (jc != null) {
                            cache.setItemOfClass(jc, ie.getValue().getAsInt());
                        }
                    }
                }
                // see _kb/09-ldes-and-cache.md (Delayed-hit retrieval, api/retrieval)
                JsonElement rsEl = cacheField(nodeObj, "retrievalSystem", "retrievalSystem");
                if (rsEl != null) {
                    JsonObject rsObj = rsEl.getAsJsonObject();
                    int nItems = cache.getNumberOfItems();
                    JsonObject byClass = rsObj.getAsJsonObject("byClass");
                    for (Map.Entry<String, JsonElement> ce : byClass.entrySet()) {
                        JobClass jobinClass = classMap.get(ce.getKey());
                        if (jobinClass == null) continue;
                        JsonObject entry = ce.getValue().getAsJsonObject();
                        List<Integer> queueIndices = new ArrayList<Integer>();
                        if (entry.has("queues")) {
                            for (JsonElement qe : entry.getAsJsonArray("queues")) {
                                Node qn = nodeMap.get(qe.getAsString());
                                if (qn != null) queueIndices.add(qn.getNodeIndex());
                            }
                        }
                        JobClass[] retrievalClassByItem = new JobClass[nItems];
                        if (entry.has("items")) {
                            JsonObject itemsObj = entry.getAsJsonObject("items");
                            for (Map.Entry<String, JsonElement> ie : itemsObj.entrySet()) {
                                int it = Integer.parseInt(ie.getKey());
                                if (it >= 0 && it < nItems) {
                                    retrievalClassByItem[it] = classMap.get(ie.getValue().getAsString());
                                }
                            }
                        }
                        cache.attachRetrievalSystem(jobinClass, queueIndices, retrievalClassByItem);
                    }
                }
                // Access-cost (list-move) structure: per-item graph (shared by all
                // classes, as in Network.sanitize) or full per-class accessProb.
                JsonElement agEl = cacheField(nodeObj, "accessGraph", "accessGraph");
                JsonElement apbEl = cacheField(nodeObj, "accessProb", "accessProb");
                if (agEl != null) {
                    JsonArray gArr = agEl.getAsJsonArray();
                    int K = model.getNumberOfClasses();
                    int nItemsAc = cache.getNumberOfItems();
                    Matrix[][] ap = new Matrix[K][nItemsAc];
                    for (int it = 0; it < nItemsAc; it++) {
                        Matrix g = null;
                        if (it < gArr.size() && gArr.get(it).isJsonArray()
                                && gArr.get(it).getAsJsonArray().size() > 0) {
                            g = jsonToMatrix2D(gArr.get(it).getAsJsonArray());
                        }
                        for (int v = 0; v < K; v++) {
                            ap[v][it] = g;
                        }
                    }
                    cache.accessProb = ap;
                } else if (apbEl != null) {
                    JsonArray apArr = apbEl.getAsJsonArray();
                    Matrix[][] ap = new Matrix[apArr.size()][];
                    for (int v = 0; v < apArr.size(); v++) {
                        JsonArray rowArr = apArr.get(v).getAsJsonArray();
                        ap[v] = new Matrix[rowArr.size()];
                        for (int it = 0; it < rowArr.size(); it++) {
                            JsonArray mArr = rowArr.get(it).getAsJsonArray();
                            ap[v][it] = (mArr.size() > 0) ? jsonToMatrix2D(mArr) : null;
                        }
                    }
                    cache.accessProb = ap;
                }
                // Initial cache state [class counts | contents | retrieval bitmap]
                JsonElement isEl = cacheField(nodeObj, "initialState", "initialState");
                if (isEl != null) {
                    JsonArray stArr = isEl.getAsJsonArray();
                    Matrix st = new Matrix(1, stArr.size());
                    for (int k = 0; k < stArr.size(); k++) {
                        st.set(0, k, stArr.get(k).getAsDouble());
                    }
                    cache.setState(st);
                }
            }
        }

        // Phase 3b: Configure Transition modes (after all classes exist)
        for (Map.Entry<String, JsonObject> njEntry : nodeJsonMap.entrySet()) {
            String nodeName = njEntry.getKey();
            JsonObject nodeObj = njEntry.getValue();
            String nodeType = nodeObj.get("type").getAsString();
            if (!"Transition".equals(nodeType)) continue;
            if (!nodeObj.has("modes")) continue;

            Node node = nodeMap.get(nodeName);
            if (!(node instanceof Transition)) continue;
            Transition tnode = (Transition) node;

            JsonArray modesArr = nodeObj.getAsJsonArray("modes");
            for (JsonElement modeEl : modesArr) {
                JsonObject md = modeEl.getAsJsonObject();
                String modeName = md.has("name") ? md.get("name").getAsString() : "Mode";
                Mode mode = tnode.addMode(modeName);

                // Timing strategy (must be set before distribution, since
                // setTimingStrategy(IMMEDIATE) overwrites distribution with Immediate)
                if (md.has("timingStrategy")) {
                    String tsStr = md.get("timingStrategy").getAsString();
                    if ("IMMEDIATE".equals(tsStr)) {
                        tnode.setTimingStrategy(mode, TimingStrategy.IMMEDIATE);
                    } else {
                        tnode.setTimingStrategy(mode, TimingStrategy.TIMED);
                    }
                }
                // Distribution
                if (md.has("distribution")) {
                    Distribution dist = deserializeDistribution(md.getAsJsonObject("distribution"));
                    if (dist != null) {
                        tnode.setDistribution(mode, dist);
                    }
                }
                // Number of servers
                if (md.has("numServers")) {
                    JsonElement nsEl = md.get("numServers");
                    int ns;
                    if (nsEl.isJsonPrimitive() && nsEl.getAsJsonPrimitive().isString()) {
                        ns = "Infinity".equalsIgnoreCase(nsEl.getAsString()) ? Integer.MAX_VALUE : Integer.parseInt(nsEl.getAsString());
                    } else {
                        ns = nsEl.getAsInt();
                    }
                    tnode.setNumberOfServers(mode, ns);
                }
                // Firing priority
                if (md.has("firingPriority")) {
                    tnode.setFiringPriorities(mode, (int) md.get("firingPriority").getAsDouble());
                }
                // Firing weight
                if (md.has("firingWeight")) {
                    tnode.setFiringWeights(mode, md.get("firingWeight").getAsDouble());
                }
                // Enabling conditions
                if (md.has("enablingConditions")) {
                    for (JsonElement ecEl : md.getAsJsonArray("enablingConditions")) {
                        JsonObject ec = ecEl.getAsJsonObject();
                        Node ecNode = nodeMap.get(ec.get("node").getAsString());
                        JobClass ecClass = classMap.get(ec.get("class").getAsString());
                        if (ecNode instanceof Place && ecClass != null) {
                            tnode.setEnablingConditions(mode, ecClass, (Place) ecNode, (int) ec.get("count").getAsDouble());
                        }
                    }
                }
                // Inhibiting conditions
                if (md.has("inhibitingConditions")) {
                    for (JsonElement icEl : md.getAsJsonArray("inhibitingConditions")) {
                        JsonObject ic = icEl.getAsJsonObject();
                        Node icNode = nodeMap.get(ic.get("node").getAsString());
                        JobClass icClass = classMap.get(ic.get("class").getAsString());
                        if (icNode instanceof Place && icClass != null) {
                            tnode.setInhibitingConditions(mode, icClass, (Place) icNode, (int) ic.get("count").getAsDouble());
                        }
                    }
                }
                // Firing outcomes
                if (md.has("firingOutcomes")) {
                    for (JsonElement foEl : md.getAsJsonArray("firingOutcomes")) {
                        JsonObject fo = foEl.getAsJsonObject();
                        Node foNode = nodeMap.get(fo.get("node").getAsString());
                        JobClass foClass = classMap.get(fo.get("class").getAsString());
                        if (foNode != null && foClass != null) {
                            tnode.setFiringOutcome(mode, foClass, foNode, (int) fo.get("count").getAsDouble());
                        }
                    }
                }
                // Marking-dependent firing-rate multiplier (after enabling/timing/
                // distribution are set so the setter guard sees the final state).
                if (md.has("firingRateDependence")) {
                    SerializableFunction<Matrix, Double> g =
                        firingDepTableToHandle(md.getAsJsonObject("firingRateDependence"), nodeMap, classMap);
                    if (g != null) {
                        tnode.setFiringRateDependence(mode, g);
                    }
                }
            }
        }

        // Phase 3c: Restore the initial state of EVERY stateful node the document
        // names. A Place's token counts and a closed pass-and-swap queue's ordered
        // job placement are two readings of one field, and `hasInitState` tests
        // every stateful node, so restoring only those two left the model reading
        // as uninitialized: getState() then ran initDefault() and overwrote the
        // state, state space and prior this reader had just installed.
        for (Map.Entry<String, JsonObject> njEntry : nodeJsonMap.entrySet()) {
            JsonObject nodeObj = njEntry.getValue();
            if (!nodeObj.has("initialState")) continue;
            Node node = nodeMap.get(njEntry.getKey());
            if (!(node instanceof StatefulNode)) continue;

            JsonElement isElem = nodeObj.get("initialState");
            Matrix stateVec;
            if (isElem.isJsonArray()) {
                JsonArray isArr = isElem.getAsJsonArray();
                stateVec = new Matrix(1, isArr.size());
                for (int k = 0; k < isArr.size(); k++) {
                    stateVec.set(0, k, isArr.get(k).getAsDouble());
                }
            } else {
                stateVec = Matrix.singleton(isElem.getAsDouble());
            }
            StatefulNode sfNode = (StatefulNode) node;
            sfNode.setState(stateVec);
            // AND ITS ONE-ROW STATE SPACE AND TRIVIAL PRIOR. The writer omits a
            // [1] prior over one row because `initialState` already carries that
            // row, so restoring the row alone leaves a node whose state, space
            // and prior disagree -- which is not what initDefault or any
            // initFromMarginal* builds: all of them set the TRIO together, and a
            // solver that indexes the state space finds it empty. A space or
            // prior the document DID carry is installed by the loop below and
            // overwrites this pair.
            if (sfNode.getStateSpace().isEmpty()) {
                sfNode.setStateSpace(stateVec);
                sfNode.setStatePrior(Matrix.singleton(1.0));
            }
        }

        // see _kb/09-ldes-and-cache.md (statePrior is emitted only as a pair with stateSpace)
        for (Map.Entry<String, JsonObject> njEntry : nodeJsonMap.entrySet()) {
            JsonObject nodeObj = njEntry.getValue();
            if (!nodeObj.has("statePrior")) continue;
            Node node = nodeMap.get(njEntry.getKey());
            if (!(node instanceof StatefulNode)) continue;
            StatefulNode sfNode = (StatefulNode) node;
            if (nodeObj.has("stateSpace")) {
                sfNode.setStateSpace(jsonToMatrix2D(nodeObj.getAsJsonArray("stateSpace")));
            }
            JsonArray priorArr = nodeObj.getAsJsonArray("statePrior");
            // The prior indexes the rows of the state space, so it is a column.
            Matrix prior = new Matrix(priorArr.size(), 1);
            for (int k = 0; k < priorArr.size(); k++) {
                prior.set(k, 0, priorArr.get(k).getAsDouble());
            }
            sfNode.setStatePrior(prior);
        }

        // Phase 4: Build routing matrix and link
        if (routingObj != null && routingObj.size() > 0) {
            List<Node> nodeList = model.getNodes();
            RoutingMatrix P = new RoutingMatrix(model, jobClasses, nodeList);

            // Build reverse index maps
            Map<String, Node> nodeNameMap = new HashMap<String, Node>();
            for (Node n : nodeList) {
                nodeNameMap.put(n.getName(), n);
            }
            Map<String, JobClass> classNameMap = new HashMap<String, JobClass>();
            for (JobClass jc : jobClasses) {
                classNameMap.put(jc.getName(), jc);
            }

            for (Map.Entry<String, JsonElement> rtEntry : routingObj.entrySet()) {
                String key = rtEntry.getKey();
                String[] parts = key.split(",", 2);
                if (parts.length != 2) {
                    continue;
                }
                JobClass originClass = classNameMap.get(parts[0].trim());
                JobClass targetClass = classNameMap.get(parts[1].trim());
                if (originClass == null || targetClass == null) {
                    continue;
                }

                JsonObject fromToObj = rtEntry.getValue().getAsJsonObject();
                for (Map.Entry<String, JsonElement> fromEntry : fromToObj.entrySet()) {
                    Node srcNode = nodeNameMap.get(fromEntry.getKey());
                    if (srcNode == null) {
                        continue;
                    }
                    JsonObject destsObj = fromEntry.getValue().getAsJsonObject();
                    for (Map.Entry<String, JsonElement> destEntry : destsObj.entrySet()) {
                        Node destNode = nodeNameMap.get(destEntry.getKey());
                        if (destNode == null) {
                            continue;
                        }
                        double prob = destEntry.getValue().getAsDouble();
                        P.addConnection(srcNode, destNode, originClass, targetClass, prob);
                    }
                }
            }

            model.link(P);
        }

        // Phase 5: Restore non-PROB routing strategies and WRROBIN weights
        if (routingRaw != null) {
            JsonObject rsObj = null;
            JsonObject rwObj = null;
            JsonObject rpObj = null;
            // routingStrategies/routingWeights may be inside the routing object or at model level
            if (routingRaw.has("routingStrategies")) {
                rsObj = routingRaw.getAsJsonObject("routingStrategies");
            }
            if (routingRaw.has("routingWeights")) {
                rwObj = routingRaw.getAsJsonObject("routingWeights");
            }
            if (routingRaw.has("routingParams")) {
                rpObj = routingRaw.getAsJsonObject("routingParams");
            }
            // Also check model level (MATLAB/Python save format)
            if (rsObj == null && modelObj.has("routingStrategies")) {
                rsObj = modelObj.getAsJsonObject("routingStrategies");
            }
            if (rwObj == null && modelObj.has("routingWeights")) {
                rwObj = modelObj.getAsJsonObject("routingWeights");
            }
            if (rpObj == null && modelObj.has("routingParams")) {
                rpObj = modelObj.getAsJsonObject("routingParams");
            }

            if (rsObj != null) {
                for (Map.Entry<String, JsonElement> nodeEntry : rsObj.entrySet()) {
                    Node node = nodeMap.get(nodeEntry.getKey());
                    if (node == null) continue;
                    JsonObject classStrats = nodeEntry.getValue().getAsJsonObject();
                    for (Map.Entry<String, JsonElement> classEntry : classStrats.entrySet()) {
                        JobClass jc = classMap.get(classEntry.getKey());
                        if (jc == null) continue;
                        String stratName = classEntry.getValue().getAsString();
                        RoutingStrategy rs = parseRoutingStrategy(stratName);
                        if (rs == null) {
                            line_warning(mfilename(new Object() {
                        }),
                                    "Unrecognized routing strategy '%s' at node %s for class %s; "
                                            + "the strategy is not restored.",
                                    stratName, node.getName(), jc.getName());
                            continue;
                        }
                        // see _kb/09-ldes-and-cache.md for the routing-strategy apply-ordering rationale
                        // RAND is restored here, after link(P): it takes destinations from the connection matrix, so the single null-dest entry setRouting leaves is correct.
                        if (rs != RoutingStrategy.WRROBIN
                                && rs != RoutingStrategy.PROB
                                && rs != RoutingStrategy.SQ) {
                            node.setRouting(jc, rs);
                        }
                    }
                }
                // SQ: applied after the plain strategies so its parameter and
                // the strategy are installed by one call.
                for (Map.Entry<String, JsonElement> nodeEntry : rsObj.entrySet()) {
                    Node node = nodeMap.get(nodeEntry.getKey());
                    if (node == null) continue;
                    JsonObject classStrats = nodeEntry.getValue().getAsJsonObject();
                    for (Map.Entry<String, JsonElement> classEntry : classStrats.entrySet()) {
                        JobClass jc = classMap.get(classEntry.getKey());
                        if (jc == null) continue;
                        RoutingStrategy rs = parseRoutingStrategy(classEntry.getValue().getAsString());
                        if (rs != RoutingStrategy.SQ) {
                            continue;
                        }
                        JsonObject params = null;
                        if (rpObj != null && rpObj.has(nodeEntry.getKey())) {
                            JsonObject nodeParams = rpObj.getAsJsonObject(nodeEntry.getKey());
                            if (nodeParams.has(classEntry.getKey())) {
                                params = nodeParams.getAsJsonObject(classEntry.getKey());
                            }
                        }
                        if (rs == RoutingStrategy.SQ) {
                            if (params == null || !params.has("d")) {
                                line_warning(mfilename(new Object() {
                        }),
                                        "Node %s routes class %s by SQ but the file carries no "
                                                + "routingParams.d; falling back to d=2.",
                                        node.getName(), jc.getName());
                                node.setSQRouting(jc, 2);
                            } else {
                                node.setSQRouting(jc, params.get("d").getAsInt());
                            }
                        }
                    }
                }
            }

            if (rwObj != null) {
                for (Map.Entry<String, JsonElement> nodeEntry : rwObj.entrySet()) {
                    Node node = nodeMap.get(nodeEntry.getKey());
                    if (node == null) continue;
                    JsonObject classWeights = nodeEntry.getValue().getAsJsonObject();
                    for (Map.Entry<String, JsonElement> classEntry : classWeights.entrySet()) {
                        JobClass jc = classMap.get(classEntry.getKey());
                        if (jc == null) continue;
                        JsonObject destWeights = classEntry.getValue().getAsJsonObject();
                        for (Map.Entry<String, JsonElement> destEntry : destWeights.entrySet()) {
                            Node dest = nodeMap.get(destEntry.getKey());
                            if (dest == null) continue;
                            double weight = destEntry.getValue().getAsDouble();
                            node.setRouting(jc, RoutingStrategy.WRROBIN, dest, weight);
                        }
                    }
                }
            }
        }

        // Krzesinski state-dependent routing, restored after link(P) and after the
        // non-PROB strategies above: link writes the uniform placeholder into the
        // entry row and routingStrategies names that row SDR, neither of which says
        // what the subnetwork looks like. See _kb/16-state-dependent-routing.md
        if (modelObj.has("stateDepRouting")) {
            restoreStateDepRouting(modelObj.getAsJsonObject("stateDepRouting"), nodeMap, classMap);
        }

        // Global (Whittle) dependence phi(n), rebuilt from the materialized slot lattice
        if (modelObj.has("globalDependence")) {
            restoreGlobalDependence(model, modelObj.getAsJsonObject("globalDependence"),
                    nodeMap, classMap);
        }

        // Finite capacity regions
        if (modelObj.has("finiteCapacityRegions")) {
            JsonArray fcrArr = modelObj.getAsJsonArray("finiteCapacityRegions");
            for (JsonElement fcrEl : fcrArr) {
                JsonObject rj = fcrEl.getAsJsonObject();
                // Collect region nodes from "stations" array or legacy "nodes" array
                List<Node> regionNodes = new ArrayList<Node>();
                if (rj.has("stations")) {
                    for (JsonElement stEl : rj.getAsJsonArray("stations")) {
                        JsonObject sj = stEl.getAsJsonObject();
                        String nodeName = sj.get("node").getAsString();
                        Node n = nodeMap.get(nodeName);
                        if (n != null) regionNodes.add(n);
                    }
                } else if (rj.has("nodes")) {
                    for (JsonElement nEl : rj.getAsJsonArray("nodes")) {
                        String nodeName = nEl.getAsString();
                        Node n = nodeMap.get(nodeName);
                        if (n != null) regionNodes.add(n);
                    }
                }
                if (!regionNodes.isEmpty()) {
                    Region region = model.addRegion(regionNodes);
                    if (rj.has("name")) {
                        region.setName(rj.get("name").getAsString());
                    }
                    if (rj.has("globalMaxJobs")) {
                        region.setGlobalMaxJobs(rj.get("globalMaxJobs").getAsInt());
                    }
                    if (rj.has("globalMaxMemory")) {
                        region.setGlobalMaxMemory(rj.get("globalMaxMemory").getAsDouble());
                    }
                    if (rj.has("classMaxJobs")) {
                        JsonObject cmjObj = rj.getAsJsonObject("classMaxJobs");
                        for (Map.Entry<String, JsonElement> entry : cmjObj.entrySet()) {
                            JobClass jc = classMap.get(entry.getKey());
                            if (jc != null) {
                                region.setClassMaxJobs(jc, entry.getValue().getAsInt());
                            }
                        }
                    }
                    if (rj.has("classMaxMemory")) {
                        JsonObject cmmObj = rj.getAsJsonObject("classMaxMemory");
                        for (Map.Entry<String, JsonElement> entry : cmmObj.entrySet()) {
                            JobClass jc = classMap.get(entry.getKey());
                            if (jc != null) {
                                region.setClassMaxMemory(jc, entry.getValue().getAsInt());
                            }
                        }
                    }
                    if (rj.has("dropRule")) {
                        JsonObject drObj = rj.getAsJsonObject("dropRule");
                        for (Map.Entry<String, JsonElement> entry : drObj.entrySet()) {
                            JobClass jc = classMap.get(entry.getKey());
                            if (jc != null) {
                                region.setDropRule(jc, parseDropStrategy(entry.getValue().getAsString()));
                            }
                        }
                    }
                    // General linear admission constraints A*n <= b
                    if (rj.has("constraintA") && rj.has("constraintB")) {
                        Matrix A = jsonToMatrix2DLenient(rj.get("constraintA"));
                        Matrix bFlat = jsonToMatrix2DLenient(rj.get("constraintB"));
                        // b is a capacity vector, one entry per constraint row
                        Matrix b = new Matrix(bFlat.getNumElements(), 1);
                        for (int k = 0; k < bFlat.getNumElements(); k++) {
                            b.set(k, 0, bFlat.get(k));
                        }
                        region.setLinearConstraints(A, b);
                    }
                    // Per-station classWeight and classSize
                    if (rj.has("stations")) {
                        for (JsonElement stEl : rj.getAsJsonArray("stations")) {
                            JsonObject sj = stEl.getAsJsonObject();
                            if (sj.has("classWeight")) {
                                JsonObject cwObj = sj.getAsJsonObject("classWeight");
                                for (Map.Entry<String, JsonElement> entry : cwObj.entrySet()) {
                                    JobClass jc = classMap.get(entry.getKey());
                                    if (jc != null) {
                                        region.setClassWeight(jc, entry.getValue().getAsDouble());
                                    }
                                }
                            }
                            if (sj.has("classSize")) {
                                JsonObject csObj = sj.getAsJsonObject("classSize");
                                for (Map.Entry<String, JsonElement> entry : csObj.entrySet()) {
                                    JobClass jc = classMap.get(entry.getKey());
                                    if (jc != null) {
                                        region.setClassSize(jc, entry.getValue().getAsDouble());
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        loadNetworkRewards(model, modelObj, nodeMap, classMap);

        return model;
    }

    // ========================================================================
    // LAYERED NETWORK LOAD
    // ========================================================================

    private static LayeredNetwork loadLayeredNetwork(JsonObject modelObj) throws IOException {
        String name = modelObj.has("name") ? modelObj.get("name").getAsString() : "model";
        LayeredNetwork model = new LayeredNetwork(name);

        // Create hosts (using Processor, which extends Host, so Task.on() accepts them)
        Map<String, Processor> hostMap = new LinkedHashMap<String, Processor>();
        if (modelObj.has("hosts")) {
            JsonArray hostsArr = modelObj.getAsJsonArray("hosts");
            for (JsonElement hostEl : hostsArr) {
                JsonObject hostObj = hostEl.getAsJsonObject();
                String hostName = hostObj.get("name").getAsString();
                int mult = hostObj.has("multiplicity") ? hostObj.get("multiplicity").getAsInt() : 1;
                String schedStr = hostObj.has("scheduling") ? hostObj.get("scheduling").getAsString() : "PS";
                SchedStrategy sched = parseSchedStrategy(schedStr);
                double quantum = hostObj.has("quantum") ? hostObj.get("quantum").getAsDouble() : 0.01;
                double speedFactor = hostObj.has("speedFactor") ? hostObj.get("speedFactor").getAsDouble() : 1.0;

                Processor host = new Processor(model, hostName, mult, sched, quantum, speedFactor);
                if (hostObj.has("replication")) {
                    host.setReplication(hostObj.get("replication").getAsInt());
                }
                if (hostObj.has("admissionConstraints")) {
                    applyLincon(host, hostObj.getAsJsonArray("admissionConstraints"));
                }
                hostMap.put(hostName, host);
            }
        }

        // Create tasks
        Map<String, Task> taskMap = new LinkedHashMap<String, Task>();
        if (modelObj.has("tasks")) {
            JsonArray tasksArr = modelObj.getAsJsonArray("tasks");
            for (JsonElement taskEl : tasksArr) {
                JsonObject taskObj = taskEl.getAsJsonObject();
                String taskName = taskObj.get("name").getAsString();
                int mult = taskObj.has("multiplicity") ? taskObj.get("multiplicity").getAsInt() : 1;
                String schedStr = taskObj.has("scheduling") ? taskObj.get("scheduling").getAsString() : "INF";
                SchedStrategy sched = parseSchedStrategy(schedStr);

                Task task;
                String taskType = taskObj.has("taskType") ? taskObj.get("taskType").getAsString() : "";
                if ("CacheTask".equals(taskType)) {
                    int totalItems = taskObj.has("totalItems") ? taskObj.get("totalItems").getAsInt() : 1;
                    // MATLAB and Python write the per-level capacity verbatim, so a
                    // single-level cache arrives as a scalar rather than a 1-element array
                    int[] itemCap = jsonToIntArray(taskObj, "cacheCapacity");
                    if (itemCap.length == 0) {
                        itemCap = new int[]{1};
                    }
                    String replStr = taskObj.has("replacementStrategy") ? taskObj.get("replacementStrategy").getAsString() : "FIFO";
                    ReplacementStrategy repl = parseReplacementStrategy(replStr);
                    task = new CacheTask(model, taskName, totalItems, itemCap, repl, mult, sched);
                } else if ("SetupTask".equals(taskType) || "FunctionTask".equals(taskType)) {
                    // FunctionTask is the legacy name of SetupTask on the wire
                    task = new SetupTask(model, taskName, mult, sched);
                } else {
                    task = new Task(model, taskName, mult, sched);
                }
                if (taskObj.has("host")) {
                    String hostName = taskObj.get("host").getAsString();
                    Processor host = hostMap.get(hostName);
                    if (host != null) {
                        task.on(host);
                    }
                }
                if (taskObj.has("replication")) {
                    task.setReplication(taskObj.get("replication").getAsInt());
                }
                // Both spellings are accepted; the object form wins because it
                // carries the whole distribution, where the mean alone can only
                // be rebuilt as an exponential.
                if (taskObj.has("thinkTime") && taskObj.get("thinkTime").isJsonObject()) {
                    task.setThinkTime(deserializeDistribution(taskObj.getAsJsonObject("thinkTime")));
                } else if (taskObj.has("thinkTimeMean")) {
                    double thinkMean = taskObj.get("thinkTimeMean").getAsDouble();
                    task.setThinkTime(thinkMean);
                }
                if (taskObj.has("priority")) {
                    task.setPriority(taskObj.get("priority").getAsInt());
                }
                if (taskObj.has("fanIn")) {
                    JsonObject fanInObj = taskObj.getAsJsonObject("fanIn");
                    for (Map.Entry<String, JsonElement> fi : fanInObj.entrySet()) {
                        task.setFanIn(fi.getKey(), fi.getValue().getAsInt());
                    }
                }
                if (taskObj.has("fanOut")) {
                    JsonObject fanOutObj = taskObj.getAsJsonObject("fanOut");
                    for (Map.Entry<String, JsonElement> fo : fanOutObj.entrySet()) {
                        task.setFanOut(fo.getKey(), fo.getValue().getAsInt());
                    }
                }
                if (taskObj.has("setupTime") && taskObj.get("setupTime").isJsonObject()) {
                    task.setSetupTime(deserializeDistribution(taskObj.getAsJsonObject("setupTime")));
                } else if (taskObj.has("setupTimeMean")) {
                    double setupMean = taskObj.get("setupTimeMean").getAsDouble();
                    task.setSetupTime(setupMean);
                }
                if (taskObj.has("delayOffTime") && taskObj.get("delayOffTime").isJsonObject()) {
                    task.setDelayOffTime(deserializeDistribution(taskObj.getAsJsonObject("delayOffTime")));
                } else if (taskObj.has("delayOffTimeMean")) {
                    double delayOffMean = taskObj.get("delayOffTimeMean").getAsDouble();
                    task.setDelayOffTime(delayOffMean);
                }
                if (taskObj.has("admissionConstraints")) {
                    applyLincon(task, taskObj.getAsJsonArray("admissionConstraints"));
                }
                taskMap.put(taskName, task);
            }
        }

        // Create entries
        Map<String, Entry> entryMap = new LinkedHashMap<String, Entry>();
        if (modelObj.has("entries")) {
            JsonArray entriesArr = modelObj.getAsJsonArray("entries");
            for (JsonElement entryEl : entriesArr) {
                JsonObject entryObj = entryEl.getAsJsonObject();
                String entryName = entryObj.get("name").getAsString();
                Entry entry;
                boolean isItemEntry = entryObj.has("entryType") && "ItemEntry".equals(entryObj.get("entryType").getAsString());
                if (isItemEntry) {
                    int cardinality = entryObj.has("totalItems") ? entryObj.get("totalItems").getAsInt() : 1;
                    Distribution popDist = Immediate.getInstance();
                    if (entryObj.has("accessProb")) {
                        popDist = deserializeDistribution(entryObj.getAsJsonObject("accessProb"));
                    }
                    entry = new ItemEntry(model, entryName, cardinality, popDist);
                } else {
                    entry = new Entry(model, entryName);
                }
                if (entryObj.has("task")) {
                    String taskName = entryObj.get("task").getAsString();
                    Task task = taskMap.get(taskName);
                    if (task != null) {
                        entry.on(task);
                    }
                }
                if (entryObj.has("arrival")) {
                    Distribution arrDist = deserializeDistribution(entryObj.getAsJsonObject("arrival"));
                    entry.setArrival(arrDist);
                }
                if (entryObj.has("forwarding")) {
                    // Defer forwarding until all entries are created
                }
                entryMap.put(entryName, entry);
            }

            // Second pass: set forwarding
            for (JsonElement entryEl : entriesArr) {
                JsonObject entryObj = entryEl.getAsJsonObject();
                if (entryObj.has("forwarding")) {
                    String entryName = entryObj.get("name").getAsString();
                    Entry entry = entryMap.get(entryName);
                    if (entry != null) {
                        JsonArray fwArr = entryObj.getAsJsonArray("forwarding");
                        for (JsonElement fwEl : fwArr) {
                            JsonObject fwObj = fwEl.getAsJsonObject();
                            String destName = fwObj.get("dest").getAsString();
                            double prob = fwObj.has("prob") ? fwObj.get("prob").getAsDouble() : 1.0;
                            Entry destEntry = entryMap.get(destName);
                            if (destEntry != null) {
                                entry.forward(destEntry, prob);
                            } else {
                                entry.forward(destName, prob);
                            }
                        }
                    }
                }
            }
        }

        // Create activities
        Map<String, Activity> activityMap = new LinkedHashMap<String, Activity>();
        if (modelObj.has("activities")) {
            JsonArray actsArr = modelObj.getAsJsonArray("activities");
            for (JsonElement actEl : actsArr) {
                JsonObject actObj = actEl.getAsJsonObject();
                String actName = actObj.get("name").getAsString();

                Distribution hostDemand = Immediate.getInstance();
                if (actObj.has("hostDemand")) {
                    hostDemand = deserializeDistribution(actObj.getAsJsonObject("hostDemand"));
                }

                Activity activity = new Activity(model, actName, hostDemand);

                if (actObj.has("task")) {
                    String taskName = actObj.get("task").getAsString();
                    Task task = taskMap.get(taskName);
                    if (task != null) {
                        activity.on(task);
                    }
                }
                if (actObj.has("boundToEntry")) {
                    String entryName = actObj.get("boundToEntry").getAsString();
                    Entry entry = entryMap.get(entryName);
                    if (entry != null) {
                        activity.boundTo(entry);
                    } else {
                        activity.boundTo(entryName);
                    }
                }
                if (actObj.has("repliesTo")) {
                    String entryName = actObj.get("repliesTo").getAsString();
                    Entry entry = entryMap.get(entryName);
                    if (entry != null) {
                        activity.repliesTo(entry);
                    }
                }
                if (actObj.has("thinkTime")) {
                    Distribution thinkDist = deserializeDistribution(actObj.getAsJsonObject("thinkTime"));
                    activity.setThinkTime(thinkDist);
                }
                if (actObj.has("callOrder")) {
                    activity.setCallOrder(actObj.get("callOrder").getAsString());
                }
                if (actObj.has("synchCalls")) {
                    JsonArray syncArr = actObj.getAsJsonArray("synchCalls");
                    for (JsonElement callEl : syncArr) {
                        JsonObject callObj = callEl.getAsJsonObject();
                        String destName = callObj.get("dest").getAsString();
                        double mean = callObj.has("mean") ? callObj.get("mean").getAsDouble() : 1.0;
                        Entry destEntry = entryMap.get(destName);
                        if (destEntry != null) {
                            activity.synchCall(destEntry, mean);
                        } else {
                            activity.synchCall(destName, mean);
                        }
                    }
                }
                if (actObj.has("asynchCalls")) {
                    JsonArray asyncArr = actObj.getAsJsonArray("asynchCalls");
                    for (JsonElement callEl : asyncArr) {
                        JsonObject callObj = callEl.getAsJsonObject();
                        String destName = callObj.get("dest").getAsString();
                        double mean = callObj.has("mean") ? callObj.get("mean").getAsDouble() : 1.0;
                        Entry destEntry = entryMap.get(destName);
                        if (destEntry != null) {
                            activity.asynchCall(destEntry, mean);
                        } else {
                            activity.asynchCall(destName, mean);
                        }
                    }
                }
                activityMap.put(actName, activity);
            }
        }

        // see _kb/09-ldes-and-cache.md for the LQN precedence dual-wire-schema rationale
        if (modelObj.has("precedences")) {
            JsonArray precsArr = modelObj.getAsJsonArray("precedences");
            for (JsonElement precEl : precsArr) {
                JsonObject precObj = precEl.getAsJsonObject();
                if (!precObj.has("task")) {
                    continue;
                }
                Task task = taskMap.get(precObj.get("task").getAsString());
                if (task == null) {
                    continue;
                }
                if (precObj.has("preActs") || precObj.has("postActs")) {
                    loadPrecedenceExplicit(precObj, task);
                } else if (precObj.has("type")) {
                    loadPrecedenceTyped(precObj, task, activityMap);
                }
            }
        }

        return model;
    }

    /**
     * Load a precedence given in the explicit
     * "preActs"/"postActs"/"preType"/"postType" schema, which carries the
     * precedence types verbatim and so needs no shape inference.
     */
    private static void loadPrecedenceExplicit(JsonObject precObj, Task task) {
        List<String> preActs = jsonToStringList(precObj, "preActs");
        List<String> postActs = jsonToStringList(precObj, "postActs");

        String preType = precObj.has("preType") ? precObj.get("preType").getAsString() : ActivityPrecedenceType.PRE_SEQ;
        String postType = precObj.has("postType") ? precObj.get("postType").getAsString() : ActivityPrecedenceType.POST_SEQ;

        Matrix preParams = null;
        if (precObj.has("preParams")) {
            preParams = jsonToRowVector(precObj.getAsJsonArray("preParams"));
        }
        // Default quorum for AND-join: all predecessors required
        if (preParams == null && preType.equals(ActivityPrecedenceType.PRE_AND)) {
            preParams = Matrix.ones(1, preActs.size());
        }
        Matrix postParams = null;
        if (precObj.has("postParams")) {
            postParams = jsonToRowVector(precObj.getAsJsonArray("postParams"));
        }

        task.addPrecedence(new ActivityPrecedence(preActs, postActs, preType, postType, preParams, postParams));
    }

    /**
     * Load a precedence given in the canonical "type"/"activities" schema
     * written by MATLAB and Python. The pre/post split is implied by the type:
     * a fork names its predecessor first, a join names its successor last.
     */
    private static void loadPrecedenceTyped(JsonObject precObj, Task task, Map<String, Activity> activityMap) {
        String type = precObj.get("type").getAsString();

        List<String> acts = new ArrayList<String>();
        for (String name : jsonToStringList(precObj, "activities")) {
            if (activityMap.containsKey(name)) {
                acts.add(name);
            }
        }
        if (acts.size() < 2 && !"Loop".equals(type)) {
            return;
        }
        List<String> tail = new ArrayList<String>(acts.subList(1, acts.size()));
        List<String> head = new ArrayList<String>(acts.subList(0, acts.size() - 1));
        String first = acts.isEmpty() ? null : acts.get(0);
        String last = acts.isEmpty() ? null : acts.get(acts.size() - 1);

        if ("Serial".equals(type)) {
            task.addPrecedence(ActivityPrecedence.Serial(acts));
        } else if ("AndFork".equals(type)) {
            task.addPrecedence(ActivityPrecedence.AndFork(first, tail));
        } else if ("AndJoin".equals(type)) {
            task.addPrecedence(ActivityPrecedence.AndJoin(head, last));
        } else if ("OrFork".equals(type)) {
            Matrix probs;
            if (precObj.has("probabilities")) {
                probs = jsonToRowVector(precObj.getAsJsonArray("probabilities"));
            } else {
                // Unspecified branch probabilities default to uniform
                probs = Matrix.ones(1, tail.size()).scale(1.0 / tail.size());
            }
            task.addPrecedence(ActivityPrecedence.OrFork(first, tail, probs));
        } else if ("OrJoin".equals(type)) {
            task.addPrecedence(ActivityPrecedence.OrJoin(head, last));
        } else if ("Loop".equals(type)) {
            double count = precObj.has("loopCount") ? precObj.get("loopCount").getAsDouble() : 1.0;
            String preName = precObj.has("preActivity") ? precObj.get("preActivity").getAsString() : null;
            if (preName != null && activityMap.containsKey(preName)) {
                // The loop body is carried whole in "activities"
                task.addPrecedence(ActivityPrecedence.Loop(preName, acts, count));
            } else if (acts.size() >= 2) {
                // Legacy form: the first activity is the loop trigger
                task.addPrecedence(ActivityPrecedence.Loop(first, tail, count));
            }
        } else if ("CacheAccess".equals(type)) {
            task.addPrecedence(ActivityPrecedence.CacheAccess(first, tail));
        }
    }

    /** Read a JSON int array, tolerating an absent key or a bare number. */
    private static int[] jsonToIntArray(JsonObject obj, String key) {
        if (!obj.has(key) || obj.get(key).isJsonNull()) {
            return new int[0];
        }
        JsonElement el = obj.get(key);
        if (!el.isJsonArray()) {
            return new int[]{el.getAsInt()};
        }
        JsonArray arr = el.getAsJsonArray();
        int[] out = new int[arr.size()];
        for (int i = 0; i < arr.size(); i++) {
            out[i] = arr.get(i).getAsInt();
        }
        return out;
    }

    /** Read a JSON string array, tolerating an absent key or a bare string. */
    private static List<String> jsonToStringList(JsonObject obj, String key) {
        List<String> out = new ArrayList<String>();
        if (!obj.has(key) || obj.get(key).isJsonNull()) {
            return out;
        }
        JsonElement el = obj.get(key);
        if (el.isJsonArray()) {
            for (JsonElement e : el.getAsJsonArray()) {
                out.add(e.getAsString());
            }
        } else {
            out.add(el.getAsString());
        }
        return out;
    }

    // ========================================================================
    // WORKFLOW LOAD
    // ========================================================================

    private static Workflow loadWorkflow(JsonObject modelObj) throws IOException {
        String name = modelObj.has("name") ? modelObj.get("name").getAsString() : "workflow";
        Workflow wf = new Workflow(name);

        if (modelObj.has("activities")) {
            JsonArray actsArr = modelObj.getAsJsonArray("activities");
            for (JsonElement actEl : actsArr) {
                JsonObject actObj = actEl.getAsJsonObject();
                String actName = actObj.get("name").getAsString();
                if (actObj.has("hostDemand")) {
                    Distribution hostDemand = deserializeDistribution(actObj.getAsJsonObject("hostDemand"));
                    wf.addActivity(actName, hostDemand);
                } else {
                    wf.addActivity(actName, 1.0);
                }
            }
        }

        if (modelObj.has("precedences")) {
            JsonArray precsArr = modelObj.getAsJsonArray("precedences");
            for (JsonElement precEl : precsArr) {
                JsonObject precObj = precEl.getAsJsonObject();
                List<String> preActs = new ArrayList<String>();
                for (JsonElement e : precObj.getAsJsonArray("preActs")) {
                    preActs.add(e.getAsString());
                }
                List<String> postActs = new ArrayList<String>();
                for (JsonElement e : precObj.getAsJsonArray("postActs")) {
                    postActs.add(e.getAsString());
                }
                String preType = precObj.has("preType") ? precObj.get("preType").getAsString() : ActivityPrecedenceType.PRE_SEQ;
                String postType = precObj.has("postType") ? precObj.get("postType").getAsString() : ActivityPrecedenceType.POST_SEQ;
                Matrix preParams = null;
                if (precObj.has("preParams")) {
                    preParams = jsonToRowVector(precObj.getAsJsonArray("preParams"));
                }
                Matrix postParams = null;
                if (precObj.has("postParams")) {
                    postParams = jsonToRowVector(precObj.getAsJsonArray("postParams"));
                }
                wf.addPrecedence(new ActivityPrecedence(preActs, postActs, preType, postType, preParams, postParams));
            }
        }

        return wf;
    }

    // ========================================================================
    // ENVIRONMENT LOAD
    // ========================================================================

    private static Environment loadEnvironment(JsonObject modelObj) throws IOException {
        String name = modelObj.has("name") ? modelObj.get("name").getAsString() : "env";
        int numStages = modelObj.has("numStages") ? modelObj.get("numStages").getAsInt() : 0;

        // see _kb/09-ldes-and-cache.md (nodeFailures has two roles on load)
        JsonArray nfArr = modelObj.has("nodeFailures") ? modelObj.getAsJsonArray("nodeFailures") : new JsonArray();
        JsonArray stagesArr = modelObj.has("stages") ? modelObj.getAsJsonArray("stages") : new JsonArray();

        List<String> declaredNames = new ArrayList<String>();
        for (int i = 0; i < stagesArr.size(); i++) {
            JsonObject stageObj = stagesArr.get(i).getAsJsonObject();
            declaredNames.add(stageObj.has("name") ? stageObj.get("name").getAsString() : "Stage" + i);
        }

        boolean macroMode = nfArr.size() > 0;
        for (int k = 0; k < nfArr.size(); k++) {
            JsonObject nf = nfArr.get(k).getAsJsonObject();
            if (!nf.has("node")) {
                throw new RuntimeException("A \"nodeFailures\" entry is missing the required \"node\" field.");
            }
            if (declaredNames.contains("DOWN_" + nf.get("node").getAsString())) {
                macroMode = false;
                break;
            }
        }

        // see _kb/09-ldes-and-cache.md for the Environment stage pre-allocation rationale
        if (macroMode) {
            numStages = Math.max(numStages, 1 + nfArr.size());
        }
        Environment env = new Environment(name, numStages);

        if (macroMode) {
            if (stagesArr.size() != 1) {
                throw new RuntimeException("\"nodeFailures\" expands the base model into the UP and DOWN_<node> "
                        + "stages, so \"stages\" must declare exactly one stage, holding the base (UP) model.");
            }
            if (modelObj.has("transitions") && modelObj.getAsJsonArray("transitions").size() > 0) {
                throw new RuntimeException("\"nodeFailures\" implies the breakdown and repair transitions; "
                        + "\"transitions\" must not be declared alongside it.");
            }
            JsonObject baseStage = stagesArr.get(0).getAsJsonObject();
            if (!baseStage.has("model")) {
                throw new RuntimeException("\"nodeFailures\" requires the base stage to carry a \"model\".");
            }
            Network baseModel = loadNetwork(baseStage.getAsJsonObject("model"));
            for (int k = 0; k < nfArr.size(); k++) {
                JsonObject nf = nfArr.get(k).getAsJsonObject();
                String nodeName = nf.get("node").getAsString();
                Markovian breakdown = nodeFailureDist(nf, "breakdownRate", nodeName, true);
                Markovian downService = nodeFailureDist(nf, "downService", nodeName, true);
                Markovian repair = nodeFailureDist(nf, "repairRate", nodeName, false);
                String resetB = nodeFailureResetPolicy(nf, "breakdownResetPolicy");
                String resetR = nodeFailureResetPolicy(nf, "repairResetPolicy");
                if (repair == null) {
                    env.addNodeBreakdown(baseModel, nodeName, breakdown, downService, resetB);
                } else {
                    env.addNodeFailureRepair(baseModel, nodeName, breakdown, repair, downService, resetB, resetR);
                }
            }
        } else {
            for (int i = 0; i < stagesArr.size(); i++) {
                JsonObject stageObj = stagesArr.get(i).getAsJsonObject();
                String stageName = declaredNames.get(i);
                String stageType = stageObj.has("type") ? stageObj.get("type").getAsString() : "";
                Network stageModel = null;
                if (stageObj.has("model")) {
                    stageModel = loadNetwork(stageObj.getAsJsonObject("model"));
                }
                if (stageModel != null) {
                    env.addStage(i, stageName, stageType, stageModel);
                }
            }

            if (modelObj.has("transitions")) {
                JsonArray transArr = modelObj.getAsJsonArray("transitions");
                for (JsonElement transEl : transArr) {
                    JsonObject transObj = transEl.getAsJsonObject();
                    int from = transObj.get("from").getAsInt();
                    int to = transObj.get("to").getAsInt();
                    Distribution dist = deserializeDistribution(transObj.getAsJsonObject("distribution"));
                    if (dist instanceof Markovian) {
                        env.addTransition(from, to, (Markovian) dist);
                    }
                }
            }

            // Re-attach the node-failure descriptors and their reset policies to the
            // stages just built, so that the environment serializes back identically.
            for (int k = 0; k < nfArr.size(); k++) {
                JsonObject nf = nfArr.get(k).getAsJsonObject();
                String nodeName = nf.get("node").getAsString();
                Markovian breakdown = nodeFailureDist(nf, "breakdownRate", nodeName, true);
                Markovian downService = nodeFailureDist(nf, "downService", nodeName, true);
                Markovian repair = nodeFailureDist(nf, "repairRate", nodeName, false);
                String resetB = nodeFailureResetPolicy(nf, "breakdownResetPolicy");
                String resetR = nodeFailureResetPolicy(nf, "repairResetPolicy");
                env.registerNodeFailure(nodeName, breakdown, repair, downService, resetB, resetR);
            }
        }

        env.init();
        return env;
    }

    /**
     * Decodes one distribution field of a "nodeFailures" entry. Note that
     * breakdownRate/repairRate carry full distributions, not scalar rates.
     */
    private static Markovian nodeFailureDist(JsonObject nf, String key, String nodeName, boolean required) {
        if (!nf.has(key) || nf.get(key).isJsonNull()) {
            if (required) {
                throw new RuntimeException("Node failure on \"" + nodeName + "\" is missing the required \""
                        + key + "\" field.");
            }
            return null;
        }
        Distribution dist = deserializeDistribution(nf.getAsJsonObject(key));
        if (!(dist instanceof Markovian)) {
            throw new RuntimeException("Node failure on \"" + nodeName + "\" has a \"" + key
                    + "\" that is not a Markovian distribution.");
        }
        return (Markovian) dist;
    }

    /**
     * Decodes one reset-policy field of a "nodeFailures" entry, defaulting to "keep".
     */
    private static String nodeFailureResetPolicy(JsonObject nf, String key) {
        if (!nf.has(key) || nf.get(key).isJsonNull()) {
            return Environment.RESET_POLICY_KEEP;
        }
        String value = nf.get(key).getAsString();
        if (value.isEmpty()) {
            return Environment.RESET_POLICY_KEEP;
        }
        return value;
    }

    // ========================================================================
    // MATRIX <-> JSON CONVERSION UTILITIES
    // ========================================================================

    /**
     * Converts a Matrix to a JSON array. If the matrix is a row vector (1 row),
     * returns a flat array. Otherwise returns a 2D array.
     */
    private static JsonArray matrixToJsonArray(Matrix m) {
        if (m == null) {
            return new JsonArray();
        }
        int rows = m.getNumRows();
        int cols = m.getNumCols();

        JsonArray arr = new JsonArray();
        if (rows == 1) {
            // Row vector: flat array
            for (int j = 0; j < cols; j++) {
                arr.add(m.get(0, j));
            }
        } else if (cols == 1) {
            // Column vector: flat array
            for (int i = 0; i < rows; i++) {
                arr.add(m.get(i, 0));
            }
        } else {
            // 2D matrix
            for (int i = 0; i < rows; i++) {
                JsonArray row = new JsonArray();
                for (int j = 0; j < cols; j++) {
                    row.add(m.get(i, j));
                }
                arr.add(row);
            }
        }
        return arr;
    }

    /**
     * Converts a Matrix to a 2D JSON array (always nested, even for 1x1).
     */
    private static JsonArray matrixToJson2D(Matrix m) {
        JsonArray arr = new JsonArray();
        if (m == null) {
            return arr;
        }
        int rows = m.getNumRows();
        int cols = m.getNumCols();
        for (int i = 0; i < rows; i++) {
            JsonArray row = new JsonArray();
            for (int j = 0; j < cols; j++) {
                row.add(m.get(i, j));
            }
            arr.add(row);
        }
        return arr;
    }

    /**
     * Converts a row-vector Matrix to a flat JSON array.
     */
    private static JsonArray matrixToJsonRowVector(Matrix m) {
        JsonArray arr = new JsonArray();
        if (m == null) {
            return arr;
        }
        // Treat as 1D regardless of shape
        int len = m.getNumCols();
        if (m.getNumRows() > 1 && m.getNumCols() == 1) {
            // Column vector
            len = m.getNumRows();
            for (int i = 0; i < len; i++) {
                arr.add(m.get(i, 0));
            }
        } else {
            for (int j = 0; j < len; j++) {
                arr.add(m.get(0, j));
            }
        }
        return arr;
    }

    /**
     * Converts a 2D JSON array to a Matrix.
     */
    private static Matrix jsonToMatrix2D(JsonArray arr) {
        if (arr == null || arr.size() == 0) {
            return new Matrix(0, 0);
        }
        int rows = arr.size();
        JsonArray firstRow = arr.get(0).getAsJsonArray();
        int cols = firstRow.size();
        Matrix m = new Matrix(rows, cols);
        for (int i = 0; i < rows; i++) {
            JsonArray row = arr.get(i).getAsJsonArray();
            for (int j = 0; j < cols; j++) {
                m.set(i, j, row.get(j).getAsDouble());
            }
        }
        return m;
    }

    /**
     * Lenient 2-D matrix parser: accepts a bare scalar (1x1), a 1-D array
     * (single row) or a 2-D array-of-arrays. The MATLAB linemodel_save
     * encoder emits 1x1 matrices as bare scalars.
     */
    /**
     * A block array that may sit under "params" (what linemodel_save.m and
     * linemodel_io.py write) or at the top level (what this class wrote before
     * 2026-08-01). BMAP and MarkedMMPP were unreachable across codebases while
     * the two disagreed: a MATLAB or Python model.json carrying a BMAP service
     * failed to load with "BMAP requires D0 and at least one batch matrix".
     */
    private static JsonArray blockArray(JsonObject obj, String key) {
        JsonObject params = obj.getAsJsonObject("params");
        if (params != null && params.getAsJsonArray(key) != null) {
            return params.getAsJsonArray(key);
        }
        return obj.getAsJsonArray(key);
    }

    private static Matrix jsonToMatrix2DLenient(JsonElement el) {
        if (el == null || el.isJsonNull()) {
            return new Matrix(0, 0);
        }
        if (el.isJsonPrimitive()) {
            Matrix m = new Matrix(1, 1);
            m.set(0, 0, el.getAsDouble());
            return m;
        }
        JsonArray arr = el.getAsJsonArray();
        if (arr.size() == 0) {
            return new Matrix(0, 0);
        }
        if (arr.get(0).isJsonPrimitive()) {
            Matrix m = new Matrix(1, arr.size());
            for (int j = 0; j < arr.size(); j++) {
                m.set(0, j, arr.get(j).getAsDouble());
            }
            return m;
        }
        return jsonToMatrix2D(arr);
    }

    /**
     * Converts a flat JSON array to a row-vector Matrix (1 x n).
     */
    private static Matrix jsonToRowVector(JsonArray arr) {
        if (arr == null || arr.size() == 0) {
            return new Matrix(1, 0);
        }
        // Check if this is a nested 2D array or flat
        if (arr.get(0).isJsonArray()) {
            // 2D: return first row as vector (or full matrix if multi-row)
            return jsonToMatrix2D(arr);
        }
        int n = arr.size();
        Matrix m = new Matrix(1, n);
        for (int j = 0; j < n; j++) {
            m.set(0, j, arr.get(j).getAsDouble());
        }
        return m;
    }

    // ========================================================================
    // ENUM PARSING UTILITIES
    // ========================================================================

    /**
     * Parses a scheduling strategy string to a {@link SchedStrategy} enum value.
     */
    private static SchedStrategy parseSchedStrategy(String str) {
        if (str == null || str.isEmpty()) {
            return SchedStrategy.PS;
        }
        // Try direct enum match first
        try {
            return SchedStrategy.valueOf(str);
        } catch (IllegalArgumentException e) {
            // Fall through to manual matching
        }
        // Try common string representations
        String upper = str.toUpperCase();
        if ("FCFS".equals(upper) || "FIFO".equals(upper)) return SchedStrategy.FCFS;
        if ("LCFS".equals(upper) || "LIFO".equals(upper)) return SchedStrategy.LCFS;
        if ("PS".equals(upper)) return SchedStrategy.PS;
        if ("INF".equals(upper)) return SchedStrategy.INF;
        if ("SIRO".equals(upper) || "RAND".equals(upper)) return SchedStrategy.SIRO;
        if ("HOL".equals(upper)) return SchedStrategy.HOL;
        if ("DPS".equals(upper)) return SchedStrategy.DPS;
        if ("GPS".equals(upper)) return SchedStrategy.GPS;
        if ("SEPT".equals(upper)) return SchedStrategy.SEPT;
        if ("LEPT".equals(upper)) return SchedStrategy.LEPT;
        if ("SJF".equals(upper)) return SchedStrategy.SJF;
        if ("LJF".equals(upper)) return SchedStrategy.LJF;
        if ("LCFSPR".equals(upper)) return SchedStrategy.LCFSPR;
        if ("FCFSPRIO".equals(upper)) return SchedStrategy.FCFSPRIO;
        if ("PSPRIO".equals(upper)) return SchedStrategy.PSPRIO;
        if ("DPSPRIO".equals(upper)) return SchedStrategy.DPSPRIO;
        if ("GPSPRIO".equals(upper)) return SchedStrategy.GPSPRIO;
        if ("REF".equals(upper)) return SchedStrategy.REF;
        if ("EXT".equals(upper)) return SchedStrategy.EXT;
        if ("FORK".equals(upper)) return SchedStrategy.FORK;
        if ("POLLING".equals(upper)) return SchedStrategy.POLLING;
        if ("LPS".equals(upper)) return SchedStrategy.LPS;
        // Default
        return SchedStrategy.PS;
    }

    /**
     * Parses a replacement strategy string to a {@link ReplacementStrategy} enum value.
     */
    private static ReplacementStrategy parseReplacementStrategy(String str) {
        if (str == null || str.isEmpty()) {
            return ReplacementStrategy.LRU;
        }
        try {
            return ReplacementStrategy.valueOf(str.toUpperCase());
        } catch (IllegalArgumentException e) {
            // Fall through
        }
        String upper = str.toUpperCase();
        if ("LRU".equals(upper)) return ReplacementStrategy.LRU;
        if ("FIFO".equals(upper)) return ReplacementStrategy.FIFO;
        if ("RANDOM".equals(upper) || "RR".equals(upper)) return ReplacementStrategy.RR;
        return ReplacementStrategy.LRU;
    }

    /**
     * Parses a routing strategy string to a {@link RoutingStrategy} enum value.
     */
    private static RoutingStrategy parseRoutingStrategy(String str) {
        if (str == null || str.isEmpty()) {
            return null;
        }
        try {
            return RoutingStrategy.valueOf(str);
        } catch (IllegalArgumentException e) {
            // Fall through to manual matching
        }
        String upper = str.toUpperCase();
        if ("RAND".equals(upper)) return RoutingStrategy.RAND;
        if ("PROB".equals(upper)) return RoutingStrategy.PROB;
        if ("RROBIN".equals(upper)) return RoutingStrategy.RROBIN;
        if ("WRROBIN".equals(upper)) return RoutingStrategy.WRROBIN;
        if ("JSQ".equals(upper)) return RoutingStrategy.JSQ;
        if ("SQ".equals(upper)) return RoutingStrategy.SQ;
        if ("DISABLED".equals(upper)) return RoutingStrategy.DISABLED;
        if ("FIRING".equals(upper)) return RoutingStrategy.FIRING;
        return null;
    }

    /**
     * Converts a DropStrategy enum to a schema-compatible string.
     */
    private static String dropStrategyToStr(DropStrategy ds) {
        if (ds == DropStrategy.Drop) return "drop";
        if (ds == DropStrategy.WaitingQueue) return "waitingQueue";
        if (ds == DropStrategy.BlockingAfterService) return "blockingAfterService";
        if (ds == DropStrategy.Retrial) return "retrial";
        if (ds == DropStrategy.RetrialWithLimit) return "retrialWithLimit";
        return "waitingQueue";
    }

    /**
     * Serialize the Krzesinski (1987) state-dependent routing declaration, or null
     * when the model carries none.
     *
     * <p>Centers travel by NODE NAME, which is what makes the block language
     * independent: a node reordering on either side cannot shift a center. Branch
     * index 1 denotes the complement M-V and is written as an empty list, keeping
     * the paper's own numbering, so the array read back is index-for-index the one
     * {@link Node#setStateDepRouting} expects.</p>
     *
     * <p>Only one declaration is carried, which is the same restriction the
     * struct itself has: {@link NetworkStruct#sdr} holds a single structure.</p>
     */
    private static JsonObject serializeStateDepRouting(Network model) {
        List<Node> nodes = model.getNodes();
        List<JobClass> classes = model.getClasses();
        for (int i = 0; i < nodes.size(); i++) {
            Node node = nodes.get(i);
            for (int r = 0; r < classes.size(); r++) {
                StateDepRouting sdr = node.getStateDepRouting(classes.get(r));
                if (sdr == null || sdr.branchNodes == null) {
                    continue;
                }
                JsonObject out = new JsonObject();
                out.addProperty("entry", node.getName());
                out.addProperty("departure", sdr.departureNode.getName());
                out.addProperty("class", classes.get(r).getName());
                JsonArray branches = new JsonArray();
                for (int b = 0; b < sdr.branchNodes.size(); b++) {
                    JsonArray bn = new JsonArray();
                    List<Node> centers = sdr.branchNodes.get(b);
                    if (centers != null) {
                        for (int q = 0; q < centers.size(); q++) {
                            bn.add(centers.get(q).getName());
                        }
                    }
                    branches.add(bn);
                }
                out.add("branches", branches);
                JsonArray level = new JsonArray();
                for (int b = 0; b < sdr.level.length; b++) {
                    level.add(sdr.level[b]);
                }
                JsonArray cArr = new JsonArray();
                for (int t = 0; t < sdr.C.length; t++) {
                    cArr.add(sdr.C[t]);
                }
                JsonArray dArr = new JsonArray();
                for (int t = 0; t < sdr.d.length; t++) {
                    JsonArray row = new JsonArray();
                    for (int b = 0; b < sdr.d[t].length; b++) {
                        row.add(sdr.d[t][b]);
                    }
                    dArr.add(row);
                }
                out.add("level", level);
                out.add("C", cArr);
                out.add("d", dArr);
                return out;
            }
        }
        return null;
    }

    /**
     * Restore a state-dependent routing declaration onto a loaded model.
     *
     * <p>Called AFTER link(P), whose uniform placeholder in the entry row the
     * declaration supersedes, and after the non-PROB strategies, which have named
     * that row SDR without saying what the subnetwork looks like.</p>
     */
    private static void restoreStateDepRouting(JsonObject sdrObj, Map<String, Node> nodeMap,
                                               Map<String, JobClass> classMap) {
        Node entry = nodeMap.get(sdrObj.get("entry").getAsString());
        Node departure = nodeMap.get(sdrObj.get("departure").getAsString());
        JobClass jobClass = classMap.get(sdrObj.get("class").getAsString());
        if (entry == null || departure == null || jobClass == null) {
            line_warning(mfilename(new Object() {
            }), "The stateDepRouting block names a node or class the model does not declare; "
                    + "the declaration is not restored.");
            return;
        }
        JsonArray branchArr = sdrObj.getAsJsonArray("branches");
        List<List<Node>> branches = new ArrayList<List<Node>>();
        for (int b = 0; b < branchArr.size(); b++) {
            JsonArray bn = branchArr.get(b).getAsJsonArray();
            List<Node> centers = new ArrayList<Node>();
            for (int q = 0; q < bn.size(); q++) {
                Node c = nodeMap.get(bn.get(q).getAsString());
                if (c == null) {
                    line_warning(mfilename(new Object() {
                    }), "The stateDepRouting block names node '%s', which the model does not "
                            + "declare; the declaration is not restored.", bn.get(q).getAsString());
                    return;
                }
                centers.add(c);
            }
            branches.add(centers);
        }
        JsonArray levelArr = sdrObj.getAsJsonArray("level");
        int[] level = new int[levelArr.size()];
        for (int b = 0; b < levelArr.size(); b++) {
            level[b] = levelArr.get(b).getAsInt();
        }
        JsonArray cArr = sdrObj.getAsJsonArray("C");
        double[] C = new double[cArr.size()];
        for (int t = 0; t < cArr.size(); t++) {
            C[t] = cArr.get(t).getAsDouble();
        }
        JsonArray dArr = sdrObj.getAsJsonArray("d");
        double[][] d = new double[dArr.size()][];
        for (int t = 0; t < dArr.size(); t++) {
            JsonArray row = dArr.get(t).getAsJsonArray();
            d[t] = new double[row.size()];
            for (int b = 0; b < row.size(); b++) {
                d[t][b] = row.get(b).getAsDouble();
            }
        }
        entry.setStateDepRouting(jobClass, departure, branches, level, C, d);
    }

    /**
     * Per-class cutoffs for materializing a class-dependence handle onto a
     * bounded lattice. A closed class cannot exceed its population; an open
     * class is unbounded, so it gets the same saturation cutoff the OI/PAS rate
     * table uses, beyond which beta is taken to be constant. Mirrors the rule in
     * the MATLAB writer (linemodel_save).
     */
    private static final int GD_MAX_LATTICE = 200000;

    /**
     * Materialize the network-level global (Whittle) dependence phi(n) onto the wire.
     * n is the full (nstations x nclasses) population matrix, but only the slots a
     * class can occupy carry a coordinate: a Source holds no jobs and a class with
     * zero per-class capacity at a station never appears there, so those entries are
     * pinned to 0. The restriction is lossless, since no DEP or PHASE event ever
     * fires at a slot the class cannot occupy. Mirrors gd_block in linemodel_save.m.
     */
    private static JsonObject serializeGlobalDependence(Network model) {
        NetworkStruct sn = model.getStruct();
        int M = sn.nstations;
        int K = sn.nclasses;
        int wcut = model.getGlobalDependenceCutoff();
        SerializableFunction<Matrix, Matrix> phi = model.getGlobalDependence();
        Matrix peak = model.getGlobalDependencePeak();

        List<String> stationNames = new ArrayList<String>();
        for (int i = 0; i < M; i++) {
            stationNames.add(model.getNodes().get((int) sn.stationToNode.get(i)).getName());
        }
        List<String> classNames = new ArrayList<String>();
        for (int r = 0; r < K; r++) {
            classNames.add(model.getClasses().get(r).getName());
        }

        List<Integer> slotSt = new ArrayList<Integer>();
        List<Integer> slotCl = new ArrayList<Integer>();
        List<Integer> cuts = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            if (sn.nodetype.get((int) sn.stationToNode.get(i)) == NodeType.Source) {
                continue;
            }
            for (int r = 0; r < K; r++) {
                double cap = sn.classcap.get(i, r);
                if (!(cap > 0)) {
                    continue;
                }
                double nj = sn.njobs.get(r);
                int c = Double.isFinite(nj) ? (int) Math.round(nj) : wcut;
                if (Double.isFinite(cap)) {
                    c = Math.min(c, (int) Math.round(cap));
                }
                slotSt.add(i);
                slotCl.add(r);
                cuts.add(Math.max(c, 0));
            }
        }

        int P = cuts.size();
        long total = 1;
        for (int d = 0; d < P; d++) {
            total *= (cuts.get(d) + 1);
            if (total > GD_MAX_LATTICE) {
                throw new IllegalArgumentException(
                        "The global dependence lattice exceeds the wire limit of " + GD_MAX_LATTICE
                        + " points (" + P + " varying station-class slots). Lower the wireCutoff "
                        + "argument of setGlobalDependence, or solve the model natively.");
            }
        }

        JsonObject blk = new JsonObject();
        blk.addProperty("type", "globalDependent");
        JsonArray stArr = new JsonArray();
        for (String nm : stationNames) stArr.add(nm);
        blk.add("stations", stArr);
        JsonArray clArr = new JsonArray();
        for (String nm : classNames) clArr.add(nm);
        blk.add("classes", clArr);
        JsonArray slotArr = new JsonArray();
        JsonArray cutArr = new JsonArray();
        for (int d = 0; d < P; d++) {
            JsonObject sm = new JsonObject();
            sm.addProperty("station", stationNames.get(slotSt.get(d)));
            sm.addProperty("class", classNames.get(slotCl.get(d)));
            slotArr.add(sm);
            cutArr.add(cuts.get(d));
        }
        blk.add("slots", slotArr);
        blk.add("cutoffs", cutArr);
        blk.addProperty("cutoff", wcut);

        JsonObject tbl = new JsonObject();
        int[] c = new int[P];
        for (long li = 0; li < total; li++) {
            long rem = li;
            for (int d = 0; d < P; d++) {
                c[d] = (int) (rem % (cuts.get(d) + 1));
                rem /= (cuts.get(d) + 1);
            }
            Matrix n = new Matrix(M, K);
            for (int d = 0; d < P; d++) {
                n.set(slotSt.get(d), slotCl.get(d), c[d]);
            }
            Matrix v = phi.apply(n);
            JsonArray va = new JsonArray();
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) {
                    double x;
                    if (v.getNumElements() == 1) {
                        x = v.get(0);
                    } else if (v.getNumCols() == 1) {
                        x = v.get(i, 0);
                    } else {
                        x = v.get(i, r);
                    }
                    va.add(Double.isFinite(x) ? x : 0.0);
                }
            }
            StringBuilder sb = new StringBuilder();
            if (P == 0) {
                sb.append("0");
            } else {
                for (int d = 0; d < P; d++) {
                    if (d > 0) sb.append(",");
                    sb.append(c[d]);
                }
            }
            tbl.add(sb.toString(), va);
        }
        blk.add("scaling", tbl);

        JsonArray pkArr = new JsonArray();
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                pkArr.add(peak == null || peak.isEmpty() ? 1.0 : peak.get(i, r));
            }
        }
        blk.add("peak", pkArr);
        return blk;
    }

    /**
     * Rebuild the network-level global (Whittle) dependence from the slot lattice
     * written by {@link #serializeGlobalDependence(Network)} and install it on the
     * model. Slots carry station and class NAMES, resolved through the model's own
     * index spaces, so a reordering on the writing side cannot shift a coordinate.
     */
    private static void restoreGlobalDependence(Network model, JsonObject blk,
                                                Map<String, Node> nodeMap,
                                                Map<String, JobClass> classMap) {
        if (blk == null || !blk.has("scaling")) {
            return;
        }
        String type = blk.has("type") ? blk.get("type").getAsString() : "";
        if (!"globalDependent".equals(type)) {
            return;
        }
        NetworkStruct sn = model.getStruct();
        final int M = sn.nstations;
        final int K = sn.nclasses;

        List<Integer> slotStList = new ArrayList<Integer>();
        List<Integer> slotClList = new ArrayList<Integer>();
        if (blk.has("slots")) {
            for (JsonElement el : blk.getAsJsonArray("slots")) {
                JsonObject sm = el.getAsJsonObject();
                Node node = nodeMap.get(sm.get("station").getAsString());
                JobClass jc = classMap.get(sm.get("class").getAsString());
                if (node == null || jc == null) {
                    return;
                }
                slotStList.add((int) sn.nodeToStation.get(node.getNodeIndex()));
                slotClList.add(jc.getIndex() - 1);
            }
        }
        final int P = slotStList.size();
        final int[] slotSt = new int[P];
        final int[] slotCl = new int[P];
        for (int d = 0; d < P; d++) {
            slotSt[d] = slotStList.get(d);
            slotCl[d] = slotClList.get(d);
        }
        final int[] cuts = new int[P];
        if (blk.has("cutoffs")) {
            JsonArray ca = blk.getAsJsonArray("cutoffs");
            for (int d = 0; d < P && d < ca.size(); d++) {
                cuts[d] = ca.get(d).getAsInt();
            }
        }
        int wcut = blk.has("cutoff") ? blk.get("cutoff").getAsInt() : 10;

        final HashMap<String, double[]> tbl = new HashMap<String, double[]>();
        for (Map.Entry<String, JsonElement> e : blk.getAsJsonObject("scaling").entrySet()) {
            JsonArray va = e.getValue().getAsJsonArray();
            double[] v = new double[va.size()];
            for (int i = 0; i < va.size(); i++) {
                v[i] = va.get(i).getAsDouble();
            }
            tbl.put(e.getKey(), v);
        }

        Matrix peak = new Matrix(M, K);
        if (blk.has("peak")) {
            JsonArray pa = blk.getAsJsonArray("peak");
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) {
                    int flat = i * K + r;
                    peak.set(i, r, flat < pa.size() ? pa.get(flat).getAsDouble() : 1.0);
                }
            }
        } else {
            peak = Matrix.ones(M, K);
        }

        SerializableFunction<Matrix, Matrix> phi = (Matrix n) -> {
            StringBuilder sb = new StringBuilder();
            if (P == 0) {
                sb.append("0");
            } else {
                for (int d = 0; d < P; d++) {
                    int x = (int) Math.round(n.get(slotSt[d], slotCl[d]));
                    if (x < 0) x = 0;
                    if (x > cuts[d]) x = cuts[d];
                    if (d > 0) sb.append(",");
                    sb.append(x);
                }
            }
            double[] v = tbl.get(sb.toString());
            Matrix out = new Matrix(M, K);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) {
                    out.set(i, r, v == null ? 1.0 : v[i * K + r]);
                }
            }
            return out;
        };
        model.setGlobalDependence(phi, peak, wcut);
    }

    private static int[] classDependenceCutoffs(Network model, int K) {
        int[] maxc = new int[K];
        for (int r = 0; r < K; r++) {
            JobClass jc = model.getClasses().get(r);
            if (jc instanceof ClosedClass
                    && Double.isFinite(((ClosedClass) jc).getPopulation())) {
                maxc[r] = (int) Math.round(((ClosedClass) jc).getPopulation());
            } else {
                maxc[r] = 10;
            }
        }
        return maxc;
    }

    /**
     * Rebuild a marking-dependent firing-rate dependence handle g(M) from the
     * materialized lattice written by the transition-mode serializer. M is the
     * node-indexed marking matrix; the enabling (place,class) slots are read off,
     * clamped to the cutoffs, and the tabulated scalar multiplier is looked up
     * (default 1 outside the tabulated range). Mirrors firingdep_table_to_handle
     * in linemodel_load.m and the Python reader.
     */
    private static SerializableFunction<Matrix, Double> firingDepTableToHandle(
            JsonObject frm, Map<String, Node> nodeMap, Map<String, JobClass> classMap) {
        if (frm == null || !frm.has("slots") || !frm.has("scaling")) return null;
        JsonArray slotsArr = frm.getAsJsonArray("slots");
        int P = slotsArr.size();
        final int[] slotRow = new int[P];
        final int[] slotCol = new int[P];
        for (int s = 0; s < P; s++) {
            JsonObject sm = slotsArr.get(s).getAsJsonObject();
            Node nd = nodeMap.get(sm.get("node").getAsString());
            JobClass jc = classMap.get(sm.get("class").getAsString());
            if (nd == null || jc == null) return null;
            slotRow[s] = nd.getNodeIndex();
            slotCol[s] = jc.getIndex() - 1;
        }
        final int[] cutoffs = new int[P];
        if (frm.has("cutoffs")) {
            JsonArray cutArr = frm.getAsJsonArray("cutoffs");
            for (int s = 0; s < P && s < cutArr.size(); s++) cutoffs[s] = cutArr.get(s).getAsInt();
        } else {
            for (int s = 0; s < P; s++) cutoffs[s] = Integer.MAX_VALUE;
        }
        final Map<String, Double> table = new HashMap<>();
        for (Map.Entry<String, JsonElement> e : frm.getAsJsonObject("scaling").entrySet()) {
            table.put(e.getKey(), e.getValue().getAsDouble());
        }
        return (Matrix m) -> {
            StringBuilder sb = new StringBuilder();
            for (int s = 0; s < P; s++) {
                int c = (int) Math.round(m.get(slotRow[s], slotCol[s]));
                if (c < 0) c = 0;
                if (cutoffs[s] < Integer.MAX_VALUE && c > cutoffs[s]) c = cutoffs[s];
                if (s > 0) sb.append(",");
                sb.append(c);
            }
            Double v = table.get(sb.toString());
            return v != null ? v : 1.0;   // marking outside the table is neutral
        };
    }

    /**
     * Resolves a Cache node field from either wire form. The canonical form is
     * the flat node-level key that every writer in the project emits; the
     * schema, however, also defines a nested {@code cache} object
     * (CacheConfig), which a hand-written or spec-conformant file may use and
     * which was previously ignored, silently loading such a cache with default
     * parameters. The nested object is consulted first and the flat key is the
     * fallback, matching the reference reader in linemodel_io.py.
     *
     * @param nodeObj   the node JSON object
     * @param nestedKey the key inside the nested {@code cache} object
     * @param flatKey   the equivalent node-level flat key
     * @return the resolved element, or null when neither form declares it
     */
    private static JsonElement cacheField(JsonObject nodeObj, String nestedKey, String flatKey) {
        if (nodeObj.has("cache") && nodeObj.get("cache").isJsonObject()) {
            JsonObject cacheObj = nodeObj.getAsJsonObject("cache");
            if (cacheObj.has(nestedKey) && !cacheObj.get(nestedKey).isJsonNull()) {
                return cacheObj.get(nestedKey);
            }
        }
        if (nodeObj.has(flatKey) && !nodeObj.get(flatKey).isJsonNull()) {
            return nodeObj.get(flatKey);
        }
        return null;
    }

    /**
     * Parses a drop strategy string to a DropStrategy enum value.
     */
    private static DropStrategy parseDropStrategy(String str) {
        if (str == null || str.isEmpty()) return DropStrategy.WaitingQueue;
        if ("drop".equals(str)) return DropStrategy.Drop;
        if ("waitingQueue".equals(str)) return DropStrategy.WaitingQueue;
        if ("blockingAfterService".equals(str)) return DropStrategy.BlockingAfterService;
        if ("retrial".equals(str)) return DropStrategy.Retrial;
        if ("retrialWithLimit".equals(str)) return DropStrategy.RetrialWithLimit;
        return DropStrategy.WaitingQueue;
    }
}
