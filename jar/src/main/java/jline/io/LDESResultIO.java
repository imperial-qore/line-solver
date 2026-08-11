/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import com.google.gson.*;

import jline.lang.NetworkStruct;
import jline.solvers.SolverResult;
import jline.solvers.ldes.LDESResult;
import jline.util.matrix.Matrix;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Provides save/load functionality for LDES solver results to/from JSON format.
 *
 * <p>Serializes {@link LDESResult} (and base {@link SolverResult}) fields to a JSON file
 * that can be parsed by the Python wrapper or other tools without needing JVM access.</p>
 *
 * <p>Usage:
 * <pre>
 *   // Save result
 *   LDESResultIO.save(result, sn, "result.json");
 *
 *   // Load result
 *   LDESResult result = LDESResultIO.load("result.json");
 * </pre>
 * </p>
 */
public class LDESResultIO {

    private static final String FORMAT_NAME = "ldes-result";
    private static final String FORMAT_VERSION = "1.0";

    // ========================================================================
    // SAVE
    // ========================================================================

    /**
     * Saves an {@link LDESResult} and associated {@link NetworkStruct} to a JSON file.
     *
     * @param result   the LDES solver result
     * @param sn       the network structure (for station/class names)
     * @param filename the output file path
     * @throws IOException if the file cannot be written
     */
    public static void save(LDESResult result, NetworkStruct sn, String filename) throws IOException {
        save(result, sn, filename, false);
    }

    /**
     * Saves an {@link LDESResult} and associated {@link NetworkStruct} to a JSON file.
     *
     * @param result      the LDES solver result
     * @param sn          the network structure (for station/class names)
     * @param filename    the output file path
     * @param trajectory  if true, include transient trajectory data (QNt, UNt, TNt, t)
     * @throws IOException if the file cannot be written
     */
    public static void save(LDESResult result, NetworkStruct sn, String filename, boolean trajectory) throws IOException {
        save(result, sn, filename, trajectory, false);
    }

    /**
     * Saves an {@link LDESResult} and associated {@link NetworkStruct} to a JSON file.
     *
     * @param result          the LDES solver result
     * @param sn              the network structure (for station/class names)
     * @param filename        the output file path
     * @param trajectory      if true, include transient trajectory data (QNt, UNt, TNt, t)
     * @param resptSamples    if true, include the per-job response time samples, the
     *                        empirical CDF input of the subprocess clients. They are
     *                        emitted as their own top-level block rather than inside
     *                        the trajectory data, because they are a distributional
     *                        output and not a time series, and they are gated because
     *                        a long run holds one sample per completion.
     * @throws IOException if the file cannot be written
     */
    public static void save(LDESResult result, NetworkStruct sn, String filename,
                            boolean trajectory, boolean resptSamples) throws IOException {
        JsonObject doc = new JsonObject();
        doc.addProperty("format", FORMAT_NAME);
        doc.addProperty("version", FORMAT_VERSION);

        // Solver metadata
        doc.addProperty("solver", result.solver != null ? result.solver : "LDES");
        doc.addProperty("method", result.method != null ? result.method : "default");
        doc.addProperty("runtime", result.runtime);
        doc.addProperty("converged", result.converged);
        doc.addProperty("stoppingReason", result.stoppingReason != null ? result.stoppingReason : "");
        doc.addProperty("convergenceBatches", result.convergenceBatches);
        doc.addProperty("totalSimulatedEvents", result.totalSimulatedEvents);

        // Dimensions and names
        JsonObject dimensions = new JsonObject();
        dimensions.addProperty("nstations", sn.nstations);
        dimensions.addProperty("nclasses", sn.nclasses);
        dimensions.addProperty("nchains", sn.nchains);

        // Station names (from nodenames, filtered to stations only)
        JsonArray stationNames = new JsonArray();
        for (int i = 0; i < sn.nstations; i++) {
            int nodeIdx = (int) sn.stationToNode.get(i);
            if (nodeIdx >= 0 && nodeIdx < sn.nodenames.size()) {
                stationNames.add(sn.nodenames.get(nodeIdx));
            } else {
                stationNames.add("Station" + i);
            }
        }
        dimensions.add("stationNames", stationNames);

        // Class names
        JsonArray classNames = new JsonArray();
        for (int r = 0; r < sn.nclasses; r++) {
            if (r < sn.classnames.size()) {
                classNames.add(sn.classnames.get(r));
            } else {
                classNames.add("Class" + r);
            }
        }
        dimensions.add("classNames", classNames);
        doc.add("dimensions", dimensions);

        // Performance metrics
        JsonObject metrics = new JsonObject();
        addMatrix(metrics, "QN", result.QN);
        addMatrix(metrics, "UN", result.UN);
        addMatrix(metrics, "RN", result.RN);
        addMatrix(metrics, "TN", result.TN);
        addMatrix(metrics, "AN", result.AN);
        addMatrix(metrics, "WN", result.WN);
        addMatrix(metrics, "CN", result.CN);
        addMatrix(metrics, "XN", result.XN);
        addMatrix(metrics, "DropRateJoin", result.DropRateJoin);
        doc.add("metrics", metrics);

        // Exact joint-state residence-time histogram (when --export-histogram is set),
        // used by MATLAB/native-Python to evaluate arbitrary Markov rewards host-side.
        if (result.stateHistogramSpace != null && result.stateHistogramTime != null
                && result.stateHistogramSpace.getNumRows() > 0) {
            JsonObject hist = new JsonObject();
            addMatrix(hist, "space", result.stateHistogramSpace);
            addMatrix(hist, "time", result.stateHistogramTime);
            if (result.stateTrajectorySpace != null && result.stateTrajectoryTime != null
                    && result.stateTrajectorySpace.getNumRows() > 0) {
                addMatrix(hist, "trajSpace", result.stateTrajectorySpace);
                addMatrix(hist, "trajTime", result.stateTrajectoryTime);
            }
            doc.add("stateHistogram", hist);
        }

        // Per-cache hit/miss ratios and expected latency (set on the Cache nodes
        // by the simulator). Indexed by job class. Lets the native bridge restore
        // cache.get_hit_ratio()/get_miss_ratio() and node-level hit/miss tables.
        JsonObject cacheMetrics = new JsonObject();
        if (sn.nodes != null) {
            for (jline.lang.nodes.Node node : sn.nodes) {
                if (node instanceof jline.lang.nodes.Cache) {
                    jline.lang.nodes.Cache cache = (jline.lang.nodes.Cache) node;
                    JsonObject cm = new JsonObject();
                    addMatrix(cm, "hit", cache.getHitRatio());
                    addMatrix(cm, "delayed", cache.getDelayedHitRatio());
                    addMatrix(cm, "miss", cache.getMissRatio());
                    addMatrix(cm, "latency", cache.getResidT());
                    addMatrix(cm, "hitList", cache.getHitRatioByList());
                    addMatrix(cm, "itemProb", cache.getItemProb());
                    cacheMetrics.add(node.getName(), cm);
                }
            }
        }
        if (cacheMetrics.size() > 0) {
            doc.add("cacheMetrics", cacheMetrics);
        }

        // Finite capacity region (FCR) metrics: per region (rows) x class (cols).
        // WeightNfcr is the weighted occupation (Total Weight) and MemOccNfcr the
        // memory occupation; row sums give the region aggregates.
        if (sn.nregions > 0 && result.QNfcr != null) {
            JsonObject fcr = new JsonObject();
            fcr.addProperty("nregions", sn.nregions);
            addMatrix(fcr, "QNfcr", result.QNfcr);
            addMatrix(fcr, "UNfcr", result.UNfcr);
            addMatrix(fcr, "RNfcr", result.RNfcr);
            addMatrix(fcr, "TNfcr", result.TNfcr);
            addMatrix(fcr, "ANfcr", result.ANfcr);
            addMatrix(fcr, "WNfcr", result.WNfcr);
            addMatrix(fcr, "WeightNfcr", result.WeightNfcr);
            addMatrix(fcr, "MemOccNfcr", result.MemOccNfcr);
            addMatrix(fcr, "DropRateNfcr", result.DropRateNfcr);
            doc.add("fcr", fcr);
        }

        // Sample counts per metric
        JsonObject sampleCounts = new JsonObject();
        addMatrix(sampleCounts, "QNSamples", result.QNSamples);
        addMatrix(sampleCounts, "UNSamples", result.UNSamples);
        addMatrix(sampleCounts, "RNSamples", result.RNSamples);
        addMatrix(sampleCounts, "TNSamples", result.TNSamples);
        doc.add("sampleCounts", sampleCounts);

        // Confidence intervals
        JsonObject ci = new JsonObject();
        addMatrix(ci, "QNCI", result.QNCI);
        addMatrix(ci, "UNCI", result.UNCI);
        addMatrix(ci, "RNCI", result.RNCI);
        addMatrix(ci, "TNCI", result.TNCI);
        addMatrix(ci, "ANCI", result.ANCI);
        addMatrix(ci, "WNCI", result.WNCI);
        doc.add("confidenceIntervals", ci);

        // Relative precision
        JsonObject relPrec = new JsonObject();
        addMatrix(relPrec, "QNRelPrec", result.QNRelPrec);
        addMatrix(relPrec, "UNRelPrec", result.UNRelPrec);
        addMatrix(relPrec, "RNRelPrec", result.RNRelPrec);
        addMatrix(relPrec, "TNRelPrec", result.TNRelPrec);
        doc.add("relativePrecision", relPrec);

        // Impatience statistics (only if any non-null)
        if (result.renegedCustomers != null || result.balkedCustomers != null
                || result.retriedCustomers != null || result.avgOrbitSize != null) {
            JsonObject impatience = new JsonObject();
            addMatrix(impatience, "renegedCustomers", result.renegedCustomers);
            addMatrix(impatience, "avgRenegingWaitTime", result.avgRenegingWaitTime);
            addMatrix(impatience, "renegingRate", result.renegingRate);
            addMatrix(impatience, "balkedCustomers", result.balkedCustomers);
            addMatrix(impatience, "balkingProbability", result.balkingProbability);
            addMatrix(impatience, "retriedCustomers", result.retriedCustomers);
            addMatrix(impatience, "retrialDropped", result.retrialDropped);
            addMatrix(impatience, "avgOrbitSize", result.avgOrbitSize);
            doc.add("impatience", impatience);
        }

        // Per-job response time samples, the empirical CDF input
        if (resptSamples && result.respTimeSamples != null) {
            JsonArray rtsArr = new JsonArray();
            for (int i = 0; i < result.respTimeSamples.length; i++) {
                JsonArray stationArr = new JsonArray();
                for (int k = 0; k < result.respTimeSamples[i].length; k++) {
                    JsonArray samples = new JsonArray();
                    if (result.respTimeSamples[i][k] != null) {
                        for (Double v : result.respTimeSamples[i][k]) {
                            if (v != null && !Double.isNaN(v)) {
                                samples.add(v);
                            }
                        }
                    }
                    stationArr.add(samples);
                }
                rtsArr.add(stationArr);
            }
            doc.add("respTimeSamples", rtsArr);
        }

        // Trajectory data (transient time series)
        if (trajectory) {
            JsonObject transientData = new JsonObject();

            // Time vector
            if (result.t != null) {
                JsonArray tArr = new JsonArray();
                for (int i = 0; i < result.t.getNumRows(); i++) {
                    tArr.add(result.t.get(i, 0));
                }
                transientData.add("t", tArr);
            }

            // QNt[station][class] -> Matrix(numTimePoints, 2) with cols [value, time]
            if (result.QNt != null) {
                JsonArray qntArr = new JsonArray();
                for (int i = 0; i < result.QNt.length; i++) {
                    JsonArray stationArr = new JsonArray();
                    for (int k = 0; k < result.QNt[i].length; k++) {
                        addMatrix(stationArr, result.QNt[i][k]);
                    }
                    qntArr.add(stationArr);
                }
                transientData.add("QNt", qntArr);
            }

            // UNt[station][class]
            if (result.UNt != null) {
                JsonArray untArr = new JsonArray();
                for (int i = 0; i < result.UNt.length; i++) {
                    JsonArray stationArr = new JsonArray();
                    for (int k = 0; k < result.UNt[i].length; k++) {
                        addMatrix(stationArr, result.UNt[i][k]);
                    }
                    untArr.add(stationArr);
                }
                transientData.add("UNt", untArr);
            }

            // TNt[station][class]
            if (result.TNt != null) {
                JsonArray tntArr = new JsonArray();
                for (int i = 0; i < result.TNt.length; i++) {
                    JsonArray stationArr = new JsonArray();
                    for (int k = 0; k < result.TNt[i].length; k++) {
                        addMatrix(stationArr, result.TNt[i][k]);
                    }
                    tntArr.add(stationArr);
                }
                transientData.add("TNt", tntArr);
            }

            // Response time samples [station][class] -> list of doubles
            if (result.respTimeSamples != null) {
                JsonArray rtsArr = new JsonArray();
                for (int i = 0; i < result.respTimeSamples.length; i++) {
                    JsonArray stationArr = new JsonArray();
                    for (int k = 0; k < result.respTimeSamples[i].length; k++) {
                        JsonArray samples = new JsonArray();
                        if (result.respTimeSamples[i][k] != null) {
                            for (Double v : result.respTimeSamples[i][k]) {
                                if (v != null && !Double.isNaN(v)) {
                                    samples.add(v);
                                }
                            }
                        }
                        stationArr.add(samples);
                    }
                    rtsArr.add(stationArr);
                }
                transientData.add("respTimeSamples", rtsArr);
            }

            doc.add("transient", transientData);
        }

        Gson gson = new GsonBuilder().setPrettyPrinting().create();
        Writer writer = new BufferedWriter(new FileWriter(filename));
        try {
            gson.toJson(doc, writer);
        } finally {
            writer.close();
        }
    }

    // ========================================================================
    // LOAD
    // ========================================================================

    /**
     * Loads an {@link LDESResult} from a JSON file.
     *
     * @param filename the input file path
     * @return the deserialized LDESResult
     * @throws IOException if the file cannot be read or the format is invalid
     */
    public static LDESResult load(String filename) throws IOException {
        Reader reader = new BufferedReader(new FileReader(filename));
        JsonObject doc;
        try {
            doc = JsonParser.parseReader(reader).getAsJsonObject();
        } finally {
            reader.close();
        }

        String format = doc.has("format") ? doc.get("format").getAsString() : "";
        if (!FORMAT_NAME.equals(format)) {
            throw new IOException("Invalid format: expected '" + FORMAT_NAME + "', got '" + format + "'");
        }

        LDESResult result = new LDESResult();

        // Solver metadata
        result.solver = getStringOrNull(doc, "solver");
        result.method = getStringOrNull(doc, "method");
        result.runtime = doc.has("runtime") ? doc.get("runtime").getAsDouble() : 0.0;
        result.converged = doc.has("converged") && doc.get("converged").getAsBoolean();
        result.stoppingReason = getStringOrNull(doc, "stoppingReason");
        result.convergenceBatches = doc.has("convergenceBatches") ? doc.get("convergenceBatches").getAsInt() : 0;
        result.totalSimulatedEvents = doc.has("totalSimulatedEvents") ? doc.get("totalSimulatedEvents").getAsLong() : 0;

        // Performance metrics
        if (doc.has("metrics")) {
            JsonObject metrics = doc.getAsJsonObject("metrics");
            result.QN = loadMatrix(metrics, "QN");
            result.UN = loadMatrix(metrics, "UN");
            result.RN = loadMatrix(metrics, "RN");
            result.TN = loadMatrix(metrics, "TN");
            result.AN = loadMatrix(metrics, "AN");
            result.WN = loadMatrix(metrics, "WN");
            result.CN = loadMatrix(metrics, "CN");
            result.XN = loadMatrix(metrics, "XN");
            result.DropRateJoin = loadMatrix(metrics, "DropRateJoin");
        }

        // Finite capacity region (FCR) metrics
        if (doc.has("fcr")) {
            JsonObject fcr = doc.getAsJsonObject("fcr");
            result.QNfcr = loadMatrix(fcr, "QNfcr");
            result.UNfcr = loadMatrix(fcr, "UNfcr");
            result.RNfcr = loadMatrix(fcr, "RNfcr");
            result.TNfcr = loadMatrix(fcr, "TNfcr");
            result.ANfcr = loadMatrix(fcr, "ANfcr");
            result.WNfcr = loadMatrix(fcr, "WNfcr");
            result.WeightNfcr = loadMatrix(fcr, "WeightNfcr");
            result.MemOccNfcr = loadMatrix(fcr, "MemOccNfcr");
            result.DropRateNfcr = loadMatrix(fcr, "DropRateNfcr");
        }

        // Sample counts
        if (doc.has("sampleCounts")) {
            JsonObject sc = doc.getAsJsonObject("sampleCounts");
            result.QNSamples = loadMatrix(sc, "QNSamples");
            result.UNSamples = loadMatrix(sc, "UNSamples");
            result.RNSamples = loadMatrix(sc, "RNSamples");
            result.TNSamples = loadMatrix(sc, "TNSamples");
        }

        // Confidence intervals
        if (doc.has("confidenceIntervals")) {
            JsonObject ci = doc.getAsJsonObject("confidenceIntervals");
            result.QNCI = loadMatrix(ci, "QNCI");
            result.UNCI = loadMatrix(ci, "UNCI");
            result.RNCI = loadMatrix(ci, "RNCI");
            result.TNCI = loadMatrix(ci, "TNCI");
            result.ANCI = loadMatrix(ci, "ANCI");
            result.WNCI = loadMatrix(ci, "WNCI");
        }

        // Relative precision
        if (doc.has("relativePrecision")) {
            JsonObject relPrec = doc.getAsJsonObject("relativePrecision");
            result.QNRelPrec = loadMatrix(relPrec, "QNRelPrec");
            result.UNRelPrec = loadMatrix(relPrec, "UNRelPrec");
            result.RNRelPrec = loadMatrix(relPrec, "RNRelPrec");
            result.TNRelPrec = loadMatrix(relPrec, "TNRelPrec");
        }

        // Impatience
        if (doc.has("impatience")) {
            JsonObject impatience = doc.getAsJsonObject("impatience");
            result.renegedCustomers = loadMatrix(impatience, "renegedCustomers");
            result.avgRenegingWaitTime = loadMatrix(impatience, "avgRenegingWaitTime");
            result.renegingRate = loadMatrix(impatience, "renegingRate");
            result.balkedCustomers = loadMatrix(impatience, "balkedCustomers");
            result.balkingProbability = loadMatrix(impatience, "balkingProbability");
            result.retriedCustomers = loadMatrix(impatience, "retriedCustomers");
            result.retrialDropped = loadMatrix(impatience, "retrialDropped");
            result.avgOrbitSize = loadMatrix(impatience, "avgOrbitSize");
        }

        return result;
    }

    /**
     * Returns station names from a result JSON file, or null if unavailable.
     *
     * @param filename the result JSON file path
     * @return list of station names
     * @throws IOException if the file cannot be read
     */
    public static List<String> loadStationNames(String filename) throws IOException {
        Reader reader = new BufferedReader(new FileReader(filename));
        JsonObject doc;
        try {
            doc = JsonParser.parseReader(reader).getAsJsonObject();
        } finally {
            reader.close();
        }
        if (!doc.has("dimensions")) {
            return null;
        }
        JsonObject dims = doc.getAsJsonObject("dimensions");
        if (!dims.has("stationNames")) {
            return null;
        }
        JsonArray arr = dims.getAsJsonArray("stationNames");
        List<String> names = new ArrayList<String>();
        for (int i = 0; i < arr.size(); i++) {
            names.add(arr.get(i).getAsString());
        }
        return names;
    }

    /**
     * Returns class names from a result JSON file, or null if unavailable.
     *
     * @param filename the result JSON file path
     * @return list of class names
     * @throws IOException if the file cannot be read
     */
    public static List<String> loadClassNames(String filename) throws IOException {
        Reader reader = new BufferedReader(new FileReader(filename));
        JsonObject doc;
        try {
            doc = JsonParser.parseReader(reader).getAsJsonObject();
        } finally {
            reader.close();
        }
        if (!doc.has("dimensions")) {
            return null;
        }
        JsonObject dims = doc.getAsJsonObject("dimensions");
        if (!dims.has("classNames")) {
            return null;
        }
        JsonArray arr = dims.getAsJsonArray("classNames");
        List<String> names = new ArrayList<String>();
        for (int i = 0; i < arr.size(); i++) {
            names.add(arr.get(i).getAsString());
        }
        return names;
    }

    // ========================================================================
    // Helpers
    // ========================================================================

    /**
     * Adds a Matrix to a JsonObject as a JSON array of arrays.
     * NaN values become JSON null, Infinity becomes 1e308.
     */
    private static void addMatrix(JsonObject parent, String key, Matrix matrix) {
        if (matrix == null) {
            parent.add(key, JsonNull.INSTANCE);
            return;
        }
        int rows = matrix.getNumRows();
        int cols = matrix.getNumCols();
        JsonArray outer = new JsonArray();
        for (int i = 0; i < rows; i++) {
            JsonArray row = new JsonArray();
            for (int j = 0; j < cols; j++) {
                double v = matrix.get(i, j);
                if (Double.isNaN(v)) {
                    row.add(JsonNull.INSTANCE);
                } else if (Double.isInfinite(v)) {
                    row.add(v > 0 ? 1e308 : -1e308);
                } else {
                    row.add(v);
                }
            }
            outer.add(row);
        }
        parent.add(key, outer);
    }

    /**
     * Adds a Matrix to a JsonArray as a JSON array of arrays.
     * NaN values become JSON null, Infinity becomes 1e308.
     */
    private static void addMatrix(JsonArray parent, Matrix matrix) {
        if (matrix == null) {
            parent.add(JsonNull.INSTANCE);
            return;
        }
        int rows = matrix.getNumRows();
        int cols = matrix.getNumCols();
        JsonArray outer = new JsonArray();
        for (int i = 0; i < rows; i++) {
            JsonArray row = new JsonArray();
            for (int j = 0; j < cols; j++) {
                double v = matrix.get(i, j);
                if (Double.isNaN(v)) {
                    row.add(JsonNull.INSTANCE);
                } else if (Double.isInfinite(v)) {
                    row.add(v > 0 ? 1e308 : -1e308);
                } else {
                    row.add(v);
                }
            }
            outer.add(row);
        }
        parent.add(outer);
    }

    /**
     * Loads a Matrix from a JsonObject field containing a JSON array of arrays.
     * JSON null values become NaN.
     */
    private static Matrix loadMatrix(JsonObject parent, String key) {
        if (!parent.has(key) || parent.get(key).isJsonNull()) {
            return null;
        }
        JsonArray outer = parent.getAsJsonArray(key);
        int rows = outer.size();
        if (rows == 0) {
            return new Matrix(0, 0);
        }
        int cols = outer.get(0).getAsJsonArray().size();
        Matrix m = new Matrix(rows, cols);
        for (int i = 0; i < rows; i++) {
            JsonArray row = outer.get(i).getAsJsonArray();
            for (int j = 0; j < cols; j++) {
                JsonElement elem = row.get(j);
                if (elem.isJsonNull()) {
                    m.set(i, j, Double.NaN);
                } else {
                    m.set(i, j, elem.getAsDouble());
                }
            }
        }
        return m;
    }

    private static String getStringOrNull(JsonObject obj, String key) {
        if (!obj.has(key) || obj.get(key).isJsonNull()) {
            return null;
        }
        return obj.get(key).getAsString();
    }
}
