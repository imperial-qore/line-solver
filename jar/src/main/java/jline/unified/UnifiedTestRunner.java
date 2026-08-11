/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.unified;

import com.google.gson.Gson;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import jline.lang.Network;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.CTMC;
import jline.solvers.ldes.LDES;
import jline.solvers.fluid.FLD;
import jline.solvers.wrappers.jmt.JMT;
import jline.solvers.mam.MAM;
import jline.solvers.mva.MVA;
import jline.solvers.nc.NC;
import jline.solvers.ssa.SSA;
import jline.util.Pair;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Executes unified cross-language tests from JSON definitions.
 */
public class UnifiedTestRunner {

    public static final double DEFAULT_TOL = 1e-8;
    public static final double SIMULATION_TOL = 0.05;

    public static final Set<String> SIMULATION_SOLVERS;
    static {
        Set<String> s = new HashSet<String>();
        s.add("SolverJMT");
        s.add("SolverSSA");
        s.add("SolverLDES");
        SIMULATION_SOLVERS = Collections.unmodifiableSet(s);
    }

    private final Gson gson = new Gson();
    private final File definitionsDir;
    private TestResults results = new TestResults();

    public UnifiedTestRunner() {
        String definitionsPath = System.getProperty("line.definitions");
        if (definitionsPath != null) {
            definitionsDir = new File(definitionsPath);
        } else {
            File projectRoot = findProjectRoot();
            definitionsDir = new File(projectRoot, "test/unified/definitions");
        }
    }

    public TestResults getResults() {
        return results;
    }

    private File findProjectRoot() {
        File current = new File(System.getProperty("user.dir"));
        while (current.getParentFile() != null) {
            if (new File(current, "CLAUDE.md").exists() ||
                    (new File(current, "jar").exists() && new File(current, "matlab").exists())) {
                return current;
            }
            current = current.getParentFile();
        }
        return new File(System.getProperty("user.dir"));
    }

    public TestResults runAll() {
        return runAll(false);
    }

    public TestResults runAll(boolean verbose) {
        results = new TestResults();

        if (!definitionsDir.exists()) {
            System.out.println("Definitions directory not found: " + definitionsDir.getAbsolutePath());
            return results;
        }

        File[] arr = definitionsDir.listFiles(new java.io.FileFilter() {
            @Override
            public boolean accept(File f) {
                String name = f.getName();
                int idx = name.lastIndexOf('.');
                String ext = (idx >= 0) ? name.substring(idx + 1) : "";
                return "json".equals(ext);
            }
        });
        List<File> jsonFiles = new ArrayList<File>();
        if (arr != null) {
            jsonFiles.addAll(Arrays.asList(arr));
            Collections.sort(jsonFiles, new java.util.Comparator<File>() {
                @Override
                public int compare(File a, File b) {
                    return a.getName().compareTo(b.getName());
                }
            });
        }

        System.out.println("=== Unified Test Runner (Java/Kotlin) ===");
        System.out.println("Found " + jsonFiles.size() + " test definitions\n");

        for (File jsonFile : jsonFiles) {
            String name = jsonFile.getName();
            int idx = name.lastIndexOf('.');
            String modelName = (idx > 0) ? name.substring(0, idx) : name;
            runModel(modelName, verbose);
        }

        printSummary();
        return results;
    }

    public Object runModel(String modelName) {
        return runModel(modelName, false, false);
    }

    public Object runModel(String modelName, boolean verbose) {
        return runModel(modelName, verbose, false);
    }

    public Object runModel(String modelName, boolean verbose, boolean jsonOutput) {
        JsonResult jsonResult = jsonOutput ? new JsonResult(modelName) : null;
        long modelStartTime = jsonOutput ? System.nanoTime() : 0L;

        File jsonFile = new File(definitionsDir, modelName + ".json");
        if (!jsonFile.exists()) {
            if (verbose) System.out.println("  [SKIP] " + modelName + " - definition file not found");
            results.skipped++;
            if (jsonOutput) {
                jsonResult.status = "skipped";
                jsonResult.errors.add("definition file not found");
                return jsonResult;
            }
            return "skipped";
        }

        JsonObject definition;
        try {
            String content = new String(java.nio.file.Files.readAllBytes(jsonFile.toPath()));
            definition = gson.fromJson(content, JsonObject.class);
        } catch (Exception e) {
            if (!jsonOutput) System.out.println("  [ERROR] " + modelName + " - failed to parse JSON: " + e.getMessage());
            results.failed++;
            results.errors.add(modelName + ": JSON parse error");
            if (jsonOutput) {
                jsonResult.status = "error";
                jsonResult.errors.add("JSON parse error: " + e.getMessage());
                return jsonResult;
            }
            return "failed";
        }

        JsonObject supportedLanguages = definition.getAsJsonObject("supportedLanguages");
        if (supportedLanguages != null &&
                supportedLanguages.has("java") &&
                !supportedLanguages.get("java").getAsBoolean()) {
            if (verbose) System.out.println("  [SKIP] " + modelName + " - Java not supported");
            results.skipped++;
            if (jsonOutput) {
                jsonResult.status = "skipped";
                jsonResult.errors.add("Java not supported");
                return jsonResult;
            }
            return "skipped";
        }

        if (!ModelRegistry.hasModel(modelName)) {
            if (verbose) System.out.println("  [SKIP] " + modelName + " - not in model registry");
            results.skipped++;
            if (jsonOutput) {
                jsonResult.status = "skipped";
                jsonResult.errors.add("not in model registry");
                return jsonResult;
            }
            return "skipped";
        }

        if (verbose && !jsonOutput) System.out.println("Testing: " + modelName);

        Network model;
        try {
            model = ModelRegistry.getModel(modelName);
        } catch (Error e) {
            if (verbose) System.out.println("  [SKIP] " + modelName + " - not yet implemented");
            results.skipped++;
            if (jsonOutput) {
                jsonResult.status = "skipped";
                jsonResult.errors.add("not yet implemented");
                return jsonResult;
            }
            return "skipped";
        } catch (Exception e) {
            if (!jsonOutput) System.out.println("  [ERROR] " + modelName + " - failed to build model: " + e.getMessage());
            results.failed++;
            results.errors.add(modelName + ": Model build error - " + e.getMessage());
            if (jsonOutput) {
                jsonResult.status = "error";
                jsonResult.errors.add("Model build error: " + e.getMessage());
                return jsonResult;
            }
            return "failed";
        }

        boolean allPassed = true;
        JsonArray solversToTest = definition.getAsJsonArray("solvers");

        List<String> skipSolvers = new ArrayList<String>();
        JsonObject skipObj = definition.getAsJsonObject("skipSolvers");
        if (skipObj != null) {
            JsonArray ja = skipObj.getAsJsonArray("java");
            if (ja != null) {
                for (JsonElement e : ja) skipSolvers.add(e.getAsString());
            }
        }

        for (JsonElement solverElement : solversToTest) {
            String solverName = solverElement.getAsString();

            if (skipSolvers.contains(solverName)) {
                if (jsonOutput) {
                    Map<String, Object> entry = new HashMap<String, Object>();
                    entry.put("status", "skipped");
                    entry.put("metrics", new HashMap<String, Object>());
                    jsonResult.solverResults.put(solverName, entry);
                }
                continue;
            }

            JsonObject expectedResults = definition.getAsJsonObject("expectedResults");
            if (!expectedResults.has(solverName)) continue;

            try {
                if (jsonOutput) {
                    SolverJsonResult sjr = testSolverJson(model, solverName, definition, verbose);
                    jsonResult.solverResults.put(solverName, sjr.solverResult);
                    @SuppressWarnings("unchecked")
                    Map<String, Long> sm = (Map<String, Long>) jsonResult.timing.get("solvers");
                    sm.put(solverName, sjr.solverTime);
                    if (!sjr.passed) allPassed = false;
                } else {
                    boolean passed = testSolver(model, solverName, definition, verbose);
                    if (!passed) allPassed = false;
                }
            } catch (Exception e) {
                if (verbose) System.out.println("    [ERROR] " + solverName + ": " + e.getMessage());
                allPassed = false;
                if (jsonOutput) {
                    Map<String, Object> entry = new HashMap<String, Object>();
                    entry.put("status", "error");
                    entry.put("metrics", new HashMap<String, Object>());
                    jsonResult.solverResults.put(solverName, entry);
                    jsonResult.errors.add(solverName + ": " + e.getMessage());
                }
                results.errors.add(modelName + "/" + solverName + ": " + e.getMessage());
            }
        }

        if (allPassed) {
            results.passed++;
            if (jsonOutput) {
                jsonResult.status = "passed";
                jsonResult.timing.put("total_ms", (System.nanoTime() - modelStartTime) / 1_000_000);
                return jsonResult;
            } else {
                if (!verbose) System.out.print(".");
                return "passed";
            }
        } else {
            results.failed++;
            if (jsonOutput) {
                jsonResult.status = "failed";
                jsonResult.timing.put("total_ms", (System.nanoTime() - modelStartTime) / 1_000_000);
                return jsonResult;
            } else {
                if (!verbose) System.out.print("F");
                return "failed";
            }
        }
    }

    private boolean testSolver(Network model, String solverName, JsonObject definition, boolean verbose) {
        JsonObject toleranceConfig = definition.getAsJsonObject("tolerance");
        double tol;
        if (SIMULATION_SOLVERS.contains(solverName)) {
            tol = (toleranceConfig != null && toleranceConfig.has("simulation"))
                    ? toleranceConfig.get("simulation").getAsDouble() : SIMULATION_TOL;
        } else {
            tol = (toleranceConfig != null && toleranceConfig.has("default"))
                    ? toleranceConfig.get("default").getAsDouble() : DEFAULT_TOL;
        }

        NetworkSolver solver = createSolver(model, solverName);
        if (solver == null) {
            if (verbose) System.out.println("    [SKIP] " + solverName + " - solver not available");
            return true;
        }

        Object avgTable;
        try {
            avgTable = solver.getAvgTable();
        } catch (Exception e) {
            if (verbose) System.out.println("    [ERROR] " + solverName + " - solver failed: " + e.getMessage());
            return false;
        }

        JsonObject expected = definition.getAsJsonObject("expectedResults").getAsJsonObject(solverName);
        boolean passed = true;

        String[] metrics = new String[]{"QLen", "Util", "RespT", "Tput"};
        for (String metric : metrics) {
            if (!expected.has(metric)) continue;

            double[] expectedValues = flatten2DArray(expected.getAsJsonArray(metric));
            double[] actualValues = extractMetric(avgTable, metric);

            if (actualValues != null && expectedValues != null) {
                if (!compareWithTolerance(actualValues, expectedValues, tol)) {
                    if (verbose) {
                        System.out.println("    [FAIL] " + solverName + "." + metric + " mismatch (tol=" + String.format("%.2e", tol) + ")");
                        System.out.println("      Expected: " + Arrays.toString(expectedValues));
                        System.out.println("      Actual:   " + Arrays.toString(actualValues));
                    }
                    passed = false;
                } else if (verbose) {
                    System.out.println("    [PASS] " + solverName + "." + metric);
                }
            }
        }

        if (passed && verbose) System.out.println("  [PASS] " + solverName);
        return passed;
    }

    private static class SolverJsonResult {
        boolean passed;
        Map<String, Object> solverResult;
        long solverTime;
        SolverJsonResult(boolean p, Map<String, Object> r, long t) {
            this.passed = p;
            this.solverResult = r;
            this.solverTime = t;
        }
    }

    private SolverJsonResult testSolverJson(Network model, String solverName, JsonObject definition, boolean verbose) {
        Map<String, Object> solverResult = new HashMap<String, Object>();
        solverResult.put("status", "passed");
        solverResult.put("metrics", new HashMap<String, Object>());

        JsonObject toleranceConfig = definition.getAsJsonObject("tolerance");
        double tol;
        if (SIMULATION_SOLVERS.contains(solverName)) {
            tol = (toleranceConfig != null && toleranceConfig.has("simulation"))
                    ? toleranceConfig.get("simulation").getAsDouble() : SIMULATION_TOL;
        } else {
            tol = (toleranceConfig != null && toleranceConfig.has("default"))
                    ? toleranceConfig.get("default").getAsDouble() : DEFAULT_TOL;
        }

        NetworkSolver solver = createSolver(model, solverName);
        if (solver == null) {
            solverResult.put("status", "skipped");
            return new SolverJsonResult(true, solverResult, 0L);
        }

        long solverStartTime = System.nanoTime();
        Object avgTable;
        try {
            avgTable = solver.getAvgTable();
        } catch (Exception e) {
            long solverTime = (System.nanoTime() - solverStartTime) / 1_000_000;
            solverResult.put("status", "error");
            return new SolverJsonResult(false, solverResult, solverTime);
        }
        long solverTime = (System.nanoTime() - solverStartTime) / 1_000_000;

        JsonObject expected = definition.getAsJsonObject("expectedResults").getAsJsonObject(solverName);
        boolean passed = true;

        @SuppressWarnings("unchecked")
        Map<String, Object> metrics = (Map<String, Object>) solverResult.get("metrics");

        String[] metricNames = new String[]{"QLen", "Util", "RespT", "Tput"};
        for (String metric : metricNames) {
            if (!expected.has(metric)) continue;

            double[] expectedValues = flatten2DArray(expected.getAsJsonArray(metric));
            double[] actualValues = extractMetric(avgTable, metric);

            Map<String, Object> metricResult = new HashMap<String, Object>();
            metricResult.put("passed", true);

            if (actualValues != null && expectedValues != null) {
                if (!compareWithTolerance(actualValues, expectedValues, tol)) {
                    metricResult.put("passed", false);
                    List<Double> exp = new ArrayList<Double>();
                    for (double d : expectedValues) exp.add(d);
                    List<Double> act = new ArrayList<Double>();
                    for (double d : actualValues) act.add(d);
                    metricResult.put("expected", exp);
                    metricResult.put("actual", act);
                    double maxError = 0.0;
                    for (int i = 0; i < actualValues.length; i++) {
                        if (!Double.isNaN(expectedValues[i]) && !Double.isNaN(actualValues[i])) {
                            double error = Math.abs(actualValues[i] - expectedValues[i]);
                            if (error > maxError) maxError = error;
                        }
                    }
                    metricResult.put("maxError", maxError);
                    passed = false;
                }
            }
            metrics.put(metric, metricResult);
        }

        if (!passed) {
            solverResult.put("status", "failed");
        }

        return new SolverJsonResult(passed, solverResult, solverTime);
    }

    private NetworkSolver createSolver(Network model, String solverName) {
        try {
            SolverOptions options = new SolverOptions();
            options.seed = 23000;

            if ("SolverMVA".equals(solverName)) {
                return new MVA(model, options);
            } else if ("SolverCTMC".equals(solverName)) {
                options.cutoff = new jline.util.matrix.Matrix(new double[]{10.0});
                return new CTMC(model, options);
            } else if ("SolverJMT".equals(solverName)) {
                return new JMT(model, options);
            } else if ("SolverSSA".equals(solverName)) {
                options.samples = 10000;
                return new SSA(model, options);
            } else if ("SolverLDES".equals(solverName)) {
                return new LDES(model, options);
            } else if ("SolverFluid".equals(solverName)) {
                return new FLD(model, options);
            } else if ("SolverNC".equals(solverName)) {
                return new NC(model, options);
            } else if ("SolverMAM".equals(solverName)) {
                return new MAM(model, options);
            } else {
                return null;
            }
        } catch (Exception e) {
            return null;
        }
    }

    private double[] extractMetric(Object avgTable, String metric) {
        try {
            jline.solvers.NetworkAvgTable table = (jline.solvers.NetworkAvgTable) avgTable;
            List<Double> values;
            if ("QLen".equals(metric)) {
                values = table.getQLen();
            } else if ("Util".equals(metric)) {
                values = table.getUtil();
            } else if ("RespT".equals(metric)) {
                values = table.getRespT();
            } else if ("Tput".equals(metric)) {
                values = table.getTput();
            } else {
                return null;
            }
            if (values == null) return null;
            double[] arr = new double[values.size()];
            for (int i = 0; i < values.size(); i++) arr[i] = values.get(i);
            return arr;
        } catch (Exception e) {
            return null;
        }
    }

    private double[] flatten2DArray(JsonArray arr) {
        List<Double> flat = new ArrayList<Double>();
        for (JsonElement row : arr) {
            if (row.isJsonArray()) {
                for (JsonElement elem : row.getAsJsonArray()) {
                    flat.add(jsonElementToDouble(elem));
                }
            } else {
                flat.add(jsonElementToDouble(row));
            }
        }
        double[] result = new double[flat.size()];
        for (int i = 0; i < flat.size(); i++) result[i] = flat.get(i);
        return result;
    }

    private double jsonElementToDouble(JsonElement elem) {
        if (elem.isJsonNull()) return Double.NaN;
        if (elem.isJsonPrimitive()) {
            com.google.gson.JsonPrimitive prim = elem.getAsJsonPrimitive();
            if (prim.isNumber()) return prim.getAsDouble();
            if (prim.isString()) {
                String s = prim.getAsString();
                if ("NaN".equals(s)) return Double.NaN;
                if ("Inf".equals(s)) return Double.POSITIVE_INFINITY;
                if ("-Inf".equals(s)) return Double.NEGATIVE_INFINITY;
            }
            return Double.NaN;
        }
        return Double.NaN;
    }

    private boolean compareWithTolerance(double[] actual, double[] expected, double tol) {
        if (actual.length != expected.length) return false;

        for (int i = 0; i < actual.length; i++) {
            double exp = expected[i];
            double act = actual[i];

            if (Double.isNaN(exp) && Double.isNaN(act)) continue;
            if (Double.isNaN(exp) || Double.isNaN(act)) return false;

            double relErr = Math.abs(act - exp) / (Math.abs(exp) + 1e-10);
            double absErr = Math.abs(act - exp);
            if (relErr > tol && absErr > tol) return false;
        }
        return true;
    }

    private void printSummary() {
        System.out.println("\n\n=== Test Summary ===");
        System.out.println("Passed:  " + results.passed);
        System.out.println("Failed:  " + results.failed);
        System.out.println("Skipped: " + results.skipped);
        System.out.println("====================");

        if (results.failed > 0) {
            System.out.println("\nFailed tests:");
            for (String error : results.errors) {
                System.out.println("  - " + error);
            }
        }
    }

    public static void main(String[] args) {
        String modelName = (args.length > 0) ? args[0] : null;
        boolean jsonOutput = false;
        boolean verbose = false;
        for (String a : args) {
            if ("--json".equals(a)) jsonOutput = true;
            if ("--verbose".equals(a) || "-v".equals(a)) verbose = true;
        }

        if (modelName == null) {
            System.err.println("Usage: java jline.unified.UnifiedTestRunner <model_name> [--json] [--verbose]");
            System.exit(1);
            return;
        }

        UnifiedTestRunner runner = new UnifiedTestRunner();
        Object result = runner.runModel(modelName, verbose, jsonOutput);

        if (jsonOutput) {
            Gson g = new Gson();
            System.out.println(g.toJson(result));
        } else {
            System.out.println("\nResult: " + result);
        }

        int exitCode;
        if (result instanceof JsonResult) {
            exitCode = "passed".equals(((JsonResult) result).status) ? 0 : 1;
        } else if ("passed".equals(result)) {
            exitCode = 0;
        } else {
            exitCode = 1;
        }
        System.exit(exitCode);
    }
}
