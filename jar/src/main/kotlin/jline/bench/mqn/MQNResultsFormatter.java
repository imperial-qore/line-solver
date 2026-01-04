/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench.mqn;

import java.util.List;
import java.util.Map;
import java.util.ArrayList;

/**
 * Formats MQN benchmark results in MATLAB-style output format
 */
public class MQNResultsFormatter {

    private static final String[] SOLVER_ORDER = {"Fluid", "MVA", "NC", "Auto", "QNS", "MAM"};
    private static final int BENCHMARK_NAME_WIDTH = 45;

    /**
     * Benchmark result data structure
     */
    public static class MQNBenchmarkResult {
        public String benchmarkName;
        public Map<String, Double> errQ;
        public Map<String, Double> errW;
        public Map<String, Double> errU;
        public Map<String, Double> errT;

        @SuppressWarnings("unchecked")
        public MQNBenchmarkResult(String benchmarkName, Map<String, Object> results) {
            this.benchmarkName = benchmarkName;
            this.errQ = (Map<String, Double>) results.get("ERRQ");
            this.errW = (Map<String, Double>) results.get("ERRR");
            this.errU = (Map<String, Double>) results.get("ERRU");
            this.errT = (Map<String, Double>) results.get("ERRT");
        }
    }

    /**
     * Accumulates benchmark results for batch formatting
     */
    public static class MQNResultsAccumulator {
        private List<MQNBenchmarkResult> results = new ArrayList<>();
        private boolean headerPrinted = false;

        public void addResult(String benchmarkName, Map<String, Object> benchmarkResults) {
            results.add(new MQNBenchmarkResult(benchmarkName, benchmarkResults));
        }

        public void printResults() {
            if (!headerPrinted) {
                headerPrinted = true;
            }

            // Group results by benchmark name and average across iterations
            Map<String, List<MQNBenchmarkResult>> groupedResults = new java.util.LinkedHashMap<>();
            for (MQNBenchmarkResult result : results) {
                // Remove iteration suffix (e.g., "_1", "_2", etc.) from benchmark name
                String baseName = result.benchmarkName.replaceAll("_\\d+$", "");
                groupedResults.computeIfAbsent(baseName, k -> new ArrayList<>()).add(result);
            }

            // Print averaged results for each benchmark
            for (Map.Entry<String, List<MQNBenchmarkResult>> entry : groupedResults.entrySet()) {
                String benchmarkName = entry.getKey();
                List<MQNBenchmarkResult> iterations = entry.getValue();
                printAveragedBenchmarkResult(benchmarkName, iterations);
            }

            if (!results.isEmpty()) {
                printSummaryTable();
                printLegend();
            }
        }

        public void clear() {
            results.clear();
        }

        private void printAveragedBenchmarkResult(String benchmarkName, List<MQNBenchmarkResult> iterations) {
            // Average error values across all iterations
            Map<String, Double> avgErrQ = averageErrors(iterations, r -> r.errQ);
            Map<String, Double> avgErrW = averageErrors(iterations, r -> r.errW);
            Map<String, Double> avgErrU = averageErrors(iterations, r -> r.errU);
            Map<String, Double> avgErrT = averageErrors(iterations, r -> r.errT);

            StringBuilder line = new StringBuilder();
            String paddedName = String.format("%-" + BENCHMARK_NAME_WIDTH + "s", benchmarkName);
            line.append(paddedName).append(": ");

            line.append("ErrQ: ").append(formatErrorValues(avgErrQ)).append(" ");
            line.append("ErrW: ").append(formatErrorValues(avgErrW)).append(" ");
            line.append("ErrU: ").append(formatErrorValues(avgErrU)).append(" ");
            line.append("ErrT: ").append(formatErrorValues(avgErrT));

            System.out.println(line.toString());
        }

        private Map<String, Double> averageErrors(List<MQNBenchmarkResult> iterations,
                                                   java.util.function.Function<MQNBenchmarkResult, Map<String, Double>> errorGetter) {
            Map<String, Double> avgErrors = new java.util.HashMap<>();
            if (iterations.isEmpty()) return avgErrors;

            // For each solver, compute the mean across all iterations
            for (String solverName : SOLVER_ORDER) {
                String solverKey = mapSolverName(solverName);
                double sum = 0.0;
                int count = 0;

                for (MQNBenchmarkResult result : iterations) {
                    Map<String, Double> errors = errorGetter.apply(result);
                    if (errors != null && errors.containsKey(solverKey)) {
                        Double error = errors.get(solverKey);
                        if (error != null && error != Double.MAX_VALUE) {
                            sum += error;
                            count++;
                        }
                    }
                }

                if (count > 0) {
                    avgErrors.put(solverKey, sum / count);
                } else {
                    avgErrors.put(solverKey, Double.NaN);
                }
            }

            return avgErrors;
        }

        private String formatErrorValues(Map<String, Double> errors) {
            StringBuilder values = new StringBuilder();

            for (int i = 0; i < SOLVER_ORDER.length; i++) {
                String solverKey = mapSolverName(SOLVER_ORDER[i]);
                Double error = errors != null ? errors.get(solverKey) : null;

                if (error != null && error != Double.MAX_VALUE && !error.isNaN()) {
                    values.append(String.format("%.3f", error));
                } else {
                    values.append("  NaN");
                }

                // Always add trailing space after each value (MATLAB format)
                values.append(" ");
            }

            return values.toString();
        }

        private String mapSolverName(String displayName) {
            switch (displayName) {
                case "Fluid": return "SolverFluid";
                case "MVA": return "SolverMVA";
                case "NC": return "SolverNC";
                case "Auto": return "SolverAUTO";
                case "QNS": return "SolverQNS";
                case "MAM": return "SolverMAM";
                default: return displayName;
            }
        }

        private void printSummaryTable() {
            if (results.isEmpty()) return;

            // Compute average errors per solver across all benchmarks
            Map<String, double[]> solverErrors = new java.util.LinkedHashMap<>();
            for (String solver : SOLVER_ORDER) {
                solverErrors.put(solver, new double[]{0.0, 0.0, 0.0, 0.0}); // QLen, RespT, Util, Tput
            }

            // Count valid entries for averaging
            Map<String, int[]> solverCounts = new java.util.LinkedHashMap<>();
            for (String solver : SOLVER_ORDER) {
                solverCounts.put(solver, new int[]{0, 0, 0, 0});
            }

            // Accumulate errors from all results
            for (MQNBenchmarkResult result : results) {
                for (String solver : SOLVER_ORDER) {
                    String solverKey = mapSolverName(solver);
                    double[] errors = solverErrors.get(solver);
                    int[] counts = solverCounts.get(solver);

                    // ErrQ -> QLen (index 0)
                    if (result.errQ != null && result.errQ.containsKey(solverKey)) {
                        Double val = result.errQ.get(solverKey);
                        if (val != null && val != Double.MAX_VALUE && !val.isNaN()) {
                            errors[0] += val;
                            counts[0]++;
                        }
                    }
                    // ErrW -> RespT (index 1)
                    if (result.errW != null && result.errW.containsKey(solverKey)) {
                        Double val = result.errW.get(solverKey);
                        if (val != null && val != Double.MAX_VALUE && !val.isNaN()) {
                            errors[1] += val;
                            counts[1]++;
                        }
                    }
                    // ErrU -> Util (index 2)
                    if (result.errU != null && result.errU.containsKey(solverKey)) {
                        Double val = result.errU.get(solverKey);
                        if (val != null && val != Double.MAX_VALUE && !val.isNaN()) {
                            errors[2] += val;
                            counts[2]++;
                        }
                    }
                    // ErrT -> Tput (index 3)
                    if (result.errT != null && result.errT.containsKey(solverKey)) {
                        Double val = result.errT.get(solverKey);
                        if (val != null && val != Double.MAX_VALUE && !val.isNaN()) {
                            errors[3] += val;
                            counts[3]++;
                        }
                    }
                }
            }

            // Compute averages
            for (String solver : SOLVER_ORDER) {
                double[] errors = solverErrors.get(solver);
                int[] counts = solverCounts.get(solver);
                for (int i = 0; i < 4; i++) {
                    if (counts[i] > 0) {
                        errors[i] /= counts[i];
                    } else {
                        errors[i] = Double.NaN;
                    }
                }
            }

            // Print summary table
            System.out.println("\n--- Summary Table (Average Errors) ---");
            System.out.println(String.format("%-10s %10s %10s %10s %10s", "", "QLen", "RespT", "Util", "Tput"));
            for (String solver : SOLVER_ORDER) {
                double[] errors = solverErrors.get(solver);
                System.out.println(String.format("%-10s %10s %10s %10s %10s",
                    solver, formatSummaryValue(errors[0]), formatSummaryValue(errors[1]),
                    formatSummaryValue(errors[2]), formatSummaryValue(errors[3])));
            }
        }

        private String formatSummaryValue(double value) {
            if (Double.isNaN(value)) {
                return "NaN";
            }
            return String.format("%.4f", value);
        }

        private void printLegend() {
            System.out.println("\nLegend: Fluid, MVA, NC, Auto, QNS, MAM");
        }
    }

    private static final MQNResultsAccumulator globalAccumulator = new MQNResultsAccumulator();

    public static void addGlobalResult(String benchmarkName, Map<String, Object> results) {
        globalAccumulator.addResult(benchmarkName, results);
    }

    public static void printAndClearGlobalResults() {
        globalAccumulator.printResults();
        globalAccumulator.clear();
    }

    public static void printImmediateResult(String benchmarkName, Map<String, Object> results) {
        MQNResultsAccumulator tempAccumulator = new MQNResultsAccumulator();
        tempAccumulator.addResult(benchmarkName, results);
        tempAccumulator.printResults();
    }

    public static MQNResultsAccumulator createAccumulator() {
        return new MQNResultsAccumulator();
    }
}
