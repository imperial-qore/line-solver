/**
 * @file Workflow Management System
 *
 * @since LINE 3.0
 */
package jline.api.wf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.lang.Network;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

/**
 * Main facade for workflow management and optimization in LINE.
 */
public class WorkflowManager {

    private final Network network;
    private final SolverOptions options;
    private final Wf_analyzer workflowAnalyzer;
    private final Wf_auto_integration autoIntegration;

    public WorkflowManager(Network network, SolverOptions options) {
        this.network = network;
        this.options = options;
        this.workflowAnalyzer = new Wf_analyzer(network);
        this.autoIntegration = new Wf_auto_integration(network, options);
    }

    public WorkflowManager(Network network) {
        this(network, new SolverOptions());
    }

    /** Comprehensive workflow analysis result. */
    public static final class WorkflowAnalysisResult {
        public final Wf_analyzer.WorkflowAnalysis patternAnalysis;
        public final Wf_auto_integration.ExtendedSolverRecommendation solverRecommendation;
        public final Map<String, Object> optimizationInsights;
        public final Map<String, Double> performanceMetrics;

        public WorkflowAnalysisResult(Wf_analyzer.WorkflowAnalysis patternAnalysis,
                                      Wf_auto_integration.ExtendedSolverRecommendation solverRecommendation,
                                      Map<String, Object> optimizationInsights,
                                      Map<String, Double> performanceMetrics) {
            this.patternAnalysis = patternAnalysis;
            this.solverRecommendation = solverRecommendation;
            this.optimizationInsights = optimizationInsights;
            this.performanceMetrics = performanceMetrics;
        }

        public Wf_analyzer.WorkflowAnalysis getPatternAnalysis() { return patternAnalysis; }
        public Wf_auto_integration.ExtendedSolverRecommendation getSolverRecommendation() { return solverRecommendation; }
        public Map<String, Object> getOptimizationInsights() { return optimizationInsights; }
        public Map<String, Double> getPerformanceMetrics() { return performanceMetrics; }
    }

    public WorkflowAnalysisResult analyzeWorkflow() {
        Wf_analyzer.WorkflowAnalysis patternAnalysis = workflowAnalyzer.analyzeWorkflow();
        Wf_auto_integration.ExtendedSolverRecommendation solverRecommendation =
                autoIntegration.recommendSolverWithWorkflowAnalysis();
        Map<String, Object> optimizationInsights = autoIntegration.getOptimizationInsights();
        Map<String, Double> performanceMetrics = calculatePerformanceMetrics(patternAnalysis, solverRecommendation);
        return new WorkflowAnalysisResult(patternAnalysis, solverRecommendation,
                optimizationInsights, performanceMetrics);
    }

    public NetworkSolver createOptimizedSolver() {
        return autoIntegration.createOptimalSolver();
    }

    public Wf_analyzer.DetectedPatterns getPatternAnalysis() {
        Wf_analyzer.WorkflowAnalysis analysis = workflowAnalyzer.analyzeWorkflow();
        return analysis.detectedPatterns;
    }

    public List<String> getOptimizationRecommendations() {
        Wf_analyzer.WorkflowAnalysis analysis = workflowAnalyzer.analyzeWorkflow();
        return workflowAnalyzer.getOptimizationRecommendations(analysis);
    }

    public Map<String, Object> generateComplexityReport() {
        Wf_analyzer.WorkflowAnalysis analysis = workflowAnalyzer.analyzeWorkflow();
        Map<String, Object> report = new HashMap<String, Object>();

        Object originalComplexityObj = analysis.statistics.get("originalComplexity");
        Object optimizedComplexityObj = analysis.statistics.get("optimizedComplexity");
        @SuppressWarnings("unchecked")
        Map<String, Object> originalComplexity = originalComplexityObj instanceof Map
                ? (Map<String, Object>) originalComplexityObj : null;
        @SuppressWarnings("unchecked")
        Map<String, Object> optimizedComplexity = optimizedComplexityObj instanceof Map
                ? (Map<String, Object>) optimizedComplexityObj : null;

        report.put("originalMetrics", originalComplexity != null ? originalComplexity : new HashMap<String, Object>());
        report.put("optimizedMetrics", optimizedComplexity != null ? optimizedComplexity : new HashMap<String, Object>());

        Wf_analyzer.DetectedPatterns patterns = analysis.detectedPatterns;
        Map<String, Object> patternComplexity = new HashMap<String, Object>();

        Map<String, Object> seqMap = new HashMap<String, Object>();
        seqMap.put("count", patterns.sequences.size());
        int seqTotalNodes = 0;
        int seqComplexity = 0;
        for (List<Integer> s : patterns.sequences) {
            seqTotalNodes += s.size();
            seqComplexity += s.size() * s.size();
        }
        seqMap.put("totalNodes", seqTotalNodes);
        seqMap.put("complexity", seqComplexity);
        patternComplexity.put("sequences", seqMap);

        Map<String, Object> parMap = new HashMap<String, Object>();
        parMap.put("count", patterns.parallels.size());
        int parTotalNodes = 0;
        int parComplexity = 0;
        for (List<Integer> p : patterns.parallels) {
            parTotalNodes += p.size();
            parComplexity += factorial(p.size());
        }
        parMap.put("totalNodes", parTotalNodes);
        parMap.put("complexity", parComplexity);
        patternComplexity.put("parallels", parMap);

        Map<String, Object> loopMap = new HashMap<String, Object>();
        loopMap.put("count", patterns.loops.size());
        loopMap.put("complexity", patterns.loops.size() * 10);
        patternComplexity.put("loops", loopMap);

        Map<String, Object> branchMap = new HashMap<String, Object>();
        branchMap.put("count", patterns.branches.size());
        int totalBranches = 0;
        int branchComplexity = 0;
        for (Wf_branch_detector.BranchPattern bp : patterns.branches) {
            totalBranches += bp.getBranchNodes().size();
            branchComplexity += bp.getBranchNodes().size();
        }
        branchMap.put("totalBranches", totalBranches);
        branchMap.put("complexity", branchComplexity);
        patternComplexity.put("branches", branchMap);

        report.put("patternComplexity", patternComplexity);

        int originalNodes = 0;
        int originalLinks = 0;
        if (originalComplexity != null) {
            Object n = originalComplexity.get("totalNodes");
            Object l = originalComplexity.get("totalLinks");
            if (n instanceof Number) originalNodes = ((Number) n).intValue();
            if (l instanceof Number) originalLinks = ((Number) l).intValue();
        }
        double complexityScore = calculateComplexityScore(originalNodes, originalLinks, patterns);
        report.put("overallComplexityScore", complexityScore);

        String level;
        if (complexityScore < 10) level = "LOW";
        else if (complexityScore < 50) level = "MEDIUM";
        else if (complexityScore < 100) level = "HIGH";
        else level = "VERY_HIGH";
        report.put("complexityLevel", level);

        return report;
    }

    public Map<String, Map<String, Object>> benchmarkSolvers(List<String> solvers) {
        Map<String, Map<String, Object>> results = new HashMap<String, Map<String, Object>>();
        for (String solverName : solvers) {
            try {
                NetworkSolver solver = createSolver(solverName);
                long startTime = System.currentTimeMillis();
                SolverResult solverResult = solver.getAvg();
                long endTime = System.currentTimeMillis();

                Map<String, Object> metrics = new HashMap<String, Object>();
                metrics.put("success", true);
                metrics.put("solveTime", endTime - startTime);
                metrics.put("hasResults", solverResult != null);

                if (solverResult != null && solverResult.QN != null) {
                    metrics.put("totalQueueLength", solverResult.QN.elementSum());
                    metrics.put("avgQueueLength", solverResult.QN.elementSum() / solverResult.QN.length());
                }

                results.put(solverName, metrics);
            } catch (Exception e) {
                Map<String, Object> metrics = new HashMap<String, Object>();
                metrics.put("success", false);
                metrics.put("error", e.getMessage() != null ? e.getMessage() : "Unknown error");
                results.put(solverName, metrics);
            }
        }
        return results;
    }

    public Map<String, Map<String, Object>> benchmarkSolvers() {
        return benchmarkSolvers(Arrays.asList("MVA", "NC", "SSA", "FLUID"));
    }

    public String exportAnalysis(String format) {
        WorkflowAnalysisResult analysis = analyzeWorkflow();
        String fmt = format.toUpperCase();
        if (fmt.equals("JSON")) return exportToJson(analysis);
        if (fmt.equals("CSV")) return exportToCsv(analysis);
        if (fmt.equals("SUMMARY")) return exportToSummary(analysis);
        throw new IllegalArgumentException("Unsupported export format: " + format);
    }

    public String exportAnalysis() {
        return exportAnalysis("SUMMARY");
    }

    public Map<String, Object> validateWorkflow() {
        Map<String, Object> validation = new HashMap<String, Object>();
        List<String> issues = new ArrayList<String>();
        try {
            if (network.getNodes().isEmpty()) {
                issues.add("Network has no nodes");
            }
            Wf_analyzer.WorkflowAnalysis analysis = workflowAnalyzer.analyzeWorkflow();
            if (!workflowAnalyzer.validateAnalysis(analysis)) {
                issues.add("Workflow analysis validation failed");
            }
            if (!autoIntegration.validateWorkflowEnhancement()) {
                issues.add("Workflow-enhanced solver selection validation failed");
            }
            validation.put("isValid", issues.isEmpty());
            validation.put("issues", issues);
            validation.put("nodeCount", network.getNodes().size());
            validation.put("hasRouting", !network.getLinkedRoutingMatrix().isEmpty());
        } catch (Exception e) {
            validation.put("isValid", false);
            List<String> errs = new ArrayList<String>();
            errs.add("Validation error: " + e.getMessage());
            validation.put("issues", errs);
        }
        return validation;
    }

    private Map<String, Double> calculatePerformanceMetrics(
            Wf_analyzer.WorkflowAnalysis analysis,
            Wf_auto_integration.ExtendedSolverRecommendation recommendation) {
        Map<String, Double> metrics = new HashMap<String, Double>();
        Object updateStatsObj = analysis.statistics.get("updateStats");
        @SuppressWarnings("unchecked")
        Map<String, Object> updateStats = updateStatsObj instanceof Map ? (Map<String, Object>) updateStatsObj : null;
        double reductionRatio = 0.0;
        if (updateStats != null) {
            Object rr = updateStats.get("reductionRatio");
            if (rr instanceof Number) reductionRatio = ((Number) rr).doubleValue();
        }
        metrics.put("complexityReduction", reductionRatio);
        metrics.put("solverConfidence", recommendation.confidence);

        Wf_analyzer.DetectedPatterns patterns = analysis.detectedPatterns;
        metrics.put("sequenceEfficiency", calculateSequenceEfficiency(patterns.sequences));
        metrics.put("parallelEfficiency", calculateParallelEfficiency(patterns.parallels));
        metrics.put("loopEfficiency", calculateLoopEfficiency(patterns.loops, analysis.originalWorkflow.linkMatrix));
        metrics.put("branchEfficiency", calculateBranchEfficiency(patterns.branches));
        return metrics;
    }

    private double calculateSequenceEfficiency(List<List<Integer>> sequences) {
        if (sequences.isEmpty()) return 1.0;
        int totalLength = 0;
        for (List<Integer> s : sequences) totalLength += s.size();
        double avgLength = (double) totalLength / sequences.size();
        return Math.min(1.0, avgLength / 5.0);
    }

    private double calculateParallelEfficiency(List<List<Integer>> parallels) {
        if (parallels.isEmpty()) return 1.0;
        int totalParallelism = 0;
        for (List<Integer> p : parallels) totalParallelism += p.size();
        double avgParallelism = (double) totalParallelism / parallels.size();
        return Math.max(0.1, 1.0 - (avgParallelism - 2.0) / 10.0);
    }

    private double calculateLoopEfficiency(List<Integer> loops, Matrix linkMatrix) {
        if (loops.isEmpty()) return 1.0;
        // Simplified: average probability of 0.5 per loop
        double avgProb = 0.5;
        return 1.0 - avgProb;
    }

    private double calculateBranchEfficiency(List<Wf_branch_detector.BranchPattern> branches) {
        if (branches.isEmpty()) return 1.0;
        double sum = 0.0;
        for (Wf_branch_detector.BranchPattern bp : branches) {
            Map<String, Double> diversity = Wf_branch_detector.calculateBranchDiversity(bp);
            Double e = diversity.get("normalizedEntropy");
            sum += e != null ? e.doubleValue() : 0.0;
        }
        return sum / branches.size();
    }

    private double calculateComplexityScore(int nodes, int links, Wf_analyzer.DetectedPatterns patterns) {
        double score = (double) nodes + links * 0.5;
        for (List<Integer> s : patterns.sequences) score += s.size() * s.size() * 0.1;
        for (List<Integer> p : patterns.parallels) score += factorial(p.size()) * 0.2;
        score += patterns.loops.size() * 10;
        for (Wf_branch_detector.BranchPattern bp : patterns.branches) score += bp.getBranchNodes().size() * 2;
        return score;
    }

    private static int factorial(int n) {
        if (n <= 1) return 1;
        return n * factorial(n - 1);
    }

    private NetworkSolver createSolver(String solverName) {
        String s = solverName.toUpperCase();
        if (s.equals("MVA")) return new jline.solvers.mva.SolverMVA(network, options);
        if (s.equals("NC")) return new jline.solvers.nc.SolverNC(network, options);
        if (s.equals("SSA")) return new jline.solvers.ssa.SolverSSA(network, options);
        if (s.equals("FLUID")) return new jline.solvers.fluid.SolverFluid(network, options);
        if (s.equals("JMT")) return new jline.solvers.wrappers.jmt.SolverJMT(network, options);
        if (s.equals("CTMC")) return new jline.solvers.ctmc.SolverCTMC(network, options);
        throw new IllegalArgumentException("Unknown solver: " + solverName);
    }

    private String exportToJson(WorkflowAnalysisResult analysis) {
        StringBuilder json = new StringBuilder();
        json.append("{\n");
        json.append("  \"solver_recommendation\": \"").append(analysis.solverRecommendation.recommendedSolver).append("\",\n");
        json.append("  \"confidence\": ").append(analysis.solverRecommendation.confidence).append(",\n");
        json.append("  \"patterns\": {\n");
        json.append("    \"sequences\": ").append(analysis.patternAnalysis.detectedPatterns.sequences.size()).append(",\n");
        json.append("    \"parallels\": ").append(analysis.patternAnalysis.detectedPatterns.parallels.size()).append(",\n");
        json.append("    \"loops\": ").append(analysis.patternAnalysis.detectedPatterns.loops.size()).append(",\n");
        json.append("    \"branches\": ").append(analysis.patternAnalysis.detectedPatterns.branches.size()).append("\n");
        json.append("  }\n");
        json.append("}");
        return json.toString();
    }

    private String exportToCsv(WorkflowAnalysisResult analysis) {
        StringBuilder csv = new StringBuilder();
        csv.append("Metric,Value\n");
        csv.append("Recommended Solver,").append(analysis.solverRecommendation.recommendedSolver).append("\n");
        csv.append("Confidence,").append(analysis.solverRecommendation.confidence).append("\n");
        csv.append("Sequences,").append(analysis.patternAnalysis.detectedPatterns.sequences.size()).append("\n");
        csv.append("Parallels,").append(analysis.patternAnalysis.detectedPatterns.parallels.size()).append("\n");
        csv.append("Loops,").append(analysis.patternAnalysis.detectedPatterns.loops.size()).append("\n");
        csv.append("Branches,").append(analysis.patternAnalysis.detectedPatterns.branches.size()).append("\n");
        return csv.toString();
    }

    private String exportToSummary(WorkflowAnalysisResult analysis) {
        StringBuilder summary = new StringBuilder();
        summary.append("=== Workflow Analysis Summary ===\n\n");
        summary.append("Recommended Solver: ").append(analysis.solverRecommendation.recommendedSolver).append("\n");
        summary.append("Confidence: ").append(String.format("%.2f", analysis.solverRecommendation.confidence)).append("\n\n");
        summary.append("Detected Patterns:\n");
        summary.append("- Sequences: ").append(analysis.patternAnalysis.detectedPatterns.sequences.size()).append("\n");
        summary.append("- Parallels: ").append(analysis.patternAnalysis.detectedPatterns.parallels.size()).append("\n");
        summary.append("- Loops: ").append(analysis.patternAnalysis.detectedPatterns.loops.size()).append("\n");
        summary.append("- Branches: ").append(analysis.patternAnalysis.detectedPatterns.branches.size()).append("\n\n");
        summary.append("Reasoning:\n");
        int idx = 0;
        for (String reason : analysis.solverRecommendation.reasoning) {
            summary.append(idx + 1).append(". ").append(reason).append("\n");
            idx++;
        }
        return summary.toString();
    }

    /** Quick analysis method for simple workflow inspection. */
    public static String quickAnalysis(Network network) {
        WorkflowManager manager = new WorkflowManager(network);
        return manager.exportAnalysis("SUMMARY");
    }

    /** Get optimal solver for a network without detailed analysis. */
    public static NetworkSolver getOptimalSolver(Network network, SolverOptions options) {
        WorkflowManager manager = new WorkflowManager(network, options);
        return manager.createOptimizedSolver();
    }

    public static NetworkSolver getOptimalSolver(Network network) {
        return getOptimalSolver(network, new SolverOptions());
    }
}
