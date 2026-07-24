package jline.api.wf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.lang.Network;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.auto.ModelAnalyzer;
import jline.solvers.auto.SolverAUTO;

/**
 * Integration class connecting workflow analysis to the AUTO solver.
 */
public class Wf_auto_integration {
    private final Network network;
    private final SolverOptions options;
    private final Wf_analyzer workflowAnalyzer;
    private final ModelAnalyzer modelAnalyzer;

    public Wf_auto_integration(Network network) { this(network, new SolverOptions()); }

    public Wf_auto_integration(Network network, SolverOptions options) {
        this.network = network;
        this.options = options;
        this.workflowAnalyzer = new Wf_analyzer(network);
        this.modelAnalyzer = new ModelAnalyzer(network);
    }

    public static final class ExtendedSolverRecommendation {
        public final String recommendedSolver;
        public final double confidence;
        public final List<String> reasoning;
        public final Map<String, Object> workflowFeatures;
        public final List<String> alternativeSolvers;
        public ExtendedSolverRecommendation(String recommendedSolver, double confidence,
                                            List<String> reasoning, Map<String, Object> workflowFeatures,
                                            List<String> alternativeSolvers) {
            this.recommendedSolver = recommendedSolver;
            this.confidence = confidence;
            this.reasoning = reasoning;
            this.workflowFeatures = workflowFeatures;
            this.alternativeSolvers = alternativeSolvers;
        }
        public String getRecommendedSolver() { return recommendedSolver; }
        public double getConfidence() { return confidence; }
        public List<String> getReasoning() { return reasoning; }
        public Map<String, Object> getWorkflowFeatures() { return workflowFeatures; }
        public List<String> getAlternativeSolvers() { return alternativeSolvers; }
    }

    public ExtendedSolverRecommendation recommendSolverWithWorkflowAnalysis() {
        Wf_analyzer.WorkflowAnalysis workflowAnalysis = workflowAnalyzer.analyzeWorkflow();
        Map<String, Object> workflowFeatures = extractWorkflowFeatures(workflowAnalysis);
        String baseRecommendation = getBaseRecommendation();
        return enhanceRecommendationWithWorkflow(baseRecommendation, workflowFeatures, workflowAnalysis);
    }

    private Map<String, Object> extractWorkflowFeatures(Wf_analyzer.WorkflowAnalysis analysis) {
        Map<String, Object> features = new HashMap<String, Object>();
        Wf_analyzer.DetectedPatterns patterns = analysis.detectedPatterns;
        Map<?, ?> originalComplexity = (Map<?, ?>) analysis.statistics.get("originalComplexity");
        Map<?, ?> optimizedComplexity = (Map<?, ?>) analysis.statistics.get("optimizedComplexity");

        features.put("hasSequencePatterns", !patterns.sequences.isEmpty());
        features.put("hasParallelPatterns", !patterns.parallels.isEmpty());
        features.put("hasLoopPatterns", !patterns.loops.isEmpty());
        features.put("hasBranchPatterns", !patterns.branches.isEmpty());
        features.put("numSequences", patterns.sequences.size());
        features.put("numParallels", patterns.parallels.size());
        features.put("numLoops", patterns.loops.size());
        features.put("numBranches", patterns.branches.size());

        features.put("originalNodeCount", originalComplexity != null ? originalComplexity.get("totalNodes") : 0);
        features.put("originalLinkCount", originalComplexity != null ? originalComplexity.get("totalLinks") : 0);
        features.put("optimizedNodeCount", optimizedComplexity != null ? optimizedComplexity.get("totalNodes") : 0);
        features.put("optimizedLinkCount", optimizedComplexity != null ? optimizedComplexity.get("totalLinks") : 0);

        if (!patterns.sequences.isEmpty()) {
            double total = 0.0;
            int max = 0;
            for (List<?> seq : patterns.sequences) {
                total += seq.size();
                if (seq.size() > max) max = seq.size();
            }
            features.put("avgSequenceLength", total / patterns.sequences.size());
            features.put("maxSequenceLength", max);
        }
        if (!patterns.parallels.isEmpty()) {
            double total = 0.0;
            int max = 0;
            for (List<?> par : patterns.parallels) {
                total += par.size();
                if (par.size() > max) max = par.size();
            }
            features.put("avgParallelism", total / patterns.parallels.size());
            features.put("maxParallelism", max);
        }
        if (!patterns.loops.isEmpty()) {
            Map<?, ?> loopStats = (Map<?, ?>) analysis.statistics.get("loopStats");
            features.put("avgLoopProbability", loopStats != null ? loopStats.get("avgLoopProbability") : 0.0);
            features.put("maxLoopProbability", loopStats != null ? loopStats.get("maxLoopProbability") : 0.0);
        }
        if (!patterns.branches.isEmpty()) {
            Map<?, ?> branchStats = (Map<?, ?>) analysis.statistics.get("branchStats");
            features.put("avgBranches", branchStats != null ? branchStats.get("avgBranches") : 0.0);
            features.put("maxBranches", branchStats != null ? branchStats.get("maxBranches") : 0);

            double totalEntropy = 0.0;
            for (Wf_branch_detector.BranchPattern branch : patterns.branches) {
                Map<String, Double> diversity = Wf_branch_detector.calculateBranchDiversity(branch);
                Object e = diversity.get("entropy");
                totalEntropy += (e instanceof Double) ? (Double) e : 0.0;
            }
            features.put("avgBranchEntropy", totalEntropy / patterns.branches.size());
        }
        return features;
    }

    private String getBaseRecommendation() {
        boolean hasProductForm = modelAnalyzer.hasProductForm();
        boolean hasSingleChain = modelAnalyzer.hasSingleChain();
        boolean hasMultiChain = modelAnalyzer.hasMultiChain();
        int totalJobs = modelAnalyzer.getTotalJobs();
        double avgJobsPerChain = modelAnalyzer.getAvgJobsPerChain();
        if (hasSingleChain) return "NC";
        if (hasMultiChain && hasProductForm && totalJobs < 10) return "NC";
        if (hasMultiChain && hasProductForm && avgJobsPerChain < 30) return "MVA";
        if (hasMultiChain && avgJobsPerChain > 30) return "FLUID";
        return "MVA";
    }

    private ExtendedSolverRecommendation enhanceRecommendationWithWorkflow(String baseRecommendation,
                                                                           Map<String, Object> workflowFeatures,
                                                                           Wf_analyzer.WorkflowAnalysis analysis) {
        List<String> reasoning = new ArrayList<String>();
        List<String> alternatives = new ArrayList<String>();
        String finalRecommendation = baseRecommendation;
        double confidence = 0.7;

        boolean hasSequences = (Boolean) workflowFeatures.get("hasSequencePatterns");
        boolean hasParallels = (Boolean) workflowFeatures.get("hasParallelPatterns");
        boolean hasLoops = (Boolean) workflowFeatures.get("hasLoopPatterns");
        boolean hasBranches = (Boolean) workflowFeatures.get("hasBranchPatterns");

        if (hasSequences) {
            reasoning.add("Detected sequence patterns - suitable for analytical methods");
            confidence += 0.1;
            Object maxSeqLengthObj = workflowFeatures.get("maxSequenceLength");
            int maxSeqLength = (maxSeqLengthObj instanceof Integer) ? (Integer) maxSeqLengthObj : 0;
            if (maxSeqLength > 10) {
                reasoning.add("Long sequences detected - consider FLUID approximation");
                if ("MVA".equals(finalRecommendation)) alternatives.add("FLUID");
            }
        }
        if (hasParallels) {
            reasoning.add("Detected parallel patterns - fork-join structures present");
            Object maxParObj = workflowFeatures.get("maxParallelism");
            int maxParallelism = (maxParObj instanceof Integer) ? (Integer) maxParObj : 0;
            if (maxParallelism > 5) {
                reasoning.add("High parallelism detected - exact methods may be computationally expensive");
                if ("NC".equals(finalRecommendation) || "MVA".equals(finalRecommendation)) {
                    finalRecommendation = "SSA";
                    reasoning.add("Switching to SSA for high-parallelism workflow");
                }
                alternatives.add("JMT");
            } else {
                confidence += 0.05;
            }
        }
        if (hasLoops) {
            Object avgLoopProbObj = workflowFeatures.get("avgLoopProbability");
            double avgLoopProb = (avgLoopProbObj instanceof Double) ? (Double) avgLoopProbObj : 0.0;
            Object maxLoopProbObj = workflowFeatures.get("maxLoopProbability");
            double maxLoopProb = (maxLoopProbObj instanceof Double) ? (Double) maxLoopProbObj : 0.0;
            reasoning.add("Detected loop patterns with avg probability " + avgLoopProb);
            if (maxLoopProb > 0.8) {
                reasoning.add("High loop probability detected - may cause numerical instability");
                confidence -= 0.1;
                if ("NC".equals(finalRecommendation) || "MVA".equals(finalRecommendation)) {
                    alternatives.add("SSA");
                    alternatives.add("FLUID");
                }
            } else if (maxLoopProb > 0.5) {
                reasoning.add("Moderate loop probability - analytical methods suitable");
                confidence += 0.05;
            }
        }
        if (hasBranches) {
            Object avgEntropyObj = workflowFeatures.get("avgBranchEntropy");
            double avgEntropy = (avgEntropyObj instanceof Double) ? (Double) avgEntropyObj : 0.0;
            Object maxBranchesObj = workflowFeatures.get("maxBranches");
            int maxBranches = (maxBranchesObj instanceof Integer) ? (Integer) maxBranchesObj : 0;
            reasoning.add("Detected branch patterns with avg entropy " + avgEntropy);
            if (avgEntropy > 1.5) {
                reasoning.add("High branching entropy - complex decision structure");
                if (maxBranches > 5) {
                    reasoning.add("Many branches detected - consider simulation methods");
                    alternatives.add("JMT");
                    alternatives.add("SSA");
                }
            }
            if (avgEntropy < 0.5) {
                reasoning.add("Low branching entropy - deterministic-like behavior");
                confidence += 0.1;
            }
        }

        Object originalNodesObj = workflowFeatures.get("originalNodeCount");
        int originalNodes = (originalNodesObj instanceof Integer) ? (Integer) originalNodesObj : 0;
        Object optimizedNodesObj = workflowFeatures.get("optimizedNodeCount");
        int optimizedNodes = (optimizedNodesObj instanceof Integer) ? (Integer) optimizedNodesObj : 0;

        if (originalNodes > optimizedNodes) {
            double reduction = (double) (originalNodes - optimizedNodes) / originalNodes;
            reasoning.add("Workflow complexity reduced by " + ((int) (reduction * 100)) + "% through pattern optimization");
            confidence += 0.1;
        }
        if (optimizedNodes > 50) {
            reasoning.add("Large optimized workflow - consider approximation methods");
            if ("NC".equals(finalRecommendation)) {
                finalRecommendation = "MVA";
                reasoning.add("Switching from NC to MVA for large workflow");
            }
            alternatives.add("FLUID");
        }
        confidence = Math.min(1.0, Math.max(0.1, confidence));

        List<String> standardSolvers = Arrays.asList("MVA", "NC", "SSA", "FLUID", "JMT", "CTMC");
        for (String solver : standardSolvers) {
            if (!solver.equals(finalRecommendation) && !alternatives.contains(solver)) {
                alternatives.add(solver);
            }
        }
        List<String> top3 = new ArrayList<String>(alternatives.subList(0, Math.min(3, alternatives.size())));
        return new ExtendedSolverRecommendation(finalRecommendation, confidence, reasoning, workflowFeatures, top3);
    }

    public NetworkSolver createOptimalSolver() {
        ExtendedSolverRecommendation recommendation = recommendSolverWithWorkflowAnalysis();
        String r = recommendation.recommendedSolver.toUpperCase();
        if ("MVA".equals(r)) return new jline.solvers.mva.SolverMVA(network, options);
        if ("NC".equals(r)) return new jline.solvers.nc.SolverNC(network, options);
        if ("SSA".equals(r)) return new jline.solvers.ssa.SolverSSA(network, options);
        if ("FLUID".equals(r)) return new jline.solvers.fluid.SolverFluid(network, options);
        if ("JMT".equals(r)) return new jline.solvers.wrappers.jmt.SolverJMT(network, options);
        if ("CTMC".equals(r)) return new jline.solvers.ctmc.SolverCTMC(network, options);
        return new SolverAUTO(network, options);
    }

    public Map<String, Object> getOptimizationInsights() {
        Wf_analyzer.WorkflowAnalysis analysis = workflowAnalyzer.analyzeWorkflow();
        Map<String, Object> insights = new HashMap<String, Object>();
        insights.put("recommendations", workflowAnalyzer.getOptimizationRecommendations(analysis));
        insights.put("patternInsights", generatePatternInsights(analysis.detectedPatterns));
        insights.put("performancePredictions", generatePerformancePredictions(analysis));
        return insights;
    }

    private Map<String, Object> generatePatternInsights(Wf_analyzer.DetectedPatterns patterns) {
        Map<String, Object> insights = new HashMap<String, Object>();
        if (!patterns.sequences.isEmpty()) {
            insights.put("sequenceOptimization", "Consider merging sequential services to reduce overhead");
        }
        if (!patterns.parallels.isEmpty()) {
            insights.put("parallelOptimization", "Parallel patterns can benefit from resource pooling strategies");
        }
        if (!patterns.loops.isEmpty()) {
            int highProbLoops = patterns.loops.size();
            if (highProbLoops > 0) {
                insights.put("loopOptimization", "High-probability loops may benefit from caching or memoization");
            }
        }
        if (!patterns.branches.isEmpty()) {
            double total = 0.0;
            for (Object b : patterns.branches) {
                try {
                    java.lang.reflect.Field f = b.getClass().getDeclaredField("branchNodes");
                    f.setAccessible(true);
                    Object branchNodes = f.get(b);
                    if (branchNodes instanceof List<?>) total += ((List<?>) branchNodes).size();
                } catch (Exception e) {
                    // ignore
                }
            }
            double avgBranches = total / patterns.branches.size();
            if (avgBranches > 3) {
                insights.put("branchOptimization", "Complex branching patterns - consider load balancing strategies");
            }
        }
        return insights;
    }

    private Map<String, Object> generatePerformancePredictions(Wf_analyzer.WorkflowAnalysis analysis) {
        Map<String, Object> predictions = new HashMap<String, Object>();
        Map<?, ?> updateStats = (Map<?, ?>) analysis.statistics.get("updateStats");
        Object reductionObj = (updateStats != null) ? updateStats.get("reductionRatio") : 0.0;
        double reductionRatio = (reductionObj instanceof Double) ? (Double) reductionObj : 0.0;
        if (reductionRatio > 0.1) {
            predictions.put("complexityReduction", "Expected " + ((int) (reductionRatio * 100)) + "% reduction in solve time");
        }
        Wf_analyzer.DetectedPatterns patterns = analysis.detectedPatterns;
        if (!patterns.parallels.isEmpty()) {
            predictions.put("parallelizationPotential", "High potential for parallel execution optimization");
        }
        if (!patterns.loops.isEmpty()) {
            predictions.put("convergenceConsiderations", "Loop patterns may affect solver convergence rates");
        }
        return predictions;
    }

    public boolean validateWorkflowEnhancement() {
        try {
            ExtendedSolverRecommendation recommendation = recommendSolverWithWorkflowAnalysis();
            Wf_analyzer.WorkflowAnalysis analysis = workflowAnalyzer.analyzeWorkflow();
            return recommendation.confidence > 0.1
                    && !recommendation.reasoning.isEmpty()
                    && workflowAnalyzer.validateAnalysis(analysis);
        } catch (Exception e) {
            return false;
        }
    }
}
