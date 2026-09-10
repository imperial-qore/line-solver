/**
 * @file Workflow Analysis Engine
 *
 * @since LINE 3.0
 */
package jline.api.wf;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jline.GlobalConstants;
import jline.lang.Network;
import jline.lang.nodes.ClassSwitch;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Fork;
import jline.lang.nodes.Join;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Router;
import jline.util.matrix.Matrix;

/**
 * Main workflow analyzer that coordinates pattern detection and optimization.
 */
public class Wf_analyzer {

    private final Network network;

    public Wf_analyzer(Network network) {
        this.network = network;
    }

    /** Comprehensive workflow analysis result. */
    public static final class WorkflowAnalysis {
        public final WorkflowRepresentation originalWorkflow;
        public final DetectedPatterns detectedPatterns;
        public final Wf_pattern_updater.UpdatedWorkflow optimizedWorkflow;
        public final Map<String, Object> statistics;

        public WorkflowAnalysis(WorkflowRepresentation originalWorkflow,
                                DetectedPatterns detectedPatterns,
                                Wf_pattern_updater.UpdatedWorkflow optimizedWorkflow,
                                Map<String, Object> statistics) {
            this.originalWorkflow = originalWorkflow;
            this.detectedPatterns = detectedPatterns;
            this.optimizedWorkflow = optimizedWorkflow;
            this.statistics = statistics;
        }
        public WorkflowRepresentation getOriginalWorkflow() { return originalWorkflow; }
        public DetectedPatterns getDetectedPatterns() { return detectedPatterns; }
        public Wf_pattern_updater.UpdatedWorkflow getOptimizedWorkflow() { return optimizedWorkflow; }
        public Map<String, Object> getStatistics() { return statistics; }
    }

    /** Workflow in matrix form. */
    public static final class WorkflowRepresentation {
        public final Matrix linkMatrix;
        public final List<Integer> serviceNodes;
        public final List<Integer> forkNodes;
        public final List<Integer> joinNodes;
        public final List<Integer> routerNodes;
        public final Map<Integer, Wf_pattern_updater.ServiceParameters> serviceParameters;

        public WorkflowRepresentation(Matrix linkMatrix,
                                      List<Integer> serviceNodes,
                                      List<Integer> forkNodes,
                                      List<Integer> joinNodes,
                                      List<Integer> routerNodes,
                                      Map<Integer, Wf_pattern_updater.ServiceParameters> serviceParameters) {
            this.linkMatrix = linkMatrix;
            this.serviceNodes = serviceNodes;
            this.forkNodes = forkNodes;
            this.joinNodes = joinNodes;
            this.routerNodes = routerNodes;
            this.serviceParameters = serviceParameters;
        }

        public Matrix getLinkMatrix() { return linkMatrix; }
        public List<Integer> getServiceNodes() { return serviceNodes; }
        public List<Integer> getForkNodes() { return forkNodes; }
        public List<Integer> getJoinNodes() { return joinNodes; }
        public List<Integer> getRouterNodes() { return routerNodes; }
        public Map<Integer, Wf_pattern_updater.ServiceParameters> getServiceParameters() { return serviceParameters; }
    }

    /** All detected patterns. */
    public static final class DetectedPatterns {
        public final List<List<Integer>> sequences;
        public final List<List<Integer>> parallels;
        public final List<Integer> loops;
        public final List<Wf_branch_detector.BranchPattern> branches;

        public DetectedPatterns(List<List<Integer>> sequences,
                                List<List<Integer>> parallels,
                                List<Integer> loops,
                                List<Wf_branch_detector.BranchPattern> branches) {
            this.sequences = sequences;
            this.parallels = parallels;
            this.loops = loops;
            this.branches = branches;
        }

        public List<List<Integer>> getSequences() { return sequences; }
        public List<List<Integer>> getParallels() { return parallels; }
        public List<Integer> getLoops() { return loops; }
        public List<Wf_branch_detector.BranchPattern> getBranches() { return branches; }
    }

    /** Node classification result. */
    public static final class NodeClassification {
        public final List<Integer> serviceNodes;
        public final List<Integer> forkNodes;
        public final List<Integer> joinNodes;
        public final List<Integer> routerNodes;

        public NodeClassification(List<Integer> serviceNodes, List<Integer> forkNodes,
                                  List<Integer> joinNodes, List<Integer> routerNodes) {
            this.serviceNodes = serviceNodes;
            this.forkNodes = forkNodes;
            this.joinNodes = joinNodes;
            this.routerNodes = routerNodes;
        }
    }

    public WorkflowAnalysis analyzeWorkflow() {
        WorkflowRepresentation workflowRep = convertNetworkToWorkflow();
        DetectedPatterns patterns = detectAllPatterns(workflowRep);
        Wf_pattern_updater.UpdatedWorkflow optimizedWorkflow = Wf_pattern_updater.updatePatterns(
                workflowRep.linkMatrix,
                workflowRep.serviceNodes,
                workflowRep.forkNodes,
                workflowRep.joinNodes,
                workflowRep.routerNodes,
                workflowRep.serviceParameters);
        Map<String, Object> statistics = generateAnalysisStatistics(workflowRep, patterns, optimizedWorkflow);
        return new WorkflowAnalysis(workflowRep, patterns, optimizedWorkflow, statistics);
    }

    private WorkflowRepresentation convertNetworkToWorkflow() {
        List<Node> nodes = network.getNodes();
        Matrix linkMatrix = buildLinkMatrix();
        NodeClassification nodeClassification = classifyNodes(nodes);
        Map<Integer, Wf_pattern_updater.ServiceParameters> serviceParams =
                extractServiceParameters(nodeClassification.serviceNodes);
        return new WorkflowRepresentation(
                linkMatrix,
                nodeClassification.serviceNodes,
                nodeClassification.forkNodes,
                nodeClassification.joinNodes,
                nodeClassification.routerNodes,
                serviceParams);
    }

    private Matrix buildLinkMatrix() {
        Map<jline.lang.JobClass, Map<jline.lang.JobClass, Matrix>> routingMatrix = network.getLinkedRoutingMatrix();
        List<Node> nodes = network.getNodes();
        List<double[]> links = new ArrayList<double[]>();

        if (routingMatrix != null && !routingMatrix.isEmpty()) {
            for (Map.Entry<jline.lang.JobClass, Map<jline.lang.JobClass, Matrix>> e1 : routingMatrix.entrySet()) {
                for (Map.Entry<jline.lang.JobClass, Matrix> e2 : e1.getValue().entrySet()) {
                    Matrix matrix = e2.getValue();
                    for (int i = 0; i < matrix.getNumRows(); i++) {
                        for (int j = 0; j < matrix.getNumCols(); j++) {
                            double prob = matrix.get(i, j);
                            if (prob > GlobalConstants.Zero) {
                                links.add(new double[]{(double) i, (double) j, prob});
                            }
                        }
                    }
                }
            }
        }

        if (links.isEmpty() && nodes.size() > 1) {
            for (int i = 0; i < nodes.size() - 1; i++) {
                links.add(new double[]{(double) i, (double) (i + 1), 1.0});
            }
        }

        Matrix linkMatrix = Matrix.zeros(Math.max(1, links.size()), 3);
        for (int idx = 0; idx < links.size(); idx++) {
            double[] link = links.get(idx);
            linkMatrix.set(idx, 0, link[0]);
            linkMatrix.set(idx, 1, link[1]);
            linkMatrix.set(idx, 2, link[2]);
        }
        return linkMatrix;
    }

    private NodeClassification classifyNodes(List<Node> nodes) {
        List<Integer> serviceNodes = new ArrayList<Integer>();
        List<Integer> forkNodes = new ArrayList<Integer>();
        List<Integer> joinNodes = new ArrayList<Integer>();
        List<Integer> routerNodes = new ArrayList<Integer>();

        for (int index = 0; index < nodes.size(); index++) {
            Node node = nodes.get(index);
            if (node instanceof Queue) {
                serviceNodes.add(index);
            } else if (node instanceof Delay) {
                serviceNodes.add(index);
            } else if (node instanceof Fork) {
                forkNodes.add(index);
            } else if (node instanceof Join) {
                joinNodes.add(index);
            } else if (node instanceof Router) {
                routerNodes.add(index);
            } else if (node instanceof ClassSwitch) {
                routerNodes.add(index);
            }
        }
        return new NodeClassification(serviceNodes, forkNodes, joinNodes, routerNodes);
    }

    private Map<Integer, Wf_pattern_updater.ServiceParameters> extractServiceParameters(List<Integer> serviceNodes) {
        Map<Integer, Wf_pattern_updater.ServiceParameters> params = new HashMap<Integer, Wf_pattern_updater.ServiceParameters>();
        List<Node> nodes = network.getNodes();
        for (Integer nodeIndex : serviceNodes) {
            if (nodeIndex < nodes.size()) {
                Node node = nodes.get(nodeIndex);
                Wf_pattern_updater.ServiceParameters serviceParam;
                if (node instanceof Queue) {
                    serviceParam = extractPHFromQueue((Queue) node);
                } else if (node instanceof Delay) {
                    serviceParam = extractPHFromDelay((Delay) node);
                } else {
                    serviceParam = createDefaultPH();
                }
                params.put(nodeIndex, serviceParam);
            }
        }
        return params;
    }

    private Wf_pattern_updater.ServiceParameters extractPHFromQueue(Queue queue) {
        return createDefaultPH();
    }

    private Wf_pattern_updater.ServiceParameters extractPHFromDelay(Delay delay) {
        return createDefaultPH();
    }

    private Wf_pattern_updater.ServiceParameters createDefaultPH() {
        Matrix alpha = Matrix.ones(1, 1);
        Matrix T = Matrix.ones(1, 1);
        T.set(0, 0, -1.0);
        return new Wf_pattern_updater.ServiceParameters(alpha, T);
    }

    private DetectedPatterns detectAllPatterns(WorkflowRepresentation workflow) {
        List<List<Integer>> sequences = Wf_sequence_detector.detectSequences(
                workflow.linkMatrix, workflow.serviceNodes);
        List<List<Integer>> parallels = Wf_parallel_detector.detectParallel(
                workflow.linkMatrix, workflow.serviceNodes, workflow.forkNodes, workflow.joinNodes);
        List<Integer> loops = Wf_loop_detector.detectLoops(
                workflow.linkMatrix, workflow.serviceNodes, workflow.routerNodes, workflow.joinNodes);
        List<Wf_branch_detector.BranchPattern> branches = Wf_branch_detector.detectBranches(
                workflow.linkMatrix, workflow.serviceNodes, workflow.joinNodes);
        return new DetectedPatterns(sequences, parallels, loops, branches);
    }

    private Map<String, Object> generateAnalysisStatistics(
            WorkflowRepresentation originalWorkflow,
            DetectedPatterns patterns,
            Wf_pattern_updater.UpdatedWorkflow optimizedWorkflow) {
        Map<String, Object> stats = new HashMap<String, Object>();
        stats.put("sequenceStats", Wf_sequence_detector.getSequenceStats(patterns.sequences));
        stats.put("parallelStats", Wf_parallel_detector.getParallelStats(patterns.parallels));
        stats.put("loopStats", Wf_loop_detector.getLoopStats(
                patterns.loops, originalWorkflow.linkMatrix, originalWorkflow.routerNodes));
        stats.put("branchStats", Wf_branch_detector.getBranchStats(patterns.branches));
        stats.put("updateStats", Wf_pattern_updater.getUpdateStats(
                originalWorkflow.linkMatrix, optimizedWorkflow));
        stats.put("originalComplexity", calculateWorkflowComplexity(originalWorkflow));
        stats.put("optimizedComplexity", calculateOptimizedComplexity(optimizedWorkflow));
        return stats;
    }

    private Map<String, Object> calculateWorkflowComplexity(WorkflowRepresentation workflow) {
        Map<String, Object> complexity = new HashMap<String, Object>();
        complexity.put("totalNodes", workflow.serviceNodes.size() + workflow.forkNodes.size()
                + workflow.joinNodes.size() + workflow.routerNodes.size());
        complexity.put("totalLinks", workflow.linkMatrix.getNumRows());
        complexity.put("serviceNodes", workflow.serviceNodes.size());
        complexity.put("controlNodes", workflow.forkNodes.size() + workflow.joinNodes.size() + workflow.routerNodes.size());

        Set<Integer> nodeSet = new HashSet<Integer>();
        for (int i = 0; i < workflow.linkMatrix.getNumRows(); i++) {
            nodeSet.add((int) workflow.linkMatrix.get(i, 0));
            nodeSet.add((int) workflow.linkMatrix.get(i, 1));
        }
        complexity.put("connectedNodes", nodeSet.size());

        Map<Integer, Integer> degrees = new HashMap<Integer, Integer>();
        for (int i = 0; i < workflow.linkMatrix.getNumRows(); i++) {
            int start = (int) workflow.linkMatrix.get(i, 0);
            int end = (int) workflow.linkMatrix.get(i, 1);
            Integer cur = degrees.get(start);
            degrees.put(start, (cur == null ? 0 : cur) + 1);
            cur = degrees.get(end);
            degrees.put(end, (cur == null ? 0 : cur) + 1);
        }
        double avgDegree = 0.0;
        if (!degrees.isEmpty()) {
            double sum = 0.0;
            for (Integer v : degrees.values()) sum += v;
            avgDegree = sum / degrees.size();
        }
        complexity.put("avgDegree", avgDegree);
        return complexity;
    }

    private Map<String, Object> calculateOptimizedComplexity(Wf_pattern_updater.UpdatedWorkflow optimizedWorkflow) {
        Map<String, Object> complexity = new HashMap<String, Object>();
        complexity.put("totalNodes", optimizedWorkflow.getServiceParameters().size());
        complexity.put("totalLinks", optimizedWorkflow.getLinkMatrix().getNumRows());

        Set<Integer> nodeSet = new HashSet<Integer>();
        for (int i = 0; i < optimizedWorkflow.getLinkMatrix().getNumRows(); i++) {
            nodeSet.add((int) optimizedWorkflow.getLinkMatrix().get(i, 0));
            nodeSet.add((int) optimizedWorkflow.getLinkMatrix().get(i, 1));
        }
        complexity.put("connectedNodes", nodeSet.size());
        return complexity;
    }

    public List<String> getOptimizationRecommendations(WorkflowAnalysis analysis) {
        List<String> recommendations = new ArrayList<String>();
        DetectedPatterns patterns = analysis.detectedPatterns;
        Object updateStatsObj = analysis.statistics.get("updateStats");
        @SuppressWarnings("unchecked")
        Map<String, Object> updateStats = updateStatsObj instanceof Map ? (Map<String, Object>) updateStatsObj : null;

        if (!patterns.sequences.isEmpty()) {
            recommendations.add("Found " + patterns.sequences.size() + " sequence patterns that can be simplified");
        }
        if (!patterns.parallels.isEmpty()) {
            recommendations.add("Found " + patterns.parallels.size() + " parallel patterns for potential optimization");
        }
        if (!patterns.loops.isEmpty()) {
            recommendations.add("Found " + patterns.loops.size() + " loop patterns - consider loop unrolling for performance");
        }
        if (!patterns.branches.isEmpty()) {
            recommendations.add("Found " + patterns.branches.size() + " branch patterns - analyze probability distributions");
            int highEntropy = 0;
            for (Wf_branch_detector.BranchPattern bp : patterns.branches) {
                Map<String, Double> diversity = Wf_branch_detector.calculateBranchDiversity(bp);
                Double e = diversity.get("entropy");
                double entropy = e != null ? e.doubleValue() : 0.0;
                if (entropy > 1.0) highEntropy++;
            }
            if (highEntropy > 0) {
                recommendations.add(highEntropy + " branches have high entropy - consider load balancing");
            }
        }

        if (updateStats != null) {
            Object rrObj = updateStats.get("reductionRatio");
            double reductionRatio = rrObj instanceof Number ? ((Number) rrObj).doubleValue() : 0.0;
            if (reductionRatio > 0.1) {
                recommendations.add("Workflow complexity reduced by " + (int) (reductionRatio * 100)
                        + "% through pattern optimization");
            }
        }
        return recommendations;
    }

    public boolean validateAnalysis(WorkflowAnalysis analysis) {
        if (!Wf_pattern_updater.validateUpdatedWorkflow(analysis.optimizedWorkflow)) {
            return false;
        }
        for (List<Integer> sequence : analysis.detectedPatterns.sequences) {
            if (!Wf_sequence_detector.validateSequence(sequence, analysis.originalWorkflow.linkMatrix)) {
                return false;
            }
        }
        for (List<Integer> parallel : analysis.detectedPatterns.parallels) {
            if (!Wf_parallel_detector.validateParallelPattern(
                    parallel,
                    analysis.originalWorkflow.linkMatrix,
                    analysis.originalWorkflow.forkNodes,
                    analysis.originalWorkflow.joinNodes)) {
                return false;
            }
        }
        for (Integer loop : analysis.detectedPatterns.loops) {
            if (!Wf_loop_detector.validateLoopPattern(
                    loop,
                    analysis.originalWorkflow.linkMatrix,
                    analysis.originalWorkflow.routerNodes)) {
                return false;
            }
        }
        for (Wf_branch_detector.BranchPattern branch : analysis.detectedPatterns.branches) {
            if (!Wf_branch_detector.validateBranchPattern(branch, analysis.originalWorkflow.linkMatrix)) {
                return false;
            }
        }
        return true;
    }
}
