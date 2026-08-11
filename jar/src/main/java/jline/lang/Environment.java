/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

// Copyright (c) 2012-2026, Imperial College London
// All rights reserved.

package jline.lang;

import jline.api.mam.Map_mean;
import jline.api.mam.Mmap_count_lambda;
import jline.api.mam.Mmap_normalize;
import jline.lang.layered.LayeredNetwork;
import jline.lang.nodes.Node;
import jline.lang.nodes.ServiceNode;
import jline.lang.processes.Markovian;
import jline.lang.processes.ContinuousDistribution;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import static jline.api.mc.Ctmc_solve.ctmc_solve;


/**
 * An environment model defined by a collection of network sub-models coupled with an environment transition rule
 * that selects the active sub-model.
 */
public class Environment extends Ensemble {

    public final ContinuousDistribution[][] env;
    private final String[] names;
    private final String[] types;
    private final Model[] models;
    // Markovian representation of each stage transition
    public MatrixCell[][] proc;
    public MatrixCell[] holdTime; // Holding times
    public Matrix probEnv; // Steady-stage probability of the environment
    public Matrix probOrig; // Probability that a request originated from phase
    public ResetQueueLengthsFunction[][] resetQLFun; // Function implementing the reset policy
    private final List<NodeFailure> nodeFailures = new ArrayList<NodeFailure>(); // Breakdown/repair descriptors
    public ResetStateFunction[][] resetStateFun; // State-vector reset policy (statevec analyzer)
    public ResetEnvRatesFunction[][] resetEnvRatesFun;

    @SuppressWarnings("unchecked")
    /**
     * Creates a new environment model with the specified name.
     * Defaults to 10 stages. Each stage can hold a different network sub-model with transitions between them.
     *
     * @param name the name of this environment model
     */
    public Environment(String name) {
        this(name, 10);
    }

    @SuppressWarnings("unchecked")
    /**
     * Creates a new environment model with the specified number of stages.
     * Each stage can hold a different network sub-model with transitions between them.
     *
     * @param name the name of this environment model
     * @param numStages the number of stages in this environment
     */
    public Environment(String name, int numStages) {
        super(new ArrayList<>());
        this.setName(name);
        this.env = new ContinuousDistribution[numStages][numStages];
        this.names = new String[numStages];
        this.types = new String[numStages];
        this.models = new Model[numStages];
        this.proc = new MatrixCell[numStages][numStages];
        this.holdTime = new MatrixCell[numStages];
        this.probEnv = new Matrix(0, 0);
        this.probOrig = new Matrix(0, 0);
        this.resetQLFun = new ResetQueueLengthsFunction[numStages][numStages];
        this.resetEnvRatesFun = new ResetEnvRatesFunction[numStages][numStages];
        // Default state-vector reset is the identity map (statevec analyzer).
        this.resetStateFun = new ResetStateFunction[numStages][numStages];
        for (int i = 0; i < numStages; i++) {
            for (int j = 0; j < numStages; j++) {
                this.resetStateFun[i][j] = (Matrix piExit) -> piExit;
            }
        }
    }

    /**
     * Adds a network model to a specific stage of the environment.
     * All stages must have networks with the same number of stateful nodes.
     * 
     * @param stageIdx the index of the stage (0-based)
     * @param name the name of this stage
     * @param type the type classification for this stage
     * @param model the network model to associate with this stage
     * @throws RuntimeException if the model has a different number of stateful nodes than other stages
     */
    public void addStage(int stageIdx, String name, String type, Model model) {
        this.names[stageIdx] = name;
        this.types[stageIdx] = type;
        this.models[stageIdx] = model;
        if (stageIdx > 0 && statefulNodesOf(model) != statefulNodesOf(models[0])) {
            throw new RuntimeException(
                    "Unsupported feature. Random environment stages must map to networks with identical number of stateful nodes.");
        }
        // The inherited (List<Network>) ensemble holds only flat-network stages,
        // for backward compatibility with getEnsemble(). LayeredNetwork stages
        // are reached through getStageModels() instead.
        if (model instanceof Network) {
            this.ensemble.add(stageIdx, (Network) model);
        }
    }

    /** Stage models by absolute index; may hold Network or LayeredNetwork. */
    public Model[] getStageModels() {
        return models;
    }

    private static int statefulNodesOf(Model m) {
        if (m instanceof Network) return ((Network) m).getNumberOfStatefulNodes();
        if (m instanceof LayeredNetwork) return ((LayeredNetwork) m).getNumberOfStatefulNodes();
        return 0;
    }

    private static int nodesOf(Model m) {
        if (m instanceof Network) return ((Network) m).getNumberOfNodes();
        if (m instanceof LayeredNetwork) return ((LayeredNetwork) m).getNumberOfNodes();
        return 0;
    }

    private static int classesOf(Model m) {
        if (m instanceof Network) return ((Network) m).getNumberOfClasses();
        if (m instanceof LayeredNetwork) return ((LayeredNetwork) m).getNumberOfClasses();
        return 0;
    }

    /**
     * Adds a transition between two stages with default reset function (identity).
     * 
     * @param fromStageIdx the source stage index
     * @param toStageIdx the destination stage index
     * @param distrib the Markovian distribution governing this transition
     */
    public void addTransition(int fromStageIdx, int toStageIdx, Markovian distrib) {
        this.addTransition(fromStageIdx, toStageIdx, distrib, input -> input);
    }

    /**
     * Adds a transition between two stages with a custom reset function.
     *
     * @param fromStageIdx the source stage index
     * @param toStageIdx the destination stage index
     * @param distrib the Markovian distribution governing this transition
     * @param resetFun function to apply when transitioning to reset queue lengths
     */
    public void addTransition(
            int fromStageIdx,
            int toStageIdx,
            Markovian distrib,
            ResetQueueLengthsFunction resetFun) {
        this.env[fromStageIdx][toStageIdx] = distrib;
        this.resetQLFun[fromStageIdx][toStageIdx] = resetFun;
        // Don't set resetEnvRatesFun - leave it null for state-independent transitions
    }

    /**
     * Adds a transition between two stages with custom reset functions for both
     * queue lengths and environment rates.
     *
     * @param fromStageIdx the source stage index
     * @param toStageIdx the destination stage index
     * @param distrib the Markovian distribution governing this transition
     * @param resetFun function to apply when transitioning to reset queue lengths
     * @param resetEnvRatesFun function to apply when transitioning to reset environment rates
     */
    public void addTransition(
            int fromStageIdx,
            int toStageIdx,
            Markovian distrib,
            ResetQueueLengthsFunction resetFun,
            ResetEnvRatesFunction resetEnvRatesFun) {
        this.env[fromStageIdx][toStageIdx] = distrib;
        this.resetQLFun[fromStageIdx][toStageIdx] = resetFun;
        this.resetEnvRatesFun[fromStageIdx][toStageIdx] = resetEnvRatesFun;
    }

    @SuppressWarnings("unchecked")
    /**
     * Initializes the environment by computing stage transition rates, holding times,
     * steady-state probabilities, and stage embedding probabilities.
     * This method must be called after all stages and transitions have been added.
     */
    public void init() {
        int E = this.models.length;
        Matrix Pemb = new Matrix(E, E);

        // Analyse holding times
        MatrixCell[][] emmap = new MatrixCell[E][E];
        for (int e = 0; e < E; e++) {
            for (int h = 0; h < E; h++) {
                // Multiclass MMAP representation
                if (this.env[e][h] == null) {
                    Matrix zero = new Matrix(1, 1);
                    MatrixCell zeroMap = new MatrixCell();
                    zeroMap.set(0, zero);
                    zeroMap.set(1, zero);
                    emmap[e][h] = zeroMap;
                } else {
                    emmap[e][h] = this.env[e][h].getProcess();
                }
                for (int j = 0; j < E; j++) {
                    emmap[e][h].set(j + 2, emmap[e][h].get(1).copy());
                    if (j != h) {
                        emmap[e][h].get(j + 2).zero();
                    }
                }
            }

            holdTime[e] = new MatrixCell();
            for (int i = 0; i < emmap[e][e].size(); i++) {
                holdTime[e].set(i, emmap[e][e].get(i).copy());
            }

            for (int h = 0; h < E; h++) {
                if (h != e) {
                    this.holdTime[e].set(0, this.holdTime[e].get(0).krons(emmap[e][h].get(0)));
                    for (int j = 1; j < E + 2; j++) {
                        this.holdTime[e].set(j, this.holdTime[e].get(j).krons(emmap[e][h].get(j)));
                        Matrix ones = new Matrix(holdTime[e].get(j).length(), 1);
                        ones.ones();
                        Matrix completionRates = holdTime[e].get(j).mult(ones, new Matrix(0, 0));
                        holdTime[e].get(j).zero();
                        for (int row = 0; row < completionRates.getNumRows(); row++) {
                            holdTime[e].get(j).set(row, 0, completionRates.get(row, 0));
                        }
                    }
                    holdTime[e] = Mmap_normalize.mmap_normalize(holdTime[e]);
                }
            }
            // Completion rates for the different transitions
            Matrix countLambda = Mmap_count_lambda.mmap_count_lambda(holdTime[e]);
            double sumCountLambda = countLambda.sumRows(0);
            for (int col = 0; col < Pemb.getNumCols(); col++) {
                Pemb.set(e, col, countLambda.get(0, col) / sumCountLambda);
            }
        }
        this.proc = emmap;

        Matrix lambda = new Matrix(1, E);
        Matrix A = new Matrix(E, E);
        Matrix I = Matrix.eye(E);
        for (int e = 0; e < E; e++) {
            lambda.set(
                    0,
                    e,
                    1
                            / Map_mean.map_mean(holdTime[e].get(0), holdTime[e].get(1)));
            for (int h = 0; h < E; h++) {
                A.set(e, h, -lambda.get(0, e) * (I.get(e, h) - Pemb.get(e, h)));
            }
        }

        int countLambdaValuesLEQZero = 0;
        for (int col = 0; col < lambda.length(); col++) {
            if (lambda.get(0, col) <= 0) {
                countLambdaValuesLEQZero++;
            }
        }
        if (countLambdaValuesLEQZero == 0) {
            this.probEnv = ctmc_solve(A);
            this.probOrig = new Matrix(E, E);
            for (int e = 0; e < E; e++) {
                for (int h = 0; h < E; h++) {
                    probOrig.set(h, e, probEnv.get(0, h) * lambda.get(0, h) * Pemb.get(h, e));
                }
                if (probEnv.get(0, e) > 0) {
                    double probOrigSumCol = probOrig.sumCols(e);
                    for (int row = 0; row < E; row++) {
                        probOrig.set(row, e, probOrig.get(row, e) / probOrigSumCol);
                    }
                }
            }
        }
    }

    /**
     * Prints a formatted table showing all stages, their properties, and transitions.
     * Displays stage names, types, associated networks, and transition rates.
     */
    public void printStageTable() {
        int numStages = names.length;
        System.out.println("Stage Table:");
        System.out.println("============");
        for (int i = 0; i < numStages; i++) {
            System.out.printf("Stage %d: %s (Type: %s)%n", i + 1, names[i], types[i]);
            if (models[i] != null) {
                System.out.printf("  - Network: %s%n", models[i].getName());
                System.out.printf("  - Nodes: %d%n", nodesOf(models[i]));
                System.out.printf("  - Classes: %d%n", classesOf(models[i]));
            }
        }
        System.out.println("\nTransitions:");
        for (int i = 0; i < numStages; i++) {
            for (int j = 0; j < numStages; j++) {
                if (env[i][j] != null) {
                    double rate = env[i][j].getMean();
                    if (rate > 0) {
                        System.out.printf("  %s -> %s: rate = %.4f%n", names[i], names[j], 1.0 / rate);
                    }
                }
            }
        }
    }

    /**
     * Finds the stage index by name.
     *
     * @param stageName the name of the stage to find
     * @return the index of the stage, or -1 if not found
     */
    public int findStageByName(String stageName) {
        for (int i = 0; i < names.length; i++) {
            if (names[i] != null && names[i].equals(stageName)) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Gets the name of the stage at the given index.
     *
     * @param stageIdx the index of the stage (0-based)
     * @return the name of the stage
     * @throws ArrayIndexOutOfBoundsException if stageIdx is out of range
     */
    public String getStageName(int stageIdx) {
        return names[stageIdx];
    }

    /**
     * Gets the number of stages in this environment.
     *
     * @return the number of stages
     */
    public int getNumberOfStages() {
        return names.length;
    }

    /**
     * Adds a breakdown stage for a specific node in the network.
     * This creates an UP stage (if not already present) and a DOWN stage where the specified
     * node has a reduced service rate. A transition from UP to DOWN is added with the given
     * breakdown distribution.
     *
     * @param baseModel the base network model with normal (UP) service rates
     * @param nodeName the name of the node that can break down
     * @param breakdownDist the distribution for time until breakdown (UP->DOWN transition)
     * @param downServiceDist the service distribution when the node is down
     * @throws RuntimeException if the node is not found in the base model
     *
     * @example
     * <pre>
     * Network model = new Network("MyNetwork");
     * Queue queue = new Queue(model, "Server1", SchedStrategy.FCFS);
     * ClosedClass jobClass = new ClosedClass(model, "Jobs", 10, queue, 0);
     * queue.setService(jobClass, new Exp(2.0)); // UP service rate
     *
     * Environment env = new Environment("ServerEnv", 2);
     * env.addNodeBreakdown(model, "Server1", new Exp(0.1), new Exp(0.5));
     * </pre>
     */
    public void addNodeBreakdown(Network baseModel, String nodeName, Markovian breakdownDist, Markovian downServiceDist) {
        addNodeBreakdown(baseModel, nodeName, breakdownDist, downServiceDist, input -> input);
    }

    /**
     * Adds a breakdown stage for a specific node in the network.
     * Overload that accepts a Node object instead of a node name.
     *
     * @param baseModel the base network model with normal (UP) service rates
     * @param node the node that can break down
     * @param breakdownDist the distribution for time until breakdown (UP->DOWN transition)
     * @param downServiceDist the service distribution when the node is down
     */
    public void addNodeBreakdown(Network baseModel, Node node, Markovian breakdownDist, Markovian downServiceDist) {
        addNodeBreakdown(baseModel, node.getName(), breakdownDist, downServiceDist);
    }

    /**
     * Adds a breakdown stage for a specific node with a custom reset policy.
     *
     * @param baseModel the base network model with normal (UP) service rates
     * @param nodeName the name of the node that can break down
     * @param breakdownDist the distribution for time until breakdown
     * @param downServiceDist the service distribution when the node is down
     * @param resetFun function to reset queue lengths on breakdown
     */
    public void addNodeBreakdown(Network baseModel, String nodeName, Markovian breakdownDist,
                                  Markovian downServiceDist, ResetQueueLengthsFunction resetFun) {
        addNodeBreakdownInternal(baseModel, nodeName, breakdownDist, downServiceDist, resetFun);
    }

    /**
     * Adds a breakdown stage for a specific node with a custom reset policy.
     * Overload that accepts a Node object instead of a node name.
     *
     * @param baseModel the base network model with normal (UP) service rates
     * @param node the node that can break down
     * @param breakdownDist the distribution for time until breakdown
     * @param downServiceDist the service distribution when the node is down
     * @param resetFun function to reset queue lengths on breakdown
     */
    public void addNodeBreakdown(Network baseModel, Node node, Markovian breakdownDist,
                                  Markovian downServiceDist, ResetQueueLengthsFunction resetFun) {
        addNodeBreakdownInternal(baseModel, node.getName(), breakdownDist, downServiceDist, resetFun);
    }

    /**
     * Adds a breakdown with a named reset policy.
     *
     * @param baseModel the base network model with normal (UP) service rates
     * @param nodeName the name of the node that can break down
     * @param breakdownDist the distribution for time until breakdown
     * @param downServiceDist the service distribution when the node is down
     * @param resetPolicy named reset policy, "keep" or "clear"
     */
    public void addNodeBreakdown(Network baseModel, String nodeName, Markovian breakdownDist,
                                 Markovian downServiceDist, String resetPolicy) {
        String name = normalizeResetPolicyName(resetPolicy);
        addNodeBreakdownInternal(baseModel, nodeName, breakdownDist, downServiceDist,
                resolveResetPolicy(name), name);
    }

    /**
     * Adds a breakdown with a named reset policy.
     * Overload that accepts a Node object instead of a node name.
     *
     * @param baseModel the base network model with normal (UP) service rates
     * @param node the node that can break down
     * @param breakdownDist the distribution for time until breakdown
     * @param downServiceDist the service distribution when the node is down
     * @param resetPolicy named reset policy, "keep" or "clear"
     */
    public void addNodeBreakdown(Network baseModel, Node node, Markovian breakdownDist,
                                 Markovian downServiceDist, String resetPolicy) {
        addNodeBreakdown(baseModel, node.getName(), breakdownDist, downServiceDist, resetPolicy);
    }

    private void addNodeBreakdownInternal(Network baseModel, String nodeName, Markovian breakdownDist,
                                  Markovian downServiceDist, ResetQueueLengthsFunction resetFun) {
        addNodeBreakdownInternal(baseModel, nodeName, breakdownDist, downServiceDist, resetFun,
                RESET_POLICY_CUSTOM);
    }

    private void addNodeBreakdownInternal(Network baseModel, String nodeName, Markovian breakdownDist,
                                  Markovian downServiceDist, ResetQueueLengthsFunction resetFun,
                                  String resetName) {
        // Create UP stage if this is the first call
        int upStageIdx = findStageByName("UP");
        if (upStageIdx == -1) {
            // Find first empty stage
            for (int i = 0; i < models.length; i++) {
                if (models[i] == null) {
                    Network upModel = baseModel.copy();
                    addStage(i, "UP", "operational", upModel);
                    upStageIdx = i;
                    break;
                }
            }
        }

        // Create DOWN stage with modified service rate for the specified node
        Network downModel = baseModel.copy();

        // Find the node to modify
        int nodeIdx = -1;
        for (int i = 0; i < downModel.getNodes().size(); i++) {
            if (downModel.getNodes().get(i).getName().equals(nodeName)) {
                nodeIdx = i;
                break;
            }
        }

        if (nodeIdx == -1) {
            throw new RuntimeException("Node \"" + nodeName + "\" not found in the base model.");
        }

        // Update service distribution for the down node
        Node node = downModel.getNodes().get(nodeIdx);
        if (node instanceof ServiceNode) {
            ServiceNode serviceNode = (ServiceNode) node;
            for (int c = 0; c < downModel.getNumberOfClasses(); c++) {
                JobClass jobClass = downModel.getClasses().get(c);
                if (serviceNode.getServiceProcess(jobClass) != null) {
                    serviceNode.setService(jobClass, downServiceDist);
                }
            }
        }

        // Add DOWN stage
        String downStageName = "DOWN_" + nodeName;
        int downStageIdx = -1;
        for (int i = 0; i < models.length; i++) {
            if (models[i] == null) {
                addStage(i, downStageName, "failed", downModel);
                downStageIdx = i;
                break;
            }
        }

        // Add breakdown transition (UP -> DOWN)
        addTransition(upStageIdx, downStageIdx, breakdownDist, resetFun);

        // Record the descriptor so that the breakdown can be serialized declaratively.
        recordNodeFailure(nodeName, breakdownDist, downServiceDist, resetName);
    }

    /**
     * Adds a repair transition from DOWN to UP stage for a previously added breakdown.
     *
     * @param nodeName the name of the node that can be repaired
     * @param repairDist the distribution for repair time (DOWN->UP transition)
     * @throws RuntimeException if the DOWN stage for this node is not found
     *
     * @example
     * <pre>
     * env.addNodeRepair("Server1", new Exp(1.0));
     * </pre>
     */
    public void addNodeRepair(String nodeName, Markovian repairDist) {
        addNodeRepairInternal(nodeName, repairDist, keepResetPolicy(), RESET_POLICY_KEEP);
    }

    /**
     * Adds a repair transition with a named reset policy.
     *
     * @param nodeName the name of the node that can be repaired
     * @param repairDist the distribution for repair time
     * @param resetPolicy named reset policy, "keep" or "clear"
     */
    public void addNodeRepair(String nodeName, Markovian repairDist, String resetPolicy) {
        String name = normalizeResetPolicyName(resetPolicy);
        addNodeRepairInternal(nodeName, repairDist, resolveResetPolicy(name), name);
    }

    /**
     * Adds a repair transition with a named reset policy.
     * Overload that accepts a Node object instead of a node name.
     *
     * @param node the node that can be repaired
     * @param repairDist the distribution for repair time
     * @param resetPolicy named reset policy, "keep" or "clear"
     */
    public void addNodeRepair(Node node, Markovian repairDist, String resetPolicy) {
        addNodeRepair(node.getName(), repairDist, resetPolicy);
    }

    /**
     * Adds a repair transition from DOWN to UP stage for a previously added breakdown.
     * Overload that accepts a Node object instead of a node name.
     *
     * @param node the node that can be repaired
     * @param repairDist the distribution for repair time (DOWN->UP transition)
     */
    public void addNodeRepair(Node node, Markovian repairDist) {
        addNodeRepair(node.getName(), repairDist);
    }

    /**
     * Adds a repair transition with a custom reset policy.
     *
     * @param nodeName the name of the node that can be repaired
     * @param repairDist the distribution for repair time
     * @param resetFun function to reset queue lengths on repair
     */
    public void addNodeRepair(String nodeName, Markovian repairDist, ResetQueueLengthsFunction resetFun) {
        addNodeRepairInternal(nodeName, repairDist, resetFun);
    }

    /**
     * Adds a repair transition with a custom reset policy.
     * Overload that accepts a Node object instead of a node name.
     *
     * @param node the node that can be repaired
     * @param repairDist the distribution for repair time
     * @param resetFun function to reset queue lengths on repair
     */
    public void addNodeRepair(Node node, Markovian repairDist, ResetQueueLengthsFunction resetFun) {
        addNodeRepairInternal(node.getName(), repairDist, resetFun);
    }

    private void addNodeRepairInternal(String nodeName, Markovian repairDist, ResetQueueLengthsFunction resetFun) {
        addNodeRepairInternal(nodeName, repairDist, resetFun, RESET_POLICY_CUSTOM);
    }

    private void addNodeRepairInternal(String nodeName, Markovian repairDist, ResetQueueLengthsFunction resetFun,
                                       String resetName) {
        String downStageName = "DOWN_" + nodeName;
        int downStageIdx = findStageByName(downStageName);
        int upStageIdx = findStageByName("UP");

        if (downStageIdx == -1) {
            throw new RuntimeException("DOWN stage for node \"" + nodeName + "\" not found. Call addNodeBreakdown first.");
        }

        if (upStageIdx == -1) {
            throw new RuntimeException("UP stage not found. Call addNodeBreakdown first.");
        }

        // Add repair transition (DOWN -> UP)
        addTransition(downStageIdx, upStageIdx, repairDist, resetFun);

        // Complete the descriptor recorded by addNodeBreakdown.
        int idx = findNodeFailure(nodeName);
        if (idx >= 0) {
            this.nodeFailures.get(idx).repair = repairDist;
            this.nodeFailures.get(idx).repairResetPolicy = resetName;
        }
    }

    /**
     * Convenience method to add both breakdown and repair for a node.
     * This is the most common use case where a node can fail and be repaired.
     *
     * @param baseModel the base network model with normal (UP) service rates
     * @param nodeName the name of the node that can break down and repair
     * @param breakdownDist the distribution for time until breakdown
     * @param repairDist the distribution for repair time
     * @param downServiceDist the service distribution when the node is down
     *
     * @example
     * <pre>
     * Environment env = new Environment("ServerEnv", 2);
     * env.addNodeFailureRepair(model, "Server1", new Exp(0.1), new Exp(1.0), new Exp(0.5));
     * env.init();
     * </pre>
     */
    public void addNodeFailureRepair(Network baseModel, String nodeName,
                                     Markovian breakdownDist, Markovian repairDist,
                                     Markovian downServiceDist) {
        addNodeFailureRepair(baseModel, nodeName, breakdownDist, repairDist, downServiceDist,
                           RESET_POLICY_KEEP, RESET_POLICY_KEEP);
    }

    /**
     * Adds both breakdown and repair for a node with named reset policies.
     *
     * @param baseModel the base network model with normal (UP) service rates
     * @param nodeName the name of the node that can break down and repair
     * @param breakdownDist the distribution for time until breakdown
     * @param repairDist the distribution for repair time
     * @param downServiceDist the service distribution when the node is down
     * @param resetBreakdown named breakdown reset policy, "keep" or "clear"
     * @param resetRepair named repair reset policy, "keep" or "clear"
     */
    public void addNodeFailureRepair(Network baseModel, String nodeName,
                                     Markovian breakdownDist, Markovian repairDist,
                                     Markovian downServiceDist,
                                     String resetBreakdown, String resetRepair) {
        addNodeBreakdown(baseModel, nodeName, breakdownDist, downServiceDist, resetBreakdown);
        addNodeRepair(nodeName, repairDist, resetRepair);
    }

    /**
     * Adds both breakdown and repair for a node with named reset policies.
     * Overload that accepts a Node object instead of a node name.
     *
     * @param baseModel the base network model with normal (UP) service rates
     * @param node the node that can break down and repair
     * @param breakdownDist the distribution for time until breakdown
     * @param repairDist the distribution for repair time
     * @param downServiceDist the service distribution when the node is down
     * @param resetBreakdown named breakdown reset policy, "keep" or "clear"
     * @param resetRepair named repair reset policy, "keep" or "clear"
     */
    public void addNodeFailureRepair(Network baseModel, Node node,
                                     Markovian breakdownDist, Markovian repairDist,
                                     Markovian downServiceDist,
                                     String resetBreakdown, String resetRepair) {
        addNodeFailureRepair(baseModel, node.getName(), breakdownDist, repairDist, downServiceDist,
                           resetBreakdown, resetRepair);
    }

    /**
     * Convenience method to add both breakdown and repair for a node.
     * Overload that accepts a Node object instead of a node name.
     *
     * @param baseModel the base network model with normal (UP) service rates
     * @param node the node that can break down and repair
     * @param breakdownDist the distribution for time until breakdown
     * @param repairDist the distribution for repair time
     * @param downServiceDist the service distribution when the node is down
     */
    public void addNodeFailureRepair(Network baseModel, Node node,
                                     Markovian breakdownDist, Markovian repairDist,
                                     Markovian downServiceDist) {
        addNodeFailureRepair(baseModel, node.getName(), breakdownDist, repairDist, downServiceDist);
    }

    /**
     * Adds both breakdown and repair for a node with custom reset policies.
     *
     * @param baseModel the base network model with normal (UP) service rates
     * @param nodeName the name of the node that can break down and repair
     * @param breakdownDist the distribution for time until breakdown
     * @param repairDist the distribution for repair time
     * @param downServiceDist the service distribution when the node is down
     * @param resetBreakdown reset function for breakdown transition
     * @param resetRepair reset function for repair transition
     */
    public void addNodeFailureRepair(Network baseModel, String nodeName,
                                     Markovian breakdownDist, Markovian repairDist,
                                     Markovian downServiceDist,
                                     ResetQueueLengthsFunction resetBreakdown,
                                     ResetQueueLengthsFunction resetRepair) {
        addNodeBreakdown(baseModel, nodeName, breakdownDist, downServiceDist, resetBreakdown);
        addNodeRepair(nodeName, repairDist, resetRepair);
    }

    /**
     * Adds both breakdown and repair for a node with custom reset policies.
     * Overload that accepts a Node object instead of a node name.
     *
     * @param baseModel the base network model with normal (UP) service rates
     * @param node the node that can break down and repair
     * @param breakdownDist the distribution for time until breakdown
     * @param repairDist the distribution for repair time
     * @param downServiceDist the service distribution when the node is down
     * @param resetBreakdown reset function for breakdown transition
     * @param resetRepair reset function for repair transition
     */
    public void addNodeFailureRepair(Network baseModel, Node node,
                                     Markovian breakdownDist, Markovian repairDist,
                                     Markovian downServiceDist,
                                     ResetQueueLengthsFunction resetBreakdown,
                                     ResetQueueLengthsFunction resetRepair) {
        addNodeFailureRepair(baseModel, node.getName(), breakdownDist, repairDist, downServiceDist,
                           resetBreakdown, resetRepair);
    }

    /**
     * Sets the reset policy for queue lengths when a node breaks down.
     *
     * @param nodeName the name of the node
     * @param resetFun function to reset queue lengths: resetFun(q) -> q_new
     *                 Common policies:
     *                 - input -> input (keep all jobs, default)
     *                 - input -> input.mult(0, null) (clear all queues)
     * @throws RuntimeException if the breakdown transition is not found
     */
    public void setBreakdownResetPolicy(String nodeName, ResetQueueLengthsFunction resetFun) {
        setBreakdownResetPolicyInternal(nodeName, resetFun);
    }

    /**
     * Sets the reset policy for queue lengths when a node breaks down.
     * Overload that accepts a Node object instead of a node name.
     *
     * @param node the node
     * @param resetFun function to reset queue lengths: resetFun(q) -> q_new
     */
    public void setBreakdownResetPolicy(Node node, ResetQueueLengthsFunction resetFun) {
        setBreakdownResetPolicyInternal(node.getName(), resetFun);
    }

    /**
     * Sets a named reset policy for queue lengths when a node breaks down.
     *
     * @param nodeName the name of the node
     * @param resetPolicy named reset policy, "keep" or "clear"
     */
    public void setBreakdownResetPolicy(String nodeName, String resetPolicy) {
        String name = normalizeResetPolicyName(resetPolicy);
        setBreakdownResetPolicyInternal(nodeName, resolveResetPolicy(name), name);
    }

    /**
     * Sets a named reset policy for queue lengths when a node breaks down.
     * Overload that accepts a Node object instead of a node name.
     *
     * @param node the node
     * @param resetPolicy named reset policy, "keep" or "clear"
     */
    public void setBreakdownResetPolicy(Node node, String resetPolicy) {
        setBreakdownResetPolicy(node.getName(), resetPolicy);
    }

    private void setBreakdownResetPolicyInternal(String nodeName, ResetQueueLengthsFunction resetFun) {
        setBreakdownResetPolicyInternal(nodeName, resetFun, RESET_POLICY_CUSTOM);
    }

    private void setBreakdownResetPolicyInternal(String nodeName, ResetQueueLengthsFunction resetFun,
                                                 String resetName) {
        String downStageName = "DOWN_" + nodeName;
        int upStageIdx = findStageByName("UP");
        int downStageIdx = findStageByName(downStageName);

        if (upStageIdx == -1 || downStageIdx == -1) {
            throw new RuntimeException("Breakdown transition for node \"" + nodeName + "\" not found.");
        }

        resetQLFun[upStageIdx][downStageIdx] = resetFun;
        int idx = findNodeFailure(nodeName);
        if (idx >= 0) {
            this.nodeFailures.get(idx).breakdownResetPolicy = resetName;
        }
    }

    /**
     * Sets the reset policy for queue lengths when a node is repaired.
     *
     * @param nodeName the name of the node
     * @param resetFun function to reset queue lengths: resetFun(q) -> q_new
     * @throws RuntimeException if the repair transition is not found
     */
    public void setRepairResetPolicy(String nodeName, ResetQueueLengthsFunction resetFun) {
        setRepairResetPolicyInternal(nodeName, resetFun);
    }

    /**
     * Sets the reset policy for queue lengths when a node is repaired.
     * Overload that accepts a Node object instead of a node name.
     *
     * @param node the node
     * @param resetFun function to reset queue lengths: resetFun(q) -> q_new
     */
    public void setRepairResetPolicy(Node node, ResetQueueLengthsFunction resetFun) {
        setRepairResetPolicyInternal(node.getName(), resetFun);
    }

    /**
     * Sets a named reset policy for queue lengths when a node is repaired.
     *
     * @param nodeName the name of the node
     * @param resetPolicy named reset policy, "keep" or "clear"
     */
    public void setRepairResetPolicy(String nodeName, String resetPolicy) {
        String name = normalizeResetPolicyName(resetPolicy);
        setRepairResetPolicyInternal(nodeName, resolveResetPolicy(name), name);
    }

    /**
     * Sets a named reset policy for queue lengths when a node is repaired.
     * Overload that accepts a Node object instead of a node name.
     *
     * @param node the node
     * @param resetPolicy named reset policy, "keep" or "clear"
     */
    public void setRepairResetPolicy(Node node, String resetPolicy) {
        setRepairResetPolicy(node.getName(), resetPolicy);
    }

    private void setRepairResetPolicyInternal(String nodeName, ResetQueueLengthsFunction resetFun) {
        setRepairResetPolicyInternal(nodeName, resetFun, RESET_POLICY_CUSTOM);
    }

    private void setRepairResetPolicyInternal(String nodeName, ResetQueueLengthsFunction resetFun,
                                              String resetName) {
        String downStageName = "DOWN_" + nodeName;
        int upStageIdx = findStageByName("UP");
        int downStageIdx = findStageByName(downStageName);

        if (upStageIdx == -1 || downStageIdx == -1) {
            throw new RuntimeException("Repair transition for node \"" + nodeName + "\" not found.");
        }

        resetQLFun[downStageIdx][upStageIdx] = resetFun;
        int idx = findNodeFailure(nodeName);
        if (idx >= 0) {
            this.nodeFailures.get(idx).repairResetPolicy = resetName;
        }
    }

    /**
     * Computes system-wide reliability metrics (MTTF, MTTR, MTBF, Availability).
     * This method analyzes the breakdown/repair transitions configured using
     * addNodeBreakdown() and addNodeRepair() to compute reliability metrics.
     *
     * @return a Map containing reliability metrics with keys:
     *         "MTTF" - Mean Time To Failure (time from UP to DOWN)
     *         "MTTR" - Mean Time To Repair (time from DOWN to UP)
     *         "MTBF" - Mean Time Between Failures (MTTF + MTTR)
     *         "Availability" - Steady-state probability of being in UP state
     * @throws RuntimeException if the environment has not been initialized or
     *                          if no breakdown/repair transitions are configured
     *
     * @example
     * <pre>
     * Environment env = new Environment("ServerEnv", 2);
     * env.addNodeFailureRepair(model, "Server", new Exp(0.1), new Exp(1.0), new Exp(0.5));
     * env.init();
     * Map&lt;String, Double&gt; reliabilityMetrics = env.getReliabilityTable();
     * System.out.println("MTTF: " + reliabilityMetrics.get("MTTF"));
     * System.out.println("MTTR: " + reliabilityMetrics.get("MTTR"));
     * System.out.println("MTBF: " + reliabilityMetrics.get("MTBF"));
     * System.out.println("Availability: " + reliabilityMetrics.get("Availability"));
     * </pre>
     */
    public java.util.Map<String, Double> getReliabilityTable() {
        // Step 1: Initialize and validate
        if (probEnv == null || probEnv.isEmpty()) {
            init();
        }

        int E = names.length;
        if (E == 0) {
            throw new RuntimeException("Environment has no stages. Add stages before computing reliability metrics.");
        }

        // Step 2: Identify stage types
        int upIdx = findStageByName("UP");
        if (upIdx == -1) {
            throw new RuntimeException("No UP stage found. Use addNodeBreakdown/addNodeRepair to configure breakdown/repair transitions.");
        }

        // Find all DOWN stages
        java.util.ArrayList<Integer> downIndices = new java.util.ArrayList<Integer>();
        for (int i = 0; i < names.length; i++) {
            if (names[i] != null && names[i].startsWith("DOWN_")) {
                downIndices.add(i);
            }
        }

        if (downIndices.isEmpty()) {
            throw new RuntimeException("No DOWN stages found. Use addNodeBreakdown/addNodeRepair to configure breakdown/repair transitions.");
        }

        // Step 3: Extract breakdown rates (UP -> DOWN_*)
        java.util.ArrayList<Double> breakdownRates = new java.util.ArrayList<Double>();
        for (Integer h : downIndices) {
            if (env[upIdx][h] != null) {
                double lambda_h = 1.0 / env[upIdx][h].getMean();
                breakdownRates.add(lambda_h);
            }
        }

        if (breakdownRates.isEmpty()) {
            throw new RuntimeException("No breakdown transitions found (UP -> DOWN_*).");
        }

        // Total failure rate (competing risks)
        double lambda_total = 0.0;
        for (Double rate : breakdownRates) {
            lambda_total += rate;
        }
        double MTTF = 1.0 / lambda_total;

        // Step 4: Extract repair rates (DOWN_* -> UP)
        java.util.ArrayList<Double> repairRates = new java.util.ArrayList<Double>();
        java.util.ArrayList<Double> downProbs = new java.util.ArrayList<Double>();

        for (Integer e : downIndices) {
            if (env[e][upIdx] != null) {
                double mu_e = 1.0 / env[e][upIdx].getMean();
                repairRates.add(mu_e);
                downProbs.add(probEnv.get(0, e));
            }
        }

        if (repairRates.isEmpty()) {
            throw new RuntimeException("No repair transitions found (DOWN_* -> UP).");
        }

        // Normalize probabilities over DOWN states only
        double totalDownProb = 0.0;
        for (Double prob : downProbs) {
            totalDownProb += prob;
        }

        double MTTR;
        if (totalDownProb > 0) {
            // Weighted average repair time
            MTTR = 0.0;
            for (int i = 0; i < repairRates.size(); i++) {
                double downProbNorm = downProbs.get(i) / totalDownProb;
                MTTR += downProbNorm / repairRates.get(i);
            }
        } else {
            // Fallback: simple average if no steady-state probability
            double sum = 0.0;
            for (Double rate : repairRates) {
                sum += 1.0 / rate;
            }
            MTTR = sum / repairRates.size();
        }

        // Step 5: Compute derived metrics
        double MTBF = MTTF + MTTR;

        // Availability from steady-state probabilities
        double availUp = probEnv.get(0, upIdx);
        double availDown = 0.0;
        for (Integer idx : downIndices) {
            availDown += probEnv.get(0, idx);
        }
        double Availability = availUp / (availUp + availDown);

        // Step 6: Create output map
        java.util.LinkedHashMap<String, Double> result = new java.util.LinkedHashMap<String, Double>();
        result.put("MTTF", MTTF);
        result.put("MTTR", MTTR);
        result.put("MTBF", MTBF);
        result.put("Availability", Availability);

        return result;
    }

    /**
     * Short alias for getReliabilityTable.
     *
     * @return a Map containing reliability metrics (MTTF, MTTR, MTBF, Availability)
     * @see #getReliabilityTable()
     */
    public java.util.Map<String, Double> relT() {
        return getReliabilityTable();
    }

    /**
     * Short alias for getReliabilityTable.
     *
     * @return a Map containing reliability metrics (MTTF, MTTR, MTBF, Availability)
     * @see #getReliabilityTable()
     */
    public java.util.Map<String, Double> getRelT() {
        return getReliabilityTable();
    }

    /**
     * Short alias for getReliabilityTable.
     *
     * @return a Map containing reliability metrics (MTTF, MTTR, MTBF, Availability)
     * @see #getReliabilityTable()
     */
    public java.util.Map<String, Double> relTable() {
        return getReliabilityTable();
    }

    /**
     * Short alias for getReliabilityTable.
     *
     * @return a Map containing reliability metrics (MTTF, MTTR, MTBF, Availability)
     * @see #getReliabilityTable()
     */
    public java.util.Map<String, Double> getRelTable() {
        return getReliabilityTable();
    }

    public interface ResetQueueLengthsFunction {
        Matrix reset(Matrix input);
    }

    /**
     * Descriptor of a node breakdown/repair macro applied through
     * {@code addNodeBreakdown} / {@code addNodeRepair}.
     *
     * The stages and transitions of the environment already carry the structure of a
     * breakdown losslessly. This descriptor additionally records the queue-length reset
     * policies, which are functions and are otherwise unrecoverable, so that the
     * breakdown can be serialized declaratively (the {@code nodeFailures} JSON key).
     */
    public static class NodeFailure implements Serializable {
        private static final long serialVersionUID = 1L;

        /** Name of the node that breaks down. */
        public String node;
        /** Distribution of the time until breakdown (UP -> DOWN). */
        public Markovian breakdown;
        /** Distribution of the repair time (DOWN -> UP), or null when no repair was added. */
        public Markovian repair;
        /** Service distribution of the node while it is down. */
        public Markovian downService;
        /** Named breakdown reset policy: "keep", "clear", or "custom" (not serializable). */
        public String breakdownResetPolicy;
        /** Named repair reset policy: "keep", "clear", "custom", or "" when no repair was added. */
        public String repairResetPolicy;
    }

    /** Named reset policy carrying the queue lengths across a transition unchanged. */
    public static final String RESET_POLICY_KEEP = "keep";
    /** Named reset policy emptying the queues on a transition. */
    public static final String RESET_POLICY_CLEAR = "clear";
    /** Marks a reset policy given as an arbitrary function, which cannot be serialized. */
    public static final String RESET_POLICY_CUSTOM = "custom";

    /**
     * Resolves a named queue-length reset policy into its function.
     *
     * @param policyName "keep" (identity) or "clear" (empty the queues), case-insensitive
     * @return the reset function implementing the policy
     * @throws RuntimeException if the policy name is not recognised
     */
    public static ResetQueueLengthsFunction resolveResetPolicy(String policyName) {
        if (policyName == null) {
            return keepResetPolicy();
        }
        String name = policyName.toLowerCase();
        if (RESET_POLICY_KEEP.equals(name)) {
            return keepResetPolicy();
        }
        if (RESET_POLICY_CLEAR.equals(name)) {
            return clearResetPolicy();
        }
        throw new RuntimeException("Unknown reset policy \"" + policyName
                + "\". Use \"keep\", \"clear\", or a ResetQueueLengthsFunction.");
    }

    /**
     * @return the reset function of the "keep" policy: the queue lengths are unchanged
     */
    public static ResetQueueLengthsFunction keepResetPolicy() {
        return new ResetQueueLengthsFunction() {
            @Override
            public Matrix reset(Matrix input) {
                return input;
            }
        };
    }

    /**
     * @return the reset function of the "clear" policy: all queues are emptied
     */
    public static ResetQueueLengthsFunction clearResetPolicy() {
        return new ResetQueueLengthsFunction() {
            @Override
            public Matrix reset(Matrix input) {
                if (input == null) {
                    return null;
                }
                return new Matrix(input.getNumRows(), input.getNumCols());
            }
        };
    }

    /**
     * Normalises a named reset policy, rejecting unknown names.
     */
    private static String normalizeResetPolicyName(String policyName) {
        if (policyName == null) {
            return RESET_POLICY_KEEP;
        }
        String name = policyName.toLowerCase();
        if (RESET_POLICY_KEEP.equals(name) || RESET_POLICY_CLEAR.equals(name)) {
            return name;
        }
        throw new RuntimeException("Unknown reset policy \"" + policyName
                + "\". Use \"keep\", \"clear\", or a ResetQueueLengthsFunction.");
    }

    /**
     * @return the node breakdown/repair descriptors recorded on this environment
     */
    public List<NodeFailure> getNodeFailures() {
        return this.nodeFailures;
    }

    /**
     * Finds the node-failure descriptor for a node.
     *
     * @param nodeName the node name
     * @return the index of the descriptor, or -1 when absent
     */
    public int findNodeFailure(String nodeName) {
        for (int i = 0; i < this.nodeFailures.size(); i++) {
            if (this.nodeFailures.get(i).node.equals(nodeName)) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Attaches a node breakdown/repair descriptor to stages that already exist.
     *
     * This is the counterpart of addNodeBreakdown/addNodeRepair for the case where the UP
     * and DOWN_&lt;node&gt; stages and their transitions have already been built (for
     * instance by LineModelIO reading the expanded stages/transitions form). It records
     * the descriptor and applies the queue-length reset policies, which the expanded form
     * cannot carry.
     *
     * @param nodeName the node that breaks down
     * @param breakdownDist the time-to-breakdown distribution
     * @param repairDist the repair-time distribution, or null when no repair is defined
     * @param downServiceDist the service distribution while down
     * @param breakdownPolicy named breakdown reset policy ("keep" or "clear")
     * @param repairPolicy named repair reset policy ("keep" or "clear")
     */
    public void registerNodeFailure(String nodeName, Markovian breakdownDist, Markovian repairDist,
                                    Markovian downServiceDist, String breakdownPolicy, String repairPolicy) {
        String downStageName = "DOWN_" + nodeName;
        int upStageIdx = findStageByName("UP");
        int downStageIdx = findStageByName(downStageName);
        if (upStageIdx == -1) {
            throw new RuntimeException("Cannot register a node failure on \"" + nodeName
                    + "\": no UP stage is defined in this environment.");
        }
        if (downStageIdx == -1) {
            throw new RuntimeException("Cannot register a node failure on \"" + nodeName
                    + "\": no \"" + downStageName + "\" stage is defined in this environment.");
        }

        String breakdownName = normalizeResetPolicyName(breakdownPolicy);
        resetQLFun[upStageIdx][downStageIdx] = resolveResetPolicy(breakdownName);
        String repairName = "";
        if (repairDist != null) {
            repairName = normalizeResetPolicyName(repairPolicy);
            resetQLFun[downStageIdx][upStageIdx] = resolveResetPolicy(repairName);
        }

        recordNodeFailure(nodeName, breakdownDist, downServiceDist, breakdownName);
        int idx = findNodeFailure(nodeName);
        this.nodeFailures.get(idx).repair = repairDist;
        this.nodeFailures.get(idx).repairResetPolicy = repairName;
    }

    /**
     * Records (or replaces) the breakdown half of a node-failure descriptor.
     */
    private void recordNodeFailure(String nodeName, Markovian breakdownDist, Markovian downServiceDist,
                                   String breakdownPolicyName) {
        NodeFailure nf = new NodeFailure();
        nf.node = nodeName;
        nf.breakdown = breakdownDist;
        nf.downService = downServiceDist;
        nf.repair = null;
        nf.breakdownResetPolicy = breakdownPolicyName;
        nf.repairResetPolicy = "";
        int idx = findNodeFailure(nodeName);
        if (idx >= 0) {
            this.nodeFailures.set(idx, nf);
        } else {
            this.nodeFailures.add(nf);
        }
    }

    /**
     * State-vector reset policy for the SolverENV state-vector analyzer
     * (options.method='statevec'). Maps the exit state distribution of the
     * origin stage onto the state space of the destination stage. The default
     * is the identity (valid when both stages share the same enumerated state
     * space); supply a custom map when the per-stage state spaces differ.
     */
    public interface ResetStateFunction {
        Matrix reset(Matrix piExit);
    }

    public interface ResetEnvRatesFunction {
        Markovian reset(
                Markovian originalDist,
                Matrix QExit,
                Matrix UExit,
                Matrix TExit);
    }

    /**
     * Returns a table of stage information for this environment.
     * Matches MATLAB Environment.getStageTable() which returns a Table with columns:
     * Stage, Name, Type, Prob, HoldT, Model.
     *
     * @return a formatted string table of stage information
     */
    public String getStageTable() {
        int E = this.names.length;
        if (probEnv == null || probEnv.isEmpty()) {
            init();
        }

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("%-8s %-20s %-15s %-12s %-15s %-20s%n",
                "Stage", "Name", "Type", "Prob", "HoldT", "Model"));
        sb.append(String.format("%-8s %-20s %-15s %-12s %-15s %-20s%n",
                "-----", "----", "----", "----", "-----", "-----"));

        for (int e = 0; e < E; e++) {
            if (names[e] == null) continue;
            String stageName = names[e] != null ? names[e] : "";
            String stageType = types[e] != null ? types[e] : "";
            double prob = (probEnv != null && probEnv.length() > e) ? probEnv.get(0, e) : Double.NaN;
            String holdTimeStr = (holdTime[e] != null) ? String.format("MMAP(%d)", holdTime[e].size()) : "N/A";
            String modelName = (models[e] != null) ? models[e].getName() : "N/A";

            sb.append(String.format("%-8d %-20s %-15s %-12.6f %-15s %-20s%n",
                    e + 1, stageName, stageType, prob, holdTimeStr, modelName));
        }

        return sb.toString();
    }

    /**
     * Short alias for getStageTable.
     *
     * @return a formatted string table of stage information
     * @see #getStageTable()
     */
    public String getStageT() {
        return getStageTable();
    }

    // NOTE: the following LINE methods have not been migrated to JLINE
    // a) getEnv() - env has been made public instead, therefore no need for getter
    // b) setEnv() - env has been made public instead, therefore no need for setter
    // c) setStageName() - appears to be legacy code
    // d) setStageType() - appears to be legacy code
    // e) copyElement - overrides an unimplemented method in "Copyable", and is unused
}
