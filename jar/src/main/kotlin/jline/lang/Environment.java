/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

// Copyright (c) 2012-2026, Imperial College London
// All rights reserved.

package jline.lang;

import jline.api.mam.Map_meanKt;
import jline.api.mam.Mmap_count_lambdaKt;
import jline.api.mam.Mmap_normalizeKt;
import jline.lang.nodes.Node;
import jline.lang.nodes.ServiceNode;
import jline.lang.processes.Markovian;
import jline.lang.processes.ContinuousDistribution;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.ArrayList;

import static jline.api.mc.Ctmc_solveKt.ctmc_solve;


/**
 * An environment model defined by a collection of network sub-models coupled with an environment transition rule
 * that selects the active sub-model.
 */
public class Environment extends Ensemble {

    public final ContinuousDistribution[][] env;
    private final String[] names;
    private final String[] types;
    private final Network[] models;
    // Markovian representation of each stage transition
    public MatrixCell[][] proc;
    public MatrixCell[] holdTime; // Holding times
    public Matrix probEnv; // Steady-stage probability of the environment
    public Matrix probOrig; // Probability that a request originated from phase
    public ResetQueueLengthsFunction[][] resetQLFun; // Function implementing the reset policy
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
        this.models = new Network[numStages];
        this.proc = new MatrixCell[numStages][numStages];
        this.holdTime = new MatrixCell[numStages];
        this.probEnv = new Matrix(0, 0);
        this.probOrig = new Matrix(0, 0);
        this.resetQLFun = new ResetQueueLengthsFunction[numStages][numStages];
        this.resetEnvRatesFun = new ResetEnvRatesFunction[numStages][numStages];
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
    public void addStage(int stageIdx, String name, String type, Network model) {
        this.names[stageIdx] = name;
        this.types[stageIdx] = type;
        this.models[stageIdx] = model;
        if (stageIdx > 0 && model.getNumberOfStatefulNodes() != models[0].getNumberOfStatefulNodes()) {
            throw new RuntimeException(
                    "Unsupported feature. Random environment stages must map to networks with identical number of stateful nodes.");
        }
        this.ensemble.add(stageIdx, model);
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
                    holdTime[e] = Mmap_normalizeKt.mmap_normalize(holdTime[e]);
                }
            }
            // Completion rates for the different transitions
            Matrix countLambda = Mmap_count_lambdaKt.mmap_count_lambda(holdTime[e]);
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
                            / Map_meanKt.map_mean(holdTime[e].get(0), holdTime[e].get(1)));
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
                System.out.printf("  - Nodes: %d%n", models[i].getNumberOfNodes());
                System.out.printf("  - Classes: %d%n", models[i].getNumberOfClasses());
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

    private void addNodeBreakdownInternal(Network baseModel, String nodeName, Markovian breakdownDist,
                                  Markovian downServiceDist, ResetQueueLengthsFunction resetFun) {
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
        addNodeRepair(nodeName, repairDist, input -> input);
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
                           input -> input, input -> input);
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

    private void setBreakdownResetPolicyInternal(String nodeName, ResetQueueLengthsFunction resetFun) {
        String downStageName = "DOWN_" + nodeName;
        int upStageIdx = findStageByName("UP");
        int downStageIdx = findStageByName(downStageName);

        if (upStageIdx == -1 || downStageIdx == -1) {
            throw new RuntimeException("Breakdown transition for node \"" + nodeName + "\" not found.");
        }

        resetQLFun[upStageIdx][downStageIdx] = resetFun;
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

    private void setRepairResetPolicyInternal(String nodeName, ResetQueueLengthsFunction resetFun) {
        String downStageName = "DOWN_" + nodeName;
        int upStageIdx = findStageByName("UP");
        int downStageIdx = findStageByName(downStageName);

        if (upStageIdx == -1 || downStageIdx == -1) {
            throw new RuntimeException("Repair transition for node \"" + nodeName + "\" not found.");
        }

        resetQLFun[downStageIdx][upStageIdx] = resetFun;
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

    public interface ResetEnvRatesFunction {
        Markovian reset(
                Markovian originalDist,
                Matrix QExit,
                Matrix UExit,
                Matrix TExit);
    }

    // NOTE: the following LINE methods have not been migrated to JLINE
    // a) getEnv() - env has been made public instead, therefore no need for getter
    // b) setEnv() - env has been made public instead, therefore no need for setter
    // c) setStageName() - appears to be legacy code
    // d) setStageType() - appears to be legacy code
    // e) copyElement - overrides an unimplemented method in "Copyable", and is unused
}
