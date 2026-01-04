/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import jline.lang.JobClass;
import jline.lang.OutputStrategy;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategyType;
import jline.lang.nodes.Node;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Output section of a node
 */
public class OutputSection extends Section implements Serializable {
    protected SchedStrategyType schedPolicy;
    protected List<OutputStrategy> outputStrategies;

    public OutputSection(String className) {
        super(className);
        outputStrategies = new ArrayList<OutputStrategy>();
    }

    public final List<OutputStrategy> getOutputStrategies() {
        return this.outputStrategies;
    }

    public List<OutputStrategy> getOutputStrategyByClass(JobClass jobClass) {
        List<OutputStrategy> res = new ArrayList<OutputStrategy>();
        for (OutputStrategy outputStrategy : outputStrategies) {
            // Use index comparison to handle Signal resolution (Signal -> OpenSignal/ClosedSignal)
            if (outputStrategy.getJobClass().getIndex() == jobClass.getIndex()) {
                res.add(outputStrategy);
            }
        }
        return res;
    }

    public void printSummary() {
        System.out.println("Outputs:");
        for (OutputStrategy outputStrategies : outputStrategies) {
            if (outputStrategies.getDestination() != null) {
                System.out.format("-%s (%s)\n", outputStrategies.getDestination().getName(),
                        outputStrategies.getDestination().getClass());
            }
        }
    }

    protected void probabilityUpdate() {
        Map<JobClass, Double> totalNonRandProbability = new HashMap<JobClass, Double>();
        Map<JobClass, Integer> totalProbServers = new HashMap<JobClass, Integer>();
        for (OutputStrategy outputStrategy : this.outputStrategies) {
            JobClass jobClassIter = outputStrategy.getJobClass();
            if (outputStrategy.getRoutingStrategy() == RoutingStrategy.PROB) {
                double cProb = 0;
                if (totalNonRandProbability.containsKey(jobClassIter)) {
                    cProb = totalNonRandProbability.get(jobClassIter);
                }
                totalNonRandProbability.put(jobClassIter, cProb + outputStrategy.getProbability());
            } else if (outputStrategy.getRoutingStrategy() == RoutingStrategy.RAND) {
                int serverCt = 0;
                if (totalProbServers.containsKey(jobClassIter)) {
                    serverCt = totalProbServers.get(jobClassIter);
                }
                totalProbServers.put(jobClassIter, serverCt + 1);
            } else if (outputStrategy.getRoutingStrategy() == RoutingStrategy.DISABLED) {
                outputStrategy.setProbability(0);
            }
        }

        for (OutputStrategy outputStrategy : this.outputStrategies) {
            if (outputStrategy.getRoutingStrategy() == RoutingStrategy.RAND) {
                JobClass jobClassIter = outputStrategy.getJobClass();
                double randProb = 1;
                if (totalNonRandProbability.containsKey(jobClassIter)) {
                    randProb = 1 - totalNonRandProbability.get(jobClassIter);
                }
                if (totalProbServers.containsKey(jobClassIter)) {
                    randProb /= totalProbServers.get(jobClassIter);
                }
                outputStrategy.setProbability(randProb);
            }
        }
    }

    public void resetRouting() {
        List<OutputStrategy> newOutputStrategies = new ArrayList<OutputStrategy>();
        for (OutputStrategy outputStrategy : this.outputStrategies) {
            if (outputStrategy.getDestination() == null) {
                newOutputStrategies.add(outputStrategy);
            }
        }
        this.outputStrategies = newOutputStrategies;
    }

    public void setOutputStrategy(JobClass jobClass, RoutingStrategy routingStrategy) {
        // Remove all existing entries for this class first to ensure consistency
        // This is needed because RAND/DISABLED routing should have only one entry (with null dest)
        // while PROB routing may have multiple entries (with different dests)
        java.util.Iterator<OutputStrategy> iter = outputStrategies.iterator();
        while (iter.hasNext()) {
            OutputStrategy os = iter.next();
            if (os.getJobClass().getIndex() == jobClass.getIndex()) {
                iter.remove();
            }
        }

        // Add a single new entry for this class
        OutputStrategy outputStrategy = new OutputStrategy(jobClass, routingStrategy);
        outputStrategies.add(outputStrategy);
    }

    public void setOutputStrategy(JobClass jobClass, RoutingStrategy routingStrategy, Node destination, double probability) {
        for (OutputStrategy outputStrategy : this.outputStrategies) {
            // Use index comparison for JobClass to handle dynamically created classes
            if ((outputStrategy.getJobClass().getIndex() == jobClass.getIndex()) &&
                (outputStrategy.getDestination() == null || outputStrategy.getDestination() == destination)) {
                outputStrategy.setRoutingStrategy(routingStrategy);
                outputStrategy.setProbability(probability);
                if (outputStrategy.getDestination() == null) {
                    outputStrategy.setDestination(destination);
                }
                this.probabilityUpdate();
                return;
            }
        }

        OutputStrategy outputStrategy = new OutputStrategy(jobClass, routingStrategy, destination, probability);
        outputStrategies.add(outputStrategy);
        this.probabilityUpdate();
    }
}
