package jline.lang.sections;

import java.io.Serializable;
import java.util.*;

//import jline.solvers.ssa.events.SynchedEvent;
import jline.lang.distributions.CumulativeDistribution;
import jline.lang.*;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategyType;
import jline.lang.nodes.*;
import jline.util.Pair;

/**
 * Output section of a node
 */
public class OutputSection extends Section implements Serializable {
    protected SchedStrategyType schedPolicy;
    protected List<OutputStrategy> outputStrategies;
//    protected Map<OutputStrategy, SynchedEvent> outputEvents;
    protected boolean isClassSwitch;

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
                totalNonRandProbability.put(jobClassIter, cProb+outputStrategy.getProbability());
            } else if (outputStrategy.getRoutingStrategy() == RoutingStrategy.RAND) {
                int serverCt = 0;
                if (totalProbServers.containsKey(jobClassIter)) {
                    serverCt = totalProbServers.get(jobClassIter);
                }
                totalProbServers.put(jobClassIter, serverCt+1);
            } else if (outputStrategy.getRoutingStrategy() == RoutingStrategy.DISABLED) {
                outputStrategy.setProbability(0);
            }
        }

        for (OutputStrategy outputStrategy : this.outputStrategies) {
            if (outputStrategy.getRoutingStrategy() == RoutingStrategy.RAND) {
                JobClass jobClassIter = outputStrategy.getJobClass();
                double randProb = 1;
                if (totalNonRandProbability.containsKey(jobClassIter)) {
                    randProb = 1-totalNonRandProbability.get(jobClassIter);
                }
                if (totalProbServers.containsKey(jobClassIter)) {
                    randProb /= totalProbServers.get(jobClassIter);
                }
                outputStrategy.setProbability(randProb);
            }
        }
    }

    public OutputSection(String className) {
        super(className);
        this.isClassSwitch = false;
        outputStrategies = new ArrayList<OutputStrategy>();
//        outputEvents = new HashMap<OutputStrategy, SynchedEvent>();
    }

    public void setOutputStrategy(JobClass jobClass, RoutingStrategy routingStrategy) {
        for (OutputStrategy outputStrategy : this.outputStrategies) {
            if ((outputStrategy.getJobClass() == jobClass)) {
                outputStrategy.setRoutingStrategy(routingStrategy);
                return;
            }
        }

        OutputStrategy outputStrategy = new OutputStrategy(jobClass, routingStrategy);
        outputStrategies.add(outputStrategy);
    }

    public void setOutputStrategy(JobClass jobClass, RoutingStrategy routingStrategy, Node destination, double probability) {
        for (OutputStrategy outputStrategy : this.outputStrategies) {
            if ((outputStrategy.getJobClass() == jobClass) && (outputStrategy.getDestination() == null || outputStrategy.getDestination() == destination))  {
                outputStrategy.setRoutingStrategy(routingStrategy);
                outputStrategy.setProbability(probability);
                if(outputStrategy.getDestination() == null){
                    outputStrategy.setDestination(destination);
                }
                this.probabilityUpdate();
                return;
            }
        }

        OutputStrategy outputStrategy = new OutputStrategy(jobClass, routingStrategy, destination, probability);
        outputStrategies.add(outputStrategy);
//        outputEvents.put(outputStrategy, new SynchedEvent(this, destination, jobClass, this.isClassSwitch));
        this.probabilityUpdate();
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

    public void printSummary() {
        System.out.println("Outputs:");
        for (OutputStrategy outputStrategies : outputStrategies) {
            if (outputStrategies.getDestination() != null) {
                System.out.format("-%s (%s)\n", outputStrategies.getDestination().getName(),
                        outputStrategies.getDestination().getClass());
            }
        }
    }

    public final List<OutputStrategy> getOutputStrategies() {
        return this.outputStrategies;
    }

//    public SynchedEvent getOutputEvent(JobClass jobClass, Random random) {
//        CumulativeDistribution<OutputStrategy> outputStrategyCumulativeDistribution = new CumulativeDistribution<OutputStrategy>(random);
//
//        if (this.outputStrategies.size() == 0) {
//            throw new RuntimeException("No output strategies found!");
//        }
//
//        for (OutputStrategy outputStrategy : this.outputStrategies) {
//            if (outputStrategy.getDestination() == null) {
//                continue;
//            } else if (outputStrategy.getJobClass() != jobClass) {
//                continue;
//            }
//
//            outputStrategyCumulativeDistribution.addElement(outputStrategy, outputStrategy.getProbability());
//        }
//
//        return this.outputEvents.get(outputStrategyCumulativeDistribution.sample(random));
//    }

    @SuppressWarnings("unchecked")
//    public ArrayList<Pair<SynchedEvent,Double>>  getOutputEvents(JobClass jobClass, Random random) {
//        CumulativeDistribution<OutputStrategy> outputStrategyCumulativeDistribution = new CumulativeDistribution<OutputStrategy>(random);
//
//        if (this.outputStrategies.size() == 0) {
//            throw new RuntimeException("No output strategies found!");
//        }
//
//        for (OutputStrategy outputStrategy : this.outputStrategies) {
//            if (outputStrategy.getDestination() == null) {
//                continue;
//            } else if (outputStrategy.getJobClass() != jobClass) {
//                continue;
//            }
//
//            outputStrategyCumulativeDistribution.addElement(outputStrategy, outputStrategy.getProbability());
//        }
//        ArrayList<Pair<Double,OutputStrategy>> outputStrategies = outputStrategyCumulativeDistribution.getPossibleEventProbability();
//        ArrayList<Pair<SynchedEvent,Double>> outputEvents = new ArrayList<>();
//        for(Pair<Double,OutputStrategy> pair : outputStrategies){
//            outputEvents.add(new Pair(this.outputEvents.get(pair.getRight()),pair.getLeft()));
//        }
//        return outputEvents;
//    }

//    @SuppressWarnings("unchecked")
//    public SynchedEvent getOutputEvent(JobClass jobClass) {
//        return this.getOutputEvent(jobClass, new Random());
//    }

    public List<OutputStrategy> getOutputStrategyByClass(JobClass jobClass) {
    	List<OutputStrategy> res = new ArrayList<OutputStrategy>();
    	for(OutputStrategy outputStrategy : outputStrategies) {
    		if (outputStrategy.getJobClass().equals(jobClass)) {
    			res.add(outputStrategy);
    		}
    	}
    	return res;
    }
}
