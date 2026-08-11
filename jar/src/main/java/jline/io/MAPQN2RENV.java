/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

package jline.io;

import jline.lang.*;
import jline.lang.nodes.*;
import jline.lang.processes.*;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Transforms a queueing network with MMPP2 service into a random environment model.
 *
 * The transformation converts MMPP2 (2-phase Markov Modulated Poisson Process) service
 * distributions into an equivalent random environment model with exponential services.
 * The environment has two stages (one per MMPP phase) with transitions defined by the
 * MMPP D0 matrix.
 */
public class MAPQN2RENV {

    /**
     * Transform a queueing network with MMPP2 service into a random environment model.
     *
     * @param model Network with MMPP2 service distributions
     * @return Environment model with exponential services modulated by MMPP phases
     * @throws RuntimeException if no MMPP2 service distribution is found
     */
    public static Environment mapqn2renv(Network model) {
        // Phase 1: Validate and extract MMPP parameters
        MMPP2Params params = validateAndExtractMMPP(model);

        if (params == null) {
            throw new RuntimeException("Network must contain at least one MMPP2 service distribution");
        }

        // Phase 2: Extract MMPP Parameters
        Matrix D0 = params.D0;
        Matrix D1 = params.D1;

        // Extract service rates from D1 diagonal
        double lambda0 = D1.get(0, 0);
        double lambda1 = D1.get(1, 1);

        // Extract transition rates from D0 off-diagonal
        double sigma01 = D0.get(0, 1);
        double sigma10 = D0.get(1, 0);

        // Validate rates are non-negative
        if (lambda0 < 0 || lambda1 < 0 || sigma01 < 0 || sigma10 < 0) {
            throw new RuntimeException("All extracted rates must be non-negative");
        }

        // Phase 3: Create Environment Model
        Environment envModel = new Environment("MAPQN_Env", 2);

        // Phase 4: Build Stage Networks
        // Create stage network for Phase 0
        Network stageNet0 = buildStageNetwork(model, "Phase0", lambda0);
        envModel.addStage(0, "Phase0", "item", stageNet0);

        // Create stage network for Phase 1
        Network stageNet1 = buildStageNetwork(model, "Phase1", lambda1);
        envModel.addStage(1, "Phase1", "item", stageNet1);

        // Phase 5: Add Environment Transitions
        if (sigma01 > 0) {
            envModel.addTransition(0, 1, new Exp(sigma01));
        }

        if (sigma10 > 0) {
            envModel.addTransition(1, 0, new Exp(sigma10));
        }

        return envModel;
    }

    /**
     * Holds MMPP2 D0 and D1 matrices.
     */
    private static class MMPP2Params {
        Matrix D0;
        Matrix D1;

        MMPP2Params(Matrix D0, Matrix D1) {
            this.D0 = D0;
            this.D1 = D1;
        }
    }

    /**
     * Validate the network has MMPP2 service and extract parameters.
     */
    private static MMPP2Params validateAndExtractMMPP(Network model) {
        MMPP2 firstMMPP2 = null;

        // Iterate through all nodes to find MMPP2 distributions
        for (Node node : model.getNodes()) {
            // Check if node is a Queue or Delay
            if (!(node instanceof Queue) && !(node instanceof Delay)) {
                continue;
            }

            if (node instanceof Queue) {
                Queue queue = (Queue) node;
                // Check service distributions for all classes
                for (JobClass jobClass : model.getClasses()) {
                    Distribution dist = queue.getService(jobClass);
                    if (dist instanceof MMPP2) {
                        if (firstMMPP2 == null) {
                            firstMMPP2 = (MMPP2) dist;
                        }
                    }
                }
            }
        }

        if (firstMMPP2 != null) {
            Matrix D0 = firstMMPP2.D(0);
            Matrix D1 = firstMMPP2.D(1);
            return new MMPP2Params(D0, D1);
        }

        return null;
    }

    /**
     * Build a stage network by cloning the original and replacing MMPP2 with Exp.
     */
    private static Network buildStageNetwork(Network originalModel, String stageName, double expRate) {
        // Clone the network structure and replace MMPP2 with Exp
        Network stageNet = new Network(originalModel.getName() + "_" + stageName);

        // Build mapping from original nodes to cloned nodes
        Map<String, Node> nodeMap = new HashMap<String, Node>();
        List<JobClass> origClasses = originalModel.getClasses();

        // PASS 1: Create all nodes first
        for (Node origNode : originalModel.getNodes()) {
            Node newNode = null;
            if (origNode instanceof Source) {
                newNode = new Source(stageNet, origNode.getName());
            } else if (origNode instanceof Sink) {
                newNode = new Sink(stageNet, origNode.getName());
            } else if (origNode instanceof Delay) {
                newNode = new Delay(stageNet, origNode.getName());
            } else if (origNode instanceof Queue) {
                Queue origQueue = (Queue) origNode;
                newNode = new Queue(stageNet, origNode.getName(), origQueue.getSchedStrategy());
                if (origQueue.getNumberOfServers() > 1) {
                    ((Queue) newNode).setNumberOfServers(origQueue.getNumberOfServers());
                }
            }
            if (newNode != null) {
                nodeMap.put(origNode.getName(), newNode);
            }
        }

        // PASS 2: Create all job classes
        for (JobClass origClass : origClasses) {
            if (origClass instanceof OpenClass) {
                new OpenClass(stageNet, origClass.getName());
            } else if (origClass instanceof ClosedClass) {
                ClosedClass closedClass = (ClosedClass) origClass;
                Station refStat = closedClass.getReferenceStation();
                Node newRefStat = nodeMap.get(refStat.getName());
                if (newRefStat instanceof Station) {
                    new ClosedClass(stageNet, origClass.getName(),
                                   (int) closedClass.getPopulation(), (Station) newRefStat, 0);
                }
            }
        }

        // PASS 3: Set arrivals from Source nodes
        for (Node origNode : originalModel.getNodes()) {
            if (origNode instanceof Source && nodeMap.containsKey(origNode.getName())) {
                Source origSource = (Source) origNode;
                Source newSource = (Source) nodeMap.get(origNode.getName());

                for (int r = 0; r < origClasses.size(); r++) {
                    JobClass origJobClass = origClasses.get(r);
                    Distribution arrDist = origSource.getArrivalDistribution(origJobClass);
                    if (arrDist != null && r < stageNet.getClasses().size()) {
                        JobClass newJobClass = stageNet.getClasses().get(r);
                        newSource.setArrival(newJobClass, arrDist);
                    }
                }
            }
        }

        // PASS 4: Set services on Queue/Delay nodes, replacing MMPP2 with Exp
        for (Node origNode : originalModel.getNodes()) {
            if ((origNode instanceof Queue || origNode instanceof Delay)
                && nodeMap.containsKey(origNode.getName())) {

                Node newNode = nodeMap.get(origNode.getName());

                for (int r = 0; r < origClasses.size(); r++) {
                    JobClass origJobClass = origClasses.get(r);
                    Distribution origDist = null;

                    if (origNode instanceof Queue) {
                        origDist = ((Queue) origNode).getService(origJobClass);
                    } else if (origNode instanceof Delay) {
                        origDist = ((Delay) origNode).getService(origJobClass);
                    }

                    if (origDist != null && r < stageNet.getClasses().size()) {
                        JobClass newJobClass = stageNet.getClasses().get(r);
                        Distribution newDist;

                        if (origDist instanceof MMPP2) {
                            // Replace with exponential
                            newDist = new Exp(expRate);
                        } else {
                            // Keep other distributions as-is
                            newDist = origDist;
                        }

                        if (newNode instanceof Queue) {
                            ((Queue) newNode).setService(newJobClass, newDist);
                        } else if (newNode instanceof Delay) {
                            ((Delay) newNode).setService(newJobClass, newDist);
                        }
                    }
                }
            }
        }

        // PASS 5: Setup routing
        boolean isClosedNetwork = false;
        for (JobClass c : origClasses) {
            if (c instanceof ClosedClass) {
                isClosedNetwork = true;
                break;
            }
        }

        // Collect route nodes in order
        List<Node> routeNodes = new ArrayList<Node>();
        for (Node origNode : originalModel.getNodes()) {
            if (nodeMap.containsKey(origNode.getName())) {
                routeNodes.add(nodeMap.get(origNode.getName()));
            }
        }

        // Setup routing based on network type
        if (routeNodes.size() >= 2) {
            try {
                stageNet.link(Network.serialRouting(routeNodes.toArray(new Node[0])));
            } catch (Exception e) {
                // If serial routing fails, skip
            }
        }

        return stageNet;
    }
}
