/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

// Copyright (c) 2012-2026, Imperial College London
// All rights reserved.

package jline.lang;

import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.lang.processes.Disabled;
import jline.lang.processes.Distribution;
import jline.lang.sections.Forker;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static jline.io.InputOutputKt.line_warning;
import static jline.io.InputOutputKt.mfilename;

/**
 * A model defined by a collection of sub-models
 */
public class Ensemble extends Model {

    protected List<Network> ensemble;

    /**
     * Creates a new ensemble containing the specified network models.
     * 
     * @param models the list of network models to include in this ensemble
     */
    public Ensemble(List<Network> models) {
        super("Ensemble");
        this.ensemble = models;
    }

    /**
     * Creates a new ensemble with the specified name.
     * The ensemble's model list must be set separately.
     * 
     * @param name the name for this ensemble
     */
    public Ensemble(String name) {
        super(name);
    }

    /**
     * Gets the list of network models in this ensemble.
     * 
     * @return the list of network models
     */
    public List<Network> getEnsemble() {
        return this.ensemble;
    }

    /**
     * Sets the list of network models for this ensemble.
     * 
     * @param ensemble the list of network models to set
     */
    public void setEnsemble(List<Network> ensemble) {
        this.ensemble = ensemble;
    }

    /**
     * Gets a specific network model from the ensemble by index.
     * 
     * @param modelIdx the index of the model to retrieve (0-based)
     * @return the network model at the specified index
     * @throws IndexOutOfBoundsException if the index is out of range
     */
    public Network getModel(int modelIdx) {
        return this.ensemble.get(modelIdx);
    }

    /**
     * Gets the number of network models in this ensemble.
     * 
     * @return the number of models in the ensemble
     */
    public int size() {
        return this.ensemble.size();
    }

    // NOTE: the following LINE methods have not been migrated to JLINE
    // a) copyElement - overrides an unimplemented method in "Copyable", and is unused

    /**
     * Creates a union Network from all Networks in an ensemble.
     * Returns a single Network containing all nodes and classes from each Network
     * in the ensemble as disconnected subnetworks. Node and class names are prefixed
     * with their originating model name to avoid collisions.
     *
     * @param ensemble the Ensemble object to merge
     * @return Network object containing merged subnetworks
     */
    public static Network merge(Ensemble ensemble) {
        return mergeNetworks(ensemble.getEnsemble());
    }

    /**
     * Creates a union Network from a list of Networks.
     * Returns a single Network containing all nodes and classes from each Network
     * as disconnected subnetworks. Node and class names are prefixed with their
     * originating model name to avoid collisions.
     *
     * @param models the list of Network objects to merge
     * @return Network object containing merged subnetworks
     */
    public static Network merge(List<Network> models) {
        return mergeNetworks(models);
    }

    /**
     * Internal implementation of merge for a list of Networks.
     */
    private static Network mergeNetworks(List<Network> models) {
        // Edge case: empty ensemble
        if (models == null || models.isEmpty()) {
            line_warning(mfilename(new Object() {}), "Empty ensemble provided, returning empty Network");
            return new Network("EmptyMergedNetwork");
        }

        // Create union network with combined name
        StringBuilder nameBuilder = new StringBuilder();
        for (int i = 0; i < models.size(); i++) {
            if (i > 0) {
                nameBuilder.append("_");
            }
            nameBuilder.append(models.get(i).getName());
        }
        String unionName = nameBuilder.toString();
        if (unionName.length() > 50) {
            unionName = "MergedNetwork_" + models.size();
        }
        Network unionNetwork = new Network(unionName);

        // Check if any model has open classes (needs Source/Sink)
        boolean hasOpenClasses = false;
        for (Network model : models) {
            if (model.hasOpenClasses()) {
                hasOpenClasses = true;
                break;
            }
        }

        // Create single Source and Sink if needed
        Source unionSource = null;
        Sink unionSink = null;
        if (hasOpenClasses) {
            unionSource = new Source(unionNetwork, "MergedSource");
            unionSink = new Sink(unionNetwork, "MergedSink");
        }

        // Storage for mapping old nodes/classes to new ones
        List<Map<String, Node>> nodeMap = new ArrayList<Map<String, Node>>();
        List<Map<String, JobClass>> classMap = new ArrayList<Map<String, JobClass>>();
        List<JoinPendingItem> joinPendingList = new ArrayList<JoinPendingItem>();

        // Initialize maps for each model
        for (int m = 0; m < models.size(); m++) {
            nodeMap.add(new HashMap<String, Node>());
            classMap.add(new HashMap<String, JobClass>());
        }

        // Phase 1: Create nodes for each model
        for (int m = 0; m < models.size(); m++) {
            Network model = models.get(m);
            String modelName = model.getName();
            String prefix = modelName + "_";

            List<Node> nodes = model.getNodes();
            Map<String, Node> currentNodeMap = nodeMap.get(m);

            for (Node oldNode : nodes) {
                String newName = prefix + oldNode.getName();
                Node newNode = null;

                if (oldNode instanceof Source) {
                    // Map to merged source
                    currentNodeMap.put(oldNode.getName(), unionSource);
                    continue;
                } else if (oldNode instanceof Sink) {
                    // Map to merged sink
                    currentNodeMap.put(oldNode.getName(), unionSink);
                    continue;
                } else if (oldNode instanceof Delay) {
                    newNode = new Delay(unionNetwork, newName);
                } else if (oldNode instanceof Queue) {
                    Queue oldQueue = (Queue) oldNode;
                    newNode = new Queue(unionNetwork, newName, oldQueue.getSchedStrategy());
                    Queue newQueue = (Queue) newNode;
                    if (oldQueue.getNumberOfServers() != Integer.MAX_VALUE && oldQueue.getNumberOfServers() > 1) {
                        newQueue.setNumberOfServers(oldQueue.getNumberOfServers());
                    }
                    Matrix lldScaling = oldQueue.getLimitedLoadDependence();
                    if (lldScaling != null && !lldScaling.isEmpty()) {
                        newQueue.setLoadDependence(lldScaling);
                    }
                    if (oldQueue.getCap() < Integer.MAX_VALUE) {
                        newQueue.setCapacity((int) oldQueue.getCap());
                    }
                } else if (oldNode instanceof Router) {
                    newNode = new Router(unionNetwork, newName);
                } else if (oldNode instanceof ClassSwitch) {
                    ClassSwitch oldCS = (ClassSwitch) oldNode;
                    // Skip auto-added ClassSwitch nodes (recreated by link())
                    if (oldCS.autoAdded) {
                        continue;
                    }
                    newNode = new ClassSwitch(unionNetwork, newName);
                } else if (oldNode instanceof Fork) {
                    Fork oldFork = (Fork) oldNode;
                    newNode = new Fork(unionNetwork, newName);
                    Fork newFork = (Fork) newNode;
                    if (oldFork.getOutput() instanceof Forker) {
                        Forker forker = (Forker) oldFork.getOutput();
                        if (forker.tasksPerLink > 1) {
                            newFork.setTasksPerLink((int) forker.tasksPerLink);
                        }
                    }
                } else if (oldNode instanceof Join) {
                    // Queue for later processing after Forks exist
                    joinPendingList.add(new JoinPendingItem(m, (Join) oldNode, newName));
                    continue;
                } else if (oldNode instanceof Logger) {
                    Logger oldLogger = (Logger) oldNode;
                    newNode = new Logger(unionNetwork, newName, oldLogger.getFilePath() + oldLogger.getFileName());
                } else if (oldNode instanceof Cache) {
                    Cache oldCache = (Cache) oldNode;
                    newNode = new Cache(unionNetwork, newName,
                            oldCache.getNumberOfItems(),
                            oldCache.getItemLevelCap(),
                            oldCache.getReplacementStrategy());
                } else if (oldNode instanceof Place) {
                    newNode = new Place(unionNetwork, newName);
                } else if (oldNode instanceof Transition) {
                    newNode = new Transition(unionNetwork, newName);
                } else {
                    line_warning(mfilename(new Object() {}),
                            "Unsupported node type: " + oldNode.getClass().getSimpleName() + ", skipping");
                    continue;
                }

                currentNodeMap.put(oldNode.getName(), newNode);
            }
        }

        // Phase 2: Process pending Join nodes (after Forks exist)
        for (JoinPendingItem item : joinPendingList) {
            int m = item.modelIdx;
            Join oldJoin = item.oldNode;
            String newName = item.newName;
            Map<String, Node> currentNodeMap = nodeMap.get(m);

            Node newNode;
            if (oldJoin.joinOf != null) {
                String forkName = oldJoin.joinOf.getName();
                if (currentNodeMap.containsKey(forkName)) {
                    Node forkNode = currentNodeMap.get(forkName);
                    newNode = new Join(unionNetwork, newName, forkNode);
                } else {
                    newNode = new Join(unionNetwork, newName);
                }
            } else {
                newNode = new Join(unionNetwork, newName);
            }
            currentNodeMap.put(oldJoin.getName(), newNode);
        }

        // Phase 3: Create job classes for each model
        for (int m = 0; m < models.size(); m++) {
            Network model = models.get(m);
            String modelName = model.getName();
            String prefix = modelName + "_";

            List<JobClass> classes = model.getClasses();
            Map<String, JobClass> currentClassMap = classMap.get(m);
            Map<String, Node> currentNodeMap = nodeMap.get(m);

            for (JobClass oldClass : classes) {
                String newName = prefix + oldClass.getName();
                JobClass newClass = null;

                if (oldClass instanceof OpenClass) {
                    OpenClass oldOpen = (OpenClass) oldClass;
                    newClass = new OpenClass(unionNetwork, newName, oldOpen.getPriority());
                } else if (oldClass instanceof SelfLoopingClass) {
                    SelfLoopingClass oldSelf = (SelfLoopingClass) oldClass;
                    // Map reference station to new node
                    String oldRefStatName = oldSelf.getReferenceStation().getName();
                    if (currentNodeMap.containsKey(oldRefStatName)) {
                        Node newRefStat = currentNodeMap.get(oldRefStatName);
                        if (newRefStat instanceof Station) {
                            newClass = new SelfLoopingClass(unionNetwork, newName,
                                    (long) oldSelf.getPopulation(), (Station) newRefStat, oldSelf.getPriority());
                        }
                    }
                    if (newClass == null) {
                        line_warning(mfilename(new Object() {}),
                                "Reference station " + oldRefStatName + " not found for class " + oldClass.getName());
                        continue;
                    }
                } else if (oldClass instanceof ClosedClass) {
                    ClosedClass oldClosed = (ClosedClass) oldClass;
                    // Map reference station to new node
                    String oldRefStatName = oldClosed.getReferenceStation().getName();
                    if (currentNodeMap.containsKey(oldRefStatName)) {
                        Node newRefStat = currentNodeMap.get(oldRefStatName);
                        if (newRefStat instanceof Station) {
                            newClass = new ClosedClass(unionNetwork, newName,
                                    oldClosed.getPopulation(), (Station) newRefStat, oldClosed.getPriority());
                        }
                    }
                    if (newClass == null) {
                        line_warning(mfilename(new Object() {}),
                                "Reference station " + oldRefStatName + " not found for class " + oldClass.getName());
                        continue;
                    }
                } else {
                    line_warning(mfilename(new Object() {}),
                            "Unknown class type: " + oldClass.getClass().getSimpleName() + ", skipping");
                    continue;
                }

                currentClassMap.put(oldClass.getName(), newClass);
            }
        }

        // Phase 4: Set service and arrival distributions
        for (int m = 0; m < models.size(); m++) {
            Network model = models.get(m);
            List<Node> nodes = model.getNodes();
            List<JobClass> classes = model.getClasses();
            Map<String, Node> currentNodeMap = nodeMap.get(m);
            Map<String, JobClass> currentClassMap = classMap.get(m);

            for (Node oldNode : nodes) {
                // Handle Source arrivals
                if (oldNode instanceof Source) {
                    Source oldSource = (Source) oldNode;
                    for (JobClass oldClass : classes) {
                        if (oldClass instanceof OpenClass) {
                            if (currentClassMap.containsKey(oldClass.getName())) {
                                JobClass newClass = currentClassMap.get(oldClass.getName());
                                Distribution dist = oldSource.getArrivalDistribution(oldClass);
                                if (dist != null && !(dist instanceof Disabled)) {
                                    unionSource.setArrival(newClass, dist);
                                }
                            }
                        }
                    }
                    continue;
                }

                // Skip Sink
                if (oldNode instanceof Sink) {
                    continue;
                }

                // Skip auto-added ClassSwitch
                if (oldNode instanceof ClassSwitch && ((ClassSwitch) oldNode).autoAdded) {
                    continue;
                }

                // Get mapped node
                if (!currentNodeMap.containsKey(oldNode.getName())) {
                    continue;
                }
                Node newNode = currentNodeMap.get(oldNode.getName());

                // Set service distributions for Queue/Delay
                if (oldNode instanceof ServiceStation) {
                    ServiceStation oldStation = (ServiceStation) oldNode;
                    ServiceStation newStation = (ServiceStation) newNode;
                    for (JobClass oldClass : classes) {
                        if (!currentClassMap.containsKey(oldClass.getName())) {
                            continue;
                        }
                        JobClass newClass = currentClassMap.get(oldClass.getName());

                        try {
                            Distribution dist = oldStation.getServiceProcess(oldClass);
                            if (dist != null && !(dist instanceof Disabled)) {
                                newStation.setService(newClass, dist);
                            }
                        } catch (Exception e) {
                            // Service not defined for this class, skip
                        }
                    }
                }
            }
        }

        // Phase 5: Build routing matrix from linked routing matrices
        RoutingMatrix P = unionNetwork.initRoutingMatrix();
        int nUnionClasses = unionNetwork.getNumberOfClasses();
        int nUnionNodes = unionNetwork.getNumberOfNodes();

        for (int m = 0; m < models.size(); m++) {
            Network thisModel = models.get(m);
            List<Node> nodes = thisModel.getNodes();
            List<JobClass> classes = thisModel.getClasses();
            Map<String, Node> currentNodeMap = nodeMap.get(m);
            Map<String, JobClass> currentClassMap = classMap.get(m);

            // Get the linked routing matrix from this model
            Map<JobClass, Map<JobClass, Matrix>> Pm = null;
            try {
                Pm = thisModel.getLinkedRoutingMatrix();
            } catch (Exception e) {
                Pm = null;
            }

            if (Pm == null || Pm.isEmpty()) {
                // Fall back to outputStrategy-based extraction if no linked matrix
                // This is a simplified fallback - the linked matrix approach is preferred
                continue;
            }

            // Use the linked routing matrix (includes inter-class routing)
            int nModelClasses = classes.size();
            int nModelNodes = nodes.size();

            for (JobClass oldClassR : classes) {
                if (!currentClassMap.containsKey(oldClassR.getName())) {
                    continue;
                }
                JobClass newClassR = currentClassMap.get(oldClassR.getName());

                Map<JobClass, Matrix> PmR = Pm.get(oldClassR);
                if (PmR == null) {
                    continue;
                }

                for (JobClass oldClassS : classes) {
                    if (!currentClassMap.containsKey(oldClassS.getName())) {
                        continue;
                    }
                    JobClass newClassS = currentClassMap.get(oldClassS.getName());

                    Matrix Prs = PmR.get(oldClassS);
                    if (Prs == null || Prs.isEmpty()) {
                        continue;
                    }

                    // Copy routing probabilities, mapping old node indices to new
                    for (int i = 0; i < Math.min(nModelNodes, Prs.getNumRows()); i++) {
                        Node oldNodeI = nodes.get(i);

                        // Skip auto-added ClassSwitch nodes
                        if (oldNodeI instanceof ClassSwitch && ((ClassSwitch) oldNodeI).autoAdded) {
                            continue;
                        }

                        if (!currentNodeMap.containsKey(oldNodeI.getName())) {
                            continue;
                        }
                        Node newNodeI = currentNodeMap.get(oldNodeI.getName());

                        for (int j = 0; j < Math.min(nModelNodes, Prs.getNumCols()); j++) {
                            double prob = Prs.get(i, j);
                            if (prob == 0) {
                                continue;
                            }

                            Node oldNodeJ = nodes.get(j);

                            // Skip auto-added ClassSwitch nodes
                            if (oldNodeJ instanceof ClassSwitch && ((ClassSwitch) oldNodeJ).autoAdded) {
                                continue;
                            }

                            if (!currentNodeMap.containsKey(oldNodeJ.getName())) {
                                continue;
                            }
                            Node newNodeJ = currentNodeMap.get(oldNodeJ.getName());

                            P.set(newClassR, newClassS, newNodeI, newNodeJ, prob);
                        }
                    }
                }
            }
        }

        // Link the network with the routing matrix
        unionNetwork.link(P);

        return unionNetwork;
    }

    /**
     * Alias for merge - creates a union Network from this ensemble.
     *
     * @param ensemble the Ensemble object to convert
     * @return Network object containing merged subnetworks
     */
    public static Network toNetwork(Ensemble ensemble) {
        return merge(ensemble);
    }

    /**
     * Helper class for storing pending Join node creation info.
     */
    private static class JoinPendingItem {
        int modelIdx;
        Join oldNode;
        String newName;

        JoinPendingItem(int modelIdx, Join oldNode, String newName) {
            this.modelIdx = modelIdx;
            this.oldNode = oldNode;
            this.newName = newName;
        }
    }
}
