/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.gen;

import jline.lang.*;
import jline.lang.constant.*;
import jline.lang.nodes.*;
import jline.lang.processes.*;
import jline.util.matrix.Matrix;
import jline.util.Maths;
import jline.util.graph.DirectedGraph;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * A generator object that generates queueing network models
 * based on user specification. Characteristics of generated
 * models can be configured via the generator's properties.
 */
public class NetworkGenerator {
    
    private String schedStrat;
    private String routingStrat;
    private String distribution;
    private String cclassJobLoad;
    private boolean hasVaryingServiceRates;
    private boolean hasMultiServerQueues;
    private boolean hasRandomCSNodes;
    private boolean hasMultiChainCS;
    private boolean initializeStates;
    private Function<Integer, Matrix> topologyFcn;
    
    private static final int MAX_SERVERS = 40;
    private static final int[] HIGH_JOB_LOAD_RANGE = {31, 40};
    private static final int[] MED_JOB_LOAD_RANGE = {11, 20};
    private static final int[] LOW_JOB_LOAD_RANGE = {1, 5};
    
    private final Random random = new Random();
    
    /**
     * Constructor with default settings
     */
    public NetworkGenerator() {
        this("randomize", "randomize", "randomize", "randomize", 
             true, true, true, true, true, NetworkGenerator::randGraph);
    }
    
    /**
     * Constructor with custom settings
     */
    public NetworkGenerator(String schedStrat, String routingStrat, String distribution,
                          String cclassJobLoad, boolean hasVaryingServiceRates,
                          boolean hasMultiServerQueues, boolean hasRandomCSNodes,
                          boolean hasMultiChainCS, Function<Integer, Matrix> topologyFcn) {
        this(schedStrat, routingStrat, distribution, cclassJobLoad,
             hasVaryingServiceRates, hasMultiServerQueues, hasRandomCSNodes,
             hasMultiChainCS, true, topologyFcn);
    }
    
    /**
     * Constructor with all options including state initialization
     */
    public NetworkGenerator(String schedStrat, String routingStrat, String distribution,
                          String cclassJobLoad, boolean hasVaryingServiceRates,
                          boolean hasMultiServerQueues, boolean hasRandomCSNodes,
                          boolean hasMultiChainCS, boolean initializeStates,
                          Function<Integer, Matrix> topologyFcn) {
        setSchedStrat(schedStrat);
        setRoutingStrat(routingStrat);
        setDistribution(distribution);
        setCclassJobLoad(cclassJobLoad);
        this.hasVaryingServiceRates = hasVaryingServiceRates;
        this.hasMultiServerQueues = hasMultiServerQueues;
        this.hasRandomCSNodes = hasRandomCSNodes;
        this.hasMultiChainCS = hasMultiChainCS;
        this.initializeStates = initializeStates;
        setTopologyFcn(topologyFcn);
    }
    
    /**
     * Main function to call. Returns a generated QN model according to
     * specified properties of the NetworkGenerator object
     */
    public Network generate() {
        return generate(random.nextInt(8) + 1, null, 0, random.nextInt(4) + 1);
    }
    
    public Network generate(int numQueues) {
        return generate(numQueues, null, 0, random.nextInt(4) + 1);
    }
    
    public Network generate(int numQueues, Integer numDelays) {
        return generate(numQueues, numDelays, 0, random.nextInt(4) + 1);
    }
    
    public Network generate(int numQueues, Integer numDelays, int numOClass) {
        return generate(numQueues, numDelays, numOClass, random.nextInt(4) + 1);
    }
    
    public Network generate(int numQueues, Integer numDelays, int numOClass, int numCClass) {
        if (numDelays == null) {
            if (numQueues > 1) {
                numDelays = random.nextInt(2);
            } else {
                numDelays = 1;
            }
        }
        
        validateArgs(numQueues, numDelays, numOClass, numCClass);
        
        Network model = new Network("nw");
        List<Station> stations = createStations(model, numQueues, numDelays, numOClass > 0);
        List<JobClass> classes = createClasses(model, numOClass, numCClass);
        setServiceProcesses(stations, classes);
        defineTopology(model, topologyFcn.apply(stations.size()));
        
        // Initialize default states if requested
        if (initializeStates) {
            initDefaultStates(model);
        }
        
        return model;
    }
    
    // Property setters with validation
    public void setSchedStrat(String strat) {
        if (strat.equals("fcfs") || strat.equals("ps") || strat.equals("inf") ||
            strat.equals("lcfs") || strat.equals("lcfspr") || strat.equals("siro") ||
            strat.equals("sjf") || strat.equals("ljf") || strat.equals("sept") ||
            strat.equals("lept") || strat.equals("randomize")) {
            this.schedStrat = strat;
        } else {
            throw new IllegalArgumentException("Scheduling strategy does not exist or is not supported");
        }
    }
    
    public void setRoutingStrat(String strat) {
        if (strat.equals("Probabilities") || strat.equals("Random") || strat.equals("randomize")) {
            this.routingStrat = strat;
        } else {
            throw new IllegalArgumentException("Routing strategy does not exist or is not supported");
        }
    }
    
    public void setDistribution(String distrib) {
        if (distrib.equalsIgnoreCase("Exp") || distrib.equalsIgnoreCase("HyperExp") ||
            distrib.equalsIgnoreCase("Erlang") || distrib.equalsIgnoreCase("randomize")) {
            this.distribution = distrib;
        } else {
            throw new IllegalArgumentException("Distribution does not exist or is not supported");
        }
    }
    
    public void setCclassJobLoad(String load) {
        if (load.equalsIgnoreCase("high") || load.equalsIgnoreCase("medium") ||
            load.equalsIgnoreCase("low") || load.equalsIgnoreCase("randomize")) {
            this.cclassJobLoad = load;
        } else {
            throw new IllegalArgumentException("Model load can only be high, medium, or low");
        }
    }
    
    public void setTopologyFcn(Function<Integer, Matrix> fcn) {
        try {
            Matrix adj = fcn.apply(2);
            if (adj == null || adj.getNumRows() != 2 || adj.getNumCols() != 2) {
                throw new IllegalArgumentException("topologyFcn should take a positive integer and return an adjacency matrix");
            }
        } catch (Exception e) {
            throw new IllegalArgumentException("topologyFcn should take a positive integer and return an adjacency matrix");
        }
        this.topologyFcn = fcn;
    }
    
    // Getters
    public String getSchedStrat() { return schedStrat; }
    public String getRoutingStrat() { return routingStrat; }
    public String getDistribution() { return distribution; }
    public String getCclassJobLoad() { return cclassJobLoad; }
    public boolean isHasVaryingServiceRates() { return hasVaryingServiceRates; }
    public boolean isHasMultiServerQueues() { return hasMultiServerQueues; }
    public boolean isHasRandomCSNodes() { return hasRandomCSNodes; }
    public boolean isHasMultiChainCS() { return hasMultiChainCS; }
    public boolean isInitializeStates() { return initializeStates; }
    
    public void setInitializeStates(boolean initializeStates) {
        this.initializeStates = initializeStates;
    }
    
    // Private methods
    
    private void validateArgs(int numQueues, int numDelays, int numOClass, int numCClass) {
        if (numQueues < 0 || numDelays < 0 || numOClass < 0 || numCClass < 0) {
            throw new IllegalArgumentException("Arguments must be non-negative");
        }
        if (numQueues + numDelays <= 0 || numOClass + numCClass <= 0) {
            throw new IllegalArgumentException("At least one station and one job class required");
        }
    }
    
    private List<Station> createStations(Network model, int numQueues, int numDelays, boolean hasOClass) {
        List<Station> stations = new ArrayList<>();
        
        // Create queues
        for (int i = 0; i < numQueues; i++) {
            jline.lang.nodes.Queue queue = new jline.lang.nodes.Queue(model, name("queue", i + 1), chooseSchedStrat());
            queue.setNumberOfServers(chooseNumServers());
            stations.add(queue);
        }
        
        // Create delays
        for (int i = 0; i < numDelays; i++) {
            Delay delay = new Delay(model, name("delay", i + 1));
            stations.add(delay);
        }
        
        // Create source and sink if open classes exist
        if (hasOClass) {
            new Source(model, "source");
            new Sink(model, "sink");
        }
        
        return stations;
    }
    
    private List<JobClass> createClasses(Network model, int numOClass, int numCClass) {
        List<JobClass> classes = new ArrayList<>();
        
        // Create open classes
        for (int i = 0; i < numOClass; i++) {
            OpenClass openClass = new OpenClass(model, name("OClass", i + 1));
            model.getSource().setArrival(openClass, chooseDistribution());
            classes.add(openClass);
        }
        
        // Create closed classes
        if (numCClass > 0) {
            List<Station> stations = model.getStations();
            Station refStat = stations.get(random.nextInt(stations.size()));
            
            for (int i = 0; i < numCClass; i++) {
                ClosedClass closedClass = new ClosedClass(model, name("CClass", i + 1), 
                                                         chooseNumJobs(), refStat);
                classes.add(closedClass);
            }
        }
        
        return classes;
    }
    
    private void setServiceProcesses(List<Station> stations, List<JobClass> classes) {
        for (JobClass jobClass : classes) {
            for (Station station : stations) {
                ((ServiceStation)station).setService(jobClass, chooseDistribution());
            }
        }
    }
    
    private void defineTopology(Network model, Matrix topology) {
        int numStations = model.getNumberOfStations();
        Matrix csMask = genCSMask(model);
        
        // Add source and sink to topology graph if needed
        if (model.hasOpenClasses()) {
            Matrix newTopology = new Matrix(topology.getNumRows() + 2, topology.getNumCols() + 2);
            // Copy existing edges
            for (int i = 0; i < topology.getNumRows(); i++) {
                for (int j = 0; j < topology.getNumCols(); j++) {
                    newTopology.set(i, j, topology.get(i, j));
                }
            }
            // Add source and sink connections
            int sourceIdx = model.getIndexSourceNode();
            int sinkIdx = model.getIndexSinkNode();
            newTopology.set(sourceIdx, random.nextInt(numStations - 1), 1.0);
            newTopology.set(random.nextInt(numStations - 1), sinkIdx, 1.0);
            topology = newTopology;
        }
        
        // Set up routing for each station
        for (int i = 0; i < numStations; i++) {
            List<Integer> destIDs = new ArrayList<>();
            for (int j = 0; j < topology.getNumCols(); j++) {
                if (topology.get(i, j) > 0) {
                    destIDs.add(j);
                }
            }
            Collections.sort(destIDs);
            List<Node> outgoingNodes = addOutgoingLinks(model, i, destIDs, csMask);
            setRoutingStrategies(model, model.getStationByIndex(i), outgoingNodes);
        }
    }
    
    private Matrix genCSMask(Network model) {
        int numOClasses = 0;
        int numCClasses = 0;
        for (JobClass jc : model.getClasses()) {
            if (jc instanceof OpenClass) {
                numOClasses++;
            } else if (jc instanceof ClosedClass) {
                numCClasses++;
            }
        }
        int totalClasses = numOClasses + numCClasses;
        Matrix mask = new Matrix(totalClasses, totalClasses);
        
        if (!hasMultiChainCS) {
            // Open classes can switch among themselves
            for (int i = 0; i < numOClasses; i++) {
                for (int j = 0; j < numOClasses; j++) {
                    mask.set(i, j, 1.0);
                }
            }
            // Closed classes can switch among themselves
            for (int i = numOClasses; i < totalClasses; i++) {
                for (int j = numOClasses; j < totalClasses; j++) {
                    mask.set(i, j, 1.0);
                }
            }
            return mask;
        }
        
        // Multi-chain class switching logic
        List<Integer> allChains = new ArrayList<>();
        if (model.hasOpenClasses()) {
            allChains.addAll(assignChains(numOClasses, random.nextInt(numOClasses) + 1));
        }
        if (model.hasClosedClasses()) {
            allChains.addAll(assignChains(numCClasses, random.nextInt(numCClasses) + 1));
        }
        
        int startIdx = 0;
        for (int chainSize : allChains) {
            int endIdx = startIdx + chainSize;
            for (int i = startIdx; i < endIdx; i++) {
                for (int j = startIdx; j < endIdx; j++) {
                    mask.set(i, j, 1.0);
                }
            }
            startIdx = endIdx;
        }
        
        return mask;
    }
    
    private List<Integer> assignChains(int numClasses, int numChains) {
        return randintfixedsum(numClasses, numChains);
    }
    
    private List<Integer> randintfixedsum(int s, int n) {
        if (n == 1) {
            return Arrays.asList(s);
        } else if (s == n) {
            return Collections.nCopies(n, 1);
        }
        
        int first = random.nextInt(s - n) + 1;
        List<Integer> rest = randintfixedsum(s - first, n - 1);
        List<Integer> result = new ArrayList<>();
        result.add(first);
        result.addAll(rest);
        
        // Shuffle the result
        Collections.shuffle(result);
        return result;
    }
    
    private List<Node> addOutgoingLinks(Network model, int sourceID, List<Integer> destIDs, Matrix csMask) {
        Node sourceNode = model.getNodeByIndex(sourceID);
        List<Node> outgoingNodes = new ArrayList<>();
        
        for (int destID : destIDs) {
            Node destNode = model.getNodeByIndex(destID);
            
            if (hasRandomCSNodes && random.nextBoolean() && 
                !(sourceNode instanceof Source) && !(destNode instanceof Sink)) {
                ClassSwitch cs = randClassSwitchNode(model, csMask, sourceNode, destNode);
                outgoingNodes.add(cs);
                model.addLink(sourceNode, cs);
                model.addLink(cs, destNode);
                
                // Set routing from class switch to destination
                for (JobClass jobClass : model.getClasses()) {
                    cs.setProbRouting(jobClass, destNode, 1.0);
                }
            } else {
                outgoingNodes.add(destNode);
                model.addLink(sourceNode, destNode);
            }
        }
        
        return outgoingNodes;
    }
    
    private ClassSwitch randClassSwitchNode(Network model, Matrix mask, Node sourceNode, Node destNode) {
        String name = name("cs", sourceNode.getName() + "_" + destNode.getName());
        Matrix csMatrix = randClassSwitchMatrix(mask);
        return new ClassSwitch(model, name, csMatrix);
    }
    
    private Matrix randClassSwitchMatrix(Matrix mask) {
        int size = mask.getNumRows();
        Matrix matrix = new Matrix(size, size);
        
        for (int i = 0; i < size; i++) {
            List<Integer> validIndices = new ArrayList<>();
            for (int j = 0; j < size; j++) {
                if (mask.get(i, j) > 0) {
                    validIndices.add(j);
                }
            }
            
            if (!validIndices.isEmpty()) {
                double[] probs = randfixedsumone(validIndices.size());
                for (int k = 0; k < validIndices.size(); k++) {
                    matrix.set(i, validIndices.get(k), probs[k]);
                }
            }
        }
        
        return matrix;
    }
    
    private void setRoutingStrategies(Network model, Station station, List<Node> outgoingNodes) {
        if (station instanceof Source) {
            // Source nodes route all open classes with probability 1
            for (JobClass jobClass : model.getClasses()) {
                if (jobClass instanceof OpenClass) {
                    for (Node node : outgoingNodes) {
                        station.setProbRouting(jobClass, node, 1.0);
                    }
                }
            }
            return;
        }
        
        for (JobClass jobClass : model.getClasses()) {
            String strat = chooseRoutingStrat();
            // Normalize routing strategies for compatibility with OutputStrategy
            String stratEnum = strat.toUpperCase();
            if (stratEnum.equals("PROBABILITIES")) {
                stratEnum = "PROB";
            } else if (stratEnum.equals("RANDOM")) {
                stratEnum = "RAND";
            }
            station.setRouting(jobClass, RoutingStrategy.valueOf(stratEnum));
            
            if (strat.equals("Probabilities")) {
                double[] probs;
                
                // Closed classes cannot route to sink
                if (jobClass instanceof ClosedClass && 
                    outgoingNodes.get(outgoingNodes.size() - 1) instanceof Sink) {
                    station.setProbRouting(jobClass, outgoingNodes.get(outgoingNodes.size() - 1), 0.0);
                    probs = randfixedsumone(outgoingNodes.size() - 1);
                } else {
                    probs = randfixedsumone(outgoingNodes.size());
                }
                
                for (int j = 0; j < probs.length; j++) {
                    station.setProbRouting(jobClass, outgoingNodes.get(j), probs[j]);
                }
            }
        }
    }
    
    private double[] randfixedsumone(int numElems) {
        if (numElems == 0) return new double[0];
        if (numElems == 1) return new double[]{1.0};
        
        // Generate random values and normalize
        double[] values = new double[numElems];
        double sum = 0;
        for (int i = 0; i < numElems; i++) {
            values[i] = random.nextDouble();
            sum += values[i];
        }
        
        // Normalize to sum to 1
        for (int i = 0; i < numElems; i++) {
            values[i] = Math.ceil(values[i] / sum * 1000) / 1000;
        }
        
        // Adjust largest element to ensure exact sum of 1
        int maxIdx = 0;
        for (int i = 1; i < numElems; i++) {
            if (values[i] > values[maxIdx]) {
                maxIdx = i;
            }
        }
        
        double currentSum = 0;
        for (double v : values) {
            currentSum += v;
        }
        values[maxIdx] -= (currentSum - 1.0);
        
        return values;
    }
    
    private String name(String prefix, Object suffix) {
        if (prefix.equalsIgnoreCase("cs")) {
            return prefix + "_" + suffix;
        } else {
            return prefix + suffix;
        }
    }
    
    private int chooseNumServers() {
        if (hasMultiServerQueues) {
            return random.nextInt(MAX_SERVERS) + 1;
        } else {
            return 1;
        }
    }
    
    private int chooseNumJobs() {
        switch (cclassJobLoad.toLowerCase()) {
            case "high":
                return random.nextInt(HIGH_JOB_LOAD_RANGE[1] - HIGH_JOB_LOAD_RANGE[0] + 1) + HIGH_JOB_LOAD_RANGE[0];
            case "medium":
                return random.nextInt(MED_JOB_LOAD_RANGE[1] - MED_JOB_LOAD_RANGE[0] + 1) + MED_JOB_LOAD_RANGE[0];
            case "low":
                return random.nextInt(LOW_JOB_LOAD_RANGE[1] - LOW_JOB_LOAD_RANGE[0] + 1) + LOW_JOB_LOAD_RANGE[0];
            case "randomize":
            default:
                return random.nextInt(HIGH_JOB_LOAD_RANGE[1]) + 1;
        }
    }
    
    private SchedStrategy chooseSchedStrat() {
        if (schedStrat.equalsIgnoreCase("randomize")) {
            int id = random.nextInt(2) + 1;
            switch (id) {
                case 1:
                    return SchedStrategy.FCFS;
                case 2:
                default:
                    return SchedStrategy.PS;
            }
        } else {
            return SchedStrategy.valueOf(schedStrat.toUpperCase());
        }
    }
    
    private String chooseRoutingStrat() {
        if (routingStrat.equalsIgnoreCase("randomize")) {
            int id = random.nextInt(2) + 1;
            switch (id) {
                case 1:
                    return "Random";
                case 2:
                default:
                    return "Probabilities";
            }
        } else {
            return routingStrat;
        }
    }
    
    private Distribution chooseDistribution() {
        int id;
        switch (distribution.toLowerCase()) {
            case "exp":
                id = 1;
                break;
            case "erlang":
                id = 2;
                break;
            case "hyperexp":
                id = 3;
                break;
            case "randomize":
            default:
                id = random.nextInt(3) + 1;
                break;
        }
        
        double mean;
        if (hasVaryingServiceRates) {
            mean = Math.pow(2, random.nextInt(13) - 6);
        } else {
            mean = 1.0;
        }
        
        switch (id) {
            case 1:
                return new Exp(1.0 / mean);
            case 2:
                int k = (int)Math.pow(2, random.nextInt(7));
                return new Erlang(k / mean, k);
            case 3:
            default:
                // Create a 2-phase hyperexponential with high variance
                double scv = Math.pow(2, random.nextInt(7));
                return HyperExp.fitMeanAndSCV(mean, scv);
        }
    }
    
    /**
     * Generate a random strongly connected graph topology
     * This implements the algorithm from MATLAB's randGraph function
     */
    public static Matrix randGraph(int numVertices) {
        if (numVertices <= 0) {
            throw new IllegalArgumentException("Number of vertices must be positive");
        }
        
        if (numVertices == 1) {
            Matrix adj = new Matrix(1, 1);
            adj.set(0, 0, 1.0);
            return adj;
        }
        
        // Generate random spanning tree
        Matrix tree = randSpanningTree(numVertices);
        
        // Apply strongConnect algorithm to ensure strong connectivity
        RandGraphState state = new RandGraphState(numVertices);
        Matrix strongGraph = strongConnect(tree, 0, state);
        
        // Randomly permute nodes
        List<Integer> perm = IntStream.range(0, numVertices).boxed().collect(Collectors.toList());
        Collections.shuffle(perm);
        Matrix permuted = new Matrix(numVertices, numVertices);
        for (int i = 0; i < numVertices; i++) {
            for (int j = 0; j < numVertices; j++) {
                if (strongGraph.get(i, j) > 0) {
                    permuted.set(perm.get(i), perm.get(j), 1.0);
                }
            }
        }
        
        return permuted;
    }
    
    /**
     * Helper class to maintain state for randGraph algorithm
     */
    private static class RandGraphState {
        int globalStartTime = 0;
        int[] startTime;
        int[] lowestLink;
        int[] invStartTime;
        Random rand = new Random();
        
        RandGraphState(int n) {
            startTime = new int[n];
            lowestLink = new int[n];
            invStartTime = new int[n];
            Arrays.fill(startTime, -1);
            Arrays.fill(lowestLink, -1);
            Arrays.fill(invStartTime, -1);
        }
    }
    
    /**
     * DFS-based algorithm to ensure strong connectivity
     */
    private static Matrix strongConnect(Matrix g, int v, RandGraphState state) {
        state.startTime[v] = state.globalStartTime;
        state.invStartTime[state.startTime[v]] = v;
        state.lowestLink[v] = state.startTime[v];
        state.globalStartTime++;
        
        // Find outgoing edges from v
        List<Integer> outVertices = new ArrayList<>();
        for (int w = 0; w < g.getNumCols(); w++) {
            if (g.get(v, w) > 0) {
                outVertices.add(w);
            }
        }
        
        for (int w : outVertices) {
            if (state.startTime[w] == -1) {
                g = strongConnect(g, w, state);
                state.lowestLink[v] = Math.min(state.lowestLink[v], state.lowestLink[w]);
            } else {
                state.lowestLink[v] = Math.min(state.lowestLink[v], state.startTime[w]);
            }
        }
        
        // If v is root of SCC but not the entire graph's root
        if (state.lowestLink[v] == state.startTime[v] && state.startTime[v] > 0) {
            int descendantST = state.rand.nextInt(state.globalStartTime - state.startTime[v]) + state.startTime[v];
            int ancestorST = state.rand.nextInt(state.startTime[v]);
            // Add edge to ensure strong connectivity
            g.set(state.invStartTime[descendantST], state.invStartTime[ancestorST], 1.0);
            state.lowestLink[v] = ancestorST;
        }
        
        return g;
    }
    
    /**
     * Generate a random spanning tree
     */
    private static Matrix randSpanningTree(int numVertices) {
        Matrix tree = new Matrix(numVertices, numVertices);
        Random rand = new Random();
        
        for (int i = 1; i < numVertices; i++) {
            int parent = rand.nextInt(i);
            tree.set(parent, i, 1.0);
        }
        
        return tree;
    }
    
    /**
     * Generate a cyclic graph topology
     */
    public static Matrix cyclicGraph(int numVertices) {
        if (numVertices <= 0) {
            throw new IllegalArgumentException("Number of vertices must be positive");
        }
        
        Matrix adj = new Matrix(numVertices, numVertices);
        
        // Create a cycle
        for (int i = 0; i < numVertices - 1; i++) {
            adj.set(i, i + 1, 1.0);
        }
        adj.set(numVertices - 1, 0, 1.0);
        
        return adj;
    }
    
    /**
     * Initialize default states for all nodes in the network
     * Based on MATLAB's initDefaultCustom function
     */
    public static void initDefaultStates(Network model) {
        NetworkStruct sn = model.getStruct(false);
        int R = sn.nclasses;
        Matrix N = sn.njobs.transpose();
        
        for (int i = 0; i < model.getNumberOfNodes(); i++) {
            Node node = model.getNodeByIndex(i);
            
            if (node instanceof Station) {
                // Initialize station states
                double[] n0 = new double[R];
                double[] s0 = new double[R];
                
                // Find the station index for this node
                int stationIdx = -1;
                for (int j = 0; j < sn.nstations; j++) {
                    if ((int) sn.stationToNode.get(0, j) == i) {
                        stationIdx = j;
                        break;
                    }
                }
                
                if (stationIdx == -1) {
                    continue; // Skip if station mapping not found
                }
                
                double s = sn.nservers.get(stationIdx, 0);
                
                // Allocate jobs for closed classes
                for (int r = 0; r < R; r++) {
                    if (Double.isFinite(N.get(r, 0))) {
                        if (stationIdx == (int)sn.refstat.get(r, 0)) {
                            n0[r] = N.get(r, 0);
                        }
                        s0[r] = Math.min(n0[r], s);
                        s = s - s0[r];
                    }
                }
                
                // Create state from marginal
                int[] state = new int[2 * R];
                for (int r = 0; r < R; r++) {
                    state[r] = (int) n0[r];
                    state[R + r] = (int) s0[r];
                }
                
                if (node instanceof StatefulNode) {
                    Matrix stateMatrix = new Matrix(1, state.length);
                    for (int j = 0; j < state.length; j++) {
                        stateMatrix.set(0, j, state[j]);
                    }
                    ((StatefulNode) node).setState(stateMatrix);
                }
            } else if (node instanceof StatefulNode) {
                // Initialize non-station stateful nodes
                if (node instanceof Cache) {
                    // Initialize cache state
                    int[] state = new int[R];
                    Matrix stateMatrix = new Matrix(1, state.length);
                    for (int k = 0; k < state.length; k++) {
                        stateMatrix.set(0, k, state[k]);
                    }
                    ((Cache) node).setState(stateMatrix);
                } else if (node instanceof Router) {
                    // Initialize router state
                    ((Router) node).setState(0);
                }
            }
        }
    }
}