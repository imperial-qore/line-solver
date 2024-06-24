package jline.lang;

import jline.api.DTMC;
import jline.api.SN;
import jline.lang.processes.MAP;
import jline.lang.processes.Replayer;
import jline.util.Maths;
import jline.lang.constant.*;
import jline.lang.distributions.*;
import jline.lang.nodes.Queue;
import jline.lang.nodes.*;
import jline.lang.sections.*;
import jline.lang.state.State;
import jline.solvers.SolverHandles;
import jline.solvers.jmt.SolverJMT;
import jline.util.*;
import org.apache.commons.lang3.NotImplementedException;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

import static java.lang.Double.isFinite;
import static jline.lib.KPCToolbox.*;

/**
 * A queueing network model
 */
public class Network extends Model implements Serializable {
    private boolean doChecks;
    private boolean hasState;
    private String logPath;
    private boolean usedFeatures;

    private List<Node> nodes;
    private final List<JobClass> jobClasses;
    private final List<Station> stations;

    private boolean hasStruct;
    private NetworkStruct sn;
    private Matrix csMatrix;
    private Matrix connections;

    private List<Matrix> handles;
    private List<FiniteCapacityRegion> regions;

    // caches
    private Map<Node, Map<JobClass, List<Node>>> classLinks;

    private final NetworkAttribute attribute;
    private final List<ItemSet> items;

    public Network(String modelName) {
        super(modelName);

        this.hasState = false;
        this.doChecks = true;

        this.nodes = new ArrayList<Node>();
        this.jobClasses = new ArrayList<JobClass>();
        this.stations = new ArrayList<Station>();

        this.classLinks = new HashMap<Node, Map<JobClass, List<Node>>>();
        this.items = new ArrayList<>();
        this.regions = new ArrayList<>();

        this.hasStruct = false;
        this.csMatrix = null;
        this.sn = null;
        this.connections = null;
        this.attribute = new NetworkAttribute();
    }

    public void setDoChecks(boolean doChecks) {
        this.doChecks = doChecks;
    }

    public int[] getSize() {
        int[] outInt = new int[2];
        outInt[0] = this.getNumberOfNodes();
        outInt[1] = this.getNumberOfClasses();
        return outInt;
    }

    public boolean hasOpenClasses() {
        for (JobClass temp : this.jobClasses) {
            if (temp instanceof OpenClass) {
                return true;
            }
        }

        return false;
    }

    public boolean hasProductFormSolution() {
        return SN.snHasProductForm(this.getStruct(false));
    }

    public int getJobClassIndex(JobClass jobClass) {
        return this.jobClasses.indexOf(jobClass);
    }

    public JobClass getJobClassFromIndex(int inIdx) {
        return this.jobClasses.get(inIdx);
    }

    public List<Integer> getIndexOpenClasses() {
        List<Integer> outList = new ArrayList<Integer>();
        for (int i = 0; i < this.jobClasses.size(); i++) {
            if (this.jobClasses.get(i) instanceof OpenClass) {
                outList.add(this.getJobClassIndex(this.jobClasses.get(i)));
            }
        }
        return outList;
    }

    public List<Integer> getIndexStatefulNodes() {
        List<Integer> outList = new ArrayList<Integer>();
        for (int i = 0; i < this.nodes.size(); i++) {
            if (this.nodes.get(i) instanceof StatefulNode) {
                outList.add(i);
            }
        }
        return outList;
    }

    public boolean hasClosedClasses() {
        for (JobClass temp : this.jobClasses) {
            if (temp instanceof ClosedClass) {
                return true;
            }
        }
        return false;
    }

    public List<Integer> getIndexClosedClasses() {
        List<Integer> outList = new ArrayList<Integer>();
        for (int i = 0; i < this.jobClasses.size(); i++) {
            if (this.jobClasses.get(i) instanceof ClosedClass) {
                outList.add(this.getJobClassIndex(this.jobClasses.get(i)));
            }
        }
        return outList;
    }

    public boolean hasClasses() {
        return !this.jobClasses.isEmpty();
    }

    public List<JobClass> getClasses() {
        return this.jobClasses;
    }

    public List<Node> getNodes() {
        return this.nodes;
    }

    public List<Station> getStations() {
        return this.stations;
    }

    public void addNode(Node node) {
        node.setNodeIdx(node.getNodeIdx()); // searches within the model nodes
        nodes.add(node);
        if (node instanceof Station) {
            node.setStationIdx(node.getStationIdx()); // searches within the model stations
            stations.add((Station) node);
        }
    }

    public void setInitialized(boolean initStatus) {
        this.hasState = initStatus;
    }

    public int getNumberOfNodes() {
        return this.nodes.size();
    }

    public int getNumberOfStations() {
        return this.stations.size();
    }

    public int getNumberOfStatefulNodes() {
        int ct = 0;
        for (Node node : this.nodes) {
            if (node instanceof StatefulNode) {
                ct++;
            }
        }
        return ct;
    }

    public int getNumberOfClasses() {
        return this.jobClasses.size();
    }

    public JobClass getClassByName(String name) {
        for (JobClass jobClass : this.jobClasses) {
            if (jobClass.getName().equals(name)) {
                return jobClass;
            }
        }
        return null;
    }

    public JobClass getClassByIndex(int index) {
        for (JobClass jobClass : this.jobClasses) {
            if (this.getJobClassIndex(jobClass) == index) {
                return jobClass;
            }
        }
        return null;
    }

    public void addJobClass(JobClass jobClass) {
        this.jobClasses.add(jobClass);
    }

    public Node getNodeByName(String name) {
        for (Node node : this.nodes) {
            if (node.getName().equals(name)) {
                return node;
            }
        }

        return null;
    }

    public Node getNodeByIndex(int idx) {
        int nodesPassed = 0;
        for (Node nodeIter : this.nodes) {
            if (nodesPassed == idx) {
                return nodeIter;
            }
            nodesPassed++;
        }
        return null;
    }

    public Node getNodeByStatefulIndex(int idx) {
        int nodesPassed = 0;
        for (Node nodeIter : this.nodes) {
            if (nodeIter instanceof StatefulNode) {
                if (nodesPassed == idx) {
                    return nodeIter;
                }
                nodesPassed++;
            }
        }

        return null;
    }

    public int getNodeIndex(Node node) {
        return this.nodes.indexOf(node);
    }

    public int getNodeIndex(String name) {
        for (Node node : this.nodes) {
            if (node.getName().equals(name))
                return getNodeIndex(node);
        }
        return -1;
    }

    public int getStatefulNodeIndex(Node node) {
        if (!(node instanceof StatefulNode))
            return -1;

        int outIdx = 0;
        for (Node nodeIter : this.nodes) {
            if (nodeIter == node) {
                return outIdx;
            } else if (nodeIter instanceof StatefulNode) {
                outIdx++;
            }
        }

        return -1;
    }

    public Node getStatefulNodeFromIndex(int inIdx) {
        int outIdx = inIdx;
        for (Node nodeIter : this.nodes) {
            if (nodeIter instanceof StatefulNode) {
                if (outIdx == 0) {
                    return nodeIter;
                }
                outIdx--;
            }
        }

        return null;
    }

    public int getStationIndex(Node node) {
        return this.stations.indexOf(node);
    }

    public Node getStationFromIndex(int inIdx) {
        return this.stations.get(inIdx);
    }

    public static jline.lang.RoutingMatrix serialRouting(List<JobClass> jobClasses, Node... nodes) {
        if (nodes.length == 0) {
            return new jline.lang.RoutingMatrix();
        }

        Network network = nodes[0].model;
        jline.lang.RoutingMatrix outMatrix = new jline.lang.RoutingMatrix(network, jobClasses, network.nodes);

        for (int i = 1; i < nodes.length; i++) {
            //System.out.format("Loading connection %s->%s\n", nodes[i-1].getName(), nodes[i].getName());
            outMatrix.addConnection(nodes[i - 1], nodes[i], 1.0);
        }

        if (!(nodes[nodes.length - 1] instanceof Sink)) {
            outMatrix.addConnection(nodes[nodes.length - 1], nodes[0], 1.0);
        }

        return outMatrix;
    }

    public static jline.lang.RoutingMatrix serialRouting(List<JobClass> jobClasses, List<Node> nodes) {
        if (nodes.isEmpty()) {
            return new jline.lang.RoutingMatrix();
        }

        Network network = nodes.get(0).model;
        jline.lang.RoutingMatrix outMatrix = new jline.lang.RoutingMatrix(network, jobClasses, nodes);

        for (int i = 1; i < nodes.size(); i++) {
            //System.out.format("Loading connection %s->%s\n", nodes[i-1].getName(), nodes[i].getName());
            outMatrix.addConnection(nodes.get(i - 1), nodes.get(i), 1.0);
        }

        if (!(nodes.get(nodes.size() - 1) instanceof Sink)) {
            outMatrix.addConnection(nodes.get(nodes.size() - 1), nodes.get(0), 1.0);
        }

        return outMatrix;
    }

    public static jline.lang.RoutingMatrix serialRouting(JobClass jobClass, Node... nodes) {
        List<JobClass> jobClasses = new ArrayList<JobClass>();
        jobClasses.add(jobClass);

        return Network.serialRouting(jobClasses, nodes);
    }

    public static jline.lang.RoutingMatrix serialRouting(JobClass jobClass, List<Node> nodes) {
        List<JobClass> jobClasses = new ArrayList<JobClass>();
        jobClasses.add(jobClass);

        return Network.serialRouting(jobClasses, nodes);
    }

    public static jline.lang.RoutingMatrix serialRouting(Node... nodes) {
        if (nodes.length == 0) {
            return new jline.lang.RoutingMatrix();
        }
        Network network = nodes[0].model;
        return Network.serialRouting(network.jobClasses, nodes);
    }

    public static jline.lang.RoutingMatrix serialRouting(java.util.List<Node> nodes) {
        if (nodes.size() == 0) {
            return new jline.lang.RoutingMatrix();
        }
        Network network = nodes.get(0).model;
        return Network.serialRouting(network.jobClasses, nodes);
    }

    public void link(jline.lang.RoutingMatrix routing) {
        /*
             Input:
                routing: row: source, column: dest
         */
        sanitize();
        routing.setRouting(this);
    }

    public void unLink() {
        for (Node node : this.nodes) {
            node.resetRouting();
        }
    }

    public void addLink(Node sourceNode, Node destNode) {
        int sourceNodeIdx = this.getNodeIndex(sourceNode);
        int destNodeIdx = this.getNodeIndex(destNode);
        this.addLink(sourceNodeIdx, destNodeIdx);
    }

    public void addLink(int sourceNodeIdx, int destNodeIdx) {
        if (this.connections == null || this.connections.isEmpty())
            this.connections = new Matrix(nodes.size(), nodes.size());

        if (this.connections.getNumRows() != this.nodes.size())
            this.connections.expandMatrix(this.nodes.size(), this.nodes.size(), this.nodes.size() * this.nodes.size());

        this.connections.set(sourceNodeIdx, destNodeIdx, 1.0);
    }

    public void addLinks(Node[][] links) {
        for (Node[] linkPair : links) {
            if (linkPair.length == 2) {
                addLink(linkPair[0], linkPair[1]);
            } else {
                throw new IllegalArgumentException("Each link pair must contain exactly two objects.");
            }
        }
    }

    // Get initial state
    public State getState() {
        if (!this.hasInitState()) {
            this.initDefault();
        }

        for (int i = 0; i < this.getNumberOfNodes(); i++) {
            if (this.nodes.get(i).isStateful()) {
                Matrix initialState = ((StatefulNode) this.nodes.get(i)).getState();
                Matrix priorInitialState = ((StatefulNode) this.nodes.get(i)).getStatePrior();
                this.sn.state.put(((StatefulNode) this.nodes.get(i)), initialState);
                this.sn.stateprior.put(((StatefulNode) this.nodes.get(i)), priorInitialState);
            }
        }

        return new State(this.sn.state, this.sn.stateprior);
    }

    public boolean hasInitState() {
        boolean output = true;
        if (!this.hasState) { // check if all stations are initialized
            for (int i = 0; i < this.getNumberOfNodes(); i++) {
                if (this.nodes.get(i) instanceof StatefulNode) {
                    if (((StatefulNode) this.nodes.get(i)).getState().isEmpty()) {
                        output = false;
                    }
                }
            }
        }
        return output;
    }

    public void initDefault() {
        // TODO: is it necessary to have a version where, per LINE, nodes can be passed in as a parameter?
        // open classes are empty
        // closed classes are initialized at reference station
        // running jobs are allocated in class id order until all servers are busy

        NetworkStruct sn = this.getStruct(false);
        int R = sn.nclasses;
        Matrix N = sn.njobs.transpose();

        for (int i = 0; i < this.getNumberOfNodes(); i++) {
            if (sn.isstation.get(i, 0) == 1) {
                Matrix n0 = new Matrix(1, N.length());
                n0.zero();
                Matrix s0 = new Matrix(1, N.length());
                s0.zero();
                double s = sn.nservers.get((int) sn.nodeToStation.get(0, i), 0); // allocate

                for (int r = 0; r < N.getNumRows(); r++) {
                    if (isFinite(N.get(r, 0))) { // for all closed classes
                        if (sn.nodeToStation.get(0, i) == sn.refstat.get(r, 0)) {
                            n0.set(0, r, N.get(r, 0));
                        }
                    }
                    s0.set(0, r, Math.min(n0.get(0, r), s));
                    s -= s0.get(0, r);
                }

                Matrix state_i = State.fromMarginalAndStarted(sn, i, n0, s0);

                switch (sn.nodetypes.get(i)) {
                    case Cache:
                        int nVarsInt = (int) sn.nvars.get(i, 2 * R);
                        Matrix newState_i = new Matrix(1, state_i.getNumCols() + nVarsInt);
                        for (int p = 0; p < state_i.length(); p++) {
                            newState_i.set(0, p, state_i.get(0, p));
                        }
                        int addition = 0;
                        for (int p = state_i.length(); p < newState_i.length(); p++) {
                            newState_i.set(0, p, addition);
                            addition++;
                        }
                        state_i = newState_i.clone();
                        break;

                    case Place:
                        state_i.zero(); // for now PNs are single class
                        break;

                    default:
                        if (sn.isstation.get(i, 0) == 1) {
                            for (int r = 0; r < sn.nclasses; r++) {
                                if (sn.proctype.get(sn.nodes.get(i)).get(sn.jobclasses.get(r)) == ProcessType.MAP) {
                                    Matrix one = new Matrix(1, 1, 1);
                                    one.set(0, 0, 1);
                                    state_i = Matrix.decorate(state_i, one);
                                }
                            }
                        }
                }

                for (int r = 0; r < sn.nclasses; r++) {
                    if ((sn.routing.get(sn.nodes.get(i)).get(sn.jobclasses.get(r))
                            == RoutingStrategy.RROBIN)
                            || (sn.routing.get(sn.nodes.get(i)).get(sn.jobclasses.get(r))
                            == RoutingStrategy.WRROBIN)) {
                        // Start from first connected queue
                        List<Integer> findSnRt = new ArrayList<>();
                        for (int p = 0; p < sn.rt.getNumCols(); p++) {
                            if (sn.rt.get(i, p) == 1) {
                                findSnRt.add(p);
                            }
                        }
                        Matrix newState_i = new Matrix(1, state_i.getNumCols() + findSnRt.size());
                        for (int p = 0; p < state_i.length(); p++) {
                            newState_i.set(0, p, state_i.get(0, p));
                        }
                        for (int p = state_i.length(); p < newState_i.length(); p++) {
                            newState_i.set(0, p, findSnRt.get(p - state_i.length()));
                        }
                        state_i = newState_i.clone();
                        break;
                    }
                }

                if (state_i.isEmpty()) {
                    System.err.format("Default initialisation failed on station %d", i);
                } else {
                    ((StatefulNode) this.nodes.get(i)).setState(state_i);

                    Matrix prior_state_i = new Matrix(state_i.getNumRows(), 1);
                    prior_state_i.zero();
                    prior_state_i.set(0, 0, 1);
                    ((StatefulNode) this.nodes.get(i)).setStatePrior(prior_state_i);
                }
            } else if (sn.isstateful.get(i, 0) == 1) { // Not a station
                if (this.nodes.get(i) instanceof Cache) {
                    Cache cacheNode = (Cache) this.nodes.get(i);
                    Matrix state_i = new Matrix(1, this.getNumberOfClasses() + (int) cacheNode.getItemLevelCap().elementSum());
                    for (int idx = this.getNumberOfClasses(); idx < state_i.getNumCols(); idx++) {
                        state_i.set(idx, idx + 1);
                    }
                    cacheNode.setState(state_i);
                } else if (this.nodes.get(i) instanceof Router) {
                    Matrix one = new Matrix(1, 1, 1);
                    one.set(0, 0, 1);
                    //Router temporarily changed from StatefulNode to Node as SSA does not work otherwise
                    //((Router)this.nodes.get(i)).setState(one);
                } else {
                    ((StatefulNode) this.nodes.get(i)).setState(new Matrix(0, 0));
                }
            }
        }

        if (this.isStateValid()) {
            this.hasState = true;
        } else {
            System.err.println("Default initialisation failed.");
        }
    }

    private boolean isStateValid() {
        // This code in LINE is found in api/sn/snIsStateValid.m
        // It is only ever called by isStateValid method in Network
        // For that reason, and to avoid unnecessary complication, I've consolidated into one method in JLINE

        // Modified so not using this.sn
        NetworkStruct snTmp = this.getStruct(true);
        Matrix nir = new Matrix(snTmp.nstations, snTmp.nclasses);
        Matrix sir = new Matrix(snTmp.nstations, snTmp.nclasses);

        for (int ist = 0; ist < snTmp.nstations; ist++) {
            int isf = (int) snTmp.stationToStateful.get(0, ist);
            if (snTmp.state.get(snTmp.stations.get(ist)).getNumRows() > 1) {
                if (snTmp.stateprior.get(snTmp.stations.get(ist)).elementMax() < 1 - GlobalConstants.FineTol) {
                    System.err.format("isStateValid will ignore some states of station %d, define a unique initial state to address this problem.\n", ist);
                }
                Matrix initialState = new Matrix(1, snTmp.state.get(snTmp.stations.get(ist)).getNumCols());
                Matrix.extractRows(snTmp.state.get(snTmp.stations.get(ist)), 0, 1, initialState);
                snTmp.state.put(snTmp.stations.get(ist), initialState);
            }

            State.StateMarginalStatistics stats = State.toMarginal(snTmp,
                    (int) snTmp.stationToNode.get(0, ist),
                    snTmp.state.get(snTmp.stations.get(ist)),
                    null, null, null, null, null);

            for (int i = 0; i < snTmp.nclasses; i++) {
                nir.set(ist, i, stats.nir.get(0, i));
                sir.set(ist, i, stats.sir.get(0, i));
            }
        }

        return State.isValid(snTmp, nir, sir);
    }

    public void printSummary() {
        System.out.format("jline.Network model: %s\n", this.getName());
        System.out.format("--------------------------------------------------------\n");
        System.out.format("Job classes: \n");
        for (JobClass jobClass : this.jobClasses) {
            jobClass.printSummary();
        }
        System.out.format("--------------------------------------------------------\n");
        System.out.format("Nodes: \n");
        for (Node node : this.nodes) {
            node.printSummary();
            System.out.format("--------\n");
        }
    }

    public void clearCaches() {
        this.classLinks = new HashMap<Node, Map<JobClass, List<Node>>>();
    }

    protected void generateClassLinks() {
        this.classLinks = new HashMap<Node, Map<JobClass, List<Node>>>();
        for (Node node : this.nodes) {
            Map<JobClass, List<Node>> nodeMap = new HashMap<JobClass, List<Node>>();

            for (JobClass jobClass : this.jobClasses) {
                nodeMap.put(jobClass, new ArrayList<Node>());
            }
            classLinks.put(node, nodeMap);
        }

        for (Node node : this.nodes) {
            for (final OutputStrategy outputStrategy : node.getOutputStrategies()) {
                Node destNode = outputStrategy.getDestination();
                if (destNode == null) {
                    continue;
                }
                JobClass jobClass = outputStrategy.getJobClass();
                this.classLinks.get(destNode).get(jobClass).add(node);
            }
        }
    }

    public int getClassLinks(Node node, JobClass jobClass) {
        if (this.classLinks.isEmpty()) {
            this.generateClassLinks();
        }
        return this.classLinks.get(node).get(jobClass).size();
    }

    public double minRate() {
        double acc = Double.POSITIVE_INFINITY;
        for (Node node : this.nodes) {
            if (node instanceof HasSchedStrategy) {
                double accVal = ((HasSchedStrategy) node).minRate();
                if (accVal != 0) {
                    acc = Math.min(acc, accVal);
                }
            }
        }
        return acc;
    }

    public double maxRate() {
        double acc = 0;
        for (Node node : this.nodes) {
            if (node instanceof HasSchedStrategy) {
                double accVal = ((HasSchedStrategy) node).maxRate();
                if (accVal != Double.POSITIVE_INFINITY) {
                    acc = Math.max(acc, accVal);
                }
            }
        }
        return acc;
    }

    public double avgRate() {
        double acc = 0;
        int accCt = 0;
        for (Node node : this.nodes) {
            if (node instanceof HasSchedStrategy) {
                double accVal = ((HasSchedStrategy) node).avgRate();
                int valCt = ((HasSchedStrategy) node).rateCt();
                if ((accVal != Double.POSITIVE_INFINITY) && (accVal != 0)) {
                    acc += accVal;
                    accCt += valCt;
                }
            }
        }
        return acc / accCt;
    }

    public Matrix getCsMatrix() {
        return this.csMatrix;
    }

    public void setCsMatrix(Matrix csMatrix) {
        this.csMatrix = csMatrix;
    }

    public Matrix getConnectionMatrix() {
        if (this.connections.getNumCols() < this.getNumberOfNodes() ||
                this.connections.getNumRows() < this.getNumberOfNodes())
            this.connections.expandMatrix(this.getNumberOfNodes(), this.getIndexSourceNode(), this.getNumberOfNodes() * this.getNumberOfNodes());
        return this.connections;
    }

    public void setConnectionMatrix(Matrix connection) {
        this.connections = connection;
    }

    public Matrix getForkJoins() {
        int I = this.getNumberOfNodes();
        Matrix fjPairs = new Matrix(I, I);

        for (int i = 0; i < I; i++) {
            Node node = this.nodes.get(i);
            if (node instanceof Fork) {
                //no-op
            } else if (node instanceof Join) {
                fjPairs.set(((Join) node).joinOf.getNodeIdx(), node.getNodeIdx(), 1.0);
            }
        }
        return fjPairs;
    }

    public boolean getHasStruct() {
        return this.hasStruct;
    }

    public void setHasStruct(boolean hasStruct) {
        this.hasStruct = hasStruct;
    }

    public void setStruct(NetworkStruct sn) {
        this.sn = sn;
    }

    public NetworkStruct getStruct() {
        return this.getStruct(true);
    }

    public NetworkStruct getStruct(boolean wantInitialState) {
        if (!this.hasStruct)
            refreshStruct(true);

        if (wantInitialState)
            getState();

        return this.sn;
    }

    public void refreshStruct(boolean hardRefresh) {

        sanitize();

        List<NodeType> nodetypes;
        List<String> classnames;
        List<String> nodenames;
        Matrix refstat;
        Matrix conn;
        Matrix njobs;
        Matrix numservers;
        Matrix lldscaling;
        Map<Station, SerializableFunction<Matrix, Double>> cdscaling;
        Map<Node, Map<JobClass, RoutingStrategy>> routing;


        if (this.hasStruct && !hardRefresh) {
            nodetypes = sn.nodetypes;
            classnames = sn.classnames;
            nodenames = sn.nodenames;
            refstat = sn.refstat;
        } else {
            nodetypes = getNodeTypes();
            classnames = getClassNames();
            nodenames = getNodeNames();
            refstat = getReferenceStations();
        }

        conn = getConnectionMatrix();
        njobs = getNumberOfJobs();
        numservers = getStationServers();
        lldscaling = getLimitedLoadDependence();
        cdscaling = getLimitedClassDependence();

        if (sn == null)
            sn = new NetworkStruct();

        sn.nnodes = nodenames.size();
        sn.nclasses = classnames.size();
        sn.stations = this.stations;
        sn.stateful = new ArrayList<>();
        sn.jobclasses = this.jobClasses;
        sn.nodes = this.nodes;

        routing = new HashMap<Node, Map<JobClass, RoutingStrategy>>();
        for (Node node : this.nodes) {
            if (node.isStateful()) {
                sn.stateful.add((StatefulNode) node);
            }
            Map<JobClass, RoutingStrategy> map = new HashMap<JobClass, RoutingStrategy>();
            for (JobClass jobclass : this.jobClasses) {
                map.put(jobclass, getRoutingStrategyFromNodeAndClassPair(node, jobclass));
            }
            routing.put(node, map);
        }

        sn.nclosedjobs = DoubleStream.of(njobs.getNonZeroValues()).boxed().filter(val -> !Double.isInfinite(val)).reduce(0.0, (a, b) -> a + b);
        sn.nservers = numservers;
        sn.isstation = getIsStationArray();
        sn.nstations = stations.size();
        sn.nodetypes = nodetypes;
        sn.scv = new Matrix(sn.nstations, sn.nclasses, sn.nstations * sn.nclasses);
        sn.scv.fill(1.0);
        sn.njobs = njobs.transpose();
        sn.refstat = refstat;
        sn.space = new HashMap<Station, Matrix>();
        for (int i = 0; i < sn.nstations; i++)
            sn.space.put(stations.get(i), new Matrix(0, 0, 0));
        sn.routing = routing;
        sn.chains = new Matrix(0, 0);
        sn.lst = null;
        sn.lldscaling = lldscaling;
        sn.cdscaling = cdscaling;
        sn.nodetypes = nodetypes;
        sn.isstateful = getIsStatefulArray();
        sn.isstatedep = new Matrix(sn.nnodes, 3, 3 * sn.nnodes);
        for (int i = 0; i < sn.nnodes; i++) {
            //Line 72-80 is ignored since JLine not support cache node
            Node node = this.nodes.get(i);
            for (int j = 0; j < sn.nclasses; j++) {
                JobClass jobclass = this.jobClasses.get(j);
                switch (sn.routing.get(node).get(jobclass)) {
                    case RROBIN:
                    case WRROBIN:
                    case JSQ:
                        sn.isstatedep.set(i, 2, 1.0);
                    default:
                        continue;
                }
            }
        }
        sn.nstateful = getNumberOfStatefulNodes();
        sn.state = new HashMap<StatefulNode, Matrix>(sn.nstations);
        sn.stateprior = new HashMap<StatefulNode, Matrix>(sn.nstations);
        sn.space = new HashMap<Station, Matrix>(sn.nstations);
        for (int i = 0; i < sn.nstations; i++)
            sn.state.put(stations.get(i), new Matrix(0, 0, 0));
        sn.nodenames = nodenames;
        sn.classnames = classnames;
        sn.connmatrix = conn;

        sn.eventCache = new HashMap<>();

        //line 97-108 is ignored since for transition node

        sn.nodeToStateful = new Matrix(1, nodes.size(), nodes.size());
        sn.nodeToStation = new Matrix(1, nodes.size(), nodes.size());
        sn.stationToNode = new Matrix(1, stations.size(), stations.size());
        sn.stationToStateful = new Matrix(1, stations.size(), stations.size());
        sn.statefulToNode = new Matrix(1, sn.nstateful, sn.nstateful);
        for (int i = 0; i < nodes.size(); i++) {
            sn.nodeToStateful.set(0, i, nodes.get(i).getStatefulIdx());
            sn.nodeToStation.set(0, i, nodes.get(i).getStationIdx());
        }
        for (int i = 0; i < stations.size(); i++) {
            sn.stationToNode.set(0, i, stations.get(i).getNodeIdx());
            sn.stationToStateful.set(0, i, stations.get(i).getStatefulIdx());
        }
        for (int isf = 0; isf < sn.nstateful; isf++) {
            sn.statefulToNode.set(0, isf, getStatefulNodeFromIndex(isf).getNodeIdx());
        }

        refreshPriorities();
        refreshService(null, null);

        refreshChains(!sn.nodetypes.contains(NodeType.Cache));

        for (int c = 0; c < sn.nclasses; c++) {
            if (this.jobClasses.get(c) instanceof SelfLoopingClass) {
                sn.isslc.set(c, 0, 1.0);
            }
        }

        Matrix refclasses = this.getReferenceClasses();
        Matrix refclass = new Matrix(1, sn.nchains);
        for (int c = 0; c < sn.nchains; c++) {
            Matrix inchain_c = sn.inchain.get(c);
            Matrix find_refclasses = refclasses.find();
            List<Double> isect = Matrix.intersect(inchain_c, find_refclasses); // can have a single element
            if (!isect.isEmpty()) {
                refclass.set(0, c, isect.get(0));
            } else {
                refclass.set(0, c, -1);
            }
        }

        this.sn.refclass = refclass;
        this.sn.fj = this.getForkJoins();

        refreshLocalVars();
        refreshSync();
        //refreshPetriNetNodes()
        this.hasStruct = true;
    }

    public void sanitize() {
        //THIS FUNCTION IS NOT TESTED
        if (this.sn == null) {
            int M = this.stations.size();
            int K = this.jobClasses.size();
            for (int i = 0; i < this.nodes.size(); i++) {
                Node node = this.nodes.get(i);
                if (node instanceof Cache) {
                    Cache cache = (Cache) node;
                    for (int k = 0; k < K; k++) {
                        if (k >= cache.popularityLength() || cache.popularityGet(k) == null) {
                            cache.popularitySet(k, Disabled.getInstance());
                        }
                    }
                    if (cache.accessProb == null || cache.accessProb.length == 0) {
                        cache.accessProb = new Matrix[K][cache.getItems().getNumberOfItems()];
                        for (int v = 0; v < K; v++) {
                            for (int k = 0; k < cache.getItems().getNumberOfItems(); k++) {
                                // accessProb[v][k](l,p) is the cost (probability) for a user-v request to item k in list l to access list p
                                if (cache.getGraph() == null) {
                                    Matrix diag = new Matrix(cache.getnLevels() + 1, cache.getnLevels() + 1);
                                    for (int j = 0; j < cache.getnLevels(); j++) {
                                        diag.set(j, j + 1, 1);
                                    }
                                    cache.accessProb[v][k] = diag;
                                    cache.accessProb[v][k].set(cache.getnLevels(), cache.getnLevels(), 1);
                                } else {
                                    cache.accessProb[v][k] = cache.getGraph()[k];
                                }
                            }
                        }
                    }
                    for (int r = 0; r < cache.getCacheServer().hitClass.getNumRows(); r++) {
                        for (int c = 0; c < cache.getCacheServer().hitClass.getNumCols(); c++) {
                            cache.getCacheServer().hitClass.set(r, c, Math.round(cache.getCacheServer().hitClass.get(r, c)));
                        }
                    }
                    for (int r = 0; r < cache.getCacheServer().missClass.getNumRows(); r++) {
                        for (int c = 0; c < cache.getCacheServer().missClass.getNumCols(); c++) {
                            cache.getCacheServer().missClass.set(r, c, Math.round(cache.getCacheServer().missClass.get(r, c)));
                        }
                    }
                } else if (node instanceof Logger) {
                    //do nothing
                } else if (node instanceof ClassSwitch) {
                    //do nothing
                } else if (node instanceof Join) {
                    Join join = (Join) node;
                    for (int k = 0; k < K; k++) {
                        JobClass jobclass = this.jobClasses.get(k);
                        join.setClassCap(jobclass, Double.POSITIVE_INFINITY);
                        join.setDropRule(jobclass, DropStrategy.WaitingQueue);
                    }
                } else if (node instanceof Delay) {
                    Delay delay = (Delay) node;
                    for (int k = 0; k < K; k++) {
                        JobClass jobclass = this.jobClasses.get(k);
                        if (!delay.getServer().containsJobClass(jobclass)) {
                            delay.setService(jobclass, new Disabled(), 0);
                            delay.setClassCap(jobclass, 0);
                            delay.getInput().setInputJobProcess(new InputBinding(jobclass, SchedStrategyType.NP, DropStrategy.WaitingQueue));
                        }
                    }

                    if (Objects.requireNonNull(delay.getSchedStrategy()) == SchedStrategy.SEPT) {
                        ArrayList<Double> svcTime = new ArrayList<Double>();
                        for (int k = 0; k < K; k++)
                            svcTime.add(delay.getServiceProcess(this.jobClasses.get(k)).getMean());
                        Collections.sort(svcTime);

                        for (int k = 0; k < K; k++)
                            delay.setSchedStrategyPar(this.jobClasses.get(k), svcTime.get(k));
                    } else {
                        continue;
                    }

                } else if (node instanceof Queue) {
                    Queue queue = (Queue) node;
                    for (int k = 0; k < K; k++) {
                        JobClass jobclass = this.jobClasses.get(k);
                        if (!queue.getServer().containsJobClass(jobclass)) {
                            queue.setService(jobclass, new Disabled(), 0);
                            queue.setClassCap(jobclass, 0);
                            queue.getInput().setInputJobProcess(new InputBinding(jobclass, SchedStrategyType.NP, DropStrategy.WaitingQueue));
                        }
                    }

                    switch (queue.getSchedStrategy()) {
                        case SEPT:
                            ArrayList<Double> svcTime = new ArrayList<Double>();
                            for (int k = 0; k < K; k++)
                                svcTime.add(queue.getServiceProcess(this.jobClasses.get(k)).getMean());

                            if (svcTime.stream().distinct().collect(Collectors.toList()).size() != K)
                                throw new RuntimeException("SEPT does not support identical service time means.");

                            ArrayList<Double> svcTimeSorted = new ArrayList<Double>(svcTime);
                            Collections.sort(svcTimeSorted);
                            for (int k = 0; k < K; k++)
                                queue.setSchedStrategyPar(this.jobClasses.get(k), svcTimeSorted.indexOf(svcTime.get(k)) + 1);
                            break;
                        case LEPT:
                            svcTime = new ArrayList<Double>();
                            for (int k = 0; k < K; k++)
                                svcTime.add(queue.getServiceProcess(this.jobClasses.get(k)).getMean());

                            if (svcTime.stream().distinct().collect(Collectors.toList()).size() != K)
                                throw new RuntimeException("SEPT does not support identical service time means.");

                            svcTimeSorted = new ArrayList<Double>(svcTime);
                            Collections.sort(svcTimeSorted, Collections.reverseOrder());
                            for (int k = 0; k < K; k++)
                                queue.setSchedStrategyPar(this.jobClasses.get(k), svcTimeSorted.indexOf(svcTime.get(k)) + 1);
                            break;
                        default:
                            continue;
                    }
                } else if (node instanceof Sink) {
                    //do nothing
                } else if (node instanceof Source) {
                    Source source = (Source) node;
                    for (int k = 0; k < K; k++) {
                        JobClass jobclass = this.jobClasses.get(k);
                        if (!source.containsJobClass(jobclass))
                            source.setArrival(jobclass, new Disabled());
                    }
                } else if (node instanceof Router) {
                    // Do nothing
                } else if (node instanceof Place) {
                    throw new NotImplementedException("Place node not yet supported in JLINE");
                } else if (node instanceof Transition) {
                    throw new NotImplementedException("Transition node not yet supported in JLINE");
                }
            }

            int sourceIdx = this.getIndexSourceNode();
            for (int i = 0; i < M; i++) {
                if ((sourceIdx == -1) || (i != sourceIdx)) {
                    for (int r = 0; r < K; r++) {
                        ServiceSection server = this.stations.get(i).getServer();
                        if (server instanceof ServiceTunnel) {
                            //do nothing
                        } //else if (server instanceof CacheClassSwitcher) {}
                        else {
                            if (!this.stations.get(i).getServer().containsJobClass(this.jobClasses.get(r)))
                                this.stations.get(i).getServer().setServiceProcesses(new ServiceBinding(this.jobClasses.get(r), ServiceStrategy.LI, new Disabled()));
                        }
                    }
                }
            }
        }
    }

    public List<NodeType> getNodeTypes() {
        int M = getNumberOfNodes();
        List<NodeType> nodetypes = new ArrayList<NodeType>(M);

        try {
            for (int i = 0; i < M; i++) {
                Node nodeIter = this.nodes.get(i);
                if (nodeIter instanceof Logger)
                    nodetypes.add(NodeType.Logger);
                else if (nodeIter instanceof ClassSwitch)
                    nodetypes.add(NodeType.ClassSwitch);
                else if (nodeIter instanceof Join)
                    nodetypes.add(NodeType.Join);
                else if (nodeIter instanceof Sink)
                    nodetypes.add(NodeType.Sink);
                else if (nodeIter instanceof Router)
                    nodetypes.add(NodeType.Router);
                else if (nodeIter instanceof Delay)
                    nodetypes.add(NodeType.Delay);
                else if (nodeIter instanceof Fork)
                    nodetypes.add(NodeType.Fork);
                else if (nodeIter instanceof Queue)
                    nodetypes.add(NodeType.Queue);
                else if (nodeIter instanceof Source)
                    nodetypes.add(NodeType.Source);
//				Below node types are not supported in JLine
//        		else if (nodeIter instanceof Place)
//        			nodetypes.add(NodeType.Place);
//        		else if (nodeIter instanceof Transition)
//        			nodetypes.add(NodeType.Transition);
                else if (nodeIter instanceof Cache)
                    nodetypes.add(NodeType.Cache);
                else
                    throw new Exception("Unknown node type.");
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        return nodetypes;
    }

    public List<String> getClassNames() {
        if (hasStruct && sn.classnames != null)
            return sn.classnames;

        int K = getNumberOfClasses();
        List<String> classnames = new ArrayList<String>();
        for (int i = 0; i < K; i++)
            classnames.add(jobClasses.get(i).getName());

        return classnames;
    }

    public List<String> getNodeNames() {
        if (hasStruct && sn.classnames != null)
            return sn.nodenames;

        int M = getNumberOfNodes();
        List<String> nodenames = new ArrayList<String>();
        for (int i = 0; i < M; i++)
            nodenames.add(nodes.get(i).getName());

        return nodenames;
    }

    public Matrix getReferenceStations() {
        int K = getNumberOfClasses();
        Matrix refstat = new Matrix(K, 1, K);

        for (int i = 0; i < K; i++) {
            if (jobClasses.get(i).type == JobClassType.Open) {
                refstat.set(i, 0, getIndexSourceStation());
            } else {
                ClosedClass cc = (ClosedClass) jobClasses.get(i);
                refstat.set(i, 0, getStationIndex(cc.getRefstat()));
            }
        }


        return refstat;
    }

    public Matrix getReferenceClasses() {
        int K = this.jobClasses.size();
        Matrix refclass = new Matrix(K, 1);
        for (int i = 0; i < K; i++) {
            if (this.jobClasses.get(i).isReferenceClass())
                refclass.set(i, 0, 1.0);
        }
        return refclass;
    }

    public int getIndexSourceNode() {
        int res = 0;
        for (Node nodeIter : this.nodes) {
            if (nodeIter instanceof Source)
                return res;
            res++;
        }

        return -1;
    }

    public Source getSource() {
        int index = this.getIndexSourceNode();
        if (index == -1) {
            System.err.println("The given model does not have a source");
            return null;
        }
        return (Source) this.nodes.get(index);
    }

    public Sink getSink() {
        int index = this.getIndexSinkNode();
        if (index == -1) {
            System.err.println("The given model does not have a sink");
            return null;
        }
        return (Sink) this.nodes.get(index);
    }

    public int getIndexSinkNode() {
        int res = 0;
        for (Node nodeIter : this.nodes) {
            if (nodeIter instanceof Sink)
                return res;
            res++;
        }

        return -1;
    }

    public Matrix getNumberOfJobs() {
        int K = getNumberOfClasses();
        Matrix njobs = new Matrix(K, 1, K);
        for (int i = 0; i < K; i++) {
            if (jobClasses.get(i).type == JobClassType.Open)
                njobs.set(i, 0, Double.POSITIVE_INFINITY);
            else if (jobClasses.get(i).type == JobClassType.Closed)
                njobs.set(i, 0, jobClasses.get(i).getNumberOfJobs());
        }
        return njobs;
    }

    public Matrix getStationServers() {
        int I = stations.size();
        Matrix numservers = new Matrix(I, 1, I);
        for (int i = 0; i < I; i++) {
            if (stations.get(i).getNumberOfServers() == Integer.MAX_VALUE)
                numservers.set(i, 0, Double.POSITIVE_INFINITY);
            else
                numservers.set(i, 0, stations.get(i).getNumberOfServers());
        }

        return numservers;
    }

    public Matrix getLimitedLoadDependence() {
        List<Matrix> mus = new ArrayList<Matrix>();
        int maxsize = 0;

        for (Station station : this.stations) {
            mus.add(station.getLimitedLoadDependence());
            maxsize = Math.max(maxsize, station.getLimitedLoadDependence().length());
        }

        int M = this.stations.size();
        Matrix alpha = new Matrix(M, maxsize);
        alpha.fill(1.0);
        for (int i = 0; i < M; i++) {
            Matrix mu = mus.get(i);
            if (mu.length() > 0) {
                Matrix.extract(mu, 0, 1, 0, mu.length(), alpha, i, 0);
                for (int j = 0; j < mu.length(); j++) {
                    if (alpha.get(i, j) == 0)
                        alpha.set(i, j, 1.0);
                }
            }
        }
        return alpha;
    }

    public Map<Station, SerializableFunction<Matrix, Double>> getLimitedClassDependence() {
        Map<Station, SerializableFunction<Matrix, Double>> gamma = new HashMap<Station, SerializableFunction<Matrix, Double>>();

        for (Station station : this.stations) {
            if (station.getLimitedClassDependence() != null)
                gamma.put(station, station.getLimitedClassDependence());
        }

        if (gamma.size() > 0) {
            for (Station station : this.stations) {
                if (gamma.getOrDefault(station, null) == null) {
                    gamma.put(station, (nvec) -> {
                        return 1.0;
                    });
                }
            }
            return gamma;
        } else {
            //Set to null which represents the cell(nstations, 0) in matlab
            return null;
        }
    }

    public RoutingStrategy getRoutingStrategyFromNodeAndClassPair(Node node, JobClass c) {
        //Another approach is to use the last routing strategy
        if (node.getOutput().getOutputStrategyByClass(c).size() == 0)
            if (node instanceof Sink) {
                return RoutingStrategy.DISABLED;
            } else {
                return RoutingStrategy.RAND;
            }
        RoutingStrategy res = null;
        try {
            for (OutputStrategy outputStrategy : node.getOutputStrategies()) {
                if (outputStrategy.getJobClass().equals(c)) {
                    if (res == null)
                        res = outputStrategy.getRoutingStrategy();
                    else if (!res.equals(outputStrategy.getRoutingStrategy()))
                        throw new Exception("Inconsistent routing strategy.");
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        if (res != null)
            return res;
        else if (c instanceof OpenClass)
            return RoutingStrategy.RAND;
        else
            return RoutingStrategy.DISABLED;
    }

    private Matrix getIsStationArray() {
        Matrix isStation = new Matrix(nodes.size(), 1, nodes.size());
        for (int i = 0; i < nodes.size(); i++) {
            if (nodes.get(i) instanceof Station)
                isStation.set(i, 0, 1);
        }
        return isStation;
    }

    private Matrix getIsStatefulArray() {
        Matrix isStateful = new Matrix(nodes.size(), 1, nodes.size());
        for (int i = 0; i < nodes.size(); i++) {
            if (nodes.get(i) instanceof StatefulNode)
                isStateful.set(i, 0, 1);
        }
        return isStateful;
    }

    public void setSn(NetworkStruct sn) {
        this.sn = sn;
    }

    public void refreshPriorities() {
        int K = this.jobClasses.size();
        Matrix classprio = new Matrix(1, K, K);
        for (int i = 0; i < K; i++) {
            classprio.set(0, i, jobClasses.get(i).priority);
        }

        if (this.sn != null)
            sn.classprio = classprio;
    }

    public void refreshService(List<Integer> statSet, List<Integer> classSet) {
        boolean[] status = refreshRates(statSet, classSet);
        boolean hasSCVChanged = status[1];
        boolean hasRateChanged = status[0];

        if (hasSCVChanged) {
            refreshServiceTypes(statSet, classSet);
            refreshServicePhases(statSet, classSet);
            refreshLST(statSet, classSet);
        }

        if (this.sn.sched == null) {
            refreshScheduling();
        } else {
            for (Station station : this.stations) {
                SchedStrategy schedStrategy = this.sn.sched.getOrDefault(station, null);
                if (schedStrategy == SchedStrategy.SEPT || schedStrategy == SchedStrategy.LEPT) {
                    refreshScheduling();
                    break;
                }
            }
        }
    }

    public void refreshJobs() {
        refreshStruct(true);
        Matrix njobs = getNumberOfJobs();
        double njobsSum = 0;
        for (int j = 0; j < njobs.length(); j++) {
            List<Double> njobsList = njobs.toList1D();
            if (!njobsList.get(j).isInfinite()) {
                njobsSum += njobsList.get(j);
            }
        }
        this.sn.nclosedjobs = njobsSum;
        this.sn.njobs = njobs.transpose();

    }

    public boolean[] refreshRates(List<Integer> statSet, List<Integer> classSet) {
        boolean hasRateChanged = false;
        boolean hasSCVChanged = false;
        int M = this.stations.size();
        int K = this.jobClasses.size();
        Matrix rates = null;
        Matrix scv = null;
        Matrix rates_orig = null;
        Matrix scv_orig = null;

        if (statSet == null && classSet == null) {
            statSet = new ArrayList<Integer>();
            for (int i = 0; i < M; i++)
                statSet.add(i);

            classSet = new ArrayList<Integer>();
            for (int i = 0; i < K; i++)
                classSet.add(i);

            rates = new Matrix(M, K, M * K);
            scv = new Matrix(M, K, M * K);
            scv.fill(Double.NaN);
            hasRateChanged = true;
            hasSCVChanged = true;
        } else {
            if (statSet == null) {
                statSet = new ArrayList<Integer>();
                for (int i = 0; i < M; i++)
                    statSet.add(i);
            }

            if (classSet == null) {
                classSet = new ArrayList<Integer>();
                for (int i = 0; i < K; i++)
                    classSet.add(i);
            }

            rates = this.sn.rates.clone();
            scv = this.sn.scv.clone();
            rates_orig = this.sn.rates.clone();
            scv_orig = this.sn.scv.clone();
        }
        boolean hasOpenClasses = this.hasOpenClasses();
        int sourceIdx = getIndexSourceNode();

        for (Integer i : statSet) {
            Station station = stations.get(i);
            for (Integer r : classSet) {
                if (station.getServer() instanceof ServiceTunnel) {
                    if (station instanceof Source) {
                        if (!((Source) station).containsJobClass(this.jobClasses.get(r))) {
                            rates.set(i, r, Double.NaN);
                            scv.set(i, r, Double.NaN);
                        } else {
                            Distribution distr = ((Source) station).getArrivalDistribution(this.jobClasses.get(r));
                            rates.set(i, r, distr.getRate());
                            scv.set(i, r, distr.getSCV());
                        }
                    } else if (station instanceof Join) {
                        rates.set(i, r, Double.POSITIVE_INFINITY);
                        scv.set(i, r, 0.0);
                    }
                    //Line 55-57 is ignored since no support of Place in JLINE
                } else {
                    if (!hasOpenClasses || i != sourceIdx) {
                        if (!station.getServer().containsJobClass(this.jobClasses.get(r))) {
                            rates.set(i, r, Double.NaN);
                            scv.set(i, r, Double.NaN);
                        } else {
                            Distribution distr = station.getServer().getServiceDistribution(this.jobClasses.get(r));
                            rates.set(i, r, distr.getRate());
                            scv.set(i, r, distr.getSCV());
                        }
                    }
                }
            }
        }

        if (!hasRateChanged) {
            Matrix tmp = rates.sub(1, rates_orig);
            tmp.abs();
            if (tmp.elementSum() > 0)
                hasRateChanged = true;
        }

        if (!hasSCVChanged) {
            Matrix tmp = scv.sub(1, scv_orig);
            tmp.abs();
            if (tmp.elementSum() > 0)
                hasSCVChanged = true;
        }

        if (hasRateChanged) {
            this.sn.rates = rates;
        }

        if (hasSCVChanged) {
            this.sn.scv = scv;
        }

        return new boolean[]{hasRateChanged, hasSCVChanged};
    }

    private void refreshServiceTypes(List<Integer> statSet, List<Integer> classSet) {
        int M = this.stations.size();
        int K = this.jobClasses.size();
        Map<Station, Map<JobClass, ProcessType>> proctype = null;

        if (statSet == null && classSet == null) {
            statSet = new ArrayList<Integer>();
            for (int i = 0; i < M; i++)
                statSet.add(i);

            classSet = new ArrayList<Integer>();
            for (int i = 0; i < K; i++)
                classSet.add(i);

            proctype = new HashMap<Station, Map<JobClass, ProcessType>>();
        } else if (statSet == null || classSet == null) {
            try {
                throw new Exception("refreshServiceTypes requires either both null or not null parameters");
            } catch (Exception e) {
                e.printStackTrace();
            }
        } else {
            proctype = this.sn.proctype;
        }
        boolean hasOpenClass = this.hasOpenClasses();
        int sourceIdx = getIndexSourceNode();

        for (Integer i : statSet) {
            Station station = this.stations.get(i);
            Map<JobClass, ProcessType> map = new HashMap<JobClass, ProcessType>();
            for (Integer r : classSet) {
                JobClass jobclass = this.jobClasses.get(r);
                if (station.getServer() instanceof ServiceTunnel) {
                    if (station instanceof Source) {
                        Distribution distr = ((Source) station).getArrivalDistribution(jobclass);
                        map.put(jobclass, getProcessType(distr));
                    } else if (station instanceof Join) {
                        map.put(jobclass, ProcessType.IMMEDIATE);
                    }
                } else {
                    if (!hasOpenClass || i != sourceIdx) {
                        if (!station.getServer().containsJobClass(jobclass)) {
                            map.put(jobclass, ProcessType.DISABLED);
                        } else {
                            Distribution distr = station.getServer().getServiceDistribution(jobclass);
                            map.put(jobclass, getProcessType(distr));
                        }
                    }
                }
            }
            proctype.put(station, map);
        }

        if (this.sn != null)
            this.sn.proctype = proctype;
    }

    public ProcessType getProcessType(Distribution distr) {
        if (distr instanceof Erlang) {
            return ProcessType.ERLANG;
        } else if (distr instanceof Exp) {
            return ProcessType.EXP;
        } else if (distr instanceof HyperExp) {
            return ProcessType.HYPEREXP;
        } else if (distr instanceof APH) {
            return ProcessType.APH;
        } else if (distr instanceof Det) {
            return ProcessType.DET;
        } else if (distr instanceof Coxian) {
            return ProcessType.COXIAN;
        } else if (distr instanceof Poisson) {
            return ProcessType.POISSON;
        } else if (distr instanceof Replayer) {
            return ProcessType.REPLAYER;
        } else if (distr instanceof Binomial) {
            return ProcessType.BINOMIAL;
        } else if (distr instanceof MAP) {
            return ProcessType.MAP;
        } else if (distr instanceof Immediate) {
            return ProcessType.IMMEDIATE;
        } else {
            return ProcessType.DISABLED;
        }
    }

    @SuppressWarnings("unchecked")
    public void refreshServicePhases(List<Integer> statSet, List<Integer> classSet) {
        int M = this.stations.size();
        int K = this.jobClasses.size();
        Map<Station, Map<JobClass, Matrix>> mu;
        Map<Station, Map<JobClass, Matrix>> phi;
        Matrix phases;

        if (statSet != null && classSet != null && this.sn.mu != null && this.sn.phi != null && this.sn.phases != null) {
            mu = this.sn.mu;
            phi = this.sn.phi;
            phases = this.sn.phases;
        } else {
            mu = new HashMap<Station, Map<JobClass, Matrix>>();
            phi = new HashMap<Station, Map<JobClass, Matrix>>();
            phases = new Matrix(stations.size(), jobClasses.size(), stations.size() * jobClasses.size());
            for (Station station : this.stations) {
                mu.put(station, new HashMap<JobClass, Matrix>());
                phi.put(station, new HashMap<JobClass, Matrix>());
            }
        }

        if (statSet == null) {
            statSet = new ArrayList<Integer>();
            for (int i = 0; i < M; i++)
                statSet.add(i);
        }
        if (classSet == null) {
            classSet = new ArrayList<Integer>();
            for (int i = 0; i < K; i++)
                classSet.add(i);
        }
        int sourceIdx = this.getIndexSourceStation();

        for (Integer i : statSet) {
            Station station = stations.get(i);
            Map<JobClass, Matrix> mu_i = null;
            Map<JobClass, Matrix> phi_i = null;
            if (i == sourceIdx) {
                List<Object> res = station.getSourceRates();
                mu_i = (Map<JobClass, Matrix>) res.get(1);
                phi_i = (Map<JobClass, Matrix>) res.get(2);
            } else {
                //Line 56 - 63 is ignored since fork is not station
                if (station instanceof Join) {
                    mu_i = new HashMap<JobClass, Matrix>();
                    phi_i = new HashMap<JobClass, Matrix>();
                    for (Integer r : classSet) {
                        Matrix mu_i_val = new Matrix(1, 1, 1);
                        Matrix phi_i_val = new Matrix(1, 1, 1);
                        mu_i_val.set(0, 0, Double.NaN);
                        phi_i_val.set(0, 0, Double.NaN);
                        mu_i.put(this.jobClasses.get(r), mu_i_val);
                        phi_i.put(this.jobClasses.get(r), phi_i_val);
                    }
                } else {
                    List<Object> res = station.getServiceRates();
                    mu_i = (Map<JobClass, Matrix>) res.get(1);
                    phi_i = (Map<JobClass, Matrix>) res.get(2);
                }
            }

            mu.put(station, mu_i);
            phi.put(station, phi_i);
            for (Integer r : classSet) {
                double[] mu_val = mu_i.get(this.jobClasses.get(r)).getNonZeroValues();
                boolean flag = true;
                for (int idx = 0; idx < mu_val.length; idx++)
                    flag = flag && Double.isNaN(mu_val[idx]);

                if (!flag)
                    phases.set(i, r, mu_val.length);
            }
        }

        if (this.sn != null) {
            this.sn.mu = mu;
            this.sn.phi = phi;
            this.sn.phases = phases;
            this.sn.phasessz = new Matrix(stations.size(), jobClasses.size(), stations.size() * jobClasses.size());
            this.sn.phaseshift = new Matrix(0, 0);
            for (int i = 0; i < stations.size(); i++) {
                for (int j = 0; j < jobClasses.size(); j++) {
                    this.sn.phasessz.set(i, j, Math.max(1.0, phases.get(i, j)));
                }
            }
            Matrix.concatColumns(new Matrix(this.sn.phases.getNumRows(), 1), this.sn.phasessz.cumsumViaRow(), this.sn.phaseshift);
        }
        refreshServiceRepresentations();
    }

    public void refreshLST(List<Integer> statSet, List<Integer> classSet) {
        int M = this.stations.size();
        int K = this.jobClasses.size();
        Map<Station, Map<JobClass, SerializableFunction<Double, Double>>> lst;

        if (statSet == null) {
            statSet = new ArrayList<Integer>();
            for (int i = 0; i < M; i++)
                statSet.add(i);
        }
        if (classSet == null) {
            classSet = new ArrayList<Integer>();
            for (int i = 0; i < K; i++)
                classSet.add(i);
        }

        if (this.sn.lst != null) {
            lst = this.sn.lst;
        } else {
            lst = new HashMap<Station, Map<JobClass, SerializableFunction<Double, Double>>>();
            for (Station station : stations) {
                lst.put(station, new HashMap<JobClass, SerializableFunction<Double, Double>>());
            }
        }
        int sourceIdx = this.getIndexSourceNode();

        for (Integer i : statSet) {
            Station station = this.stations.get(i);
            Map<JobClass, SerializableFunction<Double, Double>> map = new HashMap<JobClass, SerializableFunction<Double, Double>>();
            for (Integer r : classSet) {
                JobClass jobclass = this.jobClasses.get(r);
                if (i == sourceIdx) {
                    Distribution distr = ((Source) station).getArrivalDistribution(jobclass);
                    if (distr instanceof Disabled)
                        map.put(jobclass, null);
                    else
                        map.put(jobclass, e -> distr.evalLST(e));
                } else {
                    //line 45-46 is ignored since Fork is not station
                    if (station instanceof Join)
                        map.put(jobclass, null);
                    else
                        map.put(jobclass, e -> station.getServer().getServiceDistribution(jobclass).evalLST(e));
                }
            }
            lst.put(station, map);
        }

        if (this.sn != null)
            this.sn.lst = lst;
    }

    @SuppressWarnings("unchecked")
    //Current this function not implement the return value
    public void refreshServiceRepresentations() {
        int M = this.stations.size();
        int K = this.jobClasses.size();
        Map<Station, Map<JobClass, Map<Integer, Matrix>>> ph = new HashMap<Station, Map<JobClass, Map<Integer, Matrix>>>();
        for (int i = 0; i < M; i++)
            ph.put(this.stations.get(i), new HashMap<JobClass, Map<Integer, Matrix>>());
        Matrix phases = new Matrix(M, K, M * K);
        int sourceIdx = this.getIndexSourceStation();

        for (int i = 0; i < M; i++) {
            Station station = this.stations.get(i);
            Map<JobClass, Map<Integer, Matrix>> ph_i = new HashMap<JobClass, Map<Integer, Matrix>>();
            if (i == sourceIdx) {
                ph_i = (Map<JobClass, Map<Integer, Matrix>>) station.getSourceRates().get(0);
            } else {
                if (station instanceof Join) {
                    Coxian coxian = new Coxian(new ArrayList<Double>(Collections.singletonList(Double.NaN)), new ArrayList<Double>(Collections.singletonList(Double.NaN)));
                    for (JobClass jobclass : this.jobClasses)
                        ph_i.put(jobclass, coxian.getRepres());
                } else {
                    ph_i = (Map<JobClass, Map<Integer, Matrix>>) station.getServiceRates().get(0);
                }
            }
            ph.put(station, ph_i);

            for (int r = 0; r < K; r++) {
                Map<Integer, Matrix> ph_i_r = ph_i.get(this.jobClasses.get(r));
                if (ph_i_r == null)
                    phases.set(i, r, 1.0);
                else if (!(ph_i_r.get(0).hasNaN() || ph_i_r.get(1).hasNaN()))
                    phases.set(i, r, ph_i_r.get(0).getNumCols());
                //Other situation set to 0 (The matrix initial value is 0)
            }
        }

        if (this.sn != null) {
            Map<Station, Map<JobClass, Matrix>> pie = new HashMap<Station, Map<JobClass, Matrix>>();
            for (int i = 0; i < M; i++) {
                Station station = this.stations.get(i);
                Map<JobClass, Matrix> pie_i = new HashMap<JobClass, Matrix>();
                for (int r = 0; r < K; r++) {
                    JobClass jobclass = this.jobClasses.get(r);
                    Map<Integer, Matrix> map_ir = ph.get(station).get(jobclass);
                    if (map_ir != null) {
                        pie_i.put(jobclass, map_pie(map_ir.get(0), map_ir.get(1)));
                    } else {
                        Matrix tmp = new Matrix(1, 1, 1);
                        tmp.set(0, 0, Double.NaN);
                        pie_i.put(jobclass, tmp);
                    }
                }
                pie.put(station, pie_i);
            }

            this.sn.proc = ph;
            this.sn.pie = pie;
            this.sn.phases = phases;
            this.sn.phasessz = new Matrix(M, K, M * K);
            this.sn.phaseshift = new Matrix(0, 0);
            for (int i = 0; i < stations.size(); i++) {
                for (int j = 0; j < jobClasses.size(); j++) {
                    this.sn.phasessz.set(i, j, Math.max(1.0, phases.get(i, j)));
                }
            }
            //self.sn.phasessz(self.sn.nodeToStation(self.sn.nodetype == NodeType.Join),:)=phases(self.sn.nodeToStation(self.sn.nodetype == NodeType.Join),:);
            //Not tested, since current JLine does not support Join Node
            if (this.sn.nodeToStation != null) {
                for (int i = 0; i < this.sn.nodetypes.size(); i++) {
                    if (this.sn.nodetypes.get(i) != NodeType.Join) {
                        continue;
                    }
                    for (int j = 0; j < this.sn.phases.getNumCols(); j++) {
                        this.sn.phasessz.set((int) this.sn.nodeToStation.get(i), j, this.sn.phases.get((int) this.sn.nodeToStation.get(i), j));
                    }
                }
            }
            Matrix.concatColumns(new Matrix(this.sn.phases.getNumRows(), 1), this.sn.phasessz.cumsumViaRow(), this.sn.phaseshift);
        }
    }

    public void refreshScheduling() {
        int M = this.stations.size();
        int K = this.jobClasses.size();
        Map<Station, SchedStrategy> sched = getStationScheduling();
        Matrix schedparam = new Matrix(M, K, M * K);
        int sourceIdx = this.getIndexSourceNode();

        for (int i = 0; i < M; i++) {
            Station station = this.stations.get(i);
            if (sourceIdx == -1 || i != sourceIdx) {
                if (!(station.getServer() instanceof ServiceTunnel)) {
                    Queue queue = (Queue) station;
                    boolean NaNFlag = false;
                    for (int r = 0; r < K; r++) {
                        double val = queue.getSchedStrategyPar(this.jobClasses.get(r));
                        if (Double.isNaN(val)) {
                            NaNFlag = true;
                            for (int idx = 0; idx < r; idx++)
                                schedparam.remove(i, idx);
                            break;
                        } else {
                            schedparam.set(i, r, val);
                        }
                    }

                    if (NaNFlag) {
                        //Current JLine not support SEPT and LEPT schedule strategy.
                        continue;
                    }
                }
            }
        }

        if (this.sn != null) {
            this.sn.sched = sched;
            this.sn.schedparam = schedparam;
            for (int i = 0; i < M; i++) {
                this.sn.sched.put(sn.stations.get(i), sched.getOrDefault(this.stations.get(i), null));
            }
            //No schedid in JLine
        }
    }

    public Map<Station, SchedStrategy> getStationScheduling() {
        Map<Station, SchedStrategy> res = new HashMap<Station, SchedStrategy>();
        for (Station station : this.stations) {
            if (station.getNumberOfServers() == Integer.MAX_VALUE) {
                res.put(station, SchedStrategy.INF);
            } else {
                if (station instanceof Source) {
                    res.put(station, SchedStrategy.EXT);
                } else {
                    HasSchedStrategy s = (HasSchedStrategy) station;
                    res.put(station, s.getSchedStrategy());
                }
            }
        }
        return res;
    }

    public void refreshChains(boolean propagate) {
        propagate = true;
        refreshRoutingMatrix(this.sn.rates);
        Matrix rt = this.sn.rt;
        Matrix rtnodes = this.sn.rtnodes;

        Matrix stateful = this.sn.isstateful.find();
        int K = this.sn.nclasses;

        if (this.csMatrix == null) {
            Matrix csmask = new Matrix(K, K);
            for (int r = 0; r < K; r++) {
                for (int s = 0; s < K; s++) {
                    for (int isf = 0; isf < stateful.getNumRows(); isf++) {
                        for (int jsf = 0; jsf < stateful.getNumRows(); jsf++) {
                            if (rt.get(isf * K + r, jsf * K + s) > 0)
                                csmask.set(r, s, 1.0);
                        }
                    }
                }
            }

            for (int isf = 0; isf < stateful.getNumRows(); isf++) {
                int ind = (int) this.sn.statefulToNode.get(0, isf);
                boolean isCS = (this.sn.nodetypes.get(ind) == NodeType.Cache) || (this.sn.nodetypes.get(ind) == NodeType.ClassSwitch);
                for (int r = 0; r < K; r++) {
                    csmask.set(r, r, 1.0);
                    for (int s = 0; s < K; s++) {
                        if (r != s) {
                            if (isCS) {
                                ClassSwitcher classSwitcher = (ClassSwitcher) this.nodes.get(ind).getServer();
                                if (classSwitcher.applyCsFun(r, s) > 0)
                                    csmask.set(r, s, 1.0);
                            }
                        }
                    }
                }
            }
            this.sn.csmask = csmask;
        } else {
            this.sn.csmask = this.csMatrix;
        }

        if (((sn.refclass != null) && (!sn.refclass.isEmpty())) && (sn.refclass.length() < sn.nchains)) {
            sn.refclass.expandMatrix(1, sn.nchains, sn.nchains);
        }

        if (propagate) {
            //Compute visits
            this.sn = SN.snRefreshVisits(this.sn, this.sn.chains, rt, rtnodes);

            //Call dependent capacity refresh
            refreshCapacity();
        }
    }

    public void refreshRoutingMatrix(Matrix rates) {
        if (rates == null)
            throw new RuntimeException("refreshRoutingMatrix cannot retrieve station rates, pass them as an input parameters.");

        int M = this.getNumberOfNodes();
        int K = this.getNumberOfClasses();
        Matrix arvRates = new Matrix(1, K, K);
        List<Integer> stateful = this.getIndexStatefulNodes();
        int indSourceStation = this.getIndexSourceStation();
        for (Integer i : this.getIndexOpenClasses()) {
            arvRates.set(0, i, rates.get(indSourceStation, i));
        }

        routingMatrixReturn res = getRoutingMatrix(arvRates, 4);
        Matrix rt = res.rt;
        Matrix rtnodes = res.rtnodes;
        Matrix linksmat = res.linksmat;
        Matrix chains = res.chains;

        if (this.doChecks) {
            outerloop:
            for (JobClass jobclass : this.jobClasses) {
                for (Map<JobClass, RoutingStrategy> nodeRoutingMap : this.sn.routing.values()) {
                    if (nodeRoutingMap.get(jobclass) != RoutingStrategy.DISABLED)
                        continue outerloop;
                }
                throw new RuntimeException("Routing strategy is unspecified at all nodes");
            }
        }

        boolean isStateDep = (Matrix.extractColumn(this.sn.isstatedep, 2, null).getNonZeroLength() > 0);
        Map<Integer, Map<Integer, SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Double>>> rtnodefuncell =
                new HashMap<Integer, Map<Integer, SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Double>>>();

        if (isStateDep) {
            for (int ind = 0; ind < M; ind++) {
                final int ind_final = ind;
                for (int jnd = 0; jnd < M; jnd++) {
                    final int jnd_final = jnd;
                    for (int r = 0; r < K; r++) {
                        final int r_final = r;
                        for (int s = 0; s < K; s++) {
                            final int s_final = s;
                            Map<Integer, SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Double>> map =
                                    rtnodefuncell.getOrDefault(ind * K + r, new HashMap<Integer, SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Double>>());
                            if (this.sn.isstatedep.get(ind, 2) > 0) {
                                switch (this.sn.routing.get(this.nodes.get(ind)).get(this.jobClasses.get(r))) {
                                    case RROBIN:
                                        map.put(jnd * K + s, (pair) -> sub_rr_wrr(ind_final, jnd_final, r_final, s_final, linksmat, pair.getLeft(), pair.getRight()));
                                    case WRROBIN:
                                        map.put(jnd * K + s, (pair) -> sub_rr_wrr(ind_final, jnd_final, r_final, s_final, linksmat, pair.getLeft(), pair.getRight()));
                                    case JSQ:
                                        map.put(jnd * K + s, (pair) -> sub_jsq(ind_final, jnd_final, r_final, s_final, linksmat, pair.getLeft(), pair.getRight()));
                                    default:
                                        map.put(jnd * K + s, (pair) -> rtnodes.get(ind_final * K + r_final, jnd_final * K + s_final));
                                }
                            } else {
                                map.put(jnd * K + s, (pair) -> rtnodes.get(ind_final * K + r_final, jnd_final * K + s_final));
                            }
                            rtnodefuncell.put(ind * K + r, map);
                        }
                    }
                }
            }
        }

        /* we now generate the node routing matrix for the given state and then
         * lump the states for non-stateful nodes so that run gives the routing
         *  table for stateful nodes only */
        List<Integer> statefulNodeClasses = new ArrayList<Integer>(); //Not using JLineMatrix for performance consideration
        for (int i = 0; i < stateful.size(); i++) {
            for (int j = 0; j < K; j++) {
                statefulNodeClasses.add(stateful.get(i) * K + j);
            }
        }

        SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Matrix> rtfun = null;
        if (isStateDep) {
            rtfun = (pair) -> {
                Matrix cellfunnodes = new Matrix(M * K, M * K);
                for (int ind = 0; ind < M; ind++) {
                    for (int jnd = 0; jnd < M; jnd++) {
                        for (int r = 0; r < K; r++) {
                            for (int s = 0; s < K; s++) {
                                int row = ind * K + r;
                                int col = jnd * K + s;
                                double val = rtnodefuncell.get(row).get(col).apply(pair);
                                if (val != 0)
                                    cellfunnodes.set(row, col, val);
                            }
                        }
                    }
                }
                return DTMC.dtmc_stochcomp(cellfunnodes, statefulNodeClasses);
            };
        } else {
            rtfun = ((pair) -> DTMC.dtmc_stochcomp(rtnodes, statefulNodeClasses));
        }

        int nchains = chains.getNumRows();
        Map<Integer, Matrix> inchain = new HashMap<Integer, Matrix>();
        for (int c = 0; c < nchains; c++) {
            Matrix chains_c = new Matrix(1, chains.getNumCols());
            Matrix.extract(chains, c, c + 1, 0, chains.getNumCols(), chains_c, 0, 0);
            Matrix chains_c_t = chains_c.find().transpose();
            inchain.put(c, chains_c_t);
        }

        this.sn.rt = rt;
        this.sn.rtnodes = rtnodes;
        this.sn.rtfun = rtfun;
        this.sn.chains = chains;
        this.sn.nchains = nchains;
        this.sn.inchain = inchain;
        for (int c = 0; c < nchains; c++) {
            Matrix inchain_c = inchain.get(c);
            double val = this.sn.refstat.get((int) inchain_c.get(0, 0), 0);
            for (int col = 1; col < inchain_c.getNumCols(); col++) {
                if (val != this.sn.refstat.get((int) inchain_c.get(0, col), 0))
                    throw new RuntimeException("Classes within chain have different reference station");
            }
        }
    }

    public void printRoutingMatrix() {
        this.getStruct(false);
        for (int ind = 0; ind < sn.nnodes; ind++) {
            for (int jnd = 0; jnd < sn.nnodes; jnd++) {
                for (int r = 0; r < sn.nclasses; r++) {
                    for (int s = 0; s < sn.nclasses; s++) {
                        double pr = sn.rtnodes.get(ind * sn.nclasses + r, jnd * sn.nclasses + s);
                        if (pr > 0)
                            System.out.println(sn.nodenames.get(ind) + " [" + sn.classnames.get(r) + "] => " + sn.nodenames.get(jnd) + " [" + sn.classnames.get(s) + "] : Pr=" + pr);
                    }
                }
            }
        }
    }

    public routingMatrixReturn getRoutingMatrix(Matrix arvRates, int returnVal) {

        int idxSource, idxSink, I, K;
        List<Integer> idxOpenClasses;
        boolean hasOpen;
        Matrix conn, NK;

        if (this.hasStruct) {
            idxSource = this.sn.nodetypes.indexOf(NodeType.Source);
            idxSink = this.sn.nodetypes.indexOf(NodeType.Sink);
            idxOpenClasses = new ArrayList<Integer>();
            for (int col = 0; col < this.sn.njobs.getNumCols(); col++) {
                if (Double.isInfinite((this.sn.njobs.get(0, col))))
                    idxOpenClasses.add(col);
            }
            hasOpen = !idxOpenClasses.isEmpty();
            if ((arvRates == null) || arvRates.isEmpty()) {
                arvRates = new Matrix(1, idxOpenClasses.size(), idxOpenClasses.size());
                for (int i = 0; i < idxOpenClasses.size(); i++) {
                    arvRates.set(0, i, this.sn.rates.get(idxSource, idxOpenClasses.get(i)));
                }
            }
            conn = this.sn.connmatrix;
            I = this.sn.nnodes;
            K = this.sn.nclasses;
            NK = this.sn.njobs;
        } else {
            idxSource = this.getIndexSourceNode();
            idxSink = this.getIndexSinkNode();
            idxOpenClasses = this.getIndexOpenClasses();
            conn = this.getConnectionMatrix();
            hasOpen = this.hasOpenClasses();
            I = this.getNumberOfNodes();
            K = this.getNumberOfClasses();
            NK = this.getNumberOfJobs().transpose();

            if (this.sn == null)
                this.sn = new NetworkStruct();
            this.sn.connmatrix = conn;
            sn.routing = new HashMap<Node, Map<JobClass, RoutingStrategy>>();
            for (Node node : this.nodes) {
                Map<JobClass, RoutingStrategy> map = new HashMap<JobClass, RoutingStrategy>();
                for (JobClass jobclass : this.jobClasses) {
                    map.put(jobclass, getRoutingStrategyFromNodeAndClassPair(node, jobclass));
                }
                sn.routing.put(node, map);
            }
        }

        Matrix rtnodes = new Matrix(I * K, I * K, 0);
        Matrix chains = null;

        // The first loop considers the class at which a job enters the
        for (int i = 0; i < I; i++) {
            Node node = this.nodes.get(i);
            if (node.getOutput() instanceof Forker) {
                //Line 1294-1314 not test
                for (int j = 0; j < I; j++) {
                    for (int k = 0; k < K; k++) {
                        if (conn.get(i, j) > 0) {
                            JobClass jobclass = this.jobClasses.get(k);
                            List<OutputStrategy> outputStrategy_k = node.getOutput().getOutputStrategyByClass(jobclass);
                            if (outputStrategy_k.size() > 0)
                                rtnodes.set(i * K + k, j * K + k, 1.0 / outputStrategy_k.size());
                            else
                                rtnodes.set(i * K + k, j * K + k, 1.0);
                            if (this.sn.routing.get(node).get(jobclass) == RoutingStrategy.PROB) {
                                //Check the number of outgoing links
                                int sum = (int) conn.sumRows(i);
                                //Fork must have all the output strategy towards all outgoing links.
                                if (outputStrategy_k.size() != sum)
                                    throw new RuntimeException("Fork must have all the output strategy towards all outgoing links.");
                                //Fork must have 1.0 routing probability towards all outgoing links.
                                for (OutputStrategy ops : outputStrategy_k) {
                                    if (ops.getProbability() != 1.0)
                                        throw new RuntimeException("Fork must have 1.0 routing probability towards all outgoing links.");
                                }
                            }
                        }
                    }
                }
            } else {
                boolean isSink_i = (i == idxSink);
                boolean isSource_i = (i == idxSource);
                for (int k = 0; k < K; k++) {
                    JobClass jobclass = this.jobClasses.get(k);
                    List<OutputStrategy> outputStrategy_k = node.getOutput().getOutputStrategyByClass(jobclass);
                    switch (this.sn.routing.get(node).get(jobclass)) {
                        case PROB:
                            if (Double.isInfinite((NK.get(0, k))) || !isSink_i) {
                                for (OutputStrategy ops : outputStrategy_k) {
                                    int j = ops.getDestination().getNodeIdx();
                                    rtnodes.set(i * K + k, j * K + k, ops.getProbability());
                                }
                            }
                            break;
                        //Not tested the following situation
                        case DISABLED:
                            double sum = conn.sumRows(i);
                            for (int j = 0; j < I; j++) {
                                if (conn.get(i, j) > 0)
                                    rtnodes.set(i * K + k, j * K + k, 1.0 / sum);
                            }
                            break;
                        case RAND:
                        case RROBIN:
                        case WRROBIN:
                        case JSQ:
                            if (Double.isInfinite((NK.get(0, k)))) {
                                sum = conn.sumRows(i);
                                for (int j = 0; j < I; j++) {
                                    if (conn.get(i, j) > 0)
                                        rtnodes.set(i * K + k, j * K + k, 1.0 / sum);
                                }
                            } else if (!isSource_i && !isSink_i) {
                                Matrix connectionClosed = conn.clone();
                                if (idxSink >= 0) { // this is empty set in MATLAB
                                    if (connectionClosed.get(i, idxSink) > 0)
                                        connectionClosed.remove(i, idxSink);
                                }
                                sum = connectionClosed.sumRows(i);
                                for (int j = 0; j < I; j++) {
                                    if (connectionClosed.get(i, j) > 0)
                                        rtnodes.set(i * K + k, j * K + k, 1.0 / sum);
                                }
                            }
                            break;
                        default:
                            for (int j = 0; j < I; j++) {
                                if (conn.get(i, j) > 0)
                                    rtnodes.set(i * K + k, j * K + k, GlobalConstants.Zero);
                            }
                    }

                }
            }
        }

        // The second loop corrects the first one at nodes that change the class of the job in the service section.
        for (int i = 0; i < I; i++) {
            Node node = this.nodes.get(i);
            if (node.getServer() instanceof StatelessClassSwitcher) {
                Matrix Pi = new Matrix(K, rtnodes.getNumCols(), 0);
                Matrix.extract(rtnodes, i * K, (i + 1) * K, 0, rtnodes.getNumCols(), Pi, 0, 0);
                Matrix Pcs = new Matrix(K, K);
                ClassSwitcher classSwitcher = (ClassSwitcher) node.getServer();
                for (int r = 0; r < K; r++) {
                    for (int s = 0; s < K; s++) {
                        Pcs.set(r, s, classSwitcher.applyCsFun(r, s));
                    }
                }

                for (int jnd = 0; jnd < I; jnd++) {
                    Matrix Pij = new Matrix(K, K, K * K);
                    Matrix.extract(Pi, 0, K, K * jnd, K * (jnd + 1), Pij, 0, 0);
                    //diag(Pij)'
                    Matrix diagPij = new Matrix(1, K, K);
                    Matrix.extractDiag(Pij, diagPij);
                    //repmat(diag(Pij)', K, 1)
                    Matrix repmatPij = diagPij.repmat(K, 1);
                    //rtnodes(((ind-1)*K+1) : ((ind-1)*K+K),(jnd-1)*K+(1:K)) = Pcs.*repmat(diag(Pij)',K,1);
                    for (int row = 0; row < K; row++) {
                        for (int col = 0; col < K; col++) {
                            double val = repmatPij.get(row, col) * Pcs.get(row, col);
                            if (val != 0)
                                rtnodes.set(i * K + row, jnd * K + col, val);
                            else
                                rtnodes.remove(i * K + row, jnd * K + col);
                        }
                    }
                }
            } else if (node.getServer() instanceof StatefulClassSwitcher) {
                Matrix Pi = new Matrix(K, rtnodes.getNumCols(), 0);
                Matrix.extract(rtnodes, i * K, (i + 1) * K, 0, rtnodes.getNumCols(), Pi, 0, 0);
                Pi = new Matrix(Pi);
                Matrix Pcs = new Matrix(K, K);
                for (int r = 0; r < K; r++) {
                    for (int s = 0; s < K; s++) {
                        Pcs.set(r, s, ((StatefulClassSwitcher) node.getServer()).applyCsFun(r, s));
                    }
                }
                for (int row = i * K; row < (i + 1) * K; row++) {
                    for (int col = 0; col < rtnodes.getNumCols(); col++) {
                        rtnodes.set(row, col, 0);
                    }
                }
                if (node.getServer() instanceof CacheClassSwitcher) {
                    CacheClassSwitcher cacheClassSwitcher = (CacheClassSwitcher) node.getServer();
                    for (int r = 0; r < K; r++) {
                        boolean flagHit = false, flagMiss = false;
                        for (int j = 0; j < cacheClassSwitcher.hitClass.length(); j++) {
                            if (cacheClassSwitcher.hitClass.get(j) == r) {
                                flagHit = true;
                                break;
                            }
                        }
                        for (int j = 0; j < cacheClassSwitcher.missClass.length(); j++) {
                            if (cacheClassSwitcher.missClass.get(j) == r) {
                                flagMiss = true;
                                break;
                            }
                        }
                        if (!flagHit && !flagMiss) {
                            double rowsum = Matrix.extractRows(Pcs, r, r + 1, null).elementSum();
                            for (int j = 0; j < Pcs.getNumCols(); j++) {
                                Pcs.set(r, j, Pcs.get(r, j) / rowsum);
                            }
                        }
                    }
                    for (int r = 0; r < K; r++) {
                        boolean flagHit = false, flagMiss = false;
                        for (int j = 0; j < cacheClassSwitcher.hitClass.length(); j++) {
                            if (cacheClassSwitcher.hitClass.get(j) == r) {
                                flagHit = true;
                                break;
                            }
                        }
                        for (int j = 0; j < cacheClassSwitcher.missClass.length(); j++) {
                            if (cacheClassSwitcher.missClass.get(j) == r) {
                                flagMiss = true;
                                break;
                            }
                        }
                        if (!flagHit && !flagMiss) {
                            for (int jnd = 0; jnd < I; jnd++) {
                                for (int s = 0; s < K; s++) {
                                    if (i * K + r >= Pi.getNumRows()) {
                                        Matrix newPi = new Matrix(i * K + r + 1, Pi.getNumCols());
                                        for (int x = 0; x < Pi.getNumRows(); x++) {
                                            for (int y = 0; y < Pi.getNumCols(); y++) {
                                                newPi.set(x, y, Pi.get(x, y));
                                            }
                                        }
                                        Pi = newPi;
                                    }
                                    Pi.set(i * K + r, jnd * K + s, 0);
                                }
                            }
                        }
                    }
                    for (int r = 0; r < K; r++) {
                        if (cacheClassSwitcher.actualHitProb.length() > r && cacheClassSwitcher.hitClass.length() > r &&
                                cacheClassSwitcher.hitClass.get(r) != -1) {
                            double ph = cacheClassSwitcher.actualHitProb.get(r);
                            double pm = cacheClassSwitcher.actualMissProb.get(r);

                            int h = (int) cacheClassSwitcher.hitClass.get(r);
                            int m = (int) cacheClassSwitcher.missClass.get(r);

                            rtnodes.set(i * K + r, i * K + h, ph);
                            rtnodes.set(i * K + r, i * K + m, pm);
                        } else {
                            if (cacheClassSwitcher.hitClass.length() > r && cacheClassSwitcher.hitClass.get(r) != -1) {
                                int h = (int) cacheClassSwitcher.hitClass.get(r);
                                int m = (int) cacheClassSwitcher.missClass.get(r);
                                rtnodes.set(i * K + r, i * K + h, Double.NaN);
                                rtnodes.set(i * K + r, i * K + m, Double.NaN);
                            }
                        }
                    }
                    for (int jnd = 0; jnd < I; jnd++) {
                        Matrix Pij = new Matrix(K, K);
                        Matrix.extract(Pi, 0, K, jnd * K, (jnd + 1) * K, Pij, 0, 0);
                        for (int r = 0; r < K; r++) {
                            boolean flagHit = false, flagMiss = false;
                            for (int j = 0; j < cacheClassSwitcher.hitClass.length(); j++) {
                                if (cacheClassSwitcher.hitClass.get(j) == r) {
                                    flagHit = true;
                                    break;
                                }
                            }
                            for (int j = 0; j < cacheClassSwitcher.missClass.length(); j++) {
                                if (cacheClassSwitcher.missClass.get(j) == r) {
                                    flagMiss = true;
                                    break;
                                }
                            }
                            if (flagHit || flagMiss) {
                                for (int s = 0; s < K; s++) {
                                    rtnodes.set(i * K + r, jnd * K + s, Pcs.get(r, s) * Pij.get(s, s));
                                }
                            }
                        }
                    }
                }

            }
        }
        // ignore all chains containing a Pnodes column that sums to 0, since these are classes that cannot arrive to the node unless this column belongs to the source
        Matrix sumRtnodesCols = rtnodes.sumCols();
        Set<Integer> colsToIgnore = new HashSet<Integer>();
        for (int col = 0; col < sumRtnodesCols.getNumCols(); col++) {
            if (sumRtnodesCols.get(col) == 0)
                colsToIgnore.add(col);
        }
        if (hasOpen) {
            for (int i = idxSource * K; i < (idxSource + 1) * K; i++)
                colsToIgnore.remove(i);
        }

        /* We route back from the sink to the source. Since open classes
         * have an infinite population, if there is a class switch QN
         * with the following chains
         * Source -> (A or B) -> C -> Sink
         * Source -> D -> Sink
         * We can re-route class C into the source either as A or B or C.
         * We here re-route back as C and leave for the chain analyzer
         * to detect that C is in a chain with A and B and change this
         * part.
         */
        if (this.csMatrix == null) {
            Matrix param = rtnodes.add(1, rtnodes.transpose());
            Set<Set<Integer>> chainCandidates = Matrix.weaklyConnect(param, colsToIgnore);

            Matrix chainstmp = new Matrix(chainCandidates.size(), K);
            Iterator<Set<Integer>> it = chainCandidates.iterator();
            int tmax = 0;
            while (it.hasNext()) {
                Set<Integer> set = it.next();
                if (set.size() > 1) {
                    for (Integer num : set)
                        chainstmp.set(tmax, num % K, 1.0);
                    tmax++;
                }
            }
            chains = new Matrix(tmax, K);
            Matrix.extract(chainstmp, 0, tmax, 0, chainstmp.getNumCols(), chains, 0, 0);
        } else {
            Set<Set<Integer>> chainCandidates = Matrix.weaklyConnect(this.csMatrix, new HashSet<Integer>());
            chains = new Matrix(chainCandidates.size(), K);
            Iterator<Set<Integer>> it = chainCandidates.iterator();
            int tmax = 0;
            while (it.hasNext()) {
                Set<Integer> set = it.next();
                for (Integer num : set)
                    chains.set(tmax, num, 1.0);
                tmax++;
            }
        }
        this.sn.chains = chains;

        //Split chains block
        List<Integer> splitChains = new ArrayList<Integer>();
        Matrix sumCol = chains.sumCols();
        for (int i = 0; i < sumCol.getNumCols(); i++) {
            if (sumCol.get(i) > 1) {
                //rows = find(chains(:,col));
                List<Integer> rows = new ArrayList<Integer>();
                for (int j = 0; j < chains.getNumRows(); j++) {
                    if (chains.get(j, i) == 1)
                        rows.add(j);
                }

                if (rows.size() > 1) {
                    int row = rows.get(0);
                    for (int j = 1; j < rows.size(); j++) {
                        //chains(rows(1),:) = chains(row(1),:) | chains(r,:);
                        for (int k = 0; k < chains.getNumCols(); k++) {
                            if (chains.get(row, k) == 1 || chains.get(rows.get(j), k) == 1)
                                chains.set(row, k, 1);
                            else
                                chains.set(row, k, 0);
                        }
                        splitChains.add(j);
                    }
                }
            }
        }

        if (splitChains.size() > 0) {
            Matrix newChains = new Matrix(chains.getNumRows() - splitChains.size(), chains.getNumCols());
            for (int i = 0; i < chains.getNumRows(); i++) {
                if (!splitChains.contains(i)) {
                    for (int j = 0; j < chains.getNumCols(); j++)
                        newChains.set(i, j, chains.get(i, j));
                }
            }
            chains = newChains;
        }

        /* We now obtain the routing matrix P by ignoring the non-stateful
         * nodes and calculating by the stochastic complement method the
         * correct transition probabilities, that includes the effects
         * of the non-stateful nodes (e.g., ClassSwitch)
         */
        List<Integer> statefulNodes = this.getIndexStatefulNodes();
        List<Integer> statefulNodeClasses = new ArrayList<Integer>(); //Not using JLineMatrix for performance consideration
        for (int i = 0; i < statefulNodes.size(); i++) {
            for (int j = 0; j < K; j++) {
                statefulNodeClasses.add(statefulNodes.get(i) * K + j);
            }
        }

        /* this routes open classes back from the sink into the source
         * it will not work with non-renewal arrivals as it chooses in which open
         * class to reroute a job with probability depending on the arrival rates
         */
        if (hasOpen) {
            arvRates.removeNaN();
            for (int i = 0; i < idxOpenClasses.size(); i++) {
                //s_chain = find(chains(:,s));
                Matrix s_chain = new Matrix(chains.getNumRows(), 1, 0);
                Matrix.extract(chains, 0, chains.getNumRows(), idxOpenClasses.get(i), idxOpenClasses.get(i) + 1, s_chain, 0, 0);
                s_chain = s_chain.find();
                //others_in_chain = find(chains(s_chain,:));
                Matrix others_in_chain = new Matrix(s_chain.getNumRows(), chains.getNumCols(), 0);
                for (int row = 0; row < s_chain.getNumRows(); row++)
                    Matrix.extract(chains, (int) s_chain.get(row, 0), (int) s_chain.get(row, 0) + 1, 0, chains.getNumCols(), others_in_chain, row, 0);
                others_in_chain = others_in_chain.find();
                //arvRates(others_in_chain)/sum(arvRates(others_in_chain))
                Matrix arv_rates_others_in_chain = new Matrix(1, others_in_chain.getNumRows(), 0);
                for (int row = 0; row < others_in_chain.getNumRows(); row++)
                    arv_rates_others_in_chain.set(0, row, arvRates.get(0, (int) others_in_chain.get(row, 0)));
                arv_rates_others_in_chain.divide(arv_rates_others_in_chain.sumRows(0), arv_rates_others_in_chain, true);
                //repmat(arvRates(others_in_chain)/sum(arvRates(others_in_chain)),length(others_in_chain),1);
                Matrix rep_res = arv_rates_others_in_chain.repmat(others_in_chain.getNumRows(), 1);
                //rtnodes((idxSink-1)*K+others_in_chain,(idxSource-1)*K+others_in_chain) = rep_res
                for (int row1 = 0; row1 < others_in_chain.getNumRows(); row1++) {
                    for (int row2 = 0; row2 < others_in_chain.getNumRows(); row2++) {
                        rtnodes.set(idxSink * K + (int) others_in_chain.get(row1, 0), idxSource * K + (int) others_in_chain.get(row2, 0), rep_res.get(row1, row2));
                    }
                }
            }
        }

        /* Hide the nodes that are not stateful */
        Matrix rt = DTMC.dtmc_stochcomp(rtnodes, statefulNodeClasses);
        this.sn.rt = rt;

        /* Compute the optional outputs */
        Map<JobClass, Map<JobClass, Matrix>> rtNodesByClass = null;
        if (returnVal >= 5) {
            rtNodesByClass = new HashMap<JobClass, Map<JobClass, Matrix>>();
            for (int r = 0; r < K; r++) {
                Map<JobClass, Matrix> map = new HashMap<JobClass, Matrix>();
                for (int s = 0; s < K; s++) {
                    Matrix matrix = new Matrix(I, I, I * I);
                    for (int i = 0; i < I; i++) {
                        for (int j = 0; j < I; j++) {
                            matrix.set(i, j, rtnodes.get(i * K + s, j * K + r));
                        }
                    }
                    map.put(this.jobClasses.get(s), matrix);
                }
                rtNodesByClass.put(this.jobClasses.get(r), map);
            }
        }

        Map<Node, Map<Node, Matrix>> rtNodesByStation = null;
        if (returnVal >= 6) {
            rtNodesByStation = new HashMap<Node, Map<Node, Matrix>>();
            for (int i = 0; i < I; i++) {
                Map<Node, Matrix> map = new HashMap<Node, Matrix>();
                for (int j = 0; j < I; j++) {
                    Matrix matrix = new Matrix(K, K, K * K);
                    for (int r = 0; r < K; r++) {
                        for (int s = 0; s < K; s++) {
                            matrix.set(r, s, rtnodes.get(i * K + s, j * K + r));
                        }
                    }
                    map.put(this.nodes.get(j), matrix);
                }
                rtNodesByStation.put(this.nodes.get(i), map);
            }
        }

        //Return
        return new routingMatrixReturn(rt, rtnodes, this.sn.connmatrix, chains, rtNodesByClass, rtNodesByStation);
    }

    public void refreshCapacity() {
        int M = this.stations.size();
        int K = this.jobClasses.size();
        int C = this.sn.nchains;

        Matrix classcap = new Matrix(M, K);
        classcap.fill(Double.POSITIVE_INFINITY);
        Matrix chaincap = new Matrix(M, K);
        chaincap.fill(Double.POSITIVE_INFINITY);
        Matrix capacity = new Matrix(M, 1);
        //Something wrong with dropRule in LINE
        Map<Station, Map<JobClass, DropStrategy>> dropRule = new HashMap<Station, Map<JobClass, DropStrategy>>();
        for (Station station : this.stations) {
            Map<JobClass, DropStrategy> dropRule_station = new HashMap<JobClass, DropStrategy>();
            for (JobClass jobclass : this.jobClasses)
                dropRule_station.put(jobclass, DropStrategy.WaitingQueue);
            dropRule.put(station, dropRule_station);
        }

        Matrix njobs = this.sn.njobs;
        Matrix rates = this.sn.rates;

        for (int c = 0; c < C; c++) {
            Matrix inchain_c = this.sn.inchain.get(c);
            //chainCap = sum(njobs(inchain));
            double chainCap = 0;
            for (int idx = 0; idx < inchain_c.length(); idx++)
                chainCap += njobs.get(0, (int) inchain_c.get(0, idx));

            for (int idx = 0; idx < inchain_c.length(); idx++) {
                int r = (int) inchain_c.get(idx);
                JobClass jobclass = this.jobClasses.get(r);
                for (int i = 0; i < M; i++) {
                    Station station = this.stations.get(i);
                    if (!(station instanceof Source)) {
                        dropRule.get(station).put(jobclass, station.getDropRule(jobclass));
                    }
                    if (Double.isNaN(rates.get(i, r))) {
                        classcap.remove(i, r);
                        chaincap.remove(i, c);
                    } else {
                        classcap.set(i, r, chainCap);
                        chaincap.set(i, c, chainCap);
                        if (station.getClassCap(jobclass) >= 0)
                            classcap.set(i, r, Math.min(classcap.get(i, r), station.getClassCap(jobclass)));
                        if (station.getCap() >= 0)
                            classcap.set(i, r, Math.min(classcap.get(i, r), station.getCap()));
                    }
                }
            }
        }

        for (int i = 0; i < M; i++)
            capacity.set(i, 0, Math.min(chaincap.sumRows(i), classcap.sumRows(i)));

        this.sn.cap = capacity;
        this.sn.classcap = classcap;
        this.sn.droprule = dropRule;
    }

    public void refreshLocalVars() {
        int R = this.jobClasses.size();
        int I = this.nodes.size();
        Matrix nvars = new Matrix(I, 2 * R + 1);
        Map<Node, NodeParam> nodeparam = new HashMap<Node, NodeParam>();

        for (int i = 0; i < I; i++) {
            Node node = this.nodes.get(i);
            NodeParam param = new NodeParam();
            switch (this.sn.nodetypes.get(i)) {
                case Cache:
                    Cache cache = (Cache) node;
                    nvars.set(i, 2 * R, cache.getItemLevelCap().elementSum());
                    param.nitems = 0;
                    param.accost = cache.accessProb;
                    for (int r = 0; r < this.getNumberOfClasses(); r++) {
                        if (!cache.popularityGet(r).isDisabled()) {
                            param.nitems = (int) Maths.max(param.nitems, cache.popularityGet(r).getSupport().getRight());
                        }
                    }
                    param.itemcap = cache.getItemLevelCap();
                    param.pread = new HashMap<>();
                    for (int r = 0; r < this.getNumberOfClasses(); r++) {
                        if (cache.popularityGet(r).isDisabled()) {
                            param.pread.put(r, null);
                        } else {
                            List<Double> t = new ArrayList<>();
                            for (int j = 1; j <= param.nitems; j++) {
                                t.add((double) j);
                            }
                            param.pread.put(r, cache.popularityGet(r).evalPMF(t).toList1D());
                        }
                    }
                    param.rpolicy = cache.getReplacementPolicy();
                    param.hitclass = new Matrix(cache.getCacheServer().hitClass.getNumRows(), cache.getCacheServer().hitClass.getNumCols());
                    for (int r = 0; r < cache.getCacheServer().hitClass.getNumRows(); r++) {
                        for (int c = 0; c < cache.getCacheServer().hitClass.getNumCols(); c++) {
                            param.hitclass.set(r, c, Math.round(cache.getCacheServer().hitClass.get(r, c)));
                        }
                    }
                    param.missclass = new Matrix(cache.getCacheServer().missClass.getNumRows(), cache.getCacheServer().missClass.getNumCols());
                    for (int r = 0; r < cache.getCacheServer().missClass.getNumRows(); r++) {
                        for (int c = 0; c < cache.getCacheServer().missClass.getNumCols(); c++) {
                            param.missclass.set(r, c, Math.round(cache.getCacheServer().missClass.get(r, c)));
                        }
                    }
                    break;
                case Fork:
                    param.fanOut = ((Forker) node.getOutput()).taskPerLink;
                    break;
                case Join:
                    Joiner joiner = (Joiner) node.getInput();
                    param.joinStrategy = joiner.joinStrategy;
                    param.fanIn = joiner.joinRequired;
                    break;
                case Logger:
                    // TODOgetAvgSysRespT
                    throw new RuntimeException("Logger node is not supported in JLINE");
                case Source:
                    // TODO, refactor getServiceProcess? in MATLAB this is .arrivalProcess
                    for (int r = 0; r < R; r++) {
                        ServiceBinding serviceProcess = node.getServer().getServiceProcess(this.getClassByIndex(r));
                        if (serviceProcess != null) {
                            Distribution serviceDistrib = serviceProcess.getDistribution();
                            param.fileName = new ArrayList<String>(Collections.nCopies(R, null));
                            if (serviceDistrib instanceof Replayer) {
                                param.fileName.set(r, ((Replayer) serviceDistrib).getFileName());
                            }
                        }
                    }
                    break;
                case Queue:
                    // TODO
                    for (int r = 0; r < R; r++) {
                        ServiceBinding serviceProcess = node.getServer().getServiceProcess(this.getClassByIndex(r));
                        if (serviceProcess != null) {
                            Distribution serviceDistrib = serviceProcess.getDistribution();
                            param.fileName = new ArrayList<String>(Collections.nCopies(R, null));
                            if (serviceDistrib instanceof Replayer) {
                                param.fileName.set(r, ((Replayer) serviceDistrib).getFileName());
                            }
                        }
                    }
                    break;
                case Delay:
                    // TODO
                    for (int r = 0; r < R; r++) {
                        ServiceBinding serviceProcess = node.getServer().getServiceProcess(this.getClassByIndex(r));
                        Distribution serviceDistrib = serviceProcess.getDistribution();
                        param.fileName = new ArrayList<String>(Collections.nCopies(R, null));
                        if (serviceDistrib instanceof Replayer) {
                            param.fileName.set(r, ((Replayer) serviceDistrib).getFileName());
                        }
                    }
                    break;
                case Transition:
                    // TODO
//	    			for(int r = 0; r < R; r++) {
//	    				Distribution distrb = node.getServer().getServiceDistribution(this.jobClasses.get(r));
//	    				if (distrb instanceof MAP)
//	    					throw new RuntimeException("MAP is not supported in JLINE");
//	    				else if (distrb instanceof Replayer)
//	    					throw new RuntimeException("Replayer is not supported in JLINE");
//	    			}
                    break;
                default:
                    break;
            }

            for (int r = 0; r < R; r++) {
                JobClass jobclass = this.jobClasses.get(r);
                switch (this.sn.routing.get(node).get(jobclass)) {
                    case KCHOICES:
                        throw new RuntimeException("Routing Strategy KCHOICES is not supported in JLINE");
                    case WRROBIN:
                        param.weights = new HashMap<JobClass, Matrix>();
                        param.outlinks = new HashMap<JobClass, Matrix>();
                        nvars.set(i, R + r, nvars.get(i, R + r) + 1);

                        //varsparam{ind}{r}.weights = zeros(1,self.sn.nnodes);
                        param.weights.put(jobclass, new Matrix(1, this.sn.nnodes));
                        //varsparam{ind}{r}.outlinks = find(self.sn.connmatrix(ind,:));
                        Matrix conn_i = new Matrix(0, 0);
                        Matrix.extractRows(this.sn.connmatrix, i, i + 1, conn_i);
                        Matrix conn_i_transpose = conn_i.find().transpose();
                        param.outlinks.put(jobclass, conn_i_transpose);

                        List<OutputStrategy> outputStrategy_r = node.getOutput().getOutputStrategyByClass(jobclass);
                        for (int c = 0; c < outputStrategy_r.size(); c++) {
                            Node destination = outputStrategy_r.get(c).getDestination();
                            Double weight = outputStrategy_r.get(c).getProbability();
                            param.weights.get(jobclass).set(0, destination.getNodeIdx(), weight);
                        }
                        break;
                    case RROBIN:
                        param.outlinks = new HashMap<JobClass, Matrix>();
                        nvars.set(i, R + r, nvars.get(i, R + r) + 1);

                        //varsparam{ind}{r}.outlinks = find(self.sn.connmatrix(ind,:));
                        conn_i = new Matrix(0, 0);
                        Matrix.extractRows(this.sn.connmatrix, i, i + 1, conn_i);
                        conn_i_transpose = conn_i.find().transpose();
                        param.outlinks.put(jobclass, conn_i_transpose);
                        break;
                    default:
                        break;
                }
            }
            nodeparam.put(node, param);
        }

        if (this.sn != null) {
            this.sn.nvars = nvars;
            this.sn.nodeparam = nodeparam;
        }
    }

    public void refreshSync() {
        int local = this.nodes.size();
        int nclasses = this.sn.nclasses;
        Map<Integer, Sync> sync = new HashMap<Integer, Sync>();    //Index starts from 0
        Map<Node, Matrix> emptystate = new HashMap<Node, Matrix>();
        for (Node node : this.nodes)
            emptystate.put(node, new Matrix(0, 0));

        Matrix rtmask;
        if (this.sn.isstatedep.getNonZeroLength() > 0) {
            rtmask = this.sn.rtfun.apply(new Pair<Map<Node, Matrix>, Map<Node, Matrix>>(emptystate, emptystate));
        } else {
            //ceil(self.sn.rt);
            Matrix rt = this.sn.rt;
            rtmask = new Matrix(rt.getNumRows(), rt.getNumCols());
            for (int colIdx = 0; colIdx < rt.getNumCols(); colIdx++) {
                int col1 = rt.getColIndexes()[colIdx];
                int col2 = rt.getColIndexes()[colIdx + 1];

                for (int i = col1; i < col2; i++) {
                    int rowIdx = rt.getNonZeroRows()[i];
                    double value = rt.getNonZeroValues()[i];
                    rtmask.set(rowIdx, colIdx, Math.ceil(value));
                }
            }
        }

        for (int i = 0; i < sn.nnodes; i++) {
            for (int r = 0; r < nclasses; r++) {
                if (sn.isstation.get(i, 0) > 0 && sn.phases.get((int) sn.nodeToStation.get(0, i), r) > 1) {
                    Sync synct = new Sync();
                    synct.active.put(0, new NetworkEvent(EventType.PHASE, i, r, Double.NaN, new Matrix(0, 0), Double.NaN, Double.NaN));
                    synct.passive.put(0, new NetworkEvent(EventType.LOCAL, local, r, 1.0, new Matrix(0, 0), Double.NaN, Double.NaN));
                    sync.put(sync.size(), synct);
                }
                if (sn.isstateful.get(i, 0) > 0) {
                    //Line 24 - 29 is ignored since cache node
                    int isf = (int) sn.nodeToStateful.get(0, i);
                    for (int j = 0; j < sn.nnodes; j++) {
                        if (sn.isstateful.get(j, 0) > 0) {
                            int jsf = (int) sn.nodeToStateful.get(0, j);
                            for (int s = 0; s < nclasses; s++) {
                                double p = rtmask.get(isf * nclasses + r, jsf * nclasses + s);
                                if (p > 0) {
                                    Sync synct = new Sync();
                                    synct.active.put(0, new NetworkEvent(EventType.DEP, i, r, Double.NaN, new Matrix(0, 0), Double.NaN, Double.NaN));
                                    switch (sn.routing.get(this.nodes.get(i)).get(this.jobClasses.get(s))) {
                                        case RROBIN:
                                        case WRROBIN:
                                        case JSQ:
                                            final int i_final = i, j_final = j, r_final = r, s_final = s;
                                            synct.passive.put(0, new NetworkEvent(EventType.ARV, j, s,
                                                    ((pair) -> sn.rtfun.apply(pair).get(i_final * nclasses + r_final, j_final * nclasses + s_final)),
                                                    new Matrix(0, 0), Double.NaN, Double.NaN));
                                            break;
                                        default:
                                            synct.passive.put(0, new NetworkEvent(EventType.ARV, j, s,
                                                    sn.rt.get(isf * nclasses + r, jsf * nclasses + s),
                                                    new Matrix(0, 0), Double.NaN, Double.NaN));
                                    }
                                    sync.put(sync.size(), synct);
                                }
                            }
                        }
                    }
                }
            }
        }

        if (this.sn != null)
            this.sn.sync = sync;
    }

    public double sub_rr_wrr(int ind, int jnd, int r, int s, Matrix linksmat, Map<Node, Matrix> state_before, Map<Node, Matrix> state_after) {
        int R = this.sn.nclasses;
        int isf = (int) this.sn.nodeToStateful.get(0, ind);
        Node statefulNode = this.getStatefulNodeFromIndex(isf);
        if (!state_before.containsKey(statefulNode)) {
            return Math.min(linksmat.get(ind, jnd), 1.0);
        } else {
            if (r == s) {
                Matrix jm = state_after.get(statefulNode);
                return ((int) jm.get(jm.getNumCols() - 1 - R + r) == jnd) ? 1.0 : 0.0;
            } else {
                return 0.0;
            }
        }
    }

    public double sub_jsq(int ind, int jnd, int r, int s, Matrix linksmat, Map<Node, Matrix> state_before, Map<Node, Matrix> state_after) {
        int isf = (int) this.sn.nodeToStateful.get(0, ind);
        Node statefulNode = this.getStatefulNodeFromIndex(isf);
        if (!state_before.containsKey(statefulNode)) {
            return Math.min(linksmat.get(ind, jnd), 1.0);
        } else {
            if (r == s) {
                Matrix n = new Matrix(1, this.sn.nnodes);
                n.fill(Double.POSITIVE_INFINITY);
                for (int knd = 0; knd < this.sn.nnodes; knd++) {
                    if (linksmat.get(ind, knd) > 0) {
                        Node statefulNode_knd = this.getStatefulNodeFromIndex((int) this.sn.nodeToStateful.get(0, knd));
                        n.set(0, knd, State.toMarginal(this.sn, knd, state_before.get(statefulNode), null, null, null, null, null).ni.get(0, 0));
                    }
                }
                double min = n.elementMin();
                if (n.get(jnd) == min)
                    return 1.0 / n.count(min);
                else
                    return 0.0;
            } else {
                return 0.0;
            }
        }
    }

    public SolverHandles getAvgHandles() {
        int M = this.stations.size();
        int K = this.jobClasses.size();

        Matrix isSource = new Matrix(M, 1);
        Matrix isSink = new Matrix(M, 1);
        Matrix hasServiceTunnel = new Matrix(M, 1);
        Matrix isServiceDefined = Matrix.ones(M, K);

        for (int i = 0; i < M; i++) {
            if (this.stations.get(i) instanceof Source)
                isSource.set(i, 0, 1);
            if (((Node) this.stations.get(i)) instanceof Sink)
                isSink.set(i, 0, 1);

            if (this.stations.get(i).getServer() instanceof ServiceTunnel)
                hasServiceTunnel.set(i, 0, 1);
            else {
                for (int r = 0; r < K; r++) {
//                    if (!this.stations.get(i).getServer().containsJobClass(this.jobClasses.get(r)))
                    if (this.stations.get(i).getServer().getServiceDistribution(this.jobClasses.get(r)).isDisabled())
                        isServiceDefined.remove(i, r);
                }
            }
        }

        //Calculate Q
        Map<Station, Map<JobClass, SolverHandles.Metric>> Q = new HashMap<>();
        for (int i = 0; i < M; i++) {
            Map<JobClass, SolverHandles.Metric> map = new HashMap<>();
            for (int r = 0; r < K; r++) {
                SolverHandles.Metric Qir = new SolverHandles.Metric();
                Qir.type = "Number of Customers";
                Qir.jobClass = this.jobClasses.get(r);
                Qir.station = this.stations.get(i);
                if (isSource.get(i, 0) > 0)
                    Qir.isDisabled = true;
                else if (isSink.get(i, 0) > 0)
                    Qir.isDisabled = true;
                else Qir.isDisabled = hasServiceTunnel.get(i, 0) == 0 && isServiceDefined.get(i, r) == 0;
                map.put(this.jobClasses.get(r), Qir);
            }
            Q.put(this.stations.get(i), map);
        }

        //Calculate U
        Map<Station, Map<JobClass, SolverHandles.Metric>> U = new HashMap<>();
        for (int i = 0; i < M; i++) {
            Map<JobClass, SolverHandles.Metric> map = new HashMap<>();
            for (int r = 0; r < K; r++) {
                SolverHandles.Metric Uir = new SolverHandles.Metric();
                Uir.type = "Utilization";
                Uir.jobClass = this.jobClasses.get(r);
                Uir.station = this.stations.get(i);
                if (isSource.get(i, 0) > 0)
                    Uir.isDisabled = true;
                else if (isSink.get(i, 0) > 0)
                    Uir.isDisabled = true;
                else if (this.stations.get(i) instanceof Join)
                    Uir.isDisabled = true;
                else Uir.isDisabled = hasServiceTunnel.get(i, 0) == 0 && isServiceDefined.get(i, r) == 0;
                map.put(this.jobClasses.get(r), Uir);
            }
            U.put(this.stations.get(i), map);
        }

        //Calculate R
        Map<Station, Map<JobClass, SolverHandles.Metric>> R = new HashMap<>();
        for (int i = 0; i < M; i++) {
            Map<JobClass, SolverHandles.Metric> map = new HashMap<>();
            for (int r = 0; r < K; r++) {
                SolverHandles.Metric Rir = new SolverHandles.Metric();
                Rir.type = "Response Time";
                Rir.jobClass = this.jobClasses.get(r);
                Rir.station = this.stations.get(i);
                if (isSource.get(i, 0) > 0)
                    Rir.isDisabled = true;
                else if (isSink.get(i, 0) > 0)
                    Rir.isDisabled = true;
                else Rir.isDisabled = hasServiceTunnel.get(i, 0) == 0 && isServiceDefined.get(i, r) == 0;
                map.put(this.jobClasses.get(r), Rir);
            }
            R.put(this.stations.get(i), map);
        }

        //Calculate T
        Map<Station, Map<JobClass, SolverHandles.Metric>> T = new HashMap<>();
        for (int i = 0; i < M; i++) {
            Map<JobClass, SolverHandles.Metric> map = new HashMap<>();
            for (int r = 0; r < K; r++) {
                SolverHandles.Metric Tir = new SolverHandles.Metric();
                Tir.type = "Throughput";
                Tir.jobClass = this.jobClasses.get(r);
                Tir.station = this.stations.get(i);
                Tir.isDisabled = hasServiceTunnel.get(i, 0) == 0 && isServiceDefined.get(i, r) == 0;
                map.put(this.jobClasses.get(r), Tir);
            }
            T.put(this.stations.get(i), map);
        }

        //Calculate A
        Map<Station, Map<JobClass, SolverHandles.Metric>> A = new HashMap<>();
        for (int i = 0; i < M; i++) {
            Map<JobClass, SolverHandles.Metric> map = new HashMap<>();
            for (int r = 0; r < K; r++) {
                SolverHandles.Metric Air = new SolverHandles.Metric();
                Air.type = "Arrival Rate";
                Air.jobClass = this.jobClasses.get(r);
                Air.station = this.stations.get(i);
                Air.isDisabled = hasServiceTunnel.get(i, 0) == 0 && isServiceDefined.get(i, r) == 0;
                map.put(this.jobClasses.get(r), Air);
            }
            A.put(this.stations.get(i), map);
        }

        return new SolverHandles(Q, U, R, T, A);
    }

    public void initFromMarginal(Matrix n) {

        if (!hasStruct) {
            refreshStruct(true);
        }
        getStruct(true);

        if (!State.isValid(sn, n, new Matrix(0, 0))) {
            System.err.println("Initial state not contained in the state spac. Trying to recover.");
            for (int row = 0; row < n.getNumRows(); row++) {
                for (int col = 0; col < n.getNumCols(); col++) {
                    n.set(row, col, Math.round(n.get(row, col)));
                }
            }
            if (!State.isValid(sn, n, new Matrix(0, 0))) {
                throw new RuntimeException("Cannot recover - stopping.");
            }
        }

        for (int i = 0; i < sn.nnodes; i++) {
            if (sn.isstateful.get(i) == 1) {
                int ist = (int) sn.nodeToStation.get(i);
                if (sn.nodetypes.get(i) == NodeType.Place) {
                    // Must be single class token
                    stations.get(i).setState(n.sumRows(0, n.getNumCols()));
                } else {
                    stations.get(i).setState(State.fromMarginal(
                            sn,
                            i,
                            new Matrix(Matrix.extractRows(n, ist, ist + 1, null))));
                }
                if (stations.get(i).getState().isEmpty()) {
                    throw new RuntimeException("Invalid state assignment for a station.");
                }
            }
        }

        hasState = true;
    }

    public boolean hasFork() {
        for (NodeType type : this.getNodeTypes()) {
            if (type == NodeType.Fork) {
                return true;
            }
        }
        return false;
    }

    public boolean hasJoin() {
        for (NodeType type : this.getNodeTypes()) {
            if (type == NodeType.Join) {
                return true;
            }
        }
        return false;
    }

    public List<JobClass> getJobClass() {
        return this.jobClasses;
    }

    public RoutingMatrix initRoutingMatrix() {
        RoutingMatrix rt = new RoutingMatrix(this, jobClasses, nodes);
        for (JobClass j : jobClasses) {
            if (j instanceof SelfLoopingClass) {
                for (Node i : nodes) {
                    //rt.set(j, j, ((SelfLoopingClass) j).getRefstat(), ((SelfLoopingClass) j).getRefstat(), 1.0);
                    rt.set(j, j, i, i, 1.0);
                }
            }
        }
        return rt;
    }

    public NetworkAttribute getAttribute() {
        return attribute;
    }

    public String getLogPath() {
        return logPath;
    }

    public List<FiniteCapacityRegion> getRegions() {
        return regions;
    }


    /**
     * Adds an item set to the current list of items
     *
     * @param itemSet - the item set to be added
     */
    public void addItemSet(ItemSet itemSet) {
        this.items.forEach(s -> {
            if (s.name.equals(itemSet.name)) {
                throw new RuntimeException("An item type with name " + s.name + " already exists.");
            }
        });
        int nItemSet = this.items.size();
        itemSet.setIndex(nItemSet);
        this.items.add(itemSet);
    }

    public void resetHandles() {
        this.handles = new ArrayList<>();
    }

    public void resetModel(boolean resetState) {
        this.resetHandles();

        if (this.hasStruct) {
            Map<JobClass, Map<JobClass, Matrix>> rtorig;
            if (this.sn == null) {
                rtorig = null;
            } else {
                if (this.sn.rtorig != null) {
                    rtorig = new HashMap<>(this.sn.rtorig); // save linked routing table
                } else {
                    rtorig = null;
                }
            }
            this.sn = new NetworkStruct();
            this.sn.rtorig = rtorig;
            this.hasStruct = false;
        }

        if (resetState) {
            this.hasStruct = false;
        }

        for (int ind = 0; ind < this.getNumberOfNodes(); ind++) {
            this.nodes.get(ind).reset();
        }
    }

    public void reset() {
        this.resetModel(false);
        this.hasState = true;
    }

    public void reset(boolean resetState) {
        this.resetModel(resetState);
        this.hasState = true;
    }

    public void resetNetwork() {
        this.resetNetwork(true);
    }

    /**
     * Resets the topology of the current network
     *
     * @param deleteCSNodes - flag to indicate whether to delete the class switch nodes
     */
    public List<Node> resetNetwork(boolean deleteCSNodes) {
        int M = this.getNumberOfStations();

        // Remove class switch nodes
        if (deleteCSNodes) {
            List<Node> oldNodes = this.nodes;
            this.nodes = new ArrayList<>();
            for (Node n : oldNodes) {
                if (!(n instanceof ClassSwitch)) {
                    this.nodes.add(n);
                }
            }
        }

        for (int i = 0; i < M; i++) {
            ((Dispatcher) this.stations.get(i).getOutput()).initDispatcherJobClasses(this.getClasses());
        }

        this.handles = new ArrayList<>();
        this.connections = new Matrix(this.getNumberOfNodes(), this.getNumberOfNodes());
        return this.getNodes();
    }

    /**
     * Resets the struct of a given network
     */
    public void resetStruct() {
        this.sn = null;
        this.hasStruct = false;
    }

    /**
     * Gets the station index of the source
     *
     * @return -
     */
    public int getIndexSourceStation() {
        for (int i = 0; i < stations.size(); i++) {
            if (stations.get(i) instanceof Source) {
                return i;
            }
        }
        return -1;
    }

    public static class routingMatrixReturn {

        public Matrix rt;
        public Matrix rtnodes;
        public Matrix linksmat;
        public Matrix chains;
        public Map<JobClass, Map<JobClass, Matrix>> rtNodesByClass;
        public Map<Node, Map<Node, Matrix>> rtNodesByStation;

        public routingMatrixReturn(Matrix rt, Matrix rtnodes, Matrix linksmat,
                                   Matrix chains, Map<JobClass, Map<JobClass, Matrix>> rtNodesByClass, Map<Node, Map<Node, Matrix>> rtNodesByStation) {
            this.rt = rt;
            this.rtnodes = rtnodes;
            this.linksmat = linksmat;
            this.chains = chains;
            this.rtNodesByClass = rtNodesByClass;
            this.rtNodesByStation = rtNodesByStation;
        }
    }

    /**
     * Returns the language features used by the given network
     *
     * @return - the language features used by the given network
     */
    public FeatureSet getUsedLangFeatures() {
        FeatureSet s = new FeatureSet();
        HashMap<String, Boolean> f = new HashMap<>();
        if (!this.getIndexClosedClasses().isEmpty()) {
            f.put("ClosedClass", true);
        }
        if (!this.getIndexOpenClasses().isEmpty()) {
            f.put("OpenClass", true);
        }
        for (int i = 0; i < getNumberOfNodes(); i++) {
            for (int r = 0; r < getNumberOfClasses(); r++) {
                Node n = this.nodes.get(i);
                if (n instanceof Queue || n instanceof Delay) {
                    ServiceBinding serviceProcess = n.getServer().getServiceProcess(this.getClassByIndex(r));
                    if (serviceProcess != null) {
                        if (!(serviceProcess.getDistribution() instanceof Disabled) && !(serviceProcess.getDistribution() instanceof Immediate)) {
                            f.put(serviceProcess.getDistribution().getName(), true);
                        }
                        String sched = "";
                        if (n instanceof Delay) {
                            f.put("Delay", true);
                            sched = SchedStrategy.toFeature(((Delay) n).getSchedStrategy());
                        } else {
                            f.put("Queue", true);
                            sched = SchedStrategy.toFeature(((Queue) n).getSchedStrategy());
                        }
                        if (sched.length() > 0) {
                            f.put(sched, true);
                        }
                        if (r < n.getOutput().getOutputStrategies().size()) {
                            String routing = RoutingStrategy.toFeature(n.getOutput().getOutputStrategies().get(r).getRoutingStrategy());
                            if (routing.length() > 0) {
                                f.put(routing, true);
                            }
                        }
                    }
                } else if (n instanceof Router) {
                    if (r < n.getOutput().getOutputStrategies().size()) {
                        String routing = RoutingStrategy.toFeature(n.getOutput().getOutputStrategies().get(r).getRoutingStrategy());
                        if (routing.length() > 0) {
                            f.put(routing, true);
                        }
                    }
                } else if (n instanceof Source) {
                    Distribution serviceProcess = ((Source) n).getServiceProcess(this.getClassByIndex(r));
                    if (!(serviceProcess instanceof Disabled) && !(serviceProcess instanceof Immediate)) {
                        f.put(serviceProcess.getName(), true);
                    }
                    f.put("Source", true);
                } else if (n instanceof ClassSwitch) {
                    f.put("StatelessClassSwitcher", true);
                    f.put("ClassSwitch", true);
                } else if (n instanceof Fork) {
                    f.put("Fork", true);
                    f.put("Forker", true);
                } else if (n instanceof Join) {
                    f.put("Join", true);
                    f.put("Joiner", true);
                } else if (n instanceof Sink) {
                    f.put("Sink", true);
                } else if (n instanceof Cache) {
                    f.put("CacheClassSwitcher", true);
                    f.put("Cache", true);
                } /*else if(n instanceof Transition){
                TODO: implement transition
                } else if(n instaneof Place){
                TODO: implement place
                }*/
            }
        }
        s.setTrue(f.keySet().toArray(new String[0]));
        return s;
    }

    public void setNodeScheduling(int nodeIdx, JobClass jobClass, SchedStrategy schedStrategy) {
        this.nodes.get(nodeIdx).setScheduling(jobClass, schedStrategy);
    }

    public void setJoinNodeStrategy(int nodeIdx, JobClass jobClass, JoinStrategy joinStrategy) {
        Join join = (Join) this.nodes.get(nodeIdx);
        join.setStrategy(jobClass, joinStrategy);
        this.nodes.set(nodeIdx, join);
    }

    public void setJoinNodeRequired(int nodeIdx, JobClass jobClass, int njobs) {
        Join join = (Join) this.nodes.get(nodeIdx);
        join.setRequired(jobClass, njobs);
        this.nodes.set(nodeIdx, join);
    }

    public void setNodeRouting(int nodeIdx, JobClass jobClass, RoutingStrategy routingStrategy) {
        Node node = this.nodes.get(nodeIdx);
        node.setRouting(jobClass, routingStrategy);
        this.nodes.set(nodeIdx, node);
    }

    public void jsimgView() {
        SolverJMT jmt = new SolverJMT(this);
        jmt.jsimgView();
    }

    public boolean hasFCFS() {
        return SN.snHasFCFS(getStruct());
    }

    public boolean hasHomogeneousScheduling(SchedStrategy strategy) {
        return SN.snHasHomogeneousScheduling(getStruct(), strategy);
    }

    public boolean hasDPS() {
        return SN.snHasDPS(getStruct());
    }

    public boolean hasGPS() {
        return SN.snHasGPS(getStruct());
    }

    public boolean hasINF() {
        return SN.snHasINF(getStruct());
    }

    public boolean hasPS() {
        return SN.snHasPS(getStruct());
    }

    public boolean hasRAND() {
        return SN.snHasRAND(getStruct());
    }

    public boolean hasHOL() {
        return SN.snHasHOL(getStruct());
    }

    public boolean hasLCFS() {
        return SN.snHasLCFS(getStruct());
    }

    public boolean hasSEPT() {
        return SN.snHasSEPT(getStruct());
    }

    public boolean hasLEPT() {
        return SN.snHasLEPT(getStruct());
    }

    public boolean hasSJF() {
        return SN.snHasSJF(getStruct());
    }

    public boolean hasLJF() {
        return SN.snHasLJF(getStruct());
    }

    public boolean hasMultiClassFCFS() {
        return SN.snHasMultiClassFCFS(getStruct());
    }

    public boolean hasMultiClassHeterFCFS() {
        return SN.snHasMultiClassHeterFCFS(getStruct());
    }

    public boolean hasMultiServer() {
        return SN.snHasMultiServer(getStruct());
    }

    public boolean hasSingleChain() {
        return SN.snHasSingleChain(getStruct());
    }

    public boolean hasMultiChain() {
        return SN.snHasMultiChain(getStruct());
    }

    public boolean hasSingleClass() {
        return SN.snHasSingleClass(getStruct());
    }

    public boolean hasMultiClass() {
        return SN.snHasMultiClass(getStruct());
    }

    public boolean hasClassSwitching() {
        return SN.snHasClassSwitching(getStruct());
    }

}


