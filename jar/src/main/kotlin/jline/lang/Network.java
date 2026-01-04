/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import static jline.GlobalConstants.Inf;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.io.Ret;
import jline.io.ModelVisualizer;
import jline.lang.constant.*;
import jline.lang.nodeparam.*;
import jline.lang.nodes.*;
import jline.lang.nodes.Queue;
import jline.lang.processes.*;
import jline.lang.sections.*;
import jline.lang.state.FromMarginal;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.solvers.AvgHandle;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverAvgHandles;
import jline.solvers.SolverTranHandles;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.mva.SolverMVA;
import jline.util.Maths;
import jline.util.Pair;
import jline.lang.reward.RewardFunction;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

import javax.xml.parsers.ParserConfigurationException;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

import static java.lang.Double.NaN;
import static java.lang.Double.isFinite;
import static jline.io.InputOutputKt.line_debug;
import static jline.api.mam.Map_pieKt.map_pie;
import static jline.api.mam.Map_erlangKt.map_erlang;
import static jline.api.mc.Dtmc_stochcompKt.dtmc_stochcomp;
import static jline.api.sn.SnGetDemandsChainKt.snGetDemandsChain;
import static jline.api.sn.SnGetProductFormChainParamsKt.snGetProductFormChainParams;
import static jline.api.sn.SnGetProductFormParamsKt.snGetProductFormParams;
import static jline.api.sn.SnHasClassSwitchingKt.snHasClassSwitching;
import static jline.api.sn.SnHasDPSKt.snHasDPS;
import static jline.api.sn.SnHasDPSPRIOKt.snHasDPSPRIO;
import static jline.api.sn.SnHasFCFSKt.snHasFCFS;
import static jline.api.sn.SnHasGPSKt.snHasGPS;
import static jline.api.sn.SnHasGPSPRIOKt.snHasGPSPRIO;
import static jline.api.sn.SnHasHOLKt.snHasHOL;
import static jline.api.sn.SnHasHomogeneousSchedulingKt.snHasHomogeneousScheduling;
import static jline.api.sn.SnHasINFKt.snHasINF;
import static jline.api.sn.SnHasLCFSKt.snHasLCFS;
import static jline.api.sn.SnHasLCFSPRKt.snHasLCFSPR;
import static jline.api.sn.SnHasLEPTKt.snHasLEPT;
import static jline.api.sn.SnHasLJFKt.snHasLJF;
import static jline.api.sn.SnHasMultiChainKt.snHasMultiChain;
import static jline.api.sn.SnHasMultiClassFCFSKt.snHasMultiClassFCFS;
import static jline.api.sn.SnHasMultiClassHeterFCFSKt.snHasMultiClassHeterFCFS;
import static jline.api.sn.SnHasMultiClassKt.snHasMultiClass;
import static jline.api.sn.SnHasMultiServerKt.snHasMultiServer;
import static jline.api.sn.SnHasPSKt.snHasPS;
import static jline.api.sn.SnHasPSPRIOKt.snHasPSPRIO;
import static jline.api.sn.SnHasProductFormKt.snHasProductForm;
import static jline.api.sn.SnHasSEPTKt.snHasSEPT;
import static jline.api.sn.SnHasSIROKt.snHasSIRO;
import static jline.api.sn.SnHasSJFKt.snHasSJF;
import static jline.api.sn.SnHasSingleChainKt.snHasSingleChain;
import static jline.api.sn.SnHasSingleClassKt.snHasSingleClass;
import static jline.api.sn.SnIsStateValidKt.snIsStateValid;
import static jline.api.sn.SnRefreshVisitsKt.snRefreshVisits;
import static jline.api.sn.SnRtnodesToRtorigKt.snRtnodesToRtorig;
import static jline.io.InputOutputKt.*;
import static jline.lang.ModelAdapter.mmt;
import static jline.util.Maths.circul;
import static jline.util.Utils.isInf;

/**
 * A queueing network model
 * 
 * TABLE OF CONTENTS:
 * 1. FIELDS AND INITIALIZATION
 * 2. CONSTRUCTORS
 * 3. FACTORY METHODS
 * 4. NODE AND COMPONENT MANAGEMENT
 * 5. STATE AND CACHE MANAGEMENT
 * 6. GETTER METHODS - HANDLES
 * 7. GETTER METHODS - CLASSES AND CHAINS
 * 8. GETTER METHODS - CONFIGURATION
 * 9. GETTER METHODS - INDEXES
 * 10. GETTER METHODS - ROUTING
 * 11. GETTER METHODS - NODES
 * 12. GETTER METHODS - COUNTS
 * 13. GETTER METHODS - PROCESS AND PRODUCT FORM
 * 14. ROUTING MATRIX MANAGEMENT
 * 15. GETTER METHODS - SOURCES AND SINKS
 * 16. QUERY METHODS (HAS/IS)
 * 17. INITIALIZATION METHODS
 * 18. REFRESH METHODS
 * 19. RESET METHODS
 * 20. VALIDATION AND CONFIGURATION
 * 21. UTILITY METHODS
 * 22. INNER CLASSES
 */
public class Network extends Model implements Copyable {
    
    // ================================================================================
    // SECTION 1: FIELDS AND INITIALIZATION
    // ================================================================================
    // Private fields and member variables
    private final List<JobClass> jobClasses;
    private final List<Station> stations;
    private final NetworkAttribute attribute;
    private final List<ItemSet> items;
    private final List<Chain> chains;
    private final List<Region> regions;
    private List<StatefulNode> stateful;
    private boolean enableChecks;
    private boolean hasState;
    private String logPath;
    private FeatureSet usedFeatures;
    private List<Node> nodes;
    private boolean hasStruct;
    private NetworkStruct sn;
    private Matrix csMatrix;
    private Matrix connections;
    private boolean allowReplace;
    private boolean initializingState;

    private List<Object> handles;
    // caches
    private Map<Node, Map<JobClass, List<Node>>> classLinks;

    // ================================================================================
    // SECTION 2: CONSTRUCTORS
    // ================================================================================
    // Network initialization and construction
    
    /**
     * Creates a new queueing network model with the specified name.
     * Initializes all internal data structures and sets default configuration.
     *
     * @param modelName the name for this network model
     */
    public Network(String modelName) {
        super(modelName);

        this.hasState = false;
        this.enableChecks = true;

        this.nodes = new ArrayList<Node>();
        this.jobClasses = new ArrayList<JobClass>();
        this.stations = new ArrayList<Station>();
        this.stateful = new ArrayList<StatefulNode>();
        this.chains = new ArrayList<Chain>();

        this.classLinks = new HashMap<Node, Map<JobClass, List<Node>>>();
        this.items = new ArrayList<>();
        this.regions = new ArrayList<>();

        this.hasStruct = false;
        this.csMatrix = null;
        this.sn = null;
        this.connections = null;
        this.allowReplace = false;
        this.attribute = new NetworkAttribute();
    }

    // ========================================================================
    // SECTION 3: FACTORY METHODS
    // Static methods for creating network instances with specific configurations
    // ========================================================================

    /**
     * Creates a cyclic queueing network model with specified job populations, service demands, scheduling strategies, and server counts.
     * 
     * @param N        matrix of job populations for each class [1 x R] where R is the number of job classes
     * @param D        matrix of service demands [M x R] where M is the number of stations and R is the number of classes
     * @param strategy array of scheduling strategies for each station [M x 1]
     * @param S        matrix or vector specifying number of servers for each station [M x 1]
     * @return a configured closed queueing network model
     * @throws RuntimeException if server count matrix dimensions are invalid
     */
    public static Network cyclic(Matrix N, Matrix D, SchedStrategy[] strategy, Matrix S) {
        Network model = new Network("Model");
        int M = D.getNumRows();
        int R = D.getNumCols();

        if (S.getNumRows() == 1 && S.getNumCols() == M) {
            S = S.transpose();
        } else if (S.getNumRows() != M || S.getNumCols() != 1) {
            line_error(mfilename(new Object() {
            }), "The vector specifying the number of servers must be of size Mx1, where M is the number of stations.");
        }


        List<Node> nodes = new ArrayList<>(M);
        List<ClosedClass> jobclasses = new ArrayList<>(R);

        int nqueues = 0;
        int ndelays = 0;

        for (int i = 0; i < M; i++) {
            if (strategy[i] == SchedStrategy.INF) {
                ndelays++;
                nodes.add(new Delay(model, "Delay" + ndelays));
            } else {
                nqueues++;
                Queue queueNode = new Queue(model, "Queue" + nqueues, strategy[i]);
                queueNode.setNumberOfServers((int) S.get(i, 0));
                nodes.add(queueNode);
            }
        }

        for (int r = 0; r < R; r++) {
            ClosedClass newclass = new ClosedClass(model, "Class" + (r + 1), (int) N.get(0, r), (Station) nodes.get(0), 0);
            jobclasses.add(newclass);
        }

        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                ((Queue) nodes.get(i)).setService(jobclasses.get(r), Exp.fitMean(D.get(i, r)));
            }
        }

        RoutingMatrix P = model.initRoutingMatrix();
        for (int r = 0; r < R; r++) {
            P.set(jobclasses.get(r), circul(M));
        }

        model.link(P);
        return model;
    }

    /**
     * Creates a cyclic queueing network with First Come First Served (FCFS) scheduling at all stations.
     * 
     * @param N matrix of job populations for each class [1 x R] where R is the number of job classes
     * @param D matrix of service demands [M x R] where M is the number of stations and R is the number of classes
     * @return a configured closed queueing network model with FCFS scheduling
     */
    public static Network cyclicFcfs(Matrix N, Matrix D) {
        return Network.cyclicFcfs(N, D, Matrix.ones(D.getNumRows(), 1));
    }

    /**
     * Creates a cyclic queueing network with FCFS scheduling and specified server counts.
     * 
     * @param N matrix of job populations for each class [1 x R] where R is the number of job classes
     * @param D matrix of service demands [M x R] where M is the number of stations and R is the number of classes
     * @param S matrix specifying number of servers for each station [M x 1]
     * @return a configured closed queueing network model with FCFS scheduling
     */
    public static Network cyclicFcfs(Matrix N, Matrix D, Matrix S) {
        int M = D.getNumRows();
        SchedStrategy[] strategy = new SchedStrategy[M];
        for (int i = 0; i < M; i++) strategy[i] = SchedStrategy.FCFS;
        return Network.cyclic(N, D, strategy, S);
    }

    /**
     * Creates a cyclic network with infinite server (delay) stations followed by FCFS queue stations.
     * 
     * @param N matrix of job populations for each class [1 x R] where R is the number of job classes
     * @param D matrix of service demands at queue stations [M x R]
     * @param Z matrix of think times at delay stations [MZ x R] where MZ is the number of delay stations
     * @return a configured mixed queueing network model with delay and FCFS queue stations
     */
    public static Network cyclicFcfsInf(Matrix N, Matrix D, Matrix Z) {
        return cyclicFcfsInf(N, D, Z, Matrix.ones(D.getNumRows(), 1));
    }

    /**
     * Creates a cyclic network with infinite server stations followed by FCFS queue stations with specified server counts.
     * 
     * @param N matrix of job populations for each class [1 x R] where R is the number of job classes
     * @param D matrix of service demands at queue stations [M x R]
     * @param Z matrix of think times at delay stations [MZ x R] where MZ is the number of delay stations
     * @param S matrix specifying number of servers for queue stations [M x 1]
     * @return a configured mixed queueing network model with delay and FCFS queue stations
     */
    public static Network cyclicFcfsInf(Matrix N, Matrix D, Matrix Z, Matrix S) {
        int M = D.getNumRows();
        int MZ = Z.getNumRows();
        if (Z.elementMax() == 0) MZ = 0;
        int R = D.getNumCols();
        SchedStrategy[] strategy = new SchedStrategy[M + MZ];
        for (int i = 0; i < MZ; i++) strategy[i] = SchedStrategy.INF;
        for (int i = 0; i < M; i++) strategy[MZ + i] = SchedStrategy.FCFS;

        Matrix Dnew = new Matrix(M + MZ, R);
        if (MZ > 0) {
            Matrix.concatRows(Z, D, Dnew);
        } else {
            Dnew = D;
        }
        Matrix Snew = new Matrix(M + MZ, 1);
        Matrix Sinf = new Matrix(MZ, 1);
        Sinf.fill(Inf);
        if (MZ > 0) {
            Matrix.concatRows(Sinf, S, Snew);
        } else {
            Snew = S;
        }
        return Network.cyclic(N, Dnew, strategy, Snew);
    }

    /**
     * Creates a cyclic queueing network with Processor Sharing (PS) scheduling at all stations.
     * 
     * @param N matrix of job populations for each class [1 x R] where R is the number of job classes
     * @param D matrix of service demands [M x R] where M is the number of stations and R is the number of classes
     * @return a configured closed queueing network model with PS scheduling
     */
    public static Network cyclicPs(Matrix N, Matrix D) {
        return Network.cyclicPs(N, D, Matrix.ones(D.getNumRows(), 1));
    }

    /**
     * Creates a cyclic queueing network with PS scheduling and specified server counts.
     * 
     * @param N matrix of job populations for each class [1 x R] where R is the number of job classes
     * @param D matrix of service demands [M x R] where M is the number of stations and R is the number of classes
     * @param S matrix specifying number of servers for each station [M x 1]
     * @return a configured closed queueing network model with PS scheduling
     */
    public static Network cyclicPs(Matrix N, Matrix D, Matrix S) {
        int M = D.getNumRows();
        SchedStrategy[] strategy = new SchedStrategy[M];
        for (int i = 0; i < M; i++) strategy[i] = SchedStrategy.PS;
        return Network.cyclic(N, D, strategy, S);
    }

    /**
     * Creates a cyclic network with infinite server (delay) stations followed by PS queue stations.
     * 
     * @param N matrix of job populations for each class [1 x R] where R is the number of job classes
     * @param D matrix of service demands at queue stations [M x R]
     * @param Z matrix of think times at delay stations [MZ x R] where MZ is the number of delay stations
     * @return a configured mixed queueing network model with delay and PS queue stations
     */
    public static Network cyclicPsInf(Matrix N, Matrix D, Matrix Z) {
        return cyclicPsInf(N, D, Z, Matrix.ones(D.getNumRows(), 1));
    }

    /**
     * Creates a cyclic network with infinite server stations followed by PS queue stations with specified server counts.
     * 
     * @param N matrix of job populations for each class [1 x R] where R is the number of job classes
     * @param D matrix of service demands at queue stations [M x R]
     * @param Z matrix of think times at delay stations [MZ x R] where MZ is the number of delay stations
     * @param S matrix specifying number of servers for queue stations [M x 1]
     * @return a configured mixed queueing network model with delay and PS queue stations
     */
    public static Network cyclicPsInf(Matrix N, Matrix D, Matrix Z, Matrix S) {
        int M = D.getNumRows();
        int MZ = Z.getNumRows();
        if (Z.elementMax() == 0) MZ = 0;
        int R = D.getNumCols();
        SchedStrategy[] strategy = new SchedStrategy[M + MZ];
        for (int i = 0; i < MZ; i++) strategy[i] = SchedStrategy.INF;
        for (int i = 0; i < M; i++) strategy[MZ + i] = SchedStrategy.PS;

        Matrix Dnew = new Matrix(M + MZ, R);
        if (MZ > 0) {
            Matrix.concatRows(Z, D, Dnew);
        } else {
            Dnew = D;
        }
        Matrix Snew = new Matrix(M + MZ, 1);
        Matrix Sinf = new Matrix(MZ, 1);
        Sinf.fill(Inf);
        if (MZ > 0) {
            Matrix.concatRows(Sinf, S, Snew);
        } else {
            Snew = S;
        }
        return Network.cyclic(N, Dnew, strategy, Snew);
    }

    /**
     * Creates a serial routing matrix connecting nodes in sequence.
     * Jobs flow from each node to the next in the provided order.
     * The last node connects back to the first unless it's a Sink.
     *
     * @param jobClasses list of job classes to route
     * @param nodes      nodes to connect in serial order
     * @return routing matrix with serial connections
     */
    public static RoutingMatrix serialRouting(List<JobClass> jobClasses, Node... nodes) {
        if (nodes.length == 0) {
            return new RoutingMatrix();
        }

        Network network = nodes[0].model;
        RoutingMatrix outMatrix = new RoutingMatrix(network, jobClasses, network.nodes);

        for (int i = 1; i < nodes.length; i++) {
            //System.out.format("Loading connection %s->%s\n", nodes[i-1].getName(), nodes[i].getName());
            outMatrix.addConnection(nodes[i - 1], nodes[i], 1.0);
        }

        if (!(nodes[nodes.length - 1] instanceof Sink)) {
            outMatrix.addConnection(nodes[nodes.length - 1], nodes[0], 1.0);
        }

        return outMatrix;
    }

    /**
     * Creates a serial routing matrix connecting nodes in sequence.
     * Jobs flow from each node to the next in the provided order.
     * The last node connects back to the first unless it's a Sink.
     *
     * @param jobClasses list of job classes to route
     * @param nodes      list of nodes to connect in serial order
     * @return routing matrix with serial connections
     */
    public static RoutingMatrix serialRouting(List<JobClass> jobClasses, List<Node> nodes) {
        if (nodes.isEmpty()) {
            return new RoutingMatrix();
        }

        Network network = nodes.get(0).model;
        RoutingMatrix outMatrix = new RoutingMatrix(network, jobClasses, nodes);

        for (int i = 1; i < nodes.size(); i++) {
            //System.out.format("Loading connection %s->%s\n", nodes[i-1].getName(), nodes[i].getName());
            outMatrix.addConnection(nodes.get(i - 1), nodes.get(i), 1.0);
        }

        if (!(nodes.get(nodes.size() - 1) instanceof Sink)) {
            outMatrix.addConnection(nodes.get(nodes.size() - 1), nodes.get(0), 1.0);
        }

        return outMatrix;
    }

    /**
     * Creates a serial routing matrix for a single job class.
     *
     * @param jobClass the job class to route
     * @param nodes    nodes to connect in serial order
     * @return routing matrix with serial connections
     */
    public static RoutingMatrix serialRouting(JobClass jobClass, Node... nodes) {
        List<JobClass> jobClasses = new ArrayList<JobClass>();
        jobClasses.add(jobClass);

        return Network.serialRouting(jobClasses, nodes);
    }

    /**
     * Creates a serial routing matrix for a single job class.
     *
     * @param jobClass the job class to route
     * @param nodes    list of nodes to connect in serial order
     * @return routing matrix with serial connections
     */
    public static RoutingMatrix serialRouting(JobClass jobClass, List<Node> nodes) {
        List<JobClass> jobClasses = new ArrayList<JobClass>();
        jobClasses.add(jobClass);

        return Network.serialRouting(jobClasses, nodes);
    }

    /**
     * Creates a serial routing matrix for all job classes in the network.
     *
     * @param nodes nodes to connect in serial order
     * @return routing matrix with serial connections
     */
    public static RoutingMatrix serialRouting(Node... nodes) {
        if (nodes.length == 0) {
            return new RoutingMatrix();
        }
        Network network = nodes[0].model;
        return Network.serialRouting(network.jobClasses, nodes);
    }

    /**
     * Creates a serial routing matrix for all job classes in the network.
     *
     * @param nodes list of nodes to connect in serial order
     * @return routing matrix with serial connections
     */
    public static RoutingMatrix serialRouting(List<Node> nodes) {
        if (nodes.size() == 0) {
            return new RoutingMatrix();
        }
        Network network = nodes.get(0).model;
        return Network.serialRouting(network.jobClasses, nodes);
    }

    /**
     * Creates a tandem queueing network with specified arrival rates and service demands.
     *
     * @param lambda   matrix of arrival rates [classes x sources]
     * @param D        matrix of service demands [stations x classes]
     * @param strategy array of scheduling strategies for each station
     * @param S        matrix of server counts [stations x classes]
     * @return configured tandem network model
     */
    public static Network tandem(Matrix lambda, Matrix D, SchedStrategy[] strategy, Matrix S) {
        Network model = new Network("Model");
        int M = D.getNumRows();
        int R = D.getNumCols();

        List<Node> nodes = new ArrayList<>(M);
        List<OpenClass> jobclasses = new ArrayList<>(R);

        int nqueues = 0;
        int ndelays = 0;
        nodes.add(new Source(model, "Source"));

        for (int i = 0; i < M; i++) {
            if (strategy[i] == SchedStrategy.INF) {
                ndelays++;
                nodes.add(new Delay(model, "Delay" + ndelays));
            } else {
                nqueues++;
                Queue queueNode = new Queue(model, "Queue" + nqueues, strategy[i]);
                queueNode.setNumberOfServers((int) S.get(i, 0));
                nodes.add(queueNode);
            }
        }

        nodes.add(new Sink(model, "Sink"));

        for (int r = 0; r < R; r++) {
            OpenClass newclass = new OpenClass(model, "Class" + (r + 1), 0);
            jobclasses.add(newclass);
        }

        for (int r = 0; r < R; r++) {
            ((Source) nodes.get(0)).setArrival(jobclasses.get(r), Exp.fitMean(1.0 / lambda.get(0, r)));
            for (int i = 0; i < M; i++) {
                ((Queue) nodes.get(1 + i)).setService(jobclasses.get(r), Exp.fitMean(D.get(i, r)));
            }
        }

        RoutingMatrix P = model.initRoutingMatrix();
        for (int r = 0; r < R; r++) {
            Matrix Pr = circul(nodes.size());
            P.set(jobclasses.get(r), Pr);
        }

        model.link(P);
        return model;
    }

    public static Network tandemFcfs(Matrix lambda, Matrix D, Matrix S) {
        int M = D.getNumRows();
        SchedStrategy[] strategy = new SchedStrategy[M];
        for (int i = 0; i < M; i++) strategy[i] = SchedStrategy.FCFS;

        return Network.tandem(lambda, D, strategy, S);
    }

    /**
     * Creates a tandem network with FCFS infinite servers.
     *
     * @param lambda arrival rate matrix
     * @param D      service demand matrix
     * @return configured tandem FCFS infinite server network
     */
    public static Network tandemFcfsInf(Matrix lambda, Matrix D) {
        return tandemFcfsInf(lambda, D, new Matrix(""));
    }

    /**
     * Creates a tandem network with FCFS infinite servers and delay centers.
     *
     * @param lambda arrival rate matrix
     * @param D      service demand matrix
     * @param Z      delay time matrix
     * @return configured tandem FCFS infinite server network with delays
     */
    public static Network tandemFcfsInf(Matrix lambda, Matrix D, Matrix Z) {
        Matrix S = new Matrix(D.getNumRows(), 1, D.getNumRows());
        S.fill(1.0);
        return tandemFcfsInf(lambda, D, Z, S);
    }

    /**
     * Creates a tandem network with FCFS infinite servers, delays, and specified server counts.
     *
     * @param lambda arrival rate matrix
     * @param D      service demand matrix
     * @param Z      delay time matrix
     * @param S      server count matrix
     * @return configured tandem FCFS infinite server network
     */
    public static Network tandemFcfsInf(Matrix lambda, Matrix D, Matrix Z, Matrix S) {
        int M = D.getNumRows();
        int MZ = Z.getNumRows();
        if (Z.elementMax() == 0) MZ = 0;
        int R = D.getNumCols();
        SchedStrategy[] strategy = new SchedStrategy[M + MZ];
        for (int i = 0; i < MZ; i++) strategy[i] = SchedStrategy.INF;
        for (int i = 0; i < M; i++) strategy[MZ + i] = SchedStrategy.FCFS;

        Matrix Dnew = new Matrix(M + MZ, R);
        if (MZ > 0) {
            Matrix.concatRows(Z, D, Dnew);
        } else {
            Dnew = D;
        }
        Matrix Snew = new Matrix(M + MZ, 1);
        Matrix Sinf = new Matrix(MZ, 1);
        Sinf.fill(Inf);
        if (MZ > 0) {
            Matrix.concatRows(Sinf, S, Snew);
        } else {
            Snew = S;
        }
        return Network.tandem(lambda, Dnew, strategy, Snew);
    }

    public static Network tandemPs(Matrix lambda, Matrix D, Matrix S) {
        int M = D.getNumRows();
        SchedStrategy[] strategy = new SchedStrategy[M];
        for (int i = 0; i < M; i++) strategy[i] = SchedStrategy.PS;

        return Network.tandem(lambda, D, strategy, S);
    }

    /**
     * Creates a tandem network with processor sharing infinite servers.
     *
     * @param lambda arrival rate matrix
     * @param D      service demand matrix
     * @return configured tandem PS infinite server network
     */
    public static Network tandemPsInf(Matrix lambda, Matrix D) {
        return tandemPsInf(lambda, D, new Matrix(""));
    }

    /**
     * Creates a tandem network with processor sharing infinite servers and delays.
     *
     * @param lambda arrival rate matrix
     * @param D      service demand matrix
     * @param Z      delay time matrix
     * @return configured tandem PS infinite server network with delays
     */
    public static Network tandemPsInf(Matrix lambda, Matrix D, Matrix Z) {
        Matrix S = new Matrix(D.getNumRows(), 1, D.getNumRows());
        S.fill(1.0);
        return tandemPsInf(lambda, D, Z, S);
    }

    public static Network tandemPsInf(Matrix lambda, Matrix D, Matrix Z, Matrix S) {
        int M = D.getNumRows();
        int MZ = Z.getNumRows();
        if (Z.elementMax() == 0) MZ = 0;
        int R = D.getNumCols();
        SchedStrategy[] strategy = new SchedStrategy[M + MZ];
        for (int i = 0; i < MZ; i++) strategy[i] = SchedStrategy.INF;
        for (int i = 0; i < M; i++) strategy[MZ + i] = SchedStrategy.PS;

        Matrix Dnew = new Matrix(M + MZ, R);
        if (MZ > 0) {
            Matrix.concatRows(Z, D, Dnew);
        } else {
            Dnew = D;
        }
        Matrix Snew = new Matrix(M + MZ, 1);
        Matrix Sinf = new Matrix(MZ, 1);
        Sinf.fill(Inf);
        if (MZ > 0) {
            Matrix.concatRows(Sinf, S, Snew);
        } else {
            Snew = S;
        }
        return Network.tandem(lambda, Dnew, strategy, Snew);
    }

    /**
     * Adds an item set to the current list of items
     *
     * @param itemSet - the item set to be added
     */
    // ========================================================================
    // SECTION 4: NODE AND COMPONENT MANAGEMENT
    // Methods for adding and managing network nodes and components
    // ========================================================================

    /**
     * Adds an item set to the network model. Item sets define the types of resources or items 
     * that can be processed by the network nodes.
     * 
     * @param itemSet the item set to add to the network
     * @throws RuntimeException if an item set with the same name already exists
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

    /**
     * Adds a job class to the network model. Job classes define different types of jobs 
     * that traverse the network with potentially different service requirements and routing.
     * 
     * @param jobClass the job class to add to the network
     * @throws RuntimeException if a class with the same name already exists (when validation is enabled)
     */
    public void addJobClass(JobClass jobClass) {
        if (this.enableChecks) {
            if (!(this.getClassByName(jobClass.getName()) == null)) {
                line_error(mfilename(new Object() {
                }), "A class with name " + jobClass.getName() + " already exists.\n");
            }
        }
        this.jobClasses.add(jobClass);
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

    /**
     * Adds a node to this network.
     * If the node is a station, it's also added to the stations list.
     * If allowReplace is true and a node with the same name already exists, it will be replaced.
     *
     * @param node the node to add to the network
     * @return true if the node replaced an existing node, false otherwise
     */
    public boolean addNode(Node node) {
        // Check if a node with the same name already exists
        Node existingNode = null;
        int existingIndex = -1;
        for (int i = 0; i < nodes.size(); i++) {
            if (nodes.get(i).getName().equals(node.getName())) {
                existingNode = nodes.get(i);
                existingIndex = i;
                break;
            }
        }
        
        if (existingNode != null) {
            if (allowReplace) {
                // Replace the existing node
                node.setNodeIdx(existingNode.getNodeIndex());
                nodes.set(existingIndex, node);
                
                // If the node is a station, also replace it in the stations list
                if (node instanceof Station) {
                    if (existingNode instanceof Station) {
                        Station existingStation = (Station) existingNode;
                        node.setStationIdx(existingStation.getStationIdx());
                        stations.set(existingStation.getStationIdx(), (Station) node);
                    } else {
                        // The new node is a station but the existing one wasn't
                        node.setStationIdx(node.getStationIdx());
                        stations.add((Station) node);
                    }
                } else if (existingNode instanceof Station) {
                    // The existing node was a station but the new one isn't
                    // Remove the old station from the stations list
                    stations.remove((Station) existingNode);
                }
                return true; // Node was replaced
            } else {
                throw new RuntimeException("A node with an identical name already exists: " + node.getName());
            }
        } else {
            // Add the new node normally
            node.setNodeIdx(node.getNodeIndex()); // searches within the model nodes
            nodes.add(node);
            if (node instanceof Station) {
                node.setStationIdx(node.getStationIdx()); // searches within the model stations
                stations.add((Station) node);
            }
            return false; // Node was added, not replaced
        }
    }

    /**
     * Adds a finite capacity region to this network.
     *
     * @param nodes list of nodes forming the capacity region
     * @return the created finite capacity region
     */
    public Region addRegion(List<Node> nodes) {
        Region region = new Region(nodes, this.jobClasses);
        int regionIndex = this.regions.size() + 1;
        region.setName("FCR" + regionIndex);
        this.regions.add(region);
        return region;
    }

    // ========================================================================
    // SECTION 5: STATE AND CACHE MANAGEMENT
    // Methods for managing network state and clearing cached data
    // ========================================================================

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

    // ========================================================================
    // SECTION 6: GETTER METHODS - HANDLES
    // Methods for retrieving performance metrics and solver handles
    // ========================================================================

    public NetworkAttribute getAttribute() {
        return attribute;
    }

    public AvgHandle getAvgArvRHandles() {
        return this.getAvgHandles().getAvgArvRHandles();
    }

    public SolverAvgHandles getAvgHandles() {
        int M = this.stations.size();
        int K = this.jobClasses.size();

        Matrix isSource = new Matrix(M, 1);
        Matrix isSink = new Matrix(M, 1);
        Matrix hasServiceTunnel = new Matrix(M, 1);
        Matrix isServiceDefined = Matrix.ones(M, K);

        for (int i = 0; i < M; i++) {
            if (this.stations.get(i) instanceof Source) isSource.set(i, 0, 1);
            if (((Node) this.stations.get(i)) instanceof Sink) isSink.set(i, 0, 1);

            if (this.stations.get(i).getServer() instanceof ServiceTunnel) hasServiceTunnel.set(i, 0, 1);
            else {
                for (int r = 0; r < K; r++) {
//                    if (!this.stations.get(i).getServer().containsJobClass(this.jobClasses.get(r)))
                    if (this.stations.get(i) instanceof ServiceStation) {
                        // Check standard service distribution
                        boolean standardServiceDefined = !this.stations.get(i).getServer().getServiceDistribution(this.jobClasses.get(r)).isDisabled();
                        // Also check for heterogeneous service definitions
                        boolean heteroServiceDefined = false;
                        if (this.stations.get(i) instanceof Queue) {
                            Queue queue = (Queue) this.stations.get(i);
                            if (queue.isHeterogeneous()) {
                                // Check all server types for service to this class
                                Map<ServerType, Map<JobClass, Distribution>> heteroDistrs = queue.getHeteroServiceDistributions();
                                for (Map<JobClass, Distribution> classMap : heteroDistrs.values()) {
                                    Distribution distr = classMap.get(this.jobClasses.get(r));
                                    if (distr != null && !distr.isDisabled()) {
                                        heteroServiceDefined = true;
                                        break;
                                    }
                                }
                            }
                        }
                        if (!standardServiceDefined && !heteroServiceDefined) {
                            isServiceDefined.remove(i, r);
                        }
                    }
                }
            }
        }

        //Calculate Q
        AvgHandle Q = new AvgHandle();
        for (int i = 0; i < M; i++) {
            Map<JobClass, Metric> map = new HashMap<>();
            for (int r = 0; r < K; r++) {
                Metric Qir = new Metric();
                Qir.type = "Number of Customers";
                Qir.jobClass = this.jobClasses.get(r);
                Qir.station = this.stations.get(i);
                if (isSource.get(i, 0) > 0) Qir.isDisabled = true;
                else if (isSink.get(i, 0) > 0) Qir.isDisabled = true;
                else Qir.isDisabled = hasServiceTunnel.get(i, 0) == 0 && isServiceDefined.get(i, r) == 0;
                Q.put(this.stations.get(i), this.jobClasses.get(r), Qir);
            }
        }

        //Calculate U
        AvgHandle U = new AvgHandle();
        for (int i = 0; i < M; i++) {
            Map<JobClass, Metric> map = new HashMap<>();
            for (int r = 0; r < K; r++) {
                Metric Uir = new Metric();
                Uir.type = "Utilization";
                Uir.jobClass = this.jobClasses.get(r);
                Uir.station = this.stations.get(i);
                if (isSource.get(i, 0) > 0) Uir.isDisabled = true;
                else if (isSink.get(i, 0) > 0) Uir.isDisabled = true;
                else if (this.stations.get(i) instanceof Join) Uir.isDisabled = true;
                else Uir.isDisabled = hasServiceTunnel.get(i, 0) == 0 && isServiceDefined.get(i, r) == 0;
                U.put(this.stations.get(i), this.jobClasses.get(r), Uir);
            }
        }

        //Calculate R
        AvgHandle R = new AvgHandle();
        for (int i = 0; i < M; i++) {
            Map<JobClass, Metric> map = new HashMap<>();
            for (int r = 0; r < K; r++) {
                Metric Rir = new Metric();
                Rir.type = "Response Time";
                Rir.jobClass = this.jobClasses.get(r);
                Rir.station = this.stations.get(i);
                if (isSource.get(i, 0) > 0) Rir.isDisabled = true;
                else if (isSink.get(i, 0) > 0) Rir.isDisabled = true;
                else Rir.isDisabled = hasServiceTunnel.get(i, 0) == 0 && isServiceDefined.get(i, r) == 0;
                R.put(this.stations.get(i), this.jobClasses.get(r), Rir);
            }
        }

        //Calculate W
        AvgHandle W = new AvgHandle();
        for (int i = 0; i < M; i++) {
            Map<JobClass, Metric> map = new HashMap<>();
            for (int r = 0; r < K; r++) {
                Metric Wir = new Metric();
                Wir.type = "Residence Time";
                Wir.jobClass = this.jobClasses.get(r);
                Wir.station = this.stations.get(i);
                if (isSource.get(i, 0) > 0) Wir.isDisabled = true;
                else if (isSink.get(i, 0) > 0) Wir.isDisabled = true;
                else Wir.isDisabled = hasServiceTunnel.get(i, 0) == 0 && isServiceDefined.get(i, r) == 0;
                W.put(this.stations.get(i), this.jobClasses.get(r), Wir);
            }
        }

        //Calculate T
        AvgHandle T = new AvgHandle();
        for (int i = 0; i < M; i++) {
            Map<JobClass, Metric> map = new HashMap<>();
            for (int r = 0; r < K; r++) {
                Metric Tir = new Metric();
                Tir.type = "Throughput";
                Tir.jobClass = this.jobClasses.get(r);
                Tir.station = this.stations.get(i);
                Tir.isDisabled = hasServiceTunnel.get(i, 0) == 0 && isServiceDefined.get(i, r) == 0;
                T.put(this.stations.get(i), this.jobClasses.get(r), Tir);
            }
        }

        //Calculate A
        AvgHandle A = new AvgHandle();
        for (int i = 0; i < M; i++) {
            Map<JobClass, Metric> map = new HashMap<>();
            for (int r = 0; r < K; r++) {
                Metric Air = new Metric();
                Air.type = "Arrival Rate";
                Air.jobClass = this.jobClasses.get(r);
                Air.station = this.stations.get(i);
                Air.isDisabled = hasServiceTunnel.get(i, 0) == 0 && isServiceDefined.get(i, r) == 0;
                A.put(this.stations.get(i), this.jobClasses.get(r), Air);
            }
        }

        //Calculate Tard
        AvgHandle Tard = new AvgHandle();
        for (int i = 0; i < M; i++) {
            Map<JobClass, Metric> map = new HashMap<>();
            for (int r = 0; r < K; r++) {
                Metric Tardir = new Metric();
                Tardir.type = "Tardiness";
                Tardir.jobClass = this.jobClasses.get(r);
                Tardir.station = this.stations.get(i);
                if (isSource.get(i, 0) > 0) Tardir.isDisabled = true;
                else Tardir.isDisabled = hasServiceTunnel.get(i, 0) == 0 && isServiceDefined.get(i, r) == 0;
                Tard.put(this.stations.get(i), this.jobClasses.get(r), Tardir);
            }
        }

        //Calculate SysTard
        AvgHandle SysTard = new AvgHandle();
        for (int r = 0; r < K; r++) {
            Metric SysTardr = new Metric();
            SysTardr.type = "System Tardiness";
            SysTardr.jobClass = this.jobClasses.get(r);
            SysTardr.station = null;
            SysTard.put(null, this.jobClasses.get(r), SysTardr);
        }

        return new SolverAvgHandles(Q, U, R, W, T, A, Tard, SysTard);
    }

    public AvgHandle getAvgQLenHandles() {
        return this.getAvgHandles().getAvgQLenHandles();
    }

    public AvgHandle getAvgResidTHandles() {
        return this.getAvgHandles().getAvgResidTHandles();
    }

    public AvgHandle getAvgRespTHandles() {
        return this.getAvgHandles().getAvgRespTHandles();
    }

    public AvgHandle getAvgTputHandles() {
        return this.getAvgHandles().getAvgTputHandles();
    }

    public AvgHandle getAvgUtilHandles() {
        return this.getAvgHandles().getAvgUtilHandles();
    }

    /**
     * Returns the list of job chains in this network.
     * Job chains group classes that follow similar routing patterns.
     *
     * @return list of job chains
     */
    // ========================================================================
    // SECTION 7: GETTER METHODS - CLASSES AND CHAINS
    // Methods for retrieving job classes and chains information
    // ========================================================================

    public List<Chain> getChains() {
        NetworkStruct sn = this.getStruct();
        for (int c = 0; c < sn.chains.getNumElements(); c++) {
            List<JobClass> chainClasses = new ArrayList<>();
            for (int r = 0; r < sn.nclasses; r++) {
                if (sn.inchain.get(c).get(0, r) == 1) {
                    chainClasses.add(this.jobClasses.get(r));
                }
            }
            chains.add(new Chain("Chain" + c, chainClasses, this.stations));
        }
        return chains;
    }

    public JobClass getClassByIndex(int index) {
        return this.jobClasses.get(index);
    }

    public JobClass getClassByName(String name) {
        for (JobClass jobClass : this.jobClasses) {
            if (jobClass.getName().equals(name)) {
                return jobClass;
            }
        }
        return null;
    }

    public Chain getClassChain(JobClass jobClass) {
        NetworkStruct sn = this.getStruct();
        for (int c = 0; c < sn.chains.getNumElements(); c++) {
            for (int r = 0; r < sn.nclasses; r++) {
                if (sn.inchain.get(c).get(0, r) == 1) return chains.get(c);
            }
        }
        return null;
    }

    public int getClassChainIndex(JobClass jobClass) {
        NetworkStruct sn = this.getStruct();
        for (int c = 0; c < sn.chains.getNumElements(); c++) {
            for (int r = 0; r < sn.nclasses; r++) {
                if (sn.inchain.get(c).get(0, r) == 1) return c;
            }
        }
        return -1;
    }

    public int getClassIndex(JobClass jobclass) {
        int outIdx = 0;
        for (JobClass classIter : this.jobClasses) {
            if (classIter == jobclass) {
                return outIdx;
            } else {
                outIdx++;
            }
        }
        return -1;
    }

    public int getClassIndex(String name) {
        int ret = -1;
        List<String> classNames = getClassNames();
        for (int i = 0; i < this.getNumberOfClasses(); i++) {
            if (classNames.get(i).equals(name)) {
                ret = i;
                break;
            }
        }
        return ret;
    }

    public int getClassLinks(Node node, JobClass jobClass) {
        if (this.classLinks.isEmpty()) {
            this.generateClassLinks();
        }
        return this.classLinks.get(node).get(jobClass).size();
    }

    public List<String> getClassNames() {
        if (hasStruct && sn.classnames != null) return sn.classnames;

        int K = getNumberOfClasses();
        List<String> classnames = new ArrayList<String>();
        for (int i = 0; i < K; i++)
            classnames.add(jobClasses.get(i).getName());

        return classnames;
    }

    public Matrix getClassSwitchingMask() {
        return this.getStruct().csmask;
    }

    /**
     * Returns the list of job classes in this network.
     *
     * @return list of job classes
     */
    public List<JobClass> getClasses() {
        return this.jobClasses;
    }

    public Matrix getConnectionMatrix() {
        if (this.connections == null || this.connections.isEmpty())
            this.connections = new Matrix(this.getNumberOfNodes(), this.getNumberOfNodes());
        if (this.connections.getNumCols() < this.getNumberOfNodes() || this.connections.getNumRows() < this.getNumberOfNodes())
            this.connections.expandMatrix(this.getNumberOfNodes(), this.getIndexSourceNode(), this.getNumberOfNodes() * this.getNumberOfNodes());
        return this.connections;
    }

    public void setConnectionMatrix(Matrix connection) {
        this.connections = connection;
    }

    public Matrix getCsMatrix() {
        return this.csMatrix;
    }

    public void setCsMatrix(Matrix csMatrix) {
        this.csMatrix = csMatrix;
    }

    // ========================================================================
    // SECTION 8: GETTER METHODS - CONFIGURATION
    // Methods for retrieving service demands and configuration data
    // ========================================================================

    public Ret.snGetDemands getDemands() {
        Ret.snGetProductFormParams ret = snGetProductFormParams(getStruct());
        return new Ret.snGetDemands(ret.D, ret.Z);
    }

    public Matrix getDemandsChain() {
        Ret.snGetDemands ret = snGetDemandsChain(getStruct());
        return ret.Dchain;
    }

    public Matrix getForkJoins() {
        int I = this.getNumberOfNodes();
        Matrix fjPairs = new Matrix(I, I);

        int forkCount = 0, joinCount = 0, linkedJoinCount = 0;
        for (int i = 0; i < I; i++) {
            Node node = this.nodes.get(i);
            if (node instanceof Fork) {
                forkCount++;
                //no-op
            } else if (node instanceof Join) {
                joinCount++;
                Join joinNode = (Join) node;
                if (joinNode.joinOf != null) {
                    fjPairs.set(joinNode.joinOf.getNodeIndex(), node.getNodeIndex(), 1.0);
                    linkedJoinCount++;
                } else {
                    line_warning(mfilename(new Object() {
                    }), String.format("Join node '%s' at index %d has null joinOf reference", 
                        joinNode.getName(), i));
                }
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

    public boolean getAllowReplace() {
        return this.allowReplace;
    }

    public void setAllowReplace(boolean allowReplace) {
        this.allowReplace = allowReplace;
    }

    /**
     * Returns the indices of all closed job classes in this network.
     *
     * @return list of indices for closed job classes
     */
    // ========================================================================
    // SECTION 9: GETTER METHODS - INDEXES
    // Methods for retrieving various index mappings
    // ========================================================================

    public List<Integer> getIndexClosedClasses() {
        List<Integer> outList = new ArrayList<Integer>();
        for (int i = 0; i < this.jobClasses.size(); i++) {
            if (this.jobClasses.get(i) instanceof ClosedClass) {
                outList.add(this.getJobClassIndex(this.jobClasses.get(i)));
            }
        }
        return outList;
    }

    /**
     * Returns the indices of all open job classes in this network.
     *
     * @return list of indices for open job classes
     */
    public List<Integer> getIndexOpenClasses() {
        List<Integer> outList = new ArrayList<Integer>();
        for (int i = 0; i < this.jobClasses.size(); i++) {
            if (this.jobClasses.get(i) instanceof OpenClass) {
                outList.add(this.getJobClassIndex(this.jobClasses.get(i)));
            }
        }
        return outList;
    }

    public int getIndexSinkNode() {
        int res = 0;
        for (Node nodeIter : this.nodes) {
            if (nodeIter instanceof Sink) return res;
            res++;
        }

        return -1;
    }

    public int getIndexSourceNode() {
        int res = 0;
        for (Node nodeIter : this.nodes) {
            if (nodeIter instanceof Source) return res;
            res++;
        }

        return -1;
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

    public List<Integer> getIndexStatefulNodes() {
        List<Integer> outList = new ArrayList<Integer>();
        for (int i = 0; i < this.nodes.size(); i++) {
            if (this.nodes.get(i) instanceof StatefulNode) {
                outList.add(i);
            }
        }
        return outList;
    }

    /**
     * Returns the job class at the specified index.
     *
     * @param inIdx index of the job class
     * @return job class at the given index
     * @throws IndexOutOfBoundsException if index is invalid
     */
    public JobClass getJobClassFromIndex(int inIdx) {
        return this.jobClasses.get(inIdx);
    }

    /**
     * Returns the index of the specified job class in this network.
     *
     * @param jobClass the job class to find
     * @return index of the job class, or -1 if not found
     */
    public int getJobClassIndex(JobClass jobClass) {
        return this.jobClasses.indexOf(jobClass);
    }

    public List<JobClass> getJobClasses() {
        return this.jobClasses;
    }

    public Map<Station, SerializableFunction<Matrix, Double>> getLimitedClassDependence() {
        Map<Station, SerializableFunction<Matrix, Double>> gamma = new HashMap<Station, SerializableFunction<Matrix, Double>>();

        for (Station station : this.stations) {
            if (station.getLimitedClassDependence() != null) gamma.put(station, station.getLimitedClassDependence());
        }

        if (!gamma.isEmpty()) {
            for (Station station : this.stations) {
                if (gamma.getOrDefault(station, null) == null) {
                    gamma.put(station, (nvec) -> {
                        return 1.0;
                    });
                }
            }
        }
        return gamma;
    }

    public Matrix getLimitedLoadDependence() {
        List<Matrix> mus = new ArrayList<Matrix>();
        int maxsize = 0;

        for (Station station : this.stations) {
            mus.add(station.getLimitedLoadDependence());
            maxsize = FastMath.max(maxsize, station.getLimitedLoadDependence().length());
        }

        int M = this.stations.size();
        Matrix alpha = new Matrix(M, maxsize);
        alpha.fill(1.0);
        for (int i = 0; i < M; i++) {
            Matrix mu = mus.get(i);
            if (!mu.isEmpty()) {
                Matrix.extract(mu, 0, 1, 0, mu.length(), alpha, i, 0);
                for (int j = 0; j < mu.length(); j++) {
                    if (alpha.get(i, j) == 0) alpha.set(i, j, 1.0);
                }
            }
        }
        return alpha;
    }

    /**
     * Gets the limited joint-dependent scaling tables and cutoffs for all stations.
     *
     * @return a Pair containing: (1) Map of station -> linearized scaling table,
     *         (2) Map of station -> per-class cutoffs
     */
    public Pair<Map<Station, Matrix>, Map<Station, Matrix>> getLimitedJointDependence() {
        Map<Station, Matrix> ljdscaling = new HashMap<Station, Matrix>();
        Map<Station, Matrix> ljdcutoffs = new HashMap<Station, Matrix>();

        for (Station station : this.stations) {
            Matrix scaling = station.getLimitedJointDependence();
            Matrix cutoffs = station.getLimitedJointDependenceCutoffs();
            if (scaling != null && !scaling.isEmpty()) {
                ljdscaling.put(station, scaling);
                ljdcutoffs.put(station, cutoffs);
            }
        }

        return new Pair<Map<Station, Matrix>, Map<Station, Matrix>>(ljdscaling, ljdcutoffs);
    }

    /**
     * Gets the limited joint-class-dependent scaling tables and cutoffs for all stations.
     * Each class has its own scaling table indexed by the population vector.
     *
     * @return a Pair containing: (1) Map of station -> (Map of class -> scaling table),
     *         (2) Map of station -> per-class cutoffs
     */
    public Pair<Map<Station, Map<JobClass, Matrix>>, Map<Station, Matrix>> getLimitedJointClassDependence() {
        Map<Station, Map<JobClass, Matrix>> ljcdscaling = new HashMap<Station, Map<JobClass, Matrix>>();
        Map<Station, Matrix> ljcdcutoffs = new HashMap<Station, Matrix>();

        for (Station station : this.stations) {
            Map<JobClass, Matrix> scaling = station.getLimitedJointClassDependence();
            Matrix cutoffs = station.getLimitedJointClassDependenceCutoffs();
            if (scaling != null && !scaling.isEmpty()) {
                ljcdscaling.put(station, scaling);
                ljcdcutoffs.put(station, cutoffs);
            }
        }

        return new Pair<Map<Station, Map<JobClass, Matrix>>, Map<Station, Matrix>>(ljcdscaling, ljcdcutoffs);
    }

    // ========================================================================
    // SECTION 10: GETTER METHODS - ROUTING
    // Methods for retrieving routing information and matrices
    // ========================================================================

    public Map<JobClass, Map<JobClass, Matrix>> getLinkedRoutingMatrix() {
        return getStruct(false).rtorig;
    }

    // ========================================================================
    // SECTION 11: GETTER METHODS - NODES
    // Methods for retrieving node-related information
    // ========================================================================

    public String getLogPath() {
        return logPath;
    }

    public void setLogPath(String logPath) {
        this.logPath = logPath;
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

    public Node getNodeByName(String name) {
        for (Node node : this.nodes) {
            if (node.getName().equals(name)) {
                return node;
            }
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
            if (node.getName().equals(name)) return getNodeIndex(node);
        }
        return -1;
    }

    public List<String> getNodeNames() {
        if (hasStruct && sn.classnames != null) return sn.nodenames;

        int M = getNumberOfNodes();
        List<String> nodenames = new ArrayList<String>();
        for (int i = 0; i < M; i++)
            nodenames.add(nodes.get(i).getName());

        return nodenames;
    }

    public List<NodeType> getNodeTypes() {
        int M = getNumberOfNodes();
        List<NodeType> nodetypes = new ArrayList<NodeType>(M);

        try {
            for (int i = 0; i < M; i++) {
                Node nodeIter = this.nodes.get(i);
                if (nodeIter instanceof Logger) nodetypes.add(NodeType.Logger);
                else if (nodeIter instanceof ClassSwitch) nodetypes.add(NodeType.ClassSwitch);
                else if (nodeIter instanceof Join) nodetypes.add(NodeType.Join);
                else if (nodeIter instanceof Sink) nodetypes.add(NodeType.Sink);
                else if (nodeIter instanceof Router) nodetypes.add(NodeType.Router);
                else if (nodeIter instanceof Delay) nodetypes.add(NodeType.Delay);
                else if (nodeIter instanceof Fork) nodetypes.add(NodeType.Fork);
                else if (nodeIter instanceof Queue) {
                    if (((Queue) nodeIter).getNumberOfServers() == Integer.MAX_VALUE) {
                        nodetypes.add(NodeType.Delay);
                    } else {
                        nodetypes.add(NodeType.Queue);
                    }
                } else if (nodeIter instanceof Source) nodetypes.add(NodeType.Source);
                else if (nodeIter instanceof Place) nodetypes.add(NodeType.Place);
                else if (nodeIter instanceof Transition) nodetypes.add(NodeType.Transition);
                else if (nodeIter instanceof Cache) nodetypes.add(NodeType.Cache);
                else throw new Exception("Unknown node type.");
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException("Fatal error in Network.getNodeTypes() call.", e);
        }

        return nodetypes;
    }

    /**
     * Returns the list of all nodes in this network.
     *
     * @return list of nodes including stations and non-station nodes
     */
    public List<Node> getNodes() {
        return this.nodes;
    }

    // ========================================================================
    // SECTION 12: GETTER METHODS - COUNTS
    // Methods for retrieving counts of various network elements
    // ========================================================================

    public int getNumberOfChains() {
        return getStruct(false).nchains;
    }

    public int getNumberOfClasses() {
        return this.jobClasses.size();
    }

    public int getNumberOfOpenClasses() {
        int count = 0;
        for (JobClass jobClass : this.jobClasses) {
            if (jobClass.type == JobClassType.OPEN) {
                count++;
            }
        }
        return count;
    }

    public int getNumberOfClosedClasses() {
        int count = 0;
        for (JobClass jobClass : this.jobClasses) {
            if (jobClass.type == JobClassType.CLOSED) {
                count++;
            }
        }
        return count;
    }

    public Matrix getNumberOfJobs() {
        int K = getNumberOfClasses();
        Matrix njobs = new Matrix(K, 1, K);
        for (int i = 0; i < K; i++) {
            if (jobClasses.get(i).type == JobClassType.OPEN) njobs.set(i, 0, Inf);
            else if (jobClasses.get(i).type == JobClassType.CLOSED)
                njobs.set(i, 0, jobClasses.get(i).getNumberOfJobs());
            else if (jobClasses.get(i).type == JobClassType.DISABLED)
                njobs.set(i, 0, 0);
        }
        return njobs;
    }

    /**
     * Returns the total number of nodes in this network.
     *
     * @return number of nodes
     */
    public int getNumberOfNodes() {
        return this.nodes.size();
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

    /**
     * Returns the total number of stations in this network.
     *
     * @return number of service stations
     */
    public int getNumberOfStations() {
        return this.stations.size();
    }

    // ========================================================================
    // SECTION 13: GETTER METHODS - PROCESS AND PRODUCT FORM
    // Methods for retrieving process types and product form parameters
    // ========================================================================

    public ProcessType getProcessType(Distribution distr) {

        if (distr instanceof Erlang) {
            return ProcessType.ERLANG;
        } else if (distr instanceof Exp) {
            return ProcessType.EXP;
        } else if (distr instanceof HyperExp) {
            return ProcessType.HYPEREXP;
        } else if (distr instanceof APH) {
            return ProcessType.APH;
        } else if (distr instanceof Geometric) {
            return ProcessType.GEOMETRIC;
        } else if (distr instanceof Bernoulli) {
            return ProcessType.BERNOULLI;
        } else if (distr instanceof Cox2) {
            return ProcessType.COX2;
        } else if (distr instanceof PH) {
            return ProcessType.PH;
        } else if (distr instanceof MMPP2) {
            return ProcessType.MMPP2;
        } else if (distr instanceof Lognormal) {
            return ProcessType.LOGNORMAL;
        } else if (distr instanceof Pareto) {
            return ProcessType.PARETO;
        } else if (distr instanceof Weibull) {
            return ProcessType.WEIBULL;
        } else if (distr instanceof DiscreteUniform) {
            return ProcessType.DUNIFORM;
        } else if (distr instanceof Uniform) {
            return ProcessType.UNIFORM;
        } else if (distr instanceof Gamma) {
            return ProcessType.GAMMA;
        } else if (distr instanceof Det) {
            return ProcessType.DET;
        } else if (distr instanceof Coxian) {
            return ProcessType.COXIAN;
        } else if (distr instanceof Poisson) {
            return ProcessType.POISSON;
        } else if (distr instanceof Replayer || distr instanceof Trace) {
            return ProcessType.REPLAYER;
        } else if (distr instanceof Binomial) {
            return ProcessType.BINOMIAL;
        } else if (distr instanceof BMAP) {
            return ProcessType.BMAP;
        } else if (distr instanceof MarkedMAP) {
            return ProcessType.MMAP;
        } else if (distr instanceof MAP) {
            return ProcessType.MAP;
        } else if (distr instanceof ME) {
            return ProcessType.ME;
        } else if (distr instanceof RAP) {
            return ProcessType.RAP;
        } else if (distr instanceof Immediate) {
            return ProcessType.IMMEDIATE;
        } else {
            return ProcessType.DISABLED;
        }
    }

    public Ret.snGetProductFormParams getProductFormChainParameters() {
        return snGetProductFormChainParams(getStruct());
    }

    public Ret.snGetProductFormParams getProductFormParameters() {
        return snGetProductFormParams(getStruct());
    }

    public Matrix getReferenceClasses() {
        int K = this.jobClasses.size();
        Matrix refclass = new Matrix(K, 1);
        for (int i = 0; i < K; i++) {
            if (this.jobClasses.get(i).isReferenceClass()) refclass.set(i, 0, 1.0);
        }
        return refclass;
    }

    public Matrix getReferenceStations() {
        int K = getNumberOfClasses();
        Matrix refstat = new Matrix(K, 1, K);

        for (int i = 0; i < K; i++) {
            if (jobClasses.get(i).type == JobClassType.OPEN) {
                refstat.set(i, 0, getIndexSourceStation());
            } else {
                ClosedClass cc = (ClosedClass) jobClasses.get(i);
                refstat.set(i, 0, getStationIndex(cc.getReferenceStation()));
            }
        }


        return refstat;
    }

    public List<Region> getRegions() {
        return regions;
    }

    // ========================================================================
    // SECTION 14: ROUTING MATRIX MANAGEMENT
    // Methods for managing and retrieving routing matrices for stations
    // ========================================================================

    public routingMatrixReturn getRoutingMatrix(Matrix arvRates, int returnVal) {

        int idxSource, idxSink, I, K;
        List<Integer> idxOpenClasses;
        boolean hasOpen;
        Matrix conn, NK;

        // Validate that cached struct values match actual network state
        // If there's a mismatch, we cannot use cached values
        boolean useCache = this.hasStruct;
        if (useCache) {
            // Check for consistency between cached values and actual network
            if (this.sn.nnodes != this.nodes.size() || this.sn.nclasses != this.jobClasses.size()) {
                // Cached values are stale, don't use them
                useCache = false;
            }
        }

        if (useCache) {
            idxSource = this.sn.nodetype.indexOf(NodeType.Source);
            idxSink = this.sn.nodetype.indexOf(NodeType.Sink);
            idxOpenClasses = new ArrayList<Integer>();
            for (int col = 0; col < this.sn.njobs.getNumCols(); col++) {
                if (isInf((this.sn.njobs.get(0, col)))) idxOpenClasses.add(col);
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

            if (this.sn == null) this.sn = new NetworkStruct();
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
        for (int ind = 0; ind < I; ind++) {
            Node node = this.nodes.get(ind);
            if (node.getOutput() instanceof Forker) {
                for (int jnd = 0; jnd < I; jnd++) {
                    for (int k = 0; k < K; k++) {
                        if (conn.get(ind, jnd) > 0) {
                            JobClass jobclass = this.jobClasses.get(k);
                            List<OutputStrategy> outputStrategy_k = node.getOutput().getOutputStrategyByClass(jobclass);
                            // MATLAB: if outputStrategy{k} has destinations, fanout = length of that
                            // Otherwise fanout = sum(conn(ind,:)) - number of outgoing connections
                            // Count only entries with non-null destinations
                            int fanout = (int) conn.sumRows(ind); // Default: number of outgoing connections
                            if (outputStrategy_k != null) {
                                int destCount = 0;
                                for (OutputStrategy os : outputStrategy_k) {
                                    if (os.getDestination() != null) {
                                        destCount++;
                                    }
                                }
                                if (destCount > 0) {
                                    fanout = destCount;
                                }
                            }
                            rtnodes.set(ind * K + k, jnd * K + k, 1.0 / fanout);
                            if (this.sn.routing.get(node).get(jobclass) == RoutingStrategy.PROB) {
                                //Check the number of outgoing links
                                int sum = (int) conn.sumRows(ind);
                                //Fork must have all the output strategy towards all outgoing links.
                                if (outputStrategy_k.size() != sum)
                                    line_warning(mfilename(new Object() {
                                    }), "Fork must have all the output strategy towards all outgoing links.");
                                //Fork must have 1.0 routing probability towards all outgoing links.
                                for (OutputStrategy ops : outputStrategy_k) {
                                    if (ops.getProbability() != 1.0)
                                        line_warning(mfilename(new Object() {
                                        }), "Fork must have 1.0 routing probability towards all outgoing links.");
                                }
                            }
                        }
                    }
                }
            } else {
                boolean isSink_i = (ind == idxSink);
                boolean isSource_i = (ind == idxSource);
                for (int k = 0; k < K; k++) {
                    JobClass jobclass = this.jobClasses.get(k);
                    List<OutputStrategy> outputStrategy_k = node.getOutput().getOutputStrategyByClass(jobclass);
                    // Declare variables outside switch to avoid scope issues
                    double sum;
                    switch (this.sn.routing.get(node).get(jobclass)) {
                        case PROB:
                            if (isInf((NK.get(0, k))) || !isSink_i) {
                                for (OutputStrategy ops : outputStrategy_k) {
                                    // Use this network's index for the destination, not the destination's cached index
                                    // This handles cases where the destination node may belong to a different model
                                    int j = this.getNodeIndex(ops.getDestination());
                                    if (j >= 0 && j < I) {
                                        rtnodes.set(ind * K + k, j * K + k, ops.getProbability());
                                    }
                                }
                            }
                            break;
                        //Not tested the following situation
                        case DISABLED:
                            sum = conn.sumRows(ind);
                            for (int j = 0; j < I; j++) {
                                if (conn.get(ind, j) > 0) rtnodes.set(ind * K + k, j * K + k, 1.0 / sum);
                            }
                            break;
                        case RAND:
                        case RROBIN:
                        case WRROBIN:
                        case JSQ:
                            if (isInf((NK.get(0, k)))) {
                                sum = conn.sumRows(ind);
                                for (int j = 0; j < I; j++) {
                                    if (conn.get(ind, j) > 0) rtnodes.set(ind * K + k, j * K + k, 1.0 / sum);
                                }
                            } else if (!isSource_i && !isSink_i) {
                                Matrix connectionClosed = conn.copy();
                                if (idxSink >= 0) { // this is empty set in MATLAB
                                    if (connectionClosed.get(ind, idxSink) > 0) connectionClosed.remove(ind, idxSink);
                                }
                                sum = connectionClosed.sumRows(ind);
                                for (int j = 0; j < I; j++) {
                                    if (connectionClosed.get(ind, j) > 0)
                                        rtnodes.set(ind * K + k, j * K + k, 1.0 / sum);
                                }
                            }
                            break;
                        default:
                            for (int j = 0; j < I; j++) {
                                if (conn.get(ind, j) > 0) rtnodes.set(ind * K + k, j * K + k, GlobalConstants.Zero);
                            }
                    }
                }
            }
        }

        // The second loop corrects the first one at nodes that change the class of the job in the service section.
        for (int i = 0; i < I; i++) {
            Node node = this.nodes.get(i);
            if (node.getServer() instanceof StatelessClassSwitcher) {
                StatelessClassSwitcher classSwitcher = (StatelessClassSwitcher) node.getServer();
                Matrix Pi = new Matrix(K, rtnodes.getNumCols(), 0);
                Matrix.extract(rtnodes, i * K, (i + 1) * K, 0, rtnodes.getNumCols(), Pi, 0, 0);
                Matrix Pcs = new Matrix(K, K);
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
                            if (val != 0) rtnodes.set(i * K + row, jnd * K + col, val);
                            else rtnodes.remove(i * K + row, jnd * K + col);
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
                        for (int j = 0; j < cacheClassSwitcher.hitClass.getNumCols(); j++) {
                            if (cacheClassSwitcher.hitClass.get(0, j) == r) {
                                flagHit = true;
                                break;
                            }
                        }
                        for (int j = 0; j < cacheClassSwitcher.missClass.getNumCols(); j++) {
                            if (cacheClassSwitcher.missClass.get(0, j) == r) {
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
                        for (int j = 0; j < cacheClassSwitcher.hitClass.getNumCols(); j++) {
                            if (cacheClassSwitcher.hitClass.get(0, j) == r) {
                                flagHit = true;
                                break;
                            }
                        }
                        for (int j = 0; j < cacheClassSwitcher.missClass.getNumCols(); j++) {
                            if (cacheClassSwitcher.missClass.get(0, j) == r) {
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
                        if (cacheClassSwitcher.actualHitProb.length() > r && cacheClassSwitcher.hitClass.getNumCols() > r && cacheClassSwitcher.hitClass.get(0, r) != -1) {
                            double ph = cacheClassSwitcher.actualHitProb.get(r);
                            double pm = cacheClassSwitcher.actualMissProb.get(r);

                            int h = (int) cacheClassSwitcher.hitClass.get(0, r);
                            int m = (int) cacheClassSwitcher.missClass.get(0, r);

                            rtnodes.set(i * K + r, i * K + h, ph);
                            rtnodes.set(i * K + r, i * K + m, pm);
                        } else {
                            if (cacheClassSwitcher.hitClass.getNumCols() > r && cacheClassSwitcher.hitClass.get(0, r) != -1) {
                                int h = (int) cacheClassSwitcher.hitClass.get(0, r);
                                int m = (int) cacheClassSwitcher.missClass.get(0, r);
                                rtnodes.set(i * K + r, i * K + h, NaN);
                                rtnodes.set(i * K + r, i * K + m, NaN);
                            }
                        }
                    }
                    for (int jnd = 0; jnd < I; jnd++) {
                        Matrix Pij = new Matrix(K, K);
                        Matrix.extract(Pi, 0, K, jnd * K, (jnd + 1) * K, Pij, 0, 0);
                        for (int r = 0; r < K; r++) {
                            boolean flagHit = false, flagMiss = false;
                            for (int j = 0; j < cacheClassSwitcher.hitClass.getNumCols(); j++) {
                                if (cacheClassSwitcher.hitClass.get(0, j) == r) {
                                    flagHit = true;
                                    break;
                                }
                            }
                            for (int j = 0; j < cacheClassSwitcher.missClass.getNumCols(); j++) {
                                if (cacheClassSwitcher.missClass.get(0, j) == r) {
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

        // Handle self-looping classes
        for (int r = 0; r < K; r++) {
            if (this.jobClasses.get(r) instanceof SelfLoopingClass) {
                SelfLoopingClass slClass = (SelfLoopingClass) this.jobClasses.get(r);
                int slRefNode = slClass.getReferenceStation().getNodeIndex();

                for (int ind = 0; ind < I; ind++) {            // source
                    for (int jnd = 0; jnd < I; jnd++) {        // destination
                        for (int s = 0; s < K; s++) {          // class at destination
                            rtnodes.set(ind * K + r, jnd * K + s, 0.0); // reset
                        }
                    }
                    // Route the self-looping class only to its reference station
                    rtnodes.set(slRefNode * K + r, slRefNode * K + r, 1.0);
                }
            }
        }

        // ignore all chains containing a Pnodes column that sums to 0, since these are classes that cannot arrive to the node unless this column belongs to the source
        Matrix sumRtnodesCols = rtnodes.sumCols();
        Set<Integer> colsToIgnore = new HashSet<Integer>();
        for (int col = 0; col < sumRtnodesCols.getNumCols(); col++) {
            if (sumRtnodesCols.get(col) == 0) colsToIgnore.add(col);
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
            /* rtnodes+rtnodes' makes the matrix undirected, it helps grouping
               transient states with all the recurrent states their can reach
               into a single "chain" of the LINE model
            */
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
                    if (chains.get(j, i) == 1) rows.add(j);
                }

                if (rows.size() > 1) {
                    int row = rows.get(0);
                    for (int j = 1; j < rows.size(); j++) {
                        //chains(rows(1),:) = chains(row(1),:) | chains(r,:);
                        for (int k = 0; k < chains.getNumCols(); k++) {
                            if (chains.get(row, k) == 1 || chains.get(rows.get(j), k) == 1) chains.set(row, k, 1);
                            else chains.set(row, k, 0);
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
        Matrix rt = dtmc_stochcomp(rtnodes, statefulNodeClasses);
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

    public RoutingStrategy getRoutingStrategyFromNodeAndClassPair(Node node, JobClass c) {
        //Another approach is to use the last routing strategy
        if (node.getOutput().getOutputStrategyByClass(c).isEmpty()) {
            if (node instanceof Sink || node instanceof Source) {
                return RoutingStrategy.DISABLED;
            } else {
                return RoutingStrategy.RAND;
            }
        }
        RoutingStrategy res = null;
        try {
            for (OutputStrategy outputStrategy : node.getOutputStrategies()) {
                // Use index comparison for JobClass to handle dynamically created classes
                if (outputStrategy.getJobClass().getIndex() == c.getIndex()) {
                    if (res == null)
                        res = outputStrategy.getRoutingStrategy();
                    else if (!res.equals(outputStrategy.getRoutingStrategy())) {
                        line_error(mfilename(new Object() {
                        }), "Inconsistent routing strategy.");
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            //System.exit(1);
            throw new RuntimeException("Fatal error in Network.getRoutingStrategyFromNodeAndClassPair() call.", e);
        }

        if (res != null) return res;
        else if (c instanceof OpenClass) return RoutingStrategy.RAND;
        else return RoutingStrategy.DISABLED;
    }

    // ========================================================================
    // SECTION 15: GETTER METHODS - SOURCES AND SINKS
    // Methods for retrieving source and sink nodes
    // ========================================================================

    public Sink getSink() {
        int index = this.getIndexSinkNode();
        if (index == -1) {
            line_error(mfilename(new Object() {
            }), "The given model does not have a sink");
            return null;
        }
        return (Sink) this.nodes.get(index);
    }

    /**
     * Returns the dimensions of this network as [nodes, classes].
     *
     * @return array containing [number of nodes, number of job classes]
     */
    public int[] getSize() {
        int[] outInt = new int[2];
        outInt[0] = this.getNumberOfNodes();
        outInt[1] = this.getNumberOfClasses();
        return outInt;
    }

    public Source getSource() {
        int index = this.getIndexSourceNode();
        if (index == -1) {
            line_error(mfilename(new Object() {
            }), "The given model does not have a source or this is declared after OpenClass.");
            return null;
        }
        return (Source) this.nodes.get(index);
    }

    // Get initial state
    public State getState() {
        if (!this.hasInitState()) {
            // Prevent infinite recursion during state initialization
            if (this.initializingState) {
                return null;
            }
            this.initializingState = true;
            try {
                this.initDefault();
            } finally {
                this.initializingState = false;
            }
        }

        // Ensure state maps are initialized
        if (this.sn != null) {
            if (this.sn.state == null) {
                this.sn.state = new HashMap<StatefulNode, Matrix>();
            }
            if (this.sn.stateprior == null) {
                this.sn.stateprior = new HashMap<StatefulNode, Matrix>();
            }
            if (this.sn.space == null) {
                this.sn.space = new HashMap<StatefulNode, Matrix>();
            }
        }

        for (int i = 0; i < this.getNumberOfNodes(); i++) {
            if (this.nodes.get(i).isStateful()) {
                Node node_i = this.nodes.get(i);
                Matrix initialState = ((StatefulNode) node_i).getState();
                Matrix priorInitialState = ((StatefulNode) node_i).getStatePrior();
                Matrix initialStateSpace = ((StatefulNode) node_i).getStateSpace();
                
                this.sn.state.put((StatefulNode) node_i, initialState);
                this.sn.stateprior.put((StatefulNode) node_i, priorInitialState);
                this.sn.space.put((StatefulNode) node_i, initialStateSpace);
            }
        }

        return new State(this.sn.state, this.sn.stateprior, this.sn.space);
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

    public int getStatefulNodeIndex(Node node) {
        if (!(node instanceof StatefulNode)) return -1;

        int outIdx = 0;
        for (Node nodeIter : this.nodes) {
            if (nodeIter == node) {
                return outIdx;
            } else if (nodeIter instanceof StatefulNode) {
                outIdx++;
            }
        }

        // Node not found - provide debugging information
        StringBuilder debug = new StringBuilder();
        debug.append("Node '").append(node.getName()).append("' (").append(node.getClass().getSimpleName()).append(") not found in network nodes list.\n");
        debug.append("Network contains ").append(this.nodes.size()).append(" nodes:\n");
        for (int i = 0; i < this.nodes.size(); i++) {
            Node n = this.nodes.get(i);
            debug.append("  [").append(i).append("] ").append(n.getName()).append(" (").append(n.getClass().getSimpleName()).append(")");
            if (n instanceof StatefulNode) {
                debug.append(" [StatefulNode]");
            }
            if (n.getName().equals(node.getName())) {
                debug.append(" [SAME NAME, DIFFERENT REFERENCE]");
            }
            debug.append("\n");
        }
        debug.append("This suggests a reference mismatch between the node and the network's node list.");
        
        return -1;
    }

    public int getStatefulNodeIndex(String name) {
        int ret = -1;
        List<String> statefulNodeNames = getStatefulNodeNames();
        for (int i = 0; i < this.getNumberOfStatefulNodes(); i++) {
            if (statefulNodeNames.get(i).equals(name)) {
                ret = i;
                break;
            }
        }
        return ret;
    }

    public List<String> getStatefulNodeNames() {
        List<String> statefulNodeNames = new ArrayList<>();
        for (int i = 0; i < this.getNumberOfNodes(); i++) {
            if (nodes.get(i).isStateful()) {
                statefulNodeNames.add(nodes.get(i).getName());
            }
        }
        return statefulNodeNames;
    }

    public List<StatefulNode> getStatefulNodes() {
        if (stateful.isEmpty()) {
            List<StatefulNode> statefulNodes = new ArrayList<>();
            for (int i = 0; i < this.getNumberOfNodes(); i++) {
                if (nodes.get(i).isStateful()) {
                    statefulNodes.add((StatefulNode) nodes.get(i));
                }
            }
            stateful = statefulNodes;
        }
        return stateful;
    }

    public Matrix getStatefulServers() {
        int I = getStatefulNodes().size();
        Matrix numservers = new Matrix(I, 1, I);
        for (int i = 0; i < I; i++) {
            if (getStatefulNodes().get(i).getNumberOfServers() == Integer.MAX_VALUE)
                numservers.set(i, 0, Inf);
            else numservers.set(i, 0, stations.get(i).getNumberOfServers());
        }

        return numservers;
    }

    public Station getStationByIndex(int index) {
        return this.stations.get(index);
    }

    public Station getStationByName(String name) {
        for (Station stat : this.stations) {
            if (stat.getName().equals(name)) {
                return stat;
            }
        }
        return null;
    }

    public Node getStationFromIndex(int inIdx) {
        return this.stations.get(inIdx);
    }

    public int getStationIndex(Node node) {
        return this.stations.indexOf(node);
    }

    public int getStationIndex(String name) {
        int ret = -1;
        List<String> stationNodeNames = getStationNames();
        for (int i = 0; i < this.getNumberOfStatefulNodes(); i++) {
            if (stationNodeNames.get(i).equals(name)) {
                ret = i;
                break;
            }
        }
        return ret;
    }

    public List<Integer> getStationIndexes(int index) {
        List<Integer> statIndexes = new ArrayList<>();

        for (int i = 0; i < this.stations.size(); i++) {
            statIndexes.add(this.stations.get(i).getStationIdx());
        }
        return statIndexes;
    }

    public List<String> getStationNames() {
        List<String> stationNodeNames = new ArrayList<>();
        for (int i = 0; i < this.getNumberOfNodes(); i++) {
            if (nodes.get(i) instanceof Station) {
                stationNodeNames.add(nodes.get(i).getName());
            }
        }
        return stationNodeNames;
    }

    public Map<Station, SchedStrategy> getStationScheduling() {
        Map<Station, SchedStrategy> res = new HashMap<Station, SchedStrategy>();
        for (Station station : this.stations) {
            if (station.getNumberOfServers() == Integer.MAX_VALUE) {
                res.put(station, SchedStrategy.INF);
            } else {
                if (station instanceof Source) {
                    res.put(station, SchedStrategy.EXT);
                } else if (station instanceof Station) {
                    res.put(station, station.getSchedStrategy());
                }
            }
        }
        return res;
    }

    public Matrix getStationServers() {
        int I = stations.size();
        Matrix numservers = new Matrix(I, 1, I);
        for (int i = 0; i < I; i++) {
            if (stations.get(i).getNumberOfServers() == Integer.MAX_VALUE)
                numservers.set(i, 0, Inf);
            else numservers.set(i, 0, stations.get(i).getNumberOfServers());
        }

        return numservers;
    }

    /**
     * Returns the list of stations in this network.
     * Stations are nodes that can provide service to jobs.
     *
     * @return list of service stations
     */
    public List<Station> getStations() {
        return this.stations;
    }

    public NetworkStruct getStruct() {
        return this.getStruct(true);
    }

    public void setStruct(NetworkStruct sn) {
        this.sn = sn;
    }

    public NetworkStruct getStruct(boolean wantInitialState) {
        if (!this.hasStruct) {
            refreshStruct(true);
            if (GlobalConstants.getVerbose() == VerboseLevel.DEBUG) {
                this.sn.print();
            }
        }

        if (wantInitialState) {
            State state = getState();
            // If getState returns null due to recursion protection, skip state initialization for now
            if (state == null && this.initializingState) {
                // State initialization is in progress, skip for now
            }
        }

        return this.sn;
    }

    public SolverTranHandles getTranHandles() {
        int M = this.stations.size();
        int K = this.jobClasses.size();

        Matrix isSource = new Matrix(M, 1);
        Matrix isSink = new Matrix(M, 1);
        Matrix hasServiceTunnel = new Matrix(M, 1);
        Matrix isServiceDefined = Matrix.ones(M, K);

        for (int i = 0; i < M; i++) {
            if (this.stations.get(i) instanceof Source) isSource.set(i, 0, 1);
            if (((Node) this.stations.get(i)) instanceof Sink) isSink.set(i, 0, 1);

            if (this.stations.get(i).getServer() instanceof ServiceTunnel) hasServiceTunnel.set(i, 0, 1);
            else {
                for (int r = 0; r < K; r++) {
//                    if (!this.stations.get(i).getServer().containsJobClass(this.jobClasses.get(r)))
                    if (this.stations.get(i).getServer().getServiceDistribution(this.jobClasses.get(r)).isDisabled())
                        isServiceDefined.remove(i, r);
                }
            }
        }

        //Calculate Qt
        AvgHandle Qt = new AvgHandle();
        for (int i = 0; i < M; i++) {
            Map<JobClass, Metric> map = new HashMap<>();
            for (int r = 0; r < K; r++) {
                Metric Qir = new Metric();
                Qir.type = "Number of Customers";
                Qir.jobClass = this.jobClasses.get(r);
                Qir.station = this.stations.get(i);
                if (isSource.get(i, 0) > 0) Qir.isDisabled = true;
                else if (isSink.get(i, 0) > 0) Qir.isDisabled = true;
                else Qir.isDisabled = hasServiceTunnel.get(i, 0) == 0 && isServiceDefined.get(i, r) == 0;
                Qt.put(this.stations.get(i), this.jobClasses.get(r), Qir);
            }
        }

        //Calculate Ut
        AvgHandle Ut = new AvgHandle();
        for (int i = 0; i < M; i++) {
            Map<JobClass, Metric> map = new HashMap<>();
            for (int r = 0; r < K; r++) {
                Metric Uir = new Metric();
                Uir.type = "Utilization";
                Uir.jobClass = this.jobClasses.get(r);
                Uir.station = this.stations.get(i);
                if (isSource.get(i, 0) > 0) Uir.isDisabled = true;
                else if (isSink.get(i, 0) > 0) Uir.isDisabled = true;
                else if (this.stations.get(i) instanceof Join) Uir.isDisabled = true;
                else Uir.isDisabled = hasServiceTunnel.get(i, 0) == 0 && isServiceDefined.get(i, r) == 0;
                Ut.put(this.stations.get(i), this.jobClasses.get(r), Uir);
            }
        }

        //Calculate Tt
        AvgHandle Tt = new AvgHandle();
        for (int i = 0; i < M; i++) {
            Map<JobClass, Metric> map = new HashMap<>();
            for (int r = 0; r < K; r++) {
                Metric Tir = new Metric();
                Tir.type = "Throughput";
                Tir.jobClass = this.jobClasses.get(r);
                Tir.station = this.stations.get(i);
                Tir.isDisabled = hasServiceTunnel.get(i, 0) == 0 && isServiceDefined.get(i, r) == 0;
                Tt.put(this.stations.get(i), this.jobClasses.get(r), Tir);
            }
        }

        return new SolverTranHandles(Qt, Ut, Tt);
    }

    public AvgHandle getTranQLenHandles() {
        return this.getTranHandles().getTranQLenHandles();
    }

    public AvgHandle getTranTputHandles() {
        return this.getTranHandles().getTranTputHandles();
    }

    public AvgHandle getTranUtilHandles() {
        return this.getTranHandles().getTranUtilHandles();
    }

    /**
     * Returns the language features used by the given network
     *
     * @return - the language features used by the given network
     */
    public FeatureSet getUsedLangFeatures() {
        usedFeatures = new FeatureSet();
        if (!this.getIndexClosedClasses().isEmpty()) {
            setUsedLangFeature("ClosedClass");
        }
        if (!this.getIndexOpenClasses().isEmpty()) {
            setUsedLangFeature("OpenClass");
        }
        // Check for self-looping classes
        for (JobClass jc : this.jobClasses) {
            if (jc instanceof SelfLoopingClass) {
                setUsedLangFeature("SelfLoopingClass");
                break;
            }
        }
        for (int i = 0; i < getNumberOfNodes(); i++) {
            for (int r = 0; r < getNumberOfClasses(); r++) {
                Node n = this.nodes.get(i);
                if (n instanceof Queue || n instanceof Delay) {
                    ServiceBinding serviceProcess = n.getServer().getServiceProcess(this.getClassByIndex(r));
                    if (serviceProcess != null) {
                        if (!(serviceProcess.getDistribution() instanceof Disabled) && !(serviceProcess.getDistribution() instanceof Immediate)) {
                            setUsedLangFeature(serviceProcess.getDistribution().getName());
                        }
                        String sched = "";
                        if (n instanceof Delay) {
                            setUsedLangFeature("Delay");
                            sched = SchedStrategy.toFeature(((Delay) n).getSchedStrategy());
                        } else {
                            setUsedLangFeature("Queue");
                            sched = SchedStrategy.toFeature(((Queue) n).getSchedStrategy());
                        }
                        if (sched.length() > 0) {
                            setUsedLangFeature(sched);
                        }
                        if (r < n.getOutput().getOutputStrategies().size()) {
                            String routing = RoutingStrategy.toFeature(n.getOutput().getOutputStrategies().get(r).getRoutingStrategy());
                            if (routing.length() > 0) {
                                setUsedLangFeature(routing);
                            }
                        }
                    }
                } else if (n instanceof Router) {
                    if (r < n.getOutput().getOutputStrategies().size()) {
                        String routing = RoutingStrategy.toFeature(n.getOutput().getOutputStrategies().get(r).getRoutingStrategy());
                        if (routing.length() > 0) {
                            setUsedLangFeature(routing);
                        }
                    }
                } else if (n instanceof Source) {
                    Distribution serviceProcess = ((Source) n).getArrivalProcess(this.getClassByIndex(r));
                    if (!(serviceProcess instanceof Disabled) && !(serviceProcess instanceof Immediate)) {
                        setUsedLangFeature(serviceProcess.getName());
                    }
                    setUsedLangFeature("Source");
                } else if (n instanceof ClassSwitch) {
                    setUsedLangFeature("StatelessClassSwitcher");
                    setUsedLangFeature("ClassSwitch");
                } else if (n instanceof Fork) {
                    setUsedLangFeature("Fork");
                    setUsedLangFeature("Forker");
                } else if (n instanceof Join) {
                    setUsedLangFeature("Join");
                    setUsedLangFeature("Joiner");
                } else if (n instanceof Sink) {
                    setUsedLangFeature("JobSink");
                } else if (n instanceof Cache) {
                    setUsedLangFeature("CacheClassSwitcher");
                    setUsedLangFeature("Cache");
                } else if (n instanceof Transition) {
                    setUsedLangFeature("Transition");
                    setUsedLangFeature("Enabling");
                    setUsedLangFeature("Timing");
                    setUsedLangFeature("Firing");
                } else if (n instanceof Place) {
                    setUsedLangFeature("Storage");
                    setUsedLangFeature("Linkage");
                    setUsedLangFeature("Place");
                }
            }
        }
        // Check for load-dependent service (LLD, LJD, or LJCD)
        for (Station station : this.stations) {
            if ((station.getLimitedLoadDependence() != null && !station.getLimitedLoadDependence().isEmpty()) ||
                (station.getLimitedJointDependence() != null && !station.getLimitedJointDependence().isEmpty()) ||
                (station.getLimitedJointClassDependence() != null && !station.getLimitedJointClassDependence().isEmpty())) {
                setUsedLangFeature("LoadDependence");
                break;
            }
        }
        return usedFeatures;
    }

    // ========================================================================
    // SECTION 16: QUERY METHODS (HAS/IS)
    // Methods for querying network properties and characteristics
    // ========================================================================

    public boolean hasClassSwitching() {
        return snHasClassSwitching(getStruct());
    }

    /**
     * Checks if this network has any job classes defined.
     *
     * @return true if job classes exist, false otherwise
     */
    public boolean hasClasses() {
        return !this.jobClasses.isEmpty();
    }

    /**
     * Checks if this network contains any closed job classes.
     * Closed classes have fixed populations with no external arrivals.
     *
     * @return true if closed classes exist, false otherwise
     */
    public boolean hasClosedClasses() {
        for (JobClass temp : this.jobClasses) {
            if (temp instanceof ClosedClass) {
                return true;
            }
        }
        return false;
    }

    public boolean hasDPS() {
        return snHasDPS(getStruct());
    }

    public boolean hasDPSPrio() {
        return snHasDPSPRIO(getStruct());
    }

    public boolean hasFCFS() {
        return snHasFCFS(getStruct());
    }

    public boolean hasFork() {
        for (NodeType type : this.getNodeTypes()) {
            if (type == NodeType.Fork) {
                return true;
            }
        }
        return false;
    }

    public boolean hasGPS() {
        return snHasGPS(getStruct());
    }

    public boolean hasGPSPrio() {
        return snHasGPSPRIO(getStruct());
    }

    public boolean hasHOL() {
        return snHasHOL(getStruct());
    }

    public boolean hasHomogeneousScheduling(SchedStrategy strategy) {
        return snHasHomogeneousScheduling(getStruct(), strategy);
    }

    public boolean hasINF() {
        return snHasINF(getStruct());
    }

    public boolean hasInitState() {
        boolean output = true;
        if (!this.hasState) { // check if all stations are initialized
            for (int i = 0; i < this.getNumberOfNodes(); i++) {
                if (this.nodes.get(i) instanceof StatefulNode) {
                    if (((StatefulNode) this.nodes.get(i)).getState().isEmpty()) {
                        output = false;
                        break;
                    }
                }
            }
        }
        return output;
    }

    public boolean hasJoin() {
        for (NodeType type : this.getNodeTypes()) {
            if (type == NodeType.Join) {
                return true;
            }
        }
        return false;
    }

    public boolean hasLCFS() {
        return snHasLCFS(getStruct());
    }

    public boolean hasLCFSPR() {
        return snHasLCFSPR(getStruct());
    }

    public boolean hasLEPT() {
        return snHasLEPT(getStruct());
    }

    public boolean hasLJF() {
        return snHasLJF(getStruct());
    }

    public boolean hasMultiChain() {
        return snHasMultiChain(getStruct());
    }

    public boolean hasMultiClass() {
        return snHasMultiClass(getStruct());
    }

    public boolean hasMultiClassFCFS() {
        return snHasMultiClassFCFS(getStruct());
    }

    public boolean hasMultiClassHeterFCFS() {
        return snHasMultiClassHeterFCFS(getStruct());
    }

    public boolean hasMultiServer() {
        return snHasMultiServer(getStruct());
    }

    /**
     * Checks if this network contains any open job classes.
     * Open classes have external arrivals and departures.
     *
     * @return true if open classes exist, false otherwise
     */
    public boolean hasOpenClasses() {
        for (JobClass temp : this.jobClasses) {
            if (temp instanceof OpenClass) {
                return true;
            }
        }

        return false;
    }

    public boolean hasPS() {
        return snHasPS(getStruct());
    }

    public boolean hasPSPrio() {
        return snHasPSPRIO(getStruct());
    }

    /**
     * Checks if this network has a product-form solution.
     * Product-form networks can be solved efficiently using MVA methods.
     *
     * @return true if the network has product-form, false otherwise
     */
    public boolean hasProductFormSolution() {
        return snHasProductForm(this.getStruct(false));
    }

    public boolean hasSEPT() {
        return snHasSEPT(getStruct());
    }

    public boolean hasSIRO() {
        return snHasSIRO(getStruct());
    }

    public boolean hasSJF() {
        return snHasSJF(getStruct());
    }

    public boolean hasSingleChain() {
        return snHasSingleChain(getStruct());
    }

    public boolean hasSingleClass() {
        return snHasSingleClass(getStruct());
    }

    // ========================================================================
    // SECTION 17: INITIALIZATION METHODS
    // Methods for initializing network parameters and structures
    // ========================================================================

    public void initDefault() {
        // TODO: is it necessary to have a version where, per LINE, nodes can be passed in as a parameter?
        // open classes are empty
        // closed classes are initialized at reference station
        // running jobs are allocated in class id order until all servers are busy

        NetworkStruct sn = this.getStruct(false);
        int R = sn.nclasses;
        Matrix N = sn.njobs.transpose();

        for (int ind = 0; ind < this.getNumberOfNodes(); ind++) {
            // Check if user has already set a state for this node - if so, preserve it
            if (this.nodes.get(ind).isStateful()) {
                Matrix existingState = ((StatefulNode) this.nodes.get(ind)).getState();
                if (existingState != null && !existingState.isEmpty()) {
                    // User pre-loaded a state, skip default initialization for this node
                    continue;
                }
            }

            Matrix state_i = new Matrix(0, 0);
            if (sn.isstation.get(ind, 0) == 1) {
                Matrix n0 = new Matrix(1, N.length());
                n0.zero();
                Matrix s0 = new Matrix(1, N.length());
                s0.zero();
                double s = sn.nservers.get((int) sn.nodeToStation.get(ind), 0); // allocate

                for (int r = 0; r < N.getNumRows(); r++) {
                    if (isFinite(N.get(r, 0))) { // for all closed classes
                        if (sn.nodeToStation.get(ind) == sn.refstat.get(r, 0)) {
                            n0.set(0, r, N.get(r, 0));
                        }
                    }
                    s0.set(0, r, FastMath.min(n0.get(0, r), s));
                    s -= s0.get(0, r);
                }

                switch (sn.nodetype.get(ind)) {
                    case Cache:
                        state_i = FromMarginal.fromMarginalAndStarted(sn, ind, n0, s0);
                        int nVarsInt = (int) sn.nvars.get(ind, 2 * R);
                        Matrix newState_i = new Matrix(1, state_i.getNumCols() + nVarsInt);
                        for (int p = 0; p < state_i.length(); p++) {
                            newState_i.set(0, p, state_i.get(0, p));
                        }
                        int addition = 0;
                        for (int p = state_i.length(); p < newState_i.length() - 1; p++) {
                            newState_i.set(0, p, addition);
                            addition++;
                        }
                        state_i = newState_i.copy();
                        break;
                    case Place: {
                        StatefulNode placeNode = (StatefulNode) this.nodes.get(ind);
                        if (placeNode.getState() != null && placeNode.getState().elementSum() > 0) {
                            // user pre-loaded something, keep it
                            state_i = placeNode.getState().copy();
                        } else {
                            state_i = new Matrix(1, this.getNumberOfClasses());
                            for (int r = 0; r < sn.nclasses; r++) {
                                if (sn.refstat.get(r, 0) == sn.nodeToStation.get(ind)) {
                                    state_i.set(0, r, sn.njobs.get(r, 0));
                                } else {
                                    state_i.set(0, r, 0);
                                }
                            }
                        }
                        break;
                    }
                    default:
                        state_i = FromMarginal.fromMarginalAndStarted(sn, ind, n0, s0);
                        if (sn.isstation.get(ind, 0) == 1) {
                            for (int r = 0; r < sn.nclasses; r++) {
                                if (sn.procid.get(sn.nodes.get(ind)).get(sn.jobclasses.get(r)) == ProcessType.MAP) {
                                    Matrix one = new Matrix(1, 1, 1);
                                    one.set(0, 0, 1);
                                    state_i = Matrix.cartesian(state_i, one);
                                }
                            }
                        }
                }

                for (int r = 0; r < sn.nclasses; r++) {
                    if ((sn.routing.get(sn.nodes.get(ind)).get(sn.jobclasses.get(r)) == RoutingStrategy.RROBIN) || (sn.routing.get(sn.nodes.get(ind)).get(sn.jobclasses.get(r)) == RoutingStrategy.WRROBIN)) {
                        // Start from first connected queue
                        for (int p = 0; p < sn.connmatrix.getNumCols(); p++) {
                            if (sn.connmatrix.get(ind, p) == 1) {
                                state_i = state_i.concatCols(Matrix.singleton(p));
                                break;
                            }
                        }
                    }
                }

                if (state_i.isEmpty()) {
                    if (GlobalConstants.Verbose == VerboseLevel.DEBUG) {
                        line_warning(mfilename(new Object() {
                        }), "Default initialisation failed on station " + ind);
                    }
                }
            } else if (this.nodes.get(ind).isStateful()) {
                if (this.nodes.get(ind) instanceof Cache) {
                    Cache cacheNode = (Cache) this.nodes.get(ind);
                    state_i = new Matrix(1, this.getNumberOfClasses() + (int) cacheNode.getItemLevelCap().elementSum());
                    int ctr = 0;
                    for (int idx = this.getNumberOfClasses(); idx < state_i.getNumCols(); idx++) {
                        state_i.set(idx, ++ctr);
                    }
                } else if (this.nodes.get(ind) instanceof Router) {
                    state_i = Matrix.zeros(1, this.getNumberOfClasses());
                } else if (this.nodes.get(ind) instanceof Transition) {
                    Transition tr = (Transition) this.nodes.get(ind);
                    // local buffer: [disabled-servers, zeros(firing-phases), fired-servers]
                    state_i = ((TransitionNodeParam) sn.nodeparam.get(tr)).nmodeservers.copy();
                    for (int k = 0; k < ((TransitionNodeParam) sn.nodeparam.get(tr)).nmodes; k++) {
                        state_i.concatCols(Matrix.zeros(1, (int) ((TransitionNodeParam) sn.nodeparam.get(tr)).firingphases.get(k)));
                    }
                    state_i.concatCols(Matrix.zeros(1, ((TransitionNodeParam) sn.nodeparam.get(tr)).nmodes));
                    // middle region (firing phases) already all zeros
                    // last nmodes positions stay zero as well
                    break;
                } else {
                    state_i = new Matrix(0, 0);
                }
            }

            if (this.nodes.get(ind).isStateful()) {
                StatefulNode node_i = (StatefulNode) this.nodes.get(ind);
                if (state_i.getNumRows() == 1) {
                    node_i.setState(state_i.getRow(0));
                    node_i.setStateSpace(state_i);
                    node_i.setStatePrior(Matrix.singleton(1));
                } else if (state_i.getNumRows() > 1) {
                    Matrix prior_state_i = new Matrix(state_i.getNumRows(), 1);
                    prior_state_i.zero();
                    prior_state_i.set(0, 0, 1.0);
                    node_i.setStateSpace(state_i);
                    node_i.setStatePrior(prior_state_i);
                    node_i.setState(state_i.getRow(0));
                } else {
                    node_i.setState(new Matrix(0, 0));
                    node_i.setStateSpace(new Matrix(0, 0));
                    node_i.setStatePrior(new Matrix(0, 0));
                }
            }
        }

        // Fix: Remove circular dependency - isStateValid() calls getStruct() which calls getState()
        // which can lead to infinite recursion during initialization
        // TODO: Implement a lighter-weight state validation that doesn't require full NetworkStruct
        this.hasState = true;
    }

    public void initFromAvgQLen(Matrix AvgQLen) {
        int rows = AvgQLen.getNumRows();
        int cols = AvgQLen.getNumCols();
        Matrix n = new Matrix(rows, cols);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                n.set(i, j, FastMath.round(AvgQLen.get(i, j)));
            }
        }

        Matrix njobs = n.sumCols();

        for (int r = 0; r < cols; r++) {
            if (njobs.get(r) > AvgQLen.sumCols(r)) {
                int i = njobs.getColMax(r);
                n.set(i, r, n.get(i, r) - 1);
            }
        }

        try {
            this.initFromMarginal(n);
        } catch (Exception e) {
            this.initDefault();
        }
    }

    public void initFromAvgTable(NetworkAvgTable nt) {
        Matrix q = new Matrix(nt.getQLen());
        q.reshape(this.getNumberOfClasses(), this.getNumberOfStations());
        q.transpose();
        this.initFromAvgQLen(q);
    }

    public void initFromMarginal(Matrix n) {

        if (!hasStruct) {
            refreshStruct(true);
        }
        getStruct(true);

        // Convert station-based n to node-based n if necessary
        if (this.getNumberOfStations() < this.getNumberOfNodes()) {
            if (n.getNumRows() == this.getNumberOfStations()) {
                Matrix nnodes = new Matrix(sn.nnodes, n.getNumCols());
                nnodes.zero();
                for (int ist = 0; ist < sn.nstations; ist++) {
                    int ind = (int) sn.stationToNode.get(ist);
                    for (int c = 0; c < n.getNumCols(); c++) {
                        nnodes.set(ind, c, n.get(ist, c));
                    }
                }
                n = nnodes;
            } else if (n.getNumRows() == this.getNumberOfNodes()) {
                // no-op
            } else {
                line_error(mfilename(new Object() {}), 
                    "The supplied matrix of marginal states does not have the correct number of rows. One either one per station or one per node.");
            }
        }

        if (!State.isValid(sn, n, new Matrix(0, 0))) {
            if (GlobalConstants.Verbose == VerboseLevel.DEBUG) {
                line_warning(mfilename(new Object() {
                }), "Initial state not contained in the state spac. Trying to recover.");
            }
            for (int row = 0; row < n.getNumRows(); row++) {
                for (int col = 0; col < n.getNumCols(); col++) {
                    n.set(row, col, FastMath.round(n.get(row, col)));
                }
            }
            if (!State.isValid(sn, n, new Matrix(0, 0))) {
                throw new RuntimeException("Cannot recover from failed state initialization - stopping.");
            }
        }

        for (int i = 0; i < sn.nnodes; i++) {
            Matrix state_i = new Matrix(0, 0);
            if (sn.isstateful.get(i) == 1) {
                int ist = (int) sn.nodeToStation.get(i);
                if (sn.nodetype.get(i) == NodeType.Place) {
                    // Must be single class token
                    state_i = n.sumRows(0, n.getNumCols());
                } else {
                    // Pass the marginal state for this node (row i of n)
                    state_i = FromMarginal.fromMarginal(sn, i, n.getRow(i));
                }
                if (stations.get(ist).getState().isEmpty()) {
                    throw new RuntimeException("Invalid state assignment for a station.");
                }
            }
            if (this.nodes.get(i).isStateful()) {
                StatefulNode node_i = (StatefulNode) this.nodes.get(i);
                if (state_i.getNumRows() == 1) {
                    node_i.setState(state_i.getRow(0));
                    node_i.setStateSpace(state_i);
                    node_i.setStatePrior(Matrix.singleton(1));
                } else if (state_i.getNumRows() > 1) {
                    Matrix prior_state_i = new Matrix(state_i.getNumRows(), 1);
                    prior_state_i.zero();
                    prior_state_i.set(0, 0, 1.0);
                    node_i.setStateSpace(state_i);
                    node_i.setStatePrior(prior_state_i);
                    node_i.setState(state_i.getRow(0));
                } else {
                    node_i.setState(new Matrix(0, 0));
                    node_i.setStateSpace(new Matrix(0, 0));
                    node_i.setStatePrior(new Matrix(0, 0));
                }
            }
        }

        hasState = true;
    }

    public void initFromMarginalAndRunning(Matrix n, Matrix s) {
        NetworkStruct sn = getStruct();
        boolean isvalidn = State.isValid(sn, n, s);
        if (!isvalidn) line_error(mfilename(new Object() {
        }), "Initial state is not valid.");

        for (int ind = 0; ind < this.getNumberOfNodes(); ind++) {
            if (nodes.get(ind).isStateful()) {
                int ist = (int) sn.nodeToStation.get(ind);
                StatefulNode node_i = (StatefulNode) this.nodes.get(ind);
                Matrix state_i = FromMarginal.fromMarginalAndRunning(sn, ind, n.getRow(ist), s.getRow(ist));
                if (((StatefulNode) nodes.get(ind)).getState().isEmpty()) {
                    line_error(mfilename(new Object() {
                    }), "Invalid state assignment for station " + ind + "\n");
                }
                if (state_i.getNumRows() == 1) {
                    node_i.setState(state_i.getRow(0));
                    node_i.setStateSpace(state_i);
                    node_i.setStatePrior(Matrix.singleton(1));
                } else if (state_i.getNumRows() > 1) {
                    Matrix prior_state_i = new Matrix(state_i.getNumRows(), 1);
                    prior_state_i.zero();
                    prior_state_i.set(0, 0, 1.0);
                    node_i.setStateSpace(state_i);
                    node_i.setStatePrior(prior_state_i);
                    node_i.setState(state_i.getRow(0));
                } else {
                    node_i.setState(new Matrix(0, 0));
                    node_i.setStateSpace(new Matrix(0, 0));
                    node_i.setStatePrior(new Matrix(0, 0));
                }
            }
        }
        hasState = true;
    }

    public void initFromMarginalAndStarted(Matrix n, Matrix s) {
        NetworkStruct sn = getStruct();
        
        // Check if we need to validate with node-based indexing
        Matrix nValidation = n;
        Matrix sValidation = s;
        
        if (this.getNumberOfStations() < this.getNumberOfNodes()) {
            if (n.getNumRows() == this.getNumberOfStations()) {
                // Create node-based matrices for validation only
                Matrix nnodes = new Matrix(sn.nnodes, n.getNumCols());
                Matrix snodes = new Matrix(sn.nnodes, s.getNumCols());
                nnodes.zero();
                snodes.zero();
                for (int ist = 0; ist < sn.nstations; ist++) {
                    int ind = (int) sn.stationToNode.get(ist);
                    for (int c = 0; c < n.getNumCols(); c++) {
                        nnodes.set(ind, c, n.get(ist, c));
                        snodes.set(ind, c, s.get(ist, c));
                    }
                }
                nValidation = nnodes;
                sValidation = snodes;
            } else if (n.getNumRows() == this.getNumberOfNodes()) {
                // Already node-based
            } else {
                line_error(mfilename(new Object() {}), 
                    "The supplied matrix of marginal states does not have the correct number of rows. One either one per station or one per node.");
            }
        }
        
        boolean isvalidn = State.isValid(sn, nValidation, sValidation);
        if (!isvalidn) line_error(mfilename(new Object() {
        }), "Initial state is not valid.");

        for (int ind = 0; ind < this.getNumberOfNodes(); ind++) {
            if (nodes.get(ind).isStateful()) {
                int ist = (int) sn.nodeToStation.get(ind);
                StatefulNode node_i = (StatefulNode) this.nodes.get(ind);
                Matrix state_i = FromMarginal.fromMarginalAndStarted(sn, ind, n.getRow(ist), s.getRow(ist));
//                if (((StatefulNode) nodes.get(ind)).getState().isEmpty()) {
//                    line_error(mfilename(new Object() {
//                    }), "Invalid state assignment for station " + ind + "\n");
//                }
                if (state_i.getNumRows() == 1) {
                    node_i.setState(state_i.getRow(0));
                    node_i.setStatePrior(Matrix.singleton(1));
                    node_i.setStateSpace(state_i);
                } else if (state_i.getNumRows() > 1) {
                    Matrix prior_state_i = new Matrix(state_i.getNumRows(), 1);
                    prior_state_i.zero();
                    prior_state_i.set(0, 0, 1.0);
                    node_i.setStateSpace(state_i);
                    node_i.setStatePrior(prior_state_i);
                    node_i.setState(state_i.getRow(0));
                } else {
                    node_i.setState(new Matrix(0, 0));
                    node_i.setStateSpace(new Matrix(0, 0));
                    node_i.setStatePrior(new Matrix(0, 0));
                }
            }
        }

        hasState = true;
    }

    public RoutingMatrix initRoutingMatrix() {
        RoutingMatrix rt = new RoutingMatrix(this, jobClasses, nodes);
        return rt;
    }

    /**
     * Checks if this network is a Java native (JNetwork) implementation.
     * Always returns true for JNetwork implementations.
     *
     * @return true, as this is a Java implementation
     */
    public boolean isJavaNative() {
        return true;
    }

    public boolean isLimitedLoadDependent() {
        return this.getStruct().lldscaling.isEmpty();
    }

    /**
     * Checks if this network is a MATLAB native (MNetwork) implementation.
     * Always returns false for JNetwork implementations.
     *
     * @return false, as this is a Java implementation
     */
    public boolean isMatlabNative() {
        return false;
    }

    public boolean isStateValid() {
        NetworkStruct sn = this.getStruct();
        return snIsStateValid(sn);
    }

    public void jsimgView() {
        SolverJMT jmt = new SolverJMT(this);
        jmt.jsimgView();
    }

    public void jsimwView() {
        SolverJMT jmt = new SolverJMT(this);
        try {
            jmt.jsimwView();
        } catch (ParserConfigurationException e) {
            e.printStackTrace();
        }
    }

    // ========== TIKZ VISUALIZATION ==========


    /**
     * Generates TikZ code for visualizing this network.
     *
     * @return Complete LaTeX document with TikZ diagram
     */
    public String toTikZ() {
        return new jline.io.tikz.TikZExporter(this).generateTikZ();
    }

    /**
     * Generates TikZ code with custom options.
     *
     * @param options Configuration options for the visualization
     * @return Complete LaTeX document with TikZ diagram
     */
    public String toTikZ(jline.io.tikz.TikZOptions options) {
        return new jline.io.tikz.TikZExporter(this, options).generateTikZ();
    }

    /**
     * Displays this network as a TikZ diagram in a PDF viewer.
     * Requires pdflatex to be installed on the system.
     */
    public void tikzView() {
        new jline.io.tikz.TikZExporter(this).display();
    }

    /**
     * Displays this network as a TikZ diagram with custom options.
     *
     * @param options Configuration options for the visualization
     */
    public void tikzView(jline.io.tikz.TikZOptions options) {
        new jline.io.tikz.TikZExporter(this, options).display();
    }

    /**
     * Exports this network to a PNG file using TikZ.
     *
     * @param filePath The output file path (.png extension added if missing)
     * @throws java.io.IOException If file operations fail
     */
    public void tikzExportPNG(String filePath) throws java.io.IOException {
        String pngPath = filePath.endsWith(".png") ? filePath : filePath + ".png";
        new jline.io.tikz.TikZExporter(this).exportToPNG(pngPath);
    }

    /**
     * Exports this network to a PNG file using TikZ with custom DPI.
     *
     * @param filePath The output file path (.png extension added if missing)
     * @param dpi Resolution in dots per inch
     * @throws java.io.IOException If file operations fail
     */
    public void tikzExportPNG(String filePath, int dpi) throws java.io.IOException {
        String pngPath = filePath.endsWith(".png") ? filePath : filePath + ".png";
        new jline.io.tikz.TikZExporter(this).exportToPNG(pngPath, dpi);
    }

    /**
     * Exports this network to a PDF file using TikZ.
     *
     * @param filePath The output file path (without extension)
     * @return The generated PDF file
     * @throws java.io.IOException If file operations fail
     */
    public java.io.File exportTikZ(String filePath) throws java.io.IOException {
        jline.io.tikz.TikZExporter exporter = new jline.io.tikz.TikZExporter(this);
        String texPath = filePath.endsWith(".tex") ? filePath : filePath + ".tex";
        exporter.exportToFile(texPath);
        return exporter.exportToPDF();
    }

    /**
     * Exports this network's TikZ code to a .tex file.
     *
     * @param filePath The output file path
     * @throws java.io.IOException If file operations fail
     */
    public void exportTikZToFile(String filePath) throws java.io.IOException {
        jline.io.tikz.TikZExporter exporter = new jline.io.tikz.TikZExporter(this);
        exporter.exportToFile(filePath);
    }

    // ========== END TIKZ VISUALIZATION ==========

    public void link(RoutingMatrix P) {
        // MATLAB compatibility: check if network was already linked (re-link case)
        boolean isReset = false;
        if (this.sn != null) {
            isReset = true;
            resetNetwork(); // remove artificial class switch nodes, also clears connections
            // After resetNetwork, we need to ensure struct will be regenerated
            this.hasStruct = false;
        }

        // Only warn if addLink() was used BEFORE the first link() call (mixing approaches)
        // Note: Skip warning if this is a re-link case (isReset=true) since connections
        // were populated by the previous link() call, not by user calling addLink()
        if (!isReset && connections != null && !connections.isEmpty() && connections.elementSum() > 0) {
            line_warning(mfilename(new Object() {
            }), "The Network.link method cannot be used after calling the addLink() method. Use Node.setProbRouting instead to configure routing probabilities.");
        }
        
        sanitize();
        P.setRouting(this);

        for (Node node : this.nodes) {
            if (node instanceof Place) {
                ((Place) node).init();
            }
        }
        
        // Note: MATLAB also calls refreshChains at the end if isReset is true,
        // but this happens automatically in Java during struct regeneration
    }

    /**
     * Links the network with logging capability
     * Creates Logger nodes before and after specified stations and updates routing matrix
     *
     * @param P            the routing matrix
     * @param isNodeLogged boolean array indicating which nodes should be logged
     * @param logPath      path where log files will be stored (optional, uses existing logPath if null)
     * @return array containing [loggersBefore, loggersAfter] as List arrays
     */
    public List<Logger>[] linkAndLog(RoutingMatrix P, boolean[] isNodeLogged, String logPath) {
        // Reset struct to regenerate with loggers
        this.resetStruct();

        if (this.hasState) {
            // In MATLAB version this throws an error, but here we'll just warn
            line_warning(mfilename(new Object() {
            }), "The network state should be reset before calling linkAndLog.");
        }

        // Only reset if connections were actually populated (not just initialized empty by getConnectionMatrix)
        if (this.connections != null && !this.connections.isEmpty() && this.connections.elementSum() > 0) {
            line_warning(mfilename(new Object() {
            }), "Network topology already instantiated. Calling resetNetwork automatically before adding loggers.");
            this.resetNetwork();
        }

        int R = this.getNumberOfClasses();
        int Mnodes = this.getNumberOfNodes();

        if (Mnodes != isNodeLogged.length) {
            line_error(mfilename(new Object() {
            }), "The size of the isNodeLogged array does not match the number of nodes.");
        }

        // Prevent logging source and sink nodes - match MATLAB order: sink first, then source
        int sinkIndex = this.getIndexSinkNode();
        int sourceIndex = this.getIndexSourceNode();

        if (sinkIndex != -1 && isNodeLogged[sinkIndex]) {
            line_warning(mfilename(new Object() {
            }), "Sink station cannot be logged, ignoring.");
            isNodeLogged[sinkIndex] = false;
        }

        if (sourceIndex != -1 && isNodeLogged[sourceIndex]) {
            line_warning(mfilename(new Object() {
            }), "Source station cannot be logged, ignoring.");
            isNodeLogged[sourceIndex] = false;
        }

        // Set log path
        if (logPath != null) {
            this.setLogPath(logPath);
        } else {
            logPath = this.getLogPath();
        }

        if (logPath == null || logPath.isEmpty()) {
            line_error(mfilename(new Object() {
            }),"To instantiate a Logger, first use setLogPath method on the Network object to define the global path to save logs.");
        }

        // Create logger lists
        List<Logger> loggersBefore = new ArrayList<>();
        List<Logger> loggersAfter = new ArrayList<>();

        // Create departure loggers FIRST (to match expected node ordering)
        for (int ind = 0; ind < Mnodes; ind++) {
            if (isNodeLogged[ind]) {
                String nodeName = this.getNodeNames().get(ind);
                String logFileName = logPath + java.io.File.separator + nodeName + "-Dep.csv";
                String loggerName = "Dep_" + nodeName;
                Logger logger = new Logger(this, loggerName, logFileName);

                // Set default logging options
                logger.setTimestamp(true);
                logger.setJobID(true);
                logger.setJobClass(true);

                // Set routing for all job classes (matching MATLAB lines 63, 68)
                for (JobClass jobClass : this.getClasses()) {
                    logger.setRouting(jobClass, RoutingStrategy.RAND);
                }

                loggersAfter.add(logger);
            }
        }

        // Create arrival loggers SECOND
        for (int ind = 0; ind < Mnodes; ind++) {
            if (isNodeLogged[ind]) {
                String nodeName = this.getNodeNames().get(ind);
                String logFileName = logPath + java.io.File.separator + nodeName + "-Arv.csv";
                String loggerName = "Arv_" + nodeName;
                Logger logger = new Logger(this, loggerName, logFileName);

                // Set default logging options
                logger.setTimestamp(true);
                logger.setJobID(true);
                logger.setJobClass(true);

                // Set routing for all job classes (matching MATLAB lines 55, 60)
                for (JobClass jobClass : this.getClasses()) {
                    logger.setRouting(jobClass, RoutingStrategy.RAND);
                }

                loggersBefore.add(logger);
            }
        }

        // Create new routing matrix following MATLAB step-by-step approach
        // Step 1: Create expanded matrix of size 3*Mnodes (matching MATLAB line 81)
        int Mnodesnew = 3 * Mnodes;
        Matrix[][][] newP = new Matrix[R][R][1];
        for (int r = 0; r < R; r++) {
            for (int s = 0; s < R; s++) {
                newP[r][s][0] = new Matrix(Mnodesnew, Mnodesnew);
            }
        }

        // Step 2: Map original routing matrix considering logging (matching MATLAB lines 82-110)
        for (int r = 0; r < R; r++) {
            for (int s = 0; s < R; s++) {
                // class indexing starts at 1
                Matrix routingMatrix = P.get(1+r, 1+s);
                Matrix expandedMatrix = newP[r][s][0];

                for (int ind = 0; ind < Mnodes; ind++) {
                    for (int jnd = 0; jnd < Mnodes; jnd++) {
                        double prob = routingMatrix.get(ind, jnd);
                        if (prob > 0) {
                            if (isNodeLogged[ind] && isNodeLogged[jnd]) {
                                // link departure logger of source to arrival logger of destination
                                // Dep_i -> Arv_j (not Arv_i -> Dep_j)
                                expandedMatrix.set(Mnodes + ind, 2 * Mnodes + jnd, prob);
                            } else if (isNodeLogged[ind] && !isNodeLogged[jnd]) {
                                // link departure logger of source to destination
                                expandedMatrix.set(Mnodes + ind, jnd, prob);
                            } else if (!isNodeLogged[ind] && isNodeLogged[jnd]) {
                                // link source to arrival logger of destination
                                expandedMatrix.set(ind, 2 * Mnodes + jnd, prob);
                            } else {
                                // link i to j (matching MATLAB line 98)
                                expandedMatrix.set(ind, jnd, prob);
                            }
                        }
                    }
                }

                // Step 3: Add internal logger connections
                for (int ind = 0; ind < Mnodes; ind++) {
                    if (isNodeLogged[ind]) {
                        if (r == s) { // Only for same class transitions
                            // Flow: external -> Arv_Logger -> Node -> Dep_Logger -> external
                            expandedMatrix.set(2 * Mnodes + ind, ind, 1.0); // arrival logger -> original node
                            expandedMatrix.set(ind, Mnodes + ind, 1.0); // original node -> departure logger
                        }
                    }
                }
            }
        }

        // Step 4: Create index list of logged nodes (matching MATLAB line 113: idx = find(isNodeLogged))
        List<Integer> loggedIndices = new ArrayList<>();
        for (int i = 0; i < Mnodes; i++) {
            if (isNodeLogged[i]) {
                loggedIndices.add(i);
            }
        }

        // Step 5: Create final node list following MATLAB pattern [1:Mnodes,Mnodes+idx,2*Mnodes+idx]
        List<Node> finalNodes = new ArrayList<>();
        
        // Add all original nodes [1:Mnodes]
        for (int i = 0; i < Mnodes; i++) {
            finalNodes.add(this.nodes.get(i));
        }
        
        // Add departure loggers for logged nodes [Mnodes+idx]
        int loggerIdx = 0;
        for (int idx : loggedIndices) {
            finalNodes.add(loggersAfter.get(loggerIdx++));
        }
        
        // Add arrival loggers for logged nodes [2*Mnodes+idx]
        loggerIdx = 0;
        for (int idx : loggedIndices) {
            finalNodes.add(loggersBefore.get(loggerIdx++));
        }

        // Step 6: Extract final routing matrix (matching MATLAB line 114)
        RoutingMatrix finalP = new RoutingMatrix(this, this.jobClasses, finalNodes);
        int finalSize = finalNodes.size();
        
        for (int r = 0; r < R; r++) {
            for (int s = 0; s < R; s++) {
                Matrix finalMatrix = new Matrix(finalSize, finalSize);
                Matrix expandedMatrix = newP[r][s][0];
                
                // Create index mapping for final matrix extraction
                List<Integer> finalIndices = new ArrayList<>();
                // Add [1:Mnodes]
                for (int i = 0; i < Mnodes; i++) {
                    finalIndices.add(i);
                }
                // Add [Mnodes+idx]
                for (int idx : loggedIndices) {
                    finalIndices.add(Mnodes + idx);
                }
                // Add [2*Mnodes+idx]
                for (int idx : loggedIndices) {
                    finalIndices.add(2 * Mnodes + idx);
                }
                
                // Extract submatrix
                for (int i = 0; i < finalSize; i++) {
                    for (int j = 0; j < finalSize; j++) {
                        int origI = finalIndices.get(i);
                        int origJ = finalIndices.get(j);
                        finalMatrix.set(i, j, expandedMatrix.get(origI, origJ));
                    }
                }

                finalP.set(this.jobClasses.get(r), this.jobClasses.get(s), finalMatrix);
            }
        }

        // Link the network
        this.link(finalP);

        // Return loggers
        @SuppressWarnings("unchecked")
        List<Logger>[] result = new List[2];
        result[0] = loggersBefore;
        result[1] = loggersAfter;
        return result;
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

    // ========================================================================
    // SECTION 18: REFRESH METHODS
    // Methods for refreshing various network parameters and structures
    // ========================================================================

    public void refreshCapacity() {
        int M = this.stations.size();
        int K = this.jobClasses.size();
        int C = this.sn.nchains;

        Matrix classcap = new Matrix(M, K);
        classcap.fill(Inf);
        Matrix chaincap = new Matrix(M, K);
        chaincap.fill(Inf);
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
            for (int idx = 0; idx < inchain_c.length(); idx++) {
                chainCap += njobs.get(0, (int) inchain_c.get(0, idx));
            }

            if (chainCap >= Integer.MAX_VALUE) {
                chainCap = Inf;
            }

            for (int idx = 0; idx < inchain_c.length(); idx++) {
                int r = (int) inchain_c.get(idx);
                JobClass jobclass = this.jobClasses.get(r);
                for (int i = 0; i < M; i++) {
                    Station station = this.stations.get(i);
                    if (!(station instanceof Source)) {
                        DropStrategy stationDropRule = station.getDropRule(jobclass);
                        // Default to Drop for finite capacity queues if no explicit drop rule set
                        if (stationDropRule == null && !Double.isInfinite(station.getCap()) && station.getCap() < Integer.MAX_VALUE) {
                            stationDropRule = DropStrategy.Drop;
                        }
                        dropRule.get(station).put(jobclass, stationDropRule);
                    }
                    if (Double.isNaN(rates.get(i, r)) && !(station instanceof Place)) {
                        // Class doesn't visit this station - set to 0 per MATLAB convention
                        classcap.set(i, r, 0);
                        chaincap.set(i, c, 0);
                    } else {
                        classcap.set(i, r, chainCap);
                        chaincap.set(i, c, chainCap);
                        if (station.getClassCap(jobclass) >= 0)
                            classcap.set(i, r, FastMath.min(classcap.get(i, r), station.getClassCap(jobclass)));
                        if (station.getCap() >= 0)
                            classcap.set(i, r, FastMath.min(classcap.get(i, r), station.getCap()));
                    }
                }
            }
        }

        for (int i = 0; i < M; i++) {
            Station station = this.stations.get(i);
            // If station has explicit finite cap set, use it directly (Kendall K notation)
            // Otherwise use minimum of chain cap sum and class cap sum
            // Note: Integer.MAX_VALUE is the default "unset" value, not an explicit capacity
            if (station.getCap() >= 0 && station.getCap() < Integer.MAX_VALUE && !Double.isInfinite(station.getCap())) {
                // Explicit capacity set - use directly as total capacity
                capacity.set(i, 0, station.getCap());
            } else {
                capacity.set(i, 0, FastMath.min(chaincap.sumRows(i), classcap.sumRows(i)));
            }
        }

        this.sn.cap = capacity;
        this.sn.classcap = classcap;
        this.sn.droprule = dropRule;

        refreshRegions();
    }

    /**
     * Populate finite capacity region information in sn struct.
     * <p>
     * region is a CellMatrix of size F (number of regions).
     * region.get(f) is Matrix(M, K+1) where:
     *   entry (i,r) = max jobs of class r at station i in region f
     *   entry (i,K) = global max jobs at station i in region f
     *   -1 = infinite capacity
     *
     * @return the updated NetworkStruct
     */
    public NetworkStruct refreshRegions() {
        int M = this.sn.nstations;
        int K = this.sn.nclasses;
        int F = this.regions.size();
        this.sn.nregions = F;
        this.sn.region = new MatrixCell(F);
        // regionrule(f, r) = DropStrategy id for class r in region f
        this.sn.regionrule = new Matrix(F, K);
        this.sn.regionrule.fill(DropStrategy.Drop.getID());  // Default to drop
        // regionweight(f, r) = class weight for class r in region f
        this.sn.regionweight = new Matrix(F, K);
        this.sn.regionweight.fill(1.0);  // Default weight = 1.0
        // regionsz(f, r) = class size/memory for class r in region f
        this.sn.regionsz = new Matrix(F, K);
        this.sn.regionsz.fill(1.0);  // Default size = 1

        for (int f = 0; f < F; f++) {
            Region fcr = this.regions.get(f);
            // Matrix with M rows (stations) and K+1 columns (K classes + 1 global)
            Matrix regionMatrix = new Matrix(M, K + 1);
            regionMatrix.fill(-1);  // Initialize all to infinite (-1)

            // Find which stations are in this region and set their capacities
            for (Node node : fcr.getNodes()) {
                for (int i = 0; i < M; i++) {
                    if (this.stations.get(i) == node) {
                        // Set per-class max jobs for this station in this region
                        for (int r = 0; r < K; r++) {
                            int classMax = fcr.getClassMaxJobs(this.jobClasses.get(r));
                            regionMatrix.set(i, r, classMax);
                        }
                        // Set global max jobs for this station in this region (column K)
                        int globalMax = fcr.getGlobalMaxJobs();
                        regionMatrix.set(i, K, globalMax);
                        break;
                    }
                }
            }

            // Extract per-class drop rules, weights, and sizes for this region
            for (int r = 0; r < K; r++) {
                JobClass jobClass = this.jobClasses.get(r);
                DropStrategy classDropStrategy = fcr.getDropStrategy(jobClass);
                this.sn.regionrule.set(f, r, classDropStrategy.getID());
                this.sn.regionweight.set(f, r, fcr.getClassWeight(jobClass));
                this.sn.regionsz.set(f, r, fcr.getClassSize(jobClass));
            }

            this.sn.region.set(f, regionMatrix);
        }
        return this.sn;
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
                            if (rt.get(isf * K + r, jsf * K + s) > 0) csmask.set(r, s, 1.0);
                        }
                    }
                }
            }

            for (int isf = 0; isf < stateful.getNumRows(); isf++) {
                int ind = (int) this.sn.statefulToNode.get(0, isf);
                boolean isCS = (this.sn.nodetype.get(ind) == NodeType.Cache) || (this.sn.nodetype.get(ind) == NodeType.ClassSwitch);
                for (int r = 0; r < K; r++) {
                    csmask.set(r, r, 1.0);
                    for (int s = 0; s < K; s++) {
                        if (r != s) {
                            if (isCS) {
                                ClassSwitcher classSwitcher = (ClassSwitcher) this.nodes.get(ind).getServer();
                                if (classSwitcher.applyCsFun(r, s) > 0) csmask.set(r, s, 1.0);
                            }
                        }
                    }
                }
            }
            this.sn.csmask = csmask;
        } else {
            // Start with csMatrix from link(), but also include any transitions
            // that may have been created by stochastic complement through
            // non-stateful nodes (e.g., auto-added ClassSwitch nodes)
            Matrix csmask = this.csMatrix.copy();
            for (int r = 0; r < K; r++) {
                for (int s = 0; s < K; s++) {
                    for (int isf = 0; isf < stateful.getNumRows(); isf++) {
                        for (int jsf = 0; jsf < stateful.getNumRows(); jsf++) {
                            if (rt.get(isf * K + r, jsf * K + s) > 0) {
                                csmask.set(r, s, 1.0);
                            }
                        }
                    }
                }
            }
            this.sn.csmask = csmask;
        }

        if (((sn.refclass != null) && (!sn.refclass.isEmpty())) && (sn.refclass.length() < sn.nchains)) {
            sn.refclass.expandMatrix(1, sn.nchains, sn.nchains);
        }

        if (propagate) {
            //Compute visits
            this.sn = snRefreshVisits(this.sn, this.sn.chains, rt, rtnodes);

            //Call dependent capacity refresh
            refreshCapacity();
        }

        // Populate rtorig from rtnodes if not already set (when link() was not called)
        if (this.sn.rtorig == null || this.sn.rtorig.isEmpty()) {
            kotlin.Pair<List<List<Matrix>>, Matrix> rtorigResult = snRtnodesToRtorig(this.sn);
            List<List<Matrix>> rtorigcell = rtorigResult.getFirst();

            // Convert rtorigcell (List<List<Matrix>>) to Map<JobClass, Map<JobClass, Matrix>>
            Map<JobClass, Map<JobClass, Matrix>> rtorigMap = new HashMap<JobClass, Map<JobClass, Matrix>>();
            for (int r = 0; r < K; r++) {
                JobClass fromClass = this.jobClasses.get(r);
                Map<JobClass, Matrix> innerMap = new HashMap<JobClass, Matrix>();
                for (int s = 0; s < K; s++) {
                    JobClass toClass = this.jobClasses.get(s);
                    innerMap.put(toClass, rtorigcell.get(r).get(s));
                }
                rtorigMap.put(fromClass, innerMap);
            }
            this.sn.rtorig = rtorigMap;
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
        this.sn.nclosedjobs = Math.toIntExact(FastMath.round(njobsSum));
        this.sn.njobs = njobs.transpose();

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
        int sourceIdx = this.getIndexSourceStation();

        for (Integer i : statSet) {
            Station station = this.stations.get(i);
            Map<JobClass, SerializableFunction<Double, Double>> map = new HashMap<JobClass, SerializableFunction<Double, Double>>();
            for (Integer r : classSet) {
                JobClass jobclass = this.jobClasses.get(r);
                if (i == sourceIdx) {
                    Distribution distr = ((Source) station).getArrivalDistribution(jobclass);
                    if (distr instanceof Disabled) map.put(jobclass, null);
                    else map.put(jobclass, e -> ((ContinuousDistribution) distr).evalLST(e));
                } else {
                    //line 45-46 is ignored since Fork is not station
                    if (station instanceof Join) map.put(jobclass, null);
                    else
                        map.put(jobclass, e -> ((ContinuousDistribution) station.getServer().getServiceDistribution(jobclass)).evalLST(e));
                }
            }
            lst.put(station, map);
        }

        if (this.sn != null) this.sn.lst = lst;
    }

    public void refreshLocalVars() {
        int R = this.jobClasses.size();
        int I = this.nodes.size();
        Matrix nvars = new Matrix(I, 2 * R + 1);
        Map<Node, NodeParam> nodeparam = new HashMap<Node, NodeParam>();

        for (int ind = 0; ind < I; ind++) {
            Node node = this.nodes.get(ind);
            NodeParam param = null;
            switch (this.sn.nodetype.get(ind)) {
                case Cache:
                    CacheNodeParam cacheParam = new CacheNodeParam();
                    Cache cache = (Cache) node;
                    nvars.set(ind, 2 * R, cache.getItemLevelCap().elementSum());
                    cacheParam.nitems = 0;
                    cacheParam.accost = cache.accessProb;
                    for (int r = 0; r < this.getNumberOfClasses(); r++) {
                        if (!cache.popularityGet(r).isDisabled()) {
                            cacheParam.nitems = (int) Maths.max(cacheParam.nitems, cache.popularityGet(r).getSupport().getRight());
                        }
                    }
                    cacheParam.itemcap = cache.getItemLevelCap();
                    cacheParam.pread = new HashMap<>();
                    for (int r = 0; r < this.getNumberOfClasses(); r++) {
                        if (cache.popularityGet(r).isDisabled()) {
                            cacheParam.pread.put(r, null);
                        } else {
                            List<Double> t = new ArrayList<>();
                            for (int j = 1; j <= cacheParam.nitems; j++) {
                                t.add((double) j);
                            }
                            cacheParam.pread.put(r, ((DiscreteDistribution) cache.popularityGet(r)).evalPMF(t).toList1D());
                        }
                    }
                    cacheParam.replacestrat = cache.getReplacementStrategy();
                    cacheParam.hitclass = new Matrix(cache.getCacheServer().hitClass.getNumRows(), cache.getCacheServer().hitClass.getNumCols());
                    for (int r = 0; r < cache.getCacheServer().hitClass.getNumRows(); r++) {
                        for (int c = 0; c < cache.getCacheServer().hitClass.getNumCols(); c++) {
                            cacheParam.hitclass.set(r, c, FastMath.round(cache.getCacheServer().hitClass.get(r, c)));
                        }
                    }
                    cacheParam.missclass = new Matrix(cache.getCacheServer().missClass.getNumRows(), cache.getCacheServer().missClass.getNumCols());
                    for (int r = 0; r < cache.getCacheServer().missClass.getNumRows(); r++) {
                        for (int c = 0; c < cache.getCacheServer().missClass.getNumCols(); c++) {
                            cacheParam.missclass.set(r, c, FastMath.round(cache.getCacheServer().missClass.get(r, c)));
                        }
                    }
                    param = cacheParam;
                    break;
                case Fork:
                    ForkNodeParam forkParam = new ForkNodeParam();
                    forkParam.fanOut = ((Forker) node.getOutput()).tasksPerLink;
                    param = forkParam;
                    break;
                case Join:
                    JoinNodeParam joinParam = new JoinNodeParam();
                    Joiner joiner = (Joiner) node.getInput();
                    joinParam.joinStrategy = joiner.joinStrategy;
                    joinParam.fanIn = joiner.joinRequired;
                    param = joinParam;
                    break;
                case Logger:
                    LoggerNodeParam loggerParam = new LoggerNodeParam();
                    Logger logger = (Logger) node;
                    loggerParam.fileName.set(0, logger.getFileName());
                    loggerParam.filePath = logger.getFilePath();
                    loggerParam.startTime = logger.getStartTime();
                    loggerParam.loggerName = logger.getLoggerName();
                    loggerParam.timestamp = logger.getTimestamp();
                    loggerParam.jobID = logger.getJobID();
                    loggerParam.jobClass = logger.getJobClass();
                    loggerParam.timeSameClass = logger.getTimeSameClass();
                    loggerParam.timeAnyClass = logger.getTimeAnyClass();
                    param = loggerParam;
                    break;
                case Source:
                    ServiceNodeParam sourceParam = new ServiceNodeParam();
                    for (int r = 0; r < R; r++) {
                        Distribution arrivalDistrib = ((Source) node).getArrivalProcess(this.getClassByIndex(r));

                        if (arrivalDistrib != null) {
                            if (arrivalDistrib instanceof MAP) {
                                if (sourceParam.fileName == null)
                                    sourceParam.fileName = new ArrayList<>(Collections.nCopies(R, null));
                                nvars.set(ind, r, nvars.get(ind, r) + 1);
                            }

                            if (arrivalDistrib instanceof Replayer || arrivalDistrib instanceof Trace) {
                                if (sourceParam.fileName == null)
                                    sourceParam.fileName = new ArrayList<>(Collections.nCopies(R, null));
                                sourceParam.fileName.set(r, ((Replayer) arrivalDistrib).getFileName());
                            }
                        }
                    }
                    param = sourceParam;
                    break;
                case Queue:
                    QueueNodeParam queueParam = new QueueNodeParam();
                    Queue queue = (Queue) node;

                    // Handle polling parameters
                    if (queue.getSchedStrategy() == SchedStrategy.POLLING && queue.getServer() instanceof PollingServer) {
                        PollingServer pollingServer = (PollingServer) queue.getServer();
                        queueParam.pollingType = pollingServer.getPollingType();
                        if (queueParam.pollingType == PollingType.KLIMITED) {
                            queueParam.pollingPar = pollingServer.getPollingK();
                        }
                        
                        // Initialize switchover times for each job class
                        for (int r = 0; r < R; r++) {
                            JobClass jobClass = this.getClassByIndex(r);
                            Distribution switchover = pollingServer.getSwitchover(jobClass);
                            if (switchover != null) {
                                queueParam.switchoverTime.put(jobClass, switchover);
                            }
                        }
                    }

                    for (int r = 0; r < R; r++) {
                        ServiceBinding serviceProcess = node.getServer().getServiceProcess(this.getClassByIndex(r));
                        if (serviceProcess != null) {
                            Distribution serviceDistrib = serviceProcess.getDistribution();

                            if (serviceDistrib instanceof MAP) {
                                if (queueParam.fileName == null)
                                    queueParam.fileName = new ArrayList<>(Collections.nCopies(R, null));
                                nvars.set(ind, r, nvars.get(ind, r) + 1);
                            }

                            if (serviceDistrib instanceof Replayer || serviceDistrib instanceof Trace) {
                                if (queueParam.fileName == null)
                                    queueParam.fileName = new ArrayList<>(Collections.nCopies(R, null));
                                queueParam.fileName.set(r, ((Replayer) serviceDistrib).getFileName());
                            }
                        }
                    }

                    param = queueParam;
                    break;
                case Delay:
                    ServiceNodeParam delayParam = new ServiceNodeParam();

                    for (int r = 0; r < R; r++) {
                        ServiceBinding serviceProcess = node.getServer().getServiceProcess(this.getClassByIndex(r));
                        if (serviceProcess != null) {
                            Distribution serviceDistrib = serviceProcess.getDistribution();

                            if (serviceDistrib instanceof MAP) {
                                if (delayParam.fileName == null)
                                    delayParam.fileName = new ArrayList<>(Collections.nCopies(R, null));
                                nvars.set(ind, r, nvars.get(ind, r) + 1);
                            }

                            if (serviceDistrib instanceof Replayer || serviceDistrib instanceof Trace) {
                                if (delayParam.fileName == null)
                                    delayParam.fileName = new ArrayList<>(Collections.nCopies(R, null));
                                delayParam.fileName.set(r, ((Replayer) serviceDistrib).getFileName());
                            }
                        }
                    }

                    param = delayParam;
                    break;
                case Transition:
                    TransitionNodeParam transitionParam = new TransitionNodeParam();

                    for (int r = 0; r < R; r++) {
                        Distribution firingDistrib = node.getServer().getServiceDistribution(this.jobClasses.get(r));

                        if (firingDistrib instanceof MAP) {
                            if (transitionParam.fileName == null)
                                transitionParam.fileName = new ArrayList<>(Collections.nCopies(R, null));
                            nvars.set(ind, r, nvars.get(ind, r) + 1);
                        }

                        if (firingDistrib instanceof Replayer || firingDistrib instanceof Trace) {
                            if (transitionParam.fileName == null)
                                transitionParam.fileName = new ArrayList<>(Collections.nCopies(R, null));
                            transitionParam.fileName.set(r, ((Replayer) firingDistrib).getFileName());
                        }
                    }

                    param = transitionParam;
                    break;
                default:
                    break;
            }

            for (int r = 0; r < R; r++) {
                JobClass jobclass = this.jobClasses.get(r);
                if (param == null) {
                    param = new NodeParam();
                }
                // Declare variables outside switch to avoid scope issues
                Matrix conn_i;
                Matrix conn_i_transpose;
                switch (this.sn.routing.get(node).get(jobclass)) {
                    case KCHOICES:
                        throw new RuntimeException("Routing Strategy KCHOICES is not supported in JLINE");
                    case WRROBIN:
                        param.weights = new HashMap<JobClass, Matrix>();
                        param.outlinks = new HashMap<JobClass, Matrix>();
                        nvars.set(ind, R + r, nvars.get(ind, R + r) + 1);

                        //varsparam{ind}{r}.weights = zeros(1,self.sn.nnodes);
                        param.weights.put(jobclass, new Matrix(1, this.sn.nnodes));
                        //varsparam{ind}{r}.outlinks = find(self.sn.connmatrix(ind,:));
                        conn_i = new Matrix(0, 0);
                        Matrix.extractRows(this.sn.connmatrix, ind, ind + 1, conn_i);
                        conn_i_transpose = conn_i.find().transpose();
                        param.outlinks.put(jobclass, conn_i_transpose);

                        List<OutputStrategy> outputStrategy_r = node.getOutput().getOutputStrategyByClass(jobclass);
                        for (int c = 0; c < outputStrategy_r.size(); c++) {
                            Node destination = outputStrategy_r.get(c).getDestination();
                            Double weight = outputStrategy_r.get(c).getProbability();
                            param.weights.get(jobclass).set(0, destination.getNodeIndex(), weight);
                        }
                        break;
                    case RROBIN:
                        if (param.outlinks == null) { param.outlinks = new HashMap<JobClass, Matrix>(); }

                        nvars.set(ind, R + r, nvars.get(ind, R + r) + 1);

                        //varsparam{ind}{r}.outlinks = find(self.sn.connmatrix(ind,:));
                        conn_i = new Matrix(0, 0);
                        Matrix.extractRows(this.sn.connmatrix, ind, ind + 1, conn_i);
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
            // Initialize varsparam for cache state management
            this.sn.varsparam = new Matrix(I, 1);
            this.sn.varsparam.fill(-1); // -1 indicates no specific item selected
        }
    }

    public void refreshPetriNetNodes() {
        for (int ind = 0; ind < this.getNumberOfNodes(); ind++) {
            Node node = this.getNodeByIndex(ind);
            if (node instanceof Transition) {
                Transition transition = (Transition) node;
                TransitionNodeParam transitionParam = (TransitionNodeParam) this.sn.nodeparam.get(transition);
                if (transitionParam == null) {
                    // Create TransitionNodeParam if it doesn't exist (e.g., for transitions with no modes)
                    transitionParam = new TransitionNodeParam();
                    this.sn.nodeparam.put(transition, transitionParam);
                }
                transitionParam.nmodes = transition.getNumberOfModes();
                transitionParam.modenames = transition.getModeNames();
                transitionParam.enabling = new ArrayList<>();
                transitionParam.inhibiting = new ArrayList<>();
                transitionParam.firing = new ArrayList<>();
                transitionParam.timing = new ArrayList<>();
                for (Mode m : transition.getModes()) {
                    transitionParam.enabling.add(transition.enablingConditions.get(m));
                    transitionParam.inhibiting.add(transition.inhibitingConditions.get(m));
                    transitionParam.firing.add(transition.firingOutcomes.get(m));
                    transitionParam.timing.add(transition.timingStrategies.get(m));
                }
                transitionParam.nmodeservers = transition.getNumberOfModeServers();
                transitionParam.firingprio = transition.firingPriorities;
                transitionParam.fireweight = transition.firingWeights;

                transitionParam.firingproc = new HashMap<>();
                transitionParam.firingpie = new HashMap<>();
                transitionParam.firingphases = new Matrix(1, transition.getNumberOfModes());

                for (Mode m : transition.getModes()) {
                    if (transition.getFiringDistribution(m) instanceof Markovian) {
                        transitionParam.firingproc.put(m, ((Markovian) transition.getFiringDistribution(m)).getProcess());
                        transitionParam.firingpie.put(m, ((Markovian) transition.getFiringDistribution(m)).getInitProb());
                        transitionParam.firingphases.set(0, m.getIndex() - 1, (int) ((Markovian) transition.getFiringDistribution(m)).getNumberOfPhases());
                    } else {
                        // there is no simple way to pass the mean and scv since
                        // a transition is not a station so sn.scv and sn.rates do
                        // not have it, hence we save the parameters here
                        double mean = transition.getFiringDistribution(m).getMean();
                        double scv = transition.getFiringDistribution(m).getSCV();
                        MatrixCell representation = new MatrixCell();
                        representation.set(0, Matrix.singleton(mean));
                        representation.set(1, Matrix.singleton(scv));
                        transitionParam.firingproc.put(m, representation);
                        transitionParam.firingpie.put(m, Matrix.singleton(NaN));
                        transitionParam.firingphases.set(m.getIndex() - 1, NaN);
                    }
                    if (transitionParam.firingprocid == null) {
                        transitionParam.firingprocid = new HashMap<>();
                    }
                    transitionParam.firingprocid.put(m, ProcessType.fromDistribution(transition.getFiringDistribution(m)));
                }
            }
        }
    }

    public void refreshPriorities() {
        int K = this.jobClasses.size();
        Matrix classprio = new Matrix(1, K, K);
        for (int i = 0; i < K; i++) {
            classprio.set(0, i, jobClasses.get(i).priority);
        }

        if (this.sn != null) sn.classprio = classprio;
    }

    /**
     * Refreshes the deadline configuration for all job classes in the network structure.
     * Extracts deadline values from JobClass objects and populates the classdeadline matrix.
     */
    public void refreshDeadlines() {
        int K = this.jobClasses.size();
        Matrix classdeadline = new Matrix(1, K, K);
        for (int i = 0; i < K; i++) {
            classdeadline.set(0, i, jobClasses.get(i).deadline);
        }

        if (this.sn != null) sn.classdeadline = classdeadline;
    }

    @SuppressWarnings("unchecked")
    public void refreshProcessPhases(List<Integer> statSet, List<Integer> classSet) {
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
                        mu_i_val.set(0, 0, NaN);
                        phi_i_val.set(0, 0, NaN);
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

                if (!flag) phases.set(i, r, mu_val.length);
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
                    this.sn.phasessz.set(i, j, FastMath.max(1.0, phases.get(i, j)));
                }
            }
            Matrix.concatColumns(new Matrix(this.sn.phases.getNumRows(), 1), this.sn.phasessz.cumsumViaRow(), this.sn.phaseshift);
        }
        refreshProcessRepresentations();
    }

    @SuppressWarnings("unchecked")
    //Current this function not implement the return value
    public void refreshProcessRepresentations() {
        int M = this.stations.size();
        int K = this.jobClasses.size();
        Map<Station, Map<JobClass, MatrixCell>> ph = new HashMap<>();
        for (int i = 0; i < M; i++) {
            ph.put(this.stations.get(i), new HashMap<>());
        }
        Matrix phases = new Matrix(M, K, M * K);
        int sourceIdx = this.getIndexSourceStation();

        for (int i = 0; i < M; i++) {
            Station station = this.stations.get(i);
            Map<JobClass, MatrixCell> ph_i = new HashMap<JobClass, MatrixCell>();
            if (i == sourceIdx) {
                ph_i = (Map<JobClass, MatrixCell>) station.getSourceRates().get(0);
            } else {
                if (station instanceof Join) {
                    Coxian coxian = new Coxian(new ArrayList<Double>(Collections.singletonList(NaN)), new ArrayList<Double>(Collections.singletonList(NaN)));
                    for (JobClass jobclass : this.jobClasses)
                        ph_i.put(jobclass, coxian.getProcess());
                } else {
                    ph_i = (Map<JobClass, MatrixCell>) station.getServiceRates().get(0);
                }
            }
            for (JobClass jc : ph_i.keySet()) {
            }
            ph.put(station, ph_i);

            for (int r = 0; r < K; r++) {
                MatrixCell ph_i_r = ph_i.get(this.jobClasses.get(r));
                if (ph_i_r == null) {
                    phases.set(i, r, 1.0);
                } else if (ph_i_r.get(0) == null || ph_i_r.get(0).hasNaN()) {
                    // Disabled or non-Markovian with NaN parameters
                    phases.set(i, r, 0.0);
                } else if (ph_i_r.get(1) == null) {
                    // Non-Markovian distributions with single parameter (e.g., Det)
                    JobClass jobclass = this.jobClasses.get(r);
                    Distribution dist = null;
                    if (i == sourceIdx) {
                        dist = ((Source) station).getArrivalDistribution(jobclass);
                    } else if (station instanceof ServiceStation) {
                        dist = ((ServiceStation) station).getServiceProcess(jobclass);
                    }

                    // Convert non-Markovian distributions (including Det) to Erlang MAP approximation (mirrors MATLAB behavior)
                    if (dist != null && !(dist instanceof Disabled)) {
                        double targetMean = dist.getMean();
                        double scv = dist.getSCV();
                        int nPhases;
                        if (scv < GlobalConstants.CoarseTol) {
                            // Deterministic or near-deterministic: use 20 phases
                            nPhases = 20;
                        } else {
                            // Match SCV: for Erlang, SCV = 1/n, so n = 1/SCV
                            nPhases = FastMath.max(1, (int) FastMath.ceil(1.0 / scv));
                            nPhases = FastMath.min(nPhases, 100); // Cap at 100 phases
                        }
                        MatrixCell erlangMAP = map_erlang(targetMean, nPhases);
                        ph_i.put(jobclass, erlangMAP);
                        phases.set(i, r, nPhases);
                    } else {
                        phases.set(i, r, 1.0);
                    }
                } else if (ph_i_r.get(1).hasNaN()) {
                    phases.set(i, r, 0.0);
                } else {
                    phases.set(i, r, ph_i_r.get(0).getNumCols());
                }
                //Other situation set to 0 (The matrix initial value is 0)
            }
        }

        if (this.sn != null) {
            Map<Station, Map<JobClass, Matrix>> pie = new HashMap<>();
            for (int i = 0; i < M; i++) {
                Station station = this.stations.get(i);
                Map<JobClass, Matrix> pie_i = new HashMap<>();
                for (int r = 0; r < K; r++) {
                    JobClass jobclass = this.jobClasses.get(r);
                    MatrixCell map_ir = ph.get(station).get(jobclass);
                    if (map_ir != null && map_ir.get(0) != null && map_ir.get(1) != null) {
                        // PH/MAP representation has both D0 and D1
                        if (Double.isNaN(map_pie(map_ir.get(0), map_ir.get(1)).value())) {
                            map_pie(map_ir.get(0), map_ir.get(1));
                        }
                        pie_i.put(jobclass, map_pie(map_ir.get(0), map_ir.get(1)));
                    } else {
                        // Non-Markovian distributions or disabled - use NaN pie
                        Matrix tmp = new Matrix(1, 1, 0);
                        tmp.set(0, 0, NaN);
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
                    this.sn.phasessz.set(i, j, FastMath.max(1.0, phases.get(i, j)));
                }
            }
            //self.sn.phasessz(self.sn.nodeToStation(self.sn.nodetype == NodeType.Join),:)=phases(self.sn.nodeToStation(self.sn.nodetype == NodeType.Join),:);
            //Not tested, since current JLine does not support Join Node
            if (this.sn.nodeToStation != null) {
                for (int i = 0; i < this.sn.nodetype.size(); i++) {
                    if (this.sn.nodetype.get(i) != NodeType.Join) {
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

    public void refreshImpatience() {
        int M = this.stations.size();
        int K = this.jobClasses.size();

        // Initialize impatience data structures
        Map<Station, Map<JobClass, ProcessType>> impatienceType = new HashMap<>();
        Map<Station, Map<JobClass, Matrix>> impatienceMu = new HashMap<>();
        Map<Station, Map<JobClass, Matrix>> impatiencePhi = new HashMap<>();
        Map<Station, Map<JobClass, MatrixCell>> impatienceProc = new HashMap<>();
        Map<Station, Map<JobClass, Matrix>> impatiencePie = new HashMap<>();
        Map<Station, Map<JobClass, Integer>> impatiencePhases = new HashMap<>();

        for (int i = 0; i < M; i++) {
            Station station = this.stations.get(i);
            impatienceType.put(station, new HashMap<>());
            impatienceMu.put(station, new HashMap<>());
            impatiencePhi.put(station, new HashMap<>());
            impatienceProc.put(station, new HashMap<>());
            impatiencePie.put(station, new HashMap<>());
            impatiencePhases.put(station, new HashMap<>());

            for (int r = 0; r < K; r++) {
                JobClass jobclass = this.jobClasses.get(r);

                // Get patience distribution for this station-class pair
                Distribution patience = station.getPatience(jobclass);

                if (patience != null && !patience.isDisabled()) {
                    // Patience is configured - extract parameters
                    ProcessType procType = getProcessType(patience);
                    impatienceType.get(station).put(jobclass, procType);

                    // Extract rate/mean parameters
                    double mean = patience.getMean();
                    double rate = (mean > 0) ? 1.0 / mean : 0.0;
                    Matrix mu_val = new Matrix(1, 1);
                    mu_val.set(0, 0, rate);
                    impatienceMu.get(station).put(jobclass, mu_val);

                    // Extract SCV
                    double scv = patience.getSCV();
                    Matrix phi_val = new Matrix(1, 1);
                    phi_val.set(0, 0, scv);
                    impatiencePhi.get(station).put(jobclass, phi_val);

                    // Extract process representation for PH-type distributions
                    if (patience instanceof jline.lang.processes.Markovian) {
                        MatrixCell proc_val = ((jline.lang.processes.Markovian) patience).getProcess();
                        impatienceProc.get(station).put(jobclass, proc_val);

                        // Extract initial probability vector
                        if (proc_val != null && proc_val.get(0) != null && proc_val.get(1) != null) {
                            Matrix pie_val = map_pie(proc_val.get(0), proc_val.get(1));
                            impatiencePie.get(station).put(jobclass, pie_val);

                            // Number of phases
                            int nPhases = proc_val.get(0).getNumCols();
                            impatiencePhases.get(station).put(jobclass, nPhases);
                        } else {
                            // Simple distribution (Exp, Det, etc.)
                            Matrix pie_val = new Matrix(1, 1);
                            pie_val.set(0, 0, 1.0);
                            impatiencePie.get(station).put(jobclass, pie_val);
                            impatiencePhases.get(station).put(jobclass, 1);
                        }
                    } else {
                        // Non-Markovian distribution
                        Matrix pie_val = new Matrix(1, 1);
                        pie_val.set(0, 0, 1.0);
                        impatiencePie.get(station).put(jobclass, pie_val);
                        impatiencePhases.get(station).put(jobclass, 1);
                    }
                }
                // If no patience configured, don't add entries (sparse storage)
            }
        }

        // Store in network structure
        if (this.sn != null) {
            this.sn.impatienceType = impatienceType;
            this.sn.impatienceMu = impatienceMu;
            this.sn.impatiencePhi = impatiencePhi;
            this.sn.impatienceProc = impatienceProc;
            this.sn.impatiencePie = impatiencePie;
            this.sn.impatiencePhases = impatiencePhases;
        }
    }

    /**
     * Refreshes balking configuration in the network structure.
     * Extracts balking strategy and thresholds from all station-class pairs.
     */
    public void refreshBalking() {
        int M = this.stations.size();
        int K = this.jobClasses.size();

        // Initialize balking data structures
        Map<Station, Map<JobClass, BalkingStrategy>> balkingStrategy = new HashMap<>();
        Map<Station, Map<JobClass, List<BalkingThreshold>>> balkingThresholds = new HashMap<>();

        for (int i = 0; i < M; i++) {
            Station station = this.stations.get(i);
            balkingStrategy.put(station, new HashMap<>());
            balkingThresholds.put(station, new HashMap<>());

            for (int r = 0; r < K; r++) {
                JobClass jobclass = this.jobClasses.get(r);

                // Get balking configuration for this station-class pair
                BalkingStrategy strategy = station.getBalkingStrategy(jobclass);
                List<BalkingThreshold> thresholds = station.getBalkingThresholds(jobclass);

                if (strategy != null && thresholds != null && !thresholds.isEmpty()) {
                    balkingStrategy.get(station).put(jobclass, strategy);
                    balkingThresholds.get(station).put(jobclass, thresholds);
                }
            }
        }

        // Store in network structure
        if (this.sn != null) {
            this.sn.balkingStrategy = balkingStrategy;
            this.sn.balkingThresholds = balkingThresholds;
        }
    }

    /**
     * Refreshes retrial configuration in the network structure.
     * Extracts retrial delay distributions and max attempts from all station-class pairs.
     */
    public void refreshRetrial() {
        int M = this.stations.size();
        int K = this.jobClasses.size();

        // Initialize retrial data structures
        Map<Station, Map<JobClass, ProcessType>> retrialType = new HashMap<>();
        Map<Station, Map<JobClass, Matrix>> retrialMu = new HashMap<>();
        Map<Station, Map<JobClass, Matrix>> retrialPhi = new HashMap<>();
        Map<Station, Map<JobClass, MatrixCell>> retrialProc = new HashMap<>();
        Map<Station, Map<JobClass, Integer>> retrialMaxAttempts = new HashMap<>();

        for (int i = 0; i < M; i++) {
            Station station = this.stations.get(i);
            retrialType.put(station, new HashMap<>());
            retrialMu.put(station, new HashMap<>());
            retrialPhi.put(station, new HashMap<>());
            retrialProc.put(station, new HashMap<>());
            retrialMaxAttempts.put(station, new HashMap<>());

            for (int r = 0; r < K; r++) {
                JobClass jobclass = this.jobClasses.get(r);

                // Get retrial configuration for this station-class pair
                Distribution retrialDist = station.getRetrialDelayDistribution(jobclass);

                if (retrialDist != null && !retrialDist.isDisabled()) {
                    // Retrial is configured - extract parameters
                    ProcessType procType = getProcessType(retrialDist);
                    retrialType.get(station).put(jobclass, procType);

                    // Extract rate/mean parameters
                    double mean = retrialDist.getMean();
                    double rate = (mean > 0) ? 1.0 / mean : 0.0;
                    Matrix mu_val = new Matrix(1, 1);
                    mu_val.set(0, 0, rate);
                    retrialMu.get(station).put(jobclass, mu_val);

                    // Extract SCV
                    double scv = retrialDist.getSCV();
                    Matrix phi_val = new Matrix(1, 1);
                    phi_val.set(0, 0, scv);
                    retrialPhi.get(station).put(jobclass, phi_val);

                    // Extract process representation for PH-type distributions
                    if (retrialDist instanceof jline.lang.processes.Markovian) {
                        MatrixCell proc_val = ((jline.lang.processes.Markovian) retrialDist).getProcess();
                        retrialProc.get(station).put(jobclass, proc_val);
                    }

                    // Store max attempts
                    int maxAttempts = station.getMaxRetrialAttempts(jobclass);
                    retrialMaxAttempts.get(station).put(jobclass, maxAttempts);
                }
            }
        }

        // Store in network structure
        if (this.sn != null) {
            this.sn.retrialType = retrialType;
            this.sn.retrialMu = retrialMu;
            this.sn.retrialPhi = retrialPhi;
            this.sn.retrialProc = retrialProc;
            this.sn.retrialMaxAttempts = retrialMaxAttempts;
        }
    }

    /**
     * Populates heterogeneous server configuration in NetworkStruct.
     * <p>
     * This method extracts heterogeneous server type information from Queue nodes
     * and populates the corresponding fields in NetworkStruct:
     * <ul>
     * <li>nservertypes - number of server types per station</li>
     * <li>servertypenames - names of server types per station</li>
     * <li>serverspertype - number of servers per type per station</li>
     * <li>servercompat - compatibility matrix (server type x job class)</li>
     * <li>heterorates - service rates per server type per class</li>
     * <li>heteroproc - PH process representation per server type per class</li>
     * <li>heteroprocid - process type per server type per class</li>
     * <li>heteroschedpolicy - heterogeneous scheduling policy per station</li>
     * </ul>
     */
    public void refreshHeterogeneousServers() {
        if (this.sn == null) return;

        int M = this.stations.size();
        int K = this.jobClasses.size();

        // Initialize all heterogeneous fields
        sn.nservertypes = new Matrix(M, 1, M);
        sn.servertypenames = new HashMap<Station, List<String>>();
        sn.serverspertype = new HashMap<Station, Matrix>();
        sn.servercompat = new HashMap<Station, Matrix>();
        sn.heterorates = new HashMap<Station, Map<Integer, Map<Integer, Double>>>();
        sn.heteroproc = new HashMap<Station, Map<Integer, Map<Integer, MatrixCell>>>();
        sn.heteroprocid = new HashMap<Station, Map<Integer, Map<Integer, ProcessType>>>();
        sn.heteroschedpolicy = new HashMap<Station, HeteroSchedPolicy>();

        for (int i = 0; i < M; i++) {
            Station station = this.stations.get(i);

            // Only Queue stations can have heterogeneous servers
            if (!(station instanceof Queue)) {
                sn.nservertypes.set(i, 0, 0.0);
                continue;
            }

            Queue queue = (Queue) station;
            if (!queue.isHeterogeneous()) {
                sn.nservertypes.set(i, 0, 0.0);
                continue;
            }

            List<ServerType> serverTypes = queue.getServerTypes();
            int nTypes = serverTypes.size();
            sn.nservertypes.set(i, 0, nTypes);

            // Populate server type names
            List<String> names = new ArrayList<String>();
            for (ServerType st : serverTypes) {
                names.add(st.getName());
            }
            sn.servertypenames.put(station, names);

            // Populate servers per type: Matrix (nTypes x 1)
            Matrix spt = new Matrix(nTypes, 1, nTypes);
            for (int t = 0; t < nTypes; t++) {
                spt.set(t, 0, serverTypes.get(t).getNumOfServers());
            }
            sn.serverspertype.put(station, spt);

            // Populate compatibility matrix: Matrix (nTypes x K)
            Matrix compat = new Matrix(nTypes, K, nTypes * K);
            for (int t = 0; t < nTypes; t++) {
                ServerType st = serverTypes.get(t);
                for (int r = 0; r < K; r++) {
                    JobClass jc = this.jobClasses.get(r);
                    compat.set(t, r, st.isCompatible(jc) ? 1.0 : 0.0);
                }
            }
            sn.servercompat.put(station, compat);

            // Populate heterogeneous service rates and processes
            Map<Integer, Map<Integer, Double>> stationRates = new HashMap<Integer, Map<Integer, Double>>();
            Map<Integer, Map<Integer, MatrixCell>> stationProc = new HashMap<Integer, Map<Integer, MatrixCell>>();
            Map<Integer, Map<Integer, ProcessType>> stationProcId = new HashMap<Integer, Map<Integer, ProcessType>>();

            Map<ServerType, Map<JobClass, Distribution>> heteroDistrs = queue.getHeteroServiceDistributions();
            for (int t = 0; t < nTypes; t++) {
                ServerType st = serverTypes.get(t);
                Map<JobClass, Distribution> classMap = heteroDistrs.get(st);

                Map<Integer, Double> typeRates = new HashMap<Integer, Double>();
                Map<Integer, MatrixCell> typeProc = new HashMap<Integer, MatrixCell>();
                Map<Integer, ProcessType> typeProcId = new HashMap<Integer, ProcessType>();

                for (int r = 0; r < K; r++) {
                    JobClass jc = this.jobClasses.get(r);
                    if (classMap != null && classMap.containsKey(jc)) {
                        Distribution distr = classMap.get(jc);
                        typeRates.put(r, distr.getRate());
                        // Get PH representation for Markovian distributions
                        if (distr instanceof Markovian) {
                            typeProc.put(r, ((Markovian) distr).getProcess());
                        }
                        typeProcId.put(r, getProcessType(distr));
                    }
                }
                stationRates.put(t, typeRates);
                stationProc.put(t, typeProc);
                stationProcId.put(t, typeProcId);
            }
            sn.heterorates.put(station, stationRates);
            sn.heteroproc.put(station, stationProc);
            sn.heteroprocid.put(station, stationProcId);

            // Scheduling policy
            sn.heteroschedpolicy.put(station, queue.getHeteroSchedPolicy());
        }
    }

    public void refreshProcessTypes(List<Integer> statSet, List<Integer> classSet) {
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
                throw new Exception("refreshProcessTypes requires either both null or not null parameters");
            } catch (Exception e) {
                e.printStackTrace();
            }
        } else {
            proctype = this.sn.procid;
        }
        boolean hasOpenClass = this.hasOpenClasses();
        int sourceIdx = getIndexSourceStation();

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
                    } else {
                        // Fallback for any other station with ServiceTunnel (e.g., Fork if treated as station)
                        map.put(jobclass, ProcessType.DISABLED);
                    }
                } else {
                    if (!hasOpenClass || i != sourceIdx) {
                        if (station instanceof Place || !station.getServer().containsJobClass(jobclass)) {
                            map.put(jobclass, ProcessType.DISABLED);
                        } else {
                            Distribution distr = station.getServer().getServiceDistribution(jobclass);
                            map.put(jobclass, getProcessType(distr));
                        }
                    } else {
                        // For open class source stations with non-ServiceTunnel servers, leave as disabled
                        // (matching MATLAB behavior where procid stays NaN)
                        map.put(jobclass, ProcessType.DISABLED);
                    }
                }
            }
            proctype.put(station, map);
        }

        if (this.sn != null) this.sn.procid = proctype;
    }

    public void refreshProcesses(List<Integer> statSet, List<Integer> classSet) {
        boolean[] status = refreshRates(statSet, classSet);
        boolean hasSCVChanged = status[1];
        boolean hasRateChanged = status[0];

        if (hasSCVChanged) {
            refreshProcessTypes(statSet, classSet);
            refreshProcessPhases(statSet, classSet);
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

    public void refreshProcesses() {
        refreshProcesses(null, null);
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
            scv.fill(NaN);
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

            rates = this.sn.rates.copy();
            scv = this.sn.scv.copy();
            rates_orig = this.sn.rates.copy();
            scv_orig = this.sn.scv.copy();
        }
        boolean hasOpenClasses = this.hasOpenClasses();
        int sourceIdx = getIndexSourceStation();

        for (Integer i : statSet) {
            Station station = stations.get(i);
            for (Integer r : classSet) {
                if (station instanceof Place) {
                    rates.set(i, r, NaN);
                    scv.set(i, r, NaN);
                } else if (station.getServer() instanceof ServiceTunnel) {
                    if (station instanceof Source) {
                        if (!((Source) station).containsJobClass(this.jobClasses.get(r))) {
                            rates.set(i, r, NaN);
                            scv.set(i, r, NaN);
                        } else {
                            Distribution distr = ((Source) station).getArrivalDistribution(this.jobClasses.get(r));
                            rates.set(i, r, distr.getRate());
                            scv.set(i, r, distr.getSCV());
                        }
                    } else if (station instanceof Join) {
                        rates.set(i, r, Inf);
                        scv.set(i, r, 0.0);
                    }
                } else {
                    if (!hasOpenClasses || i != sourceIdx) {
                        if (!station.getServer().containsJobClass(this.jobClasses.get(r))) {
                            rates.set(i, r, NaN);
                            scv.set(i, r, NaN);
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
            tmp.absEq();
            if (tmp.elementSum() > 0) hasRateChanged = true;
        }

        if (!hasSCVChanged) {
            Matrix tmp = scv.sub(1, scv_orig);
            tmp.absEq();
            if (tmp.elementSum() > 0) hasSCVChanged = true;
        }

        if (hasRateChanged) {
            this.sn.rates = rates;
        }

        if (hasSCVChanged) {
            this.sn.scv = scv;
        }

        return new boolean[]{hasRateChanged, hasSCVChanged};
    }

    public void refreshRoutingMatrix(Matrix rates) {
        if (rates == null)
            line_error(mfilename(new Object() {
            }), "refreshRoutingMatrix cannot retrieve station rates, pass them as an input parameters.");

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

        if (this.enableChecks) {
            outerloop:
            for (JobClass jobclass : this.jobClasses) {
                for (Map<JobClass, RoutingStrategy> nodeRoutingMap : this.sn.routing.values()) {
                    if (nodeRoutingMap.get(jobclass) != RoutingStrategy.DISABLED) continue outerloop;
                }
                throw new RuntimeException("Routing strategy is unspecified at all nodes");
            }
        }

        boolean isStateDep = (Matrix.extractColumn(this.sn.isstatedep, 2, null).getNonZeroLength() > 0);
        Map<Integer, Map<Integer, SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Double>>> rtnodefuncell = new HashMap<Integer, Map<Integer, SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Double>>>();

        if (isStateDep) {
            for (int ind = 0; ind < M; ind++) {
                final int ind_final = ind;
                for (int jnd = 0; jnd < M; jnd++) {
                    final int jnd_final = jnd;
                    for (int r = 0; r < K; r++) {
                        final int r_final = r;
                        for (int s = 0; s < K; s++) {
                            final int s_final = s;
                            Map<Integer, SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Double>> map = rtnodefuncell.getOrDefault(ind * K + r, new HashMap<Integer, SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Double>>());
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
                                if (val != 0) cellfunnodes.set(row, col, val);
                            }
                        }
                    }
                }
                return dtmc_stochcomp(cellfunnodes, statefulNodeClasses);
            };
        } else {
            rtfun = ((pair) -> dtmc_stochcomp(rtnodes, statefulNodeClasses));
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
            double val = this.sn.refstat.get((int) inchain_c.value(), 0);
            for (int col = 1; col < inchain_c.getNumCols(); col++) {
                if (val != this.sn.refstat.get((int) inchain_c.get(0, col), 0))
                    throw new RuntimeException("Classes within chain have different reference station");
            }
        }
    }

    public void refreshScheduling() {
        int M = this.stations.size();
        int K = this.jobClasses.size();
        Map<Station, SchedStrategy> sched = getStationScheduling();
        Matrix schedparam = new Matrix(M, K, M * K);
        int sourceIdx = this.getIndexSourceStation();

        for (int i = 0; i < M; i++) {
            Station station = this.stations.get(i);
            if (sourceIdx == -1 || i != sourceIdx) {
                if (!(station.getServer() instanceof ServiceTunnel)) {
                    if (station instanceof Queue) {
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
                            // Handle SEPT and LEPT scheduling strategies
                            SchedStrategy schedStrategy = queue.getSchedStrategy();
                            if (schedStrategy == SchedStrategy.SEPT || schedStrategy == SchedStrategy.LEPT) {
                                // Initialize schedparam for SEPT/LEPT
                                // Collect service times for all classes at this station
                                List<Double> serviceTimes = new ArrayList<>();
                                List<Integer> classIndices = new ArrayList<>();
                                
                                for (int r = 0; r < K; r++) {
                                    Distribution serviceDistr = queue.getService(this.jobClasses.get(r));
                                    if (serviceDistr != null) {
                                        double serviceTime = serviceDistr.getMean();
                                        if (!Double.isNaN(serviceTime) && !Double.isInfinite(serviceTime)) {
                                            serviceTimes.add(serviceTime);
                                            classIndices.add(r);
                                        }
                                    }
                                }
                                
                                if (!serviceTimes.isEmpty()) {
                                    // Create a list of unique service times
                                    List<Double> uniqueTimes = new ArrayList<>(new HashSet<>(serviceTimes));
                                    
                                    // Sort based on strategy
                                    if (schedStrategy == SchedStrategy.SEPT) {
                                        // Shortest first
                                        Collections.sort(uniqueTimes);
                                    } else { // LEPT
                                        // Longest first
                                        Collections.sort(uniqueTimes, Collections.reverseOrder());
                                    }
                                    
                                    // Assign priorities based on sorted order
                                    for (int idx = 0; idx < classIndices.size(); idx++) {
                                        int classIdx = classIndices.get(idx);
                                        double serviceTime = serviceTimes.get(idx);
                                        int priority = uniqueTimes.indexOf(serviceTime) + 1;
                                        schedparam.set(i, classIdx, priority);
                                        queue.setSchedStrategyPar(this.jobClasses.get(classIdx), priority);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if (this.sn != null) {
            this.sn.sched = sched;
            this.sn.schedparam = schedparam;
            for (int i = 0; i < M; i++) {
                this.sn.sched.put(this.getStations().get(i), sched.getOrDefault(this.stations.get(i), null));
            }
            //No schedid in JLine
        }
    }

    public void refreshStruct() {
        refreshStruct(true);

    }

    public void refreshStruct(boolean hardRefresh) {

        sanitize();
        resolveSignals();  // Resolve Signal placeholders to OpenSignal or ClosedSignal

        List<NodeType> nodetypes;
        List<String> classnames;
        List<String> nodenames;
        Matrix refstat;
        Matrix conn;
        Matrix njobs;
        Matrix numservers;
        Matrix lldscaling;
        Map<Station, SerializableFunction<Matrix, Double>> cdscaling;
        Map<Station, Matrix> ljdscaling;
        Map<Station, Matrix> ljdcutoffs;
        Map<Station, Map<JobClass, Matrix>> ljcdscaling;
        Map<Station, Matrix> ljcdcutoffs;
        Map<Node, Map<JobClass, RoutingStrategy>> routing;


        if (this.hasStruct && !hardRefresh) {
            nodetypes = sn.nodetype;
            classnames = sn.classnames;
            nodenames = sn.nodenames;
            refstat = sn.refstat;
        } else {
            nodetypes = getNodeTypes();
            classnames = getClassNames();
            nodenames = getNodeNames();
            refstat = getReferenceStations();

            // Append FCR names and types to node lists (only when refreshing)
            for (Region fcr : this.regions) {
                nodenames.add(fcr.getName());
                nodetypes.add(NodeType.Region);
            }
        }

        conn = getConnectionMatrix();
        njobs = getNumberOfJobs();
        numservers = getStationServers();
        lldscaling = getLimitedLoadDependence();
        cdscaling = getLimitedClassDependence();
        Pair<Map<Station, Matrix>, Map<Station, Matrix>> ljdResult = getLimitedJointDependence();
        ljdscaling = ljdResult.getLeft();
        ljdcutoffs = ljdResult.getRight();
        Pair<Map<Station, Map<JobClass, Matrix>>, Map<Station, Matrix>> ljcdResult = getLimitedJointClassDependence();
        ljcdscaling = ljcdResult.getLeft();
        ljcdcutoffs = ljcdResult.getRight();

        if (sn == null) sn = new NetworkStruct();

        // sn.nnodes counts physical nodes only; FCRs are virtual nodes appended to nodenames/nodetypes
        sn.nnodes = this.nodes.size();
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

        sn.isslc = new Matrix(sn.nclasses, 1, sn.nclasses);
        for (int c = 0; c < sn.nclasses; c++) {
            if (this.jobClasses.get(c) instanceof SelfLoopingClass) {
                sn.isslc.set(c, 0, 1.0);
            }
        }

        // Initialize issignal - default to 0 (false) for all classes
        sn.issignal = new Matrix(sn.nclasses, 1, sn.nclasses);
        // Initialize signaltype - null for non-signal classes
        sn.signaltype = new ArrayList<SignalType>(sn.nclasses);
        for (int c = 0; c < sn.nclasses; c++) {
            sn.signaltype.add(null);
        }
        // Detect Signal, OpenSignal, and ClosedSignal classes and populate issignal/signaltype
        for (int c = 0; c < sn.nclasses; c++) {
            if (this.jobClasses.get(c) instanceof Signal) {
                Signal signalClass = (Signal) this.jobClasses.get(c);
                sn.issignal.set(c, 0, 1.0);
                sn.signaltype.set(c, signalClass.getSignalType());
            } else if (this.jobClasses.get(c) instanceof OpenSignal) {
                OpenSignal signalClass = (OpenSignal) this.jobClasses.get(c);
                sn.issignal.set(c, 0, 1.0);
                sn.signaltype.set(c, signalClass.getSignalType());
            } else if (this.jobClasses.get(c) instanceof ClosedSignal) {
                ClosedSignal signalClass = (ClosedSignal) this.jobClasses.get(c);
                sn.issignal.set(c, 0, 1.0);
                sn.signaltype.set(c, signalClass.getSignalType());
            }
        }

        // Initialize syncreply - maps each class to its expected reply signal class index (-1 if none)
        // Note: JobClass.replySignalClassIndex is 1-based (MATLAB convention), convert to 0-based for sn.syncreply
        sn.syncreply = new Matrix(sn.nclasses, 1, sn.nclasses);
        for (int c = 0; c < sn.nclasses; c++) {
            int replyIdx = this.jobClasses.get(c).getReplySignalClassIndex();
            // Convert from 1-based to 0-based if valid (>= 1 means valid signal class index)
            sn.syncreply.set(c, 0, replyIdx >= 1 ? replyIdx - 1 : -1);
        }

        // Initialize signal removal configuration fields
        sn.signalRemovalDist = new ArrayList<DiscreteDistribution>(sn.nclasses);
        sn.signalRemovalPolicy = new ArrayList<RemovalPolicy>(sn.nclasses);
        sn.isCatastrophe = new Matrix(sn.nclasses, 1, sn.nclasses);
        for (int c = 0; c < sn.nclasses; c++) {
            sn.signalRemovalDist.add(null);
            sn.signalRemovalPolicy.add(null);
        }
        // Populate removal configuration from Signal, OpenSignal, and ClosedSignal classes
        for (int c = 0; c < sn.nclasses; c++) {
            JobClass jc = this.jobClasses.get(c);
            if (jc instanceof Signal) {
                Signal signalClass = (Signal) jc;
                if (signalClass.isCatastrophe()) {
                    sn.isCatastrophe.set(c, 0, 1.0);
                }
                sn.signalRemovalDist.set(c, signalClass.getRemovalDistribution());
                sn.signalRemovalPolicy.set(c, signalClass.getRemovalPolicy());
            } else if (jc instanceof OpenSignal) {
                OpenSignal signalClass = (OpenSignal) jc;
                if (signalClass.isCatastrophe()) {
                    sn.isCatastrophe.set(c, 0, 1.0);
                }
                sn.signalRemovalDist.set(c, signalClass.getRemovalDistribution());
                sn.signalRemovalPolicy.set(c, signalClass.getRemovalPolicy());
            } else if (jc instanceof ClosedSignal) {
                ClosedSignal signalClass = (ClosedSignal) jc;
                if (signalClass.isCatastrophe()) {
                    sn.isCatastrophe.set(c, 0, 1.0);
                }
                sn.signalRemovalDist.set(c, signalClass.getRemovalDistribution());
                sn.signalRemovalPolicy.set(c, signalClass.getRemovalPolicy());
            }
        }

        sn.nclosedjobs = DoubleStream.of(njobs.getNonZeroValues())
                .filter(Double::isFinite)
                .mapToInt(d -> (int) Math.round(d))
                .sum();
        sn.nservers = numservers;
        sn.isstation = new Matrix(nodes.size(), 1, nodes.size());
        for (int i = 0; i < sn.nnodes; i++) {
            Node node = this.nodes.get(i);
            sn.isstation.set(i, 0, 0.0); //state dependent service
            if (node instanceof Source || node instanceof Delay || node instanceof Queue || node instanceof Join || node instanceof Place) {
                sn.isstation.set(i, 0, 1.0); //state dependent service
            }
        }
        sn.nstations = stations.size();
        sn.nodetype = nodetypes;
        sn.scv = new Matrix(sn.nstations, sn.nclasses, sn.nstations * sn.nclasses);
        sn.scv.fill(NaN);
        sn.njobs = njobs.transpose();
        sn.refstat = refstat;
        sn.space = new HashMap<StatefulNode, Matrix>();
        for (int i = 0; i < sn.nstations; i++)
            sn.space.put(stations.get(i), new Matrix(0, 0, 0));
        sn.routing = routing;
        sn.chains = new Matrix(0, 0);
        sn.lst = null;
        sn.lldscaling = lldscaling;
        sn.cdscaling = cdscaling;
        sn.ljdscaling = ljdscaling;
        sn.ljdcutoffs = ljdcutoffs;
        sn.ljcdscaling = ljcdscaling;
        sn.ljcdcutoffs = ljcdcutoffs;
        sn.nodetype = nodetypes;
        sn.isstateful = new Matrix(nodes.size(), 1, nodes.size());
        for (int i = 0; i < sn.nnodes; i++) {
            Node node = this.nodes.get(i);
            sn.isstateful.set(i, 0, 0.0); //state dependent service
            if (node instanceof Source || node instanceof Delay || node instanceof Queue || node instanceof Cache || node instanceof Join || node instanceof Router || node instanceof Place || node instanceof Transition) {
                sn.isstateful.set(i, 0, 1.0); //state dependent service
            }
        }
        sn.isstatedep = new Matrix(sn.nnodes, 3, 3 * sn.nnodes);
        for (int i = 0; i < sn.nnodes; i++) {
            Node node = this.nodes.get(i);
            if (node instanceof Cache) {
                sn.isstatedep.set(i, 1, 1.0); //state dependent service
            }
            for (int j = 0; j < sn.nclasses; j++) {
                JobClass jobclass = this.jobClasses.get(j);
                switch (sn.routing.get(node).get(jobclass)) {
                    case RROBIN:
                    case WRROBIN:
                    case JSQ:
                    case RL:
                    case KCHOICES:
                        sn.isstatedep.set(i, 2, 1.0); // state dependent routing
                        break;
                    default:
                        continue;
                }
            }
        }
        sn.nstateful = getNumberOfStatefulNodes();
        sn.state = new HashMap<StatefulNode, Matrix>(sn.nstations);
        sn.stateprior = new HashMap<StatefulNode, Matrix>(sn.nstations);
        sn.space = new HashMap<StatefulNode, Matrix>(sn.nstations);
        
        // Populate state maps from the actual node states
        for (int i = 0; i < sn.nstateful; i++) {
            StatefulNode node = sn.stateful.get(i);
            Matrix nodeState = node.getState();
            Matrix nodePrior = node.getStatePrior();
            Matrix nodeSpace = node.getStateSpace();
            
            // Use the actual state from the node, or empty if not initialized
            // Note: getState() returns a copy, which is what we want
            sn.state.put(node, nodeState != null && !nodeState.isEmpty() ? nodeState : new Matrix(0, 0, 0));
            sn.stateprior.put(node, nodePrior != null && !nodePrior.isEmpty() ? nodePrior : new Matrix(0, 0, 0));
            sn.space.put(node, nodeSpace != null && !nodeSpace.isEmpty() ? nodeSpace : new Matrix(0, 0, 0));
        }
        sn.nodenames = nodenames;
        sn.classnames = classnames;
        sn.connmatrix = conn;

        //line 97-108 is ignored since for transition node

        sn.nodeToStateful = new Matrix(1, nodes.size(), nodes.size());
        sn.nodeToStation = new Matrix(1, nodes.size(), nodes.size());
        sn.stationToNode = new Matrix(1, stations.size(), stations.size());
        sn.stationToStateful = new Matrix(1, stations.size(), stations.size());
        sn.statefulToStation = new Matrix(1, sn.nstateful, sn.nstateful);
        sn.statefulToNode = new Matrix(1, sn.nstateful, sn.nstateful);
        for (int i = 0; i < nodes.size(); i++) {
            sn.nodeToStateful.set(0, i, nodes.get(i).getStatefulIdx());
            sn.nodeToStation.set(0, i, nodes.get(i).getStationIdx());
        }
        for (int i = 0; i < stations.size(); i++) {
            sn.stationToNode.set(0, i, stations.get(i).getNodeIndex());
            sn.stationToStateful.set(0, i, stations.get(i).getStatefulIdx());
        }
        for (int isf = 0; isf < sn.nstateful; isf++) {
            sn.statefulToNode.set(0, isf, getStatefulNodeFromIndex(isf).getNodeIndex());
            sn.statefulToStation.set(0, isf, getStatefulNodeFromIndex(isf).getStationIdx());
        }

        // Initialize immediate feedback matrix (station x class)
        sn.immfeed = new Matrix(sn.nstations, sn.nclasses, sn.nstations * sn.nclasses);
        for (int ist = 0; ist < sn.nstations; ist++) {
            int nodeIdx = (int) sn.stationToNode.get(ist);
            Node node = this.nodes.get(nodeIdx);
            for (int c = 0; c < sn.nclasses; c++) {
                boolean stationHas = false;
                if (node instanceof jline.lang.nodes.Queue) {
                    jline.lang.nodes.Queue queue = (jline.lang.nodes.Queue) node;
                    stationHas = queue.hasImmediateFeedback(c);
                }
                boolean classHas = this.jobClasses.get(c).hasImmediateFeedback();
                if (stationHas || classHas) {
                    sn.immfeed.set(ist, c, 1.0);
                }
            }
        }

        refreshPriorities();
        refreshDeadlines();
        refreshProcesses(null, null);
        refreshImpatience();
        refreshBalking();
        refreshRetrial();
        refreshHeterogeneousServers();

        // Check if priorities are specified but no priority-aware scheduling policy is used
        // Priority check: non-uniform priorities indicate priority-based scheduling
        boolean hasPriorities = false;
        double firstPrio = sn.classprio.get(0, 0);
        for (int i = 1; i < sn.classprio.getNumCols(); i++) {
            if (sn.classprio.get(0, i) != firstPrio) {
                hasPriorities = true;
                break;
            }
        }

        if (hasPriorities && sn.sched != null) {
            // Priority classes exist, check if any station uses priority-aware scheduling
            boolean hasPriorityScheduling = false;
            for (SchedStrategy schedStrategy : sn.sched.values()) {
                if (schedStrategy == SchedStrategy.PSPRIO ||
                    schedStrategy == SchedStrategy.DPSPRIO ||
                    schedStrategy == SchedStrategy.GPSPRIO ||
                    schedStrategy == SchedStrategy.HOL ||
                    schedStrategy == SchedStrategy.FCFSPRIO ||
                    schedStrategy == SchedStrategy.LCFSPRIO ||
                    schedStrategy == SchedStrategy.LCFSPRPRIO ||
                    schedStrategy == SchedStrategy.LCFSPIPRIO ||
                    schedStrategy == SchedStrategy.FCFSPRPRIO ||
                    schedStrategy == SchedStrategy.FCFSPIPRIO ||
                    schedStrategy == SchedStrategy.SRPTPRIO) {
                    hasPriorityScheduling = true;
                    break;
                }
            }

            if (!hasPriorityScheduling) {
                line_warning(mfilename(new Object() {}),
                    "Priority classes are specified but no priority-aware scheduling policy (PSPRIO, DPSPRIO, GPSPRIO, HOL, FCFSPRIO, FCFSPRPRIO, FCFSPIPRIO, LCFSPRIO, LCFSPRPRIO, LCFSPIPRIO, SRPTPRIO) is used in the model. Priorities will be ignored.");
            } else if (GlobalConstants.Verbose != VerboseLevel.SILENT) {
                // Display priority info
                double minPrio = Double.MAX_VALUE;
                double maxPrio = Double.MIN_VALUE;
                for (int i = 0; i < sn.classprio.getNumCols(); i++) {
                    double p = sn.classprio.get(0, i);
                    if (p < minPrio) minPrio = p;
                    if (p > maxPrio) maxPrio = p;
                }
                StringBuilder highNames = new StringBuilder();
                StringBuilder lowNames = new StringBuilder();
                for (int i = 0; i < sn.classprio.getNumCols(); i++) {
                    double p = sn.classprio.get(0, i);
                    if (p == minPrio) {
                        if (highNames.length() > 0) highNames.append(",");
                        highNames.append(sn.classnames.get(i));
                    }
                    if (p == maxPrio) {
                        if (lowNames.length() > 0) lowNames.append(",");
                        lowNames.append(sn.classnames.get(i));
                    }
                }
                System.out.printf("Priority: highest=%s, lowest=%s%n", highNames, lowNames);
            }
        }

        refreshChains(!sn.nodetype.contains(NodeType.Cache));

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
        refreshPetriNetNodes();
        refreshSync();
        refreshGlobalSync();
        refreshRegions();
        this.hasStruct = true;

        if (this.sn.fj.any()) {
            Matrix forkLambda = new Matrix(1, this.sn.nclasses).fill(GlobalConstants.FineTol);
            Ret.FJApprox approxReturn = mmt(this, forkLambda);
            Network nonfjmodel = approxReturn.nonfjmodel;
            Map<Integer, Integer> fanout = approxReturn.fanout;
            Matrix forkmap = approxReturn.fjforkmap;
            Matrix fjclassmap = approxReturn.fjclassmap;
            if (fanout.values().stream().anyMatch(value -> value == 1)) { //if any(fanOut==1)
                line_warning(mfilename(new Object() {
                }), "The specified fork-join topology has partial support, only SolverJMT simulation results may be reliable.\n");
            }
            NetworkStruct fsn = nonfjmodel.getStruct();
            Matrix[] origNodeVisits = new Matrix[sn.nodevisits.size()];
            for (int oldChain = 0; oldChain < sn.nodevisits.size(); oldChain++) {
                origNodeVisits[oldChain] = sn.nodevisits.get(oldChain).copy();
            }
            for (int newChain = sn.nchains; newChain < fsn.nchains; newChain++) {
                // Find an auxiliary class in this chain (one with valid mapping in fjclassmap)
                // The chain may contain both original and auxiliary classes due to chain detection grouping them together
                int anyAuxClass = -1;
                Matrix chainClasses = fsn.inchain.get(newChain);
                for (int col = 0; col < chainClasses.getNumCols(); col++) {
                    int classIdx = (int) chainClasses.get(col);
                    // Check if this class has a valid mapping in fjclassmap (auxiliary classes have non-negative mappings)
                    if (classIdx < fjclassmap.getNumCols() && fjclassmap.get(0, classIdx) >= 0) {
                        anyAuxClass = classIdx;
                        break;
                    }
                }
                if (anyAuxClass < 0) {
                    continue; // No auxiliary class found in this chain, skip processing
                }
                int origFork = (int) forkmap.get(0, anyAuxClass);
                int origClassIdx = (int) fjclassmap.get(0, anyAuxClass);
                int origChain = (int) sn.chains.getColumn(origClassIdx).find().value();
                HashSet<Integer> rowsIdx = new HashSet<>();
                for (int row = 0; row < fsn.nodetype.size(); row++) {
                    if (fsn.nodetype.get(row) == NodeType.Source || fsn.nodetype.get(row) == NodeType.Sink || fsn.nodetype.get(row) == NodeType.Fork) {
                        for (int col = 0; col < fsn.nodevisits.get(newChain).getNumCols(); col++)
                            fsn.nodevisits.get(newChain).set(row, col, 0.0);
                    }
                }
                Matrix Vaux = fsn.nodevisits.get(newChain);
                // inchain contains 0-indexed class numbers, use directly for column extraction
                Collection<Integer> colsAuxIdx = Arrays.stream(fsn.inchain.get(newChain).toIntArray1D())
                        .boxed()
                        .collect(Collectors.toList());
                Vaux.keepCols(colsAuxIdx);

                if (fsn.nnodes != sn.nnodes) {
                    // Build mapping from fsn nodes to sn nodes by name
                    // This handles cases where ClassSwitch/Source/Sink were added in fsn
                    Map<String, Integer> snNodeNameToIdx = new HashMap<>();
                    for (int i = 0; i < sn.nnodes; i++) {
                        snNodeNameToIdx.put(sn.nodenames.get(i), i);
                    }

                    // Create a new matrix with sn.nnodes rows, mapping fsn rows to sn positions
                    Matrix VauxMapped = new Matrix(sn.nnodes, Vaux.getNumCols());
                    for (int fsnRow = 0; fsnRow < fsn.nnodes; fsnRow++) {
                        String nodeName = fsn.nodenames.get(fsnRow);
                        Integer snRow = snNodeNameToIdx.get(nodeName);
                        if (snRow != null) {
                            // This fsn node exists in sn - copy its visit data
                            for (int col = 0; col < Vaux.getNumCols(); col++) {
                                VauxMapped.set(snRow, col, Vaux.get(fsnRow, col));
                            }
                        }
                        // Nodes not in sn (ClassSwitch, Source, Sink) are skipped
                    }
                    Vaux = VauxMapped;
                }
                int fanOut = (int) ((ForkNodeParam) sn.nodeparam.get(this.getNodeByIndex(origFork))).fanOut;
                Matrix X = origNodeVisits[origChain].scale(fanOut);
                for (int jaux = 0; jaux < fsn.inchain.get(newChain).getNumCols(); jaux++) {
                    // inchain contains 0-indexed class numbers, use directly
                    int j = (int) sn.inchain.get(origChain).get(jaux);
                    Matrix Y = Vaux.getColumn(jaux).scale(fanOut);
                    sn.nodevisits.get(origChain).setColumn(j, X.getColumn(j).add(Y));
                }
            }
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
                    rtmask.set(rowIdx, colIdx, FastMath.ceil(value));
                }
            }
        }

        for (int i = 0; i < sn.nnodes; i++) {
            for (int r = 0; r < nclasses; r++) {
                if (sn.isstation.get(i, 0) > 0 && sn.phases.get((int) sn.nodeToStation.get(0, i), r) > 1) {
                    Sync synct = new Sync();
                    synct.active.put(0, new Event(EventType.PHASE, i, r, NaN, new Matrix(0, 0), NaN, NaN));
                    synct.passive.put(0, new Event(EventType.LOCAL, local, r, 1.0, new Matrix(0, 0), NaN, NaN));
                    sync.put(sync.size(), synct);
                }
                if (sn.isstateful.get(i, 0) > 0) {
                    if (sn.nodetype.get(i).equals(NodeType.Cache)) {
                        if (((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(i))).pread.get(r) != null) {
                            Sync synct = new Sync();
                            synct.active.put(0, new Event(EventType.READ, i, r, NaN, new Matrix(0, 0), NaN, NaN));
                            synct.passive.put(0, new Event(EventType.READ, local, r, 1.0, new Matrix(0, 0), NaN, NaN));
                            sync.put(sync.size(), synct);
                        }
                    } else if (sn.nodetype.get(i).equals(NodeType.Transition)) {
                        // For Transitions, create sync entries for each mode (not class)
                        // This matches MATLAB: for m=1:sn.nodeparam{ind}.nmodes
                        if (r == 0) { // Only do this once per node, not for each class
                            TransitionNodeParam transParam = (TransitionNodeParam) sn.nodeparam.get(sn.nodes.get(i));
                            int nmodes = transParam.nmodes;
                            for (int m = 0; m < nmodes; m++) {
                                Sync synct = new Sync();
                                synct.active.put(0, new Event(EventType.PHASE, i, m, NaN, new Matrix(0, 0), NaN, NaN));
                                synct.passive.put(0, new Event(EventType.LOCAL, local, m, 1.0, new Matrix(0, 0), NaN, NaN));
                                sync.put(sync.size(), synct);
                            }
                        }
                    }
                    int isf = (int) sn.nodeToStateful.get(0, i);
                    for (int j = 0; j < sn.nnodes; j++) {
                        if (sn.isstateful.get(j, 0) > 0) {
                            int jsf = (int) sn.nodeToStateful.get(0, j);
                            for (int s = 0; s < nclasses; s++) {
                                double p = rtmask.get(isf * nclasses + r, jsf * nclasses + s);
                                if (p > 0) {
                                    Sync synct = new Sync();
                                    synct.active.put(0, new Event(EventType.DEP, i, r, NaN, new Matrix(0, 0), NaN, NaN));
                                    switch (sn.routing.get(this.nodes.get(i)).get(this.jobClasses.get(s))) {
                                        case RROBIN:
                                        case WRROBIN:
                                        case JSQ:
                                            final int i_final = i, j_final = j, r_final = r, s_final = s;
                                            synct.passive.put(0, new Event(EventType.ARV, j, s, ((pair) -> sn.rtfun.apply(pair).get(i_final * nclasses + r_final, j_final * nclasses + s_final)), new Matrix(0, 0), NaN, NaN));
                                            break;
                                        default:
                                            synct.passive.put(0, new Event(EventType.ARV, j, s, sn.rt.get(isf * nclasses + r, jsf * nclasses + s), new Matrix(0, 0), NaN, NaN));
                                    }
                                    sync.put(sync.size(), synct);
                                }
                            }
                        }
                    }
                }
            }
        }

        if (this.sn != null) this.sn.sync = sync;
    }

    public void refreshGlobalSync() {
        int nclasses = this.sn.nclasses;
        Map<Integer, GlobalSync> gsync = new HashMap<Integer, GlobalSync>();
        Map<Node, Matrix> emptystate = new HashMap<Node, Matrix>();
        for (Node node : this.nodes)
            emptystate.put(node, new Matrix(0, 0));

        Matrix rtmask;
        if (this.sn.isstatedep.getNonZeroLength() > 0) {
            rtmask = this.sn.rtfun.apply(new Pair<Map<Node, Matrix>, Map<Node, Matrix>>(emptystate, emptystate));
        } else {
            Matrix rt = this.sn.rt;
            rtmask = new Matrix(rt.getNumRows(), rt.getNumCols());
            for (int colIdx = 0; colIdx < rt.getNumCols(); colIdx++) {
                int col1 = rt.getColIndexes()[colIdx];
                int col2 = rt.getColIndexes()[colIdx + 1];

                for (int i = col1; i < col2; i++) {
                    int rowIdx = rt.getNonZeroRows()[i];
                    double value = rt.getNonZeroValues()[i];
                    rtmask.set(rowIdx, colIdx, FastMath.ceil(value));
                }
            }
        }

        for (int ind = 0; ind < sn.nnodes; ind++) {
            for (int r = 0; r < nclasses; r++) {
                if (sn.isstateful.get(ind, 0) > 0) {
                    if (sn.nodetype.get(ind) == NodeType.Transition) {
                        TransitionNodeParam transParam = (TransitionNodeParam) sn.nodeparam.get(sn.nodes.get(ind));
                        for (int m = 0; m < transParam.nmodes; m++) {
                            // Mode enabling
                            Matrix enablingMatrix = transParam.enabling.get(m);
                            List<Integer> enablingPlaces = new ArrayList<>();
                            for (int ep = 0; ep < enablingMatrix.getNumRows(); ep++) {
                                if (enablingMatrix.get(ep, 0) > 0) {
                                    enablingPlaces.add(ep);
                                }
                            }

                            if (!enablingPlaces.isEmpty()) {
                                GlobalSync enableSync = new GlobalSync();
                                List<ModeEvent> activeEvents = new ArrayList<>();
                                activeEvents.add(new ModeEvent(EventType.ENABLE, ind, m, 1.0));
                                enableSync.setActive(activeEvents);

                                List<ModeEvent> passiveEvents = new ArrayList<>();
                                for (int ep : enablingPlaces) {
                                    passiveEvents.add(new ModeEvent(EventType.LOCAL, ep, m, 1.0));
                                }
                                enableSync.setPassive(passiveEvents);
                                gsync.put(gsync.size(), enableSync);
                            }

                            // Mode firing  
                            Matrix firingMatrix = transParam.firing.get(m);
                            List<Integer> firingPlaces = new ArrayList<>();
                            for (int fp = 0; fp < firingMatrix.getNumRows(); fp++) {
                                if (firingMatrix.get(fp, 0) > 0) {
                                    firingPlaces.add(fp);
                                }
                            }

                            if (!enablingPlaces.isEmpty() || !firingPlaces.isEmpty()) {
                                GlobalSync fireSync = new GlobalSync();
                                List<ModeEvent> activeEvents = new ArrayList<>();
                                activeEvents.add(new ModeEvent(EventType.FIRE, ind, m));
                                fireSync.setActive(activeEvents);

                                List<ModeEvent> passiveEvents = new ArrayList<>();
                                for (int ep : enablingPlaces) {
                                    double enablingValue = transParam.enabling.get(m).get(ep, 0);
                                    passiveEvents.add(new ModeEvent(EventType.PRE, ep, m, enablingValue));
                                }
                                for (int fp : firingPlaces) {
                                    double firingValue = transParam.firing.get(m).get(fp, 0);
                                    passiveEvents.add(new ModeEvent(EventType.POST, fp, m, firingValue));
                                }
                                fireSync.setPassive(passiveEvents);
                                gsync.put(gsync.size(), fireSync);
                            }
                        }
                    }
                }
            }
        }

        if (this.sn != null) this.sn.gsync = gsync;
    }

    // ========================================================================
    // SECTION 19: RESET METHODS
    // Methods for resetting and relinking network components
    // ========================================================================

    public void relink(RoutingMatrix P) {
        resetNetwork();
        link(P);
    }

    public void reset() {
        this.resetModel(false);
        //this.resetNetwork();
        this.hasState = false;
    }

    public void reset(boolean resetState) {
        this.resetModel(resetState);
        this.hasState = false;
    }

    public void resetHandles() {
        this.handles = new ArrayList<>();
    }

    /**
     * Check if the network contains any Logger nodes
     * @return true if Logger nodes exist in the network
     */
    public boolean hasExistingLoggers() {
        for (Node node : this.nodes) {
            if (node instanceof Logger) {
                return true;
            }
        }
        return false;
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

        // Remove class switch nodes and logger nodes
        if (deleteCSNodes) {
            List<Node> oldNodes = this.nodes;
            this.nodes = new ArrayList<>();
            for (Node n : oldNodes) {
                if (!(n instanceof ClassSwitch) && !(n instanceof Logger)) {
                    // Reset cached node index since position in list may change
                    n.setNodeIdx(-1);
                    this.nodes.add(n);
                }
            }
        }

        for (int i = 0; i < M; i++) {
            ((Dispatcher) this.stations.get(i).getOutput()).initDispatcherJobClasses(this.getClasses());
        }

        this.handles = new ArrayList<>();
        this.connections = null;
        return this.getNodes();
    }

    /**
     * Resets the struct of a given network
     */
    public void resetStruct() {
        this.sn = null;
        this.hasStruct = false;
    }

    // ========================================================================
    // SECTION 20: VALIDATION AND CONFIGURATION
    // Methods for validating and configuring network parameters
    // ========================================================================

    public void sanitize() {
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
                            cache.getCacheServer().hitClass.set(r, c, FastMath.round(cache.getCacheServer().hitClass.get(r, c)));
                        }
                    }
                    for (int r = 0; r < cache.getCacheServer().missClass.getNumRows(); r++) {
                        for (int c = 0; c < cache.getCacheServer().missClass.getNumCols(); c++) {
                            cache.getCacheServer().missClass.set(r, c, FastMath.round(cache.getCacheServer().missClass.get(r, c)));
                        }
                    }
                } else if (node instanceof Logger) {
                    //no-op
                } else if (node instanceof ClassSwitch) {
                    //no-op
                } else if (node instanceof Join) {
                    Join join = (Join) node;
                    for (int k = 0; k < K; k++) {
                        JobClass jobclass = this.jobClasses.get(k);
                        join.setClassCap(jobclass, Integer.MAX_VALUE);
                        join.setDropRule(jobclass, DropStrategy.WaitingQueue);
                    }
                } else if (node instanceof Delay) {
                    Delay delay = (Delay) node;
                    boolean hasEnabledService = false;
                    for (int k = 0; k < K; k++) {
                        JobClass jobclass = this.jobClasses.get(k);
                        if (!delay.getServer().containsJobClass(jobclass)) {
                            delay.setService(jobclass, new Disabled(), 0);
                            delay.setClassCap(jobclass, 0);
                        } else {
                            // Check if this class has an enabled service
                            if (!(delay.getServiceProcess(jobclass) instanceof Disabled)) {
                                hasEnabledService = true;
                            }
                        }
                    }
                    // Validate that at least one class has service configured
                    // Skip validation for models with Cache, Petri net elements, or Source
                    // (Source indicates an open model from async calls where Delay may only have disabled services)
                    boolean hasSpecialElements = false;
                    for (int j = 0; j < this.nodes.size(); j++) {
                        Node nodeCheck = this.nodes.get(j);
                        if (nodeCheck instanceof Cache || nodeCheck instanceof Place || nodeCheck instanceof Transition || nodeCheck instanceof Source) {
                            hasSpecialElements = true;
                            break;
                        }
                    }
                    if (!hasEnabledService && !hasSpecialElements) {
                        line_error(mfilename(new Object() {}), "Delay '" + delay.getName() + "' has no service configured for any job class. Use setService() to configure service times.");
                    }
                    for (int k = 0; k < K; k++) {
                        JobClass jobclass = this.jobClasses.get(k);
                        if (delay.getServiceProcess(jobclass) instanceof Disabled)
                            delay.setRouting(jobclass, RoutingStrategy.DISABLED);
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
                    boolean hasEnabledService = false;
                    for (int k = 0; k < K; k++) {
                        JobClass jobclass = this.jobClasses.get(k);
                        if (!queue.getServer().containsJobClass(jobclass)) {
                            queue.setService(jobclass, new Disabled(), 0);
                            queue.setClassCap(jobclass, 0);
                        } else {
                            // Check if this class has an enabled service
                            if (!(queue.getServiceProcess(jobclass) instanceof Disabled)) {
                                hasEnabledService = true;
                            }
                        }
                    }
                    // Validate that at least one class has service configured
                    // Skip validation for models with Cache, Petri net elements, or Source
                    // (Source indicates an open model from async calls where Queue may only have disabled services)
                    boolean hasSpecialElements = false;
                    for (int j = 0; j < this.nodes.size(); j++) {
                        Node nodeCheck = this.nodes.get(j);
                        if (nodeCheck instanceof Cache || nodeCheck instanceof Place || nodeCheck instanceof Transition || nodeCheck instanceof Source) {
                            hasSpecialElements = true;
                            break;
                        }
                    }
                    if (!hasEnabledService && !hasSpecialElements) {
                        line_error(mfilename(new Object() {}), "Queue '" + queue.getName() + "' has no service configured for any job class. Use setService() to configure service times.");
                    }

                    for (int k = 0; k < K; k++) {
                        JobClass jobclass = this.jobClasses.get(k);
                        if (queue.getServiceProcess(jobclass) instanceof Disabled)
                            queue.setRouting(jobclass, RoutingStrategy.DISABLED);
                    }

                    // Declare variables outside switch to avoid scope issues
                    ArrayList<Double> svcTime;
                    ArrayList<Double> svcTimeSorted;
                    
                    switch (queue.getSchedStrategy()) {
                        case SEPT:
                            svcTime = new ArrayList<Double>();
                            for (int k = 0; k < K; k++)
                                svcTime.add(queue.getServiceProcess(this.jobClasses.get(k)).getMean());

                            if (svcTime.stream().distinct().collect(Collectors.toList()).size() != K)
                                throw new RuntimeException("SEPT does not support identical service time means.");

                            svcTimeSorted = new ArrayList<Double>(svcTime);
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
                    for (int r = 0; r < this.getNumberOfClasses(); r++) {
                        node.setRouting(this.getClassByIndex(r), RoutingStrategy.DISABLED);
                    }
                } else if (node instanceof Source) {
                    Source source = (Source) node;
                    for (int k = 0; k < K; k++) {
                        JobClass jobclass = this.jobClasses.get(k);
                        if (!source.containsJobClass(jobclass)) {
                            source.setArrival(jobclass, new Disabled());
                            source.setRouting(jobclass, RoutingStrategy.DISABLED);
                        }
                    }
                } else if (node instanceof Place) {
                    Place place = (Place) node;
                    for (int k = 0; k < K; k++) {
                        JobClass jobclass = this.jobClasses.get(k);
                        place.setDropRule(jobclass, DropStrategy.Drop);
                    }
                } else if (node instanceof Transition) {
                    Transition transition = (Transition) node;
                    for (int k = 0; k < K; k++) {
                        JobClass jobclass = this.jobClasses.get(k);
                        if (jobclass != null) transition.setService(jobclass, new Disabled());
                    }
                }
            }

            int sourceIdx = this.getIndexSourceStation();
            for (int i = 0; i < M; i++) {
                if ((sourceIdx == -1) || (i != sourceIdx)) {
                    for (int r = 0; r < K; r++) {
                        ServiceSection server = this.stations.get(i).getServer();
                        if (server instanceof ServiceTunnel) {
                            //do nothing
                        } //else if (server instanceof CacheClassSwitcher) {}
                        else {
                            if (this.stations.get(i).getServer() != null) {
                                if (!this.stations.get(i).getServer().containsJobClass(this.jobClasses.get(r)))
                                    this.stations.get(i).getServer().setServiceProcesses(new ServiceBinding(this.jobClasses.get(r), ServiceStrategy.LI, new Disabled()));
                            }
                        }

                    }
                }
            }

            // Check if model has Cache, Petri net elements, or Source
            // (Source indicates an open model from async calls)
            boolean hasSpecialElements = false;
            for (int i = 0; i < this.nodes.size(); i++) {
                Node node = this.nodes.get(i);
                if (node instanceof Cache || node instanceof Place || node instanceof Transition || node instanceof Source) {
                    hasSpecialElements = true;
                    break;
                }
            }

            // Validate that each job class has service configured at least one station
            // Skip validation for models with Cache, Petri net elements, or Source
            if (!hasSpecialElements) {
                for (int k = 0; k < K; k++) {
                    boolean hasServiceAtAnyStation = false;
                    boolean isUsedInCacheSetRead = false;

                    // Check if job class is used in setRead at a Cache node
                    for (int i = 0; i < this.nodes.size(); i++) {
                        Node node = this.nodes.get(i);
                        if (node instanceof Cache) {
                            Cache cache = (Cache) node;
                            if (cache.getCacheServer().inputJobClasses.containsKey(this.jobClasses.get(k).getIndex())) {
                                isUsedInCacheSetRead = true;
                                break;
                            }
                        }
                    }

                    for (int i = 0; i < M; i++) {
                        ServiceSection server = this.stations.get(i).getServer();
                        if (server != null && server.containsJobClass(this.jobClasses.get(k))) {
                            // Check if this is a valid service configuration (not disabled)
                            hasServiceAtAnyStation = true;
                            break;
                        }
                    }

                    // Check all classes now (removed open class exemption)
                    if (!hasServiceAtAnyStation && !isUsedInCacheSetRead) {
                        line_error(mfilename(new Object() {}), "Job class '" + this.jobClasses.get(k).getName() + "' has no service configured at any station. Every job class must have service configured at least one station using setService().");
                    }
                }
            }
        }
    }

    /**
     * Resolves Signal placeholders to OpenSignal or ClosedSignal based on network structure.
     *
     * <p>This method is called during model finalization (refreshStruct) to convert
     * Signal placeholder objects to their concrete types. For open networks (with Source),
     * Signal becomes OpenSignal. For closed networks (no Source), Signal becomes ClosedSignal.
     */
    public void resolveSignals() {
        // Check if there are any Signal placeholders to resolve
        boolean hasSignals = false;
        for (int i = 0; i < this.jobClasses.size(); i++) {
            JobClass jc = this.jobClasses.get(i);
            if (jc instanceof Signal && !(jc instanceof OpenSignal) && !(jc instanceof ClosedSignal)) {
                hasSignals = true;
                break;
            }
        }
        if (!hasSignals) {
            return;  // No Signal placeholders to resolve
        }

        // Determine if the network is open (has Source node)
        boolean isOpen = (getIndexSourceNode() != -1);

        // For closed networks, find a default reference station
        Station defaultRefstat = null;
        if (!isOpen) {
            // Look for first Delay node, then first Queue node
            for (Station station : this.stations) {
                if (station instanceof Delay) {
                    defaultRefstat = station;
                    break;
                }
            }
            if (defaultRefstat == null) {
                for (Station station : this.stations) {
                    if (station instanceof Queue) {
                        defaultRefstat = station;
                        break;
                    }
                }
            }
            if (defaultRefstat == null && !this.stations.isEmpty()) {
                defaultRefstat = this.stations.get(0);
            }
        }

        // Resolve each Signal placeholder
        for (int i = 0; i < this.jobClasses.size(); i++) {
            JobClass jc = this.jobClasses.get(i);
            if (jc instanceof Signal && !(jc instanceof OpenSignal) && !(jc instanceof ClosedSignal)) {
                Signal sig = (Signal) jc;

                // Determine reference station for closed signals
                Station refstat = defaultRefstat;
                if (!isOpen && sig.getTargetJobClass() != null) {
                    // Prefer the reference station of the targetJobClass
                    refstat = sig.getTargetJobClass().getReferenceStation();
                    if (refstat == null) {
                        refstat = defaultRefstat;
                    }
                }

                // Resolve and replace
                JobClass concrete = sig.resolve(isOpen, refstat);
                this.jobClasses.set(i, concrete);
            }
        }
    }

    /**
     * Enables or disables validation checks for this network.
     *
     * @param doChecks true to enable validation checks, false to disable
     */
    public void setChecks(boolean doChecks) {
        this.enableChecks = doChecks;
    }

    /**
     * Sets the initialization status of this network.
     *
     * @param initStatus true if the network is initialized, false otherwise
     */
    public void setInitialized(boolean initStatus) {
        this.hasState = initStatus;
    }

    public void setJoinNodeRequired(int nodeIdx, JobClass jobClass, int njobs) {
        Join join = (Join) this.nodes.get(nodeIdx);
        join.setRequired(jobClass, njobs);
        this.nodes.set(nodeIdx, join);
    }

    public void setJoinNodeStrategy(int nodeIdx, JobClass jobClass, JoinStrategy joinStrategy) {
        Join join = (Join) this.nodes.get(nodeIdx);
        join.setStrategy(jobClass, joinStrategy);
        this.nodes.set(nodeIdx, join);
    }

    public void setNodeRouting(int nodeIdx, JobClass jobClass, RoutingStrategy routingStrategy) {
        Node node = this.nodes.get(nodeIdx);
        node.setRouting(jobClass, routingStrategy);
        this.nodes.set(nodeIdx, node);
    }

    public void setSn(NetworkStruct sn) {
        this.sn = sn;
    }

    public void setUsedLangFeature(String feature) {
        usedFeatures.setTrue(feature);
    }

    // ========================================================================
    // SECTION 21: UTILITY METHODS
    // Helper and utility methods for various network operations
    // ========================================================================

    public double sub_jsq(int ind, int jnd, int r, int s, Matrix
            linksmat, Map<Node, Matrix> state_before, Map<Node, Matrix> state_after) {
        int isf = (int) this.sn.nodeToStateful.get(ind);
        Node statefulNode = this.getStatefulNodeFromIndex(isf);
        if (!state_before.containsKey(statefulNode)) {
            return FastMath.min(linksmat.get(ind, jnd), 1.0);
        } else {
            if (r == s) {
                Matrix n = new Matrix(1, this.sn.nnodes);
                n.fill(Inf);
                for (int knd = 0; knd < this.sn.nnodes; knd++) {
                    if (linksmat.get(ind, knd) > 0) {
                        Node statefulNode_knd = this.getStatefulNodeFromIndex((int) this.sn.nodeToStateful.get(0, knd));
                        n.set(0, knd, ToMarginal.toMarginal(this.sn, knd, state_before.get(statefulNode), null, null, null, null, null).ni.value());
                    }
                }
                double min = n.elementMin();
                if (n.get(jnd) == min) return 1.0 / n.count(min);
                else return 0.0;
            } else {
                return 0.0;
            }
        }
    }

    public double sub_rr_wrr(int ind, int jnd, int r, int s, Matrix linksmat, Map<Node, Matrix> state_before, Map<Node, Matrix> state_after) {
        int R = this.sn.nclasses;
        int isf = (int) this.sn.nodeToStateful.get(ind);
        Node statefulNode = this.getStatefulNodeFromIndex(isf);
        if (!state_before.containsKey(statefulNode)) {
            return FastMath.min(linksmat.get(ind, jnd), 1.0);
        } else {
            if (r == s) {
                Matrix jm = state_after.get(statefulNode);
                return ((int) jm.get(jm.getNumCols() - 1 - R + r) == jnd) ? 1.0 : 0.0;
            } else {
                return 0.0;
            }
        }
    }

    public void summary() {
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

    public void unLink() {
        for (Node node : this.nodes) {
            node.resetRouting();
        }
    }

    /**
     * Validates that routing probabilities don't exceed 1.0 for any node/class combination
     *
     * @param P the routing matrix to validate
     */
    private void validateRoutingProbabilities(RoutingMatrix P) {
        // Check routing matrix for valid probability sums
        for (int r = 0; r < jobClasses.size(); r++) {
            JobClass jobClass = jobClasses.get(r);
            for (int i = 0; i < nodes.size(); i++) {
                Node node = nodes.get(i);
                // Only validate stations, and exclude FORK stations from validation
                if (node instanceof Station) {
                    Station station = (Station) node;
                    if (station.getSchedStrategy() != SchedStrategy.FORK) {
                        double totalProb = 0.0;

                        // Sum routing probabilities from this node to all other nodes for this class
                        for (int j = 0; j < nodes.size(); j++) {
                            // Sum across all possible destination classes for class switching
                            for (int s = 0; s < jobClasses.size(); s++) {
                                double prob = P.get(r, s, i, j);
                                totalProb += prob;
                            }
                        }

                        // Check if total probability exceeds 1.0 (with tolerance)
                        if (totalProb > 1.0 + GlobalConstants.FineTol) {
                            line_error(mfilename(new Object() {
                                    }),
                                    String.format("The total routing probability for jobs leaving node %s in class %s is greater than 1.0 (%.6f).",
                                            node.getName(), jobClass.getName(), totalProb));
                        }
                    }
                }
            }
        }
    }

    public void view() {
        jsimgView();
    }

    public void modelView() {
        plot();
    }

    // ========================================================================
    // SECTION 22: INNER CLASSES
    // Inner classes for specialized return types and data structures
    // ========================================================================

    public static class routingMatrixReturn {

        public Matrix rt;
        public Matrix rtnodes;
        public Matrix linksmat;
        public Matrix chains;
        public Map<JobClass, Map<JobClass, Matrix>> rtNodesByClass;
        public Map<Node, Map<Node, Matrix>> rtNodesByStation;

        public routingMatrixReturn(Matrix rt, Matrix rtnodes, Matrix linksmat, Matrix chains, Map<JobClass, Map<JobClass, Matrix>> rtNodesByClass, Map<Node, Map<Node, Matrix>> rtNodesByStation) {
            this.rt = rt;
            this.rtnodes = rtnodes;
            this.linksmat = linksmat;
            this.chains = chains;
            this.rtNodesByClass = rtNodesByClass;
            this.rtNodesByStation = rtNodesByStation;
        }
    }

    // ========================================================================
    // SECTION 22: REWARD COMPUTATION
    // Methods for defining and managing reward functions for CTMC analysis
    // ========================================================================

    /**
     * Define a reward function for CTMC reward computation.
     *
     * The reward function maps a state vector and network structure to a scalar
     * reward value. Multiple rewards can be defined with different names.
     *
     * @param name The unique name for this reward
     * @param rewardFn The reward function: (state, sn) -> double
     *
     * Example:
     * <pre>
     * // Queue length reward
     * model.setReward("qlen", (state, sn) -> state.get(0, 1));
     *
     * // Throughput reward
     * model.setReward("throughput", (state, sn) -> {
     *     double n = state.get(0, 1);
     *     return n > 0 ? n * sn.rates.get(1, 0) / n : 0;
     * });
     * </pre>
     */
    public void setReward(String name, RewardFunction rewardFn) {
        if (this.sn == null) {
            this.sn = new NetworkStruct();
        }
        if (this.sn.reward == null) {
            this.sn.reward = new HashMap<String, RewardFunction>();
        }
        this.sn.reward.put(name, rewardFn);
    }

    /**
     * Get all defined reward functions.
     *
     * @return Map from reward name to reward function, or null if no rewards defined
     */
    public Map<String, RewardFunction> getRewards() {
        if (this.sn != null) {
            return this.sn.reward;
        }
        return null;
    }

    /**
     * Get a specific reward function by name.
     *
     * @param name The reward name
     * @return The reward function, or null if not found
     */
    public RewardFunction getReward(String name) {
        if (this.sn != null && this.sn.reward != null) {
            return this.sn.reward.get(name);
        }
        return null;
    }

    /**
     * Remove all defined reward functions.
     */
    public void clearRewards() {
        if (this.sn != null) {
            this.sn.reward = null;
        }
    }

    /**
     * Check if any rewards are defined.
     *
     * @return true if at least one reward is defined
     */
    public boolean hasRewards() {
        return this.sn != null && this.sn.reward != null && !this.sn.reward.isEmpty();
    }

    // ========================================================================
    // SECTION 23: VISUALIZATION METHODS
    // Methods for graphical display of the network
    // ========================================================================

    /**
     * Displays an interactive visualization of this queueing network.
     * Uses default window title and dimensions.
     */
    public void plot() {
        plot("Network: " + this.getName());
    }

    /**
     * Displays an interactive visualization of this queueing network.
     *
     * @param title the window title
     */
    public void plot(String title) {
        plot(title, 800, 600);
    }

    /**
     * Displays an interactive visualization of this queueing network.
     *
     * @param title  the window title
     * @param width  the window width
     * @param height the window height
     */
    public void plot(String title, int width, int height) {
        ModelVisualizer visualizer = new ModelVisualizer(this);
        visualizer.buildGraph();
        visualizer.show(title, width, height);
    }

}


