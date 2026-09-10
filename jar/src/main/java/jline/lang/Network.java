/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import static jline.GlobalConstants.Inf;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.io.Ret;
// import static jline.io.SysUtils.lineViewerGetPath;
import static jline.io.SysUtils.jmtGetPath;
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
import jline.solvers.wrappers.jmt.SolverJMT;
import jline.solvers.mva.SolverMVA;
import jline.util.Maths;
import jline.util.Pair;
import jline.lang.reward.RewardFunction;
import jline.util.SerializableFunction;
import org.apache.commons.math3.complex.Complex;
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
import static jline.io.InputOutput.line_debug;
import static jline.api.mam.Map_pie.map_pie;
import static jline.api.mam.Map_erlang.map_erlang;
import static jline.api.mc.Dtmc_stochcomp.dtmc_stochcomp;
import static jline.api.sn.SnGetDemandsChain.snGetDemandsChain;
import static jline.api.sn.SnGetProductFormChainParams.snGetProductFormChainParams;
import static jline.api.sn.SnGetProductFormParams.snGetProductFormParams;
import static jline.api.sn.SnHasClassSwitching.snHasClassSwitching;
import static jline.api.sn.SnHasDPS.snHasDPS;
import static jline.api.sn.SnHasDPSPRIO.snHasDPSPRIO;
import static jline.api.sn.SnHasFCFS.snHasFCFS;
import static jline.api.sn.SnHasGPS.snHasGPS;
import static jline.api.sn.SnHasGPSPRIO.snHasGPSPRIO;
import static jline.api.sn.SnHasHOL.snHasHOL;
import static jline.api.sn.SnHasHomogeneousScheduling.snHasHomogeneousScheduling;
import static jline.api.sn.SnHasINF.snHasINF;
import static jline.api.sn.SnHasLCFS.snHasLCFS;
import static jline.api.sn.SnHasLCFSPR.snHasLCFSPR;
import static jline.api.sn.SnHasLEPT.snHasLEPT;
import static jline.api.sn.SnHasLJF.snHasLJF;
import static jline.api.sn.SnHasMultiChain.snHasMultiChain;
import static jline.api.sn.SnHasMultiClassFCFS.snHasMultiClassFCFS;
import static jline.api.sn.SnHasMultiClassHeterFCFS.snHasMultiClassHeterFCFS;
import static jline.api.sn.SnHasMultiClass.snHasMultiClass;
import static jline.api.sn.SnHasMultiServer.snHasMultiServer;
import static jline.api.sn.SnHasPS.snHasPS;
import static jline.api.sn.SnHasPSPRIO.snHasPSPRIO;
import static jline.api.sn.SnHasProductForm.snHasProductForm;
import static jline.api.sn.SnHasSEPT.snHasSEPT;
import static jline.api.sn.SnHasSIRO.snHasSIRO;
import static jline.api.sn.SnHasSJF.snHasSJF;
import static jline.api.sn.SnHasSingleChain.snHasSingleChain;
import static jline.api.sn.SnHasSingleClass.snHasSingleClass;
import static jline.api.sn.SnIsStateValid.snIsStateValid;
import static jline.api.sn.SnRefreshVisits.snRefreshVisits;
import static jline.api.sn.SnRtnodesToRtorig.snRtnodesToRtorig;
import static jline.io.InputOutput.*;
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
    public List<Node> nodes;
    private boolean hasStruct;
    public boolean isFJAugmented = false; // true on FJ tag-augmented copies (see ModelAdapter.fjtag)
    private NetworkStruct sn;
    // Network-level globally state-dependent (Whittle) rate scaling phi(n) over the
    // FULL (nstations x nclasses) population matrix; see setGlobalDependence
    private SerializableFunction<Matrix, Matrix> gdScaling;
    private Matrix gdScalingPeak;
    // Open-class truncation used to materialize gdScaling onto the JSON wire; solving ignores it.
    private int gdScalingCutoff = 10;
    // Markov reward definitions live on the model (not the transient struct) so they
    // survive resetStruct()/refreshStruct(); refreshStruct copies them into sn.reward.
    private Map<String, RewardFunction> rewardFunctions;
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
        super.setAttribute(this.attribute); // keep Model.attribute pointing at the typed container
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

        // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
        if (!(nodes[nodes.length - 1] instanceof Sink) && nodes[nodes.length - 1] != nodes[0]) {
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

        // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
        if (!(nodes.get(nodes.size() - 1) instanceof Sink) && nodes.get(nodes.size() - 1) != nodes.get(0)) {
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

        nodes.add(new Source(model, "Source"));

        for (int i = 0; i < M; i++) {
            if (strategy[i] == SchedStrategy.INF) {
                nodes.add(new Delay(model, "Station" + (i + 1)));
            } else {
                Queue queueNode = new Queue(model, "Station" + (i + 1), strategy[i]);
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
     * Creates an open cluster network: Source -> Dispatcher (Router) -> Servers -> Sink.
     *
     * <p>The dispatcher is a {@link Router} that distributes incoming jobs to the M parallel
     * server queues according to the supplied {@link RoutingStrategy} (RAND, RROBIN, JSQ, ...).
     *
     * @param lambda      arrival rate matrix [1 x R]; entry r is the per-class arrival rate
     * @param D           service time matrix [M x R]; entry (i, r) is the mean service time of class r at server i
     * @param strategy    per-server scheduling strategies (length M)
     * @param S           server count matrix [M x 1]; entry i is the multiplicity of server i (1 = single server)
     * @param dispatching dispatching policy applied at the router for every class
     * @return configured open cluster model
     */
    public static Network cluster(Matrix lambda, Matrix D, SchedStrategy[] strategy, Matrix S,
                                     RoutingStrategy dispatching) {
        Network model = new Network("Cluster");
        int M = D.getNumRows();
        int R = D.getNumCols();

        Source source = new Source(model, "Source");
        Router dispatcher = new Router(model, "Dispatcher");
        Queue[] servers = new Queue[M];
        for (int i = 0; i < M; i++) {
            servers[i] = new Queue(model, "Station" + (i + 1), strategy[i]);
            int c = (S != null && S.getNumRows() > i) ? (int) S.get(i, 0) : 1;
            if (c > 1) {
                servers[i].setNumberOfServers(c);
            }
        }
        Sink sink = new Sink(model, "Sink");

        List<OpenClass> jobclasses = new ArrayList<>(R);
        for (int r = 0; r < R; r++) {
            OpenClass cls = new OpenClass(model, "Class" + (r + 1), 0);
            jobclasses.add(cls);
            source.setArrival(cls, Exp.fitMean(1.0 / lambda.get(0, r)));
            for (int i = 0; i < M; i++) {
                servers[i].setService(cls, Exp.fitMean(D.get(i, r)));
            }
        }

        model.addLink(source, dispatcher);
        for (int i = 0; i < M; i++) {
            model.addLink(dispatcher, servers[i]);
            model.addLink(servers[i], sink);
        }
        for (int r = 0; r < R; r++) {
            dispatcher.setRouting(jobclasses.get(r), dispatching);
        }
        return model;
    }

    /**
     * Creates an open FCFS cluster with one server per queue.
     *
     * @param lambda      arrival rate matrix [1 x R]
     * @param D           service time matrix [M x R]
     * @param S           server count matrix [M x 1]
     * @param dispatching dispatching policy applied at the router
     * @return configured open cluster with FCFS servers
     */
    public static Network clusterFcfs(Matrix lambda, Matrix D, Matrix S, RoutingStrategy dispatching) {
        int M = D.getNumRows();
        SchedStrategy[] strategy = new SchedStrategy[M];
        for (int i = 0; i < M; i++) strategy[i] = SchedStrategy.FCFS;
        return Network.cluster(lambda, D, strategy, S, dispatching);
    }

    /**
     * Creates an open PS cluster.
     *
     * @param lambda      arrival rate matrix [1 x R]
     * @param D           service time matrix [M x R]
     * @param S           server count matrix [M x 1]
     * @param dispatching dispatching policy applied at the router
     * @return configured open cluster with PS servers
     */
    public static Network clusterPs(Matrix lambda, Matrix D, Matrix S, RoutingStrategy dispatching) {
        int M = D.getNumRows();
        SchedStrategy[] strategy = new SchedStrategy[M];
        for (int i = 0; i < M; i++) strategy[i] = SchedStrategy.PS;
        return Network.cluster(lambda, D, strategy, S, dispatching);
    }

    /**
     * Creates an open PS cluster with one server per queue.
     *
     * @param lambda      arrival rate matrix [1 x R]
     * @param D           service time matrix [M x R]
     * @param dispatching dispatching policy applied at the router
     * @return configured open cluster with single-server PS queues
     */
    public static Network clusterPs(Matrix lambda, Matrix D, RoutingStrategy dispatching) {
        int M = D.getNumRows();
        Matrix S = new Matrix(M, 1, M);
        S.fill(1.0);
        return Network.clusterPs(lambda, D, S, dispatching);
    }

    /**
     * Creates a closed cluster network: Think (Delay) -> Dispatcher (Router) -> Servers -> Think.
     *
     * <p>The think station is the reference station for every closed class; jobs cycle from the
     * delay to the dispatcher, are dispatched to one of the parallel servers, and return to the
     * delay on completion.
     *
     * @param N           per-class population matrix [1 x R]
     * @param Z           per-class think time matrix [1 x R]
     * @param D           service time matrix [M x R]
     * @param strategy    per-server scheduling strategies (length M)
     * @param S           server count matrix [M x 1]
     * @param dispatching dispatching policy applied at the router for every class
     * @return configured closed cluster model
     */
    public static Network clusterClosed(Matrix N, Matrix Z, Matrix D, SchedStrategy[] strategy,
                                           Matrix S, RoutingStrategy dispatching) {
        Network model = new Network("Cluster");
        int M = D.getNumRows();
        int R = D.getNumCols();

        Delay think = new Delay(model, "Think");
        Router dispatcher = new Router(model, "Dispatcher");
        Queue[] servers = new Queue[M];
        for (int i = 0; i < M; i++) {
            servers[i] = new Queue(model, "Station" + (i + 1), strategy[i]);
            int c = (S != null && S.getNumRows() > i) ? (int) S.get(i, 0) : 1;
            if (c > 1) {
                servers[i].setNumberOfServers(c);
            }
        }

        List<ClosedClass> jobclasses = new ArrayList<>(R);
        for (int r = 0; r < R; r++) {
            ClosedClass cls = new ClosedClass(model, "Class" + (r + 1), (int) N.get(0, r), think, 0);
            jobclasses.add(cls);
            think.setService(cls, Exp.fitMean(Z.get(0, r)));
            for (int i = 0; i < M; i++) {
                servers[i].setService(cls, Exp.fitMean(D.get(i, r)));
            }
        }

        model.addLink(think, dispatcher);
        for (int i = 0; i < M; i++) {
            model.addLink(dispatcher, servers[i]);
            model.addLink(servers[i], think);
        }
        for (int r = 0; r < R; r++) {
            dispatcher.setRouting(jobclasses.get(r), dispatching);
        }
        return model;
    }

    /**
     * Creates a mixed cluster network in which open and closed classes share the dispatcher
     * and the servers: open classes flow Source -&gt; Dispatcher -&gt; Servers -&gt; Sink while closed
     * classes cycle Think (Delay) -&gt; Dispatcher -&gt; Servers -&gt; Think.
     *
     * <p>Classes are ordered open first: columns 1..Ro of {@code D} refer to the open classes
     * and columns Ro+1..Ro+Rc to the closed ones.
     *
     * @param lambda      per-class arrival rate matrix [1 x Ro] of the open classes
     * @param N           per-class population matrix [1 x Rc] of the closed classes
     * @param Z           per-class think time matrix [1 x Rc] of the closed classes
     * @param D           service time matrix [M x (Ro+Rc)]; entry (i, r) is the mean service time
     *                    of class r at server i
     * @param strategy    per-server scheduling strategies (length M)
     * @param S           server count matrix [M x 1]; entry i is the multiplicity of server i
     * @param dispatching dispatching policy applied at the router for every class
     * @return configured mixed cluster model
     */
    public static Network clusterMixed(Matrix lambda, Matrix N, Matrix Z, Matrix D,
                                       SchedStrategy[] strategy, Matrix S,
                                       RoutingStrategy dispatching) {
        int M = D.getNumRows();
        int R = D.getNumCols();
        int Ro = lambda.length();
        int Rc = N.length();
        if (R != Ro + Rc) {
            throw new IllegalArgumentException("D must have lambda.length()+N.length() columns");
        }
        if (Z.length() != Rc) {
            throw new IllegalArgumentException("N and Z must have the same length");
        }

        Network model = new Network("Cluster");

        Source source = new Source(model, "Source");
        Delay think = new Delay(model, "Think");
        Router dispatcher = new Router(model, "Dispatcher");
        Queue[] servers = new Queue[M];
        for (int i = 0; i < M; i++) {
            servers[i] = new Queue(model, "Station" + (i + 1), strategy[i]);
            int c = (S != null && S.getNumRows() > i) ? (int) S.get(i, 0) : 1;
            if (c > 1) {
                servers[i].setNumberOfServers(c);
            }
        }
        Sink sink = new Sink(model, "Sink");

        List<JobClass> jobclasses = new ArrayList<>(R);
        for (int r = 0; r < Ro; r++) {
            OpenClass cls = new OpenClass(model, "Class" + (r + 1), 0);
            jobclasses.add(cls);
            source.setArrival(cls, Exp.fitMean(1.0 / lambda.get(r)));
            think.setService(cls, Disabled.getInstance());
            for (int i = 0; i < M; i++) {
                servers[i].setService(cls, Exp.fitMean(D.get(i, r)));
            }
        }
        for (int c = 0; c < Rc; c++) {
            int r = Ro + c;
            ClosedClass cls = new ClosedClass(model, "Class" + (r + 1), (int) N.get(c), think, 0);
            jobclasses.add(cls);
            think.setService(cls, Exp.fitMean(Z.get(c)));
            for (int i = 0; i < M; i++) {
                servers[i].setService(cls, Exp.fitMean(D.get(i, r)));
            }
        }

        model.addLink(source, dispatcher);
        model.addLink(think, dispatcher);
        for (int i = 0; i < M; i++) {
            model.addLink(dispatcher, servers[i]);
            model.addLink(servers[i], sink);
            model.addLink(servers[i], think);
        }

        // Class-specific exits: the servers feed the sink for open classes and the delay for
        // closed ones, so the shared arcs carry explicit per-class probabilities.
        for (int r = 0; r < R; r++) {
            JobClass cls = jobclasses.get(r);
            dispatcher.setRouting(cls, dispatching);
            double isOpen = r < Ro ? 1.0 : 0.0;
            for (int i = 0; i < M; i++) {
                servers[i].setProbRouting(cls, sink, isOpen);
                servers[i].setProbRouting(cls, think, 1.0 - isOpen);
            }
            source.setProbRouting(cls, dispatcher, isOpen);
            think.setProbRouting(cls, dispatcher, 1.0 - isOpen);
        }
        return model;
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
        region.setModel(this);
        // a cached struct predates this region, so it would carry nregions=0
        resetStruct();
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

    /**
     * Mean tardiness handles, Tard(i,r) for class r at station i.
     *
     * <p>The handles were already built by {@link #getAvgHandles()}; only the
     * named accessor was missing, so a caller had to reach through
     * SolverAvgHandles by field. Twin of the MATLAB
     * {@code MNetwork.getAvgTardHandles}.
     */
    public AvgHandle getAvgTardHandles() {
        return this.getAvgHandles().getAvgTardHandles();
    }

    /** Mean system tardiness handles, SysTard(1,r) for class r. */
    public AvgHandle getAvgSysTardHandles() {
        return this.getAvgHandles().getAvgSysTardHandles();
    }

    /** The reducibility structure plus one suggested repair per absorbing station. */
    public RoutingErgodicity.ReducibilityInfo getReducibilityInfo() {
        return RoutingErgodicity.getReducibilityInfo(this);
    }

    /** The stations that are absorbing: once a job enters, it never leaves. */
    public List<Station> getAbsorbingStations() {
        return RoutingErgodicity.getAbsorbingStations(this);
    }

    /**
     * A routing matrix that makes the network ergodic. Does NOT relink; apply it
     * with {@code model.link(P)}.
     */
    public RoutingMatrix makeErgodic() {
        return RoutingErgodicity.makeErgodic(this, null);
    }

    /** As {@link #makeErgodic()}, routing absorbing stations to TARGETNAME. */
    public RoutingMatrix makeErgodic(String targetName) {
        return RoutingErgodicity.makeErgodic(this, targetName);
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

    /**
     * The chains of this network, each carrying the job classes it contains.
     *
     * <p>{@code sn.inchain.get(c)} IS A LIST OF CLASS INDICES, not an indicator
     * vector over classes -- which is how the rest of this file reads it (see
     * the chain-capacity loops around line 5100). Testing {@code get(0,r) == 1}
     * instead returned a chain with NO classes on every single-class model,
     * because inchain[0] is then the vector [0] and 0 != 1. Everything built on
     * top silently degraded: SolverCTMC.getCdfRespT looped over zero classes and
     * its catch returned a zero matrix.
     *
     * <p>The chain list is also rebuilt rather than appended to; the old code
     * accumulated a fresh set of chains on every call.
     */
    public List<Chain> getChains() {
        NetworkStruct sn = this.getStruct();
        chains.clear();
        for (int c = 0; c < sn.chains.getNumElements(); c++) {
            List<JobClass> chainClasses = new ArrayList<>();
            Matrix inchain_c = sn.inchain.get(c);
            for (int idx = 0; idx < inchain_c.length(); idx++) {
                int r = (int) inchain_c.get(idx);
                if (r >= 0 && r < sn.nclasses) {
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

    /**
     * The chain containing the given class.
     *
     * <p>The old form ignored {@code jobClass} entirely and returned the first
     * chain whose indicator test passed, so on a multi-chain model it answered
     * with the wrong chain. It reads inchain as a list of class indices, as
     * {@link #getChains()} does.
     */
    public Chain getClassChain(JobClass jobClass) {
        int c = getClassChainIndex(jobClass);
        if (c < 0) return null;
        List<Chain> all = getChains();
        return c < all.size() ? all.get(c) : null;
    }

    public int getClassChainIndex(JobClass jobClass) {
        NetworkStruct sn = this.getStruct();
        int target = jobClass.getIndex() - 1;
        for (int c = 0; c < sn.chains.getNumElements(); c++) {
            Matrix inchain_c = sn.inchain.get(c);
            for (int idx = 0; idx < inchain_c.length(); idx++) {
                if ((int) inchain_c.get(idx) == target) return c;
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


    /**
     * Whether the cached class-switching mask still describes the current class
     * set. link() records csMatrix for the class count of the model it linked,
     * but a fork-join transformation copies the model and then adds auxiliary
     * classes, so the copy can carry a mask of a different width. Reading a
     * stale mask indexes classes that no longer exist, which surfaced as an
     * out-of-bounds write while grouping chains. When the width disagrees the
     * mask is ignored and the chains are recovered from the routing instead,
     * as MATLAB does when no mask is available.
     */
    private boolean hasUsableCsMatrix(int K) {
        return this.csMatrix != null
                && this.csMatrix.getNumRows() == K
                && this.csMatrix.getNumCols() == K;
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

    /**
     * Gets the class-dependence functions beta_i(n) of the stations that declare
     * one. Stations without class dependence are deliberately ABSENT from the map
     * rather than mapped to a constant 1: the entry is a class-dependent service
     * RATE (Sauer 1983, eq. (40)), for which a constant 1 would assert that every
     * class completes at rate 1, which is not load independence (a load-independent
     * station is beta_{i,r}(n) = mu_i * n_r/|n|, whose n_r/|n| factor is what
     * regenerates the multinomial). Consumers treat a missing entry as "no class
     * dependence": Pfqn_cdfun skips it and Pfqn_conv keeps the station on the
     * load-independent recurrence. Consumers that index every station
     * unconditionally fill the gaps with the neutral scaling themselves.
     *
     * @return map from station to its class-dependence function; empty if none
     */
    public Map<Station, SerializableFunction<Matrix, Matrix>> getLimitedClassDependence() {
        Map<Station, SerializableFunction<Matrix, Matrix>> gamma = new HashMap<Station, SerializableFunction<Matrix, Matrix>>();

        for (Station station : this.stations) {
            if (station.getLimitedClassDependence() != null) gamma.put(station, station.getLimitedClassDependence());
        }

        return gamma;
    }

    /**
     * Peak (max) class-dependent rate scaling per class for each class-dependent
     * station, as a 1xR row vector (scalar declarations broadcast to R classes).
     * Used to normalize utilization as U = T*S/peak. A station that declares a
     * class dependence without a peak is a MODEL DEFECT and is refused here, as
     * in MATLAB's getLimitedClassDependencePeak.m: the peak is not recoverable
     * from the handle, and guessing it reports a number that is not a
     * utilization.
     *
     * @return map from station to its 1xR peak vector; empty if no class dependence
     * @throws RuntimeException if a class-dependent station declared no peak
     */
    public Map<Station, Matrix> getLimitedClassDependencePeak() {
        Map<Station, Matrix> peakMap = new HashMap<Station, Matrix>();
        int K = getNumberOfClasses();
        for (Station station : this.stations) {
            SerializableFunction<Matrix, Matrix> beta = station.getLimitedClassDependence();
            if (beta == null) continue;
            Matrix declared = station.getLimitedClassDependencePeak();
            Matrix peakVec = new Matrix(1, K);
            if (declared != null && !declared.isEmpty()) {
                if (declared.length() == 1) {
                    for (int r = 0; r < K; r++) peakVec.set(0, r, declared.get(0));
                } else {
                    for (int r = 0; r < K; r++) peakVec.set(0, r, declared.get(r));
                }
            } else {
                // NOT DERIVED. Sweeping beta for max_n beta(n) needs a bound the
                // handle does not carry, and an open class has no bound at all,
                // so a "derived" peak there is just beta(0) -- a number that
                // reads as a utilization and is not one. MATLAB's
                // getLimitedClassDependencePeak.m errors on exactly this.
                throw new RuntimeException("Class-dependent station '" + station.getName()
                        + "' has no declared peak rate; use setClassDependence(beta, peakRatePerClass).");
            }
            peakMap.put(station, peakVec);
        }
        return peakMap;
    }

    /**
     * Map from station to its joint-dependence (non-product-form) function
     * eta_i(n); empty if none. Twin of {@link #getLimitedClassDependence()}.
     *
     * @return map from station to its joint-dependence function; empty if none
     */
    public Map<Station, SerializableFunction<Matrix, Matrix>> getLimitedJointDependence() {
        Map<Station, SerializableFunction<Matrix, Matrix>> eta = new HashMap<Station, SerializableFunction<Matrix, Matrix>>();

        for (Station station : this.stations) {
            if (station.getLimitedJointDependence() != null) eta.put(station, station.getLimitedJointDependence());
        }

        return eta;
    }

    /**
     * Peak joint-dependent rate scaling per class for each joint-dependent
     * station, as a 1xR row vector. Twin of
     * {@link #getLimitedClassDependencePeak()}.
     *
     * @return map from station to its 1xR peak vector; empty if no joint dependence
     * @throws RuntimeException if a joint-dependent station declared no peak
     */
    public Map<Station, Matrix> getLimitedJointDependencePeak() {
        Map<Station, Matrix> peakMap = new HashMap<Station, Matrix>();
        int K = getNumberOfClasses();
        for (Station station : this.stations) {
            SerializableFunction<Matrix, Matrix> eta = station.getLimitedJointDependence();
            if (eta == null) continue;
            Matrix declared = station.getLimitedJointDependencePeak();
            Matrix peakVec = new Matrix(1, K);
            if (declared != null && !declared.isEmpty()) {
                if (declared.length() == 1) {
                    for (int r = 0; r < K; r++) peakVec.set(0, r, declared.get(0));
                } else {
                    for (int r = 0; r < K; r++) peakVec.set(0, r, declared.get(r));
                }
            } else {
                throw new RuntimeException("Joint-dependent station '" + station.getName()
                        + "' has no declared peak rate; use setJointDependence(eta, peakRatePerClass).");
            }
            peakMap.put(station, peakVec);
        }
        return peakMap;
    }

    /**
     * Declares a globally state-dependent service-rate scaling phi(n), where n is
     * the FULL (nstations x nclasses) population matrix rather than the population
     * local to one station. This is the Whittle-network primitive: when phi satisfies
     * phi_s(n) phi_t(n-e_s) = phi_t(n) phi_s(n-e_t) the chain is reversible, has the
     * product form pi(n) ~ Phi(n) prod rho_s^n_s and is insensitive. It also expresses
     * bandwidth sharing, where one route holds several links at once and no per-station
     * scaling can reproduce the coupling.
     *
     * <p>phi returns a 1x1 scalar (broadcast), an (M x 1) column (per station) or an
     * (M x K) matrix. The effective rate of class r at station i is its base rate times
     * phi(i,r), composing multiplicatively with any load-, class- or joint-dependence.
     * Only SolverCTMC declares support for it.
     *
     * @param phi  the scaling handle over the full population matrix
     * @param peak REQUIRED peak scaling, 1x1, (M x 1) or (M x K), normalizing Util=T*S/peak
     */
    public void setGlobalDependence(SerializableFunction<Matrix, Matrix> phi, Matrix peak) {
        setGlobalDependence(phi, peak, 10);
    }

    /**
     * As {@link #setGlobalDependence(SerializableFunction, Matrix)}, with an explicit
     * per-slot OPEN-class truncation used when phi is materialized onto the JSON wire
     * (closed classes are tabulated up to their own population). It plays no part in
     * solving, and exists because a handle cannot cross a language boundary: the writer
     * needs to know how far the lattice extends. Set it to the cutoff the model is
     * solved at.
     *
     * @param phi        the scaling handle over the full population matrix
     * @param peak       REQUIRED peak scaling, 1x1, (M x 1) or (M x K)
     * @param wireCutoff positive open-class truncation for JSON serialization
     */
    public void setGlobalDependence(SerializableFunction<Matrix, Matrix> phi, Matrix peak, int wireCutoff) {
        if (phi == null) {
            throw new IllegalArgumentException("Global dependence must be specified through a function.");
        }
        if (peak == null || peak.isEmpty()) {
            throw new IllegalArgumentException("Global dependence requires an explicit peak rate: setGlobalDependence(phi, peak).");
        }
        int M = getNumberOfStations();
        int K = getNumberOfClasses();
        for (int i = 0; i < peak.length(); i++) {
            if (peak.get(i) <= 0) {
                throw new IllegalArgumentException("peak must be positive.");
            }
        }
        // Probe now so a wrong output shape is refused at declaration time rather
        // than midway through state-space generation.
        Matrix[] probes = new Matrix[]{new Matrix(M, K), Matrix.ones(M, K)};
        for (int p = 0; p < probes.length; p++) {
            Matrix v = phi.apply(probes[p]);
            if (v == null || v.isEmpty()) {
                throw new IllegalArgumentException("The global dependence handle returned no scaling.");
            }
            boolean okShape = (v.length() == 1)
                    || (v.getNumRows() == M && v.getNumCols() == 1)
                    || (v.getNumRows() == M && v.getNumCols() == K);
            if (!okShape) {
                throw new IllegalArgumentException("The global dependence handle must return a scalar, an (" + M + " x 1) column or an (" + M + " x " + K + ") matrix.");
            }
            for (int j = 0; j < v.length(); j++) {
                if (!Double.isFinite(v.get(j)) || v.get(j) < 0) {
                    throw new IllegalArgumentException("The global dependence handle must return finite nonnegative scalings.");
                }
            }
        }
        if (wireCutoff < 1) {
            throw new IllegalArgumentException("wireCutoff must be a positive integer.");
        }
        this.gdScaling = phi;
        this.gdScalingPeak = expandGlobalPeak(peak, M, K);
        this.gdScalingCutoff = wireCutoff;
        this.hasStruct = false;
        this.sn = null;
    }

    /** Network-level global dependence handle, or null if the model declares none. */
    public SerializableFunction<Matrix, Matrix> getGlobalDependence() {
        return this.gdScaling;
    }

    /** Per-slot open-class wire truncation of the global dependence (default 10). */
    public int getGlobalDependenceCutoff() {
        return this.gdScalingCutoff;
    }

    /** (nstations x nclasses) peak of the global dependence, or null if none is declared. */
    public Matrix getGlobalDependencePeak() {
        if (this.gdScaling == null) {
            return null;
        }
        return expandGlobalPeak(this.gdScalingPeak, getNumberOfStations(), getNumberOfClasses());
    }

    private static Matrix expandGlobalPeak(Matrix peak, int M, int K) {
        Matrix out = new Matrix(M, K);
        if (peak.length() == 1) {
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) out.set(i, r, peak.get(0));
            }
        } else if (peak.getNumRows() == M && peak.getNumCols() == 1) {
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) out.set(i, r, peak.get(i, 0));
            }
        } else if (peak.getNumRows() == M && peak.getNumCols() == K) {
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) out.set(i, r, peak.get(i, r));
            }
        } else {
            throw new IllegalArgumentException("peak must be a scalar, an (" + M + " x 1) column or an (" + M + " x " + K + ") matrix.");
        }
        return out;
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

    /**
     * Fills the variable-forking-level matrices of a fork's node parameters.
     *
     * <p>Twin of the {@code case 'Fork'} block of MATLAB {@code refreshLocalVars}
     * and of {@code Network._refresh_nodeparam} in native Python. The matrices
     * are (nnodes x nclasses) and indexed by DESTINATION NODE, so a relink
     * cannot permute them; connectivity is read off the fork's own output
     * strategies, which is where the destinations live before {@code connmatrix}
     * is rebuilt.</p>
     *
     * <p>When any override is present the scalar {@code fanOut} is reset to the
     * mean over the connected links, so a solver that has not been taught the
     * matrices degrades to E[tasks per link] rather than to a number the fork
     * never emits.</p>
     */
    private void buildForkFanout(ForkNodeParam forkParam, Forker forker, Node node) {
        int I = this.getNumberOfNodes();
        int K = this.getNumberOfClasses();
        Matrix fanOutLink = new Matrix(I, K);
        Matrix fanOutProb = new Matrix(I, K);
        DiscreteSampler[][] fanOutDist = new DiscreteSampler[I][K];
        boolean[][] connrc = new boolean[I][K];

        for (OutputStrategy os : forker.getOutputStrategies()) {
            if (os.getDestination() == null || os.getJobClass() == null) continue;
            int k = this.getNodeIndex(os.getDestination());
            int r = os.getJobClass().getIndex() - 1;  // getIndex() is 1-based
            if (k < 0 || k >= I || r < 0 || r >= K) continue;
            connrc[k][r] = true;
            fanOutLink.set(k, r, forker.tasksPerLink);
            fanOutProb.set(k, r, 1.0);
        }

        for (Forker.ForkOverride ov : forker.tasksPerLinkByDest)
            for (int k : forkOverrideDests(ov.dest, connrc, I))
                fanOutLink.set(k, ov.jobClass - 1, ov.value);
        for (Forker.ForkOverride ov : forker.branchProb)
            for (int k : forkOverrideDests(ov.dest, connrc, I))
                fanOutProb.set(k, ov.jobClass - 1, ov.value);
        for (Forker.ForkOverride ov : forker.tasksPerLinkDist)
            for (int k : forkOverrideDests(ov.dest, connrc, I)) {
                fanOutDist[k][ov.jobClass - 1] = ov.dist;
                // the scalar slot carries the mean, so a consumer that only
                // reads fanOutLink still sees E[tasks per link]
                fanOutLink.set(k, ov.jobClass - 1, ov.dist.getMean());
            }

        forkParam.fanOutLink = fanOutLink;
        forkParam.fanOutProb = fanOutProb;
        forkParam.fanOutDist = fanOutDist;

        if (!forker.tasksPerLinkDist.isEmpty() || !forker.tasksPerLinkByDest.isEmpty()
                || !forker.branchProb.isEmpty()) {
            // The branch probability is folded in HERE and not into fanOutLink,
            // because JMT and LDES read the two separately: fanOutLink is the
            // count GIVEN the branch fires, fanOutProb is whether it fires.
            double acc = 0.0;
            int cnt = 0;
            for (int k = 0; k < I; k++)
                for (int r = 0; r < K; r++)
                    if (connrc[k][r]) { acc += fanOutLink.get(k, r) * fanOutProb.get(k, r); cnt++; }
            if (cnt > 0) forkParam.fanOut = acc / cnt;
        }
    }

    /** Node indexes a fork override applies to; an empty name means every connected link. */
    private List<Integer> forkOverrideDests(String dest, boolean[][] connrc, int I) {
        List<Integer> out = new ArrayList<Integer>();
        if (dest == null || dest.isEmpty()) {
            for (int k = 0; k < I; k++) {
                boolean any = false;
                for (int r = 0; r < connrc[k].length; r++) any = any || connrc[k][r];
                if (any) out.add(k);
            }
            return out;
        }
        int k = this.getNodeIndex(dest);
        if (k < 0) {
            line_error(mfilename(new Object() {}),
                    "Fork override names destination \"" + dest + "\", which is not a node of this model.");
        }
        out.add(k);
        return out;
    }

    /**
     * True when {@code station} is a queue registered as a retrieval-system queue
     * of some Cache (via {@link jline.lang.nodes.Cache#setRetrievalSystem}). Such a
     * queue carries a finite per-class capacity and a WaitingQueue drop marker as
     * an internal default of the delayed-hit retrieval mechanism, not as a user
     * request, so it is excluded from the open-class WAITQ+finite-capacity guard in
     * the same way Place and Join are.
     */
    private boolean isRetrievalSystemQueue(Station station) {
        if (station == null) {
            return false;
        }
        int nodeIdx = this.getNodeIndex(station);
        if (nodeIdx < 0) {
            return false;
        }
        for (Node node : this.nodes) {
            if (node instanceof Cache) {
                Map<Integer, List<Integer>> rsq = ((Cache) node).getRetrievalSystemQueueIndices();
                if (rsq != null) {
                    for (List<Integer> idxs : rsq.values()) {
                        if (idxs != null && idxs.contains(nodeIdx)) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
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
                else if (nodeIter instanceof StatefulFork) nodetypes.add(NodeType.Fork);
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
        } else if (distr instanceof DMAP) {
            return ProcessType.DMAP;
        } else if (distr instanceof MAP) {
            return ProcessType.MAP;
        } else if (distr instanceof ME) {
            return ProcessType.ME;
        } else if (distr instanceof RAP) {
            return ProcessType.RAP;
        } else if (distr instanceof Immediate) {
            return ProcessType.IMMEDIATE;
        } else if (distr instanceof NHPP) {
            return ProcessType.NHPP;
        } else if (distr instanceof MAPt) {
            return ProcessType.MAPT;
        } else if (distr instanceof PHt) {
            return ProcessType.PHT;
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
                // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
                for (int jnd = 0; jnd < I; jnd++) {
                    for (int k = 0; k < K; k++) {
                        if (conn.get(ind, jnd) > 0) {
                            JobClass jobclass = this.jobClasses.get(k);
                            List<OutputStrategy> outputStrategy_k = node.getOutput().getOutputStrategyByClass(jobclass);
                            rtnodes.set(ind * K + k, jnd * K + k, 1.0);
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
                            // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
                            if (!(jobclass instanceof SelfLoopingClass)
                                    && !(jobclass instanceof Signal)
                                    && !(jobclass instanceof OpenSignal)
                                    && !(jobclass instanceof ClosedSignal)) {
                                sum = conn.sumRows(ind);
                                if (sum > 0) {
                                    for (int jnd = 0; jnd < I; jnd++) {
                                        if (conn.get(ind, jnd) > 0) {
                                            rtnodes.set(ind * K + k, jnd * K + k, 1.0 / sum);
                                        }
                                    }
                                }
                            }
                            break;
                        case WRROBIN:
                            // Use actual WRROBIN weights from output strategy
                            if (!outputStrategy_k.isEmpty()) {
                                double totalWeight = 0;
                                double[] destWeights = new double[I];
                                for (OutputStrategy ops : outputStrategy_k) {
                                    int j = this.getNodeIndex(ops.getDestination());
                                    if (j >= 0 && j < I) {
                                        double weight = ops.getProbability();  // Weight is stored in probability field
                                        destWeights[j] = weight;
                                        totalWeight += weight;
                                    }
                                }
                                if (totalWeight > 0) {
                                    for (int j = 0; j < I; j++) {
                                        if (destWeights[j] > 0) {
                                            rtnodes.set(ind * K + k, j * K + k, destWeights[j] / totalWeight);
                                        }
                                    }
                                }
                            } else {
                                // Fallback to uniform if no weights defined
                                if (isInf((NK.get(0, k)))) {
                                    sum = conn.sumRows(ind);
                                    for (int j = 0; j < I; j++) {
                                        if (conn.get(ind, j) > 0) rtnodes.set(ind * K + k, j * K + k, 1.0 / sum);
                                    }
                                } else if (!isSource_i && !isSink_i) {
                                    Matrix connectionClosed = conn.copy();
                                    if (idxSink >= 0) {
                                        if (connectionClosed.get(ind, idxSink) > 0) connectionClosed.remove(ind, idxSink);
                                    }
                                    sum = connectionClosed.sumRows(ind);
                                    for (int j = 0; j < I; j++) {
                                        if (connectionClosed.get(ind, j) > 0)
                                            rtnodes.set(ind * K + k, j * K + k, 1.0 / sum);
                                    }
                                }
                            }
                            break;
                        case RAND:
                        case RROBIN:
                        case JSQ:
                        case SQ:
                        case SDR:
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
                        boolean flagHit = false, flagMiss = false, flagReceived = false;
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
                        for (int item = 0; item < cacheClassSwitcher.retrievalClasses.getNumRows(); item++) {
                            for (int j = 0; j < cacheClassSwitcher.retrievalClasses.getNumCols(); j++) {
                                if (cacheClassSwitcher.retrievalClasses.get(item, j) == r) {
                                    flagReceived = true;
                                    break;
                                }
                            }
                            if (flagReceived) break;
                        }
                        if (!flagHit && !flagMiss && !flagReceived) {
                            double rowsum = Matrix.extractRows(Pcs, r, r + 1, null).elementSum();
                            for (int j = 0; j < Pcs.getNumCols(); j++) {
                                Pcs.set(r, j, Pcs.get(r, j) / rowsum);
                            }
                        }
                    }
                    for (int r = 0; r < K; r++) {
                        boolean flagHit = false, flagMiss = false, flagReceived = false;
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
                        for (int item = 0; item < cacheClassSwitcher.retrievalClasses.getNumRows(); item++) {
                            for (int j = 0; j < cacheClassSwitcher.retrievalClasses.getNumCols(); j++) {
                                if (cacheClassSwitcher.retrievalClasses.get(item, j) == r) {
                                    flagReceived = true;
                                    break;
                                }
                            }
                            if (flagReceived) break;
                        }
                        if (!flagHit && !flagMiss && !flagReceived) {
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
                            // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
                            if (cacheClassSwitcher.actualDelayedHitProb != null
                                    && cacheClassSwitcher.actualDelayedHitProb.length() > r
                                    && !Double.isNaN(cacheClassSwitcher.actualDelayedHitProb.get(r))) {
                                ph = ph + cacheClassSwitcher.actualDelayedHitProb.get(r);
                            }

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
                            boolean flagHit = false, flagMiss = false, flagReceived = false;
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
                            for (int item = 0; item < cacheClassSwitcher.retrievalClasses.getNumRows(); item++) {
                                for (int j = 0; j < cacheClassSwitcher.retrievalClasses.getNumCols(); j++) {
                                    if (cacheClassSwitcher.retrievalClasses.get(item, j) == r) {
                                        flagReceived = true;
                                        break;
                                    }
                                }
                                if (flagReceived) break;
                            }
                            if (flagHit || flagMiss || flagReceived) {
                                // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
                                for (int s = 0; s < K; s++) {
                                    double routingProb = Pij.get(r, s);
                                    if (routingProb > 0) {
                                        rtnodes.set(i * K + r, jnd * K + s, routingProb);
                                    }
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

        // ignore chains with a Pnodes column summing to 0: these classes cannot arrive to the node unless the column belongs to the source
        Matrix sumRtnodesCols = rtnodes.sumCols();
        Set<Integer> colsToIgnore = new HashSet<Integer>();
        for (int col = 0; col < sumRtnodesCols.getNumCols(); col++) {
            if (sumRtnodesCols.get(col) == 0) colsToIgnore.add(col);
        }
        if (hasOpen) {
            for (int i = idxSource * K; i < (idxSource + 1) * K; i++)
                colsToIgnore.remove(i);
        }

        // Sink-to-source reroute: see _kb/04-networkstruct.md ("rt is
        // pseudo-closed") for the rationale.
        if (!hasUsableCsMatrix(K)) {
            // symmetrized to group transient+recurrent states into chains: see _kb/04-networkstruct.md
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
                        // see _kb/04-networkstruct.md (refreshStruct.m field-population notes) for rationale
                        splitChains.add(rows.get(j));
                    }
                }
            }
        }

        if (splitChains.size() > 0) {
            Matrix newChains = new Matrix(chains.getNumRows() - splitChains.size(), chains.getNumCols());
            // see _kb/04-networkstruct.md (refreshStruct.m field-population notes) for rationale
            int newRow = 0;
            for (int i = 0; i < chains.getNumRows(); i++) {
                if (!splitChains.contains(i)) {
                    for (int j = 0; j < chains.getNumCols(); j++)
                        newChains.set(newRow, j, chains.get(i, j));
                    newRow++;
                }
            }
            chains = newChains;
        }

        // Stochastic-complement construction of rt: see _kb/04-networkstruct.md.
        List<Integer> statefulNodes = this.getIndexStatefulNodes();
        List<Integer> statefulNodeClasses = new ArrayList<Integer>(); //Not using JLineMatrix for performance consideration
        for (int i = 0; i < statefulNodes.size(); i++) {
            for (int j = 0; j < K; j++) {
                statefulNodeClasses.add(statefulNodes.get(i) * K + j);
            }
        }

        // sink-to-source reroute for open classes: see _kb/04-networkstruct.md ("Which class the Sink reroutes into")
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
            if (node instanceof Sink) {
                return RoutingStrategy.DISABLED;
            } else if (node instanceof Source) {
                // Source nodes should route open classes even without explicit arrivals
                return (c instanceof OpenClass) ? RoutingStrategy.RAND : RoutingStrategy.DISABLED;
            } else if (node instanceof jline.lang.nodes.Router) {
                // Router nodes should only route classes with explicit routing defined
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
                    SchedStrategy ss = station.getSchedStrategy();
                    // see _kb/04-networkstruct.md (State-space construction conventions) for rationale
                    if (ss == SchedStrategy.OI) {
                        ss = SchedStrategy.PAS;
                    }
                    res.put(station, ss);
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
     * Which solvers and solver methods can analyze THIS model.
     *
     * <pre>
     * model.findSolver()                  every (solver, method) pair that runs
     * model.findSolver("cdf", false)      ... that returns a passage-time law
     * model.findSolver("getCdfRespT", false)  the same question, by accessor
     * model.findSolver("", true)          also the pairs that are refused, and why
     * </pre>
     *
     * <p>One row per pair; see {@link jline.solvers.auto.SolverCandidate} for the
     * columns and {@code SolverCandidate.toTable} to print them. The method
     * column is the method name to pass as a solver method, so a row can be acted on
     * directly.
     *
     * <p>{@link #findMethod()} and {@link #help()} are aliases.
     *
     * @return one row per runnable (family, method) pair
     */
    public java.util.List<jline.solvers.auto.SolverCandidate> findSolver() {
        return findSolver("", false);
    }

    /**
     * Which solvers and solver methods can analyze this model, narrowed to one
     * measure and optionally including the refused pairs.
     *
     * @param metric  a measure group ("cdf") or the accessor that returns it
     *                ("getCdfRespT"); "" or "any" keeps every pair
     * @param showAll keep the refused pairs too, with the reason each was refused
     * @return the matching rows
     */
    public java.util.List<jline.solvers.auto.SolverCandidate> findSolver(String metric,
                                                                        boolean showAll) {
        // The gate lives in SolverAUTO, which is the class that already knows
        // every family, how to build one and what each refuses. Asking it here
        // rather than reimplementing the walk is what keeps the model's answer
        // and AUTO's own dispatch from being two opinions.
        //
        // The construction is silenced as well as the walk: SolverAUTO probes
        // every candidate with supports(model), which warns on a model one of
        // them refuses, and a report must not print.
        VerboseLevel saved = GlobalConstants.getVerbose();
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
        jline.solvers.auto.SolverAUTO auto;
        try {
            auto = new jline.solvers.auto.SolverAUTO(this);
        } finally {
            GlobalConstants.setVerbose(saved);
        }
        return auto.findSolver(metric, showAll);
    }

    /**
     * Alias of {@link #findSolver()}: which solvers and solver methods can
     * analyze this model.
     *
     * <p>The two names exist because the question is asked both ways round --
     * "which solver do I use" and "which method do I pass" -- and the answer is
     * the same table, whose method column carries the method name either caller needs.
     *
     * @return one row per runnable (family, method) pair
     */
    public java.util.List<jline.solvers.auto.SolverCandidate> findMethod() {
        return findSolver("", false);
    }

    /**
     * Alias of {@link #findSolver(String, boolean)}.
     *
     * @param metric  a measure group or the accessor that returns it
     * @param showAll keep the refused pairs too
     * @return the matching rows
     */
    public java.util.List<jline.solvers.auto.SolverCandidate> findMethod(String metric,
                                                                        boolean showAll) {
        return findSolver(metric, showAll);
    }

    /**
     * Alias of {@link #findSolver()}: what can this model be solved with?
     *
     * @return one row per runnable (family, method) pair
     */
    public java.util.List<jline.solvers.auto.SolverCandidate> help() {
        return findSolver("", false);
    }

    /**
     * Alias of {@link #findSolver(String, boolean)}.
     *
     * @param metric  a measure group or the accessor that returns it
     * @param showAll keep the refused pairs too
     * @return the matching rows
     */
    public java.util.List<jline.solvers.auto.SolverCandidate> help(String metric,
                                                                  boolean showAll) {
        return findSolver(metric, showAll);
    }

    /**
     * The answer of {@link #findBindingCapacity()}: whether a buffer binds and,
     * when it does, which one. {@code classIndex} is the 0-based class index of
     * a per-class buffer and -1 for a station-level one; {@code isOpen} says
     * whether an open class reaches the buffer, which decides the fallback
     * advice the refusal gives.
     */
    public static class BindingCapacity {
        public final boolean binds;
        public final Station station;
        public final double cap;
        public final int classIndex;
        public final JobClass jobClass;
        public final boolean isOpen;

        BindingCapacity() {
            this(null, Double.POSITIVE_INFINITY, -1, null, false);
        }

        BindingCapacity(Station station, double cap, int classIndex, JobClass jobClass, boolean isOpen) {
            this.binds = station != null;
            this.station = station;
            this.cap = cap;
            this.classIndex = classIndex;
            this.jobClass = jobClass;
            this.isOpen = isOpen;
        }
    }

    /**
     * The first station whose finite capacity can actually BIND, or a result
     * whose {@code binds} is false when no buffer in the model can refuse a job.
     * <p>
     * ONE PREDICATE, TWO CALLERS, the port of MATLAB
     * {@code MNetwork.findBindingCapacity}. {@code
     * NetworkSolver.bindingCapacityReason} turns the answer into the refusal the
     * product-form solvers raise, and {@link #getUsedLangFeatures()} marks the
     * registry name {@code FiniteCapacity} on the same answer, so a solver that
     * does not declare the name refuses exactly the models the structural gate
     * refuses.
     * <p>
     * The test reads the node-level cap / classCap the user set and the class
     * populations from the CLASS OBJECTS ({@link #getNumberOfJobs()}), never
     * sn.cap / sn.classcap: refreshCapacity derives a FINITE classcap (the chain
     * population) for every closed model, so an sn-level test would call every
     * closed model capped, and reading the struct from the recorder would
     * trigger a refresh on every feature query.
     * <p>
     * Only a capacity that can bind counts. A closed model whose station
     * capacity is at least the total population can never block a job, so the
     * declaration is a no-op (setCapacity(N) on a station of an N-job closed
     * model is a common idiom). The population of an open class is Inf, so any
     * finite capacity an open class can reach binds. A Cache model is exempt:
     * Cache builds retrieval queues that legitimately carry a per-class capacity
     * of 1, and the cache analyzers solve those rather than treating them as a
     * buffer constraint.
     *
     * @return - the binding buffer, or a result whose {@code binds} is false
     */
    public BindingCapacity findBindingCapacity() {
        for (Node node : this.nodes) {
            if (node instanceof Cache) {
                return new BindingCapacity();
            }
        }
        Matrix njobs = getNumberOfJobs();
        double totalJobs = 0; // Inf as soon as one class is open
        boolean anyOpen = false;
        for (int r = 0; r < njobs.length(); r++) {
            totalJobs += njobs.get(r);
            if (Double.isInfinite(njobs.get(r))) {
                anyOpen = true;
            }
        }
        for (Node node : this.nodes) {
            if (!(node instanceof Station) || node instanceof Source || node instanceof Sink) {
                continue;
            }
            Station station = (Station) node;
            // hasFiniteCap() decodes the three "unbounded" encodings (MAX_VALUE, Inf, and
            // JMT2LINE's negative sentinel) in one place; see Station.hasFiniteCap.
            if (station.hasFiniteCap() && station.getCap() < totalJobs) {
                return new BindingCapacity(station, station.getCap(), -1, null, anyOpen);
            }
            for (int r = 0; r < Math.min(this.jobClasses.size(), njobs.length()); r++) {
                JobClass jobClass = this.jobClasses.get(r);
                double classCap = station.getClassCap(jobClass);
                if (classCap > 0 && classCap < Integer.MAX_VALUE && classCap < njobs.get(r)) {
                    return new BindingCapacity(station, classCap, r, jobClass,
                            Double.isInfinite(njobs.get(r)));
                }
            }
        }
        return new BindingCapacity();
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
        // see _kb/04-networkstruct.md (Node/process construction notes: getUsedLangFeatures) for rationale
        if (!this.regions.isEmpty()) {
            setUsedLangFeature("Region");
        }
        // see _kb/04-networkstruct.md (Node/process construction notes: getUsedLangFeatures) for rationale
        for (JobClass jc : this.jobClasses) {
            SignalType signalType;
            DiscreteDistribution removalDist;
            RemovalPolicy removalPolicy;
            if (jc instanceof ClosedSignal) {
                setUsedLangFeature("ClosedSignal");
                signalType = ((ClosedSignal) jc).getSignalType();
                removalDist = ((ClosedSignal) jc).getRemovalDistribution();
                removalPolicy = ((ClosedSignal) jc).getRemovalPolicy();
            } else if (jc instanceof OpenSignal) {
                setUsedLangFeature("OpenSignal");
                signalType = ((OpenSignal) jc).getSignalType();
                removalDist = ((OpenSignal) jc).getRemovalDistribution();
                removalPolicy = ((OpenSignal) jc).getRemovalPolicy();
            } else if (jc instanceof Signal) {
                // unresolved placeholder: resolves by presence of a Source
                setUsedLangFeature(this.getIndexSourceNode() < 0 ? "ClosedSignal" : "OpenSignal");
                signalType = ((Signal) jc).getSignalType();
                removalDist = ((Signal) jc).getRemovalDistribution();
                removalPolicy = ((Signal) jc).getRemovalPolicy();
            } else {
                continue;
            }
            if (signalType != null) {
                switch (signalType) {
                    case NEGATIVE:
                        setUsedLangFeature("SignalType_NEGATIVE");
                        break;
                    case REPLY:
                        setUsedLangFeature("SignalType_REPLY");
                        break;
                    case CATASTROPHE:
                        setUsedLangFeature("SignalType_CATASTROPHE");
                        break;
                }
            }
            // see _kb/04-networkstruct.md (Node/process construction notes: getUsedLangFeatures) for rationale
            if (removalDist != null) {
                setUsedLangFeature("SignalBatchRemoval");
            }
            if (removalPolicy != null && removalPolicy != RemovalPolicy.RANDOM) {
                setUsedLangFeature("SignalRemovalPolicy");
            }
        }
        for (int i = 0; i < getNumberOfNodes(); i++) {
            for (int r = 0; r < getNumberOfClasses(); r++) {
                Node n = this.nodes.get(i);
                if (n instanceof Queue || n instanceof Delay) {
                    ServiceBinding serviceProcess = n.getServer().getServiceProcess(this.getClassByIndex(r));
                    if (serviceProcess != null) {
                        if (!(serviceProcess.getDistribution() instanceof Disabled) && !(serviceProcess.getDistribution() instanceof Immediate)) {
                            setUsedLangFeature(serviceProcess.getDistribution().getFeatureName());
                        }
                        // A finite-server station serving several jobs at once. A Delay
                        // carries Integer.MAX_VALUE here and is NOT one: the single-server
                        // recursions answered a c-server station as one server of the same
                        // rate, and only a structural predicate could say so. Marked once
                        // per model; setUsedLangFeature is idempotent.
                        int nserv = ((Station) n).getNumberOfServers();
                        if (nserv > 1 && nserv < Integer.MAX_VALUE) {
                            setUsedLangFeature("MultiServer");
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
                        setUsedLangFeature(serviceProcess.getFeatureName());
                    }
                    setUsedLangFeature("Source");
                    // see _kb/04-networkstruct.md (Node/process construction notes: getUsedLangFeatures) for rationale
                    if (((Source) n).getArrivalBatch(this.getClassByIndex(r)) != null) {
                        setUsedLangFeature("BatchArrival");
                    }
                } else if (n instanceof ClassSwitch) {
                    setUsedLangFeature("StatelessClassSwitcher");
                    setUsedLangFeature("ClassSwitch");
                } else if (n instanceof Fork) {
                    setUsedLangFeature("Fork");
                    setUsedLangFeature("Forker");
                    // variable forking levels: a name declared but never marked
                    // is a name no solver can ever refuse, so mark all three here
                    Forker fk = (Forker) n.getOutput();
                    if (!fk.tasksPerLinkByDest.isEmpty()) setUsedLangFeature("ForkFanoutVector");
                    if (!fk.tasksPerLinkDist.isEmpty()) setUsedLangFeature("ForkFanoutRandom");
                    if (!fk.branchProb.isEmpty()) setUsedLangFeature("ForkBranchProbability");
                } else if (n instanceof Join) {
                    setUsedLangFeature("Join");
                    setUsedLangFeature("Joiner");
                } else if (n instanceof Sink) {
                    setUsedLangFeature("JobSink");
                } else if (n instanceof Cache) {
                    setUsedLangFeature("CacheClassSwitcher");
                    setUsedLangFeature("Cache");
                    String replStrat = ReplacementStrategy.toFeature(((Cache) n).getReplacementStrategy());
                    if (replStrat != null && !replStrat.isEmpty()) {
                        setUsedLangFeature(replStrat);
                    }
                    // per-list storage cost caps with per-item sizes
                    if (((Cache) n).getCostCaps() != null && !((Cache) n).getCostCaps().isEmpty()) {
                        setUsedLangFeature("CacheItemSize");
                    }
                    // see _kb/04-networkstruct.md (Node/process construction notes: getUsedLangFeatures) for rationale
                    if (!((Cache) n).getRetrievalClassIndices().isEmpty()) {
                        setUsedLangFeature("CacheRetrieval");
                    }
                } else if (n instanceof Transition) {
                    setUsedLangFeature("Transition");
                    setUsedLangFeature("Enabling");
                    setUsedLangFeature("Timing");
                    setUsedLangFeature("Firing");
                    // Inhibitor arcs: flag only when a finite threshold is set,
                    // so plain SPNs are not gated out of solvers lacking it.
                    Map<Mode, Matrix> inhCond = ((Transition) n).inhibitingConditions;
                    if (inhCond != null) {
                        boolean hasInhibitor = false;
                        for (Matrix inhM : inhCond.values()) {
                            if (inhM == null) continue;
                            for (int ii = 0; ii < inhM.getNumRows() && !hasInhibitor; ii++) {
                                for (int jj = 0; jj < inhM.getNumCols(); jj++) {
                                    double v = inhM.get(ii, jj);
                                    if (!Double.isInfinite(v) && !Double.isNaN(v)) {
                                        hasInhibitor = true;
                                        break;
                                    }
                                }
                            }
                            if (hasInhibitor) break;
                        }
                        if (hasInhibitor) {
                            setUsedLangFeature("Inhibiting");
                        }
                    }
                } else if (n instanceof Place) {
                    setUsedLangFeature("Storage");
                    setUsedLangFeature("Linkage");
                    setUsedLangFeature("Place");
                    if (((Place) n).isQueueing()) {
                        // see _kb/04-networkstruct.md (Node/process construction notes: getUsedLangFeatures) for rationale
                        setUsedLangFeature("QueueingPlace");
                    }
                }
            }
        }
        // see _kb/04-networkstruct.md (Node/process construction notes: getUsedLangFeatures) for rationale
        for (Station station : this.stations) {
            if (station.getLimitedLoadDependence() != null && !station.getLimitedLoadDependence().isEmpty()) {
                setUsedLangFeature("LoadDependence");
                break;
            }
        }
        for (Station station : this.stations) {
            if (station.getLimitedClassDependence() != null) {
                setUsedLangFeature("ClassDependence");
                break;
            }
        }
        // Joint-dependent (non-product-form) scaling: registered so solvers not
        // plumbing eta_i(n) reject the model rather than silently ignoring it.
        for (Station station : this.stations) {
            if (station.getLimitedJointDependence() != null) {
                setUsedLangFeature("JointDependence");
                break;
            }
        }
        // Globally state-dependent scaling phi(n) over the full network state, the
        // Whittle primitive. Only SolverCTMC plumbs it, so every other solver must
        // reject the model rather than solve it unscaled.
        if (this.gdScaling != null) {
            setUsedLangFeature("GlobalDependence");
        }
        // see _kb/04-networkstruct.md (Node/process construction notes: getUsedLangFeatures) for rationale
        for (Station station : this.stations) {
            if (station instanceof Queue && ((Queue) station).isDelayOffEnabled()) {
                setUsedLangFeature("SetupDelayOff");
                break;
            }
        }
        // Server parallelism: a job seizing n>1 servers changes the effective
        // capacity of the station, so a solver that cannot honour it must reject
        // the model rather than solve it as if every job seized one server.
        for (Station station : this.stations) {
            if (station instanceof Queue && ((Queue) station).hasServerParallelism()) {
                setUsedLangFeature("ServerParallelism");
                break;
            }
        }
        // see _kb/04-networkstruct.md (Node/process construction notes: getUsedLangFeatures) for rationale
        for (Station station : this.stations) {
            boolean reneging = false;
            for (JobClass jobClass : this.jobClasses) {
                if (station.getImpatienceType(jobClass) == ImpatienceType.RENEGING) {
                    reneging = true;
                    break;
                }
            }
            if (reneging) {
                setUsedLangFeature("Reneging");
                break;
            }
        }
        for (Station station : this.stations) {
            boolean balking = false;
            for (JobClass jobClass : this.jobClasses) {
                if (station.getBalkingStrategy(jobClass) != null) {
                    balking = true;
                    break;
                }
            }
            if (balking) {
                setUsedLangFeature("Balking");
                break;
            }
        }
        // see _kb/04-networkstruct.md (Node/process construction notes: getUsedLangFeatures) for rationale
        for (Node node : this.nodes) {
            if (node instanceof Queue && ((Queue) node).hasBreakdown()) {
                setUsedLangFeature("Breakdown");
                break;
            }
        }
        // A retrial orbit (Station.setRetrial / setOrbit): the same per-class test
        // refreshRetrial makes for sn.retrialProc, a configured delay that is not
        // the Disabled placeholder, and restricted to a Queue as refreshStruct is.
        // A solver that reads no sn.retrial* field would answer the model with the
        // refused jobs simply lost, so the orbit is gated by name.
        for (Node node : this.nodes) {
            if (!(node instanceof Queue)) {
                continue;
            }
            boolean retrial = false;
            for (JobClass jobClass : this.jobClasses) {
                if (((Queue) node).hasRetrial(jobClass)) {
                    retrial = true;
                    break;
                }
            }
            if (retrial) {
                setUsedLangFeature("Retrial");
                break;
            }
        }
        // A genuine QUORUM join (PARTIAL, k of n siblings with k < n) fires on the
        // k-th branch completion, which is not their maximum, and the n-k stragglers
        // are discarded on arrival: a solver that serves it as a full join returns a
        // silent wrong answer, so the rule is gated by its own name. k >= n and
        // k <= 0 are full joins and are not flagged.
        if (!this.nodes.isEmpty()) {
            Matrix conn = getConnectionMatrix();
            for (Node node : this.nodes) {
                if (!(node instanceof Join) || !(node.getInput() instanceof Joiner)) {
                    continue;
                }
                Joiner joiner = (Joiner) node.getInput();
                // siblings are counted at the FORK, as the engines do: its out-degree
                // times tasksPerLink. The join in-degree is the fallback.
                int joinIdx = getNodeIndex(node);
                double nsib = 0;
                for (int i = 0; i < conn.getNumRows(); i++) {
                    if (conn.get(i, joinIdx) > 0) nsib++;
                }
                Node forkNode = ((Join) node).joinOf;
                if (forkNode instanceof Fork && forkNode.getOutput() instanceof Forker) {
                    int forkIdx = getNodeIndex(forkNode);
                    double w = Math.max(1.0, Math.round(((Forker) forkNode.getOutput()).tasksPerLink));
                    double outdeg = 0;
                    for (int j = 0; j < conn.getNumCols(); j++) {
                        if (conn.get(forkIdx, j) > 0) outdeg++;
                    }
                    nsib = outdeg * w;
                }
                for (JobClass jobClass : this.jobClasses) {
                    JoinStrategy js = joiner.joinStrategy.get(jobClass);
                    Double kreq = joiner.joinRequired.get(jobClass);
                    if (js == null || js == JoinStrategy.STD || kreq == null || kreq <= 0) {
                        continue;
                    }
                    if (nsib <= 0 || kreq < nsib) {
                        setUsedLangFeature("JoinPartial");
                    }
                }
            }
        }
        // Heterogeneous server pools (Queue.addServerType): the station is
        // served by several pools with their own counts, class compatibilities
        // and per-(type,class) rates. A solver that reads only sn.nservers
        // answers for a homogeneous station of the same total size, which is a
        // different system, so the pools are gated by their own name rather
        // than silently flattened.
        for (Station station : this.stations) {
            if (station instanceof Queue && !((Queue) station).getServerTypes().isEmpty()) {
                setUsedLangFeature("HeteroServers");
                break;
            }
        }
        // A FIFO depository (Place.setDepartureDiscipline) releases a served
        // token only after the tokens that entered service before it, so which
        // output transitions are enabled depends on the arrival order and not
        // only on the marking. No solver implements it, hence a clean refusal
        // instead of a Normal depository's answer under the user's name.
        for (Node node : this.nodes) {
            if (!(node instanceof Place)) {
                continue;
            }
            Place place = (Place) node;
            boolean nonNormal = false;
            for (JobClass jobClass : this.jobClasses) {
                if (place.getDepartureDiscipline(jobClass) != DepartureDiscipline.Normal) {
                    nonNormal = true;
                    break;
                }
            }
            if (nonNormal) {
                setUsedLangFeature("DepartureDiscipline");
                break;
            }
        }
        // A station or per-class buffer that can BIND: the one predicate
        // NetworkSolver.bindingCapacityReason refuses on (node-level caps against
        // the class populations, open classes always bind, Cache models exempt),
        // asked here so the refusal has a registry name and a solver method that
        // does not declare it is gated on exactly the models the structural gate
        // refuses.
        if (findBindingCapacity().binds) {
            setUsedLangFeature("FiniteCapacity");
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

    /**
     * Returns a directed-graph view of the network topology. Nodes are the
     * network node names and edges carry the per-class routing probability
     * (weight) between nodes.
     *
     * @return a {@link Graph} describing the routing topology
     */
    public Graph getGraph() {
        List<String> nodeNames = new ArrayList<String>();
        for (Node node : this.nodes) {
            nodeNames.add(node.getName());
        }
        List<Graph.Edge> edges = new ArrayList<Graph.Edge>();
        for (Node node : this.nodes) {
            for (OutputStrategy os : node.getOutputStrategies()) {
                if (os.getDestination() != null && os.getProbability() > 0) {
                    String cls = os.getJobClass() != null ? os.getJobClass().getName() : "";
                    edges.add(new Graph.Edge(node.getName(),
                        os.getDestination().getName(), cls, os.getProbability()));
                }
            }
        }
        return new Graph(nodeNames, edges);
    }

    /**
     * Removes the specified job class from this model, updating all node
     * configurations (service, capacity, routing, arrival and class-switching)
     * accordingly, and re-initialising the model.
     *
     * @param jobclass the job class to remove
     */
    public void removeClass(JobClass jobclass) {
        // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
        JobClass target = null;
        for (JobClass candidate : this.jobClasses) {
            if (candidate == jobclass) {
                target = candidate;
                break;
            }
        }
        if (target == null && jobclass != null && jobclass.getName() != null) {
            for (JobClass candidate : this.jobClasses) {
                if (jobclass.getName().equals(candidate.getName())) {
                    target = candidate;
                    break;
                }
            }
        }
        if (this.jobClasses.size() <= 1) {
            if (target != null) {
                throw new RuntimeException(
                    "The network has a single class, it cannot be removed from the model.");
            }
            return;
        }
        if (target == null) {
            return;
        }
        int r = this.jobClasses.indexOf(target);
        for (Node node : this.nodes) {
            node.removeJobClass(target);
        }
        this.jobClasses.remove(target);
        // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
        int K = this.jobClasses.size();
        if (this.csMatrix != null && this.csMatrix.getNumRows() == K + 1 && r >= 0) {
            Matrix reduced = new Matrix(K, K);
            for (int i = 0; i < K; i++) {
                int oldRow = i < r ? i : i + 1;
                for (int j = 0; j < K; j++) {
                    int oldCol = j < r ? j : j + 1;
                    reduced.set(i, j, this.csMatrix.get(oldRow, oldCol));
                }
            }
            this.csMatrix = reduced;
        }
        this.reset(true);
    }

    /**
     * Returns a copy of this model in which all classes belonging to the same
     * chain are merged into a single aggregate class, so that the aggregated
     * model has one class per chain of this model. Class switching is
     * eliminated in the process.
     *
     * <p>The result also carries the (nstations x nclasses) matrix of
     * aggregation factors alpha and the deaggregation record needed to map
     * chain-level metrics back to class-level metrics with
     * {@code SnDeaggregateChainResults.snDeaggregateChainResults}. The
     * aggregation is exact for product-form models and approximate otherwise,
     * since one aggregate service process replaces the per-class ones weighted
     * by alpha.
     *
     * @return the aggregated model, alpha, and the deaggregation record
     */
    public ModelAdapter.AggregateChainResult aggregateChains() {
        return ModelAdapter.aggregateChains(this, "");
    }

    /**
     * Chain-aggregated copy of this model, with a suffix appended to the names
     * of the aggregate classes. See {@link #aggregateChains()}.
     *
     * @param suffix suffix for the aggregate class names, may be null
     * @return the aggregated model, alpha, and the deaggregation record
     */
    public ModelAdapter.AggregateChainResult aggregateChains(String suffix) {
        return ModelAdapter.aggregateChains(this, suffix);
    }

    /**
     * Returns a copy of this model with the given job class removed, leaving
     * this model untouched. Use it when the original model must stay solvable,
     * for instance in ablation studies or a per-class decomposition.
     *
     * @param jobclass the job class to remove from the copy
     * @return a new model without the specified class
     */
    public Network withoutClass(JobClass jobclass) {
        return ModelAdapter.removeClass(this, jobclass);
    }

    /**
     * Check if the queueing network routing matrix is ergodic (irreducible).
     *
     * This checks only the routing structure, not the full CTMC state space.
     * A routing is ergodic if all stations communicate, meaning the routing
     * matrix does not create absorbing states or disconnected components.
     *
     * @return RoutingErgodicityResult containing isErgodic flag and details
     */
    public RoutingErgodicityResult isRoutingErgodic() {
        return isRoutingErgodic((RoutingMatrix) null);
    }

    /**
     * Check if the queueing network routing matrix is ergodic (irreducible).
     *
     * @param P Optional RoutingMatrix. If null, will be computed from network structure.
     * @return RoutingErgodicityResult containing isErgodic flag and details
     */
    public RoutingErgodicityResult isRoutingErgodic(RoutingMatrix P) {
        RoutingErgodicityResult info = new RoutingErgodicityResult();

        int K = this.getNumberOfClasses();
        int I = this.getNumberOfNodes();

        if (K == 0 || I == 0) {
            info.isErgodic = true;
            return info;
        }

        List<String> nodeNames = this.getNodeNames();

        // Build aggregate adjacency matrix across all classes
        Matrix adjMatrix = new Matrix(I, I);

        if (P != null) {
            // Use provided RoutingMatrix
            // Note: RoutingMatrix.get(int, int) expects 1-based indices
            for (int r = 1; r <= K; r++) {
                for (int s = 1; s <= K; s++) {
                    Matrix Prs = P.get(r, s);
                    if (Prs != null && !Prs.isEmpty()) {
                        int nRows = Math.min(I, Prs.getNumRows());
                        int nCols = Math.min(I, Prs.getNumCols());
                        for (int i = 0; i < nRows; i++) {
                            for (int j = 0; j < nCols; j++) {
                                if (Prs.get(i, j) > 0) {
                                    adjMatrix.set(i, j, 1.0);
                                }
                            }
                        }
                    }
                }
            }
        } else {
            // Use getLinkedRoutingMatrix (requires struct)
            Map<JobClass, Map<JobClass, Matrix>> rtorig = getLinkedRoutingMatrix();
            if (rtorig == null || rtorig.isEmpty()) {
                info.isErgodic = true;
                return info;
            }
            for (Map.Entry<JobClass, Map<JobClass, Matrix>> entry1 : rtorig.entrySet()) {
                for (Map.Entry<JobClass, Matrix> entry2 : entry1.getValue().entrySet()) {
                    Matrix Prs = entry2.getValue();
                    if (Prs != null && !Prs.isEmpty()) {
                        int nRows = Math.min(I, Prs.getNumRows());
                        int nCols = Math.min(I, Prs.getNumCols());
                        for (int i = 0; i < nRows; i++) {
                            for (int j = 0; j < nCols; j++) {
                                if (Prs.get(i, j) > 0) {
                                    adjMatrix.set(i, j, 1.0);
                                }
                            }
                        }
                    }
                }
            }
        }

        // Find strongly connected components using DirectedGraph
        jline.util.graph.DirectedGraph graph = new jline.util.graph.DirectedGraph(adjMatrix);
        jline.util.graph.DirectedGraph.SCCResult sccResult = graph.stronglyconncomp();
        int[] sccLabels = sccResult.I;
        boolean[] isRecurrent = sccResult.recurrent;

        int numSCCs = isRecurrent.length;
        info.numSCCs = numSCCs;

        // Get station indexes (node index for each station)
        List<Integer> stationIdxs = new ArrayList<>();
        for (int i = 0; i < this.nodes.size(); i++) {
            if (this.nodes.get(i) instanceof Station) {
                stationIdxs.add(i);
            }
        }

        // Check for absorbing stations (stations with only self-loops or no outgoing edges)
        for (int idx : stationIdxs) {
            // Check if this station has any outgoing transitions to other nodes
            boolean hasOutgoing = false;
            for (int j = 0; j < I; j++) {
                if (j != idx && adjMatrix.get(idx, j) > 0) {
                    hasOutgoing = true;
                    break;
                }
            }

            if (!hasOutgoing) {
                // No outgoing transitions except possibly self-loop
                // Check if there's routing defined for this node
                boolean hasRouting = false;
                for (int i = 0; i < I; i++) {
                    if (adjMatrix.get(idx, i) > 0) {
                        hasRouting = true;
                        break;
                    }
                }

                if (hasRouting) {
                    info.absorbingStations.add(nodeNames.get(idx));
                }
            }
        }

        // Identify transient stations (stations in transient SCCs)
        for (int i = 0; i < numSCCs; i++) {
            if (!isRecurrent[i]) {
                // This SCC is transient - find nodes in it
                for (int nodeIdx = 0; nodeIdx < I; nodeIdx++) {
                    if (sccLabels[nodeIdx] == i + 1 && stationIdxs.contains(nodeIdx)) {
                        String nodeName = nodeNames.get(nodeIdx);
                        if (!info.absorbingStations.contains(nodeName)) {
                            info.transientStations.add(nodeName);
                        }
                    }
                }
            }
        }

        // Remove Sink from absorbing list (it's expected to be absorbing in open networks)
        int sinkIdx = this.getIndexSinkNode();
        if (sinkIdx >= 0 && sinkIdx < nodeNames.size()) {
            String sinkName = nodeNames.get(sinkIdx);
            info.absorbingStations.remove(sinkName);
        }

        // Determine ergodicity
        boolean hasAbsorbingStations = !info.absorbingStations.isEmpty();
        int recurrentCount = 0;
        for (boolean rec : isRecurrent) {
            if (rec) recurrentCount++;
        }
        boolean hasMultipleRecurrentSCCs = recurrentCount > 1;

        if (hasAbsorbingStations || hasMultipleRecurrentSCCs) {
            info.isErgodic = false;
            info.isReducible = true;
        } else {
            info.isErgodic = true;
            info.isReducible = false;
        }

        return info;
    }

    /**
     * Result class for isRoutingErgodic method
     */
    public static class RoutingErgodicityResult {
        public boolean isErgodic = true;
        public boolean isReducible = false;
        public List<String> absorbingStations = new ArrayList<>();
        public List<String> transientStations = new ArrayList<>();
        public int numSCCs = 1;
    }

    // ========================================================================
    // SECTION 17: INITIALIZATION METHODS
    // Methods for initializing network parameters and structures
    // ========================================================================

    public void initDefault() {
        // see _kb/04-networkstruct.md (initDefault.m/spaceGenerator.m section) for rationale

        NetworkStruct sn = this.getStruct(false);
        int R = sn.nclasses;
        Matrix N = sn.njobs.transpose();

        // see _kb/04-networkstruct.md (initDefault.m/spaceGenerator.m section) for rationale
        Matrix nplace = new Matrix(sn.nstations, R);
        nplace.zero();
        double[] totplace = new double[sn.nstations];
        for (int r = 0; r < N.getNumRows(); r++) {
            if (!isFinite(N.get(r, 0))) {
                continue;
            }
            int refist = (int) sn.refstat.get(r, 0);
            if (sn.nodetype.get((int) sn.stationToNode.get(refist)) == NodeType.Place) {
                nplace.set(refist, r, N.get(r, 0));
                totplace[refist] += N.get(r, 0);
                continue;
            }
            double remaining = N.get(r, 0);
            for (int off = 0; off < sn.nstations && remaining > 0; off++) {
                int jst = (off == 0) ? refist : (off <= refist ? off - 1 : off);
                if (sn.sched.get(sn.stations.get(jst)) == SchedStrategy.EXT
                        || sn.nodetype.get((int) sn.stationToNode.get(jst)) == NodeType.Place) {
                    continue;
                }
                double avail = FastMath.min(sn.classcap.get(jst, r) - nplace.get(jst, r),
                        sn.cap.get(jst, 0) - totplace[jst]);
                double take = FastMath.min(remaining, FastMath.max(0, avail));
                nplace.set(jst, r, nplace.get(jst, r) + take);
                totplace[jst] += take;
                remaining -= take;
            }
            if (remaining > 0) {
                throw new RuntimeException("initDefault: Cannot place the population of class "
                        + r + ": total station capacity is insufficient.");
            }
        }

        for (int ind = 0; ind < this.getNumberOfNodes(); ind++) {
            // Check if user has already set a state for this node - if so, preserve it
            if (this.nodes.get(ind).isStateful()) {
                Matrix existingState = ((StatefulNode) this.nodes.get(ind)).getState();
                if (existingState != null && !existingState.isEmpty()) {
                    boolean isPasStation = sn.isstation.get(ind, 0) == 1
                            && sn.sched.get((Station) this.nodes.get(ind)) == SchedStrategy.PAS;
                    if (isPasStation && existingState.elementSum() > 0) {
                        // see _kb/04-networkstruct.md (initDefault.m/spaceGenerator.m section) for rationale
                        StatefulNode pasNode = (StatefulNode) this.nodes.get(ind);
                        int istp = (int) sn.nodeToStation.get(ind);
                        Matrix nUser = new Matrix(1, R);
                        nUser.zero();
                        for (int c = 0; c < existingState.getNumCols(); c++) {
                            int cls = (int) existingState.get(0, c);
                            if (cls >= 1 && cls <= R) {
                                nUser.set(0, cls - 1, nUser.get(0, cls - 1) + 1);
                            }
                        }
                        Matrix sU = new Matrix(1, R);
                        sU.zero();
                        double ssU = sn.nservers.get(istp, 0);
                        for (int r = 0; r < R; r++) {
                            double v = FastMath.min(nUser.get(0, r), ssU);
                            sU.set(0, r, v);
                            ssU -= v;
                        }
                        Matrix space_i = FromMarginal.fromMarginalAndStarted(sn, ind, nUser, sU);
                        int W = space_i.getNumCols();
                        Matrix urow = new Matrix(1, W);
                        urow.zero();
                        for (int c = 0; c < Math.min(existingState.getNumCols(), W); c++) {
                            urow.set(0, c, existingState.get(0, c));
                        }
                        int uidx = Matrix.matchrow(space_i, urow);
                        if (uidx < 0) {
                            uidx = 0;
                        }
                        Matrix prior = new Matrix(space_i.getNumRows(), 1);
                        prior.zero();
                        prior.set(uidx, 0, 1.0);
                        pasNode.setStateSpace(space_i);
                        pasNode.setStatePrior(prior);
                        pasNode.setState(urow);
                        continue;
                    }
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
                        n0.set(0, r, nplace.get((int) sn.nodeToStation.get(ind), r));
                    }
                    s0.set(0, r, FastMath.min(n0.get(0, r), s));
                    s -= s0.get(0, r);
                }

                switch (sn.nodetype.get(ind)) {
                    case Cache:
                        int totalCacheCapacity = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).totalCacheCapacity;
                        state_i = FromMarginal.fromMarginalAndStarted(sn, ind, n0, s0);
                        int nVarsInt = (int) sn.nvars.get(ind, 2 * R);
                        Matrix newState_i = new Matrix(1, state_i.getNumCols() + nVarsInt);
                        for (int p = 0; p < state_i.length(); p++) {
                            newState_i.set(0, p, state_i.get(0, p));
                        }
                        int addition = 0;
                        for (int p = state_i.length(); p < state_i.length() + totalCacheCapacity; p++) {
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
                                    state_i.set(0, r, sn.njobs.get(0, r));
                                } else {
                                    state_i.set(0, r, 0);
                                }
                            }
                        }
                        break;
                    }
                    default:
                        if (sn.isstation.get(ind, 0) == 1
                                && sn.sched.get((Station) this.nodes.get(ind)) == SchedStrategy.PAS) {
                            // see _kb/04-networkstruct.md (initDefault.m/spaceGenerator.m section) for rationale
                            QueueNodeParam qnp = (QueueNodeParam) sn.nodeparam.get(this.nodes.get(ind));
                            boolean hasSwap = qnp != null && qnp.swapGraph != null
                                    && !qnp.swapGraph.isEmpty() && qnp.swapGraph.elementMaxAbs() > 0;
                            boolean isClosed = false;
                            for (int r = 0; r < N.getNumRows(); r++) {
                                if (isFinite(N.get(r, 0))) {
                                    isClosed = true;
                                    break;
                                }
                            }
                            if (hasSwap && isClosed && n0.elementSum() > 1) {
                                throw new RuntimeException("A closed pass-and-swap station with a non-empty swapping "
                                        + "graph requires an explicit initial job placement (node " + ind + "). Call "
                                        + "setState on the station with the ordered class list (oldest first) before solving.");
                            }
                            state_i = FromMarginal.fromMarginalAndStarted(sn, ind, n0, s0);
                            break;
                        }
                        state_i = FromMarginal.fromMarginalAndStarted(sn, ind, n0, s0);
                }

                // The synchronous-call (REPLY) counters are the LAST local variables (see
                // ReplyBlock), but the modulation, routing and node blocks below are
                // appended after them. Detach the counters here and re-attach them at the
                // tail, otherwise the initial row is a column permutation of every
                // enumerated row, matchrow fails, and Solver_ctmc silently skips its
                // unreachable-state pruning.
                Matrix replyCols = new Matrix(0, 0);
                if (jline.lang.state.ReplyBlock.holds(sn, ind)) {
                    int replyw = jline.lang.state.ReplyBlock.info(sn, ind).width;
                    if (replyw > 0 && state_i.getNumCols() >= replyw) {
                        int wkeep = state_i.getNumCols() - replyw;
                        replyCols = Matrix.extract(state_i, 0, state_i.getNumRows(), wkeep, state_i.getNumCols());
                        state_i = Matrix.extract(state_i, 0, state_i.getNumRows(), 0, wkeep);
                    }
                }

                // Markov-modulated service keeps a phase-restart slot per class (nvars
                // columns 0..R-1, allocated by refreshLocalVars for MAP/DMAP/MMPP2/BMAP);
                // the enumerated space stores it 1-based, so the server starts in phase 1.
                if (sn.isstation.get(ind, 0) == 1 && sn.nvars != null && !sn.nvars.isEmpty()) {
                    for (int r = 0; r < sn.nclasses; r++) {
                        for (int v = 0; v < (int) sn.nvars.get(ind, r); v++) {
                            state_i = state_i.concatCols(Matrix.singleton(1));
                        }
                    }
                }

                for (int r = 0; r < sn.nclasses; r++) {
                    if (sn.routing.get(sn.nodes.get(ind)).get(sn.jobclasses.get(r)) == RoutingStrategy.WRROBIN) {
                        // see _kb/04-networkstruct.md (State-space construction conventions) for rationale
                        state_i = state_i.concatCols(Matrix.singleton(1));
                    } else if (sn.routing.get(sn.nodes.get(ind)).get(sn.jobclasses.get(r)) == RoutingStrategy.RROBIN) {
                        // Start from first connected queue
                        for (int p = 0; p < sn.connmatrix.getNumCols(); p++) {
                            if (sn.connmatrix.get(ind, p) == 1) {
                                state_i = state_i.concatCols(Matrix.singleton(p));
                                break;
                            }
                        }
                    }
                }

                // see _kb/04-networkstruct.md (refreshStruct.m field-population notes) for rationale
                if (sn.isstation.get(ind, 0) > 0
                        && sn.sched.get(sn.stations.get((int) sn.nodeToStation.get(0, ind))) == SchedStrategy.POLLING) {
                    jline.lang.state.Polling.Info pinfoI = jline.lang.state.Polling.info(sn, ind);
                    if (pinfoI != null && pinfoI.width > 0) {
                        int srvclassI = -1;
                        for (int r = 0; r < sn.nclasses; r++) {
                            if (s0.get(0, r) > 0) { srvclassI = r; break; }
                        }
                        int[] nbufI = new int[sn.nclasses];
                        for (int r = 0; r < sn.nclasses; r++) {
                            nbufI[r] = (int) Math.round(n0.get(0, r));
                        }
                        if (srvclassI >= 0) {
                            nbufI[srvclassI] -= 1;
                        }
                        int posI, swkI, ctrI;
                        if (srvclassI >= 0) {
                            posI = srvclassI;
                            swkI = 0;
                            ctrI = jline.lang.state.Polling.budget(pinfoI, nbufI[srvclassI] + 1);
                        } else {
                            int[] resI = jline.lang.state.Polling.next(pinfoI, 0, nbufI, sn.nclasses, true);
                            posI = resI[0];
                            ctrI = 0;
                            swkI = 0;
                            if (resI[1] == jline.lang.state.Polling.MODE_SWITCH) {
                                Matrix swpieI = pinfoI.swPie[posI];
                                for (int kk = 0; kk < swpieI.length(); kk++) {
                                    if (swpieI.get(kk) > 0) { swkI = kk + 1; break; }
                                }
                            }
                        }
                        state_i = state_i.concatCols(jline.lang.state.Polling.project(pinfoI, new int[]{posI, swkI, ctrI}));
                    }
                }

                // No call is outstanding at time zero, so the counters are zero, but the
                // columns must trail every other local variable (see ReplyBlock).
                if (!replyCols.isEmpty()) {
                    state_i = state_i.concatCols(replyCols);
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
                    int numClasses = this.getNumberOfClasses();
                    int cacheCapacity = cacheNode.getTotalCacheCapacity();
                    int retrievalSystemCapacity = cacheNode.getRetrievalSystemCapacity();
                    int nItems = ((CacheNodeParam) sn.nodeparam.get(this.nodes.get(ind))).nitems;
                    // When a retrieval system is configured the local state appends block A (a
                    // per-item occupancy bitmap) and block B (per-retrieval-class delayed-hit
                    // counts); otherwise only the cache contents are stored.
                    int retrievalBitmapWidth = (retrievalSystemCapacity > 0) ? nItems : 0;
                    int retrievalPendingWidth = (retrievalSystemCapacity > 0)
                            ? cacheNode.getRetrievalClassIndices().size() : 0;

                    // Cache state appears as a flattened vector of [Classes | cache contents | block A | block B]
                    state_i = new Matrix(1, numClasses + cacheCapacity + retrievalBitmapWidth + retrievalPendingWidth);
                    int ctr = 0;
                    for (int idx = numClasses; idx < numClasses + cacheCapacity; idx++) {
                        state_i.set(idx, ++ctr);
                    }
                } else if (this.nodes.get(ind) instanceof Router) {
                    state_i = Matrix.zeros(1, this.getNumberOfClasses());
                    // A router is stateful precisely because a cyclic strategy keeps a
                    // pointer there, and refreshLocalVars counts that pointer in nvars.
                    // Without the column the row is one short of what every reader
                    // derives from nvars, and the per-class marginal runs off its end.
                    for (int r = 0; r < sn.nclasses; r++) {
                        RoutingStrategy rstrat = sn.routing.get(sn.nodes.get(ind)).get(sn.jobclasses.get(r));
                        if (rstrat == RoutingStrategy.WRROBIN) {
                            state_i = state_i.concatCols(Matrix.singleton(1));
                        } else if (rstrat == RoutingStrategy.RROBIN) {
                            for (int p = 0; p < sn.connmatrix.getNumCols(); p++) {
                                if (sn.connmatrix.get(ind, p) == 1) {
                                    state_i = state_i.concatCols(Matrix.singleton(p));
                                    break;
                                }
                            }
                        }
                    }
                } else if (this.nodes.get(ind) instanceof StatefulFork) {
                    // Stateful Fork (FJ tag-augmented copies only): per-class
                    // count of parent jobs held before the fork firing
                    state_i = Matrix.zeros(1, this.getNumberOfClasses());
                } else if (this.nodes.get(ind) instanceof Transition) {
                    Transition tr = (Transition) this.nodes.get(ind);
                    // local buffer: [disabled-servers, zeros(firing-phases), fired-servers]
                    TransitionNodeParam transParam = (TransitionNodeParam) sn.nodeparam.get(tr);
                    state_i = transParam.nmodeservers.copy();
                    for (int k = 0; k < transParam.nmodes; k++) {
                        double fp = transParam.firingphases.get(k);
                        int fpInt = Double.isNaN(fp) ? 1 : (int) fp;
                        state_i = state_i.concatCols(Matrix.zeros(1, fpInt));
                    }
                    state_i = state_i.concatCols(Matrix.zeros(1, transParam.nmodes));
                    // middle region (firing phases) already all zeros
                    // last nmodes positions stay zero as well
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

        // see _kb/04-networkstruct.md (refreshStruct.m field-population notes) for rationale
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
                    // see _kb/04-networkstruct.md (refreshStruct.m field-population notes) for rationale
                    Matrix nodeRow = n.getRow(i);
                    double maxFracDiff = 0;
                    for (int c = 0; c < nodeRow.getNumCols(); c++) {
                        double val = nodeRow.get(0, c);
                        if (!Double.isInfinite(val) && !Double.isNaN(val)) {
                            maxFracDiff = Math.max(maxFracDiff, Math.abs(val - Math.round(val)));
                        }
                    }
                    if (maxFracDiff < GlobalConstants.CoarseTol) {
                        // Integer-like: use standard fromMarginal with rounded values
                        Matrix roundedRow = new Matrix(nodeRow);
                        for (int c = 0; c < roundedRow.getNumCols(); c++) {
                            double val = roundedRow.get(0, c);
                            if (!Double.isInfinite(val) && !Double.isNaN(val)) {
                                roundedRow.set(0, c, Math.round(val));
                            }
                        }
                        state_i = FromMarginal.fromMarginal(sn, i, roundedRow);
                    } else {
                        // Fractional: store raw values directly (for Fluid solver init)
                        state_i = nodeRow;
                    }
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
            // Match MATLAB: validate the freshly-assigned state (not the stale default)
            if (sn.isstateful.get(i) == 1 && this.nodes.get(i).getState().isEmpty()) {
                throw new RuntimeException("Invalid state assignment for station " + i + ".");
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

        // see _kb/04-networkstruct.md (Routing-matrix validation) for rationale
        if (!isReset && connections != null && !connections.isEmpty() && connections.elementSum() > 0) {
            line_warning(mfilename(new Object() {
            }), "The Network.link method cannot be used after calling the addLink() method. Use Node.setProbRouting instead to configure routing probabilities.");
        }
        
        sanitize();

        // see _kb/04-networkstruct.md (Routing-matrix validation) for rationale
        if (this.enableChecks) {
            validateRoutingProbabilities(P);
        }

        // see _kb/04-networkstruct.md (Routing-matrix validation) for rationale
        for (Node node : this.nodes) {
            if (!(node instanceof Station) && node.getOutput() instanceof Dispatcher) {
                ((Dispatcher) node.getOutput()).initDispatcherJobClassesPreservingRouted(this.getClasses());
            }
        }

        injectRetrievalRouting(P);

        P.setRouting(this);

        for (Node node : this.nodes) {
            if (node instanceof Place) {
                ((Place) node).init();
            }
        }

        // Check for reducible routing (absorbing states)
        // This matches MATLAB's link.m lines 304-312
        if (this.enableChecks) {
            RoutingErgodicityResult ergResult = isRoutingErgodic(P);
            if (!ergResult.isErgodic && !ergResult.absorbingStations.isEmpty()) {
                line_warning(mfilename(new Object() {}),
                    "Reducible network topology detected, results may be unreliable.");
            }
        }

        // Check that order-independent (OI) stations have a permutation-invariant
        // rate mu(c). Disabled by model.setChecks(false).
        boolean hasOI = false;
        for (Node node : this.nodes) {
            if (node instanceof Queue) {
                Queue q = (Queue) node;
                if (q.getSchedStrategy() == SchedStrategy.OI && q.getServiceRateFunction() != null) {
                    hasOI = true;
                    break;
                }
            }
        }
        if (this.enableChecks && hasOI) {
            int K = this.jobClasses.size();
            double[] pop = new double[K];
            double ntot = 0;
            for (int r = 0; r < K; r++) {
                pop[r] = this.jobClasses.get(r).getNumberOfJobs();
                ntot += pop[r];
            }
            // The conserved quantity under class switching is the CHAIN
            // population, not the per-class one: a class reaches counts up to
            // its chain's total, and a class declared empty and filled only by
            // a switch (a ClassSwitch target, a Cache hit/miss class) reaches
            // them from a declared population of 0. Bounding each class by its
            // OWN population leaves reachable microstates unenumerated -- every
            // mixed one, when that population is 0 -- and the check then passes
            // vacuously on plainly order-dependent rates. csMatrix is the mask
            // P.setRouting has just recorded above; its weakly connected
            // components are the chains. Without class switching each component
            // is a single class and Nvec falls back to pop. A stale mask (a
            // fork-join copy widens the class set) is ignored, as elsewhere.
            double[] Nvec = new double[K];
            System.arraycopy(pop, 0, Nvec, 0, K);
            if (hasUsableCsMatrix(K)) {
                for (Set<Integer> chain : Matrix.weaklyConnect(this.csMatrix, new HashSet<Integer>())) {
                    double chainPop = 0;
                    for (Integer r : chain) {
                        if (r >= 0 && r < K) chainPop += pop[r];
                    }
                    for (Integer r : chain) {
                        if (r >= 0 && r < K) Nvec[r] = chainPop;
                    }
                }
            }
            for (Node node : this.nodes) {
                if (node instanceof Queue) {
                    Queue q = (Queue) node;
                    if (q.getSchedStrategy() == SchedStrategy.OI && q.getServiceRateFunction() != null) {
                        // ntot caps the microstate LENGTH: with Nvec now
                        // chain-wide, summing it inside checkPermInvariance
                        // would count each chain once per class and admit
                        // microstates longer than the network can produce.
                        Queue.PermCheckResult pr = q.checkPermInvariance(Nvec, Math.min(q.getCap(), ntot));
                        if (!pr.ok) {
                            line_error(mfilename(new Object() {}),
                                "Order-independent (OI) station '" + q.getName() + "' has a service rate "
                                + "function that is not permutation-invariant: mu(c) differs for a reordering "
                                + "of the microstate " + java.util.Arrays.toString(pr.badc) + ". Use "
                                + "SchedStrategy.PAS for order-dependent service, or disable this check with "
                                + "model.setChecks(false).");
                        }
                        // (1): per-job rates must be non-negative. Runs second
                        // because permutation invariance is what makes mu a
                        // function of the count vector, which is what the
                        // increment test walks.
                        Queue.MonoCheckResult mr =
                                q.checkRateMonotonicity(Nvec, Math.min(q.getCap(), ntot));
                        if (!mr.ok) {
                            line_error(mfilename(new Object() {}),
                                "Order-independent (OI) station '" + q.getName() + "' has a service rate "
                                + "function with a negative per-job service rate: mu(c) DROPS when a class-"
                                + mr.badr + " job joins, at microstate "
                                + java.util.Arrays.toString(mr.badc) + ". An OI rate must satisfy "
                                + "mu(c[:j]) >= mu(c[:j-1]), so that every mu_j(c) is non-negative. A "
                                + "processor-sharing total rate (sum_j mu_{c_j})/n has this shape whenever "
                                + "classes have different rates -- use SchedStrategy.PS for that station, "
                                + "or disable this check with model.setChecks(false).");
                        } else if (pr.partial || mr.partial) {
                            line_warning(mfilename(new Object() {}),
                                "Order-independent (OI) station '" + q.getName() + "': the permutation-"
                                + "invariance check was only partial because the reachable population is "
                                + "large; a subset of microstates was verified. To skip this check, call "
                                + "model.setChecks(false) before link().");
                        }
                    }
                }
            }
        }

        // Note: MATLAB also calls refreshChains at the end if isReset is true,
        // but this happens automatically in Java during struct regeneration
    }

    /**
     * Inject retrieval-system routing into the routing matrix P before it is applied
     * (mirrors MATLAB link.m). For every Cache with a retrieval system: (1) lift the read
     * class's routing among the retrieval queues (edges drawn in P over {queues}u{cache})
     * into each item's auto-generated retrieval class as the per-item default, then consume
     * those read-class template edges so the read class's own flow is unchanged; (2) inject
     * the deferred explicit entries registered by the setItem methods / setRetrievalSystem,
     * which override the defaults (last-wins).
     */
    private void injectRetrievalRouting(RoutingMatrix P) {
        for (Node node : this.nodes) {
            if (!(node instanceof Cache)) {
                continue;
            }
            Cache cache = (Cache) node;
            Map<Integer, List<Integer>> qmap = cache.getRetrievalSystemQueueIndices();
            if (qmap == null || qmap.isEmpty()) {
                continue;
            }
            Matrix rclasses = cache.getRetrievalClasses();
            int nItems = cache.getNumberOfItems();

            for (Map.Entry<Integer, List<Integer>> e : qmap.entrySet()) {
                int readIdx0 = e.getKey();                       // read class index (0-based)
                JobClass readClass = this.getJobClassFromIndex(readIdx0);
                List<Integer> Q = e.getValue();                  // retrieval queue node indices

                // (1) default per-item routing inherited from the read class topology in P
                for (int it = 0; it < nItems; it++) {
                    int rcIdx = (int) rclasses.get(it, readIdx0);
                    if (rcIdx < 0) {
                        continue;
                    }
                    JobClass rClass = this.getJobClassFromIndex(rcIdx);
                    for (int qi = 0; qi < Q.size(); qi++) {
                        Node qn = this.getNodeByIndex(Q.get(qi));
                        double pEntry = P.get(readClass, readClass, cache, qn);   // cache -> queue
                        if (pEntry > 0) {
                            P.set(rClass, rClass, cache, qn, pEntry);
                        }
                        double pExit = P.get(readClass, readClass, qn, cache);    // queue -> cache
                        if (pExit > 0) {
                            P.set(rClass, rClass, qn, cache, pExit);
                        }
                        for (int qj = 0; qj < Q.size(); qj++) {                   // queue -> queue
                            Node qn2 = this.getNodeByIndex(Q.get(qj));
                            double pqq = P.get(readClass, readClass, qn, qn2);
                            if (pqq > 0) {
                                P.set(rClass, rClass, qn, qn2, pqq);
                            }
                        }
                    }
                }

                // consume the read class's template edges over the queue set
                for (int qi = 0; qi < Q.size(); qi++) {
                    Node qn = this.getNodeByIndex(Q.get(qi));
                    P.set(readClass, readClass, cache, qn, 0.0);
                    P.set(readClass, readClass, qn, cache, 0.0);
                    for (int qj = 0; qj < Q.size(); qj++) {
                        P.set(readClass, readClass, qn, this.getNodeByIndex(Q.get(qj)), 0.0);
                    }
                }
            }

            // (2) explicit deferred entries override the defaults
            for (Cache.RetrievalRoutingEntry ent : cache.getRetrievalRoutingEntries()) {
                P.set(ent.fromClass, ent.toClass, ent.srcNode, ent.destNode, ent.prob);
            }
        }
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
        jline.api.sn.SnPrintRoutingMatrix.snPrintRoutingMatrix(this.sn, null);
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

        // see _kb/04-networkstruct.md (droprule/classcap derivation) for rationale
        boolean[] spawnFedChain = new boolean[C];
        for (int r0 = 0; r0 < K; r0++) {
            int spawnIdx = this.jobClasses.get(r0).getSpawnClassIndex();
            if (spawnIdx < 1 || spawnIdx > K) continue;
            for (int c0 = 0; c0 < C; c0++) {
                Matrix inchain_c0 = this.sn.inchain.get(c0);
                for (int idx = 0; idx < inchain_c0.length(); idx++) {
                    if ((int) inchain_c0.get(idx) == spawnIdx - 1) {
                        spawnFedChain[c0] = true;
                    }
                }
            }
        }

        // A chain routed through a QUORUM join is not population-conserving either,
        // and for the same reason: the join releases the parent at the k-th of n
        // siblings and the n-k stragglers stay in the branches, so the parent forks
        // again while they are still in flight. Nothing bounds that backlog, so a
        // branch station holds no more than the class population only under a
        // STANDARD join; capping it at sum(njobs) makes a simulator drop a closed
        // job. Read off the node objects, not off sn.nodeparam: refreshCapacity
        // runs before refreshLocalVars rebuilds nodeparam.
        // see _kb/05-solvers-overview.md
        boolean[] quorumClass = new boolean[K];
        // Read this.connections DIRECTLY and never getConnectionMatrix(), which
        // ALLOCATES and EXPANDS it as a side effect. Doing that from inside a
        // refresh resized the connection matrix of a tag-augmented copy (which
        // carries more nodes than the model it came from) and broke State.fromMarginal
        // on it. A null or undersized matrix here just leaves nsib at 0, which the
        // quorum test below already treats as "sibling count unknown".
        Matrix connq = this.connections;
        if (!this.nodes.isEmpty()) {
            for (Node node : this.nodes) {
                if (!(node instanceof Join) || !(node.getInput() instanceof Joiner)) continue;
                Joiner joiner = (Joiner) node.getInput();
                int joinIdx = getNodeIndex(node);
                double nsib = 0;
                if (connq != null && joinIdx < connq.getNumCols()) {
                    for (int i = 0; i < connq.getNumRows(); i++) {
                        if (connq.get(i, joinIdx) > 0) nsib++;
                    }
                    Node forkNode = ((Join) node).joinOf;
                    if (forkNode instanceof Fork && forkNode.getOutput() instanceof Forker) {
                        int forkIdx = getNodeIndex(forkNode);
                        if (forkIdx >= 0 && forkIdx < connq.getNumRows()) {
                            double w = Math.max(1.0,
                                    Math.round(((Forker) forkNode.getOutput()).tasksPerLink));
                            double outdeg = 0;
                            for (int j = 0; j < connq.getNumCols(); j++) {
                                if (connq.get(forkIdx, j) > 0) outdeg++;
                            }
                            nsib = outdeg * w;
                        }
                    }
                }
                for (int r0 = 0; r0 < K; r0++) {
                    JobClass jc = this.jobClasses.get(r0);
                    JoinStrategy js = joiner.joinStrategy.get(jc);
                    Double kreq = joiner.joinRequired.get(jc);
                    if (js == null || js == JoinStrategy.STD || kreq == null || kreq <= 0) continue;
                    if (nsib <= 0 || kreq < nsib) quorumClass[r0] = true;
                }
            }
        }
        boolean[] quorumFedChain = new boolean[C];
        for (int c0 = 0; c0 < C; c0++) {
            Matrix inchain_c0 = this.sn.inchain.get(c0);
            for (int idx = 0; idx < inchain_c0.length(); idx++) {
                int r0 = (int) inchain_c0.get(idx);
                if (r0 >= 0 && r0 < K && quorumClass[r0]) quorumFedChain[c0] = true;
            }
        }
        // A fork with tasksPerLink = w > 1 puts w tasks of the SAME parent on one link,
        // so a branch station can hold w jobs per circulating parent and the chain
        // population is no longer its bound. The multiplier is the PRODUCT over the
        // forks, because a fork nested in another's branch multiplies again; that is an
        // upper bound for forks in series, where a cap that never binds costs nothing,
        // and exact for the single-fork case. Without it a simulator drops a closed job
        // at a branch station. see _kb/04-networkstruct.md
        double forkTaskFactor = 1.0;
        for (Node node : this.nodes) {
            if (!(node instanceof Fork) || !(node.getOutput() instanceof Forker)) continue;
            double tpl = ((Forker) node.getOutput()).tasksPerLink;
            if (!Double.isNaN(tpl) && tpl >= 1) {
                forkTaskFactor *= Math.max(1.0, Math.round(tpl));
            }
        }

        for (int c = 0; c < C; c++) {
            Matrix inchain_c = this.sn.inchain.get(c);
            //chainCap = sum(njobs(inchain));
            double chainCap = 0;
            for (int idx = 0; idx < inchain_c.length(); idx++) {
                chainCap += njobs.get(0, (int) inchain_c.get(0, idx));
            }
            chainCap *= forkTaskFactor;

            if (chainCap >= Integer.MAX_VALUE || spawnFedChain[c] || quorumFedChain[c]) {
                chainCap = Inf;
            }

            for (int idx = 0; idx < inchain_c.length(); idx++) {
                int r = (int) inchain_c.get(idx);
                JobClass jobclass = this.jobClasses.get(r);
                for (int i = 0; i < M; i++) {
                    Station station = this.stations.get(i);
                    if (!(station instanceof Source)) {
                        DropStrategy stationDropRule = station.getDropRule(jobclass);
                        // see _kb/04-networkstruct.md (droprule/classcap derivation) for rationale
                        boolean isUserRuleNode = !(station instanceof Place) && !(station instanceof Join)
                                && !isRetrievalSystemQueue(station);
                        if (isUserRuleNode && stationDropRule == DropStrategy.WaitingQueue
                                && Double.isInfinite(njobs.get(0, r))) {
                            boolean stationCapIsFinite = station.hasFiniteCap();
                            // see _kb/04-networkstruct.md (droprule/classcap derivation) for rationale
                            boolean classCapIsFinite = station.getClassCap(jobclass) >= 1
                                    && !Double.isInfinite(station.getClassCap(jobclass))
                                    && station.getClassCap(jobclass) < Integer.MAX_VALUE;
                            if (stationCapIsFinite || classCapIsFinite) {
                                throw new RuntimeException("Station '" + station.getName()
                                        + "' declares setDropRule(WAITQ) for the open class '"
                                        + jobclass.getName() + "' at a finite capacity. LINE does not"
                                        + " implement waiting-room blocking for an open arrival at a"
                                        + " plain finite buffer: no solver honours this combination"
                                        + " (SolverCTMC drops the arrival, SolverJMT blocks at the"
                                        + " Source and ignores the capacity). Use"
                                        + " DropStrategy.Drop for a loss station (M/M/1/K), or one of"
                                        + " the blocking policies LINE implements (DropStrategy.BAS,"
                                        + " DropStrategy.BBS, DropStrategy.RSRD) for blocking between"
                                        + " stations.");
                            }
                        }
                        // see _kb/04-networkstruct.md (droprule/classcap derivation) for rationale
                        if (stationDropRule == null) {
                            if (Double.isInfinite(station.getCap()) || station.getCap() >= Integer.MAX_VALUE) {
                                // No buffer: the rule is never consulted.
                                stationDropRule = DropStrategy.WaitingQueue;
                            } else if (station.getCap() >= 0 && njobs.get(0, r) < Integer.MAX_VALUE) {
                                // see _kb/04-networkstruct.md (droprule/classcap derivation) for rationale
                                stationDropRule = DropStrategy.WaitingQueue;
                            } else {
                                // see _kb/04-networkstruct.md (droprule/classcap derivation) for rationale
                                stationDropRule = DropStrategy.Drop;
                            }
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
                        // see _kb/04-networkstruct.md (droprule/classcap derivation) for rationale
                        int orbitCap = station.getOrbitMaxJobs(jobclass);
                        if (orbitCap >= 0) {
                            int nsrv = station.getNumberOfServers();
                            if (nsrv < 1 || nsrv == Integer.MAX_VALUE) {
                                nsrv = 1;
                            }
                            classcap.set(i, r, FastMath.min(classcap.get(i, r), nsrv + orbitCap));
                        }
                    }
                }
            }
        }

        for (int i = 0; i < M; i++) {
            Station station = this.stations.get(i);
            // see _kb/04-networkstruct.md (droprule/classcap derivation) for rationale
            if (station.hasFiniteCap()) {
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
        this.sn.regionlincon = new HashMap<Integer, MatrixCell>();
        // see _kb/04-networkstruct.md (refreshRegions.m section) for rationale
        this.sn.regionmaxmem = new MatrixCell(F);
        // see _kb/04-networkstruct.md (refreshRegions.m section) for rationale
        this.sn.regionmembers = new MatrixCell(F);

        for (int f = 0; f < F; f++) {
            Region fcr = this.regions.get(f);
            // Matrix with M rows (stations) and K+1 columns (K classes + 1 global)
            Matrix regionMatrix = new Matrix(M, K + 1);
            regionMatrix.fill(-1);  // Initialize all to infinite (-1)
            Matrix regionMemMatrix = new Matrix(M, 1);
            regionMemMatrix.fill(-1);  // Initialize all to unbounded (-1)
            Matrix regionMemberMask = new Matrix(M, 1);
            regionMemberMask.fill(0);  // membership, independent of the caps

            // Find which stations are in this region and set their capacities
            for (Node node : fcr.getNodes()) {
                for (int i = 0; i < M; i++) {
                    if (this.stations.get(i) == node) {
                        regionMemberMask.set(i, 0, 1);
                        // see _kb/04-networkstruct.md (refreshRegions.m section) for rationale
                        for (int r = 0; r < K; r++) {
                            JobClass jc = this.jobClasses.get(r);
                            int classMax = fcr.getClassMaxJobs(jc);
                            int memMax = fcr.getClassMaxMemory(jc);
                            double szr = fcr.getClassSize(jc);
                            if (memMax != Region.UNBOUNDED && szr > 0) {
                                int memJobs = (int) Math.floor(memMax / szr);
                                classMax = (classMax == Region.UNBOUNDED) ? memJobs : Math.min(classMax, memJobs);
                            }
                            regionMatrix.set(i, r, classMax);
                        }
                        // Set global max jobs for this station in this region (column K)
                        int globalMax = fcr.getGlobalMaxJobs();
                        regionMatrix.set(i, K, globalMax);
                        // Replicate the region-global memory budget on this member row
                        regionMemMatrix.set(i, 0, fcr.getGlobalMaxMemory());
                        break;
                    }
                }
            }
            this.sn.regionmaxmem.set(f, regionMemMatrix);
            this.sn.regionmembers.set(f, regionMemberMask);

            // Extract per-class drop rules, weights, and sizes for this region
            for (int r = 0; r < K; r++) {
                JobClass jobClass = this.jobClasses.get(r);
                DropStrategy classDropStrategy = fcr.getDropStrategy(jobClass);
                this.sn.regionrule.set(f, r, classDropStrategy.getID());
                this.sn.regionweight.set(f, r, fcr.getClassWeight(jobClass));
                this.sn.regionsz.set(f, r, fcr.getClassSize(jobClass));
            }

            this.sn.region.set(f, regionMatrix);

            // Serialize linear constraints if set: store the (A,b) pair on the
            // same cell row of regionlincon (get(0) = A, get(1) = b)
            if (fcr.hasLinearConstraints()) {
                Matrix[] linCon = fcr.getLinearConstraints();
                this.sn.regionlincon.put(f, new MatrixCell(linCon[0], linCon[1]));
            }
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

        if (!hasUsableCsMatrix(K)) {
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
            // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
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
            jline.util.Pair<List<List<Matrix>>, Matrix> rtorigResult = snRtnodesToRtorig(this.sn);
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
        Map<Station, Map<JobClass, SerializableFunction<Complex, Complex>>> lst;

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
            lst = new HashMap<Station, Map<JobClass, SerializableFunction<Complex, Complex>>>();
            for (Station station : stations) {
                lst.put(station, new HashMap<JobClass, SerializableFunction<Complex, Complex>>());
            }
        }
        int sourceIdx = this.getIndexSourceStation();

        for (Integer i : statSet) {
            Station station = this.stations.get(i);
            Map<JobClass, SerializableFunction<Complex, Complex>> map = new HashMap<JobClass, SerializableFunction<Complex, Complex>>();
            for (Integer r : classSet) {
                JobClass jobclass = this.jobClasses.get(r);
                if (i == sourceIdx) {
                    Distribution distr = ((Source) station).getArrivalDistribution(jobclass);
                    if (distr instanceof Disabled) map.put(jobclass, null);
                    // see _kb/04-networkstruct.md (refreshStruct.m field-population notes) for rationale
                    else map.put(jobclass, e -> distr.evalLST(e));  // the Complex overload
                } else {
                    //line 45-46 is ignored since Fork is not station
                    if (station instanceof Join) map.put(jobclass, null);
                    else
                        map.put(jobclass, e -> station.getServer().getServiceDistribution(jobclass).evalLST(e));
                }
            }
            lst.put(station, map);
        }

        if (this.sn != null) this.sn.lst = lst;
    }

    /**
     * Decides whether station {@code ind} must carry a true-BAS blocked marker.
     *
     * A station needs the marker when it can be BLOCKED, and blocking is caused by
     * the DESTINATION being full, not by this station's own capacity. Testing this
     * station's own capacity (as this predicate once did) silently downgraded the
     * canonical BAS shape -- upstream unbounded, downstream finite -- to repetitive
     * service, and made a NON-binding capacity on the upstream change the answer by
     * 29%. So the test is reachability: is there a directly reachable destination
     * station with a finite capacity, for which BAS is declared?
     *
     * BAS must be declared HERE, on the station that does the blocking (the
     * convention of LINE's own examples/basic/closedQN/cqn_bas_blocking.m). The
     * marker column is shared with the polling controller and the cache width, so
     * FromMarginal disambiguates it by re-testing this station's own BAS rule;
     * declaration and enumeration must agree on that same test, which is why a
     * destination-declared BAS (the JMT convention, used by line-test
     * testsAdvFeatures/des/test_des_bas_closed.m) is NOT honoured here -- see
     * BUG-83. Note the complementary half of the mechanism -- the destination
     * refusing the arrival in AfterEventStation -- reads the DESTINATION's rule, so
     * the two halves only ever agreed when BAS was set on every station, which is
     * what masked the guard defect.
     *
     * @param ind      the node index of the candidate station
     * @param R        the number of job classes
     * @param destmask (nstations x nclasses) accumulator, marked at the (destination
     *                 station, class) pairs whose refusals must block {@code ind}
     *                 rather than be lost; may be null when only the flag is wanted
     * @return true when the station must reserve the blocked marker
     */
    private boolean declaresBlockedMarker(int ind, int R, Matrix destmask) {
        Node bnode = this.nodes.get(ind);
        if (!(bnode instanceof Station) || bnode instanceof Source || bnode instanceof Cache) {
            return false;
        }
        Station bst = (Station) bnode;
        boolean tf = false;
        List<Integer> dests = downstreamStations(ind);
        for (int r = 0; r < R; r++) {
            JobClass jobclass = this.jobClasses.get(r);
            boolean hereBAS = bst.getDropRule(jobclass) == DropStrategy.BlockingAfterService;
            for (int d = 0; d < dests.size(); d++) {
                int dind = dests.get(d).intValue();
                Node dnode = this.nodes.get(dind);
                if (!(dnode instanceof Station) || dnode instanceof Source) {
                    continue;
                }
                double dcap = ((Station) dnode).getCap();
                if (!(dcap > 0 && dcap < Integer.MAX_VALUE && !Double.isInfinite(dcap))) {
                    continue; // the destination must be able to fill
                }
                boolean thereBAS = ((Station) dnode).getDropRule(jobclass) == DropStrategy.BlockingAfterService;
                if (hereBAS || thereBAS) {
                    // see _kb/04-networkstruct.md (BAS blocking marker section) for rationale
                    tf = true;
                    if (destmask != null && this.sn != null && this.sn.nodeToStation != null) {
                        int jst = (int) this.sn.nodeToStation.get(dind);
                        if (jst >= 0 && jst < destmask.getNumRows()) {
                            destmask.set(jst, r, 1.0);
                        }
                    }
                }
            }
        }
        return tf;
    }

    /**
     * Node indices of the stations directly downstream of {@code ind}, walking through
     * intermediate stateless nodes (Router, ClassSwitch, ...) but stopping at the first
     * station on each path, since a blocked job is held for its immediate destination.
     *
     * @param ind the source node index
     * @return the node indices of the immediately downstream stations
     */
    private List<Integer> downstreamStations(int ind) {
        List<Integer> out = new ArrayList<Integer>();
        if (this.sn == null || this.sn.connmatrix == null || this.sn.connmatrix.isEmpty()) {
            return out;
        }
        int n = this.sn.connmatrix.getNumCols();
        boolean[] seen = new boolean[this.nodes.size()];
        Deque<Integer> queue = new ArrayDeque<Integer>();
        queue.add(ind);
        seen[ind] = true;
        while (!queue.isEmpty()) {
            int cur = queue.poll();
            for (int j = 0; j < n && j < this.nodes.size(); j++) {
                if (this.sn.connmatrix.get(cur, j) != 1 || seen[j]) {
                    continue;
                }
                seen[j] = true;
                Node jn = this.nodes.get(j);
                if (jn instanceof Station) {
                    out.add(j); // stop here: this is an immediate destination
                } else {
                    queue.add(j); // stateless hop, keep walking
                }
            }
        }
        return out;
    }

    public void refreshLocalVars() {
        int R = this.jobClasses.size();
        int I = this.nodes.size();
        // Columns 0..R-1 modulation phases, R..2R-1 routing vars, 2R the shared node
        // block (cache width / BAS marker / polling controller). Columns 2R+1+r are the
        // synchronous-call (REPLY) blocked-server counters, appended so that every
        // existing nvars reader keeps its indices; see ReplyBlock. They stay zero unless
        // the model declares a REPLY signal, so no other model changes state width.
        Matrix nvars = new Matrix(I, 3 * R + 1);
        // see _kb/04-networkstruct.md (BAS blocking marker section) for rationale
        Matrix isbasblocking = new Matrix(this.nodes.size(), 1);
        isbasblocking.zero();
        // see _kb/04-networkstruct.md (BAS blocking marker section) for rationale
        Matrix isbasdestination = new Matrix(this.stations.size(), R);
        isbasdestination.zero();
        Map<Node, NodeParam> nodeparam = new HashMap<Node, NodeParam>();
        // Polling stations, sized in a second pass below: Polling.info derives the
        // controller width from nodeparam, which is only complete once this loop ends.
        List<Integer> pollingParamNodes = new ArrayList<Integer>();
        List<QueueNodeParam> pollingParamValues = new ArrayList<QueueNodeParam>();

        for (int ind = 0; ind < I; ind++) {
            Node node = this.nodes.get(ind);
            NodeParam param = null;
            switch (this.sn.nodetype.get(ind)) {
                case Cache:
                    CacheNodeParam cacheParam = new CacheNodeParam();
                    Cache cache = (Cache) node;
                    cacheParam.nitems = 0;
                    cacheParam.accost = cache.accessProb;
                    // popularity is keyed (itemSetIndex, class): address this cache's own
                    // row explicitly, since the linear form aliases to it only when the
                    // cache's item set is the first in the model, i.e. only when the model
                    // holds a single cache.
                    int itemRow = cache.getItems().getIndex();
                    for (int r = 0; r < this.getNumberOfClasses(); r++) {
                        Distribution pop = cache.popularityGet(itemRow, r);
                        if (pop != null && !pop.isDisabled()) {
                            cacheParam.nitems = (int) Maths.max(cacheParam.nitems, pop.getSupport().getRight());
                        }
                    }
                    // see _kb/04-networkstruct.md (State-space construction conventions) for rationale
                    int retrievalBitmapWidth = (cache.getRetrievalSystemCapacity() > 0) ? cacheParam.nitems : 0;
                    int retrievalPendingWidth = (cache.getRetrievalSystemCapacity() > 0)
                            ? cache.getRetrievalClassIndices().size() : 0;
                    nvars.set(ind, 2 * R, cache.getTotalCacheCapacity() + retrievalBitmapWidth + retrievalPendingWidth);
                    cacheParam.itemcap = cache.getItemLevelCap();
                    // per-item storage costs and per-list cost caps (ton21cache Sec. IX)
                    cacheParam.itemsize = cache.getItemSizes();
                    cacheParam.costcap = cache.getCostCaps();
                    cacheParam.costcapglobal = cache.isCostCapGlobal();
                    cacheParam.totalCacheCapacity = cache.getTotalCacheCapacity();
                    cacheParam.retrievalSystemCapacity = cache.getRetrievalSystemCapacity();
                    cacheParam.pread = new HashMap<>();
                    for (int r = 0; r < this.getNumberOfClasses(); r++) {
                        Distribution pop = cache.popularityGet(itemRow, r);
                        if (pop == null || pop.isDisabled()) {
                            cacheParam.pread.put(r, null);
                        } else {
                            List<Double> t = new ArrayList<>();
                            for (int j = 1; j <= cacheParam.nitems; j++) {
                                t.add((double) j);
                            }
                            cacheParam.pread.put(r, ((DiscreteDistribution) pop).evalPMF(t).toList1D());
                        }
                    }
                    cacheParam.classitem = new HashMap<>();
                    for (int r = 0; r < this.getNumberOfClasses(); r++) {
                        cacheParam.classitem.put(r, cache.getItemOfClass(r));
                    }
                    cacheParam.replacestrat = cache.getReplacementStrategy();
                    cacheParam.qlru = cache.getAdmissionProb();
                    // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
                    if (cacheParam.replacestrat == ReplacementStrategy.CLIMB) {
                        int Cclimb = (int) cacheParam.itemcap.elementSum();
                        Matrix unitCap = new Matrix(1, Cclimb);
                        for (int lc = 0; lc < Cclimb; lc++) unitCap.set(0, lc, 1.0);
                        cacheParam.itemcap = unitCap;
                        cacheParam.replacestrat = ReplacementStrategy.FIFO;
                        Matrix chain = new Matrix(Cclimb + 1, Cclimb + 1);
                        chain.zero();
                        for (int rc = 0; rc < Cclimb; rc++) chain.set(rc, rc + 1, 1.0);
                        chain.set(Cclimb, Cclimb, 1.0);
                        int Kclimb = this.getNumberOfClasses();
                        Matrix[][] acClimb = new Matrix[Kclimb][cacheParam.nitems];
                        for (int vc = 0; vc < Kclimb; vc++) {
                            for (int kc = 0; kc < cacheParam.nitems; kc++) {
                                acClimb[vc][kc] = chain;
                            }
                        }
                        cacheParam.accost = acClimb;
                    }
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
                    cacheParam.retrievalClasses = new Matrix(cache.getCacheServer().retrievalClasses.getNumRows(),
                            cache.getCacheServer().retrievalClasses.getNumCols());
                    for (int r = 0; r < cache.getCacheServer().retrievalClasses.getNumRows(); r++) {
                        for (int c = 0; c < cache.getCacheServer().retrievalClasses.getNumCols(); c++) {
                            cacheParam.retrievalClasses.set(r, c, FastMath.round(cache.getCacheServer().retrievalClasses.get(r, c)));
                        }
                    }
                    cacheParam.retrievalClassIndices = new HashSet<>();
                    cacheParam.retrievalClassIndices.addAll(cache.getRetrievalClassIndices());
                    cacheParam.retrievalSystemQueueIndices = new HashMap<>();
                    for (Map.Entry<Integer, List<Integer>> entry : cache.getRetrievalSystemQueueIndices().entrySet()) {
                        Integer jobClass = entry.getKey();
                        List<Integer> queueIndices = entry.getValue();
                        List<Integer> newRetrievalSystemQueueIndices = new ArrayList<>();
                        for (Integer index : queueIndices) {
                            newRetrievalSystemQueueIndices.add(index);
                        }
                        cacheParam.retrievalSystemQueueIndices.put(jobClass, newRetrievalSystemQueueIndices);
                    }
                    param = cacheParam;
                    break;
                case Fork:
                    ForkNodeParam forkParam = new ForkNodeParam();
                    Forker forker = (Forker) node.getOutput();
                    forkParam.fanOut = forker.tasksPerLink;
                    buildForkFanout(forkParam, forker, node);
                    param = forkParam;
                    break;
                case Join:
                    JoinNodeParam joinParam = new JoinNodeParam();
                    Joiner joiner = (Joiner) node.getInput();
                    joinParam.joinStrategy = joiner.joinStrategy;
                    // fanIn is the JMT numRequired of a STANDARD join (-1 = every
                    // sibling), joinRequired the quorum k of a PARTIAL one: the two
                    // read the same field but are written to different JMT elements,
                    // so both are carried
                    joinParam.fanIn = joiner.joinRequired;
                    joinParam.joinRequired = joiner.joinRequired;
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
                            // see _kb/04-networkstruct.md (State-space construction conventions) for rationale
                            if (arrivalDistrib instanceof MAP || arrivalDistrib instanceof DMAP
                                    || arrivalDistrib instanceof MMPP2 || arrivalDistrib instanceof BMAP) {
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

                    // see _kb/04-networkstruct.md (refreshStruct.m field-population notes) for rationale
                    if (queue.isDelayOffEnabled()) {
                        for (int r = 0; r < R; r++) {
                            JobClass jobClass = this.getClassByIndex(r);
                            Distribution setup = queue.getSetupTime(jobClass);
                            Distribution delayoff = queue.getDelayOffTime(jobClass);
                            if (setup != null && delayoff != null) {
                                queueParam.setupTime.put(jobClass, setup);
                                queueParam.delayoffTime.put(jobClass, delayoff);
                            }
                        }
                    }

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
                        // see _kb/04-networkstruct.md (refreshStruct.m field-population notes) for rationale
                        pollingParamNodes.add(Integer.valueOf(ind));
                        pollingParamValues.add(queueParam);
                    }

                    for (int r = 0; r < R; r++) {
                        ServiceBinding serviceProcess = node.getServer().getServiceProcess(this.getClassByIndex(r));
                        if (serviceProcess != null) {
                            Distribution serviceDistrib = serviceProcess.getDistribution();

                            // Markov-modulated service needs the phase-restart state slot
                            // (must stay in sync with ismkvmodclass in State.afterEventInit)
                            if (serviceDistrib instanceof MAP || serviceDistrib instanceof DMAP
                                    || serviceDistrib instanceof MMPP2 || serviceDistrib instanceof BMAP) {
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

                    // see _kb/04-networkstruct.md (State-space construction conventions) for rationale
                    if (queue.getSchedStrategy() == SchedStrategy.PAS
                            || queue.getSchedStrategy() == SchedStrategy.OI) {
                        Matrix sg;
                        if (queue.getSchedStrategy() == SchedStrategy.OI) {
                            // Order-independent: swap graph is always zero (empty), so
                            // class order is preserved on completion (plain OI).
                            sg = new Matrix(R, R);
                        } else {
                            sg = queue.getSwapGraph();
                            if (sg == null) {
                                // PAS default: complete compatibility graph (no self-loops).
                                sg = new Matrix(R, R);
                                for (int a = 0; a < R; a++) {
                                    for (int b = 0; b < R; b++) {
                                        sg.set(a, b, a == b ? 0 : 1);
                                    }
                                }
                            }
                        }
                        queueParam.swapGraph = sg;
                        queueParam.svcRateFun = queue.getServiceRateFunction();
                    }

                    param = queueParam;
                    break;
                case Delay:
                    ServiceNodeParam delayParam = new ServiceNodeParam();

                    for (int r = 0; r < R; r++) {
                        ServiceBinding serviceProcess = node.getServer().getServiceProcess(this.getClassByIndex(r));
                        if (serviceProcess != null) {
                            Distribution serviceDistrib = serviceProcess.getDistribution();

                            // Markov-modulated service needs the phase-restart state slot
                            // (must stay in sync with ismkvmodclass in State.afterEventInit)
                            if (serviceDistrib instanceof MAP || serviceDistrib instanceof DMAP
                                    || serviceDistrib instanceof MMPP2 || serviceDistrib instanceof BMAP) {
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

                        if (firingDistrib instanceof MAP || firingDistrib instanceof DMAP) {
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
                    case SDR: {
                        // Krzesinski SDR: resolve the declared branch topology from
                        // node objects to node indices. The product-form solvers read
                        // the station-indexed twin from sn.sdr, built by
                        // refreshStateDepRouting.
                        jline.lang.StateDepRouting decl = node.getStateDepRouting(jobclass);
                        if (decl == null) {
                            throw new RuntimeException("Node " + node.getName() + " declares state-dependent routing "
                                    + "without a structure; use Node.setStateDepRouting.");
                        }
                        if (param.sdr == null) param.sdr = new HashMap<JobClass, jline.lang.StateDepRouting>();
                        int Bsdr = decl.branchNodes.size();
                        jline.lang.StateDepRouting resolved = new jline.lang.StateDepRouting();
                        resolved.entry = ind;
                        resolved.departure = decl.departureNode.getNodeIndex();
                        resolved.branch = new int[Bsdr][];
                        resolved.entryOf = new int[Bsdr];
                        resolved.departureOf = new int[Bsdr];
                        for (int b = 1; b < Bsdr; b++) {
                            List<Node> bn = decl.branchNodes.get(b);
                            int[] idx = new int[bn.size()];
                            for (int k = 0; k < bn.size(); k++) {
                                idx[k] = bn.get(k).getNodeIndex();
                            }
                            resolved.branch[b] = idx;
                            resolved.entryOf[b] = bn.get(0).getNodeIndex();
                            resolved.departureOf[b] = bn.get(bn.size() - 1).getNodeIndex();
                        }
                        resolved.level = decl.level.clone();
                        resolved.C = decl.C.clone();
                        resolved.d = new double[decl.d.length][];
                        for (int t = 0; t < decl.d.length; t++) {
                            resolved.d[t] = decl.d[t].clone();
                        }
                        resolved.entryNode = decl.entryNode;
                        resolved.departureNode = decl.departureNode;
                        resolved.branchNodes = decl.branchNodes;
                        param.sdr.put(jobclass, resolved);
                        break;
                    }
                    case SQ: {
                        if (param.d == null) param.d = new HashMap<JobClass, Integer>();
                        if (param.outlinks == null) param.outlinks = new HashMap<JobClass, Matrix>();

                        // Read d from the matching OutputStrategy entry.
                        int dVal = 2;
                        List<OutputStrategy> kchOs = node.getOutput().getOutputStrategyByClass(jobclass);
                        if (kchOs != null && !kchOs.isEmpty()) {
                            for (OutputStrategy os : kchOs) {
                                if (os.getRoutingStrategy() == jline.lang.constant.RoutingStrategy.SQ) {
                                    dVal = os.getSqD();
                                    break;
                                }
                            }
                        }
                        param.d.put(jobclass, dVal);

                        // outlinks for completeness (mirrors WRROBIN/RROBIN bookkeeping).
                        conn_i = new Matrix(0, 0);
                        Matrix.extractRows(this.sn.connmatrix, ind, ind + 1, conn_i);
                        conn_i_transpose = conn_i.find().transpose();
                        param.outlinks.put(jobclass, conn_i_transpose);
                        break;
                    }
                    case WRROBIN:
                        param.weights = new HashMap<JobClass, Matrix>();
                        param.outlinks = new HashMap<JobClass, Matrix>();
                        param.weightedOutlinks = new HashMap<JobClass, Matrix>();
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
                            // see _kb/04-networkstruct.md (State-space construction conventions) for rationale
                            if (destination == null) {
                                continue;
                            }
                            Double weight = outputStrategy_r.get(c).getProbability();
                            param.weights.get(jobclass).set(0, destination.getNodeIndex(), weight);
                        }
                        // see _kb/04-networkstruct.md (State-space construction conventions) for rationale
                        Matrix wrrOutlinks = param.outlinks.get(jobclass);
                        Matrix wrrWeights = param.weights.get(jobclass);
                        List<Integer> wrrCycle = new ArrayList<Integer>();
                        for (int oi = 0; oi < wrrOutlinks.length(); oi++) {
                            int destNode = (int) wrrOutlinks.get(oi);
                            int w = (int) Math.round(wrrWeights.get(0, destNode));
                            if (w < 1) {
                                w = 1;
                            }
                            for (int wc = 0; wc < w; wc++) {
                                wrrCycle.add(destNode);
                            }
                        }
                        Matrix wrrWol = new Matrix(wrrCycle.size(), 1);
                        for (int ci = 0; ci < wrrCycle.size(); ci++) {
                            wrrWol.set(ci, 0, wrrCycle.get(ci));
                        }
                        param.weightedOutlinks.put(jobclass, wrrWol);
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
            // see _kb/04-networkstruct.md (BAS blocking marker section) for rationale
            if (declaresBlockedMarker(ind, R, isbasdestination)) {
                // see _kb/04-networkstruct.md (BAS blocking marker section) for rationale
                if (((Station) node).getSchedStrategy() == SchedStrategy.POLLING) {
                    throw new RuntimeException("True BAS blocking is not supported at a polling station ("
                            + node.getName() + "): the polling controller and the BAS blocked marker share one "
                            + "local-state column. Use a non-polling scheduling strategy at the blocking station, "
                            + "or remove the BAS drop rule.");
                }
                nvars.set(ind, 2 * R, 1);
                isbasblocking.set(ind, 0, 1);
            }
            // see _kb/04-networkstruct.md (BAS blocking marker section) for rationale
            if (this.sn != null && this.sn.hasbreakdown != null && ind < this.sn.hasbreakdown.length()
                    && this.sn.hasbreakdown.get(ind) == 1) {
                if (isbasblocking.get(ind, 0) == 1) {
                    throw new RuntimeException("Station '" + node.getName() + "' combines server breakdowns with "
                            + "true-BAS blocking: the breakdown status and the BAS blocked marker share one "
                            + "local-state column. Remove the BAS drop rule or the breakdown.");
                }
                if (((Station) node).getSchedStrategy() == SchedStrategy.POLLING) {
                    throw new RuntimeException("Server breakdowns are not supported at a polling station ("
                            + node.getName() + "): the polling controller and the breakdown status share one "
                            + "local-state column.");
                }
                nvars.set(ind, 2 * R, 1);
            }
            nodeparam.put(node, param);
        }

        // Synchronous call (REPLY signal): a job of a class that expects a reply leaves
        // this station for the callee but KEEPS its server, which is released only when
        // the matching REPLY signal class arrives back here. Reserve one counter column
        // per (station, calling class) so the state can carry the held servers; the job
        // itself is at the callee and therefore absent from this station's marginal.
        //
        // The holding station is identified structurally, since a CTMC has no job
        // identity to key on as LDES does: it is a station that the REPLY class is
        // routed INTO. In the canonical shape Client -> Server -> (switch to Reply) ->
        // Client, that selects the client and NOT the server -- the server's departure
        // is the one that CREATES the reply, and LDES likewise does not block it
        // (Solver_ssj: !classSwitchedToReply). Stations that never receive the reply
        // class carry no counter and are untouched.
        Matrix replyblock = new Matrix(I, R);
        replyblock.zero();
        if (this.sn != null && this.sn.syncreply != null && !this.sn.syncreply.isEmpty()
                && this.sn.rtnodes != null && !this.sn.rtnodes.isEmpty()) {
            Matrix rtnodes = this.sn.rtnodes;
            for (int r = 0; r < R; r++) {
                int s = (int) this.sn.syncreply.get(r, 0); // stored 0-based
                if (s < 0 || s >= R) {
                    continue;
                }
                for (int ind = 0; ind < I; ind++) {
                    Node node = this.nodes.get(ind);
                    if (!(node instanceof Station) || node instanceof Source) {
                        continue;
                    }
                    // An INF station has a server for every job, so holding one is
                    // immaterial and needs no state.
                    SchedStrategy sched = ((Station) node).getSchedStrategy();
                    if (sched == SchedStrategy.INF) {
                        continue;
                    }
                    // Does the reply class ever arrive here?
                    boolean arrivesHere = false;
                    for (int i = 0; i < I && !arrivesHere; i++) {
                        for (int q = 0; q < R; q++) {
                            if (rtnodes.get(i * R + q, ind * R + s) > 0) {
                                arrivesHere = true;
                                break;
                            }
                        }
                    }
                    if (!arrivesHere) {
                        continue;
                    }
                    if (sched != SchedStrategy.FCFS) {
                        throw new RuntimeException("Synchronous calls (REPLY signals) are supported only at FCFS "
                                + "stations, but " + node.getName() + " uses " + SchedStrategy.toText(sched)
                                + ". A held server is encoded as a per-class counter, which is exact only where "
                                + "servers are interchangeable (FCFS) or unlimited (INF); the other disciplines "
                                + "are not yet encoded rather than infeasible. Set this station to FCFS or INF, "
                                + "or simulate the layered model directly with SolverLDES.");
                    }
                    nvars.set(ind, 2 * R + 1 + r, 1);
                    replyblock.set(ind, r, 1);
                }
            }
        }

        if (this.sn != null) {
            this.sn.nvars = nvars;
            this.sn.nodeparam = nodeparam;
            refreshStateDepRouting();
            // Size the polling controller now that nodeparam is complete: Polling.info
            // reads sn.nodeparam/sn.sched/sn.proc, and memoizes itself on the param.
            for (int pi = 0; pi < pollingParamNodes.size(); pi++) {
                int ind = pollingParamNodes.get(pi).intValue();
                pollingParamValues.get(pi).pollinfo = null; // rebuild against the fresh struct
                jline.lang.state.Polling.Info pinfo = jline.lang.state.Polling.info(this.sn, ind);
                if (pinfo != null) {
                    nvars.set(ind, 2 * R, pinfo.width);
                }
            }
            this.sn.nvars = nvars;
            this.sn.isbasblocking = isbasblocking;
            this.sn.isbasdestination = isbasdestination;
            this.sn.replyblock = replyblock;
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
                transitionParam.firingdep = new ArrayList<>();
                int nnodes = this.getNumberOfNodes();
                int nclasses = this.getNumberOfClasses();
                for (Mode m : transition.getModes()) {
                    // Pad enabling/inhibiting/firing to (nnodes x nclasses) since
                    // addMode initializes with node count at creation time
                    Matrix en = transition.enablingConditions.get(m);
                    if (en.getNumRows() < nnodes || en.getNumCols() < nclasses) {
                        Matrix enFull = new Matrix(nnodes, nclasses);
                        enFull.zero();
                        for (int pi = 0; pi < en.getNumRows(); pi++) {
                            for (int pj = 0; pj < en.getNumCols(); pj++) {
                                enFull.set(pi, pj, en.get(pi, pj));
                            }
                        }
                        en = enFull;
                    }
                    Matrix inh = transition.inhibitingConditions.get(m);
                    if (inh.getNumRows() < nnodes || inh.getNumCols() < nclasses) {
                        Matrix inhFull = new Matrix(nnodes, nclasses);
                        inhFull.fill(Double.POSITIVE_INFINITY);
                        for (int pi = 0; pi < inh.getNumRows(); pi++) {
                            for (int pj = 0; pj < inh.getNumCols(); pj++) {
                                inhFull.set(pi, pj, inh.get(pi, pj));
                            }
                        }
                        inh = inhFull;
                    }
                    Matrix fir = transition.firingOutcomes.get(m);
                    if (fir.getNumRows() < nnodes || fir.getNumCols() < nclasses) {
                        Matrix firFull = new Matrix(nnodes, nclasses);
                        firFull.zero();
                        for (int pi = 0; pi < fir.getNumRows(); pi++) {
                            for (int pj = 0; pj < fir.getNumCols(); pj++) {
                                firFull.set(pi, pj, fir.get(pi, pj));
                            }
                        }
                        fir = firFull;
                    }
                    transitionParam.enabling.add(en);
                    transitionParam.inhibiting.add(inh);
                    transitionParam.firing.add(fir);
                    transitionParam.timing.add(transition.timingStrategies.get(m));
                    // Marking-dependent firing-rate multiplier (null == unit).
                    transitionParam.firingdep.add(transition.getFiringRateDependence(m));
                }
                transitionParam.nmodeservers = transition.getNumberOfModeServers();
                transitionParam.firingprio = transition.firingPriorities;
                transitionParam.fireweight = transition.firingWeights;

                transitionParam.firingproc = new LinkedHashMap<>();
                transitionParam.firingpie = new LinkedHashMap<>();
                transitionParam.firingphases = new Matrix(1, transition.getNumberOfModes());

                for (Mode m : transition.getModes()) {
                    if (transition.getFiringDistribution(m) instanceof Markovian) {
                        transitionParam.firingproc.put(m, ((Markovian) transition.getFiringDistribution(m)).getProcess());
                        transitionParam.firingpie.put(m, ((Markovian) transition.getFiringDistribution(m)).getInitProb());
                        transitionParam.firingphases.set(0, m.getIndex() - 1, (int) ((Markovian) transition.getFiringDistribution(m)).getNumberOfPhases());
                    } else if (transition.getFiringDistribution(m) instanceof Pareto) {
                        // Pareto has its own getProcess() method returning {shape, scale}
                        transitionParam.firingproc.put(m, ((Pareto) transition.getFiringDistribution(m)).getProcess());
                        transitionParam.firingpie.put(m, Matrix.singleton(NaN));
                        transitionParam.firingphases.set(m.getIndex() - 1, NaN);
                    } else {
                        // see _kb/04-networkstruct.md (refreshStruct.m field-population notes) for rationale
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
                        transitionParam.firingprocid = new LinkedHashMap<>();
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
            applyMarkedPhaseAccounting();
            Matrix.concatColumns(new Matrix(this.sn.phases.getNumRows(), 1), this.sn.phasessz.cumsumViaRow(), this.sn.phaseshift);
        }
        refreshProcessRepresentations();
    }

    /**
     * Marked (MMAP) source classes share the carrier's modulating chain: the
     * non-carrier classes (mark index > 1) contribute a single always-zero
     * state column rather than their own phase block (mirrors MATLAB
     * refreshProcessRepresentations).
     */
    private void applyMarkedPhaseAccounting() {
        if (this.sn == null || this.sn.markidx == null || this.sn.phasessz == null) {
            return;
        }
        int M = Math.min(this.sn.markidx.getNumRows(), this.sn.phasessz.getNumRows());
        int K = Math.min(this.sn.markidx.getNumCols(), this.sn.phasessz.getNumCols());
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < K; j++) {
                if (this.sn.markidx.get(i, j) > 1) {
                    this.sn.phasessz.set(i, j, 1.0);
                }
            }
        }
    }

    /**
     * True when the station's process representation for this class holds raw
     * distribution parameters rather than a (D0,D1) pair. Gamma, Weibull,
     * Lognormal, Pareto and Uniform return two 1x1 parameters from getProcess,
     * which is shape-indistinguishable from a genuine single-phase MAP, so the
     * process type has to be asked directly. Mirrors isRawParam_ir in
     * refreshProcessRepresentations.m.
     */
    private boolean isRawParameterProcess(Station station, boolean isSource, JobClass jobclass) {
        Distribution dist = null;
        if (isSource && station instanceof Source) {
            dist = ((Source) station).getArrivalDistribution(jobclass);
        } else if (station instanceof ServiceStation) {
            dist = ((ServiceStation) station).getServiceProcess(jobclass);
        }
        return dist instanceof jline.lang.processes.Gamma
                || dist instanceof jline.lang.processes.Weibull
                || dist instanceof jline.lang.processes.Lognormal
                || dist instanceof jline.lang.processes.Pareto
                || dist instanceof jline.lang.processes.Uniform;
    }

    @SuppressWarnings("unchecked")
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
            ph.put(station, ph_i);

            for (int r = 0; r < K; r++) {
                MatrixCell ph_i_r = ph_i.get(this.jobClasses.get(r));
                if (ph_i_r == null) {
                    phases.set(i, r, 1.0);
                } else if (ph_i_r.get(0) == null || ph_i_r.get(0).hasNaN()) {
                    // Disabled or non-Markovian with NaN parameters
                    phases.set(i, r, 0.0);
                } else if (ph_i_r.get(1) == null || isRawParameterProcess(station, i == sourceIdx, this.jobClasses.get(r))) {
                    // see _kb/04-networkstruct.md (refreshStruct.m field-population notes) for rationale
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
                    // NHPP carries a rate schedule, not D0/D1: always 1 active phase
                    JobClass jc4phases = this.jobClasses.get(r);
                    Distribution d4phases = (i == sourceIdx && station instanceof Source)
                            ? ((Source) station).getArrivalDistribution(jc4phases)
                            : (station instanceof ServiceStation
                                ? ((ServiceStation) station).getServiceProcess(jc4phases) : null);
                    if (d4phases instanceof MAPt) {
                        // a matrix schedule modulates a phase structure, so the phase
                        // count is the order of the segment matrices, not 1
                        phases.set(i, r, ((MAPt) d4phases).getNumberOfPhases());
                    } else if (d4phases instanceof PHt) {
                        phases.set(i, r, ((PHt) d4phases).getNumberOfPhases());
                    } else if (d4phases instanceof NHPP) {
                        phases.set(i, r, 1.0);
                    } else {
                        phases.set(i, r, ph_i_r.get(0).getNumCols());
                    }
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
                    if (map_ir != null && map_ir.get(0) != null && map_ir.get(1) != null
                            && !map_ir.get(0).hasNaN() && !map_ir.get(1).hasNaN()) {
                        // NHPP carries a rate schedule, not D0/D1: skip map_pie
                        Distribution d4pie = (i == sourceIdx && station instanceof Source)
                                ? ((Source) station).getArrivalDistribution(jobclass)
                                : (station instanceof ServiceStation
                                    ? ((ServiceStation) station).getServiceProcess(jobclass) : null);
                        if (d4pie instanceof MAPt || d4pie instanceof PHt) {
                            // pie of the time-averaged nominal, which the fluid carrier
                            // uses; NaN would strand the phase structure
                            MatrixCell nominal = (d4pie instanceof MAPt)
                                    ? ((MAPt) d4pie).getTimeAverageProcess()
                                    : ((PHt) d4pie).getTimeAverageProcessMAP();
                            pie_i.put(jobclass, map_pie(nominal.get(0), nominal.get(1)));
                        } else if (d4pie instanceof NHPP) {
                            Matrix nanPie = new Matrix(1, 1, 0);
                            nanPie.set(0, 0, NaN);
                            pie_i.put(jobclass, nanPie);
                        } else {
                            // PH/MAP representation has both D0 and D1
                            pie_i.put(jobclass, map_pie(map_ir.get(0), map_ir.get(1)));
                        }
                    } else {
                        // Non-Markovian distributions or disabled - use NaN pie
                        Matrix tmp = new Matrix(1, 1, 0);
                        tmp.set(0, 0, NaN);
                        pie_i.put(jobclass, tmp);
                    }
                }
                pie.put(station, pie_i);
            }

            // Record, per station-class, whether the representation admits a phase-type
            // reading. Consumers of mu/phi/pie (CTMC state space, SSA, fluid ODEs) treat
            // those as rates and probabilities, which only holds when the pair is
            // Markovian; a matrix-exponential process gives signed per-phase values.
            Map<Station, Map<JobClass, Boolean>> isph = new HashMap<>();
            for (int i = 0; i < M; i++) {
                Station station = this.stations.get(i);
                Map<JobClass, Boolean> isph_i = new HashMap<>();
                for (int r = 0; r < K; r++) {
                    JobClass jobclass = this.jobClasses.get(r);
                    isph_i.put(jobclass, jline.api.sn.SnIsPhaseType.snIsPhaseType(
                            ph.get(station).get(jobclass), pie.get(station).get(jobclass)));
                }
                isph.put(station, isph_i);
            }
            this.sn.isph = isph;

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
            applyMarkedPhaseAccounting();
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
        Map<Station, Map<JobClass, ImpatienceType>> impatienceClass = new HashMap<>();

        for (int i = 0; i < M; i++) {
            Station station = this.stations.get(i);
            impatienceType.put(station, new HashMap<>());
            impatienceMu.put(station, new HashMap<>());
            impatiencePhi.put(station, new HashMap<>());
            impatienceProc.put(station, new HashMap<>());
            impatiencePie.put(station, new HashMap<>());
            impatiencePhases.put(station, new HashMap<>());
            impatienceClass.put(station, new HashMap<>());

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

                // see _kb/04-networkstruct.md (refreshStruct.m field-population notes) for rationale
                ImpatienceType impClass = station.getImpatienceType(jobclass);
                if (impClass != null) {
                    impatienceClass.get(station).put(jobclass, impClass);
                }
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
            this.sn.impatienceClass = impatienceClass;
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
        // see _kb/04-networkstruct.md (refreshStruct.m field-population notes) for rationale
        Map<Station, Map<JobClass, Integer>> retrialPolicy = new HashMap<>();
        Map<Station, Map<JobClass, Integer>> orbitMaxJobs = new HashMap<>();

        for (int i = 0; i < M; i++) {
            Station station = this.stations.get(i);
            retrialType.put(station, new HashMap<>());
            retrialMu.put(station, new HashMap<>());
            retrialPhi.put(station, new HashMap<>());
            retrialProc.put(station, new HashMap<>());
            retrialMaxAttempts.put(station, new HashMap<>());
            retrialPolicy.put(station, new HashMap<>());
            orbitMaxJobs.put(station, new HashMap<>());

            for (int r = 0; r < K; r++) {
                JobClass jobclass = this.jobClasses.get(r);

                retrialPolicy.get(station).put(jobclass, station.getRetrialPolicy(jobclass));
                orbitMaxJobs.get(station).put(jobclass, station.getOrbitMaxJobs(jobclass));

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
            this.sn.retrialPolicy = retrialPolicy;
            this.sn.orbitMaxJobs = orbitMaxJobs;
        }
    }

    /**
     * Refreshes the server breakdown / repair configuration in the network structure.
     *
     * <p>Failure and repair are properties of the SERVER, so they are stored per
     * station; the optional degraded service used while the server is down is a
     * service distribution and is therefore per class. Mirrors the sn.hasbreakdown /
     * breakdownMu / repairMu / breakdownProc / repairProc / downServiceRates block
     * built by MATLAB refreshStruct.m.</p>
     */
    public void refreshBreakdown() {
        int M = this.stations.size();
        int K = this.jobClasses.size();
        int I = this.nodes.size();

        Matrix hasbreakdown = new Matrix(I, 1);
        hasbreakdown.zero();
        Matrix breakdownMu = new Matrix(M, 1);
        breakdownMu.zero();
        Matrix repairMu = new Matrix(M, 1);
        repairMu.zero();
        Matrix downServiceRates = new Matrix(M, K);
        downServiceRates.zero();
        Map<Station, MatrixCell> breakdownProc = new HashMap<Station, MatrixCell>();
        Map<Station, MatrixCell> repairProc = new HashMap<Station, MatrixCell>();

        for (int ist = 0; ist < M; ist++) {
            Station station = this.stations.get(ist);
            if (!(station instanceof Queue)) {
                continue;
            }
            Queue queue = (Queue) station;
            if (!queue.hasBreakdown()) {
                continue;
            }
            int ind = this.nodes.indexOf(station);
            if (ind >= 0) {
                hasbreakdown.set(ind, 0, 1);
            }
            breakdownMu.set(ist, 0, 1.0 / queue.getBreakdownFailure().getMean());
            repairMu.set(ist, 0, 1.0 / queue.getBreakdownRepair().getMean());
            if (queue.getBreakdownFailure() instanceof jline.lang.processes.Markovian) {
                breakdownProc.put(station,
                        ((jline.lang.processes.Markovian) queue.getBreakdownFailure()).getProcess());
            }
            if (queue.getBreakdownRepair() instanceof jline.lang.processes.Markovian) {
                repairProc.put(station,
                        ((jline.lang.processes.Markovian) queue.getBreakdownRepair()).getProcess());
            }
            for (int r = 0; r < K; r++) {
                JobClass jobclass = this.jobClasses.get(r);
                Distribution dsvc = queue.getDownService(jobclass);
                if (dsvc == null || dsvc.isDisabled()) {
                    continue;
                }
                if (!(dsvc instanceof Exp)) {
                    throw new RuntimeException("Station '" + station.getName() + "': the down-server service "
                            + "distribution must be exponential; a phase-type degraded service would need its own "
                            + "phase block in the joint chain.");
                }
                double dmean = dsvc.getMean();
                if (dmean > 0) {
                    downServiceRates.set(ist, r, 1.0 / dmean);
                }
            }
        }

        if (this.sn != null) {
            this.sn.hasbreakdown = hasbreakdown;
            this.sn.breakdownMu = breakdownMu;
            this.sn.repairMu = repairMu;
            this.sn.breakdownProc = breakdownProc;
            this.sn.repairProc = repairProc;
            this.sn.downServiceRates = downServiceRates;
        }
    }

    /**
     * Refreshes orbit impatience configuration in the network structure.
     * Extracts the (D0,D1) process representation of the orbit abandonment
     * distribution from all station-class pairs into sn.orbitImpatience.
     */
    public void refreshOrbitImpatience() {
        int M = this.stations.size();
        int K = this.jobClasses.size();

        Map<Station, Map<JobClass, MatrixCell>> orbitImpatience = new HashMap<>();

        for (int i = 0; i < M; i++) {
            Station station = this.stations.get(i);
            orbitImpatience.put(station, new HashMap<>());

            for (int r = 0; r < K; r++) {
                JobClass jobclass = this.jobClasses.get(r);

                Distribution orbitDist = station.getOrbitImpatience(jobclass);

                if (orbitDist != null && !orbitDist.isDisabled()
                        && orbitDist instanceof jline.lang.processes.Markovian) {
                    MatrixCell proc_val = ((jline.lang.processes.Markovian) orbitDist).getProcess();
                    orbitImpatience.get(station).put(jobclass, proc_val);
                }
            }
        }

        // Store in network structure
        if (this.sn != null) {
            this.sn.orbitImpatience = orbitImpatience;
        }
    }

    /**
     * Refreshes the batch rejection probabilities in the network structure.
     * Exports the per station-class batch reject probability configured on the
     * stations into sn.batchRejectProb. Mirrors the (nstations x nclasses)
     * sn.batchRejectProb matrix built by MATLAB refreshStruct.m and by the
     * native Python _refresh_balking_retrial(); entries default to 0, meaning
     * that partial admission of an arriving batch is allowed.
     */
    public void refreshBatchRejectProb() {
        int M = this.stations.size();
        int K = this.jobClasses.size();

        Map<Station, Map<JobClass, Double>> batchRejectProb = new HashMap<Station, Map<JobClass, Double>>();

        for (int i = 0; i < M; i++) {
            Station station = this.stations.get(i);
            Map<JobClass, Double> stationProb = new HashMap<JobClass, Double>();
            batchRejectProb.put(station, stationProb);

            for (int r = 0; r < K; r++) {
                JobClass jobclass = this.jobClasses.get(r);
                stationProb.put(jobclass, station.getBatchRejectProbability(jobclass));
            }
        }

        // Store in network structure
        if (this.sn != null) {
            this.sn.batchRejectProb = batchRejectProb;
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

        for (int i = 0; i < M; i++) {
            Station station = this.stations.get(i);

            // Heterogeneous-server parameters live on the station's ServiceNodeParam
            jline.lang.nodeparam.ServiceNodeParam snp = sn.getServiceParam(station);

            // Only Queue stations can have heterogeneous servers
            if (!(station instanceof Queue)) {
                if (snp != null) {
                    snp.nservertypes = 0;
                }
                continue;
            }

            Queue queue = (Queue) station;

            // Job parallelism is independent of the pools: a homogeneous station
            // may declare it, and a heterogeneous one may not.
            if (snp != null) {
                if (queue.hasServerParallelism()) {
                    Matrix par = new Matrix(1, K, K);
                    for (int r = 0; r < K; r++) {
                        par.set(0, r, queue.getServerParallelism(this.jobClasses.get(r)));
                    }
                    snp.serverparallelism = par;
                } else {
                    snp.serverparallelism = null;
                }
            }

            if (!queue.isHeterogeneous()) {
                if (snp != null) {
                    snp.nservertypes = 0;
                }
                continue;
            }

            // A heterogeneous Queue must carry a ServiceNodeParam (QueueNodeParam)
            if (snp == null) {
                continue;
            }

            List<ServerType> serverTypes = queue.getServerTypes();
            int nTypes = serverTypes.size();
            snp.nservertypes = nTypes;

            // Populate server type names
            List<String> names = new ArrayList<String>();
            for (ServerType st : serverTypes) {
                names.add(st.getName());
            }
            snp.servertypenames = names;

            // Populate servers per type: Matrix (nTypes x 1)
            Matrix spt = new Matrix(nTypes, 1, nTypes);
            for (int t = 0; t < nTypes; t++) {
                spt.set(t, 0, serverTypes.get(t).getNumOfServers());
            }
            snp.serverspertype = spt;

            // Populate compatibility matrix: Matrix (nTypes x K)
            Matrix compat = new Matrix(nTypes, K, nTypes * K);
            for (int t = 0; t < nTypes; t++) {
                ServerType st = serverTypes.get(t);
                for (int r = 0; r < K; r++) {
                    JobClass jc = this.jobClasses.get(r);
                    compat.set(t, r, st.isCompatible(jc) ? 1.0 : 0.0);
                }
            }
            snp.servercompat = compat;

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
            snp.heterorates = stationRates;
            snp.heteroproc = stationProc;
            snp.heteroprocid = stationProcId;

            // Scheduling policy
            snp.heteroschedpolicy = queue.getHeteroSchedPolicy();
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
                        if ((station instanceof Place && !((Place) station).isQueueing())
                                || !station.getServer().containsJobClass(jobclass)) {
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

        // Marked (MMAP) source arrivals: mark index (1-based) of class r at
        // source station i; -1 = not a marked class (mirrors MATLAB sn.markidx)
        Matrix markidx = new Matrix(M, K);
        markidx.fill(-1);
        for (int i = 0; i < M; i++) {
            Station station = this.stations.get(i);
            if (station instanceof Source && ((Source) station).getMarkedClasses() != null) {
                List<JobClass> markedClasses = ((Source) station).getMarkedClasses();
                for (int k = 0; k < markedClasses.size(); k++) {
                    int r = this.getJobClassIndex(markedClasses.get(k));
                    if (r >= 0) {
                        markidx.set(i, r, k + 1);
                    }
                }
            }
        }

        if (this.sn != null) {
            this.sn.procid = proctype;
            this.sn.markidx = markidx;
        }
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
                    Place place = (Place) station;
                    Distribution svc = place.isQueueing() ? place.getService(this.jobClasses.get(r)) : null;
                    if (svc != null && !(svc instanceof Disabled)) {
                        rates.set(i, r, svc.getRate());
                        scv.set(i, r, svc.getSCV());
                    } else {
                        rates.set(i, r, NaN);
                        scv.set(i, r, NaN);
                    }
                } else if (station.getServer() instanceof ServiceTunnel) {
                    if (station instanceof Source) {
                        if (!((Source) station).containsJobClass(this.jobClasses.get(r))) {
                            rates.set(i, r, NaN);
                            scv.set(i, r, NaN);
                        } else {
                            Source source = (Source) station;
                            Distribution distr = source.getArrivalDistribution(this.jobClasses.get(r));
                            int markOfClass = -1;
                            if (distr instanceof jline.lang.processes.MarkedMAP && source.getMarkedClasses() != null) {
                                markOfClass = source.getMarkedClasses().indexOf(this.jobClasses.get(r));
                            }
                            if (markOfClass >= 0) {
                                // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
                                jline.lang.processes.MAP marginal =
                                        ((jline.lang.processes.MarkedMAP) distr).toMAPs(markOfClass + 1);
                                rates.set(i, r, marginal.getRate());
                                scv.set(i, r, marginal.getSCV());
                            } else {
                                rates.set(i, r, distr.getRate());
                                scv.set(i, r, distr.getSCV());
                            }
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
                throw new RuntimeException("Routing strategy is unspecified at all nodes for class " + jobclass.getName() + " in model " + this.getName());
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
                                        break;
                                    case WRROBIN:
                                        map.put(jnd * K + s, (pair) -> sub_rr_wrr(ind_final, jnd_final, r_final, s_final, linksmat, pair.getLeft(), pair.getRight()));
                                        break;
                                    case JSQ:
                                        map.put(jnd * K + s, (pair) -> sub_jsq(ind_final, jnd_final, r_final, s_final, linksmat, pair.getLeft(), pair.getRight()));
                                        break;
                                    case SQ:
                                        map.put(jnd * K + s, (pair) -> sub_sq(ind_final, jnd_final, r_final, s_final, linksmat, pair.getLeft(), pair.getRight()));
                                        break;
                                    case SDR:
                                        map.put(jnd * K + s, (pair) -> sub_sdr(ind_final, jnd_final, r_final, s_final, linksmat, pair.getLeft(), pair.getRight()));
                                        break;
                                    default:
                                        map.put(jnd * K + s, (pair) -> rtnodes.get(ind_final * K + r_final, jnd_final * K + s_final));
                                        break;
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

        // stochastic complementation over non-stateful nodes: see _kb/04-networkstruct.md
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
        Map<Station, SerializableFunction<Matrix, Matrix>> cdscaling;
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
        Map<Station, Matrix> cdscalingpeak = getLimitedClassDependencePeak();
        Map<Station, SerializableFunction<Matrix, Matrix>> jdscaling = getLimitedJointDependence();
        Map<Station, Matrix> jdscalingpeak = getLimitedJointDependencePeak();
        SerializableFunction<Matrix, Matrix> gdscaling = getGlobalDependence();
        Matrix gdscalingpeak = getGlobalDependencePeak();
        int gdscalingcutoff = getGlobalDependenceCutoff();

        if (sn == null) sn = new NetworkStruct();

        jline.io.LineConsole.compiling(this.getName());
        jline.io.LineConsole.compileDetail("reading the routing strategies of %s and %s",
                jline.io.LineConsole.plural(this.nodes.size(), "node", "nodes"),
                jline.io.LineConsole.plural(this.jobClasses.size(), "class", "classes"));

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
        // signaltarget(c) = 0-based index of the (positive) class whose jobs signal
        // class c removes; -1 for non-signals or signals with no target.
        sn.signaltarget = new Matrix(sn.nclasses, 1, sn.nclasses);
        for (int c = 0; c < sn.nclasses; c++) sn.signaltarget.set(c, 0, -1);
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
                int tgt = signalClass.getTargetJobClassIndex();
                if (tgt >= 1) sn.signaltarget.set(c, 0, tgt - 1);
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

        // Initialize classspawn - class injected at the same station on each
        // service completion of class c (-1 if none), 1-based to 0-based
        sn.classspawn = new Matrix(sn.nclasses, 1, sn.nclasses);
        for (int c = 0; c < sn.nclasses; c++) {
            int spawnIdx = this.jobClasses.get(c).getSpawnClassIndex();
            sn.classspawn.set(c, 0, spawnIdx >= 1 ? spawnIdx - 1 : -1);
        }

        // see _kb/04-networkstruct.md (refreshStruct.m field-population notes) for rationale
        sn.arrivalbatch = new ArrayList<DiscreteDistribution>(sn.nclasses);
        for (int c = 0; c < sn.nclasses; c++) {
            sn.arrivalbatch.add(null);
        }
        int batchSourceIdx = this.getIndexSourceNode();
        if (batchSourceIdx >= 0 && this.nodes.get(batchSourceIdx) instanceof Source) {
            Source batchSource = (Source) this.nodes.get(batchSourceIdx);
            for (int c = 0; c < sn.nclasses; c++) {
                sn.arrivalbatch.set(c, batchSource.getArrivalBatch(this.jobClasses.get(c)));
            }
        }

        // Initialize signal removal configuration fields
        sn.signalremdist = new ArrayList<DiscreteDistribution>(sn.nclasses);
        sn.signalrempolicy = new ArrayList<RemovalPolicy>(sn.nclasses);
        sn.iscatastrophe = new Matrix(sn.nclasses, 1, sn.nclasses);
        for (int c = 0; c < sn.nclasses; c++) {
            sn.signalremdist.add(null);
            sn.signalrempolicy.add(null);
        }
        // Populate removal configuration from Signal, OpenSignal, and ClosedSignal classes
        for (int c = 0; c < sn.nclasses; c++) {
            JobClass jc = this.jobClasses.get(c);
            if (jc instanceof Signal) {
                Signal signalClass = (Signal) jc;
                if (signalClass.isCatastrophe()) {
                    sn.iscatastrophe.set(c, 0, 1.0);
                }
                sn.signalremdist.set(c, signalClass.getRemovalDistribution());
                sn.signalrempolicy.set(c, signalClass.getRemovalPolicy());
            } else if (jc instanceof OpenSignal) {
                OpenSignal signalClass = (OpenSignal) jc;
                if (signalClass.isCatastrophe()) {
                    sn.iscatastrophe.set(c, 0, 1.0);
                }
                sn.signalremdist.set(c, signalClass.getRemovalDistribution());
                sn.signalrempolicy.set(c, signalClass.getRemovalPolicy());
            } else if (jc instanceof ClosedSignal) {
                ClosedSignal signalClass = (ClosedSignal) jc;
                if (signalClass.isCatastrophe()) {
                    sn.iscatastrophe.set(c, 0, 1.0);
                }
                sn.signalremdist.set(c, signalClass.getRemovalDistribution());
                sn.signalrempolicy.set(c, signalClass.getRemovalPolicy());
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
        sn.cdscalingpeak = cdscalingpeak;
        sn.jdscaling = jdscaling;
        sn.jdscalingpeak = jdscalingpeak;
        sn.gdscaling = gdscaling;
        sn.gdscalingpeak = gdscalingpeak;
        sn.gdscalingcutoff = gdscalingcutoff;
        sn.nodetype = nodetypes;
        sn.isstateful = new Matrix(nodes.size(), 1, nodes.size());
        for (int i = 0; i < sn.nnodes; i++) {
            Node node = this.nodes.get(i);
            sn.isstateful.set(i, 0, 0.0); //state dependent service
            if (node instanceof Source || node instanceof Delay || node instanceof Queue || node instanceof Cache || node instanceof Join || node instanceof Router || node instanceof Place || node instanceof Transition || node instanceof StatefulFork) {
                sn.isstateful.set(i, 0, 1.0); //state dependent service
            }
        }
        // Station-indexed mask of queue stations carrying setup/delay-off
        // times, as in MATLAB refreshStruct.m.
        sn.hassetup = new Matrix(sn.nstations, 1, sn.nstations);
        for (int i = 0; i < sn.nstations; i++) {
            Station station = this.stations.get(i);
            boolean hasSetup = station instanceof Queue && ((Queue) station).isDelayOffEnabled();
            sn.hassetup.set(i, 0, hasSetup ? 1.0 : 0.0);
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
                    case SQ:
                    case SDR:
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

        jline.io.LineConsole.compileDetail("refreshing class priorities");
        refreshPriorities();
        refreshDeadlines();
        refreshProcesses(null, null);
        refreshImpatience();
        refreshBalking();
        refreshRetrial();
        refreshBreakdown();
        refreshOrbitImpatience();
        refreshBatchRejectProb();

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

        jline.io.LineConsole.compileDetail("computing the routing table and the chains");
        refreshChains(!sn.nodetype.contains(NodeType.Cache));
        jline.io.LineConsole.compileDetail("found %s over %s",
                jline.io.LineConsole.plural(sn.nchains, "chain", "chains"),
                jline.io.LineConsole.plural(sn.nclasses, "class", "classes"));

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

        jline.io.LineConsole.compileDetail("refreshing node parameters and state-dependent routing");
        refreshLocalVars();
        // Must run after refreshLocalVars(), which populates the per-station ServiceNodeParam
        // (QueueNodeParam) that getServiceParam() reads to store nservertypes/heterorates.
        refreshHeterogeneousServers();
        refreshPetriNetNodes();
        jline.io.LineConsole.compileDetail("building the synchronization events");
        refreshSync();
        refreshGlobalSync();
        refreshRegions();
        this.hasStruct = true;
        // Reward definitions live on the model; copy them into the (re)built struct so
        // they survive resetStruct()/refreshStruct() and are visible to reward analyzers.
        this.sn.reward = this.rewardFunctions;
        this.sn.isfjaugmented = this.isFJAugmented;

        if (this.sn.fj.any() && !this.isFJAugmented) {
            // skipped on FJ tag-augmented copies, where ModelAdapter.fjtag
            // overwrites the auxiliary-class visits explicitly
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
                    // see _kb/04-networkstruct.md (refreshStruct.m: fjclassmap pairing) for rationale
                    int a = (int) fsn.inchain.get(newChain).get(jaux);
                    if (a >= fjclassmap.getNumCols() || fjclassmap.get(0, a) < 0) {
                        continue; // not an auxiliary class
                    }
                    int j = (int) fjclassmap.get(0, a);
                    Matrix Y = Vaux.getColumn(jaux).scale(fanOut);
                    sn.nodevisits.get(origChain).setColumn(j, X.getColumn(j).add(Y));
                }
            }

            // see _kb/04-networkstruct.md (refreshStruct.m field-population notes) for rationale
        }

        // Needs the visit ratios, so it runs after refreshChains.
        checkServiceReachable();
    }

    /**
     * Refuses a class that is ROUTED TO a station which cannot serve it.
     *
     * <p>sanitize() disables the OUTGOING routing of a class a station cannot
     * serve, which is what keeps it out of that station's visit ratios -- but
     * nothing stopped the class being routed IN, and a class that arrives where
     * it cannot be served is a flow sink: it enters and never leaves. The
     * station-level guard next to it cannot see this, because it asks whether the
     * station serves ANY class, not whether it serves the classes that reach
     * it.</p>
     *
     * <p>One such model gave three different wrong answers, none flagged, on a
     * closed cycle D &lt;-&gt; Q whose class C2 has no service at Q: MVA reported
     * Q/C2 with ArvR 1 against Tput 0, CTMC dropped class C2 entirely, and SSA
     * returned D/C2 QLen 2e-06 with the Q rows absent.</p>
     *
     * <p>AN ABSENT SERVICE AND AN EXPLICIT Disabled ARE TREATED ALIKE. sanitize
     * fills an absent slot with Disabled, so the two reach every solver as the
     * identical struct and produce identical wrong numbers; the spelling cannot
     * decide. What decides is whether anything ROUTES THE CLASS IN, which is what
     * this tests, so the class-switching idiom -- each queue serving one class and
     * marking the rest Disabled -- has no incoming flow there and is untouched.</p>
     *
     * <p>READS sn.rtnodes, WALKED FORWARD FROM THE FEED POINTS -- not sn.nodevisits,
     * which this guard read until 2026-09-02 and which 94d5570f3 had made blind to
     * the very case it exists for. That commit extended the `served` mask of
     * sn_refresh_visits from the station chain to the NODE chain, and it had to: on
     * a materialised LQN replica the unserved states close into a spurious cycle.
     * But the mask zeroes exactly the (station, class) cell a flow sink shows up in.
     * On Source -&gt; Q -&gt; Sink with class B unservable at Q, B's chain went from
     * Q = 1 to Q = 0 and the guard fell silent, while the Sink still read 1 -- flow
     * arriving downstream of a node it never visited. A MASKED VISIT VECTOR CANNOT
     * ANSWER THIS QUESTION, because the mask IS the answer being looked for. Do not
     * route this guard back through nodevisits or visits; both carry that mask.</p>
     *
     * <p>rtnodes on its own over-approximates -- it says where a class WOULD go if
     * one existed -- and the WALK is what removes the slack. It starts only at
     * (Source, class) pairs whose arrival is not Disabled, and at the reference
     * station of each closed class with a positive population, so the
     * Disabled-arrival row a class-switching Source carries is never entered. Three
     * rules keep it honest: a Sink is ABSORBING (rtnodes wraps it back to the Source
     * to close the kernel, and following that wrap re-enters every Source row,
     * including the Disabled ones the seeding just excluded); a Source is expanded
     * ONLY AS A SEED, for the same reason; and an unservable (station, class) is
     * REACHED BUT NOT EXPANDED, since nothing leaves it -- that is the whole
     * complaint -- so nothing downstream of it is evidence of anything.</p>
     *
     * <p>This SUBSUMES the fed-chain precondition the guard used to carry
     * separately: a chain no job can enter has no seed, so its rows are never walked
     * at all. That is strictly finer than the per-chain test it replaces, which
     * admitted every class of a chain any one of whose classes was fed.</p>
     *
     * <p>A SYNCHRONOUS REPLY IS NOT SERVED BY THE STATION IT RETURNS TO: it releases
     * the server that station held across the call, which is the whole content of
     * setSyncReply. Its Disabled service there is the marker of the feature, not a
     * flow sink, so the (station, reply class) pairs sn.replyblock marks are exempt.</p>
     */
    /**
     * Whether a job of class r can LEAVE node ind again.
     *
     * <p>True for anything that is not a service station, and for a station that
     * serves r, declares heterogeneous server types (its per-class slot is empty by
     * construction) or holds a server across a synchronous call whose reply class is
     * r. False only for the flow sink itself, which is what stops the walk in
     * reachedNodeClasses -- the same three exemptions the guard applies, kept in one
     * place so the walk and the verdict cannot drift apart.</p>
     *
     * @param ind 0-based node index
     * @param r   0-based class index
     * @return true when the pair is not a flow sink
     */
    private boolean servesClass(int ind, int r) {
        Node nd = this.nodes.get(ind);
        if (!(nd instanceof Queue) && !(nd instanceof Delay)) {
            return true;
        }
        if (nd instanceof Queue && !((Queue) nd).getServerTypes().isEmpty()) {
            return true;
        }
        JobClass jobclass = this.jobClasses.get(r);
        Station st = (Station) nd;
        if (st.getServer().containsJobClass(jobclass)
                && !(st.getServiceProcess(jobclass) instanceof Disabled)) {
            return true;
        }
        return holdsReplyFor(ind, r);
    }

    /**
     * The (node, class) pairs a job can actually ARRIVE at.
     *
     * <p>A forward walk of sn.rtnodes from the feed points. See
     * checkServiceReachable for why the evidence is the UNMASKED routing kernel
     * rather than sn.nodevisits, and for the three rules -- absorbing Sink, Source
     * expanded only as a seed, unservable pair reached but not expanded -- that keep
     * the walk from over-approximating.</p>
     *
     * @return an [nodes][classes] table, or null when rtnodes is unreadable or is
     *         not the expected (N*R) square, which leaves the caller checking
     *         nothing exactly as before
     */
    private boolean[][] reachedNodeClasses() {
        Matrix rt = this.sn == null ? null : this.sn.rtnodes;
        int N = this.nodes.size();
        int R = this.jobClasses.size();
        if (rt == null || N < 1 || R < 1
                || rt.getNumRows() < N * R || rt.getNumCols() < N * R) {
            return null;
        }
        boolean[][] reached = new boolean[N][R];
        boolean[][] seed = new boolean[N][R];
        Deque<int[]> stack = new ArrayDeque<int[]>();
        for (int ind = 0; ind < N; ind++) {
            Node nd = this.nodes.get(ind);
            if (!(nd instanceof Source)) {
                continue;
            }
            for (int r = 0; r < R; r++) {
                Distribution arv = ((Source) nd).getArrivalDistribution(this.jobClasses.get(r));
                if (arv != null && !(arv instanceof Disabled)) {
                    seed[ind][r] = true;
                }
            }
        }
        for (int r = 0; r < R; r++) {
            JobClass jobclass = this.jobClasses.get(r);
            if (!(jobclass instanceof ClosedClass)
                    || ((ClosedClass) jobclass).getPopulation() <= 0) {
                continue;
            }
            if (this.sn.refstat == null || r >= this.sn.refstat.getNumRows()
                    || this.sn.stationToNode == null) {
                continue;
            }
            int ist = (int) this.sn.refstat.get(r, 0);
            if (ist < 0 || ist >= this.sn.stationToNode.getNumRows()) {
                continue;
            }
            int ind = (int) this.sn.stationToNode.get(ist, 0);
            if (ind >= 0 && ind < N) {
                seed[ind][r] = true;
            }
        }
        for (int ind = 0; ind < N; ind++) {
            for (int r = 0; r < R; r++) {
                if (seed[ind][r] && !reached[ind][r]) {
                    reached[ind][r] = true;
                    stack.push(new int[]{ind, r});
                }
            }
        }
        while (!stack.isEmpty()) {
            int[] cur = stack.pop();
            int ind = cur[0];
            int r = cur[1];
            Node nd = this.nodes.get(ind);
            if (nd instanceof Sink) {
                continue;
            }
            if (nd instanceof Source && !seed[ind][r]) {
                continue;
            }
            if (!servesClass(ind, r)) {
                continue;
            }
            int row = ind * R + r;
            for (int col = 0; col < N * R; col++) {
                if (rt.get(row, col) <= GlobalConstants.Zero) {
                    continue;
                }
                int j = col / R;
                int sIdx = col % R;
                if (!reached[j][sIdx]) {
                    reached[j][sIdx] = true;
                    stack.push(new int[]{j, sIdx});
                }
            }
        }
        return reached;
    }

    /**
     * Whether node ind holds a server across a synchronous call whose reply class
     * is r, i.e. r returns there to RELEASE a server rather than to be served by
     * one.
     *
     * <p>sn.syncreply is indexed by the CALLING class and holds the 0-based reply
     * class, -1 where no reply is expected; sn.replyblock marks the (node, calling
     * class) pairs that hold a server across the call.</p>
     *
     * @param ind 0-based node index
     * @param r   0-based class index, tested as a reply class
     * @return true when a Disabled service for r at ind is the marker of a
     *         synchronous call rather than a flow sink
     */
    private boolean holdsReplyFor(int ind, int r) {
        if (this.sn.replyblock == null || this.sn.replyblock.isEmpty()
                || this.sn.syncreply == null || this.sn.syncreply.isEmpty()
                || ind >= this.sn.replyblock.getNumRows()) {
            return false;
        }
        int nk = FastMath.min(this.sn.syncreply.getNumRows(), this.sn.replyblock.getNumCols());
        for (int k = 0; k < nk; k++) {
            if ((int) this.sn.syncreply.get(k, 0) == r && this.sn.replyblock.get(ind, k) > 0) {
                return true;
            }
        }
        return false;
    }

    private void checkServiceReachable() {
        if (!this.enableChecks || this.sn == null) {
            return;
        }
        // Same exemption as the sanitize checks: a station of a cache, Petri-net
        // or fork-join model legitimately carries no per-class service.
        for (int ind = 0; ind < this.nodes.size(); ind++) {
            Node nd = this.nodes.get(ind);
            if (nd instanceof Cache || nd instanceof Place || nd instanceof Transition
                    || nd instanceof Fork || nd instanceof Join) {
                return;
            }
        }
        boolean[][] reached = reachedNodeClasses();
        if (reached == null) {
            return;
        }
        int K = this.jobClasses.size();
        for (int ind = 0; ind < this.nodes.size(); ind++) {
            Node nd = this.nodes.get(ind);
            if (!(nd instanceof Queue) && !(nd instanceof Delay)) {
                continue;
            }
            if (nd instanceof Queue && !((Queue) nd).getServerTypes().isEmpty()) {
                continue;
            }
            if (ind >= reached.length) {
                continue;
            }
            Station st = (Station) nd;
            for (int r = 0; r < K; r++) {
                if (r >= reached[ind].length || !reached[ind][r]) {
                    continue;
                }
                JobClass jobclass = this.jobClasses.get(r);
                if (st.getServer().containsJobClass(jobclass)
                        && !(st.getServiceProcess(jobclass) instanceof Disabled)) {
                    continue;
                }
                // A SYNCHRONOUS REPLY is not served by the station it returns to:
                // it releases the server that station held across the call, which
                // is the whole content of setSyncReply.
                if (holdsReplyFor(ind, r)) {
                    continue;
                }
                String kind = (nd instanceof Delay) ? "Delay" : "Queue";
                line_error(mfilename(new Object() {
                }), kind + " '" + nd.getName() + "' has no service configured for job class '"
                        + jobclass.getName() + "', but the class is routed to it. Jobs would arrive "
                        + "and never leave. Call setService() for that class, route it elsewhere, "
                        + "or disable this check with model.setChecks(false).");
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
                // see _kb/04-networkstruct.md (refreshGlobalSync/refreshSync section) for rationale
                if (sn.isstation.get(i, 0) > 0 && sn.impatienceType != null) {
                    jline.lang.nodes.Station rst = sn.stations.get((int) sn.nodeToStation.get(0, i));
                    jline.lang.JobClass rjc = sn.jobclasses.get(r);
                    Map<JobClass, ProcessType> imap = sn.impatienceType.get(rst);
                    ProcessType itype = (imap != null) ? imap.get(rjc) : null;
                    if (itype == ProcessType.EXP) {
                        Sync synct = new Sync();
                        synct.active.put(0, new Event(EventType.RENEGE, i, r, NaN, new Matrix(0, 0), NaN, NaN));
                        synct.passive.put(0, new Event(EventType.LOCAL, local, r, 1.0, new Matrix(0, 0), NaN, NaN));
                        sync.put(sync.size(), synct);
                    }
                }
                // Retrial action (exponential retrial delay): an orbiting class-r
                // job retries entry; succeeds only when a server is free.
                if (sn.isstation.get(i, 0) > 0 && sn.retrialProc != null && sn.retrialType != null) {
                    jline.lang.nodes.Station qst = sn.stations.get((int) sn.nodeToStation.get(0, i));
                    jline.lang.JobClass qjc = sn.jobclasses.get(r);
                    Map<JobClass, jline.util.matrix.MatrixCell> rpmap = sn.retrialProc.get(qst);
                    jline.util.matrix.MatrixCell rproc = (rpmap != null) ? rpmap.get(qjc) : null;
                    Map<JobClass, ProcessType> rtmap = sn.retrialType.get(qst);
                    ProcessType rtype = (rtmap != null) ? rtmap.get(qjc) : null;
                    if (rproc != null && rproc.size() > 0 && rtype == ProcessType.EXP) {
                        Sync synct = new Sync();
                        synct.active.put(0, new Event(EventType.RETRY, i, r, NaN, new Matrix(0, 0), NaN, NaN));
                        synct.passive.put(0, new Event(EventType.LOCAL, local, r, 1.0, new Matrix(0, 0), NaN, NaN));
                        sync.put(sync.size(), synct);
                    }
                }
                // see _kb/04-networkstruct.md (refreshGlobalSync/refreshSync section) for rationale
                if (r == 0 && sn.isstation.get(i, 0) > 0 && sn.hasbreakdown != null
                        && i < sn.hasbreakdown.length() && sn.hasbreakdown.get(i) == 1) {
                    Sync failSync = new Sync();
                    failSync.active.put(0, new Event(EventType.FAILURE, i, r, NaN, new Matrix(0, 0), NaN, NaN));
                    failSync.passive.put(0, new Event(EventType.LOCAL, local, r, 1.0, new Matrix(0, 0), NaN, NaN));
                    sync.put(sync.size(), failSync);
                    Sync repSync = new Sync();
                    repSync.active.put(0, new Event(EventType.REPAIR, i, r, NaN, new Matrix(0, 0), NaN, NaN));
                    repSync.passive.put(0, new Event(EventType.LOCAL, local, r, 1.0, new Matrix(0, 0), NaN, NaN));
                    sync.put(sync.size(), repSync);
                }
                // see _kb/04-networkstruct.md (refreshGlobalSync/refreshSync section) for rationale
                if (sn.isstation.get(i, 0) > 0
                        && sn.sched.get(sn.stations.get((int) sn.nodeToStation.get(0, i))) == SchedStrategy.POLLING) {
                    jline.lang.state.Polling.Info pinfoSync = jline.lang.state.Polling.info(sn, i);
                    if (pinfoSync != null && pinfoSync.hasSw[r]) {
                        Sync synct = new Sync();
                        synct.active.put(0, new Event(EventType.SWITCH, i, r, NaN, new Matrix(0, 0), NaN, NaN));
                        synct.passive.put(0, new Event(EventType.LOCAL, local, r, 1.0, new Matrix(0, 0), NaN, NaN));
                        sync.put(sync.size(), synct);
                    }
                }
                if (sn.isstateful.get(i, 0) > 0) {
                    if (sn.nodetype.get(i).equals(NodeType.Fork)) {
                        // see _kb/04-networkstruct.md (refreshGlobalSync/refreshSync section) for rationale
                        continue;
                    }
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
                                        case SQ:
                                        case SDR:
                                            // see _kb/04-networkstruct.md (refreshGlobalSync/refreshSync section) for rationale
                                            final int isf_final = isf, jsf_final = jsf, r_final = r, s_final = s;
                                            synct.passive.put(0, new Event(EventType.ARV, j, s, ((pair) -> sn.rtfun.apply(pair).get(isf_final * nclasses + r_final, jsf_final * nclasses + s_final)), new Matrix(0, 0), NaN, NaN));
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
                    // see _kb/04-networkstruct.md (refreshGlobalSync/refreshSync section) for rationale
                    if (r > 0) continue;
                    if (sn.nodetype.get(ind) == NodeType.Transition) {
                        TransitionNodeParam transParam = (TransitionNodeParam) sn.nodeparam.get(sn.nodes.get(ind));
                        // First loop: all mode enabling events (matching MATLAB refreshGlobalSync ordering)
                        for (int m = 0; m < transParam.nmodes; m++) {
                            Matrix enablingMatrix = transParam.enabling.get(m);
                            List<Integer> enablingPlaces = new ArrayList<>();
                        // see _kb/04-networkstruct.md (SPN/GSPN state generation) for rationale
                        for (int ep = 0; ep < enablingMatrix.getNumRows(); ep++) {
                                boolean anyClass = false;
                                for (int rr = 0; rr < enablingMatrix.getNumCols(); rr++) {
                                    if (enablingMatrix.get(ep, rr) > 0) { anyClass = true; break; }
                                }
                                if (anyClass) {
                                    enablingPlaces.add(ep);
                                }
                            }
                            // see _kb/04-networkstruct.md (refreshGlobalSync/refreshSync section) for rationale
                            Matrix inhibitingMatrix = transParam.inhibiting.get(m);
                            List<Integer> inhibitingPlaces = new ArrayList<>();
                            for (int ip = 0; ip < inhibitingMatrix.getNumRows(); ip++) {
                                if (inhibitingMatrix.get(ip, 0) < Double.POSITIVE_INFINITY && !enablingPlaces.contains(ip)) {
                                    inhibitingPlaces.add(ip);
                                }
                            }

                            if (!enablingPlaces.isEmpty() || !inhibitingPlaces.isEmpty()) {
                                GlobalSync enableSync = new GlobalSync();
                                List<ModeEvent> activeEvents = new ArrayList<>();
                                activeEvents.add(new ModeEvent(EventType.ENABLE, ind, m, 1.0));
                                enableSync.setActive(activeEvents);

                                List<ModeEvent> passiveEvents = new ArrayList<>();
                                for (int ep : enablingPlaces) {
                                    passiveEvents.add(new ModeEvent(EventType.LOCAL, ep, m, 1.0));
                                }
                                for (int ip : inhibitingPlaces) {
                                    passiveEvents.add(new ModeEvent(EventType.LOCAL, ip, m, 1.0));
                                }
                                enableSync.setPassive(passiveEvents);
                                gsync.put(gsync.size(), enableSync);
                            }
                        }
                        // Second loop: all mode firing events (matching MATLAB refreshGlobalSync ordering)
                        for (int m = 0; m < transParam.nmodes; m++) {
                            Matrix enablingMatrix = transParam.enabling.get(m);
                            List<Integer> enablingPlaces = new ArrayList<>();
                        // see _kb/04-networkstruct.md (SPN/GSPN state generation) for rationale
                        for (int ep = 0; ep < enablingMatrix.getNumRows(); ep++) {
                                boolean anyClass = false;
                                for (int rr = 0; rr < enablingMatrix.getNumCols(); rr++) {
                                    if (enablingMatrix.get(ep, rr) > 0) { anyClass = true; break; }
                                }
                                if (anyClass) {
                                    enablingPlaces.add(ep);
                                }
                            }

                            Matrix firingMatrix = transParam.firing.get(m);
                            List<Integer> firingPlaces = new ArrayList<>();
                            for (int fp = 0; fp < firingMatrix.getNumRows(); fp++) {
                                boolean anyClassF = false;
                                for (int rr = 0; rr < firingMatrix.getNumCols(); rr++) {
                                    if (firingMatrix.get(fp, rr) > 0) { anyClassF = true; break; }
                                }
                                if (anyClassF) {
                                    firingPlaces.add(fp);
                                }
                            }
                            // see _kb/04-networkstruct.md (refreshGlobalSync/refreshSync section) for rationale
                            Matrix firingInhibitingMatrix = transParam.inhibiting.get(m);
                            List<Integer> firingInhibitingPlaces = new ArrayList<>();
                            for (int ip = 0; ip < firingInhibitingMatrix.getNumRows(); ip++) {
                                if (firingInhibitingMatrix.get(ip, 0) < Double.POSITIVE_INFINITY
                                        && !enablingPlaces.contains(ip) && !firingPlaces.contains(ip)) {
                                    firingInhibitingPlaces.add(ip);
                                }
                            }

                            if (!enablingPlaces.isEmpty() || !firingPlaces.isEmpty() || !firingInhibitingPlaces.isEmpty()) {
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
                                for (int ip : firingInhibitingPlaces) {
                                    passiveEvents.add(new ModeEvent(EventType.LOCAL, ip, m, 1.0));
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

    /**
     * Relink the network from a modified rtorig map.
     * Converts the map to a RoutingMatrix and calls link().
     * This matches MATLAB's relink(P) behavior when P comes from getLinkedRoutingMatrix.
     */
    public void relinkFromRtorig(Map<JobClass, Map<JobClass, Matrix>> rtorig) {
        resetNetwork();
        List<Node> nodes = this.getNodes();
        List<JobClass> classes = this.getJobClasses();
        RoutingMatrix P = new RoutingMatrix(this, classes, nodes);
        for (Map.Entry<JobClass, Map<JobClass, Matrix>> e1 : rtorig.entrySet()) {
            JobClass fromClass = e1.getKey();
            for (Map.Entry<JobClass, Matrix> e2 : e1.getValue().entrySet()) {
                JobClass toClass = e2.getKey();
                Matrix mat = e2.getValue();
                P.set(fromClass, toClass, mat);
            }
        }
        link(P);
    }

    /**
     * Updates sn.rtorig without replacing the entire NetworkStruct.
     * Used by RoutingMatrix.setRouting() to preserve computed fields (rates, isstatedep, etc.)
     * when re-linking during solver iteration.
     */
    public void updateRtorig(Map<JobClass, Map<JobClass, Matrix>> rtorig) {
        if (this.sn != null) {
            this.sn.rtorig = rtorig;
        } else {
            NetworkStruct newSn = new NetworkStruct();
            newSn.rtorig = rtorig;
            this.sn = newSn;
        }
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
                    // address this cache's own item-set row, not a linear index into the
                    // (itemSetIndex, class) map, which aliases wrongly once a second
                    // cache gives the map more than one row
                    int itemRow = cache.getItems().getIndex();
                    for (int k = 0; k < K; k++) {
                        if (cache.popularityGet(itemRow, k) == null) {
                            cache.popularitySet(itemRow, k, Disabled.getInstance());
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
                    // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
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
                    // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
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
                            // Open classes without arrivals should still route (for class switching pass-through)
                            source.setRouting(jobclass, (jobclass instanceof OpenClass) ? RoutingStrategy.RAND : RoutingStrategy.DISABLED);
                        }
                    }
                } else if (node instanceof Place) {
                    Place place = (Place) node;
                    for (int k = 0; k < K; k++) {
                        JobClass jobclass = this.jobClasses.get(k);
                        // MATLAB refreshCapacity.m defaults all droprule to WAITQ
                        place.setDropRule(jobclass, DropStrategy.WaitingQueue);
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
                rekeyRtorig(sig, concrete);
            }
        }
    }

    /**
     * Replaces a resolved signal's placeholder in the rtorig keys. rtorig is
     * populated at link() time and keyed by the class objects then in the
     * model, so once resolveSignals substitutes a placeholder with its
     * concrete OpenSignal/ClosedSignal every lookup through the new object
     * silently misses: in ModelAdapter.mmt this dropped all routes touching a
     * signal class, leaving its fork-join auxiliary class with no routing.
     */
    private void rekeyRtorig(JobClass oldKey, JobClass newKey) {
        if (this.sn == null || this.sn.rtorig == null) {
            return;
        }
        Map<JobClass, Map<JobClass, Matrix>> rtorig = this.sn.rtorig;
        Map<JobClass, Matrix> outer = rtorig.remove(oldKey);
        if (outer != null) {
            rtorig.put(newKey, outer);
        }
        for (Map<JobClass, Matrix> row : rtorig.values()) {
            Matrix inner = row.remove(oldKey);
            if (inner != null) {
                row.put(newKey, inner);
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
     * Returns whether validation checks are enabled for this network.
     *
     * @return true if validation checks are enabled
     */
    public boolean getChecks() {
        return this.enableChecks;
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
        // see _kb/04-networkstruct.md (Node/process construction notes: getUsedLangFeatures) for rationale
        if (feature.equals("MarkedMAP") || feature.equals("MarkedMMPP")) {
            feature = "MMAP";
        }
        usedFeatures.setTrue(feature);
    }

    // ========================================================================
    // SECTION 21: UTILITY METHODS
    // Helper and utility methods for various network operations
    // ========================================================================

    /**
     * Builds sn.sdr, the station-indexed twin of the Krzesinski state-dependent
     * routing structure declared on the entry center. The node-indexed copy stays
     * in sn.nodeparam.get(entry).sdr, where the routing function sub_sdr reads it.
     *
     * <p>A network admits one subnetwork Q(V,V) and every class routed by it must
     * declare the same one: the routing probabilities of Krzesinski (1987) are
     * chain independent, so a per-class topology has no product form.</p>
     */
    public void refreshStateDepRouting() {
        if (this.sn == null) {
            return;
        }
        this.sn.sdr = null;
        if (this.sn.nodeparam == null) {
            return;
        }
        jline.lang.StateDepRouting decl = null;
        String declName = null;
        for (int ind = 0; ind < this.sn.nnodes; ind++) {
            Node node = this.nodes.get(ind);
            NodeParam np = this.sn.nodeparam.get(node);
            if (np == null || np.sdr == null) {
                continue;
            }
            for (JobClass jobclass : np.sdr.keySet()) {
                jline.lang.StateDepRouting cand = np.sdr.get(jobclass);
                if (decl == null) {
                    decl = cand;
                    declName = node.getName();
                } else if (!decl.sameAs(cand)) {
                    throw new RuntimeException("Two different state-dependent routing structures are declared (nodes "
                            + declName + " and " + node.getName() + "). The routing probabilities of Krzesinski (1987) "
                            + "are chain independent, so a network admits one subnetwork Q(V,V) and every class routed "
                            + "by it must declare the same branches, nesting and coefficients.");
                }
            }
        }
        if (decl == null) {
            return;
        }
        jline.lang.StateDepRouting stationSdr = decl.toStationIndices(this.sn.nodeToStation, this.sn.nodenames);
        jline.api.pfqn.Pfqn_sdr.pfqn_sdrcoeff(stationSdr); // validates the declaration and its bounds
        stationSdr.entryNode = decl.entryNode;
        stationSdr.departureNode = decl.departureNode;
        stationSdr.branchNodes = decl.branchNodes;
        this.sn.sdr = stationSdr;
    }

    /**
     * Krzesinski (1987) product-form state-dependent routing, eq. (10).
     *
     * <p>ind is the entry center e of Q(V,V). The probability of proceeding to a
     * branch entry is a function of the total branch and subnetwork populations,
     * and the residual mass returns the customer to the departure center d, which
     * is the busy form of waiting of Section 2.5.</p>
     *
     * @param ind entry node index
     * @param jnd destination node index
     * @param r active class
     * @param s passive class
     * @param linksmat connection matrix
     * @param state_before state before the transition
     * @param state_after state after the transition
     * @return the routing probability from ind to jnd
     */
    public double sub_sdr(int ind, int jnd, int r, int s, Matrix linksmat,
                          Map<Node, Matrix> state_before, Map<Node, Matrix> state_after) {
        int isf = (int) this.sn.nodeToStateful.get(ind);
        Node statefulNode = this.getStatefulNodeFromIndex(isf);
        if (!state_before.containsKey(statefulNode)) {
            return FastMath.min(linksmat.get(ind, jnd), 1.0);
        }
        if (r != s) {
            return 0.0;
        }
        NodeParam np = this.sn.nodeparam.get(this.nodes.get(ind));
        jline.lang.StateDepRouting sdr = np.sdr.get(this.jobClasses.get(r));
        double[] n = new double[this.sn.nnodes];
        for (int knd = 0; knd < this.sn.nnodes; knd++) {
            int ksf = (int) this.sn.nodeToStateful.get(0, knd);
            if (ksf < 0) {
                continue;
            }
            Node kNode = this.getStatefulNodeFromIndex(ksf);
            if (!state_before.containsKey(kNode)) {
                continue;
            }
            n[knd] = ToMarginal.toMarginal(this.sn, knd, state_before.get(kNode), null, null, null, null, null).ni.value();
        }
        jline.api.pfqn.Pfqn_sdr.Coeff c = jline.api.pfqn.Pfqn_sdr.pfqn_sdrcoeff(sdr);
        double[] P = jline.api.pfqn.Pfqn_sdr.pfqn_sdrprob(c, n);
        double p = 0.0;
        for (int b = 1; b < sdr.branch.length; b++) {
            if (sdr.entryOf[b] == jnd) {
                p += P[b];
            }
        }
        if (sdr.departure == jnd) {
            p += jline.api.pfqn.Pfqn_sdr.pfqn_sdrped(P);
        }
        return p;
    }

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
                        n.set(0, knd, ToMarginal.toMarginal(this.sn, knd, state_before.get(statefulNode_knd), null, null, null, null, null).ni.value());
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

    /**
     * Cache for the marginal-probability vector returned by sub_sq.
     * Key: serialized (k, m, memPos, n[]); value: per-eligible-destination
     * marginal probabilities. Same n vector reused across all (ind, jnd) cells
     * within a single rtfun invocation, so a few thousand entries amortize
     * across an entire long simulation.
     */
    private final java.util.concurrent.ConcurrentHashMap<String, double[]> sqProbCache =
            new java.util.concurrent.ConcurrentHashMap<String, double[]>();

    /**
     * Power-of-K choices marginal routing probability. Matches LDES semantics
     * in Solver_ssj.kt:selectSQDestination: enumerate all m^k
     * ordered tuples of eligible destinations sampled WITH replacement, break
     * ties by first occurrence in the tuple, and return the fraction of tuples
     * for which jnd is the JSQ winner. With memory, the prior pick is forced
     * as the last candidate and only m^(k-1) tuples are enumerated.
     */
    public double sub_sq(int ind, int jnd, int r, int s, Matrix linksmat,
                               Map<Node, Matrix> state_before, Map<Node, Matrix> state_after) {
        int isf = (int) this.sn.nodeToStateful.get(ind);
        Node statefulNode = this.getStatefulNodeFromIndex(isf);
        Matrix stateBefore = state_before.get(statefulNode);
        if (stateBefore == null || stateBefore.isEmpty()) {
            return FastMath.min(linksmat.get(ind, jnd), 1.0);
        }
        if (r != s) {
            return 0.0;
        }
        List<Integer> eligible = new ArrayList<Integer>();
        for (int knd = 0; knd < this.sn.nnodes; knd++) {
            if (linksmat.get(ind, knd) > 0) {
                eligible.add(knd);
            }
        }
        int m = eligible.size();
        if (m == 0 || !eligible.contains(jnd)) {
            return 0.0;
        }
        int kval = 2;
        Node fromNode = this.nodes.get(ind);
        NodeParam np = (this.sn.nodeparam != null) ? this.sn.nodeparam.get(fromNode) : null;
        JobClass jc = this.jobClasses.get(r);
        if (np != null && np.d != null && np.d.containsKey(jc)) {
            kval = np.d.get(jc);
        }
        if (kval < 1) kval = 1;
        if (kval > m) kval = m;

        double[] n = new double[m];
        for (int i = 0; i < m; i++) {
            int knd = eligible.get(i);
            Node statefulKnd = this.getStatefulNodeFromIndex((int) this.sn.nodeToStateful.get(0, knd));
            n[i] = ToMarginal.toMarginal(this.sn, knd, state_before.get(statefulKnd),
                    null, null, null, null, null).ni.value();
        }
        int jnd_pos = eligible.indexOf(jnd);

        // Build cache key from (d, ndest, n[]). SQ(d) carries no dispatcher
        // memory, so the marginal is a pure function of the queue lengths.
        StringBuilder kb = new StringBuilder(64);
        kb.append('k').append(kval).append('|').append('m').append(m).append('|');
        for (int i = 0; i < m; i++) {
            kb.append((int) n[i]).append(',');
        }
        String cacheKey = kb.toString();
        double[] pVec = sqProbCache.get(cacheKey);
        if (pVec == null) {
            pVec = new double[m];
            {
                long nTuples = 1;
                for (int c = 0; c < kval; c++) nTuples *= m;
                int[] tuple = new int[kval];
                for (long ti = 0; ti < nTuples; ti++) {
                    long rem = ti;
                    for (int c = 0; c < kval; c++) {
                        tuple[c] = (int) (rem % m);
                        rem /= m;
                    }
                    double minval = n[tuple[0]];
                    int winnerPos = 0;
                    for (int c = 1; c < kval; c++) {
                        if (n[tuple[c]] < minval) {
                            minval = n[tuple[c]];
                            winnerPos = c;
                        }
                    }
                    pVec[tuple[winnerPos]] += 1.0;
                }
                for (int i = 0; i < m; i++) pVec[i] /= (double) nTuples;
            }
            // Bound cache size — avoid unbounded growth in pathological cases.
            if (sqProbCache.size() > 50000) {
                sqProbCache.clear();
            }
            sqProbCache.put(cacheKey, pVec);
        }
        return pVec[jnd_pos];
    }

    /**
     * Column of a node state holding the round-robin pointer of class r at node ind.
     *
     * <p>The per-class routing pointers are appended in class order, and only for
     * the classes that actually route round-robin, followed by the node-level block
     * (nvars column 2R) if the node has one. The pointer is therefore located by
     * counting back from the end of the state: over the node-level block, then over
     * the pointers of the classes after r.</p>
     *
     * <p>Reading it as {@code state(end - R + r)} instead, as this did, assumes the
     * pointers are exactly the last R columns. That silently reads the wrong column
     * whenever the node carries a trailing block of its own (the polling controller
     * of a POLLING station, a cache, or a BAS blocked marker), and also whenever
     * only some of the classes route round-robin, in which case a pointer is
     * confused with a service phase.</p>
     *
     * @param ind    node index
     * @param r      class index
     * @param stateLen number of columns of the node state
     * @return the pointer column index
     */
    private int subRouteSlot(int ind, int r, int stateLen) {
        int R = this.sn.nclasses;
        int nodeblock = (int) this.sn.nvars.get(ind, 2 * R);
        int after = 0;
        for (int rr = r + 1; rr < R; rr++) {
            after += (int) this.sn.nvars.get(ind, R + rr);
        }
        return stateLen - nodeblock - after - 1;
    }

    public double sub_rr_wrr(int ind, int jnd, int r, int s, Matrix linksmat, Map<Node, Matrix> state_before, Map<Node, Matrix> state_after) {
        int R = this.sn.nclasses;
        int isf = (int) this.sn.nodeToStateful.get(ind);
        Node statefulNode = this.getStatefulNodeFromIndex(isf);
        // Check if state is empty (matches MATLAB: isempty(state_before{isf}))
        Matrix stateBefore = state_before.get(statefulNode);
        if (stateBefore == null || stateBefore.isEmpty()) {
            return FastMath.min(linksmat.get(ind, jnd), 1.0);
        } else {
            if (r == s) {
                // see _kb/04-networkstruct.md (State-space construction conventions) for rationale
                Matrix stateAfter = state_after.get(statefulNode);
                if (stateAfter == null || stateAfter.isEmpty()) {
                    return FastMath.min(linksmat.get(ind, jnd), 1.0);
                }
                int rrIdx = subRouteSlot(ind, r, (int) stateAfter.length());
                int slot = (int) stateAfter.get(rrIdx);
                // see _kb/04-networkstruct.md (State-space construction conventions) for rationale
                if (this.sn.routing.get(this.nodes.get(ind)).get(this.jobClasses.get(r)) == RoutingStrategy.WRROBIN) {
                    Matrix wol = this.sn.nodeparam.get(this.nodes.get(ind)).weightedOutlinks.get(this.jobClasses.get(r));
                    if (wol == null || wol.isEmpty()) {
                        return 0.0;
                    }
                    int cyc = (int) wol.length();
                    if (slot < 1 || slot > cyc) {
                        return 0.0;
                    }
                    return ((int) wol.get(slot - 1) == jnd) ? 1.0 : 0.0;
                }
                return (slot == jnd) ? 1.0 : 0.0;
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
     * Validates the routing probabilities supplied to link(): they must be
     * nonnegative, and the total leaving a node in a given class may not exceed 1.
     *
     * <p>Mirrors the equivalent block in the MATLAB {@code @MNetwork/link.m}.
     *
     * @param P the routing matrix to validate
     */
    private void validateRoutingProbabilities(RoutingMatrix P) {
        int K = this.jobClasses.size();
        int I = this.nodes.size();

        // see _kb/04-networkstruct.md (Routing-matrix validation) for rationale
        for (int r = 1; r <= K; r++) {
            for (int s = 1; s <= K; s++) {
                Matrix Prs = P.get(r, s);
                if (Prs == null) {
                    continue;
                }
                int nr = Math.min(I, Prs.getNumRows());
                int nc = Math.min(I, Prs.getNumCols());
                for (int i = 0; i < nr; i++) {
                    for (int j = 0; j < nc; j++) {
                        double prob = Prs.get(i, j);
                        if (prob < -GlobalConstants.FineTol) {
                            if (K > 1) {
                                line_error(mfilename(new Object() {
                                        }),
                                        String.format("Negative routing probability %g from node %s to node %s (class %s to class %s). Routing probabilities must be nonnegative.",
                                                prob, this.nodes.get(i).getName(), this.nodes.get(j).getName(),
                                                this.jobClasses.get(r - 1).getName(), this.jobClasses.get(s - 1).getName()));
                            } else {
                                line_error(mfilename(new Object() {
                                        }),
                                        String.format("Negative routing probability %g from node %s to node %s. Routing probabilities must be nonnegative.",
                                                prob, this.nodes.get(i).getName(), this.nodes.get(j).getName()));
                            }
                        }
                    }
                }
            }
        }

        // see _kb/04-networkstruct.md (Routing-matrix validation) for rationale
        for (int i = 0; i < I; i++) {
            Node node = this.nodes.get(i);
            // see _kb/04-networkstruct.md (Routing-matrix validation) for rationale
            if (node instanceof Fork || node instanceof Place
                    || node instanceof Transition || node instanceof Router) {
                continue;
            }
            // see _kb/04-networkstruct.md (Routing-matrix validation) for rationale
            if (node.getOutput() instanceof Forker) {
                continue;
            }
            if (node instanceof Station && ((Station) node).getSchedStrategy() == SchedStrategy.FORK) {
                continue;
            }
            for (int r = 1; r <= K; r++) {
                double totalProb = 0.0;
                for (int s = 1; s <= K; s++) {
                    Matrix Prs = P.get(r, s);
                    if (Prs == null) {
                        continue;
                    }
                    int nc = Math.min(I, Prs.getNumCols());
                    for (int j = 0; j < nc; j++) {
                        totalProb += Prs.get(i, j);
                    }
                }
                if (totalProb > 1.0 + GlobalConstants.FineTol) {
                    if (K > 1) {
                        line_error(mfilename(new Object() {
                                }),
                                String.format("The total routing probability for jobs leaving node %s in class %s is %g, which is greater than 1.0.",
                                        node.getName(), this.jobClasses.get(r - 1).getName(), totalProb));
                    } else {
                        line_error(mfilename(new Object() {
                                }),
                                String.format("The total routing probability for jobs leaving node %s is %g, which is greater than 1.0.",
                                        node.getName(), totalProb));
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
        if (this.rewardFunctions == null) {
            // LinkedHashMap, not HashMap: getRewardNames() and getAvgReward()
            // report in this map's iteration order, and the MATLAB and Python
            // references report in DECLARATION order. A HashMap here made the
            // JAR permute an otherwise identical reward vector, which is a
            // difference no value comparison catches and every transcript diff
            // does.
            this.rewardFunctions = new LinkedHashMap<String, RewardFunction>();
        }
        this.rewardFunctions.put(name, rewardFn);
        // Keep the current struct (if any) in sync so callers holding it see the reward.
        if (this.sn != null) {
            this.sn.reward = this.rewardFunctions;
        }
    }

    /**
     * Get all defined reward functions.
     *
     * @return Map from reward name to reward function, or null if no rewards defined
     */
    public Map<String, RewardFunction> getRewards() {
        return this.rewardFunctions;
    }

    /**
     * Get a specific reward function by name.
     *
     * @param name The reward name
     * @return The reward function, or null if not found
     */
    public RewardFunction getReward(String name) {
        if (this.rewardFunctions != null) {
            return this.rewardFunctions.get(name);
        }
        return null;
    }

    /**
     * Remove all defined reward functions.
     */
    public void clearRewards() {
        this.rewardFunctions = null;
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
        return this.rewardFunctions != null && !this.rewardFunctions.isEmpty();
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
        try {
            SolverJMT jmt = new SolverJMT(this);
            NetworkStruct sn = this.getStruct();
            jmt.writeJSIM(sn);
            String jsimFile = jmt.getFilePath() + java.io.File.separator + jmt.getFileName() + ".jsim";

            String jmtPath = jmtGetPath();
            String jmtJar = jmtPath + java.io.File.separator + "JMT.jar";
            if (!new java.io.File(jmtJar).exists()) {
                System.err.println("JMT.jar not found at: " + jmtJar);
                return;
            }
            ProcessBuilder pb = new ProcessBuilder("java", "-cp", jmtJar, "jmt.commandline.Jmt", "jsimg", jsimFile);
            pb.inheritIO();
            pb.start();
            // String viewerPath = lineViewerGetPath();
            // ProcessBuilder pb = new ProcessBuilder("java", "-jar", viewerPath, jsimFile);
            // pb.inheritIO();
            // pb.start();
        } catch (Exception e) {
            System.err.println("Failed to launch JSIMgraph viewer: " + e.getMessage());
        }
    }

}


