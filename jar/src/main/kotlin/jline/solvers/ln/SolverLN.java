/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ln;

import jline.VerboseLevel;
import jline.lang.*;
import jline.GlobalConstants;
import jline.lang.constant.*;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.LayeredNetworkElement;
import jline.lang.NetworkStruct;
import jline.lang.layered.LayeredNetworkStruct;
import jline.lang.nodes.*;
import jline.lang.nodes.Queue;
import jline.lang.processes.APH;
import jline.lang.processes.Disabled;
import jline.lang.processes.Distribution;
import jline.lang.processes.EmpiricalCDF;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.solvers.*;
import jline.solvers.auto.SolverAUTO;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.mam.SolverMAM;
import jline.solvers.mva.MVAOptions;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.solvers.qns.SolverQNS;
import jline.solvers.ssa.SolverSSA;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import java.util.*;

import static jline.api.mam.Aph_convseqKt.aph_convseq;
import static jline.api.mam.Aph_simplifyKt.aph_simplify;
import static jline.api.mc.Dtmc_makestochasticKt.dtmc_makestochastic;
import jline.io.Ret.DistributionResult;
import jline.util.Pair;
import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.line_warning;
import static jline.io.InputOutputKt.mfilename;
import static jline.util.Utils.isInf;
import static jline.util.matrix.Matrix.weaklyConnect;

/**
 * Solver for Layered Queueing Networks (LQN) using ensemble-based iterative methods.
 * 
 * <p>SolverLN implements layered queueing network analysis through decomposition into
 * simpler queueing models. LQNs extend traditional queueing networks by modeling
 * software systems with nested service requests, where servers can act as clients
 * to other services, creating layered dependencies.</p>
 * 
 * <p>Key LQN solver capabilities:
 * <ul>
 *   <li>Multi-layer model decomposition and iteration</li>
 *   <li>Software system modeling with nested service calls</li>
 *   <li>Client-server interaction patterns</li>
 *   <li>Convergence detection across model layers</li>
 *   <li>Ensemble-based performance analysis</li>
 * </ul>
 * </p>
 * 
 * <p>The solver iterates between layers, updating service demands and arrival rates
 * until convergence is achieved across all layers. This enables analysis of complex
 * distributed software architectures and service-oriented systems.</p>
 * 
 * @see jline.lang.layered.LayeredNetwork
 * @see EnsembleSolver
 * @see LNOptions
 * @since 1.0
 */
public class SolverLN extends EnsembleSolver {
    // registries of quantities to update at every iteration
    public int nlayers; // number of model layers
    public LayeredNetworkStruct lqn; // lqn data structure
    public boolean hasconverged; // true if last iteration converged, false otherwise
    public Integer averagingstart; // iteration at which result averaging started
    public List<Double> idxhash; // ensemble model associated to host or task
    public Matrix servtmatrix; // auxiliary matrix to determine entry servt
    public Matrix ptaskcallers; // probability that a task is called by a given task, directly or indirectly (remotely)
    public Map<Integer, Matrix> ptaskcallers_step; // probability that a task is called by a given task, directly or indirectly (remotely) up to a given step distance
    public Matrix ilscaling; // probability that a task is called by a given task, directly or indirectly (remotely) up to a given step distance
    public Matrix njobs; // number of jobs for each caller in a given submodel
    public Matrix njobsorig; // number of jobs for each caller at layer build time
    public List<Integer> routereset; // models that require hard reset of service chains
    public List<Integer> svcreset; // models that require hard reset of service process
    public List<Double> maxitererr; // maximum error at current iteration over all layers

    // performance metrics and related processes
    public Matrix util;
    //private Matrix util_ilock; // interlock matrix, element (i,j) says how much the utilization of task i is imputed to task j
    public Matrix tput;
    public Map<Integer, Distribution> tputproc;
    public Matrix servt; // this is the mean service time of an activity, which is the response time at the lower layer (if applicable)
    public Matrix residt; // this is the residence time at the lower layer (if applicable)
    public Map<Integer, Distribution> servtproc; // this is the service time process with mean fitted to the servt value
    public Map<Integer, Matrix> servtcdf; // CDF of the service time process (keyed by activity index)
    public Matrix thinkt;
    public Map<Integer, Distribution> thinkproc;
    public Map<Integer, Distribution> thinktproc;
    public Matrix callresidt;
    public Matrix callservt;
    public Map<Integer, Distribution> callservtproc;
    public Map<Integer, Matrix> callservtcdf; // CDF of the call service time (keyed by call index)
    public Matrix ignore;
    // registries of quantities to update at every iteration
    public Matrix arvproc_classes_updmap;
    public Matrix thinkt_classes_updmap;
    public Matrix servt_classes_updmap;
    public Matrix call_classes_updmap;
    public Matrix route_prob_updmap;
    public Matrix unique_route_prob_updmap; // auxiliary cache of unique route_prob_updmap rows
    // Temporary variables
    public Map<Integer, List<Integer[]>> cell_arvproc_classes_updmap; // [modelidx, actidx, node, class]
    public Map<Integer, List<Integer[]>> cell_thinkt_classes_updmap; // [modelidx, actidx, node, class]
    public Map<Integer, List<Integer[]>> cell_servt_classes_updmap; // [modelidx, actidx, node, class]
    public Map<Integer, List<Integer[]>> cell_call_classes_updmap; // [modelidx, callidx, node, class]
    public Map<Integer, List<Integer[]>> cell_route_prob_updmap; // [modelidx, actidxfrom, actidxto, nodefrom, nodeto, classfrom, classto]
    // Temporary variables used in buildLayers recursion
    public List<Network> temp_ensemble;
    public JobClass curClassC;
    public Map<Integer, APH> entryproc; // APH process for entries (keyed by entry index)
    private Matrix entrycdfrespt;

    // Under-relaxation state for convergence improvement
    public double relax_omega;        // Current relaxation factor
    public List<Double> relax_err_history;  // Error history for adaptive mode
    public Matrix servt_prev;         // Previous service times for relaxation
    public Matrix residt_prev;        // Previous residence times for relaxation
    public Matrix tput_prev;          // Previous throughputs for relaxation
    public Matrix thinkt_prev;        // Previous think times for relaxation
    public Matrix callservt_prev;     // Previous call service times for relaxation

    // MOL (Method of Layers) properties for hierarchical iteration
    public List<Integer> hostLayerIndices;   // Indices of host (processor) layers in ensemble
    public List<Integer> taskLayerIndices;   // Indices of task layers in ensemble
    public int molItHostOuter;               // MOL host outer iteration counter
    public int molItTaskInner;               // MOL task inner iteration counter

    // Solver factory for update_solver functionality
    public SolverFactory solverFactory;

    // Phase-2 support properties
    public boolean hasPhase2;           // Flag: model has phase-2 activities
    public Matrix servt_ph1;            // Phase-1 service time per activity (nidx x 1)
    public Matrix servt_ph2;            // Phase-2 service time per activity (nidx x 1)
    public Matrix util_ph1;             // Phase-1 utilization per entry
    public Matrix util_ph2;             // Phase-2 utilization per entry
    public Matrix prOvertake;           // Overtaking probability per entry (nentries x 1)

    public SolverLN(LayeredNetwork lqnmodel) {
        this(lqnmodel, new SolverOptions(SolverType.LN));
    }

    public SolverLN(LayeredNetwork lqnmodel, SolverFactory solverFactory) {
        this(lqnmodel, solverFactory, new SolverOptions(SolverType.LN));
    }

    public SolverLN(LayeredNetwork lqnmodel, SolverOptions options) {
        this(lqnmodel, new DefaultSolverFactory(), options);
    }

    public SolverLN(LayeredNetwork lqnmodel, SolverType solverType) {
        this(lqnmodel, solverType, new SolverOptions(SolverType.LN));
    }

    public SolverLN(LayeredNetwork lqnmodel, SolverType solverType, SolverOptions options) {
        this(lqnmodel, createSolverFactory(solverType), options);
    }

    public SolverLN(LayeredNetwork lqnmodel, SolverType solverType, LNOptions lnOptions, SolverOptions solverOptions) {
        this(lqnmodel, createSolverFactory(solverType, solverOptions), lnOptions);
    }

    public SolverLN(LayeredNetwork lqnmodel, SolverFactory solverFactory, SolverOptions options) {
        super(null, "SolverLN", options); // first argument is null as the ensemble cannot be built yet
        this.lqn = lqnmodel.getStruct();

        // Detect and initialize phase-2 support
        this.hasPhase2 = false;
        if (lqn.actphase != null) {
            for (int a = 1; a <= lqn.nacts; a++) {
                if (lqn.actphase.get(0, a) > 1) {
                    this.hasPhase2 = true;
                    break;
                }
            }
        }
        if (this.hasPhase2) {
            this.servt_ph1 = new Matrix(1, lqn.nidx, lqn.nidx);
            this.servt_ph2 = new Matrix(1, lqn.nidx, lqn.nidx);
            this.util_ph1 = new Matrix(1, lqn.nidx, lqn.nidx);
            this.util_ph2 = new Matrix(1, lqn.nidx, lqn.nidx);
            this.prOvertake = new Matrix(1, lqn.nentries, lqn.nentries);
        }

        construct();
        solvers = new NetworkSolver[nlayers];

        // Check if model has FunctionTasks (matches MATLAB SolverLN.m lines 138-148)
        boolean hasFunctionTasks = false;
        if (lqn.isfunction != null) {
            for (int i = 0; i < lqn.isfunction.getNumCols(); i++) {
                if (lqn.isfunction.get(0, i) == 1) {
                    hasFunctionTasks = true;
                    break;
                }
            }
        }

        for (int i = 0; i < nlayers; i++) {
            // For FunctionTask layers with setupTime, use SolverMAM with dec.poisson method
            // This is required because standard MVA doesn't handle setup/delayoff queue semantics
            // (matches MATLAB SolverLN.m lines 138-148: checks stations{2}.setupTime)
            if (hasFunctionTasks && ensemble[i].getStations().size() > 1) {
                // MATLAB checks stations{2} (1-indexed) = station index 1 (0-indexed)
                Station station2 = ensemble[i].getStations().get(1);
                if (station2 instanceof Queue) {
                    Queue serverQueue = (Queue) station2;
                    if (serverQueue.isDelayOffEnabled()) {
                        // Use SolverMAM with dec.poisson for FunctionTask layers
                        SolverOptions mamOptions = new SolverOptions(SolverType.MAM);
                        mamOptions.verbose = VerboseLevel.SILENT;
                        mamOptions.method = "dec.poisson";
                        solvers[i] = new SolverMAM(ensemble[i], mamOptions);
                        continue;
                    }
                }
            }
            solvers[i] = solverFactory.at(ensemble[i]);
        }
        this.solverFactory = solverFactory; // Store for later use
    }

    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.LN);
    }

    @Override
    public SolverResult analyze(int it, int e) {
        SolverResult result1 = new SolverResult();
        System.out.flush();
        solvers[e].getAvg();

        if (e == 1 && solvers[e].options.verbose != VerboseLevel.SILENT) {
            System.out.println();
        }

        result1.QN = solvers[e].result.QN;
        result1.UN = solvers[e].result.UN;
        result1.RN = solvers[e].result.RN;
        result1.TN = solvers[e].result.TN;
        result1.AN = solvers[e].result.AN;
        result1.WN = solvers[e].result.WN;
        result1.CN = solvers[e].result.CN; // not in MATLAB
        result1.XN = solvers[e].result.XN; // not in MATLAB

        return result1;
    }

    public void buildLayers() {

        this.temp_ensemble = new ArrayList<>();
        this.cell_servt_classes_updmap = new HashMap<>(lqn.nhosts + lqn.ntasks);
        this.cell_call_classes_updmap = new HashMap<>(lqn.nhosts + lqn.ntasks);
        this.cell_arvproc_classes_updmap = new HashMap<>(lqn.nhosts + lqn.ntasks);
        this.cell_thinkt_classes_updmap = new HashMap<>(lqn.nhosts + lqn.ntasks);
        this.cell_route_prob_updmap = new HashMap<>(lqn.nhosts + lqn.ntasks);

        double temp_idxhash = 1;

        // build one subnetwork for every processor
        for (int hidx = 1; hidx <= lqn.nhosts; hidx++) {
            if (this.ignore.get(hidx) == 0) {
                List<Integer> callers = lqn.tasksof.get(hidx);
                buildLayersRecursive(hidx, callers, true);
                this.idxhash.add(temp_idxhash);
                temp_idxhash++;
            } else {
                this.idxhash.add(Double.NaN);
            }
        }

        // build one subnetwork for every task
        for (int t = 1; t <= lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            // MATLAB: ~self.ignore(tidx) & ~lqn.isref(tidx) & ~(isempty(find(self.lqn.iscaller(tidx,:), 1)) & isempty(find(self.lqn.iscaller(:,tidx), 1)))
            boolean isolated_task = true;
            for (int i = 1; i <= lqn.nidx; i++) {
                if (lqn.iscaller.isAssigned(tidx, i) || lqn.iscaller.isAssigned(i, tidx)) {
                    isolated_task = false;
                    break;
                }
            }
            if (this.ignore.get(tidx) == 0 && (int) lqn.isref.get(tidx) == 0 && !isolated_task) {
                // obtain the activity graph of each task that calls some entry in t
                // [calling_idx, called_entries] = find(lqn.iscaller(:, lqn.entriesof{tidx}));
                List<Integer> calling_idx = new ArrayList<>();
                for (int eidx : lqn.entriesof.get(tidx)) {
                    for (int i = 0; i < lqn.iscaller.getNumRows(); i++) {
                        if (lqn.iscaller.get(i, eidx) != 0) {
                            calling_idx.add(i);
                        }
                    }
                }
                // callers = intersect(lqn.tshift+(1:lqn.ntasks), unique(calling_idx)');
                List<Integer> taskRange = new ArrayList<>();
                for (int i = lqn.tshift + 1; i <= lqn.tshift + lqn.ntasks; i++) {
                    taskRange.add(i);
                }
                List<Integer> callers = new ArrayList<>();
                for (int idx : calling_idx) {
                    if (taskRange.contains(idx) && !callers.contains(idx)) {
                        callers.add(idx);
                    }
                }
                if (!callers.isEmpty()) {
                    buildLayersRecursive(tidx, callers, false);
                    idxhash.add(temp_idxhash);
                    temp_idxhash++;

                } else {
                    idxhash.add(Double.NaN);
                }

            } else {
                idxhash.add(Double.NaN);
            }
        }
        thinkt_classes_updmap = integerMapToMatrix(cell_thinkt_classes_updmap);
        call_classes_updmap = integerMapToMatrix(cell_call_classes_updmap);
        servt_classes_updmap = integerMapToMatrix(cell_servt_classes_updmap);
        arvproc_classes_updmap = integerMapToMatrix(cell_arvproc_classes_updmap);
        route_prob_updmap = integerMapToMatrix(cell_route_prob_updmap);

        this.ensemble = new Network[nlayers];

        for (int i = 0; i < nlayers; i++) {
            ensemble[i] = temp_ensemble.get(i);
        }

        // Classify layers as host (processor) or task for MOL iteration
        this.hostLayerIndices = new ArrayList<>();
        this.taskLayerIndices = new ArrayList<>();

        // Host layers: indices 1:nhosts (before idxhash remapping)
        // Note: idxhash uses 1-based indexing (NaN at index 0), so use hidx directly
        for (int hidx = 1; hidx <= lqn.nhosts; hidx++) {
            if (!Double.isNaN(idxhash.get(hidx))) {
                hostLayerIndices.add(idxhash.get(hidx).intValue() - 1); // Convert to 0-indexed
            }
        }

        // Task layers: indices tshift+1:tshift+ntasks
        // Note: idxhash uses 1-based indexing (NaN at index 0), so use tidx directly
        for (int t = 1; t <= lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            if (tidx < idxhash.size() && !Double.isNaN(idxhash.get(tidx))) {
                taskLayerIndices.add(idxhash.get(tidx).intValue() - 1); // Convert to 0-indexed
            }
        }
    }

    public void buildLayersRecursive(int idx, List<Integer> callers, boolean ishostlayer) {
        nlayers++;
        Matrix jobPosKey = new Matrix(1, lqn.nidx + 1);
        Map<Integer, JobClass> curClassKey = new HashMap<>(lqn.nidx);
        int nreplicas = (int) lqn.repl.get(0, idx);
        Matrix mult = lqn.maxmult.copy(); // this removes spare capacity that cannot be used
        lqn.mult = mult.copy(); // using maxmult as in MATLAB version
        Network model = new Network(lqn.hashnames.get(idx));
        model.setChecks(false);

        boolean hasSynccaller = false;
        if (!ishostlayer) {
            for (int caller : callers) {
                for (int entries : lqn.entriesof.get(idx)) {
                    if (lqn.issynccaller.isAssigned(caller, entries)) {
                        hasSynccaller = true;
                    }
                }
            }
        }
        boolean hasAsynccaller = false;
        if (!ishostlayer) {
            for (int caller : callers) {
                for (int entries : lqn.entriesof.get(idx)) {
                    if (lqn.isasynccaller.isAssigned(caller, entries)) {
                        hasAsynccaller = true;
                    }
                }
            }
        }
        Delay clientDelay = null;
        if (ishostlayer || hasSynccaller || hasAsynccaller) {
            clientDelay = new Delay(model, "Clients");
            model.getAttribute().setClientIdx(1);
            model.getAttribute().setServerIdx(2);
            model.getAttribute().setSourceIdx(-1);
        } else {
            model.getAttribute().setSourceIdx(-1);
            model.getAttribute().setServerIdx(1);
            model.getAttribute().setClientIdx(-1);
        }

        Map<Integer, Queue> serverStation = new HashMap<Integer, Queue>(nreplicas);
        for (int i = 1; i <= nreplicas; i++) {
            if (i == 1) {
                serverStation.put(i, new Queue(model, lqn.hashnames.get(idx), lqn.sched.get(idx)));
            } else {
                String name = lqn.hashnames.get(idx).concat(".");
                serverStation.put(i, new Queue(model, name.concat(Integer.toString(i)), lqn.sched.get(idx)));
            }
            serverStation.get(i).setNumberOfServers((int) mult.get(0, idx));
            serverStation.get(i).getAttribute().setIsHost(ishostlayer);
            serverStation.get(i).getAttribute().setIdx(idx);
        }

        Cache cacheNode = null;
        
        // Detect cache layer: all callers must be cache tasks and this must be a host layer
        boolean iscachelayer = false;
        if (ishostlayer) {
            boolean tempIsCacheLayer = true;
            for (int caller : callers) {
                if (lqn.iscache.get(0, caller) == 0) {
                    tempIsCacheLayer = false;
                    break;
                }
            }
            iscachelayer = tempIsCacheLayer;
        }
        
        // Detect function layer: all callers must be function tasks and this must be a host layer
        boolean isfunctionlayer = false;
        if (ishostlayer) {
            boolean tempIsFunctionLayer = true;
            for (int caller : callers) {
                if (lqn.isfunction != null && lqn.isfunction.get(0, caller) == 0) {
                    tempIsFunctionLayer = false;
                    break;
                }
            }
            isfunctionlayer = tempIsFunctionLayer;
        }
        
        if (iscachelayer) {
            // Create cache node - use first caller to get cache parameters
            int callerIdx = callers.get(0);
            ReplacementStrategy replStrategy;
            int strategyConstant = (int) lqn.replacestrat.get(0, callerIdx);
            switch (strategyConstant) {
                case 0: replStrategy = ReplacementStrategy.RR; break;
                case 1: replStrategy = ReplacementStrategy.FIFO; break;
                case 2: replStrategy = ReplacementStrategy.SFIFO; break;
                case 3: replStrategy = ReplacementStrategy.LRU; break;
                default: replStrategy = ReplacementStrategy.RR; break;
            }
            // Convert multi-level cache array to Matrix
            int[] levelCaps = lqn.itemcap.get(callerIdx);
            Matrix itemLevelCapMatrix = new Matrix(1, levelCaps.length);
            for (int i = 0; i < levelCaps.length; i++) {
                itemLevelCapMatrix.set(0, i, levelCaps[i]);
            }
            cacheNode = new Cache(model, lqn.hashnames.get(callerIdx),
                                 (int) lqn.nitems.get(0, callerIdx),
                                 itemLevelCapMatrix,  // Pass Matrix for multi-level cache support
                                 replStrategy);
        }

        List<Integer> actsInCaller = new ArrayList<>();
        for (int i : callers) {
            actsInCaller.addAll(lqn.actsof.get(i));
        }


        Matrix isPostAndAct = new Matrix(lqn.actposttype.getNumRows(), lqn.actposttype.getNumCols(), lqn.nacts);
        Matrix isPreAndAct = new Matrix(lqn.actpretype.getNumRows(), lqn.actpretype.getNumCols(), lqn.actpretype.getNonZeroLength());
        for (int i = lqn.nhosts + lqn.ntasks + lqn.nentries + 1; i <= lqn.nidx; i++) {
            if (lqn.actposttype.get(0, i) == ActivityPrecedenceType.ID_POST_AND) {
                isPostAndAct.set(0, i, 1);
            }
            if (lqn.actpretype.get(0, i) == ActivityPrecedenceType.ID_PRE_AND) {
                isPreAndAct.set(0, i, 1);
            }
        }

        boolean hasFork = false;
        for (int i : actsInCaller) {
            if (isPostAndAct.get(0, i) != 0) {
                hasFork = true;
                break;
            }
        }
        int maxfanout = 1;
        for (int aidx : actsInCaller) {
            List<Integer> successors = new ArrayList<>();
            for (int i = 1; i < lqn.graph.getNumCols(); i++) {
                if (lqn.graph.get(aidx, i) != 0) {
                    successors.add(i);
                }
            }
            int postAndCount = 0;
            for (int i : successors) {
                if (isPostAndAct.get(0, i) == 1) {
                    postAndCount++;
                }
            }
            if (postAndCount > 0) {
                maxfanout = FastMath.max(maxfanout, postAndCount);
            }
        }


        Fork forkNode = null;
        Stack<JobClass> forkClassStack = null;
        Map<Integer, Router> forkOutputRouter = new HashMap<>(maxfanout);
        if (hasFork) {
            forkNode = new Fork(model, "Fork_PostAnd");
            for (int f = 1; f <= maxfanout; f++) {
                forkOutputRouter.put(f, new Router(model, "Fork_PostAnd_" + f));
            }
            forkClassStack = new Stack<>();
        }


        boolean hasJoin = false;
        for (int i : actsInCaller) {
            if (isPreAndAct.get(0, i) != 0) {
                hasJoin = true;
                break;
            }
        }

        Join joinNode = null;
        if (hasJoin) {
            joinNode = new Join(model, "Join_PreAnd", forkNode);
        }

        Map<Integer, JobClass> aidxclass = new HashMap<>(lqn.nentries + lqn.nacts);
        Map<Integer, JobClass> cidxclass = new HashMap<>();
        Map<Integer, JobClass> cidxauxclass = new HashMap<>();

        if (ishostlayer) {
            model.getAttribute().addHosts(new Integer[]{null, model.getAttribute().getServerIdx()});
        } else {
            model.getAttribute().addTasks(new Integer[]{null, model.getAttribute().getServerIdx()});
        }

        Source sourceStation = null;
        Sink sinkStation = null;
        Matrix openClasses = new Matrix(lqn.ncalls + 1, 4, lqn.ncalls * 3);
        final int[] openClassesAssignedLine = {0}; // Use array to make it effectively final
        Matrix entryOpenClasses = new Matrix(lqn.nentries + 1, 3, lqn.nentries); // track entry-level open arrivals
        final int[] entryOpenClassesAssignedLine = {0}; // counter for entry open classes
        //  first pass: create the classes
        double njobs;
        Map<Integer, Double> callmean = new HashMap<>();
        for (int tidx_caller : callers) {
            // Per-caller: check if this caller is a sync caller to entries of idx (task layer path)
            // For host layers, idx is a host index with no entries, so isSyncCallerToEntries is always false
            boolean isSyncCallerToEntries = false;
            if (!ishostlayer && lqn.entriesof.get(idx) != null) {
                for (int entry_idx : lqn.entriesof.get(idx)) {
                    if (lqn.issynccaller.get(tidx_caller, entry_idx) != 0) {
                        isSyncCallerToEntries = true;
                        break;
                    }
                }
            }
            // Per-caller: for host layers, check if this caller task has direct callers
            // (matching MATLAB buildLayersRecursive hasDirectCallers logic)
            boolean hasDirectCallers = false;
            if (ishostlayer) {
                if (lqn.isref.get(tidx_caller) != 0) {
                    hasDirectCallers = true;
                } else {
                    for (int eidx : lqn.entriesof.get(tidx_caller)) {
                        // Check if any task is a sync or async caller to this entry
                        for (int row = 1; row <= lqn.ntasks; row++) {
                            int tidx_row = row + lqn.tshift;
                            if (lqn.issynccaller.get(tidx_row, eidx) != 0 || lqn.isasynccaller.get(tidx_row, eidx) != 0) {
                                hasDirectCallers = true;
                                break;
                            }
                        }
                        if (hasDirectCallers) break;
                        // Check for open arrivals on this entry
                        if (lqn.arrival != null && lqn.arrival.containsKey(eidx) && lqn.arrival.get(eidx) != null) {
                            hasDirectCallers = true;
                            break;
                        }
                    }
                }
            }
            if ((ishostlayer && hasDirectCallers) || isSyncCallerToEntries) {
                if (this.njobs.get(tidx_caller, idx) == 0) {
                    njobs = mult.get(0, tidx_caller) * lqn.repl.get(0, tidx_caller);
                    if (isInf(njobs)) {
                        njobs = 0;
                        for (int row = 0; row < lqn.taskgraph.getNumRows(); row++) {
                            if ((int) lqn.taskgraph.get(row, tidx_caller) != 0) {
                                int caller_of_tidx_caller = row;
                                njobs += mult.get(0, caller_of_tidx_caller);
                            }
                        }
                        if (isInf(njobs)) {
                            // If the callers of tidx_caller are inf servers, then use a heuristic
                            njobs = 0;
                            for (int i = 1; i < mult.getNumCols(); i++) {
                                if (!isInf(mult.get(0, i)) && !Double.isNaN(mult.get(0, i))) {
                                    if (i < lqn.repl.getNumCols()) {
                                        njobs = njobs + mult.get(0, i) * lqn.repl.get(0, i);
                                    }
                                }
                            }
                            njobs = FastMath.min(njobs, 1000);  // Match MATLAB cap
                        }
                    }
                    this.njobs.set(tidx_caller, idx, njobs);
                } else {
                    njobs = this.njobs.get(tidx_caller, idx);
                }
                String caller_name = lqn.hashnames.get(tidx_caller);
                aidxclass.put(tidx_caller, new ClosedClass(model, caller_name, njobs, clientDelay));
                clientDelay.setService(aidxclass.get(tidx_caller), Disabled.getInstance());
                for (int m = 1; m <= nreplicas; m++) {
                    serverStation.get(m).setService(aidxclass.get(tidx_caller), Disabled.getInstance());
                }
                aidxclass.get(tidx_caller).setReferenceClass(true);
                aidxclass.get(tidx_caller).setAttribute(new Integer[]{LayeredNetworkElement.TASK, tidx_caller});
                aidxclass.get(tidx_caller).setCompletes(false);
                model.getAttribute().addTasks(new Integer[]{aidxclass.get(tidx_caller).getIndex(), tidx_caller});
                assert clientDelay != null;
                clientDelay.setService(aidxclass.get(tidx_caller), thinkproc.get(tidx_caller));
                if (lqn.isref.get(tidx_caller) == 0) {
                    if (!cell_thinkt_classes_updmap.containsKey(idx)) {
                        cell_thinkt_classes_updmap.put(idx, new ArrayList<>());
                    }
                    cell_thinkt_classes_updmap.get(idx).add(new Integer[]{idx, tidx_caller, 1, aidxclass.get(tidx_caller).getIndex()});
                }

                for (int eidx : lqn.entriesof.get(tidx_caller)) {
                    aidxclass.put(eidx, new ClosedClass(model, lqn.hashnames.get(eidx), 0, clientDelay));
                    clientDelay.setService(aidxclass.get(eidx), Disabled.getInstance());
                    for (int m = 1; m <= nreplicas; m++) {
                        serverStation.get(m).setService(aidxclass.get(eidx), Disabled.getInstance());
                    }
                    aidxclass.get(eidx).setCompletes(false);
                    aidxclass.get(eidx).setAttribute(new Integer[]{LayeredNetworkElement.ENTRY, eidx});
                    model.getAttribute().addEntries(new Integer[]{aidxclass.get(eidx).getIndex(), eidx});
                    clientDelay.setService(aidxclass.get(eidx), Immediate.getInstance());

                    // Check for open arrival distribution on this entry
                    if (lqn.arrival != null && lqn.arrival.containsKey(eidx) &&
                        lqn.arrival.get(eidx) != null) {

                        if (sourceStation == null) {
                            model.getAttribute().setSourceIdx(model.getNumberOfNodes() + 1);
                            sourceStation = new jline.lang.nodes.Source(model, "Source");
                            sinkStation = new Sink(model, "Sink");
                        }

                        // Create open class for this entry
                        OpenClass openClassForEntry = new OpenClass(model, lqn.hashnames.get(eidx) + "_Open", 0);
                        sourceStation.setArrival(openClassForEntry, lqn.arrival.get(eidx));
                        clientDelay.setService(openClassForEntry, Disabled.getInstance());
                        for (int m = 1; m <= nreplicas; m++) {
                            serverStation.get(m).setService(openClassForEntry, servtproc.get(eidx));
                        }

                        // Track for routing setup later: [class_index, entry_index]
                        entryOpenClassesAssignedLine[0]++;
                        entryOpenClasses.set(entryOpenClassesAssignedLine[0], 1, openClassForEntry.getIndex());
                        entryOpenClasses.set(entryOpenClassesAssignedLine[0], 2, eidx);

                        // Track: Use negative entry index to distinguish from call arrivals
                        if (!cell_arvproc_classes_updmap.containsKey(idx)) {
                            cell_arvproc_classes_updmap.put(idx, new ArrayList<>());
                        }
                        cell_arvproc_classes_updmap.get(idx).add(new Integer[]{idx, -eidx,
                            model.getNodeIndex(sourceStation), openClassForEntry.getIndex()});

                        openClassForEntry.setCompletes(false);
                        openClassForEntry.setAttribute(new Integer[]{LayeredNetworkElement.ENTRY, eidx});
                    }
                }
            }
            // for each activity of the calling task
            for (int aidx : lqn.actsof.get(tidx_caller)) {
                if (ishostlayer || isSyncCallerToEntries) {
                    aidxclass.put(aidx, new ClosedClass(model, lqn.hashnames.get(aidx), 0, clientDelay));
                    clientDelay.setService(aidxclass.get(aidx), Disabled.getInstance());
                    for (int m = 1; m <= nreplicas; m++) {
                        serverStation.get(m).setService(aidxclass.get(aidx), Disabled.getInstance());
                    }
                    aidxclass.get(aidx).setCompletes(false);
                    aidxclass.get(aidx).setAttribute(new Integer[]{LayeredNetworkElement.ACTIVITY, aidx});
                    model.getAttribute().addActivities(new Integer[]{aidxclass.get(aidx).getIndex(), aidx});
                    if (!(ishostlayer && (lqn.parent.get(0, (int) lqn.parent.get(0, aidx)) == idx))) {
                        clientDelay.setService(aidxclass.get(aidx), servtproc.get(aidx));
                    }
                    if (iscachelayer) {
                        // Set cache read item entry for activities in cache layer
                        for (int eidx = 1; eidx <= lqn.nentries; eidx++) {
                            if (lqn.graph.get(lqn.eshift + eidx, aidx) != 0) {
                                clientDelay.setService(aidxclass.get(aidx), servtproc.get(aidx));
                                break;
                            }
                        }
                    }
                }
                // add a class for each outgoing call from this activity
                for (int cidx : lqn.callsof.get(aidx)) {
                    callmean.put(cidx, lqn.callproc_mean.getOrDefault(cidx, Double.NaN));
                    if (lqn.calltype.get(cidx) == CallType.ASYNC) {
                        int serverIdx = (int) lqn.callpair.get(cidx, 2);
                        if (serverIdx >= 0 && lqn.parent.get(0, serverIdx) == idx) {
                            if (sourceStation == null) {
                                model.getAttribute().setSourceIdx(model.getNumberOfNodes() + 1);
                                sourceStation = new jline.lang.nodes.Source(model, "Source");
                                sinkStation = new Sink(model, "Sink");
                            }
                            cidxclass.put(cidx, new OpenClass(model, lqn.callhashnames.get(cidx), 0));
                            sourceStation.setArrival(cidxclass.get(cidx), Immediate.getInstance());
                            clientDelay.setService(cidxclass.get(cidx), Disabled.getInstance());
                            for (int m = 1; m <= nreplicas; m++) {
                                serverStation.get(m).setService(cidxclass.get(cidx), Immediate.getInstance());
                            }
                            openClassesAssignedLine[0]++;
                            openClasses.set(openClassesAssignedLine[0], 1, cidxclass.get(cidx).getIndex());
                            openClasses.set(openClassesAssignedLine[0], 2, callmean.get(cidx));
                            openClasses.set(openClassesAssignedLine[0], 3, cidx);
                            model.getAttribute().addCalls(new Integer[]{cidxclass.get(cidx).getIndex(), cidx, (int) lqn.callpair.get(cidx, 1), (int) lqn.callpair.get(cidx, 2)});
                            cidxclass.get(cidx).setCompletes(false);
                            cidxclass.get(cidx).setAttribute(new Integer[]{LayeredNetworkElement.CALL, cidx});
                            double minRespT = 0;
                            if (lqn.actsof.get(idx) != null) {
                                for (int tidx_act : lqn.actsof.get(idx)) {
                                    minRespT = minRespT + lqn.hostdem_mean.getOrDefault(tidx_act, 0.0); // upper bound, uses all activities not just the ones reachable by this entry
                                }
                            }
                            for (int m = 1; m <= nreplicas; m++) {
                                if (minRespT == 0) {
                                    serverStation.get(m).setService(cidxclass.get(cidx), Immediate.getInstance());
                                } else {
                                    serverStation.get(m).setService(cidxclass.get(cidx), Exp.fitMean(minRespT));
                                }
                            }
                        }
                    } else if (lqn.calltype.get(cidx) == CallType.SYNC) {
                        cidxclass.put(cidx, new ClosedClass(model, lqn.callhashnames.get(cidx), 0, clientDelay));
                        clientDelay.setService(cidxclass.get(cidx), Disabled.getInstance());
                        for (int m = 1; m <= nreplicas; m++) {
                            serverStation.get(m).setService(cidxclass.get(cidx), Disabled.getInstance());
                        }
                        model.getAttribute().addCalls(new Integer[]{cidxclass.get(cidx).getIndex(), cidx, (int) lqn.callpair.get(cidx, 1), (int) lqn.callpair.get(cidx, 2)});
                        cidxclass.get(cidx).setCompletes(false);
                        cidxclass.get(cidx).setAttribute(new Integer[]{LayeredNetworkElement.CALL, cidx});
                        double minRespT = 0;
                        if (lqn.actsof.get(idx) != null) {
                            for (int tidx_act : lqn.actsof.get(idx)) {
                                minRespT = minRespT + lqn.hostdem_mean.get(tidx_act); // upper bound, uses all activities not just the ones reachable by this entry
                            }
                        }
                        for (int m = 1; m <= nreplicas; m++) {
                            if (minRespT == 0) {
                                serverStation.get(m).setService(cidxclass.get(cidx), Immediate.getInstance());
                            } else {
                                serverStation.get(m).setService(cidxclass.get(cidx), Exp.fitMean(minRespT));
                            }
                        }
                    }

                    if (callmean.get(cidx) != nreplicas) {
                        if (lqn.calltype.get(cidx) == CallType.SYNC) {
                            cidxauxclass.put(cidx, new ClosedClass(model, lqn.callhashnames.get(cidx) + ".Aux", 0, clientDelay));
                            cidxauxclass.get(cidx).setCompletes(false);
                            cidxauxclass.get(cidx).setAttribute(new Integer[]{LayeredNetworkElement.CALL, cidx});
                            assert clientDelay != null; // safety check not in MATLAB
                            clientDelay.setService(cidxauxclass.get(cidx), Immediate.getInstance());
                            for (int m = 1; m <= nreplicas; m++) {
                                serverStation.get(m).setService(cidxauxclass.get(cidx), Disabled.getInstance());
                            }
                        }
                    }
                }
            }
        }

        RoutingMatrix P = model.initRoutingMatrix();
        if (sourceStation != null) {
            for (int o = 1; o <= openClassesAssignedLine[0]; o++) {
                int oidx = (int) openClasses.get(o, 1) - 1;
                double p = 1.0 / openClasses.get(o, 2);
                for (int m = 1; m <= nreplicas; m++) {
                    P.addConnection(model.getClasses().get(oidx), model.getClasses().get(oidx), sourceStation, serverStation.get(m), 1.0 / (double) nreplicas);
                    for (int n = 1; n <= nreplicas; n++) {
                        P.addConnection(model.getClasses().get(oidx), model.getClasses().get(oidx), serverStation.get(m), serverStation.get(n), (1.0 - p) / (double) nreplicas);
                    }
                    P.addConnection(model.getClasses().get(oidx), model.getClasses().get(oidx), serverStation.get(m), sinkStation, p);
                }
                int cidx = (int) openClasses.get(o, 3); // 3 = source
                if (!cell_arvproc_classes_updmap.containsKey(idx)) {
                    cell_arvproc_classes_updmap.put(idx, new ArrayList<>());
                }
                cell_arvproc_classes_updmap.get(idx).add(new Integer[]{idx, cidx, model.getNodeIndex(sourceStation) + 1, oidx + 1});
                for (int m = 1; m <= nreplicas; m++) {
                    if (!cell_call_classes_updmap.containsKey(idx)) {
                        cell_call_classes_updmap.put(idx, new ArrayList<>());
                    }
                    cell_call_classes_updmap.get(idx).add(new Integer[]{idx, cidx, model.getNodeIndex(serverStation.get(m)) + 1, oidx + 1});
                }
            }

            // Setup routing for entry-level open arrivals
            for (int e = 1; e <= entryOpenClassesAssignedLine[0]; e++) {
                int eoidx = (int) entryOpenClasses.get(e, 1) - 1; // class index (0-based)
                for (int m = 1; m <= nreplicas; m++) {
                    // Route: Source → ServerStation → Sink
                    P.addConnection(model.getClasses().get(eoidx), model.getClasses().get(eoidx),
                        sourceStation, serverStation.get(m), 1.0 / (double) nreplicas);
                    P.addConnection(model.getClasses().get(eoidx), model.getClasses().get(eoidx),
                        serverStation.get(m), sinkStation, 1.0);
                }
            }
        }

        int atClient = 1;
        int atServer = 2;
        int atCache = 3;

        // Create final copies for inner class access
        final boolean finalIsCacheLayer = iscachelayer;
        final boolean finalIsFunctionLayer = isfunctionlayer;
        final Cache finalCacheNode = cacheNode;

        class InnerRecurActGraph {

            recurActGraphReturnType recurActGraph(RoutingMatrix P, int tidx_caller, int aidx, JobClass curClass, int jobPos, Source sourceStation, Delay clientDelay, Sink sinkStation, Join joinNode, Fork forkNode, Stack<JobClass> forkClassStack) {
                jobPosKey.set(0, aidx, jobPos);
                curClassKey.put(aidx, curClass);
                List<Integer> nextaidxs = new ArrayList<>();
                for (int i = 1; i < lqn.graph.getNumCols(); i++) {
                    if (lqn.graph.isAssigned(aidx, i)) {
                        nextaidxs.add(i);
                    }
                }


                Matrix isNextPrecFork = new Matrix(1, lqn.nidx + 1, lqn.nidx);
                if (!nextaidxs.isEmpty()) {
                    isNextPrecFork.set(0, aidx, 0);
                    for (int i : nextaidxs) {
                        if (isPostAndAct.isAssigned(0, i)) {
                            isNextPrecFork.set(0, aidx, 1);
                            break;
                        }
                    }
                }
                if (!nextaidxs.isEmpty()) {
                    for (int nextaidx : nextaidxs) {
                        boolean isLoop = lqn.graph.get(aidx, nextaidx) != lqn.dag.get(aidx, nextaidx);
                        // in the activity graph, the following if is entered only
                        // by an edge that is the return from a LOOP activity
                        if (lqn.parent.get(0, aidx) != lqn.parent.get(0, nextaidx)) { // if different parent task
                            int cidx = 0;
                            for (int i = 1; i <= lqn.ncalls; i++) { // find the call index
                                if (lqn.callpair.get(i, 1) == aidx && lqn.callpair.get(i, 2) == nextaidx) {
                                    cidx = i;
                                    break;
                                }
                            }
                            if (lqn.calltype.get(cidx) == CallType.ASYNC) {
                                // add only if the target is serverStation
                                int serverIdx = (int) lqn.callpair.get(cidx, 2);
                                if (serverIdx >= 0 && lqn.parent.get(serverIdx) == idx) {
                                    if (sourceStation == null) { // we need to add source and sink to the model
                                        model.getAttribute().setSourceIdx(model.getNumberOfNodes() + 1);
                                        sourceStation = new jline.lang.nodes.Source(model, "Source");
                                        sinkStation = new Sink(model, "Sink");
                                    }
                                    JobClass cidxClass = new OpenClass(model, lqn.callhashnames.get(cidx), 0);
                                    cidxclass.put(cidx, cidxClass);
                                    sourceStation.setArrival(cidxClass, Immediate.getInstance());
                                    clientDelay.setService(cidxClass, Disabled.getInstance());
                                    for (int m = 1; m <= nreplicas; m++) {
                                        serverStation.get(m).setService(cidxClass, Immediate.getInstance());
                                    }
                                    openClassesAssignedLine[0]++;
                                    openClasses.set(openClassesAssignedLine[0], 1, cidxClass.getIndex());
                                    openClasses.set(openClassesAssignedLine[0], 2, callmean.get(cidx));
                                    openClasses.set(openClassesAssignedLine[0], 3, cidx);
                                    model.getAttribute().addCalls(new Integer[]{cidxClass.getIndex(), cidx, (int) lqn.callpair.get(cidx, 1), (int) lqn.callpair.get(cidx, 2)});
                                    cidxClass.setCompletes(false);
                                    cidxClass.setAttribute(new Integer[]{LayeredNetworkElement.CALL, cidx});
                                    double minRespT = 0;
                                    List<Integer> actsof_idx = lqn.actsof.get(idx);
                                    for (int tidx_act : actsof_idx) {
                                        minRespT += lqn.hostdem_mean.getOrDefault(tidx_act, 0.0); // upper bound, uses all activities not just the ones reachable by this entry
                                    }
                                    for (int m = 1; m <= nreplicas; m++) {
                                        serverStation.get(m).setService(cidxClass, Exp.fitMean(minRespT));
                                    }
                                }
                            } else if (lqn.calltype.get(cidx) == CallType.SYNC) {
                                // START routeSynchCall in MATLAB
                                if (jobPos == atClient) {
                                    int serverIdx = (int) lqn.callpair.get(cidx, 2);
                                    if (serverIdx >= 0 && lqn.parent.get(serverIdx) == idx) {
                                        if (callmean.get(cidx) < nreplicas) {
                                            P.addConnection(curClass, cidxauxclass.get(cidx), clientDelay, clientDelay, 1.0 - callmean.get(cidx));
                                            for (int m = 1; m <= nreplicas; m++) {
                                                P.addConnection(curClass, cidxclass.get(cidx), clientDelay, serverStation.get(m), callmean.get(cidx) / (double) nreplicas);
                                                P.addConnection(cidxclass.get(cidx), cidxclass.get(cidx), serverStation.get(m), clientDelay, 1.0);
                                            }
                                            P.addConnection(cidxauxclass.get(cidx), cidxclass.get(cidx), clientDelay, clientDelay, 1.0);
                                        } else if (callmean.get(cidx) == nreplicas) {
                                            for (int m = 1; m <= nreplicas; m++) {
                                                P.addConnection(curClass, cidxclass.get(cidx), clientDelay, serverStation.get(m), 1.0 / (double) nreplicas);
                                                P.addConnection(cidxclass.get(cidx), cidxclass.get(cidx), serverStation.get(m), clientDelay, 1.0);
                                            }
                                        } else { // callmean.get(cidx) > nreplicas
                                            for (int m = 1; m <= nreplicas; m++) {
                                                P.addConnection(curClass, cidxclass.get(cidx), clientDelay, serverStation.get(m), 1.0 / (double) nreplicas);
                                                P.addConnection(cidxclass.get(cidx), cidxauxclass.get(cidx), serverStation.get(m), clientDelay, 1.0);
                                                P.addConnection(cidxauxclass.get(cidx), cidxclass.get(cidx), clientDelay, serverStation.get(m), 1.0 - 1.0 / (callmean.get(cidx) / nreplicas));
                                            }
                                            P.addConnection(cidxauxclass.get(cidx), cidxclass.get(cidx), clientDelay, clientDelay, 1.0 / callmean.get(cidx));
                                        }
                                        jobPos = atClient;
                                        clientDelay.setService(cidxclass.get(cidx), Immediate.getInstance());
                                        if (!cell_call_classes_updmap.containsKey(idx)) {
                                            cell_call_classes_updmap.put(idx, new ArrayList<>());
                                        }
                                        for (int m = 1; m <= nreplicas; m++) {
                                            serverStation.get(m).setService(cidxclass.get(cidx), callservtproc.get(cidx));
                                            cell_call_classes_updmap.get(idx).add(new Integer[]{idx, cidx, model.getNodeIndex(serverStation.get(m)) + 1, cidxclass.get(cidx).getIndex()});
                                        }
                                        curClass = cidxclass.get(cidx);
                                    } else { // if it is not a call to an entry of the server
                                        if (callmean.get(cidx) < nreplicas) {
                                            P.addConnection(curClass, cidxclass.get(cidx), clientDelay, clientDelay, callmean.get(cidx) / nreplicas); // the mean number of calls is now embedded in the demand
                                            P.addConnection(cidxclass.get(cidx), cidxauxclass.get(cidx), clientDelay, clientDelay, 1.0); // the mean number of calls is now embedded in the demand
                                            P.addConnection(curClass, cidxauxclass.get(cidx), clientDelay, clientDelay, 1.0 - callmean.get(cidx) / nreplicas); // the mean number of calls is now embedded in the demand
                                            curClass = cidxauxclass.get(cidx);
                                        } else if (callmean.get(cidx) == nreplicas) {
                                            P.addConnection(curClass, cidxclass.get(cidx), clientDelay, clientDelay, 1.0);
                                            curClass = cidxclass.get(cidx);
                                        } else {  // callmean.get(cidx) > nreplicas
                                            P.addConnection(curClass, cidxclass.get(cidx), clientDelay, clientDelay, 1.0); // the mean number of calls is now embedded in the demand
                                            P.addConnection(cidxclass.get(cidx), cidxauxclass.get(cidx), clientDelay, clientDelay, 1.0);  // the mean number of calls is now embedded in the demand
                                            curClass = cidxauxclass.get(cidx);
                                        }
                                        jobPos = atClient;
                                        clientDelay.setService(cidxclass.get(cidx), callservtproc.get(cidx));
                                        if (!cell_call_classes_updmap.containsKey(idx)) {
                                            cell_call_classes_updmap.put(idx, new ArrayList<>());
                                        }
                                        cell_call_classes_updmap.get(idx).add(new Integer[]{idx, cidx, 1, cidxclass.get(cidx).getIndex()});
                                    }
                                } else if (jobPos == atServer) {
                                    int serverIdx = (int) lqn.callpair.get(cidx, 2);
                                    if (serverIdx >= 0 && lqn.parent.get(serverIdx) == idx) {
                                        if (callmean.get(cidx) < nreplicas) {
                                            for (int m = 1; m <= nreplicas; m++) {
                                                P.addConnection(curClass, cidxclass.get(cidx), serverStation.get(m), clientDelay, 1.0 - callmean.get(cidx));
                                                P.addConnection(curClass, cidxclass.get(cidx), serverStation.get(m), serverStation.get(m), callmean.get(cidx));
                                                serverStation.get(m).setService(cidxclass.get(cidx), callservtproc.get(cidx));
                                            }
                                            jobPos = atClient;
                                            curClass = cidxauxclass.get(cidx);
                                        } else if (callmean.get(cidx) == nreplicas) {
                                            for (int m = 1; m <= nreplicas; m++) {
                                                P.addConnection(curClass, cidxclass.get(cidx), serverStation.get(m), serverStation.get(m), 1.0);
                                            }
                                            jobPos = atServer;
                                            curClass = cidxclass.get(cidx);
                                        } else {
                                            for (int m = 1; m <= nreplicas; m++) {
                                                P.addConnection(curClass, cidxclass.get(cidx), serverStation.get(m), serverStation.get(m), 1.0);
                                                P.addConnection(cidxclass.get(cidx), cidxclass.get(cidx), serverStation.get(m), serverStation.get(m), 1.0 - 1.0 / callmean.get(cidx));
                                                P.addConnection(cidxclass.get(cidx), cidxauxclass.get(cidx), serverStation.get(m), clientDelay, 1.0 / callmean.get(cidx));
                                            }
                                            jobPos = atClient;
                                            curClass = cidxauxclass.get(cidx);
                                        }
                                        if (!cell_call_classes_updmap.containsKey(idx)) {
                                            cell_call_classes_updmap.put(idx, new ArrayList<>());
                                        }
                                        for (int m = 1; m <= nreplicas; m++) {
                                            serverStation.get(m).setService(cidxclass.get(cidx), callservtproc.get(cidx));
                                            cell_call_classes_updmap.get(idx).add(new Integer[]{idx, cidx, model.getNodeIndex(serverStation.get(m)) + 1, cidxclass.get(cidx).getIndex()}); // check if getIndex+1
                                        }
                                    } else {
                                        // if it is not a call to an entry of the server
                                        // callmean not needed since we switched
                                        // to ResidT to model service time at client
                                        if (callmean.get(cidx) < nreplicas) {
                                            for (int m = 1; m <= nreplicas; m++) {
                                                P.addConnection(curClass, cidxclass.get(cidx), serverStation.get(m), clientDelay, 1.0);
                                            }
                                            P.addConnection(cidxclass.get(cidx), cidxauxclass.get(cidx), clientDelay, clientDelay, 1.0);
                                            curClass = cidxauxclass.get(cidx);
                                        } else if (callmean.get(cidx) == nreplicas) {
                                            for (int m = 1; m <= nreplicas; m++) {
                                                P.addConnection(curClass, cidxclass.get(cidx), serverStation.get(m), clientDelay, 1.0);
                                            }
                                            curClass = cidxclass.get(cidx);
                                        } else { // callmean.get(cidx) > nreplicas
                                            for (int m = 1; m <= nreplicas; m++) {
                                                P.addConnection(curClass, cidxclass.get(cidx), serverStation.get(m), clientDelay, 1.0);
                                            }
                                            P.addConnection(cidxclass.get(cidx), cidxauxclass.get(cidx), clientDelay, clientDelay, 1.0);
                                            curClass = cidxauxclass.get(cidx);
                                        }
                                        jobPos = atClient;
                                        clientDelay.setService(cidxclass.get(cidx), callservtproc.get(cidx));
                                        if (!cell_call_classes_updmap.containsKey(idx)) {
                                            cell_call_classes_updmap.put(idx, new ArrayList<>());
                                        }
                                        cell_call_classes_updmap.get(idx).add(new Integer[]{idx, cidx, 1, cidxclass.get(cidx).getIndex()}); // check if getIndex+1
                                    }
                                    // END routeSynchCall
                                }
                            }
                        } else {
                            boolean intersects;
                            intersects = false;
                            Set<Integer> range = new HashSet<>();
                            for (int i = 1; i <= lqn.nentries; i++) {
                                range.add(lqn.eshift + i);
                            }
                            for (int num : nextaidxs) {
                                if (range.contains(num)) {
                                    intersects = true;
                                    break;
                                }
                            }
                            if (!intersects) {
                                jobPos = (int) jobPosKey.get(0, aidx);
                                curClass = curClassKey.get(aidx);
                            } else {
                                // Find the index of nextaidx in nextaidxs
                                int index = nextaidxs.indexOf(nextaidx);

                                // Find the previous index if the current index is greater than 0
                                boolean isMember = false;
                                if (index > 0) {
                                    int previousValue = nextaidxs.get(index - 1);
                                    // Check if previousValue is in the range
                                    isMember = range.contains(previousValue);
                                }
                                if (isMember) {
                                    curClassC = curClass;
                                }
                                jobPos = atClient;
                                curClass = curClassC;
                            }

                            if (jobPos == atClient) {
                                if (ishostlayer) {
                                    if (!finalIsCacheLayer) {
                                        for (int m = 1; m <= nreplicas; m++) {
                                        if (isNextPrecFork.get(0, aidx) != 0) {
                                            P.addConnection(curClass, curClass, clientDelay, forkNode, 1.0);
                                            // Find the index of nextaidx within the subset of Post-And activities (matching MATLAB)
                                            List<Integer> postAndActivities = new ArrayList<>();
                                            for (int i = 0; i < nextaidxs.size(); i++) {
                                                if (isPostAndAct.get(0, nextaidxs.get(i)) == 1) {
                                                    postAndActivities.add(nextaidxs.get(i));
                                                }
                                            }
                                            int fIdx = postAndActivities.indexOf(nextaidx);
                                            if (fIdx >= 0) {
                                                int f = fIdx + 1; // MATLAB uses 1-based indexing
                                                forkClassStack.add(curClass);
                                                P.addConnection(curClass, curClass, forkNode, forkOutputRouter.get(f), 1.0);
                                                P.addConnection(curClass, aidxclass.get(nextaidx), forkOutputRouter.get(f), serverStation.get(m), 1.0);
                                            } else {
                                                // If nextaidx is not a post-and activity, use default routing without fork
                                                P.addConnection(curClass, aidxclass.get(nextaidx), clientDelay, serverStation.get(m), lqn.graph.get(aidx, nextaidx));
                                            }
                                        } else {
                                            if (isPreAndAct.get(0, aidx) != 0) {
                                                JobClass forkClass = (JobClass) forkClassStack.pop();
                                                P.addConnection(curClass, forkClass, clientDelay, joinNode, 1.0);
                                                P.addConnection(forkClass, aidxclass.get(nextaidx), joinNode, serverStation.get(m), 1.0);
                                            } else {
                                                P.addConnection(curClass, aidxclass.get(nextaidx), clientDelay, serverStation.get(m), lqn.graph.get(aidx, nextaidx));
                                            }
                                        }
                                        serverStation.get(m).setService(aidxclass.get(nextaidx), lqn.hostdem.get(nextaidx));
                                        // Check if the task owning this activity is a FunctionTask (has setup/delayoff)
                                        // This is different from finalIsFunctionLayer which checks all callers
                                        int nextaidxParentIdx = (int) lqn.parent.get(0, nextaidx);
                                        boolean isNextaidxParentFunctionTask = nextaidxParentIdx >= 0 && lqn.isfunction != null &&
                                                                        lqn.isfunction.getNumCols() > nextaidxParentIdx &&
                                                                        lqn.isfunction.get(0, nextaidxParentIdx) == 1;
                                        if (isNextaidxParentFunctionTask && lqn.setuptime != null && lqn.delayofftime != null) {
                                            Distribution setupTime = lqn.setuptime.get(nextaidxParentIdx);
                                            Distribution delayoffTime = lqn.delayofftime.get(nextaidxParentIdx);
                                            if (setupTime != null && delayoffTime != null) {
                                                serverStation.get(m).setDelayOff(aidxclass.get(nextaidx), setupTime, delayoffTime);
                                            }
                                        }
                                        }
                                        jobPos = atServer;
                                        curClass = aidxclass.get(nextaidx);
                                        if (!cell_servt_classes_updmap.containsKey(idx)) {
                                            cell_servt_classes_updmap.put(idx, new ArrayList<>());
                                        }
                                        cell_servt_classes_updmap.get(idx).add(new Integer[]{idx, nextaidx, 2, aidxclass.get(nextaidx).getIndex()});
                                    } else {
                                        // Cache layer routing: client -> cache -> server
                                        P.addConnection(curClass, aidxclass.get(nextaidx), clientDelay, finalCacheNode, lqn.graph.get(aidx, nextaidx));

                                        // Setup cache read item entry
                                        finalCacheNode.setReadItemEntry(aidxclass.get(nextaidx), lqn.itemproc.get(aidx), (int) lqn.nitems.get(0, aidx));
                                        
                                        // Find hit and miss activities
                                        List<Integer> hitmissaidx = new ArrayList<>();
                                        for (int i = 1; i <= lqn.nidx; i++) {
                                            if (lqn.graph.get(nextaidx, i) != 0) {
                                                hitmissaidx.add(i);
                                            }
                                        }
                                        
                                        if (hitmissaidx.size() >= 2) {
                                            int hitaidx = hitmissaidx.get(0);
                                            int missaidx = hitmissaidx.get(1);
                                            
                                            lqn.hitmissaidx = new ArrayList<>(hitmissaidx);
                                            lqn.hitaidx = hitaidx;
                                            lqn.missaidx = missaidx;
                                            
                                            // Set hit and miss classes
                                            finalCacheNode.setHitClass(aidxclass.get(nextaidx), aidxclass.get(hitaidx));
                                            finalCacheNode.setMissClass(aidxclass.get(nextaidx), aidxclass.get(missaidx));
                                        }

                                        jobPos = atCache;
                                        curClass = aidxclass.get(nextaidx);
                                    }
                                } else {
                                    if (isNextPrecFork.get(0, aidx) != 0) {
                                        P.addConnection(curClass, curClass, clientDelay, forkNode, 1.0);
                                        int f = 0;
                                        boolean foundNextaidx = false;
                                        for (int i = 0; i < nextaidxs.size(); i++) {
                                            if (isPostAndAct.get(0, nextaidxs.get(i)) == 1) {
                                                f++;
                                                if (nextaidxs.get(i) == nextaidx) {
                                                    foundNextaidx = true;
                                                    break;
                                                }
                                            }
                                        }
                                        if (foundNextaidx && f > 0 && forkOutputRouter.containsKey(f)) {
                                            forkClassStack.add(curClass);
                                            P.addConnection(curClass, curClass, forkNode, forkOutputRouter.get(f), 1.0);
                                            P.addConnection(curClass, aidxclass.get(nextaidx), forkOutputRouter.get(f), clientDelay, 1.0);
                                        } else {
                                            // If nextaidx is not a post-and activity or routing is unavailable, use default routing
                                            P.addConnection(curClass, aidxclass.get(nextaidx), clientDelay, clientDelay, lqn.graph.get(aidx, nextaidx));
                                        }
                                    } else {
                                        if (isPreAndAct.get(0, aidx) != 0) {
                                            JobClass forkClass = (JobClass) forkClassStack.pop();
                                            P.addConnection(curClass, forkClass, clientDelay, joinNode, 1.0);
                                            P.addConnection(forkClass, aidxclass.get(nextaidx), joinNode, clientDelay, 1.0);
                                        } else {
                                            P.addConnection(curClass, aidxclass.get(nextaidx), clientDelay, clientDelay, lqn.graph.get(aidx, nextaidx));
                                        }
                                    }
                                    jobPos = atClient;
                                    curClass = aidxclass.get(nextaidx);
                                    clientDelay.setService(aidxclass.get(nextaidx), servtproc.get(nextaidx));
                                    if (!cell_thinkt_classes_updmap.containsKey(idx)) {
                                        cell_thinkt_classes_updmap.put(idx, new ArrayList<>());
                                    }
                                    cell_thinkt_classes_updmap.get(idx).add(new Integer[]{idx, nextaidx, 1, aidxclass.get(nextaidx).getIndex()});
                                }
                            } else if (jobPos == atServer || jobPos == atCache) {
                                if (ishostlayer) {
                                    if (finalIsCacheLayer) {
                                        // Cache layer routing: cache -> server
                                        curClass = aidxclass.get(nextaidx);
                                        for (int m = 1; m <= nreplicas; m++) {
                                            if (isNextPrecFork.get(0, aidx) != 0) {
                                                P.addConnection(curClass, curClass, finalCacheNode, forkNode, 1.0);
                                                int f = 0;
                                                boolean foundNextaidx = false;
                                                for (int i = 0; i < nextaidxs.size(); i++) {
                                                    if (isPostAndAct.get(0, nextaidxs.get(i)) == 1) {
                                                        f++;
                                                        if (nextaidxs.get(i) == nextaidx) {
                                                            foundNextaidx = true;
                                                            break;
                                                        }
                                                    }
                                                }
                                                if (foundNextaidx && f > 0 && forkOutputRouter.containsKey(f)) {
                                                    forkClassStack.add(curClass);
                                                    P.addConnection(curClass, curClass, forkNode, forkOutputRouter.get(f), 1.0);
                                                    P.addConnection(curClass, aidxclass.get(nextaidx), forkOutputRouter.get(f), serverStation.get(m), 1.0);
                                                } else {
                                                    // If nextaidx is not a post-and activity or routing is unavailable, use default routing
                                                    P.addConnection(curClass, aidxclass.get(nextaidx), finalCacheNode, serverStation.get(m), lqn.graph.get(aidx, nextaidx));
                                                }
                                            } else {
                                                if (isPreAndAct.get(0, aidx) != 0) {
                                                    JobClass forkClass = (JobClass) forkClassStack.pop();
                                                    P.addConnection(curClass, forkClass, finalCacheNode, joinNode, 1.0);
                                                    P.addConnection(forkClass, aidxclass.get(nextaidx), joinNode, serverStation.get(m), 1.0);
                                                } else {
                                                    P.addConnection(curClass, aidxclass.get(nextaidx), finalCacheNode, serverStation.get(m), lqn.graph.get(aidx, nextaidx));
                                                }
                                            }
                                            Distribution baseService = lqn.hostdem.get(nextaidx);
                                            // Check if the task owning this activity is a FunctionTask (has setup/delayoff)
                                            int nextaidxParentId = (int) lqn.parent.get(0, nextaidx);
                                            boolean isNextaidxParentFunc = nextaidxParentId >= 0 && lqn.isfunction != null &&
                                                                            lqn.isfunction.getNumCols() > nextaidxParentId &&
                                                                            lqn.isfunction.get(0, nextaidxParentId) == 1;
                                            if (isNextaidxParentFunc && lqn.setuptime != null && lqn.delayofftime != null) {
                                                Distribution setupTime = lqn.setuptime.get(nextaidxParentId);
                                                Distribution delayoffTime = lqn.delayofftime.get(nextaidxParentId);
                                                if (setupTime != null && delayoffTime != null) {
                                                    // For function tasks, the effective service time includes setup and delay-off times
                                                    // Create a combined distribution that represents setup + service + delayoff
                                                    double totalMean = baseService.getMean() + setupTime.getMean() + delayoffTime.getMean();
                                                    // For now, use exponential approximation with combined mean
                                                    Distribution effectiveService = new jline.lang.processes.Exp(1.0 / totalMean);
                                                    serverStation.get(m).setService(aidxclass.get(nextaidx), effectiveService);
                                                    serverStation.get(m).setDelayOff(aidxclass.get(nextaidx), setupTime, delayoffTime);
                                                } else {
                                                    serverStation.get(m).setService(aidxclass.get(nextaidx), baseService);
                                                }
                                            } else {
                                                serverStation.get(m).setService(aidxclass.get(nextaidx), baseService);
                                            }
                                        }
                                    } else {
                                        for (int m = 1; m <= nreplicas; m++) {
                                        if (isNextPrecFork.get(0, aidx) != 0) {
                                            P.addConnection(curClass, curClass, serverStation.get(m), forkNode, 1.0);
                                            // Find the index of nextaidx within the subset of Post-And activities (matching MATLAB)
                                            List<Integer> postAndActivities = new ArrayList<>();
                                            for (int i = 0; i < nextaidxs.size(); i++) {
                                                if (isPostAndAct.get(0, nextaidxs.get(i)) == 1) {
                                                    postAndActivities.add(nextaidxs.get(i));
                                                }
                                            }
                                            int fIdx = postAndActivities.indexOf(nextaidx);
                                            if (fIdx >= 0) {
                                                int f = fIdx + 1; // MATLAB uses 1-based indexing
                                                forkClassStack.add(curClass);
                                                P.addConnection(curClass, curClass, forkNode, forkOutputRouter.get(f), 1.0);
                                                P.addConnection(curClass, aidxclass.get(nextaidx), forkOutputRouter.get(f), serverStation.get(m), 1.0);
                                            } else {
                                                // If nextaidx is not a post-and activity, use default routing without fork
                                                P.addConnection(curClass, aidxclass.get(nextaidx), serverStation.get(m), serverStation.get(m), lqn.graph.get(aidx, nextaidx));
                                            }
                                        } else {
                                            if (isPreAndAct.get(0, aidx) != 0) {
                                                JobClass forkClass = (JobClass) forkClassStack.pop();
                                                P.addConnection(curClass, forkClass, serverStation.get(m), joinNode, 1.0);
                                                P.addConnection(forkClass, aidxclass.get(nextaidx), joinNode, serverStation.get(m), 1.0);
                                            } else {
                                                P.addConnection(curClass, aidxclass.get(nextaidx), serverStation.get(m), serverStation.get(m), lqn.graph.get(aidx, nextaidx));
                                            }
                                        }
                                        serverStation.get(m).setService(aidxclass.get(nextaidx), lqn.hostdem.get(nextaidx));
                                        // Check if the task owning this activity is a FunctionTask (has setup/delayoff)
                                        int nextaidxParent = (int) lqn.parent.get(0, nextaidx);
                                        boolean isNextaidxParentFunction = nextaidxParent >= 0 && lqn.isfunction != null &&
                                                                        lqn.isfunction.getNumCols() > nextaidxParent &&
                                                                        lqn.isfunction.get(0, nextaidxParent) == 1;
                                        if (isNextaidxParentFunction && lqn.setuptime != null && lqn.delayofftime != null) {
                                            Distribution setupTime = lqn.setuptime.get(nextaidxParent);
                                            Distribution delayoffTime = lqn.delayofftime.get(nextaidxParent);
                                            if (setupTime != null && delayoffTime != null) {
                                                serverStation.get(m).setDelayOff(aidxclass.get(nextaidx), setupTime, delayoffTime);
                                            }
                                        }
                                        }
                                    }

                                    if (finalIsCacheLayer) {
                                        jobPos = atServer;
                                        curClass = aidxclass.get(nextaidx);
                                        if (!cell_servt_classes_updmap.containsKey(idx)) {
                                            cell_servt_classes_updmap.put(idx, new ArrayList<>());
                                        }
                                        cell_servt_classes_updmap.get(idx).add(new Integer[]{idx, nextaidx, 2, aidxclass.get(nextaidx).getIndex()});
                                    } else {
                                        jobPos = atServer;
                                        curClass = aidxclass.get(nextaidx);
                                        if (!cell_servt_classes_updmap.containsKey(idx)) {
                                            cell_servt_classes_updmap.put(idx, new ArrayList<>());
                                        }
                                        cell_servt_classes_updmap.get(idx).add(new Integer[]{idx, nextaidx, 2, aidxclass.get(nextaidx).getIndex()});
                                    }
                                } else {
                                    for (int m = 1; m <= nreplicas; m++) {
                                        if (isNextPrecFork.get(0, aidx) != 0) {
                                            P.addConnection(curClass, curClass, serverStation.get(m), forkNode, 1.0);
                                            // Find the index of nextaidx within the subset of Post-And activities (matching MATLAB)
                                            List<Integer> postAndActivities = new ArrayList<>();
                                            for (int i = 0; i < nextaidxs.size(); i++) {
                                                if (isPostAndAct.get(0, nextaidxs.get(i)) == 1) {
                                                    postAndActivities.add(nextaidxs.get(i));
                                                }
                                            }
                                            int fIdx = postAndActivities.indexOf(nextaidx);
                                            if (fIdx >= 0) {
                                                int f = fIdx + 1; // MATLAB uses 1-based indexing
                                                forkClassStack.add(curClass);
                                                P.addConnection(curClass, curClass, forkNode, forkOutputRouter.get(f), 1.0);
                                                P.addConnection(curClass, aidxclass.get(nextaidx), forkOutputRouter.get(f), clientDelay, 1.0);
                                            } else {
                                                // If nextaidx is not a post-and activity, use default routing without fork
                                                P.addConnection(curClass, aidxclass.get(nextaidx), serverStation.get(m), clientDelay, lqn.graph.get(aidx, nextaidx));
                                            }
                                        } else {
                                            if (isPreAndAct.get(0, aidx) != 0) {
                                                JobClass forkClass = (JobClass) forkClassStack.pop();
                                                P.addConnection(curClass, forkClass, serverStation.get(m), joinNode, 1.0);
                                                P.addConnection(forkClass, aidxclass.get(nextaidx), joinNode, clientDelay, 1.0);
                                            } else {
                                                P.addConnection(curClass, aidxclass.get(nextaidx), serverStation.get(m), clientDelay, lqn.graph.get(aidx, nextaidx));
                                            }
                                        }
                                        jobPos = atClient;
                                        curClass = aidxclass.get(nextaidx);
                                        clientDelay.setService(aidxclass.get(nextaidx), servtproc.get(nextaidx));
                                        if (!cell_thinkt_classes_updmap.containsKey(idx)) {
                                            cell_thinkt_classes_updmap.put(idx, new ArrayList<>());
                                        }
                                        cell_thinkt_classes_updmap.get(idx).add(new Integer[]{idx, nextaidx, 1, aidxclass.get(nextaidx).getIndex()});
                                    }
                                }
                            }
                            if (aidx != nextaidx && !isLoop) {
                                recurActGraphReturnType returnType = recurActGraph(P, tidx_caller, nextaidx, curClass, jobPos, sourceStation, clientDelay, sinkStation, joinNode, forkNode, forkClassStack);
                                P = returnType.P;
                                curClass = returnType.curClass;
                                jobPos = returnType.jobPos;

                                if (jobPos == atClient) {
                                    P.addConnection(curClass, aidxclass.get(tidx_caller), clientDelay, clientDelay, 1.0);
                                    if (!curClass.getName().endsWith(".Aux")) {
                                        curClass.setCompletes(true);
                                    }
                                } else {
                                    for (int m = 1; m <= nreplicas; m++) {
                                        P.addConnection(curClass, aidxclass.get(tidx_caller), serverStation.get(m), clientDelay, 1.0);
                                    }
                                    if (!curClass.getName().endsWith(".Aux")) {
                                        curClass.setCompletes(true);
                                    }
                                }
                            }
                        }
                    }
                }
                return new

                        recurActGraphReturnType(curClass, jobPos, P);
            }
        }

        int jobPos = atClient; // start at client
        // second pass: setup the routing out of entries
        for (int tidx_caller : callers) {
            // Per-caller: check if this caller is a sync caller to entries of idx
            // For host layers, idx is a host index with no entries, so this is always false
            boolean isSyncCallerToTargetEntries = false;
            if (!ishostlayer && lqn.entriesof.get(idx) != null) {
                for (int entry_idx : lqn.entriesof.get(idx)) {
                    if (lqn.issynccaller.get(tidx_caller, entry_idx) != 0) {
                        isSyncCallerToTargetEntries = true;
                        break;
                    }
                }
            }
            // Per-caller: recompute hasDirectCallers for host layers (same as first pass)
            boolean hasDirectCallers2 = false;
            if (ishostlayer) {
                if (lqn.isref.get(tidx_caller) != 0) {
                    hasDirectCallers2 = true;
                } else {
                    for (int eidx : lqn.entriesof.get(tidx_caller)) {
                        for (int row = 1; row <= lqn.ntasks; row++) {
                            int tidx_row = row + lqn.tshift;
                            if (lqn.issynccaller.get(tidx_row, eidx) != 0 || lqn.isasynccaller.get(tidx_row, eidx) != 0) {
                                hasDirectCallers2 = true;
                                break;
                            }
                        }
                        if (hasDirectCallers2) break;
                        if (lqn.arrival != null && lqn.arrival.containsKey(eidx) && lqn.arrival.get(eidx) != null) {
                            hasDirectCallers2 = true;
                            break;
                        }
                    }
                }
            }
            if ((ishostlayer && hasDirectCallers2) || isSyncCallerToTargetEntries) { // if it is only an asynch caller the closed classes are not needed
                int ncaller_entries = lqn.entriesof.get(tidx_caller).size();
                for (int eidx : lqn.entriesof.get(tidx_caller)) {
                    JobClass aidxClass_eidx = aidxclass.get(eidx);
                    JobClass aidxClass_tidx_caller = aidxclass.get(tidx_caller);
                    P.addConnection(aidxClass_tidx_caller, aidxClass_eidx, clientDelay, clientDelay, 1.0 / (double) ncaller_entries);
                    if (ncaller_entries > 1) {
                        // at successive iterations make sure to replace this with throughput ratio
                        if (!cell_route_prob_updmap.containsKey(idx)) {
                            cell_route_prob_updmap.put(idx, new ArrayList<>());
                        }
                        cell_route_prob_updmap.get(idx).add(new Integer[]{idx, tidx_caller, eidx, 1, 1, aidxClass_tidx_caller.getIndex(), aidxClass_eidx.getIndex()});
                    }
                    P = new InnerRecurActGraph().recurActGraph(P, tidx_caller, eidx, aidxClass_eidx, jobPos, sourceStation, clientDelay, sinkStation, joinNode, forkNode, forkClassStack).P;
                }
            }
        }
        model.link(P);
        temp_ensemble.add(model);
    }

    public void construct() {
        // mark down to ignore unreachable disconnected components
        this.ignore = new Matrix(lqn.nidx + 1, 1);
        Set<Set<Integer>> wccs = weaklyConnect(lqn.graph, null);
        // Remove the connected component consisting of a singleton {0} introduced by 0-padding
        wccs.removeIf(component -> component.size() == 1 && component.contains(0));
        if (wccs.size() > 1) {
            // The model has disconnected submodels
            for (Set<Integer> component : wccs) {
                boolean hasREF = false;
                // Find if the wcc has a reference task
                for (Integer component_idx : component) {
                    if (component_idx < lqn.eshift) {
                        if (lqn.sched.get(component_idx) == SchedStrategy.REF) {
                            hasREF = true;
                            break;
                        }
                    }
                }
                if (!hasREF) {
                    for (Integer component_idx : component) {
                        this.ignore.set(component_idx, 0, 1.0); // true
                    }
                }
            }
        }

        // initialize internal data structures
        this.nlayers = 0;
        this.entrycdfrespt = new Matrix(lqn.nentries, 1);
        this.hasconverged = false;

        // initialize CDF and process maps for moment3 method
        this.servtcdf = new HashMap<Integer, Matrix>();
        this.callservtcdf = new HashMap<Integer, Matrix>();
        this.entryproc = new HashMap<Integer, APH>();

        // initialize svc and think times
        this.servtproc = new HashMap<Integer, Distribution>();
        this.tputproc = new HashMap<Integer, Distribution>();
        this.servtproc.putAll(lqn.hostdem);
        this.thinkproc = new HashMap<Integer, Distribution>();
        this.thinkproc.putAll(lqn.think);
        this.callservtproc = new HashMap<Integer, Distribution>();
        for (int cidx = 1; cidx <= lqn.ncalls; cidx++) {
            int serverIdx = (int) lqn.callpair.get(cidx, 2);
            if (serverIdx >= 0 && lqn.hostdem.containsKey(serverIdx)) {
                callservtproc.put(cidx, lqn.hostdem.get(serverIdx));
            }
        }

        // perform layering
        this.njobs = new Matrix(lqn.tshift + lqn.ntasks + 1, lqn.tshift + lqn.ntasks + 1, lqn.nidx * lqn.nidx);
        this.idxhash = new ArrayList<>();
        idxhash.add(Double.NaN);
        buildLayers();
        this.solvers = new NetworkSolver[nlayers + 1];
        this.njobsorig = new Matrix(this.njobs);

        // initialize data structures for interlock correction
        this.ptaskcallers = new Matrix(lqn.nhosts + lqn.ntasks + 1, lqn.nhosts + lqn.ntasks + 1, lqn.nidx * lqn.nidx);
        this.ptaskcallers_step = new HashMap<Integer, Matrix>(nlayers);
        for (int i = 1; i <= this.nlayers; i++) {
            this.ptaskcallers_step.put(i, new Matrix(lqn.nhosts + lqn.ntasks, lqn.nhosts + lqn.ntasks, lqn.nidx * lqn.nidx));
        }

        // layering generates update maps that we use here to cache the elements that need reset
        this.routereset = new ArrayList<>();
        for (int i = 1; i < route_prob_updmap.getNumRows(); i++) {
            int idx = (int) route_prob_updmap.get(i, 1);
            if (idx >= 0 && idx < idxhash.size()) {
                int buffer = idxhash.get(idx).intValue();
                if (!routereset.contains(buffer)) { // unique
                    routereset.add(buffer);
                }
            }
        }

        this.svcreset = new ArrayList<>();
        for (int i = 1; i < thinkt_classes_updmap.getNumRows(); i++) {
            int idx = (int) thinkt_classes_updmap.get(i, 1);
            if (idx >= 0 && idx < idxhash.size()) {
                int buffer = idxhash.get(idx).intValue();
                if (!svcreset.contains(buffer)) { // unique
                    svcreset.add(buffer);
                }
            }
        }
        for (int i = 1; i < call_classes_updmap.getNumRows(); i++) {
            int idx = (int) call_classes_updmap.get(i, 1);
            if (idx >= 0 && idx < idxhash.size()) {
                int buffer = idxhash.get(idx).intValue();
                if (!svcreset.contains(buffer)) { // unique
                    svcreset.add(buffer);
                }
            }
        }
        Collections.sort(svcreset);
    }

    public boolean converged(int it) {
        /* Apply convergence test to SolverLN iterations. As the solver keeps iterating, this method maintains a
         * moving avg of the recent results based on which it averages across the layer the maximum queue-length
         * error. Convergence is tested by resetting all layers (to avoid caching) and doing an extra iteration.
         * If the iteration keeps fulfilling the error requirements for convergence, the solver completes.*/

        boolean bool = false;
        int iter_min = FastMath.max(2 * this.ensemble.length, (int) FastMath.ceil(this.options.iter_max / 4.0));
        int E = this.nlayers;
        Map<Integer, Map<Integer, SolverResult>> results = this.results;

        // Start moving average to help convergence

        if (false) {
            if (averagingstart != null) {
                int wnd_size = it - this.averagingstart + 1;
                double mov_avg_weight = 1.0 / (double) wnd_size;
                // assume ready state
                if (it >= iter_min) {
                    for (int e = 0; e < E; e++) {
                        results.get(results.size()).get(e).QN.add(mov_avg_weight - 1, results.get(results.size()).get(e).QN);
                        results.get(results.size()).get(e).UN.add(mov_avg_weight - 1, results.get(results.size()).get(e).UN);
                        results.get(results.size()).get(e).RN.add(mov_avg_weight - 1, results.get(results.size()).get(e).RN);
                        results.get(results.size()).get(e).TN.add(mov_avg_weight - 1, results.get(results.size()).get(e).TN);
                        results.get(results.size()).get(e).AN.add(mov_avg_weight - 1, results.get(results.size()).get(e).AN);
                        results.get(results.size()).get(e).WN.add(mov_avg_weight - 1, results.get(results.size()).get(e).WN);

                        for (int k = 1; k < wnd_size; k++) {
                            results.get(results.size()).get(e).QN.add(mov_avg_weight, results.get(results.size() - k).get(e).QN);
                            results.get(results.size()).get(e).UN.add(mov_avg_weight, results.get(results.size() - k).get(e).UN);
                            results.get(results.size()).get(e).RN.add(mov_avg_weight, results.get(results.size() - k).get(e).RN);
                            results.get(results.size()).get(e).TN.add(mov_avg_weight, results.get(results.size() - k).get(e).TN);
                            results.get(results.size()).get(e).AN.add(mov_avg_weight, results.get(results.size() - k).get(e).AN);
                            results.get(results.size()).get(e).WN.add(mov_avg_weight, results.get(results.size() - k).get(e).WN);

                        }
                    }
                }
            }
        } else {
            int wnd_size = Integer.max(5, (int) FastMath.ceil(iter_min / 5.0));
            double mov_avg_weight = 1.0 / (double) wnd_size;
            results = this.results;
            // assume ready state
            if (it >= iter_min) {
                for (int e = 0; e < E; e++) {
                    Matrix filler = results.get(results.size()).get(e).QN.copy();
                    filler.fill(mov_avg_weight);
                    results.get(results.size()).get(e).QN = results.get(results.size()).get(e).QN.elementMult(filler, null);
                    results.get(results.size()).get(e).UN = results.get(results.size()).get(e).UN.elementMult(filler, null);
                    results.get(results.size()).get(e).RN = results.get(results.size()).get(e).RN.elementMult(filler, null);
                    results.get(results.size()).get(e).TN = results.get(results.size()).get(e).TN.elementMult(filler, null);
                    results.get(results.size()).get(e).AN = results.get(results.size()).get(e).AN.elementMult(filler, null);
                    results.get(results.size()).get(e).WN = results.get(results.size()).get(e).WN.elementMult(filler, null);
                    for (int k = 1; k < wnd_size; k++) {
                        results.get(results.size()).get(e).QN = results.get(results.size()).get(e).QN.add(1, results.get(results.size() - k).get(e).QN.elementMult(filler, null));
                        results.get(results.size()).get(e).UN = results.get(results.size()).get(e).UN.add(1, results.get(results.size() - k).get(e).UN.elementMult(filler, null));
                        results.get(results.size()).get(e).RN = results.get(results.size()).get(e).RN.add(1, results.get(results.size() - k).get(e).RN.elementMult(filler, null));
                        results.get(results.size()).get(e).TN = results.get(results.size()).get(e).TN.add(1, results.get(results.size() - k).get(e).TN.elementMult(filler, null));
                        results.get(results.size()).get(e).AN = results.get(results.size()).get(e).AN.add(1, results.get(results.size() - k).get(e).AN.elementMult(filler, null));
                        results.get(results.size()).get(e).WN = results.get(results.size()).get(e).WN.add(1, results.get(results.size() - k).get(e).WN.elementMult(filler, null));

                    }
                }
            }
        }

        this.results = results;

        // Take as error metric the max qlen-error averaged across layers
        if (it > 1) {
            if (it == 2) { // initialize
                this.maxitererr = new ArrayList<Double>();
                this.maxitererr.add(0.0);
                this.maxitererr.add(0.0);
            }

            this.maxitererr.add(0.0);

            for (int e = 0; e < E; e++) {
                Matrix metric = results.get(results.size()).get(e).QN;
                Matrix metric_1 = results.get(results.size() - 1).get(e).QN;
                double N = this.ensemble[e].getNumberOfJobs().elementSum();
                if (N > 0) {
                    double IterErr;
                    try {
                        Matrix difference01 = metric.sub(1, metric_1);
                        difference01.removeNaN();
                        difference01.absEq();
                        IterErr = difference01.elementMax() / N;
                    } catch (Exception exception) {
                        IterErr = 0.0;
                    }
                    this.maxitererr.set(it, this.maxitererr.get(it) + IterErr);
                }
//                if (this.options.verbose != VerboseLevel.SILENT) {
//                    if (this.solvers[e].options.verbose != VerboseLevel.SILENT) {
//                        String msg = String.format("\bQLen change: %.5f.\n", this.maxitererr.get(it) / E);
//                        System.out.print(msg);
//                    }
//                }
                if (it == iter_min) {
                    if (this.options.verbose != VerboseLevel.SILENT) {
                        System.out.print(". Started averaging to aid convergence.");
                    }
                    this.averagingstart = it;
                }
            }

            // Update relaxation factor for adaptive/auto modes
            String relaxMode = this.options.config.relax;
            if (relaxMode == null) relaxMode = "none";
            if (relaxMode.equalsIgnoreCase("adaptive") || relaxMode.equalsIgnoreCase("auto")) {
                // Track error history
                this.relax_err_history.add(this.maxitererr.get(it));
                int wnd = this.options.config.relax_history;
                while (this.relax_err_history.size() > wnd) {
                    this.relax_err_history.remove(0);
                }

                if (this.relax_err_history.size() >= 3) {
                    // Detect oscillation by counting sign changes in error differences
                    int signChanges = 0;
                    for (int i = 0; i < this.relax_err_history.size() - 2; i++) {
                        double diff1 = this.relax_err_history.get(i + 1) - this.relax_err_history.get(i);
                        double diff2 = this.relax_err_history.get(i + 2) - this.relax_err_history.get(i + 1);
                        if (diff1 * diff2 < 0) {
                            signChanges++;
                        }
                    }

                    int numDiffs = this.relax_err_history.size() - 1;
                    if (relaxMode.equalsIgnoreCase("auto") && this.relax_omega == 1.0) {
                        // For 'auto' mode: enable relaxation when oscillation detected
                        if (signChanges >= numDiffs * 0.5) {
                            this.relax_omega = this.options.config.relax_factor;
                            // Debug output removed
                            if (this.options.verbose != VerboseLevel.SILENT) {
                                // System.out.printf(" [enabling relaxation, omega=%.2f]", this.relax_omega);
                            }
                        }
                    } else if (relaxMode.equalsIgnoreCase("adaptive")) {
                        // For 'adaptive' mode: adjust omega based on error trajectory
                        if (signChanges >= numDiffs * 0.5) {
                            // Oscillating - reduce omega
                            this.relax_omega = FastMath.max(this.options.config.relax_min, this.relax_omega * 0.8);
                            // Debug output removed
                            if (this.options.verbose != VerboseLevel.SILENT) {
                                // System.out.printf(" [omega=%.2f]", this.relax_omega);
                            }
                        } else if (signChanges == 0 && this.maxitererr.get(it) < this.maxitererr.get(it - 1)) {
                            // Monotonically decreasing - can increase omega slightly
                            this.relax_omega = FastMath.min(1.0, this.relax_omega * 1.05);
                        }
                    }
                }
            }
        }

        // Check convergence. Do not allow to converge in less than 2 iterations.
        if (it == 0 && (this.options.verbose != VerboseLevel.SILENT)) {
            // Debug output removed
        } else if ((it > 2) && (this.maxitererr.get(it) < this.options.iter_tol) && (this.maxitererr.get(it - 1) < this.options.iter_tol) && (this.maxitererr.get(it - 1) < this.options.iter_tol)) {
            // if potential convergence has just been detected, do a hard reset of every layer to check that this is
            // really the fixed point
            if (!this.hasconverged) {
                for (int e = 0; e < E; e++) {
                    this.ensemble[e].reset(false);
                }
                if (this.options.verbose != VerboseLevel.SILENT) {
                    //if (this.solvers[this.solvers.length - 1].options.verbose != VerboseLevel.SILENT) {
                    // Debug output removed
                    //}
                }
                //If it passes the change again next time then complete
                this.hasconverged = true;
            } else {
                if (this.options.verbose != VerboseLevel.SILENT) {
                    //if (this.solvers[this.solvers.length - 1].options.verbose != VerboseLevel.SILENT) {
                    //    String msg = String.format("SolverLN completed in %d iterations.", results.get(1).size());
                    //    System.out.println(msg);
                    //} else {
                    // Debug output removed
                    //}
                }
                bool = true;
            }
        } else {
            this.hasconverged = false;
        }
        return bool;
    }

    public void finish() {
        for (int e = 0; e < this.ensemble.length; e++) {
            solvers[e].getAvg();
        }
        //this.model.ensemble = this.ensemble; // not included as this comes through Diamond inheritance
    }

    public Matrix getArvproc_classes_updmap() {
        return arvproc_classes_updmap;
    }

    public AvgTable getAvgTable() {
        return getEnsembleAvg();
    }

    public AvgTable avgTable() {
        return getAvgTable();
    }

    public AvgTable avgT() {
        return getAvgTable();
    }

    public AvgTable aT() {
        return getAvgTable();
    }

    public Matrix getCall_classes_updmap() {
        return call_classes_updmap;
    }

    public List<Network> getEnsemble() {
        List<Network> myEnsemble = new ArrayList<>();
        Collections.addAll(myEnsemble, ensemble);
        return myEnsemble;
    }

    //GC
    @Override
    public AvgTable getEnsembleAvg() {
        return getEnsembleAvgInternal();
    }

    private AvgTable getEnsembleAvgInternal() {
        // Check if solver was properly constructed (may have returned early due to unsupported features)
        if (this.ensemble == null || this.ensemble.length == 0) {
            return null;
        }
        try {
            this.iterate();
        } catch (Exception e) {
            line_warning(mfilename(new Object(){}.getClass().getEnclosingMethod()),
                "SolverLN iteration failed: %s at %s", e.getMessage(),
                e.getStackTrace().length > 0 ? e.getStackTrace()[0].toString() : "unknown");
            e.printStackTrace();
            return null;
        }
        // NOTE: TestSolverLN, TestSolverLN2, TestSolverLN3, had problems here due to
        // different values returned by getAvg() in MATLAB and LINE on these examples

        // Check if results were populated during iteration
        if (this.results == null || this.results.isEmpty()) {
            line_warning(mfilename(new Object(){}.getClass().getEnclosingMethod()),
                "SolverLN iteration did not produce results. Returning null AvgTable.");
            return null;
        }

        Matrix QN = new Matrix(1, this.lqn.nidx + 1, this.lqn.nidx);
        for (int i = 1; i <= this.lqn.nidx; i++)
            QN.set(i, Double.NaN);

        Matrix UN = new Matrix(1, this.lqn.nidx + 1, this.lqn.nidx);
        for (int i = 1; i <= this.lqn.nidx; i++)
            UN.set(i, Double.NaN);

        Matrix RN = new Matrix(1, this.lqn.nidx + 1, this.lqn.nidx);
        for (int i = 1; i <= this.lqn.nidx; i++)
            RN.set(i, Double.NaN);

        Matrix TN = new Matrix(1, this.lqn.nidx + 1, this.lqn.nidx);
        for (int i = 1; i <= this.lqn.nidx; i++)
            TN.set(i, Double.NaN);

        // utilization will be first stored here
        Matrix PN = new Matrix(1, this.lqn.nidx + 1, this.lqn.nidx);
        for (int i = 1; i <= this.lqn.nidx; i++)
            PN.set(i, Double.NaN);

        // response time will be first stored here
        Matrix SN = new Matrix(1, this.lqn.nidx + 1, this.lqn.nidx);
        for (int i = 1; i <= this.lqn.nidx; i++)
            SN.set(i, Double.NaN);

        //residence time
        Matrix WN = new Matrix(1, this.lqn.nidx + 1, this.lqn.nidx);
        for (int i = 1; i <= this.lqn.nidx; i++)
            WN.set(i, Double.NaN);

        // not available yet
        Matrix AN = new Matrix(1, this.lqn.nidx + 1, this.lqn.nidx);
        for (int i = 1; i <= this.lqn.nidx; i++)
            AN.set(i, Double.NaN);

        // Track which activities have already been accumulated into task WN
        // to prevent double-counting when an activity appears in multiple layers
        boolean[] wnProcessed = new boolean[this.lqn.nidx + 1];

        int E = this.nlayers;

        for (int e = 0; e < E; e++) {
            int clientIdx = this.ensemble[e].getAttribute().getClientIdx();
            int serverIdx = this.ensemble[e].getAttribute().getServerIdx();
            int sourceIdx = this.ensemble[e].getAttribute().getSourceIdx();

            // determine processor metrics
            if (serverIdx != -1) {
            Station s = this.ensemble[e].getStations().get(serverIdx - 1);
            Queue q = (Queue) s;
            if (q.getAttribute().getIsHost()) {
                int hidx = q.getAttribute().getIdx();
                TN.set(0, hidx, 0);
                PN.set(0, hidx, 0);
                for (int c = 0; c < this.ensemble[e].getNumberOfClasses(); c++) {
                    if (this.ensemble[e].getClasses().get(c).getCompletes()) {
                        double t = 0;
                        double u = 0;
                        if (clientIdx != -1) {
                            t = FastMath.max(t, this.results.get(this.results.size()).get(e).TN.get(clientIdx - 1, c));
                        }
                        if (sourceIdx != -1) {
                            t = FastMath.max(t, this.results.get(this.results.size()).get(e).TN.get(sourceIdx - 1, c));
                        }
                        TN.set(0, hidx, TN.get(hidx) + FastMath.max(t, this.results.get(this.results.size()).get(e).TN.get(serverIdx - 1, c)));
                    }
                    int type = this.ensemble[e].getClasses().get(c).getAttribute()[0];
                    if (type == LayeredNetworkElement.ACTIVITY) {
                        int aidx = this.ensemble[e].getClasses().get(c).getAttribute()[1];
                        int tidx = (int) this.lqn.parent.get(0, aidx);
                        if (Double.isNaN(PN.get(aidx))) PN.set(0, aidx, 0);
                        if (tidx >= 0) {
                            if (Double.isNaN(PN.get(tidx))) PN.set(0, tidx, 0);
                            PN.set(0, aidx, PN.get(aidx) + this.results.get(this.results.size()).get(e).UN.get(serverIdx - 1, c));
                            PN.set(0, tidx, PN.get(tidx) + this.results.get(this.results.size()).get(e).UN.get(serverIdx - 1, c));
                        } else {
                            PN.set(0, aidx, PN.get(aidx) + this.results.get(this.results.size()).get(e).UN.get(serverIdx - 1, c));
                        }
                        PN.set(0, hidx, PN.get(hidx) + this.results.get(this.results.size()).get(e).UN.get(serverIdx - 1, c));
                    }
                }
                TN.set(0, hidx, Double.NaN); // Added for consistency with LQNS
            }
            }

            //determine remaining metrics
            // Declare variables outside switch to avoid scope issues
            int tidx, eidx, cidx, aidx;
            Station task_s, entry_s;
            Queue task_q, entry_q;

            if (serverIdx != -1) {
            for (int c = 0; c < this.ensemble[e].getNumberOfClasses(); c++) {
                int type = this.ensemble[e].getClasses().get(c).getAttribute()[0];
                switch (type) {
                    case LayeredNetworkElement.TASK:
                        tidx = this.ensemble[e].getClasses().get(c).getAttribute()[1];
                        task_s = this.ensemble[e].getStations().get(serverIdx - 1);
                        task_q = (Queue) task_s;
                        if (task_q.getAttribute().getIsHost()) {
                            if (Double.isNaN(TN.get(tidx)) && clientIdx != -1) {
                                // Read task throughput from layer results (like MATLAB)
                                // this.tput is not populated for REF tasks
                                TN.set(tidx, this.results.get(this.results.size()).get(e).TN.get(clientIdx - 1, c));
                            }
                        }
                        break;

                    case LayeredNetworkElement.ENTRY:
                        eidx = this.ensemble[e].getClasses().get(c).getAttribute()[1];
//                        tidx = (int) this.lqn.parent.get(eidx);  //unused parameter
                        // For phase-2 models, use residt (caller's view with overtaking)
                        // Otherwise use servt (total service time = response time)
                        // Use 1D access (idx) since servt/residt can be row or column vectors
                        if (this.hasPhase2 && this.servt_ph2 != null && this.servt_ph2.get(eidx - 1) > GlobalConstants.FineTol) {
                            SN.set(eidx, this.residt.get(eidx - 1));  // Phase-1 + overtaking correction
                        } else {
                            SN.set(eidx, this.servt.get(eidx - 1));
                        }
                        entry_s = this.ensemble[e].getStations().get(serverIdx - 1);
                        entry_q = (Queue) entry_s;
                        if (entry_q.getAttribute().getIsHost()) {
                            if (Double.isNaN(TN.get(eidx)) && clientIdx != -1) {
                                // store the result in th eprocessor model
                                TN.set(eidx, this.results.get(this.results.size()).get(e).TN.get(clientIdx - 1, c));
                            }
                        }
                        break;

                    case LayeredNetworkElement.CALL:
                        cidx = this.ensemble[e].getClasses().get(c).getAttribute()[1];
                        aidx = (int) this.lqn.callpair.get(cidx, 1);
                        // Only sync calls contribute to caller's response time
                        if (this.lqn.calltype.get(cidx) == CallType.SYNC) {
                            SN.set(aidx, SN.get(aidx) + this.results.get(this.results.size()).get(e).RN.get(serverIdx - 1, c) * this.lqn.callproc_mean.getOrDefault(cidx, 1.0));
                        }
                        if (Double.isNaN(QN.get(aidx))) {
                            QN.set(aidx, 0);
                        }
                        QN.set(aidx, QN.get(aidx) + this.results.get(this.results.size()).get(e).QN.get(serverIdx - 1, c));
                        break;

                    case LayeredNetworkElement.ACTIVITY:
                        aidx = this.ensemble[e].getClasses().get(c).getAttribute()[1];

                        // For activities with Immediate host demand, set RespT to 0 explicitly.
                        // These activities have zero service time by definition.
                        // The NC solver may return non-zero RN due to how it handles Immediate
                        // service, so we override with 0 and skip RN accumulation.
                        Distribution hostDemand = this.lqn.hostdem.get(aidx);
                        boolean skipRNAccumulation = (hostDemand instanceof Immediate);

                        // Also skip RN accumulation if this class has Disabled service at the server station.
                        // This happens for RefTask activities in non-host layers where they are callers,
                        // not the ones being served. Accumulating RN from these layers would incorrectly
                        // add non-zero response times to activities that have Immediate host demand.
                        // Note: We only skip RN/SN, not other metrics like QN, TN, WN.
                        Station serverStation_check = this.ensemble[e].getStations().get(serverIdx - 1);
                        JobClass activityClass_check = this.ensemble[e].getClasses().get(c);
                        if (serverStation_check instanceof ServiceStation) {
                            Distribution serverService_check = ((ServiceStation) serverStation_check).getServiceProcess(activityClass_check);
                            if (serverService_check instanceof Disabled) {
                                skipRNAccumulation = true;  // Skip RN/SN but continue with other metrics
                            }
                        }
                        if (skipRNAccumulation) {
                            if (Double.isNaN(SN.get(aidx))) {
                                SN.set(aidx, 0);  // Set RespT = 0 for Immediate host demand
                            }
                            if (Double.isNaN(RN.get(aidx))) {
                                RN.set(aidx, 0);  // Set RN = 0 for Immediate host demand
                            }
                            // Continue to process other metrics but skip RN accumulation below
                        }
                        tidx = (int) this.lqn.parent.get(0, aidx);
                        if (Double.isNaN(QN.get(tidx))) QN.set(tidx, 0);
                        QN.set(tidx, QN.get(tidx) + this.results.get(this.results.size()).get(e).QN.get(serverIdx - 1, c));

                        if (Double.isNaN(TN.get(aidx))) {
                            TN.set(aidx, 0);
                        }
                        if (Double.isNaN(QN.get(aidx))) {
                            QN.set(aidx, 0);
                        }

                        // For forwarding targets: propagate activity metrics to task
                        // Check if task has its own class (non-forwarding targets do)
                        boolean hasTaskClass = false;
                        Map<Integer, Integer[]> tasks = this.ensemble[e].getAttribute().getTasks();
                        if (tasks != null) {
                            for (Integer[] taskAttr : tasks.values()) {
                                if (taskAttr.length > 1 && taskAttr[1] != null && taskAttr[1].intValue() == tidx) {
                                    hasTaskClass = true;
                                    break;
                                }
                            }
                        }
                        if (!hasTaskClass) {
                            // Propagate activity throughput to task for forwarding targets
                            if (Double.isNaN(TN.get(tidx))) {
                                TN.set(tidx, 0);
                            }
                            switch (this.ensemble[e].getClasses().get(c).getJobClassType()) {
                                case CLOSED:
                                    TN.set(tidx, TN.get(tidx) + this.results.get(this.results.size()).get(e).TN.get(serverIdx - 1, c));
                                    break;
                                case OPEN:
                                    if (sourceIdx != -1) {
                                        TN.set(tidx, TN.get(tidx) + this.results.get(this.results.size()).get(e).TN.get(sourceIdx - 1, c));
                                    }
                                    break;
                                default:
                            }
                        }

                        // Find the entry this activity is bound to (matches MATLAB lines 130-153)
                        List<Integer> entriesOfTask = this.lqn.entriesof.get(tidx);
                        if (entriesOfTask != null) {
                            for (int eidx_check : entriesOfTask) {
                                // Check if activity is bound to this entry via graph
                                if (this.lqn.graph.get(eidx_check, aidx) > 0) {
                                    if (Double.isNaN(TN.get(eidx_check))) TN.set(eidx_check, 0);
                                    if (Double.isNaN(QN.get(eidx_check))) QN.set(eidx_check, 0);
                                    if (Double.isNaN(SN.get(eidx_check))) SN.set(eidx_check, 0);

                                    // Get activity throughput
                                    double actTput = 0;
                                    switch (this.ensemble[e].getClasses().get(c).getJobClassType()) {
                                        case CLOSED:
                                            actTput = this.results.get(this.results.size()).get(e).TN.get(serverIdx - 1, c);
                                            break;
                                        case OPEN:
                                            if (sourceIdx != -1) {
                                                actTput = this.results.get(this.results.size()).get(e).TN.get(sourceIdx - 1, c);
                                            }
                                            break;
                                        default:
                                    }

                                    // Only add if entry doesn't have its own class (forwarding target case)
                                    boolean hasEntryClass = false;
                                    Map<Integer, Integer[]> entries = this.ensemble[e].getAttribute().getEntries();
                                    if (entries != null) {
                                        for (Integer[] entryAttr : entries.values()) {
                                            if (entryAttr.length > 1 && entryAttr[1] != null && entryAttr[1].intValue() == eidx_check) {
                                                hasEntryClass = true;
                                                break;
                                            }
                                        }
                                    }
                                    if (!hasEntryClass) {
                                        TN.set(eidx_check, TN.get(eidx_check) + actTput);
                                        QN.set(eidx_check, QN.get(eidx_check) + this.results.get(this.results.size()).get(e).QN.get(serverIdx - 1, c));
                                        SN.set(eidx_check, SN.get(eidx_check) + this.results.get(this.results.size()).get(e).RN.get(serverIdx - 1, c));
                                    }
                                    break;
                                }
                            }
                        }

                        switch (this.ensemble[e].getClasses().get(c).getJobClassType()) {
                            case CLOSED:
                                TN.set(aidx, TN.get(aidx) + this.results.get(this.results.size()).get(e).TN.get(serverIdx - 1, c));
                                break;

                            case OPEN:
                                if (sourceIdx != -1) {
                                    TN.set(aidx, TN.get(aidx) + this.results.get(this.results.size()).get(e).TN.get(sourceIdx - 1, c));
                                }
                                break;

                            default:
                        }
                        // Skip RN accumulation for activities with Immediate host demand
                        // (skipRNAccumulation flag was set earlier in this case block)
                        if (!skipRNAccumulation) {
                            if (Double.isNaN(SN.get(aidx))) {
                                SN.set(aidx, 0);
                            }
                            SN.set(aidx, SN.get(aidx) + this.results.get(this.results.size()).get(e).RN.get(serverIdx - 1, c));

                            if (Double.isNaN(RN.get(aidx))) {
                                RN.set(aidx, 0);
                            }
                            RN.set(aidx, RN.get(aidx) + this.results.get(this.results.size()).get(e).RN.get(serverIdx - 1, c));
                        }

                        if (Double.isNaN(WN.get(aidx))) {
                            WN.set(aidx, 0);
                        }

                        if (Double.isNaN(WN.get(tidx))) {
                            WN.set(tidx, 0);
                        }
                        // Use this.residt (computed via QN/TN_ref in updateMetricsDefault)
                        // instead of layer WN to avoid fork+loop visit distortion
                        // (matches MATLAB getEnsembleAvg.m)
                        WN.set(aidx, this.residt.get(aidx - 1));
                        if (!wnProcessed[aidx]) {
                            WN.set(tidx, WN.get(tidx) + this.residt.get(aidx - 1));
                            wnProcessed[aidx] = true;
                        }
                        if (Double.isNaN(QN.get(aidx))) {
                            QN.set(aidx, 0);
                        }
                        QN.set(aidx, QN.get(aidx) + this.results.get(this.results.size()).get(e).QN.get(serverIdx - 1, c));
                        break;
                    default:
                }
            }
            }
        }
        for (int e = 1; e <= this.lqn.nentries; e++) {
            int eidx = this.lqn.eshift + e;
            int tidx = (int) this.lqn.parent.get(0, eidx);
            if (Double.isNaN(UN.get(tidx))) UN.set(tidx, 0);

            // Phase-2 support: utilization includes both phases
            if (this.hasPhase2 && this.servt_ph2 != null && this.servt_ph2.get(eidx - 1) > GlobalConstants.FineTol) {
                // Phase-1 utilization
                this.util_ph1.set(eidx - 1, TN.get(eidx) * this.servt_ph1.get(eidx - 1));
                // Phase-2 utilization
                this.util_ph2.set(eidx - 1, TN.get(eidx) * this.servt_ph2.get(eidx - 1));
                // Total utilization = both phases (server is busy during both)
                UN.set(eidx, this.util_ph1.get(eidx - 1) + this.util_ph2.get(eidx - 1));
            } else {
                // Standard calculation for entries without phase-2
                UN.set(eidx, TN.get(eidx) * SN.get(eidx));
            }

            if (tidx >= 0) {
                for (int i = 0; i < this.lqn.actsof.get(tidx).size(); i++) {
                    int aidx = this.lqn.actsof.get(tidx).get(i);
                    UN.set(aidx, TN.get(aidx) * SN.get(aidx));
                }
                UN.set(tidx, UN.get(tidx) + UN.get(eidx));
            }
        }

        for (double idx : this.ignore.find().toList1D()) {
            int idxInt = (int) idx;
            if (idxInt >= 0 && idxInt <= this.lqn.nidx) {
                QN.set(idxInt, 0.0);
                UN.set(idxInt, 0.0);
                RN.set(idxInt, 0.0);
                TN.set(idxInt, 0.0);
                PN.set(idxInt, 0.0);
                SN.set(idxInt, 0.0);
                WN.set(idxInt, 0.0);
                AN.set(idxInt, 0.0);
            }
        }

        // if LN standard naming
        QN = UN.copy();
        UN = PN.copy();
        RN = SN.copy();

        int maxnamelength = 12;
        for (int i = 1; i <= lqn.names.size(); i++) {
            maxnamelength = FastMath.max(maxnamelength, lqn.names.get(i).length());
        }
        /* LayeredNetworkAvgTable Generation Boilerplate */
        List<String> nodeNames = new ArrayList<>(lqn.names.values());
        List<String> nodeTypes = new ArrayList<>();

        for (int o = 0; o < nodeNames.size(); o++) {
            switch ((int) lqn.type.get(1 + o)) {
                case LayeredNetworkElement.PROCESSOR:
                    nodeTypes.add("Processor");
                    break;
                case LayeredNetworkElement.TASK:
                    if (lqn.sched.get(1 + o) == SchedStrategy.REF) {
                        nodeTypes.add("RefTask");
                    } else {
                        nodeTypes.add("Task");
                    }
                    break;
                case LayeredNetworkElement.ENTRY:
                    nodeTypes.add("Entry");
                    break;
                case LayeredNetworkElement.ACTIVITY:
                    nodeTypes.add("Activity");
                    break;
                case LayeredNetworkElement.CALL:
                    nodeTypes.add("Call");
                    break;
            }
        }

        List<Double> Qval = QN.toList1D();
        List<Double> Uval = UN.toList1D();
        List<Double> Rval = RN.toList1D();
        List<Double> Residval = WN.toList1D();
        List<Double> Aval = AN.toList1D();
        List<Double> Tval = TN.toList1D();

        Qval.remove(0);
        Uval.remove(0);
        Rval.remove(0);
        Residval.remove(0);
        Aval.remove(0);
        Tval.remove(0);
        LayeredNetworkAvgTable AvgTable = new LayeredNetworkAvgTable(Qval, Uval, Rval, Residval, Aval, Tval);
        AvgTable.setNodeNames(nodeNames);
        AvgTable.setNodeTypes(nodeTypes);
        AvgTable.setOptions(this.options);
        if (this.options.verbose == VerboseLevel.DEBUG) {
            AvgTable.print(true);
        }
        return AvgTable;
    }

    public Matrix getEntryServiceMatrix() {
        //task19:getEntryServiceMatrix function to be written
        //matrix that returns the entry servt after multiplication with residt of entries and activities
        int eshift = this.lqn.eshift;
        int sidelengthU = this.lqn.nidx + this.lqn.ncalls;
        //U starts from (0,0)
        Matrix U = new Matrix(sidelengthU, sidelengthU, sidelengthU * sidelengthU);
        int eidx;
        for (int e = 1; e <= this.lqn.nentries; e++) {
            eidx = eshift + e;
            U = getEntryServiceMatrixRecursion(this.lqn, eidx, eidx, U);
        }

        U.apply(0, 1.0, "great");
        U.apply(0, 0.0, "lessequal");
        return U;
    }

    public Matrix getEntryServiceMatrixRecursion(LayeredNetworkStruct lqn, int aidx, int eidx, Matrix U) {
        //auxiliary function to getServiceMatrix
        Matrix aidxrow = new Matrix(1, lqn.graph.getNumCols(), lqn.graph.getNumCols());
        aidxrow = Matrix.extractRows(lqn.graph, aidx, aidx + 1, aidxrow);
        Matrix nextaidxs = aidxrow.find();
        for (int i = 0; i < nextaidxs.getNumRows(); i++) {
            int nextaidx = (int) nextaidxs.get(i);
            boolean isLoop = lqn.graph.get(aidx, nextaidx) != lqn.dag.get(aidx, nextaidx);
            // in the activity graph, the following if is entered only
            // by an edge that is the return from a LOOP activity
            if (lqn.parent.get(0, aidx) != lqn.parent.get(0, nextaidx)) {
                //if the successor activity is  a call
                for (int j = 0; j < lqn.callsof.get(aidx).size(); j++) {
                    int cidx = lqn.callsof.get(aidx).get(j);
                    if (lqn.calltype.get(cidx) == CallType.SYNC) {
                        // mean number of calls alrady factored in
                        U.set(eidx - 1, lqn.nidx + cidx - 1, 1);
                    } else if (lqn.calltype.get(cidx) == CallType.ASYNC) {
                        // nop - doesn't contribute to respt
                    }
                }
            }

            //here we have processed all calls, let us do the activities now
            // if the successor activity is not a call
            if (lqn.parent.get(0, aidx) == lqn.parent.get(0, nextaidx)) {
                if (nextaidx != aidx && !isLoop) {
                    double Gvalue = lqn.graph.get(aidx, nextaidx) > 0 ? lqn.graph.get(aidx, nextaidx) : 0;
                    U.set(eidx - 1, nextaidx - 1, U.get(eidx - 1, nextaidx - 1) + Gvalue);
                    U = getEntryServiceMatrixRecursion(lqn, nextaidx, eidx, U);
                }
            }
        }
        return U;
    }

    public List<Double> getIdxhash() {
        return idxhash;
    }

    public Matrix getRoute_prob_updmap() {
        return route_prob_updmap;
    }

    public Matrix getServt_classes_updmap() {
        return servt_classes_updmap;
    }

    public Matrix getThinkt_classes_updmap() {
        return thinkt_classes_updmap;
    }

    public void init() {
        //operation before starting to iterate

        List<Double> numSet = new ArrayList<Double>();
        if (this.route_prob_updmap.getNonZeroLength() == 0) {
            this.unique_route_prob_updmap = this.route_prob_updmap;
        } else {
            this.unique_route_prob_updmap = this.route_prob_updmap.uniqueInCol(1);
//            for (int i = 1; i < this.route_prob_updmap.getNumRows(); i++) {
//                boolean unique = true;
//                // check if this.route_prob_updmap.get(i, 1) is already in numSet
//                for (double j : numSet) {
//                    if (this.route_prob_updmap.get(i, 1) == j) {
//                        unique = false;
//                        break;
//                    }
//                }
//                // add it if not
//                if (unique)
//                    numSet.add(this.route_prob_updmap.get(i, 1));
//            }
//            unique_route_prob_updmap = new Matrix(1, numSet.size(), numSet.size());
//            for (int k = 0; k < numSet.size(); k++)
//                this.unique_route_prob_updmap.set(k, numSet.get(k));
        }

        // Use row vectors (1 x nidx) consistently for all metric arrays
        this.tput = new Matrix(1, this.lqn.nidx, this.lqn.nidx);
        this.util = new Matrix(1, this.lqn.nidx, this.lqn.nidx);
        this.servt = new Matrix(1, this.lqn.nidx, this.lqn.nidx);
        this.servtmatrix = this.getEntryServiceMatrix();
        for (int e = 0; e < this.nlayers; e++)
            this.solvers[e].enableChecks = false;

        // Initialize under-relaxation state
        String relaxMode = this.options.config.relax;
        if (relaxMode == null) relaxMode = "none";
        switch (relaxMode.toLowerCase()) {
            case "auto":
                this.relax_omega = 1.0; // Start without relaxation
                break;
            case "fixed":
            case "adaptive":
                this.relax_omega = this.options.config.relax_factor;
                break;
            default: // 'none' or unrecognized
                this.relax_omega = 1.0; // No relaxation
                break;
        }
        this.relax_err_history = new ArrayList<>();
        this.servt_prev = new Matrix(1, this.lqn.nidx, this.lqn.nidx);
        this.servt_prev.fill(Double.NaN);
        this.residt_prev = new Matrix(1, this.lqn.nidx, this.lqn.nidx);
        this.residt_prev.fill(Double.NaN);
        this.tput_prev = new Matrix(1, this.lqn.nidx, this.lqn.nidx);
        this.tput_prev.fill(Double.NaN);
        this.thinkt_prev = new Matrix(1, this.lqn.ntasks + this.lqn.tshift, this.lqn.ntasks + this.lqn.ntasks - 1);
        this.thinkt_prev.fill(Double.NaN);
        this.callservt_prev = new Matrix(1, this.lqn.ncalls, this.lqn.ncalls);
        this.callservt_prev.fill(Double.NaN);
    }

    public Matrix integerMapToMatrix(Map<Integer, List<Integer[]>> cell) {
        Set<Integer> keys = cell.keySet();
        int lines = 1;
        int columns = 0;
        for (int i : keys) {
            if (!cell.get(i).isEmpty()) {
                lines = lines + cell.get(i).size();
                columns = 1 + cell.get(i).get(0).length;
            }
        }
        Matrix matrix = new Matrix(lines, columns, (lines - 1) * (columns - 1));
        int lineToAssign = 1;
        List<Integer> keysList = new ArrayList<>(keys);
        Collections.sort(keysList);
        for (int i : keysList) {
            for (int j = 0; j < cell.get(i).size(); j++) {
                for (int c = 1; c < columns; c++) {
                    matrix.set(lineToAssign, c, cell.get(i).get(j)[c - 1]);
                }
                lineToAssign++;
            }
        }
        return matrix;
    }

    public void post(int it) {

        updateMetrics(it);
        updateThinkTimes(it);

        if (this.options.config.interlocking) {
            updatePopulations(it);
        }
        updateLayers(it);
        updateRoutingProbabilities(it);

        for (int e : routereset) {
            ensemble[e - 1].refreshChains(true);
            // Note: refreshChains(true) already calls snRefreshVisits() internally,
            // matching MATLAB's refreshChains() behavior. No additional visit refresh needed.
            solvers[e - 1].reset();
        }

        for (int e : svcreset) {
            List<Integer> statSet = new ArrayList<>();
            List<Integer> classSet = new ArrayList<>();
            for (int i = 0; i < ensemble[e - 1].getNumberOfClasses(); i++) {
                classSet.add(i);
            }
            for (int i = 0; i < ensemble[e - 1].getNumberOfStations(); i++) {
                statSet.add(i);
            }
            String solverName = solvers[e - 1].name;
            if ("SolverMVA".equals(solverName) || "SolverNC".equals(solverName)) {
                // Leaner than refreshProcesses, no need to refresh phases
                if ("default".equals(this.options.method)) {
                    ensemble[e - 1].refreshRates(null, null);
                } else if ("moment3".equals(this.options.method)) {
                    ensemble[e - 1].refreshProcesses();
                }
            } else {
                ensemble[e - 1].refreshProcesses();
            }
            solvers[e - 1].reset();
        }

        if (this.options.config.interlocking) {
            for (int e = 0; e < this.nlayers; e++) {
                ensemble[e].refreshJobs();
            }
        }

        if (it == 1) {
            for (int e = 0; e < ensemble.length; e++) {
                solvers[e].enableChecks = false;
            }
        }
    }

    public void pre(int it) {
        //no op
    }

    @Override
    public void runAnalyzer() throws IllegalAccessException {


    }

    public boolean supports(Ensemble ensemble) {
        boolean bool = true;
        for (int e = 0; e < ensemble.size(); e++) {
            bool = bool && solvers[e].supports(ensemble.getModel(e));
        }
        return bool;
    }

    public void updateLayers(int it) {
        //task14: updateLayers function to be written
        // reassign service times
        for (int r = 1; r < this.thinkt_classes_updmap.getNumRows(); r++) {
            int ri;
            if (it % 2 == 1) { // elevator - alternate direction on odd iterations
                ri = this.thinkt_classes_updmap.getNumRows() - r;
            } else {
                ri = r;
            }

            double idx = this.thinkt_classes_updmap.get(ri, 1);
            double aidx = this.thinkt_classes_updmap.get(ri, 2);
            double nodeidx = this.thinkt_classes_updmap.get(ri, 3);
            double classidx = this.thinkt_classes_updmap.get(ri, 4);
            int idxInt = (int) idx;
            if (idxInt < 0 || idxInt >= idxhash.size()) {
                continue;  // Skip if index is out of bounds
            }
            JobClass tmp_class = this.ensemble[this.idxhash.get(idxInt).intValue() - 1].getClassByIndex((int) classidx - 1);
            // here update the number of jobs in the task chain
            if (aidx <= (this.lqn.tshift + this.lqn.ntasks)) {
                // aidx here is actually set to tidx in buildLayersRecursive
                if (tmp_class.getJobClassType() == JobClassType.CLOSED) {
                    if (this.options.config.interlocking) {
                        ClosedClass tmp_class_c = (ClosedClass) tmp_class;
                        tmp_class_c.setPopulation(this.njobs.get((int) aidx, (int) idx));
                    }
                }
            }
            ServiceStation node = (ServiceStation) this.ensemble[idxhash.get((int) idx).intValue() - 1].getNodeByStatefulIndex((int) nodeidx - 1);
            // Case 1
            if ((int) nodeidx == this.ensemble[(idxhash.get((int) idx).intValue()) - 1].getAttribute().getClientIdx()) {
                if (this.lqn.type.get((int) aidx) == LayeredNetworkElement.TASK) {
                    if (this.lqn.sched.get((int) aidx) != SchedStrategy.REF) {
                        if (this.thinktproc.get((int) aidx) != null) {
                            node.setService(tmp_class, this.thinktproc.get((int) aidx));
                        }
                    } else {
                        node.setService(tmp_class, this.servtproc.get((int) aidx));
                    }
                } else {
                    node.setService(tmp_class, this.servtproc.get((int) aidx));
                }
            }
            // Case 2
            if ((int) nodeidx == this.ensemble[this.idxhash.get((int) idx).intValue() - 1].getAttribute().getServerIdx()) {
                node.setService(tmp_class, this.servtproc.get((int) aidx));
            }
        }

        // reassign arrival rates
        for (int r = 1; r < this.arvproc_classes_updmap.getNumRows(); r++) {
            int ri;
            if (it % 2 == 1) { // elevator - alternate direction on odd iterations
                ri = this.arvproc_classes_updmap.getNumRows() - r;
            } else {
                ri = r;
            }
            double idx = this.arvproc_classes_updmap.get(ri, 1);
            double eidx_or_cidx = this.arvproc_classes_updmap.get(ri, 2);
            double nodeidx = this.arvproc_classes_updmap.get(ri, 3);
            double classidx = this.arvproc_classes_updmap.get(ri, 4);
            JobClass tmp_class = this.ensemble[this.idxhash.get((int) idx).intValue() - 1].getClassByIndex((int) classidx - 1);
            Source node = (Source) this.ensemble[this.idxhash.get((int) idx).intValue() - 1].getNodeByStatefulIndex((int) nodeidx - 1);

            if (eidx_or_cidx < 0) {  // Entry-level arrival (negative index)
                int eidx = -(int) eidx_or_cidx;
                node.setArrival(tmp_class, this.lqn.arrival.get(eidx));
            } else {  // Async call arrival (positive index)
                int cidx = (int) eidx_or_cidx;
                node.setArrival(tmp_class, this.tputproc.get((int) this.lqn.callpair.get(cidx, 1)));
            }
        }

        // reassign call service time / response time
        for (int c = 1; c < this.call_classes_updmap.getNumRows(); c++) {
            int ci;
            if (it % 2 == 1) { // elevator - alternate direction on odd iterations
                ci = this.call_classes_updmap.getNumRows() - c;
            } else {
                ci = c;
            }
            double idx = this.call_classes_updmap.get(ci, 1);
            double cidx = this.call_classes_updmap.get(ci, 2);
            double nodeidx = this.call_classes_updmap.get(ci, 3);
            double classidx = this.call_classes_updmap.get(ci, 4);
            JobClass tmp_class = this.ensemble[this.idxhash.get((int) idx).intValue() - 1].getClassByIndex((int) classidx - 1);
            Queue node = (Queue) this.ensemble[this.idxhash.get((int) idx).intValue() - 1].getNodeByStatefulIndex((int) nodeidx - 1);

            // Case 1 - client
            if ((int) nodeidx == this.ensemble[this.idxhash.get((int) idx).intValue() - 1].getAttribute().getClientIdx()) {
                node.setService(tmp_class, this.callservtproc.get((int) cidx));
            }

            // Case 2 - the call is processed by the server, then replace with the svc time
            if ((int) nodeidx == this.ensemble[this.idxhash.get((int) idx).intValue() - 1].getAttribute().getServerIdx()) {
                double eidx = this.lqn.callpair.get((int) cidx, 2);
                node.setService(tmp_class, this.servtproc.get((int) eidx));
            }
        }
    }

    public void updateMetrics(int it) {
        switch (options.method) {
            case "default":
                this.updateMetricsDefault(it);
                break;
            case "moment3":
                this.updateMetricsMomentBased(it);
                break;
        }
    }

    public void updateMetricsDefault(int it) {
        LayeredNetworkStruct lqn = this.lqn;
        int rLen = this.results.size();

        // Safety check: if results is empty, skip metric updates
        if (this.results == null || this.results.isEmpty() || rLen == 0) {
            return;
        }

        // obtain the activity service times
        this.servt = new Matrix(1, lqn.nidx, lqn.nidx);
        this.residt = new Matrix(1, lqn.nidx, lqn.nidx);
        for (int r = 1; r < this.servt_classes_updmap.getNumRows(); r++) {
            int idx = (int) this.servt_classes_updmap.get(r, 1);     //layer
            int aidx = (int) this.servt_classes_updmap.get(r, 2);    //activity
            int nodeidx = (int) this.servt_classes_updmap.get(r, 3); //node
            int classidx = (int) this.servt_classes_updmap.get(r, 4); //jobclass

            // store the residence times and tput at this layer to become
            // the servt / tputs of aidx in another layer, as needed
            // this.servt starts from 0 for JLineMatrix Multiplication

            // Compute residt from QN/TN_ref instead of WN to avoid
            // fork+loop visit distortion (WN uses visits from DTMC solve
            // which are distorted when Fork non-stochastic rows coexist
            // with loop back-edges in the routing matrix)
            // (matches MATLAB updateMetricsDefault.m)
            int layerIdx_0 = this.idxhash.get(idx).intValue() - 1;
            NetworkStruct layerSn = this.ensemble[layerIdx_0].getStruct(false);
            int chainIdx = -1;
            if (layerSn != null && layerSn.chains != null && !layerSn.chains.isEmpty()) {
                for (int ch = 0; ch < layerSn.nchains; ch++) {
                    if (layerSn.chains.get(ch, classidx - 1) > 0) {
                        chainIdx = ch;
                        break;
                    }
                }
            }
            int refclass_c = (chainIdx >= 0 && layerSn.refclass != null) ? (int) layerSn.refclass.get(chainIdx) : classidx - 1;
            int refstat_k = (layerSn != null && layerSn.refstat != null) ? (int) layerSn.refstat.get(classidx - 1) : 0;

            int iter_min = (int) FastMath.min(30, FastMath.ceil(this.options.iter_max / 4.0));
            if (this.averagingstart != null && it >= iter_min) {
                int wnd_size = it - this.averagingstart + 1;
                this.servt.set(aidx - 1, 0);
                this.residt.set(aidx - 1, 0);
                this.tput.set(aidx - 1, 0);
                for (int w = 1; w < wnd_size; w++) {
                    this.servt.set(aidx - 1, this.servt.get(aidx - 1) + this.results.get(rLen - w).get(layerIdx_0).RN.get(nodeidx - 1, classidx - 1) / wnd_size);
                    double TN_ref_w = this.results.get(rLen - w).get(layerIdx_0).TN.get(refstat_k, refclass_c);
                    if (TN_ref_w > GlobalConstants.FineTol) {
                        this.residt.set(aidx - 1, this.residt.get(aidx - 1) + this.results.get(rLen - w).get(layerIdx_0).QN.get(nodeidx - 1, classidx - 1) / TN_ref_w / wnd_size);
                    } else {
                        this.residt.set(aidx - 1, this.residt.get(aidx - 1) + this.results.get(rLen - w).get(layerIdx_0).WN.get(nodeidx - 1, classidx - 1) / wnd_size);
                    }
                    this.tput.set(aidx - 1, this.tput.get(aidx - 1) + this.results.get(rLen - w).get(layerIdx_0).TN.get(nodeidx - 1, classidx - 1) / wnd_size);
                }
            } else {
                this.servt.set(aidx - 1, this.results.get(results.size()).get(layerIdx_0).RN.get(nodeidx - 1, classidx - 1));
                double TN_ref = this.results.get(results.size()).get(layerIdx_0).TN.get(refstat_k, refclass_c);
                if (TN_ref > GlobalConstants.FineTol) {
                    this.residt.set(aidx - 1, this.results.get(results.size()).get(layerIdx_0).QN.get(nodeidx - 1, classidx - 1) / TN_ref);
                } else {
                    this.residt.set(aidx - 1, this.results.get(results.size()).get(layerIdx_0).WN.get(nodeidx - 1, classidx - 1));
                }
                this.tput.set(aidx - 1, this.results.get(results.size()).get(layerIdx_0).TN.get(nodeidx - 1, classidx - 1));
            }

            // Force servt/residt to 0 for activities with Immediate service time
            // The NC solver may return non-zero RN due to numerical issues, but
            // Immediate activities have zero service time by definition
            if (lqn.hostdem.get(aidx) instanceof Immediate) {
                this.servt.set(aidx - 1, 0);
                this.residt.set(aidx - 1, 0);
            }

            // Fix for async-only entry targets: use RN (response time per visit) for residt
            // The host layer closed model incorrectly splits residence time (WN) between
            // activities when an entry only receives async calls (no sync callers).
            // For async-only entries, use RN instead of WN since the async arrivals
            // don't share the closed chain's visit ratio - each async arrival gets
            // the full response time per visit.
            if (aidx > lqn.ashift && aidx <= lqn.ashift + lqn.nacts) {
                // This is an activity - find its bound entry
                for (int eidx = lqn.eshift + 1; eidx <= lqn.eshift + lqn.nentries; eidx++) {
                    if (lqn.graph.get(eidx, aidx) > 0) {
                        // Found bound entry - check if async-only
                        boolean hasSyncCallers = false;
                        boolean hasAsyncCallers = false;
                        for (int i = 1; i <= lqn.nidx; i++) {
                            if (lqn.issynccaller.get(i, eidx) > 0) {
                                hasSyncCallers = true;
                            }
                            if (lqn.isasynccaller.get(i, eidx) > 0) {
                                hasAsyncCallers = true;
                            }
                        }
                        if (hasAsyncCallers && !hasSyncCallers) {
                            // Async-only target: use RN (response time per visit)
                            // instead of WN (residence time with visit ratio)
                            this.residt.set(aidx - 1, this.servt.get(aidx - 1)); // servt already has RN
                        }
                        break;
                    }
                }
            }

            // Apply under-relaxation if enabled and not first iteration
            double omega = this.relax_omega;
            if (omega < 1.0 && it > 1) {
                if (!Double.isNaN(this.servt_prev.get(aidx - 1))) {
                    this.servt.set(aidx - 1, omega * this.servt.get(aidx - 1) + (1 - omega) * this.servt_prev.get(aidx - 1));
                }
                if (!Double.isNaN(this.residt_prev.get(aidx - 1))) {
                    this.residt.set(aidx - 1, omega * this.residt.get(aidx - 1) + (1 - omega) * this.residt_prev.get(aidx - 1));
                }
                if (!Double.isNaN(this.tput_prev.get(aidx - 1))) {
                    this.tput.set(aidx - 1, omega * this.tput.get(aidx - 1) + (1 - omega) * this.tput_prev.get(aidx - 1));
                }
            }
            // Store current values for next iteration
            this.servt_prev.set(aidx - 1, this.servt.get(aidx - 1));
            this.residt_prev.set(aidx - 1, this.residt.get(aidx - 1));
            this.tput_prev.set(aidx - 1, this.tput.get(aidx - 1));

            // Preserve Immediate type for activities that originally had Immediate service times
            if (lqn.hostdem.get(aidx) instanceof Immediate) {
                this.servtproc.put(aidx, Immediate.getInstance());
            } else {
                this.servtproc.put(aidx, Exp.fitMean(this.servt.get(aidx - 1)));
            }
            this.tputproc.put(aidx, Exp.fitRate(this.tput.get(aidx - 1)));
        }

        // Phase-2 support: split activity service times by phase
        if (this.hasPhase2) {
            // Reset phase-specific arrays
            this.servt_ph1 = new Matrix(1, lqn.nidx, lqn.nidx);
            this.servt_ph2 = new Matrix(1, lqn.nidx, lqn.nidx);

            // Split activity service times by phase
            for (int a = 1; a <= lqn.nacts; a++) {
                int aidx = lqn.ashift + a;
                if (lqn.actphase.get(0, a) == 1) {
                    this.servt_ph1.set(aidx - 1, this.servt.get(aidx - 1));
                } else {
                    this.servt_ph2.set(aidx - 1, this.servt.get(aidx - 1));
                }
            }

            // Aggregate phase service times to entry level
            for (int e = 1; e <= lqn.nentries; e++) {
                int eidx = lqn.eshift + e;
                List<Integer> acts = lqn.actsof.get(eidx);
                if (acts != null) {
                    for (int aidx : acts) {
                        int a = aidx - lqn.ashift;
                        if (a > 0 && a <= lqn.nacts) {
                            if (lqn.actphase.get(0, a) == 1) {
                                this.servt_ph1.set(eidx - 1, this.servt_ph1.get(eidx - 1) + this.servt_ph1.get(aidx - 1));
                            } else {
                                this.servt_ph2.set(eidx - 1, this.servt_ph2.get(eidx - 1) + this.servt_ph2.get(aidx - 1));
                            }
                        }
                    }
                }
            }
        }

        // obtain throughput for activities in thinkt_classes_updmap (needed for async calls)
        // this ensures tputproc is set for activities that make async calls from client nodes
        for (int r = 1; r < this.thinkt_classes_updmap.getNumRows(); r++) {
            int idx = (int) this.thinkt_classes_updmap.get(r, 1);     // layer
            int aidx = (int) this.thinkt_classes_updmap.get(r, 2);    // activity
            int nodeidx = (int) this.thinkt_classes_updmap.get(r, 3); // node
            int classidx = (int) this.thinkt_classes_updmap.get(r, 4); // jobclass

            // only update if not already set by servt_classes_updmap processing
            if (!this.tputproc.containsKey(aidx)) {
                int iter_min = (int) FastMath.min(30, FastMath.ceil(this.options.iter_max / 4.0));
                if (this.averagingstart != null && it >= iter_min) {
                    int wnd_size = it - this.averagingstart + 1;
                    this.tput.set(aidx - 1, 0);
                    for (int w = 1; w < wnd_size; w++) {
                        this.tput.set(aidx - 1, this.tput.get(aidx - 1) + this.results.get(rLen - w).get(this.idxhash.get(idx).intValue() - 1).TN.get(nodeidx - 1, classidx - 1) / wnd_size);
                    }
                } else {
                    this.tput.set(aidx - 1, this.results.get(results.size()).get(this.idxhash.get(idx).intValue() - 1).TN.get(nodeidx - 1, classidx - 1));
                }
                this.tputproc.put(aidx, Exp.fitRate(this.tput.get(aidx - 1)));
            }
        }

        // obtain the call residence time
        this.callservt = new Matrix(1, lqn.ncalls, lqn.ncalls);
        this.callresidt = new Matrix(1, lqn.ncalls, lqn.ncalls);
        for (int r = 1; r < this.call_classes_updmap.getNumRows(); r++) {
            int idx = (int) this.call_classes_updmap.get(r, 1);     // layer
            int cidx = (int) this.call_classes_updmap.get(r, 2);    // call
            int nodeidx = (int) this.call_classes_updmap.get(r, 3);// node
            int classidx = (int) this.call_classes_updmap.get(r, 4);// jobclass

            if (this.call_classes_updmap.get(r, 3) > 1) {
                if (nodeidx == 1) {
                    this.callservt.set(cidx - 1, 0.0);
                } else {
                    this.callservt.set(cidx - 1, this.results.get(results.size()).get(this.idxhash.get(idx).intValue() - 1).RN.get(nodeidx - 1, classidx - 1) * this.lqn.callproc_mean.getOrDefault(cidx, 1.0));
                    this.callresidt.set(cidx - 1, this.results.get(results.size()).get(this.idxhash.get(idx).intValue() - 1).WN.get(nodeidx - 1, classidx - 1));
                }
                // Apply under-relaxation to call service times
                double omega = this.relax_omega;
                if (omega < 1.0 && it > 1 && !Double.isNaN(this.callservt_prev.get(cidx - 1))) {
                    this.callservt.set(cidx - 1, omega * this.callservt.get(cidx - 1) + (1 - omega) * this.callservt_prev.get(cidx - 1));
                }
                this.callservt_prev.set(cidx - 1, this.callservt.get(cidx - 1));
            }
        }

        //then resolve the entry servt summming up these contributions
        Matrix out = new Matrix(1, lqn.nidx + lqn.ncalls + 1, lqn.nidx + lqn.ncalls);
        Matrix entry_servt = new Matrix(this.servtmatrix.getNumRows(), 1, this.servtmatrix.getNumRows());
        Matrix.concatColumns(this.residt, this.callresidt, out);
        this.servtmatrix.mult(out.transpose(), entry_servt);

        for (int i = 0; i < lqn.eshift; i++) {
            entry_servt.set(i, 0, 0);
        }


        // this block fixes the problem that ResidT is scaled so that the task as Vtask = 1,
        // but in call servt the entries need to have Ventry = 1
        for (int eidx = lqn.eshift + 1; eidx < lqn.eshift + lqn.nentries + 1; eidx++) {
            int tidx = (int) lqn.parent.get(0, eidx); //  task of entry
            if (tidx < 0) {
                continue;  // Skip if task index is invalid
            }
            int hidx = (int) lqn.parent.get(0, tidx); // host of entry
            // Skip ignored tasks and hosts
            if (this.ignore.get(tidx) != 0 || this.ignore.get(hidx) != 0) {
                continue;
            }

            // Check if this entry has sync callers (which create closed classes)
            boolean hasSyncCallers = false;
            for (int i = 1; i <= lqn.nidx; i++) {
                if (lqn.issynccaller.get(i, eidx) > 0) {
                    hasSyncCallers = true;
                    break;
                }
            }

            if (hasSyncCallers) {
                // Original logic for entries with sync callers
                // get class in host layer of task and entry
                List<Integer> tidxclass = new ArrayList<Integer>();
                List<Integer> eidxclass = new ArrayList<Integer>();

                for (int i = 1; i <= ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getTasks().size(); i++) {
                    if (tidx == ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getTasks().get(i)[1]) {
                        if (ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getTasks().get(i)[0] != null)
                            tidxclass.add(ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getTasks().get(i)[0]);
                    }
                }

                Map<Integer, Integer[]> m = ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getEntries();
                for (int i = 1; i <= ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getEntries().size(); i++) {
                    if (eidx == ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getEntries().get(i)[1]) {
                        eidxclass.add(ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getEntries().get(i)[0]);
                    }
                }

                double task_tput = 0;
                double entry_tput = 0;

                for (int i = 0; i < tidxclass.size(); i++) {
                    task_tput += this.results.get(results.size()).get(this.idxhash.get(hidx).intValue() - 1).TN.get(ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getClientIdx() - 1, tidxclass.get(i) - 1);
                }


                for (int i = 0; i < eidxclass.size(); i++) {
                    entry_tput += this.results.get(results.size()).get(this.idxhash.get(hidx).intValue() - 1).TN.get(ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getClientIdx() - 1, eidxclass.get(i) - 1);
                }

                //ServiceStation entry_refstat = ((ServiceStation) ensemble[this.idxhash.get(hidx).intValue() - 1].getClassByIndex(tidxclass.get(0) - 1).getReferenceStation());
                //double entry_servt_z = entry_refstat.getServiceProcess(ensemble[this.idxhash.get(hidx).intValue() - 1].getClassByIndex(tidxclass.get(0) - 1)).getMean();
                //entry_servt.set(eidx - 1, ensemble[this.idxhash.get(hidx).intValue() - 1].getClassByIndex(tidxclass.get(0) - 1).getNumberOfJobs() / entry_tput - entry_servt_z);
                //System.out.println("value: "+entry_servt.get(eidx - 1) * task_tput / FastMath.max(GlobalConstants.Zero, entry_tput));
                this.servt.set(eidx - 1, entry_servt.get(eidx - 1) * task_tput / FastMath.max(GlobalConstants.Zero, entry_tput));
                this.residt.set(eidx - 1, entry_servt.get(eidx - 1) * task_tput / FastMath.max(GlobalConstants.Zero, entry_tput));
            } else {
                // For async-only targets, use entry_servt directly
                // No throughput ratio scaling needed since there are no closed classes
                this.servt.set(eidx - 1, entry_servt.get(eidx - 1));
                this.residt.set(eidx - 1, entry_servt.get(eidx - 1));
            }
        }

        // Phase-2 support: compute overtaking probability and apply correction
        // This must happen AFTER entry throughput is available (computed above)
        if (this.hasPhase2) {
            for (int e = 1; e <= lqn.nentries; e++) {
                int eidx = lqn.eshift + e;
                int tidx = (int) lqn.parent.get(0, eidx);

                if (this.servt_ph2.get(eidx - 1) > GlobalConstants.FineTol) {
                    // Get entry throughput (use task throughput as approximation if entry not available)
                    double entry_tput;
                    if (this.tput.get(eidx - 1) > GlobalConstants.FineTol) {
                        entry_tput = this.tput.get(eidx - 1);
                    } else if (tidx > 0 && this.tput.get(tidx - 1) > GlobalConstants.FineTol) {
                        entry_tput = this.tput.get(tidx - 1);
                    } else {
                        entry_tput = 0;
                    }

                    // Compute overtaking probability now that throughput is available
                    if (entry_tput > GlobalConstants.FineTol) {
                        this.prOvertake.set(e - 1, this.overtakeProb(eidx));
                    } else {
                        this.prOvertake.set(e - 1, 0);
                    }

                    // Caller's response time = phase-1 only + P(overtake) * phase-2
                    double overtake_delay = this.prOvertake.get(e - 1) * this.servt_ph2.get(eidx - 1);

                    // The caller sees phase-1 + overtaking correction (not full phase-2)
                    this.residt.set(eidx - 1, this.servt_ph1.get(eidx - 1) + overtake_delay);
                }
            }
        }

        for (int r = 1; r < this.call_classes_updmap.getNumRows(); r++) {
            int cidx = (int) this.call_classes_updmap.get(r, 2);
            int eidx = (int) lqn.callpair.get(cidx, 2);
            if (this.call_classes_updmap.get(r, 3) > 1) {
                this.servtproc.put(eidx, Exp.fitMean(this.servt.get(eidx - 1)));
            }
        }

        // determine call response time processes
        // this.callresidt starts from 0
        for (int r = 1; r < this.call_classes_updmap.getNumRows(); r++) {
            int cidx = (int) this.call_classes_updmap.get(r, 2);
            int eidx = (int) lqn.callpair.get(cidx, 2);
            if (this.call_classes_updmap.get(r, 3) > 1) {
                if (it == 1) {
                    // note that respt is per visit, so number of calls is 1
                    this.callservt.set(cidx - 1, this.servt.get(eidx - 1));
                    this.callservtproc.put(cidx, this.servtproc.get(eidx));
                } else {
                    // note that respt is per visit, so number of calls is 1
                    this.callservtproc.put(cidx, Exp.fitMean(this.callservt.get(cidx - 1)));
                }
            }
        }

        this.ptaskcallers = new Matrix(this.ptaskcallers.getNumRows(), this.ptaskcallers.getNumCols(), this.ptaskcallers.getNumRows() * this.ptaskcallers.getNumCols() - 1);

        for (int i = 0; i < this.ptaskcallers.getNumRows(); i++) {
            for (int j = 0; j < this.ptaskcallers.getNumCols(); j++) {
                this.ptaskcallers.set(i, j, 0);
            }
        }

        // determine ptaskcallers for direct callers to tasks
        for (int t = 1; t <= lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            if ((int) lqn.isref.get(tidx) == 0) {
                List<Integer> calling_idx = new ArrayList<>();
                for (int entry : lqn.entriesof.get(tidx)) {
                    for (int i = 0; i < lqn.iscaller.getNumRows(); i++) {
                        if (lqn.iscaller.get(i, entry) != 0) {
                            calling_idx.add(i);
                        }
                    }
                }

                List<Integer> uniqueCallingIdx = new ArrayList<>();
                for (int idx : calling_idx) {
                    if (!uniqueCallingIdx.contains(idx)) {
                        uniqueCallingIdx.add(idx);
                    }
                }
                Collections.sort(uniqueCallingIdx);

                SortedSet<Integer> callers = new TreeSet<>();
                for (int row = lqn.tshift + 1; row < lqn.tshift + lqn.ntasks + 1; row++) {
                    if (uniqueCallingIdx.contains(row)) {
                        callers.add(row);
                    }
                }

                Matrix caller_tput = new Matrix(lqn.ntasks + 1, 1, lqn.ntasks);
                // Skip tasks without a valid layer mapping (ignored, isolated, etc.)
                if (tidx >= idxhash.size() || Double.isNaN(idxhash.get(tidx))) {
                    continue;
                }
                for (int caller_idx : callers) {
                    List<Double> caller_idxclass = new ArrayList<Double>();
                    List<Integer> keys = new ArrayList<Integer>();
                    Map<Integer, Integer[]> taskmap = this.ensemble[idxhash.get(tidx).intValue() - 1].getAttribute().getTasks();
                    for (int i = 2; i <= taskmap.size(); i++) {
                        keys.add(taskmap.get(i)[1]);
                    }
                    Integer caller_idx_found_index = 0;
                    for (int i = 0; i < keys.size(); i++) {
                        if (keys.get(i) == caller_idx) {
                            caller_idx_found_index = i + 1;
                        }
                    }
                    // Skip callers that are not in the taskmap (e.g., async-only callers)
                    // In MATLAB, find() returns empty and the loop body is effectively skipped
                    if (caller_idx_found_index == 0) {
                        continue;
                    }
                    caller_idxclass.add((double) taskmap.get(1 + caller_idx_found_index)[0]);

                    double sum = 0;
                    for (double j : caller_idxclass) {
                        Matrix Tn = results.get(results.size()).get(idxhash.get(tidx).intValue() - 1).TN;
                        sum += Tn.get(this.ensemble[idxhash.get(tidx).intValue() - 1].getAttribute().getClientIdx() - 1, (int) j - 1);
                    }
                    caller_tput.set(caller_idx - lqn.tshift, sum);
                }
                double task_tput = 0;
                for (int i = 1; i < caller_tput.getNumRows(); i++) {
                    task_tput += caller_tput.get(i);
                }
                for (int i = 1; i <= lqn.ntasks; i++) {
                    this.ptaskcallers.set(tidx, lqn.tshift + i, caller_tput.get(i) / FastMath.max(GlobalConstants.Zero, task_tput));
                }
            }
        }

        // determine ptaskcallers for direct callers to hosts
        for (int hidx = 1; hidx <= lqn.nhosts; hidx++) {
            // Skip ignored hosts
            if (this.ignore.get(hidx) != 0) {
                continue;
            }
            Matrix caller_tput = new Matrix(1, lqn.ntasks + 1, lqn.ntasks);
            List<Integer> callers = lqn.tasksof.get(hidx);

            for (int caller_idx : callers) {
                List<Integer> caller_idxclass = new ArrayList<Integer>();
                for (int i = 1; i <= ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getTasks().size(); i++) {
                    if (caller_idx == ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getTasks().get(i)[1]) {
                        caller_idxclass.add(ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getTasks().get(i)[0]);
                    }
                }
                double sum = 0;
                for (int i : caller_idxclass) {
                    sum += this.results.get(this.results.size()).get(this.idxhash.get(hidx).intValue() - 1).TN.get(this.ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getClientIdx() - 1, i - 1);
                }
                caller_tput.set(caller_idx - lqn.tshift, caller_tput.get(caller_idx - lqn.tshift) + sum);
            }
            double host_tput = caller_tput.elementSum();
            for (int i = 1; i <= lqn.ntasks; i++) {
                this.ptaskcallers.set(hidx, lqn.tshift + i, caller_tput.get(i) / FastMath.max(GlobalConstants.Zero, host_tput));
            }
        }

        // impute call probability using a DTMC random walk on the taskcaller graph
        // for matrix multiplication, let P to be the same size of that in Matlab
        Matrix P = new Matrix(this.ptaskcallers.getNumRows() - 1, this.ptaskcallers.getNumCols() - 1, (this.ptaskcallers.getNumRows() - 1) * (this.ptaskcallers.getNumCols() - 1));
        for (int i = 1; i < this.ptaskcallers.getNumRows(); i++) {
            for (int j = 1; j < this.ptaskcallers.getNumCols(); j++) {
                P.set(i - 1, j - 1, this.ptaskcallers.get(i, j));
            }
        }
        P = dtmc_makestochastic(P); // hold mass at reference stations when there

        this.ptaskcallers_step.put(1, P.copy()); // configure step = 1
        // Random walk block validated: indexing matches MATLAB updateMetricsDefault.m lines 319-351
        for (int h = 1; h <= lqn.nhosts; h++) {
            int hidx = h;
            for (int i = 0; i < lqn.tasksof.get(hidx).size(); i++) {
                int tidx = lqn.tasksof.get(hidx).get(i);
                // initialize the probability mass on tidx
                Matrix x0 = new Matrix(1, P.length(), P.length());
                x0.set(hidx - 1, 1);
                // start the walk backward to impute probability of indirect callers
                Matrix x = new Matrix(x0.getNumRows(), P.getNumCols(), x0.getNumRows() * P.getNumCols());
                x0.mult(P, x);
                for (int step = 2; step <= this.nlayers; step++) {
                    Matrix xret = Matrix.createLike(x0);
                    x.mult(P, xret);
                    x.setTo(xret);

                    // ptaskcallers_step.get(step) is a matrix that is not 0-padded -- index starts at 0
                    for (int remidx = 0; remidx < this.ptaskcallers_step.get(step).getNumCols(); remidx++) {
                        this.ptaskcallers_step.get(step).set(tidx - 1, remidx, x.get(remidx));
                        double scaled = this.ptaskcallers.get(hidx, tidx) * x.get(remidx);
                        this.ptaskcallers_step.get(step).set(hidx - 1, remidx, scaled);
                    }

                    double sum = 0;
                    Matrix ref_nonzero = lqn.isref.find();
                    for (int k = 0; k < ref_nonzero.getNumElements(); k++) {
                        sum += x.get((int) ref_nonzero.get(k) - 1);
                    }
                    if (sum > 1.0 - this.options.tol) break;

                    // ptaskcallers is a matrix that is 0-padded -- index starts at 1
                    for (int index = 1; index < this.ptaskcallers.getNumRows(); index++) {
                        double max = FastMath.max(this.ptaskcallers.get(index, tidx), x.get(0, index - 1));
                        max = max >= 0 ? max : 0;
                        this.ptaskcallers.set(index, tidx, max);
                    }
                }
            }
        }
    }

    public void updateMetricsMomentBased(int it) {
        LayeredNetworkStruct lqn = this.lqn;

        // This method propagates through the layers 3 moments of the
        // response time distribution computed from the CDF obtained by the
        // solvers of the individual layers. In the present implementation,
        // calls are still assumed to be exponentially distributed.

        if (!this.hasconverged) {
            // ===== PRE-CONVERGENCE: Mean-based propagation using exponential fits =====

            // First obtain servt of activities at hostlayers
            this.servt = new Matrix(1, lqn.nidx, lqn.nidx);
            this.residt = new Matrix(1, lqn.nidx, lqn.nidx);
            for (int r = 1; r < this.servt_classes_updmap.getNumRows(); r++) {
                int idx = (int) this.servt_classes_updmap.get(r, 1);     // layer
                int aidx = (int) this.servt_classes_updmap.get(r, 2);    // activity
                int nodeidx = (int) this.servt_classes_updmap.get(r, 3); // node
                int classidx = (int) this.servt_classes_updmap.get(r, 4); // jobclass

                // Use RN as indicated in MATLAB version (with debugging note)
                this.servt.set(aidx - 1, this.results.get(results.size()).get(this.idxhash.get(idx).intValue() - 1).RN.get(nodeidx - 1, classidx - 1));
                this.tput.set(aidx - 1, this.results.get(results.size()).get(this.idxhash.get(idx).intValue() - 1).TN.get(nodeidx - 1, classidx - 1));
                this.servtproc.put(aidx, Exp.fitMean(this.servt.get(aidx - 1)));

                // Compute residt from QN/TN_ref (matching MATLAB updateMetricsMomentBased)
                int layerIdx_0 = this.idxhash.get(idx).intValue() - 1;
                NetworkStruct layerSn = this.ensemble[layerIdx_0].getStruct(false);
                int chainIdx = -1;
                if (layerSn != null && layerSn.chains != null && !layerSn.chains.isEmpty()) {
                    for (int ch = 0; ch < layerSn.nchains; ch++) {
                        if (layerSn.chains.get(ch, classidx - 1) > 0) {
                            chainIdx = ch;
                            break;
                        }
                    }
                }
                int refclass_c = (chainIdx >= 0 && layerSn.refclass != null) ? (int) layerSn.refclass.get(chainIdx) : classidx - 1;
                int refstat_k = (layerSn != null && layerSn.refstat != null) ? (int) layerSn.refstat.get(classidx - 1) : 0;
                double TN_ref = this.results.get(results.size()).get(layerIdx_0).TN.get(refstat_k, refclass_c);
                if (TN_ref > GlobalConstants.FineTol) {
                    this.residt.set(aidx - 1, this.results.get(results.size()).get(layerIdx_0).QN.get(nodeidx - 1, classidx - 1) / TN_ref);
                } else {
                    this.residt.set(aidx - 1, this.results.get(results.size()).get(layerIdx_0).WN.get(nodeidx - 1, classidx - 1));
                }
            }

            // Estimate call response times at hostlayers
            this.callservt = new Matrix(1, lqn.ncalls, lqn.ncalls);
            this.callresidt = new Matrix(1, lqn.ncalls, lqn.ncalls);
            for (int r = 1; r < this.call_classes_updmap.getNumRows(); r++) {
                int idx = (int) this.call_classes_updmap.get(r, 1);     // layer
                int cidx = (int) this.call_classes_updmap.get(r, 2);    // call
                int nodeidx = (int) this.call_classes_updmap.get(r, 3); // node
                int classidx = (int) this.call_classes_updmap.get(r, 4); // jobclass

                if (this.call_classes_updmap.get(r, 3) > 1) {
                    if (nodeidx == 1) {
                        this.callservt.set(cidx - 1, 0.0);
                        this.callresidt.set(cidx - 1, 0.0);
                    } else {
                        // Include call multiplicity in callservt (matching MATLAB)
                        this.callservt.set(cidx - 1, this.results.get(results.size()).get(this.idxhash.get(idx).intValue() - 1).RN.get(nodeidx - 1, classidx - 1) * this.lqn.callproc_mean.getOrDefault(cidx, 1.0));
                        // callresidt uses WN which already includes visit multiplicity
                        this.callresidt.set(cidx - 1, this.results.get(results.size()).get(this.idxhash.get(idx).intValue() - 1).WN.get(nodeidx - 1, classidx - 1));
                    }
                }
            }

            // Then resolve the entry servt summing up these contributions
            Matrix out = new Matrix(1, lqn.nidx + lqn.ncalls, lqn.nidx + lqn.ncalls);
            Matrix entry_servt = new Matrix(this.servtmatrix.getNumRows(), 1, this.servtmatrix.getNumRows());
            Matrix.concatColumns(this.servt, this.callservt, out);

            // Solve the system: entry_servt = (I - servtmatrix)^(-1) * [servt; callservt]
            Matrix identity = Matrix.eye(lqn.nidx + lqn.ncalls);
            Matrix system = identity.sub(this.servtmatrix);
            Matrix rhs = out.transpose();
            entry_servt = system.inv().mult(rhs);

            // Clear entries up to eshift
            for (int i = 0; i < lqn.eshift; i++) {
                entry_servt.set(i, 0, 0);
            }

            // Propagate forwarding calls: add target entry's service time to source entry
            for (int cidx = 1; cidx <= lqn.ncalls; cidx++) {
                if (lqn.calltype.get(cidx) == CallType.FWD) {
                    int source_eidx = (int) lqn.callpair.get(cidx, 1);
                    int target_eidx = (int) lqn.callpair.get(cidx, 2);
                    double fwd_prob = this.lqn.callproc_mean.getOrDefault(cidx, 1.0);
                    entry_servt.set(source_eidx - 1, 0, entry_servt.get(source_eidx - 1, 0) + fwd_prob * entry_servt.get(target_eidx - 1, 0));
                }
            }

            // Update servt for entries
            for (int i = lqn.eshift; i < lqn.eshift + lqn.nentries; i++) {
                this.servt.set(i, entry_servt.get(i, 0));
            }

            // Clear activities after ashift
            for (int i = lqn.ashift; i < entry_servt.getNumRows(); i++) {
                entry_servt.set(i, 0, 0);
            }

            // Compute entry-level residt using servtmatrix and activity residt
            Matrix residt_out = new Matrix(1, lqn.nidx + lqn.ncalls + 1, lqn.nidx + lqn.ncalls);
            Matrix.concatColumns(this.residt, this.callresidt, residt_out);
            Matrix entry_residt = new Matrix(this.servtmatrix.getNumRows(), 1, this.servtmatrix.getNumRows());
            this.servtmatrix.mult(residt_out.transpose(), entry_residt);
            for (int i = 0; i < lqn.eshift; i++) {
                entry_residt.set(i, 0, 0);
            }

            // Scale entry residt (and servt) by task/entry throughput ratio
            for (int eidx = lqn.eshift + 1; eidx < lqn.eshift + lqn.nentries + 1; eidx++) {
                int tidx = (int) lqn.parent.get(0, eidx);
                if (tidx < 0) continue;
                int hidx = (int) lqn.parent.get(0, tidx);
                if (this.ignore.get(tidx) != 0 || this.ignore.get(hidx) != 0) continue;

                boolean hasSyncCallers = false;
                for (int ii = 1; ii <= lqn.nidx; ii++) {
                    if (lqn.issynccaller.get(ii, eidx) > 0) {
                        hasSyncCallers = true;
                        break;
                    }
                }

                if (hasSyncCallers) {
                    List<Integer> tidxclass = new ArrayList<Integer>();
                    List<Integer> eidxclass = new ArrayList<Integer>();
                    for (int ii = 1; ii <= ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getTasks().size(); ii++) {
                        if (tidx == ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getTasks().get(ii)[1]) {
                            if (ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getTasks().get(ii)[0] != null)
                                tidxclass.add(ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getTasks().get(ii)[0]);
                        }
                    }
                    for (int ii = 1; ii <= ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getEntries().size(); ii++) {
                        if (eidx == ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getEntries().get(ii)[1]) {
                            eidxclass.add(ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getEntries().get(ii)[0]);
                        }
                    }
                    double task_tput = 0;
                    double entry_tput = 0;
                    for (int ii = 0; ii < tidxclass.size(); ii++) {
                        task_tput += this.results.get(results.size()).get(this.idxhash.get(hidx).intValue() - 1).TN.get(ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getClientIdx() - 1, tidxclass.get(ii) - 1);
                    }
                    for (int ii = 0; ii < eidxclass.size(); ii++) {
                        entry_tput += this.results.get(results.size()).get(this.idxhash.get(hidx).intValue() - 1).TN.get(ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getClientIdx() - 1, eidxclass.get(ii) - 1);
                    }
                    this.servt.set(eidx - 1, entry_servt.get(eidx - 1, 0) * task_tput / FastMath.max(GlobalConstants.Zero, entry_tput));
                    this.residt.set(eidx - 1, entry_residt.get(eidx - 1, 0) * task_tput / FastMath.max(GlobalConstants.Zero, entry_tput));
                } else {
                    this.residt.set(eidx - 1, entry_residt.get(eidx - 1, 0));
                }
            }

            // Update servtproc for entries based on call classes
            for (int r = 1; r < this.call_classes_updmap.getNumRows(); r++) {
                int cidx = (int) this.call_classes_updmap.get(r, 2);    // call
                int eidx = (int) lqn.callpair.get(cidx, 2);             // entry index (1-indexed)
                if (this.call_classes_updmap.get(r, 3) > 1) {
                    this.servtproc.put(eidx, Exp.fitMean(this.servt.get(eidx - 1)));
                }
            }

            // Determine call response times processes
            for (int r = 1; r < this.call_classes_updmap.getNumRows(); r++) {
                int cidx = (int) this.call_classes_updmap.get(r, 2);    // call
                int eidx = (int) lqn.callpair.get(cidx, 2);             // entry index (1-indexed)
                if (this.call_classes_updmap.get(r, 3) > 1) {
                    if (it == 1) {
                        // Note that respt is per visit, so number of calls is 1
                        this.callservt.set(cidx - 1, this.servt.get(eidx - 1));
                        this.callservtproc.put(cidx, this.servtproc.get(eidx));
                    } else {
                        // Note that respt is per visit, so number of calls is 1
                        this.callservtproc.put(cidx, Exp.fitMean(this.callservt.get(cidx - 1)));
                    }
                }
            }
        } else {
            // ===== POST-CONVERGENCE: Full CDF-based 3-moment APH fitting =====

            // Initialize CDF storage
            this.servtcdf = new HashMap<Integer, Matrix>();
            Map<Integer, DistributionResult> repo = new HashMap<Integer, DistributionResult>();

            // First obtain servt of activities at hostlayers
            this.servt = new Matrix(1, lqn.nidx, lqn.nidx);
            this.residt = new Matrix(1, lqn.nidx, lqn.nidx);
            for (int r = 1; r < this.servt_classes_updmap.getNumRows(); r++) {
                int idx = (int) this.servt_classes_updmap.get(r, 1);
                int aidx = (int) this.servt_classes_updmap.get(r, 2);
                int nodeidx = (int) this.servt_classes_updmap.get(r, 3);
                int classidx = (int) this.servt_classes_updmap.get(r, 4);

                this.tput.set(aidx - 1, this.results.get(results.size()).get(this.idxhash.get(idx).intValue() - 1).TN.get(nodeidx - 1, classidx - 1));

                // Compute residt from QN/TN_ref (matching updateMetricsDefault)
                int layerIdx_0 = this.idxhash.get(idx).intValue() - 1;
                NetworkStruct layerSn = this.ensemble[layerIdx_0].getStruct(false);
                int chainIdx = -1;
                if (layerSn != null && layerSn.chains != null && !layerSn.chains.isEmpty()) {
                    for (int ch = 0; ch < layerSn.nchains; ch++) {
                        if (layerSn.chains.get(ch, classidx - 1) > 0) {
                            chainIdx = ch;
                            break;
                        }
                    }
                }
                int refclass_c = (chainIdx >= 0 && layerSn.refclass != null) ? (int) layerSn.refclass.get(chainIdx) : classidx - 1;
                int refstat_k = (layerSn != null && layerSn.refstat != null) ? (int) layerSn.refstat.get(classidx - 1) : 0;
                double TN_ref = this.results.get(results.size()).get(layerIdx_0).TN.get(refstat_k, refclass_c);
                if (TN_ref > GlobalConstants.FineTol) {
                    this.residt.set(aidx - 1, this.results.get(results.size()).get(layerIdx_0).QN.get(nodeidx - 1, classidx - 1) / TN_ref);
                } else {
                    this.residt.set(aidx - 1, this.results.get(results.size()).get(layerIdx_0).WN.get(nodeidx - 1, classidx - 1));
                }

                int submodelidx = this.idxhash.get(idx).intValue();
                if (!repo.containsKey(submodelidx)) {
                    // Try SolverFluid first for actual response time CDFs with
                    // higher-moment information; fall back to layer solver's
                    // exponential approximation if Fluid fails on the layer model
                    try {
                        SolverOptions fluidOpts = SolverFluid.defaultOptions();
                        fluidOpts.iter_max = 1000;
                        SolverFluid fluidSolver = new SolverFluid(ensemble[submodelidx - 1], fluidOpts);
                        repo.put(submodelidx, fluidSolver.getCdfRespT());
                    } catch (Exception e) {
                        try {
                            DistributionResult cdfResult = this.solvers[submodelidx - 1].getCdfRespT();
                            repo.put(submodelidx, cdfResult);
                        } catch (Exception e2) {
                            // Skip if cannot get CDF
                            continue;
                        }
                    }
                }

                // Store CDF for this activity
                DistributionResult cdfResult = repo.get(submodelidx);
                if (cdfResult != null && cdfResult.cdfData != null &&
                    nodeidx - 1 < cdfResult.cdfData.size() &&
                    classidx - 1 < cdfResult.cdfData.get(nodeidx - 1).size()) {
                    this.servtcdf.put(aidx, cdfResult.cdfData.get(nodeidx - 1).get(classidx - 1));
                }
            }

            // Initialize callservtcdf
            this.callservtcdf = new HashMap<Integer, Matrix>();
            this.callservt = new Matrix(1, lqn.ncalls, lqn.ncalls);
            this.callresidt = new Matrix(1, lqn.ncalls, lqn.ncalls);

            // Estimate call response times at hostlayers
            for (int r = 1; r < this.call_classes_updmap.getNumRows(); r++) {
                int idx = (int) this.call_classes_updmap.get(r, 1);
                int cidx = (int) this.call_classes_updmap.get(r, 2);
                int nodeidx = (int) this.call_classes_updmap.get(r, 3);
                int classidx = (int) this.call_classes_updmap.get(r, 4);

                if (this.call_classes_updmap.get(r, 3) > 1) {
                    int submodelidx = this.idxhash.get(idx).intValue();
                    if (!repo.containsKey(submodelidx)) {
                        try {
                            SolverOptions fluidOpts2 = SolverFluid.defaultOptions();
                            fluidOpts2.iter_max = 1000;
                            SolverFluid fluidSolver = new SolverFluid(ensemble[submodelidx - 1], fluidOpts2);
                            repo.put(submodelidx, fluidSolver.getCdfRespT());
                        } catch (Exception e) {
                            try {
                                DistributionResult cdfResult = this.solvers[submodelidx - 1].getCdfRespT();
                                repo.put(submodelidx, cdfResult);
                            } catch (Exception e2) {
                                continue;
                            }
                        }
                    }

                    try {
                        DistributionResult cdfResult = repo.get(submodelidx);
                        if (cdfResult != null && cdfResult.cdfData != null &&
                            nodeidx - 1 < cdfResult.cdfData.size() &&
                            classidx - 1 < cdfResult.cdfData.get(nodeidx - 1).size()) {
                            this.callservtcdf.put(cidx, cdfResult.cdfData.get(nodeidx - 1).get(classidx - 1));
                        }
                    } catch (Exception e) {
                        // Fallback: create default CDF matrix
                        Matrix defaultCdf = new Matrix(3, 2);
                        defaultCdf.set(0, 0, 0); defaultCdf.set(0, 1, 0);
                        defaultCdf.set(1, 0, 0.5); defaultCdf.set(1, 1, 0);
                        defaultCdf.set(2, 0, 1); defaultCdf.set(2, 1, 0);
                        this.callservtcdf.put(cidx, defaultCdf);
                    }

                    // Also set callresidt from WN (includes visit multiplicity)
                    if (nodeidx > 1) {
                        this.callresidt.set(cidx - 1, this.results.get(results.size()).get(this.idxhash.get(idx).intValue() - 1).WN.get(nodeidx - 1, classidx - 1));
                    }
                }
            }

            // Build combined CDF map (combining servtcdf and callservtcdf)
            Map<Integer, Matrix> cdf = new HashMap<Integer, Matrix>();
            for (Map.Entry<Integer, Matrix> entry : this.servtcdf.entrySet()) {
                cdf.put(entry.getKey(), entry.getValue());
            }
            for (Map.Entry<Integer, Matrix> entry : this.callservtcdf.entrySet()) {
                cdf.put(lqn.nidx + entry.getKey(), entry.getValue());
            }

            // Resolve entry service times by summing contributions
            Matrix identity = Matrix.eye(lqn.nidx + lqn.ncalls);
            Matrix system = identity.sub(this.servtmatrix);
            Matrix matrix = system.inv();  // (I - servtmatrix)^(-1)

            // Process each entry
            for (int i = 1; i <= lqn.nentries; i++) {
                int eidx = lqn.eshift + i;

                // Find contributing indices (where matrix(eidx,:) > 0)
                List<Integer> convolidx = new ArrayList<Integer>();
                for (int j = 0; j < matrix.getNumCols(); j++) {
                    if (matrix.get(eidx - 1, j) > 0) {
                        // Skip entries (indices <= eshift + nentries)
                        if (j + 1 > lqn.eshift + lqn.nentries) {
                            convolidx.add(j + 1);  // 1-indexed
                        }
                    }
                }

                // Build APH convolution list
                List<Pair<Matrix, Matrix>> paramList = new ArrayList<Pair<Matrix, Matrix>>();

                for (int fitidx : convolidx) {
                    Matrix cdfMatrix = cdf.get(fitidx);
                    if (cdfMatrix == null || cdfMatrix.isEmpty()) {
                        continue;
                    }

                    // Get moments from empirical CDF
                    EmpiricalCDF empiricalCdf = new EmpiricalCDF(cdfMatrix);
                    double[] moments = empiricalCdf.getMoments();
                    double m1 = moments[0];
                    double m2 = moments[1];
                    double m3 = moments[2];

                    if (m1 > GlobalConstants.CoarseTol) {
                        // Fit APH from raw moments
                        APH fitdist = APH.fitRawMoments(m1, m2, m3);
                        Matrix alpha = fitdist.getInitProb();
                        Matrix T = (Matrix) fitdist.getParam(3).getValue();

                        // For call indices, multiply repetitions by mean number of calls
                        // servtmatrix has 1.0 for calls, but we need callproc_mean repetitions
                        double repetitions = matrix.get(eidx - 1, fitidx - 1);
                        if (fitidx > lqn.nidx) {
                            int cidx_local = fitidx - lqn.nidx;
                            repetitions = repetitions * this.lqn.callproc_mean.getOrDefault(cidx_local, 1.0);
                        }
                        int integerRepetitions = (int) Math.floor(repetitions);
                        double fractionalPart = repetitions - integerRepetitions;

                        if (fractionalPart == 0) {
                            // Only integer repetitions
                            for (int rep = 0; rep < integerRepetitions; rep++) {
                                paramList.add(new Pair<Matrix, Matrix>(alpha, T));
                            }
                        } else if (integerRepetitions > 0 && fractionalPart > 0) {
                            // Integer repetitions + fractional part
                            for (int rep = 0; rep < integerRepetitions; rep++) {
                                paramList.add(new Pair<Matrix, Matrix>(alpha, T));
                            }
                            // Add fractional part using branch structure (pattern=3)
                            APH zeroDist = APH.fitMeanAndSCV(GlobalConstants.FineTol, 0.99);
                            Matrix zeroAlpha = zeroDist.getInitProb();
                            Matrix zeroT = (Matrix) zeroDist.getParam(3).getValue();

                            kotlin.Pair<Matrix, Matrix> branchResult = aph_simplify(
                                alpha, T, zeroAlpha, zeroT,
                                fractionalPart, 1.0 - fractionalPart, 3);
                            paramList.add(new Pair<Matrix, Matrix>(branchResult.getFirst(), branchResult.getSecond()));
                        } else {
                            // Only fractional part
                            APH zeroDist = APH.fitMeanAndSCV(GlobalConstants.FineTol, 0.99);
                            Matrix zeroAlpha = zeroDist.getInitProb();
                            Matrix zeroT = (Matrix) zeroDist.getParam(3).getValue();

                            kotlin.Pair<Matrix, Matrix> branchResult = aph_simplify(
                                alpha, T, zeroAlpha, zeroT,
                                fractionalPart, 1.0 - fractionalPart, 3);
                            paramList.add(new Pair<Matrix, Matrix>(branchResult.getFirst(), branchResult.getSecond()));
                        }

                        // Update servtproc and callservtproc based on fitidx
                        if (fitidx <= lqn.nidx) {
                            this.servtproc.put(fitidx, Exp.fitMean(m1));
                            this.servt.set(fitidx - 1, m1);
                        } else {
                            this.callservtproc.put(fitidx - lqn.nidx, Exp.fitMean(m1));
                            this.callservt.set(fitidx - lqn.nidx - 1, m1);
                        }
                    }
                }

                // Convolve all contributions
                if (paramList.isEmpty()) {
                    this.servt.set(eidx - 1, 0);
                } else {
                    // Convert to Kotlin Pair list for aph_convseq
                    List<kotlin.Pair<Matrix, Matrix>> kotlinParamList = new ArrayList<kotlin.Pair<Matrix, Matrix>>();
                    for (Pair<Matrix, Matrix> p : paramList) {
                        kotlinParamList.add(new kotlin.Pair<Matrix, Matrix>(p.getLeft(), p.getRight()));
                    }

                    kotlin.Pair<Matrix, Matrix> convResult = aph_convseq(kotlinParamList);
                    APH entryDist = new APH(convResult.getFirst(), convResult.getSecond());

                    // Store results
                    int entryIndex = eidx - (lqn.nhosts + lqn.ntasks);
                    this.entryproc.put(entryIndex, entryDist);
                    this.servt.set(eidx - 1, entryDist.getMean());
                    this.servtproc.put(eidx, Exp.fitMean(this.servt.get(eidx - 1)));

                    // Evaluate and store CDF
                    Matrix evalCdf = entryDist.evalCDFMatrix();
                    this.entrycdfrespt.set(entryIndex - 1, 0, evalCdf.get(0, 0));
                }
            }

            // Propagate forwarding calls: add target entry's service time to source entry
            for (int cidx = 1; cidx <= lqn.ncalls; cidx++) {
                if (lqn.calltype.get(cidx) == CallType.FWD) {
                    int source_eidx = (int) lqn.callpair.get(cidx, 1);
                    int target_eidx = (int) lqn.callpair.get(cidx, 2);
                    double fwd_prob = this.lqn.callproc_mean.getOrDefault(cidx, 1.0);
                    this.servt.set(source_eidx - 1, this.servt.get(source_eidx - 1) + fwd_prob * this.servt.get(target_eidx - 1));
                    this.servtproc.put(source_eidx, Exp.fitMean(this.servt.get(source_eidx - 1)));
                }
            }

            // Compute entry-level residt using servtmatrix and activity residt
            Matrix residt_out_pc = new Matrix(1, lqn.nidx + lqn.ncalls + 1, lqn.nidx + lqn.ncalls);
            Matrix.concatColumns(this.residt, this.callresidt, residt_out_pc);
            Matrix entry_residt_pc = new Matrix(this.servtmatrix.getNumRows(), 1, this.servtmatrix.getNumRows());
            this.servtmatrix.mult(residt_out_pc.transpose(), entry_residt_pc);
            for (int i = 0; i < lqn.eshift; i++) {
                entry_residt_pc.set(i, 0, 0);
            }

            // Scale entry residt by task/entry throughput ratio
            for (int eidx = lqn.eshift + 1; eidx < lqn.eshift + lqn.nentries + 1; eidx++) {
                int tidx = (int) lqn.parent.get(0, eidx);
                if (tidx < 0) continue;
                int hidx = (int) lqn.parent.get(0, tidx);
                if (this.ignore.get(tidx) != 0 || this.ignore.get(hidx) != 0) continue;

                boolean hasSyncCallers = false;
                for (int ii = 1; ii <= lqn.nidx; ii++) {
                    if (lqn.issynccaller.get(ii, eidx) > 0) {
                        hasSyncCallers = true;
                        break;
                    }
                }

                if (hasSyncCallers) {
                    List<Integer> tidxclass = new ArrayList<Integer>();
                    List<Integer> eidxclass = new ArrayList<Integer>();
                    for (int ii = 1; ii <= ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getTasks().size(); ii++) {
                        if (tidx == ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getTasks().get(ii)[1]) {
                            if (ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getTasks().get(ii)[0] != null)
                                tidxclass.add(ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getTasks().get(ii)[0]);
                        }
                    }
                    for (int ii = 1; ii <= ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getEntries().size(); ii++) {
                        if (eidx == ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getEntries().get(ii)[1]) {
                            eidxclass.add(ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getEntries().get(ii)[0]);
                        }
                    }
                    double task_tput = 0;
                    double entry_tput = 0;
                    for (int ii = 0; ii < tidxclass.size(); ii++) {
                        task_tput += this.results.get(results.size()).get(this.idxhash.get(hidx).intValue() - 1).TN.get(ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getClientIdx() - 1, tidxclass.get(ii) - 1);
                    }
                    for (int ii = 0; ii < eidxclass.size(); ii++) {
                        entry_tput += this.results.get(results.size()).get(this.idxhash.get(hidx).intValue() - 1).TN.get(ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getClientIdx() - 1, eidxclass.get(ii) - 1);
                    }
                    this.residt.set(eidx - 1, entry_residt_pc.get(eidx - 1, 0) * task_tput / FastMath.max(GlobalConstants.Zero, entry_tput));
                } else {
                    this.residt.set(eidx - 1, entry_residt_pc.get(eidx - 1, 0));
                }
            }

            // Determine call response times processes (final loop)
            for (int r = 1; r < this.call_classes_updmap.getNumRows(); r++) {
                int cidx = (int) this.call_classes_updmap.get(r, 2);
                int eidx = (int) lqn.callpair.get(cidx, 2);
                if (this.call_classes_updmap.get(r, 3) > 1) {
                    if (it == 1) {
                        this.callservt.set(cidx - 1, this.servt.get(eidx - 1));
                        this.callservtproc.put(cidx, Exp.fitMean(this.servt.get(eidx - 1)));
                    }
                }
            }
        }
    }

    public void updatePopulations(int it) {
        //task15: updatePopulations function to be written
        LayeredNetworkStruct lqn = this.lqn;
        Matrix call_mult_count = this.njobsorig;
        // interlock scaling factors (not 0-padded)
        Matrix ilscaling = Matrix.ones(lqn.nhosts + lqn.ntasks, lqn.nhosts + lqn.ntasks);
        double minremote;
        for (int h = 1; h <= lqn.nhosts; h++) {
            minremote = Integer.MAX_VALUE;
            for (int hops = 1; hops <= this.nlayers; hops++) {
                int hidx = h;
                ilscaling.set(hidx - 1, 1.0);
                if (lqn.isref.get(hidx) == 0) {
                    //the following are remote (indirect) callers that are certain to be callers
                    //of task t, hence if they have multiplicity m ten task t cannot have as
                    //a matter of fact multiplicity more than m
                    List<Integer> callers = lqn.tasksof.get(hidx);
                    // caller_conn_components = lqn.conntasks(callers-lqn.tshift);
                    List<Integer> caller_conn_components = new ArrayList<>();
                    for (int caller : callers) {
                        caller_conn_components.add((int) lqn.conntasks.get(0, caller - lqn.tshift - 1));
                    }
                    int multcallers = 0;
                    for (int i : callers) {
                        multcallers += call_mult_count.get(i, hidx);
                    }
                    Matrix indirect_callers = this.ptaskcallers_step.get(hops).getRow(hidx - 1).find();
                    double multremote = 0;
                    for (Double remidxDouble : indirect_callers.toList1D()) {
                        // first we consider the update where the remote caller is an infinite server
                        // but since the ref task has finite multiplicity it is treated
                        // as normal

                        int remidx = remidxDouble.intValue() + 1;
                        if (lqn.sched.get(remidx) == SchedStrategy.INF && lqn.isref.get(remidx) == 0) {
                            multremote = Integer.MAX_VALUE;
                        } else {
                            // now we multiply the probability that a request to hidx
                            // originates from remidx
                            multremote += this.ptaskcallers_step.get(hops).get(hidx - 1, remidx - 1) * call_mult_count.get(remidx, hidx);
                        }
                    }
                    if ((multcallers > multremote && multremote > GlobalConstants.CoarseTol) && (multremote != Integer.MAX_VALUE) && multremote < minremote) {
                        minremote = multremote;
                        // we spread the scaling proportionally to the direct caller probabilities
                        List<Double> caller_spreading_ratio = new ArrayList<Double>();
                        for (int i = 0; i < callers.size(); i++) {
                            caller_spreading_ratio.add(i, this.ptaskcallers.get(hidx, callers.get(i)));
                        }
                        List<Integer> uniqueCallerConnComponents = new ArrayList<>(new HashSet<>(caller_conn_components));
                        for (int u : uniqueCallerConnComponents) {
                            double caller_ratio_sum = 0;
                            for (int i = 0; i < callers.size(); i++) {
                                if (caller_conn_components.get(i) == u) {
                                    caller_ratio_sum += caller_spreading_ratio.get(i);
                                }
                            }
                            for (int i = 0; i < callers.size(); i++) {
                                if (caller_conn_components.get(i) == u) {
                                    caller_spreading_ratio.set(i, caller_spreading_ratio.get(i) / caller_ratio_sum);
                                }
                            }
                        }
                        for (int k = 0; k < callers.size(); k++) {
                            int c = callers.get(k);
                            double num = FastMath.min(1, multremote / (double) multcallers * caller_spreading_ratio.get(k));
                            ilscaling.set(c - 1, hidx - 1, num);
                        }
                    }
                }
            }
        }


        int maxhops = this.nlayers;
        for (int t = 1; t <= lqn.ntasks; t++) {
            minremote = Integer.MAX_VALUE;
            for (int hops = 1; hops <= maxhops; hops++) { // hops <= this.nlayers
                int tidx = lqn.tshift + t;
                if (lqn.isref.get(tidx) == 0) {
                    // the following are remote (indirect) callers that certain to be
                    // callers of task t, hence if they have multiplicity m then task t
                    // cannot have as a matter of fact multiplicity more than m
                    // [calling_idx, called_entries] = find(lqn.iscaller(:, lqn.entriesof{tidx}));
                    List<Integer> calling_idx = new ArrayList<>();
                    for (int eidx : lqn.entriesof.get(tidx)) {
                        for (int i = 0; i < lqn.iscaller.getNumRows(); i++) {
                            if (lqn.iscaller.get(i, eidx) != 0) {
                                calling_idx.add(i);
                            }
                        }
                    }
                    // callers = intersect(lqn.tshift+(1:lqn.ntasks), unique(calling_idx)');
                    List<Integer> taskRange = new ArrayList<>();
                    for (int i = lqn.tshift + 1; i <= lqn.tshift + lqn.ntasks; i++) {
                        taskRange.add(i);
                    }
                    List<Integer> callers = new ArrayList<>();
                    for (int idx : calling_idx) {
                        if (taskRange.contains(idx) && !callers.contains(idx)) {
                            callers.add(idx);
                        }
                    }
                    // caller_conn_components = lqn.conntasks(callers-lqn.tshift);
                    List<Integer> caller_conn_components = new ArrayList<>();
                    for (int caller : callers) {
                        caller_conn_components.add((int) lqn.conntasks.get(0, caller - lqn.tshift - 1));
                    }
                    int multcallers = 0;
                    for (int i : callers) {
                        multcallers += this.njobsorig.get(i, tidx);
                    }
                    Matrix rowhidx = new Matrix(1, this.ptaskcallers_step.get(hops).getNumCols(), this.ptaskcallers_step.get(hops).getNumCols());
                    Matrix.extractRows(this.ptaskcallers_step.get(hops), tidx - 1, tidx, rowhidx);
                    Matrix indirect_callers = rowhidx.find(); // not 0-padded
                    double multremote = 0;
                    for (Double remidxDouble : indirect_callers.toList1D()) {
                        int remidx = remidxDouble.intValue() + 1;
                        if (lqn.sched.get(remidx) == SchedStrategy.INF) {
                            multremote = Integer.MAX_VALUE;
                        } else {
                            multremote += this.ptaskcallers_step.get(hops).get(tidx - 1, remidx - 1) * call_mult_count.get(remidx, tidx);
                        }
                    }
                    if ((multcallers > multremote && multremote > GlobalConstants.CoarseTol) && (multremote != Integer.MAX_VALUE) && multremote < minremote) {
                        minremote = multremote;
                        // we spread the scaling proportionally to the direct caller probabilities
                        List<Double> caller_spreading_ratio = new ArrayList<Double>();
                        for (int i = 0; i < callers.size(); i++) {
                            caller_spreading_ratio.add(i, this.ptaskcallers.get(tidx, callers.get(i)));
                        }
                        List<Integer> uniqueCallerConnComponents = new ArrayList<>(new HashSet<>(caller_conn_components));
                        for (int u : uniqueCallerConnComponents) {
                            double caller_ratio_sum = 0;
                            for (int i = 0; i < callers.size(); i++) {
                                if (caller_conn_components.get(i) == u) {
                                    caller_ratio_sum += caller_spreading_ratio.get(i);
                                }
                            }
                            for (int i = 0; i < callers.size(); i++) {
                                if (caller_conn_components.get(i) == u) {
                                    caller_spreading_ratio.set(i, caller_spreading_ratio.get(i) / caller_ratio_sum);
                                }
                            }
                        }
                        for (int k = 0; k < callers.size(); k++) {
                            int c = callers.get(k);
                            double num = FastMath.min(1, multremote / (double) multcallers * caller_spreading_ratio.get(k));
                            ilscaling.set(c - 1, tidx - 1, num);
                        }
                    }
                }
            }
        }
        // this.ilscaling starting from (0,0)
        // this.njobs starting from (1,1)
        this.ilscaling = ilscaling.copy();
        for (int i = 0; i < this.ilscaling.getNumRows(); i++) {
            for (int j = 0; j < this.ilscaling.getNumCols(); j++) {
                this.njobs.set(i + 1, j + 1, this.njobsorig.get(i + 1, j + 1) * this.ilscaling.get(i, j));
            }
        }
    }

    public void updateRoutingProbabilities(int it) {
        int map_length = 0;
        if (unique_route_prob_updmap.getNumRows() != 0 && unique_route_prob_updmap.getNumCols() != 0) {
            map_length = unique_route_prob_updmap.length() - 1;
        }

        for (int u = 1; u <= map_length; u++) {
            int idx;
            if (it != 0) { // same as mod(it,0)
                idx = (int) this.unique_route_prob_updmap.get(u);
            } else {
                idx = (int) this.unique_route_prob_updmap.get(this.unique_route_prob_updmap.length() - u + 1);
            }
            boolean idx_updated = false;
            java.util.Set<String> updatedArcs = new java.util.HashSet<>();

            Network nt = this.ensemble[this.idxhash.get(idx).intValue() - 1];
            List<jline.lang.nodes.Node> ntNodes = nt.getNodes();
            // MATLAB: P = self.ensemble{self.idxhash(idx)}.getLinkedRoutingMatrix;
            // Get existing routing matrix instead of creating new empty one
            Map<JobClass, Map<JobClass, Matrix>> rtorig = nt.getLinkedRoutingMatrix();

            Matrix tmp_rpu = Matrix.extractColumn(this.route_prob_updmap, 1, null);
            Matrix tmp_rpu_find = tmp_rpu.countEachRow(idx).find();
            for (int i = 0; i < tmp_rpu_find.length(); i++) {
                int r = (int) tmp_rpu_find.get(i);

                double host = this.route_prob_updmap.get(r, 1);
                double tidx_caller = this.route_prob_updmap.get(r, 2);
                double eidx = this.route_prob_updmap.get(r, 3);
                double nodefrom = this.route_prob_updmap.get(r, 4);
                double nodeto = this.route_prob_updmap.get(r, 5);
                double classidxfrom = this.route_prob_updmap.get(r, 6);
                double classidxto = this.route_prob_updmap.get(r, 7);
                // Cache routing logic - implemented following MATLAB SolverLN/updateRoutingProbabilities.m
                try {
                    int idxInt = (int) idx;
                    if (idxInt >= 0 && idxInt < this.idxhash.size()) {
                        Double idxValue = this.idxhash.get(idxInt);
                        if (idxValue != null && !Double.isNaN(idxValue)) {
                            Network network = this.ensemble[idxValue.intValue() - 1];
                            
                            // Check if idx is a cache node - following MATLAB: ~isempty(self.ensemble{self.idxhash(idx)}.items)
                            boolean isCacheNode = false;
                            for (jline.lang.nodes.Node node : network.getNodes()) {
                                if (node instanceof jline.lang.nodes.Cache) {
                                    isCacheNode = true;
                                    break;
                                }
                            }
                            
                            if (isCacheNode) { // if idx is a cache node - use host throughput
                                int hostInt = (int) host;
                                if (hostInt >= 0 && hostInt < this.idxhash.size()) {
                                    Double hostValue = this.idxhash.get(hostInt);
                                    if (hostValue != null && !Double.isNaN(hostValue)) {
                                        Network hostNetwork = this.ensemble[hostValue.intValue() - 1];
                                        // MATLAB: Xtot = sum(self.results{end,self.idxhash(host)}.TN(self.ensemble{self.idxhash(host)}.attribute.serverIdx,:))
                                        Matrix TN_copy = this.results.get(this.results.size() - 1).get(hostValue.intValue() - 1).TN;
                                        double Xtot = TN_copy.sumRows(hostNetwork.getAttribute().getServerIdx() - 1);
                                        if (Xtot > 0) {
                                            // MATLAB: hm_tput = sum(self.results{end,self.idxhash(host)}.TN(self.ensemble{self.idxhash(host)}.attribute.serverIdx,classidxto))
                                            double hm_tput = TN_copy.get(hostNetwork.getAttribute().getServerIdx() - 1, (int) classidxto - 1);
                                            double prob = hm_tput / Xtot;
                                            // MATLAB: P{classidxfrom,classidxto}(nodefrom, nodeto) = hm_tput / Xtot;
                                            // Directly modify rtorig instead of creating new RoutingMatrix
                                            JobClass classFrom = nt.getJobClassFromIndex((int) classidxfrom - 1);
                                            JobClass classTo = nt.getJobClassFromIndex((int) classidxto - 1);
                                            if (rtorig.containsKey(classFrom) && rtorig.get(classFrom).containsKey(classTo)) {
                                                Matrix routingMatrix = rtorig.get(classFrom).get(classTo);
                                                routingMatrix.set((int) nodefrom - 1, (int) nodeto - 1, prob);
                                                idx_updated = true;
                                                // Track arc for ClassSwitch node update
                                                int nodeFromIdx = (int) nodefrom - 1;
                                                int nodeToIdx = (int) nodeto - 1;
                                                if (nodeFromIdx >= 0 && nodeFromIdx < ntNodes.size() &&
                                                    nodeToIdx >= 0 && nodeToIdx < ntNodes.size()) {
                                                    String arcKey = ntNodes.get(nodeFromIdx).getName() + "_TO_" + ntNodes.get(nodeToIdx).getName();
                                                    updatedArcs.add(arcKey);
                                                }
                                            }
                                        }
                                    }
                                }
                            } else { // if idx is not a cache - use caller throughput
                                int tidxCallerInt = (int) tidx_caller;
                                if (tidxCallerInt >= 0 && tidxCallerInt < this.idxhash.size()) {
                                    Double tidxCallerValue = this.idxhash.get(tidxCallerInt);
                                    if (tidxCallerValue != null && !Double.isNaN(tidxCallerValue)) {
                                        Network callerNetwork = this.ensemble[tidxCallerValue.intValue() - 1];
                                        // MATLAB: Xtot = sum(self.results{end,self.idxhash(tidx_caller)}.TN(self.ensemble{self.idxhash(tidx_caller)}.attribute.serverIdx,:))
                                        Matrix TN_copy = this.results.get(this.results.size() - 1).get(tidxCallerValue.intValue() - 1).TN;
                                        double Xtot = TN_copy.sumRows(callerNetwork.getAttribute().getServerIdx() - 1);
                                        if (Xtot > 0) {
                                            Map<Integer, Integer[]> call_map = callerNetwork.getAttribute().getCalls();
                                            // Find the entry class for this call - MATLAB: calls(:,4) == eidx
                                            // MATLAB column 4 is 1-indexed, so in Java 0-indexed array it's index 3
                                            Integer foundEidxclass = null;
                                            for (Map.Entry<Integer, Integer[]> entry : call_map.entrySet()) {
                                                Integer[] callInfo = entry.getValue();
                                                if (callInfo.length > 3 && callInfo[3] == (int) eidx) {
                                                    // MATLAB: calls(..., 1) - column 1 = index 0
                                                    foundEidxclass = callInfo[0];
                                                    break;
                                                }
                                            }
                                            if (foundEidxclass != null) {
                                                int eidxclass = foundEidxclass;
                                                // MATLAB: entry_tput = sum(self.results{end,self.idxhash(tidx_caller)}.TN(self.ensemble{self.idxhash(tidx_caller)}.attribute.serverIdx,eidxclass))
                                                double entry_tput = TN_copy.get(callerNetwork.getAttribute().getServerIdx() - 1, eidxclass - 1);
                                                double prob = entry_tput / Xtot;
                                                // MATLAB: P{classidxfrom,classidxto}(nodefrom, nodeto) = entry_tput / Xtot;
                                                // Directly modify rtorig instead of creating new RoutingMatrix
                                                JobClass classFrom = nt.getJobClassFromIndex((int) classidxfrom - 1);
                                                JobClass classTo = nt.getJobClassFromIndex((int) classidxto - 1);
                                                if (rtorig.containsKey(classFrom) && rtorig.get(classFrom).containsKey(classTo)) {
                                                    Matrix routingMatrix = rtorig.get(classFrom).get(classTo);
                                                    // nodefrom and nodeto are stateful indices (1-indexed from MATLAB)
                                                    if (GlobalConstants.getVerbose() == VerboseLevel.DEBUG) {
                                                        System.out.println("[DEBUG] updateRoutingProbabilities: setting rtorig[" + classFrom.getName() +
                                                            "][" + classTo.getName() + "][" + ((int) nodefrom - 1) + "," + ((int) nodeto - 1) +
                                                            "] = " + prob + " (entry_tput=" + entry_tput + ", Xtot=" + Xtot + ")");
                                                    }
                                                    routingMatrix.set((int) nodefrom - 1, (int) nodeto - 1, prob);
                                                    idx_updated = true;
                                                    // Track arc for ClassSwitch node update
                                                    int nodeFromIdx = (int) nodefrom - 1;
                                                    int nodeToIdx = (int) nodeto - 1;
                                                    if (nodeFromIdx >= 0 && nodeFromIdx < ntNodes.size() &&
                                                        nodeToIdx >= 0 && nodeToIdx < ntNodes.size()) {
                                                        String arcKey = ntNodes.get(nodeFromIdx).getName() + "_TO_" + ntNodes.get(nodeToIdx).getName();
                                                        updatedArcs.add(arcKey);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                } catch (Exception e) {
                    // Fallback to avoid breaking existing functionality
                    // Cache routing logic failed, continue with normal processing
                }
            }

            // MATLAB: self.ensemble{self.idxhash(idx)}.relink(P)
            if (idx_updated) {
                nt.relinkFromRtorig(rtorig);
            }
        }
    }

    /**
     * Updates ClassSwitch nodes' switching matrices based on rtorig.
     * This is needed because updating rtorig alone doesn't propagate to ClassSwitch nodes,
     * which are used when refreshChains recomputes rtnodes from network topology.
     * Matches Python-native behavior in solver_ln.py lines 4200-4210.
     *
     * @param nt The network containing ClassSwitch nodes
     * @param updatedArcs Set of (nodefrom, nodeto) pairs that had routing updates
     */
    private void updateClassSwitchNodes(Network nt, java.util.Set<String> updatedArcs) {
        if (updatedArcs == null || updatedArcs.isEmpty()) {
            return;
        }

        NetworkStruct ntSn = nt.getStruct(false);
        if (ntSn == null || ntSn.rtorig == null || ntSn.rtorig.isEmpty()) {
            return;
        }

        List<jline.lang.nodes.Node> nodes = nt.getNodes();
        int K = ntSn.nclasses;

        // For each updated arc, find and update the corresponding ClassSwitch node
        for (String arcKey : updatedArcs) {
            String[] parts = arcKey.split("_TO_");
            if (parts.length != 2) continue;
            String nodeFromName = parts[0];
            String nodeToName = parts[1];

            // Find ClassSwitch node with name pattern CS_{nodeFrom}_to_{nodeTo}
            String csName = "CS_" + nodeFromName + "_to_" + nodeToName;
            jline.lang.nodes.ClassSwitch csNode = null;
            for (jline.lang.nodes.Node node : nodes) {
                if (node instanceof jline.lang.nodes.ClassSwitch && node.getName().equals(csName)) {
                    csNode = (jline.lang.nodes.ClassSwitch) node;
                    break;
                }
            }

            if (csNode == null) {
                continue;
            }

            // Find node indices
            int nodeFromIdx = -1;
            int nodeToIdx = -1;
            for (int i = 0; i < nodes.size(); i++) {
                if (nodes.get(i).getName().equals(nodeFromName)) {
                    nodeFromIdx = i;
                }
                if (nodes.get(i).getName().equals(nodeToName)) {
                    nodeToIdx = i;
                }
            }
            if (nodeFromIdx < 0 || nodeToIdx < 0) continue;

            // Build ClassSwitch matrix from rtorig for this arc
            jline.lang.ClassSwitchMatrix csMatrix = csNode.initClassSwitchMatrix();

            // First, compute row sums for normalization
            double[] rowSums = new double[K];
            for (int r = 0; r < K; r++) {
                JobClass fromClass = nt.getJobClassFromIndex(r);
                if (fromClass == null || !ntSn.rtorig.containsKey(fromClass)) continue;
                Map<JobClass, Matrix> fromMap = ntSn.rtorig.get(fromClass);
                for (int s = 0; s < K; s++) {
                    JobClass toClass = nt.getJobClassFromIndex(s);
                    if (toClass == null || !fromMap.containsKey(toClass)) continue;
                    Matrix routeMat = fromMap.get(toClass);
                    if (routeMat == null || nodeFromIdx >= routeMat.getNumRows() || nodeToIdx >= routeMat.getNumCols()) continue;
                    double prob = routeMat.get(nodeFromIdx, nodeToIdx);
                    if (prob > 0) {
                        rowSums[r] += prob;
                    }
                }
            }

            // Build the ClassSwitch matrix with normalized probabilities
            for (int r = 0; r < K; r++) {
                JobClass fromClass = nt.getJobClassFromIndex(r);
                if (fromClass == null || !ntSn.rtorig.containsKey(fromClass)) continue;
                Map<JobClass, Matrix> fromMap = ntSn.rtorig.get(fromClass);
                for (int s = 0; s < K; s++) {
                    JobClass toClass = nt.getJobClassFromIndex(s);
                    if (toClass == null || !fromMap.containsKey(toClass)) continue;
                    Matrix routeMat = fromMap.get(toClass);
                    if (routeMat == null || nodeFromIdx >= routeMat.getNumRows() || nodeToIdx >= routeMat.getNumCols()) continue;
                    double prob = routeMat.get(nodeFromIdx, nodeToIdx);
                    if (prob > 0 && rowSums[r] > 0) {
                        // Normalize to ensure row sums to 1
                        csMatrix.set(r, s, prob / rowSums[r]);
                    }
                }
            }

            // Update the ClassSwitch node
            csNode.setClassSwitchingMatrix(csMatrix);
        }
    }

    public void updateThinkTimes(int it) {
        //task16:updateThinkTimes function to be written
        if (this.lqn.iscaller.getNumCols() > 0) { // ignore models without callers
            Matrix torder = new Matrix(1, lqn.ntasks + 1, lqn.ntasks);
            for (int i = 1; i < lqn.ntasks + 1; i++)
                torder.set(i, i);

            this.thinkt = new Matrix(1, this.lqn.ntasks + this.lqn.tshift, this.lqn.ntasks + this.lqn.ntasks - 1);
            this.thinktproc = new HashMap<Integer, Distribution>();

            // solve all task models
            for (int t = 1; t <= this.lqn.ntasks; t++) {
                int tidx = this.lqn.tshift + t;
                double tidx_thinktime = this.lqn.think_mean.getOrDefault(tidx, 0.0); // user specified think time
                if (!Double.isNaN(this.idxhash.get(tidx) - 1)) { // this skips all REF tasks
                    // obtain total self.tput of task t
                    // mean throughput of task t in the model where it is a server, summed across replicas

                    // we use njobs to adapt to interlocking corrections
                    double njobs = Matrix.extractRows(this.njobs, tidx, tidx + 1, null).elementMax();

                    Matrix matrixExtracted = this.results.get(this.results.size()).get(this.idxhash.get(tidx).intValue() - 1).TN;
//                    Matrix serverIdxRow = new Matrix(1, matrixExtracted.getNumCols(), matrixExtracted.getNumCols());
//                    int extractRowIndex = (int) this.results.get(this.results.size()).get(this.idxhash.get(tidx).intValue() - 1).TN.
//                            get(this.ensemble[this.idxhash.get(tidx).intValue() - 1].getAttribute().getServerIdx());
//                    Matrix.extractRows(matrixExtracted, extractRowIndex, extractRowIndex + 1, serverIdxRow);
                    // lqn.repl has padded 0 index, while tput does not.
                    this.tput.set(tidx - 1, this.lqn.repl.get(0, tidx) * matrixExtracted.sumRows(this.ensemble[this.idxhash.get(tidx).intValue() - 1].getAttribute().getServerIdx() - 1));

                    // obtain total self.utilization of task t
                    Matrix UmatrixExtracted = this.results.get(this.results.size()).get(this.idxhash.get(tidx).intValue() - 1).UN;
                    this.util.set(tidx - 1, UmatrixExtracted.sumRows(this.ensemble[this.idxhash.get(tidx).intValue() - 1].getAttribute().getServerIdx() - 1));

                    if (this.lqn.sched.get(tidx) == SchedStrategy.INF) { // first we consider the update where t is an infinite server
                        // key think time update formula for LQNs, this accounts for the fact that in LINE infinite server self.utilization is dimensionally a mean number of jobs
                        this.thinkt.set(tidx - 1, FastMath.max(GlobalConstants.Zero, (njobs - this.util.get(tidx - 1)) / this.tput.get(tidx - 1) - tidx_thinktime));
                    } else { // otherwise we consider the case where t is a regular queueing station (other than an infinite server)
                        // key think time update formula for LQNs, this accounts that in LINE self.utilization is scaled in [0,1] for all queueing stations irrespectively of the number of servers
                        this.thinkt.set(tidx - 1, FastMath.max(GlobalConstants.Zero, njobs * FastMath.abs(1 - this.util.get(tidx - 1)) / this.tput.get(tidx - 1) - tidx_thinktime));
                    }
                    // Apply under-relaxation to think time if enabled
                    double omega = this.relax_omega;
                    if (omega < 1.0 && it > 1 && !Double.isNaN(this.thinkt_prev.get(tidx - 1))) {
                        this.thinkt.set(tidx - 1, omega * this.thinkt.get(tidx - 1) + (1 - omega) * this.thinkt_prev.get(tidx - 1));
                    }
                    this.thinkt_prev.set(tidx - 1, this.thinkt.get(tidx - 1));
                    Exp exponential = Exp.fitMean(this.thinkt.get(tidx - 1) + tidx_thinktime);
                    this.thinktproc.put(tidx, exponential);
                } else { // set to zero if this is a ref task
                    this.thinkt.set(tidx - 1, GlobalConstants.FineTol);
                    this.thinktproc.put(tidx, Immediate.getInstance());
                }
            }
        }
    }

    private static SolverFactory createSolverFactory(SolverType solverType) {
        return new SolverFactory() {
            @Override
            public NetworkSolver at(Network model) {
                SolverOptions options = new SolverOptions(solverType);
                options.verbose = VerboseLevel.SILENT;
                
                switch (solverType) {
                    case AUTO:
                        return new SolverAUTO(model, options);
                    case CTMC:
                        return new SolverCTMC(model, options);
                    case FLUID:
                        return new SolverFluid(model, options);
                    case JMT:
                        return new SolverJMT(model, options);
                    case MAM:
                        return new SolverMAM(model, options);
                    case MVA:
                        return new SolverMVA(model, options);
                    case NC:
                        return new SolverNC(model, options);
                    case QNS:
                        return new SolverQNS(model, options);
                    case SSA:
                        return new SolverSSA(model, options);
                    default:
                        line_error(mfilename(new Object() {}), "Unsupported solver type: " + solverType);
                        return null; // This line will never be reached due to line_error
                }
            }
        };
    }

    private static SolverFactory createSolverFactory(SolverType solverType, SolverOptions solverOptions) {
        return new SolverFactory() {
            @Override
            public NetworkSolver at(Network model) {
                switch (solverType) {
                    case AUTO:
                        return new SolverAUTO(model, solverOptions);
                    case CTMC:
                        return new SolverCTMC(model, solverOptions);
                    case FLUID:
                        return new SolverFluid(model, solverOptions);
                    case JMT:
                        return new SolverJMT(model, solverOptions);
                    case MAM:
                        return new SolverMAM(model, solverOptions);
                    case MVA:
                        return new SolverMVA(model, solverOptions);
                    case NC:
                        return new SolverNC(model, solverOptions);
                    case QNS:
                        return new SolverQNS(model, solverOptions);
                    case SSA:
                        return new SolverSSA(model, solverOptions);
                    default:
                        line_error(mfilename(new Object() {}), "Unsupported solver type: " + solverType);
                        return null; // This line will never be reached due to line_error
                }
            }
        };
    }

    static class DefaultSolverFactory implements SolverFactory {
        SolverOptions defaultOptions;

        public DefaultSolverFactory() {
            defaultOptions = new MVAOptions();
            defaultOptions.verbose = VerboseLevel.SILENT;
            //defaultOptions.config.fork_join="ht";
        }

        public NetworkSolver at(Network model) {
            return new SolverMVA(model, defaultOptions);
        }
    }

    protected static class recurActGraphReturnType {
        public JobClass curClass;
        public int jobPos;

        public RoutingMatrix P;

        public recurActGraphReturnType(JobClass curClass, int jobPos, RoutingMatrix P) {
            this.curClass = curClass;
            this.jobPos = jobPos;
            this.P = P;
        }
    }

    /**
     * Export current solver state for continuation.
     * <p>
     * Returns a LNState object containing the current solution state,
     * which can be used to continue iteration with a different solver via setState().
     * </p>
     *
     * @return LNState object containing exported state
     */
    public LNState getState() {
        LNState state = new LNState();

        // Service/think time processes
        state.servtproc = this.servtproc != null ? new HashMap<>(this.servtproc) : null;
        state.thinktproc = this.thinktproc != null ? new HashMap<>(this.thinktproc) : null;
        state.callservtproc = this.callservtproc != null ? new HashMap<>(this.callservtproc) : null;
        state.tputproc = this.tputproc != null ? new HashMap<>(this.tputproc) : null;
        state.entryproc = this.entryproc != null ? new HashMap<>(this.entryproc) : null;

        // Performance metrics
        state.util = this.util != null ? this.util.copy() : null;
        state.tput = this.tput != null ? this.tput.copy() : null;
        state.servt = this.servt != null ? this.servt.copy() : null;
        state.residt = this.residt != null ? this.residt.copy() : null;
        state.thinkt = this.thinkt != null ? this.thinkt.copy() : null;
        state.callresidt = this.callresidt != null ? this.callresidt.copy() : null;
        state.callservt = this.callservt != null ? this.callservt.copy() : null;

        // Relaxation state
        state.relaxOmega = this.relax_omega;
        state.servtPrev = this.servt_prev != null ? this.servt_prev.copy() : null;
        state.residtPrev = this.residt_prev != null ? this.residt_prev.copy() : null;
        state.tputPrev = this.tput_prev != null ? this.tput_prev.copy() : null;
        state.thinktPrev = this.thinkt_prev != null ? this.thinkt_prev.copy() : null;
        state.callservtPrev = this.callservt_prev != null ? this.callservt_prev.copy() : null;

        // Results from last iteration
        state.results = this.results != null ? new HashMap<>(this.results) : null;

        // Interlock data
        state.njobs = this.njobs != null ? this.njobs.copy() : null;
        state.ptaskcallers = this.ptaskcallers != null ? this.ptaskcallers.copy() : null;
        state.ilscaling = this.ilscaling != null ? this.ilscaling.copy() : null;

        return state;
    }

    /**
     * Import solution state for continuation.
     * <p>
     * Initializes the solver with a previously exported state, allowing iteration
     * to continue from where a previous solver left off.
     * </p>
     *
     * @param state LNState object to import
     */
    public void setState(LNState state) {
        // Service/think time processes
        if (state.servtproc != null) this.servtproc = new HashMap<>(state.servtproc);
        if (state.thinktproc != null) this.thinktproc = new HashMap<>(state.thinktproc);
        if (state.callservtproc != null) this.callservtproc = new HashMap<>(state.callservtproc);
        if (state.tputproc != null) this.tputproc = new HashMap<>(state.tputproc);
        if (state.entryproc != null) this.entryproc = new HashMap<>(state.entryproc);

        // Performance metrics
        if (state.util != null) this.util = state.util.copy();
        if (state.tput != null) this.tput = state.tput.copy();
        if (state.servt != null) this.servt = state.servt.copy();
        if (state.residt != null) this.residt = state.residt.copy();
        if (state.thinkt != null) this.thinkt = state.thinkt.copy();
        if (state.callresidt != null) this.callresidt = state.callresidt.copy();
        if (state.callservt != null) this.callservt = state.callservt.copy();

        // Relaxation state
        this.relax_omega = state.relaxOmega;
        if (state.servtPrev != null) this.servt_prev = state.servtPrev.copy();
        if (state.residtPrev != null) this.residt_prev = state.residtPrev.copy();
        if (state.tputPrev != null) this.tput_prev = state.tputPrev.copy();
        if (state.thinktPrev != null) this.thinkt_prev = state.thinktPrev.copy();
        if (state.callservtPrev != null) this.callservt_prev = state.callservtPrev.copy();

        // Results
        if (state.results != null) this.results = new HashMap<>(state.results);

        // Interlock data
        if (state.njobs != null) this.njobs = state.njobs.copy();
        if (state.ptaskcallers != null) this.ptaskcallers = state.ptaskcallers.copy();
        if (state.ilscaling != null) this.ilscaling = state.ilscaling.copy();

        // Update layer models with imported state
        int it = 1;
        if (this.results != null && !this.results.isEmpty()) {
            it = this.results.size();
        }
        updateLayers(it);

        // Refresh all layer solvers with new parameters
        for (int e = 0; e < nlayers; e++) {
            ensemble[e].refreshChains(true);
            solvers[e].reset();
        }
    }

    /**
     * Change the solver for all layers.
     * <p>
     * Replaces all layer solvers with new solvers created by the given factory function.
     * This allows switching between different solving methods (e.g., from MVA to JMT)
     * while preserving the current solution state.
     * </p>
     *
     * @param newSolverFactory Factory to create new layer solvers
     */
    public void updateSolver(SolverFactory newSolverFactory) {
        this.solverFactory = newSolverFactory;

        // Replace all layer solvers
        for (int e = 0; e < nlayers; e++) {
            solvers[e] = newSolverFactory.at(ensemble[e]);
        }
    }

    /**
     * Compute overtaking probability using transient Markov chain.
     * <p>
     * This computes the probability that a new arrival to entry eidx finds
     * the server in phase-2 (post-reply processing).
     * </p>
     * <p>
     * Uses a 3-state Continuous Time Markov Chain (CTMC):
     * <ul>
     *   <li>State 0: Server idle</li>
     *   <li>State 1: Server in phase-1 (caller is blocked)</li>
     *   <li>State 2: Server in phase-2 (caller has been released)</li>
     * </ul>
     * By PASTA (Poisson Arrivals See Time Averages), the overtaking probability
     * equals the steady-state probability of being in phase-2.
     * </p>
     *
     * @param eidx Entry index
     * @return Overtaking probability (0 to 1)
     */
    public double overtakeProb(int eidx) {
        int e = eidx - lqn.eshift;
        int tidx = (int) lqn.parent.get(0, eidx);

        double S1 = this.servt_ph1.get(eidx - 1);
        double S2 = this.servt_ph2.get(eidx - 1);

        // Get throughput - use entry if available, otherwise use task
        double lambda;
        if (this.tput.get(eidx - 1) > GlobalConstants.FineTol) {
            lambda = this.tput.get(eidx - 1);
        } else if (tidx > 0 && this.tput.get(tidx - 1) > GlobalConstants.FineTol) {
            lambda = this.tput.get(tidx - 1);
        } else {
            lambda = 0;
        }

        int c = (int) lqn.mult.get(tidx);  // number of servers

        // Handle degenerate cases
        if (S2 < GlobalConstants.FineTol || lambda < GlobalConstants.FineTol || S1 < GlobalConstants.FineTol) {
            return 0;
        }

        double mu1 = 1.0 / S1;
        double mu2 = 1.0 / S2;
        double prOt;

        if (c == 1) {
            // Single server: exact CTMC solution
            // Generator matrix Q for states {idle, phase-1, phase-2}
            // Q = | -lambda   lambda      0     |
            //     |    0      -mu1       mu1    |
            //     |   mu2       0       -mu2    |
            // Solve pi*Q = 0, sum(pi) = 1 using augmented system

            Matrix Q = new Matrix(3, 3);
            Q.set(0, 0, -lambda);
            Q.set(0, 1, lambda);
            Q.set(0, 2, 0);
            Q.set(1, 0, 0);
            Q.set(1, 1, -mu1);
            Q.set(1, 2, mu1);
            Q.set(2, 0, mu2);
            Q.set(2, 1, 0);
            Q.set(2, 2, -mu2);

            // Augmented system: [Q'; 1 1 1] * pi = [0; 0; 0; 1]
            Matrix A = new Matrix(4, 3);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    A.set(i, j, Q.get(j, i));  // Q transpose
                }
            }
            A.set(3, 0, 1);
            A.set(3, 1, 1);
            A.set(3, 2, 1);

            Matrix b = new Matrix(4, 1);
            b.set(0, 0, 0);
            b.set(1, 0, 0);
            b.set(2, 0, 0);
            b.set(3, 0, 1);

            // Solve using least squares (A' * A * x = A' * b)
            Matrix At = A.transpose();
            Matrix AtA = At.mult(A);
            Matrix Atb = At.mult(b);

            // Simple 3x3 solve using Cramer's rule or direct solve
            // For numerical stability, use the closed-form solution for 3-state birth-death
            // pi0 = mu1*mu2 / (lambda*mu2 + mu1*mu2 + lambda*mu1)
            // pi1 = lambda*mu2 / (lambda*mu2 + mu1*mu2 + lambda*mu1)
            // pi2 = lambda*mu1 / (lambda*mu2 + mu1*mu2 + lambda*mu1)
            double denom = lambda * mu2 + mu1 * mu2 + lambda * mu1;
            prOt = lambda * mu1 / denom;

            // PASTA: P(find in phase-2) = pi(2)
            prOt = FastMath.max(0, FastMath.min(1, prOt));
        } else {
            // Multi-server approximation
            double rho = lambda * (S1 + S2) / c;

            if (rho >= 1) {
                // Saturated system: probability proportional to phase-2 fraction
                prOt = S2 / (S1 + S2);
            } else {
                // Probability a random server is in phase-2
                prOt = (S2 / (S1 + S2)) * rho;
            }

            // Bound the result
            prOt = FastMath.max(0, FastMath.min(1, prOt));
        }

        return prOt;
    }

    /**
     * Get the fork fanout correction factor for an activity.
     *
     * For activities that are in a fork branch (after fork, before join),
     * returns the number of parallel branches so throughput can be corrected.
     * For fork sources, join targets, and non-fork activities, returns 1.
     *
     * Without Fork/Join nodes, probabilistic routing divides throughput by
     * the number of branches. This function identifies activities that need
     * correction (multiplication by fanout) to recover the correct throughput.
     *
     * Activities needing correction:
     * 1. Fork sources - routing normalization divides their throughput
     * 2. POST_AND activities (fork branch targets)
     * 3. Activities in fork branch chains (successors of POST_AND before join)
     *
     * Activities NOT needing correction:
     * - Join targets - receive sum from all branches
     *
     * @param aidx The activity index
     * @return The fork fanout correction factor (1 if no correction needed)
     */
    private int getForkFanout(int aidx) {
        if (this.lqn == null || this.lqn.graph == null) {
            return 1;
        }

        Matrix graph = this.lqn.graph;
        if (graph.getNumRows() == 0 || graph.getNumCols() == 0) {
            return 1;
        }

        return getForkFanoutRecursive(aidx, new HashSet<Integer>());
    }

    /**
     * Recursive helper for getForkFanout.
     */
    private int getForkFanoutRecursive(int aidx, Set<Integer> visited) {
        if (visited.contains(aidx)) {
            return 1;
        }
        visited.add(aidx);

        Matrix graph = this.lqn.graph;

        // Check if this activity is a join target (has PRE_AND predecessors)
        // Join targets don't need correction - they receive from all branches
        if (isJoinTarget(aidx)) {
            return 1;
        }

        // Check if this activity is a fork source (has POST_AND successors)
        int postAndSuccessors = countPostAndSuccessors(aidx);
        if (postAndSuccessors > 1 && !isPostAndActivity(aidx)) {
            return postAndSuccessors;
        }

        // Check if this activity is POST_AND (fork branch target)
        if (isPostAndActivity(aidx)) {
            // Find predecessor (fork source) and count its POST_AND successors
            for (int predIdx = 1; predIdx < graph.getNumRows(); predIdx++) {
                if (predIdx != aidx && graph.get(predIdx, aidx) != 0) {
                    int fanout = countPostAndSuccessors(predIdx);
                    if (fanout > 1) {
                        return fanout;
                    }
                }
            }
            return 1;
        }

        // Check if this activity is in a fork branch chain (predecessor has fanout)
        for (int predIdx = 1; predIdx < graph.getNumRows(); predIdx++) {
            if (predIdx != aidx && graph.get(predIdx, aidx) != 0) {
                int predFanout = getForkFanoutRecursive(predIdx, visited);
                if (predFanout > 1) {
                    return predFanout;
                }
            }
        }

        return 1;
    }

    /**
     * Check if an activity is POST_AND (fork branch target).
     */
    private boolean isPostAndActivity(int aidx) {
        if (this.lqn.actposttype == null) {
            return false;
        }
        if (aidx < 0 || aidx >= this.lqn.actposttype.getNumCols()) {
            return false;
        }
        return this.lqn.actposttype.get(0, aidx) == ActivityPrecedenceType.ID_POST_AND;
    }

    /**
     * Check if an activity is a join target (has PRE_AND predecessors).
     */
    private boolean isJoinTarget(int aidx) {
        if (this.lqn.actpretype == null || this.lqn.graph == null) {
            return false;
        }
        Matrix graph = this.lqn.graph;

        // Check if ANY predecessor of this activity is PRE_AND
        for (int predIdx = 1; predIdx < graph.getNumRows(); predIdx++) {
            if (predIdx != aidx && graph.get(predIdx, aidx) != 0) {
                if (predIdx < this.lqn.actpretype.getNumCols() &&
                    this.lqn.actpretype.get(0, predIdx) == ActivityPrecedenceType.ID_PRE_AND) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Count POST_AND successors of an activity.
     */
    private int countPostAndSuccessors(int aidx) {
        if (this.lqn.actposttype == null || this.lqn.graph == null) {
            return 0;
        }
        Matrix graph = this.lqn.graph;
        int count = 0;
        for (int succIdx = 1; succIdx < graph.getNumCols(); succIdx++) {
            if (graph.get(aidx, succIdx) != 0 && isPostAndActivity(succIdx)) {
                count++;
            }
        }
        return count;
    }

    /**
     * State class for exporting/importing SolverLN solution state.
     */
    public static class LNState {
        public Map<Integer, Distribution> servtproc;
        public Map<Integer, Distribution> thinktproc;
        public Map<Integer, Distribution> callservtproc;
        public Map<Integer, Distribution> tputproc;
        public Map<Integer, APH> entryproc;

        public Matrix util;
        public Matrix tput;
        public Matrix servt;
        public Matrix residt;
        public Matrix thinkt;
        public Matrix callresidt;
        public Matrix callservt;

        public double relaxOmega;
        public Matrix servtPrev;
        public Matrix residtPrev;
        public Matrix tputPrev;
        public Matrix thinktPrev;
        public Matrix callservtPrev;

        public Map<Integer, Map<Integer, SolverResult>> results;

        public Matrix njobs;
        public Matrix ptaskcallers;
        public Matrix ilscaling;
    }
}


