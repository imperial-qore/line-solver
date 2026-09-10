/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ln;

import jline.VerboseLevel;
import jline.api.fj.FJ_quorum;
import jline.api.lqn.LqnPh;
import jline.lang.processes.Markovian;
import jline.lang.workflow.Workflow;
import jline.lang.processes.Geometric;
import jline.lang.*;
import jline.GlobalConstants;
import jline.lang.constant.*;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.LayeredNetworkElement;
import jline.lang.NetworkStruct;
import jline.lang.layered.LayeredNetworkStruct;
import jline.lang.state.State;
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
import jline.solvers.wrappers.jmt.SolverJMT;
import jline.solvers.mam.SolverMAM;
import jline.solvers.mva.MVAOptions;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.solvers.wrappers.qns.SolverQNS;
import jline.solvers.ssa.SolverSSA;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import jline.api.sn.SnCompatRate;

import java.util.*;

import static jline.api.mam.Aph_convseq.aph_convseq;
import static jline.api.mam.Aph_simplify.aph_simplify;
import static jline.api.mc.Dtmc_makestochastic.dtmc_makestochastic;
import jline.io.Ret.DistributionResult;
import jline.util.Pair;
import static jline.io.InputOutput.line_debug;
import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.line_warning;
import static jline.io.InputOutput.mfilename;
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
    public Matrix joint; // join times at AND-Join activities (synchronization delay)
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
    /** True per ensemble index if the layer carries an admission constraint. */
    public boolean[] layerHasRegion;
    /** Cached sn.chains of a constrained layer, per ensemble index. */
    public Matrix[] layerChains;
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
    /** method='moment3': the entry response time law tabulated as an (n x 2) [F(t), t] matrix, keyed by the 0-based entry-local index. */
    private Map<Integer, Matrix> entrycdfrespt;
    /** method='moment3': true once the moment-based entry-law pass has run. */
    private boolean momentPassDone;
    // The LQN this solver was built on. The inherited Solver.model stays null
    // (the ensemble does not exist yet when super() runs), so the layered model
    // is kept here: it is what lets the layer submodels be shared back with the
    // LayeredNetwork and what carries a warm start into the layer transients.
    private LayeredNetwork lqnModel;


    // Under-relaxation state for convergence improvement
    public double relax_omega;        // Current relaxation factor
    public List<Double> relax_err_history;  // Error history for adaptive mode
    public Matrix servt_prev;         // Previous service times for relaxation
    public Matrix residt_prev;        // Previous residence times for relaxation
    public Matrix tput_prev;          // Previous throughputs for relaxation
    public Matrix thinkt_prev;        // Previous think times for relaxation
    public java.util.Set<Integer> singleReplicaTasks; // Task indices modeled as single representative replica (fan-out)
    public Matrix callservt_prev;     // Previous call service times for relaxation
    public Matrix callresidt_prev;    // Previous call residence times for growth rate capping

    // Stochastic iteration (Robbins-Monro / Polyak-Ruppert) state, used when
    // one or more layer solvers return noisy estimates
    public String stochiterMode;      // resolved mode: "rm" | "crn" | "off"
    public boolean stochiterAuto;     // true if mode was resolved from 'auto'
    public Integer stochiterStart;    // iteration at which RM averaging started
    public boolean[] stochlayers;     // true if the layer solver is stochastic
    public Map<Integer, SolverResult> stochAvg; // Polyak-Ruppert averages of layer results
    public int stochAvgCount;         // iterations accumulated into stochAvg
    public Matrix stochServtAvg;      // Polyak-Ruppert average of the servt iterate
    public Matrix stochResidtAvg;     // Polyak-Ruppert average of the residt iterate

    // MOL (Method of Layers) properties for hierarchical iteration
    public List<Integer> hostLayerIndices;   // Indices of host (processor) layers in ensemble
    public List<Integer> taskLayerIndices;   // Indices of task layers in ensemble

    // Solver factory for update_solver functionality
    public SolverFactory solverFactory;
    /** layer solver kind when built from a SolverType; null when a bare factory was given */
    private SolverType layerSolverType;

    // Phase-2 support properties
    public boolean hasPhase2;           // Flag: model has phase-2 activities
    public Matrix servt_ph1;            // Phase-1 service time per activity (nidx x 1)
    public Matrix servt_ph2;            // Phase-2 service time per activity (nidx x 1)
    public Matrix util_ph1;             // Phase-1 utilization per entry
    public Matrix util_ph2;             // Phase-2 utilization per entry
    public Matrix prOvertake;           // Overtaking probability per entry (nentries x 1)

    // Interlock path tables of Franks (1999), Ch. 4 (built once at init)
    public double[][] il_table_all;              // (nentries x nentries) reachability, all phases
    public double[][] il_table_ph1;              // (nentries x nentries) reachability, phase-1 only
    public List<Integer>[] il_common_entries;     // common parent entry abs-indices per server
    public List<Integer>[] il_source_tasks_all;   // all-phase source tasks per server
    public List<Integer>[] il_source_tasks_ph2;   // phase-2 source tasks per server
    public double[] il_num_sources;              // total source multiplicity per server

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
        this(lqnmodel, createSolverFactory(solverType), options, solverType);
    }

    public SolverLN(LayeredNetwork lqnmodel, SolverType solverType, LNOptions lnOptions, SolverOptions solverOptions) {
        this(lqnmodel, createSolverFactory(solverType, solverOptions), lnOptions, solverType);
    }

    /**
     * Forwarding transformation of Franks (1999), Sec. 3.3.1 and Fig. 3.8:
     * each forwarding chain reachable from a synchronous call is reconnected
     * to the client that issued the original rendezvous, as a pseudo
     * rendezvous (SYNC) call whose mean is the original call mean times the
     * product of the forwarding probabilities along the path. One level of
     * servers disappears from the layering and the forwarded workload is
     * carried by ordinary SYNC call classes, so layer construction, think
     * times, populations and the interlock analysis all see plain rendezvous
     * arcs. As the thesis notes, the pseudo arcs are excluded from the slice
     * times and from the overtaking and interlock probabilities. FWD calls
     * remain in the struct but no longer contribute blocking anywhere in
     * SolverLN. Asynchronous calls into a forwarding chain are left
     * untouched, since a send-no-reply terminates the chain of blocking.
     */
    private void applyForwardingRendezvous() {
        boolean hasFwd = false;
        for (int c = 0; c < lqn.ncalls; c++) {
            if (lqn.calltype.get(c) == CallType.FWD) {
                hasFwd = true;
                break;
            }
        }
        if (!hasFwd) return;

        int ncalls0 = lqn.ncalls;
        for (int cidx = 0; cidx < ncalls0; cidx++) {
            if (lqn.calltype.get(cidx) != CallType.SYNC) continue;
            int aidx = (int) lqn.callpair.get(cidx, 0);
            int tidx = (int) lqn.parent.get(0, aidx);
            double base_mean = lqn.callproc_mean.getOrDefault(cidx, 0.0);
            if (base_mean <= 0) continue;
            // BFS through the forwarding chain of the sync target
            List<Integer> frontier = new ArrayList<>();
            List<Double> probs = new ArrayList<>();
            List<Integer> visitedE = new ArrayList<>();
            frontier.add((int) lqn.callpair.get(cidx, 1));
            probs.add(1.0);
            while (!frontier.isEmpty()) {
                int eidx = frontier.remove(0);
                double p_path = probs.remove(0);
                if (visitedE.contains(eidx)) continue;
                visitedE.add(eidx);
                for (int fcidx = 0; fcidx < ncalls0; fcidx++) {
                    if (lqn.calltype.get(fcidx) != CallType.FWD || (int) lqn.callpair.get(fcidx, 0) != eidx) continue;
                    double fprob = lqn.callproc_mean.getOrDefault(fcidx, 0.0);
                    int tgt = (int) lqn.callpair.get(fcidx, 1);
                    double pseudo_mean = base_mean * p_path * fprob;
                    int target_tidx = (int) lqn.parent.get(0, tgt);
                    if (pseudo_mean > 0 && target_tidx != tidx) {
                        // Merge into an existing SYNC call with the same
                        // (activity, target) pair if any; otherwise append a
                        // new pseudo SYNC call.
                        // Call indices are 0-based, so 0 is a REAL call and
                        // cannot double as the not-found sentinel.
                        int mrow = -1;
                        for (int scan = 0; scan < lqn.ncalls; scan++) {
                            if (lqn.calltype.get(scan) == CallType.SYNC
                                    && (int) lqn.callpair.get(scan, 0) == aidx
                                    && (int) lqn.callpair.get(scan, 1) == tgt) {
                                mrow = scan;
                                break;
                            }
                        }
                        if (mrow >= 0) {
                            double newmean = lqn.callproc_mean.get(mrow) + pseudo_mean;
                            Distribution d = LayeredNetwork.callCountDist(newmean);
                            lqn.callproc.put(mrow, d);
                            lqn.callproc_mean.put(mrow, newmean);
                            lqn.callproc_scv.put(mrow, d.getSCV());
                        } else {
                            // The appended call takes the next FREE 0-based
                            // index, which is the old count; ncalls is the
                            // count, not the last index.
                            int ncall = lqn.ncalls;
                            lqn.ncalls = ncall + 1;
                            Distribution d = LayeredNetwork.callCountDist(pseudo_mean);
                            if (ncall >= lqn.callpair.getNumRows()) {
                                Matrix newPair = new Matrix(ncall + 1, lqn.callpair.getNumCols());
                                for (int r = 0; r < lqn.callpair.getNumRows(); r++) {
                                    for (int cc = 0; cc < lqn.callpair.getNumCols(); cc++) {
                                        newPair.set(r, cc, lqn.callpair.get(r, cc));
                                    }
                                }
                                lqn.callpair = newPair;
                            }
                            lqn.calltype.put(ncall, CallType.SYNC);
                            lqn.callpair.set(ncall, 0, aidx);
                            lqn.callpair.set(ncall, 1, tgt);
                            lqn.callnames.put(ncall, lqn.names.get(aidx) + "=>" + lqn.names.get(tgt));
                            lqn.callhashnames.put(ncall, lqn.hashnames.get(aidx) + "=>" + lqn.hashnames.get(tgt));
                            lqn.callproc.put(ncall, d);
                            // Only the mean (and process object) are consumed
                            // by SolverLN; mirror the base call for the rest
                            lqn.callproc_type.put(ncall, lqn.callproc_type.get(cidx));
                            lqn.callproc_params.put(ncall, lqn.callproc_params.get(cidx));
                            lqn.callproc_mean.put(ncall, pseudo_mean);
                            lqn.callproc_scv.put(ncall, d.getSCV());
                            lqn.callproc_proc.put(ncall, lqn.callproc_proc.get(cidx));
                            if (lqn.callsof.get(aidx) == null) {
                                lqn.callsof.put(aidx, new ArrayList<Integer>());
                            }
                            lqn.callsof.get(aidx).add(ncall);
                            lqn.iscaller.set(tidx, target_tidx, 1.0);
                            lqn.iscaller.set(aidx, target_tidx, 1.0);
                            lqn.iscaller.set(tidx, tgt, 1.0);
                            lqn.iscaller.set(aidx, tgt, 1.0);
                            lqn.issynccaller.set(tidx, target_tidx, 1.0);
                            lqn.issynccaller.set(aidx, target_tidx, 1.0);
                            lqn.issynccaller.set(tidx, tgt, 1.0);
                            lqn.issynccaller.set(aidx, tgt, 1.0);
                            lqn.graph.set(aidx, tgt, 1.0);
                            lqn.taskgraph.set(tidx, target_tidx, 1.0);
                        }
                    }
                    // Follow the chain
                    if (!visitedE.contains(tgt) && !frontier.contains(tgt)) {
                        frontier.add(tgt);
                        probs.add(p_path * fprob);
                    }
                }
            }
        }
    }

    public SolverLN(LayeredNetwork lqnmodel, SolverFactory solverFactory, SolverOptions options) {
        this(lqnmodel, solverFactory, options, null);
    }

    /**
     * layerSolverType is taken before construct(): buildLayers gates routed call groups on
     * the layer solver, so it has to be readable by then.
     */
    public SolverLN(LayeredNetwork lqnmodel, SolverFactory solverFactory, SolverOptions options,
                    SolverType layerSolverType) {
        super(null, "SolverLN", options); // first argument is null as the ensemble cannot be built yet
        this.layerSolverType = layerSolverType;
        this.lqnModel = lqnmodel;
        this.lqn = lqnmodel.getStruct();
        // Rewrite forwarding chains as caller-side pseudo rendezvous calls,
        // see applyForwardingRendezvous
        applyForwardingRendezvous();

        // Detect and initialize phase-2 support
        this.hasPhase2 = false;
        if (lqn.actphase != null) {
            for (int a = 0; a < lqn.nacts; a++) {
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

        for (int i = 0; i < nlayers; i++) {
            // A setup no longer forces the MAM decomposition on the layer. The open
            // M/G/1-with-setup QBD reads the idle period from the Poisson rate 1/X,
            // and in a CLOSED layer the idle period a thread sees is the rest of the
            // cycle, 1/X - S: on lqn_setup that is 1.0 against the 2.29 the open
            // reading gives, so the thread was powered down far more often than it
            // is, and the answer landed 12.67% below LDES. The cold start is charged
            // to the ENTRY instead, with the probability that the thread was
            // actually found down -- see setupCharge.
            solvers[i] = silenced(solverFactory.at(ensemble[i]));
            assertLayerSolverSupportsModel(solvers[i], ensemble[i], i);
        }
        this.solverFactory = solverFactory; // Store for later use
        line_debug(options.verbose, String.format("LN: constructed %d layers, nhosts=%d, ntasks=%d, nentries=%d, nacts=%d",
            nlayers, lqn.nhosts, lqn.ntasks, lqn.nentries, lqn.nacts));
        for (int i = 0; i < nlayers; i++) {
            line_debug(options.verbose, String.format("LN: layer %d solver=%s, nstations=%d, nclasses=%d",
                i, solvers[i].getName(), ensemble[i].getNumberOfStations(), ensemble[i].getNumberOfClasses()));
        }
    }

    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.LN);
    }

    @Override
    public SolverResult analyze(int it, int e) {
        SolverResult result1 = new SolverResult();
        System.out.flush();
        solvers[e].getAvg();

        // A layer solver is SILENT (see silenced), so it prints no banner of its
        // own and there is nothing left here to separate: the blank line that
        // used to precede layer 1 is gone with the banner.

        result1.QN = solvers[e].result.QN;
        result1.UN = solvers[e].result.UN;
        result1.RN = solvers[e].result.RN;
        result1.TN = solvers[e].result.TN;
        result1.AN = solvers[e].result.AN;
        result1.WN = solvers[e].result.WN;
        result1.CN = solvers[e].result.CN; // not in MATLAB
        result1.XN = solvers[e].result.XN; // not in MATLAB

        // Refresh the stochastic classification from the method the layer
        // solver actually resolved at runtime (e.g. an NC layer with method
        // 'default' falling back to Monte Carlo integration); post() resets
        // the layer solvers, so this must be captured here while results are
        // still attached.
        if (it == 1 && this.stochlayers != null) {
            this.stochlayers[e] = solvers[e].isStochastic();
        }

        // Warm-start the next AMVA solve of this layer from the current
        // solution (chain-aggregated queue lengths). The layer AMVA can
        // admit multiple fixed points (e.g., multiserver FCFS layers near
        // saturation), so a cold restart may jump between solution
        // branches under infinitesimal input changes, which prevents
        // outer-loop convergence.
        if ("SolverMVA".equals(solvers[e].name) && result1.QN != null) {
            NetworkStruct sne = ensemble[e].getStruct(false);
            if (result1.QN.getNumRows() == sne.nstations && result1.QN.getNumCols() == sne.nclasses) {
                Matrix Qch = new Matrix(sne.nstations, sne.nchains);
                for (int c = 0; c < sne.nchains; c++) {
                    for (int r = 0; r < sne.nclasses; r++) {
                        if (sne.chains.get(c, r) > 0) {
                            for (int i = 0; i < sne.nstations; i++) {
                                double q = result1.QN.get(i, r);
                                if (!Double.isNaN(q) && !Double.isInfinite(q)) {
                                    Qch.set(i, c, Qch.get(i, c) + q);
                                }
                            }
                        }
                    }
                }
                solvers[e].options.init_sol = Qch;
            }
        }

        return result1;
    }

    /**
     * The method the layers were actually BUILT for: 'srvn.ph', 'srvn.cs',
     * 'flat.cs' or 'moment3'. Resolved once in buildLayers, because the alias
     * 'srvn' may fall back; every dispatch reads this and not options.method, so
     * a reconstruction can never disagree with the layers it is reading.
     */
    public String lnmethod;

    /**
     * The layering the METHOD asked for, or null when the method names none and
     * options.config.layering decides. Set once in buildLayers.
     */
    private String lnlayering;

    /**
     * Normalise a method name onto one the solver dispatches on. A method name
     * carries TWO decisions: the LAYERING, which fixes what a submodel is, and
     * the ENCODING, which fixes how an activity graph is written into it.
     *
     * 'srvn.cs' encodes the activity graph as ROUTING, 'srvn.ph' as a composed
     * phase-type server law, 'srvn' is the alias that takes 'srvn.ph' where it
     * can serve the model and 'srvn.cs' otherwise, 'flat.cs' squashes every
     * server into one submodel with the routing encoding, 'flat.ph' squashes
     * them with the composed one ('flat' is the alias of 'flat.cs' and resolves
     * unconditionally rather than probing 'flat.ph', because a model is squashed
     * in order to express what only the routing encoding carries), and 'moment3'
     * is the three-moment distribution pass over the routing layers. 'default'
     * is the srvn alias; an unrecognised method name takes 'srvn.cs'.
     *
     * @param method the requested method, may be null
     * @return one of "srvn", "srvn.ph", "srvn.cs", "flat.cs", "moment3"
     */
    public static String lnRequestedMethod(String method) {
        if (method == null || method.isEmpty()) {
            return "srvn";
        }
        String m = method.toLowerCase();
        if ("srvn.ph".equals(m) || "ph".equals(m)) {
            return "srvn.ph";
        }
        if ("srvn.cs".equals(m) || "srvncs".equals(m) || "cs".equals(m)) {
            return "srvn.cs";
        }
        if ("srvn".equals(m) || "default".equals(m) || "auto".equals(m)) {
            return "srvn";
        }
        if ("flat.cs".equals(m) || "flatcs".equals(m) || "flat".equals(m)
                || "squashed".equals(m)) {
            return "flat.cs";
        }
        if ("flat.ph".equals(m) || "flatph".equals(m) || "squashed.ph".equals(m)) {
            return "flat.ph";
        }
        if ("moment3".equals(m)) {
            return "moment3";
        }
        // An unrecognised token takes the routing encoding, which is what every
        // name other than 'moment3' resolved to before the alias existed.
        return "srvn.cs";
    }

    /**
     * Valid methods for this solver, SolverLN.m verbatim.
     *
     * Each name states the LAYERING and the ENCODING; lnRequestedMethod
     * normalises the alias spellings ('ph', 'cs', 'srvncs', 'flatcs',
     * 'squashed', 'squashed.ph') onto these, and they are left out here to keep
     * the list unambiguous, exactly as the reference does.
     *
     * @return the method names SolverLN accepts
     */
    public String[] listValidMethods() {
        return new String[]{"srvn", "srvn.ph", "srvn.cs", "flat", "flat.cs", "flat.ph",
                "moment3", "default"};
    }

    /**
     * The encoding rules the layer builders enforce at solve time, stated here so
     * that a CALLER can see them before running.
     *
     * <p>{@code srvn.ph} and {@code flat.ph} compose each entry into ONE
     * phase-type law, and several constructs have nowhere to go in that law: a
     * forwarding call whose target is not in the caller's activity graph, a routed
     * call group whose dispatch order the composition folds away, a cache task, an
     * admission constraint, a queue-dependent rate on a station the composition
     * replaces. {@code flat.ph} additionally squashes every layer into one
     * network, which per-layer state (a replica, a powered-down setup thread)
     * cannot survive.</p>
     *
     * <p>None of these is a feature name, so none can be a feature-set delta: they
     * are properties of what the METHOD does to the model. Left only in
     * {@code phFlatServerSet} and the composer they were invisible to every gate
     * above them, and {@code listValidMethods} returns the same eight names for
     * every model, so a report offered every encoding on every layered model.</p>
     *
     * <p>Phase 2 is deliberately NOT tested: that refusal reads
     * {@code this.hasPhase2}, which is built during layering rather than being a
     * property of the model, so a gate cannot ask it without doing the layering it
     * precedes. Mirrors MATLAB {@code ln_method_refusal}.</p>
     *
     * @param method the concrete method name
     * @return empty string when the method can encode this model, else the reason
     */
    @Override
    public String supportsModelMethod(String method) {
        String reason = super.supportsModelMethod(method);
        if (!reason.isEmpty()) {
            return reason;
        }
        if (lqn == null || !("srvn.ph".equals(method) || "flat.ph".equals(method))) {
            return "";
        }
        // The squashing refusals, 'flat.ph' only: each carries PER-LAYER state
        // that one submodel cannot hold.
        if ("flat.ph".equals(method)) {
            int nelem = lqn.nhosts + lqn.ntasks;
            for (int i = 0; i < nelem; i++) {
                if (lqn.repl.get(0, i) > 1) {
                    return "method='flat.ph' does not support replicated processors or tasks, "
                            + "whose replicas need a submodel each. Use method='srvn.ph'.";
                }
                if (lqn.hassetup != null && i < lqn.hassetup.getNumCols()
                        && lqn.hassetup.get(0, i) != 0) {
                    return "method='flat.ph' does not support setup tasks, whose powered-down "
                            + "threads are per-layer state. Use method='srvn.ph'.";
                }
            }
        }
        // The composed-entry-law refusals, both PH encodings.
        if (lqn.iscache != null) {
            for (int i = 0; i < lqn.iscache.getNumCols(); i++) {
                if (lqn.iscache.get(0, i) != 0) {
                    return "method='" + method + "' does not support cache tasks. "
                            + "Use method='default'.";
                }
            }
        }
        for (int cidx = 0; cidx < lqn.ncalls; cidx++) {
            if (lqn.calltype.get(cidx) == CallType.FWD) {
                return "method='" + method + "' does not support forwarding calls, whose target "
                        + "is not part of the caller's activity graph. Use method='default'.";
            }
        }
        if (lqn.hassetup != null) {
            for (int i = 0; i < lqn.hassetup.getNumCols(); i++) {
                if (lqn.hassetup.get(0, i) == 0) {
                    continue;
                }
                if (lqn.sched.get(i) == SchedStrategy.INF
                        || Double.isInfinite(lqn.mult.get(0, i))) {
                    return "method='" + method + "': task '" + lqn.names.get(i)
                            + "' declares a setup time on an infinite-server task, which holds no "
                            + "thread to power down; give it a finite multiplicity.";
                }
            }
        }
        if (lqn.callgroups != null && !lqn.callgroups.isEmpty()) {
            return "method='" + method + "' does not support routed call groups, whose dispatch "
                    + "order is a routing property. Use method='flat.cs'.";
        }
        if (lqn.lincon != null && !lqn.lincon.isEmpty()) {
            return "method='" + method + "' does not support admission constraints on a layer "
                    + "station. Use method='default'.";
        }
        return "";
    }

    /** True when the layers are the collapsed phase-type ones of method 'srvn.ph'. */
    public boolean isSrvnPH() {
        return "srvn.ph".equals(lnmethod);
    }

    /**
     * True when the layers carry the COMPOSED phase-type server law rather than
     * the routing encoding of the activity graph, under either layering. The
     * encoding, not the layering, decides which update and reconstruction passes
     * run, so every such dispatch asks this and not for one method name.
     *
     * @return true for 'srvn.ph' and for 'flat.ph'
     */
    public boolean isPHEncoding() {
        return "srvn.ph".equals(lnmethod) || "flat.ph".equals(lnmethod);
    }

    public void buildLayers() {

        // Method resolution. A method name carries both the LAYERING and the
        // ENCODING: 'srvn.ph' replaces the routing encoding of the activity graph
        // by a composed phase-type server law, 'srvn' is the alias that takes it
        // where it can serve the model and 'srvn.cs' otherwise, and 'flat.cs'
        // squashes every server into one submodel. The choice is made ONCE, here,
        // and every later dispatch reads lnmethod.
        // See _kb/06-solver-catalog.md (LN section).
        String requested = lnRequestedMethod(options == null ? null : options.method);
        // the method names the layering, so it sets it
        this.lnlayering = ("flat.cs".equals(requested) || "flat.ph".equals(requested)) ? "flat" : null;
        if ("flat.ph".equals(requested)) {
            // the squashed layering with the composed law: ONE submodel holding
            // every server, and a caller visiting each of them once per
            // invocation. The feature gate is the srvn.ph one plus the refusals a
            // single submodel carries -- see phFlatServerSet.
            this.phLawsReady = false;
            this.lnmethod = "flat.ph";
            buildLayersPH(true);
            return;
        }
        if ("srvn.ph".equals(requested) || "srvn".equals(requested)) {
            boolean hard = "srvn.ph".equals(requested);
            if (!flatServerSet().isEmpty()) {
                if (hard) {
                    throw new IllegalStateException("method='srvn.ph' requires the srvn layering, "
                            + "because it replaces each server by a submodel of its own. Use "
                            + "method='srvn.cs' for that layering.");
                }
                line_debug(options.verbose, "LN: method=srvn cannot use srvn.ph under the flat layering");
            } else {
                this.phLawsReady = false;
                if (hard || probeSrvnPH()) {
                    this.lnmethod = "srvn.ph";
                    buildLayersPH();
                    return;
                }
            }
        }
        // The label reports what was BUILT, so a model squashed through
        // options.config.layering reads back as 'flat.cs' even when no method
        // named it.
        if ("moment3".equals(requested)) {
            this.lnmethod = "moment3";
        } else {
            this.lnmethod = flatServerSet().isEmpty() ? "srvn.cs" : "flat.cs";
        }

        this.temp_ensemble = new ArrayList<>();
        this.cell_servt_classes_updmap = new HashMap<>(lqn.nhosts + lqn.ntasks);
        this.cell_call_classes_updmap = new HashMap<>(lqn.nhosts + lqn.ntasks);
        this.cell_arvproc_classes_updmap = new HashMap<>(lqn.nhosts + lqn.ntasks);
        this.cell_thinkt_classes_updmap = new HashMap<>(lqn.nhosts + lqn.ntasks);
        this.cell_route_prob_updmap = new HashMap<>(lqn.nhosts + lqn.ntasks);

        double temp_idxhash = 1;

        // see _kb/06-solver-catalog.md (LN section) for the layering taxonomy
        List<Integer> flatServers = flatServerSet();
        if (!flatServers.isEmpty()) {
            buildFlatLayer(flatServers);
            return;
        }

        // build one subnetwork for every processor
        for (int hidx = 0; hidx < lqn.nhosts; hidx++) {
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
        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            // MATLAB: ~self.ignore(tidx) & ~lqn.isref(tidx) & ~(isempty(find(self.lqn.iscaller(tidx,:), 1)) & isempty(find(self.lqn.iscaller(:,tidx), 1)))
            boolean isolated_task = true;
            for (int i = 0; i < lqn.nidx; i++) {
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
                for (int i = lqn.tshift; i < lqn.tshift + lqn.ntasks; i++) {
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

        // Share the layer submodels with the LayeredNetwork so the model and
        // this solver operate on the same layer objects. SolverENV over an LQN
        // stage relies on this: model.initFromMarginal (warm start) and the
        // aggregate getEnvStruct/getStations must reach the very layers whose
        // transient this solver later reads. Mirrors the MATLAB/Python
        // model.ensemble = solver-layers assignment.
        if (this.lqnModel != null) {
            this.lqnModel.setEnsemble(
                    new java.util.ArrayList<Network>(java.util.Arrays.asList(this.ensemble)));
        }

        // Layers carrying an admission constraint need the region wait recovered in
        // updateMetricsDefault -- see _kb/06-solver-catalog.md (LN section)
        this.layerHasRegion = new boolean[nlayers];
        this.layerChains = new Matrix[nlayers];
        for (int e = 0; e < nlayers; e++) {
            this.layerHasRegion[e] = this.ensemble[e].getRegions() != null && !this.ensemble[e].getRegions().isEmpty();
            if (this.layerHasRegion[e]) {
                // layer structure is iteration-invariant, so cache the chain matrix
                this.layerChains[e] = this.ensemble[e].getStruct().chains;
            }
        }

        // Classify layers as host (processor) or task for MOL iteration
        this.hostLayerIndices = new ArrayList<>();
        this.taskLayerIndices = new ArrayList<>();

        // Host layers: indices 1:nhosts (before idxhash remapping)
        // Note: idxhash uses 1-based indexing (NaN at index 0), so use hidx directly
        for (int hidx = 0; hidx < lqn.nhosts; hidx++) {
            if (!Double.isNaN(idxhash.get(hidx))) {
                hostLayerIndices.add(idxhash.get(hidx).intValue() - 1); // Convert to 0-indexed
            }
        }

        // Task layers: indices tshift+1:tshift+ntasks
        // Note: idxhash uses 1-based indexing (NaN at index 0), so use tidx directly
        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            if (tidx < idxhash.size() && !Double.isNaN(idxhash.get(tidx))) {
                taskLayerIndices.add(idxhash.get(tidx).intValue() - 1); // Convert to 0-indexed
            }
        }
    }

    /**
     * Processors and called tasks that become stations of the flat layer, empty
     * under the default 'srvn' layering. The features a single submodel cannot
     * carry are rejected here rather than silently dropped.
     */
    /**
     * Rejects routed call groups under any layering that cannot carry them.
     *
     * A group states the order in which one caller visits several callees. The
     * srvn layering puts every callee in a submodel of its own and replaces it,
     * in the caller's submodel, by a surrogate delay, so the callees are never
     * co-resident and no node has arcs to more than one of them: the order has
     * nowhere to be expressed and would be silently degraded to the aggregate
     * call means. The squashed layering keeps all of them as stations of one
     * model, which is what makes the strategy representable.
     */
    private void assertCallGroups() {
        if (lqn.callgroups == null || lqn.callgroups.isEmpty()) {
            return;
        }
        String layering = this.options != null && this.options.config != null
                ? this.options.config.layering : null;
        boolean flat = layering != null
                && (layering.equalsIgnoreCase("flat") || layering.equalsIgnoreCase("squashed"));
        if (!flat) {
            throw new RuntimeException(
                    "Call groups routed by a routing strategy require the squashed layering; "
                            + "set options.config.layering='flat'. Under srvn the targets never "
                            + "share a submodel, so the dispatch order cannot be represented.");
        }
        // Only a layer solver with state-dependent routing honours the strategy;
        // MVA, NC and FLD would silently return the probabilistic split instead.
        String st = this.layerSolverType != null ? this.layerSolverType.name().toUpperCase() : "";
        if (!(st.contains("CTMC") || st.contains("SSA") || st.contains("LDES") || st.contains("JMT"))) {
            throw new RuntimeException(
                    "Routed call groups need a layer solver with state-dependent routing "
                            + "(CTMC or SSA); MVA, NC and FLD would silently return the "
                            + "probabilistic split under a round-robin or JSQ label.");
        }
    }

    private List<Integer> flatServerSet() {
        assertCallGroups();
        List<Integer> servers = new ArrayList<>();
        // lnlayering is what the METHOD asked for and wins where it is set;
        // otherwise the layering is the one named in the options.
        String layering = this.lnlayering != null ? this.lnlayering
                : (this.options != null && this.options.config != null
                        ? this.options.config.layering : null);
        if (layering == null || layering.equalsIgnoreCase("srvn")) {
            return servers;
        }
        if (!layering.equalsIgnoreCase("flat") && !layering.equalsIgnoreCase("squashed")) {
            throw new RuntimeException("Unknown layering strategy '" + layering
                    + "', use 'srvn' or 'flat'.");
        }
        int nelem = lqn.nhosts + lqn.ntasks;
        for (int i = 0; i < nelem; i++) {
            if (lqn.repl.get(0, i) > 1) {
                throw new RuntimeException("Flat layering does not support replicated processors or tasks, use the default 'srvn' layering.");
            }
            if (lqn.iscache != null && lqn.iscache.get(0, i) != 0) {
                throw new RuntimeException("Flat layering does not support cache tasks, use the default 'srvn' layering.");
            }
            if (lqn.hassetup != null && lqn.hassetup.get(0, i) != 0) {
                throw new RuntimeException("Flat layering does not support setup tasks, use the default 'srvn' layering.");
            }
        }
        for (int hidx = 0; hidx < lqn.nhosts; hidx++) {
            if (this.ignore.get(hidx) == 0 && lqn.tasksof.get(hidx) != null
                    && !lqn.tasksof.get(hidx).isEmpty()) {
                servers.add(hidx);
            }
        }
        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            if (this.ignore.get(tidx) != 0 || (int) lqn.isref.get(tidx) != 0) {
                continue;
            }
            boolean isolated = true;
            for (int i = 0; i < lqn.nidx; i++) {
                if (lqn.iscaller.isAssigned(tidx, i) || lqn.iscaller.isAssigned(i, tidx)) {
                    isolated = false;
                    break;
                }
            }
            if (isolated) {
                continue;
            }
            boolean hasTaskCaller = false;
            for (int eidx : lqn.entriesof.get(tidx)) {
                for (int i = lqn.tshift; i < lqn.tshift + lqn.ntasks; i++) {
                    if (lqn.iscaller.get(i, eidx) != 0) {
                        hasTaskCaller = true;
                        break;
                    }
                }
                if (hasTaskCaller) {
                    break;
                }
            }
            if (hasTaskCaller) {
                servers.add(tidx);
            }
        }
        if (servers.isEmpty()) {
            throw new RuntimeException("Flat layering found no server: the model has no processor with tasks.");
        }
        return servers;
    }

    /** Build the single flat layer holding every processor and task. */
    private void buildFlatLayer(List<Integer> flatServers) {
        List<Integer> flatCallers = new ArrayList<>();
        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            if (this.ignore.get(tidx) == 0) {
                flatCallers.add(tidx);
            }
        }
        buildLayersRecursive(flatServers, flatCallers, false, true);

        thinkt_classes_updmap = integerMapToMatrix(cell_thinkt_classes_updmap);
        call_classes_updmap = integerMapToMatrix(cell_call_classes_updmap);
        servt_classes_updmap = integerMapToMatrix(cell_servt_classes_updmap);
        arvproc_classes_updmap = integerMapToMatrix(cell_arvproc_classes_updmap);
        route_prob_updmap = integerMapToMatrix(cell_route_prob_updmap);

        this.ensemble = new Network[nlayers];
        for (int i = 0; i < nlayers; i++) {
            ensemble[i] = temp_ensemble.get(i);
        }
        if (this.lqnModel != null) {
            this.lqnModel.setEnsemble(
                    new java.util.ArrayList<Network>(java.util.Arrays.asList(this.ensemble)));
        }

        // every server resolves to the single flat layer, which is at once the
        // host layer and the task layer
        for (int i = 0; i < lqn.nhosts + lqn.ntasks; i++) {
            this.idxhash.add(Double.NaN);
        }
        for (int sidx : flatServers) {
            this.idxhash.set(sidx, 1.0);
        }
        this.hostLayerIndices = new ArrayList<>();
        this.taskLayerIndices = new ArrayList<>();
        this.hostLayerIndices.add(0);
        this.taskLayerIndices.add(0);
    }

    /**
     * Station index (1-based) of ELEMIDX inside layer E, falling back to the
     * layer's own server when ELEMIDX is not a server there.
     */
    public int stationIdxOf(int e, int elemIdx) {
        NetworkAttribute attr = this.ensemble[e].getAttribute();
        Integer stn = attr.getServerIdxOf().get(elemIdx);
        return stn == null ? attr.getServerIdx() : stn;
    }

    /**
     * Station of layer E that class C (0-based) is served at: the processor of an
     * activity, the called task of a call, the layer's server otherwise.
     */
    public int stationIdxOfClass(int e, int c) {
        Integer[] attr = this.ensemble[e].getClasses().get(c).getAttribute();
        int elem = -1;
        if (attr != null && attr.length > 1 && attr[1] != null) {
            if (attr[0] == LayeredNetworkElement.ACTIVITY) {
                elem = (int) lqn.parent.get(0, (int) lqn.parent.get(0, attr[1]));
            } else if (attr[0] == LayeredNetworkElement.CALL) {
                elem = (int) lqn.parent.get(0, (int) lqn.callpair.get(attr[1], 1));
            }
        }
        return stationIdxOf(e, elem);
    }

    /** Station indices of the host (ISHOST true) or task servers of layer E. */
    public List<Integer> serverStationsOf(int e, boolean ishost) {
        NetworkAttribute attr = this.ensemble[e].getAttribute();
        return ishost ? attr.getHostStations() : attr.getTaskStations();
    }

    /** Stations of ELEMIDX when it is a server of this layer, null otherwise. */
    private static Map<Integer, Queue> serversFor(Map<Integer, Map<Integer, Queue>> srvStation, int elemIdx) {
        Map<Integer, Queue> st = srvStation.get(elemIdx);
        return (st == null || st.isEmpty()) ? null : st;
    }

    /** Declare JOBCLASS at every server station of this layer. */
    private static void setAllServers(Map<Integer, Map<Integer, Queue>> srvStation,
                                      JobClass jobclass, Distribution dist) {
        for (Map<Integer, Queue> stns : srvStation.values()) {
            for (Queue q : stns.values()) {
                q.setService(jobclass, dist);
            }
        }
    }

    /**
     * Initial service time of a call class at the stations of SEEDIDX, refined at
     * every iteration through call_classes_updmap.
     */
    /**
     * Resolves lqn.callgroups from target entries to call indices. ofCidx[cidx] is the
     * 1-based group index of that call, 0 when it is dispatched on its own; members.get(g-1)
     * lists the call indices of group g in declaration order and strategy.get(g-1) is its
     * RoutingStrategy. A group that does not resolve to at least two calls is dropped, so a
     * stale group cannot silently rewrite a single call's routing.
     */
    private void callGroupsByCidx(int[] ofCidx, List<List<Integer>> members,
                                  List<RoutingStrategy> strategy) {
        if (lqn.callgroups == null || lqn.callgroups.isEmpty()) {
            return;
        }
        for (LayeredNetworkStruct.CallGroupStruct grp : lqn.callgroups) {
            List<Integer> found = new ArrayList<>();
            List<Integer> callsOfCaller = lqn.callsof.get(grp.caller);
            if (callsOfCaller == null) {
                continue;
            }
            for (int cidx : callsOfCaller) {
                if (grp.targets.contains((int) lqn.callpair.get(cidx, 1))) {
                    found.add(cidx);
                }
            }
            if (found.size() >= 2) {
                members.add(found);
                strategy.add(grp.strategy);
                for (int cidx : found) {
                    ofCidx[cidx] = members.size();
                }
            }
        }
    }

    private void seedCallService(Map<Integer, Map<Integer, Queue>> srvStation, int seedIdx, JobClass cidxClass) {
        double minRespT = 0;
        if (lqn.actsof.get(seedIdx) != null) {
            for (int tidx_act : lqn.actsof.get(seedIdx)) {
                // upper bound, uses all activities not just the ones reachable by this entry
                minRespT = minRespT + lqn.hostdem_mean.getOrDefault(tidx_act, 0.0);
            }
        }
        Map<Integer, Queue> stns = serversFor(srvStation, seedIdx);
        if (stns == null) {
            return;
        }
        for (Queue q : stns.values()) {
            q.setService(cidxClass, minRespT == 0 ? Immediate.getInstance() : Exp.fitMean(minRespT));
        }
    }

    public void buildLayersRecursive(int idx, List<Integer> callers, boolean ishostlayer) {
        List<Integer> idxSet = new ArrayList<>();
        idxSet.add(idx);
        buildLayersRecursive(idxSet, callers, ishostlayer, false);
    }

    public void buildLayersRecursive(List<Integer> idxSet, List<Integer> callers, boolean ishostlayer, final boolean flat) {
        final int idx = idxSet.get(0); // layer key: model name, ensemble slot and update-map column
        nlayers++;
        Matrix jobPosKey = new Matrix(1, lqn.nidx);
        Map<Integer, JobClass> curClassKey = new HashMap<>(lqn.nidx);
        final Map<Integer, Map<Integer, Queue>> curStationKey = new HashMap<>(lqn.nidx);
        // Fan-out check: when all callers have fan-out >= nreplicas for this task,
        // each replica sees the full caller traffic (fork-join semantics).
        // Model a single representative replica; updateThinkTimes multiplies by K.
        // For host layers: if all caller tasks have the same replication as the host,
        // the host is co-replicated with the task, so also use single-replica modeling.
        final int nreplicas;
        {
            int rawReplicas = (int) lqn.repl.get(0, idx);
            boolean reduceFanout = false;
            if (!flat && rawReplicas > 1 && !callers.isEmpty()) {
                if (!ishostlayer && lqn.fanout != null) {
                    // Task layer: check if all callers fan-out to all replicas
                    reduceFanout = true;
                    for (int caller : callers) {
                        double fo = lqn.fanout.get(caller, idx);
                        if (fo < rawReplicas) {  // 0 means not set (default=1)
                            reduceFanout = false;
                            break;
                        }
                    }
                } else if (ishostlayer) {
                    // Host layer: if all caller tasks have the same replication as the host,
                    // model a single representative host replica (co-replicated with tasks)
                    reduceFanout = true;
                    for (int caller : callers) {
                        if ((int) lqn.repl.get(0, caller) != rawReplicas) {
                            reduceFanout = false;
                            break;
                        }
                    }
                }
            }
            // flat layering rejects replicated elements upstream
            nreplicas = (reduceFanout || flat) ? 1 : rawReplicas;
        }
        final boolean singleReplicaMode = (nreplicas == 1 && (int) lqn.repl.get(0, idx) > 1);
        // If this is a task layer with fan-out reduction, record the server task as single-replica
        if (singleReplicaMode && !ishostlayer) {
            this.singleReplicaTasks.add(idx);
        }
        Matrix mult = lqn.maxmult.copy(); // this removes spare capacity that cannot be used
        lqn.mult = mult.copy(); // using maxmult as in MATLAB version
        Network model = new Network(flat ? lqn.hashnames.get(idx) + ".Flat" : lqn.hashnames.get(idx));
        model.setChecks(false);

        // entries of every server element of this layer
        final List<Integer> entriesOfSet = new ArrayList<>();
        for (int sidx : idxSet) {
            if (lqn.entriesof.get(sidx) != null) {
                entriesOfSet.addAll(lqn.entriesof.get(sidx));
            }
        }

        boolean hasSynccaller = false;
        if (!ishostlayer) {
            for (int caller : callers) {
                for (int entries : entriesOfSet) {
                    if (lqn.issynccaller.isAssigned(caller, entries)) {
                        hasSynccaller = true;
                    }
                }
            }
        }
        boolean hasAsynccaller = false;
        if (!ishostlayer) {
            for (int caller : callers) {
                for (int entries : entriesOfSet) {
                    if (lqn.isasynccaller.isAssigned(caller, entries)) {
                        hasAsynccaller = true;
                    }
                }
            }
        }
        Delay clientDelay = null;
        if (flat || ishostlayer || hasSynccaller || hasAsynccaller) {
            clientDelay = new Delay(model, "Clients");
            model.getAttribute().setClientIdx(1);
            model.getAttribute().setServerIdx(2);
            model.getAttribute().setSourceIdx(-1);
        } else {
            model.getAttribute().setSourceIdx(-1);
            model.getAttribute().setServerIdx(1);
            model.getAttribute().setClientIdx(-1);
        }

        // One station (times its replicas) per server element of the layer
        final Map<Integer, Map<Integer, Queue>> srvStation = new HashMap<Integer, Map<Integer, Queue>>();
        for (int sidx : idxSet) {
            boolean sishost = sidx < lqn.nhosts;
            Map<Integer, Queue> stns = new HashMap<Integer, Queue>(nreplicas);
            for (int i = 1; i <= nreplicas; i++) {
                if (i == 1) {
                    stns.put(i, new Queue(model, lqn.hashnames.get(sidx), lqn.sched.get(sidx)));
                } else {
                    String name = lqn.hashnames.get(sidx).concat(".");
                    stns.put(i, new Queue(model, name.concat(Integer.toString(i)), lqn.sched.get(sidx)));
                }
                stns.get(i).setNumberOfServers((int) mult.get(0, sidx));
                stns.get(i).getAttribute().setIsHost(sishost);
                stns.get(i).getAttribute().setIdx(sidx);
                // LQN successive activities on the same host retain the server:
                // mark the layer queue with immediate feedback so that simulators
                // do not re-queue the job behind other waiting jobs on self-loops.
                stns.get(i).setImmediateFeedback(true);
            }
            srvStation.put(sidx, stns);
            model.getAttribute().putServerIdxOf(sidx, model.getNodeIndex(stns.get(1)) + 1);
        }
        // the layer's own server, sole server under 'srvn' layering
        Map<Integer, Queue> serverStation = srvStation.get(idx);

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
        
        // Detect function layer: all callers must be setup tasks and this must be a host layer
        boolean issetuplayer = false;
        if (ishostlayer) {
            boolean tempIsSetupLayer = true;
            for (int caller : callers) {
                if (lqn.hassetup != null && lqn.hassetup.get(0, caller) == 0) {
                    tempIsSetupLayer = false;
                    break;
                }
            }
            issetuplayer = tempIsSetupLayer;
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
        for (int i = lqn.ashift; i < lqn.nidx; i++) {
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
            for (int i = 0; i < lqn.graph.getNumCols(); i++) {
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

        // Routed call groups: one dispatch with n destinations, taken at a router
        // whose only links are those destinations. See _kb/06-solver-catalog.md
        int[] cgroupOfCidx = new int[lqn.ncalls + 1];
        List<List<Integer>> cgroupMembers = new ArrayList<>();
        List<RoutingStrategy> cgroupStrategy = new ArrayList<>();
        callGroupsByCidx(cgroupOfCidx, cgroupMembers, cgroupStrategy);
        Map<Integer, JobClass> grpDispatchClass = new HashMap<>();
        Map<Integer, JobClass> grpReturnClass = new HashMap<>();
        Map<Integer, Router> grpRouter = new HashMap<>();
        Set<Integer> groupRouted = new HashSet<>();
        List<Integer[]> routedGroupSites = new ArrayList<>(); // [router node, dispatch class, group]

        // Station indices of the layer's servers, kept apart from the hosts /
        // tasks rows whose entries are [class index, LQN element] pairs
        for (int sidx : idxSet) {
            int stn = model.getAttribute().getServerIdxOf().get(sidx);
            if (sidx < lqn.nhosts) {
                model.getAttribute().getHostStations().add(stn);
                if (!flat) {
                    model.getAttribute().addHosts(new Integer[]{null, stn});
                }
            } else {
                model.getAttribute().getTaskStations().add(stn);
                if (!flat) {
                    model.getAttribute().addTasks(new Integer[]{null, stn});
                }
            }
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
            if (!ishostlayer || flat) {
                for (int entry_idx : entriesOfSet) {
                    if (lqn.issynccaller.get(tidx_caller, entry_idx) != 0) {
                        isSyncCallerToEntries = true;
                        break;
                    }
                }
            }
            // Per-caller: for host layers, check if this caller task has direct callers
            // or any entry that is a forwarding target
            // (matching MATLAB buildLayersRecursive hasDirectCallers logic)
            boolean hasDirectCallers = false;
            boolean isForwardingTarget = false;
            boolean hasOpenArrival = false;
            boolean hostIsServer = serversFor(srvStation, (int) lqn.parent.get(0, tidx_caller)) != null;
            if (hostIsServer) {
                if (lqn.isref.get(tidx_caller) != 0) {
                    hasDirectCallers = true;
                } else {
                    for (int eidx : lqn.entriesof.get(tidx_caller)) {
                        // Check if any task is a sync or async caller to this entry
                        for (int row = 0; row < lqn.ntasks; row++) {
                            int tidx_row = row + lqn.tshift;
                            if (lqn.issynccaller.get(tidx_row, eidx) != 0 || lqn.isasynccaller.get(tidx_row, eidx) != 0) {
                                hasDirectCallers = true;
                                break;
                            }
                        }
                        // Check for open arrivals on this entry
                        if (lqn.arrival != null && lqn.arrival.containsKey(eidx) && lqn.arrival.get(eidx) != null) {
                            hasOpenArrival = true;
                        }
                        // Check if this entry is a forwarding target
                        for (int cidx_fwd = 0; cidx_fwd < lqn.ncalls; cidx_fwd++) {
                            if (lqn.calltype.get(cidx_fwd) == CallType.FWD && (int) lqn.callpair.get(cidx_fwd, 1) == eidx) {
                                isForwardingTarget = true;
                                break;
                            }
                        }
                    }
                }
            }
            // A task reached ONLY by an entry arrival has no task layer, because no task
            // calls it, so updateThinkTimes never gives its caller class a surrogate delay
            // and the class cycles against an Immediate one. Adding an open stream on top
            // of that unthrottled chain saturated lqn_open_arrival: the processor at 0.68
            // against 0.32 from lqns, lqsim and LDES alike. The chain is the better of the
            // two representations here, so the stream is dropped and the chain is closed
            // on the known arrival rate instead, exactly as a forwarding target is.
            boolean openArrivalOnly = hasOpenArrival && !hasDirectCallers && !isForwardingTarget
                    && !isSyncCallerToEntries && lqn.isref.get(tidx_caller) == 0;
            if ((hostIsServer && (hasDirectCallers || hasOpenArrival || isForwardingTarget)) || isSyncCallerToEntries) {
                if (this.njobs.get(tidx_caller, idx) == 0) {
                    // Use single-replica njobs if either: (1) this layer is in single-replica mode,
                    // or (2) the caller task itself is in single-replica mode
                    boolean callerIsSingleReplica = singleReplicaMode || this.singleReplicaTasks.contains(tidx_caller);
                    njobs = mult.get(0, tidx_caller) * (callerIsSingleReplica ? 1.0 : lqn.repl.get(0, tidx_caller));
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
                            for (int i = 0; i < mult.getNumCols(); i++) {
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
                setAllServers(srvStation, aidxclass.get(tidx_caller), Disabled.getInstance());
                aidxclass.get(tidx_caller).setReferenceClass(true);
                aidxclass.get(tidx_caller).setAttribute(new Integer[]{LayeredNetworkElement.TASK, tidx_caller});
                aidxclass.get(tidx_caller).setCompletes(false);
                model.getAttribute().addTasks(new Integer[]{aidxclass.get(tidx_caller).getIndex(), tidx_caller});
                assert clientDelay != null;
                if (lqn.isref.get(tidx_caller) != 0) {
                    clientDelay.setService(aidxclass.get(tidx_caller), thinkproc.get(tidx_caller));
                } else {
                    // a served task's declared think time is not a per-request
                    // delay, so the seed carries none either; updateThinkTimes
                    // replaces this from the first iteration on
                    clientDelay.setService(aidxclass.get(tidx_caller), Immediate.getInstance());
                }
                if (lqn.isref.get(tidx_caller) == 0) {
                    if (!cell_thinkt_classes_updmap.containsKey(idx)) {
                        cell_thinkt_classes_updmap.put(idx, new ArrayList<>());
                    }
                    cell_thinkt_classes_updmap.get(idx).add(new Integer[]{idx, tidx_caller, 1, aidxclass.get(tidx_caller).getIndex()});
                }

                for (int eidx : lqn.entriesof.get(tidx_caller)) {
                    aidxclass.put(eidx, new ClosedClass(model, lqn.hashnames.get(eidx), 0, clientDelay));
                    clientDelay.setService(aidxclass.get(eidx), Disabled.getInstance());
                    setAllServers(srvStation, aidxclass.get(eidx), Disabled.getInstance());
                    aidxclass.get(eidx).setCompletes(false);
                    aidxclass.get(eidx).setAttribute(new Integer[]{LayeredNetworkElement.ENTRY, eidx});
                    model.getAttribute().addEntries(new Integer[]{aidxclass.get(eidx).getIndex(), eidx});
                    clientDelay.setService(aidxclass.get(eidx), Immediate.getInstance());

                    // Check for open arrival distribution on this entry
                    if (!openArrivalOnly && lqn.arrival != null && lqn.arrival.containsKey(eidx) &&
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
                        // entries are Immediate, so the work is the bound activity's host
                        // demand (buildLayersRecursive.m:266-276); servtproc.get(eidx) is 0
                        // -1, not 0: the element space is 0-based, so 0 is the first
                        // host and cannot double as "no bound activity".
                        int boundAidx = -1;
                        for (int cand = 0; cand < lqn.nidx; cand++) {
                            if (lqn.graph.get(eidx, cand) > 0) {
                                boundAidx = cand;
                                break;
                            }
                        }
                        setAllServers(srvStation, openClassForEntry,
                                servtproc.get(boundAidx >= 0 ? boundAidx : eidx));

                        // Track for routing setup later: [class_index, entry_index]
                        entryOpenClassesAssignedLine[0]++;
                        entryOpenClasses.set(entryOpenClassesAssignedLine[0], 1, openClassForEntry.getIndex());
                        entryOpenClasses.set(entryOpenClassesAssignedLine[0], 2, eidx);

                        // Track: Use negative entry index to distinguish from call arrivals
                        if (!cell_arvproc_classes_updmap.containsKey(idx)) {
                            cell_arvproc_classes_updmap.put(idx, new ArrayList<>());
                        }
                        cell_arvproc_classes_updmap.get(idx).add(new Integer[]{idx, -eidx,
                            model.getNodeIndex(sourceStation) + 1, openClassForEntry.getIndex() + 1});

                        openClassForEntry.setCompletes(false);
                        openClassForEntry.setAttribute(new Integer[]{LayeredNetworkElement.ENTRY, eidx});
                    }
                }
            }
            // for each activity of the calling task
            for (int aidx : lqn.actsof.get(tidx_caller)) {
                if (hostIsServer || isSyncCallerToEntries) {
                    aidxclass.put(aidx, new ClosedClass(model, lqn.hashnames.get(aidx), 0, clientDelay));
                    clientDelay.setService(aidxclass.get(aidx), Disabled.getInstance());
                    setAllServers(srvStation, aidxclass.get(aidx), Disabled.getInstance());
                    aidxclass.get(aidx).setCompletes(false);
                    aidxclass.get(aidx).setAttribute(new Integer[]{LayeredNetworkElement.ACTIVITY, aidx});
                    model.getAttribute().addActivities(new Integer[]{aidxclass.get(aidx).getIndex(), aidx});
                    if (serversFor(srvStation, (int) lqn.parent.get(0, (int) lqn.parent.get(0, aidx))) == null) {
                        clientDelay.setService(aidxclass.get(aidx), servtproc.get(aidx));
                    }
                    if (iscachelayer) {
                        // Set cache read item entry for activities in cache layer
                        for (int eidx = 0; eidx < lqn.nentries; eidx++) {
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
                        int serverIdx = (int) lqn.callpair.get(cidx, 1);
                        if (serverIdx >= 0 && serversFor(srvStation, (int) lqn.parent.get(0, serverIdx)) != null) {
                            if (sourceStation == null) {
                                model.getAttribute().setSourceIdx(model.getNumberOfNodes() + 1);
                                sourceStation = new jline.lang.nodes.Source(model, "Source");
                                sinkStation = new Sink(model, "Sink");
                            }
                            cidxclass.put(cidx, new OpenClass(model, lqn.callhashnames.get(cidx), 0));
                            sourceStation.setArrival(cidxclass.get(cidx), Immediate.getInstance());
                            clientDelay.setService(cidxclass.get(cidx), Disabled.getInstance());
                            setAllServers(srvStation, cidxclass.get(cidx), Immediate.getInstance());
                            openClassesAssignedLine[0]++;
                            openClasses.set(openClassesAssignedLine[0], 1, cidxclass.get(cidx).getIndex());
                            openClasses.set(openClassesAssignedLine[0], 2, callmean.get(cidx));
                            openClasses.set(openClassesAssignedLine[0], 3, cidx);
                            model.getAttribute().addCalls(new Integer[]{cidxclass.get(cidx).getIndex(), cidx, (int) lqn.callpair.get(cidx, 0), (int) lqn.callpair.get(cidx, 1)});
                            cidxclass.get(cidx).setCompletes(false);
                            cidxclass.get(cidx).setAttribute(new Integer[]{LayeredNetworkElement.CALL, cidx});
                            seedCallService(srvStation, flat ? (int) lqn.parent.get(0, (int) lqn.callpair.get(cidx, 1)) : idx,
                                    cidxclass.get(cidx));
                        }
                    } else if (lqn.calltype.get(cidx) == CallType.SYNC) {
                        int gid = cgroupOfCidx[cidx];
                        if (gid > 0) {
                            // The members share the dispatch class, which is the one the
                            // strategy routes and the one that visits the targets. The hop
                            // must not switch class, because a state-dependent routing
                            // function is evaluated at zero off the class diagonal; the
                            // switch goes on the return arc into the group class.
                            if (!grpDispatchClass.containsKey(gid)) {
                                grpRouter.put(gid, new Router(model, lqn.hashnames.get(aidx) + ".Dispatch" + gid + ".Router"));
                                JobClass dispatchClass = new ClosedClass(model, lqn.hashnames.get(aidx) + ".Dispatch" + gid, 0, clientDelay);
                                clientDelay.setService(dispatchClass, Immediate.getInstance());
                                setAllServers(srvStation, dispatchClass, Disabled.getInstance());
                                dispatchClass.setCompletes(false);
                                dispatchClass.setAttribute(new Integer[]{LayeredNetworkElement.CALL, cidx});
                                JobClass returnClass = new ClosedClass(model, lqn.callhashnames.get(cidx) + ".Group" + gid, 0, clientDelay);
                                clientDelay.setService(returnClass, Immediate.getInstance());
                                setAllServers(srvStation, returnClass, Disabled.getInstance());
                                returnClass.setCompletes(false);
                                returnClass.setAttribute(new Integer[]{LayeredNetworkElement.CALL, cidx});
                                grpDispatchClass.put(gid, dispatchClass);
                                grpReturnClass.put(gid, returnClass);
                                routedGroupSites.add(new Integer[]{model.getNodeIndex(grpRouter.get(gid)), dispatchClass.getIndex(), gid});
                            }
                            cidxclass.put(cidx, grpDispatchClass.get(gid));
                        } else {
                            cidxclass.put(cidx, new ClosedClass(model, lqn.callhashnames.get(cidx), 0, clientDelay));
                            clientDelay.setService(cidxclass.get(cidx), Disabled.getInstance());
                            setAllServers(srvStation, cidxclass.get(cidx), Disabled.getInstance());
                            cidxclass.get(cidx).setCompletes(false);
                            cidxclass.get(cidx).setAttribute(new Integer[]{LayeredNetworkElement.CALL, cidx});
                        }
                        model.getAttribute().addCalls(new Integer[]{cidxclass.get(cidx).getIndex(), cidx, (int) lqn.callpair.get(cidx, 0), (int) lqn.callpair.get(cidx, 1)});
                        seedCallService(srvStation, flat ? (int) lqn.parent.get(0, (int) lqn.callpair.get(cidx, 1)) : idx,
                                cidxclass.get(cidx));
                    }

                    // an Aux class is needed whenever the call does not happen exactly once
                    // per activity execution, which is a property of callmean alone. A group
                    // member's mean is the 1/n share the dispatch already carries, so an Aux
                    // there would charge the skip path twice
                    if (callmean.get(cidx) != 1 && cgroupOfCidx[cidx] == 0) {
                        if (lqn.calltype.get(cidx) == CallType.SYNC) {
                            cidxauxclass.put(cidx, new ClosedClass(model, lqn.callhashnames.get(cidx) + ".Aux", 0, clientDelay));
                            cidxauxclass.get(cidx).setCompletes(false);
                            cidxauxclass.get(cidx).setAttribute(new Integer[]{LayeredNetworkElement.CALL, cidx});
                            assert clientDelay != null; // safety check not in MATLAB
                            clientDelay.setService(cidxauxclass.get(cidx), Immediate.getInstance());
                            setAllServers(srvStation, cidxauxclass.get(cidx), Disabled.getInstance());
                        }
                    }

                    // For SYNC calls in task layers, create classes for forwarding
                    // calls from the target entry. This implements synthetic
                    // synchronization: the caller blocks until the forwarding
                    // target completes, modeling contention correctly.
                    // FWD calls have the source ENTRY (not activity) in callpair(:,1),
                    // so we scan all calls to find FWD calls from the target entry.
                    // For chain forwarding (e0->e1->e2), recursively follow the chain.
                    // Forwarding chains are represented by caller-side pseudo
                    // rendezvous calls added by applyForwardingRendezvous, which
                    // are ordinary SYNC calls handled above; FWD calls need no
                    // classes in the layers.
                }
            }
        }

        // The fork-join transform mints its own Source/Sink pair, detaching the open
        // stream already routed through this one: see _kb/06-solver-catalog.md (LN section)
        if (sourceStation != null && hasFork) {
            throw new RuntimeException("SolverLN: layer '" + model.getName()
                    + "' carries both an AND fork and an open stream (an async call or an entry "
                    + "arrival); the fork-join transform needs a Source of its own");
        }

        RoutingMatrix P = model.initRoutingMatrix();
        if (sourceStation != null) {
            for (int o = 1; o <= openClassesAssignedLine[0]; o++) {
                int oidx = (int) openClasses.get(o, 1) - 1;
                double p = 1.0 / openClasses.get(o, 2);
                int cidx = (int) openClasses.get(o, 3); // 3 = source
                Map<Integer, Queue> tgtStation = serversFor(srvStation, (int) lqn.parent.get(0, (int) lqn.callpair.get(cidx, 1)));
                if (tgtStation == null) {
                    tgtStation = serverStation;
                }
                int ntgt = tgtStation.size();
                double callmean_o = openClasses.get(o, 2);
                if (callmean_o < 1) {
                    // fewer than one call per arrival: a single Bernoulli pass, the
                    // geometric loop below would need a negative repeat probability
                    P.addConnection(model.getClasses().get(oidx), model.getClasses().get(oidx), sourceStation, sinkStation, 1.0 - callmean_o);
                    for (int m = 1; m <= ntgt; m++) {
                        P.addConnection(model.getClasses().get(oidx), model.getClasses().get(oidx), sourceStation, tgtStation.get(m), callmean_o / (double) ntgt);
                        P.addConnection(model.getClasses().get(oidx), model.getClasses().get(oidx), tgtStation.get(m), sinkStation, 1.0);
                    }
                } else {
                    for (int m = 1; m <= ntgt; m++) {
                        P.addConnection(model.getClasses().get(oidx), model.getClasses().get(oidx), sourceStation, tgtStation.get(m), 1.0 / (double) ntgt);
                        for (int n = 1; n <= ntgt; n++) {
                            P.addConnection(model.getClasses().get(oidx), model.getClasses().get(oidx), tgtStation.get(m), tgtStation.get(n), (1.0 - p) / (double) ntgt);
                        }
                        P.addConnection(model.getClasses().get(oidx), model.getClasses().get(oidx), tgtStation.get(m), sinkStation, p);
                    }
                }
                if (!cell_arvproc_classes_updmap.containsKey(idx)) {
                    cell_arvproc_classes_updmap.put(idx, new ArrayList<>());
                }
                cell_arvproc_classes_updmap.get(idx).add(new Integer[]{idx, cidx, model.getNodeIndex(sourceStation) + 1, oidx + 1});
                for (int m = 1; m <= ntgt; m++) {
                    if (!cell_call_classes_updmap.containsKey(idx)) {
                        cell_call_classes_updmap.put(idx, new ArrayList<>());
                    }
                    cell_call_classes_updmap.get(idx).add(new Integer[]{idx, cidx, model.getNodeIndex(tgtStation.get(m)) + 1, oidx + 1});
                }
            }

        }

        int atClient = 1;
        int atServer = 2;
        int atCache = 3;

        // Create final copies for inner class access
        final boolean finalIsCacheLayer = iscachelayer;
        final boolean finalIsSetupLayer = issetuplayer;
        final Cache finalCacheNode = cacheNode;

        // Delayed-hit retrieval: a dedicated fetch station in the cache sublayer so the
        // closed AMVA (da_cacheqn_retrieval) captures the finite-population coalescing.
        // EXPERIMENTAL - FURTHER WORK NEEDED. This wiring makes LN(MVA)/LN(NC) produce
        // the correct cache hit/miss PROBABILITIES for LCQ models with retrieval and
        // captures the coalescing throughput benefit in DIRECTION only (understated
        // magnitude vs LDES; the delayed-hit fraction is not recovered - see
        // Da_cacheqn_retrieval LIMITATIONS). Java port of the retrieval block in
        // matlab/src/solvers/LN/@SolverLN/buildLayersRecursive.m.
        boolean hasRetrievalCache = false;
        if (iscachelayer && lqn.hasretrieval != null) {
            for (int c : callers) if (lqn.hasretrieval.get(0, c) != 0) { hasRetrievalCache = true; break; }
        }
        Queue retrievalStation = null;
        if (hasRetrievalCache) {
            retrievalStation = new Queue(model, lqn.hashnames.get(callers.get(0)) + ".Fetch", SchedStrategy.PS);
        }
        final boolean finalHasRetrievalCache = hasRetrievalCache;
        final JobClass[] retrievalReadClass = new JobClass[1];
        final JobClass[] retrievalMissClass = new JobClass[1];
        final int[] retrievalMissAidx = new int[]{-1};

        class InnerRecurActGraph {

            recurActGraphReturnType recurActGraph(RoutingMatrix P, int tidx_caller, int aidx, JobClass curClass, int jobPos, Source sourceStation, Delay clientDelay, Sink sinkStation, Join joinNode, Fork forkNode, Stack<JobClass> forkClassStack) {
                return recurActGraph(P, tidx_caller, aidx, curClass, jobPos, sourceStation, clientDelay, sinkStation, joinNode, forkNode, forkClassStack, new HashMap<Integer, Queue>());
            }

            recurActGraphReturnType recurActGraph(RoutingMatrix P, int tidx_caller, int aidx, JobClass curClass, int jobPos, Source sourceStation, Delay clientDelay, Sink sinkStation, Join joinNode, Fork forkNode, Stack<JobClass> forkClassStack, Map<Integer, Queue> curStations) {
                jobPosKey.set(0, aidx, jobPos);
                curClassKey.put(aidx, curClass);
                curStationKey.put(aidx, curStations);
                List<Integer> nextaidxs = new ArrayList<>();
                for (int i = 0; i < lqn.graph.getNumCols(); i++) {
                    if (lqn.graph.isAssigned(aidx, i)) {
                        nextaidxs.add(i);
                    }
                }


                Matrix isNextPrecFork = new Matrix(1, lqn.nidx, lqn.nidx);
                if (!nextaidxs.isEmpty()) {
                    isNextPrecFork.set(0, aidx, 0);
                    for (int i : nextaidxs) {
                        if (isPostAndAct.isAssigned(0, i)) {
                            isNextPrecFork.set(0, aidx, 1);
                            break;
                        }
                    }
                }
                // Pre-fork state, captured at the first branch so that calls this
                // activity issues before the fork stay sequential
                boolean forkSaved = false;
                JobClass forkSaveCurClass = null;
                int forkSaveJobPos = 0;
                Map<Integer, Queue> forkSaveStations = null;
                if (!nextaidxs.isEmpty()) {
                    for (int nextaidx : nextaidxs) {
                        // Restore pre-fork state for each branch iteration
                        if (isNextPrecFork.get(0, aidx) != 0) {
                            if (!forkSaved) {
                                if (isPostAndAct.isAssigned(0, nextaidx)) {
                                    forkSaved = true;
                                    forkSaveCurClass = curClass;
                                    forkSaveJobPos = jobPos;
                                    forkSaveStations = curStations;
                                }
                            } else {
                                curClass = forkSaveCurClass;
                                jobPos = forkSaveJobPos;
                                curStations = forkSaveStations;
                            }
                        }
                        boolean isLoop = lqn.graph.get(aidx, nextaidx) != lqn.dag.get(aidx, nextaidx);
                        // in the activity graph, the following if is entered only
                        // by an edge that is the return from a LOOP activity
                        if (lqn.parent.get(0, aidx) != lqn.parent.get(0, nextaidx)) { // if different parent task
                            int cidx = 0;
                            for (int i = 0; i < lqn.ncalls; i++) { // find the call index
                                if (lqn.callpair.get(i, 0) == aidx && lqn.callpair.get(i, 1) == nextaidx) {
                                    cidx = i;
                                    break;
                                }
                            }
                            if (lqn.calltype.get(cidx) == CallType.ASYNC) {
                                // an async call does not block the caller, so it adds no
                                // routing here; its open class is declared once in the
                                // class-declaration pass above (MATLAB: empty ASYNC case)
                            } else if (lqn.calltype.get(cidx) == CallType.SYNC) {
                                // START routeSynchCall in MATLAB
                                int gidRoute = cgroupOfCidx[cidx];
                                if (gidRoute > 0) {
                                    if (groupRouted.contains(gidRoute)) {
                                        continue; // the group is one dispatch, already wired
                                    }
                                    List<Queue> grpTgtStn = new ArrayList<>();
                                    List<Integer> grpTgtCidx = new ArrayList<>();
                                    for (int mcidx : cgroupMembers.get(gidRoute - 1)) {
                                        int meidx = (int) lqn.callpair.get(mcidx, 1);
                                        Map<Integer, Queue> mstn = meidx >= 0
                                                ? serversFor(srvStation, (int) lqn.parent.get(meidx)) : null;
                                        if (mstn != null && !mstn.isEmpty()) {
                                            grpTgtStn.add(mstn.get(1));
                                            grpTgtCidx.add(mcidx);
                                        }
                                    }
                                    if (grpTgtStn.size() >= 2) {
                                        // A routed group is ONE hop with n destinations, taken
                                        // at a router whose only links are those destinations:
                                        // the strategy routes over a NODE's links, not over one
                                        // class's arcs. The 1/n split laid down here is the
                                        // probabilistic reading a solver without state-dependent
                                        // routing would see; it is replaced after link().
                                        Node fromNode = jobPos == atClient ? clientDelay
                                                : curStations.get(curStations.keySet().iterator().next());
                                        JobClass dispCls = grpDispatchClass.get(gidRoute);
                                        JobClass retCls = grpReturnClass.get(gidRoute);
                                        double share = 1.0 / grpTgtStn.size();
                                        P.addConnection(curClass, dispCls, fromNode, grpRouter.get(gidRoute), 1.0);
                                        for (int m = 0; m < grpTgtStn.size(); m++) {
                                            P.addConnection(dispCls, dispCls, grpRouter.get(gidRoute), grpTgtStn.get(m), share);
                                            P.addConnection(dispCls, retCls, grpTgtStn.get(m), clientDelay, 1.0);
                                            grpTgtStn.get(m).setService(dispCls, callservtproc.get(grpTgtCidx.get(m)));
                                            if (!cell_call_classes_updmap.containsKey(idx)) {
                                                cell_call_classes_updmap.put(idx, new ArrayList<>());
                                            }
                                            cell_call_classes_updmap.get(idx).add(new Integer[]{idx,
                                                    grpTgtCidx.get(m), model.getNodeIndex(grpTgtStn.get(m)) + 1,
                                                    dispCls.getIndex()});
                                        }
                                        groupRouted.add(gidRoute);
                                        curClass = retCls;
                                        jobPos = atClient;
                                        curStations = new HashMap<Integer, Queue>();
                                        continue;
                                    }
                                }
                                if (jobPos == atClient) {
                                    int serverIdx = (int) lqn.callpair.get(cidx, 1);
                                    // stations of the called task, null if it is not a server of this layer
                                    Map<Integer, Queue> tstn = serverIdx >= 0
                                            ? serversFor(srvStation, (int) lqn.parent.get(serverIdx)) : null;
                                    int ntgt = tstn == null ? 0 : tstn.size();
                                    if (tstn != null) {
                                        if (callmean.get(cidx) < 1) {
                                            P.addConnection(curClass, cidxauxclass.get(cidx), clientDelay, clientDelay, 1.0 - callmean.get(cidx));
                                            for (int m = 1; m <= ntgt; m++) {
                                                P.addConnection(curClass, cidxclass.get(cidx), clientDelay, tstn.get(m), callmean.get(cidx) / (double) ntgt);
                                                P.addConnection(cidxclass.get(cidx), cidxclass.get(cidx), tstn.get(m), clientDelay, 1.0);
                                            }
                                            P.addConnection(cidxauxclass.get(cidx), cidxclass.get(cidx), clientDelay, clientDelay, 1.0);
                                        } else if (callmean.get(cidx) == 1) {
                                            for (int m = 1; m <= ntgt; m++) {
                                                P.addConnection(curClass, cidxclass.get(cidx), clientDelay, tstn.get(m), 1.0 / (double) ntgt);
                                                P.addConnection(cidxclass.get(cidx), cidxclass.get(cidx), tstn.get(m), clientDelay, 1.0);
                                            }
                                        } else { // callmean.get(cidx) > 1
                                            for (int m = 1; m <= ntgt; m++) {
                                                P.addConnection(curClass, cidxclass.get(cidx), clientDelay, tstn.get(m), 1.0 / (double) ntgt);
                                                P.addConnection(cidxclass.get(cidx), cidxauxclass.get(cidx), tstn.get(m), clientDelay, 1.0);
                                                P.addConnection(cidxauxclass.get(cidx), cidxclass.get(cidx), clientDelay, tstn.get(m), (1.0 - 1.0 / callmean.get(cidx)) / ntgt);
                                            }
                                            P.addConnection(cidxauxclass.get(cidx), cidxclass.get(cidx), clientDelay, clientDelay, 1.0 / callmean.get(cidx));
                                        }
                                        jobPos = atClient;
                                        curStations = new HashMap<Integer, Queue>();
                                        clientDelay.setService(cidxclass.get(cidx), Immediate.getInstance());
                                        if (!cell_call_classes_updmap.containsKey(idx)) {
                                            cell_call_classes_updmap.put(idx, new ArrayList<>());
                                        }
                                        for (int m = 1; m <= ntgt; m++) {
                                            tstn.get(m).setService(cidxclass.get(cidx), callservtproc.get(cidx));
                                            cell_call_classes_updmap.get(idx).add(new Integer[]{idx, cidx, model.getNodeIndex(tstn.get(m)) + 1, cidxclass.get(cidx).getIndex()});
                                        }
                                        curClass = cidxclass.get(cidx);
                                    } else { // if it is not a call to an entry of a server in this layer
                                        if (callmean.get(cidx) < 1) {
                                            // The mean number of calls is embedded in the demand
                                            // (callservt = callmean * W), so the class must be
                                            // visited deterministically; a Bernoulli(callmean)
                                            // visit would discount the call time twice.
                                            P.addConnection(curClass, cidxclass.get(cidx), clientDelay, clientDelay, 1.0);
                                            P.addConnection(cidxclass.get(cidx), cidxauxclass.get(cidx), clientDelay, clientDelay, 1.0);
                                            curClass = cidxauxclass.get(cidx);
                                        } else if (callmean.get(cidx) == 1) {
                                            P.addConnection(curClass, cidxclass.get(cidx), clientDelay, clientDelay, 1.0);
                                            curClass = cidxclass.get(cidx);
                                        } else {  // callmean.get(cidx) > 1
                                            P.addConnection(curClass, cidxclass.get(cidx), clientDelay, clientDelay, 1.0); // the mean number of calls is now embedded in the demand
                                            P.addConnection(cidxclass.get(cidx), cidxauxclass.get(cidx), clientDelay, clientDelay, 1.0);  // the mean number of calls is now embedded in the demand
                                            curClass = cidxauxclass.get(cidx);
                                        }
                                        jobPos = atClient;
                                        curStations = new HashMap<Integer, Queue>();
                                        clientDelay.setService(cidxclass.get(cidx), callservtproc.get(cidx));
                                        if (!cell_call_classes_updmap.containsKey(idx)) {
                                            cell_call_classes_updmap.put(idx, new ArrayList<>());
                                        }
                                        cell_call_classes_updmap.get(idx).add(new Integer[]{idx, cidx, 1, cidxclass.get(cidx).getIndex()});
                                    }
                                } else if (jobPos == atServer) {
                                    int serverIdx = (int) lqn.callpair.get(cidx, 1);
                                    Map<Integer, Queue> tstn = serverIdx >= 0
                                            ? serversFor(srvStation, (int) lqn.parent.get(serverIdx)) : null;
                                    int ntgt = tstn == null ? 0 : tstn.size();
                                    if (tstn != null) {
                                        if (callmean.get(cidx) < 1) {
                                            // The call is skipped with probability 1-callmean; the two
                                            // flows merge back in the call class at the client, as in
                                            // the atClient branch. Routing the skip into the call class
                                            // instead would leave the Aux class with no inbound arc and
                                            // its chain without a reference class.
                                            for (int m = 1; m <= ntgt; m++) {
                                                Queue fromStn = curStations.get(Math.min(m, curStations.size()));
                                                P.addConnection(curClass, cidxauxclass.get(cidx), fromStn, clientDelay, 1.0 - callmean.get(cidx));
                                                P.addConnection(curClass, cidxclass.get(cidx), fromStn, tstn.get(m), callmean.get(cidx) / (double) ntgt);
                                                P.addConnection(cidxclass.get(cidx), cidxclass.get(cidx), tstn.get(m), clientDelay, 1.0);
                                                tstn.get(m).setService(cidxclass.get(cidx), callservtproc.get(cidx));
                                            }
                                            P.addConnection(cidxauxclass.get(cidx), cidxclass.get(cidx), clientDelay, clientDelay, 1.0);
                                            // The reply transits the client in the call class, and
                                            // sn_refresh_visits drops any (station, class) state whose
                                            // rate is NaN, so the class must be declared there
                                            clientDelay.setService(cidxclass.get(cidx), Immediate.getInstance());
                                            jobPos = atClient;
                                            curStations = new HashMap<Integer, Queue>();
                                            curClass = cidxclass.get(cidx);
                                        } else if (callmean.get(cidx) == 1) {
                                            for (int m = 1; m <= ntgt; m++) {
                                                Queue fromStn = curStations.get(Math.min(m, curStations.size()));
                                                P.addConnection(curClass, cidxclass.get(cidx), fromStn, tstn.get(m), 1.0);
                                            }
                                            if (isFlatLayering()) {
                                                // the reply returns the job to the client, which is
                                                // where the successor restoration expects it
                                                for (int m = 1; m <= ntgt; m++) {
                                                    P.addConnection(cidxclass.get(cidx), cidxclass.get(cidx), tstn.get(m), clientDelay, 1.0);
                                                }
                                                clientDelay.setService(cidxclass.get(cidx), Immediate.getInstance());
                                                jobPos = atClient;
                                                curStations = new HashMap<Integer, Queue>();
                                            } else {
                                                jobPos = atServer;
                                                curStations = tstn;
                                            }
                                            curClass = cidxclass.get(cidx);
                                        } else {
                                            for (int m = 1; m <= ntgt; m++) {
                                                Queue fromStn = curStations.get(Math.min(m, curStations.size()));
                                                P.addConnection(curClass, cidxclass.get(cidx), fromStn, tstn.get(m), 1.0);
                                            }
                                            if (isFlatLayering()) {
                                                // the geometric repeat transits the client between
                                                // visits; a self-loop would merge them into one
                                                for (int m = 1; m <= ntgt; m++) {
                                                    P.addConnection(cidxclass.get(cidx), cidxauxclass.get(cidx), tstn.get(m), clientDelay, 1.0);
                                                    P.addConnection(cidxauxclass.get(cidx), cidxclass.get(cidx), clientDelay, tstn.get(m), (1.0 - 1.0 / callmean.get(cidx)) / (double) ntgt);
                                                }
                                                P.addConnection(cidxauxclass.get(cidx), cidxclass.get(cidx), clientDelay, clientDelay, 1.0 / callmean.get(cidx));
                                                clientDelay.setService(cidxclass.get(cidx), Immediate.getInstance());
                                                curClass = cidxclass.get(cidx);
                                            } else {
                                                for (int m = 1; m <= ntgt; m++) {
                                                    P.addConnection(cidxclass.get(cidx), cidxclass.get(cidx), tstn.get(m), tstn.get(m), 1.0 - 1.0 / callmean.get(cidx));
                                                    P.addConnection(cidxclass.get(cidx), cidxauxclass.get(cidx), tstn.get(m), clientDelay, 1.0 / callmean.get(cidx));
                                                }
                                                curClass = cidxauxclass.get(cidx);
                                            }
                                            jobPos = atClient;
                                            curStations = new HashMap<Integer, Queue>();
                                        }
                                        if (!cell_call_classes_updmap.containsKey(idx)) {
                                            cell_call_classes_updmap.put(idx, new ArrayList<>());
                                        }
                                        for (int m = 1; m <= ntgt; m++) {
                                            tstn.get(m).setService(cidxclass.get(cidx), callservtproc.get(cidx));
                                            cell_call_classes_updmap.get(idx).add(new Integer[]{idx, cidx, model.getNodeIndex(tstn.get(m)) + 1, cidxclass.get(cidx).getIndex()}); // check if getIndex+1
                                        }
                                    } else {
                                        // if it is not a call to an entry of a server in this layer
                                        // callmean not needed since we switched
                                        // to ResidT to model service time at client
                                        if (callmean.get(cidx) < 1) {
                                            for (Queue fromStn : curStations.values()) {
                                                P.addConnection(curClass, cidxclass.get(cidx), fromStn, clientDelay, 1.0);
                                            }
                                            P.addConnection(cidxclass.get(cidx), cidxauxclass.get(cidx), clientDelay, clientDelay, 1.0);
                                            curClass = cidxauxclass.get(cidx);
                                        } else if (callmean.get(cidx) == 1) {
                                            for (Queue fromStn : curStations.values()) {
                                                P.addConnection(curClass, cidxclass.get(cidx), fromStn, clientDelay, 1.0);
                                            }
                                            curClass = cidxclass.get(cidx);
                                        } else { // callmean.get(cidx) > 1
                                            for (Queue fromStn : curStations.values()) {
                                                P.addConnection(curClass, cidxclass.get(cidx), fromStn, clientDelay, 1.0);
                                            }
                                            P.addConnection(cidxclass.get(cidx), cidxauxclass.get(cidx), clientDelay, clientDelay, 1.0);
                                            curClass = cidxauxclass.get(cidx);
                                        }
                                        jobPos = atClient;
                                        curStations = new HashMap<Integer, Queue>();
                                        clientDelay.setService(cidxclass.get(cidx), callservtproc.get(cidx));
                                        if (!cell_call_classes_updmap.containsKey(idx)) {
                                            cell_call_classes_updmap.put(idx, new ArrayList<>());
                                        }
                                        cell_call_classes_updmap.get(idx).add(new Integer[]{idx, cidx, 1, cidxclass.get(cidx).getIndex()}); // check if getIndex+1
                                    }
                                    // END routeSynchCall
                                }
                                // (forwarding handled via pseudo rendezvous calls)
                            }
                        } else {
                            boolean intersects;
                            intersects = false;
                            Set<Integer> range = new HashSet<>();
                            for (int i = 0; i < lqn.nentries; i++) {
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
                                curStations = curStationKey.get(aidx);
                                if (curStations == null) {
                                    curStations = new HashMap<Integer, Queue>();
                                }
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
                                curStations = new HashMap<Integer, Queue>();
                                curClass = curClassC;
                            }

                            // stations of the processor the next activity runs on, null if
                            // that processor is not a server of this layer
                            Map<Integer, Queue> hstn = serversFor(srvStation, (int) lqn.parent.get(0, (int) lqn.parent.get(0, nextaidx)));
                            int nhstn = hstn == null ? 0 : hstn.size();
                            if (jobPos == atClient) {
                                if (hstn != null) {
                                    if (!finalIsCacheLayer) {
                                        for (int m = 1; m <= nhstn; m++) {
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
                                                P.addConnection(curClass, aidxclass.get(nextaidx), forkOutputRouter.get(f), hstn.get(m), 1.0);
                                            } else {
                                                // If nextaidx is not a post-and activity, use default routing without fork
                                                P.addConnection(curClass, aidxclass.get(nextaidx), clientDelay, hstn.get(m), lqn.graph.get(aidx, nextaidx));
                                            }
                                        } else {
                                            if (isPreAndAct.get(0, aidx) != 0) {
                                                JobClass forkClass = (JobClass) forkClassStack.pop();
                                                applyJoinQuorum(joinNode, forkClass, nextaidx);
                                                P.addConnection(curClass, forkClass, clientDelay, joinNode, 1.0);
                                                P.addConnection(forkClass, aidxclass.get(nextaidx), joinNode, hstn.get(m), 1.0);
                                            } else {
                                                P.addConnection(curClass, aidxclass.get(nextaidx), clientDelay, hstn.get(m), lqn.graph.get(aidx, nextaidx));
                                            }
                                        }
                                        hstn.get(m).setService(aidxclass.get(nextaidx), lqn.hostdem.get(nextaidx));
                                        // A SetupTask's cold start is NOT wired into the layer
                                        // station any more: it is charged to the entry with
                                        // probability p, see setupCharge.
                                        }
                                        jobPos = atServer;
                                        curStations = hstn;
                                        curClass = aidxclass.get(nextaidx);
                                        if (!cell_servt_classes_updmap.containsKey(idx)) {
                                            cell_servt_classes_updmap.put(idx, new ArrayList<>());
                                        }
                                        cell_servt_classes_updmap.get(idx).add(new Integer[]{idx, nextaidx, model.getNodeIndex(hstn.get(1)) + 1, aidxclass.get(nextaidx).getIndex()});
                                    } else {
                                        // Cache layer routing: client -> cache -> server
                                        P.addConnection(curClass, aidxclass.get(nextaidx), clientDelay, finalCacheNode, lqn.graph.get(aidx, nextaidx));

                                        // Setup cache read item entry
                                        finalCacheNode.setReadItemEntry(aidxclass.get(nextaidx), lqn.itemproc.get(aidx), (int) lqn.nitems.get(0, aidx));
                                        
                                        // Find hit and miss activities
                                        List<Integer> hitmissaidx = new ArrayList<>();
                                        for (int i = 0; i < lqn.nidx; i++) {
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

                                            if (finalHasRetrievalCache) {
                                                // Record the fetch-station wiring; setRetrievalSystem
                                                // mints per-item retrieval classes, so it is deferred
                                                // to after the activity-graph traversal (before link)
                                                // to avoid perturbing the class indexing mid-build.
                                                retrievalReadClass[0] = aidxclass.get(nextaidx);
                                                retrievalMissClass[0] = aidxclass.get(missaidx);
                                                retrievalMissAidx[0] = missaidx;
                                            }
                                        }

                                        jobPos = atCache;
                                        curStations = new HashMap<Integer, Queue>();
                                        curClass = aidxclass.get(nextaidx);
                                    }
                                } else { // the processor of the next activity is not a server of this layer
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
                                            applyJoinQuorum(joinNode, forkClass, nextaidx);
                                            P.addConnection(curClass, forkClass, clientDelay, joinNode, 1.0);
                                            P.addConnection(forkClass, aidxclass.get(nextaidx), joinNode, clientDelay, 1.0);
                                        } else {
                                            P.addConnection(curClass, aidxclass.get(nextaidx), clientDelay, clientDelay, lqn.graph.get(aidx, nextaidx));
                                        }
                                    }
                                    jobPos = atClient;
                                    curStations = new HashMap<Integer, Queue>();
                                    curClass = aidxclass.get(nextaidx);
                                    clientDelay.setService(aidxclass.get(nextaidx), servtproc.get(nextaidx));
                                    if (!cell_thinkt_classes_updmap.containsKey(idx)) {
                                        cell_thinkt_classes_updmap.put(idx, new ArrayList<>());
                                    }
                                    cell_thinkt_classes_updmap.get(idx).add(new Integer[]{idx, nextaidx, 1, aidxclass.get(nextaidx).getIndex()});
                                }
                            } else if (jobPos == atServer || jobPos == atCache) {
                                if (hstn != null) {
                                    if (jobPos == atCache) {
                                        // Cache layer routing: cache -> server
                                        curClass = aidxclass.get(nextaidx);
                                        for (int m = 1; m <= nhstn; m++) {
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
                                                    P.addConnection(curClass, aidxclass.get(nextaidx), forkOutputRouter.get(f), hstn.get(m), 1.0);
                                                } else {
                                                    // If nextaidx is not a post-and activity or routing is unavailable, use default routing
                                                    P.addConnection(curClass, aidxclass.get(nextaidx), finalCacheNode, hstn.get(m), lqn.graph.get(aidx, nextaidx));
                                                }
                                            } else {
                                                if (isPreAndAct.get(0, aidx) != 0) {
                                                    JobClass forkClass = (JobClass) forkClassStack.pop();
                                                    applyJoinQuorum(joinNode, forkClass, nextaidx);
                                                    P.addConnection(curClass, forkClass, finalCacheNode, joinNode, 1.0);
                                                    P.addConnection(forkClass, aidxclass.get(nextaidx), joinNode, hstn.get(m), 1.0);
                                                } else {
                                                    P.addConnection(curClass, aidxclass.get(nextaidx), finalCacheNode, hstn.get(m), lqn.graph.get(aidx, nextaidx));
                                                }
                                            }
                                            Distribution baseService = lqn.hostdem.get(nextaidx);
                                            // A SetupTask's cold start is NOT wired into the layer
                                            // station any more, and it is not folded into the host
                                            // demand either: it is charged to the entry with
                                            // probability p, see setupCharge.
                                            hstn.get(m).setService(aidxclass.get(nextaidx), baseService);
                                        }
                                    } else {
                                        for (int m = 1; m <= nhstn; m++) {
                                        Queue fromStn = curStations.get(Math.min(m, curStations.size()));
                                        if (isNextPrecFork.get(0, aidx) != 0) {
                                            P.addConnection(curClass, curClass, fromStn, forkNode, 1.0);
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
                                                P.addConnection(curClass, aidxclass.get(nextaidx), forkOutputRouter.get(f), hstn.get(m), 1.0);
                                            } else {
                                                // If nextaidx is not a post-and activity, use default routing without fork
                                                P.addConnection(curClass, aidxclass.get(nextaidx), fromStn, hstn.get(m), lqn.graph.get(aidx, nextaidx));
                                            }
                                        } else {
                                            if (isPreAndAct.get(0, aidx) != 0) {
                                                JobClass forkClass = (JobClass) forkClassStack.pop();
                                                applyJoinQuorum(joinNode, forkClass, nextaidx);
                                                P.addConnection(curClass, forkClass, fromStn, joinNode, 1.0);
                                                P.addConnection(forkClass, aidxclass.get(nextaidx), joinNode, hstn.get(m), 1.0);
                                            } else {
                                                P.addConnection(curClass, aidxclass.get(nextaidx), fromStn, hstn.get(m), lqn.graph.get(aidx, nextaidx));
                                            }
                                        }
                                        hstn.get(m).setService(aidxclass.get(nextaidx), lqn.hostdem.get(nextaidx));
                                        // A SetupTask's cold start is NOT wired into the layer
                                        // station any more: it is charged to the entry with
                                        // probability p, see setupCharge.
                                        }
                                    }

                                    jobPos = atServer;
                                    curStations = hstn;
                                    curClass = aidxclass.get(nextaidx);
                                    if (!cell_servt_classes_updmap.containsKey(idx)) {
                                        cell_servt_classes_updmap.put(idx, new ArrayList<>());
                                    }
                                    cell_servt_classes_updmap.get(idx).add(new Integer[]{idx, nextaidx, model.getNodeIndex(hstn.get(1)) + 1, aidxclass.get(nextaidx).getIndex()});
                                } else {
                                    for (int m = 1; m <= curStations.size(); m++) {
                                        Queue fromStn = curStations.get(m);
                                        if (isNextPrecFork.get(0, aidx) != 0) {
                                            P.addConnection(curClass, curClass, fromStn, forkNode, 1.0);
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
                                                P.addConnection(curClass, aidxclass.get(nextaidx), fromStn, clientDelay, lqn.graph.get(aidx, nextaidx));
                                            }
                                        } else {
                                            if (isPreAndAct.get(0, aidx) != 0) {
                                                JobClass forkClass = (JobClass) forkClassStack.pop();
                                                applyJoinQuorum(joinNode, forkClass, nextaidx);
                                                P.addConnection(curClass, forkClass, fromStn, joinNode, 1.0);
                                                P.addConnection(forkClass, aidxclass.get(nextaidx), joinNode, clientDelay, 1.0);
                                            } else {
                                                P.addConnection(curClass, aidxclass.get(nextaidx), fromStn, clientDelay, lqn.graph.get(aidx, nextaidx));
                                            }
                                        }
                                        jobPos = atClient;
                                        curStations = new HashMap<Integer, Queue>();
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
                                // Save curClassC before recursion - it's a class-level field
                                // but in MATLAB it's local to each recurActGraph call scope
                                JobClass savedCurClassC = curClassC;
                                recurActGraphReturnType returnType = recurActGraph(P, tidx_caller, nextaidx, curClass, jobPos, sourceStation, clientDelay, sinkStation, joinNode, forkNode, forkClassStack, curStations);
                                curClassC = savedCurClassC;
                                P = returnType.P;
                                curClass = returnType.curClass;
                                jobPos = returnType.jobPos;
                                curStations = returnType.curStations != null ? returnType.curStations : new HashMap<Integer, Queue>();

                                if (jobPos == atClient) {
                                    P.addConnection(curClass, aidxclass.get(tidx_caller), clientDelay, clientDelay, 1.0);
                                    if (!curClass.getName().endsWith(".Aux")) {
                                        curClass.setCompletes(true);
                                    }
                                } else {
                                    for (Queue fromStn : curStations.values()) {
                                        P.addConnection(curClass, aidxclass.get(tidx_caller), fromStn, clientDelay, 1.0);
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

                        recurActGraphReturnType(curClass, jobPos, P, curStations);
            }
        }

        int jobPos = atClient; // start at client
        // second pass: setup the routing out of entries
        for (int tidx_caller : callers) {
            // Per-caller: check if this caller is a sync caller to entries of idx
            // For host layers, idx is a host index with no entries, so this is always false
            boolean isSyncCallerToTargetEntries = false;
            if (!ishostlayer || flat) {
                for (int entry_idx : entriesOfSet) {
                    if (lqn.issynccaller.get(tidx_caller, entry_idx) != 0) {
                        isSyncCallerToTargetEntries = true;
                        break;
                    }
                }
            }
            // Per-caller: recompute hasDirectCallers for host layers (same as first pass)
            boolean hasDirectCallers2 = false;
            boolean isForwardingTarget2 = false;
            boolean hostIsServer2 = serversFor(srvStation, (int) lqn.parent.get(0, tidx_caller)) != null;
            if (hostIsServer2) {
                if (lqn.isref.get(tidx_caller) != 0) {
                    hasDirectCallers2 = true;
                } else {
                    for (int eidx : lqn.entriesof.get(tidx_caller)) {
                        for (int row = 0; row < lqn.ntasks; row++) {
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
                        // Check if this entry is a forwarding target
                        for (int cidx_fwd = 0; cidx_fwd < lqn.ncalls; cidx_fwd++) {
                            if (lqn.calltype.get(cidx_fwd) == CallType.FWD && (int) lqn.callpair.get(cidx_fwd, 1) == eidx) {
                                isForwardingTarget2 = true;
                                break;
                            }
                        }
                    }
                }
            }
            if ((hostIsServer2 && (hasDirectCallers2 || isForwardingTarget2)) || isSyncCallerToTargetEntries) { // if it is only an asynch caller the closed classes are not needed
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

        // Setup routing for entry-level open arrivals (AFTER recurActGraph to avoid being overwritten)
        if (sourceStation != null && entryOpenClassesAssignedLine[0] > 0) {
            for (int e = 1; e <= entryOpenClassesAssignedLine[0]; e++) {
                int eoidx = (int) entryOpenClasses.get(e, 1) - 1; // class index (0-based)
                JobClass openClass = model.getClasses().get(eoidx);

                // Explicitly set routing: ONLY Source → Server → Sink
                // Zero out all routing for this class first
                for (jline.lang.nodes.Node node1 : model.getNodes()) {
                    for (jline.lang.nodes.Node node2 : model.getNodes()) {
                        P.addConnection(openClass, openClass, node1, node2, 0.0);
                    }
                }

                // Station an open arrival at this entry enters, the processor of its
                // task under host layering and the task itself under flat layering
                Map<Integer, Queue> eoStation = flat
                        ? serversFor(srvStation, (int) lqn.parent.get(0, (int) entryOpenClasses.get(e, 2)))
                        : serverStation;
                if (eoStation == null) {
                    eoStation = serverStation;
                }
                int neo = eoStation.size();
                for (int m = 1; m <= neo; m++) {
                    // Route: Source -> station of the entry -> Sink
                    P.addConnection(openClass, openClass, sourceStation, eoStation.get(m), 1.0 / (double) neo);
                    P.addConnection(openClass, openClass, eoStation.get(m), sinkStation, 1.0);
                }
            }
        }

        // Apply the deferred delayed-hit retrieval wiring (mirrors the retrieval block
        // in matlab buildLayersRecursive.m). setRetrievalSystem mints per-item retrieval
        // classes that violate the LN sublayer-class == LQN-activity invariant, so we
        // (set fetch service, grow P, pad hitClass/services, tag the retrieval classes
        // auxiliary) to keep the LN bookkeeping consistent.
        if (hasRetrievalCache && retrievalReadClass[0] != null) {
            // Fetch service = the miss activity's full service (host demand + backend call).
            retrievalStation.setService(retrievalReadClass[0], this.servtproc.get(retrievalMissAidx[0]));
            P.addConnection(retrievalReadClass[0], retrievalReadClass[0], finalCacheNode, retrievalStation, 1.0);
            P.addConnection(retrievalReadClass[0], retrievalReadClass[0], retrievalStation, finalCacheNode, 1.0);
            finalCacheNode.setRetrievalSystem(retrievalReadClass[0], retrievalMissClass[0], new Queue[]{retrievalStation});

            // Tag the auto-generated retrieval classes as auxiliary: non-LQN attribute
            // and completes=false so the LN result aggregation (getEnsembleAvg) skips
            // them (they map to no LQN activity). JAR retrievalClasses are 0-indexed.
            Matrix rc = finalCacheNode.getRetrievalClasses();
            java.util.Set<Integer> rcvals = new java.util.HashSet<Integer>();
            for (int a = 0; a < rc.getNumRows(); a++) {
                for (int b = 0; b < rc.getNumCols(); b++) {
                    int v = (int) rc.get(a, b);
                    if (v >= 0) rcvals.add(v);
                }
            }
            for (Integer rci : rcvals) {
                JobClass jc = model.getClasses().get(rci);
                jc.setAttribute(new Integer[]{-1, -1});
                jc.setCompletes(false);
            }

            // Pad hitClass/missClass to equal length (setRetrievalSystem extends missClass
            // for the new retrieval classes but leaves hitClass shorter).
            Matrix hc = finalCacheNode.getHitClass();
            Matrix mc = finalCacheNode.getMissClass();
            int Lhm = Math.max(hc.length(), mc.length());
            if (hc.length() < Lhm) hc.expandMatrix(1, Lhm, Lhm);
            if (mc.length() < Lhm) mc.expandMatrix(1, Lhm, Lhm);

            // Pad every service station's service to the full class count with Disabled:
            // the retrieval classes are served only at the fetch station.
            int KfullSvc = model.getNumberOfClasses();
            for (Station stp : model.getStations()) {
                if (stp instanceof Queue) {
                    Queue ss = (Queue) stp;
                    for (int rp = 0; rp < KfullSvc; rp++) {
                        JobClass jc = model.getClasses().get(rp);
                        if (ss.getService(jc) == null) ss.setService(jc, Disabled.getInstance());
                    }
                }
            }

            // Grow P to the full class count (empty routing for the retrieval classes -
            // their cache<->fetch circulation is applied by getStruct from the read
            // class's template) so link does not index past P.
            for (JobClass jc : model.getClasses()) P.addClass(jc);
        }

        if (flat) {
            // link() installs RAND routing for every (node, class) pair left without
            // an outgoing arc. With one station per server in a single layer those
            // spurious uniform arcs let a class wander to stations it never visits,
            // trapping the flow in a sub-cycle and leaving the reference class with
            // zero visits. Disable them, as LQN2QN does for its signal classes.
            int nnodes = model.getNumberOfNodes();
            int nclasses = model.getNumberOfClasses();
            // one pass per class pair, row sums taken in bulk: fetching the block
            // once per node instead made this step O(nnodes^2 nclasses^2) and
            // dominated construction on layers with a thousand classes
            for (int rcls = 0; rcls < nclasses; rcls++) {
                double[] outflow = new double[nnodes];
                for (int scls = 0; scls < nclasses; scls++) {
                    Matrix block = P.get(rcls + 1, scls + 1);
                    if (block == null) {
                        continue;
                    }
                    Matrix rowSums = block.sumRows();
                    int rows = FastMath.min(nnodes, rowSums.getNumRows());
                    for (int inode = 0; inode < rows; inode++) {
                        outflow[inode] += rowSums.get(inode, 0);
                    }
                }
                for (int inode = 0; inode < nnodes; inode++) {
                    if (model.getNodes().get(inode) instanceof Sink) {
                        continue;
                    }
                    if (outflow[inode] < GlobalConstants.FineTol) {
                        model.getNodes().get(inode).setRouting(model.getClasses().get(rcls), RoutingStrategy.DISABLED);
                    }
                }
            }
        }

        model.link(P);
        // link() installs the probabilistic split; the declared strategy replaces it on
        // the router (node, class) whose only links are the group's targets
        for (Integer[] site : routedGroupSites) {
            // getIndex() is 1-based on a JobClass, the node list is 0-based
            model.getNodes().get(site[0]).setRouting(model.getClasses().get(site[1] - 1),
                    cgroupStrategy.get(site[2] - 1));
        }

        // Admission constraint on the server station -- see _kb/06-solver-catalog.md (LN section)
        if (lqn.lincon != null && lqn.lincon.containsKey(idx) && lqn.lincon.get(idx)[0] != null) {
            Matrix Aelem = lqn.lincon.get(idx)[0];
            Matrix Alayer = new Matrix(Aelem.getNumRows(), model.getClasses().size());
            // column j of Aelem is task constrainedIdx(j) on a host layer, entry constrainedIdx(j) otherwise
            List<Integer> constrainedIdx = ishostlayer ? lqn.tasksof.get(idx) : lqn.entriesof.get(idx);
            if (constrainedIdx == null) {
                constrainedIdx = new java.util.ArrayList<Integer>();
            }
            for (int j = 0; j < constrainedIdx.size(); j++) {
                List<JobClass> layerClasses = new java.util.ArrayList<JobClass>();
                if (ishostlayer) {
                    // a task occupies the host through the classes of its activities
                    List<Integer> acts = lqn.actsof.get(constrainedIdx.get(j));
                    if (acts != null) {
                        for (int k = 0; k < acts.size(); k++) {
                            JobClass jc = aidxclass.get(acts.get(k));
                            if (jc != null) {
                                layerClasses.add(jc);
                            }
                        }
                    }
                } else {
                    // an entry is occupied by the classes of the calls that target it
                    for (int c = 0; c < lqn.ncalls; c++) {
                        if (lqn.callpair.get(c, 1) == constrainedIdx.get(j)) {
                            JobClass jc = cidxclass.get(c);
                            if (jc != null) {
                                layerClasses.add(jc);
                            }
                        }
                    }
                }
                for (int k = 0; k < layerClasses.size(); k++) {
                    int col = layerClasses.get(k).getIndex() - 1;
                    for (int r = 0; r < Aelem.getNumRows(); r++) {
                        Alayer.set(r, col, Alayer.get(r, col) + Aelem.get(r, j));
                    }
                }
            }
            if (Alayer.elementSum() > 0) {
                // One region spanning every replica: the constraint models a passive
                // resource of the server as a whole (a semaphore, a connection pool),
                // so replicas share the tokens rather than each holding a private copy
                List<Node> regionNodes = new java.util.ArrayList<Node>();
                for (int m = 1; m <= nreplicas; m++) {
                    regionNodes.add(serverStation.get(m));
                }
                Region fcr = model.addRegion(regionNodes);
                fcr.setLinearConstraints(Alayer, lqn.lincon.get(idx)[1]);
            }
        }

        // Service-rate dependence on a server station -- see _kb/06-solver-catalog.md (LN section)
        if (lqn.lldscaling != null || lqn.cdscaling != null || lqn.jdscaling != null
                || lqn.pools != null) {
            int R = model.getClasses().size();
            for (int sidx : idxSet) {
                boolean hasld = lqn.lldscaling != null && lqn.lldscaling.containsKey(sidx);
                boolean hascd = lqn.cdscaling != null && lqn.cdscaling.containsKey(sidx);
                boolean hasjd = lqn.jdscaling != null && lqn.jdscaling.containsKey(sidx);
                boolean haspools = lqn.pools != null && lqn.pools.containsKey(sidx);
                // A compatibility declaration is a rate law in its own right, and
                // addServerType refuses to coexist with a joint dependence, so a
                // pooled server carries NO lld/cd/jd handle: testing only those
                // three skipped every pools-only station and dropped the
                // compatibility silently, leaving the layer as the plain
                // multiserver the declaration exists to say it is not.
                if (!hasld && !hascd && !hasjd && !haspools) {
                    continue;
                }
                boolean sishost = sidx < lqn.nhosts;
                List<Integer> operandIdx = sishost ? lqn.tasksof.get(sidx) : lqn.entriesof.get(sidx);
                if (operandIdx == null) {
                    operandIdx = new java.util.ArrayList<Integer>();
                }
                List<List<Integer>> cols = new java.util.ArrayList<List<Integer>>();
                for (int j = 0; j < operandIdx.size(); j++) {
                    List<Integer> colsOfOperand = new java.util.ArrayList<Integer>();
                    if (sishost) {
                        // a task occupies the host through the classes of its activities
                        List<Integer> acts = lqn.actsof.get(operandIdx.get(j));
                        if (acts != null) {
                            for (int k = 0; k < acts.size(); k++) {
                                JobClass jc = aidxclass.get(acts.get(k));
                                if (jc != null) {
                                    colsOfOperand.add(jc.getIndex());
                                }
                            }
                        }
                    } else {
                        // an entry is occupied by the classes of the calls that target it
                        for (int c = 0; c < lqn.ncalls; c++) {
                            if (lqn.callpair.get(c, 1) == operandIdx.get(j)) {
                                JobClass jc = cidxclass.get(c);
                                if (jc != null) {
                                    colsOfOperand.add(jc.getIndex());
                                }
                            }
                        }
                    }
                    cols.add(colsOfOperand);
                }
                boolean oneClassPerOperand = true;
                for (int j = 0; j < cols.size(); j++) {
                    if (cols.get(j).size() > 1) {
                        oneClassPerOperand = false;
                    }
                }
                Map<Integer, Queue> stns = srvStation.get(sidx);
                for (int m = 1; m <= nreplicas; m++) {
                    if (hasld) {
                        stns.get(m).setLimitedLoadDependence(lqn.lldscaling.get(sidx));
                    }
                    if (hascd) {
                        // beta_{i,r} is product-form only while an operand maps to a single
                        // class; where it aggregates several, the same scaling is emitted as
                        // a joint dependence, which is numerically identical but not exact
                        SerializableFunction<Matrix, Matrix> cdh = layerDepHandle(lqn.cdscaling.get(sidx), cols, R, model);
                        Matrix cdpeak = layerPeak(lqn.cdscalingpeak.get(sidx), cols, R);
                        if (oneClassPerOperand) {
                            stns.get(m).setLimitedClassDependence(cdh, cdpeak);
                        } else {
                            stns.get(m).setLimitedJointDependence(cdh, cdpeak);
                        }
                    }
                    if (hasjd) {
                        stns.get(m).setLimitedJointDependence(layerDepHandle(lqn.jdscaling.get(sidx), cols, R, model),
                                layerPeak(lqn.jdscalingpeak.get(sidx), cols, R));
                    }
                    if (haspools) {
                        // A compatibility declaration IS a rate law: the pools clear mu(n)
                        // of SnCompatRate, order independent at every integer state.
                        // snCompatScaling normalises it against the rate the SAME
                        // population would get under full compatibility, so eta isolates
                        // the compatibility GRAPH and a fully-compatible pool is the
                        // neutral eta == 1; the low-occupancy loss stays with the solver's
                        // own multiserver term.
                        // The lowering is to a JOINT dependence, hence an approximation in
                        // the layer: see _kb/06-solver-catalog.md (LN section) for why the
                        // exact OI analyzer cannot serve a class-switching layer.
                        final LayeredNetworkStruct.ServerPools pl = lqn.pools.get(sidx);
                        final double[] pcounts = new double[pl.counts.getNumCols()];
                        final double[] prates = new double[pl.rates.getNumCols()];
                        for (int t = 0; t < pcounts.length; t++) {
                            pcounts[t] = pl.counts.get(0, t);
                            prates[t] = pl.rates.get(0, t);
                        }
                        SerializableFunction<Matrix, Matrix> etaPool =
                                new SerializableFunction<Matrix, Matrix>() {
                                    @Override
                                    public Matrix apply(Matrix nop) {
                                        double[] n = new double[nop.getNumCols()];
                                        for (int j = 0; j < n.length; j++) {
                                            n[j] = nop.get(0, j);
                                        }
                                        Matrix out = new Matrix(1, 1);
                                        out.set(0, 0, SnCompatRate.snCompatScaling(
                                                pl.compat, pcounts, prates, n));
                                        return out;
                                    }
                                };
                        // A rate-scaled station reports U = T*S/peak, so the peak
                        // is the rate the pools clear with EVERY server active --
                        // sum_t counts*rates, which for unit-rate pools is just the
                        // server count. Leaving it at 1 drops the division by the
                        // multiplicity and reports a fully-compatible pool at S
                        // times the utilization of the plain multiserver it is
                        // supposed to reproduce.
                        double poolPeak = SnCompatRate.snCompatPeak(pcounts, prates);
                        Matrix poolPeakRow = new Matrix(1, cols.size());
                        for (int j = 0; j < cols.size(); j++) {
                            poolPeakRow.set(0, j, poolPeak);
                        }
                        stns.get(m).setLimitedJointDependence(layerDepHandle(etaPool, cols, R, model),
                                layerPeak(poolPeakRow, cols, R));
                    }
                }
            }
        }

        temp_ensemble.add(model);
    }

    /**
     * <p>Lifts a service-rate dependence handle declared on a LayeredNetwork server
     * to the layer station that represents it. F maps the per-operand population
     * vector of that server to a 1x1 scaling shared by every operand or to a
     * per-operand row vector; COLS.get(j) lists the layer classes (1-based)
     * through which operand j occupies the station.</p>
     *
     * <p>Solvers evaluate the handle in two different index spaces: CTMC and the
     * exact recursions pass a per-class vector, while the AMVA and NC chain
     * recursions pass a per-chain vector. The handle therefore reads the length
     * of its argument to pick the space, aggregates the operand populations in
     * it, and answers a vector of the SAME length, since the caller indexes the
     * answer with the index it passed in. An index belonging to no operand keeps
     * the neutral scaling 1.</p>
     */
    private static SerializableFunction<Matrix, Matrix> layerDepHandle(final SerializableFunction<Matrix, Matrix> f,
                                                                      final List<List<Integer>> cols, final int R,
                                                                      final Network model) {
        return new SerializableFunction<Matrix, Matrix>() {
            private final Map<Integer, List<List<Integer>>> chainCols = new HashMap<Integer, List<List<Integer>>>();

            @Override
            public Matrix apply(Matrix n) {
                int L = n.length();
                List<List<Integer>> idx = cols;
                if (L != R) {
                    if (!chainCols.containsKey(L)) {
                        chainCols.put(L, layerChainCols(cols, model, L));
                    }
                    idx = chainCols.get(L);
                }
                int K = idx.size();
                Matrix nop = new Matrix(1, K);
                for (int j = 0; j < K; j++) {
                    double s = 0;
                    for (int k = 0; k < idx.get(j).size(); k++) {
                        s += n.get(idx.get(j).get(k) - 1);
                    }
                    nop.set(0, j, s);
                }
                Matrix w = f.apply(nop);
                Matrix v = new Matrix(1, L);
                for (int r = 0; r < L; r++) {
                    v.set(0, r, 1.0);
                }
                for (int j = 0; j < K; j++) {
                    double wj = w.get(Math.min(j, w.length() - 1));
                    for (int k = 0; k < idx.get(j).size(); k++) {
                        v.set(0, idx.get(j).get(k) - 1, wj);
                    }
                }
                return v;
            }
        };
    }

    /** Operand columns of a layer station in the chain index space. */
    private static List<List<Integer>> layerChainCols(List<List<Integer>> cols, Network model, int nchains) {
        Matrix chains = model.getStruct().chains;
        List<List<Integer>> out = new ArrayList<List<Integer>>();
        for (int j = 0; j < cols.size(); j++) {
            List<Integer> ch = new ArrayList<Integer>();
            for (int k = 0; k < cols.get(j).size(); k++) {
                int c = cols.get(j).get(k) - 1;
                for (int r = 0; r < chains.getNumRows(); r++) {
                    if (chains.get(r, c) > 0 && r + 1 <= nchains && !ch.contains(r + 1)) {
                        ch.add(r + 1);
                    }
                }
            }
            java.util.Collections.sort(ch);
            out.add(ch);
        }
        return out;
    }

    /** Spreads a per-operand peak rate scaling onto the classes of the layer station. */
    private static Matrix layerPeak(Matrix peakPerOperand, List<List<Integer>> cols, int R) {
        Matrix peak = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            peak.set(0, r, 1.0);
        }
        for (int j = 0; j < cols.size(); j++) {
            double pj = peakPerOperand.get(Math.min(j, peakPerOperand.length() - 1));
            for (int k = 0; k < cols.get(j).size(); k++) {
                peak.set(0, cols.get(j).get(k) - 1, pj);
            }
        }
        return peak;
    }

    /**
     * Rejects activity graphs whose AND forks and joins are not properly nested.
     * The traversal in buildLayersRecursive pairs a join with the most recent
     * fork through a LIFO stack of fork classes, so it can only represent
     * series-parallel graphs. A join whose inputs come from different forks pops
     * a class that was never pushed; failing here names the model instead.
     */
    public void assertSeriesParallelForks() {
        int ashift = lqn.nhosts + lqn.ntasks + lqn.nentries;
        // actpretype marks the join INPUTS; the join itself is their successor
        for (int j = ashift; j < lqn.nidx; j++) {
            java.util.Set<Integer> forks = new java.util.HashSet<Integer>();
            for (int i = ashift; i < lqn.nidx; i++) {
                if (lqn.graph.get(i, j) == 0
                        || lqn.actpretype.get(0, i) != ActivityPrecedenceType.ID_PRE_AND) {
                    continue;
                }
                forks.add(enclosingFork(i, ashift));
            }
            if (forks.size() > 1 || forks.contains(-1)) {
                line_error(mfilename(new Object() {
                }), "Activity '" + lqn.hashnames.get(j) + "' joins branches of different AND forks; "
                        + "SolverLN supports only properly nested (series-parallel) fork-join graphs.");
            }
        }
    }

    /** Fork whose branch subtree contains activity A, -1 when A is outside every fork. */
    private int enclosingFork(int a, int ashift) {
        java.util.Set<Integer> seen = new java.util.HashSet<Integer>();
        java.util.ArrayDeque<Integer> queue = new java.util.ArrayDeque<Integer>();
        queue.add(a);
        while (!queue.isEmpty()) {
            int cur = queue.poll();
            if (!seen.add(cur)) {
                continue;
            }
            for (int p = ashift; p < lqn.nidx; p++) {
                if (lqn.graph.get(p, cur) == 0) {
                    continue;
                }
                if (isAndFork(p, ashift)) {
                    return p;
                }
                queue.add(p);
            }
        }
        return -1;
    }

    /** True when activity F is the pre activity of an AND fork precedence. */
    private boolean isAndFork(int f, int ashift) {
        for (int b = ashift; b < lqn.nidx; b++) {
            if (lqn.graph.get(f, b) != 0 && lqn.actposttype.get(0, b) == ActivityPrecedenceType.ID_POST_AND) {
                return true;
            }
        }
        return false;
    }

    public void construct() {
        // mark down to ignore unreachable disconnected components
        this.ignore = new Matrix(lqn.nidx, 1);
        Set<Set<Integer>> wccs = weaklyConnect(lqn.graph, null);
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
        this.singleReplicaTasks = new java.util.HashSet<>();
        this.entrycdfrespt = new HashMap<Integer, Matrix>();
        this.hasconverged = false;
        this.momentPassDone = false;

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
        for (int cidx = 0; cidx < lqn.ncalls; cidx++) {
            int serverIdx = (int) lqn.callpair.get(cidx, 1);
            if (serverIdx >= 0 && lqn.hostdem.containsKey(serverIdx)) {
                callservtproc.put(cidx, lqn.hostdem.get(serverIdx));
            }
        }

        assertSeriesParallelForks();

        // perform layering
        this.njobs = new Matrix(lqn.tshift + lqn.ntasks, lqn.tshift + lqn.ntasks, lqn.nidx * lqn.nidx);
        this.idxhash = new ArrayList<>();
        buildLayers();
        this.solvers = new NetworkSolver[nlayers + 1];
        this.njobsorig = new Matrix(this.njobs);

        // Build the interlock path tables of Sec. 4.2
        if (this.options.config.interlocking) {
            initInterlock();
        }

        // initialize data structures for interlock correction
        this.ptaskcallers = new Matrix(lqn.nhosts + lqn.ntasks, lqn.nhosts + lqn.ntasks, lqn.nidx * lqn.nidx);
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
        // convergence test dispatch: see _kb/06-solver-catalog.md ("Convergence test: stochastic iteration dispatch")

        boolean bool = false;

        // Stochastic iteration dispatch: see _kb/06-solver-catalog.md (LN
        // section) for the rationale.
        if (this.stochiterMode != null) {
            if (this.stochiterAuto && "off".equals(this.stochiterMode) && it >= 1 && anyStochasticLayer()) {
                // a layer with method 'default' resolved at runtime to a
                // stochastic method (captured in analyze() at iteration 1)
                this.stochiterMode = "rm";
                line_debug(options.verbose, "LN: stochastic layer method detected at runtime, switching to Robbins-Monro iteration");
            }
            if ("rm".equals(this.stochiterMode)) {
                return convergedStoch(it);
            }
        }

        // The moment3 pass is terminal: it runs once hasconverged is set, and its
        // own output perturbs the error test below, which would clear hasconverged
        // and hand the next iteration back to the residence-time branch. See
        // BUGS.md BUG-97 and the note where momentPassDone is set.
        if ("moment3".equals(this.lnmethod) && this.momentPassDone) {
            return true;
        }

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
            // assume ready state.
            // The window averages results[size] ... results[size-wnd_size+1] and
            // `results` is keyed 1..it, so a window longer than the history reaches
            // key <= 0 and dereferences null. `it >= iter_min` alone is a sufficient
            // guard only while iter_min >= wnd_size, i.e. iter_min >= 5, and
            // iter_min = max(2*nlayers, ceil(iter_max/4)) drops below 5 on a SMALL
            // ensemble with a SMALL iter_max -- 2 layers and iter_max <= 19 give
            // iter_min = 4 -- which is what a bounded-iteration probe sets. The
            // second clause is the C++ twin's (solver_ln.h:2411): with too little
            // history, do not average at all rather than average over a shortened
            // window, so no run that already had enough history changes. MATLAB's
            // twin (@SolverLN/converged.m:95) is still unclamped and errors here.
            if (it >= iter_min && results.size() >= wnd_size) {
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
                jline.io.LineConsole.iter(it,
                        "layer iteration %d: max queue-length change %.3e (tolerance %.3e)",
                        it, this.maxitererr.get(it), this.options.iter_tol);
//                if (this.options.verbose != VerboseLevel.SILENT) {
//                    if (this.solvers[e].options.verbose != VerboseLevel.SILENT) {
//                        String msg = String.format("\bQLen change: %.5f.\n", this.maxitererr.get(it) / E);
//                        System.out.print(msg);
//                    }
//                }
                if (it == iter_min) {
                    if (jline.io.LineConsole.isActive()) {
                        jline.io.LineConsole.step("started averaging the iterates to aid convergence");
                    } else if (this.options.verbose != VerboseLevel.SILENT) {
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
        if (it > 1 && this.maxitererr != null) {
            line_debug(options.verbose, String.format("LN convergence check: it=%d, maxerr=%e, tol=%e",
                it, this.maxitererr.get(it), this.options.iter_tol));
        }
        if (it == 0 && (this.options.verbose != VerboseLevel.SILENT)) {
            // Debug output removed
        } else if ((it > iter_min) && (this.maxitererr.get(it) < this.options.iter_tol) && (this.maxitererr.get(it - 1) < this.options.iter_tol) && (this.maxitererr.get(it - 2) < this.options.iter_tol)) {
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

    /**
     * Convergence controller for stochastic layer solvers (Robbins-Monro mode).
     *
     * <p>When one or more layer solvers return noisy estimates (simulation,
     * e.g. JMT/SSA/LDES, or Monte Carlo integration, e.g. NC with mci/imci/ls),
     * the deterministic Picard iteration in converged() cannot terminate: the
     * successive-difference error is bounded below by the standard error of
     * the layer estimates, and the layer-reset confirmation step merely
     * resamples the noise. This routine implements a stochastic approximation
     * iteration instead:</p>
     * <ol>
     * <li>Burn-in: for the first stochiter_burnin iterations the plain Picard
     * iteration runs with the relaxation factor configured at init.</li>
     * <li>Robbins-Monro step: afterwards the relaxation factor applied by
     * updateMetrics to the fed-forward iterate (servt, residt, tput,
     * callservt) decays as omega_k = a0/k^alpha with alpha in (0.5,1]. Under
     * the contraction assumption already made by the deterministic iteration,
     * and zero-mean noise with bounded variance, the iterate converges almost
     * surely to the true fixed point (Robbins and Monro, 1951). Layer seeds
     * are rotated per iteration in pre() so successive evaluations observe
     * independent noise.</li>
     * <li>Polyak-Ruppert averaging: running averages of the layer results and
     * of the reported iterates are maintained and installed as the final
     * solution in finish(), giving the optimal O(1/sqrt(k)) rate and
     * robustness to the choice of a0 (Polyak and Juditsky, 1992).</li>
     * <li>Stopping: iteration stops when the drift of the averaged results
     * stays below iter_tol for stochiter_conseq consecutive iterations. The
     * drift of a running average decays like 1/k even under persistent noise,
     * so the test terminates, and it self-calibrates: larger noise keeps the
     * drift above tolerance longer, forcing more averaging.</li>
     * </ol>
     *
     * <p>Note: the Robbins-Monro step acts through relax_omega, which is
     * applied by the default metric update path; the moment3 update path does
     * not use relaxation, so this controller is primarily intended for method
     * 'default'.</p>
     *
     * @param it the completed iteration count
     * @return true when the averaged iterate has converged
     */
    public boolean convergedStoch(int it) {
        if (it < 1) {
            return false;
        }
        int E = this.nlayers;
        int burnin = this.options.config.stochiter_burnin;
        double a0 = this.options.config.stochiter_a0;
        double alpha = this.options.config.stochiter_alpha;

        if (this.maxitererr == null) {
            this.maxitererr = new ArrayList<Double>();
        }
        while (this.maxitererr.size() <= it) {
            this.maxitererr.add(0.0);
        }

        // Schedule the Robbins-Monro step used by updateMetrics at the next iteration
        if (it >= burnin) {
            this.relax_omega = FastMath.min(1.0, a0 / FastMath.pow(FastMath.max(1, it - burnin + 1), alpha));
        }

        if (it <= burnin) {
            // pure Picard burn-in; no averaging or convergence testing yet
            this.maxitererr.set(it, Double.POSITIVE_INFINITY);
            if (this.options.verbose != VerboseLevel.SILENT) {
                System.out.printf("Stochastic iteration burn-in %d/%d.", it, burnin);
            }
            return false;
        }

        if (this.stochiterStart == null) {
            this.stochiterStart = it;
            if (this.options.verbose != VerboseLevel.SILENT) {
                System.out.print(" Started Robbins-Monro averaging (stochastic layer solvers detected).");
            }
        }

        // Polyak-Ruppert update of the layer result averages and drift metric
        int k = this.stochAvgCount + 1;
        double err = 0.0;
        for (int e = 0; e < E; e++) {
            SolverResult raw = this.results.get(this.results.size()).get(e);
            SolverResult avg = new SolverResult();
            if (k == 1) {
                avg.QN = raw.QN == null ? null : raw.QN.copy();
                avg.UN = raw.UN == null ? null : raw.UN.copy();
                avg.RN = raw.RN == null ? null : raw.RN.copy();
                avg.TN = raw.TN == null ? null : raw.TN.copy();
                avg.AN = raw.AN == null ? null : raw.AN.copy();
                avg.WN = raw.WN == null ? null : raw.WN.copy();
            } else {
                SolverResult prev = this.stochAvg.get(e);
                avg.QN = polyakAvg(prev.QN, raw.QN, k);
                avg.UN = polyakAvg(prev.UN, raw.UN, k);
                avg.RN = polyakAvg(prev.RN, raw.RN, k);
                avg.TN = polyakAvg(prev.TN, raw.TN, k);
                avg.AN = polyakAvg(prev.AN, raw.AN, k);
                avg.WN = polyakAvg(prev.WN, raw.WN, k);
                // drift of the averaged queue lengths, normalized by population
                double N = this.ensemble[e].getNumberOfJobs().elementSum();
                if (N > 0 && avg.QN != null && prev.QN != null) {
                    try {
                        Matrix drift = avg.QN.sub(1, prev.QN);
                        drift.removeNaN();
                        drift.absEq();
                        err += drift.elementMax() / N;
                    } catch (Exception exception) {
                        // dimension mismatch across iterations; skip layer
                    }
                }
            }
            this.stochAvg.put(e, avg);
        }
        this.stochAvgCount = k;

        // Polyak-Ruppert averages of the fed-forward iterates used in reporting
        if (k == 1) {
            this.stochServtAvg = this.servt == null ? null : this.servt.copy();
            this.stochResidtAvg = this.residt == null ? null : this.residt.copy();
        } else {
            this.stochServtAvg = polyakAvg(this.stochServtAvg, this.servt, k);
            this.stochResidtAvg = polyakAvg(this.stochResidtAvg, this.residt, k);
        }

        this.maxitererr.set(it, err);
        if (this.options.verbose != VerboseLevel.SILENT) {
            System.out.printf("RMIterErr=%.6e (tol=%.6e, omega=%.3f, k=%d)", err, this.options.iter_tol, this.relax_omega, k);
        }

        // Stop when the averaged-iterate drift stays below tolerance
        int conseq = this.options.config.stochiter_conseq;
        if (k > conseq) {
            boolean below = true;
            for (int j = it - conseq + 1; j <= it; j++) {
                if (!(this.maxitererr.get(j) < this.options.iter_tol)) {
                    below = false;
                    break;
                }
            }
            this.hasconverged = below;
            return below;
        }
        return false;
    }

    /**
     * Running-mean update m = prev + (raw - prev)/k, robust to NaN entries in
     * either operand (a NaN sample leaves the average untouched).
     */
    private static Matrix polyakAvg(Matrix prev, Matrix raw, int k) {
        if (prev == null) {
            return raw == null ? null : raw.copy();
        }
        if (raw == null || raw.getNumRows() != prev.getNumRows() || raw.getNumCols() != prev.getNumCols()) {
            return prev;
        }
        Matrix m = prev.copy();
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                double p = prev.get(i, j);
                double r = raw.get(i, j);
                double v = p + (r - p) / k;
                if (Double.isNaN(v)) {
                    v = Double.isNaN(r) ? p : r;
                }
                m.set(i, j, v);
            }
        }
        return m;
    }

    public void finish() {
        line_debug(options.verbose, String.format("LN finish: collecting final results from %d layers", this.ensemble.length));
        // In Robbins-Monro mode, report the Polyak-Ruppert averaged results
        // rather than the last (noisy) iterate
        if ("rm".equals(this.stochiterMode) && this.stochAvgCount > 0) {
            Map<Integer, SolverResult> lastRow = this.results.get(this.results.size());
            for (int e = 0; e < this.ensemble.length; e++) {
                SolverResult avg = this.stochAvg.get(e);
                SolverResult res = lastRow.get(e);
                if (avg != null && res != null) {
                    res.QN = avg.QN;
                    res.UN = avg.UN;
                    res.RN = avg.RN;
                    res.TN = avg.TN;
                    res.AN = avg.AN;
                    res.WN = avg.WN;
                }
            }
            if (this.stochServtAvg != null) {
                this.servt = this.stochServtAvg;
            }
            if (this.stochResidtAvg != null) {
                this.residt = this.stochResidtAvg;
            }
        }
        for (int e = 0; e < this.ensemble.length; e++) {
            solvers[e].getAvg();
        }
        //this.model.ensemble = this.ensemble; // not included as this comes through Diamond inheritance
    }

    public Matrix getArvproc_classes_updmap() {
        return arvproc_classes_updmap;
    }

    public AvgTable getAvgTable() {
        return jline.io.LineResultRecorder.around(this, "layered", () -> getAvgTableImpl());
    }

    /**
     * Body of {@link #getAvgTable()}, split out so {@link jline.io.LineResultRecorder}
     * sees what the getter RETURNED. The JAVA cross-codebase parity row is
     * measured from that rather than from what an example printed.
     */
    protected AvgTable getAvgTableImpl() {
        if (options != null && options.method != null
                && ("mwba.upper".equals(options.method) || "mwba.lower".equals(options.method))) {
            return getBoxBoundsTable(options.method);
        }
        return getEnsembleAvg();
    }

    /**
     * Response time distribution of every entry of the layered network.
     *
     * <p>Counterpart of MATLAB {@code @SolverLN/getCdfRespT.m}. The distribution is
     * formed by the {@code moment3} pass alone -- the mean-based update builds no law
     * at all -- so a solver constructed with any other method re-runs the ensemble
     * under {@code moment3} here and restores the caller's method afterwards. The
     * routing layers already built serve {@code moment3} unchanged, so only the
     * update pass changes.</p>
     *
     * @return one entry per entry of the LQN, in the entry-local index space
     *         ({@code lqn.eshift + i}). Each element is an (n x 2) matrix whose
     *         columns are [F(t), t], the column order every CDF getter in LINE uses,
     *         or null for an entry the pass fitted no law to.
     * @throws RuntimeException if the layers were built for a phase-type encoding,
     *         which carries no activity-graph routing to re-run over
     */
    public List<Matrix> getCdfRespT() {
        if (this.entrycdfrespt == null || !this.entrycdfrespt.containsKey(0)) {
            // The distribution pass reads the routing encoding of the activity graph,
            // which srvn.ph / flat.ph layers do not carry: re-running getEnsembleAvg
            // over them would reconstruct the wrong topology rather than a coarser
            // answer. Refuse by name.
            if (isPHEncoding()) {
                throw new RuntimeException("getCdfRespT needs the routing encoding of the activity "
                        + "graph, which method='" + this.lnmethod + "' does not build. Rebuild the "
                        + "solver with method='srvn.cs' or method='moment3'.");
            }
            String curMethod = this.options == null ? null : this.options.method;
            String curLnMethod = this.lnmethod;
            // BOTH the option and the RESOLVED method have to move: updateMetrics
            // dispatches on lnmethod, which buildLayers resolved once, so setting
            // options.method alone leaves the mean-based update in place and returns
            // an EMPTY table.
            if (this.options != null) {
                this.options.method = "moment3";
            }
            this.lnmethod = "moment3";
            try {
                getEnsembleAvg();
            } finally {
                if (this.options != null) {
                    this.options.method = curMethod;
                }
                this.lnmethod = curLnMethod;
            }
        }
        List<Matrix> cdfRespT = new ArrayList<Matrix>();
        for (int e = 0; e < this.lqn.nentries; e++) {
            cdfRespT.add(this.entrycdfrespt.get(e));
        }
        return cdfRespT;
    }

    /**
     * Layer-wise performance sensitivities of the layered network with respect to
     * service rates. Counterpart of MATLAB {@code @SolverLN/getSensitivityTable.m}.
     *
     * <p>Solves the layered model and then delegates to each layer solver, returning
     * the concatenation of the layer tables with a leading Layer column. Every row is
     * therefore a (Layer, Station, JobClass) triple carrying the derivative of that
     * row's mean measures with respect to that station-class service RATE in that
     * layer: dTput_dRate, dRespT_dRate, dQLen_dRate, dUtil_dRate.</p>
     *
     * <p>IMPORTANT, on what these derivatives mean. Each entry is a derivative WITHIN
     * ITS LAYER, taken with the layer parameters that the fixed point produced held
     * fixed. It is a partial derivative of the layer submodel, not the total
     * derivative of the layered model: perturbing a host demand in one layer moves
     * the think times, populations and service rates of the other layers through the
     * fixed-point map, and that indirect term is not included here. The layer table
     * is the right object for attributing a bottleneck inside a layer, and the wrong
     * one for predicting the effect of a parameter change on the solved layered
     * model. For the latter, finite-difference the LayeredNetwork itself.</p>
     *
     * @return the layer-wise sensitivity table
     * @see #getSensitivityTable(String, double, String)
     * @see LayeredNetworkSensitivityTable
     */
    public LayeredNetworkSensitivityTable getSensitivityTable() {
        return getSensitivityTable("auto", Double.NaN, "forward");
    }

    /**
     * Layer-wise performance sensitivities with an explicit branch, step and
     * difference scheme. The options are passed through to the layer solvers
     * unchanged, with the same meaning as in
     * {@link NetworkSolver#getSensitivityTable(String, double, String)}: each layer
     * independently takes the analytic branch where its own solver supports it and
     * the layer model is in scope, and finite differences otherwise.
     *
     * @param method one of "auto", "exact", "fd"
     * @param step   relative step of the rate perturbation; NaN selects the default
     * @param scheme "forward" or "central"
     * @return the layer-wise sensitivity table, whose
     *         {@link LayeredNetworkSensitivityTable#getMethod()} is "mixed" when the
     *         layers did not all take the same branch
     * @see #getSensitivityTable()
     */
    public LayeredNetworkSensitivityTable getSensitivityTable(String method, double step, String scheme) {
        // The layer solvers must sit at the converged fixed point, so the ensemble is
        // solved first if it has not been already.
        if (this.results == null || this.results.isEmpty()) {
            getEnsembleAvg();
        }

        int E = getNumberOfModels();
        List<String> layerName = new ArrayList<>();
        List<String> stationName = new ArrayList<>();
        List<String> className = new ArrayList<>();
        List<Double> dTput = new ArrayList<>();
        List<Double> dRespT = new ArrayList<>();
        List<Double> dQLen = new ArrayList<>();
        List<Double> dUtil = new ArrayList<>();
        List<String> layerMethods = new ArrayList<>();
        List<jline.io.Ret.pfqnSens> layerSens = new ArrayList<>();

        for (int e = 0; e < E; e++) {
            NetworkSolver solver = this.solvers[e];
            if (solver == null) {
                layerMethods.add(null);
                layerSens.add(null);
                continue;
            }
            NetworkSensitivityTable T = solver.getSensitivityTable(method, step, scheme);
            layerMethods.add(T.getMethod());
            layerSens.add(T.getSens());
            String thisLayer = this.ensemble[e].getName();
            List<String> st = T.getStationNames();
            List<String> cl = T.getClassNames();
            List<Double> tT = T.getDTput();
            List<Double> tR = T.getDRespT();
            List<Double> tQ = T.getDQLen();
            List<Double> tU = T.getDUtil();
            for (int i = 0; i < st.size(); i++) {
                layerName.add(thisLayer);
                stationName.add(st.get(i));
                className.add(cl.get(i));
                dTput.add(tT.get(i));
                dRespT.add(tR.get(i));
                dQLen.add(tQ.get(i));
                dUtil.add(tU.get(i));
            }
        }

        // One branch label per layer, plus a summary that is "mixed" when the layers
        // did not all take the same branch.
        String summary = "";
        for (int e = 0; e < layerMethods.size(); e++) {
            String m = layerMethods.get(e);
            if (m == null) {
                continue;
            }
            if (summary.isEmpty()) {
                summary = m;
            } else if (!summary.equals(m)) {
                summary = "mixed";
                break;
            }
        }

        LayeredNetworkSensitivityTable table =
                new LayeredNetworkSensitivityTable(dTput, dRespT, dQLen, dUtil);
        table.setOptions(this.options);
        table.setLayerNames(layerName);
        table.setStationNames(stationName);
        table.setClassNames(className);
        table.setLayerMethods(layerMethods);
        table.setLayerSens(layerSens);
        table.setMethod(summary);
        return table;
    }

    /**
     * Majumdar-Woodside robust box bounds table for the LQN (processor-contention
     * model). Mirrors the MATLAB SolverLN mwba.upper/mwba.lower path.
     */
    private AvgTable getBoxBoundsTable(String method) {
        jline.api.lqn.Lqn_boxbounds.Result bnd = jline.api.lqn.Lqn_boxbounds.compute(this.lqn);
        boolean upper = "mwba.upper".equals(method);
        int nidx = this.lqn.nidx;

        List<Double> Qval = new ArrayList<>();
        List<Double> Uval = new ArrayList<>();
        List<Double> Rval = new ArrayList<>();
        List<Double> Residval = new ArrayList<>();
        List<Double> Aval = new ArrayList<>();
        List<Double> Tval = new ArrayList<>();
        for (int i = 0; i < nidx; i++) {
            double t = upper ? bnd.TN_up[i] : bnd.TN_lo[i];
            double u = upper ? bnd.UN_up[i] : bnd.UN_lo[i];
            Tval.add(Double.isNaN(t) ? 0.0 : t);
            Uval.add(Double.isNaN(u) ? 0.0 : u);
            Qval.add(0.0);
            Rval.add(0.0);
            Residval.add(0.0);
            Aval.add(0.0);
        }

        List<String> nodeNames = new ArrayList<>(lqn.names.values());
        List<String> nodeTypes = new ArrayList<>();
        for (int o = 0; o < nodeNames.size(); o++) {
            switch ((int) lqn.type.get(o)) {
                case LayeredNetworkElement.PROCESSOR: nodeTypes.add("Processor"); break;
                case LayeredNetworkElement.TASK:
                    nodeTypes.add(lqn.sched.get(o) == SchedStrategy.REF ? "RefTask" : "Task");
                    break;
                case LayeredNetworkElement.ENTRY: nodeTypes.add("Entry"); break;
                case LayeredNetworkElement.ACTIVITY: nodeTypes.add("Activity"); break;
                case LayeredNetworkElement.CALL: nodeTypes.add("Call"); break;
            }
        }
        LayeredNetworkAvgTable AvgTable = new LayeredNetworkAvgTable(Qval, Uval, Rval, Residval, Aval, Tval);
        AvgTable.setNodeNames(nodeNames);
        AvgTable.setNodeTypes(nodeTypes);
        AvgTable.setOptions(this.options);
        return AvgTable;
    }

    /**
     * Transient average station metrics of the layered network.
     *
     * <p>Mirrors MATLAB SolverLN.getTranAvg: runs the ensemble fixed-point solve,
     * then delegates the transient analysis to each layer solver and assembles the
     * per-layer station x class traces block-diagonally (layer e in a disjoint
     * row/column block). Off-block cells are left null.</p>
     *
     * <p>Transient traces are only produced by transient-capable layer solvers
     * (Fluid, CTMC, SSA); with steady-state-only layers (MVA, NC) the delegated
     * getTranAvg throws, matching the MATLAB behaviour.</p>
     *
     * @return block-diagonal transient queue lengths, utilizations, throughputs
     *         and per-cell time vectors
     */
    public LNTranAvgResult getTranAvg() {
        String mode = (this.options != null && this.options.config != null
                && this.options.config.ln_transient != null)
                ? this.options.config.ln_transient : "coupled";
        if ("coupled".equalsIgnoreCase(mode)) {
            return getTranAvgCoupled();
        }
        if ("decoupled".equalsIgnoreCase(mode)) {
            return getTranAvgDecoupled();
        }
        line_error(mfilename(new Object(){}),
                "Unknown ln_transient mode '" + mode + "' (use 'coupled' or 'decoupled').");
        return null;
    }

    /**
     * Decoupled (frozen-demand) layered transient: the inter-layer demands stay
     * pinned at the converged fixed point and each layer's transient runs in
     * isolation. Mirrors MATLAB {@code SolverLN.getTranAvgDecoupled}.
     */
    public LNTranAvgResult getTranAvgDecoupled() {
        // Run the ensemble fixed point (mirrors self.getAvg in MATLAB).
        getEnsembleAvg();
        return assembleLayerTransients(null);
    }

    /**
     * Runs each layer's transient (optionally with an injected per-layer rate
     * schedule) and assembles the block-diagonal aggregate. The ensemble fixed
     * point is assumed to have been solved by the caller.
     *
     * @param schedByLayer per-layer {@code rate_sched} injections, or null for
     *                     the frozen-demand (decoupled) run
     */
    private LNTranAvgResult assembleLayerTransients(
            java.util.List<jline.solvers.fluid.handlers.FluidRateMultiplier.RateEntry>[] schedByLayer) {
        // Collect the per-layer transient traces from each layer solver.
        Matrix[][][] layerQ = new Matrix[nlayers][][];
        Matrix[][][] layerU = new Matrix[nlayers][][];
        Matrix[][][] layerT = new Matrix[nlayers][][];
        Matrix[] layerTime = new Matrix[nlayers];
        int totalRows = 0;
        int totalCols = 0;
        // The layered fixed-point solve above runs each layer in steady state
        // (layer getAvg rejects a timespan). The transient window therefore
        // lives on the SolverLN options and is applied to each layer solver only
        // around its transient getTranAvg call (mirrors MATLAB/Python).
        boolean hasTs = this.options != null && this.options.timespan != null
                && this.options.timespan.length >= 2
                && !isInf(this.options.timespan[0]) && !isInf(this.options.timespan[1]);
        // The fixed-point solve above hard-resets every layer state when it
        // detects convergence (converged() calls ensemble[e].reset), and the
        // layer solvers then re-initialize to the default state. That discards
        // any warm start installed with LayeredNetwork.initFromMarginal, so
        // replay it here: the transients below, unlike the steady-state solve,
        // depend on the initial state. This is what carries the queue lengths
        // of a SolverENV stage across an environment switch.
        boolean warmStarted = false;
        if (this.lqnModel != null && this.lqnModel.getInitMarginalBlocks() != null
                && this.lqnModel.getInitMarginalBlocks().size() == nlayers) {
            java.util.List<Matrix> blocks = this.lqnModel.getInitMarginalBlocks();
            for (int e = 0; e < nlayers; e++) {
                if (ensemble[e] == null) continue;
                Matrix block = new Matrix(blocks.get(e));
                if (!(solvers[e] instanceof SolverFluid)) {
                    // A discrete layer solver needs an integer state. Round it
                    // chain-wise: element-wise rounding would perturb the layer
                    // populations and thus the fixed point the transient relaxes to.
                    State.roundMarginalPreservingChains(block, ensemble[e].getStruct(false));
                }
                ensemble[e].initFromMarginal(block);
            }
            warmStarted = true;
        }
        for (int e = 0; e < nlayers; e++) {
            NetworkSolver s = solvers[e];
            if (s == null) {
                layerQ[e] = new Matrix[0][0];
                layerU[e] = new Matrix[0][0];
                layerT[e] = new Matrix[0][0];
                continue;
            }
            double[] savedTs = null;
            if (hasTs && s.options != null) {
                savedTs = s.options.timespan;
                s.options.timespan = new double[]{this.options.timespan[0], this.options.timespan[1]};
            }
            java.util.List<jline.solvers.fluid.handlers.FluidRateMultiplier.RateEntry> savedSched = null;
            boolean injected = false;
            if (schedByLayer != null && s.options != null && s.options.config != null
                    && schedByLayer[e] != null && !schedByLayer[e].isEmpty()) {
                savedSched = s.options.config.rate_sched;
                s.options.config.rate_sched = schedByLayer[e];
                injected = true;
                // Drop any cached transient: a layer solver returns its stored
                // result on repeat calls and does NOT re-solve on an
                // options.config change, so without this the injected
                // rate_sched would be silently ignored across
                // waveform-relaxation iterations. reset() only clears the
                // cached result, not the model's equilibrium service priming.
                s.reset();
            }
            if (warmStarted && !injected) {
                // Same reason as above: the cached steady-state result would
                // otherwise shadow the transient run from the replayed state.
                s.reset();
            }
            try {
                s.getTranAvg(); // populates s.result.QNt/UNt/TNt and s.result.t
            } finally {
                if (savedTs != null) {
                    s.options.timespan = savedTs;
                }
                if (injected) {
                    s.options.config.rate_sched = savedSched;
                }
            }
            SolverResult r = s.result;
            layerQ[e] = (r != null && r.QNt != null) ? r.QNt : new Matrix[0][0];
            layerU[e] = (r != null && r.UNt != null) ? r.UNt : new Matrix[0][0];
            layerT[e] = (r != null && r.TNt != null) ? r.TNt : new Matrix[0][0];
            layerTime[e] = (r != null) ? r.t : null;
            // QN drives the block extent, matching the MATLAB assembly.
            totalRows += layerQ[e].length;
            totalCols += (layerQ[e].length > 0) ? layerQ[e][0].length : 0;
        }

        Matrix[][] QNt = new Matrix[totalRows][totalCols];
        Matrix[][] UNt = new Matrix[totalRows][totalCols];
        Matrix[][] TNt = new Matrix[totalRows][totalCols];
        Matrix[][] t = new Matrix[totalRows][totalCols];

        int r0 = 0;
        int c0 = 0;
        for (int e = 0; e < nlayers; e++) {
            int nr = layerQ[e].length;
            int nc = (nr > 0) ? layerQ[e][0].length : 0;
            for (int i = 0; i < nr; i++) {
                for (int r = 0; r < nc; r++) {
                    QNt[r0 + i][c0 + r] = layerQ[e][i][r];
                    if (i < layerU[e].length && layerU[e][i] != null && r < layerU[e][i].length) {
                        UNt[r0 + i][c0 + r] = layerU[e][i][r];
                    }
                    if (i < layerT[e].length && layerT[e][i] != null && r < layerT[e][i].length) {
                        TNt[r0 + i][c0 + r] = layerT[e][i][r];
                    }
                    t[r0 + i][c0 + r] = layerTime[e];
                }
            }
            r0 += nr;
            c0 += nc;
        }

        // Resample every cell onto a single common time grid and store the
        // dense result into this.result, so SolverENV can consume an LN stage
        // through the same result.QNt/UNt/TNt + global result.t path it uses
        // for flat networks. Layer solvers use adaptive, per-layer grids;
        // without this the block-diagonal cells each carry a different time
        // vector and the off-block cells are null (which SolverENV.post cannot
        // handle). The zero-filled off-block cells match the MATLAB/Python
        // disabled-cell semantics.
        double tEnd = 0.0;
        int npts = 0;
        for (int e = 0; e < nlayers; e++) {
            if (layerTime[e] != null && layerTime[e].getNumRows() > 0) {
                tEnd = Math.max(tEnd, layerTime[e].get(layerTime[e].getNumRows() - 1, 0));
                npts = Math.max(npts, layerTime[e].getNumRows());
            }
        }
        if (npts > 1 && tEnd > 0.0) {
            npts = Math.min(npts, 200);
            Matrix tc = new Matrix(npts, 1);
            for (int p = 0; p < npts; p++) tc.set(p, 0, tEnd * p / (npts - 1));
            Matrix[][] rQ = new Matrix[totalRows][totalCols];
            Matrix[][] rU = new Matrix[totalRows][totalCols];
            Matrix[][] rT = new Matrix[totalRows][totalCols];
            for (int i = 0; i < totalRows; i++) {
                for (int j = 0; j < totalCols; j++) {
                    rQ[i][j] = interpColumn(t[i][j], QNt[i][j], tc);
                    rU[i][j] = interpColumn(t[i][j], UNt[i][j], tc);
                    rT[i][j] = interpColumn(t[i][j], TNt[i][j], tc);
                }
            }
            if (this.result == null) this.result = new SolverResult();
            this.result.QNt = rQ;
            this.result.UNt = rU;
            this.result.TNt = rT;
            this.result.t = tc;
        }

        return new LNTranAvgResult(QNt, UNt, TNt, t);
    }

    /**
     * Coupled layered transient by waveform relaxation over the LQN ensemble.
     *
     * <p>Port of MATLAB {@code @SolverLN/getTranAvgCoupled.m}. Unlike
     * {@link #getTranAvgDecoupled()}, which freezes inter-layer demands at the
     * converged fixed point, this reconciles the per-layer transients
     * iteratively: each layer's transient is driven by TIME-VARYING inter-layer
     * demand trajectories taken from the other layers' latest transients, and
     * the loop repeats until the trajectories stop changing (sup-norm gap over
     * time). The time-varying demands are injected into each layer solver
     * through the per-(station,class) rate schedule
     * ({@code options.config.rate_sched}), honoured by the fluid rate
     * multiplier and by the CTMC time-varying transient.</p>
     *
     * <p>Iteration 0 uses the frozen equilibrium demands, so it reproduces
     * {@link #getTranAvgDecoupled()} exactly; at convergence every layer relaxes
     * to its fixed point, so the endpoint equals {@code getEnsembleAvg}. The
     * return layout is the same block-diagonal (station x class per layer).</p>
     *
     * <p>Coupled channels: task think times (client delay) and
     * synchronous-call service demands (caller client station). Both are the
     * dominant inter-layer couplings; intra-layer host service stays at its
     * equilibrium value.</p>
     */
    @SuppressWarnings("unchecked")
    public LNTranAvgResult getTranAvgCoupled() {
        // Capture the transient horizon BEFORE the fixed-point solve: the
        // layered solve strips options.timespan (each layer getAvg rejects a
        // timespan), so reading it afterwards would spuriously trigger the
        // no-horizon fallback.
        boolean hasTs = this.options != null && this.options.timespan != null
                && this.options.timespan.length >= 2
                && Double.isFinite(this.options.timespan[0])
                && Double.isFinite(this.options.timespan[1]);
        double t0 = hasTs ? this.options.timespan[0] : 0.0;
        double t1 = hasTs ? this.options.timespan[1] : 0.0;

        getEnsembleAvg();

        if (!hasTs) {
            // No finite transient horizon: nothing to co-evolve, defer to decoupled.
            return assembleLayerTransients(null);
        }

        int maxit = (this.options.config != null) ? this.options.config.ln_transient_iter_max : 20;
        double tol = (this.options.config != null) ? this.options.config.ln_transient_tol : 1e-2;
        final int ngrid = 100;
        Matrix tgrid = new Matrix(ngrid, 1);
        for (int p = 0; p < ngrid; p++) {
            tgrid.set(p, 0, t0 + (t1 - t0) * p / (double) (ngrid - 1));
        }

        // Per-layer sn (for node->station bookkeeping), cached once.
        NetworkStruct[] layerSn = new NetworkStruct[nlayers];
        for (int e = 0; e < nlayers; e++) {
            layerSn[e] = (ensemble[e] != null) ? ensemble[e].getStruct(false) : null;
        }

        // Iteration 0: decoupled transients (frozen equilibrium demands already
        // set by the fixed-point solve).
        LNTranAvgResult out = assembleLayerTransients(null);
        double[][][][] traj = collectLayerTrajectories(tgrid);

        for (int it = 1; it <= maxit; it++) {
            double[][][][] trajPrev = traj;
            // 1) recompute inter-layer demand trajectories from the latest traj
            Map<Integer, double[]> thinkt = new HashMap<Integer, double[]>();
            Map<Integer, double[]> callservt = new HashMap<Integer, double[]>();
            recomputeCoupledDemand(ngrid, thinkt, callservt);
            // 2) build the per-layer rate_sched injections from those demands
            java.util.List<jline.solvers.fluid.handlers.FluidRateMultiplier.RateEntry>[] schedByLayer =
                    (java.util.List<jline.solvers.fluid.handlers.FluidRateMultiplier.RateEntry>[]) new java.util.List[nlayers];
            buildCoupledRateSched(thinkt, callservt, tgrid, layerSn, schedByLayer);
            // 3) re-run each layer with its injected time-varying demand
            out = assembleLayerTransients(schedByLayer);
            traj = collectLayerTrajectories(tgrid);
            // 4) convergence: sup-norm gap of the queue-length trajectories
            double gap = 0.0;
            for (int e = 0; e < nlayers; e++) {
                if (traj[e] == null || trajPrev[e] == null) continue;
                for (int i = 0; i < traj[e].length && i < trajPrev[e].length; i++) {
                    for (int r = 0; r < traj[e][i].length && r < trajPrev[e][i].length; r++) {
                        for (int p = 0; p < ngrid; p++) {
                            gap = Math.max(gap, Math.abs(traj[e][i][r][p] - trajPrev[e][i][r][p]));
                        }
                    }
                }
            }
            if (this.options.verbose == VerboseLevel.STD || this.options.verbose == VerboseLevel.DEBUG) {
                System.out.format("\nLN coupled transient: iter %d, sup-norm gap %.3e", it, gap);
            }
            if (gap < tol) {
                break;
            }
        }
        return out;
    }

    /**
     * Resamples every layer's last transient onto {@code tgrid}.
     *
     * <p>Returns the queue-length trajectories {@code traj[e][i][r][p]} (the
     * quantity the relaxation convergence test uses) and caches the matching
     * utilization, throughput and residence-time trajectories in
     * {@link #coupledU}, {@link #coupledT} and {@link #coupledR}, so the demand
     * recomputation reads exactly the traces the last layer run produced.</p>
     */
    private double[][][][] collectLayerTrajectories(Matrix tgrid) {
        coupledU = collectLayerMetric(tgrid, 'U');
        coupledT = collectLayerMetric(tgrid, 'T');
        coupledR = collectLayerMetric(tgrid, 'R');
        return collectLayerMetric(tgrid, 'Q');
    }

    /**
     * Resamples one transient metric of every layer onto {@code tgrid}.
     *
     * @param which 'Q' queue length, 'U' utilization, 'T' throughput, 'R'
     *              residence time (Q/T by Little's law)
     * @return {@code out[layer][station][class][gridpoint]}
     */
    private double[][][][] collectLayerMetric(Matrix tgrid, char which) {
        int ngrid = tgrid.getNumRows();
        double[][][][] out = new double[nlayers][][][];
        for (int e = 0; e < nlayers; e++) {
            NetworkSolver s = (e < solvers.length) ? solvers[e] : null;
            SolverResult r = (s != null) ? s.result : null;
            if (r == null || r.QNt == null) {
                out[e] = new double[0][][];
                continue;
            }
            int M = r.QNt.length;
            int K = (M > 0 && r.QNt[0] != null) ? r.QNt[0].length : 0;
            out[e] = new double[M][K][];
            for (int i = 0; i < M; i++) {
                for (int c = 0; c < K; c++) {
                    Matrix src;
                    switch (which) {
                        case 'U': src = (r.UNt != null && i < r.UNt.length && r.UNt[i] != null
                                && c < r.UNt[i].length) ? r.UNt[i][c] : null; break;
                        case 'T': src = (r.TNt != null && i < r.TNt.length && r.TNt[i] != null
                                && c < r.TNt[i].length) ? r.TNt[i][c] : null; break;
                        default: src = r.QNt[i][c]; break;
                    }
                    Matrix col;
                    if (which == 'R') {
                        Matrix q = interpColumn(r.t, r.QNt[i][c], tgrid);
                        Matrix tp = interpColumn(r.t, (r.TNt != null && i < r.TNt.length
                                && r.TNt[i] != null && c < r.TNt[i].length) ? r.TNt[i][c] : null, tgrid);
                        col = new Matrix(ngrid, 1);
                        for (int p = 0; p < ngrid; p++) {
                            col.set(p, 0, q.get(p, 0) / Math.max(tp.get(p, 0), GlobalConstants.FineTol));
                        }
                    } else {
                        col = interpColumn(r.t, src, tgrid);
                    }
                    double[] v = new double[ngrid];
                    for (int p = 0; p < ngrid; p++) {
                        v[p] = col.get(p, 0);
                    }
                    out[e][i][c] = v;
                }
            }
        }
        return out;
    }

    /**
     * Recomputes the time-varying inter-layer demands from the layer
     * trajectories, pointwise in t, mirroring the scalar updateThinkTimes /
     * updateMetricsDefault formulas.
     *
     * @param ngrid     number of grid points
     * @param thinkt    filled with tidx -&gt; think-time trajectory
     * @param callservt filled with cidx -&gt; call service-time trajectory
     */
    private void recomputeCoupledDemand(int ngrid,
                                        Map<Integer, double[]> thinkt,
                                        Map<Integer, double[]> callservt) {
        double[][][][] trajU = coupledU;
        double[][][][] trajT = coupledT;
        double[][][][] trajR = coupledR;

        // Task think times: from the task's own server-layer utilization/throughput.
        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            if (tidx >= idxhash.size() || Double.isNaN(idxhash.get(tidx))) {
                continue;
            }
            if (lqn.isref != null && lqn.isref.get(tidx) > 0) {
                continue;
            }
            int e = idxhash.get(tidx).intValue() - 1;
            if (e < 0 || e >= nlayers || trajU[e].length == 0) {
                continue;
            }
            int sIdx = stationIdxOf(e, tidx) - 1;
            if (sIdx < 0 || sIdx >= trajU[e].length) {
                continue;
            }
            double njobs = Matrix.extractRows(this.njobs, tidx, tidx + 1, null).elementMax();
            // same closure as updateThinkTimes, so the same gate -- see refThinkTime
            double userthink = refThinkTime(lqn, tidx);
            boolean isInfServer = lqn.sched.get(tidx) == SchedStrategy.INF;
            double[] tk = new double[ngrid];
            for (int p = 0; p < ngrid; p++) {
                double u = 0.0;
                double tp = 0.0;
                for (int r = 0; r < trajU[e][sIdx].length; r++) {
                    u += trajU[e][sIdx][r][p];
                    tp += trajT[e][sIdx][r][p];
                }
                double tsafe = Math.max(tp, GlobalConstants.FineTol);
                double v = isInfServer ? ((njobs - u) / tsafe - userthink)
                        : (njobs * Math.abs(1 - u) / tsafe - userthink);
                // total mean including the user think time
                tk[p] = Math.max(GlobalConstants.Zero, v) + userthink;
            }
            thinkt.put(tidx, tk);
        }

        // Synchronous-call service demands: callee entry response time * call mean.
        for (int cidx = 0; cidx < lqn.ncalls; cidx++) {
            if (lqn.calltype.get(cidx) != CallType.SYNC) {
                continue;
            }
            int eidx = (int) lqn.callpair.get(cidx, 1); // callee entry (callpair carries a padded leading column)
            int tidx = (int) lqn.parent.get(eidx);      // callee task
            if (tidx < 0 || tidx >= idxhash.size() || Double.isNaN(idxhash.get(tidx))) {
                continue;
            }
            int e = idxhash.get(tidx).intValue() - 1;
            if (e < 0 || e >= nlayers || trajR[e].length == 0) {
                continue;
            }
            int sIdx = stationIdxOf(e, tidx) - 1;
            if (sIdx < 0 || sIdx >= trajR[e].length) {
                continue;
            }
            // response time of the callee at its server, summed over the entry classes
            double[] rc = new double[ngrid];
            boolean any = false;
            List<JobClass> classes = ensemble[e].getClasses();
            for (int r = 0; r < trajR[e][sIdx].length && r < classes.size(); r++) {
                Integer[] attr = classes.get(r).getAttribute();
                if (attr != null && attr.length >= 2 && attr[0] != null && attr[1] != null
                        && attr[0].intValue() == LayeredNetworkElement.ENTRY && attr[1].intValue() == eidx) {
                    for (int p = 0; p < ngrid; p++) {
                        rc[p] += trajR[e][sIdx][r][p];
                        if (rc[p] != 0.0) any = true;
                    }
                }
            }
            if (!any) {
                // fall back to the entry's activities response time
                for (int r = 0; r < trajR[e][sIdx].length; r++) {
                    for (int p = 0; p < ngrid; p++) {
                        rc[p] += trajR[e][sIdx][r][p];
                    }
                }
            }
            double callmean = lqn.callproc_mean.getOrDefault(cidx, 1.0);
            double[] d = new double[ngrid];
            for (int p = 0; p < ngrid; p++) {
                d[p] = rc[p] * callmean;
            }
            callservt.put(cidx, d);
        }
    }

    // Per-layer utilization / throughput / residence trajectories captured
    // alongside the queue lengths by collectLayerTrajectories; kept as fields so
    // the demand recomputation reads exactly the traces the last layer run
    // produced.
    private double[][][][] coupledU;
    private double[][][][] coupledT;
    private double[][][][] coupledR;

    /**
     * Maps the recomputed demand trajectories to per-layer {@code rate_sched}
     * injections, using the same update maps updateLayers uses to place
     * setService calls.
     *
     * <p>{@code options.config.ln_transient_channels} selects which inter-layer
     * coupling channels are injected: {@code "both"} (default), {@code "thinkt"}
     * (client-delay only), or {@code "callservt"} (synchronous-call service
     * only). Used to isolate each channel's contribution.</p>
     */
    private void buildCoupledRateSched(Map<Integer, double[]> thinkt,
                                       Map<Integer, double[]> callservt,
                                       Matrix tgrid, NetworkStruct[] layerSn,
                                       java.util.List<jline.solvers.fluid.handlers.FluidRateMultiplier.RateEntry>[] schedByLayer) {
        String channels = (this.options.config != null && this.options.config.ln_transient_channels != null)
                ? this.options.config.ln_transient_channels.toLowerCase() : "both";

        // think-time channel (client delay of caller tasks)
        if (("both".equals(channels) || "thinkt".equals(channels)) && thinkt_classes_updmap != null) {
            for (int row = 1; row < thinkt_classes_updmap.getNumRows(); row++) {
                // the updmap Matrix carries a padded leading row AND column, so
                // the fields live at rows 1..n, columns 1..4 (as updateLayers)
                int idx = (int) thinkt_classes_updmap.get(row, 1);
                int aidx = (int) thinkt_classes_updmap.get(row, 2);
                int nodeidx = (int) thinkt_classes_updmap.get(row, 3);
                int classidx = (int) thinkt_classes_updmap.get(row, 4);
                if (idx >= idxhash.size() || Double.isNaN(idxhash.get(idx))) {
                    continue;
                }
                int e = idxhash.get(idx).intValue() - 1;
                if (e < 0 || e >= nlayers || nodeidx != ensemble[e].getAttribute().getClientIdx()) {
                    continue;
                }
                if (lqn.type.get(aidx) == LayeredNetworkElement.TASK
                        && lqn.sched.get(aidx) != SchedStrategy.REF && thinkt.containsKey(aidx)) {
                    addCoupledSched(schedByLayer, e, layerSn[e], nodeidx, classidx, tgrid, thinkt.get(aidx));
                }
            }
        }

        // call-service channel (client station of caller for each sync call)
        if (("both".equals(channels) || "callservt".equals(channels)) && call_classes_updmap != null) {
            for (int row = 1; row < call_classes_updmap.getNumRows(); row++) {
                int idx = (int) call_classes_updmap.get(row, 1);
                int cidx = (int) call_classes_updmap.get(row, 2);
                int nodeidx = (int) call_classes_updmap.get(row, 3);
                int classidx = (int) call_classes_updmap.get(row, 4);
                if (idx >= idxhash.size() || Double.isNaN(idxhash.get(idx))) {
                    continue;
                }
                int e = idxhash.get(idx).intValue() - 1;
                if (e < 0 || e >= nlayers || nodeidx != ensemble[e].getAttribute().getClientIdx()) {
                    continue;
                }
                if (callservt.containsKey(cidx)) {
                    addCoupledSched(schedByLayer, e, layerSn[e], nodeidx, classidx, tgrid, callservt.get(cidx));
                }
            }
        }
    }

    /**
     * Appends one {@code rate_sched} entry that MODULATES the layer's
     * equilibrium rate by the ratio of the transient demand to its steady-state
     * (end-of-horizon) value: {@code effective_rate(t) = nominal*demand(end)/demand(t)}.
     *
     * <p>Passing {@code rates = 1/demand(t)} and {@code nominal = 1/demand(end)}
     * makes the multiplier {@code demand(end)/demand(t)}, which is exactly 1 at
     * the horizon end, so the layer relaxes to its unmodified fixed point
     * (endpoint == getEnsembleAvg) regardless of any small mismatch between the
     * transient residence Q/T and the scalar equilibrium demand.</p>
     */
    private static void addCoupledSched(java.util.List<jline.solvers.fluid.handlers.FluidRateMultiplier.RateEntry>[] schedByLayer,
                                        int e, NetworkStruct sn, int nodeidx, int classidx,
                                        Matrix tgrid, double[] demand) {
        if (sn == null || sn.nodeToStation == null || nodeidx - 1 >= sn.nodeToStation.length()) {
            return;
        }
        int ist = (int) sn.nodeToStation.get(nodeidx - 1);
        if (ist < 0) {
            return;
        }
        int n = demand.length;
        double dend = demand[n - 1];
        if (!(dend > GlobalConstants.FineTol)) {
            return; // degenerate steady-state demand; skip this channel
        }
        // Bound the transient demand to a physical band around its steady-state
        // value. The Little's-law demand estimates (thinkt = njobs*(1-util)/tput,
        // R = Q/T) diverge at early t when the layer transient throughput is
        // still ~0 (no jobs have completed yet), producing unbounded multipliers.
        // An inter-layer demand cannot realistically deviate from its equilibrium
        // by an unbounded factor; the fluid ODE tolerates the spikes but the CTMC
        // propagation does not (an enormous generator scaling ejects all mass
        // from a state within one segment). Clamping d to [dend/Cap, dend*Cap]
        // keeps the multiplier m = dend/d in [1/Cap, Cap] and is exactly 1 at the
        // horizon end (endpoint preserved).
        final double cap = 20.0;
        double[] tg = new double[n];
        double[] rates = new double[n];
        for (int p = 0; p < n; p++) {
            tg[p] = tgrid.get(p, 0);
            double d = Math.min(Math.max(demand[p], dend / cap), dend * cap);
            rates[p] = 1.0 / d;
        }
        if (schedByLayer[e] == null) {
            schedByLayer[e] = new ArrayList<jline.solvers.fluid.handlers.FluidRateMultiplier.RateEntry>();
        }
        schedByLayer[e].add(new jline.solvers.fluid.handlers.FluidRateMultiplier.RateEntry(
                ist, classidx - 1, tg, rates, Double.valueOf(1.0 / dend)));
    }

    /** Linear interpolation of the (tsrc, ysrc) series onto grid tc (clamped at
     *  the endpoints); returns a zero column when the source is null/empty. */
    private static Matrix interpColumn(Matrix tsrc, Matrix ysrc, Matrix tc) {
        int n = tc.getNumRows();
        Matrix out = new Matrix(n, 1);
        if (tsrc == null || ysrc == null) return out;
        int L = Math.min(tsrc.getNumRows(), ysrc.getNumRows());
        if (L == 0) return out;
        for (int p = 0; p < n; p++) {
            double x = tc.get(p, 0);
            double val;
            if (x <= tsrc.get(0, 0)) {
                val = ysrc.get(0, 0);
            } else if (x >= tsrc.get(L - 1, 0)) {
                val = ysrc.get(L - 1, 0);
            } else {
                int lo = 0;
                while (lo < L - 1 && tsrc.get(lo + 1, 0) < x) lo++;
                double x0 = tsrc.get(lo, 0), x1 = tsrc.get(lo + 1, 0);
                double y0 = ysrc.get(lo, 0), y1 = ysrc.get(lo + 1, 0);
                double f = (x1 > x0) ? (x - x0) / (x1 - x0) : 0.0;
                val = y0 + f * (y1 - y0);
            }
            out.set(p, 0, val);
        }
        return out;
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

    /**
     * Snaps near-tenth entries onto the tenth and flattens near-zero entries,
     * mirroring the sanitization MATLAB SolverLN applies before formatting
     * getAvgTable. NaN entries are left untouched, as are negative ones (the
     * MATLAB relative test is vacuous there).
     */
    private static void sanitizeAvg(Matrix m) {
        if (m == null) {
            return;
        }
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                double v = m.get(i, j);
                if (Double.isNaN(v)) {
                    continue;
                }
                double scaled = v * 10.0;
                double rounded = FastMath.round(scaled);
                if (FastMath.abs(scaled - rounded) < GlobalConstants.CoarseTol * scaled) {
                    v = rounded / 10.0;
                    m.set(i, j, v);
                }
                if (v <= GlobalConstants.FineTol) {
                    m.set(i, j, 0);
                }
            }
        }
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
        line_debug(options.verbose, String.format("LN: starting ensemble iteration with %d layers", nlayers));
        try {
            this.iterate();
        } catch (Exception e) {
            line_warning(mfilename(new Object(){}.getClass().getEnclosingMethod()),
                "SolverLN iteration failed: %s at %s", e.getMessage(),
                e.getStackTrace().length > 0 ? e.getStackTrace()[0].toString() : "unknown");
            e.printStackTrace();
            return null;
        }
        // the layers of method 'srvn.ph' carry one class per caller task, so the
        // per-element results are rebuilt analytically -- see getEnsembleAvgPH
        if (isPHEncoding()) {
            Matrix[] r = getEnsembleAvgPH();
            return assembleAvgTable(r[0], r[1], r[2], r[3], r[4], r[5]);
        }

        // NOTE: TestSolverLN, TestSolverLN2, TestSolverLN3, had problems here due to
        // different values returned by getAvg() in MATLAB and LINE on these examples

        // Check if results were populated during iteration
        if (this.results == null || this.results.isEmpty()) {
            line_warning(mfilename(new Object(){}.getClass().getEnclosingMethod()),
                "SolverLN iteration did not produce results. Returning null AvgTable.");
            return null;
        }

        Matrix QN = new Matrix(1, this.lqn.nidx, this.lqn.nidx);
        for (int i = 0; i < this.lqn.nidx; i++)
            QN.set(i, Double.NaN);

        Matrix UN = new Matrix(1, this.lqn.nidx, this.lqn.nidx);
        for (int i = 0; i < this.lqn.nidx; i++)
            UN.set(i, Double.NaN);

        Matrix RN = new Matrix(1, this.lqn.nidx, this.lqn.nidx);
        for (int i = 0; i < this.lqn.nidx; i++)
            RN.set(i, Double.NaN);

        Matrix TN = new Matrix(1, this.lqn.nidx, this.lqn.nidx);
        for (int i = 0; i < this.lqn.nidx; i++)
            TN.set(i, Double.NaN);

        // utilization will be first stored here
        Matrix PN = new Matrix(1, this.lqn.nidx, this.lqn.nidx);
        for (int i = 0; i < this.lqn.nidx; i++)
            PN.set(i, Double.NaN);

        // response time will be first stored here
        Matrix SN = new Matrix(1, this.lqn.nidx, this.lqn.nidx);
        for (int i = 0; i < this.lqn.nidx; i++)
            SN.set(i, Double.NaN);

        //residence time
        Matrix WN = new Matrix(1, this.lqn.nidx, this.lqn.nidx);
        for (int i = 0; i < this.lqn.nidx; i++)
            WN.set(i, Double.NaN);

        // not available yet
        Matrix AN = new Matrix(1, this.lqn.nidx, this.lqn.nidx);
        for (int i = 0; i < this.lqn.nidx; i++)
            AN.set(i, Double.NaN);

        // Track which activities have already been accumulated into task WN
        // to prevent double-counting when an activity appears in multiple layers
        boolean[] wnProcessed = new boolean[this.lqn.nidx + 1];

        int E = this.nlayers;

        for (int e = 0; e < E; e++) {
            int clientIdx = this.ensemble[e].getAttribute().getClientIdx();
            int serverIdx = this.ensemble[e].getAttribute().getServerIdx();
            int sourceIdx = this.ensemble[e].getAttribute().getSourceIdx();

            // determine processor metrics, one processor at a time under flat layering
            List<Integer> hostStations = serverStationsOf(e, true);
            boolean hasHostServer = !hostStations.isEmpty();
            for (int hs : hostStations) {
                Queue q = (Queue) this.ensemble[e].getStations().get(hs - 1);
                int hidx = q.getAttribute().getIdx();
                TN.set(0, hidx, 0);
                PN.set(0, hidx, 0);
                for (int c = 0; c < this.ensemble[e].getNumberOfClasses(); c++) {
                    if (this.ensemble[e].getClasses().get(c).getCompletes()) {
                        double t = 0;
                        if (clientIdx != -1) {
                            t = FastMath.max(t, this.results.get(this.results.size()).get(e).TN.get(clientIdx - 1, c));
                        }
                        if (sourceIdx != -1) {
                            t = FastMath.max(t, this.results.get(this.results.size()).get(e).TN.get(sourceIdx - 1, c));
                        }
                        TN.set(0, hidx, TN.get(hidx) + FastMath.max(t, this.results.get(this.results.size()).get(e).TN.get(hs - 1, c)));
                    }
                    int type = this.ensemble[e].getClasses().get(c).getAttribute()[0];
                    if (type == LayeredNetworkElement.ACTIVITY) {
                        if (stationIdxOfClass(e, c) != hs) {
                            continue; // the activity does not run on this processor
                        }
                        int aidx = this.ensemble[e].getClasses().get(c).getAttribute()[1];
                        int tidx = (int) this.lqn.parent.get(0, aidx);
                        if (Double.isNaN(PN.get(aidx))) PN.set(0, aidx, 0);
                        if (tidx >= 0) {
                            if (Double.isNaN(PN.get(tidx))) PN.set(0, tidx, 0);
                            PN.set(0, aidx, PN.get(aidx) + this.results.get(this.results.size()).get(e).UN.get(hs - 1, c));
                            PN.set(0, tidx, PN.get(tidx) + this.results.get(this.results.size()).get(e).UN.get(hs - 1, c));
                        } else {
                            PN.set(0, aidx, PN.get(aidx) + this.results.get(this.results.size()).get(e).UN.get(hs - 1, c));
                        }
                        PN.set(0, hidx, PN.get(hidx) + this.results.get(this.results.size()).get(e).UN.get(hs - 1, c));
                    }
                }
                TN.set(0, hidx, Double.NaN); // Added for consistency with LQNS
            }

            //determine remaining metrics
            // Declare variables outside switch to avoid scope issues
            int tidx, eidx, cidx, aidx;
            Station task_s, entry_s;
            Queue task_q, entry_q;

            if (serverIdx != -1) {
            for (int c = 0; c < this.ensemble[e].getNumberOfClasses(); c++) {
                int type = this.ensemble[e].getClasses().get(c).getAttribute()[0];
                // under flat layering each class is served at its own station, so
                // read the layer result there rather than at the layer's serverIdx
                serverIdx = stationIdxOfClass(e, c);
                switch (type) {
                    case LayeredNetworkElement.TASK:
                        tidx = this.ensemble[e].getClasses().get(c).getAttribute()[1];
                        if (hasHostServer) {
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
                        if (this.hasPhase2 && this.servt_ph2 != null && this.servt_ph2.get(eidx) > GlobalConstants.FineTol) {
                            SN.set(eidx, this.residt.get(eidx));  // Phase-1 + overtaking correction
                        } else {
                            SN.set(eidx, this.servt.get(eidx));
                        }
                        if (hasHostServer) {
                            if (Double.isNaN(TN.get(eidx)) && clientIdx != -1) {
                                // store the result in th eprocessor model
                                TN.set(eidx, this.results.get(this.results.size()).get(e).TN.get(clientIdx - 1, c));
                            }
                        }
                        break;

                    case LayeredNetworkElement.CALL:
                        cidx = this.ensemble[e].getClasses().get(c).getAttribute()[1];
                        aidx = (int) this.lqn.callpair.get(cidx, 0);
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
                        WN.set(aidx, this.residt.get(aidx));
                        if (!wnProcessed[aidx]) {
                            WN.set(tidx, WN.get(tidx) + this.residt.get(aidx));
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
        for (int e = 0; e < this.lqn.nentries; e++) {
            int eidx = this.lqn.eshift + e;
            int tidx = (int) this.lqn.parent.get(0, eidx);
            if (Double.isNaN(UN.get(tidx))) UN.set(tidx, 0);

            // Phase-2 support: utilization includes both phases
            if (this.hasPhase2 && this.servt_ph2 != null && this.servt_ph2.get(eidx) > GlobalConstants.FineTol) {
                // Phase-1 utilization
                this.util_ph1.set(eidx, TN.get(eidx) * this.servt_ph1.get(eidx));
                // Phase-2 utilization
                this.util_ph2.set(eidx, TN.get(eidx) * this.servt_ph2.get(eidx));
                // Total utilization = both phases (server is busy during both)
                UN.set(eidx, this.util_ph1.get(eidx) + this.util_ph2.get(eidx));
            } else {
                // Standard calculation for entries without phase-2
                UN.set(eidx, TN.get(eidx) * SN.get(eidx));
            }

            // Entry utilization = sum of activity processor utilizations for that entry
            List<Integer> entryActs = this.lqn.actsof.get(eidx);
            if (entryActs != null) {
                double entryUtil = 0;
                for (int aidx : entryActs) {
                    if (!Double.isNaN(PN.get(aidx))) {
                        entryUtil += PN.get(aidx);
                    }
                }
                PN.set(eidx, entryUtil);
            }

            if (tidx >= 0) {
                for (int i = 0; i < this.lqn.actsof.get(tidx).size(); i++) {
                    int aidx = this.lqn.actsof.get(tidx).get(i);
                    UN.set(aidx, TN.get(aidx) * SN.get(aidx));
                }
                UN.set(tidx, UN.get(tidx) + UN.get(eidx));
            }
        }

        // AN IGNORED ELEMENT IS IDLE, NOT UNDEFINED, and the two are different
        // cells. Its component holds no reference task, so nothing reaches it and
        // every measure it HAS is zero -- but the measures its kind never has stay
        // NaN, exactly as they do for a reachable element. A flat zero over all six
        // columns broke the table's NaN mask (a processor with a queue length of 0,
        // an arrival rate reported where no solver reports one), and the mask is
        // part of the answer: see _kb/06-solver-catalog.md. Reported columns are
        // QLen=UN, Util=PN, RespT=SN, ResidT=WN, ArvR=AN, Tput=TN, so the pre-swap
        // QN and RN are discarded below and are not written here.
        for (double idx : this.ignore.find().toList1D()) {
            int idxInt = (int) idx;
            if (idxInt >= 0 && idxInt < this.lqn.nidx) {
                PN.set(idxInt, 0.0);            // every kind reports a utilization
                AN.set(idxInt, Double.NaN);     // nothing reports an arrival rate on an LQN
                switch ((int) this.lqn.type.get(idxInt)) {
                    case LayeredNetworkElement.PROCESSOR:
                        UN.set(idxInt, Double.NaN);
                        SN.set(idxInt, Double.NaN);
                        WN.set(idxInt, Double.NaN);
                        TN.set(idxInt, Double.NaN);
                        break;
                    case LayeredNetworkElement.TASK:
                        UN.set(idxInt, 0.0);
                        SN.set(idxInt, Double.NaN);
                        WN.set(idxInt, 0.0);
                        TN.set(idxInt, 0.0);
                        break;
                    case LayeredNetworkElement.ENTRY:
                        UN.set(idxInt, 0.0);
                        SN.set(idxInt, 0.0);
                        WN.set(idxInt, Double.NaN);
                        TN.set(idxInt, 0.0);
                        break;
                    case LayeredNetworkElement.ACTIVITY:
                        UN.set(idxInt, 0.0);
                        SN.set(idxInt, 0.0);
                        WN.set(idxInt, 0.0);
                        TN.set(idxInt, 0.0);
                        break;
                    default:
                        break;
                }
            }
        }

        // if LN standard naming
        QN = UN.copy();
        UN = PN.copy();
        RN = SN.copy();
        return assembleAvgTable(QN, UN, RN, TN, AN, WN);
    }

    /**
     * Build the LQN average table from the six per-element rows, shared by the
     * default reconstruction and by method 'srvn.ph'.
     *
     * @param QN queue lengths (task and entry utilizations)
     * @param UN utilizations (processor utilizations)
     * @param RN response times
     * @param TN throughputs
     * @param AN arrival rates
     * @param WN residence times
     * @return the average table
     */
    private AvgTable assembleAvgTable(Matrix QN, Matrix UN, Matrix RN, Matrix TN, Matrix AN, Matrix WN) {
        // Sanitize small numerical perturbations, matching MATLAB SolverLN
        // getAvgTable: a value within CoarseTol of a tenth is snapped onto it,
        // and anything at or below FineTol is flattened to zero. Without this the
        // LN fixed-point residual leaks into the printed table (e.g. a residence
        // time of 0.20009 where the model demand is exactly 0.2), which reads as
        // a cross-language divergence even though the iterates agree.
        sanitizeAvg(QN);
        sanitizeAvg(UN);
        sanitizeAvg(RN);
        sanitizeAvg(TN);
        sanitizeAvg(AN);
        sanitizeAvg(WN);

        int maxnamelength = 12;
        for (int i = 0; i < lqn.names.size(); i++) {
            maxnamelength = FastMath.max(maxnamelength, lqn.names.get(i).length());
        }
        /* LayeredNetworkAvgTable Generation Boilerplate */
        List<String> nodeNames = new ArrayList<>(lqn.names.values());
        List<String> nodeTypes = new ArrayList<>();

        for (int o = 0; o < nodeNames.size(); o++) {
            switch ((int) lqn.type.get(o)) {
                case LayeredNetworkElement.PROCESSOR:
                    nodeTypes.add("Processor");
                    break;
                case LayeredNetworkElement.TASK:
                    if (lqn.sched.get(o) == SchedStrategy.REF) {
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

        LayeredNetworkAvgTable AvgTable = new LayeredNetworkAvgTable(Qval, Uval, Rval, Residval, Aval, Tval);
        AvgTable.setNodeNames(nodeNames);
        AvgTable.setNodeTypes(nodeTypes);
        AvgTable.setOptions(this.options);
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
        for (int e = 0; e < this.lqn.nentries; e++) {
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
                //if the successor activity is a call
                //entries (e.g., forwarding sources) have no callsof entry (MATLAB: empty cell)
                List<Integer> callsofAidx = lqn.callsof.get(aidx) != null ? lqn.callsof.get(aidx) : new ArrayList<Integer>();
                for (int j = 0; j < callsofAidx.size(); j++) {
                    int cidx = callsofAidx.get(j);
                    if (lqn.calltype.get(cidx) == CallType.SYNC) {
                        // mean number of calls alrady factored in
                        U.set(eidx, lqn.nidx + cidx, 1);
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
                    U.set(eidx, nextaidx, U.get(eidx, nextaidx) + Gvalue);
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
        // The moment3 pass is terminal WITHIN ONE SOLVE, so the flag is scoped to
        // one iterate(): left standing, the terminal test in converged() fires at
        // it=0 on the NEXT solve, the loop body never runs and every metric comes
        // back zero. See BUGS.md BUG-97.
        this.momentPassDone = false;
        line_debug(options.verbose, String.format("LN init: nlayers=%d, iter_max=%d, iter_tol=%e",
            nlayers, options.iter_max, options.iter_tol));

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
        // Keep the feature-set gate armed on the layer solvers: a layer whose model a
        // solver cannot represent must be rejected, not solved into silently wrong
        // numbers (e.g. MVA has no notion of a Signal class and would return a
        // product-form answer in which no job is ever removed).
        for (int e = 0; e < this.nlayers; e++)
            this.solvers[e].enableChecks = true;

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
        this.callresidt_prev = new Matrix(1, this.lqn.ncalls, this.lqn.ncalls);
        this.callresidt_prev.fill(Double.NaN);

        // Optional initialization of the layer throughput state from the
        // Majumdar-Woodside robust box bounds (geometric-mean point estimate
        // sqrt(Xlo*Xup)). Selected with options.config.layer_init = "bound". For
        // closed-chain LQNs the outer iteration re-derives throughputs from the
        // first solve, so this does not alter the converged result.
        String initMode = this.options.config.layer_init;
        if (initMode != null && (initMode.equalsIgnoreCase("bound")
                || initMode.equalsIgnoreCase("boxbound") || initMode.equalsIgnoreCase("mwba"))) {
            try {
                jline.api.lqn.Lqn_boxbounds.Result bnd = jline.api.lqn.Lqn_boxbounds.compute(this.lqn);
                for (int idx = 0; idx < this.lqn.nidx; idx++) {
                    double u = bnd.TN_up[idx];
                    double l = bnd.TN_lo[idx];
                    double x;
                    if (!Double.isNaN(u) && !Double.isNaN(l) && u > 0 && l > 0) {
                        x = Math.sqrt(u * l);
                    } else if (!Double.isNaN(u)) {
                        x = u;
                    } else if (!Double.isNaN(l)) {
                        x = l;
                    } else {
                        x = 0;
                    }
                    if (x > 0) {
                        this.tput.set(idx, x);
                        this.tputproc.put(idx, Exp.fitRate(x));
                    }
                }
                line_debug(options.verbose, "LN init: throughputs initialized from robust box bounds");
            } catch (Exception ex) {
                line_debug(options.verbose, "LN box-bound initialization skipped: " + ex.getMessage());
            }
        }

        // Resolve the stochastic iteration mode. Simulation-based or Monte
        // Carlo based layer solvers observe the layer map only up to noise,
        // for which the deterministic Picard iteration and its successive-
        // difference test are inadequate (see convergedStoch). The static
        // classification below uses options.method; a layer running method
        // 'default' may still resolve to a stochastic method at runtime, so
        // analyze() refreshes stochlayers after the first iteration and
        // converged() upgrades an 'auto' mode accordingly.
        this.stochlayers = new boolean[this.nlayers];
        for (int e = 0; e < this.nlayers; e++) {
            this.stochlayers[e] = this.solvers[e] != null && this.solvers[e].isStochastic();
        }
        String stochMode = this.options.config.stochiter;
        if (stochMode == null) stochMode = "auto";
        this.stochiterAuto = stochMode.equalsIgnoreCase("auto");
        if (this.stochiterAuto) {
            stochMode = anyStochasticLayer() ? "rm" : "off";
        }
        this.stochiterMode = stochMode.toLowerCase();
        this.stochiterStart = null;
        this.stochAvg = new HashMap<>();
        this.stochAvgCount = 0;
        this.stochServtAvg = null;
        this.stochResidtAvg = null;
    }

    private boolean anyStochasticLayer() {
        if (this.stochlayers == null) {
            return false;
        }
        for (boolean isStoch : this.stochlayers) {
            if (isStoch) {
                return true;
            }
        }
        return false;
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
        line_debug(options.verbose, String.format("LN post: iteration %d", it));

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
            // refreshChains can change the chain basis, invalidating the
            // warm-start solution cached by analyze()
            solvers[e - 1].options.init_sol = new Matrix(0, 0);
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
                // Mirrors MATLAB post(): only 'moment3' needs the full process
                // refresh; every other method takes the leaner rate refresh (no
                // need to refresh phases). Naming the methods explicitly instead
                // left an unrecognized one with no refresh at all, so the layers
                // never updated and the iteration converged on wrong numbers.
                if ("moment3".equals(this.lnmethod) || isPHEncoding()) {
                    // both carry a phase-type service law, whose phases a
                    // rate-only refresh would drop
                    ensemble[e - 1].refreshProcesses();
                } else {
                    ensemble[e - 1].refreshRates(null, null);
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
            // Keep the gate armed: see the note in the init() loop above.
            for (int e = 0; e < ensemble.length; e++) {
                solvers[e].enableChecks = true;
            }
        }
    }

    public void pre(int it) {
        // Seed control for stochastic layer solvers
        if (this.stochiterMode == null || this.stochlayers == null) {
            return;
        }
        if ("rm".equals(this.stochiterMode)) {
            // Rotate seeds so successive iterations observe the layer map
            // under independent noise, as required for Robbins-Monro
            // averaging to reduce variance
            for (int e = 0; e < this.nlayers; e++) {
                if (this.stochlayers[e]) {
                    this.solvers[e].options.seed = this.options.seed + (it - 1) * this.nlayers + e + 1;
                }
            }
        } else if ("crn".equals(this.stochiterMode)) {
            // Common random numbers: pin a constant per-layer seed so each
            // layer map is deterministic given its seed (sample average
            // approximation). The standard convergence test then applies to
            // the sample-average fixed point, which carries an
            // O(1/sqrt(samples)) bias with respect to the true fixed point.
            for (int e = 0; e < this.nlayers; e++) {
                if (this.stochlayers[e]) {
                    this.solvers[e].options.seed = this.options.seed + e + 1;
                }
            }
        }
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
        // Under 'srvn.ph' the layer classes are one per caller task and their laws
        // are composed, not read off the update maps -- see _kb/06-solver-catalog.md
        if (isPHEncoding()) {
            updateLayersPH(it);
            return;
        }
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
            if (aidx < (this.lqn.tshift + this.lqn.ntasks)) {
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
            // Case 2 - server replica (any of them)
            else {
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
                node.setArrival(tmp_class, this.tputproc.get((int) this.lqn.callpair.get(cidx, 0)));
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
            // Case 2 - server replica (any of them)
            else {
                int eidx = (int) this.lqn.callpair.get((int) cidx, 1);
                // A phase-2 entry replies before phase 2 runs, so the caller is held
                // for residt, not servt. Charging the caller servt here while its own
                // layer charges residt makes the two layers settle at different rates
                // and breaks flow conservation across the call. Under flat layering
                // both live in one model, where the correction would be applied twice.
                if (this.hasPhase2 && this.servt_ph2 != null
                        && this.servt_ph2.get(eidx) > GlobalConstants.FineTol
                        && !isFlatLayering()) {
                    node.setService(tmp_class, Exp.fitMean(this.residt.get(eidx)));
                } else {
                    node.setService(tmp_class, this.servtproc.get(eidx));
                }
            }
        }
    }

    public void updateMetrics(int it) {
        // Mirrors MATLAB updateMetrics.m: only 'moment3' propagates three moments
        // of the response time distribution; every other method name ('default',
        // 'mva', 'nc', or an unrecognized one) takes the default update.
        // Whitelisting names here instead left the metrics unset for any other
        // method, so residt stayed null and getEnsembleAvg died with an NPE.
        if (isPHEncoding()) {
            // see _kb/06-solver-catalog.md (LN section) for rationale
            updateMetricsPH(it);
        } else if ("moment3".equals(lnmethod)) {
            this.updateMetricsMomentBased(it);
        } else {
            this.updateMetricsDefault(it);
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
            int iter_min = (int) FastMath.min(30, FastMath.ceil(this.options.iter_max / 4.0));
            if (this.averagingstart != null && it >= iter_min) {
                int wnd_size = it - this.averagingstart + 1;
                this.servt.set(aidx, 0);
                this.residt.set(aidx, 0);
                this.tput.set(aidx, 0);
                for (int w = 0; w < wnd_size; w++) {
                    this.servt.set(aidx, this.servt.get(aidx) + this.results.get(rLen - w).get(layerIdx_0).RN.get(nodeidx - 1, classidx - 1) / wnd_size);
                    double TN_ref_w = chainRefTput(layerIdx_0, classidx - 1, this.results.get(rLen - w).get(layerIdx_0));
                    if (TN_ref_w > GlobalConstants.FineTol) {
                        this.residt.set(aidx, this.residt.get(aidx) + this.results.get(rLen - w).get(layerIdx_0).QN.get(nodeidx - 1, classidx - 1) / TN_ref_w / wnd_size);
                    } else {
                        this.residt.set(aidx, this.residt.get(aidx) + this.results.get(rLen - w).get(layerIdx_0).WN.get(nodeidx - 1, classidx - 1) / wnd_size);
                    }
                    this.tput.set(aidx, this.tput.get(aidx) + this.results.get(rLen - w).get(layerIdx_0).TN.get(nodeidx - 1, classidx - 1) / wnd_size);
                }
            } else {
                this.servt.set(aidx, this.results.get(results.size()).get(layerIdx_0).RN.get(nodeidx - 1, classidx - 1));
                double TN_ref = chainRefTput(layerIdx_0, classidx - 1, this.results.get(results.size()).get(layerIdx_0));
                double QN_val = this.results.get(results.size()).get(layerIdx_0).QN.get(nodeidx - 1, classidx - 1);
                if (TN_ref > GlobalConstants.FineTol) {
                    this.residt.set(aidx, QN_val / TN_ref);
                } else {
                    this.residt.set(aidx, this.results.get(results.size()).get(layerIdx_0).WN.get(nodeidx - 1, classidx - 1));
                }
                this.tput.set(aidx, this.results.get(results.size()).get(layerIdx_0).TN.get(nodeidx - 1, classidx - 1));
            }

            // Force servt/residt to 0 for activities with Immediate service time
            // The NC solver may return non-zero RN due to numerical issues, but
            // Immediate activities have zero service time by definition
            if (lqn.hostdem.get(aidx) instanceof Immediate) {
                this.servt.set(aidx, 0);
                this.residt.set(aidx, 0);
            }

            // An activity think time is a delay in series with its host demand,
            // held at the activity's own task: the task keeps its thread for the
            // whole hostdem+thinktime interval, but the host processor is released
            // for it, so the processor layer is left alone and only the activity's
            // service and residence grow. entry_servt below is summed from residt,
            // so adding it to servt alone would never reach the entry. Placed after
            // the Immediate reset above, which would otherwise wipe it. Mirrors
            // MATLAB lqn_act_thinktime / updateMetricsDefault.
            double zt_act = actThinkTime(lqn, aidx);
            if (zt_act > 0) {
                this.servt.set(aidx, this.servt.get(aidx) + zt_act);
                this.residt.set(aidx, this.residt.get(aidx) + zt_act);
            }

            // Fix for async-only entry targets: use RN (response time per visit) for residt
            // The host layer closed model incorrectly splits residence time (WN) between
            // activities when an entry only receives async calls (no sync callers).
            // For async-only entries, use RN instead of WN since the async arrivals
            // don't share the closed chain's visit ratio - each async arrival gets
            // the full response time per visit.
            if (aidx >= lqn.ashift && aidx < lqn.ashift + lqn.nacts) {
                // This is an activity - find its bound entry
                for (int eidx = lqn.eshift; eidx < lqn.eshift + lqn.nentries; eidx++) {
                    if (lqn.graph.get(eidx, aidx) > 0) {
                        // Found bound entry - check if async-only
                        boolean hasSyncCallers = false;
                        boolean hasAsyncCallers = false;
                        for (int i = 0; i < lqn.nidx; i++) {
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
                            this.residt.set(aidx, this.servt.get(aidx)); // servt already has RN
                        }
                        break;
                    }
                }
            }

            // Recover from Inf/NaN: snap back to previous iteration's value
            if (it > 1) {
                if ((Double.isInfinite(this.servt.get(aidx)) || Double.isNaN(this.servt.get(aidx))) && !Double.isNaN(this.servt_prev.get(aidx))) {
                    this.servt.set(aidx, this.servt_prev.get(aidx));
                }
                if ((Double.isInfinite(this.residt.get(aidx)) || Double.isNaN(this.residt.get(aidx))) && !Double.isNaN(this.residt_prev.get(aidx))) {
                    this.residt.set(aidx, this.residt_prev.get(aidx));
                }
                if ((Double.isInfinite(this.tput.get(aidx)) || Double.isNaN(this.tput.get(aidx))) && !Double.isNaN(this.tput_prev.get(aidx))) {
                    this.tput.set(aidx, this.tput_prev.get(aidx));
                }
            }
            // Apply under-relaxation if enabled and not first iteration
            double omega = this.relax_omega;
            if (omega < 1.0 && it > 1) {
                if (!Double.isNaN(this.servt_prev.get(aidx))) {
                    this.servt.set(aidx, omega * this.servt.get(aidx) + (1 - omega) * this.servt_prev.get(aidx));
                }
                if (!Double.isNaN(this.residt_prev.get(aidx))) {
                    this.residt.set(aidx, omega * this.residt.get(aidx) + (1 - omega) * this.residt_prev.get(aidx));
                }
                if (!Double.isNaN(this.tput_prev.get(aidx))) {
                    this.tput.set(aidx, omega * this.tput.get(aidx) + (1 - omega) * this.tput_prev.get(aidx));
                }
            }
            // Store current values for next iteration
            this.servt_prev.set(aidx, this.servt.get(aidx));
            this.residt_prev.set(aidx, this.residt.get(aidx));
            this.tput_prev.set(aidx, this.tput.get(aidx));

            // Preserve Immediate type for activities that originally had Immediate service times
            // Safeguard against MVA numerical instability producing extreme values
            double max_servt = 1e10;
            if (lqn.hostdem.get(aidx) instanceof Immediate) {
                this.servtproc.put(aidx, Immediate.getInstance());
            } else if (this.servt.get(aidx) > 0 && this.servt.get(aidx) <= max_servt) {
                this.servtproc.put(aidx, Exp.fitMean(this.servt.get(aidx)));
            }
            this.tputproc.put(aidx, Exp.fitRate(this.tput.get(aidx)));
        }

        // Phase-2 support: split activity service times by phase
        if (this.hasPhase2) {
            // Reset phase-specific arrays
            this.servt_ph1 = new Matrix(1, lqn.nidx, lqn.nidx);
            this.servt_ph2 = new Matrix(1, lqn.nidx, lqn.nidx);

            // Split activity service times by phase
            for (int a = 0; a < lqn.nacts; a++) {
                int aidx = lqn.ashift + a;
                if (lqn.actphase.get(0, a) == 1) {
                    this.servt_ph1.set(aidx, this.servt.get(aidx));
                } else {
                    this.servt_ph2.set(aidx, this.servt.get(aidx));
                }
            }

            // Aggregate phase service times to entry level
            for (int e = 0; e < lqn.nentries; e++) {
                int eidx = lqn.eshift + e;
                List<Integer> acts = lqn.actsof.get(eidx);
                if (acts != null) {
                    for (int aidx : acts) {
                        int a = aidx - lqn.ashift;
                        if (a >= 0 && a < lqn.nacts) {
                            if (lqn.actphase.get(0, a) == 1) {
                                this.servt_ph1.set(eidx, this.servt_ph1.get(eidx) + this.servt_ph1.get(aidx));
                            } else {
                                this.servt_ph2.set(eidx, this.servt_ph2.get(eidx) + this.servt_ph2.get(aidx));
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
                    this.tput.set(aidx, 0);
                    for (int w = 0; w < wnd_size; w++) {
                        this.tput.set(aidx, this.tput.get(aidx) + this.results.get(rLen - w).get(this.idxhash.get(idx).intValue() - 1).TN.get(nodeidx - 1, classidx - 1) / wnd_size);
                    }
                } else {
                    this.tput.set(aidx, this.results.get(results.size()).get(this.idxhash.get(idx).intValue() - 1).TN.get(nodeidx - 1, classidx - 1));
                }
                this.tputproc.put(aidx, Exp.fitRate(this.tput.get(aidx)));
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
                    this.callservt.set(cidx, 0.0);
                } else {
                    int eidxLayer = this.idxhash.get(idx).intValue() - 1;
                    SolverResult layerRes = this.results.get(results.size()).get(eidxLayer);
                    double fcrWait = regionWait(eidxLayer, nodeidx - 1, classidx - 1, layerRes);
                    this.callservt.set(cidx, (layerRes.RN.get(nodeidx - 1, classidx - 1) + fcrWait) * this.lqn.callproc_mean.getOrDefault(cidx, 1.0));
                    // Normalise per chain-reference visit, as residt does above. WN divides
                    // by the class's own reference rate when the layer is open (an INF client
                    // task), which is per-ENTRY visit, and the entry rescaling below would
                    // then count the call once per entry.
                    double callTNref = chainRefTput(eidxLayer, classidx - 1, layerRes);
                    if (callTNref > GlobalConstants.FineTol) {
                        this.callresidt.set(cidx, layerRes.QN.get(nodeidx - 1, classidx - 1) / callTNref + fcrWait);
                    } else {
                        this.callresidt.set(cidx, layerRes.WN.get(nodeidx - 1, classidx - 1) + fcrWait);
                    }
                }
                // Growth rate capping removed - it prevents callservt from converging
                // to the correct value when initial values are near-zero (Immediate)
                // (matches MATLAB updateMetricsDefault.m)
                // Recover from Inf/NaN: snap back to previous iteration's value
                if (it > 1) {
                    if ((Double.isInfinite(this.callservt.get(cidx)) || Double.isNaN(this.callservt.get(cidx))) && !Double.isNaN(this.callservt_prev.get(cidx))) {
                        this.callservt.set(cidx, this.callservt_prev.get(cidx));
                    }
                    if ((Double.isInfinite(this.callresidt.get(cidx)) || Double.isNaN(this.callresidt.get(cidx))) && !Double.isNaN(this.callresidt_prev.get(cidx))) {
                        this.callresidt.set(cidx, this.callresidt_prev.get(cidx));
                    }
                }
                // Apply under-relaxation to call service times
                double omega = this.relax_omega;
                if (omega < 1.0 && it > 1 && !Double.isNaN(this.callservt_prev.get(cidx))) {
                    this.callservt.set(cidx, omega * this.callservt.get(cidx) + (1 - omega) * this.callservt_prev.get(cidx));
                }
                this.callservt_prev.set(cidx, this.callservt.get(cidx));
                this.callresidt_prev.set(cidx, this.callresidt.get(cidx));
            }
        }

        //then resolve the entry servt summming up these contributions
        Matrix out = new Matrix(1, lqn.nidx + lqn.ncalls, lqn.nidx + lqn.ncalls);
        Matrix entry_servt = new Matrix(this.servtmatrix.getNumRows(), 1, this.servtmatrix.getNumRows());
        Matrix.concatColumns(this.residt, this.callresidt, out);
        this.servtmatrix.mult(out.transpose(), entry_servt);

        for (int i = 0; i < lqn.eshift; i++) {
            entry_servt.set(i, 0, 0);
        }

        // Forwarding is represented by caller-side pseudo rendezvous calls
        // added by applyForwardingRendezvous, so FWD calls carry no blocking
        // here: callservt/callresidt of FWD calls remain zero and no chain
        // delay is charged into the SYNC calls.

        // Recompute entry_servt with forwarding-adjusted callresidt
        Matrix.concatColumns(this.residt, this.callresidt, out);
        this.servtmatrix.mult(out.transpose(), entry_servt);
        for (int i = 0; i < lqn.eshift; i++) {
            entry_servt.set(i, 0, 0);
        }

        // servtmatrix is a reachability matrix, so it charges every activity of every branch
        // of an AND-fork to the entry, i.e. it serialises branches that in fact run
        // concurrently. Replace that sum by the join completion time, for each join reachable
        // from the entry.
        Matrix jointExcess = this.updateJoinDelays();
        for (int eidx = lqn.eshift; eidx < lqn.eshift + lqn.nentries; eidx++) {
            double corrected = entry_servt.get(eidx, 0);
            for (int aidx = 1; aidx < jointExcess.getNumCols(); aidx++) {
                if (jointExcess.get(0, aidx) != 0 && this.servtmatrix.get(eidx, aidx) > 0) {
                    corrected += jointExcess.get(0, aidx);
                }
            }
            entry_servt.set(eidx, 0, Math.max(corrected, 0));
        }

        // A SetupTask's cold start is charged HERE, to the entry, and with the
        // probability that the thread was actually found powered down. It is not
        // host demand, so it does not belong to any activity's residence:
        // reporting it there put RespT(A2) at 1.29479 on lqn_setup against the
        // 0.333178 LDES measures, which is the bare demand. See setupCharge.
        for (int eidx = lqn.eshift; eidx < lqn.eshift + lqn.nentries; eidx++) {
            int tidx_su = (int) lqn.parent.get(0, eidx);
            entry_servt.set(eidx, 0, entry_servt.get(eidx, 0) + setupCharge(tidx_su));
        }

        // this block fixes the problem that ResidT is scaled so that the task as Vtask = 1,
        // but in call servt the entries need to have Ventry = 1
        for (int eidx = lqn.eshift; eidx < lqn.eshift + lqn.nentries; eidx++) {
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
            for (int i = 0; i < lqn.nidx; i++) {
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
                //entry_servt.set(eidx, ensemble[this.idxhash.get(hidx).intValue() - 1].getClassByIndex(tidxclass.get(0) - 1).getNumberOfJobs() / entry_tput - entry_servt_z);
                //System.out.println("value: "+entry_servt.get(eidx) * task_tput / FastMath.max(GlobalConstants.Zero, entry_tput));
                if (entry_tput > GlobalConstants.Zero) {
                    this.servt.set(eidx, entry_servt.get(eidx) * task_tput / entry_tput);
                    this.residt.set(eidx, entry_servt.get(eidx) * task_tput / entry_tput);
                } else {
                    this.servt.set(eidx, entry_servt.get(eidx));
                    this.residt.set(eidx, entry_servt.get(eidx));
                }
            } else {
                // For async-only targets, use entry_servt directly
                // No throughput ratio scaling needed since there are no closed classes
                this.servt.set(eidx, entry_servt.get(eidx));
                this.residt.set(eidx, entry_servt.get(eidx));
            }
        }

        // Phase-2 support: compute overtaking probability and apply correction
        // This must happen AFTER entry throughput is available (computed above)
        if (this.hasPhase2) {
            for (int e = 0; e < lqn.nentries; e++) {
                int eidx = lqn.eshift + e;
                int tidx = (int) lqn.parent.get(0, eidx);

                if (this.servt_ph2.get(eidx) > GlobalConstants.FineTol) {
                    // REF tasks and entries without sync callers see the full service time
                    boolean hasSyncCallersP2 = false;
                    for (int i = 0; i < lqn.nidx; i++) {
                        if (lqn.issynccaller.get(i, eidx) > 0) {
                            hasSyncCallersP2 = true;
                            break;
                        }
                    }
                    if ((tidx > 0 && lqn.isref.get(tidx) != 0) || !hasSyncCallersP2) {
                        this.residt.set(eidx, this.servt.get(eidx));
                        continue;
                    }
                    // Get entry throughput (use task throughput as approximation if entry not available)
                    double entry_tput;
                    if (this.tput.get(eidx) > GlobalConstants.FineTol) {
                        entry_tput = this.tput.get(eidx);
                    } else if (tidx > 0 && this.tput.get(tidx) > GlobalConstants.FineTol) {
                        entry_tput = this.tput.get(tidx);
                    } else {
                        entry_tput = 0;
                    }

                    // Compute overtaking probability now that throughput is available
                    if (entry_tput > GlobalConstants.FineTol) {
                        this.prOvertake.set(e, this.overtakeProb(eidx));
                    } else {
                        this.prOvertake.set(e, 0);
                    }

                    // Caller's response time = phase-1 + P(overtake) * phase-2
                    double overtake_delay = this.prOvertake.get(e) * this.servt_ph2.get(eidx);
                    this.residt.set(eidx, this.servt_ph1.get(eidx) + overtake_delay);
                }
            }
        }

        for (int r = 1; r < this.call_classes_updmap.getNumRows(); r++) {
            int cidx = (int) this.call_classes_updmap.get(r, 2);
            int eidx = (int) lqn.callpair.get(cidx, 1);
            if (this.call_classes_updmap.get(r, 3) > 1) {
                if (this.servt.get(eidx) > 0) {
                    this.servtproc.put(eidx, Exp.fitMean(this.servt.get(eidx)));
                }
            }
        }

        // determine call response time processes
        // this.callresidt starts from 0
        for (int r = 1; r < this.call_classes_updmap.getNumRows(); r++) {
            int cidx = (int) this.call_classes_updmap.get(r, 2);
            int eidx = (int) lqn.callpair.get(cidx, 1);
            if (this.call_classes_updmap.get(r, 3) > 1) {
                if (it == 1) {
                    // note that respt is per visit, so number of calls is 1
                    this.callservt.set(cidx, this.servt.get(eidx));
                    this.callservtproc.put(cidx, this.servtproc.get(eidx));
                } else {
                    // note that respt is per visit, so number of calls is 1
                    if (this.callservt.get(cidx) > 0) {
                        this.callservtproc.put(cidx, Exp.fitMean(this.callservt.get(cidx)));
                    }
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
        for (int t = 0; t < lqn.ntasks; t++) {
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
                for (int row = lqn.tshift; row < lqn.tshift + lqn.ntasks; row++) {
                    if (uniqueCallingIdx.contains(row)) {
                        callers.add(row);
                    }
                }

                Matrix caller_tput = new Matrix(lqn.ntasks, 1, lqn.ntasks);
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
                for (int i = 0; i < caller_tput.getNumRows(); i++) {
                    task_tput += caller_tput.get(i);
                }
                for (int i = 0; i < lqn.ntasks; i++) {
                    this.ptaskcallers.set(tidx, lqn.tshift + i, caller_tput.get(i) / FastMath.max(GlobalConstants.Zero, task_tput));
                }
            }
        }

        // determine ptaskcallers for direct callers to hosts
        for (int hidx = 0; hidx < lqn.nhosts; hidx++) {
            // Skip ignored hosts
            if (this.ignore.get(hidx) != 0) {
                continue;
            }
            Matrix caller_tput = new Matrix(1, lqn.ntasks, lqn.ntasks);
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
            for (int i = 0; i < lqn.ntasks; i++) {
                this.ptaskcallers.set(hidx, lqn.tshift + i, caller_tput.get(i) / FastMath.max(GlobalConstants.Zero, host_tput));
            }
        }

        // impute call probability using a DTMC random walk on the taskcaller graph
        // for matrix multiplication, let P to be the same size of that in Matlab
        Matrix P = new Matrix(this.ptaskcallers.getNumRows(), this.ptaskcallers.getNumCols(), this.ptaskcallers.getNumRows() * this.ptaskcallers.getNumCols());
        for (int i = 0; i < this.ptaskcallers.getNumRows(); i++) {
            for (int j = 0; j < this.ptaskcallers.getNumCols(); j++) {
                P.set(i, j, this.ptaskcallers.get(i, j));
            }
        }
        P = dtmc_makestochastic(P); // hold mass at reference stations when there

        this.ptaskcallers_step.put(1, P.copy()); // configure step = 1
        // Random walk block validated: indexing matches MATLAB updateMetricsDefault.m lines 319-351
        for (int h = 0; h < lqn.nhosts; h++) {
            int hidx = h;
            for (int i = 0; i < lqn.tasksof.get(hidx).size(); i++) {
                int tidx = lqn.tasksof.get(hidx).get(i);
                // initialize the probability mass on tidx
                Matrix x0 = new Matrix(1, P.length(), P.length());
                x0.set(hidx, 1);
                // start the walk backward to impute probability of indirect callers
                Matrix x = new Matrix(x0.getNumRows(), P.getNumCols(), x0.getNumRows() * P.getNumCols());
                x0.mult(P, x);
                for (int step = 2; step <= this.nlayers; step++) {
                    Matrix xret = Matrix.createLike(x0);
                    x.mult(P, xret);
                    x.setTo(xret);

                    for (int remidx = 0; remidx < this.ptaskcallers_step.get(step).getNumCols(); remidx++) {
                        this.ptaskcallers_step.get(step).set(tidx, remidx, x.get(remidx));
                        double scaled = this.ptaskcallers.get(hidx, tidx) * x.get(remidx);
                        this.ptaskcallers_step.get(step).set(hidx, remidx, scaled);
                    }

                    double sum = 0;
                    Matrix ref_nonzero = lqn.isref.find();
                    for (int k = 0; k < ref_nonzero.getNumElements(); k++) {
                        sum += x.get((int) ref_nonzero.get(k));
                    }
                    if (sum > 1.0 - this.options.tol) break;

                    for (int index = 0; index < this.ptaskcallers.getNumRows(); index++) {
                        double max = FastMath.max(this.ptaskcallers.get(index, tidx), x.get(0, index));
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
                this.servt.set(aidx, this.results.get(results.size()).get(this.idxhash.get(idx).intValue() - 1).RN.get(nodeidx - 1, classidx - 1));
                this.tput.set(aidx, this.results.get(results.size()).get(this.idxhash.get(idx).intValue() - 1).TN.get(nodeidx - 1, classidx - 1));
                // Safeguard against MVA numerical instability producing extreme values
                double max_servt_mb = 1e10;
                if (this.servt.get(aidx) > 0 && this.servt.get(aidx) <= max_servt_mb) {
                    this.servtproc.put(aidx, Exp.fitMean(this.servt.get(aidx)));
                }

                // Compute residt from QN/TN_ref (matching MATLAB updateMetricsMomentBased)
                int layerIdx_0 = this.idxhash.get(idx).intValue() - 1;
                double TN_ref = chainRefTput(layerIdx_0, classidx - 1, this.results.get(results.size()).get(layerIdx_0));
                if (TN_ref > GlobalConstants.FineTol) {
                    this.residt.set(aidx, this.results.get(results.size()).get(layerIdx_0).QN.get(nodeidx - 1, classidx - 1) / TN_ref);
                } else {
                    this.residt.set(aidx, this.results.get(results.size()).get(layerIdx_0).WN.get(nodeidx - 1, classidx - 1));
                }

                // An activity think time is in series with the host demand
                double zt_act = actThinkTime(lqn, aidx);
                if (zt_act > 0) {
                    this.servt.set(aidx, this.servt.get(aidx) + zt_act);
                    this.residt.set(aidx, this.residt.get(aidx) + zt_act);
                    this.servtproc.put(aidx, Exp.fitMean(this.servt.get(aidx)));
                }

                // async-only targets carry no visit-ratio scaling (matching updateMetricsDefault)
                if (aidx >= lqn.ashift && aidx < lqn.ashift + lqn.nacts) {
                    for (int eidx = lqn.eshift; eidx < lqn.eshift + lqn.nentries; eidx++) {
                        if (lqn.graph.get(eidx, aidx) > 0) {
                            boolean hasSyncCallers = false;
                            boolean hasAsyncCallers = false;
                            for (int i = 0; i < lqn.nidx; i++) {
                                if (lqn.issynccaller.get(i, eidx) > 0) {
                                    hasSyncCallers = true;
                                }
                                if (lqn.isasynccaller.get(i, eidx) > 0) {
                                    hasAsyncCallers = true;
                                }
                            }
                            if (hasAsyncCallers && !hasSyncCallers) {
                                this.residt.set(aidx, this.servt.get(aidx));
                            }
                            break;
                        }
                    }
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
                        this.callservt.set(cidx, 0.0);
                        this.callresidt.set(cidx, 0.0);
                    } else {
                        int eidxLayer = this.idxhash.get(idx).intValue() - 1;
                        SolverResult layerRes = this.results.get(results.size()).get(eidxLayer);
                        double fcrWait = regionWait(eidxLayer, nodeidx - 1, classidx - 1, layerRes);
                        // Include call multiplicity in callservt (matching MATLAB)
                        this.callservt.set(cidx, (layerRes.RN.get(nodeidx - 1, classidx - 1) + fcrWait) * this.lqn.callproc_mean.getOrDefault(cidx, 1.0));
                        // callresidt uses WN which already includes visit multiplicity
                        this.callresidt.set(cidx, layerRes.WN.get(nodeidx - 1, classidx - 1) + fcrWait);
                    }
                    // Growth rate capping removed - it prevents callservt from converging
                    // to the correct value when initial values are near-zero (Immediate)
                    // (matches MATLAB updateMetricsDefault.m)
                    this.callservt_prev.set(cidx, this.callservt.get(cidx));
                    this.callresidt_prev.set(cidx, this.callresidt.get(cidx));
                }
            }

            // Then resolve the entry servt summing up these contributions; the terms are
            // residence times (Vtask=1), which the task/entry tput ratio below rescales to Ventry=1
            Matrix out = new Matrix(1, lqn.nidx + lqn.ncalls, lqn.nidx + lqn.ncalls);
            Matrix entry_servt = new Matrix(this.servtmatrix.getNumRows(), 1, this.servtmatrix.getNumRows());
            Matrix.concatColumns(this.residt, this.callresidt, out);

            // Solve the system: entry_servt = (I - servtmatrix)^(-1) * [residt; callresidt]
            Matrix identity = Matrix.eye(lqn.nidx + lqn.ncalls);
            Matrix system = identity.sub(this.servtmatrix);
            Matrix rhs = out.transpose();
            entry_servt = system.inv().mult(rhs);

            // Clear entries up to eshift
            for (int i = 0; i < lqn.eshift; i++) {
                entry_servt.set(i, 0, 0);
            }

            // NO forwarding propagation here. lqnFwdRendezvous has already reconnected
            // every forwarding chain reachable from a synchronous call to the client
            // that issued the rendezvous (Franks 1999, Sec. 3.3.1), so the forwarded
            // service is in the caller's chain before this runs; adding it again
            // inflated the caller by exactly the forwarded entry's mean. An
            // asynchronous call into a chain is left untouched there by design -- a
            // send-no-reply does not block -- so it must not accumulate the forwarded
            // service either. See BUGS.md BUG-91.

            // Update servt for entries
            for (int i = lqn.eshift; i < lqn.eshift + lqn.nentries; i++) {
                this.servt.set(i, entry_servt.get(i, 0));
            }

            // Clear activities after ashift
            for (int i = lqn.ashift; i < entry_servt.getNumRows(); i++) {
                entry_servt.set(i, 0, 0);
            }

            // Compute entry-level residt using servtmatrix and forwarding-augmented callresidt
            Matrix residt_out = new Matrix(1, lqn.nidx + lqn.ncalls, lqn.nidx + lqn.ncalls);
            Matrix.concatColumns(this.residt, this.callresidt, residt_out);
            Matrix entry_residt = new Matrix(this.servtmatrix.getNumRows(), 1, this.servtmatrix.getNumRows());
            this.servtmatrix.mult(residt_out.transpose(), entry_residt);
            for (int i = 0; i < lqn.eshift; i++) {
                entry_residt.set(i, 0, 0);
            }

            // Scale entry residt (and servt) by task/entry throughput ratio
            for (int eidx = lqn.eshift; eidx < lqn.eshift + lqn.nentries; eidx++) {
                int tidx = (int) lqn.parent.get(0, eidx);
                if (tidx < 0) continue;
                int hidx = (int) lqn.parent.get(0, tidx);
                if (this.ignore.get(tidx) != 0 || this.ignore.get(hidx) != 0) continue;

                boolean hasSyncCallers = false;
                for (int ii = 0; ii < lqn.nidx; ii++) {
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
                    if (entry_tput > GlobalConstants.Zero) {
                        this.servt.set(eidx, entry_servt.get(eidx, 0) * task_tput / entry_tput);
                        this.residt.set(eidx, entry_residt.get(eidx, 0) * task_tput / entry_tput);
                    } else {
                        this.residt.set(eidx, entry_residt.get(eidx, 0));
                    }
                } else {
                    this.residt.set(eidx, entry_residt.get(eidx, 0));
                }
            }

            // Update servtproc for entries based on call classes
            for (int r = 1; r < this.call_classes_updmap.getNumRows(); r++) {
                int cidx = (int) this.call_classes_updmap.get(r, 2);    // call
                int eidx = (int) lqn.callpair.get(cidx, 1);             // entry index (1-indexed)
                if (this.call_classes_updmap.get(r, 3) > 1) {
                    if (this.servt.get(eidx) > 0) {
                        this.servtproc.put(eidx, Exp.fitMean(this.servt.get(eidx)));
                    }
                }
            }

            // Determine call response times processes
            for (int r = 1; r < this.call_classes_updmap.getNumRows(); r++) {
                int cidx = (int) this.call_classes_updmap.get(r, 2);    // call
                int eidx = (int) lqn.callpair.get(cidx, 1);             // entry index (1-indexed)
                if (this.call_classes_updmap.get(r, 3) > 1) {
                    if (it == 1) {
                        // Note that respt is per visit, so number of calls is 1
                        this.callservt.set(cidx, this.servt.get(eidx));
                        this.callservtproc.put(cidx, this.servtproc.get(eidx));
                    } else {
                        // Note that respt is per visit, so number of calls is 1
                        if (this.callservt.get(cidx) > 0) {
                            this.callservtproc.put(cidx, Exp.fitMean(this.callservt.get(cidx)));
                        }
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

                this.tput.set(aidx, this.results.get(results.size()).get(this.idxhash.get(idx).intValue() - 1).TN.get(nodeidx - 1, classidx - 1));

                // Compute residt from QN/TN_ref (matching updateMetricsDefault)
                int layerIdx_0 = this.idxhash.get(idx).intValue() - 1;
                double TN_ref = chainRefTput(layerIdx_0, classidx - 1, this.results.get(results.size()).get(layerIdx_0));
                if (TN_ref > GlobalConstants.FineTol) {
                    this.residt.set(aidx, this.results.get(results.size()).get(layerIdx_0).QN.get(nodeidx - 1, classidx - 1) / TN_ref);
                } else {
                    this.residt.set(aidx, this.results.get(results.size()).get(layerIdx_0).WN.get(nodeidx - 1, classidx - 1));
                }

                int submodelidx = this.idxhash.get(idx).intValue();
                if (!repo.containsKey(submodelidx)) {
                    // Try SolverFluid first for actual response time CDFs with
                    // higher-moment information; fall back to layer solver's
                    // exponential approximation if Fluid fails on the layer model
                    try {
                        // Fluid iter_max counts ODE horizon extensions, each retained in full,
                        // so it must stay at the Fluid default as in MATLAB updateMetricsMomentBased
                        SolverOptions fluidOpts = SolverFluid.defaultOptions();
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
                        this.callresidt.set(cidx, this.results.get(results.size()).get(this.idxhash.get(idx).intValue() - 1).WN.get(nodeidx - 1, classidx - 1));
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
            for (int i = 0; i < lqn.nentries; i++) {
                int eidx = lqn.eshift + i;

                // Find contributing indices (where matrix(eidx,:) > 0)
                List<Integer> convolidx = new ArrayList<Integer>();
                for (int j = 0; j < matrix.getNumCols(); j++) {
                    if (matrix.get(eidx, j) > 0) {
                        // Skip hosts, tasks and entries: only the activity block convolves
                        if (j >= lqn.ashift) {
                            convolidx.add(j);
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

                    // An activity think time is in series with the host demand,
                    // so its raw moments convolve with the measured ones before
                    // the APH fit, as the reference's lqn_act_thinktime block does
                    if (fitidx < lqn.nidx && actThinkTime(lqn, fitidx) > 0
                            && lqn.actthink != null && lqn.actthink.get(fitidx) != null) {
                        Distribution ztd = lqn.actthink.get(fitidx);
                        double t1 = ztd.getMean();
                        double sig2 = ztd.getSCV() * t1 * t1;
                        double t2 = sig2 + t1 * t1;
                        double t3 = ztd.getSkewness() * Math.pow(sig2, 1.5) + 3 * t1 * t2 - 2 * t1 * t1 * t1;
                        m3 = m3 + 3 * m2 * t1 + 3 * m1 * t2 + t3;
                        m2 = m2 + 2 * m1 * t1 + t2;
                        m1 = m1 + t1;
                    }

                    if (m1 > GlobalConstants.CoarseTol) {
                        // Fit APH from raw moments
                        APH fitdist = APH.fitRawMoments(m1, m2, m3);
                        Matrix alpha = fitdist.getInitProb();
                        Matrix T = (Matrix) fitdist.getParam(3).getValue();

                        // For call indices, multiply repetitions by mean number of calls
                        // servtmatrix has 1.0 for calls, but we need callproc_mean repetitions
                        double repetitions = matrix.get(eidx, fitidx);
                        if (fitidx >= lqn.nidx) {
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

                            jline.util.Pair<Matrix, Matrix> branchResult = aph_simplify(
                                alpha, T, zeroAlpha, zeroT,
                                fractionalPart, 1.0 - fractionalPart, 3);
                            paramList.add(new Pair<Matrix, Matrix>(branchResult.getFirst(), branchResult.getSecond()));
                        } else {
                            // Only fractional part
                            APH zeroDist = APH.fitMeanAndSCV(GlobalConstants.FineTol, 0.99);
                            Matrix zeroAlpha = zeroDist.getInitProb();
                            Matrix zeroT = (Matrix) zeroDist.getParam(3).getValue();

                            jline.util.Pair<Matrix, Matrix> branchResult = aph_simplify(
                                alpha, T, zeroAlpha, zeroT,
                                fractionalPart, 1.0 - fractionalPart, 3);
                            paramList.add(new Pair<Matrix, Matrix>(branchResult.getFirst(), branchResult.getSecond()));
                        }

                        // Update servtproc and callservtproc based on fitidx
                        if (fitidx < lqn.nidx) {
                            this.servtproc.put(fitidx, Exp.fitMean(m1));
                            this.servt.set(fitidx, m1);
                        } else {
                            this.callservtproc.put(fitidx - lqn.nidx, Exp.fitMean(m1));
                            this.callservt.set(fitidx - lqn.nidx, m1);
                        }
                    }
                }

                // Convolve all contributions
                if (paramList.isEmpty()) {
                    this.servt.set(eidx, 0);
                } else {
                    // Build jline.util.Pair list for aph_convseq
                    List<jline.util.Pair<Matrix, Matrix>> jlineParamList = new ArrayList<jline.util.Pair<Matrix, Matrix>>();
                    for (Pair<Matrix, Matrix> p : paramList) {
                        jlineParamList.add(new jline.util.Pair<Matrix, Matrix>(p.getLeft(), p.getRight()));
                    }

                    jline.util.Pair<Matrix, Matrix> convResult = aph_convseq(jlineParamList);
                    APH entryDist = new APH(convResult.getFirst(), convResult.getSecond());

                    // Store results. entryIndex is the ENTRY-LOCAL index, 0-based,
                    // the same one entryproc is keyed by.
                    int entryIndex = eidx - lqn.eshift;
                    this.entryproc.put(entryIndex, entryDist);
                    this.servt.set(eidx, entryDist.getMean());
                    this.servtproc.put(eidx, Exp.fitMean(this.servt.get(eidx)));

                    // The WHOLE law, as MATLAB stores it: an (n x 2) table whose
                    // columns are [F(t), t]. Keeping F(0) alone -- zero by
                    // construction -- left getCdfRespT with nothing to report, and
                    // the 1-based slot it went to was out of range for entry 0.
                    this.entrycdfrespt.put(entryIndex, entryDist.evalCDFMatrix());
                }
            }

            // NO forwarding propagation here, for the reason given at the entry_servt
            // assembly above: lqnFwdRendezvous has already charged the forwarded
            // service to the caller. See BUGS.md BUG-91.

            // Compute entry-level residt using servtmatrix and activity residt
            Matrix residt_out_pc = new Matrix(1, lqn.nidx + lqn.ncalls, lqn.nidx + lqn.ncalls);
            Matrix.concatColumns(this.residt, this.callresidt, residt_out_pc);
            Matrix entry_residt_pc = new Matrix(this.servtmatrix.getNumRows(), 1, this.servtmatrix.getNumRows());
            this.servtmatrix.mult(residt_out_pc.transpose(), entry_residt_pc);
            for (int i = 0; i < lqn.eshift; i++) {
                entry_residt_pc.set(i, 0, 0);
            }

            // Scale entry residt by task/entry throughput ratio
            for (int eidx = lqn.eshift; eidx < lqn.eshift + lqn.nentries; eidx++) {
                int tidx = (int) lqn.parent.get(0, eidx);
                if (tidx < 0) continue;
                int hidx = (int) lqn.parent.get(0, tidx);
                if (this.ignore.get(tidx) != 0 || this.ignore.get(hidx) != 0) continue;

                boolean hasSyncCallers = false;
                for (int ii = 0; ii < lqn.nidx; ii++) {
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
                    if (entry_tput > GlobalConstants.Zero) {
                        this.residt.set(eidx, entry_residt_pc.get(eidx, 0) * task_tput / entry_tput);
                    } else {
                        this.residt.set(eidx, entry_residt_pc.get(eidx, 0));
                    }
                } else {
                    this.residt.set(eidx, entry_residt_pc.get(eidx, 0));
                }
            }

            // Determine call response times processes (final loop)
            for (int r = 1; r < this.call_classes_updmap.getNumRows(); r++) {
                int cidx = (int) this.call_classes_updmap.get(r, 2);
                int eidx = (int) lqn.callpair.get(cidx, 1);
                if (this.call_classes_updmap.get(r, 3) > 1) {
                    if (it == 1) {
                        this.callservt.set(cidx, this.servt.get(eidx));
                        this.callservtproc.put(cidx, Exp.fitMean(this.servt.get(eidx)));
                    }
                }
            }

            // This pass IS the moment3 answer, and it is TERMINAL. Its entry laws
            // are convolutions of the activities' own response distributions; the
            // pre-convergence branch instead reads QN/TN_ref, a residence per
            // REFERENCE cycle, which the entry assembly then treats as a
            // per-entry-visit time. The two disagree by the entry's visit ratio
            // whenever it is not 1, so letting the iteration fall back to that
            // branch after this one has run DISCARDS the moment-based laws and
            // reports the other quantity. See BUGS.md BUG-97.
            this.momentPassDone = true;
        }
    }

    // ========================================================================
    // Interlock: static analysis (built once at init)
    // ========================================================================

    /**
     * Build the interlock path table and locate the common parents.
     * Interlocking arises when requests issued by one client reach a common
     * lower-level server along two or more independent paths, so that
     * arrivals a layer decomposition treats as independent are in fact
     * correlated. Franks (1999), Ch. 4:
     *   Phase A: the path table path(a,b) of Sec. 4.2, the calls to entry b
     *            caused by one invocation of entry a, with a unit diagonal;
     *            a second table restricts the count to the phase-1 flow
     *   Phase B: the common-parent finder of Fig. 4.2, retaining only the
     *            entries at which the flow genuinely splits
     *   Phase C: the source tasks and the source count n_s of Eq. (4.7)
     * The phase-aware tables, the branch-point test and the source count are
     * refinements beyond the published algorithm, which assumes one path
     * table and counts source tasks directly.
     */
    @SuppressWarnings("unchecked")
    private void initInterlock() {
        LayeredNetworkStruct lqn = this.lqn;

        // Phase A: path table of Sec. 4.2
        int nentries = lqn.nentries;
        double[][] il_all = new double[nentries][nentries];
        double[][] il_ph1 = new double[nentries][nentries];

        for (int e = 0; e < nentries; e++) {
            int eidx = lqn.eshift + e;
            boolean[] visited = new boolean[nentries];
            traceInterlockPaths(lqn, eidx, e, 1.0, 1.0, visited, il_all, il_ph1, 0);
        }

        this.il_table_all = il_all;
        this.il_table_ph1 = il_ph1;

        // Phase B+C: Find common entries and sources per server entity
        int arraySize = lqn.tshift + lqn.ntasks;
        this.il_common_entries = new List[arraySize];
        this.il_source_tasks_all = new List[arraySize];
        this.il_source_tasks_ph2 = new List[arraySize];
        this.il_num_sources = new double[arraySize];

        // Process task servers
        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            if (lqn.isref.get(tidx) != 0 || lqn.sched.get(tidx) == SchedStrategy.INF) {
                continue;
            }
            int[][] result = findInterlockForServer(lqn, tidx, il_all, il_ph1);
            this.il_common_entries[tidx] = toIntList(result[0]);
            this.il_source_tasks_all[tidx] = toIntList(result[1]);
            this.il_source_tasks_ph2[tidx] = toIntList(result[2]);
            this.il_num_sources[tidx] = result[3].length > 0 ? result[3][0] : 0;
        }

        // Process host servers
        for (int h = 0; h < lqn.nhosts; h++) {
            int hidx = h;
            if (lqn.sched.get(hidx) == SchedStrategy.INF) {
                continue;
            }
            int[][] result = findInterlockForServer(lqn, hidx, il_all, il_ph1);
            this.il_common_entries[hidx] = toIntList(result[0]);
            this.il_source_tasks_all[hidx] = toIntList(result[1]);
            this.il_source_tasks_ph2[hidx] = toIntList(result[2]);
            this.il_num_sources[hidx] = result[3].length > 0 ? result[3][0] : 0;
        }
    }

    private static List<Integer> toIntList(int[] arr) {
        List<Integer> list = new ArrayList<Integer>();
        for (int v : arr) {
            list.add(v);
        }
        return list;
    }

    /**
     * Phase A: Recursive path tracing for interlock reachability.
     */
    private void traceInterlockPaths(LayeredNetworkStruct lqn, int eidx, int root_e,
                                      double prob_all, double prob_ph1, boolean[] visited,
                                      double[][] il_all, double[][] il_ph1, int depth) {
        int e = eidx - lqn.eshift;
        if (e < 0 || e >= lqn.nentries) return;
        if (visited[e]) return;
        visited[e] = true;

        // Record reachability from root to this entry
        il_all[root_e][e] += prob_all;
        il_ph1[root_e][e] += prob_ph1;

        // Follow synchronous calls from activities of this entry
        List<Integer> acts = lqn.actsof.get(eidx);
        if (acts != null) {
            for (int aidx : acts) {
                if (aidx < lqn.ashift || aidx >= lqn.ashift + lqn.nacts) continue;
                int a = aidx - lqn.ashift;

                // Pruning: at non-root entries (depth > 0), skip phase-2+ activities
                if (depth > 0 && lqn.actphase != null && lqn.actphase.get(0, a) > 1) {
                    continue;
                }

                boolean is_ph1 = true;
                if (lqn.actphase != null && lqn.actphase.get(0, a) > 1) {
                    is_ph1 = false;
                }

                // Follow calls from this activity
                List<Integer> calls_from_act = lqn.callsof.get(aidx);
                if (calls_from_act != null) {
                    for (int cidx : calls_from_act) {
                        if (cidx < 0 || cidx >= lqn.ncalls) continue;
                        if (lqn.calltype.get(cidx) != CallType.SYNC) continue;
                        double call_mean = lqn.callproc_mean.getOrDefault(cidx, 0.0);
                        if (call_mean <= 0) continue;
                        int dst_eidx = (int) lqn.callpair.get(cidx, 1);
                        int dst_e = dst_eidx - lqn.eshift;
                        if (dst_e < 0 || dst_e >= lqn.nentries) continue;

                        double next_all = prob_all * call_mean;
                        double next_ph1 = is_ph1 ? prob_ph1 * call_mean : 0;

                        traceInterlockPaths(lqn, dst_eidx, root_e, next_all, next_ph1, visited, il_all, il_ph1, depth + 1);
                    }
                }
            }
        }

        visited[e] = false;
    }

    /**
     * Phase B+C: Find interlock for a single server.
     * Returns int[4][] where [0]=commonEntries, [1]=srcAll, [2]=srcPh2, [3]={numSources}.
     */
    private int[][] findInterlockForServer(LayeredNetworkStruct lqn, int serverIdx,
                                            double[][] il_all, double[][] il_ph1) {
        int[][] empty = new int[][] { new int[0], new int[0], new int[0], new int[0] };

        // Get server entry numbers
        List<Integer> serverEntryNums = getServerEntryNums(lqn, serverIdx);
        if (serverEntryNums.isEmpty()) return empty;

        // Get client tasks
        List<Integer> clientTasks = getClientTasks(lqn, serverIdx);
        if (clientTasks.size() < 1) return empty;

        // Get client entries that reach the server
        List<int[]> clientEntryPairs = new ArrayList<int[]>(); // [taskIdx, entryNum]
        for (int ct : clientTasks) {
            List<Integer> entries = lqn.entriesof.get(ct);
            if (entries == null) continue;
            for (int ce : entries) {
                int ce_num = ce - lqn.eshift;
                if (ce_num < 0 || ce_num >= lqn.nentries) continue;
                for (int se_num : serverEntryNums) {
                    if (il_all[ce_num][se_num] > 0) {
                        clientEntryPairs.add(new int[] { ct, ce_num });
                        break;
                    }
                }
            }
        }

        if (clientEntryPairs.size() < 2) return empty;

        // Find common parent entries (branch points)
        Set<Integer> commonEntriesSet = new LinkedHashSet<Integer>();
        int nPairs = clientEntryPairs.size();
        for (int i = 0; i < nPairs; i++) {
            for (int j = i + 1; j < nPairs; j++) {
                if (clientEntryPairs.get(i)[0] == clientEntryPairs.get(j)[0]) continue; // Same task
                int entryA_num = clientEntryPairs.get(i)[1];
                int entryC_num = clientEntryPairs.get(j)[1];

                // Search all tasks for common parents
                for (int t = 0; t < lqn.ntasks; t++) {
                    int tidx = lqn.tshift + t;
                    List<Integer> entries_of_task = lqn.entriesof.get(tidx);
                    if (entries_of_task == null) continue;
                    for (int ex : entries_of_task) {
                        for (int ey : entries_of_task) {
                            int ex_num = ex - lqn.eshift;
                            int ey_num = ey - lqn.eshift;
                            if (ex_num < 0 || ey_num < 0 || ex_num >= lqn.nentries || ey_num >= lqn.nentries) continue;
                            if (il_all[ex_num][entryA_num] > 0 && il_all[ey_num][entryC_num] > 0) {
                                if (isBranchPointCheck(lqn, ex, entryA_num + lqn.eshift, ey, entryC_num + lqn.eshift, il_all)) {
                                    commonEntriesSet.add(ex);
                                }
                            }
                        }
                    }
                }
            }
        }

        if (commonEntriesSet.isEmpty()) return empty;

        List<Integer> commonEntries = new ArrayList<Integer>(commonEntriesSet);

        // Phase C: Find source tasks
        Set<Integer> interlockedTasksSet = new LinkedHashSet<Integer>();
        for (int ce_eidx : commonEntries) {
            List<Integer> itasks = findInterlockedTasks(lqn, ce_eidx, serverIdx, il_all);
            interlockedTasksSet.addAll(itasks);
        }
        List<Integer> interlockedTasks = new ArrayList<Integer>(interlockedTasksSet);

        // All source tasks = tasks owning common entries
        Set<Integer> allSrcSet = new LinkedHashSet<Integer>();
        for (int ce_eidx : commonEntries) {
            int owner_tidx = (int) lqn.parent.get(0, ce_eidx);
            allSrcSet.add(owner_tidx);
        }

        // Remove interlocked tasks from allSrcTasks
        allSrcSet.removeAll(interlockedTasksSet);

        // Ph2 sources: interlocked tasks with phase-2 activities reaching server
        Set<Integer> ph2SrcSet = new LinkedHashSet<Integer>();
        for (int it : interlockedTasks) {
            List<Integer> itEntries = lqn.entriesof.get(it);
            if (itEntries == null) continue;
            for (int ie : itEntries) {
                if (hasPhase2Activities(lqn, ie)) {
                    int ie_num = ie - lqn.eshift;
                    if (ie_num >= 0 && ie_num < lqn.nentries) {
                        boolean found = false;
                        for (int se_num : serverEntryNums) {
                            if (il_all[ie_num][se_num] - il_ph1[ie_num][se_num] > 0) {
                                found = true;
                                break;
                            }
                        }
                        if (found) {
                            ph2SrcSet.add(it);
                        }
                    }
                }
            }
        }

        // Add external sources (tasks calling into interlocked paths from outside)
        for (int it : interlockedTasks) {
            List<Integer> itEntries = lqn.entriesof.get(it);
            if (itEntries == null) continue;
            for (int ie : itEntries) {
                for (int ci = 0; ci < lqn.iscaller.getNumRows(); ci++) {
                    if (lqn.iscaller.get(ci, ie) != 0) {
                        if (ci >= lqn.tshift && ci < lqn.tshift + lqn.ntasks) {
                            if (!interlockedTasksSet.contains(ci)) {
                                allSrcSet.add(ci);
                            }
                        }
                    }
                }
            }
        }

        // Count total source multiplicity
        double nsrc = 0;
        for (int st : allSrcSet) {
            nsrc += lqn.mult.get(st);
        }

        List<Integer> allSrcTasks = new ArrayList<Integer>(allSrcSet);
        List<Integer> ph2SrcTasks = new ArrayList<Integer>(ph2SrcSet);

        return new int[][] {
            listToIntArray(commonEntries),
            listToIntArray(allSrcTasks),
            listToIntArray(ph2SrcTasks),
            new int[] { (int) nsrc }
        };
    }

    private static int[] listToIntArray(List<Integer> list) {
        int[] arr = new int[list.size()];
        for (int i = 0; i < list.size(); i++) {
            arr[i] = list.get(i);
        }
        return arr;
    }

    /** Get server entry numbers (1-based entry numbers, not absolute indices). */
    private List<Integer> getServerEntryNums(LayeredNetworkStruct lqn, int serverIdx) {
        List<Integer> nums = new ArrayList<Integer>();
        if (serverIdx < lqn.nhosts) {
            List<Integer> tasks = lqn.tasksof.get(serverIdx);
            if (tasks != null) {
                for (int tidx : tasks) {
                    List<Integer> entries = lqn.entriesof.get(tidx);
                    if (entries != null) {
                        for (int se : entries) {
                            nums.add(se - lqn.eshift);
                        }
                    }
                }
            }
        } else {
            List<Integer> entries = lqn.entriesof.get(serverIdx);
            if (entries != null) {
                for (int se : entries) {
                    nums.add(se - lqn.eshift);
                }
            }
        }
        return nums;
    }

    /** Get client tasks for a server. */
    private List<Integer> getClientTasks(LayeredNetworkStruct lqn, int serverIdx) {
        if (serverIdx < lqn.nhosts) {
            List<Integer> tasks = lqn.tasksof.get(serverIdx);
            return tasks != null ? tasks : new ArrayList<Integer>();
        } else {
            Set<Integer> clientSet = new LinkedHashSet<Integer>();
            List<Integer> server_entries = lqn.entriesof.get(serverIdx);
            if (server_entries != null) {
                for (int se : server_entries) {
                    for (int ci = 0; ci < lqn.iscaller.getNumRows(); ci++) {
                        if (lqn.iscaller.get(ci, se) != 0) {
                            if (ci >= lqn.tshift && ci < lqn.tshift + lqn.ntasks) {
                                clientSet.add(ci);
                            }
                        }
                    }
                }
            }
            return new ArrayList<Integer>(clientSet);
        }
    }

    /** Branch point check. */
    private boolean isBranchPointCheck(LayeredNetworkStruct lqn, int srcX_eidx, int entryA_eidx,
                                        int srcY_eidx, int entryB_eidx, double[][] il_all) {
        int taskA = (int) lqn.parent.get(0, entryA_eidx);
        int taskB = (int) lqn.parent.get(0, entryB_eidx);
        int taskX = (int) lqn.parent.get(0, srcX_eidx);

        // Multiserver client: if X, A, B same task => not branch point
        if (taskX == taskA && taskX == taskB) return false;

        // Quick check: direct call
        if (srcX_eidx == entryA_eidx || srcY_eidx == entryB_eidx) return true;

        int entryA_num = entryA_eidx - lqn.eshift;
        int entryB_num = entryB_eidx - lqn.eshift;

        // Check downstream calls diverge to different tasks
        List<Integer> dstTasks_X = getCallDstTasks(lqn, srcX_eidx, entryA_num, il_all);
        List<Integer> dstTasks_Y = getCallDstTasks(lqn, srcY_eidx, entryB_num, il_all);

        for (int dx : dstTasks_X) {
            for (int dy : dstTasks_Y) {
                if (dx != dy) return true;
            }
        }
        return false;
    }

    /** Get destination tasks of sync calls from an entry reaching a target. */
    private List<Integer> getCallDstTasks(LayeredNetworkStruct lqn, int src_eidx, int target_e_num, double[][] il_all) {
        Set<Integer> dstTasks = new LinkedHashSet<Integer>();
        List<Integer> acts = lqn.actsof.get(src_eidx);
        if (acts == null) return new ArrayList<Integer>(dstTasks);
        for (int aidx : acts) {
            if (aidx < lqn.ashift || aidx >= lqn.ashift + lqn.nacts) continue;
            List<Integer> calls = lqn.callsof.get(aidx);
            if (calls == null) continue;
            for (int cidx : calls) {
                if (cidx < 0 || cidx >= lqn.ncalls) continue;
                if (lqn.calltype.get(cidx) != CallType.SYNC) continue;
                int dst_eidx = (int) lqn.callpair.get(cidx, 1);
                int dst_e = dst_eidx - lqn.eshift;
                if (dst_e >= 0 && dst_e < lqn.nentries && il_all[dst_e][target_e_num] > 0) {
                    dstTasks.add((int) lqn.parent.get(0, dst_eidx));
                }
            }
        }
        return new ArrayList<Integer>(dstTasks);
    }

    /** Get interlocked tasks on paths from an entry to a server. */
    private List<Integer> findInterlockedTasks(LayeredNetworkStruct lqn, int src_eidx, int serverIdx, double[][] il_all) {
        boolean[] visited = new boolean[lqn.nentries];
        List<Integer> itasks = new ArrayList<Integer>();
        traceToServerRec(lqn, src_eidx, serverIdx, il_all, visited, itasks, true);
        return itasks;
    }

    private void traceToServerRec(LayeredNetworkStruct lqn, int eidx, int serverIdx,
                                   double[][] il_all, boolean[] visited, List<Integer> itasks, boolean isHead) {
        int e = eidx - lqn.eshift;
        if (e < 0 || e >= lqn.nentries || visited[e]) return;

        int ownerTask = (int) lqn.parent.get(0, eidx);

        // Check if we reached the server
        if (ownerTask == serverIdx) return;
        if (serverIdx < lqn.nhosts && (int) lqn.parent.get(0, ownerTask) == serverIdx) return;

        visited[e] = true;

        // Follow synchronous calls from ALL phases
        List<Integer> acts = lqn.actsof.get(eidx);
        boolean found = false;
        if (acts != null) {
            for (int aidx : acts) {
                if (aidx < lqn.ashift || aidx >= lqn.ashift + lqn.nacts) continue;
                List<Integer> calls = lqn.callsof.get(aidx);
                if (calls == null) continue;
                for (int cidx : calls) {
                    if (cidx < 0 || cidx >= lqn.ncalls || lqn.calltype.get(cidx) != CallType.SYNC) continue;
                    int dst_eidx = (int) lqn.callpair.get(cidx, 1);
                    int dst_task = (int) lqn.parent.get(0, dst_eidx);

                    // Check if destination reaches server
                    boolean reachesServer = false;
                    if (dst_task == serverIdx) {
                        reachesServer = true;
                    } else if (serverIdx < lqn.nhosts && (int) lqn.parent.get(0, dst_task) == serverIdx) {
                        reachesServer = true;
                    } else {
                        int dst_e = dst_eidx - lqn.eshift;
                        List<Integer> serverEntryNums = getServerEntryNums(lqn, serverIdx);
                        for (int se_num : serverEntryNums) {
                            if (dst_e >= 0 && dst_e < lqn.nentries && il_all[dst_e][se_num] > 0) {
                                reachesServer = true;
                                break;
                            }
                        }
                    }

                    if (reachesServer) {
                        traceToServerRec(lqn, dst_eidx, serverIdx, il_all, visited, itasks, false);
                        found = true;
                    }
                }
            }
        }

        if (found && !isHead) {
            if (!itasks.contains(ownerTask)) {
                itasks.add(ownerTask);
            }
        }

        visited[e] = false;
    }

    /** Check if entry has phase-2 activities. */
    private boolean hasPhase2Activities(LayeredNetworkStruct lqn, int eidx) {
        if (lqn.actphase == null) return false;
        List<Integer> acts = lqn.actsof.get(eidx);
        if (acts == null) return false;
        for (int aidx : acts) {
            int a = aidx - lqn.ashift;
            if (a >= 0 && a < lqn.nacts && lqn.actphase.get(0, a) > 1) {
                return true;
            }
        }
        return false;
    }

    // ========================================================================
    // Interlock: flow computation (called each iteration)
    // ========================================================================

    /**
     * Compute interlock probability for a (client, server) pair.
     */
    /** True when the layer solver applies Eq. (4.7) inside its own MVA. */
    private boolean layerTakesInterlock(int e) {
        // Only the MVA layer solver reads options.config.interlock, and only a layer whose sole
        // queueing stations are the host's own tasks can take a matrix built for that host:
        // under flat layering one layer holds every server, so the correction stays on the
        // residence times there.
        if (e < 0 || this.solvers == null || e >= this.solvers.length || this.solvers[e] == null) {
            return false;
        }
        if (!(this.solvers[e] instanceof jline.solvers.mva.SolverMVA)) {
            return false;
        }
        String layering = this.options.config.layering;
        if (layering != null && ("flat".equalsIgnoreCase(layering) || "squashed".equalsIgnoreCase(layering))) {
            return false;
        }
        // A layer whose MVA path has no interlock term would be moved to another algorithm by
        // the matrix alone: exact multiserver MVA would become AMVA, the linearizer would become
        // the load-dependent forward step. That swap is worth far more than the correction it
        // carries, and on a layer sitting near a bifurcation it turns the LN iteration into a
        // limit cycle. Such a layer keeps the residt scaling below instead.
        return jline.solvers.mva.analyzers.Solver_mva_analyzer.mvaCarriesInterlock(
                this.ensemble[e].getStruct(), this.solvers[e].options);
    }

    /**
     * Class-level interlock matrix of one host layer. IL(r,s) is the share of the class-s queue
     * that a class-r arrival must not see at the host. The matrix is CLASS-indexed, not
     * chain-indexed, so that a later refreshChains cannot leave it stale; the layer solver
     * aggregates it to chains against the struct it is about to solve. Two classes are
     * interlocked only if BOTH their tasks are, which is the 0/1 relation ir_mkj of Eq. (5);
     * the diagonal stays zero, since a request always sees its own class in full. The entry is
     * the Eq. (5) product Pr(IL_ms)*IR_ms*IR_mr, asymmetric in (r,s) because Pr(IL) is taken
     * from the QUEUED class s, so that the layer's ILw(r,s) = 1-IL(r,s) is the lower-level
     * adjustment rate r_lower of Li and Franks (2015).
     */
    private Matrix buildLayerInterlock(int e, List<Integer> host_tasks, double[] task_prIL, double[] task_PrIL) {
        int nclasses = this.ensemble[e].getClasses().size();
        double[] class_prIL = new double[nclasses];   // IR
        double[] class_PrIL = new double[nclasses];   // Pr(IL)
        for (int r = 0; r < nclasses; r++) {
            int tidx = clientTaskOfClass(e, r);
            if (tidx < 0) continue;
            int ti = host_tasks.indexOf(tidx);
            if (ti >= 0) {
                class_prIL[r] = task_prIL[ti];
                class_PrIL[r] = task_PrIL[ti];
            }
        }
        Matrix IL = new Matrix(nclasses, nclasses);
        boolean any = false;
        for (int r = 0; r < nclasses; r++) {
            if (class_prIL[r] <= GlobalConstants.FineTol) continue;
            for (int sIl = 0; sIl < nclasses; sIl++) {
                if (sIl == r || class_prIL[sIl] <= GlobalConstants.FineTol) continue;
                IL.set(r, sIl, class_PrIL[sIl] * class_prIL[sIl] * class_prIL[r]);
                any = true;
            }
        }
        if (any) {
            any = false;
            for (int r = 0; r < nclasses && !any; r++)
                for (int sIl = 0; sIl < nclasses && !any; sIl++)
                    if (IL.get(r, sIl) > GlobalConstants.FineTol) any = true;
        }
        return any ? IL : null; // nothing interlocked, keep the layer on the plain MVA path
    }

    /** Task that a layer class belongs to, -1 when the class names no task. */
    private int clientTaskOfClass(int e, int c) {
        Integer[] attr = this.ensemble[e].getClasses().get(c).getAttribute();
        if (attr == null || attr.length < 2 || attr[1] == null) {
            return -1;
        }
        int tidx = -1;
        if (attr[0] == LayeredNetworkElement.TASK) {
            tidx = attr[1];
        } else if (attr[0] == LayeredNetworkElement.ENTRY || attr[0] == LayeredNetworkElement.ACTIVITY) {
            tidx = (int) lqn.parent.get(0, attr[1]);
        } else if (attr[0] == LayeredNetworkElement.CALL) {
            tidx = (int) lqn.parent.get(0, (int) lqn.callpair.get(attr[1], 0));
        }
        if (tidx < lqn.tshift || tidx >= lqn.tshift + lqn.ntasks) {
            return -1;
        }
        return tidx;
    }

    /**
     * Interlock probability for one (client, server) pair, as {IR, Pr(IL)} of Li and Franks,
     * "An improved interlocking correction for decomposition of layered queueing networks",
     * CCECE 2015, Eqs. (3) and (4). isProcessorHost selects the m' rule of lqns
     * Interlock::ilrate_pril_flow: at a PROCESSOR the common-source population is doubled above
     * 3 customers and squared at or below it, which is what turns m = 4 into the pril = 1/8 its
     * trace reports. The two factors are multiplied into the Eq. (5) rate by buildLayerInterlock,
     * so neither carries the source count on its own -- that lives in m'. This replaces the
     * superseded (n_s-1)/n_s discount of Franks (1999), Eq. (4.7).
     */
    private double[] computeInterlockProb(LayeredNetworkStruct lqn, int client_tidx, int server_idx,
                                          boolean isProcessorHost) {
        List<Integer> commonEntries = this.il_common_entries[server_idx];
        double numSources = this.il_num_sources[server_idx];
        List<Integer> allSrcTasks = this.il_source_tasks_all[server_idx];
        List<Integer> ph2SrcTasks = this.il_source_tasks_ph2[server_idx];

        double[] none = new double[]{0, 0};
        if (numSources == 0 || commonEntries == null || commonEntries.isEmpty()) return none;

        // Get client entries
        List<Integer> client_entries = lqn.entriesof.get(client_tidx);
        if (client_entries == null) return none;

        // Interlocked flow lambda^IL of Eq. (4), and alongside it the flow weighted by 1/m',
        // which gives Pr(IL) of Eq. (3).
        double sum_flow = 0;
        double sum_pril = 0;
        for (int ce_eidx : commonEntries) {
            int srcTask = (int) lqn.parent.get(0, ce_eidx);
            int ce_num = ce_eidx - lqn.eshift;
            // population of this common source, in customer copies
            double m_src = lqn.mult.get(srcTask);
            if (Double.isNaN(m_src) || Double.isInfinite(m_src) || m_src < 1) m_src = 1;
            double m_eff;
            if (isProcessorHost) {
                m_eff = (m_src > 3) ? (m_src + m_src) : (m_src * m_src);
            } else {
                m_eff = m_src;
            }

            for (int dstA_eidx : client_entries) {
                int dstA_num = dstA_eidx - lqn.eshift;
                if (dstA_num < 0 || dstA_num >= lqn.nentries) continue;
                if (this.il_table_all[ce_num][dstA_num] <= 0) continue;

                // Get source entry throughput
                double ce_tput = getEntryTput(lqn, ce_eidx, srcTask);
                if (ce_tput <= GlobalConstants.FineTol) continue;

                // Deferred flow is scored separately from the phase-1 flow
                boolean hasP2 = hasPhase2Activities(lqn, ce_eidx);

                if (!hasP2 && allSrcTasks.contains(srcTask)) {
                    double contrib = ce_tput * this.il_table_all[ce_num][dstA_num];
                    sum_flow += contrib;
                    sum_pril += contrib / m_eff;
                } else if (hasP2 && allSrcTasks.contains(srcTask)) {
                    double contrib = ce_tput * this.il_table_ph1[ce_num][dstA_num];
                    sum_flow += contrib;
                    sum_pril += contrib / m_eff;
                }

                double ph2 = this.il_table_all[ce_num][dstA_num] - this.il_table_ph1[ce_num][dstA_num];
                if (ph2 > 0 && ph2SrcTasks.contains(srcTask)) {
                    double contrib = ce_tput * ph2;
                    sum_flow += contrib;
                    sum_pril += contrib / m_eff;
                }
            }
        }

        // Get client throughput
        double client_tput = getTaskTput(lqn, client_tidx);
        if (client_tput <= GlobalConstants.FineTol) return none;

        double IR = Math.min(sum_flow, client_tput) / client_tput;
        IR = Math.min(1.0, Math.max(0.0, IR));
        double prIL = (sum_flow <= GlobalConstants.FineTol) ? 0 : (sum_pril / sum_flow);
        prIL = Math.min(1.0, Math.max(0.0, prIL));
        return new double[]{IR, prIL};
    }

    /** Helper: get entry throughput. */
    private double getEntryTput(LayeredNetworkStruct lqn, int eidx, int taskIdx) {
        double tputVal = this.tput.get(eidx);
        if (tputVal <= GlobalConstants.FineTol) {
            List<Integer> acts = lqn.actsof.get(eidx);
            if (acts != null && !acts.isEmpty()) {
                tputVal = this.tput.get(acts.get(0));
            }
        }
        if (tputVal <= GlobalConstants.FineTol) {
            tputVal = this.tput.get(taskIdx);
        }
        return tputVal;
    }

    /** Helper: get task throughput. */
    private double getTaskTput(LayeredNetworkStruct lqn, int tidx) {
        double tputVal = this.tput.get(tidx);
        if (tputVal <= GlobalConstants.FineTol) {
            List<Integer> entries = lqn.entriesof.get(tidx);
            if (entries != null) {
                for (int eidx : entries) {
                    double et = this.tput.get(eidx);
                    if (et <= GlobalConstants.FineTol) {
                        List<Integer> acts = lqn.actsof.get(eidx);
                        if (acts != null && !acts.isEmpty()) {
                            et = this.tput.get(acts.get(0));
                        }
                    }
                    tputVal += et;
                }
            }
        }
        return tputVal;
    }

    public void updatePopulations(int it) {
        // Eq. (4.7) removes one source in n_s from the queue length inside
        // MVA; the equivalent correction is applied here to the residence
        // times returned by the layer, leaving service and utilization
        // untouched. See Franks (1999), Ch. 4.
        LayeredNetworkStruct lqn = this.lqn;

        if (!this.options.config.interlocking || this.il_common_entries == null) {
            return;
        }

        // Save originals for proportional entry_servt update
        Matrix callresidt_orig = this.callresidt.copy();
        Matrix residt_orig = this.residt.copy();
        boolean adjusted = false;

        // Pass 1: For each sync call, check if destination server has interlock
        for (int cidx = 0; cidx < lqn.ncalls; cidx++) {
            if (lqn.calltype.get(cidx) != CallType.SYNC) continue;

            int dst_eidx = (int) lqn.callpair.get(cidx, 1);
            int server_tidx = (int) lqn.parent.get(0, dst_eidx);

            // Find the server entity with interlock data
            int server_for_il = -1;
            if (server_tidx < this.il_common_entries.length && this.il_common_entries[server_tidx] != null && !this.il_common_entries[server_tidx].isEmpty()) {
                server_for_il = server_tidx;
            } else {
                // Check host server
                if (server_tidx >= lqn.tshift) {
                    int host_idx = (int) lqn.parent.get(0, server_tidx);
                    if (host_idx >= 0 && host_idx < this.il_common_entries.length && this.il_common_entries[host_idx] != null && !this.il_common_entries[host_idx].isEmpty()) {
                        server_for_il = host_idx;
                    }
                }
            }
            if (server_for_il < 0) continue;

            // Get client task (activity -> task via parent)
            int src_aidx = (int) lqn.callpair.get(cidx, 0);
            int client_tidx = (int) lqn.parent.get(0, src_aidx);

            // Interlock probability for this client and server. This path serves a TASK, not a
            // processor, so the m' rule of Li/lqns leaves the source population alone; the
            // product IR*Pr(IL) reproduces the scalar this branch used before.
            double[] ilp = computeInterlockProb(lqn, client_tidx, server_for_il, false);
            double prIL = ilp[0] * ilp[1];
            if (prIL <= GlobalConstants.FineTol) continue;

            // Compute waiting time reduction
            double S = this.servt.get(dst_eidx);  // service time at destination entry
            double call_mean = lqn.callproc_mean.getOrDefault(cidx, 0.0);
            if (call_mean <= 0 || this.callservt.get(cidx) <= 0) continue;

            double RN = this.callservt.get(cidx) / call_mean;  // response time per visit
            double W = Math.max(0, RN - S);  // waiting time per visit

            if (W > GlobalConstants.FineTol) {
                double RN_adj = S + (1 - prIL) * W;
                double scale = RN_adj / RN;
                this.callservt.set(cidx, this.callservt.get(cidx) * scale);
                this.callresidt.set(cidx, this.callresidt.get(cidx) * scale);
                if (this.callservt.get(cidx) > 0) {
                    this.callservtproc.put(cidx, Exp.fitMean(this.callservt.get(cidx)));
                }
                adjusted = true;
            }
        }

        // Pass 2: Host-level interlock -- reduce processor queueing in residt.
        // Every layer starts the pass without a matrix, so a host that stops being
        // interlocked does not keep the previous iteration's correction alive.
        for (int e = 0; e < this.nlayers; e++) {
            if (this.solvers[e] != null && this.solvers[e].options != null && this.solvers[e].options.config != null) {
                this.solvers[e].options.config.interlock = null;
            }
        }
        for (int h = 0; h < lqn.nhosts; h++) {
            int hidx = h;
            if (this.il_common_entries[hidx] == null || this.il_common_entries[hidx].isEmpty()) continue;

            // Compute prIL and processor utilization for each task on this host
            List<Integer> host_tasks = lqn.tasksof.get(hidx);
            if (host_tasks == null) continue;
            double[] task_prIL = new double[host_tasks.size()];   // IR, Eq. (4)
            double[] task_PrIL = new double[host_tasks.size()];   // Pr(IL), Eq. (3)
            double[] task_util = new double[host_tasks.size()];
            for (int ti = 0; ti < host_tasks.size(); ti++) {
                int tidx = host_tasks.get(ti);
                // The host of a task layer is a PROCESSOR, which is what selects the m' rule.
                double[] ilp = computeInterlockProb(lqn, tidx, hidx, true);
                task_prIL[ti] = ilp[0];
                task_PrIL[ti] = ilp[1];
                // Compute task's processor utilization
                List<Integer> entries = lqn.entriesof.get(tidx);
                if (entries != null) {
                    for (int eidx : entries) {
                        List<Integer> acts = lqn.actsof.get(eidx);
                        if (acts != null) {
                            for (int aidx : acts) {
                                double tputVal = this.tput.get(aidx);
                                double hostdem = lqn.hostdem_mean.containsKey(aidx) ? lqn.hostdem_mean.get(aidx) : 0;
                                task_util[ti] += tputVal * hostdem;
                            }
                        }
                    }
                }
            }

            double U_total = 0;
            double U_interlocked = 0;
            for (int ti = 0; ti < host_tasks.size(); ti++) {
                U_total += task_util[ti];
                if (task_prIL[ti] > GlobalConstants.FineTol) {
                    U_interlocked += task_util[ti];
                }
            }
            if (U_total <= GlobalConstants.FineTol || U_interlocked <= GlobalConstants.FineTol) continue;
            double il_fraction = U_interlocked / U_total;

            // When the layer solver carries Eq. (4.7) inside its own MVA, the interlock goes
            // to the layer as a class-level matrix and the residence times are left untouched.
            // Scaling them here as well would remove the same waiting twice, and would still
            // leave the layer's own THROUGHPUT uncorrected, which is what breaks flow balance
            // across a call: the reported task rate then comes from a cycle time the correction
            // has already shortened elsewhere.
            int layerOfHost = (hidx < idxhash.size() && !Double.isNaN(idxhash.get(hidx)))
                    ? idxhash.get(hidx).intValue() - 1 : -1;
            if (layerOfHost >= 0 && layerTakesInterlock(layerOfHost)) {
                Matrix ILmat = buildLayerInterlock(layerOfHost, host_tasks, task_prIL, task_PrIL);
                this.solvers[layerOfHost].options.config.interlock = ILmat;
                continue;
            }

            for (int ti = 0; ti < host_tasks.size(); ti++) {
                if (task_prIL[ti] <= GlobalConstants.FineTol) continue;
                int tidx = host_tasks.get(ti);
                // Weight by the share of host utilization that is interlocked. The rate is the
                // SAME Eq. (5) product IR*Pr(IL) that pass 1 applies to a call and that
                // buildLayerInterlock puts in the layer matrix -- IR alone is a flow SHARE, ~1
                // whenever a layer has a single common source, and using it here removed the
                // whole processor queueing rather than the interlocked part of it, which broke
                // flow balance across a call.
                double effective_prIL = task_prIL[ti] * task_PrIL[ti] * il_fraction;
                List<Integer> entries = lqn.entriesof.get(tidx);
                if (entries != null) {
                    for (int eidx : entries) {
                        List<Integer> acts = lqn.actsof.get(eidx);
                        if (acts != null) {
                            for (int aidx : acts) {
                                double D = lqn.hostdem_mean.containsKey(aidx) ? lqn.hostdem_mean.get(aidx) : 0;
                                if (D > 0 && this.residt.get(aidx) > D + GlobalConstants.FineTol) {
                                    double W_proc = this.residt.get(aidx) - D;
                                    this.residt.set(aidx, D + (1 - effective_prIL) * W_proc);
                                    adjusted = true;
                                }
                            }
                        }
                    }
                }
            }
        }

        if (!adjusted) return;

        // Recompute entry service times from adjusted callresidt/residt
        // Use proportional scaling to preserve visit ratio adjustments
        Matrix out_old = new Matrix(1, lqn.nidx + lqn.ncalls, lqn.nidx + lqn.ncalls);
        Matrix.concatColumns(residt_orig, callresidt_orig, out_old);
        Matrix entry_servt_old = new Matrix(this.servtmatrix.getNumRows(), 1, this.servtmatrix.getNumRows());
        this.servtmatrix.mult(out_old.transpose(), entry_servt_old);

        Matrix out_new = new Matrix(1, lqn.nidx + lqn.ncalls, lqn.nidx + lqn.ncalls);
        Matrix.concatColumns(this.residt, this.callresidt, out_new);
        Matrix entry_servt_new = new Matrix(this.servtmatrix.getNumRows(), 1, this.servtmatrix.getNumRows());
        this.servtmatrix.mult(out_new.transpose(), entry_servt_new);

        // The entry servt is rescaled only when it was itself assembled from these
        // residence times, which is the default path. After the moment3 pass it is
        // the MEAN OF AN APH CONVOLUTION of the activities' own response laws, and a
        // ratio of residence-time sums is not a correction to it: applying it
        // multiplies the entry law by the entry's visit ratio and reports a service
        // time BELOW that of the single activity the entry contains. The residence
        // times keep their correction either way. See BUGS.md BUG-97.
        boolean momentLaws = "moment3".equals(this.lnmethod) && this.momentPassDone;

        for (int eidx = lqn.eshift; eidx < lqn.eshift + lqn.nentries; eidx++) {
            if (entry_servt_old.get(eidx, 0) > GlobalConstants.FineTol) {
                double ratio = entry_servt_new.get(eidx, 0) / entry_servt_old.get(eidx, 0);
                if (!momentLaws) {
                    this.servt.set(eidx, this.servt.get(eidx) * ratio);
                    if (this.servt.get(eidx) > 0) {
                        this.servtproc.put(eidx, Exp.fitMean(this.servt.get(eidx)));
                    }
                }
                this.residt.set(eidx, this.residt.get(eidx) * ratio);
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
                                        Matrix TN_copy = this.results.get(this.results.size()).get(hostValue.intValue() - 1).TN;
                                        double Xtot = TN_copy.sumRows(stationIdxOf(hostValue.intValue() - 1, hostInt) - 1);
                                        if (Xtot > 0) {
                                            // MATLAB: hm_tput = sum(self.results{end,self.idxhash(host)}.TN(self.ensemble{self.idxhash(host)}.attribute.serverIdx,classidxto))
                                            double hm_tput = TN_copy.get(stationIdxOf(hostValue.intValue() - 1, hostInt) - 1, (int) classidxto - 1);
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
                                        Matrix TN_copy = this.results.get(this.results.size()).get(tidxCallerValue.intValue() - 1).TN;
                                        double Xtot = TN_copy.sumRows(stationIdxOf(tidxCallerValue.intValue() - 1, tidxCallerInt) - 1);
                                        if (Xtot > 0) {
                                            Map<Integer, Integer[]> call_map = callerNetwork.getAttribute().getCalls();
                                            // Find ALL entry classes for this call - MATLAB: calls(find(calls(:,4) == eidx), 1)
                                            // MATLAB column 4 is 1-indexed, so in Java 0-indexed array it's index 3
                                            // Must collect ALL matching classes (not just the first) for multi-entry callers
                                            List<Integer> matchingEidxClasses = new ArrayList<>();
                                            for (Map.Entry<Integer, Integer[]> entry : call_map.entrySet()) {
                                                Integer[] callInfo = entry.getValue();
                                                if (callInfo.length > 3 && callInfo[3] == (int) eidx) {
                                                    // MATLAB: calls(..., 1) - column 1 = index 0
                                                    matchingEidxClasses.add(callInfo[0]);
                                                }
                                            }
                                            {
                                                // MATLAB: entry_tput = sum(self.results{...}.TN(serverIdx, eidxclass))
                                                // Sum throughput across ALL matching entry classes; an entry with
                                                // no matching call classes gets entry_tput = 0 (MATLAB sum of empty),
                                                // so its routing probability is driven to zero rather than left at
                                                // the initial uniform value
                                                double entry_tput = 0;
                                                for (int eidxclass : matchingEidxClasses) {
                                                    entry_tput += TN_copy.get(stationIdxOf(tidxCallerValue.intValue() - 1, tidxCallerInt) - 1, eidxclass - 1);
                                                }
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
                    // Log the error for debugging
                    System.err.println("[WARN] updateRoutingProbabilities exception at row " + r + ": " + e.getMessage());
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

    /** True when the ensemble is a single flat layer. */
    private boolean isFlatLayering() {
        String layering = this.options != null && this.options.config != null
                ? this.options.config.layering : null;
        return layering != null
                && (layering.equalsIgnoreCase("flat") || layering.equalsIgnoreCase("squashed"));
    }

    /**
     * Mean cold start one request of task TIDX pays, 0 when it declares none.
     *
     * <p>A SetupTask powers a thread down when it goes idle and pays a setup before
     * it can serve again. The thread is released at a reply and starts a delay-off
     * countdown D of mean d; it powers off only if D expires before the next
     * request arrives, and a request arriving first cancels the countdown and pays
     * nothing. With the idle interval I seen by one thread and exponential D,</p>
     *
     * <pre>  p = P(D &lt; I) = E[I] / (E[I] + d),   and the charge is  p * s.</pre>
     *
     * <p>E[I] comes from the current iterate. Admission takes an ACTIVE idle thread
     * before it wakes a sleeping one, so the pool that actually cycles is only as
     * large as the load needs: with offered load b = X*S = rho*mult threads, about
     * max(1,b) stay hot, each seeing arrivals at rate X/max(1,b) and busy S per
     * arrival, so E[I] = (max(1,b) - b) / X. At mult = 1 this is (1-rho)/X and is
     * EXACT given p, returning p = a/(a+d) for the one-customer model that
     * SolverLDESLayeredSetupTest checks against. Above one thread it is an
     * approximation, the exact answer for c servers with setup being
     * matrix-analytic (Gandhi, Harchol-Balter and Adan, Performance Evaluation
     * 67(11), 2010). Twin of MATLAB lqn_setup_charge.m.</p>
     */
    private double setupCharge(int tidx) {
        if (tidx < 0 || lqn.hassetup == null || lqn.hassetup.getNumCols() <= tidx
                || lqn.hassetup.get(0, tidx) == 0) {
            return 0;
        }
        double s = setupMeanOf(lqn.setuptime, tidx);
        double d = setupMeanOf(lqn.delayofftime, tidx);
        if (!(s > GlobalConstants.FineTol) || !(d > GlobalConstants.FineTol)) {
            return 0;
        }
        double mult = lqn.mult.get(0, tidx);
        if (Double.isInfinite(mult) || Double.isNaN(mult) || mult <= 0) {
            return 0; // an infinite-server task holds no thread to power down
        }
        if (this.tput == null || this.util == null
                || this.tput.length() <= tidx || this.util.length() <= tidx) {
            return s; // nothing has arrived yet, so the thread is down when the first does
        }
        double X = this.tput.get(tidx);
        if (!(X > GlobalConstants.FineTol) || Double.isNaN(X) || Double.isInfinite(X)) {
            return s;
        }
        double rho = this.util.get(tidx);
        if (Double.isNaN(rho) || Double.isInfinite(rho) || rho < 0) {
            rho = 0;
        }
        rho = FastMath.min(rho, 1 - GlobalConstants.FineTol);
        double b = rho * mult;                    // offered load, in threads
        double EI = (FastMath.max(1, b) - b) / X; // idle interval of a hot thread
        return s * EI / (EI + d);
    }

    /** Mean of a setup or delay-off process of task TIDX, 0 when it declares none. */
    private double setupMeanOf(Map<Integer, Distribution> procs, int tidx) {
        if (procs == null) {
            return 0;
        }
        Distribution p = procs.get(tidx);
        if (p == null) {
            return 0;
        }
        double m = p.getMean();
        return (Double.isNaN(m) || Double.isInfinite(m)) ? 0 : m;
    }

    /**
     * Mean time per request that a thread of task TIDX stays busy after replying,
     * i.e. servt - residt averaged over the entries by their share of the requests.
     */
    private double phase2Tail(int tidx) {
        if (!this.hasPhase2 || this.servt_ph2 == null || this.residt == null) {
            return 0;
        }
        double tail = 0;
        double wtot = 0;
        List<Integer> entries = lqn.entriesof.get(tidx);
        if (entries == null) {
            return 0;
        }
        for (int eidx : entries) {
            if (this.servt_ph2.get(eidx) <= GlobalConstants.FineTol) {
                continue;
            }
            double w = this.tput.get(eidx);
            if (w <= GlobalConstants.FineTol) {
                w = 1;
            }
            tail += w * FastMath.max(0, this.servt.get(eidx) - this.residt.get(eidx));
            wtot += w;
        }
        if (wtot > GlobalConstants.FineTol) {
            tail = tail / wtot;
        }
        return tail;
    }

    /**
     * Total exogenous rate into the entries of task TIDX, zero unless the arrival is the
     * only way in. Matches the openArrivalOnly predicate the layer builder drops the open
     * class on: with a caller or a forwarding source the stream rides a class of its own
     * and closing the caller chain on the rate as well would load the layer twice.
     */
    private double openArrivalRateOf(int tidx) {
        if (this.lqn.isref.get(tidx) != 0 || this.lqn.arrival == null) {
            return 0;
        }
        List<Integer> entries = this.lqn.entriesof.get(tidx);
        if (entries == null) {
            return 0;
        }
        for (int eidx : entries) {
            for (int row = 0; row < this.lqn.ntasks; row++) {
                int tidx_row = row + this.lqn.tshift;
                if (this.lqn.issynccaller.get(tidx_row, eidx) != 0
                        || this.lqn.isasynccaller.get(tidx_row, eidx) != 0) {
                    return 0;
                }
            }
            for (int cidx = 0; cidx < this.lqn.ncalls; cidx++) {
                if (this.lqn.calltype.get(cidx) == CallType.FWD
                        && (int) this.lqn.callpair.get(cidx, 1) == eidx) {
                    return 0;
                }
            }
        }
        double rate = 0;
        for (int eidx : entries) {
            Distribution arv = this.lqn.arrival.get(eidx);
            if (arv == null) {
                continue;
            }
            double m = arv.getMean();
            if (!Double.isInfinite(m) && !Double.isNaN(m) && m > GlobalConstants.FineTol) {
                rate += 1 / m;
            }
        }
        return rate;
    }

    public void updateThinkTimes(int it) {
        // Under 'srvn.ph' a caller reaches the server once per invocation, so the
        // station rate is not the task's invocation rate -- see updateThinkTimesPH
        if (isPHEncoding()) {
            if (this.thinkt == null) {
                this.thinkt = new Matrix(1, this.lqn.ntasks + this.lqn.tshift,
                        this.lqn.ntasks + this.lqn.tshift);
            }
            this.thinktproc = new HashMap<Integer, Distribution>();
            updateThinkTimesPH(it);
            return;
        }
        //task16:updateThinkTimes function to be written
        if (this.lqn.iscaller.getNumCols() > 0) { // ignore models without callers
            Matrix torder = new Matrix(1, lqn.ntasks, lqn.ntasks);
            for (int i = 0; i < lqn.ntasks; i++)
                torder.set(i, i);

            this.thinkt = new Matrix(1, this.lqn.ntasks + this.lqn.tshift, this.lqn.ntasks + this.lqn.ntasks - 1);
            this.thinktproc = new HashMap<Integer, Distribution>();

            // solve all task models
            for (int t = 0; t < this.lqn.ntasks; t++) {
                int tidx = this.lqn.tshift + t;
                // only a REFERENCE task's think time separates one request from
                // the next; on a served task it is not a per-request delay and
                // charging it throttles the task -- see refThinkTime
                double tidx_thinktime = refThinkTime(this.lqn, tidx);
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
                    int layer_t = this.idxhash.get(tidx).intValue() - 1;
                    int station_t = stationIdxOf(layer_t, tidx) - 1;
                    double rawTN = matrixExtracted.sumRows(station_t);
                    this.tput.set(tidx, this.lqn.repl.get(0, tidx) * rawTN);

                    // obtain total self.utilization of task t
                    Matrix UmatrixExtracted = this.results.get(this.results.size()).get(layer_t).UN;
                    this.util.set(tidx, UmatrixExtracted.sumRows(station_t));

                    if (this.lqn.sched.get(tidx) == SchedStrategy.INF) { // first we consider the update where t is an infinite server
                        // key LQN think-time update: in LINE an infinite-server self.utilization is dimensionally a mean number of jobs
                        this.thinkt.set(tidx, FastMath.max(GlobalConstants.Zero, (njobs - this.util.get(tidx)) / this.tput.get(tidx) - tidx_thinktime));
                    } else { // otherwise we consider the case where t is a regular queueing station (other than an infinite server)
                        // key LQN think-time update: in LINE self.utilization is scaled to [0,1] for all queueing stations regardless of server count
                        this.thinkt.set(tidx, FastMath.max(GlobalConstants.Zero, njobs * FastMath.abs(1 - this.util.get(tidx)) / this.tput.get(tidx) - tidx_thinktime));
                    }
                    // Phase-2 tail: the entry replies after phase 1, so the station
                    // serves the caller for residt, but the thread stays busy for
                    // servt. What is left is busy time, not think time -- without this
                    // the task cycles slower than its callers call it, and throughput
                    // is not conserved across the call.
                    if (!isFlatLayering()) {
                        this.thinkt.set(tidx, FastMath.max(GlobalConstants.Zero,
                                this.thinkt.get(tidx) - phase2Tail(tidx)));
                    }
                    // The cold start goes the OTHER way from the phase-2 tail. A
                    // caller class cycles as delay plus station service, and the
                    // station serves only the host demand: the charge is on the
                    // entry, not on any activity's demand, so the station never sees
                    // it and the delay has to carry it. Without this the callee layer
                    // cycled at 0.529412 against the 0.5 its callers drive on
                    // lqn_setup. Zero for every task without a setup.
                    this.thinkt.set(tidx, FastMath.max(GlobalConstants.Zero,
                            this.thinkt.get(tidx) + setupCharge(tidx)));
                    // Recover from Inf/NaN: snap back to previous iteration's value
                    if (it > 1 && !Double.isNaN(this.thinkt_prev.get(tidx))) {
                        if (Double.isInfinite(this.thinkt.get(tidx)) || Double.isNaN(this.thinkt.get(tidx))) {
                            this.thinkt.set(tidx, this.thinkt_prev.get(tidx));
                        }
                    }
                    // Apply under-relaxation to think time if enabled
                    double omega = this.relax_omega;
                    if (omega < 1.0 && it > 1 && !Double.isNaN(this.thinkt_prev.get(tidx))) {
                        double rawT = this.thinkt.get(tidx);
                        double prevT = this.thinkt_prev.get(tidx);
                        // If recovering from crash (prev much larger than raw), snap to raw
                        if (prevT > 10 * rawT && rawT > GlobalConstants.FineTol) {
                            this.thinkt_prev.set(tidx, rawT); // reset prev to allow recovery
                        }
                        this.thinkt.set(tidx, omega * this.thinkt.get(tidx) + (1 - omega) * this.thinkt_prev.get(tidx));
                    }
                    this.thinkt_prev.set(tidx, this.thinkt.get(tidx));
                    Exp exponential = Exp.fitMean(this.thinkt.get(tidx) + tidx_thinktime);
                    this.thinktproc.put(tidx, exponential);
                } else { // ref task, forwarding target or open-arrival target (no task layer)
                    // A task reached only by an entry arrival has no caller, so no task
                    // layer, so nothing above would ever set its surrogate delay and the
                    // thread pool on its host layer would cycle against an Immediate one.
                    // The stream is known, so the cycle is closed on it directly: the same
                    // construction as a forwarding target below. See lqn_open_arrival.
                    double arvrate = openArrivalRateOf(tidx);
                    if (arvrate > GlobalConstants.FineTol) {
                        double njobs_arv = Matrix.extractRows(this.njobs, tidx, tidx + 1, null).elementMax();
                        if (!(njobs_arv > 0)) {
                            njobs_arv = this.lqn.maxmult.get(0, tidx);
                        }
                        this.tput.set(tidx, arvrate);
                        double host_residt = 0;
                        if (this.lqn.entriesof.get(tidx) != null) {
                            for (int eidx_arv : this.lqn.entriesof.get(tidx)) {
                                List<Integer> arvActs = this.lqn.actsof.get(eidx_arv);
                                if (arvActs != null) {
                                    for (int aidx_arv : arvActs) {
                                        if (!Double.isNaN(this.residt.get(aidx_arv))) {
                                            host_residt += this.residt.get(aidx_arv);
                                        }
                                    }
                                }
                            }
                        }
                        double zArv = FastMath.max(GlobalConstants.Zero,
                                njobs_arv / arvrate - host_residt - tidx_thinktime);
                        double omegaArv = this.relax_omega;
                        if (omegaArv < 1.0 && it > 1 && !Double.isNaN(this.thinkt_prev.get(tidx))) {
                            zArv = omegaArv * zArv + (1 - omegaArv) * this.thinkt_prev.get(tidx);
                        }
                        this.thinkt.set(tidx, zArv);
                        this.thinkt_prev.set(tidx, zArv);
                        this.thinktproc.put(tidx, Exp.fitMean(zArv + tidx_thinktime));
                        continue;
                    }
                    // Check if this is a forwarding target task
                    boolean isFwdTarget = false;
                    int fwd_cidx_found = -1;
                    int source_tidx = -1;
                    double fwd_prob = 0;
                    if (this.lqn.isref.get(tidx) == 0 && this.lqn.entriesof.get(tidx) != null) {
                        for (int eidx_fwd : this.lqn.entriesof.get(tidx)) {
                            for (int cidx_fwd = 0; cidx_fwd < this.lqn.ncalls; cidx_fwd++) {
                                if (this.lqn.calltype.get(cidx_fwd) == CallType.FWD && (int) this.lqn.callpair.get(cidx_fwd, 1) == eidx_fwd) {
                                    isFwdTarget = true;
                                    fwd_cidx_found = cidx_fwd;
                                    int source_eidx = (int) this.lqn.callpair.get(cidx_fwd, 0);
                                    source_tidx = (int) this.lqn.parent.get(0, source_eidx);
                                    fwd_prob = this.lqn.callproc_mean.getOrDefault(cidx_fwd, 1.0);
                                    break;
                                }
                            }
                            if (isFwdTarget) break;
                        }
                    }
                    if (isFwdTarget) {
                        // Forwarding target think time: derived from the source
                        // task's throughput, forwarding probability, and the
                        // processor response time for the target's activities.
                        // Formula: thinkt = njobs / arrival_rate - host_residt
                        // where host_residt is the response time at the processor.
                        double njobs_fwd = Matrix.extractRows(this.njobs, tidx, tidx + 1, null).elementMax();
                        double arrival_rate = this.tput.get(source_tidx) * fwd_prob;
                        if (arrival_rate > GlobalConstants.FineTol && njobs_fwd > 0) {
                            this.tput.set(tidx, arrival_rate);
                            // Subtract the processor response time for the target's
                            // activities (already computed by updateMetricsDefault)
                            int target_eidx = (int) this.lqn.callpair.get(fwd_cidx_found, 1);
                            double host_residt = 0;
                            List<Integer> targetActs = this.lqn.actsof.get(target_eidx);
                            if (targetActs != null) {
                                for (int aidx_fwd : targetActs) {
                                    host_residt += this.residt.get(aidx_fwd);
                                }
                            }
                            this.thinkt.set(tidx, FastMath.max(GlobalConstants.Zero, njobs_fwd / arrival_rate - host_residt - tidx_thinktime));
                        } else {
                            // Source throughput not yet available; use large think time
                            this.thinkt.set(tidx, 1000);
                        }
                        // Apply under-relaxation
                        double omega = this.relax_omega;
                        if (omega < 1.0 && it > 1 && !Double.isNaN(this.thinkt_prev.get(tidx))) {
                            this.thinkt.set(tidx, omega * this.thinkt.get(tidx) + (1 - omega) * this.thinkt_prev.get(tidx));
                        }
                        this.thinkt_prev.set(tidx, this.thinkt.get(tidx));
                        this.thinktproc.put(tidx, Exp.fitMean(this.thinkt.get(tidx) + tidx_thinktime));
                    } else { // set to zero if this is a ref task
                        this.thinkt.set(tidx, GlobalConstants.FineTol);
                        this.thinktproc.put(tidx, Immediate.getInstance());
                    }
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
            NetworkStruct sn = model.getStruct();
            if (sn.nregions > 0) {
                return regionCapableLayerSolver(model, defaultOptions.verbose);
            }
            return new SolverMVA(model, defaultOptions);
        }
    }

    /**
     * <p>Picks the first solver whose feature set covers a layer carrying an
     * admission constraint. The order is by decreasing accuracy: CTMC is exact
     * but state-space bound, LDES and SSA simulate. Selection is by
     * SolverX.supports so it self-corrects if another solver later declares
     * Region.</p>
     *
     * @param model   the layer model
     * @param verbose verbosity inherited from the layer solver options
     * @return a solver whose feature set covers the layer
     */
    /**
     * <p>Waiting time absorbed by an admission constraint in a layer, recovered by
     * Little's law from the layer population deficit. A job blocked at the
     * constraint is counted at no station (JMT WAITQ convention), so its wait is
     * absent from RN; without this the caller never sees the blocking and the
     * fixed point loses flow balance.</p>
     *
     * <p>The population is conserved per chain, not per class: a job in a layer
     * switches class along the activity graph, so the call class itself carries
     * population 0. Splitting the chain deficit by throughput gives every
     * region-visiting class the same wait.</p>
     *
     * @param eidxLayer zero-based ensemble index of the layer
     * @param nodeidx   zero-based station index of the constrained server
     * @param classidx  zero-based class index of the call
     * @param res       this iteration's layer result
     * @return the recovered waiting time, or 0 where the layer is unconstrained
     */
    /**
     * Throughput of the reference class of the chain holding CLASSIDX in a layer.
     *
     * <p>This is the normaliser {@code residt} uses, {@code TN(refstat,refclass)}, so a
     * quantity divided by it is per chain-reference visit. Feeding {@code entry_servt} a
     * mix of normalisations is what the {@code task_tput/entry_tput} rescaling then
     * mis-corrects; see _kb/06-solver-catalog.md.</p>
     *
     * @param eidxLayer zero-based ensemble index of the layer
     * @param classidx  zero-based class index
     * @param res       this iteration's layer result
     * @return the reference throughput, or 0 where it cannot be resolved
     */
    double chainRefTput(int eidxLayer, int classidx, SolverResult res) {
        if (eidxLayer < 0 || eidxLayer >= this.ensemble.length || res == null || res.TN == null) {
            return 0.0;
        }
        NetworkStruct layerSn = this.ensemble[eidxLayer].getStruct(false);
        if (layerSn == null || layerSn.refstat == null || layerSn.refclass == null) {
            return 0.0;
        }
        int chainIdx = -1;
        if (layerSn.chains != null && !layerSn.chains.isEmpty()) {
            for (int ch = 0; ch < layerSn.nchains; ch++) {
                if (layerSn.chains.get(ch, classidx) > 0) {
                    chainIdx = ch;
                    break;
                }
            }
        }
        if (chainIdx < 0) {
            return 0.0;
        }
        int refclass_c = (int) layerSn.refclass.get(chainIdx);
        int refstat_k = (int) layerSn.refstat.get(classidx);
        if (refstat_k < 0 || refstat_k >= res.TN.getNumRows() || refclass_c < 0 || refclass_c >= res.TN.getNumCols()) {
            return 0.0;
        }
        return res.TN.get(refstat_k, refclass_c);
    }

    double regionWait(int eidxLayer, int nodeidx, int classidx, SolverResult res) {
        if (this.layerHasRegion == null || eidxLayer < 0 || eidxLayer >= this.layerHasRegion.length
                || !this.layerHasRegion[eidxLayer]) {
            return 0.0;
        }
        Matrix chains = this.layerChains[eidxLayer];
        if (chains == null || chains.isEmpty()) {
            return 0.0;
        }
        int cix = -1;
        for (int c = 0; c < chains.getNumRows(); c++) {
            if (chains.get(c, classidx) > 0) {
                cix = c;
                break;
            }
        }
        if (cix < 0) {
            return 0.0;
        }
        double chainPop = 0;
        double qsum = 0;
        double xregion = 0;
        for (int k = 0; k < chains.getNumCols(); k++) {
            if (chains.get(cix, k) == 0) {
                continue;
            }
            JobClass jc = this.ensemble[eidxLayer].getClasses().get(k);
            if (jc instanceof ClosedClass) {
                chainPop += ((ClosedClass) jc).getPopulation();
            }
            for (int i = 0; i < res.QN.getNumRows(); i++) {
                qsum += res.QN.get(i, k);
            }
            xregion += res.TN.get(nodeidx, k);
        }
        double deficit = chainPop - qsum;
        if (Double.isFinite(deficit) && deficit > 0 && xregion > GlobalConstants.FineTol) {
            return deficit / xregion;
        }
        return 0.0;
    }

    static NetworkSolver regionCapableLayerSolver(Network model, VerboseLevel verbose) {
        FeatureSet featUsed = model.getUsedLangFeatures();
        if (FeatureSet.supports(SolverCTMC.getFeatureSet(), featUsed)) {
            SolverOptions o = new jline.solvers.ctmc.CTMCOptions();
            o.verbose = verbose;
            return new SolverCTMC(model, o);
        }
        if (FeatureSet.supports(jline.solvers.ldes.SolverLDES.getFeatureSet(), featUsed)) {
            SolverOptions o = new jline.solvers.ldes.LDESOptions();
            o.verbose = verbose;
            return new jline.solvers.ldes.SolverLDES(model, o);
        }
        if (FeatureSet.supports(jline.solvers.ssa.SolverSSA.getFeatureSet(), featUsed)) {
            SolverOptions o = new jline.solvers.ssa.SSAOptions();
            o.verbose = verbose;
            return new jline.solvers.ssa.SolverSSA(model, o);
        }
        throw new RuntimeException("LN layer " + model.getName()
                + " carries an admission constraint but none of SolverCTMC, SolverLDES, SolverSSA supports it. Supply a layer solver factory explicitly.");
    }

    protected static class recurActGraphReturnType {
        public JobClass curClass;
        public int jobPos;

        public RoutingMatrix P;

        // stations the job currently sits at, empty when it is at the client;
        // under flat layering successive activities may sit at different servers
        public Map<Integer, Queue> curStations;

        public recurActGraphReturnType(JobClass curClass, int jobPos, RoutingMatrix P) {
            this(curClass, jobPos, P, null);
        }

        public recurActGraphReturnType(JobClass curClass, int jobPos, RoutingMatrix P, Map<Integer, Queue> curStations) {
            this.curClass = curClass;
            this.jobPos = jobPos;
            this.P = P;
            this.curStations = curStations;
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
        state.callresidtPrev = this.callresidt_prev != null ? this.callresidt_prev.copy() : null;

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
        if (state.callresidtPrev != null) this.callresidt_prev = state.callresidtPrev.copy();

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
            // refreshChains can change the chain basis, invalidating the
            // warm-start solution cached by analyze()
            solvers[e].options.init_sol = new Matrix(0, 0);
            if (solvers[e] instanceof jline.solvers.mva.SolverMVA) {
                ((jline.solvers.mva.SolverMVA) solvers[e]).resetForkWarmStart();
            }
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
            solvers[e] = silenced(newSolverFactory.at(ensemble[e]));
            assertLayerSolverSupportsModel(solvers[e], ensemble[e], e);
        }
    }

    /**
     * Sets a layer solver to {@link VerboseLevel#SILENT} and returns it.
     *
     * <p>A LAYER SOLVER NEVER NARRATES. The fixed point runs every layer once
     * per iteration, so a layer left at the caller's verbosity prints its own
     * banner nlayers*iter_max times and buries the layered narration the caller
     * actually asked for. The level is stamped HERE rather than in the factory
     * because a factory the USER supplied never sees the LN options at all, and
     * stamping it in {@code DefaultSolverFactory} alone left exactly that case
     * loud.</p>
     *
     * <p>SolverLN's own reporting is unaffected: it reads {@code this.options.verbose},
     * not the layer's.</p>
     *
     * @param solver the freshly built layer solver, possibly null
     * @return the same solver, silenced
     */
    private static NetworkSolver silenced(NetworkSolver solver) {
        if (solver != null && solver.options != null) {
            solver.options.verbose = VerboseLevel.SILENT;
        }
        return solver;
    }

    /**
     * Reject a layer solver that cannot represent its layer model, upfront
     * with a clear message.
     * <p>
     * LN server-layer stations carry immediate feedback (sn.immfeed):
     * successive same-host activities retain the server, modelled as
     * immediate-feedback self-loops so the layer solver does not re-queue the
     * job. SolverJMT rejects any model with immfeed and returns no solution, so
     * a pure-JMT layer factory otherwise fails cryptically at the first
     * iteration. Detect it here instead.
     * </p>
     * <p>
     * The guard is CONDITIONAL: it fires only when the specific layer model
     * actually carries immfeed. A SolverJMT layer solver on an immfeed-free
     * layer is allowed, and non-JMT factories are never rejected.
     * </p>
     *
     * @param layerSolver the solver produced by the factory for this layer
     * @param layerModel  the layer network model
     * @param e           the 0-based layer index
     */
    private void assertLayerSolverSupportsModel(NetworkSolver layerSolver, Network layerModel, int e) {
        if (layerSolver instanceof SolverJMT) {
            NetworkStruct lsn = layerModel.getStruct(false);
            if (lsn.immfeed != null && lsn.immfeed.elementSum() > 0) {
                line_error(mfilename(new Object() {}), String.format(
                    "SolverJMT cannot solve LN layer %d: the layer carries immediate feedback (sn.immfeed), which SolverJMT does not support, so LN would fail at the first iteration. Use the default layer factory (MVA/NC) or another layer solver that supports immediate feedback.",
                    e));
            }
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

        double S1 = this.servt_ph1.get(eidx);
        double S2 = this.servt_ph2.get(eidx);

        // Get throughput - use entry if available, otherwise use task
        double lambda;
        if (this.tput.get(eidx) > GlobalConstants.FineTol) {
            lambda = this.tput.get(eidx);
        } else if (tidx > 0 && this.tput.get(tidx) > GlobalConstants.FineTol) {
            lambda = this.tput.get(tidx);
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
            // Single server: exact CTMC solution for states {idle, phase-1, phase-2}
            //   Q = | -lambda   lambda    0   |
            //       |    0      -mu1     mu1   |
            //       |   mu2       0     -mu2   |
            // Solving pi*Q = 0, sum(pi) = 1 gives the closed-form stationary
            // distribution of this 3-state birth-death cycle:
            //   pi0 = mu1*mu2 / D, pi1 = lambda*mu2 / D, pi2 = lambda*mu1 / D
            // where D = lambda*mu2 + mu1*mu2 + lambda*mu1. By PASTA the overtaking
            // probability is P(find in phase-2) = pi2.
            double denom = lambda * mu2 + mu1 * mu2 + lambda * mu1;
            prOt = lambda * mu1 / denom;
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
            for (int predIdx = 0; predIdx < graph.getNumRows(); predIdx++) {
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
        for (int predIdx = 0; predIdx < graph.getNumRows(); predIdx++) {
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
        for (int predIdx = 0; predIdx < graph.getNumRows(); predIdx++) {
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
        for (int succIdx = 0; succIdx < graph.getNumCols(); succIdx++) {
            if (graph.get(aidx, succIdx) != 0 && isPostAndActivity(succIdx)) {
                count++;
            }
        }
        return count;
    }

    /**
     * The activities belonging to each branch of an AND-join.
     *
     * <p>A branch is recovered by walking backwards from each immediate predecessor of the join
     * until an activity marked POST_AND is reached, that activity being the branch head spawned
     * by the AND-fork. Branches between a fork and its join are disjoint paths, so the walk is
     * unambiguous.</p>
     *
     * @param joinaidx global index of the join target activity
     * @return one list of global activity indices per branch
     */
    private List<List<Integer>> branchMembers(int joinaidx) {
        List<List<Integer>> members = new ArrayList<>();
        Matrix graph = this.lqn.graph;
        int ashift = this.lqn.ashift;
        int nacts = this.lqn.nacts;

        for (int tail = 0; tail < graph.getNumRows(); tail++) {
            if (graph.get(tail, joinaidx) <= 0) {
                continue;
            }
            if (tail <= ashift || tail > ashift + nacts) {
                continue; // not an activity
            }
            List<Integer> chain = new ArrayList<>();
            chain.add(tail);
            int cur = tail;
            for (int guard = 0; guard < nacts; guard++) {
                if (cur < this.lqn.actposttype.getNumCols() &&
                    this.lqn.actposttype.get(0, cur) == ActivityPrecedenceType.ID_POST_AND) {
                    break; // branch head
                }
                int prev = -1;
                int nprev = 0;
                for (int p = 0; p < graph.getNumRows(); p++) {
                    if (p != cur && graph.get(p, cur) > 0 && p > ashift && p <= ashift + nacts) {
                        prev = p;
                        nprev++;
                    }
                }
                if (nprev != 1) {
                    break; // a merge or the start of the graph
                }
                cur = prev;
                chain.add(cur);
            }
            members.add(chain);
        }
        return members;
    }

    /**
     * Compute the completion time of every AND-join and the correction it implies.
     *
     * <p>The branches of an AND-fork run concurrently, so the time to pass the join is the k-th
     * smallest of the branch completion times, k being the quorum of the join (k equals the
     * branch count when the join waits for all its branches). Times are taken over residt
     * because that is the quantity entry_servt aggregates.</p>
     *
     * @return a row matrix holding, per join target, the join time minus the sequential sum of
     *         its branch times, i.e. the amount by which the reachability matrix overcounts
     */
    private Matrix updateJoinDelays() {
        int nidx = this.lqn.nidx;
        this.joint = new Matrix(1, nidx, nidx);
        Matrix excess = new Matrix(1, nidx, nidx);
        if (this.lqn.actpretype == null) {
            return excess;
        }

        int ashift = this.lqn.ashift;
        int lastAct = Math.min(ashift + this.lqn.nacts, this.lqn.actpretype.getNumCols()) - 1;
        for (int aidx = ashift; aidx <= lastAct; aidx++) {
            // PRE_AND marks the branch tails, not the join target, so the join is the
            // activity whose predecessors carry that mark.
            if (!isJoinTarget(aidx)) {
                continue;
            }
            List<List<Integer>> branches = branchMembers(aidx);
            int n = branches.size();
            if (n == 0) {
                continue;
            }
            double[] branchTimes = new double[n];
            double total = 0.0;
            for (int b = 0; b < n; b++) {
                double s = 0.0;
                for (int m : branches.get(b)) {
                    s += this.residt.get(m);
                    // A branch activity with an Immediate host demand does all its work in a
                    // rendezvous, so residt alone would make this correction vanish silently.
                    if (this.lqn.callsof != null && this.callresidt != null) {
                        List<Integer> mcalls = this.lqn.callsof.get(m);
                        if (mcalls != null) {
                            for (int cidx : mcalls) {
                                if (this.lqn.calltype != null && this.lqn.calltype.get(cidx) == CallType.SYNC
                                        && cidx - 1 < this.callresidt.getNumCols()) {
                                    s += this.callresidt.get(cidx);
                                }
                            }
                        }
                    }
                }
                branchTimes[b] = s;
                total += s;
            }
            if (n == 1) {
                this.joint.set(0, aidx, branchTimes[0]);
                continue;
            }
            int quorum = n;
            if (this.lqn.actquorum != null && aidx < this.lqn.actquorum.getNumCols()) {
                int q = (int) this.lqn.actquorum.get(0, aidx);
                if (q >= 1 && q <= n) {
                    quorum = q;
                }
            }
            // Branch times are taken as exponential, so the variance is the square of the mean.
            double[] vars = new double[n];
            for (int b = 0; b < n; b++) {
                vars[b] = branchTimes[b] * branchTimes[b];
            }
            double jt = FJ_quorum.quorumMoments(branchTimes, vars, quorum)[0];
            this.joint.set(0, aidx, jt);
            excess.set(0, aidx, jt - total);
        }
        return excess;
    }

    /**
     * Count the PRE_AND predecessors of an activity, i.e. the number of branches feeding a join.
     */
    private int countPreAndPredecessors(int aidx) {
        if (this.lqn.actpretype == null || this.lqn.graph == null) {
            return 0;
        }
        Matrix graph = this.lqn.graph;
        int count = 0;
        for (int predIdx = 0; predIdx < graph.getNumRows(); predIdx++) {
            if (predIdx != aidx && graph.get(predIdx, aidx) != 0 &&
                predIdx < this.lqn.actpretype.getNumCols() &&
                this.lqn.actpretype.get(0, predIdx) == ActivityPrecedenceType.ID_PRE_AND) {
                count++;
            }
        }
        return count;
    }

    /**
     * Apply the quorum of an AND-join to the Join node for the class active at the fork.
     *
     * <p>A join whose quorum equals its branch count waits for all branches, which is already
     * the default JoinStrategy.STD, so nothing is set in that case. Only a genuine quorum
     * k &lt; n switches the node to JoinStrategy.Quorum.</p>
     *
     * @param joinNode      the Join node of the layer
     * @param forkClass     the job class active at the corresponding fork
     * @param joinTargetIdx global index of the join target activity
     */
    private void applyJoinQuorum(Join joinNode, JobClass forkClass, int joinTargetIdx) {
        if (joinNode == null || forkClass == null || this.lqn.actquorum == null) {
            return;
        }
        if (joinTargetIdx < 1 || joinTargetIdx >= this.lqn.actquorum.getNumCols()) {
            return;
        }
        int quorum = (int) this.lqn.actquorum.get(0, joinTargetIdx);
        int nbranches = countPreAndPredecessors(joinTargetIdx);
        if (quorum < 1 || nbranches < 1 || quorum >= nbranches) {
            return;
        }
        joinNode.setStrategy(forkClass, JoinStrategy.Quorum);
        joinNode.setRequired(forkClass, quorum);
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
        public Matrix callresidtPrev;

        public Map<Integer, Map<Integer, SolverResult>> results;

        public Matrix njobs;
        public Matrix ptaskcallers;
        public Matrix ilscaling;
    }
    /**
     * Think time of an activity, zero when it has none.
     * <p>
     * An activity think time is a delay in series with that activity's host
     * demand, held at the activity's own task: the task keeps its thread for the
     * whole hostdem+thinktime interval, so it serializes against the task
     * multiplicity, but the host processor is released for it. This mirrors lqns,
     * whose think-time attribute LINE already writes out.
     * <p>
     * actthink_mean is absent for an activity that was never given a think time,
     * and may hold NaN, so the value is filtered here rather than added by the
     * caller: a NaN reaching servt or residt propagates into the layer solvers,
     * whose NaN-ignoring fallbacks then launder it into a plausible-looking
     * bare-service figure instead of failing.
     *
     * @param lqn the layered network structure
     * @param aidx absolute index of the activity
     * @return the think time, or 0 when the activity has none
     */
    private static double actThinkTime(LayeredNetworkStruct lqn, int aidx) {
        if (lqn.actthink_mean == null) {
            return 0;
        }
        Double v = lqn.actthink_mean.get(aidx);
        if (v == null || Double.isNaN(v) || v <= GlobalConstants.FineTol) {
            return 0;
        }
        return v;
    }

    /**
     * Declared think time of a task as it enters the thread cycle: the value for
     * a REFERENCE task, zero for any other.
     *
     * A think time is an attribute of the closed customer population a reference
     * task stands for, and it is what separates one request of that population
     * from the next. On a served task it has no such meaning, and charging it
     * per request throttles the task: lqn_basic's T3 has 25 threads and a
     * declared think time of 4, and reading it as a per-request delay caps it at
     * 25/(4+0.02) = 6.219 completions per second. Three independent oracles put
     * the rate at five calls per caller request instead -- lqsim 66.5, LDES
     * 66.955, lqns 75.6 -- so a non-reference task's think time does not enter
     * the cycle. See _kb/06-solver-catalog.md (LN section).
     *
     * @param lqn the layered network structure
     * @param tidx absolute index of the task
     * @return the think time of a reference task, 0 otherwise
     */
    private static double refThinkTime(LayeredNetworkStruct lqn, int tidx) {
        if (lqn.isref == null || lqn.isref.get(tidx) == 0) {
            return 0;
        }
        double v = lqn.think_mean.getOrDefault(tidx, 0.0);
        if (Double.isNaN(v) || v < 0) {
            return 0;
        }
        return v;
    }

    // =====================================================================
    // Method 'srvn.ph': the activity graph of an entry as a phase-type server law
    //
    // Each layer is a two-station cycle, Delay('Clients') + Queue(server), with
    // one closed class per caller task. The sequencing the default method
    // encodes as routing -- one class per entry, per activity and per call, plus
    // Fork, Join, Router and ClassSwitch nodes -- is composed instead into a
    // single phase-type service law per (layer, caller), by the exact
    // series-parallel reduction of Workflow. The layer therefore carries only
    // the client/server back-and-forth, and the activity graph survives as a
    // distribution.
    //
    // Twin of the MATLAB @SolverLN/buildLayersPH.m, updateLayersPH.m,
    // updateMetricsPH.m, updateThinkTimesPH.m, phComposeEntryLaws.m and
    // getEnsembleAvgPH.m, of the Python solver_ln_ph routines and of the
    // C++ *_ph members of solvers/ln/solver_ln.h. See
    // _kb/06-solver-catalog.md (LN section) for the layering taxonomy.
    // =====================================================================
    /** One two-station phLayers, and the caller classes that cycle through it. */
    public static class PHLayer {
        /** Element index of the server, a processor or a task. */
        public int idx;
        /** True when the server is a processor. */
        public boolean ishost;
        /** Caller tasks with a closed class in this layer. */
        public int[] callers;
        /** Class index of each caller task, -1 when absent. */
        public int[] classOfCaller;
        /** Replicas of the server station. */
        public int nreplicas;
        /** Station indices (1-based) of the server replicas. */
        public int[] qstations;
        /** Mean of the law each class is currently served with. */
        public double[] svcmeanByClass;
        /** Open classes: {class index, entry index} or {class index, -call index}. */
        public int[][] openArrivals;
        /**
         * Closed population of the MODEL this server sits in. Under 'flat.ph'
         * that is every caller of the single network, not only the callers of
         * this one station, so it is recorded here rather than recomputed.
         */
        public double npop;
    }


    /** True once phInitLaws has composed the per-entry workflows. */
    private boolean phLawsReady;

    // --- per-entry workflows and their composed laws
    Workflow[] phWf;
    Workflow[] phWfhost;
    double[][] phExecs;      // [eidx][aidx]
    double[][] phCallexecs;  // [eidx][cidx]
    Matrix[] phHostalpha;
    Matrix[] phHostT;
    double[] phHostmean;
    Matrix[] phEntryalpha;
    Matrix[] phEntryT;
    double[] phEntrymean;
    double[] phEntryscv;
    double[] phShare;
    double[] phOverlap;
    double[] phSetupshare;

    // --- caller-side aggregates
    double[] phXdemand;
    double[][] phNcalls;    // [caller task][called entry]
    double[][] phCalltime;  // [caller task][called task]
    double[] phProcresid;
    double[] phActthinkt;
    double[] phCalltotal;

    PHLayer[] phLayers;      // indexed by element idx


    // =====================================================================
    // Layer construction
    // =====================================================================

    /**
     * Build the ensemble under method 'srvn.ph'.
     */
    /**
     * Answer whether method 'srvn.ph' can serve this model, without disturbing
     * the solver.
     *
     * Both the feature gate and the series-parallel reduction can refuse, and
     * the second only finds out by composing the per-entry workflows -- work the
     * build then reuses, since those laws do not depend on the iterate. Used by
     * the alias 'srvn' to choose between 'srvn.ph' and 'srvn.cs'.
     *
     * @return true when the layers can be built as composed phase-type ones
     */
    private boolean probeSrvnPH() {
        try {
            assertSrvnPHSupported();
            phInitLaws();
            this.phLawsReady = true;
            return true;
        } catch (RuntimeException e) {
            this.phLawsReady = false;
            line_debug(options.verbose,
                    "LN: method=srvn cannot use srvn.ph on this model (" + e.getMessage() + ")");
            return false;
        }
    }

    public void buildLayersPH() {
        buildLayersPH(false);
    }

    /**
     * Build the ensemble of a PH encoding. FLAT false is method 'srvn.ph', one
     * layer per served element; FLAT true is method 'flat.ph', ONE layer holding
     * a station for every processor and every called task, with the same one
     * closed class per caller task.
     *
     * @param flat true to squash every server into a single submodel
     */
    public void buildLayersPH(boolean flat) {
        int nelem = lqn.nhosts + lqn.ntasks;
        if (!this.phLawsReady) {
            assertSrvnPHSupported(flat);
        }

        this.cell_thinkt_classes_updmap = new java.util.HashMap<Integer, List<Integer[]>>();
        this.cell_servt_classes_updmap = new java.util.HashMap<Integer, List<Integer[]>>();
        this.cell_call_classes_updmap = new java.util.HashMap<Integer, List<Integer[]>>();
        this.cell_arvproc_classes_updmap = new java.util.HashMap<Integer, List<Integer[]>>();
        this.cell_route_prob_updmap = new java.util.HashMap<Integer, List<Integer[]>>();

        // The interlock correction rewrites the populations of the call classes,
        // which this method does not create: its callers reach the server in one
        // class each
        if (this.options.config.interlocking) {
            this.options.config.interlocking = false;
        }

        // A preceding probe has already composed the per-entry workflows; they do
        // not depend on the iterate, so they are not rebuilt here.
        if (!this.phLawsReady) {
            phInitLaws();
        }

        // Seed the fixed point with the static demands, then compose the entry laws
        this.residt = new Matrix(1, lqn.nidx, lqn.nidx);
        this.servt = new Matrix(1, lqn.nidx, lqn.nidx);
        this.callservt = new Matrix(1, FastMath.max(lqn.ncalls, 1), FastMath.max(lqn.ncalls, 1));
        this.callresidt = new Matrix(1, FastMath.max(lqn.ncalls, 1), FastMath.max(lqn.ncalls, 1));
        for (int aidx = lqn.ashift; aidx < lqn.ashift + lqn.nacts; aidx++) {
            this.residt.set(aidx, lqn.hostdem_mean.getOrDefault(aidx, 0.0));
        }
        for (int cidx = 0; cidx < lqn.ncalls; cidx++) {
            if (lqn.calltype.get(cidx) == CallType.SYNC || lqn.calltype.get(cidx) == CallType.ASYNC) {
                int eidx = (int) lqn.callpair.get(cidx, 1);
                double v = lqn.callproc_mean.getOrDefault(cidx, 0.0) * phHostmean[eidx];
                this.callservt.set(cidx, v);
                this.callresidt.set(cidx, v);
            }
        }
        phComposeEntryLaws();

        Network[] built = new Network[nelem];

        if (flat) {
            // ONE subnetwork holding every processor and every called task
            List<Integer> servers = buildPHFlatLayer();
            this.idxhash = new ArrayList<Double>();
            for (int i = 0; i < nelem; i++) {
                this.idxhash.add(Double.NaN);
            }
            for (int idx : servers) {
                this.idxhash.set(idx, 1.0);
            }
            this.nlayers = 1;
            this.ensemble = new Network[]{phFlatModel};
            if (lqnModel != null) {
                lqnModel.setEnsemble(new ArrayList<Network>(Arrays.asList(this.ensemble)));
            }
            this.hostLayerIndices = new ArrayList<Integer>();
            this.taskLayerIndices = new ArrayList<Integer>();
            this.hostLayerIndices.add(0);
            this.taskLayerIndices.add(0);
            this.layerHasRegion = new boolean[1];
            this.layerChains = new Matrix[1];
            this.thinkt = new Matrix(1, lqn.ntasks + lqn.tshift, lqn.ntasks + lqn.tshift);
            this.tput = new Matrix(1, lqn.nidx, lqn.nidx);
            this.util = new Matrix(1, lqn.nidx, lqn.nidx);
            updateLayersPH(0);
            this.thinkt_classes_updmap = integerMapToMatrix(this.cell_thinkt_classes_updmap);
            this.servt_classes_updmap = integerMapToMatrix(this.cell_servt_classes_updmap);
            this.call_classes_updmap = integerMapToMatrix(this.cell_call_classes_updmap);
            this.arvproc_classes_updmap = integerMapToMatrix(this.cell_arvproc_classes_updmap);
            this.route_prob_updmap = integerMapToMatrix(this.cell_route_prob_updmap);
            return;
        }

        // One subnetwork per processor
        for (int hidx = 0; hidx < lqn.nhosts; hidx++) {
            if (this.ignore.get(hidx) != 0) {
                continue;
            }
            List<Integer> callers = phHostLayerCallers(hidx);
            if (callers.isEmpty()) {
                continue;
            }
            built[hidx] = buildPHLayer(hidx, callers, true);
        }

        // One subnetwork per called task
        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            if (this.ignore.get(tidx) != 0 || lqn.isref.get(tidx) != 0) {
                continue;
            }
            List<Integer> callers = phTaskLayerCallers(tidx);
            List<Integer> asyncCalls = phAsyncCallsInto(tidx);
            if (callers.isEmpty() && asyncCalls.isEmpty()) {
                continue;
            }
            built[tidx] = buildPHLayer(tidx, callers, false);
        }

        // Compact the ensemble and index it
        this.idxhash = new ArrayList<Double>();
        List<Network> models = new ArrayList<Network>();
        for (int idx = 0; idx < nelem; idx++) {
            if (built[idx] == null) {
                this.idxhash.add(Double.NaN);
            } else {
                models.add(built[idx]);
                this.idxhash.add((double) models.size());
            }
        }
        this.nlayers = models.size();
        Network[] arr = new Network[this.nlayers];
        for (int i = 0; i < this.nlayers; i++) {
            arr[i] = models.get(i);
        }
        this.ensemble = arr;
        if (lqnModel != null) {
            lqnModel.setEnsemble(new ArrayList<Network>(Arrays.asList(arr)));
        }

        this.hostLayerIndices = new ArrayList<Integer>();
        this.taskLayerIndices = new ArrayList<Integer>();
        for (int hidx = 0; hidx < lqn.nhosts; hidx++) {
            if (!Double.isNaN(this.idxhash.get(hidx))) {
                this.hostLayerIndices.add(this.idxhash.get(hidx).intValue() - 1);
            }
        }
        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            if (!Double.isNaN(this.idxhash.get(tidx))) {
                this.taskLayerIndices.add(this.idxhash.get(tidx).intValue() - 1);
            }
        }

        this.layerHasRegion = new boolean[this.nlayers];
        this.layerChains = new Matrix[this.nlayers];

        // install the initial laws, so that iteration 1 sees the seeded demands
        // rather than the placeholders the stations were created with
        this.thinkt = new Matrix(1, lqn.ntasks + lqn.tshift, lqn.ntasks + lqn.tshift);
        this.tput = new Matrix(1, lqn.nidx, lqn.nidx);
        this.util = new Matrix(1, lqn.nidx, lqn.nidx);
        updateLayersPH(0);

        // The maps carry no law of this method -- updateLayersPH composes them
        // instead -- but post() resets the layers they name, so they are filled
        this.thinkt_classes_updmap = integerMapToMatrix(this.cell_thinkt_classes_updmap);
        this.servt_classes_updmap = integerMapToMatrix(this.cell_servt_classes_updmap);
        this.call_classes_updmap = integerMapToMatrix(this.cell_call_classes_updmap);
        this.arvproc_classes_updmap = integerMapToMatrix(this.cell_arvproc_classes_updmap);
        this.route_prob_updmap = integerMapToMatrix(this.cell_route_prob_updmap);
    }

    /** The single squashed network of method 'flat.ph', null under 'srvn.ph'. */
    private Network phFlatModel;

    /**
     * Build the ONE layer of method 'flat.ph': a client delay plus a station for
     * every processor and every called task.
     *
     * A caller task is one closed class, and it visits each server it uses ONCE
     * per invocation, carrying there the composed law of the demand it places on
     * that server -- the same law method 'srvn.ph' installs in the server's own
     * layer. What changes is that the servers now contend inside one network
     * instead of seeing each other through surrogate delays, so the client delay
     * keeps only the think times and whatever of the cycle this model does not
     * hold. That is the whole difference between the two encodings of the PH
     * composition, and it is why the reconstruction passes are shared verbatim.
     *
     * @return the served element indices, in station order
     */
    private List<Integer> buildPHFlatLayer() {
        List<Integer> servers = phFlatServerSet();
        int nsrv = servers.size();

        Network model = new Network("FlatPH");
        model.setChecks(false);
        model.getAttribute().setClientIdx(1);
        model.getAttribute().setSourceIdx(-1);

        Delay clientDelay = new Delay(model, "Clients");

        Queue[] srvStation = new Queue[nsrv];
        int[] stationOf = new int[lqn.nidx + 1];
        for (int si = 0; si < nsrv; si++) {
            int idx = servers.get(si);
            boolean ishost = idx < lqn.nhosts;
            Queue q = new Queue(model, lqn.hashnames.get(idx), lqn.sched.get(idx));
            q.setNumberOfServers((int) lqn.maxmult.get(0, idx));
            q.getAttribute().setIsHost(ishost);
            q.getAttribute().setIdx(idx);
            srvStation[si] = q;
            int stn = model.getNodeIndex(q) + 1;
            stationOf[idx] = stn;
            model.getAttribute().putServerIdxOf(idx, stn);
            if (ishost) {
                model.getAttribute().getHostStations().add(stn);
                model.getAttribute().addHosts(new Integer[]{null, stn});
            } else {
                model.getAttribute().getTaskStations().add(stn);
                model.getAttribute().addTasks(new Integer[]{null, stn});
            }
        }
        // the scalar fallback of stationIdxOf, which no served element reaches here
        model.getAttribute().setServerIdx(stationOf[servers.get(0)]);

        // Callers of each server, and the union of them, which becomes the class set
        Map<Integer, List<Integer>> callersOf = new HashMap<Integer, List<Integer>>();
        List<Integer> allCallers = new ArrayList<Integer>();
        for (int si = 0; si < nsrv; si++) {
            int idx = servers.get(si);
            List<Integer> cs = (idx < lqn.nhosts) ? phHostLayerCallers(idx) : phTaskLayerCallers(idx);
            callersOf.put(idx, cs);
            for (int c : cs) {
                if (!allCallers.contains(c)) {
                    allCallers.add(c);
                }
            }
        }
        Collections.sort(allCallers);

        // One closed class per caller task
        int[] classOfCaller = new int[lqn.nidx + 1];
        Arrays.fill(classOfCaller, -1);
        double npop = 0;
        for (int c : allCallers) {
            // phFlatServerSet has refused every replicated element, so the
            // per-replica reduction the srvn builder makes is the identity here
            double njobsC = phLayerPopulation(servers.get(0), c, 1);
            ClosedClass cls = new ClosedClass(model, lqn.hashnames.get(c), njobsC, clientDelay);
            cls.setReferenceClass(true);
            cls.setAttribute(new Integer[]{LayeredNetworkElement.TASK, c});
            classOfCaller[c] = cls.getIndex();
            model.getAttribute().addTasks(new Integer[]{cls.getIndex(), c});
            npop += njobsC;
            clientDelay.setService(cls, Exp.fitMean(FastMath.max(GlobalConstants.FineTol, refThinkTime(lqn, c))));
            // A station this caller never reaches must say so with Disabled, NOT
            // with a tiny placeholder law. An FCFS station carries ONE service law
            // across its classes, so a placeholder is not inert there: it is mixed
            // into the multiserver correction and invents waiting where there is
            // none. Under 'srvn.ph' the question never arises, since every class
            // of a layer visits that layer's single server.
            for (int si = 0; si < nsrv; si++) {
                srvStation[si].setService(cls, Disabled.getInstance());
            }
            for (int si = 0; si < nsrv; si++) {
                int idx = servers.get(si);
                if (!callersOf.get(idx).contains(c)) {
                    continue;
                }
                srvStation[si].setService(cls, Exp.fitMean(GlobalConstants.FineTol));
                this.njobs.set(c, idx, njobsC);
                phAddUpdMap(this.cell_thinkt_classes_updmap, idx, new Integer[]{idx, c, 1, cls.getIndex()});
                phAddUpdMap(this.cell_servt_classes_updmap, idx,
                        new Integer[]{idx, c, stationOf[idx], cls.getIndex()});
            }
        }

        // Open classes: entry arrivals on a processor station, async calls on a task one
        Map<Integer, List<int[]>> openArrivalsOf = new HashMap<Integer, List<int[]>>();
        for (int si = 0; si < nsrv; si++) {
            openArrivalsOf.put(servers.get(si), new ArrayList<int[]>());
        }
        Source sourceStation = null;
        Sink sinkStation = null;
        for (int si = 0; si < nsrv; si++) {
            int hidx = servers.get(si);
            if (hidx >= lqn.nhosts) {
                continue;
            }
            for (int c : callersOf.get(hidx)) {
                // A task no other task calls has no task station, so
                // updateThinkTimesPH never gives its caller class a surrogate
                // delay: the class cycles against an Immediate one and an open
                // stream on top of it doubles the load. The chain is the
                // representation that honours the thread pool, so it is kept and
                // closed on the arrival rate instead.
                if (phOpenArrivalOnly(c)) {
                    continue;
                }
                for (int eidx : phEntriesOf(c)) {
                    if (!phHasOpenArrival(eidx)) {
                        continue;
                    }
                    if (sourceStation == null) {
                        model.getAttribute().setSourceIdx(model.getNumberOfNodes() + 1);
                        sourceStation = new Source(model, "Source");
                        sinkStation = new Sink(model, "Sink");
                    }
                    OpenClass ocls = new OpenClass(model, lqn.hashnames.get(eidx) + ".Open", 0);
                    ocls.setAttribute(new Integer[]{LayeredNetworkElement.ENTRY, eidx});
                    sourceStation.setArrival(ocls, lqn.arrival.get(eidx));
                    clientDelay.setService(ocls, Disabled.getInstance());
                    // Disabled, not a placeholder, at every station this stream misses
                    for (int s2 = 0; s2 < nsrv; s2++) {
                        srvStation[s2].setService(ocls, Disabled.getInstance());
                    }
                    srvStation[si].setService(ocls,
                            Exp.fitMean(FastMath.max(GlobalConstants.FineTol, phHostmean[eidx])));
                    openArrivalsOf.get(hidx).add(new int[]{ocls.getIndex(), eidx});
                    model.getAttribute().addEntries(new Integer[]{ocls.getIndex(), eidx});
                    phAddUpdMap(this.cell_arvproc_classes_updmap, hidx, new Integer[]{hidx, -eidx,
                            model.getNodeIndex(sourceStation) + 1, ocls.getIndex()});
                }
            }
        }
        for (int si = 0; si < nsrv; si++) {
            int tidx = servers.get(si);
            if (tidx < lqn.nhosts) {
                continue;
            }
            for (int cidx : phAsyncCallsInto(tidx)) {
                if (sourceStation == null) {
                    model.getAttribute().setSourceIdx(model.getNumberOfNodes() + 1);
                    sourceStation = new Source(model, "Source");
                    sinkStation = new Sink(model, "Sink");
                }
                OpenClass ocls = new OpenClass(model, lqn.callhashnames.get(cidx), 0);
                ocls.setAttribute(new Integer[]{LayeredNetworkElement.CALL, cidx});
                sourceStation.setArrival(ocls, Immediate.getInstance());
                clientDelay.setService(ocls, Disabled.getInstance());
                // Disabled, not a placeholder, at every station this stream misses
                for (int s2 = 0; s2 < nsrv; s2++) {
                    srvStation[s2].setService(ocls, Disabled.getInstance());
                }
                int eidx = (int) lqn.callpair.get(cidx, 1);
                srvStation[si].setService(ocls,
                        Exp.fitMean(FastMath.max(GlobalConstants.FineTol, phEntrymean[eidx])));
                openArrivalsOf.get(tidx).add(new int[]{ocls.getIndex(), -cidx});
                model.getAttribute().addCalls(new Integer[]{ocls.getIndex(), cidx,
                        (int) lqn.callpair.get(cidx, 0), eidx});
                phAddUpdMap(this.cell_arvproc_classes_updmap, tidx, new Integer[]{tidx, cidx,
                        model.getNodeIndex(sourceStation) + 1, ocls.getIndex()});
                phAddUpdMap(this.cell_call_classes_updmap, tidx,
                        new Integer[]{tidx, cidx, stationOf[tidx], ocls.getIndex()});
            }
        }
        if (sourceStation != null) {
            for (JobClass jc : model.getClasses()) {
                if (jc instanceof ClosedClass) {
                    sourceStation.setArrival(jc, Disabled.getInstance());
                }
            }
        }

        // Routing: one visit per server the caller uses, in server order. The
        // number of calls is carried by the service law, not by a visit ratio,
        // so no arc ever moves.
        RoutingMatrix P = model.initRoutingMatrix();
        for (int c : allCallers) {
            JobClass cls = model.getClassByIndex(classOfCaller[c] - 1);
            Station prev = clientDelay;
            boolean visited = false;
            for (int si = 0; si < nsrv; si++) {
                if (!callersOf.get(servers.get(si)).contains(c)) {
                    continue;
                }
                P.set(cls, cls, prev, srvStation[si], 1.0);
                prev = srvStation[si];
                visited = true;
            }
            if (visited) {
                P.set(cls, cls, prev, clientDelay, 1.0);
            }
        }
        for (int si = 0; si < nsrv; si++) {
            for (int[] oa : openArrivalsOf.get(servers.get(si))) {
                JobClass cls = model.getClassByIndex(oa[0] - 1);
                P.set(cls, cls, sourceStation, srvStation[si], 1.0);
                P.set(cls, cls, srvStation[si], sinkStation, 1.0);
            }
        }
        model.link(P);

        int nclasses = model.getNumberOfClasses();
        for (int si = 0; si < nsrv; si++) {
            int idx = servers.get(si);
            PHLayer L = new PHLayer();
            L.idx = idx;
            L.ishost = idx < lqn.nhosts;
            L.callers = phToArray(callersOf.get(idx));
            L.classOfCaller = classOfCaller;
            L.nreplicas = 1;
            L.qstations = new int[]{stationOf[idx]};
            L.svcmeanByClass = new double[nclasses + 1];
            List<int[]> oa = openArrivalsOf.get(idx);
            L.openArrivals = oa.toArray(new int[0][]);
            L.npop = npop < 1 ? 1 : npop;
            phLayers[idx] = L;
        }
        phFlatModel = model;
        return servers;
    }

    /**
     * Processors and called tasks that become stations of the flat layer.
     *
     * The set is the elements the srvn builder would have given a layer of their
     * own, so 'flat.ph' and 'srvn.ph' place the SAME stations and differ only in
     * how many networks hold them. The refusals are those of flatServerSet, since
     * they are properties of the squashing and not of the encoding: each of these
     * carries per-layer state that one submodel cannot hold.
     *
     * @return the served element indices, processors first
     */
    private List<Integer> phFlatServerSet() {
        int nelem = lqn.nhosts + lqn.ntasks;
        for (int i = 0; i < nelem; i++) {
            if (lqn.repl.get(0, i) > 1) {
                throw new IllegalStateException("method='flat.ph' does not support replicated processors "
                        + "or tasks, whose replicas need a submodel each. Use method='srvn.ph'.");
            }
            if (lqn.iscache != null && i < lqn.iscache.getNumCols() && lqn.iscache.get(0, i) != 0) {
                throw new IllegalStateException(
                        "method='flat.ph' does not support cache tasks. Use method='default'.");
            }
            if (lqn.hassetup != null && i < lqn.hassetup.getNumCols() && lqn.hassetup.get(0, i) != 0) {
                throw new IllegalStateException("method='flat.ph' does not support setup tasks, whose "
                        + "powered-down threads are per-layer state. Use method='srvn.ph'.");
            }
        }
        List<Integer> servers = new ArrayList<Integer>();
        for (int hidx = 0; hidx < lqn.nhosts; hidx++) {
            if (this.ignore.get(hidx) != 0) {
                continue;
            }
            List<Integer> tasks = lqn.tasksof.get(hidx);
            if (tasks == null || tasks.isEmpty()) {
                continue;
            }
            if (phHostLayerCallers(hidx).isEmpty()) {
                continue;
            }
            servers.add(hidx);
        }
        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            if (this.ignore.get(tidx) != 0 || lqn.isref.get(tidx) != 0) {
                continue;
            }
            if (phTaskLayerCallers(tidx).isEmpty() && phAsyncCallsInto(tidx).isEmpty()) {
                continue;
            }
            servers.add(tidx);
        }
        if (servers.isEmpty()) {
            throw new IllegalStateException(
                    "method='flat.ph' found no server: the model has no processor with tasks.");
        }
        return servers;
    }

    /** Build the two-station phLayers of server element IDX. */
    private Network buildPHLayer(int idx, List<Integer> callers, boolean ishost) {
        Network model = new Network(lqn.hashnames.get(idx));
        model.setChecks(false);
        model.getAttribute().setClientIdx(1);
        model.getAttribute().setServerIdx(2);
        model.getAttribute().setSourceIdx(-1);

        Delay clientDelay = new Delay(model, "Clients");

        int nreplicas = phReplicaCount(idx, callers, ishost);
        Queue[] srvStation = new Queue[nreplicas + 1];
        for (int m = 1; m <= nreplicas; m++) {
            String nm = (m == 1) ? lqn.hashnames.get(idx) : (lqn.hashnames.get(idx) + "." + m);
            srvStation[m] = new Queue(model, nm, lqn.sched.get(idx));
            srvStation[m].setNumberOfServers((int) lqn.maxmult.get(0, idx));
            srvStation[m].getAttribute().setIsHost(ishost);
            srvStation[m].getAttribute().setIdx(idx);
        }
        int serverStationIdx = model.getNodeIndex(srvStation[1]) + 1;
        model.getAttribute().putServerIdxOf(idx, serverStationIdx);
        if (ishost) {
            model.getAttribute().getHostStations().add(serverStationIdx);
            model.getAttribute().addHosts(new Integer[]{null, serverStationIdx});
        } else {
            model.getAttribute().getTaskStations().add(serverStationIdx);
            model.getAttribute().addTasks(new Integer[]{null, serverStationIdx});
        }

        // --- closed class per caller task
        int[] classOfCaller = new int[lqn.nidx + 1];
        Arrays.fill(classOfCaller, -1);
        for (int c : callers) {
            double njobs = phLayerPopulation(idx, c, nreplicas);
            this.njobs.set(c, idx, njobs);
            ClosedClass cls = new ClosedClass(model, lqn.hashnames.get(c), njobs, clientDelay);
            cls.setReferenceClass(true);
            cls.setAttribute(new Integer[]{LayeredNetworkElement.TASK, c});
            classOfCaller[c] = cls.getIndex();
            model.getAttribute().addTasks(new Integer[]{cls.getIndex(), c});
            clientDelay.setService(cls, Exp.fitMean(FastMath.max(GlobalConstants.FineTol, refThinkTime(lqn, c))));
            for (int m = 1; m <= nreplicas; m++) {
                srvStation[m].setService(cls, Exp.fitMean(GlobalConstants.FineTol));
            }
            // every phLayers must be refreshed after a law change: post() resets the
            // layers named by the think-time map
            phAddUpdMap(this.cell_thinkt_classes_updmap, idx, new Integer[]{idx, c, 1, cls.getIndex()});
            phAddUpdMap(this.cell_servt_classes_updmap, idx,
                    new Integer[]{idx, c, serverStationIdx, cls.getIndex()});
        }

        // --- open classes: entry arrivals on a host phLayers, async calls on a task phLayers
        List<int[]> openArrivals = new ArrayList<int[]>();
        Source sourceStation = null;
        Sink sinkStation = null;
        if (ishost) {
            for (int c : callers) {
                // A task no other task calls has no task phLayers, so updateThinkTimesPH
                // never gives its caller class a surrogate delay: the class cycles against
                // an Immediate one and an open stream on top of it doubles the load. The
                // chain is the representation that honours the thread pool, so it is kept
                // and closed on the arrival rate instead -- see updateThinkTimesPH.
                if (phOpenArrivalOnly(c)) {
                    continue;
                }
                for (int eidx : phEntriesOf(c)) {
                    if (!phHasOpenArrival(eidx)) {
                        continue;
                    }
                    if (sourceStation == null) {
                        model.getAttribute().setSourceIdx(model.getNumberOfNodes() + 1);
                        sourceStation = new Source(model, "Source");
                        sinkStation = new Sink(model, "Sink");
                    }
                    OpenClass ocls = new OpenClass(model, lqn.hashnames.get(eidx) + ".Open", 0);
                    ocls.setAttribute(new Integer[]{LayeredNetworkElement.ENTRY, eidx});
                    sourceStation.setArrival(ocls, lqn.arrival.get(eidx));
                    clientDelay.setService(ocls, Disabled.getInstance());
                    for (int m = 1; m <= nreplicas; m++) {
                        srvStation[m].setService(ocls,
                                Exp.fitMean(FastMath.max(GlobalConstants.FineTol, phHostmean[eidx])));
                    }
                    openArrivals.add(new int[]{ocls.getIndex(), eidx});
                    model.getAttribute().addEntries(new Integer[]{ocls.getIndex(), eidx});
                    phAddUpdMap(this.cell_arvproc_classes_updmap, idx, new Integer[]{idx, -eidx,
                            model.getNodeIndex(sourceStation) + 1, ocls.getIndex()});
                }
            }
        } else {
            for (int cidx : phAsyncCallsInto(idx)) {
                if (sourceStation == null) {
                    model.getAttribute().setSourceIdx(model.getNumberOfNodes() + 1);
                    sourceStation = new Source(model, "Source");
                    sinkStation = new Sink(model, "Sink");
                }
                OpenClass ocls = new OpenClass(model, lqn.callhashnames.get(cidx), 0);
                ocls.setAttribute(new Integer[]{LayeredNetworkElement.CALL, cidx});
                sourceStation.setArrival(ocls, Immediate.getInstance());
                clientDelay.setService(ocls, Disabled.getInstance());
                int eidx = (int) lqn.callpair.get(cidx, 1);
                for (int m = 1; m <= nreplicas; m++) {
                    srvStation[m].setService(ocls,
                            Exp.fitMean(FastMath.max(GlobalConstants.FineTol, phEntrymean[eidx])));
                }
                openArrivals.add(new int[]{ocls.getIndex(), -cidx});
                model.getAttribute().addCalls(new Integer[]{ocls.getIndex(), cidx,
                        (int) lqn.callpair.get(cidx, 0), eidx});
                phAddUpdMap(this.cell_arvproc_classes_updmap, idx, new Integer[]{idx, cidx,
                        model.getNodeIndex(sourceStation) + 1, ocls.getIndex()});
                phAddUpdMap(this.cell_call_classes_updmap, idx,
                        new Integer[]{idx, cidx, serverStationIdx, ocls.getIndex()});
            }
        }

        if (sourceStation != null) {
            for (JobClass jc : model.getClasses()) {
                if (jc instanceof ClosedClass) {
                    sourceStation.setArrival(jc, Disabled.getInstance());
                }
            }
        }

        // Routing: one visit to the server per client cycle. The number of calls is
        // carried by the service law, not by a visit ratio, so no arc ever changes
        RoutingMatrix P = model.initRoutingMatrix();
        for (int c : callers) {
            JobClass cls = model.getClassByIndex(classOfCaller[c] - 1);
            for (int m = 1; m <= nreplicas; m++) {
                P.set(cls, cls, clientDelay, srvStation[m], 1.0 / nreplicas);
                P.set(cls, cls, srvStation[m], clientDelay, 1.0);
            }
        }
        for (int[] oa : openArrivals) {
            JobClass cls = model.getClassByIndex(oa[0] - 1);
            for (int m = 1; m <= nreplicas; m++) {
                P.set(cls, cls, sourceStation, srvStation[m], 1.0 / nreplicas);
                P.set(cls, cls, srvStation[m], sinkStation, 1.0);
            }
        }
        model.link(P);

        PHLayer L = new PHLayer();
        L.idx = idx;
        L.ishost = ishost;
        L.callers = phToArray(callers);
        L.classOfCaller = classOfCaller;
        L.nreplicas = nreplicas;
        L.qstations = new int[nreplicas];
        for (int m = 1; m <= nreplicas; m++) {
            L.qstations[m - 1] = 1 + m;
        }
        L.svcmeanByClass = new double[model.getNumberOfClasses() + 1];
        L.openArrivals = openArrivals.toArray(new int[0][]);
        double np = 0;
        for (int c : callers) {
            double v = this.njobs.get(c, idx);
            if (!Double.isNaN(v) && !Double.isInfinite(v) && v > 0) {
                np += v;
            }
        }
        L.npop = np < 1 ? 1 : np;
        phLayers[idx] = L;
        return model;
    }

    // =====================================================================
    // Composition of the entry laws
    // =====================================================================

    /** Build the per-entry workflows and the iteration-invariant processor law. */
    private void phInitLaws() {
        int n = lqn.nidx + 1;
        phWf = new Workflow[n];
        phWfhost = new Workflow[n];
        phExecs = new double[n][];
        phCallexecs = new double[n][];
        phHostalpha = new Matrix[n];
        phHostT = new Matrix[n];
        phHostmean = new double[n];
        phEntryalpha = new Matrix[n];
        phEntryT = new Matrix[n];
        phEntrymean = new double[n];
        phEntryscv = new double[n];
        phShare = new double[n];
        phOverlap = new double[n];
        phSetupshare = new double[n];
        phXdemand = new double[n];
        phNcalls = new double[n][n];
        phCalltime = new double[n][n];
        phProcresid = new double[n];
        phActthinkt = new double[n];
        phCalltotal = new double[n];
        phLayers = new PHLayer[lqn.nhosts + lqn.ntasks + 1];
        Arrays.fill(phEntryscv, 1.0);
        Arrays.fill(phOverlap, 1.0);

        for (int e = 0; e < lqn.nentries; e++) {
            int eidx = lqn.eshift + e;
            int tidx = (int) lqn.parent.get(0, eidx);
            if (this.ignore.get(tidx) != 0) {
                continue;
            }
            LqnPh.EntryWorkflow ew = LqnPh.entryWorkflow(lqnModel, lqn, eidx, true);
            phWf[eidx] = ew.wf;
            phExecs[eidx] = ew.execs;
            phCallexecs[eidx] = ew.callexecs;
            LqnPh.EntryWorkflow ehost = LqnPh.entryWorkflow(lqnModel, lqn, eidx, false);
            phWfhost[eidx] = ehost.wf;
            // the processor sees the WORK of concurrent branches, not their elapsed
            // time, so the host law serialises an AND fork -- see LqnPh.serialLaw
            Pair<Matrix, Matrix> hl = LqnPh.serialLaw(phWfhost[eidx]);
            phHostalpha[eidx] = hl.getLeft();
            phHostT[eidx] = hl.getRight();
            phHostmean[eidx] = LqnPh.moments(phHostalpha[eidx], phHostT[eidx])[0];
        }

        // until the first iteration reports throughputs, a task splits its
        // requests evenly over its entries
        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            List<Integer> entries = phEntriesOf(tidx);
            if (entries.isEmpty()) {
                continue;
            }
            for (int eidx : entries) {
                phShare[eidx] = 1.0 / entries.size();
            }
        }
    }

    /**
     * Recompose the entry service laws from the current fixed-point iterate.
     * <p>
     * The composed mean is NOT the sum of the leaf means when the graph forks:
     * the branches of an AND fork phOverlap, and the entry finishes with the last
     * of them. The ratio of the two, the phOverlap factor, is what the
     * caller-side aggregates are scaled by, so that the pieces of a cycle still
     * add up to the cycle.
     * </p>
     */
    public void phComposeEntryLaws() {
        double[] entrySetupShare = new double[lqn.nidx + 1];
        Arrays.fill(phOverlap, 1.0);

        for (int e = 0; e < lqn.nentries; e++) {
            int eidx = lqn.eshift + e;
            if (phWf[eidx] == null) {
                continue;
            }
            Workflow w = phWf[eidx];
            double[] ex = phExecs[eidx];
            double entrysum = 0;
            double procsum = 0;
            for (int aidx : phActsOf(eidx)) {
                double m = this.residt.get(aidx) + actThinkTime(lqn, aidx);
                procsum += ex[aidx] * m;
                w.setActivityDemandMean(lqn.names.get(aidx), FastMath.max(m, GlobalConstants.FineTol));
                for (int cidx : phCallsOf(aidx)) {
                    if (lqn.calltype.get(cidx) != CallType.SYNC) {
                        continue;
                    }
                    w.setActivityDemand(lqn.callhashnames.get(cidx), phCallBurstLaw(cidx));
                    m += this.callservt.get(cidx);
                }
                entrysum += ex[aidx] * m;
            }
            Markovian law = w.refreshPH();
            Matrix alpha = law.getInitProb();
            Matrix T = law.getSubgenerator();
            double[] mm = LqnPh.moments(alpha, T);
            double m1 = mm[0];
            double scv = mm[1];
            // All activities of an entry run on ONE processor, so the branches of an
            // AND fork cannot phOverlap the processor phResidence they request: the
            // composed maximum is a lower bound on the entry service time only above
            // that total. Where it falls below, the law is rescaled in time to it,
            // which keeps its shape, its SCV and its order.
            if (procsum > m1 + GlobalConstants.FineTol) {
                T = T.scale(m1 / procsum); // scale() returns a new matrix, it does not mutate
                m1 = procsum;
            }
            // A SetupTask powers a thread down when it goes idle, so a request may
            // find it off and pay a cold start before the entry runs at all. The
            // setup is not part of the activity graph and never enters the
            // series-parallel reduction: it is prefixed to the composed law
            // afterwards, as the mixture p*(setup THEN entry) + (1-p)*entry, which
            // is again phase-type. See phSetupProb for p.
            double p = phSetupProb(eidx);
            if (p > GlobalConstants.FineTol) {
                Pair<Matrix, Matrix> sl = phSetupLaw((int) lqn.parent.get(0, eidx));
                if (sl != null) {
                    Pair<Matrix, Matrix> c = Workflow.composeSerial(sl.getLeft(), sl.getRight(), alpha, T);
                    List<Matrix> alphas = new ArrayList<Matrix>();
                    List<Matrix> Ts = new ArrayList<Matrix>();
                    alphas.add(c.getLeft());
                    Ts.add(c.getRight());
                    alphas.add(alpha);
                    Ts.add(T);
                    Pair<Matrix, Matrix> mix = Workflow.composeMixture(alphas, Ts, new double[]{p, 1 - p});
                    alpha = mix.getLeft();
                    T = mix.getRight();
                    mm = LqnPh.moments(alpha, T);
                    m1 = mm[0];
                    scv = mm[1];
                    // The phShare of the entry law that is cold start and not work. The
                    // surrogate-delay closure measures a thread's cycle in WORK, so it
                    // must not read a station utilization that this has inflated --
                    // see updateThinkTimesPH.
                    entrySetupShare[eidx] = p * setupMeanOf(lqn.setuptime, (int) lqn.parent.get(0, eidx))
                            / FastMath.max(m1, GlobalConstants.FineTol);
                }
            }
            phEntryalpha[eidx] = alpha;
            phEntryT[eidx] = T;
            phEntrymean[eidx] = m1;
            phEntryscv[eidx] = scv;
            if (entrysum > GlobalConstants.FineTol) {
                phOverlap[eidx] = FastMath.min(1, m1 / entrysum);
            }
        }

        // Per task, the phShare-weighted fraction of its station service that is
        // cold start rather than work.
        Arrays.fill(phSetupshare, 0.0);
        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            if (this.ignore.get(tidx) != 0) {
                continue;
            }
            for (int eidx : phEntriesOf(tidx)) {
                phSetupshare[tidx] += phShare[eidx] * entrySetupShare[eidx];
            }
        }

        // Expected number of calls per invocation, and the caller-side aggregates
        for (int i = 0; i < lqn.nidx; i++) {
            Arrays.fill(phNcalls[i], 0.0);
            Arrays.fill(phCalltime[i], 0.0);
        }
        Arrays.fill(phProcresid, 0.0);
        Arrays.fill(phActthinkt, 0.0);
        Arrays.fill(phCalltotal, 0.0);

        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            if (this.ignore.get(tidx) != 0) {
                continue;
            }
            for (int eidx : phEntriesOf(tidx)) {
                double w = phShare[eidx];
                if (w <= 0 || phExecs[eidx] == null) {
                    continue;
                }
                double[] ex = phExecs[eidx];
                double r = phOverlap[eidx];
                for (int aidx : phActsOf(eidx)) {
                    phProcresid[tidx] += w * r * ex[aidx] * this.residt.get(aidx);
                    phActthinkt[tidx] += w * r * ex[aidx] * actThinkTime(lqn, aidx);
                    for (int cidx : phCallsOf(aidx)) {
                        if (lqn.calltype.get(cidx) != CallType.SYNC) {
                            continue;
                        }
                        int tgte = (int) lqn.callpair.get(cidx, 1);
                        int tgtt = (int) lqn.parent.get(0, tgte);
                        // the COUNT of calls does not change with the phOverlap, only
                        // the time the caller is held by them
                        phNcalls[tidx][tgte] += w * ex[aidx] * lqn.callproc_mean.getOrDefault(cidx, 0.0);
                        phCalltime[tidx][tgtt] += w * r * ex[aidx] * this.callservt.get(cidx);
                        phCalltotal[tidx] += w * r * ex[aidx] * this.callservt.get(cidx);
                    }
                }
            }
        }
    }

    /**
     * Law of the total time one execution of the issuing activity spends in call
     * CIDX: the geometric compound, of mean callproc_mean, of the response law of
     * the called entry. The response law is fitted to the response time reported
     * by the callee's phLayers and to the SCV of the callee's own composed law, so
     * no extra solver output is needed.
     */
    private Distribution phCallBurstLaw(int cidx) {
        double m = lqn.callproc_mean.getOrDefault(cidx, 0.0);
        int eidx = (int) lqn.callpair.get(cidx, 1);
        if (m <= GlobalConstants.FineTol) {
            return Immediate.getInstance();
        }
        double R = this.callservt.get(cidx) / m;
        double scv = phEntryscv[eidx];
        if (Double.isNaN(scv) || Double.isInfinite(scv) || scv <= GlobalConstants.FineTol) {
            scv = 1.0;
        }
        APH base = APH.fitMeanAndSCV(FastMath.max(R, GlobalConstants.FineTol), scv);
        Pair<Matrix, Matrix> loop = Workflow.composeLoopGeometric(base.getInitProb(), base.getSubgenerator(), m);
        if (Workflow.isAcyclicGenerator(loop.getRight())) {
            return new APH(loop.getLeft(), loop.getRight());
        }
        return new jline.lang.processes.PH(loop.getLeft(), loop.getRight());
    }

    // =====================================================================
    // Pushing the composed laws into the layers
    // =====================================================================

    /**
     * Push the composed laws into the layers.
     * <p>
     * A phLayers of this method carries no routing that depends on the iterate: the
     * number of calls a caller makes is folded into its service law rather than
     * into a visit ratio, so only two laws move per (phLayers, class) -- the
     * phase-type service law at the server and the mean of the surrogate delay
     * at the client.
     * </p>
     *
     * @param it iteration number
     */
    public void updateLayersPH(int it) {
        for (int idx = 0; idx < lqn.nhosts + lqn.ntasks; idx++) {
            if (Double.isNaN(this.idxhash.get(idx)) || phLayers[idx] == null) {
                continue;
            }
            PHLayer L = phLayers[idx];
            Network model = this.ensemble[this.idxhash.get(idx).intValue() - 1];
            Delay clientDelay = (Delay) model.getStations().get(0);

            for (int c : L.callers) {
                int k = L.classOfCaller[c];
                JobClass cls = model.getClassByIndex(k - 1);
                Pair<Matrix, Matrix> sl = phServiceLaw(idx, L.ishost, c);
                L.svcmeanByClass[k] = LqnPh.moments(sl.getLeft(), sl.getRight())[0];
                Distribution law = phStationLaw(sl.getLeft(), sl.getRight());
                for (int s = 0; s < L.nreplicas; s++) {
                    ((Queue) model.getStations().get(L.qstations[s] - 1)).setService(cls, law);
                }
                clientDelay.setService(cls, Exp.fitMean(
                        FastMath.max(GlobalConstants.FineTol, phDelayMean(idx, c))));
            }

            for (int[] oa : L.openArrivals) {
                int k = oa[0];
                int tag = oa[1];
                JobClass cls = model.getClassByIndex(k - 1);
                if (tag > 0) {
                    // entry arrival: the processor demand law of the entry is static
                    L.svcmeanByClass[k] = phHostmean[tag];
                    continue;
                }
                int cidx = -tag;
                int eidx = (int) lqn.callpair.get(cidx, 1);
                L.svcmeanByClass[k] = phEntrymean[eidx];
                Distribution law = phStationLaw(phEntryalpha[eidx], phEntryT[eidx]);
                for (int s = 0; s < L.nreplicas; s++) {
                    ((Queue) model.getStations().get(L.qstations[s] - 1)).setService(cls, law);
                }
                int aidx = (int) lqn.callpair.get(cidx, 0);
                double rate = this.tput.get(aidx) * lqn.callproc_mean.getOrDefault(cidx, 0.0);
                if (Double.isNaN(rate) || Double.isInfinite(rate) || rate <= GlobalConstants.FineTol) {
                    rate = GlobalConstants.FineTol;
                }
                ((Source) model.getNodes().get(model.getAttribute().getSourceIdx() - 1))
                        .setArrival(cls, Exp.fitRate(rate));
            }
        }
    }

    /** Law of the demand caller C places on the server of phLayers IDX per invocation. */
    private Pair<Matrix, Matrix> phServiceLaw(int idx, boolean ishost, int c) {
        if (ishost) {
            // mixture over the entries of C, weighted by their phShare of its requests
            List<Matrix> alphas = new ArrayList<Matrix>();
            List<Matrix> Ts = new ArrayList<Matrix>();
            List<Double> probs = new ArrayList<Double>();
            for (int eidx : phEntriesOf(c)) {
                if (phHostT[eidx] == null || phShare[eidx] <= 0) {
                    continue;
                }
                alphas.add(phHostalpha[eidx]);
                Ts.add(phHostT[eidx]);
                probs.add(phShare[eidx]);
            }
            if (alphas.isEmpty()) {
                return phImmediateLaw();
            }
            double tot = 0;
            for (double p : probs) {
                tot += p;
            }
            double[] pr = new double[probs.size()];
            for (int i = 0; i < probs.size(); i++) {
                pr[i] = probs.get(i) / tot;
            }
            return Workflow.composeMixture(alphas, Ts, pr);
        }

        // task phLayers: the total demand is the sum, over the entries of the server, of
        // a geometric compound of the entry law of mean equal to the number of calls
        Matrix alpha = null;
        Matrix T = null;
        for (int eidx : phEntriesOf(idx)) {
            double n = phNcalls[c][eidx];
            if (n <= GlobalConstants.FineTol || phEntryT[eidx] == null) {
                continue;
            }
            Pair<Matrix, Matrix> lp = Workflow.composeLoopGeometric(phEntryalpha[eidx], phEntryT[eidx], n);
            if (alpha == null) {
                alpha = lp.getLeft();
                T = lp.getRight();
            } else {
                Pair<Matrix, Matrix> s = Workflow.composeSerial(alpha, T, lp.getLeft(), lp.getRight());
                alpha = s.getLeft();
                T = s.getRight();
            }
        }
        if (alpha == null) {
            return phImmediateLaw();
        }
        return new Pair<Matrix, Matrix>(alpha, T);
    }

    private Pair<Matrix, Matrix> phImmediateLaw() {
        Matrix a = Matrix.singleton(1.0);
        Matrix T = Matrix.singleton(-GlobalConstants.Immediate);
        return new Pair<Matrix, Matrix>(a, T);
    }

    /**
     * Mean time a thread of caller C spends away from the server of phLayers IDX per
     * invocation: idle, plus whatever of its cycle the phLayers does not hold.
     */
    private double phDelayMean(int idx, int c) {
        double z = this.thinkt.get(c) + refThinkTime(lqn, c);
        if (Double.isNaN(z) || Double.isInfinite(z) || z < 0) {
            z = 0;
        }
        z += phActthinkt[c];
        // The caller's own processor residence, unless this model holds that
        // processor as a station of its own.
        int hidx = (int) lqn.parent.get(c);
        if (!phServedHere(idx, hidx)) {
            z += phProcresid[c];
        }
        // And the time spent at every callee this model does not hold. Every
        // term is SUMMED in rather than obtained by subtracting from a total.
        // That subtraction cancels catastrophically once the call time is large:
        // a caller whose only callee is this server has the two terms equal, and
        // 7 + 1.4e47 - 1.4e47 is 0, not 7, because the think time falls below the
        // ULP of the call time. The layer then sees a client delay of zero,
        // saturates, reports a residence time that inflates the very call time
        // that caused the cancellation, and the fixed point runs away --
        // lqn_sockshop reached RespT 1.4e47 this way.
        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            if (!phServedHere(idx, tidx)) {
                z += phCalltime[c][tidx];
            }
        }
        if (Double.isNaN(z) || Double.isInfinite(z) || z < 0) {
            z = GlobalConstants.FineTol;
        }
        return z;
    }

    /**
     * True when LQN element ELEM is a station of the same model that holds
     * server IDX. Under 'srvn.ph' that is ELEM == IDX, since each server has a
     * layer of its own; under 'flat.ph' it is every server of the one network.
     */
    private boolean phServedHere(int idx, int elem) {
        if (elem < 0 || elem >= this.idxhash.size() || idx < 0 || idx >= this.idxhash.size()) {
            return false;
        }
        double a = this.idxhash.get(elem);
        double b = this.idxhash.get(idx);
        return !Double.isNaN(a) && !Double.isNaN(b) && a == b;
    }

    /**
     * Station law of a composed workflow. A geometric loop over a body of two or
     * more phases closes a cycle in the phase graph, and a cyclic generator is a
     * PH and not an APH: no phLayers solver declares PH, so such a law is reduced to
     * the APH with the SAME first two moments. AMVA and NC read exactly those
     * two, so the reduction is lossless for them and is a two-moment fit for the
     * phase-aware phLayers solvers.
     */
    private Distribution phStationLaw(Matrix alpha, Matrix T) {
        if (Workflow.isAcyclicGenerator(T)) {
            return new APH(alpha, T);
        }
        double[] mm = LqnPh.moments(alpha, T);
        return APH.fitMeanAndSCV(mm[0], mm[1]);
    }

    // =====================================================================
    // Metric reconstruction
    // =====================================================================

    /**
     * Reconstruct the LQN metrics.
     * <p>
     * A phLayers of this method reports one row per caller task, not one per entry,
     * activity and call, so the per-element quantities the rest of SolverLN reads
     * -- servt, residt, callservt, callresidt, tput -- are recovered analytically
     * from the series-parallel weights of the entry workflows.
     * </p>
     * <p>
     * The split is conservative by construction. A station reports a phResidence
     * time R per visit against a service law of mean S, so the queueing inflation
     * R/S is attributed to every leaf of that visit in proportion to its own
     * mean: the pieces sum back to R exactly.
     * </p>
     *
     * @param it iteration number
     */
    public void updateMetricsPH(int it) {
        int nidx = lqn.nidx;
        this.servt = new Matrix(1, nidx, nidx);
        this.residt = new Matrix(1, nidx, nidx);
        this.callservt = new Matrix(1, FastMath.max(lqn.ncalls, 1), FastMath.max(lqn.ncalls, 1));
        this.callresidt = new Matrix(1, FastMath.max(lqn.ncalls, 1), FastMath.max(lqn.ncalls, 1));

        double[] inflNum = new double[nidx + 1];
        double[] inflDen = new double[nidx + 1];
        double[] taskTput = new double[nidx + 1];
        double[] openTput = new double[nidx + 1];

        // Host layers: the queueing inflation of the processor demand
        for (int hidx = 0; hidx < lqn.nhosts; hidx++) {
            if (Double.isNaN(this.idxhash.get(hidx)) || phLayers[hidx] == null) {
                continue;
            }
            PHLayer L = phLayers[hidx];
            SolverResult res = phLastResult(this.idxhash.get(hidx).intValue() - 1);
            double npop = phLayerPop(L, hidx);
            for (int c : L.callers) {
                int k = L.classOfCaller[c];
                int kc = k - 1; // result matrices index classes from zero
                double X = phSumOver(res.TN, L.qstations, kc);
                double R = phResidence(phSumOver(res.QN, L.qstations, kc), X, res.RN.get(L.qstations[0] - 1, kc));
                double f = phInflationOf(R, L.svcmeanByClass[k], npop);
                if (Double.isNaN(X) || Double.isInfinite(X) || X < 0) {
                    X = 0;
                }
                // TOTAL over the replicas. The processor phLayers of a replicated element
                // models ONE representative replica, so X is one replica's rate and the
                // element's own rate is REPL times it. The matching per replica quantity
                // is phXdemand, which the think-time closure divides down for the same reason.
                taskTput[c] += lqn.repl.get(0, c) * X;
                for (int eidx : phEntriesOf(c)) {
                    double w = FastMath.max(phShare[eidx], 0) * X;
                    inflNum[eidx] += w * f;
                    inflDen[eidx] += w;
                }
            }
            for (int[] oa : L.openArrivals) {
                int tag = oa[1];
                if (tag <= 0) {
                    continue; // an async call is served in the task phLayers, not here
                }
                int k = oa[0];
                int kc = k - 1;
                int eidx = tag;
                double X = phSumOver(res.TN, L.qstations, kc);
                if (Double.isNaN(X) || Double.isInfinite(X) || X <= 0) {
                    continue;
                }
                double f = phInflationOf(phResidence(phSumOver(res.QN, L.qstations, kc), X,
                        res.RN.get(L.qstations[0] - 1, kc)), L.svcmeanByClass[k], npop);
                inflNum[eidx] += X * f;
                inflDen[eidx] += X;
                openTput[eidx] += X;
                taskTput[(int) lqn.parent.get(0, eidx)] += X;
            }
        }

        for (int e = 0; e < lqn.nentries; e++) {
            int eidx = lqn.eshift + e;
            double f = 1;
            if (inflDen[eidx] > GlobalConstants.FineTol) {
                f = inflNum[eidx] / inflDen[eidx];
            }
            if (Double.isNaN(f) || Double.isInfinite(f) || f < 1) {
                f = 1; // a phResidence time cannot fall below the demand it contains
            }
            for (int aidx : phActsOf(eidx)) {
                this.residt.set(aidx, f * lqn.hostdem_mean.getOrDefault(aidx, 0.0));
            }
        }

        // Task layers: the response time of every call
        double[] relw = new double[nidx + 1];
        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            if (Double.isNaN(this.idxhash.get(tidx)) || phLayers[tidx] == null) {
                continue;
            }
            PHLayer L = phLayers[tidx];
            SolverResult res = phLastResult(this.idxhash.get(tidx).intValue() - 1);
            double npop = phLayerPop(L, tidx);
            for (int c : L.callers) {
                int k = L.classOfCaller[c];
                int kc = k - 1;
                double X = phSumOver(res.TN, L.qstations, kc);
                if (Double.isNaN(X) || Double.isInfinite(X) || X < 0) {
                    X = 0;
                }
                double g = phInflationOf(phResidence(phSumOver(res.QN, L.qstations, kc), X,
                        res.RN.get(L.qstations[0] - 1, kc)), L.svcmeanByClass[k], npop);
                for (int cidx : phSyncCallsBetween(c, tidx)) {
                    int eidx = (int) lqn.callpair.get(cidx, 1);
                    double v = lqn.callproc_mean.getOrDefault(cidx, 0.0) * g * phEntrymean[eidx];
                    this.callservt.set(cidx, v);
                    this.callresidt.set(cidx, v);
                }
                for (int eidx : phEntriesOf(tidx)) {
                    relw[eidx] += X * phNcalls[c][eidx];
                }
            }
            for (int[] oa : L.openArrivals) {
                int tag = oa[1];
                if (tag >= 0) {
                    continue;
                }
                int kc = oa[0] - 1;
                int cidx = -tag;
                int eidx = (int) lqn.callpair.get(cidx, 1);
                double X = phSumOver(res.TN, L.qstations, kc);
                double R = res.RN.get(L.qstations[0] - 1, kc);
                if (!Double.isNaN(R) && !Double.isInfinite(R) && R > 0) {
                    double v = R * lqn.callproc_mean.getOrDefault(cidx, 0.0);
                    this.callservt.set(cidx, v);
                    this.callresidt.set(cidx, v);
                }
                if (!Double.isNaN(X) && !Double.isInfinite(X) && X > 0) {
                    relw[eidx] += X;
                }
            }
        }

        // Entry shares and throughputs
        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            List<Integer> entries = phEntriesOf(tidx);
            if (entries.isEmpty()) {
                continue;
            }
            // How the requests SPLIT over the entries is a flow-balance question, and
            // is answered at the task phLayers: a caller class reaches that server once
            // per invocation of the caller, carrying its whole call burst in its
            // service law, so the station rate counts caller cycles and the per-entry
            // rate is that rate times the calls the caller makes.
            double tot = 0;
            for (int eidx : entries) {
                tot += relw[eidx] + openTput[eidx];
            }
            if (tot > GlobalConstants.FineTol) {
                for (int eidx : entries) {
                    phShare[eidx] = (relw[eidx] + openTput[eidx]) / tot;
                }
            } else {
                for (int eidx : entries) {
                    phShare[eidx] = 1.0 / entries.size();
                }
            }
            // HOW MANY requests the task completes is a different question, and the
            // flow-balance total does not answer it: that total is what the callers
            // DEMAND, not what the task's threads can deliver. A thread cycles through
            // its host demand AND then through the task think time, and only the
            // processor phLayers of the task carries both, so the rate is read there.
            if (taskTput[tidx] > GlobalConstants.FineTol) {
                this.tput.set(tidx, taskTput[tidx]);
            } else {
                this.tput.set(tidx, tot); // no processor phLayers of its own
            }
            for (int eidx : entries) {
                this.tput.set(eidx, this.tput.get(tidx) * phShare[eidx]);
            }
            // The DEMAND is kept apart because it, and not the rate just reported, is
            // what closes the surrogate delay: normalising the think time by a rate the
            // same think time produced makes the processor phLayers self-referential and it
            // settles wherever it started -- see updateThinkTimesPH. PER REPLICA,
            // because the thread count it is paired with there is per replica.
            double nrep = FastMath.max(1, lqn.repl.get(0, tidx));
            if (tot > GlobalConstants.FineTol) {
                phXdemand[tidx] = tot / nrep;
            } else {
                phXdemand[tidx] = this.tput.get(tidx) / nrep;
            }
        }

        // Recovery, under-relaxation, and the derived per-element quantities
        double omega = this.relax_omega;
        for (int aidx = lqn.ashift; aidx < lqn.ashift + lqn.nacts; aidx++) {
            double v = this.residt.get(aidx);
            if ((Double.isInfinite(v) || Double.isNaN(v)) && it > 1
                    && !Double.isNaN(this.residt_prev.get(aidx))) {
                v = this.residt_prev.get(aidx);
            }
            if (omega < 1.0 && it > 1 && !Double.isNaN(this.residt_prev.get(aidx))) {
                v = omega * v + (1 - omega) * this.residt_prev.get(aidx);
            }
            this.residt.set(aidx, v);
            this.residt_prev.set(aidx, v);
        }
        for (int cidx = 0; cidx < lqn.ncalls; cidx++) {
            double v = this.callservt.get(cidx);
            if (Double.isInfinite(v) || Double.isNaN(v)) {
                if (it > 1 && !Double.isNaN(this.callservt_prev.get(cidx))
                        && !Double.isInfinite(this.callservt_prev.get(cidx))) {
                    v = this.callservt_prev.get(cidx);
                } else {
                    v = 0;
                }
            }
            if (omega < 1.0 && it > 1 && !Double.isNaN(this.callservt_prev.get(cidx))) {
                v = omega * v + (1 - omega) * this.callservt_prev.get(cidx);
            }
            this.callservt.set(cidx, v);
            this.callresidt.set(cidx, v);
            this.callservt_prev.set(cidx, v);
            this.callresidt_prev.set(cidx, v);
            if (v > 0) {
                this.callservtproc.put(cidx, Exp.fitMean(v));
            }
        }

        // Recompose the entry laws from the iterate just computed. The entry service
        // time is then the mean of the COMPOSED law and not the sum of the parts: the
        // branches of an AND fork phOverlap, so an entry that forks finishes with the
        // last of its branches and is not charged their sum.
        phComposeEntryLaws();

        for (int e = 0; e < lqn.nentries; e++) {
            int eidx = lqn.eshift + e;
            if (phExecs[eidx] == null) {
                continue;
            }
            double[] ex = phExecs[eidx];
            for (int aidx : phActsOf(eidx)) {
                double sa = this.residt.get(aidx) + actThinkTime(lqn, aidx);
                for (int cidx : phCallsOf(aidx)) {
                    if (lqn.calltype.get(cidx) == CallType.SYNC) {
                        sa += this.callservt.get(cidx);
                    }
                }
                this.servt.set(aidx, sa);
                this.servt_prev.set(aidx, sa);
                this.tput.set(aidx, this.tput.get(eidx) * ex[aidx]);
                this.tput_prev.set(aidx, this.tput.get(aidx));
                // fitRate clamps a null rate, but a never-called activity has no arrivals: Disabled says so, as the python twin does.
                this.tputproc.put(aidx, this.tput.get(aidx) > 0
                        ? Exp.fitRate(this.tput.get(aidx)) : Disabled.getInstance());
                if (sa > 0) {
                    this.servtproc.put(aidx, Exp.fitMean(sa));
                }
            }
            this.servt.set(eidx, phEntrymean[eidx]);
            this.residt.set(eidx, phEntrymean[eidx]);
            if (this.servt.get(eidx) > 0) {
                this.servtproc.put(eidx, Exp.fitMean(this.servt.get(eidx)));
            }
        }
    }

    // =====================================================================
    // Surrogate delays
    // =====================================================================

    /**
     * Surrogate delay of every caller.
     * <p>
     * Same closure as updateThinkTimes -- a thread of the task is idle for
     * whatever of its cycle the task's own station does not hold -- but the rate
     * it is normalised by is the INVOCATION rate of the task and not the
     * throughput of its station. Under this method a caller class reaches the
     * server once per invocation of the caller, carrying its whole call burst in
     * its service law, so the station rate counts caller cycles rather than calls
     * and the two differ by the mean number of calls.
     * </p>
     *
     * @param it iteration number
     */
    public void updateThinkTimesPH(int it) {
        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            if (this.ignore.get(tidx) != 0) {
                continue;
            }
            // only a reference task's think time separates one request from the next;
            // on a served task it is not a per-request delay -- see refThinkTime
            double ztask = refThinkTime(lqn, tidx);
            if (Double.isNaN(this.idxhash.get(tidx))) {
                // A task no other task calls but whose entries carry an arrival still
                // has a cycle: its threads are driven by the stream. buildLayersPH drops
                // the open class for it precisely so this closure can set the rate.
                double arvrate = phArrivalRate(tidx);
                if (arvrate > GlobalConstants.FineTol) {
                    double njobs = lqn.maxmult.get(0, tidx);
                    if (Double.isInfinite(njobs) || Double.isNaN(njobs) || njobs <= 0) {
                        njobs = phMaxNjobs(tidx);
                    }
                    double z = FastMath.max(GlobalConstants.Zero,
                            njobs / arvrate - phHostResid(tidx) - ztask);
                    double om = this.relax_omega;
                    if (om < 1.0 && it > 1 && !Double.isNaN(this.thinkt_prev.get(tidx))) {
                        z = om * z + (1 - om) * this.thinkt_prev.get(tidx);
                    }
                    this.tput.set(tidx, arvrate);
                    this.thinkt.set(tidx, z);
                    this.thinkt_prev.set(tidx, z);
                    this.thinktproc.put(tidx, Exp.fitMean(z + ztask));
                    continue;
                }
                // a reference task, or one no other task calls: it has no station of
                // its own, so its only delay is the think time the user declared
                this.thinkt.set(tidx, GlobalConstants.FineTol);
                this.thinktproc.put(tidx, Immediate.getInstance());
                continue;
            }
            int layer_t = this.idxhash.get(tidx).intValue() - 1;
            PHLayer L = phLayers[tidx];
            SolverResult res = phLastResult(layer_t);
            double U = 0;
            for (int k = 0; k < res.UN.getNumCols(); k++) {
                double v = res.UN.get(L.qstations[0] - 1, k);
                if (!Double.isNaN(v)) {
                    U += v;
                }
            }
            this.util.set(tidx, U);
            // The closure below measures a thread's cycle in WORK: it is idle for
            // whatever of the cycle its station does not hold it working. A SetupTask's
            // station service also carries a cold start, which is time the thread is
            // unavailable but is not work, so it is taken back out of U before the
            // closure reads it. Zero for every task without a setup.
            if (phSetupshare[tidx] > 0) {
                U = U * (1 - phSetupshare[tidx]);
            }
            // the rate the CALLERS ask of the task, not the rate its processor phLayers
            // reported: the latter is itself a function of this think time
            double X = phXdemand[tidx];
            if (!(X > GlobalConstants.FineTol)) {
                X = this.tput.get(tidx);
            }
            // The thread pool of ONE replica, the convention phXdemand is kept in.
            double njobs = lqn.maxmult.get(0, tidx);
            if (Double.isInfinite(njobs) || Double.isNaN(njobs) || njobs <= 0) {
                njobs = phMaxNjobs(tidx);
            }
            double z;
            if (X > GlobalConstants.FineTol) {
                if (lqn.sched.get(tidx) == SchedStrategy.INF) {
                    // an infinite server reports a mean number of busy threads
                    z = (njobs - U) / X - ztask;
                } else {
                    z = njobs * FastMath.abs(1 - U) / X - ztask;
                }
            } else {
                z = this.thinkt.get(tidx);
            }
            z = FastMath.max(GlobalConstants.Zero, z);
            if (it > 1 && !Double.isNaN(this.thinkt_prev.get(tidx))
                    && (Double.isInfinite(z) || Double.isNaN(z))) {
                z = this.thinkt_prev.get(tidx);
            }
            double omega = this.relax_omega;
            if (omega < 1.0 && it > 1 && !Double.isNaN(this.thinkt_prev.get(tidx))) {
                z = omega * z + (1 - omega) * this.thinkt_prev.get(tidx);
            }
            this.thinkt.set(tidx, z);
            this.thinkt_prev.set(tidx, z);
            this.thinktproc.put(tidx, Exp.fitMean(z + ztask));
        }
    }

    /**
     * Total exogenous rate into the entries of task TIDX, zero unless the arrival
     * is the only way in -- the predicate buildLayersPH drops the open class on.
     */
    private double phArrivalRate(int tidx) {
        if (lqn.isref.get(tidx) != 0) {
            return 0;
        }
        for (int eidx : phEntriesOf(tidx)) {
            if (phAnyCallerOf(eidx)) {
                return 0;
            }
        }
        if (lqn.arrival == null) {
            return 0;
        }
        double rate = 0;
        for (int eidx : phEntriesOf(tidx)) {
            Distribution d = lqn.arrival.get(eidx);
            if (d != null) {
                double m = d.getMean();
                if (!Double.isNaN(m) && !Double.isInfinite(m) && m > GlobalConstants.FineTol) {
                    rate += 1 / m;
                }
            }
        }
        return rate;
    }

    /** Response time the caller class of task TIDX sees at its processor phLayers. */
    private double phHostResid(int tidx) {
        int hidx = (int) lqn.parent.get(0, tidx);
        if (hidx < 0 || hidx >= this.idxhash.size() || Double.isNaN(this.idxhash.get(hidx))) {
            return 0;
        }
        PHLayer L = phLayers[hidx];
        if (L == null || L.classOfCaller[tidx] < 0) {
            return 0;
        }
        SolverResult res = phLastResult(this.idxhash.get(hidx).intValue() - 1);
        double r = res.RN.get(L.qstations[0] - 1, L.classOfCaller[tidx] - 1);
        return Double.isNaN(r) ? 0 : r;
    }

    // =====================================================================
    // Result reconstruction
    // =====================================================================

    /**
     * LQN-level results. The layers report per caller task, so every entry,
     * activity and call figure is rebuilt from the converged fixed point rather
     * than read off a class row, in the same layout getEnsembleAvg returns.
     *
     * @return {QN, UN, RN, TN, AN, WN}, each a 1 x (nidx+1) row indexed by
     *         element index
     */
    public Matrix[] getEnsembleAvgPH() {
        int nidx = lqn.nidx;
        Matrix QN = phNanRow(nidx);
        Matrix UN = phNanRow(nidx);
        Matrix RN = phNanRow(nidx);
        Matrix TN = phNanRow(nidx);
        Matrix AN = phNanRow(nidx);
        Matrix WN = phNanRow(nidx);
        Matrix PN = phNanRow(nidx); // processor utilization
        Matrix UT = phNanRow(nidx); // task and entry utilization

        for (int a = 0; a < lqn.nacts; a++) {
            int aidx = lqn.ashift + a;
            int tidx = (int) lqn.parent.get(0, aidx);
            if (this.ignore.get(tidx) != 0) {
                continue;
            }
            int hidx = (int) lqn.parent.get(0, tidx);
            TN.set(aidx, this.tput.get(aidx));
            RN.set(aidx, this.servt.get(aidx));
            UT.set(aidx, this.tput.get(aidx) * this.servt.get(aidx));
            // LINE scales the utilization of a queueing station into [0,1] whatever its
            // multiplicity, and reports a mean number of busy servers at an infinite
            // server: the processor phShare of an activity follows the same convention
            PN.set(aidx, this.tput.get(aidx) * lqn.hostdem_mean.getOrDefault(aidx, 0.0)
                    / phHostServers(hidx));
            if (Double.isNaN(PN.get(hidx))) {
                PN.set(hidx, 0);
            }
            PN.set(hidx, PN.get(hidx) + PN.get(aidx));
        }

        for (int e = 0; e < lqn.nentries; e++) {
            int eidx = lqn.eshift + e;
            int tidx = (int) lqn.parent.get(0, eidx);
            if (this.ignore.get(tidx) != 0) {
                continue;
            }
            TN.set(eidx, this.tput.get(eidx));
            RN.set(eidx, this.servt.get(eidx));
            UT.set(eidx, this.tput.get(eidx) * this.servt.get(eidx));
            List<Integer> acts = phActsOf(eidx);
            if (!acts.isEmpty()) {
                double s = 0;
                for (int aidx : acts) {
                    s += PN.get(aidx);
                }
                PN.set(eidx, s);
            }
            // ResidT is reported per visit to the TASK, not per execution of the
            // activity: an activity of this entry runs EXECS times per invocation, and
            // the entry takes SHARE of the task's invocations. RespT stays per
            // execution.
            if (phExecs[eidx] != null) {
                double[] ex = phExecs[eidx];
                double w = phShare[eidx];
                for (int aidx : acts) {
                    WN.set(aidx, w * ex[aidx] * this.residt.get(aidx));
                }
            }
            if (Double.isNaN(UT.get(tidx))) {
                UT.set(tidx, 0);
            }
            UT.set(tidx, UT.get(tidx) + UT.get(eidx));
        }

        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            if (this.ignore.get(tidx) != 0) {
                continue;
            }
            TN.set(tidx, this.tput.get(tidx));
            List<Integer> acts = phActsOf(tidx);
            if (!acts.isEmpty()) {
                double p = 0;
                double w = 0;
                for (int aidx : acts) {
                    p += PN.get(aidx);
                    if (!Double.isNaN(WN.get(aidx))) {
                        w += WN.get(aidx);
                    }
                }
                PN.set(tidx, p);
                WN.set(tidx, w);
            }
        }

        for (int hidx = 0; hidx < lqn.nhosts; hidx++) {
            TN.set(hidx, Double.NaN); // kept NaN for consistency with LQNS
        }

        // Idle, not undefined -- the same rule getEnsembleAvg applies, and for the
        // same reason: an unreachable element reports zero for the measures its kind
        // HAS and NaN for the ones it never has, so that the table's NaN mask
        // survives a disconnected component. Reported columns here are QLen=UT,
        // Util=PN, RespT=RN, ResidT=WN, ArvR=AN, Tput=TN; the pre-swap QN and UN are
        // discarded below and are not written.
        for (int idx = 0; idx < nidx; idx++) {
            if (this.ignore.get(idx) != 0) {
                PN.set(idx, 0.0);           // every kind reports a utilization
                AN.set(idx, Double.NaN);    // nothing reports an arrival rate on an LQN
                switch ((int) lqn.type.get(idx)) {
                    case LayeredNetworkElement.PROCESSOR:
                        UT.set(idx, Double.NaN);
                        RN.set(idx, Double.NaN);
                        WN.set(idx, Double.NaN);
                        TN.set(idx, Double.NaN);
                        break;
                    case LayeredNetworkElement.TASK:
                        UT.set(idx, 0.0);
                        RN.set(idx, Double.NaN);
                        WN.set(idx, 0.0);
                        TN.set(idx, 0.0);
                        break;
                    case LayeredNetworkElement.ENTRY:
                        UT.set(idx, 0.0);
                        RN.set(idx, 0.0);
                        WN.set(idx, Double.NaN);
                        TN.set(idx, 0.0);
                        break;
                    case LayeredNetworkElement.ACTIVITY:
                        UT.set(idx, 0.0);
                        RN.set(idx, 0.0);
                        WN.set(idx, 0.0);
                        TN.set(idx, 0.0);
                        break;
                    default:
                        break;
                }
            }
        }

        return new Matrix[]{UT, PN, RN, TN, AN, WN};
    }

    // =====================================================================
    // Feature gate
    // =====================================================================

    /**
     * Features the collapsed phLayers cannot represent are refused by name rather
     * than silently degraded -- see _kb/06-solver-catalog.md (LN section).
     */
    private void assertSrvnPHSupported() {
        assertSrvnPHSupported(false);
    }

    /**
     * Features the composed law cannot represent are refused by name rather than
     * silently degraded. The list is a property of the ENCODING, so it is the
     * same under either layering; what the squashing adds on top is refused in
     * phFlatServerSet.
     *
     * @param flat true when the caller is method 'flat.ph'
     */
    private void assertSrvnPHSupported(boolean flat) {
        String mname = flat ? "flat.ph" : "srvn.ph";
        if (this.hasPhase2) {
            throw new IllegalStateException("method='" + mname + "' does not support second-phase activities: "
                    + "the composed entry law has no reply point. Use method='default'.");
        }
        for (int cidx = 0; cidx < lqn.ncalls; cidx++) {
            if (lqn.calltype.get(cidx) == CallType.FWD) {
                throw new IllegalStateException("method='" + mname + "' does not support forwarding calls, whose "
                        + "target is not part of the caller's activity graph. Use method='default'.");
            }
        }
        if (lqn.iscache != null) {
            for (int i = 0; i < lqn.iscache.getNumCols(); i++) {
                if (lqn.iscache.get(0, i) != 0) {
                    throw new IllegalStateException(
                            "method='" + mname + "' does not support cache tasks. Use method='default'.");
                }
            }
        }
        // A SetupTask IS supported: the setup is not part of the activity graph, so it
        // never enters the series-parallel reduction and is prefixed to the composed
        // entry law afterwards as the phase-type mixture. An INF task is the exception,
        // as in LDES: it holds no thread to power down, so the cycle has no meaning.
        if (lqn.hassetup != null) {
            for (int i = 0; i < lqn.hassetup.getNumCols(); i++) {
                if (lqn.hassetup.get(0, i) == 0) {
                    continue;
                }
                if (lqn.sched.get(i) == SchedStrategy.INF
                        || Double.isInfinite(lqn.mult.get(0, i))) {
                    throw new IllegalStateException("method='" + mname + "': task '" + lqn.names.get(i)
                            + "' declares a setup time on an infinite-server task, which holds no thread "
                            + "to power down; give it a finite multiplicity.");
                }
            }
        }
        if (lqn.callgroups != null && !lqn.callgroups.isEmpty()) {
            // The group states the ORDER in which one caller visits several
            // callees, and the composed law folds every call into one visit, so
            // the order has nowhere to be expressed. Squashing does not recover
            // it: 'flat.cs' is the only encoding that dispatches a group.
            throw new IllegalStateException("method='" + mname + "' does not support routed call groups, "
                    + "whose dispatch order is a routing property. Use method='flat.cs'.");
        }
        if (lqn.lincon != null && !lqn.lincon.isEmpty()) {
            throw new IllegalStateException("method='" + mname + "' does not support admission constraints on a "
                    + "phLayers station. Use method='default'.");
        }
        // A queue-dependent service rate is a property of the layer STATION, and the
        // composed law replaces that station by an entry law, so the scaling has nowhere
        // to attach. Only buildLayersRecursive emits it; this encoding used to DROP it in
        // silence, which reads as a solved model rather than a refused one.
        assertNoRateDependence(lqn.lldscaling, "lldscaling", mname);
        assertNoRateDependence(lqn.cdscaling, "cdscaling", mname);
        assertNoRateDependence(lqn.jdscaling, "jdscaling", mname);
        assertNoRateDependence(lqn.pools, "server pools", mname);
    }

    /** Refuse a rate dependence that the composed phase-type law cannot carry. */
    private void assertNoRateDependence(Map<Integer, ?> dep, String fname, String mname) {
        if (dep == null || dep.isEmpty()) {
            return;
        }
        int sidx = Collections.min(dep.keySet());
        throw new IllegalStateException("method='" + mname + "' does not support queue-dependent "
                + "service rates on a layer station ('" + lqn.names.get(sidx) + "' declares "
                + fname + "). Use method='srvn.cs'.");
    }

    // =====================================================================
    // Small helpers
    // =====================================================================

    /** Tasks that run on processor HIDX and reach it with requests. */
    private List<Integer> phHostLayerCallers(int hidx) {
        List<Integer> out = new ArrayList<Integer>();
        List<Integer> tasks = lqn.tasksof.get(hidx);
        if (tasks == null) {
            return out;
        }
        for (int tidx : tasks) {
            if (this.ignore.get(tidx) != 0) {
                continue;
            }
            if (lqn.isref.get(tidx) != 0) {
                out.add(tidx);
                continue;
            }
            boolean served = false;
            for (int eidx : phEntriesOf(tidx)) {
                if (phAnyCallerOf(eidx) || phHasOpenArrival(eidx)) {
                    served = true;
                    break;
                }
            }
            if (served) {
                out.add(tidx);
            }
        }
        return out;
    }

    /** Tasks issuing a synchronous call to an entry of TIDX. */
    private List<Integer> phTaskLayerCallers(int tidx) {
        List<Integer> out = new ArrayList<Integer>();
        for (int c = lqn.tshift; c < lqn.tshift + lqn.ntasks; c++) {
            if (c == tidx || this.ignore.get(c) != 0) {
                continue;
            }
            for (int eidx : phEntriesOf(tidx)) {
                if (lqn.issynccaller.get(c, eidx) != 0) {
                    out.add(c);
                    break;
                }
            }
        }
        return out;
    }

    /** Asynchronous calls whose target entry belongs to TIDX. */
    private List<Integer> phAsyncCallsInto(int tidx) {
        List<Integer> out = new ArrayList<Integer>();
        List<Integer> targets = phEntriesOf(tidx);
        for (int cidx = 0; cidx < lqn.ncalls; cidx++) {
            if (lqn.calltype.get(cidx) == CallType.ASYNC
                    && targets.contains((int) lqn.callpair.get(cidx, 1))) {
                out.add(cidx);
            }
        }
        return out;
    }

    /**
     * True when an entry arrival is the ONLY way requests reach task TIDX.
     * 'srvn.ph' refuses forwarding calls outright, so sync/async callers are the
     * whole test.
     */
    private boolean phOpenArrivalOnly(int tidx) {
        if (lqn.isref.get(tidx) != 0) {
            return false;
        }
        for (int eidx : phEntriesOf(tidx)) {
            if (phAnyCallerOf(eidx)) {
                return false;
            }
        }
        for (int eidx : phEntriesOf(tidx)) {
            if (phHasOpenArrival(eidx)) {
                return true;
            }
        }
        return false;
    }

    private boolean phAnyCallerOf(int eidx) {
        for (int i = 0; i < lqn.issynccaller.getNumRows(); i++) {
            if (lqn.issynccaller.get(i, eidx) != 0) {
                return true;
            }
        }
        for (int i = 0; i < lqn.isasynccaller.getNumRows(); i++) {
            if (lqn.isasynccaller.get(i, eidx) != 0) {
                return true;
            }
        }
        return false;
    }

    private boolean phHasOpenArrival(int eidx) {
        return lqn.arrival != null && lqn.arrival.get(eidx) != null;
    }

    /**
     * Replicas of the server station, with the same fan-out reduction as the
     * default builder: a caller that reaches every replica sees one
     * representative.
     */
    private int phReplicaCount(int idx, List<Integer> callers, boolean ishost) {
        int raw = (int) lqn.repl.get(0, idx);
        if (raw <= 1 || callers.isEmpty()) {
            return FastMath.max(1, raw);
        }
        boolean reduce = false;
        if (!ishost && lqn.fanout != null) {
            reduce = true;
            for (int c : callers) {
                if (lqn.fanout.get(c, idx) < raw) {
                    reduce = false;
                    break;
                }
            }
        } else if (ishost) {
            reduce = true;
            for (int c : callers) {
                if ((int) lqn.repl.get(0, c) != raw) {
                    reduce = false;
                    break;
                }
            }
        }
        if (reduce) {
            if (!ishost) {
                this.singleReplicaTasks.add(idx);
            }
            return 1;
        }
        return raw;
    }

    /** Threads of caller C present in the phLayers of IDX. */
    private double phLayerPopulation(int idx, int c, int nreplicas) {
        boolean single = (nreplicas == 1 && lqn.repl.get(0, idx) > 1)
                || this.singleReplicaTasks.contains(c);
        double njobs = single ? lqn.maxmult.get(0, c) : lqn.maxmult.get(0, c) * lqn.repl.get(0, c);
        if (Double.isInfinite(njobs)) {
            njobs = 0;
            // taskgraph spans the hosts and tasks only, not the entries and activities
            for (int i = 0; i < lqn.taskgraph.getNumRows(); i++) {
                if (lqn.taskgraph.get(i, c) != 0) {
                    njobs += lqn.maxmult.get(0, i);
                }
            }
            if (Double.isInfinite(njobs) || njobs == 0) {
                double s = 0;
                for (int i = 0; i < lqn.maxmult.getNumCols(); i++) {
                    double m = lqn.maxmult.get(0, i);
                    if (!Double.isInfinite(m) && !Double.isNaN(m)) {
                        s += m * lqn.repl.get(0, i);
                    }
                }
                njobs = FastMath.min(s, 1000);
            }
        }
        return njobs;
    }

    /** Total closed population of a phLayers, i.e. how many jobs a job can queue behind. */
    private double phLayerPop(PHLayer L, int idx) {
        if (L.npop >= 1) {
            return L.npop;
        }
        double n = 0;
        for (int c : L.callers) {
            double v = this.njobs.get(c, idx);
            if (!Double.isNaN(v) && !Double.isInfinite(v) && v > 0) {
                n += v;
            }
        }
        return n < 1 ? 1 : n;
    }

    /**
     * Residence time per visit, by Little from the queue length rather than from
     * the reported RN. A phLayers that saturates can come back from AMVA with an RN
     * that no closed model can produce, and a reconstruction that trusts it feeds
     * the impossible value straight back into the call response times.
     */
    private static double phResidence(double Q, double X, double RN) {
        if (!Double.isNaN(Q) && !Double.isInfinite(Q) && Q >= 0
                && !Double.isNaN(X) && !Double.isInfinite(X) && X > GlobalConstants.FineTol) {
            return Q / X;
        }
        return RN;
    }

    /**
     * Ratio of a phResidence time to the mean of the law it was measured against,
     * bounded above by the phLayers population: a job can wait behind at most every
     * other job in a closed phLayers.
     */
    private static double phInflationOf(double R, double S, double npop) {
        double f = 1;
        if (S > GlobalConstants.FineTol && !Double.isNaN(R) && !Double.isInfinite(R) && R > 0) {
            f = R / S;
        }
        if (Double.isNaN(f) || Double.isInfinite(f) || f < 1) {
            f = 1;
        }
        if (!Double.isNaN(npop) && !Double.isInfinite(npop) && npop >= 1 && f > npop) {
            f = npop;
        }
        return f;
    }

    /** Synchronous calls issued by task C to an entry of task TIDX. */
    private List<Integer> phSyncCallsBetween(int c, int tidx) {
        List<Integer> out = new ArrayList<Integer>();
        for (int cidx = 0; cidx < lqn.ncalls; cidx++) {
            if (lqn.calltype.get(cidx) != CallType.SYNC) {
                continue;
            }
            if ((int) lqn.parent.get(0, (int) lqn.callpair.get(cidx, 0)) == c
                    && (int) lqn.parent.get(0, (int) lqn.callpair.get(cidx, 1)) == tidx) {
                out.add(cidx);
            }
        }
        return out;
    }

    /**
     * Probability that a request for entry EIDX finds its task's thread powered
     * off. ONE closure for both methods: SolverLN.setupCharge returns p*s, so p
     * is that over s. It also answers p = 1 during construction, before the first
     * solve has sized tput or util.
     */
    private double phSetupProb(int eidx) {
        int tidx = (int) lqn.parent.get(0, eidx);
        if (lqn.hassetup == null || lqn.hassetup.getNumCols() <= tidx || lqn.hassetup.get(0, tidx) == 0) {
            return 0;
        }
        double d = setupMeanOf(lqn.delayofftime, tidx);
        double s = setupMeanOf(lqn.setuptime, tidx);
        if (!(d > GlobalConstants.FineTol) || !(s > GlobalConstants.FineTol)) {
            return 0;
        }
        return FastMath.min(1, FastMath.max(0, setupCharge(tidx) / s));
    }

    /** Phase-type law of task TIDX's setup time, null when it declares none. */
    private Pair<Matrix, Matrix> phSetupLaw(int tidx) {
        if (lqn.setuptime == null) {
            return null;
        }
        Distribution proc = lqn.setuptime.get(tidx);
        if (proc == null) {
            return null;
        }
        double m = proc.getMean();
        if (Double.isNaN(m) || Double.isInfinite(m) || m <= GlobalConstants.FineTol) {
            return null;
        }
        double scv = proc.getSCV();
        if (Double.isNaN(scv) || Double.isInfinite(scv) || scv <= GlobalConstants.FineTol) {
            scv = 1.0;
        }
        APH law = APH.fitMeanAndSCV(m, scv);
        return new Pair<Matrix, Matrix>(law.getInitProb(), law.getSubgenerator());
    }


    /**
     * Divisor that scales a processor utilization into [0,1]. An infinite server
     * reports a mean number of busy servers instead, so it divides by one.
     */
    private double phHostServers(int hidx) {
        if (lqn.sched.get(hidx) == SchedStrategy.INF) {
            return 1;
        }
        double m = lqn.maxmult.get(0, hidx);
        return (!Double.isInfinite(m) && !Double.isNaN(m) && m > 0) ? m : 1;
    }



    private double phMaxNjobs(int tidx) {
        double m = 0;
        for (int j = 0; j < this.njobs.getNumCols(); j++) {
            m = FastMath.max(m, this.njobs.get(tidx, j));
        }
        return m;
    }

    private List<Integer> phEntriesOf(int idx) {
        List<Integer> v = lqn.entriesof.get(idx);
        return v == null ? new ArrayList<Integer>() : v;
    }

    private List<Integer> phActsOf(int idx) {
        List<Integer> v = lqn.actsof.get(idx);
        return v == null ? new ArrayList<Integer>() : v;
    }

    private List<Integer> phCallsOf(int idx) {
        List<Integer> v = lqn.callsof.get(idx);
        return v == null ? new ArrayList<Integer>() : v;
    }

    private SolverResult phLastResult(int layerIdx) {
        return this.results.get(this.results.size()).get(layerIdx);
    }

    private static double phSumOver(Matrix M, int[] stations, int col) {
        double s = 0;
        for (int st : stations) {
            double v = M.get(st - 1, col);
            if (!Double.isNaN(v)) {
                s += v;
            }
        }
        return s;
    }

    private static Matrix phNanRow(int nidx) {
        Matrix M = new Matrix(1, nidx, nidx);
        for (int i = 0; i < nidx; i++) {
            M.set(i, Double.NaN);
        }
        return M;
    }

    private static void phAddUpdMap(java.util.Map<Integer, List<Integer[]>> cell, int idx, Integer[] row) {
        if (!cell.containsKey(idx)) {
            cell.put(idx, new ArrayList<Integer[]>());
        }
        cell.get(idx).add(row);
    }

    private static int[] phToArray(List<Integer> l) {
        int[] a = new int[l.size()];
        for (int i = 0; i < l.size(); i++) {
            a[i] = l.get(i);
        }
        return a;
    }
}
