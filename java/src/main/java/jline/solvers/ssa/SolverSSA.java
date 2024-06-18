package jline.solvers.ssa;

import jline.lang.*;
import jline.lang.constant.*;
import jline.lang.nodes.StatefulNode;
import jline.lang.nodes.Station;
import jline.lang.state.EventCache;
import jline.lang.state.EventCacheKey;
import jline.lang.state.EventResult;
import jline.lang.state.State;
import jline.examples.ClosedModel;
import jline.lib.KPCToolbox;
import jline.solvers.*;
import jline.util.Maths;
import jline.util.Matrix;
import jline.util.Pair;
import jline.util.UniqueRowResult;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import static java.lang.Double.POSITIVE_INFINITY;
import static java.lang.Math.ceil;
import static java.util.stream.Collectors.toMap;

public class SolverSSA extends NetworkSolver {

    private final int DEFAULT_THREADS = (int) Math.ceil(Runtime.getRuntime().availableProcessors() / 2);
    private ExecutorService threadPool;
    private int numThreads = DEFAULT_THREADS;
    private EventCache eventCache;

    public SolverSSA(Network model) {
        // If no options provided, use default options
        this(model, new SSAOptions());
        this.result = new SolverSSAResult();
    }

    public SolverSSA(Network model, Object... varargin) {
        this(model, SolverSSA.defaultOptions());
        this.options = Solver.parseOptions(this.options, varargin);
    }

    public SolverSSA(Network model, String method) {
        super(model, "SolverSSA", SolverSSA.defaultOptions().method(method));
        this.result = new SolverSSAResult();
    }

    public SolverSSA(Network model, SolverOptions options) {
        super(model, "SolverSSA", options);
        this.result = new SolverSSAResult();
    }

    public SSAValues solver_ssa(EventCache eventCache) {
        NetworkStruct sn = this.sn;
        SolverOptions options = this.options;

        if (Maths.randomMatlabStyle()) {
            Maths.setMatlabRandomSeed(options.seed);
        }



        // generate local state spaces

        int nstateful = sn.nstateful;
        int R = sn.nclasses;
        Matrix N = sn.njobs.transpose();
        Map<Integer, Sync> sync = sn.sync;
        Matrix csmask = sn.csmask;

        double cutoff_value = options.cutoff;

        Matrix cutoff = new Matrix(sn.nstations, sn.nclasses);
        cutoff.fill(cutoff_value);

        Matrix Np = N.transpose();
        Matrix capacityc = new Matrix(sn.nnodes, sn.nclasses);
        capacityc.zero();

        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.isstation.get(ind) == 1) {
                int ist = (int) sn.nodeToStation.get(ind);
                for (int r = 0; r < sn.nclasses; r++) {
                    int c = 0;
                    for (int i = 0; i < sn.chains.getNumRows(); i++) {
                        if (sn.chains.get(i, r) == 1) {
                            c = i;
                        }
                    }

                    Matrix proc_m = sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0);
                    boolean disabled = false;
                    for (int row = 0; row < proc_m.getNumRows(); row++) {
                        for (int col = 0; col < proc_m.getNumCols(); col++) {
                            if (Double.isNaN(proc_m.get(row, col))) {
                                disabled = true;
                            }
                        }
                    }

                    if (!sn.visits.get(c).isEmpty() && sn.visits.get(c).get(ist, r) == 0) {
                        capacityc.set(ind, r, 0);
                    } else if (!sn.proc.isEmpty() && !sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).isEmpty() && disabled) {
                        capacityc.set(ind, r, 0);
                    } else {
                        if (N.get(r) == Double.POSITIVE_INFINITY) {
                            capacityc.set(ind,r, Maths.min(cutoff.get(ist,r), sn.classcap.get(ist,r)));
                        } else {
                            // sum values in sn,njobs at indices where the c-th row of sn.chains is true
                            int njobs_sum = 0;
                            for (int i = 0 ; i < sn.njobs.getNumCols(); i++) {
                                if (sn.chains.get(c, i) == 1) {
                                    njobs_sum += (int) sn.njobs.get(i);
                                }
                            }
                            capacityc.set(ind,r, njobs_sum);
                        }
                    }
                }
                int capacity_sum = (int) capacityc.sumRows().get(ind);
                if (sn.nservers.get(ist) == Double.POSITIVE_INFINITY) {
                    sn.nservers.set(ist, capacity_sum);
                }
                for (int col = 0; col < sn.cap.getNumCols(); col++) {
                    sn.cap.set(ist, col, capacity_sum);
                }
                for (int col = 0; col < sn.classcap.getNumCols(); col++) {
                    sn.classcap.set(ist, col, capacityc.get(ind, col));
                }

            }
        }

        if (State.isinf(Np)) {
            // set all elements of Np where theyre Ifinity to 0
            for (int col = 0; col < Np.getNumCols(); col++) {
                if (Np.get(0, col) == Double.POSITIVE_INFINITY) {
                    Np.set(0, col, 0);
                }
            }
        }

        Matrix init_state_hashed = new Matrix(1, nstateful);
        init_state_hashed.zero(); // pick first state in space{i}

        Map<Integer, Matrix> arvRatesSamples = new HashMap<>();
        Map<Integer, Matrix> depRatesSamples = new HashMap<>();
        for (int r = 0; r < R; r++) {
            Matrix m = new Matrix((int) options.samples, nstateful);
            m.zero();
            arvRatesSamples.put(r, m);
            depRatesSamples.put(r, m.clone());
        }
        int A = sync.size();
        int samples_collected = 1;
        Matrix state = init_state_hashed.clone();
        Map<Integer, Matrix> stateCell = new HashMap<>();
        Map<Integer, Matrix> nir = new HashMap<>();

        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.isstateful.get(ind) == 1) {
                int isf = (int) sn.nodeToStateful.get(ind);
                Matrix state_space =  Matrix.extractRows(sn.space.get(sn.stations.get(isf)), (int) state.get(isf),
                        (int) state.get(isf) + 1, null);
                stateCell.put(isf, state_space);

                if (sn.isstation.get(ind) == 1) {
                    int ist = (int) sn.nodeToStation.get(ind);
                    // transpose to get nir as column i
                    nir.put(ist, State.toMarginal(sn, ind, state_space, null, null,
                            null, null, null).nir.transpose());
                }
            }
        }


        state = new Matrix(0,0);
        Matrix statelen = new Matrix(stateCell.size(), 1);
        for (int ind = 0; ind < stateCell.size(); ind++) {
            if (stateCell.containsKey(ind)) {
                Matrix row = stateCell.get(ind);
                if (state.isEmpty()) {
                    state = row;
                } else {
                    state = Matrix.concatColumns(state, row, null);
                }
                statelen.set(ind, row.getNumElements());
            } else {
                statelen.set(ind, 0);
            }
        }

        Matrix tranSync = new Matrix(samples_collected, 1);
        tranSync.zero();
        Matrix z = new Matrix(1,1);
        z.zero();
        Matrix tranState = Matrix.concatColumns(z, state, null).transpose();
        samples_collected = 1;

        Matrix SSq = new Matrix(0, 0);
        for (int ind = 0; ind < nir.size(); ind++) {
            Matrix col = nir.get(ind);
            if (SSq.isEmpty()) {
                SSq = col;
            } else {
                SSq = Matrix.concatRows(SSq, col, null);
            }
        }
        int local = sn.nnodes;

        Map<Integer, Integer> node_a = new HashMap<>();
        Map<Integer, Integer> node_p = new HashMap<>();
        Map<Integer, Integer> class_a = new HashMap<>();
        Map<Integer, Integer> class_p = new HashMap<>();
        Map<Integer, EventType> event_a = new HashMap<>();
        Map<Integer, EventType> event_p = new HashMap<>();
        Map<Integer, Double> outprob_a = new HashMap<>();
        Map<Integer, Double> outprob_p = new HashMap<>();
        for (int act = 0; act < A; act++) {
            NetworkEvent active = sync.get(act).active.get(0);
            NetworkEvent passive = sync.get(act).passive.get(0);
            node_a.put(act, active.getNodeIdx());
            node_p.put(act, passive.getNodeIdx());
            class_a.put(act, active.getJobclassIdx());
            class_p.put(act, passive.getJobclassIdx());
            event_a.put(act, active.getEvent());
            event_p.put(act, passive.getEvent());
        }

        Map<Integer, Map<Integer, Matrix>> newStateCell = new HashMap<>();
        boolean isSimulation = true; // allow state vector to grow, e.g., for FCFS buffers
        double cur_time = 0;
        Map<Integer, Double> enabled_rates = new HashMap<>();
        Map<Integer, Integer> enabled_sync = new HashMap<>();
        Map<Integer, Matrix> stateCell_1 = new HashMap<>();
        while (samples_collected < options.samples && cur_time <= options.timespan[1]) {
            if (samples_collected == 1) {
                if (options.method.equals("para") || options.method.equals ("parallel")) {
                    System.out.printf("SSA samples: %6d\n", samples_collected);
                }
            }
            int ctr = 1;
            Map<Integer, Integer> node_a_sf = new HashMap<>();
            Map<Integer, Integer> node_p_sf = new HashMap<>();
            Map<Integer, Double> prob_sync_p = new HashMap<>();
            enabled_rates.clear();
            enabled_sync.clear();

            for (int act = 0; act < A; act++) {
                Map<Integer, Matrix> rate_a = new HashMap<>();
                newStateCell.put(act, stateCell.entrySet().stream().collect(toMap(Map.Entry::getKey, e -> e.getValue().clone())));
                {
                    int isf = (int) sn.nodeToStateful.get(node_a.get(act));
                    EventResult eventResult = State.afterEvent(sn, node_a.get(act),
                            stateCell.get(isf), event_a.get(act), class_a.get(act), isSimulation, eventCache);
                    if (!eventResult.outspace.isEmpty()) {
                        newStateCell.get(act).put((int) sn.nodeToStateful.get(node_a.get(act)), eventResult.outspace);
                    } else {
                        newStateCell.get(act).remove((int) sn.nodeToStateful.get(node_a.get(act)));
                    }
                    if (!eventResult.outrate.isEmpty()) {
                        rate_a.put(act, eventResult.outrate);
                    } else {
                        rate_a.remove(act);
                    }
                    if (!eventResult.outprob.isEmpty()) {
                        outprob_a.put(act, eventResult.outprob.toDouble());
                    } else {
                        outprob_a.remove(act);
                    }
                }

                if (!newStateCell.get(act).containsKey((int) sn.nodeToStateful.get(node_a.get(act))) || !rate_a.containsKey(act)) {
                    continue;
                }

                for (int ia = 0; ia < newStateCell.get(act).
                        get((int) sn.nodeToStateful.get(node_a.get(act))).getNumRows(); ia++) {
                    if (Double.isNaN(rate_a.get(act).get(ia)) || rate_a.get(act).get(ia) == 0) {
                        // handling degenerate rate values
                        rate_a.get(act).set(ia, GlobalConstants.Zero);
                    }

                    Matrix hash_check = newStateCell.get(act).get((int) sn.nodeToStateful.get(node_a.get(act)));
                    boolean hash_found = false;
                    for (int col = 0; col < hash_check.getNumCols(); col++) {
                        if (hash_check.get(ia, col) != -1) {
                            hash_found = true;
                            break;
                        }
                    }

                    if (!hash_found) {
                        continue;
                    }

                    boolean update_cond = true;
                    if (rate_a.get(act).get(ia) > 0) {
                        if (node_p.get(act) != local) {
                            if (node_p.get(act).equals(node_a.get(act))) {

                                // self-loop
                                EventResult eventResult = State.afterEvent(sn, node_p.get(act),
                                        newStateCell.get(act).get((int) sn.nodeToStateful.get(node_a.get(act))),
                                        event_p.get(act), class_p.get(act), isSimulation, eventCache);
                                if (!eventResult.outspace.isEmpty()) {
                                    newStateCell.get(act).put((int) sn.nodeToStateful.get(node_p.get(act)), eventResult.outspace);
                                } else {
                                    newStateCell.get(act).remove((int) sn.nodeToStateful.get(node_p.get(act)));
                                }
                                if (!eventResult.outprob.isEmpty()) {
                                    outprob_p.put(act, eventResult.outprob.toDouble());
                                }

                            } else {
                                // departure
                                EventResult eventResult = State.afterEvent(sn, node_p.get(act),
                                        newStateCell.get(act).get((int) sn.nodeToStateful.get(node_p.get(act))),
                                        event_p.get(act), class_p.get(act), isSimulation, eventCache);

                                if (!eventResult.outspace.isEmpty()) {
                                    newStateCell.get(act).put((int) sn.nodeToStateful.get(node_p.get(act)), eventResult.outspace);
                                } else {
                                    newStateCell.get(act).remove((int) sn.nodeToStateful.get(node_p.get(act)));
                                }
                                if (!eventResult.outprob.isEmpty()) {
                                    outprob_p.put(act, eventResult.outprob.toDouble());
                                }
                            }

                            if (newStateCell.get(act).containsKey((int) sn.nodeToStateful.get(node_p.get(act)))) {
                                if (sn.isstatedep.get(node_p.get(act)) == 1) {
                                    // TODO: line 159 MATLAB
                                    throw new RuntimeException("UNIMPLEMENTED CASE");
                                } else {
                                    prob_sync_p.put(act, sync.get(act).passive.get(0).getProb());
                                }
                            } else {
                                prob_sync_p.put(act, 0.0);
                            }
                        }
                        if (newStateCell.get(act).containsKey((int) sn.nodeToStateful.get(node_a.get(act)))) {
                            if (node_p.get(act) == local) {
                                prob_sync_p.put(act, 1.0);
                            }
                            if (!Double.isNaN(rate_a.get(act).toDouble())) {
                                if (!newStateCell.get(act).isEmpty()) {
                                    if (event_a.get(act) == EventType.DEP) {
                                        node_a_sf.put(act, (int) sn.nodeToStateful.get(node_a.get(act)));
                                        node_p_sf.put(act, (int) sn.nodeToStateful.get(node_p.get(act)));

                                        Matrix original_departure = depRatesSamples.get(class_a.get(act));
                                        Matrix original_arrival = arvRatesSamples.get(class_p.get(act));

                                        double added_value = outprob_a.get(act) * outprob_p.get(act) * rate_a.get(act).get(ia)
                                                * prob_sync_p.get(act);

                                        int a_sf_act = node_a_sf.get(act);
                                        int p_sf_act = node_p_sf.get(act);

                                        double dep_value = original_departure.get(samples_collected - 1, a_sf_act)
                                                + added_value;
                                        double arv_val = original_arrival.get(samples_collected - 1, p_sf_act)
                                                + added_value;

                                        original_departure.set(samples_collected - 1, a_sf_act, dep_value);
                                        original_arrival.set(samples_collected - 1, p_sf_act, arv_val);

                                    }
                                    if (node_p.get(act) < local && csmask.get(class_a.get(act), class_p.get(act)) != 1
                                            && sn.nodetypes.get(node_p.get(act)) != NodeType.Source &&
                                            (rate_a.get(act).get(ia) * prob_sync_p.get(act) > 0)) {
                                        // TODO: appropriate exception here
                                        throw new RuntimeException("UNIMPLEMENTED");
                                    }

                                    enabled_rates.put(ctr-1, rate_a.get(act).get(ia) * prob_sync_p.get(act));
                                    enabled_sync.put(ctr-1, act);
                                    ctr++;
                                }
                            }
                        }
                    }
                }
            }

            Matrix enabled_rates_m = new Matrix(1, enabled_rates.size());
            for (int i = 0; i < enabled_rates.size(); i++) {
                enabled_rates_m.set(0, i, enabled_rates.get(i));
            }
            double tot_rate = enabled_rates_m.elementSum();
            Matrix cum_sum = enabled_rates_m.cumsumViaRow();
            Matrix cum_rate = Matrix.scale_mult(cum_sum, 1 / tot_rate);

            double rand = Maths.random();
            int firing_ctr = -1;
            for (int i = 0; i < cum_rate.getNumElements(); i++) {
                if (rand > cum_rate.get(i)) {
                    firing_ctr = i;
                } else {
                    break;
                }
            }
            firing_ctr++;
            if (enabled_sync.isEmpty()) {
                throw new RuntimeException("SSA simulation entered a deadlock before collecting all samples, " +
                        "no synchronization is enabled.");
            }

            // this part is needed to ensure that when the state vector grows the padding of zero is done on the left

            for (int ind = 0; ind < sn.nnodes; ind++) {
                if (sn.isstation.get(ind) == 1) {
                    int isf = (int) sn.nodeToStateful.get(ind);
                    boolean deltalen = (stateCell.get(isf).getNumElements() > statelen.get(isf));
                    if (deltalen) {
                        statelen.set(isf, stateCell.get(isf).getNumElements());
                        // padding
                        int shift = 0;
                        if (ind > 0) {
                            for (int col = 0; col < isf; col++) {
                                shift += (int) statelen.get(col);
                            }
                        }
                        Matrix pad = new Matrix(1, tranState.getNumCols());

                        Matrix top = Matrix.extractRows(tranState, 0, shift + 1, null);
                        Matrix bottom = Matrix.extractRows(tranState, shift + 1, tranState.getNumRows(), null);
                        Matrix tmp = Matrix.concatRows(top, pad, null);
                        tranState = Matrix.concatRows(tmp, bottom, null);
                    }
                }
            }

            state = new Matrix(0, 0);
            for (int ind = 0; ind < sn.nnodes; ind++) {
                if (stateCell.containsKey(ind)) {
                    Matrix row = stateCell.get(ind);
                    if (state.isEmpty()) {
                        state = row;
                    } else {
                        state = Matrix.concatColumns(state, row, null);
                    }
                }
            }

            double dt = -(Math.log(Maths.random()) / tot_rate);
            cur_time += dt;

            Matrix dt_m = new Matrix(1, 1);
            dt_m.set(0,0,dt);
            Matrix newTranState = Matrix.concatColumns(dt_m, state, null);
            // add new column to tranState if needed
            if (samples_collected - 1 >= tranState.getNumCols()) {
                Matrix tranState_new_col = new Matrix(tranState.getNumRows(), 1);
                tranState_new_col.zero();
                tranState = Matrix.concatColumns(tranState, tranState_new_col, null);
            }
            for (int row = 0; row < state.getNumElements() + 1; row++) {
                tranState.set(row, samples_collected - 1, newTranState.get(row));
            }

            // add new row to tranSync if needed
            if (samples_collected - 1 >= tranSync.getNumRows()) {
                Matrix tranSync_new_row = new Matrix(1, tranSync.getNumCols());
                tranSync_new_row.zero();
                tranSync = Matrix.concatRows(tranSync, tranSync_new_row, null);
            }
            // TODO: is this +1 needed? added to make value in tranSync align with LINE always
            tranSync.set(samples_collected - 1, 0, enabled_sync.get(firing_ctr) + 1);

            for (int ind = 0; ind < sn.nnodes; ind++) {
                if (sn.isstation.get(ind) == 1) {
                    int isf = (int) sn.nodeToStateful.get(ind);
                    int ist = (int) sn.nodeToStation.get(ind);
                    nir.put(ist, State.toMarginal(sn, ind, stateCell.get(isf), null, null, null, null, null).nir.transpose());
                }
            }

            // make one big column from all the nir matrices
            Matrix nir_col = new Matrix(0, 0);
            for (int ind = 0; ind < sn.nnodes; ind++) {
                if (nir.containsKey(ind)) {
                    Matrix col = nir.get(ind);
                    if (nir_col.isEmpty()) {
                        nir_col = col;
                    } else {
                        nir_col = Matrix.concatRows(nir_col, col, null);
                    }
                }
            }
            // assign nir_col to samples_collected column of SSq
            // add new column to SSq
            if (samples_collected - 1 >= SSq.getNumCols()) {
                Matrix SSq_new_col = new Matrix(SSq.getNumRows(), 1);
                SSq_new_col.zero();
                SSq = Matrix.concatColumns(SSq, SSq_new_col, null);
            }
            for (int row = 0; row < nir_col.getNumRows(); row++) {
                int col = samples_collected - 1;
                SSq.set(row, col, nir_col.get(row));
            }
            samples_collected++;

            stateCell_1 =
                    stateCell.entrySet().stream().collect(toMap(Map.Entry::getKey, e -> e.getValue().clone()));
            stateCell = newStateCell.get(enabled_sync.get(firing_ctr));


            if (!options.method.equals("para") && !options.method.equals ("parallel") &&
                    (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG)) {
            if (samples_collected == 2) {
                System.out.printf("\b\n SSA samples: %6d ", samples_collected);
            } else if (options.verbose == VerboseLevel.DEBUG) {
                System.out.printf("\b\b\b\b\b\b\b %6d", samples_collected);
            } else if (samples_collected % 100 == 0) {
                System.out.printf("\b\b\b\b\b\b\b %6d", samples_collected);
	            if (samples_collected == options.samples) {
	                System.out.println("\n");
	            }
            }
        }

        }

        tranState = tranState.transpose();

        UniqueRowResult uniqueRows = Matrix.uniqueRows(Matrix.extract(tranState, 0, tranState.getNumRows(),
                1, tranState.getNumCols()));
        Matrix u = uniqueRows.sortedMatrix;
        Matrix ui = uniqueRows.vi;
        Map<Integer, List<Integer>> uj = uniqueRows.vj;

        Matrix statesz = new Matrix(1, stateCell_1.size());
        for (int i = 0; i < statesz.getNumElements(); i++) {
            statesz.set(i, stateCell_1.get(i).getNumElements());
        }
        Map<Integer, Matrix> tranSysState = new HashMap<>();



        tranSysState.put(0, Matrix.extractColumn(tranState, 0, null).cumsumViaCol());
        for (int j = 0; j < statesz.getNumElements(); j++) {
            int start_index = 0;
            for (int i = 0; i <= j - 1; i++) {
                start_index += (int) statesz.get(i);

            }
            start_index++;
            int end_index = (int) (start_index + statesz.get(j));
            Matrix tmp = Matrix.extract(tranState, 0, tranState.getNumRows(), start_index, end_index);
            tranSysState.put(j+1, tmp);
        }

        Map<Integer, Matrix> arvRates = new HashMap<>();
        Map<Integer, Matrix> depRates = new HashMap<>();
        for (int i = 0; i < R; i++) {
            Matrix rate = new Matrix(u.getNumRows(), sn.nstateful);
            arvRates.put(i, rate);
            depRates.put(i, rate.clone());
        }
        Matrix pi = new Matrix(1, u.getNumRows());
        for (int s = 0; s < u.getNumRows(); s++) {
            // sum elements in tranState_first where uj == s
            double sum = 0;
            for (Integer i : uj.get(s)) {
                sum += tranState.get(i, 0);
            }

            pi.set(s, sum);
        }

        // extract columns of SSq using ui as indices
        Matrix SSq_new = new Matrix(0,0);
        for (int col_ind = 0; col_ind < ui.getNumElements(); col_ind++) {
            int col = (int) ui.get(col_ind);
            Matrix row = Matrix.extractColumn(SSq, col, null);
            if (SSq_new.isEmpty()) {
                SSq_new = row;
            } else {
                SSq_new = Matrix.concatColumns(SSq_new, row, null);
            }
        }
        SSq = SSq_new.clone().transpose();
        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.isstateful.get(ind) == 1) {
                int isf = (int) sn.nodeToStateful.get(ind);
                for (int s = 0; s < u.getNumRows(); s++) {
                    for (int r = 0; r < R; r++) {
                        arvRates.get(r).set(s, isf, arvRatesSamples.get(r).get((int) ui.get(s), isf));
                        arvRates.put(r, arvRates.get(r));
                        depRates.get(r).set(s, isf, depRatesSamples.get(r).get((int) ui.get(s), isf));
                        depRates.put(r, depRates.get(r));

                    }
                }
            }
        }
        pi = Matrix.scale_mult(pi, 1/pi.elementSum());
        if (options.method.equals("para") || options.method.equals ("parallel")) {
            System.out.printf("SSA samples: %6d\n", samples_collected);
        }
        return new SSAValues(pi, SSq, arvRates, depRates, tranSysState, tranSync, sn);
    }





    @Override
    protected void runAnalyzer() throws IllegalAccessException, ParserConfigurationException, IOException {
        threadPool = Executors.newFixedThreadPool(this.numThreads);

        long T0 = java.lang.System.currentTimeMillis();
        if (this.options == null) {
            this.options = new SolverOptions(SolverType.SSA);
        }
        if (this.enableChecks && !this.supports(this.model)) {
            throw new RuntimeException("This model is not supported by the SSA solver.");
        }
//        this.resetRandomGeneratorSeed(options.seed);

        NetworkStruct sn = getStruct();

        SolverSSAResult result = solver_ssa_analyzer();
        Matrix QN = result.QN;
        Matrix UN = result.UN;
        Matrix RN = result.RN;
        Matrix TN = result.TN;
        Matrix CN = result.CN;
        Matrix XN = result.XN;
        Map<Integer, Matrix> tranSysState = result.tranSysState;
        Matrix tranSync = result.tranSync;
        sn = result.sn;


        for (int isf=0; isf < sn.nstateful; isf++) {
            int ind = (int) sn.statefulToNode.get(isf);
            if (sn.nodetypes.get((int) sn.statefulToNode.get(isf)) == NodeType.Cache) {
                // TODO: Cache nodetype case
            }
        }
        long runtime = java.lang.System.currentTimeMillis() - T0;
        int M = sn.nstations;
        int R = sn.nclasses;
        Map<Station, Map<JobClass, SolverHandles.Metric>> T = getAvgTputHandles();
        Matrix AN = new Matrix(0,0);
        if (!T.isEmpty() && !TN.isEmpty()) {
            AN = new Matrix(M, R);
            AN.zero();
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    for (int k = 0; k < R; k++) {
                        for (int r = 0; r < R; r++) {
                            AN.set(i,k, AN.get(i,k) + TN.get(j,r) * sn.rt.get((j) * R + r, (i) * R + k));
                        }
                    }
                }
            }
        }

        SolverSSAResult res = new SolverSSAResult();
        res.QN = QN;
        res.UN = UN;
        res.RN = RN;
        res.TN = TN;
        res.AN = AN;
        res.WN = new Matrix(0,0);
        res.CN = CN;
        res.XN = XN;
        res.runtime = runtime;
        this.setAvgResults(res);
        // TODO: safe to downcast here, but indicative of bad oop design.
        //  "Result" is a SolverResult field in Solver superclass and space is specific to SSA
        SolverSSAResult ssaRes = (SolverSSAResult) this.result;
        ssaRes.space = sn.space;


    }


    private SolverSSAResult solver_ssa_analyzer() {
        long Tstart = java.lang.System.currentTimeMillis();

        sn.space.clear();
        // TODO: check conversion of Stateful node to station
        for (StatefulNode statefulNode : sn.state.keySet()) {
            int node = statefulNode.getStationIdx();
            sn.space.put(sn.stations.get(node), sn.state.get(statefulNode));
        }

        SolverSSAResult res = new SolverSSAResult();
        switch (options.method) {
            case "default":
            case "serial":
            case "ssa":
                res = solver_ssa_analyzer_serial(false);
                break;
            case "para":
            case "parallel":
                res = solver_ssa_analyzer_parallel();
                break;

        }
        // measure time
        res.runtime = java.lang.System.currentTimeMillis() - Tstart;
        return res;
    }

    private SolverSSAResult solver_ssa_analyzer_parallel() {
        int M = sn.nstations;
        int K = sn.nclasses;
        Map<Station, Map<JobClass, Map<Integer, Matrix>>> PH = sn.proc;
        Matrix S = sn.nservers;
        Matrix NK = sn.njobs.transpose();




        Map<Integer, Matrix> QNs = new ConcurrentHashMap<>();
        Map<Integer, Matrix> UNs = new ConcurrentHashMap<>();
        Map<Integer, Matrix> RNs = new ConcurrentHashMap<>();
        Map<Integer, Matrix> TNs = new ConcurrentHashMap<>();
        Map<Integer, Matrix> CNs = new ConcurrentHashMap<>();
        Map<Integer, Matrix> XNs = new ConcurrentHashMap<>();


        this.options.samples = (int) ceil(this.options.samples / (double) numThreads);

        boolean matlabRand = Maths.randomMatlabStyle();
        for (int t = 0; t < numThreads; t++) {
            final int threadIndex = t;
            threadPool.submit(() -> {

                if (matlabRand) {
                    Maths.setMatlabRandomSeed(options.seed);
                }

                Matrix probSysState = new Matrix(0,0);
                Matrix SSq = new Matrix(0,0);
                Map<Integer, Matrix> arvRates = new HashMap<>();
                Map<Integer, Matrix> depRates = new HashMap<>();

                switch (this.options.method) {
                    case "para":
                    case "parallel":
                        eventCache = new EventCache(true, this.options.config.eventcache);
                        SSAValues result = solver_ssa(eventCache);
                        probSysState = result.pi;
                        SSq = result.SSq;
                        arvRates = result.arvRates;
                        depRates = result.depRates;
                        break;
                }
                Matrix XN = new Matrix(1, K);
                XN.fill(Double.NaN);
                Matrix UN = new Matrix(M, K);
                UN.fill(Double.NaN);
                Matrix QN = new Matrix(M, K);
                QN.fill(Double.NaN);
                Matrix RN = new Matrix(M, K);
                RN.fill(Double.NaN);
                Matrix TN = new Matrix(M, K);
                TN.fill(Double.NaN);
                Matrix CN = new Matrix(1, K);
                CN.fill(Double.NaN);



                for (int k = 0; k < K; k++) {
                    double refsf = sn.stationToStateful.get((int) sn.refstat.get(k));
                    Matrix departure = depRates.get(k);
                    Matrix dep_wset_refsf = Matrix.extractColumn(departure, (int) refsf, null);
                    XN.set(k, probSysState.mult(dep_wset_refsf).toDouble());
                    for (int i = 0; i < M; i++) {
                        int isf = (int) sn.stationToStateful.get(i);

                        Matrix dep_isf = Matrix.extractColumn(departure, isf, null);
                        TN.set(i, k, probSysState.mult(dep_isf).toDouble());

                        Matrix ssq_extracted = Matrix.extractColumn(SSq, i*K + k , null);
                        QN.set(i, k, probSysState.mult(ssq_extracted).toDouble());

                        switch (sn.sched.get(sn.stations.get(i))) {
                            case INF:
                                UN.set(i, k, QN.get(i, k));
                                break;
                            default:
                                if (!PH.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).isEmpty()) {
                                    // probSysState*arvRates(:,i,k)*map_mean(PH{i}{k})/S(i);
                                    Matrix arrival = arvRates.get(k);
                                    Matrix arv_ik = Matrix.extractColumn(arrival, i, null);
                                    double map_mean = KPCToolbox.map_mean(PH.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).get(0),
                                            PH.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).get(1)) / S.get(i);
                                    UN.set(i, k, probSysState.mult(arv_ik).toDouble() * map_mean);
                                }
                        }
                    }
                }

                for (int k = 0; k < K; k++) {
                    for (int i = 0; i < M; i++) {
                        if (TN.get(i, k) > 0) {
                            RN.set(i, k, QN.get(i, k) / TN.get(i, k));
                        } else {
                            RN.set(i, k, 0);

                        }
                    }
                    CN.set(k, NK.get(k) / XN.get(k));
                }
                QN.apply(Double.NaN, 0, "equal");
                CN.apply(Double.NaN, 0, "equal");
                RN.apply(Double.NaN, 0, "equal");
                UN.apply(Double.NaN, 0, "equal");
                XN.apply(Double.NaN, 0, "equal");
                TN.apply(Double.NaN, 0, "equal");

                QNs.put(threadIndex, QN);
                CNs.put(threadIndex, CN);
                RNs.put(threadIndex, RN);
                UNs.put(threadIndex, UN);
                XNs.put(threadIndex, XN);
                TNs.put(threadIndex, TN);

            });
//            threads[t].start();
        }
        threadPool.shutdown();
        try {
            threadPool.awaitTermination(600, TimeUnit.SECONDS);
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }


        Matrix XN = new Matrix(1, K);
        for (int i = 0; i < XN.getNumRows(); i++){
            for (int j = 0; j < XN.getNumCols(); j++) {
                for (Matrix xn: XNs.values()) {
                    XN.set(i, j, XN.get(i, j) + xn.get(i, j));
                }
            }
        }
        XN = Matrix.scale_mult(XN,  1 /(double) numThreads);
        Matrix UN = new Matrix(M, K);
        for (int i = 0; i < UN.getNumRows(); i++){
            for (int j = 0; j < UN.getNumCols(); j++) {
                for (Matrix un: UNs.values()) {
                    UN.set(i, j, UN.get(i, j) + un.get(i, j));
                }
            }
        }
        UN = Matrix.scale_mult(UN,  1 /(double) numThreads);

        Matrix QN = new Matrix(M, K);
        for (int i = 0; i < QN.getNumRows(); i++){
            for (int j = 0; j < QN.getNumCols(); j++) {
                for (Matrix qn: QNs.values()) {
                    QN.set(i, j, QN.get(i, j) + qn.get(i, j));
                }
            }
        }
        QN = Matrix.scale_mult(QN,  1 /(double) numThreads);


        Matrix RN = new Matrix(M, K);
        for (int i = 0; i < RN.getNumRows(); i++){
            for (int j = 0; j < RN.getNumCols(); j++) {
                for (Matrix rn: RNs.values()) {
                    RN.set(i, j, RN.get(i, j) + rn.get(i, j));
                }
            }
        }
        RN = Matrix.scale_mult(RN,  1 /(double) numThreads);



        Matrix TN = new Matrix(M, K);
        for (int i = 0; i < TN.getNumRows(); i++){
            for (int j = 0; j < TN.getNumCols(); j++) {
                for (Matrix tn: TNs.values()) {
                    TN.set(i, j, TN.get(i, j) + tn.get(i, j));
                }
            }
        }
        TN = Matrix.scale_mult(TN,  1 /(double) numThreads);


        Matrix CN = new Matrix(1, K);
        for (int i = 0; i < CN.getNumRows(); i++){
            for (int j = 0; j < CN.getNumCols(); j++) {
                for (Matrix cn: CNs.values()) {
                    CN.set(i, j, CN.get(i, j) + cn.get(i, j));
                }
            }
        }
        CN = Matrix.scale_mult(CN,  1 /(double) numThreads);

        Map<Integer, Matrix> tranSysState = new HashMap<>();
        Matrix tranSync = new Matrix(0, 0);
        return new SolverSSAResult(QN, UN, RN, TN, CN, XN, tranSysState, tranSync, sn);

    }




    private SolverSSAResult solver_ssa_analyzer_serial(boolean hash) {
        int M = sn.nstations;
        int K = sn.nclasses;

        Matrix S = sn.nservers;
        Matrix NK = sn.njobs.transpose();
        Map<Station, SchedStrategy> schedid = sn.sched;
        long Tstart = java.lang.System.currentTimeMillis();

        Map<Station, Map<JobClass, Map<Integer, Matrix>>> PH = sn.proc;
        Map<Integer, Matrix> tranSysState = new HashMap<>();
        Matrix tranSync = new Matrix(0,0);

        Matrix XN = new Matrix(1,K);
        XN.fill(Double.NaN);
        Matrix UN = new Matrix(M, K);
        UN.fill(Double.NaN);
        Matrix QN = new Matrix(M, K);
        QN.fill(Double.NaN);
        Matrix RN = new Matrix(M, K);
        RN.fill(Double.NaN);
        Matrix TN = new Matrix(M, K);
        TN.fill(Double.NaN);
        Matrix CN = new Matrix(1, K);
        CN.fill(Double.NaN);

        this.options.samples++;
        eventCache = new EventCache(false, this.options.config.eventcache);
        SSAValues result = solver_ssa(eventCache);
        Matrix probSysState = result.pi;
        Matrix StateSpaceAggr = result.SSq;
        Map<Integer, Matrix> arvRates = result.arvRates;
        Map<Integer, Matrix> depRates = result.depRates;
        tranSysState = result.tranSysState;
        tranSync = result.tranSync;


        for (int k = 0; k < K; k++) {
            double refsf = sn.stationToStateful.get((int) sn.refstat.get(k));
            Matrix departure = depRates.get(k);
            Matrix dep_wset_refsf = new Matrix(StateSpaceAggr.getNumRows(), 1);
            for (int i = 0; i < StateSpaceAggr.getNumRows(); i++) {
                dep_wset_refsf.set(i, 0, departure.get(i, (int) refsf));
            }
            // toDouble call since 1xn mult nx1
            XN.set(k, probSysState.mult(dep_wset_refsf).toDouble());
        }

        for (int i = 0; i < M; i++) {
            int isf = (int) sn.stationToStateful.get(i);
            for (int k = 0; k < K; k++) {
                Matrix departure = depRates.get(k);
                Matrix dep_wset_isf = new Matrix(StateSpaceAggr.getNumRows(), 1);
                for (int j = 0; j < StateSpaceAggr.getNumRows(); j++) {
                    dep_wset_isf.set(j, 0, departure.get(j, isf));
                }
                TN.set(i, k, probSysState.mult(dep_wset_isf).toDouble());

                Matrix ssaggr_wset_isf = Matrix.extract(StateSpaceAggr, 0, StateSpaceAggr.getNumRows(), (i) * K + k, (i)*K+k+1);
                QN.set(i, k, probSysState.mult(ssaggr_wset_isf).toDouble());
            }
            switch (schedid.get(sn.stations.get(i))) {
                case INF:
                    for (int k = 0; k < K; k++) {
                        UN.set(i, k, QN.get(i, k));
                    }
                    break;

                default:
                    if ((sn.lldscaling == null || sn.lldscaling.isEmpty()) && (sn.cdscaling == null || sn.cdscaling.isEmpty())) {
                        for (int k = 0; k < K; k++) {
                            if (!PH.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).isEmpty()) {
                                Matrix arrival = arvRates.get(k);
                                Matrix arv_wset_isf = new Matrix(StateSpaceAggr.getNumRows(), 1);
                                for (int c = 0; c < StateSpaceAggr.getNumRows(); c++) {
                                    arv_wset_isf.set(c, 0, arrival.get(c, isf));
                                }
                                double map_mean = KPCToolbox.map_mean(PH.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).get(0),
                                        PH.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).get(1)) / S.get(i);
                                UN.set(i, k, probSysState.mult(arv_wset_isf).toDouble() * map_mean);

                            }
                        }
                    } else {
                        // lld/cd cases
                        int ind = (int) sn.stationToNode.get(i);
                        for (int col = 0; col < K; col++) {
                            UN.set(i, col, Double.NaN);
                        }
                    }

                    break;
            }
        }

        for (int k = 0; k < K; k++) {
            for (int i = 0; i < M; i++) {
                if (TN.get(i, k) > 0) {
                    RN.set(i, k, QN.get(i, k) / TN.get(i, k));
                } else {
                    RN.set(i, k, 0);

                }
            }
            CN.set(k, NK.get(k) / XN.get(k));
        }

        // update routing probabilities in nodes with state-dependent routing
        // TODO: Cache nodetype case

        // updates cache actual hit and miss data
        // TODO: Cache nodetype case

        // matrices QN, CN, RN, UN, XN, TN, where they are Double.isNan set to 0
        QN.apply(Double.NaN, 0, "equal");
        CN.apply(Double.NaN, 0,  "equal");
        RN.apply(Double.NaN, 0, "equal");
        UN.apply(Double.NaN, 0, "equal");
        XN.apply(Double.NaN, 0, "equal");
        TN.apply(Double.NaN, 0, "equal");




        return new SolverSSAResult(QN, UN, RN, TN, CN, XN, tranSysState, tranSync, sn);
    }

    // set number of threads for para SolverSSA
    public void setParallelism(int numThreads) {
        this.numThreads = numThreads;
    }



    private boolean supports (Network model) {
        FeatureSet usedLangFeatures = model.getUsedLangFeatures();
        // TODO: complete - needs SolverFeatureSet class implementation
        return true;
    }

    private NetworkStruct getStruct() {
        return this.model.getStruct(true);
    }


}
