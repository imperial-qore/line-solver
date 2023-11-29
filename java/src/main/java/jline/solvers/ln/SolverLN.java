package jline.solvers.ln;


import jline.lang.*;

import jline.lang.constant.*;
import jline.lang.distributions.Disabled;
import jline.lang.distributions.Distribution;
import jline.lang.distributions.Exp;
import jline.lang.distributions.Immediate;
import jline.lang.layerednetworks.LayeredNetwork;
import jline.lang.layerednetworks.LayeredNetworkElement;
import jline.lang.layerednetworks.LayeredNetworkStruct;
import jline.lang.nodes.*;

import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;

import jline.solvers.*;
import jline.solvers.mva.SolverMVA;
import jline.util.Matrix;
import org.apache.commons.lang3.NotImplementedException;

import java.util.*;

public class SolverLN extends EnsembleSolver {
    // TODO: unlike MATLAB LayeredNetworkSolver is not a superclass here, implement and add as an interface?
    // registries of quantities to update at every iteration
    private int nlayers;
    private LayeredNetworkStruct lqn;
    private boolean hasconverged;
    private Integer averagingstart;
    private List<Double> idxhash;
    private Matrix servtmatrix;
    private Network[] ensemble;
    private Matrix ptaskcallers;
    private Map<Integer, Matrix> ptaskcallers_step;
    private Matrix ilscaling;
    private Matrix njobs;
    private Matrix njobsorig;
    private List<Integer> routereset;
    private List<Integer> svcreset;
    private List<Double> maxitererr;
    private Matrix unique_route_prob_updmap;

    // performance metrics and related processes
    private Matrix util;
    private Matrix util_ilock;
    private Matrix tput;
    private Map<Integer, Distribution> tputproc;
    private Matrix servt;
    private Map<Integer, Distribution> servtproc;
    private Matrix servtcdf;
    private Matrix thinkt;
    private Map<Integer, Distribution> thinkproc;
    private Map<Integer, Distribution> thinktproc;
    private Matrix entryproc;
    private Matrix entrycdfrespt;
    private Matrix callresidt;
    private Map<Integer, Distribution> callresidtproc;
    private Matrix callresidtcdf;

    // registries of quantities to update at every iteration
    private Matrix arvproc_classes_updmap;
    private Matrix thinkt_classes_updmap;
    private Matrix servt_classes_updmap;
    private Matrix call_classes_updmap;
    private Matrix route_prob_updmap;

    private Map<Integer, List<Integer[]>> cell_arvproc_classes_updmap;
    private Map<Integer, List<Integer[]>> cell_thinkt_classes_updmap;
    private Map<Integer, List<Integer[]>> cell_servt_classes_updmap;
    private Map<Integer, List<Integer[]>> cell_call_classes_updmap;
    private Map<Integer, List<Integer[]>> cell_route_prob_updmap;

    // Temporary variables
    private List<Network> temp_ensemble;
    private JobClass curClassC;

    static class DefaultSolverFactory implements SolverFactory {
        public NetworkSolver at(Network model) {
            return new SolverMVA(model);
        }
    }

    public SolverLN(LayeredNetwork lqnmodel) {
        this(lqnmodel, new SolverOptions(SolverType.LN));
    }

    public SolverLN(LayeredNetwork lqnmodel, SolverFactory solverFactory) {
        this(lqnmodel, solverFactory, new SolverOptions(SolverType.LN));
    }

    public SolverLN(LayeredNetwork lqnmodel, SolverOptions options) {
        this(lqnmodel, new DefaultSolverFactory(), options);
    }

    public SolverLN(LayeredNetwork lqnmodel, SolverFactory solverFactory, SolverOptions options) {
        super(null, "SolverLN", options); // first argument is null as the ensemble cannot be built yet
        this.lqn = lqnmodel.getStruct();
        construct();
        solvers = new NetworkSolver[nlayers];
        for (int i = 0; i < nlayers; i++) {
            solvers[i] = solverFactory.at(ensemble[i]);
        }
    }

    private void construct() {
        // initialize internal data structures
        this.nlayers = 0;
        this.entrycdfrespt = new Matrix(lqn.nentries, 1);
        this.hasconverged = false;

        // initialize svc and think times
        this.servtproc = new HashMap<Integer, Distribution>();
        this.servtproc.putAll(lqn.hostdem);
        this.thinkproc = new HashMap<Integer, Distribution>();
        this.thinkproc.putAll(lqn.think);
        this.callresidtproc = new HashMap<Integer, Distribution>();
        for (int cidx = 1; cidx <= lqn.ncalls; cidx++) {
            callresidtproc.put(cidx, lqn.hostdem.get((int) lqn.callpair.get(cidx, 2)));
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
        for (int i = 1; i < route_prob_updmap.numRows; i++) {
            int buffer = idxhash.get((int) route_prob_updmap.get(i, 1)).intValue();
            if (!routereset.contains(buffer)) { // unique
                routereset.add(buffer);
            }
        }

        this.svcreset = new ArrayList<>();
        for (int i = 1; i < thinkt_classes_updmap.numRows; i++) {
            int buffer = idxhash.get((int) thinkt_classes_updmap.get(i, 1)).intValue();
            if (!svcreset.contains(buffer)) { // unique
                svcreset.add(buffer);
            }
        }
        for (int i = 1; i < call_classes_updmap.numRows; i++) {
            int buffer = idxhash.get((int) call_classes_updmap.get(i, 1)).intValue();
            if (!svcreset.contains(buffer)) { // unique
                svcreset.add(buffer);
            }
        }
        Collections.sort(svcreset);
    }

    public void initFromRawAvgTables() {
        throw new NotImplementedException("initFromRawAvgTables not yet available");
    }

    public void paramFromRawAvgTables() {
        throw new NotImplementedException("initFromRawAvgTables not yet available");
    }

    public SolverLN reset() {
        return this;
    }

    public boolean converged(int it) {
        /* Apply convergence test to SolverLN iterations. As the solver keeps iterating, this method maintains a
         * moving avg of the recent results based on which it averages across the layer the maximum queue-length
         * error. Convergence is tested by resetting all layers (to avoid caching) and doing an extra iteration.
         * If the iteration keeps fulfilling the error requirements for convergence, the solver completes.*/

        boolean bool = false;
        int iter_min = Math.min(2 * this.ensemble.length, (int) Math.ceil(this.options.iter_max / 4.0));
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
                        results.get(results.size()).get(e).QN.add(mov_avg_weight - 1,
                                results.get(results.size()).get(e).QN);
                        results.get(results.size()).get(e).UN.add(mov_avg_weight - 1,
                                results.get(results.size()).get(e).UN);
                        results.get(results.size()).get(e).RN.add(mov_avg_weight - 1,
                                results.get(results.size()).get(e).RN);
                        results.get(results.size()).get(e).TN.add(mov_avg_weight - 1,
                                results.get(results.size()).get(e).TN);
                        results.get(results.size()).get(e).AN.add(mov_avg_weight - 1,
                                results.get(results.size()).get(e).AN);
                        results.get(results.size()).get(e).WN.add(mov_avg_weight - 1,
                                results.get(results.size()).get(e).WN);

                        for (int k = 1; k < wnd_size; k++) {
                            results.get(results.size()).get(e).QN.add(mov_avg_weight,
                                    results.get(results.size() - k).get(e).QN);
                            results.get(results.size()).get(e).UN.add(mov_avg_weight,
                                    results.get(results.size() - k).get(e).UN);
                            results.get(results.size()).get(e).RN.add(mov_avg_weight,
                                    results.get(results.size() - k).get(e).RN);
                            results.get(results.size()).get(e).TN.add(mov_avg_weight,
                                    results.get(results.size() - k).get(e).TN);
                            results.get(results.size()).get(e).AN.add(mov_avg_weight,
                                    results.get(results.size() - k).get(e).AN);
                            results.get(results.size()).get(e).WN.add(mov_avg_weight,
                                    results.get(results.size() - k).get(e).WN);

                        }
                    }
                }
            }
        } else {
            int wnd_size = Integer.max(5, (int) Math.ceil(iter_min / 5));
            double mov_avg_weight = 1.0 / (double) wnd_size;
            results = this.results;
            // assume ready state
            if (it >= iter_min) {
                for (int e = 0; e < E; e++) {
                    results.get(results.size()).get(e).QN.add(mov_avg_weight - 1,
                            results.get(results.size()).get(e).QN);
                    results.get(results.size()).get(e).UN.add(mov_avg_weight - 1,
                            results.get(results.size()).get(e).UN);
                    results.get(results.size()).get(e).RN.add(mov_avg_weight - 1,
                            results.get(results.size()).get(e).RN);
                    results.get(results.size()).get(e).TN.add(mov_avg_weight - 1,
                            results.get(results.size()).get(e).TN);
                    results.get(results.size()).get(e).AN.add(mov_avg_weight - 1,
                            results.get(results.size()).get(e).AN);
                    results.get(results.size()).get(e).WN.add(mov_avg_weight - 1,
                            results.get(results.size()).get(e).WN);

                    for (int k = 1; k < wnd_size; k++) {
                        results.get(results.size()).get(e).QN.add(mov_avg_weight,
                                results.get(results.size() - k).get(e).QN);
                        results.get(results.size()).get(e).UN.add(mov_avg_weight,
                                results.get(results.size() - k).get(e).UN);
                        results.get(results.size()).get(e).RN.add(mov_avg_weight,
                                results.get(results.size() - k).get(e).RN);
                        results.get(results.size()).get(e).TN.add(mov_avg_weight,
                                results.get(results.size() - k).get(e).TN);
                        results.get(results.size()).get(e).AN.add(mov_avg_weight,
                                results.get(results.size() - k).get(e).AN);
                        results.get(results.size()).get(e).WN.add(mov_avg_weight,
                                results.get(results.size() - k).get(e).WN);

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
                int N = (int) this.ensemble[e].getNumberOfJobs().elementSum();
                if (N > 0) {
                    double IterErr;
                    try {
                        Matrix difference01 = metric.sub(1, metric_1);
                        difference01.removeNaN();
                        difference01.abs();
                        IterErr = (double) difference01.elementMax() / N;
                    } catch (Exception exception) {
                        IterErr = 0.0;
                    }
                    this.maxitererr.set(it, this.maxitererr.get(it) + IterErr);
                }
                if (this.options.verbose != VerboseLevel.SILENT) {
                    if (this.solvers[e].options.verbose != VerboseLevel.SILENT) {
                        String msg = String.format("\bQLen change: %.5f.\n", this.maxitererr.get(it) / E);
                        System.out.print(msg);
                    } else {
                        String msg = String.format("\bQLen change: %.5f.", this.maxitererr.get(it) / E);
                        System.out.print(msg);
                    }
                }
                if (it == iter_min) {
                    if (this.options.verbose != VerboseLevel.SILENT) {
                        System.out.print(" Starting averaging to help convergence.");
                    }
                    this.averagingstart = it;
                }
            }
        } else {
            if (it == 1 && (this.options.verbose != VerboseLevel.SILENT)) {
                System.out.println("");
            }
        }

        // Check convergence. Do not allow to converge in less than 2 iterations.
        if (it == 0 && (this.options.verbose != VerboseLevel.SILENT)) {
            System.out.println("SolverLN initialization completed. Starting iteration on ensemble models.");
        } else if ((it > 2) && (this.maxitererr.get(it) < this.options.iter_tol) && (this.maxitererr.get(it - 1)
                < this.options.iter_tol) && (this.maxitererr.get(it - 1) < this.options.iter_tol)) {
            // if potential convergence has just been detected, do a hard reset of every layer to check that this is
            // really the fixed point
            if (!this.hasconverged) {
                for (int e = 0; e < E; e++) {
                    this.ensemble[e].reset(false);
                }
                if (this.options.verbose != VerboseLevel.SILENT) {
                    if (this.solvers[this.solvers.length - 1].options.verbose != VerboseLevel.SILENT) {
                        String msg = " Testing convergence.";
                        System.out.print(msg);
                    } else {
                        String msg = "\b Testing convergence.";
                        System.out.println(msg);
                    }
                }
                //If it passes the change again next time then complete
                this.hasconverged = true;
            } else {
                if (this.options.verbose != VerboseLevel.SILENT) {
                    if (this.solvers[this.solvers.length - 1].options.verbose != VerboseLevel.SILENT) {
                        String msg = String.format(" SolverLN completed in %d iterations.", results.get(1).size());
                        System.out.println(msg);
                    } else {
                        String msg = String.format("\n SolverLN completed in %d iterations.", results.get(1).size());
                        System.out.println(msg);
                    }
                }
                bool = true;
            }
        } else {
            this.hasconverged = false;
        }
        return bool;
    }

    public void init() {
        //operation before starting to iterate

        List<Double> numSet = new ArrayList<Double>();
        if (this.route_prob_updmap.nz_length == 0) {
            this.unique_route_prob_updmap = this.route_prob_updmap;
        } else {
            this.unique_route_prob_updmap = this.route_prob_updmap.uniqueInCol(0);
//            for (int i = 1; i < this.route_prob_updmap.numRows; i++) {
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

        this.tput = new Matrix(this.lqn.nidx, 1, this.lqn.nidx);
        this.util = new Matrix(this.lqn.nidx, 1, this.lqn.nidx);
        this.servt = new Matrix(this.lqn.nidx, 1, this.lqn.nidx);
        this.servtmatrix = this.getEntryServiceMatrix();
        for (int e = 0; e < this.nlayers; e++)
            this.solvers[e].enableChecks = false;
    }

    public void pre(int it) {
        //no op
    }

    @Override
    public SolverResult analyze(int it, int e) {
        SolverResult result1 = new SolverResult();
        solvers[e].getAvg();

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
            solvers[e - 1].resetResults();
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
            ensemble[e - 1].refreshRates(statSet, classSet);
            solvers[e - 1].resetResults();
        }

        if (it == 1) {
            for (int e = 0; e < ensemble.length; e++) {
                solvers[e].enableChecks = false;
            }
        }

    }

    public void finish() {
        if (this.options.verbose != VerboseLevel.SILENT) {
            System.out.println("");
        }
        for(int e=0;e < this.ensemble.length;e++){
            solvers[e].getAvgTable();
        }
        //this.model.ensemble = this.ensemble; // not included as this comes through Diamond inheritance
    }

    @Override
    public boolean supports(Ensemble ensemble) {
        boolean bool = true;
        for(int e=0;e < this.ensemble.length;e++){
            //bool = bool && ((NetworkSolver)solvers[e]).supports((LayeredNetwork) ensemble.getModel(e));
        }
        return bool;
    }

    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.LN);
    }

    //GC
    @Override
    public void getEnsembleAvg() {

        this.iterate();
        // NOTE: TestSolverLN, TestSolverLN2, TestSolverLN3, had problems here due to
        // different values returned by getAvg() in MATLAB and LINE on these examples

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

        int E = this.nlayers;

        for (int e = 0; e < E; e++) {
            int clientIdx = this.ensemble[e].getAttribute().getClientIdx();
            int serverIdx = this.ensemble[e].getAttribute().getServerIdx();
            int sourceIdx = this.ensemble[e].getAttribute().getSourceIdx();

            // determine processor metrics
            Station s = this.ensemble[e].getStations().get(serverIdx - 1);
            Queue q = (Queue) s;
            if (q.getAttribute().getIsHost()) {
                int hidx = q.getAttribute().getIdx();
                TN.set(0, hidx, 0);
                PN.set(0, hidx, 0);
                for (int c = 0; c < this.ensemble[e].getNumberOfClasses(); c++) {
                    if (this.ensemble[e].getClasses().get(c).getCompletes()) {
                        int t = 0;
                        int u = 0;
                        if (clientIdx != -1) {
                            t = (int) (t > this.results.get(this.results.size()).get(e).TN.get(clientIdx - 1, c) ? t :
                                    this.results.get(results.size()).get(e).TN.get(clientIdx - 1, c));
                        }
                        if (sourceIdx != -1) {
                            t = (int) (t > this.results.get(this.results.size()).get(e).TN.get(clientIdx - 1, c) ? t :
                                    this.results.get(results.size()).get(e).TN.get(sourceIdx - 1, c));
                        }
                        double alter = this.results.get(this.results.size()).get(e).TN.get(serverIdx - 1, c);
                        if (!Double.isNaN(alter)) {
                            alter = alter > t ? alter : t;
                            TN.set(0, hidx, TN.get(hidx) + alter);
                        } else {
                            TN.set(0, hidx, TN.get(hidx) + t);
                        }
                    }
                    int type = this.ensemble[e].getClasses().get(c).getAttribute()[0];
                    if (type == LayeredNetworkElement.ACTIVITY) {
                        int aidx = this.ensemble[e].getClasses().get(c).getAttribute()[1];
                        int tidx = (int) this.lqn.parent.get(aidx);
                        if (Double.isNaN(PN.get(aidx)))
                            PN.set(0, aidx, 0);
                        if (Double.isNaN(PN.get(tidx)))
                            PN.set(0, tidx, 0);
                        PN.set(0, aidx, PN.get(aidx) + this.results.get(this.results.size()).get(e)
                                .UN.get(serverIdx - 1, c));
                        PN.set(0, tidx, PN.get(tidx) + this.results.get(this.results.size()).get(e)
                                .UN.get(serverIdx - 1, c));
                        PN.set(0, hidx, PN.get(hidx) + this.results.get(this.results.size()).get(e)
                                .UN.get(serverIdx - 1, c));
                    }
                }
            }

            //determine remaining metrics
            for (int c = 0; c < this.ensemble[e].getNumberOfClasses(); c++) {
                int type = this.ensemble[e].getClasses().get(c).getAttribute()[0];
                switch (type) {
                    case LayeredNetworkElement.TASK:
                        int tidx = this.ensemble[e].getClasses().get(c).getAttribute()[1];
                        Station task_s = this.ensemble[e].getStations().get(serverIdx - 1);
                        Queue task_q = (Queue) task_s;
                        if (task_q.getAttribute().getIsHost()) {
                            if (Double.isNaN(TN.get(tidx))) {
                                // store the result in th eprocessor model
                                TN.set(tidx, this.results.get(this.results.size()).get(e).TN.get(clientIdx - 1, c));
                            }
                        }
                        break;

                    case LayeredNetworkElement.ENTRY:
                        int eidx = this.ensemble[e].getClasses().get(c).getAttribute()[1];
                        tidx = (int) this.lqn.parent.get(eidx);  //unused parameter
                        Station entry_s = this.ensemble[e].getStations().get(serverIdx - 1);
                        Queue entry_q = (Queue) entry_s;
                        if (entry_q.getAttribute().getIsHost()) {
                            if (Double.isNaN(TN.get(eidx))) {
                                // store the result in th eprocessor model
                                TN.set(eidx, this.results.get(this.results.size()).get(e).TN.get(clientIdx - 1, c));
                            }
                        }
                        break;

                    case LayeredNetworkElement.CALL:
                        int cidx = this.ensemble[e].getClasses().get(c).getAttribute()[1];
                        int aidx = (int) this.lqn.callpair.get(cidx, 1);
                        SN.set(aidx, SN.get(aidx) + this.results.get(this.results.size()).get(e).RN.get(serverIdx - 1, c)
                                * this.lqn.callproc.get(cidx).getMean());
                        if (Double.isNaN(QN.get(aidx))) {
                            QN.set(aidx, 0);
                        }
                        QN.set(aidx, QN.get(aidx) + this.results.get(this.results.size()).get(e).QN.get(serverIdx - 1, c));
                        break;

                    case LayeredNetworkElement.ACTIVITY:

                        aidx = this.ensemble[e].getClasses().get(c).getAttribute()[1];
                        tidx = (int) this.lqn.parent.get(aidx);
                        QN.set(tidx, QN.get(tidx) + this.results.get(this.results.size()).get(e).QN.get(serverIdx - 1, c));

                        if (Double.isNaN(TN.get(aidx))) {
                            TN.set(aidx, 0);
                        }
                        if (Double.isNaN(QN.get(aidx))) {
                            QN.set(aidx, 0);
                        }

                        switch (this.ensemble[e].getClasses().get(c).getJobClassType()) {
                            case Closed:
                                TN.set(aidx, TN.get(aidx) + this.results.get(this.results.size()).get(e).TN.get(serverIdx - 1, c));
                                break;

                            case Open:
                                TN.set(aidx, TN.get(aidx) + this.results.get(this.results.size()).get(e).TN.get(sourceIdx - 1, c));
                                break;

                            default:
                        }
                        if (Double.isNaN(SN.get(aidx))) {
                            SN.set(aidx, 0);
                        }
                        SN.set(aidx, SN.get(aidx) + this.results.get(this.results.size()).get(e).RN.get(serverIdx - 1, c));

                        if (Double.isNaN(RN.get(aidx))) {
                            RN.set(aidx, 0);
                        }
                        RN.set(aidx, RN.get(aidx) + this.results.get(this.results.size()).get(e).RN.get(serverIdx - 1, c));

                        if (Double.isNaN(WN.get(aidx))) {
                            WN.set(aidx, 0);
                        }

                        if (Double.isNaN(WN.get(tidx))) {
                            WN.set(tidx, 0);
                        }
                        // WN be dealed as RN
                        WN.set(aidx, WN.get(aidx) + this.results.get(this.results.size()).get(e).RN.get(serverIdx - 1, c));
                        WN.set(tidx, WN.get(tidx) + this.results.get(this.results.size()).get(e).RN.get(serverIdx - 1, c));
                        if (Double.isNaN(QN.get(aidx))) {
                            QN.set(aidx, 0);
                        }
                        QN.set(aidx, QN.get(aidx) + this.results.get(this.results.size()).get(e).RN.get(serverIdx - 1, c));
                        break;
                    default:
                }
            }
        }
        for (int e = 1; e <= this.lqn.nentries; e++) {
            int eidx = this.lqn.eshift + e;
            int tidx = (int) this.lqn.parent.get(eidx);
            if (Double.isNaN(UN.get(tidx)))
                UN.set(tidx, 0);
            UN.set(eidx, TN.get(eidx) * SN.get(eidx));

            for (int i = 0; i < this.lqn.actsof.get(tidx).size(); i++) {
                int aidx = this.lqn.actsof.get(tidx).get(i);
                UN.set(aidx, TN.get(aidx) * SN.get(aidx));
            }
            UN.set(tidx, UN.get(tidx) + UN.get(eidx));
        }

        QN = UN.clone();
        UN = PN.clone();
        RN = SN.clone();

        int maxnamelength = 12;
        for (int i = 1; i <= lqn.names.size(); i++) {
            maxnamelength = Math.max(maxnamelength, lqn.names.get(i).length());
        }
        maxnamelength += 4;
        String Node = "Node";
        String NodeType = "NodeType";
        String nodeType = "NodeType";
        String Qlen = "QLen";
        String Util = "Util";
        String RespT = "RespT";
        String ResidT = "ResidT";
        String Tput = "Tput";
        String Processor = "Processor";
        String Task = "Task";
        String Entry = "Entry";
        String Activity = "Activitity";
        String format = "%-" + maxnamelength + "s%-16s%-16s%-16s%-16s%-16s%-16s\n";
        System.out.println("-------------------------------------------------------------------------------------------------------");
        System.out.format(format, Node, NodeType, Qlen, Util, RespT, ResidT, Tput);
        System.out.println("-------------------------------------------------------------------------------------------------------");
        for (int h = 1; h <= lqn.nhosts; h++) {
            System.out.format(format, lqn.names.get(h), Processor, String.format("%.4f", QN.get(h)), String.format("%.4f", UN.get(h)),
                    String.format("%.4f", RN.get(h)), String.format("%.4f", WN.get(h)),
                    String.format("%.4f", TN.get(h)));
        }
        for (int t = 1; t <= lqn.ntasks; t++) {
            System.out.format(format, lqn.names.get(t + lqn.tshift), Task, String.format("%.4f", QN.get(t + lqn.tshift)), String.format("%.4f", UN.get(t + lqn.tshift)),
                    String.format("%.4f", RN.get(t + lqn.tshift)), String.format("%.4f", WN.get(t + lqn.tshift)),
                    String.format("%.4f", TN.get(t + lqn.tshift)));
        }
        for (int e = 1; e <= lqn.nentries; e++) {
            System.out.format(format, lqn.names.get(e + lqn.eshift), Entry, String.format("%.4f", QN.get(e + lqn.eshift)), String.format("%.4f", UN.get(e + lqn.eshift)),
                    String.format("%.4f", RN.get(e + lqn.eshift)), String.format("%.4f", WN.get(e + lqn.eshift)),
                    String.format("%.4f", TN.get(e + lqn.eshift)));
        }
        for (int a = 1; a <= lqn.nacts; a++) {
            System.out.format(format, lqn.names.get(a + lqn.ashift), Activity, String.format("%.4f", QN.get(a + lqn.ashift)), String.format("%.4f", UN.get(a + lqn.ashift)),
                    String.format("%.4f", RN.get(a + lqn.ashift)), String.format("%.4f", WN.get(a + lqn.ashift)),
                    String.format("%.4f", TN.get(a + lqn.ashift)));
        }
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
            List<Integer> callers = lqn.tasksof.get(hidx);
            buildLayersRecursive(hidx, callers, true);
            this.idxhash.add(temp_idxhash);
            temp_idxhash++;
        }

        // build one subnetwork for every task
        for (int t = 1; t <= lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            boolean isolated_task = true;
            for (int i = 1; i <= lqn.nidx; i++) {
                if (lqn.iscaller.isAssigned(tidx, i) || lqn.iscaller.isAssigned(i, tidx)) {
                    isolated_task = false;
                    break;
                }
            }
            if ((int) lqn.isref.get(tidx) == 0 && !isolated_task) {
                // obtain the activity graph of each task that calls some entry in t
                List<Integer> callers = new ArrayList<>();
                for (int eidx : lqn.entriesof.get(tidx)) {
                    for (int row = lqn.tshift + 1; row < lqn.tshift + lqn.ntasks + 1; row++) {
                        if (lqn.iscaller.get(row, eidx) > 0) {
                            callers.add(row);
                        }
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
        thinkt_classes_updmap = cellToMatrix(cell_thinkt_classes_updmap);
        call_classes_updmap = cellToMatrix(cell_call_classes_updmap);
        servt_classes_updmap = cellToMatrix(cell_servt_classes_updmap);
        arvproc_classes_updmap = cellToMatrix(cell_arvproc_classes_updmap);
        route_prob_updmap = cellToMatrix(cell_route_prob_updmap);

        this.ensemble = new Network[nlayers];

        for (int i = 0; i < nlayers; i++) {
            ensemble[i] = temp_ensemble.get(i);
        }


    }

    public void buildLayersRecursive(int receiver_index, List<Integer> caller_index, boolean is_processor_layer) {
        nlayers++;
        Matrix jobPosKey = new Matrix(1, lqn.nidx + 1);
        Map<Integer, JobClass> curClassKey = new HashMap<>(lqn.nidx);
        int nreplicas = (int) lqn.repl.get(0, receiver_index);
        Network model = new Network(lqn.hashnames.get(receiver_index));
        model.setDoChecks(false);

        boolean hasSynccaller = false;
        if (!is_processor_layer) {
            for (int callers : caller_index) {
                for (int entries : lqn.entriesof.get(receiver_index)) {
                    if (lqn.issynccaller.isAssigned(callers, entries)) {
                        hasSynccaller = true;
                    }
                }
            }
        }
        Delay clientDelay = null;
        if (is_processor_layer || hasSynccaller) {
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
                serverStation.put(i, new Queue(model, lqn.hashnames.get(receiver_index), lqn.sched.get(receiver_index)));
            } else {
                String name = lqn.hashnames.get(receiver_index).concat(".");
                serverStation.put(i, new Queue(model, name.concat(Integer.toString(i)), lqn.sched.get(receiver_index)));
            }
            serverStation.get(i).setNumberOfServers((int) lqn.mult.get(0, receiver_index));
            serverStation.get(i).getAttribute().setIsHost(is_processor_layer);
            serverStation.get(i).getAttribute().setIdx(receiver_index);
        }

        boolean iscachelayer = false;
        // todo cache line37-40
        /*if (is_processor_layer) {
            for (int callers : caller_index) {
                if (lqn.iscache.get(0, callers) != 0) {
                    iscachelayer = true;
                    break;
                }
            }
        }
        if (iscachelayer) {

            // Cache cacheNode = new Cache(model,lqn.hashnames.get(receiver_index),lqn.itemlevelcap.get(caller_index));

        }*/

        List<Integer> actsInCaller = new ArrayList<>();
        for (int i : caller_index) {
            actsInCaller.addAll(lqn.actsof.get(i));
        }


        Matrix isPostAndAct = new Matrix(lqn.actposttype.numRows, lqn.actposttype.numCols, lqn.nacts);
        Matrix isPreAndAct = new Matrix(lqn.actpretype.numRows, lqn.actpretype.numCols, lqn.actpretype.nz_length);
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
            for (int i = 1; i < lqn.graph.numCols; i++) {
                if (lqn.graph.get(aidx, i) != 0) {
                    successors.add(i);
                }
            }
            boolean flag = true;
            for (int i : successors) {
                if (isPostAndAct.get(0, i) == 0) {
                    flag = false;
                }
            }
            if (flag) {
                maxfanout = Math.max(maxfanout, successors.size());
            }
        }


        Fork forkNode = null;
        Map<Integer, Router> forkOutputRouter = new HashMap<>(maxfanout);
        if (hasFork) {
            forkNode = new Fork(model);
            for (int f = 1; f <= maxfanout; f++) {
                forkOutputRouter.put(f, new Router(model, "forkNode_router" + f));
            }
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
            joinNode = new Join(model);
        }

        Map<Integer, JobClass> aidxclass = new HashMap<>(lqn.nentries + lqn.nacts);
        Map<Integer, JobClass> cidxclass = new HashMap<>();
        Map<Integer, JobClass> cidxauxclass = new HashMap<>();

        if (is_processor_layer) {
            model.getAttribute().addHosts(new Integer[]{null, model.getAttribute().getServerIdx()});
        } else {
            model.getAttribute().addTasks(new Integer[]{null, model.getAttribute().getServerIdx()});
        }

        Source sourceStation = null;
        Sink sinkStation = null;
        Matrix openClasses = new Matrix(lqn.ncalls + 1, 4, lqn.ncalls * 3);
        int openClassesAssignedLine = 0;

        boolean issynccaller = false;
        if (!is_processor_layer) {
            for (int tidx_caller : caller_index) {
                for (int idx : lqn.entriesof.get(receiver_index)) {
                    if (lqn.issynccaller.get(tidx_caller, idx) != 0) {
                        issynccaller = true;
                        break;
                    }
                }
            }
        }
        double curnjobs;
        Map<Integer, Double> callmean = new HashMap<>();
        for (int tidx_caller : caller_index) {
            if (is_processor_layer || issynccaller) {
                if (njobs.get(tidx_caller, receiver_index) == 0) {
                    curnjobs = lqn.mult.get(0, tidx_caller) * lqn.repl.get(0, tidx_caller);
                    if (Double.isInfinite(curnjobs)) {
                        curnjobs = 0;
                        for (int i = 1; i < lqn.mult.numCols; i++) {
                            if (Double.isFinite(lqn.mult.get(0, i))) {
                                curnjobs = curnjobs + lqn.mult.get(0, i) * lqn.repl.get(0, i);
                            }
                        }
                        curnjobs = Math.min(curnjobs, 1000000);
                    }
                    njobs.set(tidx_caller, receiver_index, curnjobs);
                } else {
                    curnjobs = njobs.get(tidx_caller, receiver_index);
                }
                String caller_name = lqn.hashnames.get(tidx_caller);
                aidxclass.put(tidx_caller, new ClosedClass(model, caller_name, (long) curnjobs, clientDelay));
                aidxclass.get(tidx_caller).setReferenceClass(true);
                aidxclass.get(tidx_caller).setAttribute(new Integer[]{LayeredNetworkElement.TASK, tidx_caller});
                aidxclass.get(tidx_caller).setCompletes(false);
                model.getAttribute().addTasks(new Integer[]{aidxclass.get(tidx_caller).getIndex(), tidx_caller});
                assert clientDelay != null;
                clientDelay.setService(aidxclass.get(tidx_caller), thinkproc.get(tidx_caller));
                if (lqn.isref.get(tidx_caller) == 0) {
                    if (!cell_thinkt_classes_updmap.containsKey(receiver_index)) {
                        cell_thinkt_classes_updmap.put(receiver_index, new ArrayList<>());
                    }
                    cell_thinkt_classes_updmap.get(receiver_index).add(new Integer[]{receiver_index, tidx_caller, 1, aidxclass.get(tidx_caller).getIndex()});
                }

                for (int eidx : lqn.entriesof.get(tidx_caller)) {
                    aidxclass.put(eidx, new ClosedClass(model, lqn.hashnames.get(eidx), 0, clientDelay));
                    aidxclass.get(eidx).setCompletes(false);
                    aidxclass.get(eidx).setAttribute(new Integer[]{LayeredNetworkElement.ENTRY, eidx});
                    model.getAttribute().addEntries(new Integer[]{aidxclass.get(eidx).getIndex(), eidx});
                    clientDelay.setService(aidxclass.get(eidx), new Immediate());
                }
            }

            for (int aidx : lqn.actsof.get(tidx_caller)) {
                if (is_processor_layer || issynccaller) {
                    aidxclass.put(aidx, new ClosedClass(model, lqn.hashnames.get(aidx), 0, clientDelay));
                    aidxclass.get(aidx).setCompletes(false);
                    aidxclass.get(aidx).setAttribute(new Integer[]{LayeredNetworkElement.ACTIVITY, aidx});
                    model.getAttribute().addActivities(new Integer[]{aidxclass.get(aidx).getIndex(), aidx});
                    if (!(is_processor_layer && (lqn.parent.get(0, (int) lqn.parent.get(0, aidx)) == receiver_index))) {
                        clientDelay.setService(aidxclass.get(aidx), servtproc.get(aidx));
                    }
                    if (iscachelayer) {
                        //todo cache
                    }
                }

                for (int cidx : lqn.callsof.get(aidx)) {
                    callmean.put(cidx, lqn.callproc.get(cidx).getMean());
                    if (lqn.calltype.get(cidx) == CallType.ASYNC) {
                        if (lqn.parent.get(0, (int) lqn.callpair.get(cidx, 2)) == receiver_index) {
                            if (sourceStation == null) {
                                model.getAttribute().setClientIdx(model.getNumberOfNodes() + 1);
                                sourceStation = new Source(model, "Source");
                                sinkStation = new Sink(model, "Sink");
                            }
                            cidxclass.put(cidx, new OpenClass(model, lqn.callhashnames.get(cidx), 0));
                            sourceStation.setArrival(cidxclass.get(cidx), new Exp(GlobalConstants.CoarseTol));
                            for (int m = 1; m <= nreplicas; m++) {
                                serverStation.get(m).setService(cidxclass.get(cidx), new Immediate());
                            }
                            openClassesAssignedLine++;
                            openClasses.set(openClassesAssignedLine, 1, cidxauxclass.get(cidx).getIndex());
                            openClasses.set(openClassesAssignedLine, 2, lqn.callproc.get(cidx).getMean());
                            openClasses.set(openClassesAssignedLine, 3, cidx);
                            model.getAttribute().addCalls(new Integer[]{cidxclass.get(cidx).getIndex(), cidx, (int) lqn.callpair.get(cidx, 1), (int) lqn.callpair.get(cidx, 2)});
                            //cidxclass.get(cidx).setCompletes(false);//todo check cidx or aidx
                            cidxclass.get(cidx).setAttribute(new Integer[]{LayeredNetworkElement.CALL, cidx});
                            double minRespT = 0;
                            if (!is_processor_layer) {
                                for (int tidx_act : lqn.actsof.get(receiver_index)) {
                                    minRespT = minRespT + lqn.hostdem.get(tidx_act).getMean();
                                }
                            }
                            for (int m = 1; m <= nreplicas; m++) {
                                serverStation.get(m).setService(cidxclass.get(cidx), new Exp(1 / minRespT));
                            }
                        }
                    } else if (lqn.calltype.get(cidx) == CallType.SYNC) {
                        cidxclass.put(cidx, new ClosedClass(model, lqn.callhashnames.get(cidx), 0, clientDelay));
                        //cidxclass.get(cidx).setCompletes(false);//todo check cidx or aidx
                        cidxclass.get(cidx).setAttribute(new Integer[]{LayeredNetworkElement.CALL, cidx});
                        double minRespT = 0;
                        if (!is_processor_layer) {
                            for (int tidx_act : lqn.actsof.get(receiver_index)) {
                                minRespT = minRespT + lqn.hostdem.get(tidx_act).getMean();
                            }
                        }
                        for (int m = 1; m <= nreplicas; m++) {
                            serverStation.get(m).setService(cidxclass.get(cidx), new Exp(1 / minRespT));
                        }
                    }

                    if (callmean.get(cidx) != nreplicas) {
                        if (lqn.calltype.get(cidx) == CallType.SYNC) {
                            cidxauxclass.put(cidx, new ClosedClass(model, lqn.callhashnames.get(cidx) + ".Aux", 0, clientDelay));
                            cidxauxclass.get(cidx).setCompletes(false);
                            cidxauxclass.get(cidx).setAttribute(new Integer[]{LayeredNetworkElement.CALL, cidx});
                            assert clientDelay != null;
                            clientDelay.setService(cidxauxclass.get(cidx), new Immediate());
                            for (int m = 1; m <= nreplicas; m++) {
                                serverStation.get(m).setService(cidxauxclass.get(cidx), new Disabled());
                            }
                        }
                    }
                }
            }
        }

        RoutingMatrix P = new RoutingMatrix(model, model.getJobClass(), model.getNodes());
        if (sourceStation != null) {
            for (int o = 1; o <= openClassesAssignedLine; o++) {
                int oidx = (int) openClasses.get(o, 1);
                double p = 1 / openClasses.get(o, 2);
                for (int m = 1; m <= nreplicas; m++) {
                    P.addConnection(model.getClasses().get(oidx), model.getClasses().get(oidx), sourceStation, serverStation.get(m), 1.0 / (double) nreplicas);
                    for (int n = 1; n <= nreplicas; n++) {
                        P.addConnection(model.getClasses().get(oidx), model.getClasses().get(oidx), serverStation.get(m), serverStation.get(n), (1.0 - p) / (double) nreplicas);
                    }
                    P.addConnection(model.getClasses().get(oidx), model.getClasses().get(oidx), serverStation.get(m), sinkStation, p);
                }
                int cidx = (int) openClasses.get(o, 3);
                if (!cell_arvproc_classes_updmap.containsKey(receiver_index)) {
                    cell_arvproc_classes_updmap.put(receiver_index, new ArrayList<>());
                }
                cell_arvproc_classes_updmap.get(receiver_index).add(new Integer[]{receiver_index, cidx, model.getNodeIndex(sourceStation) + 1, oidx});
                for (int m = 1; m <= nreplicas; m++) {
                    if (!cell_call_classes_updmap.containsKey(receiver_index)) {
                        cell_call_classes_updmap.put(receiver_index, new ArrayList<>());
                    }
                    cell_call_classes_updmap.get(receiver_index).add(new Integer[]{receiver_index, cidx, model.getNodeIndex(serverStation.get(m)) + 1, oidx});
                }
            }
        }

        int atClient = 1;
        int atServer = 2;
        //int atCache = 3; todo cache

        int jobPos = atClient;
        class Inner {

            recurActGraphReturnType recirActGraph(RoutingMatrix P, int tidx_caller, int aidx, JobClass curClass, int jobPos, Source sourceStation, Delay clientDelay, Sink sinkStation, Join joinNode, Fork forkNode) {
                jobPosKey.set(0, aidx, jobPos);
                curClassKey.put(aidx, curClass);
                List<Integer> nextaidxs = new ArrayList<>();
                for (int i = 1; i < lqn.graph.numCols; i++) {
                    if (lqn.graph.isAssigned(aidx, i)) {
                        nextaidxs.add(i);
                    }
                }
                Matrix isNextPrecFork = new Matrix(1, lqn.nidx + 1, lqn.nidx);

                if (!nextaidxs.isEmpty()) {
                    isNextPrecFork.set(0, aidx, 1);
                    for (int i : nextaidxs) {
                        if (!isPostAndAct.isAssigned(0, i)) {
                            isNextPrecFork.set(0, aidx, 0);
                            break;
                        }
                    }
                }
                if (!nextaidxs.isEmpty()) {
                    for (int nextaidx : nextaidxs) {
                        if (lqn.parent.get(0, aidx) != lqn.parent.get(0, nextaidx)) {
                            int cidx = 0;
                            for (int i = 1; i <= lqn.ncalls; i++) {
                                if (lqn.callpair.get(i, 1) == aidx && lqn.callpair.get(i, 2) == nextaidx) {
                                    cidx = i;
                                    break;
                                }
                            }
                            if (lqn.calltype.get(cidx) == CallType.SYNC) {
                                if (jobPos == atClient) {
                                    if (lqn.parent.get((int) lqn.callpair.get(cidx, 2)) == receiver_index) {
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
                                        } else {
                                            for (int m = 1; m < nreplicas; m++) {
                                                P.addConnection(curClass, cidxclass.get(cidx), clientDelay, serverStation.get(m), 1.0 / (double) nreplicas);
                                                P.addConnection(cidxclass.get(cidx), cidxauxclass.get(cidx), serverStation.get(m), clientDelay, 1.0);
                                                P.addConnection(cidxauxclass.get(cidx), cidxclass.get(cidx), clientDelay, serverStation.get(m), 1.0 - 1.0 / (callmean.get(cidx) / nreplicas));
                                            }
                                            P.addConnection(clientDelay, clientDelay, cidxauxclass.get(cidx), cidxclass.get(cidx));
                                        }
                                        clientDelay.setService(cidxclass.get(cidx), new Immediate());
                                        if (!cell_call_classes_updmap.containsKey(receiver_index)) {
                                            cell_call_classes_updmap.put(receiver_index, new ArrayList<>());
                                        }
                                        for (int m = 1; m <= nreplicas; m++) {
                                            serverStation.get(m).setService(cidxclass.get(cidx), callresidtproc.get(cidx));
                                            cell_call_classes_updmap.get(receiver_index).add(new Integer[]{receiver_index, cidx, model.getNodeIndex(serverStation.get(m)) + 1, cidxclass.get(cidx).getIndex()});
                                        }
                                        curClass = cidxclass.get(cidx);
                                    } else {
                                        if (callmean.get(cidx) < nreplicas) {
                                            P.addConnection(curClass, cidxclass.get(cidx), clientDelay, clientDelay, callmean.get(cidx) / nreplicas);
                                            P.addConnection(cidxclass.get(cidx), cidxauxclass.get(cidx), clientDelay, clientDelay, 1.0);
                                            P.addConnection(cidxclass.get(cidx), cidxclass.get(cidx), clientDelay, clientDelay, 1.0);
                                            curClass = cidxauxclass.get(cidx);
                                        } else if (callmean.get(cidx) == nreplicas) {
                                            P.addConnection(curClass, cidxclass.get(cidx), clientDelay, clientDelay, 1.0);
                                            curClass = cidxclass.get(cidx);
                                        } else {
                                            P.addConnection(curClass, cidxclass.get(cidx), clientDelay, clientDelay, 1.0);
                                            P.addConnection(cidxclass.get(cidx), cidxauxclass.get(cidx), clientDelay, clientDelay, 1.0);
                                            curClass = cidxauxclass.get(cidx);
                                        }
                                        jobPos = atClient;
                                        clientDelay.setService(cidxclass.get(cidx), callresidtproc.get(cidx));
                                        if (!cell_call_classes_updmap.containsKey(receiver_index)) {
                                            cell_call_classes_updmap.put(receiver_index, new ArrayList<>());
                                        }
                                        cell_call_classes_updmap.get(receiver_index).add(new Integer[]{receiver_index, cidx, 1, cidxclass.get(cidx).getIndex()});
                                    }
                                } else if (jobPos == atServer) {
                                    if (lqn.parent.get((int) lqn.callpair.get(cidx, 2)) == receiver_index) {
                                        if (callmean.get(cidx) < nreplicas) {
                                            for (int m = 1; m <= nreplicas; m++) {
                                                P.addConnection(curClass, cidxclass.get(cidx), serverStation.get(m), clientDelay, 1.0 - callmean.get(cidx));
                                                P.addConnection(curClass, cidxclass.get(cidx), serverStation.get(m), serverStation.get(m), callmean.get(cidx));
                                                serverStation.get(m).setService(cidxclass.get(cidx), callresidtproc.get(cidx));
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
                                        if (!cell_call_classes_updmap.containsKey(receiver_index)) {
                                            cell_call_classes_updmap.put(receiver_index, new ArrayList<>());
                                        }
                                        for (int m = 1; m <= nreplicas; m++) {
                                            serverStation.get(m).setService(cidxclass.get(cidx), callresidtproc.get(cidx));
                                            cell_call_classes_updmap.get(receiver_index).add(new Integer[]{receiver_index, cidx, model.getNodeIndex(serverStation.get(m)) + 1, cidxclass.get(cidx).getIndex()});
                                        }
                                    } else {
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
                                        } else {
                                            for (int m = 1; m <= nreplicas; m++) {
                                                P.addConnection(curClass, cidxclass.get(cidx), serverStation.get(m), clientDelay, 1.0);
                                            }
                                            P.addConnection(cidxclass.get(cidx), cidxauxclass.get(cidx), clientDelay, clientDelay, 1.0);
                                            curClass = cidxauxclass.get(cidx);
                                        }
                                        jobPos = atClient;
                                        clientDelay.setService(cidxclass.get(cidx), callresidtproc.get(cidx));
                                        if (!cell_call_classes_updmap.containsKey(receiver_index)) {
                                            cell_call_classes_updmap.put(receiver_index, new ArrayList<>());
                                        }
                                        cell_call_classes_updmap.get(receiver_index).add(new Integer[]{receiver_index, cidx, 1, cidxclass.get(cidx).getIndex()});
                                    }
                                }
                            }
                        } else {
                            if (nextaidx <= lqn.eshift || nextaidx > lqn.ashift) {
                                jobPos = (int) jobPosKey.get(0, aidx);
                                curClass = curClassKey.get(aidx);
                            } else {
                                for (int i = 0; i < nextaidxs.size(); i++) {
                                    if (nextaidxs.get(i) == nextaidx) {
                                        if ((nextaidxs.get(i - 1) >= (lqn.eshift + 1)) && (nextaidxs.get(i - 1) <= (lqn.eshift + lqn.nentries))) {
                                            curClassC = curClass;
                                        }
                                        break;
                                    }
                                }
                                jobPos = atClient;
                                curClass = curClassC;
                            }

                            if (jobPos == atClient) {
                                if (is_processor_layer) {
                                    //todo cache line 317-332
                                    for (int m = 1; m <= nreplicas; m++) {
                                        if (isNextPrecFork.get(0, aidx) != 0) {
                                            P.addConnection(curClass, curClass, clientDelay, forkNode, 1.0);
                                            int f = 0;
                                            for (int i = 0; i < nextaidxs.size(); i++) {
                                                if (nextaidxs.get(i) == nextaidx) {
                                                    f = i;
                                                }
                                            }
                                            P.addConnection(curClass, curClass, forkNode, forkOutputRouter.get(f), 1.0);
                                            P.addConnection(curClass, aidxclass.get(nextaidx), forkOutputRouter.get(f), serverStation.get(m), 1.0);
                                        } else {
                                            if (isPostAndAct.get(0, aidx) != 0) {
                                                P.addConnection(curClass, curClass, clientDelay, joinNode, 1.0);
                                                P.addConnection(curClass, aidxclass.get(nextaidx), joinNode, serverStation.get(m), 1.0);
                                            } else {
                                                P.addConnection(curClass, aidxclass.get(nextaidx), clientDelay, serverStation.get(m), lqn.graph.get(aidx, nextaidx));
                                            }
                                        }
                                        serverStation.get(m).setService(aidxclass.get(nextaidx), lqn.hostdem.get(nextaidx));
                                    }
                                    jobPos = atServer;
                                    curClass = aidxclass.get(nextaidx);
                                    if (!cell_servt_classes_updmap.containsKey(receiver_index)) {
                                        cell_servt_classes_updmap.put(receiver_index, new ArrayList<>());
                                    }
                                    cell_servt_classes_updmap.get(receiver_index).add(new Integer[]{receiver_index, nextaidx, 2, aidxclass.get(nextaidx).getIndex()});

                                } else {
                                    if (isNextPrecFork.get(0, aidx) != 0) {
                                        P.addConnection(curClass, curClass, clientDelay, forkNode, 1.0);
                                        int f = 0;
                                        for (int i = 0; i < nextaidxs.size(); i++) {
                                            if (nextaidxs.get(i) == nextaidx) {
                                                f = i;
                                            }
                                        }
                                        P.addConnection(curClass, curClass, forkNode, forkOutputRouter.get(f), 1.0);
                                        P.addConnection(curClass, aidxclass.get(nextaidx), forkOutputRouter.get(f), clientDelay, 1.0);
                                    } else {
                                        if (isPreAndAct.get(0, aidx) != 0) {
                                            P.addConnection(curClass, curClass, clientDelay, joinNode, 1.0);
                                            P.addConnection(joinNode, clientDelay, curClass, aidxclass.get(nextaidx));
                                        } else {
                                            P.addConnection(curClass, aidxclass.get(nextaidx), clientDelay, clientDelay, lqn.graph.get(aidx, nextaidx));
                                        }
                                    }
                                    jobPos = atClient;
                                    curClass = aidxclass.get(nextaidx);
                                    clientDelay.setService(aidxclass.get(nextaidx), servtproc.get(nextaidx));
                                    if (!cell_thinkt_classes_updmap.containsKey(receiver_index)) {
                                        cell_thinkt_classes_updmap.put(receiver_index, new ArrayList<>());
                                    }
                                    cell_thinkt_classes_updmap.get(receiver_index).add(new Integer[]{receiver_index, nextaidx, 1, aidxclass.get(nextaidx).getIndex()});
                                }
                            } else {
                                if (is_processor_layer) {
                                    //todo Cache line 355-375
                                    for (int m = 1; m <= nreplicas; m++) {
                                        if (isNextPrecFork.get(0, aidx) != 0) {
                                            P.addConnection(curClass, curClass, serverStation.get(m), forkNode, 1.0);
                                            int f = 0;
                                            for (int i = 0; i < nextaidxs.size(); i++) {
                                                if (nextaidxs.get(i) == nextaidx) {
                                                    f = i;
                                                }
                                            }
                                            P.addConnection(curClass, curClass, forkNode, forkOutputRouter.get(f), 1.0);
                                            P.addConnection(curClass, aidxclass.get(nextaidx), forkOutputRouter.get(f), serverStation.get(m), 1.0);
                                        } else {
                                            if (isPreAndAct.get(0, aidx) != 0) {
                                                P.addConnection(curClass, curClass, serverStation.get(m), joinNode, 1.0);
                                                P.addConnection(curClass, aidxclass.get(nextaidx), joinNode, serverStation.get(m), 1.0);
                                            } else {
                                                P.addConnection(curClass, aidxclass.get(nextaidx), serverStation.get(m), serverStation.get(m), lqn.graph.get(aidx, nextaidx));
                                            }
                                        }
                                        serverStation.get(m).setService(aidxclass.get(nextaidx), lqn.hostdem.get(nextaidx));
                                    }

                                    jobPos = atServer;
                                    curClass = aidxclass.get(nextaidx);
                                    if (!cell_servt_classes_updmap.containsKey(receiver_index)) {
                                        cell_servt_classes_updmap.put(receiver_index, new ArrayList<>());
                                    }
                                    cell_servt_classes_updmap.get(receiver_index).add(new Integer[]{receiver_index, nextaidx, 2, aidxclass.get(nextaidx).getIndex()});
                                } else {
                                    for (int m = 1; m < nreplicas; m++) {
                                        if (isNextPrecFork.get(0, aidx) != 0) {
                                            P.addConnection(curClass, curClass, serverStation.get(m), forkNode, 1.0);
                                            int f = 0;
                                            for (int i = 0; i < nextaidxs.size(); i++) {
                                                if (nextaidxs.get(i) == nextaidx) {
                                                    f = i;
                                                }
                                            }
                                            P.addConnection(curClass, curClass, forkNode, forkOutputRouter.get(f), 1.0);
                                            P.addConnection(curClass, aidxclass.get(nextaidx), forkOutputRouter.get(f), clientDelay, 1.0);
                                        } else {
                                            if (isPreAndAct.get(0, aidx) != 0) {
                                                P.addConnection(curClass, curClass, serverStation.get(m), joinNode, 1.0);
                                                P.addConnection(curClass, aidxclass.get(nextaidx), joinNode, clientDelay, 1.0);
                                            } else {
                                                P.addConnection(curClass, aidxclass.get(nextaidx), serverStation.get(m), clientDelay, lqn.graph.get(aidx, nextaidx));
                                            }
                                        }
                                        jobPos = atClient;
                                        curClass = aidxclass.get(nextaidx);
                                        clientDelay.setService(aidxclass.get(nextaidx), servtproc.get(nextaidx));
                                        if (!cell_thinkt_classes_updmap.containsKey(receiver_index)) {
                                            cell_thinkt_classes_updmap.put(receiver_index, new ArrayList<>());
                                        }
                                        cell_thinkt_classes_updmap.get(receiver_index).add(new Integer[]{receiver_index, nextaidx, 1, aidxclass.get(nextaidx).getIndex()});
                                    }
                                }
                            }
                            if (aidx != nextaidx) {
                                recurActGraphReturnType returnType = recirActGraph(P, tidx_caller, nextaidx, curClass, jobPos, sourceStation, clientDelay, sinkStation, joinNode, forkNode);
                                P = returnType.P;
                                curClass = returnType.curClass;
                                jobPos = returnType.jobPos;

                                if (jobPos == atClient) {
                                    P.addConnection(curClass, aidxclass.get(tidx_caller), clientDelay, clientDelay, 1.0);
                                    if (!curClass.getName().contains(".Aux")) {
                                        curClass.setCompletes(true);
                                    }
                                } else {
                                    for (int m = 1; m <= nreplicas; m++) {
                                        P.addConnection(curClass, aidxclass.get(tidx_caller), serverStation.get(m), clientDelay, 1.0);
                                    }
                                    if (!curClass.getName().contains(".Aux")) {
                                        curClass.setCompletes(true);
                                    }
                                }
                            }
                        }
                    }
                }
                return new recurActGraphReturnType(curClass, jobPos, P);
            }

        }


        for (int tidx_caller : caller_index) {
            if (lqn.issynccaller.get(tidx_caller, receiver_index) == 1 || is_processor_layer) {
                int ncaller_entries = lqn.entriesof.get(tidx_caller).size();
                for (int eidx : lqn.entriesof.get(tidx_caller)) {
                    JobClass aidxClass_eidx = aidxclass.get(eidx);
                    JobClass aidxClass_tidx_caller = aidxclass.get(tidx_caller);
                    P.addConnection(aidxClass_tidx_caller, aidxClass_eidx, clientDelay, clientDelay, 1.0 / (double) ncaller_entries);

                    if (ncaller_entries > 1) {
                        if (!cell_route_prob_updmap.containsKey(receiver_index)) {
                            cell_route_prob_updmap.put(receiver_index, new ArrayList<>());
                        }
                        cell_route_prob_updmap.get(receiver_index).add(new Integer[]{receiver_index, tidx_caller, eidx, 1, 1, aidxClass_tidx_caller.getIndex(), aidxClass_eidx.getIndex()});
                    }
                    P = new Inner().recirActGraph(P, tidx_caller, eidx, aidxClass_eidx, jobPos, sourceStation, clientDelay, sinkStation, joinNode, forkNode).P;
                }
            }
        }
        model.link(P);
        temp_ensemble.add(model);

    }


    public void updateLayers(int it) {
        //task14: updateLayers function to be written
        // reassign service times
        for (int r = 1; r < this.thinkt_classes_updmap.numRows; r++) {
            int ri;
            if (it != 0) {
                ri = this.thinkt_classes_updmap.numRows - r;
            } else {
                ri = r;
            }

            double idx = this.thinkt_classes_updmap.get(ri, 1);
            double aidx = this.thinkt_classes_updmap.get(ri, 2);
            double nodeidx = this.thinkt_classes_updmap.get(ri, 3);
            double classidx = this.thinkt_classes_updmap.get(ri, 4);
            JobClass tmp_class = this.ensemble[this.idxhash.get((int) idx).intValue() - 1].getClassByIndex((int) classidx - 1);
            // here update the number of jobs in the task chain
            if (aidx < (this.lqn.tshift + this.lqn.ntasks)) {
                // aidx here is actually set to tidx in buildLayersRecursive
                if (tmp_class.getJobClassType() == JobClassType.Closed) {
                    ClosedClass tmp_class_c = (ClosedClass) tmp_class;
                    tmp_class_c.setPopulation((long) this.njobs.get((int) aidx, (int) idx));
                }
            }
            Queue node = (Queue) this.ensemble[idxhash.get((int) idx).intValue() - 1].getNodeByStatefulIndex((int) nodeidx - 1);
            // Case 1
            if ((int) nodeidx == this.ensemble[(idxhash.get((int) idx).intValue()) - 1].getAttribute().getClientIdx()) {
                if (this.lqn.type.get((int) aidx) == LayeredNetworkElement.TASK) {
                    if (this.lqn.schedid.get((int) aidx) != SchedStrategy.toID(SchedStrategy.REF)) {
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
        for (int r = 1; r < this.arvproc_classes_updmap.numRows; r++) {
            int ri;
            if (it != 0) {
                ri = this.arvproc_classes_updmap.numRows - r;
            } else {
                ri = r;
            }
            double idx = this.arvproc_classes_updmap.get(ri, 1);
            double cidx = this.arvproc_classes_updmap.get(ri, 2);
            double nodeidx = this.arvproc_classes_updmap.get(ri, 3);
            double classidx = this.arvproc_classes_updmap.get(ri, 4);
            JobClass tmp_class = this.ensemble[this.idxhash.get((int) idx).intValue() - 1].getClassByIndex((int) classidx - 1);
            Source node = (Source) this.ensemble[this.idxhash.get((int) idx).intValue() - 1].getNodeByStatefulIndex((int) nodeidx - 1);
            node.setArrival(tmp_class, this.tputproc.get((int) this.lqn.callpair.get((int) cidx, 1)));
        }

        // reassign call service time / response time
        for (int c = 1; c < this.call_classes_updmap.numRows; c++) {
            int ci;
            if (it != 0) {
                ci = this.call_classes_updmap.numRows - c;
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
                node.setService(tmp_class, this.callresidtproc.get((int) cidx));
            }

            // Case 2 - the call is processed by the server, then replace with the svc time
            if ((int) nodeidx == this.ensemble[this.idxhash.get((int) idx).intValue() - 1].getAttribute().getServerIdx()) {
                double eidx = this.lqn.callpair.get((int) cidx, 2);
                node.setService(tmp_class, this.servtproc.get((int) eidx)); /* TODO update for correct assignment line 121 MATLAB */
            }
        }
    }

    public void updatePopulations(int it) {
        //task15: updatePopulations function to be written
        LayeredNetworkStruct lqn = this.lqn;
        // interlock scaling factors
        int ilscalingside = lqn.nhosts + lqn.ntasks;
        Matrix ilscaling = new Matrix(ilscalingside, ilscalingside, ilscalingside * ilscalingside);
        for (int i = 0; i < ilscalingside; i++) {
            for (int j = 0; j < ilscalingside; j++) {
                ilscaling.set(i, j, 1);
            }
        }

        for (int step = 1; step <= this.nlayers; step++) {
            for (int h = 1; h <= lqn.nhosts; h++) {
                int hidx = h;
                ilscaling.set(hidx, 1.0);
                if (lqn.isref.get(hidx) == 0) {
                    //the following are remote (indirect) callers that certain to be callers
                    //of task t, hence if they have multiplicity m ten task t cannot have as
                    //a matter of fact multiplicity more than m
                    List<Integer> callers = lqn.tasksof.get(hidx);
                    int multcallers = 0;
                    for (int i : callers) {
                        multcallers += this.njobsorig.get(i, hidx);
                    }
                    Matrix rowhidx = new Matrix(1, this.ptaskcallers_step.get(step).numCols + 1,
                            this.ptaskcallers_step.get(step).numCols);
                    Matrix.extractRows(this.ptaskcallers_step.get(step), hidx, hidx + 1, rowhidx);
                    Matrix step_callers = rowhidx.find();
                    int multremote = 0;
                    for (int i = 1; i < step_callers.length(); i++) {
                        int remidx = (int) step_callers.get(i);
                        if (lqn.schedid.get(remidx) == SchedStrategy.toID(SchedStrategy.INF)) {
                            multremote += this.ptaskcallers_step.get(step).get(hidx, remidx) * this.util.get(remidx);
                        } else {
                            multremote += this.ptaskcallers_step.get(step).get(hidx, remidx) * this.util.get(remidx) *
                                    lqn.mult.get(remidx);
                        }
                    }
                    if ((multcallers > multremote && multremote > 0) && !Double.isInfinite(multremote)) {
                        // we spread the scaling proportionally to the direct caller probabilities
                        List<Double> caller_spreading_ratio = new ArrayList<Double>();
                        for (int i = 0; i < callers.size(); i++) {
                            caller_spreading_ratio.set(i, this.ptaskcallers.get(hidx, callers.get(i)));
                        }
                        double caller_ratio_sum = 0;
                        for (int i = 0; i < callers.size(); i++) {
                            caller_ratio_sum += caller_spreading_ratio.get(i);
                        }
                        for (int i = 0; i < callers.size(); i++) {
                            caller_spreading_ratio.set(i, caller_spreading_ratio.get(i) / caller_ratio_sum);
                        }
                        for (int k = 0; k < callers.size(); k++) {
                            int c = callers.get(k);
                            double num1 = ilscaling.get(c, hidx);
                            double num2 = (double) multremote / (double) multcallers * caller_spreading_ratio.get(k);
                            double num = Math.min(num1, num2);
                            ilscaling.set(c, hidx, num);
                        }
                    }
                }
            }


            for (int t = 1; t <= lqn.ntasks; t++) {
                int tidx = lqn.tshift + t;
                if (lqn.isref.get(tidx) == 0) {
                    // the following are remote (indirect) callers that certain to be
                    // callers of task t, hence if they have multiplicity m then task t
                    // cannot have as a matter of fact multiplicity more than m
                    boolean isolated_task = true;
                    List<Integer> callers = new ArrayList<>();
                    for (int i = 1; i <= lqn.nidx; i++) {
                        if (lqn.iscaller.isAssigned(tidx, i) || lqn.iscaller.isAssigned(i, tidx)) {
                            isolated_task = false;
                            break;
                        }
                    }
                    if ((int) lqn.isref.get(tidx) == 0 && !isolated_task) {

                        for (int eidx : lqn.entriesof.get(tidx)) {
                            for (int row = lqn.tshift + 1; row < lqn.tshift + lqn.ntasks + 1; row++) {
                                if (lqn.iscaller.isAssigned(row, eidx)) {
                                    callers.add(row);
                                }
                            }
                        }
                    }
                    int multcallers = 0;
                    for (int i : callers) {
                        multcallers += this.njobsorig.get(i, tidx);
                    }
                    Matrix rowhidx = new Matrix(1, this.ptaskcallers_step.get(step).numCols,
                            this.ptaskcallers_step.get(step).numCols);
                    Matrix.extractRows(this.ptaskcallers_step.get(step), tidx, tidx, rowhidx);
                    Matrix step_callers = rowhidx.find();
                    int multremote = 0;
                    for (int i = 1; i < step_callers.length(); i++) {
                        int remidx = (int) step_callers.get(i);
                        if (lqn.schedid.get(remidx) == SchedStrategy.toID(SchedStrategy.INF)) {
                            multremote += this.ptaskcallers_step.get(step).get(tidx, remidx) * this.util.get(remidx);
                        } else {
                            multremote += this.ptaskcallers_step.get(step).get(tidx, remidx) * this.util.get(remidx) *
                                    lqn.mult.get(remidx);
                        }
                    }
                    if ((multcallers > multremote && multremote > 0) && !Double.isInfinite(multremote)) {
                        // we spread the scaling proportionally to the direct caller probabilities
                        List<Double> caller_spreading_ratio = new ArrayList<Double>();
                        for (int i = 0; i < callers.size(); i++) {
                            caller_spreading_ratio.set(i, this.ptaskcallers.get(tidx, callers.get(i)));
                        }
                        double caller_ratio_sum = 0;
                        for (int i = 0; i < callers.size(); i++) {
                            caller_ratio_sum += caller_spreading_ratio.get(i);
                        }
                        for (int i = 0; i < callers.size(); i++) {
                            caller_spreading_ratio.set(i, caller_spreading_ratio.get(i) / caller_ratio_sum);
                        }
                        for (int k = 0; k < callers.size(); k++) {
                            int c = callers.get(k);
                            double num1 = ilscaling.get(c, tidx);
                            double num2 = (double) multremote / (double) multcallers * caller_spreading_ratio.get(k);
                            double num = Math.min(num1, num2);
                            ilscaling.set(c, tidx, num);
                        }
                    }
                }
            }
        }
        // this.ilscaling starting from (0,0)
        // this.njobs starting from (1,1)
        this.ilscaling = ilscaling.clone();
        for (int i = 0; i < this.ilscaling.numRows; i++) {
            for (int j = 0; j < this.ilscaling.numCols; j++) {
                this.njobs.set(i, j, this.njobsorig.get(i, j) * this.ilscaling.get(i, j));
            }
        }
    }


    public void updateThinkTimes(int it) {
        //task16:updateThinkTimes function to be written
        if (this.lqn.iscaller.numCols > 0) { // ignore models without callers
            Matrix torder = new Matrix(1, lqn.ntasks + 1, lqn.ntasks);
            for (int i = 1; i < lqn.ntasks + 1; i++)
                torder.set(i, i);

            this.thinkt = new Matrix(1, this.lqn.ntasks + this.lqn.tshift,
                    this.lqn.ntasks + this.lqn.ntasks - 1);
            this.thinktproc = new HashMap<Integer, Distribution>();

            // solve all task models
            for (int t = 1; t <= this.lqn.ntasks; t++) {
                int tidx = this.lqn.tshift + t;
                double tidx_thinktime = this.lqn.think.get(tidx).getMean(); // user specified think time
                if (!Double.isNaN(this.idxhash.get(tidx) - 1)) { // this skips all REF tasks
                    // obtain total self.tput of task t
                    // mean throughput of task t in the model where it is a server, summed across replicas
                    double njobs = Matrix.extractRows(this.njobsorig, tidx, tidx + 1, null).elementMax();

                    Matrix matrixExtracted = this.results.get(this.results.size()).get(this.idxhash.get(tidx).intValue() - 1).TN;
                    Matrix serverIdxRow = new Matrix(1, matrixExtracted.numCols, matrixExtracted.numCols);
                    int extractRowIndex = (int) this.results.get(this.results.size()).get(this.idxhash.get(tidx).intValue() - 1).TN.
                            get(this.ensemble[this.idxhash.get(tidx).intValue() - 1].getAttribute().getServerIdx());
                    Matrix.extractRows(matrixExtracted, extractRowIndex, extractRowIndex + 1, serverIdxRow);
                    this.tput.set(tidx, this.lqn.repl.get(tidx) * serverIdxRow.elementSum());

                    // obtain total self.utilization of task t
                    Matrix UmatrixExtracted = this.results.get(this.results.size()).get(this.idxhash.get(tidx).intValue() - 1).UN;
                    Matrix UserverIdxRow = new Matrix(1, matrixExtracted.numCols, matrixExtracted.numCols);
                    int UextractRowIndex = (int) this.results.get(this.results.size()).get(this.idxhash.get(tidx).intValue() - 1).UN.
                            get(this.ensemble[this.idxhash.get(tidx).intValue() - 1].getAttribute().getServerIdx());
                    Matrix.extractRows(UmatrixExtracted, UextractRowIndex, UextractRowIndex + 1, UserverIdxRow);
                    this.util.set(tidx, UserverIdxRow.elementSum());

                    if (this.lqn.schedid.get(tidx) == SchedStrategy.toID(SchedStrategy.INF)) { // first we consider the update where t is an infinite server
                        // key think time update formula for LQNs, this accounts for the fact that in LINE infinite server self.utilization is dimensionally a mean number of jobs
                        this.thinkt.set(tidx - 1, (njobs - this.util.get(tidx)) / this.tput.get(tidx) - tidx_thinktime);
                    } else { // otherwise we consider the case where t is a regular queueing station (other than an infinite server)
                        // key think time update formula for LQNs, this accounts that in LINE self.utilization is scaled in [0,1] for all queueing stations irrespectively of the number of servers
                        this.thinkt.set(tidx - 1, njobs * Math.abs(1 - this.util.get(tidx)) / this.tput.get(tidx) - tidx_thinktime);
                    }
                    Exp exponential = new Exp(this.thinkt.get(tidx - 1) + tidx_thinktime);
                    this.thinktproc.put(tidx, exponential);
                } else { // set to zero if this is a ref task
                    this.thinkt.set(tidx - 1, GlobalConstants.Zero);
                    this.thinktproc.put(tidx, new Immediate());
                }
            }
        }
    }

    public void updateMetrics(int it) {
        switch(options.method) {
            case "default":
                this.updateMetricsDefault(it);
            case "moment3":
                // TODO
        }
    }

    public void updateMetricsDefault(int it) {
        //task17:updateMetrics function to be written
        LayeredNetworkStruct lqn = this.lqn;

        // obtain the activity service times
        this.servt = new Matrix(1, lqn.nidx, lqn.nidx);
        for (int r = 1; r < this.servt_classes_updmap.numRows; r++) {
            int idx = (int) this.servt_classes_updmap.get(r, 1);     //layer
            int aidx = (int) this.servt_classes_updmap.get(r, 2);    //activity
            int nodeidx = (int) this.servt_classes_updmap.get(r, 3); //node
            int classidx = (int) this.servt_classes_updmap.get(r, 4); //jobclass

            // store the residence times and tput at this layer to become
            // the servt / tputs of aidx in another layer, as needed
            // this.servt starts from 0 for JLineMatrix Multiplication
            this.servt.set(aidx - 1, this.results.get(results.size()).get(this.idxhash.get(idx).intValue() - 1).RN.
                    get(nodeidx - 1, classidx - 1));
            this.tput.set(aidx, this.results.get(results.size()).get(this.idxhash.get(idx).intValue() - 1).TN.
                    get(nodeidx - 1, classidx - 1));
            Exp mean = new Exp(1 / this.servt.get(aidx - 1));
            Exp rate = new Exp(this.servt.get(aidx - 1));
            this.servtproc.put(aidx, mean);
            //this.tputproc.put(aidx, rate);
        }

        // obtain the call residence time
        this.callresidt = new Matrix(1, lqn.ncalls, lqn.ncalls);
        for (int r = 1; r < this.call_classes_updmap.numRows; r++) {
            int idx = (int) this.call_classes_updmap.get(r, 1);     // layer
            int cidx = (int) this.call_classes_updmap.get(r, 2);    // call
            int nodeidx = (int) this.call_classes_updmap.get(r, 3);// node
            int classidx = (int) this.call_classes_updmap.get(r, 4);// jobclass

            if (this.call_classes_updmap.get(r, 3) > 1) {
                if (nodeidx == 1) {
                    this.callresidt.set(cidx - 1, 0);
                } else {
                    this.callresidt.set(cidx - 1, this.results.get(results.size()).get(this.idxhash.get(idx).intValue() - 1).RN.get(nodeidx - 1, classidx - 1));
                }
            }
        }

        //then resolve the entry servt summming up these contributions
        Matrix out = new Matrix(1, lqn.nidx + lqn.ncalls + 1, lqn.nidx + lqn.ncalls);
        Matrix entry_servt = new Matrix(this.servtmatrix.numRows, 1, this.servtmatrix.numRows);
        Matrix.concatColumns(this.servt, this.callresidt, out);
        this.servtmatrix.mult(out.transpose(), entry_servt);

        for (int i = 0; i < lqn.eshift; i++) {
            entry_servt.set(i, 0, 0);
        }


        // this block fixes the problem that ResidT is scaled so that the task as Vtask = 1,
        // but in call servt the entries need to have Ventry = 1
        for (int eidx = lqn.eshift + 1; eidx < lqn.eshift + lqn.nentries + 1; eidx++) {
            int tidx = (int) lqn.parent.get(eidx); //  task of entry
            int hidx = (int) lqn.parent.get(tidx); // host of entry
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
                task_tput += this.results.get(results.size()).get(this.idxhash.get(hidx).intValue() - 1).TN.
                        get(ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getClientIdx() - 1, tidxclass.get(i) - 1);
            }


            for (int i = 0; i < eidxclass.size(); i++) {
                entry_tput += this.results.get(results.size()).get(this.idxhash.get(hidx).intValue() - 1).TN.
                        get(ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getClientIdx() - 1, eidxclass.get(i) - 1);
            }

            this.servt.set(eidx - 1, entry_servt.get(eidx) * task_tput / entry_tput);
        }

        for (int i = 1; i < this.call_classes_updmap.numRows; i++) {
            int cidx = (int) this.call_classes_updmap.get(i, 2);
            int eidx = (int) lqn.callpair.get(cidx, 2);
            if (this.call_classes_updmap.get(i, 3) > 1) {
                this.servtproc.put(eidx, new Exp(1 / this.servt.get(eidx - 1)));
            }
        }

        // determine call response time processes
        // this.callresidt starts from 0
        for (int i = 1; i < this.call_classes_updmap.numRows; i++) {
            int cidx = (int) this.call_classes_updmap.get(i, 2);
            int eidx = (int) lqn.callpair.get(cidx, 2);
            if (this.call_classes_updmap.get(i, 3) > 1) {
                if (it == 1) {
                    // note that respt is per visit, so number of calls is 1
                    this.callresidt.set(cidx - 1, this.servt.get(eidx - 1));
                    this.callresidtproc.put(cidx, this.servtproc.get(eidx));
                } else {
                    // note that respt is per visit, so number of calls is 1
                    this.callresidtproc.put(cidx, new Exp(1 / this.callresidt.get(cidx - 1)));
                }
            }
        }

        this.ptaskcallers = new Matrix(this.ptaskcallers.numRows, this.ptaskcallers.numCols,
                this.ptaskcallers.numRows * this.ptaskcallers.numCols - 1);

        for (int i = 0; i < this.ptaskcallers.numRows; i++) {
            for (int j = 0; j < this.ptaskcallers.numCols; j++) {
                this.ptaskcallers.set(i, j, 0);
            }
        }

        // determine ptaskcallers for direct callers to tasks
        for (int t = 1; t <= lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            boolean isolated_task = true;
            for (int i = 1; i <= lqn.nidx; i++) {
                if (lqn.iscaller.isAssigned(tidx, i) || lqn.iscaller.isAssigned(i, tidx)) {
                    isolated_task = false;
                    break;
                }
            }
            if ((int) lqn.isref.get(tidx) == 0 && !isolated_task) {
                List<Integer> callers = new ArrayList<>();
                for (int eidx : lqn.entriesof.get(tidx)) {
                    for (int row = lqn.tshift + 1; row < lqn.tshift + lqn.ntasks + 1; row++) {
                        if (lqn.iscaller.isAssigned(row, eidx)) {
                            callers.add(row);
                        }
                    }
                }
                Matrix caller_tput = new Matrix(lqn.ntasks + 1, 1, lqn.ntasks);
                List<Double> caller_idxclass = new ArrayList<Double>();
                for (int caller_idx : callers) {
                    List<Integer> keys = new ArrayList<Integer>();
                    Map<Integer, Integer[]> taskmap = this.ensemble[idxhash.get(tidx).intValue() - 1].getAttribute().getTasks();
                    for (int i = 1; i < taskmap.size(); i++) {
                        if (taskmap.get(i)[1] == caller_idx)
                            keys.add(i);
                    }
                    for (int i : keys) {
                        caller_idxclass.add((double) taskmap.get(1 + i)[0]);
                    }

                    double sum = 0;
                    for (double j : caller_idxclass) {
                        Matrix Tn = results.get(results.size()).get(idxhash.get(tidx).intValue() - 1).TN;
                        sum += Tn.get(this.ensemble[idxhash.get(tidx).intValue() - 1].getAttribute().getClientIdx() - 1
                                , (int) j - 1);
                    }
                    caller_tput.set(caller_idx - lqn.tshift - 1, sum);
                }
                double task_tput = 0;
                for (int i = 1; i < caller_tput.numRows; i++) {
                    task_tput += caller_tput.get(i - 1);
                }
                for (int i = 1; i <= lqn.ntasks; i++) {
                    this.ptaskcallers.set(tidx, lqn.tshift + i, caller_tput.get(i - 1) / task_tput);
                }
            }
        }


        // determine ptaskcallers for direct callers to hosts
        for (int hidx = 1; hidx <= lqn.nhosts; hidx++) {
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
                    sum += this.results.get(this.results.size()).get(this.idxhash.get(hidx).intValue() - 1).TN.
                            get(this.ensemble[this.idxhash.get(hidx).intValue() - 1].getAttribute().getClientIdx() - 1,
                                    i - 1);
                }
                caller_tput.set(caller_idx - lqn.tshift, caller_tput.get(caller_idx - lqn.tshift) + sum);
            }
            double host_tput = caller_tput.elementSum();
            for (int i = 1; i <= lqn.ntasks; i++) {
                this.ptaskcallers.set(hidx, lqn.tshift + i, caller_tput.get(i) / host_tput);
            }
        }


        // impute call probability using a DTMC random walk on the taskcaller graph
        // for matrix multiplication, let P to be the same size of that in Matlab
        Matrix P = new Matrix(this.ptaskcallers.numRows - 1, this.ptaskcallers.numCols - 1,
                (this.ptaskcallers.numRows - 1) * (this.ptaskcallers.numCols - 1));
        for (int i = 1; i < this.ptaskcallers.numRows; i++) {
            for (int j = 1; j < this.ptaskcallers.numCols; j++) {
                P.set(i - 1, j - 1, this.ptaskcallers.get(i, j));
            }
        }

        for (int i = 0; i < P.numRows; i++) {
            if (P.sumRows(i) > 0) {
                for (int j = 0; j < P.numCols; j++) {
                    P.set(i, j, P.get(i, j) / P.sumRows(i));
                }
                P.set(i, i, 1 - P.sumRows(i) - P.get(i, i));
            } else {
                for (int j = 0; j < P.numCols; j++) {
                    P.set(i, j, 0);
                }
                P.set(i, i, 1);
            }
        } //P = dtmc_makestochastic(P); % hold mass at reference stations when there


        this.ptaskcallers_step.put(1, P);

        for (int h = 1; h <= lqn.nhosts; h++) {
            int hidx = h;
            for (int i = 0; i < lqn.tasksof.get(hidx).size(); i++) {
                int tidx = lqn.tasksof.get(hidx).get(i);
                // initialize the probability mass on tidx
                Matrix x0 = new Matrix(1, P.length(), P.length());
                x0.set(hidx - 1, 1);
                int step = 1;
                // start the walk backward to impute probability of indirect callers
                Matrix x = new Matrix(x0.numRows, P.numCols, x0.numRows * P.numCols);
                x0.mult(P, x);

                for (int e = 1; e <= this.nlayers; e++) {
                    step += 1;
                    Matrix y = new Matrix(x.numRows, P.length(), P.length());
                    x.mult(P, y);
                    x = y.clone();
                    double sum = 0;
                    Matrix ref_nonzero = lqn.isref.find();
                    for (int k = 1; k <= ref_nonzero.numCols; k++) {
                        sum += x.get(k);
                    }

                    if (sum > 1.0 - GlobalConstants.CoarseTol)
                        break;
                    for (int index = 1; index < x.numRows; index++) {
                        this.ptaskcallers_step.get(step).set(index, tidx, x.get(index));
                    }
                    for (int index = 1; index < x.numRows; index++) {
                        Matrix out1 = new Matrix(ptaskcallers.numRows, 1, ptaskcallers.numRows);
                        out1 = Matrix.extractColumn(this.ptaskcallers, tidx, out1);
                        double x1 = out1.elementMax();
                        double x2 = x.elementMax();
                        double max = Math.max(x1, x2);
                        max = max >= 0 ? max : 0;
                        this.ptaskcallers.set(index, tidx, max);
                    }
                }
            }
        }
    }

    public void updateRoutingProbabilities(int it) {
        //task18:updateMetrics function to be written
        int map_length = 0;
        if (unique_route_prob_updmap.numRows != 0 && unique_route_prob_updmap.numCols != 0) {
            map_length = unique_route_prob_updmap.length();
        }

        for (int u = 1; u <= map_length; u++) {
            int idx;
            if (it != 0) {
                idx = (int) this.unique_route_prob_updmap.get(u - 1);
            } else {
                idx = (int) this.unique_route_prob_updmap.get(this.unique_route_prob_updmap.length() - u + 1);
            }
            boolean idx_updated = false;

            Network nt = this.ensemble[this.idxhash.get(idx).intValue() - 1];
            RoutingMatrix P = new RoutingMatrix(nt, nt.getJobClass(), nt.getNodes());

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
                //TODO: Cache
//                if (this.ensemble[this.idxhash.get(idx).intValue() - 1].items) { // if idx is a cache node
//                    JLineMatrix TN_copy = this.results.get(this.results.size()).get(this.idxhash.get(((int) host)).intValue() - 1).TN;
//                    double Xtot = TN_copy.sumCols(this.ensemble[this.idxhash.get((int) host).intValue() - 1].getAttribute().getServerIdx() - 1);
//                    if (Xtot > 0) {
//                        JLineMatrix hm_tput_copy = this.results.get(this.results.size()).get(this.idxhash.
//                                get((int) host).intValue() - 1).TN;
//                        double hm_tput = hm_tput_copy.get(this.ensemble[this.idxhash.get((int) host).intValue() - 1].
//                                getAttribute().getServerIdx() -1 , (int) classidxto - 1);
//                        P.addConnection(nt.getNodeByStatefulIndex((int) nodefrom), nt.getNodeByStatefulIndex((int) nodeto),
//                                nt.getJobClassFromIndex((int) classidxfrom - 1), nt.getClassByIndex((int) classidxto - 1),
//                                hm_tput / Xtot);
//                        idx_updated = true;
//                    }
//                } else {
//                    JLineMatrix TN_copy = this.results.get(this.results.size()).get(this.idxhash.get(((int) tidx_caller)).intValue() - 1).TN;
//                    double Xtot = TN_copy.sumCols(this.ensemble[this.idxhash.get((int) tidx_caller).intValue() - 1].getAttribute().getServerIdx() - 1);
//                    if (Xtot > 0) {
//                        Map<Integer, Integer[]> call_map = this.ensemble[this.idxhash.get((int) tidx_caller).intValue() - 1].getAttribute().getCalls();
//                        int k;
//                        for (k = 0; k < call_map.size(); k++) {
//                            if (eidx == call_map.get(k)[4])
//                                break;
//                        }
//                        int eidxclass = call_map.get(k)[1];
//                        JLineMatrix matrixcopy = this.results.get(this.results.size()).get(this.idxhash.
//                                get((int) tidx_caller).intValue() - 1).TN;
//                        double entry_tput = matrixcopy.get(this.ensemble[this.idxhash.get((int) tidx_caller).intValue() - 1].
//                                getAttribute().getServerIdx() -1 , eidxclass - 1);
//                        P.addConnection(nt.getNodeByStatefulIndex((int) nodefrom), nt.getNodeByStatefulIndex((int) nodeto),
//                                nt.getJobClassFromIndex((int) classidxfrom), nt.getClassByIndex((int) classidxto),
//                                entry_tput / Xtot);
//                        idx_updated = true;
//                    }
//                }
            }

            if (idx_updated) {
                this.ensemble[this.idxhash.get(idx).intValue()].link(P);
            }
        }
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
        Matrix aidxrow = new Matrix(1, lqn.graph.numCols, lqn.graph.numCols);
        aidxrow = Matrix.extractRows(lqn.graph, aidx, aidx + 1, aidxrow);
        Matrix nextaidxs = aidxrow.find();
        for (int i = 0; i < nextaidxs.numRows; i++) {

            int nextaidx = (int) nextaidxs.get(i);
            if (lqn.parent.get(aidx) != lqn.parent.get(nextaidx)) {
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
            if (lqn.parent.get(aidx) == lqn.parent.get(nextaidx)) {
                if (nextaidx != aidx) {
                    double Gvalue = lqn.graph.get(aidx, nextaidx) > 0 ? lqn.graph.get(aidx, nextaidx) : 0;
                    U.set(eidx - 1, nextaidx - 1, U.get(eidx - 1, nextaidx - 1) + Gvalue);
                    U = getEntryServiceMatrixRecursion(lqn, nextaidx, eidx, U);
                }
            }
        }
        return U;
    }

    @Override
    protected void runAnalyzer() throws IllegalAccessException {


    }

    private Matrix cellToMatrix(Map<Integer, List<Integer[]>> cell) {
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

    public Matrix getCall_classes_updmap() {
        return call_classes_updmap;
    }

    public Matrix getArvproc_classes_updmap() {
        return arvproc_classes_updmap;
    }

    public Matrix getServt_classes_updmap() {
        return servt_classes_updmap;
    }

    public Matrix getRoute_prob_updmap() {
        return route_prob_updmap;
    }

    public Matrix getThinkt_classes_updmap() {
        return thinkt_classes_updmap;
    }

    public List<Network> getEnsemble() {
        List<Network> myEnsemble = new ArrayList<>();
        for (int e=0; e<ensemble.length; e++){
            myEnsemble.add(ensemble[e]);
        }
        return myEnsemble;
    }

    public List<Double> getIdxhash() {
        return idxhash;
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

}

