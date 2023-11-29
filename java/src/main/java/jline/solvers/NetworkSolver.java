// Copyright (c) 2012-2022, Imperial College London
// All rights reserved.

package jline.solvers;

import jline.lang.JobClass;
import jline.lang.constant.GlobalConstants;
import jline.lang.nodes.Station;
import jline.util.Matrix;
import jline.lang.Network;
import jline.lang.NetworkStruct;

import javax.xml.parsers.ParserConfigurationException;
import java.util.*;

// Abstract class for solvers applicable to Network models
public abstract class NetworkSolver extends Solver {

    public Network model; // Model to be solved
    public NetworkStruct sn; // Structure describing the model
    public SolverHandles handles; // Performance metric handles
    public SolverMetrics metrics; // Performance metrics

    // Construct a NetworkSolver with given model, name and options data structure
    protected NetworkSolver(Network model, String name, SolverOptions options) {
        super(name, options);
        this.model = model;
        if (model.getNumberOfNodes() == 0) {
            throw new RuntimeException("The model supplied in input is empty.");
        }
        this.handles = model.getAvgHandles();
        // TODO: get and set transient handles
        // [Qt,Ut,Tt] = model.getTranHandles;
        // self.setTranHandles(Qt,Ut,Tt);
        // this.setTranHandles();
        this.sn = model.getStruct(true); // Force model to refresh
    }

    protected NetworkSolver(Network model, String name) {
        this(model, name, defaultOptions());
    }

    protected void setTranHandles() {
        // TODO: implementation - note arguments should likely not be void
        throw new RuntimeException("setTranHandles() has not yet been implemented in JLINE.");
    }

    public void setAvgHandles(SolverHandles handles) {
        this.handles = handles;
    }

    protected void getTranHandles() {
        // TODO: implementation - note return type should likely not be void
        throw new RuntimeException("getTranHandles() has not yet been implemented in JLINE.");
    }

    public SolverHandles getAvgHandles() {
        return this.handles;
    }

    protected Map<Station, Map<JobClass, SolverHandles.Metric>> getAvgQLenHandles() {
        return this.handles.Q;
    }

    protected Map<Station, Map<JobClass, SolverHandles.Metric>> getAvgRespTHandles() {
        // TODO: implementation - note return type should likely not be void
        return this.handles.R;
    }

    protected Map<Station, Map<JobClass, SolverHandles.Metric>> getAvgUtilHandles() {
        // TODO: implementation - note return type should likely not be void
        return this.handles.U;
    }

    protected Map<Station, Map<JobClass, SolverHandles.Metric>> getAvgTputHandles() {
        // TODO: implementation - note return type should likely not be void
        return this.handles.T;
    }

    // Returns true if the solver has computed steady-state average metrics.
    protected boolean hasAvgResults() {
        return !((result.QN == null || result.QN.isEmpty()) &&
                (result.UN == null || result.UN.isEmpty()) &&
                (result.RN == null || result.RN.isEmpty()) &&
                (result.TN == null || result.TN.isEmpty()) &&
                (result.CN == null || result.CN.isEmpty()) &&
                (result.XN == null || result.XN.isEmpty()));
    }

    // Returns true if the solver has computed transient average metrics.
    protected boolean hasTranResults() {
        if (this.hasResults()) {
            return result.QNt.length > 0 && result.QNt[0].length > 0 && !result.QNt[0][0].isEmpty();
        }
        return false;
    }

    // Returns true if the solver has computed steady-state distribution metrics.
    protected boolean hasDistribResults() {
        // TODO: implementation
        throw new RuntimeException("hasDistribResults() has not yet been implemented in JLINE.");
    }

    // Get agent representation
    protected final void getAG() {
        // TODO: implementation - note return type should likely not be void
        throw new RuntimeException("getAG() has not yet been implemented in JLINE.");
    }

    // Compute average queue-lengths at steady-state
    protected final void getAvgQLen() {
        // TODO: implementation - note return type should likely not be void
        throw new RuntimeException("getAvgQLen() has not yet been implemented in JLINE.");
    }

    // Compute average utilizations at steady-state
    protected final void getAvgUtil() {
        // TODO: implementation - note return type should likely not be void
        throw new RuntimeException("getAvgUtil() has not yet been implemented in JLINE.");
    }

    // Compute average response times at steady-state
    protected final void getAvgRespT() {
        // TODO: implementation - note return type should likely not be void
        throw new RuntimeException("getAvgRespT() has not yet been implemented in JLINE.");
    }

    // Compute average waiting time in queue excluding service
    protected final void getAvgWaitT() {
        // TODO: implementation - note return type should likely not be void
        throw new RuntimeException("getAvgWaitT() has not yet been implemented in JLINE.");
    }

    // Compute average throughputs at steady-state
    protected final void getAvgTput() {
        // TODO: implementation - note return type should likely not be void
        throw new RuntimeException("getAvgTput() has not yet been implemented in JLINE.");
    }

    // Compute average arrival rate at steady-state
    protected final void getAvgArvR() {
        // TODO: implementation - note return type should likely not be void
        throw new RuntimeException("getAvgArvR() has not yet been implemented in JLINE.");
    }

    // Returns average station metrics at steady-state
    public final void getAvg() {

        // TODO: provide polymorphic version where handles can be passed in as parameters individually?
        // TODO: provide polymorphic version where handles can be passed in as a whole? (Lines 14-20)

        if (this.handles == null || this.handles.Q == null || this.handles.U == null || this.handles.R == null ||
                this.handles.T == null || this.handles.A == null) {
            resetResults();
            this.handles = model.getAvgHandles();
        }

        if (Double.isFinite(options.timespan[1])) {
            throw new RuntimeException(
                    "The getAvg method does not support the timespan option, use the getTranAvg method instead.");
        } else {
            this.options.timespan[0] = 0;
            this.options.timespan[1] = Double.POSITIVE_INFINITY;
        }

        if (!this.hasAvgResults() || !this.options.cache) {
            try {
                runAnalyzer();
            } catch (IllegalAccessException e) {
                System.out.println("IllegalAccessException upon running runAnalyzer()");
            } catch (ParserConfigurationException e) {
                System.err.println("ParserConfigurationException upon running runAnalyzer()");
            }
            // TODO: provide more granular error messaging (if useful) (Lines 33-49)
            if (!this.hasAvgResults()) {
                throw new RuntimeException("Line is unable to return results for this model.");
            }
        } // else return cached value

        int M = sn.nstations;
        int K = sn.nclasses;

        int Vrows = sn.visits.get(0).getNumRows();
        int Vcols = sn.visits.get(0).getNumCols();
        int Vcells = sn.visits.size();
        Matrix V = new Matrix(Vrows, Vcols);
        for (int i = 0; i < Vrows; i++) {
            for (int j = 0; j < Vcols; j++) {
                double tmpSum = 0;
                for (int k = 0; k < Vcells; k++) {
                    tmpSum += sn.visits.get(k).get(i, j);
                }
                V.set(i, j, tmpSum);
            }
        }

        Matrix QNclass = new Matrix(0, 0);
        Matrix UNclass = new Matrix(0, 0);
        Matrix RNclass = new Matrix(0, 0);
        Matrix TNclass = new Matrix(0, 0);
        Matrix ANclass = new Matrix(0, 0);
        Matrix WNclass;

        if (!this.result.QN.isEmpty()) {
            QNclass = new Matrix(M, K);
            for (int k = 0; k < K; k++) {
                for (int i = 0; i < M; i++) {
                    if (!this.handles.Q.get(this.model.getStations().get(i)).get(this.model.getClassByIndex(k)).isDisabled
                            && !this.result.QN.isEmpty()) {
                        QNclass.set(i, k, this.result.QN.get(i, k));
                        if (QNclass.get(i, k) < GlobalConstants.FineTol) { // Round to zero numerical perturbations
                            QNclass.set(i, k, 0);
                        }
                        if (Double.isNaN(QNclass.get(i, k))) { // Indicates that a metric is disabled
                            QNclass.set(i, k, 0);
                        }
                    } else {
                        QNclass.set(i, k, 0); // Indicates that a metric is disabled
                    }
                }
            }
        }

        if (!this.result.UN.isEmpty()) {
            UNclass = new Matrix(M, K);
            for (int k = 0; k < K; k++) {
                for (int i = 0; i < M; i++) {
                    if (!this.handles.U.get(this.model.getStations().get(i)).get(this.model.getClassByIndex(k)).isDisabled
                            && !this.result.UN.isEmpty()) {
                        UNclass.set(i, k, this.result.UN.get(i, k));
                        if (UNclass.get(i, k) < GlobalConstants.FineTol) { // Round to zero numerical perturbations
                            UNclass.set(i, k, 0);
                        }
                        if (Double.isNaN(UNclass.get(i, k))) { // Indicates that a metric is disabled
                            UNclass.set(i, k, 0);
                        }
                    } else {
                        UNclass.set(i, k, 0); // Indicates that a metric is disabled
                    }
                }
            }
        }
        if (!this.result.RN.isEmpty()) {
            RNclass = new Matrix(M, K);
            for (int k = 0; k < K; k++) {
                for (int i = 0; i < M; i++) {
                    if (!this.handles.R.get(this.model.getStations().get(i)).get(this.model.getClassByIndex(k)).isDisabled
                            && !this.result.RN.isEmpty()) {
                        RNclass.set(i, k, this.result.RN.get(i, k));
                        if (RNclass.get(i, k) < GlobalConstants.FineTol) { // Round to zero numerical perturbations
                            RNclass.set(i, k, 0);
                        }
                        if (Double.isNaN(RNclass.get(i, k))) { // Indicates that a metric is disabled
                            RNclass.set(i, k, 0);
                        }
                    } else {
                        RNclass.set(i, k, 0); // Indicates that a metric is disabled
                    }
                }
            }
        }

        if (!this.result.TN.isEmpty()) {
            TNclass = new Matrix(M, K);
            for (int k = 0; k < K; k++) {
                for (int i = 0; i < M; i++) {
                    if (!this.handles.T.get(this.model.getStations().get(i)).get(this.model.getClassByIndex(k)).isDisabled
                            && !this.result.TN.isEmpty()) {
                        TNclass.set(i, k, this.result.TN.get(i, k));
                        if (TNclass.get(i, k) < GlobalConstants.FineTol) { // Round to zero numerical perturbations
                            TNclass.set(i, k, 0);
                        }
                        if (Double.isNaN(TNclass.get(i, k))) { // Indicates that a metric is disabled
                            TNclass.set(i, k, 0);
                        }
                    } else {
                        TNclass.set(i, k, 0); // Indicates that a metric is disabled
                    }
                }
            }
        }

        if (!this.result.AN.isEmpty()) {
            ANclass = new Matrix(M, K);
            for (int k = 0; k < K; k++) {
                for (int i = 0; i < M; i++) {
                    if (!this.handles.A.get(this.model.getStations().get(i)).get(this.model.getClassByIndex(k)).isDisabled
                            && !this.result.AN.isEmpty()) {
                        ANclass.set(i, k, this.result.TN.get(i, k));
                        if (ANclass.get(i, k) < GlobalConstants.FineTol) { // Round to zero numerical perturbations
                            ANclass.set(i, k, 0);
                        }
                        if (Double.isNaN(ANclass.get(i, k))) { // Indicates that a metric is disabled
                            ANclass.set(i, k, 0);
                        }
                    } else {
                        ANclass.set(i, k, 0); // Indicates that a metric is disabled
                    }
                }
            }
        }

        // Set to zero entries that are associated to immediate transitions
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < K; j++) {
                if (RNclass.get(i, j) < 10 * GlobalConstants.FineTol) {
                    QNclass.set(i, j, 0);
                    UNclass.set(i, j, 0);
                    RNclass.set(i, j, 0);
                }
            }
        }

        WNclass = RNclass.clone();
        WNclass.zero();
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                if (!RNclass.isEmpty() && RNclass.get(i, k) > 0) {
                    int c = -1;
                    for (int chain = 0; chain < sn.chains.getNumRows(); chain++) {
                        if (sn.chains.get(chain, k) > 0) {
                            c = chain;
                            break;
                        }
                    }
                    if (RNclass.get(i, k) < GlobalConstants.FineTol) {
                        WNclass.set(i, k, RNclass.get(i, k));
                    } else {
                        int refClass = (int) sn.refclass.get(0, c);
                        if (refClass > 0) {
                            // If there is a reference class, use this:
                            WNclass.set(
                                    i,
                                    k,
                                    RNclass.get(i, k) * V.get(i, k) / V.get((int) sn.refstat.get(k, 0), refClass));
                        } else {
                            int Vrow = (int) sn.refstat.get(k, 0);
                            double Vsum = 0;
                            for (int col = 0; col < sn.inchain.get(c).getNumCols(); col++) {
                                Vsum += V.get(Vrow, (int) sn.inchain.get(c).get(0, col));
                            }
                            WNclass.set(i, k, RNclass.get(i, k) * V.get(i, k) / Vsum);
                        }
                    }
                }
            }
        }

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < K; j++) {
                if (Double.isNaN(WNclass.get(i, j))
                        || WNclass.get(i, j) < 10 * GlobalConstants.FineTol
                        || WNclass.get(i, j) < GlobalConstants.Zero) {
                    WNclass.set(i, j, 0);
                }
            }
        }

        if (!UNclass.isEmpty()) {
            boolean unstableQueueFlag = false;
            for (int i = 0; i < M; i++) {
                if (UNclass.sumRows(i) > 0.99 * sn.nservers.get(i, 0)) {
                    unstableQueueFlag = true;
                    break;
                }
            }

            boolean infJobsFlag = false;
            for (int i = 0; i < sn.njobs.length(); i++) {
                if (Double.isInfinite(sn.njobs.get(0, i))) {
                    infJobsFlag = true;
                    break;
                }
            }

            if (unstableQueueFlag && infJobsFlag) {
                System.err.println(
                        "The model has unstable queues, performance metrics may grow unbounded.");
            }
        }

        this.metrics = new SolverMetrics();
        this.metrics.QNclass = QNclass;
        this.metrics.UNclass = UNclass;
        this.metrics.RNclass = RNclass;
        this.metrics.ANclass = ANclass;
        this.metrics.TNclass = TNclass;
        this.metrics.WNclass = WNclass;
    }

    // Compute average utilizations at steady-state for all nodes
    protected final void getAvgNode() {
        // TODO: implementation - note return type and parameters should likely not be void
        throw new RuntimeException("getAvgNode() has not yet been implemented in JLINE.");
    }

    // TODO: getStageTable missing

    // Return table of average station metrics
    public final NetworkAvgTable getAvgTable() {

        // TODO: provide polymorphic version where handles can be passed in as parameters individually?
        // TODO: provide polymorphic version where handles can be passed in as a whole?
        // TODO: provide polymorphic version where keepDisabled can be passed in as a parameter?

        this.sn = model.getStruct(true);

        boolean keepDisabled = false;
        // TODO: [Q,U,R,T,~] = getAvgHandles(self);

        int M = sn.nstations;
        int K = sn.nclasses;

        try {
            if (Double.isFinite(options.timespan[1])) {
                // TODO: [Qt,Ut,Tt] = getTranHandles(self);
                getTranAvg();
            } else {
                getAvg();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        Matrix QN = this.metrics.QNclass;
        Matrix UN = this.metrics.UNclass;
        Matrix RN = this.metrics.RNclass;
        Matrix TN = this.metrics.TNclass;
        Matrix AN = this.metrics.ANclass;

        if (QN.isEmpty()) {
            throw new RuntimeException(
                    "Unable to compute results and therefore unable to print AvgTable.");
        }

        if (!keepDisabled) {
            Matrix V = new Matrix(sn.nstateful, K);
            for (int i = 0; i < sn.visits.size(); i++) {
                V = V.add(1, sn.visits.get(i));
            }
            if (V.isEmpty()) { // SSA
                // TODO: implementation for SSA
                System.out.println("Warning: unimplemented code reached in NetworkSolver.getAvgTable 1.");
            }

            List<Double> Qval = new ArrayList<>();
            List<Double> Uval = new ArrayList<>();
            List<Double> Rval = new ArrayList<>();
            List<Double> Tval = new ArrayList<>();
            List<Double> ArvR = new ArrayList<>();
            List<Double> Residval = new ArrayList<>();
            List<String> className = new ArrayList<>();
            List<String> stationName = new ArrayList<>();
            for (int i = 0; i < M; i++) {
                for (int k = 0; k < K; k++) {
                    // if (QN.get(i, k) + UN.get(i, k) + RN.get(i, k) + TN.get(i, k) > 0) {
                    int c = -1;
                    for (int row = 0; row < sn.chains.numRows; row++) {
                        if (sn.chains.get(row, k) > 0) {
                            c = row;
                            break;
                        }
                    }
                    Qval.add(QN.get(i, k));
                    Uval.add(UN.get(i, k));
                    Rval.add(RN.get(i, k));
                    ArvR.add(AN.get(i, k));
                    Tval.add(TN.get(i, k));
                    className.add(sn.jobclasses.get(k).getName());
                    stationName.add(sn.stations.get(i).getName());
                    if (RN.get(i, k) < GlobalConstants.Zero) {
                        Residval.add(RN.get(i, k));
                    } else {
                        if (sn.refclass.get(c) > 0) {
                            Residval.add(
                                    RN.get(i, k)
                                            * V.get(i, k)
                                            / V.get((int) sn.refstat.get(k), (int) sn.refclass.get(c)));
                        } else {
                            double sum = 0;
                            int row = (int) sn.refstat.get(k);
                            for (int col = 0; col < sn.chains.numCols; col++) {
                                if (sn.chains.get(c, col) > 0) sum += V.get(row, col);
                            }
                            Residval.add(RN.get(i, k) * V.get(i, k) / sum);
                        }
                    }
                    //  }
                }
            }
            NetworkAvgTable avgTable = new NetworkAvgTable(Qval, Uval, Rval, Residval, ArvR, Tval);
            avgTable.setOptions(this.options);
            avgTable.setClassNames(className);
            avgTable.setStationNames(stationName);
            return avgTable;
        } else {
            // TODO: implementation if keepDisabled is set to true
            System.out.println("Warning: unimplemented code reached in NetworkSolver.getAvgTable 2.");
        }
        return null;
    }

    // Return table of average station metrics
    protected final void getAvgQLenTable() {
        // TODO: implementation - note return type and parameters should likely not be void
        throw new RuntimeException("getAvgQLenTable() has not yet been implemented in JLINE.");
    }

    // Return table of average station metrics
    protected final void getAvgUtilTable() {
        // TODO: implementation - note return type and parameters should likely not be void
        throw new RuntimeException("getAvgUtilTable() has not yet been implemented in JLINE.");
    }

    // Return table of average station metrics
    protected final void getAvgRespTTable() {
        // TODO: implementation - note return type and parameters should likely not be void
        throw new RuntimeException("getAvgRespTTable() has not yet been implemented in JLINE.");
    }

    // Return table of average station metrics
    protected final void getAvgTputTable() {
        // TODO: implementation - note return type and parameters should likely not be void
        throw new RuntimeException("getAvgTputTable() has not yet been implemented in JLINE.");
    }

    // Return table of average node metrics
    protected final void getAvgNodeTable() {
        // TODO: implementation - note return type and parameters should likely not be void
        throw new RuntimeException("getAvgNodeTable() has not yet been implemented in JLINE.");
    }

    protected final void getAvgChainTable() {
        // TODO: implementation - note return type and parameters should likely not be void
        throw new RuntimeException("getAvgChainTable() has not yet been implemented in JLINE.");
    }

    // Return average station metrics aggregated by chain
    protected final void getAvgChain() {
        // TODO: implementation - note return type and parameters should likely not be void
        throw new RuntimeException("getAvgChain() has not yet been implemented in JLINE.");
    }

    // Return average arrival rates aggregated by chain
    protected final void getAvgArvRChain() {
        // TODO: implementation - note return type and parameters should likely not be void
        throw new RuntimeException("getAvgArvRChain() has not yet been implemented in JLINE.");
    }

    // Return average queue-lengths aggregated by chain
    protected final void getAvgQLenChain() {
        // TODO: implementation - note return type and parameters should likely not be void
        throw new RuntimeException("getAvgQLenChain() has not yet been implemented in JLINE.");
    }

    // Return average utilization aggregated by chain
    protected final void getAvgUtilChain() {
        // TODO: implementation - note return type and parameters should likely not be void
        throw new RuntimeException("getAvgUtilChain() has not yet been implemented in JLINE.");
    }

    // Return average response time aggregated by chain
    protected final void getAvgRespTChain() {
        // TODO: implementation - note return type and parameters should likely not be void
        throw new RuntimeException("getAvgRespTChain() has not yet been implemented in JLINE.");
    }

    // Return average throughputs aggregated by chain
    protected final void getAvgTputChain() {
        // TODO: implementation - note return type and parameters should likely not be void
        throw new RuntimeException("getAvgTputChain() has not yet been implemented in JLINE.");
    }

    // Return average system metrics at steady state
    public final void getAvgSys() {

        // TODO: provide polymorphic version where R and T handles can be passed in as parameters
        // individually?
        // TODO: provide polymorphic version where handles can be passed in as a whole? (Lines 14-20)

        this.sn = model.getStruct(true);

        this.getAvg();

        if (this.model.hasFork()) {
            if (this.model.hasOpenClasses()) {
                System.err.println(
                        "System response time computation not yet supported with open classes in the presence of fork nodes.");
                this.result.RN.fill(Double.NaN);
            }
        }
        if (this.model.hasJoin()) {
            if (this.model.hasOpenClasses()) {
                System.err.println(
                        "System response time computation not yet supported with open classes in the presence of join nodes.");
                this.result.RN.fill(Double.NaN);
            }
        }

        boolean[] completes = new boolean[sn.nclasses];
        for (int idx = 0; idx < sn.nclasses; idx++) {
            completes[idx] = model.getClasses().get(idx).getCompletes();
        }

        // TODO: if any(isinf(sn.njobs')) // If the model has any open class
        // This could be optimised by computing the statistics only for open chains

        // Compute chain visits
        Matrix alpha = new Matrix(sn.nstations, sn.nclasses);
        Matrix CNclass = new Matrix(1, sn.nclasses);
        if (!this.model.hasJoin() && !this.model.hasFork()) {
            for (int c = 0; c < sn.nchains; c++) {
                Matrix inchain = sn.inchain.get(c);
                for (int i = 0; i < inchain.length(); i++) {
                    int r = (int) inchain.get(i);
                    for (int j = 0; j < sn.nstations; j++) {
                        // Not empty and not a source
                        if (!this.result.RN.isEmpty() &&
                                (!(Double.isInfinite(this.sn.njobs.get(r)) && j == sn.refstat.get(r)))) {
                            CNclass.set(
                                    0,
                                    r,
                                    CNclass.get(0, r)
                                            + sn.visits.get(c).get((int) sn.stationToStateful.get(j), r)
                                            * this.result.RN.get(j, r)
                                            / sn.visits
                                            .get(c)
                                            .get((int) sn.stationToStateful.get((int) sn.refstat.get(r)), r));
                        }
                    }
                }
            }
        }


        for (int c = 0; c < sn.nchains; c++) {
            Matrix inchain = sn.inchain.get(c);
            Matrix completingClasses = Matrix.extractRows(sn.chains, c, c + 1, null);
            for (int i = 0; i < completingClasses.length(); i++) {
                if (!completes[i]) {
                    completingClasses.set(0, i, Double.NaN);
                }
            }

            for (int i = 0; i < sn.nstations; i++) {
                if (sn.refclass.get(c) >= 0) {
                    // For all classes within the chain (a class belongs to a single chain, the reference
                    // station must be identical for all classes within a chain)
                    List<Double> intersection = Matrix.intersect(sn.refclass.findNonNegative(), inchain);
                    for (double value : intersection) {
                        int k = (int) value;
                        double sumVisits = 0.0;
                        for (int idx = 0; idx < completingClasses.length(); idx++) {
                            if (completingClasses.get(idx) == 1) {
                                sumVisits +=
                                        sn.visits
                                                .get(c)
                                                .get(
                                                        (int) sn.stationToStateful.get((int) sn.refstat.get(k)),
                                                        idx);
                            }
                        }
                        alpha.set(
                                i,
                                k,
                                alpha.get(i, k)
                                        + sn.visits.get(c).get((int) sn.stationToStateful.get(i), k) / sumVisits);
                    }
                } else {
                    // For all classes within the chain (a class belongs to a single chain, the reference
                    // station must be identical for all classes within a chain)
                    for (int j = 0; j < inchain.length(); j++) {
                        int k = (int) inchain.get(j);
                        double sumVisits = 0.0;
                        for (int idx = 0; idx < completingClasses.length(); idx++) {
                            if (completingClasses.get(idx) == 1) {
                                sumVisits += sn.visits.get(c).get((int) sn.stationToStateful.get((int) sn.refstat.get(k)), idx);
                            }
                        }
                        alpha.set(i, k, alpha.get(i, k) + sn.visits.get(c).get(i, k) / sumVisits);
                    }
                }
            }
        }
        for (int i = 0; i < sn.nstations; i++) {
            for (int k = 0; k < sn.nclasses; k++) {
                if (!Double.isFinite(alpha.get(i, k))) {
                    alpha.set(i, k, 0.0);
                }
            }
        }

        // Compute average chain metrics
        this.metrics.CNchain = new Matrix(1, sn.nchains);
        this.metrics.XNchain = new Matrix(1, sn.nchains);

        for (int c = 0; c < sn.nchains; c++) {
            Matrix inchain = sn.inchain.get(c);
            Matrix completingClasses = Matrix.extractRows(sn.chains, c, c + 1, null).find();
            completingClasses = completingClasses.transpose();
            for (int i = 0; i < inchain.length(); i++) {
                if (!this.model.getClasses().get(i).getCompletes()) {
                    completingClasses.set(0, i, Double.NaN);
                }
            }

            if (!result.TN.isEmpty()) {
                // All classes in same chain must share the same refstation, so we use the first one
                int ref = (int) sn.refstat.get((int) inchain.get(0));
                // We now compute the incoming system throughput to the reference station from completing
                // classes
                for (int i = 0; i < sn.nstations; i++) {
                    for (int j = 0; j < completingClasses.length(); j++) {
                        int r = (int) completingClasses.get(j);
                        if (completingClasses.get(j) >=0) {
                            List<Double> intersection = Matrix.intersect(sn.refclass.findNonNegative(), inchain);
                            for (double value : intersection) {
                                int s = (int) value;
                                if (!Double.isNaN(this.result.TN.get(i, r))) {
                                    this.metrics.XNchain.set(
                                            0,
                                            c,
                                            this.metrics.XNchain.get(0, c)
                                                    + sn.rt.get(i * sn.nclasses + r, ref * sn.nclasses + s)
                                                    * this.result.TN.get(i, r));
                                }
                            }
                            for (int k = 0; k < inchain.length(); k++) {
                                int s = (int) inchain.get(k);
                                if (!Double.isNaN(this.result.TN.get(i, r))) {
                                    this.metrics.XNchain.set(
                                            0,
                                            c,
                                            this.metrics.XNchain.get(0, c)
                                                    + sn.rt.get(i * sn.nclasses + r, ref * sn.nclasses + s)
                                                    * this.result.TN.get(i, r));
                                }
                            }
                        }
                    }
                }
            }

            // If this is a closed chain we simply apply Little's law
            int nJobsChain = 0;
            for (int i = 0; i < sn.chains.getNumCols(); i++) {
                if (sn.chains.get(c, i) > 0) {
                    nJobsChain += sn.njobs.get(i);
                }
            }

            if (this.model.hasFork() && this.model.hasJoin()) {
                // In this case, CNclass is unreliable as it sums the contribution across all stations,
                // which would include also forked tasks, we use Little's law instead
                this.metrics.CNchain.set(0, c, nJobsChain / this.metrics.XNchain.get(0, c));
            } else {
                // TODO: implementation - note return type and parameters should likely not be void
                if (Double.isInfinite(nJobsChain)) {
                    if (inchain.length() != completingClasses.length()) {
                        throw new RuntimeException(
                                "Edge-based chain definition not yet supported for open queueing networks.");
                    }
                }
                double sumFinite = 0;
                for (int i = 0; i < inchain.length(); i++) {
                    double value = alpha.get((int) sn.refstat.get((int) inchain.get(0)), (int) inchain.get(i))
                            * CNclass.get((int) inchain.get(i));
                    if (Double.isFinite(value))
                        sumFinite += value;
                }

                this.metrics.CNchain.set(0, c, sumFinite);
            }
        }
    }

    // Return table of average system metrics
    public final NetworkAvgSysTable getAvgSysTable() {

        // TODO: provide polymorphic version where handles can be passed in as parameters individually?
        // TODO: provide polymorphic version where handles can be passed in as a whole?

        this.getAvgSys();

        NetworkAvgSysTable avgSysTable = new NetworkAvgSysTable(this.metrics.CNchain.toList1D(), this.metrics.XNchain.toList1D(), this.options);

        java.util.List<String> chainNames = new ArrayList<>();
        java.util.List<String> inChainNames = new ArrayList<>();
        for (int c = 0; c < this.sn.nchains; c++) {
            chainNames.add("Chain" + Integer.toString(c));
            Matrix inchain = sn.inchain.get(c);
            String chainMembers = new String("(");
            for (int i = 0; i < inchain.length(); i++) {
                int r = (int) inchain.get(i);
                if (i == 0) {
                    chainMembers = chainMembers.concat("" + sn.classnames.get(r));
                } else {
                    chainMembers = chainMembers.concat(" " + sn.classnames.get(r));
                }
            }
            inChainNames.add(chainMembers + ")");
        }

        avgSysTable.setChainNames(chainNames);
        avgSysTable.setInChainNames(inChainNames);

        return avgSysTable;

//    // TODO: implementation - note return type and parameters should likely not be void
//    System.out.printf(
//        "%-12s\t %-12s\t %-10s\t %-10s\t", "Chain", "JobClasses", "SysRespT", "SysTput");
//    System.out.println("\n-----------------------------------------------------");
//    NumberFormat nf = NumberFormat.getNumberInstance();
//    nf.setMinimumFractionDigits(5);
//    for (int i = 0; i < sn.nchains; i++) {
//      System.out.format(
//          "%-12s\t %-12s\t %-10s\t %-10s\n",
//          i + 1,
//          sn.jobclasses.get(i).getName(),
//          nf.format(this.metrics.CNchain.get(i)),
//          nf.format(this.metrics.XNchain.get(i)));
//    }
//    System.out.println("-----------------------------------------------------");
    }

    // Return average system response times at steady state
    protected final void getAvgSysRespT() {
        // TODO: implementation - note return type and parameters should likely not be void
        throw new RuntimeException("getAvgSysRespT() has not yet been implemented in JLINE.");
    }

    // Return average system throughputs at steady state
    protected final void getAvgSysTput() {
        // TODO: implementation - note return type and parameters should likely not be void
        throw new RuntimeException("getAvgSysTput() has not yet been implemented in JLINE.");
    }

    // Return transient average station metrics
    public final void getTranAvg() {

        // TODO: getTranHandles (lines 9 to 23)

        // NOTE: This was in LINE, but I believe is legacy as 'matrix' method can provide tran results
    /*    if (!Objects.equals(options.method, "default")) {
      System.err.println(
          "getTranAvg is not offered by the specified method. Setting the solution method to ''closing''.");
      resetResults();
    }
    options.method = "closing";*/

        sn = model.getStruct(true);
        double minRate = sn.rates.elementMin();
        if (!hasTranResults()) {
            if (Double.isInfinite(options.timespan[0]) && Double.isInfinite(options.timespan[1])) {
                options.timespan[0] = 0;
                options.timespan[1] = 30 / minRate;
                System.out.format(
                        "Timespan of transient analysis unspecified, setting the timespan option to [0, %f].\n",
                        options.timespan[1]);
            } else if (Double.isInfinite(options.timespan[0])) {
                options.timespan[0] = 0;
                System.out.format(
                        "Start time of transient analysis unspecified, setting the timespan option to [0, %f].\n",
                        options.timespan[1]);
            } else if (Double.isInfinite(options.timespan[1])) {
                options.timespan[1] = 30 / minRate;
                System.out.format(
                        "End time of transient analysis unspecified, setting the timespan option to [%f, %f].\n",
                        options.timespan[0], options.timespan[1]);
            }
            // TODO: remove init_sol re-initialising
            options.init_sol = new Matrix(0, 0);
            try {
                runAnalyzer();
            } catch (IllegalAccessException e) {
                System.err.println("IllegalAccessException upon running runAnalyzer()");
            } catch (ParserConfigurationException e) {
                System.err.println("ParserConfigurationException upon running runAnalyzer()");
            }
        }

        // TODO: store in Metrics (Lines 52-102)
    }

    // Store average metrics at steady-state
    protected final void setAvgResults() {
        // TODO: implementation - note parameters should likely not be void
        throw new RuntimeException("setAvgResults() has not yet been implemented in JLINE.");
    }

    // Store distribution metrics at steady-state
    protected final void setDistribResults() {
        // TODO: implementation - note parameters should likely not be void
        throw new RuntimeException("setDistribResults() has not yet been implemented in JLINE.");
    }

    // Store transient average metrics
    protected final void setTranProb() {
        // TODO: implementation - note parameters should likely not be void
        throw new RuntimeException("setTranProb() has not yet been implemented in JLINE.");
    }

    // Store transient average metrics
    protected final void setTranAvgResults() {
        // TODO: implementation - note parameters should likely not be void
        throw new RuntimeException("setTranAvgResults() has not yet been implemented in JLINE.");
    }

    // Return a cell array with all Network solvers
    protected static void getAllSolvers() {
        // TODO: implementation - note return type and parameters should likely not be void
        throw new RuntimeException("getAllSolvers() has not yet been implemented in JLINE.");
    }

    // Return a cell array with all Network solvers feasible for this model
    protected static void getAllFeasibleSolvers() {
        // TODO: implementation - note return type and parameters should likely not be void
        throw new RuntimeException("getAllFeasibleSolvers() has not yet been implemented in JLINE.");
    }

    // NOTE: the following LINE methods have not been migrated to JLINE
    // a) updateModel() - model is public and therefore no need for setter
    // b) all methods that are "not supported by this solver" - lack of existing method will suffice
    // rather than dedicated warning that method is not applicable


}
