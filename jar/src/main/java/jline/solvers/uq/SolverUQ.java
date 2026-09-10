/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.uq;

import jline.api.pfqn.mva.Pfqn_mva_interval;
import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.*;
import jline.lang.processes.Distribution;
import jline.lang.processes.Prior;
import jline.solvers.*;
import jline.util.matrix.Matrix;

import java.io.*;
import java.util.*;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

/**
 * UQ solver for Bayesian-style parameter uncertainty analysis.
 * <p>
 * This solver wraps another solver and handles Prior distributions by expanding
 * the model into a family of networks, one for each alternative in the Prior.
 * Results are aggregated using prior-weighted expectations.
 * <p>
 * Usage:
 * <pre>
 * Prior prior = new Prior(Arrays.asList(new Exp(1.0), new Exp(2.0)), new double[]{0.6, 0.4});
 * queue.setService(jobClass, prior);
 * SolverUQ solver = new SolverUQ(model, m -> new SolverMVA(m));
 * AvgTable avgTable = solver.getAvgTable();
 * </pre>
 */
public class SolverUQ extends EnsembleSolver {

    /**
     * Functional interface for creating solvers.
     */
    @FunctionalInterface
    public interface SolverFactory {
        NetworkSolver create(Network model);
    }

    /**
     * Information about a detected Prior distribution.
     */
    public static class PriorInfo {
        public final int nodeIdx;
        public final int classIdx;
        public final String type;  // "service" or "arrival"
        public final Prior prior;

        public PriorInfo(int nodeIdx, int classIdx, String type, Prior prior) {
            this.nodeIdx = nodeIdx;
            this.classIdx = classIdx;
            this.type = type;
            this.prior = prior;
        }
    }

    protected Network originalModel;
    protected SolverFactory solverFactory;
    protected PriorInfo priorInfo;

    /**
     * The DESIGN: which alternatives are solved and with what weight.
     *
     * The reference reduces the detected Priors to weighted design points before
     * anything is solved (UQ.buildDesign), and options.method chooses HOW:
     * 'default'/'discrete'/'quadrature' expand a discrete Prior exactly, while
     * 'montecarlo' draws options.samples of them against the prior probabilities
     * and weights each 1/n. This class used to read prior.getAlternative(i) and
     * prior.getProbabilities() directly, so it had one design and options.method
     * selected nothing.
     */
    protected Prior.Design design;
    protected SolverResult aggregatedResult;

    /**
     * Creates a SolverUQ with the given model and solver factory.
     *
     * @param model the network model containing Prior distributions
     * @param solverFactory function to create solvers for each alternative
     */
    public SolverUQ(Network model, SolverFactory solverFactory) {
        super("SolverUQ");
        this.options = new SolverOptions(SolverType.UQ);
        this.originalModel = model;
        this.solverFactory = solverFactory;
        this.results = new HashMap<>();

        // Detect Prior distributions
        this.priorInfo = detectPrior();
        if (this.priorInfo == null) {
            line_error(mfilename(new Object(){}), "No Prior distribution found in model");
        }

        // Initialize ensemble arrays
        int n = priorInfo.prior.getNumAlternatives();
        this.ensemble = new Network[n];
        this.solvers = new NetworkSolver[n];
    }

    /**
     * Creates a SolverUQ with solver options.
     *
     * @param model the network model
     * @param solverFactory function to create solvers
     * @param options solver options
     */
    public SolverUQ(Network model, SolverFactory solverFactory, SolverOptions options) {
        this(model, solverFactory);
        this.options = options;
    }

    /**
     * Returns default solver options.
     */
    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.UQ);
    }

    /**
     * Detects Prior distributions in the model.
     * Currently supports only a single Prior per model.
     *
     * @return PriorInfo for the detected Prior, or null if none found
     */
    protected PriorInfo detectPrior() {
        List<Node> nodes = originalModel.getNodes();
        List<JobClass> classes = originalModel.getClasses();

        // Check service distributions at ServiceStations
        for (int nodeIdx = 0; nodeIdx < nodes.size(); nodeIdx++) {
            Node node = nodes.get(nodeIdx);
            if (node instanceof ServiceStation) {
                ServiceStation station = (ServiceStation) node;
                for (int classIdx = 0; classIdx < classes.size(); classIdx++) {
                    JobClass jobClass = classes.get(classIdx);
                    Distribution dist = station.getServiceProcess(jobClass);
                    if (dist instanceof Prior) {
                        return new PriorInfo(nodeIdx, classIdx, "service", (Prior) dist);
                    }
                }
            }
        }

        // Check arrival distributions at Sources
        for (int nodeIdx = 0; nodeIdx < nodes.size(); nodeIdx++) {
            Node node = nodes.get(nodeIdx);
            if (node instanceof Source) {
                Source source = (Source) node;
                for (int classIdx = 0; classIdx < classes.size(); classIdx++) {
                    JobClass jobClass = classes.get(classIdx);
                    Distribution dist = source.getArrivalDistribution(jobClass);
                    if (dist instanceof Prior) {
                        return new PriorInfo(nodeIdx, classIdx, "arrival", (Prior) dist);
                    }
                }
            }
        }

        return null;
    }

    /**
     * Checks if the model has a Prior distribution.
     */
    public boolean hasPriorDistribution() {
        return priorInfo != null;
    }

    /**
     * Returns the number of alternatives in the Prior.
     */
    public int getNumAlternatives() {
        if (design != null) {
            return design.dists.size();
        }
        return priorInfo != null ? priorInfo.prior.getNumAlternatives() : 0;
    }

    /**
     * Valid methods for this solver, UQ.listValidMethods in MATLAB verbatim.
     *
     * They name the DESIGN, i.e. how the Prior is reduced to weighted design
     * points, and not the inner solver, which is chosen by the SolverFactory
     * this class is constructed with.
     *
     * @return the design methods SolverUQ accepts
     */
    public String[] listValidMethods() {
        return new String[]{"default", "discrete", "quadrature", "montecarlo"};
    }

    /**
     * Resolve the discretization method from the solver options.
     *
     * 'default' keeps the historical behaviour: a discrete Prior is expanded as
     * given. Mirrors UQ.getUQMethod.
     *
     * @return "quadrature" or "montecarlo"
     */
    public String getUQMethod() {
        String m = (options == null || options.method == null) ? "default" : options.method;
        if ("default".equals(m) || "discrete".equals(m) || "quadrature".equals(m)) {
            return "quadrature";
        }
        if ("montecarlo".equals(m)) {
            return "montecarlo";
        }
        line_error(mfilename(new Object() {
        }), "Unknown UQ method: " + m);
        return "quadrature";
    }

    /**
     * Number of nodes per Prior, from options.samples (UQ.getUQNodes).
     *
     * @return the design size for the Monte Carlo method
     */
    public int getUQNodes() {
        if (options == null || options.samples <= 0) {
            return 11;
        }
        return (int) Math.round(options.samples);
    }

    /**
     * Reduce the detected Prior to weighted design points.
     *
     * Each design point assigns one concrete Distribution to the Prior and
     * carries the weight of that assignment. Mirrors UQ.buildDesign; this solver
     * detects a single Prior, so there is no tensor product to take.
     */
    protected void buildDesign() {
        if (priorInfo == null) {
            return;
        }
        String method = getUQMethod();
        Random rng = new Random(options == null ? 23000 : options.seed);
        this.design = priorInfo.prior.discretize(getUQNodes(), method, rng);
    }

    @Override
    public int getNumberOfModels() {
        if (ensemble != null && ensemble.length > 0 && ensemble[0] != null) {
            return ensemble.length;
        }
        return getNumAlternatives();
    }

    /**
     * Deep copies a Network using serialization.
     */
    @SuppressWarnings("unchecked")
    protected Network deepCopyNetwork(Network original) {
        try {
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            ObjectOutputStream oos = new ObjectOutputStream(baos);
            oos.writeObject(original);
            oos.close();

            ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
            ObjectInputStream ois = new ObjectInputStream(bais);
            Network copy = (Network) ois.readObject();
            ois.close();

            return copy;
        } catch (IOException | ClassNotFoundException e) {
            throw new RuntimeException("Failed to deep copy network: " + e.getMessage(), e);
        }
    }

    @Override
    protected void init() {
        if (design == null) {
            buildDesign();
        }
        int n = design.dists.size();

        for (int i = 0; i < n; i++) {
            // Deep copy the original model
            Network modelCopy = deepCopyNetwork(originalModel);
            // setService does NOT invalidate a cached struct (ServiceStation:
            // deliberate, for SolverLN's iteration cost) and serialization copies
            // the cache, so without this every design point would be solved with
            // the Prior's MIXTURE moments
            modelCopy.resetStruct();

            // Get the node and class from the copy
            List<Node> nodes = modelCopy.getNodes();
            List<JobClass> classes = modelCopy.getClasses();
            Node node = nodes.get(priorInfo.nodeIdx);
            JobClass jobClass = classes.get(priorInfo.classIdx);

            // Replace Prior with the concrete distribution of this design point
            Distribution concreteDist = design.dists.get(i);

            if ("service".equals(priorInfo.type)) {
                ((ServiceStation) node).setService(jobClass, concreteDist);
            } else if ("arrival".equals(priorInfo.type)) {
                ((Source) node).setArrival(jobClass, concreteDist);
            }

            // Store the model and create solver
            ensemble[i] = modelCopy;
            solvers[i] = solverFactory.create(modelCopy);
        }
    }

    @Override
    protected void pre(int it) {
        // No pre-processing needed for single iteration
    }

    @Override
    protected SolverResult analyze(int it, int e) {
        try {
            // Use getAvg() which is the public API for running analysis
            return solvers[e].getAvg();
        } catch (Exception ex) {
            throw new RuntimeException("Failed to analyze model " + e + ": " + ex.getMessage(), ex);
        }
    }

    @Override
    protected void post(int it) {
        // Aggregate results using prior probabilities
        aggregateResults();
    }

    @Override
    protected void finish() {
        // Cleanup if needed
    }

    @Override
    protected boolean converged(int it) {
        // Single iteration only
        return it >= 1;
    }

    /**
     * Aggregates results from all alternatives using prior weights.
     */
    protected void aggregateResults() {
        int n = design.dists.size();
        double[] probs = design.weights;

        // Get dimensions from first solver's result
        SolverResult firstResult = solvers[0].result;
        if (firstResult == null || firstResult.QN == null) {
            return;
        }

        int M = firstResult.QN.getNumRows();  // stations
        int K = firstResult.QN.getNumCols();  // classes

        // Initialize aggregated matrices
        Matrix QN = Matrix.zeros(M, K);
        Matrix UN = Matrix.zeros(M, K);
        Matrix RN = Matrix.zeros(M, K);
        Matrix TN = Matrix.zeros(M, K);
        Matrix AN = Matrix.zeros(M, K);
        Matrix WN = Matrix.zeros(M, K);

        // Aggregate with prior weights
        for (int i = 0; i < n; i++) {
            SolverResult r = solvers[i].result;
            double p = probs[i];

            if (r.QN != null) {
                for (int m = 0; m < M; m++) {
                    for (int k = 0; k < K; k++) {
                        QN.set(m, k, QN.get(m, k) + p * r.QN.get(m, k));
                    }
                }
            }
            if (r.UN != null) {
                for (int m = 0; m < M; m++) {
                    for (int k = 0; k < K; k++) {
                        UN.set(m, k, UN.get(m, k) + p * r.UN.get(m, k));
                    }
                }
            }
            if (r.RN != null) {
                for (int m = 0; m < M; m++) {
                    for (int k = 0; k < K; k++) {
                        RN.set(m, k, RN.get(m, k) + p * r.RN.get(m, k));
                    }
                }
            }
            if (r.TN != null) {
                for (int m = 0; m < M; m++) {
                    for (int k = 0; k < K; k++) {
                        TN.set(m, k, TN.get(m, k) + p * r.TN.get(m, k));
                    }
                }
            }
            if (r.AN != null) {
                for (int m = 0; m < M; m++) {
                    for (int k = 0; k < K; k++) {
                        AN.set(m, k, AN.get(m, k) + p * r.AN.get(m, k));
                    }
                }
            }
            if (r.WN != null) {
                for (int m = 0; m < M; m++) {
                    for (int k = 0; k < K; k++) {
                        WN.set(m, k, WN.get(m, k) + p * r.WN.get(m, k));
                    }
                }
            }
        }

        // Store aggregated result
        aggregatedResult = new SolverResult();
        aggregatedResult.QN = QN;
        aggregatedResult.UN = UN;
        aggregatedResult.RN = RN;
        aggregatedResult.TN = TN;
        aggregatedResult.AN = AN;
        aggregatedResult.WN = WN;
        this.result = aggregatedResult;
    }

    @Override
    protected AvgTable getEnsembleAvg() {
        if (aggregatedResult == null) {
            iterate();
        }
        // Build NetworkAvgTable from aggregated results
        List<Station> stations = originalModel.getStations();
        List<JobClass> classes = originalModel.getClasses();
        int M = stations.size();
        int K = classes.size();

        List<Double> Qval = new ArrayList<>();
        List<Double> Uval = new ArrayList<>();
        List<Double> Rval = new ArrayList<>();
        List<Double> Wval = new ArrayList<>();  // Residence time (same as RespT for most)
        List<Double> Aval = new ArrayList<>();
        List<Double> Tval = new ArrayList<>();

        for (int m = 0; m < M; m++) {
            for (int k = 0; k < K; k++) {
                Qval.add(aggregatedResult.QN != null ? aggregatedResult.QN.get(m, k) : Double.NaN);
                Uval.add(aggregatedResult.UN != null ? aggregatedResult.UN.get(m, k) : Double.NaN);
                Rval.add(aggregatedResult.RN != null ? aggregatedResult.RN.get(m, k) : Double.NaN);
                Wval.add(aggregatedResult.WN != null ? aggregatedResult.WN.get(m, k) : Double.NaN);
                Aval.add(aggregatedResult.AN != null ? aggregatedResult.AN.get(m, k) : Double.NaN);
                Tval.add(aggregatedResult.TN != null ? aggregatedResult.TN.get(m, k) : Double.NaN);
            }
        }

        NetworkAvgTable avgTable = new NetworkAvgTable(Qval, Uval, Rval, Wval, Aval, Tval);
        // Set station and class names
        List<String> stationNames = new ArrayList<>();
        List<String> classNames = new ArrayList<>();
        for (int m = 0; m < M; m++) {
            for (int k = 0; k < K; k++) {
                stationNames.add(stations.get(m).getName());
                classNames.add(classes.get(k).getName());
            }
        }
        avgTable.setStationNames(stationNames);
        avgTable.setClassNames(classNames);
        avgTable.setOptions(this.options);

        return avgTable;
    }

    /**
     * Returns the prior-weighted average results.
     */
    public AvgTable getAvgTable() {
        return jline.io.LineResultRecorder.around(this, "avg", () -> getAvgTableImpl());
    }

    /**
     * Body of {@link #getAvgTable()}, split out so {@link jline.io.LineResultRecorder}
     * sees what the getter RETURNED. The JAVA cross-codebase parity row is
     * measured from that rather than from what an example printed.
     */
    protected AvgTable getAvgTableImpl() {
        return getEnsembleAvg();
    }

    /**
     * Range of every metric over the support of the Prior, weights discarded.
     * <p>
     * This is the epistemic case in which the modeller can bound a parameter but not
     * distribute it. Two regimes, distinguished by {@link Interval#exact}:
     * <ul>
     * <li>exact: the model is a single-class closed product-form network with
     * load-independent single-server queues and delays, so
     * {@link Pfqn_mva_interval} returns the exact hull of MVA over the whole demand box
     * by the monotonicity of Luthi and Haring (1998). No ensemble run is needed and the
     * interval is attained, not sampled.</li>
     * <li>not exact: fallback to the range across the alternatives that were solved. The
     * JAR Prior is discrete, so that range is again the whole support; it is reported as
     * inexact only because the monotonicity theorems do not apply to the model.</li>
     * </ul>
     * The interval is conditional on the true parameters lying inside the Prior support.
     * It is not a bound on the exact solution of the network and must not be composed
     * with SolverBA brackets.
     *
     * @return the interval-valued metrics
     */
    public Interval getInterval() {
        String why = qualifiesForIntervalMVA();
        if (why == null) {
            return intervalByMVA();
        }
        return intervalBySampling(why);
    }

    /**
     * Whether the monotonicity theorems behind {@link Pfqn_mva_interval} hold for this
     * model.
     *
     * @return null when they do, otherwise the first violated condition
     */
    public String qualifiesForIntervalMVA() {
        if (priorInfo != null && !"service".equals(priorInfo.type)) {
            return "a Prior sits on an arrival process, so the model is open";
        }
        NetworkStruct sn = originalModel.getStruct(false);
        if (sn.nclasses != 1) {
            return "the theorems are proved for a single class only";
        }
        if (sn.nclosedjobs <= 0) {
            return "the class is not closed";
        }
        if (sn.nnodes != sn.nstations) {
            return "the model has nodes that are not stations";
        }
        for (int i = 0; i < sn.nstations; i++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(i));
            if (s == SchedStrategy.INF) {
                continue;
            }
            if (sn.nservers.get(i, 0) > 1) {
                return "a queueing station has more than one server";
            }
            if (s != SchedStrategy.PS && s != SchedStrategy.FCFS) {
                return "a station is neither delay, PS nor FCFS";
            }
        }
        return null;
    }

    /**
     * Exact hull through {@link Pfqn_mva_interval}. The demand box is the nominal demand
     * vector with the prior-carrying station widened to the range of mean service times
     * over the Prior support.
     */
    public Interval intervalByMVA() {
        NetworkStruct sn = originalModel.getStruct(false);
        int M = sn.nstations;
        Matrix visits = sn.visits.get(0);

        double[] V = new double[M];
        double[] STlo = new double[M];
        double[] STup = new double[M];
        boolean[] isDelay = new boolean[M];
        for (int i = 0; i < M; i++) {
            V[i] = visits.get(i, 0);
            double rate = sn.rates.get(i, 0);
            double st = rate > 0 ? 1.0 / rate : 0.0;
            STlo[i] = st;
            STup[i] = st;
            isDelay[i] = sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF;
        }

        if (priorInfo != null) {
            // linear index: nodeToStation is stored as a row vector here, not nnodes x 1
            int ist = (int) Math.round(sn.nodeToStation.get(priorInfo.nodeIdx));
            double lo = Double.POSITIVE_INFINITY;
            double up = Double.NEGATIVE_INFINITY;
            for (int a = 0; a < priorInfo.prior.getNumAlternatives(); a++) {
                double m = priorInfo.prior.getAlternative(a).getMean();
                lo = Math.min(lo, m);
                up = Math.max(up, m);
            }
            STlo[ist] = lo;
            STup[ist] = up;
        }

        double zlo = 0.0;
        double zup = 0.0;
        int nq = 0;
        for (int i = 0; i < M; i++) {
            if (isDelay[i]) {
                zlo += V[i] * STlo[i];
                zup += V[i] * STup[i];
            } else {
                nq++;
            }
        }
        Matrix Lint = new Matrix(nq, 2);
        int[] qidx = new int[nq];
        int c = 0;
        for (int i = 0; i < M; i++) {
            if (!isDelay[i]) {
                qidx[c] = i;
                Lint.set(c, 0, V[i] * STlo[i]);
                Lint.set(c, 1, V[i] * STup[i]);
                c++;
            }
        }
        Matrix Zint = new Matrix(1, 2);
        Zint.set(0, 0, zlo);
        Zint.set(0, 1, zup);
        Matrix Nint = new Matrix(1, 1);
        Nint.set(0, 0, sn.nclosedjobs);

        Pfqn_mva_interval.Result res = Pfqn_mva_interval.pfqn_mva_interval(Lint, Nint, Zint);

        Interval ival = new Interval(M, 1);
        ival.X = res.X;
        ival.Rtot = res.Rtot;
        ival.exact = true;
        ival.method = "mvainterval";
        for (int j = 0; j < nq; j++) {
            int i = qidx[j];
            ival.Qlo.set(i, 0, res.Q.get(j, 0));
            ival.Qup.set(i, 0, res.Q.get(j, 1));
            ival.Ulo.set(i, 0, res.U.get(j, 0));
            ival.Uup.set(i, 0, res.U.get(j, 1));
            ival.Wlo.set(i, 0, res.R.get(j, 0));
            ival.Wup.set(i, 0, res.R.get(j, 1));
            ival.Rlo.set(i, 0, V[i] > 0 ? res.R.get(j, 0) / V[i] : 0.0);
            ival.Rup.set(i, 0, V[i] > 0 ? res.R.get(j, 1) / V[i] : 0.0);
        }
        // A delay station never queues, so its residence time is its own demand interval
        // and its population is the throughput times that demand, enclosed as a product of
        // two intervals.
        for (int i = 0; i < M; i++) {
            if (!isDelay[i]) {
                continue;
            }
            ival.Wlo.set(i, 0, V[i] * STlo[i]);
            ival.Wup.set(i, 0, V[i] * STup[i]);
            ival.Rlo.set(i, 0, STlo[i]);
            ival.Rup.set(i, 0, STup[i]);
            ival.Qlo.set(i, 0, res.X.get(0, 0) * V[i] * STlo[i]);
            ival.Qup.set(i, 0, res.X.get(0, 1) * V[i] * STup[i]);
            ival.Ulo.set(i, 0, ival.Qlo.get(i, 0));
            ival.Uup.set(i, 0, ival.Qup.get(i, 0));
        }
        for (int i = 0; i < M; i++) {
            ival.Tlo.set(i, 0, V[i] * res.X.get(0, 0));
            ival.Tup.set(i, 0, V[i] * res.X.get(0, 1));
        }
        return ival;
    }

    /**
     * Range of each metric across the alternatives that were solved.
     *
     * @param why the condition that ruled out the exact path, kept on the result
     */
    public Interval intervalBySampling(String why) {
        if (solvers[0] == null || solvers[0].result == null) {
            iterate();
        }
        int n = getNumberOfModels();
        SolverResult first = solvers[0].result;
        int M = first.QN.getNumRows();
        int K = first.QN.getNumCols();

        Interval ival = new Interval(M, K);
        ival.exact = false;
        ival.method = "sampled";
        ival.reason = why;
        boolean seeded = false;
        for (int e = 0; e < n; e++) {
            SolverResult r = solvers[e].result;
            if (r == null) {
                continue;
            }
            accumulate(ival.Qlo, ival.Qup, r.QN, seeded);
            accumulate(ival.Ulo, ival.Uup, r.UN, seeded);
            accumulate(ival.Rlo, ival.Rup, r.RN, seeded);
            accumulate(ival.Tlo, ival.Tup, r.TN, seeded);
            accumulate(ival.Wlo, ival.Wup, r.WN, seeded);
            seeded = true;
        }
        return ival;
    }

    private static void accumulate(Matrix lo, Matrix up, Matrix v, boolean seeded) {
        if (v == null) {
            return;
        }
        for (int m = 0; m < lo.getNumRows(); m++) {
            for (int k = 0; k < lo.getNumCols(); k++) {
                double x = v.get(m, k);
                if (!seeded) {
                    lo.set(m, k, x);
                    up.set(m, k, x);
                } else {
                    lo.set(m, k, Math.min(lo.get(m, k), x));
                    up.set(m, k, Math.max(up.get(m, k), x));
                }
            }
        }
    }

    /**
     * Interval-valued metrics returned by {@link SolverUQ#getInterval()}. Each metric is
     * a pair of nstations x nclasses matrices holding the two endpoints.
     */
    public static class Interval {
        /** Queue-length endpoints. */
        public final Matrix Qlo, Qup;
        /** Utilization endpoints. */
        public final Matrix Ulo, Uup;
        /** Response-time endpoints. */
        public final Matrix Rlo, Rup;
        /** Throughput endpoints. */
        public final Matrix Tlo, Tup;
        /** Residence-time endpoints. */
        public final Matrix Wlo, Wup;
        /** System throughput interval (1 x 2); null off the exact path. */
        public Matrix X;
        /** Total response-time interval (1 x 2); null off the exact path. */
        public Matrix Rtot;
        /** Whether the interval is the exact hull over the whole parameter box. */
        public boolean exact;
        /** "mvainterval" or "sampled". */
        public String method;
        /** When not exact, the condition that ruled out the exact path. */
        public String reason;

        public Interval(int M, int K) {
            this.Qlo = Matrix.zeros(M, K);
            this.Qup = Matrix.zeros(M, K);
            this.Ulo = Matrix.zeros(M, K);
            this.Uup = Matrix.zeros(M, K);
            this.Rlo = Matrix.zeros(M, K);
            this.Rup = Matrix.zeros(M, K);
            this.Tlo = Matrix.zeros(M, K);
            this.Tup = Matrix.zeros(M, K);
            this.Wlo = Matrix.zeros(M, K);
            this.Wup = Matrix.zeros(M, K);
        }
    }

    /**
     * Returns a table with per-alternative results and probabilities.
     *
     * @return PosteriorTable with all alternatives
     */
    public PosteriorTable getPosteriorTable() {
        if (solvers[0] == null || solvers[0].result == null) {
            iterate();
        }

        Prior prior = priorInfo.prior;
        int n = prior.getNumAlternatives();

        List<PosteriorTableRow> rows = new ArrayList<>();
        List<Station> stations = originalModel.getStations();
        List<JobClass> classes = originalModel.getClasses();

        for (int i = 0; i < n; i++) {
            SolverResult r = solvers[i].result;
            double prob = prior.getProbability(i);

            for (int m = 0; m < stations.size(); m++) {
                for (int k = 0; k < classes.size(); k++) {
                    double Q = r.QN != null ? r.QN.get(m, k) : Double.NaN;
                    double U = r.UN != null ? r.UN.get(m, k) : Double.NaN;
                    double R = r.RN != null ? r.RN.get(m, k) : Double.NaN;
                    double T = r.TN != null ? r.TN.get(m, k) : Double.NaN;
                    double A = r.AN != null ? r.AN.get(m, k) : Double.NaN;

                    rows.add(new PosteriorTableRow(
                            i, prob, stations.get(m).getName(), classes.get(k).getName(),
                            Q, U, R, T, A
                    ));
                }
            }
        }

        return new PosteriorTable(rows);
    }

    /**
     * Returns the posterior distribution for a specific metric at a station/class.
     *
     * @param metric the metric name ("Q", "U", "R", "T", "A")
     * @param station the station
     * @param jobClass the job class
     * @return EmpiricalCDF representing the posterior distribution
     */
    public EmpiricalCDF getPosteriorDist(String metric, Station station, JobClass jobClass) {
        if (solvers[0] == null || solvers[0].result == null) {
            iterate();
        }

        Prior prior = priorInfo.prior;
        int n = prior.getNumAlternatives();

        // Find station and class indices
        int stationIdx = originalModel.getStations().indexOf(station);
        int classIdx = originalModel.getClasses().indexOf(jobClass);

        if (stationIdx < 0 || classIdx < 0) {
            throw new IllegalArgumentException("Station or class not found in model");
        }

        // Collect values and probabilities
        double[] values = new double[n];
        double[] probs = prior.getProbabilities();

        for (int i = 0; i < n; i++) {
            SolverResult r = solvers[i].result;
            switch (metric.toUpperCase()) {
                case "Q":
                    values[i] = r.QN != null ? r.QN.get(stationIdx, classIdx) : Double.NaN;
                    break;
                case "U":
                    values[i] = r.UN != null ? r.UN.get(stationIdx, classIdx) : Double.NaN;
                    break;
                case "R":
                    values[i] = r.RN != null ? r.RN.get(stationIdx, classIdx) : Double.NaN;
                    break;
                case "T":
                    values[i] = r.TN != null ? r.TN.get(stationIdx, classIdx) : Double.NaN;
                    break;
                case "A":
                    values[i] = r.AN != null ? r.AN.get(stationIdx, classIdx) : Double.NaN;
                    break;
                default:
                    throw new IllegalArgumentException("Unknown metric: " + metric);
            }
        }

        return new EmpiricalCDF(values, probs);
    }

    @Override
    public boolean supports(Network model) {
        // UQ solver supports any model that has a Prior distribution
        return detectPrior() != null;
    }

    /**
     * Runs the posterior analysis.
     */
    public void runAnalyzer() {
        iterate();
    }

    /**
     * Row in the posterior table.
     */
    public static class PosteriorTableRow {
        public final int alternativeIdx;
        public final double probability;
        public final String station;
        public final String jobClass;
        public final double Q;
        public final double U;
        public final double R;
        public final double T;
        public final double A;

        public PosteriorTableRow(int alternativeIdx, double probability, String station,
                                  String jobClass, double Q, double U, double R, double T, double A) {
            this.alternativeIdx = alternativeIdx;
            this.probability = probability;
            this.station = station;
            this.jobClass = jobClass;
            this.Q = Q;
            this.U = U;
            this.R = R;
            this.T = T;
            this.A = A;
        }
    }

    /**
     * Table containing per-alternative posterior results.
     */
    public static class PosteriorTable {
        public final List<PosteriorTableRow> rows;

        public PosteriorTable(List<PosteriorTableRow> rows) {
            this.rows = rows;
        }

        public void print() {
            System.out.println("Alternative\tProbability\tStation\t\tClass\t\tQ\t\tU\t\tR\t\tT\t\tA");
            for (PosteriorTableRow row : rows) {
                System.out.printf("%d\t\t%.4f\t\t%s\t\t%s\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\n",
                        row.alternativeIdx, row.probability, row.station, row.jobClass,
                        row.Q, row.U, row.R, row.T, row.A);
            }
        }
    }

    /**
     * Empirical CDF representing a discrete posterior distribution.
     */
    public static class EmpiricalCDF {
        public final double[] values;
        public final double[] probabilities;
        public final double[] cdf;

        public EmpiricalCDF(double[] values, double[] probabilities) {
            // Sort by values and compute CDF
            int n = values.length;
            Integer[] indices = new Integer[n];
            for (int i = 0; i < n; i++) indices[i] = i;

            final double[] v = values;
            Arrays.sort(indices, Comparator.comparingDouble(i -> v[i]));

            this.values = new double[n];
            this.probabilities = new double[n];
            this.cdf = new double[n];

            double cumProb = 0;
            for (int i = 0; i < n; i++) {
                this.values[i] = values[indices[i]];
                this.probabilities[i] = probabilities[indices[i]];
                cumProb += this.probabilities[i];
                this.cdf[i] = cumProb;
            }
        }

        /**
         * Evaluates the CDF at point x.
         */
        public double evalCDF(double x) {
            for (int i = 0; i < values.length; i++) {
                if (x < values[i]) {
                    return i > 0 ? cdf[i - 1] : 0.0;
                }
            }
            return 1.0;
        }

        /**
         * Returns the mean of the distribution.
         */
        public double getMean() {
            double mean = 0;
            for (int i = 0; i < values.length; i++) {
                mean += values[i] * probabilities[i];
            }
            return mean;
        }

        public void print() {
            System.out.println("Value\t\tProbability\tCDF");
            for (int i = 0; i < values.length; i++) {
                System.out.printf("%.4f\t\t%.4f\t\t%.4f\n", values[i], probabilities[i], cdf[i]);
            }
        }
    }
}
