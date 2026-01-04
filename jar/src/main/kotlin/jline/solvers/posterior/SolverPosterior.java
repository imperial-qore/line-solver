/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.posterior;

import jline.lang.*;
import jline.lang.constant.SolverType;
import jline.lang.nodes.*;
import jline.lang.processes.Distribution;
import jline.lang.processes.Prior;
import jline.solvers.*;
import jline.util.matrix.Matrix;

import java.io.*;
import java.util.*;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

/**
 * Posterior solver for Bayesian-style parameter uncertainty analysis.
 * <p>
 * This solver wraps another solver and handles Prior distributions by expanding
 * the model into a family of networks, one for each alternative in the Prior.
 * Results are aggregated using prior-weighted expectations.
 * <p>
 * Usage:
 * <pre>
 * Prior prior = new Prior(Arrays.asList(new Exp(1.0), new Exp(2.0)), new double[]{0.6, 0.4});
 * queue.setService(jobClass, prior);
 * SolverPosterior solver = new SolverPosterior(model, m -> new SolverMVA(m));
 * AvgTable avgTable = solver.getAvgTable();
 * </pre>
 */
public class SolverPosterior extends EnsembleSolver {

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
    protected SolverResult aggregatedResult;

    /**
     * Creates a SolverPosterior with the given model and solver factory.
     *
     * @param model the network model containing Prior distributions
     * @param solverFactory function to create solvers for each alternative
     */
    public SolverPosterior(Network model, SolverFactory solverFactory) {
        super("SolverPosterior");
        this.options = new SolverOptions(SolverType.POSTERIOR);
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
     * Creates a SolverPosterior with solver options.
     *
     * @param model the network model
     * @param solverFactory function to create solvers
     * @param options solver options
     */
    public SolverPosterior(Network model, SolverFactory solverFactory, SolverOptions options) {
        this(model, solverFactory);
        this.options = options;
    }

    /**
     * Returns default solver options.
     */
    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.POSTERIOR);
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
        return priorInfo != null ? priorInfo.prior.getNumAlternatives() : 0;
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
        Prior prior = priorInfo.prior;
        int n = prior.getNumAlternatives();

        for (int i = 0; i < n; i++) {
            // Deep copy the original model
            Network modelCopy = deepCopyNetwork(originalModel);

            // Get the node and class from the copy
            List<Node> nodes = modelCopy.getNodes();
            List<JobClass> classes = modelCopy.getClasses();
            Node node = nodes.get(priorInfo.nodeIdx);
            JobClass jobClass = classes.get(priorInfo.classIdx);

            // Replace Prior with the concrete distribution
            Distribution concreteDist = prior.getAlternative(i);

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
        Prior prior = priorInfo.prior;
        int n = prior.getNumAlternatives();
        double[] probs = prior.getProbabilities();

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

        return avgTable;
    }

    /**
     * Returns the prior-weighted average results.
     */
    public AvgTable getAvgTable() {
        return getEnsembleAvg();
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
        // Posterior solver supports any model that has a Prior distribution
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
                System.out.printf("%d\t\t%.4f\t\t%s\t\t%s\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f%n",
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
                System.out.printf("%.4f\t\t%.4f\t\t%.4f%n", values[i], probabilities[i], cdf[i]);
            }
        }
    }
}
