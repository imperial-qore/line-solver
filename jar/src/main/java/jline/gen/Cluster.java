/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.gen;

import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Source;
import jline.lang.processes.APH;
import jline.solvers.NetworkAvgTable;
import jline.solvers.NetworkSolver;
import jline.util.matrix.Matrix;

import java.lang.reflect.Constructor;
import java.util.List;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.function.Function;

/**
 * Builder for cluster models with helpers to compare dispatching/scheduling
 * policies and to sweep parameters.
 *
 * <p>The builder produces models with the same topology as
 * {@link Network#cluster(Matrix, Matrix, SchedStrategy[], Matrix, RoutingStrategy)}
 * (open) or {@link Network#clusterClosed(Matrix, Matrix, Matrix, SchedStrategy[], Matrix, RoutingStrategy)}
 * (closed), and adds:
 * <ul>
 *   <li>{@link #compareDispatching(Class, RoutingStrategy...)} — solve the model under
 *       multiple dispatching policies and return the average tables side-by-side;</li>
 *   <li>{@link #compareScheduling(Class, SchedStrategy...)} — analogous comparison over
 *       per-server scheduling disciplines;</li>
 *   <li>{@link #sweepArrivalRate(double[], Class)} — open-model arrival-rate sweep;</li>
 *   <li>{@link #sweepNumStations(int[], Class)} — sweep over the number of parallel servers.</li>
 * </ul>
 *
 * <p>Solvers are passed by class reference (e.g. {@code SolverMVA.class}) and instantiated
 * via the conventional {@code Solver(Network)} constructor.
 */
public class Cluster {

    private int numStations;
    private double[] arrivalRates;
    private double[][] serviceRates;
    private int[] stationCounts;
    private SchedStrategy scheduling;
    private RoutingStrategy dispatching;
    private boolean closed;
    private int[] population;
    private double[] thinkTimes;
    private double[][] dispatchProbs;   // PROB: rows = classes (1 broadcasts), cols = servers
    private int[] dispatchWeights;      // WRROBIN: per-server integer weights
    private Integer sqD;          // SQ: d, sampled destinations
    private double[] arrivalScvs;       // per-class arrival SCVs (1 = exponential, default)
    private double[][] serviceScvs;     // (M, R) service SCVs (1 = exponential, default)

    /**
     * Creates a default cluster: 2 servers, single open class with arrival rate 1.0, service
     * rate 1.0 at each server, PS scheduling, RAND dispatching. Configure further via the
     * chainable {@code set*} methods (e.g. {@link #setNumStations(int)},
     * {@link #setArrivalRate(double)}, {@link #setServiceRate(double)}).
     */
    public Cluster() {
        this.numStations = 2;
        this.arrivalRates = new double[]{1.0};
        this.serviceRates = fillUniform(2, 1.0);
        this.stationCounts = new int[]{1, 1};
        this.scheduling = SchedStrategy.PS;
        this.dispatching = RoutingStrategy.RAND;
        this.closed = false;
    }

    /**
     * Sets the number of parallel server queues. Replicates the current single-class
     * service rate across all servers and resets the per-server multiplicity to all-1.
     */
    public Cluster setNumStations(int M) {
        if (M <= 0) throw new IllegalArgumentException("numStations must be positive");
        double sample = serviceRates.length > 0 && serviceRates[0].length > 0
                ? serviceRates[0][0] : 1.0;
        this.numStations = M;
        this.serviceRates = fillUniform(M, sample);
        this.stationCounts = new int[M];
        for (int i = 0; i < M; i++) this.stationCounts[i] = 1;
        return this;
    }

    /** Single-class arrival rate. */
    public Cluster setArrivalRate(double lambda) {
        if (lambda <= 0) throw new IllegalArgumentException("arrival rate must be positive");
        this.arrivalRates = new double[]{lambda};
        return this;
    }

    /** Per-class arrival rates. */
    public Cluster setArrivalRates(double[] lambdas) {
        for (double l : lambdas) {
            if (l <= 0) throw new IllegalArgumentException("arrival rate must be positive");
        }
        this.arrivalRates = lambdas.clone();
        return this;
    }

    /** Single service rate, broadcast across all (server, class) pairs. */
    public Cluster setServiceRate(double mu) {
        if (mu <= 0) throw new IllegalArgumentException("service rate must be positive");
        this.serviceRates = fillUniform(numStations, mu);
        return this;
    }

    /** Per-(server, class) service rates as a {@code (numStations x R)} matrix. */
    public Cluster setServiceRates(double[][] rates) {
        if (rates.length != numStations) {
            throw new IllegalArgumentException("serviceRates outer dim must equal numStations");
        }
        for (double[] row : rates) {
            for (double r : row) {
                if (r <= 0) throw new IllegalArgumentException("service rate must be positive");
            }
        }
        this.serviceRates = deepCopy(rates);
        return this;
    }

    /** Sets the dispatching strategy applied at the router. */
    public Cluster setDispatching(RoutingStrategy d) {
        this.dispatching = d;
        return this;
    }

    /** Sets the scheduling discipline used at every server. */
    public Cluster setScheduling(SchedStrategy s) {
        this.scheduling = s;
        return this;
    }

    /**
     * Configures probabilistic dispatching with per-server probabilities (broadcast to all
     * classes). Each call also flips the dispatching strategy to {@link RoutingStrategy#PROB}.
     */
    public Cluster setProbabilities(double[] probs) {
        if (probs.length != numStations) {
            throw new IllegalArgumentException("probs length must equal numStations");
        }
        this.dispatchProbs = new double[][]{probs.clone()};
        this.dispatching = RoutingStrategy.PROB;
        return this;
    }

    /**
     * Configures probabilistic dispatching with per-class, per-server probabilities. Rows index
     * classes (length R) and columns index servers (length numStations). Flips dispatching to
     * {@link RoutingStrategy#PROB}.
     */
    public Cluster setProbabilities(double[][] probs) {
        for (double[] row : probs) {
            if (row.length != numStations) {
                throw new IllegalArgumentException("each probs row must have length numStations");
            }
        }
        this.dispatchProbs = new double[probs.length][];
        for (int i = 0; i < probs.length; i++) this.dispatchProbs[i] = probs[i].clone();
        this.dispatching = RoutingStrategy.PROB;
        return this;
    }

    /**
     * Sets the squared coefficient of variation of the arrival process for every class
     * (broadcast). When SCV != 1 the arrival distribution becomes an APH fitted to the
     * given mean and SCV instead of an exponential.
     */
    public Cluster setArrivalSCV(double scv) {
        if (scv <= 0) throw new IllegalArgumentException("scv must be positive");
        int R = arrivalRates.length;
        this.arrivalScvs = new double[R];
        for (int r = 0; r < R; r++) this.arrivalScvs[r] = scv;
        return this;
    }

    /** Per-class arrival SCV. Length must equal the number of classes. */
    public Cluster setArrivalSCV(double[] scvs) {
        if (scvs.length != arrivalRates.length) {
            throw new IllegalArgumentException("scvs length must equal number of classes");
        }
        for (double s : scvs) if (s <= 0) throw new IllegalArgumentException("scv must be positive");
        this.arrivalScvs = scvs.clone();
        return this;
    }

    /**
     * Sets the squared coefficient of variation of the service distribution for every
     * (server, class) pair (broadcast). When SCV != 1 the service distribution becomes an
     * APH fitted to the supplied mean and SCV instead of an exponential.
     */
    public Cluster setServiceSCV(double scv) {
        if (scv <= 0) throw new IllegalArgumentException("scv must be positive");
        int R = closed ? (population != null ? population.length : 1) : arrivalRates.length;
        this.serviceScvs = new double[numStations][R];
        for (int i = 0; i < numStations; i++) {
            for (int r = 0; r < R; r++) this.serviceScvs[i][r] = scv;
        }
        return this;
    }

    /** Per-(server, class) service SCV. Outer dim must equal numStations. */
    public Cluster setServiceSCV(double[][] scvs) {
        if (scvs.length != numStations) {
            throw new IllegalArgumentException("serviceScvs outer dim must equal numStations");
        }
        for (double[] row : scvs) {
            for (double s : row) {
                if (s <= 0) throw new IllegalArgumentException("scv must be positive");
            }
        }
        this.serviceScvs = new double[scvs.length][];
        for (int i = 0; i < scvs.length; i++) this.serviceScvs[i] = scvs[i].clone();
        return this;
    }

    /**
     * Configures power-of-K-choices dispatching: pick the best of {@code k} randomly chosen
     * servers. Flips dispatching to {@link RoutingStrategy#SQ}.
     */
    public Cluster setSQ(int d) {
        if (d < 1) throw new IllegalArgumentException("d must be a positive integer");
        this.sqD = d;
        this.dispatching = RoutingStrategy.SQ;
        return this;
    }

    /**
     * Configures weighted round-robin dispatching with per-server integer weights. Flips
     * dispatching to {@link RoutingStrategy#WRROBIN}.
     */
    public Cluster setWeights(int[] weights) {
        if (weights.length != numStations) {
            throw new IllegalArgumentException("weights length must equal numStations");
        }
        this.dispatchWeights = weights.clone();
        this.dispatching = RoutingStrategy.WRROBIN;
        return this;
    }

    /** Sets the number of servers per station (default: all 1). */
    public Cluster setStationServers(int[] counts) {
        if (counts.length != numStations) {
            throw new IllegalArgumentException("counts length must equal numStations");
        }
        this.stationCounts = counts.clone();
        return this;
    }

    /**
     * Switches the cluster to the closed variant with a single class.
     *
     * @param population  number of jobs in the closed network
     * @param thinkTime   per-job think time at the delay
     */
    public Cluster setClosed(int population, double thinkTime) {
        return setClosed(new int[]{population}, new double[]{thinkTime});
    }

    /**
     * Switches the cluster to the closed variant with one entry per class.
     */
    public Cluster setClosed(int[] population, double[] thinkTimes) {
        if (population.length != thinkTimes.length) {
            throw new IllegalArgumentException("population and thinkTimes must have equal length");
        }
        this.closed = true;
        this.population = population.clone();
        this.thinkTimes = thinkTimes.clone();
        return this;
    }

    /** Builds the {@link Network} model from the current configuration. */
    public Network build() {
        int M = numStations;
        int R = closed ? population.length : arrivalRates.length;

        SchedStrategy[] strategies = new SchedStrategy[M];
        for (int i = 0; i < M; i++) strategies[i] = scheduling;

        Matrix S = new Matrix(M, 1);
        for (int i = 0; i < M; i++) S.set(i, 0, stationCounts[i]);

        Matrix D = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                double rate = serviceRates[i].length == 1 ? serviceRates[i][0] : serviceRates[i][r];
                if (rate <= 0) {
                    throw new IllegalStateException("Service rate must be positive");
                }
                D.set(i, r, 1.0 / rate);
            }
        }

        // Strategies that need per-destination or per-class parameters cannot be passed to
        // the static factory (which calls setRouting(jobclass, strategy) with no extras).
        // Build with RAND and apply the real strategy in post-processing.
        boolean needsPostProcess = (dispatching == RoutingStrategy.PROB && dispatchProbs != null)
                || (dispatching == RoutingStrategy.WRROBIN && dispatchWeights != null)
                || (dispatching == RoutingStrategy.SQ);
        RoutingStrategy factoryDispatch = needsPostProcess ? RoutingStrategy.RAND : dispatching;

        Network model;
        if (closed) {
            Matrix N = new Matrix(1, R);
            Matrix Z = new Matrix(1, R);
            for (int r = 0; r < R; r++) {
                N.set(0, r, population[r]);
                Z.set(0, r, thinkTimes[r]);
            }
            model = Network.clusterClosed(N, Z, D, strategies, S, factoryDispatch);
        } else {
            Matrix lambda = new Matrix(1, R);
            for (int r = 0; r < R; r++) lambda.set(0, r, arrivalRates[r]);
            model = Network.cluster(lambda, D, strategies, S, factoryDispatch);
        }

        if (dispatching == RoutingStrategy.PROB && dispatchProbs != null) {
            applyDispatcherProbabilities(model);
        } else if (dispatching == RoutingStrategy.WRROBIN && dispatchWeights != null) {
            applyDispatcherWeights(model);
        } else if (dispatching == RoutingStrategy.SQ) {
            applyDispatcherSQ(model);
        }

        applyDistributionScvs(model, R);
        return model;
    }

    private void applyDispatcherSQ(Network model) {
        Node dispatcher = model.getNodeByName("Dispatcher");
        List<JobClass> classes = model.getClasses();
        int k = sqD != null ? sqD : 2;
        for (JobClass cls : classes) {
            dispatcher.setSQRouting(cls, k);
        }
    }

    private void applyDistributionScvs(Network model, int R) {
        List<JobClass> classes = model.getClasses();

        // Arrival SCVs (open networks only).
        if (!closed && arrivalScvs != null) {
            Source src = (Source) model.getNodeByName("Source");
            for (int r = 0; r < R; r++) {
                double scv = arrivalScvs[r];
                if (scv != 1.0) {
                    src.setArrival(classes.get(r), APH.fitMeanAndSCV(1.0 / arrivalRates[r], scv));
                }
            }
        }

        // Service SCVs.
        if (serviceScvs != null) {
            for (int i = 0; i < numStations; i++) {
                Queue server = (Queue) model.getNodeByName("Station" + (i + 1));
                for (int r = 0; r < R; r++) {
                    double scv = serviceScvs[i][r];
                    if (scv != 1.0) {
                        double rate = serviceRates[i].length == 1 ? serviceRates[i][0]
                                : serviceRates[i][r];
                        server.setService(classes.get(r), APH.fitMeanAndSCV(1.0 / rate, scv));
                    }
                }
            }
        }
    }

    private void applyDispatcherProbabilities(Network model) {
        Node dispatcher = model.getNodeByName("Dispatcher");
        List<JobClass> classes = model.getClasses();
        for (int r = 0; r < classes.size(); r++) {
            double[] probs = dispatchProbs.length == 1 ? dispatchProbs[0] : dispatchProbs[r];
            for (int i = 0; i < numStations; i++) {
                Node server = model.getNodeByName("Station" + (i + 1));
                dispatcher.setProbRouting(classes.get(r), server, probs[i]);
            }
        }
    }

    private void applyDispatcherWeights(Network model) {
        Node dispatcher = model.getNodeByName("Dispatcher");
        List<JobClass> classes = model.getClasses();
        for (int r = 0; r < classes.size(); r++) {
            for (int i = 0; i < numStations; i++) {
                Node server = model.getNodeByName("Station" + (i + 1));
                dispatcher.setRouting(classes.get(r), RoutingStrategy.WRROBIN, server,
                        dispatchWeights[i]);
            }
        }
    }

    /**
     * Solves the cluster under each provided dispatching strategy.
     *
     * @param solverClass {@link NetworkSolver} subclass with a {@code Solver(Network)} constructor
     * @param policies    dispatching strategies to compare
     * @return ordered map from policy to its average-metric table
     */
    public Map<RoutingStrategy, NetworkAvgTable> compareDispatching(
            Class<? extends NetworkSolver> solverClass, RoutingStrategy... policies) {
        return compareDispatching(networkToAvgTable(solverClass), policies);
    }

    /**
     * Solves the cluster under each provided dispatching strategy using a custom solver factory.
     * Use this overload when the default {@code Solver(Network)} constructor is not adequate
     * (for example, to pass simulation options or use an alternative solver).
     */
    public Map<RoutingStrategy, NetworkAvgTable> compareDispatching(
            Function<Network, NetworkAvgTable> solverFactory, RoutingStrategy... policies) {
        Map<RoutingStrategy, NetworkAvgTable> out = new LinkedHashMap<>();
        RoutingStrategy saved = this.dispatching;
        try {
            for (RoutingStrategy p : policies) {
                this.dispatching = p;
                out.put(p, solverFactory.apply(build()));
            }
        } finally {
            this.dispatching = saved;
        }
        return out;
    }

    /**
     * Solves the cluster under each provided scheduling discipline.
     *
     * @param solverClass {@link NetworkSolver} subclass with a {@code Solver(Network)} constructor
     * @param disciplines scheduling disciplines to compare
     * @return ordered map from discipline to its average-metric table
     */
    public Map<SchedStrategy, NetworkAvgTable> compareScheduling(
            Class<? extends NetworkSolver> solverClass, SchedStrategy... disciplines) {
        return compareScheduling(networkToAvgTable(solverClass), disciplines);
    }

    /**
     * Solves the cluster under each provided scheduling discipline using a custom solver factory.
     */
    public Map<SchedStrategy, NetworkAvgTable> compareScheduling(
            Function<Network, NetworkAvgTable> solverFactory, SchedStrategy... disciplines) {
        Map<SchedStrategy, NetworkAvgTable> out = new LinkedHashMap<>();
        SchedStrategy saved = this.scheduling;
        try {
            for (SchedStrategy s : disciplines) {
                this.scheduling = s;
                out.put(s, solverFactory.apply(build()));
            }
        } finally {
            this.scheduling = saved;
        }
        return out;
    }

    /**
     * Sweeps the (single-class) arrival rate of an open cluster and solves at each value.
     *
     * @throws IllegalStateException if the cluster is closed or has more than one class
     */
    public Map<Double, NetworkAvgTable> sweepArrivalRate(double[] rates,
                                                         Class<? extends NetworkSolver> solverClass) {
        return sweepArrivalRate(rates, networkToAvgTable(solverClass));
    }

    /**
     * Same as {@link #sweepArrivalRate(double[], Class)} but with a custom solver factory.
     */
    public Map<Double, NetworkAvgTable> sweepArrivalRate(double[] rates,
                                                         Function<Network, NetworkAvgTable> solverFactory) {
        if (closed) {
            throw new IllegalStateException("sweepArrivalRate is only defined for open clusters");
        }
        if (arrivalRates.length != 1) {
            throw new IllegalStateException("sweepArrivalRate requires a single class");
        }
        Map<Double, NetworkAvgTable> out = new LinkedHashMap<>();
        double saved = arrivalRates[0];
        try {
            for (double r : rates) {
                arrivalRates[0] = r;
                out.put(r, solverFactory.apply(build()));
            }
        } finally {
            arrivalRates[0] = saved;
        }
        return out;
    }

    /**
     * Sweeps the number of parallel servers (homogeneous service rate replicated) and solves
     * at each value.
     */
    public Map<Integer, NetworkAvgTable> sweepNumStations(int[] counts,
                                                         Class<? extends NetworkSolver> solverClass) {
        return sweepNumStations(counts, networkToAvgTable(solverClass));
    }

    /**
     * Same as {@link #sweepNumStations(int[], Class)} but with a custom solver factory.
     */
    public Map<Integer, NetworkAvgTable> sweepNumStations(int[] counts,
                                                         Function<Network, NetworkAvgTable> solverFactory) {
        Map<Integer, NetworkAvgTable> out = new LinkedHashMap<>();
        int savedNumServers = this.numStations;
        double[][] savedServiceRates = this.serviceRates;
        int[] savedServerCounts = this.stationCounts;
        try {
            double[] perClassRate = savedServiceRates[0].clone();
            for (int m : counts) {
                this.numStations = m;
                this.serviceRates = new double[m][];
                for (int i = 0; i < m; i++) this.serviceRates[i] = perClassRate.clone();
                this.stationCounts = new int[m];
                for (int i = 0; i < m; i++) this.stationCounts[i] = 1;
                out.put(m, solverFactory.apply(build()));
            }
        } finally {
            this.numStations = savedNumServers;
            this.serviceRates = savedServiceRates;
            this.stationCounts = savedServerCounts;
        }
        return out;
    }

    private static Function<Network, NetworkAvgTable> networkToAvgTable(
            Class<? extends NetworkSolver> solverClass) {
        return model -> {
            try {
                Constructor<? extends NetworkSolver> ctor = solverClass.getConstructor(Network.class);
                NetworkSolver solver = ctor.newInstance(model);
                return solver.getAvgTable();
            } catch (ReflectiveOperationException e) {
                throw new RuntimeException("Failed to instantiate " + solverClass.getName()
                        + " with (Network) constructor", e);
            }
        };
    }

    private static double[][] fillUniform(int numStations, double rate) {
        double[][] out = new double[numStations][];
        for (int i = 0; i < numStations; i++) out[i] = new double[]{rate};
        return out;
    }

    private static double[][] deepCopy(double[][] src) {
        double[][] out = new double[src.length][];
        for (int i = 0; i < src.length; i++) out[i] = src[i].clone();
        return out;
    }
}
