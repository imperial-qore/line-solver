/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.GlobalConstants;
import jline.lang.ClosedClass;
import jline.lang.JobClass;
import jline.lang.OpenClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.*;
import jline.lang.processes.*;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

import java.io.Serializable;
import java.util.*;

/**
 * A node where jobs can spend time stationing there
 */
public abstract class Station extends ServiceNode implements Serializable {
    protected int numberOfServers;
    protected int cap;
    protected Map<JobClass, Integer> classCap;
    protected Map<JobClass, DropStrategy> dropRule;
    protected Matrix lldScaling;
    // Class-dependence function beta_i(n) of the per-class population vector at
    // this station. Returns either a 1x1 matrix (a chain-independent scaling
    // shared by every class) or a length-R row vector of the per-class scalings
    // [beta_{i,1}(n), ..., beta_{i,R}(n)], the latter expressing Sauer's
    // chain-dependent service rates mu_{r,i}(n) (Sauer 1983, eq. (40)).
    protected SerializableFunction<Matrix, Matrix> lcdScaling;
    // Declared peak (maximum) class-dependent rate scaling per class, a 1xR row
    // vector (scalar 1x1 broadcast to all classes). Used to normalize
    // utilization as U = T*S/peak (the T*S/c convention of multiserver stations).
    protected Matrix lcdScalingPeak;
    // Joint-dependence function eta_i(n) of the per-class population vector at
    // this station. Like lcdScaling it returns either a 1x1 matrix (a scaling
    // shared by every class) or a length-R row vector, but UNLIKE lcdScaling
    // (product-form beta_{i,r} depending on the own-class marginal n_{i,r}),
    // eta may read the joint vector arbitrarily and is therefore
    // NON-product-form. Evaluated numerically identically to lcdScaling; the
    // separate field carries the product-form vs joint distinction.
    protected SerializableFunction<Matrix, Matrix> ljdScaling;
    // Declared peak joint-dependent rate scaling per class, twin of
    // lcdScalingPeak. Used to normalize utilization as U = T*S/peak.
    protected Matrix ljdScalingPeak;
    protected Map<JobClass, Map<JobClass, Distribution>> switchoverTimes;
    protected Map<JobClass, Distribution> patienceDistributions;
    protected Map<JobClass, ImpatienceType> impatienceTypes;

    // Balking configuration (station-specific overrides)
    protected Map<JobClass, BalkingStrategy> balkingStrategies;
    protected Map<JobClass, List<BalkingThreshold>> balkingThresholds;

    // Retrial configuration (station-specific overrides)
    protected Map<JobClass, Distribution> retrialDelayDistributions;
    protected Map<JobClass, Integer> retrialMaxAttempts;

    // Orbit shape declared by setOrbit: retrial policy and orbit capacity. Both
    // default to the classical unbounded per-customer orbit, which is what a plain
    // setRetrial declares, so an absent entry is LINEAR / -1.
    protected Map<JobClass, Integer> retrialPolicies;
    protected Map<JobClass, Integer> orbitMaxJobs;

    // Orbit impatience configuration (abandonment from the retrial orbit)
    protected Map<JobClass, Distribution> orbitImpatienceDistributions;
    protected Map<JobClass, Double> batchRejectProb;

    /**
     * Creates a new station with the specified name.
     * Initializes default capacity, scheduling, and service configurations.
     *
     * @param name the name for this station
     */
    public Station(String name) {
        super(name);

        this.classCap = new HashMap<JobClass, Integer>();

        this.cap = Integer.MAX_VALUE;
        this.dropRule = new HashMap<JobClass, DropStrategy>();

        this.lldScaling = new Matrix(0, 0);
        this.lcdScaling = null;
        this.ljdScaling = null;

        this.schedStrategyPar = new HashMap<JobClass, Double>();
        this.switchoverTimes = new HashMap<JobClass, Map<JobClass, Distribution>>();
        this.patienceDistributions = new HashMap<JobClass, Distribution>();
        this.impatienceTypes = new HashMap<JobClass, ImpatienceType>();

        // Initialize balking and retrial maps
        this.balkingStrategies = new HashMap<JobClass, BalkingStrategy>();
        this.balkingThresholds = new HashMap<JobClass, List<BalkingThreshold>>();
        this.retrialDelayDistributions = new HashMap<JobClass, Distribution>();
        this.retrialMaxAttempts = new HashMap<JobClass, Integer>();
        this.retrialPolicies = new HashMap<JobClass, Integer>();
        this.orbitMaxJobs = new HashMap<JobClass, Integer>();
        this.orbitImpatienceDistributions = new HashMap<JobClass, Distribution>();
        this.batchRejectProb = new HashMap<JobClass, Double>();
    }

    @Override
    public double getCap() {
        return cap;
    }

    /**
     * Whether this station has a real, bounded buffer.
     * <p>
     * The raw {@link #getCap()} value encodes "unbounded" in THREE different ways, so a
     * bare comparison against it is a recurring source of defects:
     * <ul>
     * <li>{@code Integer.MAX_VALUE} -- the default for a station the user never capped;</li>
     * <li>{@code Double.POSITIVE_INFINITY} -- what an explicit {@code setCapacity(Inf)} stores;</li>
     * <li>a NEGATIVE value, typically {@code -1} -- the sentinel {@code JMT2LINE} copies
     * straight from JMT's {@code size="-1"} ("unlimited") when importing a .jsimg. A
     * negative cap means UNLIMITED, never "a capacity of -1".</li>
     * </ul>
     * Testing only {@code getCap() < Integer.MAX_VALUE} therefore reports a JMT-imported
     * station as finite-capacity: that defect rejected every imported model through
     * {@code NetworkSolver.bindingCapacityReason} (23 M2MTest failures), and the same
     * pattern in {@code Ensemble} copied the sentinel into a clone. Use this accessor
     * instead of re-deriving the encoding at each call site.
     *
     * @return true iff the station has a finite, non-negative buffer bound
     */
    public boolean hasFiniteCap() {
        return cap >= 0 && cap < Integer.MAX_VALUE && !Double.isInfinite(cap);
    }

    /**
     * Sets the total capacity limit for this station.
     * Capacity follows Kendall notation where K = total system capacity (queue + in-service).
     *
     * @param cap the maximum number of jobs that can be at this station
     */
    public void setCapacity(int cap) {
        this.cap = cap;
    }

    /**
     * Alias for setCapacity() for backwards compatibility.
     *
     * @param cap the maximum number of jobs that can be at this station
     */
    public void setCap(int cap) {
        this.setCapacity(cap);
    }

    @Override
    public double getClassCap(JobClass jobClass) {
        if (classCap.containsKey(jobClass)) {
            return FastMath.min(this.classCap.get(jobClass), cap);
        }

        return cap;
    }

    /**
     * Gets the drop strategy for a specific job class when capacity is exceeded.
     * 
     * @param jobclass the job class to query
     * @return the drop strategy for the job class, or null if not set
     */
    public DropStrategy getDropRule(JobClass jobclass) {
        return this.dropRule.getOrDefault(jobclass, null);
    }

    /**
     * Gets the limited class-dependent scaling function for this station.
     * 
     * @return the class-dependent scaling function, or null if not set
     */
    public SerializableFunction<Matrix, Matrix> getLimitedClassDependence() {
        return this.lcdScaling;
    }

    /**
     * Sets the limited class-dependent scaling function for this station.
     * 
     * @param gamma the class-dependent scaling function
     */
    public void setLimitedClassDependence(SerializableFunction<Matrix, Matrix> gamma) {
        this.lcdScaling = gamma;
    }

    /**
     * Sets the limited class-dependent scaling function and its declared peak
     * rate scaling per class (used to normalize utilization as U = T*S/peak).
     *
     * @param gamma the class-dependent scaling function
     * @param peakRatePerClass 1x1 (broadcast) or 1xR peak rate scaling per class
     */
    public void setLimitedClassDependence(SerializableFunction<Matrix, Matrix> gamma, Matrix peakRatePerClass) {
        this.lcdScaling = gamma;
        this.lcdScalingPeak = peakRatePerClass;
    }

    /**
     * Gets the declared peak class-dependent rate scaling per class.
     *
     * @return the 1x1 or 1xR peak rate scaling, or null if not set
     */
    public Matrix getLimitedClassDependencePeak() {
        return this.lcdScalingPeak;
    }

    /**
     * Gets the limited joint-dependent (non-product-form) scaling function.
     *
     * @return the joint-dependent scaling function eta_i, or null if not set
     */
    public SerializableFunction<Matrix, Matrix> getLimitedJointDependence() {
        return this.ljdScaling;
    }

    /**
     * Sets the limited joint-dependent scaling function for this station.
     *
     * @param eta the joint-dependent scaling function
     */
    public void setLimitedJointDependence(SerializableFunction<Matrix, Matrix> eta) {
        this.ljdScaling = eta;
    }

    /**
     * Sets the limited joint-dependent scaling function and its declared peak
     * rate scaling per class (used to normalize utilization as U = T*S/peak).
     * Unlike setLimitedClassDependence (product-form beta_{i,r}), eta may read
     * the joint per-class population vector arbitrarily and is non-product-form.
     *
     * @param eta the joint-dependent scaling function
     * @param peakRatePerClass 1x1 (broadcast) or 1xR peak rate scaling per class
     */
    public void setLimitedJointDependence(SerializableFunction<Matrix, Matrix> eta, Matrix peakRatePerClass) {
        this.ljdScaling = eta;
        this.ljdScalingPeak = peakRatePerClass;
    }

    /**
     * Gets the declared peak joint-dependent rate scaling per class.
     *
     * @return the 1x1 or 1xR peak rate scaling, or null if not set
     */
    public Matrix getLimitedJointDependencePeak() {
        return this.ljdScalingPeak;
    }

    /**
     * Gets the limited load-dependent scaling matrix for this station.
     * 
     * @return the load-dependent scaling matrix
     */
    public Matrix getLimitedLoadDependence() {
        return this.lldScaling;
    }

    /**
     * Sets the limited load-dependent scaling matrix for this station.
     *
     * @param alpha the load-dependent scaling matrix
     */
    public void setLimitedLoadDependence(Matrix alpha) {
        this.lldScaling = alpha;
    }

    @Override
    public int getNumberOfServers() {
        return this.numberOfServers;
    }

    /**
     * Sets the number of servers at this station.
     * 
     * @param numberOfServers the number of servers (use Integer.MAX_VALUE for infinite servers)
     */
    public void setNumberOfServers(int numberOfServers) {
        if (numberOfServers == -1) { // may result of a double Inf is casted as input
            numberOfServers = Integer.MAX_VALUE;
        }
        this.numberOfServers = numberOfServers;
    }

    /**
     * Sets the scheduling-strategy parameter (e.g. the class weight used by
     * DPS/GPS scheduling) for the given job class.
     *
     * @param jobClass the job class
     * @param weight   the scheduling parameter (weight) for the class
     */
    public void setStrategyParam(JobClass jobClass, double weight) {
        this.schedStrategyPar.put(jobClass, weight);
    }

    @Override
    public void removeJobClass(JobClass jobClass) {
        super.removeJobClass(jobClass);
        this.classCap.remove(jobClass);
        this.schedStrategyPar.remove(jobClass);
        this.dropRule.remove(jobClass);
    }

    /**
     * Returns the scheduling-strategy parameter for the given job class
     * (0.0 if not set).
     *
     * @param jobClass the job class
     * @return the scheduling parameter (weight) for the class
     */
    public double getStrategyParam(JobClass jobClass) {
        return this.schedStrategyPar.getOrDefault(jobClass, 0.0);
    }

    /**
     * Returns the scheduling strategy used by this station.
     *
     * @return the scheduling strategy (FCFS, PS, etc.)
     */
    public SchedStrategy getSchedStrategy() {
        return this.schedStrategy;
    }

    /**
     * Gets the service rates configured for all job classes at this station.
     * Returns a list containing rate values or Distribution objects.
     * 
     * @return list of service rates for each job class
     */
    /**
     * The {mean, SCV} pair that non-Markovian discrete distributions store in
     * {@code sn.proc}. Bernoulli, Binomial, Poisson and Geometric each declare
     * their own {@code getProcess()} returning exactly this, but the method is
     * not on {@code DiscreteDistribution}, so the pair is rebuilt here rather
     * than casting to each concrete type.
     *
     * @param distr the distribution
     * @return a two-entry cell holding the mean and the SCV
     */
    private static MatrixCell momentPair(Distribution distr) {
        MatrixCell cell = new MatrixCell();
        cell.set(0, Matrix.singleton(distr.getMean()));
        cell.set(1, Matrix.singleton(distr.getSCV()));
        return cell;
    }

    public List<Object> getServiceRates() {
        int nClasses = this.model.getNumberOfClasses();
        Map<JobClass, MatrixCell> map = new HashMap<JobClass, MatrixCell>();
        Map<JobClass, Matrix> mu = new HashMap<JobClass, Matrix>();
        Map<JobClass, Matrix> phi = new HashMap<JobClass, Matrix>();
        for (int i = 0; i < nClasses; i++) {
            JobClass jobclass = this.model.getClassByIndex(i);
            if (this instanceof Queue) {
                Queue queue = (Queue) this;
                // Since Delay, Join are all sub-class of Queue, and only queue and source
                // are the sub-class of station, we could cast this to queue to call setService method
                if (!queue.containsJobClass(jobclass)) {
                    queue.setService(jobclass, new Disabled());
                    Matrix nan_matrix = new Matrix(1, 1, 1);
                    nan_matrix.fill(Double.NaN);
                    MatrixCell tmp = new MatrixCell();
                    tmp.set(0, nan_matrix);
                    tmp.set(1, nan_matrix.copy());
                    map.put(jobclass, tmp);
                    mu.put(jobclass, nan_matrix.copy());
                    phi.put(jobclass, nan_matrix.copy());
                } else if (queue.getServiceProcess(jobclass) instanceof Immediate) {
                    Distribution distr = this.server.getServiceDistribution(jobclass);
                    Matrix map_matrix_1 = new Matrix(1, 1, 1);
                    map_matrix_1.set(0, 0, -GlobalConstants.Immediate);
                    Matrix map_matrix_2 = new Matrix(1, 1, 1);
                    map_matrix_2.set(0, 0, GlobalConstants.Immediate);
                    Matrix mu_matrix = new Matrix(1, 1, 1);
                    mu_matrix.set(0, 0, GlobalConstants.Immediate);
                    Matrix phi_matrix = new Matrix(1, 1, 1);
                    phi_matrix.set(0, 0, 1);
                    MatrixCell tmp = new MatrixCell();
                    tmp.set(0, map_matrix_1);
                    tmp.set(1, map_matrix_2);
                    map.put(jobclass, tmp);
                    mu.put(jobclass, mu_matrix);
                    phi.put(jobclass, phi_matrix);
                } else if (queue.getServiceProcess(jobclass) instanceof Det) {
                    Det distr = (Det) this.server.getServiceDistribution(jobclass);
                    map.put(jobclass, distr.getRepresentation());
                    mu.put(jobclass, Matrix.singleton(distr.getRate()));
                    phi.put(jobclass, Matrix.singleton(1.0));
                } else if ((queue.getServiceProcess(jobclass) instanceof Uniform)
                        || (queue.getServiceProcess(jobclass) instanceof Lognormal)
                        || (queue.getServiceProcess(jobclass) instanceof Weibull)
                        || (queue.getServiceProcess(jobclass) instanceof Gamma)
                        || (queue.getServiceProcess(jobclass) instanceof Pareto)) {
                    ContinuousDistribution distr = (ContinuousDistribution) this.server.getServiceDistribution(jobclass);
                    map.put(jobclass, distr.getProcess());
                    mu.put(jobclass, Matrix.singleton(distr.getRate()));
                    phi.put(jobclass, Matrix.singleton(1.0));
                } else if (queue.getServiceProcess(jobclass) instanceof Bernoulli
                        || queue.getServiceProcess(jobclass) instanceof Binomial
                        || queue.getServiceProcess(jobclass) instanceof Poisson) {
                    // Counting distributions used as a service time; see the
                    // matching branch in the source arrival path above.
                    Distribution distr = this.server.getServiceDistribution(jobclass);
                    map.put(jobclass, momentPair(distr));
                    mu.put(jobclass, Matrix.singleton(distr.getRate()));
                    phi.put(jobclass, Matrix.singleton(1.0));
                } else if (queue.getServiceProcess(jobclass) instanceof Geometric) {
                    // Lattice-valued service time with support {1,2,...}; see the
                    // matching branch in the source arrival path above.
                    Geometric distr = (Geometric) this.server.getServiceDistribution(jobclass);
                    map.put(jobclass, distr.getProcess());
                    mu.put(jobclass, Matrix.singleton(distr.getRate()));
                    phi.put(jobclass, Matrix.singleton(1.0));
                } else if ((queue.getServiceProcess(jobclass) instanceof Replayer) || (queue.getServiceProcess(jobclass) instanceof Trace)) {
                    Replayer distr = (Replayer) this.server.getServiceDistribution(jobclass);
                    Markovian aph = distr.fitAPH();
                    map.put(jobclass, aph.getProcess());
                    mu.put(jobclass, aph.getMu());
                    phi.put(jobclass, aph.getPhi());
                } else if (queue.getServiceProcess(jobclass) instanceof NHPP) {
                    // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
                    NHPP distr = (NHPP) this.server.getServiceDistribution(jobclass);
                    map.put(jobclass, distr.getProcess());
                    mu.put(jobclass, Matrix.singleton(distr.getTimeAverageRate()));
                    phi.put(jobclass, Matrix.singleton(1.0));
                } else if (queue.getServiceProcess(jobclass) instanceof Prior) {
                    // Prior distributions are handled by the UQ solver;
                    // use the prior-weighted mean rate as a placeholder for getStruct()
                    Prior priorDist = (Prior) this.server.getServiceDistribution(jobclass);
                    double rate = 1.0 / priorDist.getMean();
                    Matrix map_matrix_1 = new Matrix(1, 1, 1);
                    map_matrix_1.set(0, 0, -rate);
                    Matrix map_matrix_2 = new Matrix(1, 1, 1);
                    map_matrix_2.set(0, 0, rate);
                    MatrixCell tmp = new MatrixCell();
                    tmp.set(0, map_matrix_1);
                    tmp.set(1, map_matrix_2);
                    map.put(jobclass, tmp);
                    mu.put(jobclass, Matrix.singleton(rate));
                    phi.put(jobclass, Matrix.singleton(1.0));
                } else if ((queue.getServiceProcess(jobclass) instanceof ME)
                        || (queue.getServiceProcess(jobclass) instanceof RAP)) {
                    putMatrixExponential(map, mu, phi, jobclass,
                            (Markovian) this.server.getServiceDistribution(jobclass));
                } else if (!(queue.getServiceProcess(jobclass) instanceof Disabled)) {
                    Distribution distr = this.server.getServiceDistribution(jobclass);
                    if (distr instanceof Markovian) {
                        map.put(jobclass, ((Markovian) distr).getProcess());
                        mu.put(jobclass, ((Markovian) distr).getMu());
                        phi.put(jobclass, ((Markovian) distr).getPhi());
                    } else {
                        // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
                        Matrix nan_matrix = new Matrix(1, 1, 1);
                        nan_matrix.fill(Double.NaN);
                        MatrixCell tmp = new MatrixCell();
                        tmp.set(0, nan_matrix);
                        tmp.set(1, nan_matrix.copy());
                        map.put(jobclass, tmp);
                        mu.put(jobclass, nan_matrix.copy());
                        phi.put(jobclass, nan_matrix.copy());
                    }
                } else { // Disabled
                    Matrix nan_matrix = new Matrix(1, 1, 1);
                    nan_matrix.fill(Double.NaN);
                    MatrixCell tmp = new MatrixCell();
                    tmp.set(0, nan_matrix);
                    tmp.set(1, nan_matrix.copy());
                    map.put(jobclass, tmp);
                    mu.put(jobclass, nan_matrix.copy());
                    phi.put(jobclass, nan_matrix.copy());
                }
            } else if (this instanceof Place && ((Place) this).isQueueing()) {
                // Queueing place (QPN): the embedded queue has a real Server section, so its
                // service process is represented exactly as for a Queue station.
                Distribution svc = ((Place) this).getService(jobclass);
                if (svc == null || svc instanceof Disabled) {
                    Matrix nan_matrix = new Matrix(1, 1, 1);
                    nan_matrix.fill(Double.NaN);
                    MatrixCell tmp = new MatrixCell();
                    tmp.set(0, nan_matrix);
                    tmp.set(1, nan_matrix.copy());
                    map.put(jobclass, tmp);
                    mu.put(jobclass, nan_matrix.copy());
                    phi.put(jobclass, nan_matrix.copy());
                } else if (svc instanceof Immediate) {
                    Matrix map_matrix_1 = new Matrix(1, 1, 1);
                    map_matrix_1.set(0, 0, -GlobalConstants.Immediate);
                    Matrix map_matrix_2 = new Matrix(1, 1, 1);
                    map_matrix_2.set(0, 0, GlobalConstants.Immediate);
                    MatrixCell tmp = new MatrixCell();
                    tmp.set(0, map_matrix_1);
                    tmp.set(1, map_matrix_2);
                    map.put(jobclass, tmp);
                    mu.put(jobclass, Matrix.singleton(GlobalConstants.Immediate));
                    phi.put(jobclass, Matrix.singleton(1.0));
                } else if (svc instanceof Det) {
                    Det distr = (Det) this.server.getServiceDistribution(jobclass);
                    map.put(jobclass, distr.getRepresentation());
                    mu.put(jobclass, Matrix.singleton(distr.getRate()));
                    phi.put(jobclass, Matrix.singleton(1.0));
                } else if ((svc instanceof Uniform) || (svc instanceof Lognormal)
                        || (svc instanceof Weibull) || (svc instanceof Gamma) || (svc instanceof Pareto)) {
                    ContinuousDistribution distr = (ContinuousDistribution) this.server.getServiceDistribution(jobclass);
                    map.put(jobclass, distr.getProcess());
                    mu.put(jobclass, Matrix.singleton(distr.getRate()));
                    phi.put(jobclass, Matrix.singleton(1.0));
                } else if ((svc instanceof Replayer) || (svc instanceof Trace)) {
                    Replayer distr = (Replayer) this.server.getServiceDistribution(jobclass);
                    Markovian aph = distr.fitAPH();
                    map.put(jobclass, aph.getProcess());
                    mu.put(jobclass, aph.getMu());
                    phi.put(jobclass, aph.getPhi());
                } else if ((svc instanceof ME) || (svc instanceof RAP)) {
                    putMatrixExponential(map, mu, phi, jobclass,
                            (Markovian) this.server.getServiceDistribution(jobclass));
                } else if (svc instanceof Markovian) {
                    Markovian distr = (Markovian) this.server.getServiceDistribution(jobclass);
                    map.put(jobclass, distr.getProcess());
                    mu.put(jobclass, distr.getMu());
                    phi.put(jobclass, distr.getPhi());
                } else {
                    Matrix nan_matrix = new Matrix(1, 1, 1);
                    nan_matrix.fill(Double.NaN);
                    MatrixCell tmp = new MatrixCell();
                    tmp.set(0, nan_matrix);
                    tmp.set(1, nan_matrix.copy());
                    map.put(jobclass, tmp);
                    mu.put(jobclass, nan_matrix.copy());
                    phi.put(jobclass, nan_matrix.copy());
                }
            } else { // Disabled (ordinary Place)
                Matrix nan_matrix = new Matrix(1, 1, 1);
                nan_matrix.fill(Double.NaN);
                MatrixCell tmp = new MatrixCell();
                tmp.set(0, nan_matrix);
                tmp.set(1, nan_matrix.copy());
                map.put(jobclass, tmp);
                mu.put(jobclass, nan_matrix.copy());
                phi.put(jobclass, nan_matrix.copy());
            }
        }
        return new ArrayList<Object>(Arrays.asList(map, mu, phi));
    }

    /**
     * Records an ME or RAP process in the struct without a phase-type decomposition.
     *
     * <p>The (D0, D1) matrices are stored exactly as for any other Markovian
     * process, because they are the complete and exact description of an ME or
     * RAP and are what LDES (the only solver that accepts these processes) reads.
     * The phase quantities mu and phi are recorded as NaN instead, which is the
     * sentinel this file already uses for "no phase-type representation exists"
     * (see the non-Markovian arm of {@link #getServiceRates()}). Computing them
     * would mean calling {@code Markovian.getMu}/{@code getPhi}, whose
     * probabilistic reading of (D0, D1) fails for these processes and yields
     * entries outside [0,1]; ME and RAP now refuse those calls outright. Storing
     * NaN makes any solver that wrongly consumes phi as a probability fail
     * loudly rather than silently propagate negative rates.</p>
     *
     * @param map      per-class process representations to populate
     * @param mu       per-class phase rates to populate (set to NaN)
     * @param phi      per-class exit probabilities to populate (set to NaN)
     * @param jobclass the class whose entries are being written
     * @param distr    the ME or RAP process
     */
    private static void putMatrixExponential(Map<JobClass, MatrixCell> map,
                                             Map<JobClass, Matrix> mu,
                                             Map<JobClass, Matrix> phi,
                                             JobClass jobclass,
                                             Markovian distr) {
        // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
        map.put(jobclass, distr.getProcess());
        mu.put(jobclass, distr.getMu());
        phi.put(jobclass, distr.getPhi());
    }

    /**
     * Gets the source (arrival) rates for all job classes at this station.
     * Only applicable for source nodes.
     *
     * @return list of source rates for each job class
     */
    public List<Object> getSourceRates() {
        int nClasses = this.model.getNumberOfClasses();
        Map<JobClass, MatrixCell> map = new HashMap<JobClass, MatrixCell>();
        Map<JobClass, Matrix> mu = new HashMap<JobClass, Matrix>();
        Map<JobClass, Matrix> phi = new HashMap<JobClass, Matrix>();
        for (int r = 0; r < nClasses; r++) {
            Source source = (Source) this;
            JobClass jobclass = this.model.getClassByIndex(r);
            if (!source.containsJobClass(jobclass)) {
                source.setArrival(jobclass, new Disabled());
                // Open classes without arrivals should still route (for class switching pass-through)
                source.setRouting(jobclass, (jobclass instanceof OpenClass) ? RoutingStrategy.RAND : RoutingStrategy.DISABLED);
                Matrix nan_matrix = new Matrix(1, 1, 1);
                nan_matrix.fill(Double.NaN);
                MatrixCell tmp = new MatrixCell();
                tmp.set(0, nan_matrix);
                tmp.set(1, nan_matrix.copy());
                map.put(jobclass, tmp);
                mu.put(jobclass, nan_matrix.copy());
                phi.put(jobclass, nan_matrix.copy());
            } else if (source.getArrivalDistribution(jobclass) instanceof Det) {
                Det distr = (Det) source.getArrivalDistribution(jobclass);
                map.put(jobclass, distr.getRepresentation());
                mu.put(jobclass, Matrix.singleton(distr.getRate()));
                phi.put(jobclass, Matrix.singleton(1.0));
            } else if ((source.getArrivalDistribution(jobclass) instanceof Uniform) || (source.getArrivalDistribution(jobclass) instanceof Pareto) || (source.getArrivalDistribution(jobclass) instanceof Lognormal) || (source.getArrivalDistribution(jobclass) instanceof Gamma) || (source.getArrivalDistribution(jobclass) instanceof Weibull)) {
                ContinuousDistribution distr = (ContinuousDistribution) source.getArrivalDistribution(jobclass);
                map.put(jobclass, distr.getProcess());
                mu.put(jobclass, Matrix.singleton(distr.getRate()));
                phi.put(jobclass, Matrix.singleton(1.0));
            } else if (source.getArrivalDistribution(jobclass) instanceof Bernoulli
                    || source.getArrivalDistribution(jobclass) instanceof Binomial
                    || source.getArrivalDistribution(jobclass) instanceof Poisson) {
                // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
                Distribution distr = source.getArrivalDistribution(jobclass);
                map.put(jobclass, momentPair(distr));
                mu.put(jobclass, Matrix.singleton(distr.getRate()));
                phi.put(jobclass, Matrix.singleton(1.0));
            } else if (source.getArrivalDistribution(jobclass) instanceof Geometric) {
                // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
                Geometric distr = (Geometric) source.getArrivalDistribution(jobclass);
                map.put(jobclass, distr.getProcess());
                mu.put(jobclass, Matrix.singleton(distr.getRate()));
                phi.put(jobclass, Matrix.singleton(1.0));
            } else if ((source.getArrivalDistribution(jobclass) instanceof Replayer) || (source.getArrivalDistribution(jobclass) instanceof Trace)) {
                Replayer distr = (Replayer) source.getArrivalDistribution(jobclass);
                Markovian aph = distr.fitAPH();
                map.put(jobclass, aph.getProcess());
                mu.put(jobclass, aph.getMu());
                phi.put(jobclass, aph.getPhi());
            } else if (source.getArrivalDistribution(jobclass) instanceof NHPP) {
                NHPP distr = (NHPP) source.getArrivalDistribution(jobclass);
                map.put(jobclass, distr.getProcess());
                mu.put(jobclass, Matrix.singleton(distr.getTimeAverageRate()));
                phi.put(jobclass, Matrix.singleton(1.0));
            } else if ((source.getArrivalDistribution(jobclass) instanceof ME)
                    || (source.getArrivalDistribution(jobclass) instanceof RAP)) {
                putMatrixExponential(map, mu, phi, jobclass,
                        (Markovian) source.getArrivalDistribution(jobclass));
            } else if (!(source.getArrivalDistribution(jobclass) instanceof Disabled) && (source.getArrivalDistribution(jobclass) instanceof Markovian)) {
                Distribution distr = source.getArrivalDistribution(jobclass);
                map.put(jobclass, ((Markovian) distr).getProcess());
                mu.put(jobclass, ((Markovian) distr).getMu());
                phi.put(jobclass, ((Markovian) distr).getPhi());
            } else {
                Matrix nan_matrix = new Matrix(1, 1, 1);
                nan_matrix.fill(Double.NaN);
                MatrixCell tmp = new MatrixCell();
                tmp.set(0, nan_matrix);
                tmp.set(1, nan_matrix.copy());
                map.put(jobclass, tmp);
                mu.put(jobclass, nan_matrix.copy());
                phi.put(jobclass, nan_matrix.copy());
            }

        }
        return new ArrayList<Object>(Arrays.asList(map, mu, phi));
    }

    @Override
    public boolean isReferenceStation() {
        for (JobClass jobClass : this.model.getClasses()) {
            if (jobClass instanceof ClosedClass) {
                if (jobClass.getReferenceStation() == this) {
                    return true;
                }
            }
        }

        return false;
    }

    /**
     * Checks which job classes have defined service processes.
     *
     * @return array indicating which classes have defined service processes
     */
    public boolean[] isServiceDefined() {
        List<JobClass> jobClasses = this.model.getClasses();
        int nClasses = jobClasses.size();
        boolean[] isD = new boolean[nClasses];

        for (int r = 0; r < nClasses; r++) {
            JobClass jobClass = jobClasses.get(r);
            isD[r] = isServiceDefined(jobClass);
        }

        return isD;
    }

    /**
     * Checks if a specific job class has a defined service process.
     *
     * @param jobClass the job class to check
     * @return true if service is defined for this class, false otherwise
     */
    public boolean isServiceDefined(JobClass jobClass) {
        if (this.server == null) {
            return false;
        }

        // Check if the service process exists and is not null
        return this.server.containsJobClass(jobClass) &&
                this.server.getServiceProcess(jobClass) != null;
    }

    /**
     * Checks which job classes have disabled service processes.
     *
     * @return array indicating which classes have disabled service processes
     */
    public boolean[] isServiceDisabled() {
        List<JobClass> jobClasses = this.model.getClasses();
        int nClasses = jobClasses.size();
        boolean[] isD = new boolean[nClasses];

        for (int r = 0; r < nClasses; r++) {
            JobClass jobClass = jobClasses.get(r);
            isD[r] = isServiceDisabled(jobClass);
        }

        return isD;
    }

    /**
     * Checks if a specific job class has a disabled service process.
     *
     * @param jobClass the job class to check
     * @return true if service is disabled for this class, false otherwise
     */
    public boolean isServiceDisabled(JobClass jobClass) {
        if (this.server == null) {
            return true;
        }

        // If no service process is defined, consider it disabled
        if (!this.server.containsJobClass(jobClass)) {
            return true;
        }

        // Check if the service distribution is disabled
        return this.server.getServiceDistribution(jobClass).isDisabled();
    }

    /**
     * Sets the capacity for each chain in the network.
     * Configures per-class capacity based on chain membership and service availability.
     */
    public void setChainCapacity(double[] values) {
        if (this.model == null) {
            throw new RuntimeException("Network model not set");
        }

        // Get number of chains from the model
        int nchains = this.model.getNumberOfChains();
        if (values.length != nchains) {
            throw new IllegalArgumentException("The method requires a capacity value for each chain. Expected " + nchains + " values, got " + values.length);
        }

        // Set capacity for each class based on chain membership
        NetworkStruct sn = this.model.getStruct();
        for (int c = 0; c < nchains; c++) {
            Matrix inchain = sn.inchain.get(c);
            if (inchain != null) {
                for (int r = 0; r < inchain.getNumCols(); r++) {
                    if (inchain.get(0, r) == 1.0) { // Class r is in chain c
                        JobClass jobClass = sn.jobclasses.get(r);
                        if (!isServiceDisabled(jobClass)) {
                            setClassCap(jobClass, (int) values[c]);
                        } else {
                            setClassCap(jobClass, Integer.MAX_VALUE);
                        }
                    }
                }
            }
        }

        // Update station capacity to minimum of sum of class capacities and current capacity
        int sumClassCap = 0;
        for (Integer classCap : this.classCap.values()) {
            if (classCap > 0) {
                sumClassCap += classCap;
            }
        }
        this.cap = Math.min(sumClassCap, this.cap);
    }

    /**
     * Sets the capacity limit for a specific job class at this station.
     * 
     * @param jobClass the job class to set capacity for
     * @param cap the maximum number of jobs of this class
     */
    public void setClassCap(JobClass jobClass, int cap) {
        this.classCap.put(jobClass, cap);
    }

    /**
     * Sets the drop strategy for a specific job class when capacity is exceeded.
     * 
     * @param jobclass the job class to configure
     * @param drop the drop strategy to apply
     */
    public void setDropRule(JobClass jobclass, DropStrategy drop) {
        this.dropRule.put(jobclass, drop);
    }

    /**
     * Sets the switchover time from one job class to another.
     * This is used for class-to-class switchover functionality.
     *
     * @param fromClass the job class to switch from
     * @param toClass the job class to switch to
     * @param switchoverTime the switchover time distribution
     */
    public void setSwitchoverTime(JobClass fromClass, JobClass toClass, Distribution switchoverTime) {
        if (!this.switchoverTimes.containsKey(fromClass)) {
            this.switchoverTimes.put(fromClass, new HashMap<JobClass, Distribution>());
        }
        this.switchoverTimes.get(fromClass).put(toClass, switchoverTime);
    }

    /**
     * Gets the switchover time from one job class to another.
     *
     * @param fromClass the job class to switch from
     * @param toClass the job class to switch to
     * @return the switchover time distribution, or null if not set
     */
    public Distribution getSwitchoverTime(JobClass fromClass, JobClass toClass) {
        if (!this.switchoverTimes.containsKey(fromClass)) {
            return null;
        }
        return this.switchoverTimes.get(fromClass).get(toClass);
    }

    /**
     * Checks if a switchover time is set for the given job class pair.
     *
     * @param fromClass the job class to switch from
     * @param toClass the job class to switch to
     * @return true if switchover time is set, false otherwise
     */
    public boolean hasSwitchoverTime(JobClass fromClass, JobClass toClass) {
        return this.switchoverTimes.containsKey(fromClass) &&
               this.switchoverTimes.get(fromClass).containsKey(toClass);
    }

    /**
     * Sets the patience distribution for a specific job class at this station.
     * Jobs that wait longer than their patience time will abandon the queue.
     * This setting takes precedence over the global class patience.
     * Defaults to RENEGING patience type for backwards compatibility.
     *
     * @param jobClass the job class
     * @param distribution the patience time distribution
     * @throws IllegalArgumentException if distribution is a modulated process or arguments are null
     */
    public void setPatience(JobClass jobClass, Distribution distribution) {
        setPatience(jobClass, ImpatienceType.RENEGING, distribution);
    }

    /**
     * Sets the impatience type and distribution for a specific job class at this station.
     * Jobs that wait longer than their patience time will abandon the queue.
     * This setting takes precedence over the global class patience.
     *
     * @param jobClass the job class
     * @param impatienceType the type of impatience (RENEGING or BALKING)
     * @param distribution the patience time distribution
     * @throws IllegalArgumentException if distribution is a modulated process or arguments are null
     */
    public void setPatience(JobClass jobClass, ImpatienceType impatienceType, Distribution distribution) {
        if (jobClass == null) {
            throw new IllegalArgumentException("jobClass cannot be null");
        }
        if (impatienceType == null) {
            throw new IllegalArgumentException("impatienceType cannot be null");
        }
        if (distribution == null) {
            throw new IllegalArgumentException("distribution cannot be null");
        }

        // Validate impatience type
        if (impatienceType != ImpatienceType.RENEGING && impatienceType != ImpatienceType.BALKING) {
            throw new IllegalArgumentException("Invalid impatience type. Use ImpatienceType.RENEGING or ImpatienceType.BALKING.");
        }

        // BALKING through setPatience is not supported - use setBalking() with BalkingStrategy instead
        if (impatienceType == ImpatienceType.BALKING) {
            throw new UnsupportedOperationException(
                "BALKING is not supported via setPatience (timer-based). " +
                "Use setBalking(jobClass, BalkingStrategy, thresholds) for state-based balking, " +
                "or use ImpatienceType.RENEGING for timer-based abandonment.");
        }

        // Validate distribution type
        if (distribution instanceof BMAP || distribution instanceof MAP ||
            distribution instanceof DMAP || distribution instanceof MMPP2) {
            throw new IllegalArgumentException(
                "Modulated processes (BMAP, MAP, MMPP2) are not supported for patience distributions.");
        }

        // Check if job class is valid for this network
        if (this.model != null && !this.model.getClasses().contains(jobClass)) {
            throw new IllegalArgumentException("jobClass is not part of this network");
        }

        this.patienceDistributions.put(jobClass, distribution);
        this.impatienceTypes.put(jobClass, impatienceType);
    }

    /**
     * Gets the effective patience distribution for a specific job class.
     * Returns the station-specific setting if available, otherwise falls back
     * to the global class patience.
     *
     * @param jobClass the job class
     * @return the patience distribution, or null if not set
     */
    public Distribution getPatience(JobClass jobClass) {
        // Check station-specific patience first
        if (this.patienceDistributions.containsKey(jobClass)) {
            return this.patienceDistributions.get(jobClass);
        }

        // Fall back to global class patience
        if (jobClass != null) {
            return jobClass.getPatience();
        }

        return null;
    }

    /**
     * Gets the station-specific patience distribution (without fallback).
     *
     * @param jobClass the job class
     * @return the station-specific patience distribution, or null if not set
     */
    public Distribution getPatienceLocal(JobClass jobClass) {
        return this.patienceDistributions.get(jobClass);
    }

    /**
     * Gets the effective impatience type for a specific job class.
     * Returns the station-specific setting if available, otherwise falls back
     * to the global class impatience type.
     *
     * @param jobClass the job class
     * @return the impatience type, or null if not set
     */
    public ImpatienceType getImpatienceType(JobClass jobClass) {
        // Check station-specific impatience type first
        if (this.impatienceTypes.containsKey(jobClass)) {
            return this.impatienceTypes.get(jobClass);
        }

        // Fall back to global class impatience type
        if (jobClass != null) {
            return jobClass.getImpatienceType();
        }

        return null;
    }

    /**
     * Gets the station-specific impatience type (without fallback).
     *
     * @param jobClass the job class
     * @return the station-specific impatience type, or null if not set
     */
    public ImpatienceType getImpatienceTypeLocal(JobClass jobClass) {
        return this.impatienceTypes.get(jobClass);
    }

    /**
     * Checks if patience is configured for a job class (local or global).
     *
     * @param jobClass the job class
     * @return true if patience is configured, false otherwise
     */
    public boolean hasPatience(JobClass jobClass) {
        Distribution dist = getPatience(jobClass);
        return dist != null && !dist.isDisabled();
    }

    /**
     * Checks if patience is configured locally at this station.
     *
     * @param jobClass the job class
     * @return true if local patience is configured, false otherwise
     */
    public boolean hasPatienceLocal(JobClass jobClass) {
        return this.patienceDistributions.containsKey(jobClass);
    }

    // =================== BALKING METHODS ===================

    /**
     * Sets the balking configuration for a specific job class at this station.
     * This setting takes precedence over the global class balking configuration.
     *
     * @param jobClass the job class
     * @param strategy the balking strategy to use
     * @param thresholds list of balking thresholds mapping queue lengths to probabilities
     * @throws IllegalArgumentException if any argument is null or thresholds is empty
     */
    public void setBalking(JobClass jobClass, BalkingStrategy strategy, List<BalkingThreshold> thresholds) {
        if (jobClass == null) {
            throw new IllegalArgumentException("jobClass cannot be null");
        }
        if (strategy == null) {
            throw new IllegalArgumentException("BalkingStrategy cannot be null");
        }
        if (thresholds == null || thresholds.isEmpty()) {
            throw new IllegalArgumentException("Balking thresholds cannot be null or empty");
        }

        // Check if job class is valid for this network
        if (this.model != null && !this.model.getClasses().contains(jobClass)) {
            throw new IllegalArgumentException("jobClass is not part of this network");
        }

        this.balkingStrategies.put(jobClass, strategy);
        this.balkingThresholds.put(jobClass, new ArrayList<BalkingThreshold>(thresholds));
    }

    /**
     * Sets the balking configuration using a single threshold.
     * Convenience method for simple balking configurations.
     *
     * @param jobClass the job class
     * @param strategy the balking strategy to use
     * @param minJobs minimum queue length to trigger balking
     * @param probability balking probability when queue length >= minJobs
     */
    public void setBalking(JobClass jobClass, BalkingStrategy strategy, int minJobs, double probability) {
        List<BalkingThreshold> thresholds = new ArrayList<BalkingThreshold>();
        thresholds.add(new BalkingThreshold(minJobs, probability));
        setBalking(jobClass, strategy, thresholds);
    }

    /**
     * Gets the effective balking strategy for a specific job class.
     * Returns the station-specific setting if available, otherwise falls back
     * to the global class balking strategy.
     *
     * @param jobClass the job class
     * @return the balking strategy, or null if not configured
     */
    public BalkingStrategy getBalkingStrategy(JobClass jobClass) {
        // Check station-specific first
        if (this.balkingStrategies.containsKey(jobClass)) {
            return this.balkingStrategies.get(jobClass);
        }

        // Fall back to global class setting
        if (jobClass != null) {
            return jobClass.getBalkingStrategy();
        }

        return null;
    }

    /**
     * Gets the station-specific balking strategy (without fallback).
     *
     * @param jobClass the job class
     * @return the station-specific balking strategy, or null if not set
     */
    public BalkingStrategy getBalkingStrategyLocal(JobClass jobClass) {
        return this.balkingStrategies.get(jobClass);
    }

    /**
     * Gets the effective balking thresholds for a specific job class.
     * Returns the station-specific setting if available, otherwise falls back
     * to the global class balking thresholds.
     *
     * @param jobClass the job class
     * @return the balking thresholds, or null if not configured
     */
    public List<BalkingThreshold> getBalkingThresholds(JobClass jobClass) {
        // Check station-specific first
        if (this.balkingThresholds.containsKey(jobClass)) {
            return this.balkingThresholds.get(jobClass);
        }

        // Fall back to global class setting
        if (jobClass != null) {
            return jobClass.getBalkingThresholds();
        }

        return null;
    }

    /**
     * Gets the station-specific balking thresholds (without fallback).
     *
     * @param jobClass the job class
     * @return the station-specific balking thresholds, or null if not set
     */
    public List<BalkingThreshold> getBalkingThresholdsLocal(JobClass jobClass) {
        return this.balkingThresholds.get(jobClass);
    }

    /**
     * Checks if balking is configured for a job class (local or global).
     *
     * @param jobClass the job class
     * @return true if balking is configured, false otherwise
     */
    public boolean hasBalking(JobClass jobClass) {
        BalkingStrategy strategy = getBalkingStrategy(jobClass);
        List<BalkingThreshold> thresholds = getBalkingThresholds(jobClass);
        return strategy != null && thresholds != null && !thresholds.isEmpty();
    }

    /**
     * Checks if balking is configured locally at this station.
     *
     * @param jobClass the job class
     * @return true if local balking is configured, false otherwise
     */
    public boolean hasBalkingLocal(JobClass jobClass) {
        return this.balkingStrategies.containsKey(jobClass);
    }

    // =================== RETRIAL METHODS ===================

    /**
     * Sets the retrial configuration for a specific job class at this station.
     * This setting takes precedence over the global class retrial configuration.
     *
     * @param jobClass the job class
     * @param delayDistribution the distribution for retrial delay times
     * @param maxAttempts maximum retry attempts (-1 for unlimited, 0 for no retries)
     * @throws IllegalArgumentException if arguments are null or distribution is modulated
     */
    public void setRetrial(JobClass jobClass, Distribution delayDistribution, int maxAttempts) {
        if (jobClass == null) {
            throw new IllegalArgumentException("jobClass cannot be null");
        }
        if (delayDistribution == null) {
            throw new IllegalArgumentException("Retrial delay distribution cannot be null");
        }
        if (delayDistribution instanceof BMAP || delayDistribution instanceof MAP ||
            delayDistribution instanceof DMAP || delayDistribution instanceof MMPP2) {
            throw new IllegalArgumentException(
                "Modulated processes (BMAP, MAP, MMPP2) are not supported for retrial distributions.");
        }

        // Check if job class is valid for this network
        if (this.model != null && !this.model.getClasses().contains(jobClass)) {
            throw new IllegalArgumentException("jobClass is not part of this network");
        }

        this.retrialDelayDistributions.put(jobClass, delayDistribution);
        this.retrialMaxAttempts.put(jobClass, maxAttempts);

        // Set drop strategy to Retrial so that blocked jobs enter orbit
        if (maxAttempts == 0) {
            // maxAttempts=0 means no retries, revert to drop
            this.setDropRule(jobClass, DropStrategy.Drop);
        } else if (maxAttempts > 0) {
            this.setDropRule(jobClass, DropStrategy.RetrialWithLimit);
        } else {
            this.setDropRule(jobClass, DropStrategy.Retrial);
        }
    }

    /**
     * Sets the retrial configuration with unlimited retry attempts.
     *
     * @param jobClass the job class
     * @param delayDistribution the distribution for retrial delay times
     */
    public void setRetrial(JobClass jobClass, Distribution delayDistribution) {
        setRetrial(jobClass, delayDistribution, -1);
    }

    /**
     * Declares this station to be a retrial queue for a job class: a job that finds
     * every server busy joins an orbit and re-attempts entry after a random delay,
     * instead of waiting in a line.
     *
     * <p>This is the first-class form of the retrial idiom. It removes the waiting
     * room itself, which is what makes the station bufferless, so the caller no
     * longer has to know that setCapacity(nservers) is the way to express "no
     * waiting room, blocked jobs orbit".</p>
     *
     * <p>The capacity must nonetheless leave room for the orbiting jobs: the
     * matrix-analytic engine carries the orbit as a separate level variable, but
     * the state-space engines hold it in the station buffer, so a capacity pinned
     * to the server count would leave them nowhere to put an orbiting job and would
     * silently report an empty orbit.</p>
     *
     * @param jobClass            the job class
     * @param retrialDistribution retrial delay of an orbiting job
     * @param policy              {@link RetrialPolicy#LINEAR} (each orbiting job
     *                            retries at its own rate, aggregate rate n*nu) or
     *                            {@link RetrialPolicy#CONSTANT} (the orbit retries
     *                            as a whole at rate nu when non-empty)
     * @param maxOrbit            orbit capacity; -1 leaves the orbit unbounded. A
     *                            job that finds the orbit full is lost
     */
    public void setOrbit(JobClass jobClass, Distribution retrialDistribution, int policy, int maxOrbit) {
        if (policy != RetrialPolicy.LINEAR && policy != RetrialPolicy.CONSTANT) {
            throw new IllegalArgumentException(
                    "setOrbit requires a RetrialPolicy.LINEAR or RetrialPolicy.CONSTANT policy.");
        }
        if (maxOrbit != -1 && maxOrbit < 0) {
            throw new IllegalArgumentException(
                    "setOrbit requires maxOrbit to be -1 (unbounded) or a non-negative integer.");
        }

        int nservers = this.getNumberOfServers();
        if (nservers < 1 || nservers == Integer.MAX_VALUE) {
            nservers = 1;
            this.setNumberOfServers(nservers);
        }
        if (maxOrbit < 0) {
            this.setCapacity(Integer.MAX_VALUE);
        } else {
            this.setCapacity(nservers + maxOrbit);
        }

        this.setRetrial(jobClass, retrialDistribution, -1);
        this.retrialPolicies.put(jobClass, policy);
        this.orbitMaxJobs.put(jobClass, maxOrbit);
    }

    /**
     * Declares an unbounded orbit with the given retrial policy.
     *
     * @param jobClass            the job class
     * @param retrialDistribution retrial delay of an orbiting job
     * @param policy              the retrial policy
     */
    public void setOrbit(JobClass jobClass, Distribution retrialDistribution, int policy) {
        setOrbit(jobClass, retrialDistribution, policy, -1);
    }

    /**
     * Declares an unbounded orbit with the classical per-customer (LINEAR) policy.
     *
     * @param jobClass            the job class
     * @param retrialDistribution retrial delay of an orbiting job
     */
    public void setOrbit(JobClass jobClass, Distribution retrialDistribution) {
        setOrbit(jobClass, retrialDistribution, RetrialPolicy.LINEAR, -1);
    }

    /**
     * Retrial policy of the orbit of a job class at this station.
     *
     * @param jobClass the job class
     * @return {@link RetrialPolicy#LINEAR} or {@link RetrialPolicy#CONSTANT}
     */
    public int getRetrialPolicy(JobClass jobClass) {
        Integer p = this.retrialPolicies.get(jobClass);
        return (p == null) ? RetrialPolicy.LINEAR : p.intValue();
    }

    /**
     * Orbit capacity of a job class at this station.
     *
     * @param jobClass the job class
     * @return the orbit capacity, or -1 when the orbit is unbounded
     */
    public int getOrbitMaxJobs(JobClass jobClass) {
        Integer m = this.orbitMaxJobs.get(jobClass);
        return (m == null) ? -1 : m.intValue();
    }

    /**
     * Gets the effective retrial delay distribution for a specific job class.
     * Returns the station-specific setting if available, otherwise falls back
     * to the global class retrial distribution.
     *
     * @param jobClass the job class
     * @return the retrial delay distribution, or null if not configured
     */
    public Distribution getRetrialDelayDistribution(JobClass jobClass) {
        // Check station-specific first
        if (this.retrialDelayDistributions.containsKey(jobClass)) {
            return this.retrialDelayDistributions.get(jobClass);
        }

        // Fall back to global class setting
        if (jobClass != null) {
            return jobClass.getRetrialDelayDistribution();
        }

        return null;
    }

    /**
     * Gets the station-specific retrial delay distribution (without fallback).
     *
     * @param jobClass the job class
     * @return the station-specific retrial delay distribution, or null if not set
     */
    public Distribution getRetrialDelayDistributionLocal(JobClass jobClass) {
        return this.retrialDelayDistributions.get(jobClass);
    }

    /**
     * Gets the effective maximum retrial attempts for a specific job class.
     * Returns the station-specific setting if available, otherwise falls back
     * to the global class max attempts.
     *
     * @param jobClass the job class
     * @return maximum attempts (-1 for unlimited), or -1 if not configured
     */
    public int getMaxRetrialAttempts(JobClass jobClass) {
        // Check station-specific first
        if (this.retrialMaxAttempts.containsKey(jobClass)) {
            return this.retrialMaxAttempts.get(jobClass);
        }

        // Fall back to global class setting
        if (jobClass != null) {
            return jobClass.getMaxRetrialAttempts();
        }

        return -1;
    }

    /**
     * Gets the station-specific maximum retrial attempts (without fallback).
     *
     * @param jobClass the job class
     * @return the station-specific max attempts, or null if not set
     */
    public Integer getMaxRetrialAttemptsLocal(JobClass jobClass) {
        return this.retrialMaxAttempts.get(jobClass);
    }

    /**
     * Checks if retrial is configured for a job class (local or global).
     *
     * @param jobClass the job class
     * @return true if retrial is configured, false otherwise
     */
    public boolean hasRetrial(JobClass jobClass) {
        Distribution dist = getRetrialDelayDistribution(jobClass);
        return dist != null && !dist.isDisabled();
    }

    /**
     * Checks if retrial is configured locally at this station.
     *
     * @param jobClass the job class
     * @return true if local retrial is configured, false otherwise
     */
    public boolean hasRetrialLocal(JobClass jobClass) {
        return this.retrialDelayDistributions.containsKey(jobClass);
    }

    /**
     * Sets the orbit impatience (abandonment) distribution for a specific job class.
     * Customers waiting in the retrial orbit may abandon before successfully retrying.
     * This is separate from queue patience (reneging from the waiting queue).
     *
     * @param jobClass the job class
     * @param distribution the orbit abandonment time distribution
     * @throws IllegalArgumentException if distribution is a modulated process or arguments are null
     */
    public void setOrbitImpatience(JobClass jobClass, Distribution distribution) {
        if (jobClass == null) {
            throw new IllegalArgumentException("jobClass cannot be null");
        }
        if (distribution == null) {
            throw new IllegalArgumentException("Orbit impatience distribution cannot be null");
        }
        if (distribution instanceof BMAP || distribution instanceof MAP ||
            distribution instanceof DMAP || distribution instanceof MMPP2) {
            throw new IllegalArgumentException(
                "Modulated processes (BMAP, MAP, MMPP2) are not supported for orbit impatience distributions.");
        }

        // Check if job class is valid for this network
        if (this.model != null && !this.model.getClasses().contains(jobClass)) {
            throw new IllegalArgumentException("jobClass is not part of this network");
        }

        this.orbitImpatienceDistributions.put(jobClass, distribution);
    }

    /**
     * Gets the orbit impatience distribution for a specific job class.
     *
     * @param jobClass the job class
     * @return the orbit impatience distribution, or null if not configured
     */
    public Distribution getOrbitImpatience(JobClass jobClass) {
        return this.orbitImpatienceDistributions.get(jobClass);
    }

    /**
     * Gets the station-specific orbit impatience distribution (without fallback).
     *
     * @param jobClass the job class
     * @return the station-specific orbit impatience distribution, or null if not set
     */
    public Distribution getOrbitImpatienceLocal(JobClass jobClass) {
        return this.orbitImpatienceDistributions.get(jobClass);
    }

    /**
     * Checks if orbit impatience is configured for a job class.
     *
     * @param jobClass the job class
     * @return true if orbit impatience is configured, false otherwise
     */
    public boolean hasOrbitImpatience(JobClass jobClass) {
        Distribution dist = getOrbitImpatience(jobClass);
        return dist != null && !dist.isDisabled();
    }

    /**
     * Sets the batch rejection probability for a specific job class.
     * <p>
     * Used in BMAP/PH/N/N retrial queues. When a batch of size k arrives and only
     * m &lt; k servers are free, with probability p the entire batch is rejected to
     * the orbit, and with probability (1-p) m customers are admitted while the
     * remaining k-m go to the orbit.
     *
     * @param jobClass the job class
     * @param p probability in [0,1] that the batch is rejected rather than
     *          partially admitted; the default is 0 (partial admission allowed)
     * @throws IllegalArgumentException if jobClass is null or p is outside [0,1]
     */
    public void setBatchRejectProbability(JobClass jobClass, double p) {
        if (jobClass == null) {
            throw new IllegalArgumentException("jobClass cannot be null");
        }
        if (p < 0 || p > 1) {
            throw new IllegalArgumentException("Batch reject probability must be in [0, 1].");
        }
        if (this.model != null && !this.model.getClasses().contains(jobClass)) {
            throw new IllegalArgumentException("jobClass is not part of this network");
        }
        this.batchRejectProb.put(jobClass, p);
    }

    /**
     * Gets the batch rejection probability for a specific job class.
     *
     * @param jobClass the job class
     * @return the batch reject probability in [0,1], or 0 if not configured
     */
    public double getBatchRejectProbability(JobClass jobClass) {
        Double p = this.batchRejectProb.get(jobClass);
        if (p == null) {
            return 0.0;
        }
        return p;
    }

}
