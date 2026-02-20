/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.NodeParam;
import jline.lang.ServiceBinding;
import jline.lang.constant.ServerType;
import jline.lang.nodeparam.QueueNodeParam;
import jline.lang.constant.HeteroSchedPolicy;
import jline.lang.constant.PollingType;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SchedStrategyType;
import jline.lang.processes.Distribution;
import jline.lang.sections.*;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static jline.util.Utils.isInf;

/**
 * A queueing station that processes jobs according to various scheduling strategies.
 * 
 * <p>The Queue node is the fundamental service station in queueing networks. It models
 * a system where jobs wait in a buffer and are processed by one or more servers according
 * to a scheduling strategy. Common examples include CPU schedulers, disk I/O queues,
 * network routers, and service desks.</p>
 * 
 * <p>Key features:
 * <ul>
 *   <li>Multiple scheduling strategies: FCFS, PS, LCFS, priority-based, etc.</li>
 *   <li>Single or multi-server configurations</li>
 *   <li>Load-dependent and class-dependent service rates</li>
 *   <li>Support for preemptive and non-preemptive policies</li>
 *   <li>Polling server capabilities with various polling types</li>
 *   <li>Switchover times between job classes</li>
 * </ul>
 * </p>
 * 
 * @see SchedStrategy
 * @see ServiceStation
 * @see Server
 * @since 1.0
 */
public class Queue extends ServiceStation implements Serializable {

    /**
     * Setup time distributions for function tasks (cold start time).
     * Maps job classes to their setup time distributions.
     */
    protected Map<JobClass, Distribution> setupTimes;
    
    /**
     * Delay-off time distributions for function tasks (teardown time).
     * Maps job classes to their delay-off time distributions.
     */
    protected Map<JobClass, Distribution> delayOffTimes;

    /**
     * LPS limit: maximum number of jobs that can execute in PS mode for LPS scheduling.
     * Only used when schedStrategy is LPS.
     */
    private double lpsLimit = Double.NaN;

    /**
     * List of server types for heterogeneous multiserver queues.
     * When non-empty, this queue uses heterogeneous server configuration.
     */
    private List<ServerType> serverTypes;

    /**
     * Scheduling policy for heterogeneous servers.
     * Determines how jobs are assigned to server types when multiple compatible types exist.
     */
    private HeteroSchedPolicy heteroSchedPolicy;

    /**
     * Service distributions per server type and job class for heterogeneous queues.
     * Maps ServerType -> JobClass -> Distribution.
     */
    private Map<ServerType, Map<JobClass, Distribution>> heteroServiceDistributions;

    /**
     * Set of class indices with immediate feedback enabled.
     * When a job self-loops at this station, it stays in service instead of rejoining queue.
     * null means no classes have immediate feedback, "all" placeholder means all classes.
     */
    private java.util.Set<Integer> immediateFeedbackClasses;

    /**
     * Flag indicating if immediate feedback is enabled for all classes.
     */
    private boolean immediateFeedbackAll = false;

    /**
     * Creates a new queueing station with the specified scheduling strategy.
     * Configures the appropriate server type based on the scheduling strategy.
     * The queue is automatically added to the network model.
     * 
     * @param model the network model to add this queue to
     * @param name the name for this queueing station
     * @param schedStrategy the scheduling strategy to use (e.g., FCFS, PS, LCFS)
     * @throws RuntimeException if the scheduling strategy is not supported
     */
    public Queue(Network model, String name, SchedStrategy schedStrategy) {
        super(name);

        this.setModel(model);
        this.model.addNode(this);
        this.schedStrategy = schedStrategy;
        this.input = new Buffer(model.getClasses());
        this.output = new Dispatcher(model.getClasses());

        this.schedStrategy = schedStrategy;
        this.numberOfServers = 1;
        this.setupTimes = new HashMap<JobClass, Distribution>();
        this.delayOffTimes = new HashMap<JobClass, Distribution>();
        this.serverTypes = new ArrayList<ServerType>();
        this.heteroSchedPolicy = HeteroSchedPolicy.ORDER;
        this.heteroServiceDistributions = new HashMap<ServerType, Map<JobClass, Distribution>>();
        this.immediateFeedbackClasses = new java.util.HashSet<Integer>();
        this.immediateFeedbackAll = false;

        switch (this.schedStrategy) {
            case FCFS:
            case LCFS:
            case SIRO:
            case SEPT:
            case LEPT:
            case SJF:
            case LJF:
            case EDD:
            case SRPT:
            case SRPTPRIO:
            case SETF:
                this.schedPolicy = SchedStrategyType.NP;
                this.server = new Server(model.getClasses());
                break;
            case INF:
                this.schedPolicy = SchedStrategyType.NP;
                this.server = new InfiniteServer(model.getClasses());
                this.numberOfServers = Integer.MAX_VALUE;
                break;
            case HOL:
            case FCFSPRIO:
                this.schedPolicy = SchedStrategyType.NP;
                this.server = new Server(model.getClasses());
                break;
            case PS:
            case DPS:
            case GPS:
            case PSPRIO:
            case DPSPRIO:
            case GPSPRIO:
            case LPS:
                this.schedPolicy = SchedStrategyType.PR;
                this.server = new SharedServer(model.getClasses());
                break;
            case LCFSPR:
            case LCFSPRPRIO:
            case FCFSPR:
            case FCFSPRPRIO:
            case EDF:
            case FB:
            case PSJF:
            case LRPT:
                this.schedPolicy = SchedStrategyType.PR;
                this.server = new PreemptiveServer(model.getClasses());
                break;
            case LCFSPI:
            case LCFSPIPRIO:
                this.schedPolicy = SchedStrategyType.PNR;
                this.server = new PreemptiveServer(model.getClasses());
                break;
            case LCFSPRIO:
                this.schedPolicy = SchedStrategyType.NPPrio;
                this.server = new PreemptiveServer(model.getClasses());
                break;
            case POLLING:
                this.schedPolicy = SchedStrategyType.NP;
                this.server = new PollingServer(model.getClasses());
                break;
            default:
                throw new RuntimeException("Routing Strategy is not supported in JLINE");
        }
    }

    /**
     * Creates a new queueing station with default processor sharing (PS) scheduling.
     * 
     * @param model The network model to add this queue to
     * @param name The name for this queueing station
     */
    public Queue(Network model, String name) {
        this(model, name, SchedStrategy.PS);
    }

    /**
     * Checks whether this queue has a service process configured for the specified job class.
     * 
     * @param jobClass The job class to check
     * @return true if a service process exists for this job class, false otherwise
     */
    public boolean containsJobClass(JobClass jobClass) {
        for (ServiceBinding serviceProcess : this.serviceProcesses) {
            // Use index comparison to handle Signal resolution (Signal -> OpenSignal/ClosedSignal)
            if (serviceProcess.getJobClass().getIndex() == jobClass.getIndex()) {
                return true;
            }
        }
        return false;
    }

    /**
     * Gets the scheduling policy type (preemptive or non-preemptive).
     * 
     * @return The scheduling policy type (PR for preemptive, NP for non-preemptive)
     */
    public SchedStrategyType getSchedPolicy() {
        return this.schedPolicy;
    }

    /**
     * Gets the scheduling strategy used by this queue.
     * 
     * @return The scheduling strategy (e.g., FCFS, PS, LCFS, etc.)
     */
    public SchedStrategy getSchedStrategy() {
        return this.schedStrategy;
    }

    /**
     * Gets the scheduling strategy parameter for a specific job class.
     * 
     * <p>For strategies like DPS (Discriminatory Processor Sharing) and GPS
     * (Generalized Processor Sharing), this returns the weight or priority
     * parameter for the job class.</p>
     * 
     * @param jobClass The job class to get the parameter for
     * @return The scheduling parameter value, or 0.0 if not set
     */
    public double getSchedStrategyPar(JobClass jobClass) {
        // For LPS, return the station-wide limit (stored in lpsLimit field)
        // This will be used to populate schedparam Matrix for all classes
        if (this.schedStrategy == SchedStrategy.LPS) {
            return this.lpsLimit;
        }
        return this.schedStrategyPar.getOrDefault(jobClass, 0.0);
    }

    /**
     * Gets the service time distribution for a specific job class.
     * 
     * @param jobClass The job class to get the service distribution for
     * @return The service time distribution for this job class
     */
    public Distribution getService(JobClass jobClass) {
        return this.server.getServiceDistribution(jobClass);
    }

    /**
     * Prints a summary of this queue's configuration to standard output.
     * 
     * <p>The summary includes the queue name, service processes for each job class,
     * their mean service times and squared coefficients of variation, number of servers,
     * and output routing configuration.</p>
     */
    @Override
    public void printSummary() {
        System.out.format("jline.Queue:\n");
        System.out.format("--Name: %s\n", this.getName());
        System.out.format("--Service Processes:\n");
        for (JobClass jobClass : this.model.getClasses()) {
            System.out.format("----%s: %s (Mean: %g, SCV: %g)\n", jobClass.getName(), this.getServiceProcess(jobClass).toString(), this.getServiceProcess(jobClass).getMean(), this.getServiceProcess(jobClass).getSCV());
        }
        if (isInf(this.getNumberOfServers())) {
            System.out.format("--Number of Servers: Inf\n");
        } else {
            System.out.format("--Number of Servers: %d\n", this.getNumberOfServers());
        }
        this.output.printSummary();
    }

    /**
     * Sets a class-dependent scaling function for service rates.
     * 
     * <p>The function takes a matrix representing the number of jobs of each class
     * at the station and returns a scaling factor for the service rate. This enables
     * modeling of systems where service rates depend on the job mix.</p>
     * 
     * @param beta A function that maps job class populations to a service rate scaling factor
     * @throws RuntimeException if the scheduling strategy doesn't support class dependence
     */
    public void setClassDependence(SerializableFunction<Matrix, Double> beta) {
        switch (this.schedStrategy) {
            case PS:
            case FCFS:
                this.setLimitedClassDependence(beta);
                break;
            default:
                throw new RuntimeException("Class-dependence supported only for processor sharing (PS) and first-come first-serve (FCFS) stations.");

        }
    }

    /**
     * Sets load-dependent service rate scaling factors.
     *
     * <p>Each element alpha[n] specifies the service rate scaling when there are
     * n jobs at the station. This enables modeling of systems where performance
     * degrades under load.</p>
     *
     * @param alpha A matrix of scaling factors indexed by the number of jobs
     * @throws RuntimeException if the scheduling strategy doesn't support load dependence
     */
    public void setLoadDependence(Matrix alpha) {
        switch (this.schedStrategy) {
            case PS:
            case FCFS:
                this.setLimitedLoadDependence(alpha);
                break;
            default:
                throw new RuntimeException("Load-dependence supported only for processor sharing (PS) and first-come first-serve (FCFS) stations.");
        }
    }

    /**
     * Sets joint class-dependent scaling for service rates using a lookup table.
     *
     * <p>The scaling table is indexed by the per-class population vector (n1, n2, ..., nR)
     * and returns the service rate scaling factor for that state. The table is stored in
     * linearized form with index = n1 + n2*(N1+1) + n3*(N1+1)*(N2+1) + ...</p>
     *
     * @param scalingTable Linearized lookup table for (n1, n2, ..., nR) -> scaling
     * @param cutoffs Per-class cutoffs [N1, N2, ..., NR]
     * @throws RuntimeException if the scheduling strategy doesn't support joint dependence
     */
    public void setJointDependence(Matrix scalingTable, Matrix cutoffs) {
        switch (this.schedStrategy) {
            case PS:
            case FCFS:
                // Validate table size matches cutoffs
                int expectedSize = 1;
                for (int k = 0; k < cutoffs.length(); k++) {
                    expectedSize *= (int)(cutoffs.get(k) + 1);
                }
                if (scalingTable.length() != expectedSize) {
                    throw new RuntimeException(String.format(
                        "Scaling table size (%d) does not match expected size from cutoffs (%d).",
                        scalingTable.length(), expectedSize));
                }
                this.setLimitedJointDependence(scalingTable, cutoffs);
                break;
            default:
                throw new RuntimeException("Joint-dependence supported only for processor sharing (PS) and first-come first-serve (FCFS) stations.");
        }
    }

    /**
     * Sets joint class-dependent scaling with auto-computed cutoffs.
     *
     * <p>Cutoffs are automatically computed using the formula: ceil(6000^(1/(M*K)))
     * where M is the number of stations and K is the number of classes.</p>
     *
     * @param scalingTable Linearized lookup table for (n1, n2, ..., nR) -> scaling
     * @throws RuntimeException if the scheduling strategy doesn't support joint dependence
     */
    public void setJointDependence(Matrix scalingTable) {
        int M = this.model.getNumberOfStations();
        int K = this.model.getNumberOfClasses();
        int defaultCutoff = (int) Math.ceil(Math.pow(6000, 1.0 / (M * K)));
        Matrix cutoffs = new Matrix(1, K);
        cutoffs.fill(defaultCutoff);
        setJointDependence(scalingTable, cutoffs);
    }

    /**
     * Sets per-class joint class-dependent scaling for service rates.
     *
     * <p>Unlike setJointDependence which uses a single scaling factor for all classes,
     * this method allows each class to have its own scaling table indexed by the
     * per-class population vector (n1, n2, ..., nK).</p>
     *
     * <p>This is essential for Flow-Equivalent Server (FES) aggregation where the
     * service rate for class c in state (n1,...,nK) equals the throughput of class c
     * in an isolated subnetwork.</p>
     *
     * @param scalingTables Map from JobClass to linearized scaling table
     * @param cutoffs Per-class cutoffs [N1, N2, ..., NR]
     * @throws RuntimeException if the scheduling strategy doesn't support joint class dependence
     */
    public void setJointClassDependence(Map<JobClass, Matrix> scalingTables, Matrix cutoffs) {
        switch (this.schedStrategy) {
            case PS:
            case FCFS:
                // Validate table sizes match cutoffs
                int expectedSize = 1;
                for (int k = 0; k < cutoffs.length(); k++) {
                    expectedSize *= (int)(cutoffs.get(k) + 1);
                }
                for (Map.Entry<JobClass, Matrix> entry : scalingTables.entrySet()) {
                    if (entry.getValue().length() != expectedSize) {
                        throw new RuntimeException(String.format(
                            "Scaling table for class %s has size (%d) which does not match expected size from cutoffs (%d).",
                            entry.getKey().getName(), entry.getValue().length(), expectedSize));
                    }
                }
                this.setLimitedJointClassDependence(scalingTables, cutoffs);
                break;
            default:
                throw new RuntimeException("Joint-class-dependence supported only for processor sharing (PS) and first-come first-serve (FCFS) stations.");
        }
    }

    /**
     * Sets the number of servers at this queueing station.
     * 
     * <p>For infinite server (IS) queues, this method has no effect as they
     * always have unlimited servers.</p>
     * 
     * @param numberOfServers The number of parallel servers (must be positive)
     */
    public void setNumberOfServers(int numberOfServers) {
        /*if (this.schedStrategy == SchedStrategy.DPS || this.schedStrategy == SchedStrategy.GPS) {
            if (numberOfServers != 1) {
                throw new InvalidParameterException("Invalid number of servers for scheduling strategy");
            }
        }*/
        if (this.schedStrategy != SchedStrategy.INF) {
            this.numberOfServers = numberOfServers;
        }
    }

    /**
     * Sets the scheduling strategy parameter for a specific job class.
     * 
     * <p>For weighted scheduling strategies (DPS, GPS, etc.), this sets the
     * weight or priority parameter that determines the job class's share of
     * the service capacity.</p>
     * 
     * @param jobClass The job class to set the parameter for
     * @param weight The scheduling parameter value (e.g., weight, priority)
     */
    public void setSchedStrategyPar(JobClass jobClass, double weight) {
        this.schedStrategyPar.put(jobClass, weight);
    }

    /**
     * Sets the maximum number of jobs for LPS scheduling.
     * For LPS: limit is the max number of jobs that can execute in PS mode.
     *
     * @param limit the maximum number of concurrent jobs in PS mode
     * @throws RuntimeException if called on non-LPS queue
     */
    public void setLimit(int limit) {
        if (this.schedStrategy != SchedStrategy.LPS) {
            throw new RuntimeException("setLimit() can only be called on queues with LPS scheduling strategy");
        }
        // Store limit as station-wide parameter for LPS
        this.lpsLimit = (double) limit;
    }

    /**
     * Sets the polling type for this queue (only valid for POLLING scheduling strategy).
     *
     * @param pollingType the polling type (GATED, EXHAUSTIVE, or KLIMITED)
     */
    public void setPollingType(PollingType pollingType) {
        if (this.schedStrategy != SchedStrategy.POLLING) {
            throw new RuntimeException("setPollingType() can only be called on queues with POLLING scheduling strategy");
        }
        if (this.server instanceof PollingServer) {
            ((PollingServer) this.server).setPollingType(pollingType);
        }
        // Update NodeParam if network is available
        if (this.model != null) {
            NetworkStruct sn = this.model.getStruct();
            if (sn != null && sn.nodeparam != null) {
                NodeParam param = sn.nodeparam.get(this);
                if (param instanceof QueueNodeParam) {
                    ((QueueNodeParam) param).pollingType = pollingType;
                }
            }
        }
    }

    /**
     * Sets the polling type for this queue with K value for K-LIMITED (only valid for POLLING scheduling strategy).
     *
     * @param pollingType the polling type (GATED, EXHAUSTIVE, or KLIMITED)
     * @param k the K value for K-LIMITED polling (ignored for other types)
     */
    public void setPollingType(PollingType pollingType, int k) {
        if (this.schedStrategy != SchedStrategy.POLLING) {
            throw new RuntimeException("setPollingType() can only be called on queues with POLLING scheduling strategy");
        }
        if (this.server instanceof PollingServer) {
            ((PollingServer) this.server).setPollingType(pollingType, k);
        }
        // Update NodeParam if network is available
        if (this.model != null) {
            NetworkStruct sn = this.model.getStruct();
            if (sn != null && sn.nodeparam != null) {
                NodeParam param = sn.nodeparam.get(this);
                if (param instanceof QueueNodeParam) {
                    ((QueueNodeParam) param).pollingType = pollingType;
                    if (pollingType == PollingType.KLIMITED) {
                        ((QueueNodeParam) param).pollingPar = k;
                    }
                }
            }
        }
    }

    /**
     * Sets the K value for K-LIMITED polling (only valid for POLLING scheduling strategy with K-LIMITED type).
     *
     * @param k the K value (must be greater than 0)
     */
    public void setPollingK(int k) {
        if (this.schedStrategy != SchedStrategy.POLLING) {
            throw new RuntimeException("setPollingK() can only be called on queues with POLLING scheduling strategy");
        }
        if (this.server instanceof PollingServer) {
            ((PollingServer) this.server).setPollingK(k);
        }
        // Update NodeParam if network is available
        if (this.model != null) {
            NetworkStruct sn = this.model.getStruct();
            if (sn != null && sn.nodeparam != null) {
                NodeParam param = sn.nodeparam.get(this);
                if (param instanceof QueueNodeParam) {
                    ((QueueNodeParam) param).pollingPar = k;
                }
            }
        }
    }

    /**
     * Sets the switchover time for a job class (only valid for POLLING scheduling strategy).
     *
     * @param jobClass the job class
     * @param switchoverTime the switchover time distribution
     */
    public void setSwitchover(JobClass jobClass, Distribution switchoverTime) {
        // Validate input parameters
        if (jobClass == null) {
            throw new IllegalArgumentException("jobClass cannot be null");
        }
        if (switchoverTime == null) {
            throw new IllegalArgumentException("switchoverTime cannot be null");
        }
        
        // Check scheduling strategy
        if (this.schedStrategy != SchedStrategy.POLLING) {
            throw new RuntimeException("setSwitchover() can only be called on queues with POLLING scheduling strategy");
        }
        
        // Check if job class is valid for this network
        if (!this.model.getClasses().contains(jobClass)) {
            throw new IllegalArgumentException("jobClass is not part of this network");
        }
        
        if (this.server instanceof PollingServer) {
            ((PollingServer) this.server).setSwitchover(jobClass, switchoverTime);
        }
        // Update NodeParam if network is available
        if (this.model != null) {
            NetworkStruct sn = this.model.getStruct();
            if (sn != null && sn.nodeparam != null) {
                NodeParam param = sn.nodeparam.get(this);
                if (param instanceof QueueNodeParam) {
                    ((QueueNodeParam) param).switchoverTime.put(jobClass, switchoverTime);
                }
            }
        }
    }

    /**
     * Sets the switchover time from one job class to another (for general scheduling strategies).
     *
     * @param fromClass the job class to switch from
     * @param toClass the job class to switch to
     * @param switchoverTime the switchover time distribution
     */
    public void setSwitchover(JobClass fromClass, JobClass toClass, Distribution switchoverTime) {
        // Validate input parameters
        if (fromClass == null) {
            throw new IllegalArgumentException("fromClass cannot be null");
        }
        if (toClass == null) {
            throw new IllegalArgumentException("toClass cannot be null");
        }
        if (switchoverTime == null) {
            throw new IllegalArgumentException("switchoverTime cannot be null");
        }
        
        // Check if job classes are valid for this network
        if (!this.model.getClasses().contains(fromClass)) {
            throw new IllegalArgumentException("fromClass is not part of this network");
        }
        if (!this.model.getClasses().contains(toClass)) {
            throw new IllegalArgumentException("toClass is not part of this network");
        }
        
        // For POLLING strategy, delegate to PollingServer but still store in general storage
        if (this.schedStrategy == SchedStrategy.POLLING) {
            if (this.server instanceof PollingServer) {
                // For polling, set the switchover time for the fromClass (ignoring toClass)
                ((PollingServer) this.server).setSwitchover(fromClass, switchoverTime);
            }
        }
        
        // Store in general switchover time storage for all scheduling strategies
        this.setSwitchoverTime(fromClass, toClass, switchoverTime);
    }

    /**
     * Gets the switchover time from one job class to another.
     *
     * @param fromClass the job class to switch from
     * @param toClass the job class to switch to
     * @return the switchover time distribution, or null if not set
     */
    public Distribution getSwitchover(JobClass fromClass, JobClass toClass) {
        return this.getSwitchoverTime(fromClass, toClass);
    }

    /**
     * Gets the switchover time for a job class (POLLING scheduling strategy).
     *
     * @param jobClass the job class
     * @return the switchover time distribution, or null if not set
     */
    public Distribution getSwitchover(JobClass jobClass) {
        if (this.schedStrategy == SchedStrategy.POLLING && this.server instanceof PollingServer) {
            return ((PollingServer) this.server).getSwitchover(jobClass);
        }
        return null;
    }
    
    /**
     * Sets the setup time and delay off time for a job class.
     * This is typically used for function-based tasks that have initialization overhead.
     *
     * @param jobClass the job class
     * @param setupTime the setup time distribution
     * @param delayoffTime the delay off time distribution
     */
    public void setDelayOff(JobClass jobClass, Distribution setupTime, Distribution delayoffTime) {
        // Validate input parameters
        if (jobClass == null) {
            throw new IllegalArgumentException("jobClass cannot be null");
        }
        if (setupTime == null) {
            throw new IllegalArgumentException("setupTime cannot be null");
        }
        if (delayoffTime == null) {
            throw new IllegalArgumentException("delayoffTime cannot be null");
        }
        
        // Check if job class is valid for this network
        if (!this.model.getClasses().contains(jobClass)) {
            throw new IllegalArgumentException("jobClass is not part of this network");
        }
        
        // Store the setup and delay-off times
        this.setupTimes.put(jobClass, setupTime);
        this.delayOffTimes.put(jobClass, delayoffTime);
        
        // Update NodeParam if network is available
        if (this.model != null) {
            NetworkStruct sn = this.model.getStruct();
            if (sn != null && sn.nodeparam != null) {
                NodeParam param = sn.nodeparam.get(this);
                if (param != null) {
                    // Store in NodeParam for network structure access
                    // This ensures the values are available during analysis
                    if (param.withMemory == null) {
                        param.withMemory = new HashMap<JobClass, Matrix>();
                    }
                    // Note: We're storing references in withMemory as a temporary solution
                    // In a full implementation, NodeParam might need dedicated fields
                }
            }
        }
    }
    
    /**
     * Gets the setup time distribution for a job class.
     *
     * @param jobClass the job class
     * @return the setup time distribution, or null if not set
     */
    public Distribution getSetupTime(JobClass jobClass) {
        return this.setupTimes.get(jobClass);
    }
    
    /**
     * Gets the delay-off time distribution for a job class.
     *
     * @param jobClass the job class
     * @return the delay-off time distribution, or null if not set
     */
    public Distribution getDelayOffTime(JobClass jobClass) {
        return this.delayOffTimes.get(jobClass);
    }

    /**
     * Checks if this queue has delay-off times enabled.
     * Delay-off is considered enabled if any job class has both
     * a setup time and a delay-off time distribution configured.
     *
     * @return true if delay-off is enabled, false otherwise
     */
    public boolean isDelayOffEnabled() {
        return !this.setupTimes.isEmpty() && !this.delayOffTimes.isEmpty();
    }

    // ==================== Heterogeneous Server Methods ====================

    /**
     * Adds a server type to this queue for heterogeneous multiserver configuration.
     * <p>
     * When server types are added, the queue becomes a heterogeneous multiserver queue
     * where different server types can have different service rates and serve different
     * subsets of job classes.
     * <p>
     * The total number of servers at this queue becomes the sum of all server type counts.
     *
     * @param serverType the server type to add
     * @throws IllegalArgumentException if serverType is null or already added
     */
    public void addServerType(ServerType serverType) {
        if (serverType == null) {
            throw new IllegalArgumentException("Server type cannot be null");
        }
        if (this.serverTypes.contains(serverType)) {
            throw new IllegalArgumentException("Server type '" + serverType.getName() + "' is already added to this queue");
        }

        // Assign ID and parent
        serverType.setId(this.serverTypes.size());
        serverType.setParentQueue(this);
        this.serverTypes.add(serverType);

        // Initialize service distribution map for this server type
        this.heteroServiceDistributions.put(serverType, new HashMap<JobClass, Distribution>());

        // Update total number of servers
        updateTotalServerCount();
    }

    /**
     * Updates the total numberOfServers based on all server types.
     */
    private void updateTotalServerCount() {
        if (this.serverTypes.isEmpty()) {
            return;
        }
        int total = 0;
        for (ServerType st : this.serverTypes) {
            total += st.getNumOfServers();
        }
        this.numberOfServers = total;
    }

    /**
     * Gets the list of server types configured for this queue.
     *
     * @return a new list containing the server types
     */
    public List<ServerType> getServerTypes() {
        return new ArrayList<ServerType>(this.serverTypes);
    }

    /**
     * Gets the number of server types configured for this queue.
     *
     * @return the number of server types, or 0 if homogeneous
     */
    public int getNumServerTypes() {
        return this.serverTypes.size();
    }

    /**
     * Checks if this queue is configured as a heterogeneous multiserver queue.
     *
     * @return true if server types are defined, false for homogeneous queue
     */
    public boolean isHeterogeneous() {
        return !this.serverTypes.isEmpty();
    }

    /**
     * Sets the scheduling policy for heterogeneous servers.
     * <p>
     * This policy determines how jobs are assigned to server types when a job's
     * class is compatible with multiple server types.
     *
     * @param policy the heterogeneous scheduling policy
     * @throws IllegalArgumentException if policy is null
     */
    public void setHeteroSchedPolicy(HeteroSchedPolicy policy) {
        if (policy == null) {
            throw new IllegalArgumentException("Heterogeneous scheduling policy cannot be null");
        }
        this.heteroSchedPolicy = policy;
    }

    /**
     * Gets the scheduling policy for heterogeneous servers.
     *
     * @return the heterogeneous scheduling policy
     */
    public HeteroSchedPolicy getHeteroSchedPolicy() {
        return this.heteroSchedPolicy;
    }

    /**
     * Sets the service time distribution for a specific job class and server type.
     * <p>
     * This method is used for heterogeneous multiserver queues where different
     * server types may have different service rates for the same job class.
     *
     * @param jobClass the job class
     * @param serverType the server type
     * @param distribution the service time distribution
     * @throws IllegalArgumentException if any parameter is null or serverType is not in this queue
     */
    public void setService(JobClass jobClass, ServerType serverType, Distribution distribution) {
        if (jobClass == null) {
            throw new IllegalArgumentException("Job class cannot be null");
        }
        if (serverType == null) {
            throw new IllegalArgumentException("Server type cannot be null");
        }
        if (distribution == null) {
            throw new IllegalArgumentException("Distribution cannot be null");
        }
        if (!this.serverTypes.contains(serverType)) {
            throw new IllegalArgumentException("Server type '" + serverType.getName() +
                "' is not added to this queue. Call addServerType() first.");
        }

        // Store the heterogeneous service distribution
        Map<JobClass, Distribution> classMap = this.heteroServiceDistributions.get(serverType);
        if (classMap == null) {
            classMap = new HashMap<JobClass, Distribution>();
            this.heteroServiceDistributions.put(serverType, classMap);
        }
        classMap.put(jobClass, distribution);

        // Also ensure the job class is marked as compatible with this server type
        if (!serverType.isCompatible(jobClass)) {
            serverType.addCompatibleClass(jobClass);
        }
    }

    /**
     * Gets the service time distribution for a specific job class and server type.
     *
     * @param jobClass the job class
     * @param serverType the server type
     * @return the service time distribution, or null if not set
     */
    public Distribution getService(JobClass jobClass, ServerType serverType) {
        if (serverType == null || !this.heteroServiceDistributions.containsKey(serverType)) {
            return null;
        }
        Map<JobClass, Distribution> classMap = this.heteroServiceDistributions.get(serverType);
        return classMap != null ? classMap.get(jobClass) : null;
    }

    /**
     * Gets all heterogeneous service distributions for this queue.
     *
     * @return the map of server type to job class to distribution
     */
    public Map<ServerType, Map<JobClass, Distribution>> getHeteroServiceDistributions() {
        return this.heteroServiceDistributions;
    }

    /**
     * Gets a server type by its ID.
     *
     * @param id the server type ID
     * @return the server type, or null if not found
     */
    public ServerType getServerType(int id) {
        if (id >= 0 && id < this.serverTypes.size()) {
            return this.serverTypes.get(id);
        }
        return null;
    }

    /**
     * Gets a server type by its name.
     *
     * @param name the server type name
     * @return the server type, or null if not found
     */
    public ServerType getServerType(String name) {
        for (ServerType st : this.serverTypes) {
            if (st.getName().equals(name)) {
                return st;
            }
        }
        return null;
    }

    /**
     * Checks if all job classes have at least one compatible server type.
     * <p>
     * This validation is important because in heterogeneous queues, every job class
     * must be able to be served by at least one server type.
     *
     * @return true if all job classes have at least one compatible server type
     */
    public boolean validateCompatibility() {
        if (!isHeterogeneous()) {
            return true; // Homogeneous queues are always valid
        }

        for (JobClass jobClass : this.model.getClasses()) {
            boolean hasCompatible = false;
            for (ServerType st : this.serverTypes) {
                if (st.isCompatible(jobClass)) {
                    hasCompatible = true;
                    break;
                }
            }
            if (!hasCompatible) {
                return false;
            }
        }
        return true;
    }

    // ==================== Immediate Feedback Methods ====================

    /**
     * Enables or disables immediate feedback for all job classes at this queue.
     * When enabled, jobs that self-loop at this station stay in service instead
     * of rejoining the queue.
     *
     * @param enabled true to enable immediate feedback for all classes, false to disable
     */
    public void setImmediateFeedback(boolean enabled) {
        this.immediateFeedbackAll = enabled;
        if (!enabled) {
            this.immediateFeedbackClasses.clear();
        }
    }

    /**
     * Enables immediate feedback for a specific job class at this queue.
     *
     * @param jobClass the job class to enable immediate feedback for
     */
    public void setImmediateFeedback(JobClass jobClass) {
        if (jobClass != null) {
            this.immediateFeedbackClasses.add(jobClass.getIndex());
        }
    }

    /**
     * Enables immediate feedback for multiple job classes at this queue.
     *
     * @param jobClasses list of job classes to enable immediate feedback for
     */
    public void setImmediateFeedbackForClasses(java.util.List<JobClass> jobClasses) {
        if (jobClasses != null) {
            for (JobClass jc : jobClasses) {
                if (jc != null) {
                    this.immediateFeedbackClasses.add(jc.getIndex());
                }
            }
        }
    }

    /**
     * Checks if immediate feedback is enabled for any class at this queue.
     *
     * @return true if immediate feedback is enabled for at least one class
     */
    public boolean hasImmediateFeedback() {
        return this.immediateFeedbackAll || !this.immediateFeedbackClasses.isEmpty();
    }

    /**
     * Checks if immediate feedback is enabled for a specific class at this queue.
     *
     * @param classId the class index (0-based)
     * @return true if immediate feedback is enabled for the specified class
     */
    public boolean hasImmediateFeedback(int classId) {
        return this.immediateFeedbackAll || this.immediateFeedbackClasses.contains(classId);
    }

    /**
     * Gets the set of class indices with immediate feedback enabled.
     *
     * @return set of class indices, or null if none are set
     */
    public java.util.Set<Integer> getImmediateFeedbackClasses() {
        if (this.immediateFeedbackAll) {
            return null; // Indicates "all"
        }
        return this.immediateFeedbackClasses;
    }

    /**
     * Checks if immediate feedback is enabled for all classes.
     *
     * @return true if immediate feedback is enabled for all classes
     */
    public boolean isImmediateFeedbackAll() {
        return this.immediateFeedbackAll;
    }

}

