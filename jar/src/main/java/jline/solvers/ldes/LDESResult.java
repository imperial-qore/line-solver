/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ldes;

import jline.lang.NetworkStruct;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

import java.util.Map;

/**
 * @brief Result container for LINE Discrete Event Simulator (LDES) solver computations.
 *
 * @details LDESResult extends the base SolverResult to include LDES-specific simulation
 * results using the SSJ library. LDES uses discrete event simulation to generate
 * sample paths and estimate performance metrics.
 *
 * @section desresult_metrics Performance Metrics
 * All metrics are stored as matrices with dimensions [stations x classes]:
 * - **QN**: Average queue lengths (number of jobs)
 * - **UN**: Server utilization (fraction of busy time)
 * - **RN**: Response times (time from arrival to departure)
 * - **TN**: Throughput (jobs completed per time unit)
 * - **CN**: Visit counts (for residence time calculations)
 * - **XN**: System throughput per class
 *
 * @section desresult_transient Transient Analysis Results
 * For transient analysis (finite time horizon), additional results are available:
 * - **QNt**: Queue length evolution over time [stations][classes] arrays of time series
 * - **UNt**: Utilization evolution over time
 * - **TNt**: Throughput evolution over time
 * - **t**: Time points for transient measurements
 *
 * @section desresult_state State Information
 * - **tranSysState**: System state at each time point (for state-space analysis)
 * - **tranSync**: Synchronization information for multi-class systems
 * - **sn**: Network structure used during simulation
 *
 * @see SolverLDES Main LDES solver class
 * @see SolverResult Base result class
 * @since 1.0
 * @author QORE Lab, Imperial College London
 */
public class LDESResult extends SolverResult {

    /** Transient system state matrices indexed by time points */
    public Map<Integer, Matrix> tranSysState;

    /** Transient synchronization matrix for multi-class coordination */
    public Matrix tranSync;

    /** Network structure used for the analysis */
    public NetworkStruct sn;

    /** Confidence interval half-widths for queue lengths [stations x classes] */
    public Matrix QNCI;

    /** Confidence interval half-widths for utilizations [stations x classes] */
    public Matrix UNCI;

    /** Confidence interval half-widths for response times [stations x classes] */
    public Matrix RNCI;

    /** Confidence interval half-widths for throughputs [stations x classes] */
    public Matrix TNCI;

    /** Confidence interval half-widths for arrival rates [stations x classes] */
    public Matrix ANCI;

    /** Confidence interval half-widths for residence times [stations x classes] */
    public Matrix WNCI;

    /**
     * Variance reduction method used during simulation.
     * Values: "none" (default), "antithetic", "control", "both"
     * - "none": Standard simulation without variance reduction
     * - "antithetic": Antithetic variates using synchronized 1-U method
     * - "control": Control variates using mean-based correction
     * - "both": Combined antithetic and control variates
     */
    public String varianceReductionMethod;

    /** Flag indicating whether the simulation converged before reaching max events/time */
    public boolean converged = false;

    /**
     * Stopping reason: "convergence", "max_events" (service-completion budget,
     * options.samples/events, reached), "max_sim_events" (all-event cap,
     * LDESOptions.maxSimEvents, reached), or "max_time" (wall-clock budget).
     */
    public String stoppingReason = "";

    /** Number of batches used for final confidence interval estimation */
    public int convergenceBatches = 0;

    /** Total number of simulated service completion events */
    public long totalSimulatedEvents = 0;

    /** Number of samples (observations) used for queue length estimation [stations x classes] */
    public Matrix QNSamples;

    /** Number of samples (observations) used for utilization estimation [stations x classes] */
    public Matrix UNSamples;

    /** Number of samples (observations) used for response time estimation [stations x classes] */
    public Matrix RNSamples;

    /** Number of samples (observations) used for throughput estimation [stations x classes] */
    public Matrix TNSamples;

    /** Relative precision achieved for queue lengths (CI half-width / mean) [stations x classes] */
    public Matrix QNRelPrec;

    /** Relative precision achieved for utilizations (CI half-width / mean) [stations x classes] */
    public Matrix UNRelPrec;

    /** Relative precision achieved for response times (CI half-width / mean) [stations x classes] */
    public Matrix RNRelPrec;

    /** Relative precision achieved for throughputs (CI half-width / mean) [stations x classes] */
    public Matrix TNRelPrec;

    // ==================== Impatience Statistics ====================

    /** Number of customers who reneged (abandoned queue due to expired patience) [stations x classes] */
    public Matrix renegedCustomers;

    /** Average wait time before reneging [stations x classes] */
    public Matrix avgRenegingWaitTime;

    /** Reneging rate (reneged / (completed + reneged + dropped)) [stations x classes] */
    public Matrix renegingRate;

    /** Number of customers who balked (refused to join upon arrival) [stations x classes] */
    public Matrix balkedCustomers;

    /** Balking probability (balked / arrivals) [stations x classes] */
    public Matrix balkingProbability;

    /** Number of successful retrial attempts (customers who re-entered queue from orbit) [stations x classes] */
    public Matrix retriedCustomers;

    /** Number of customers dropped after exceeding max retrial attempts [stations x classes] */
    public Matrix retrialDropped;

    /** Average orbit size (time-weighted mean jobs in orbit) [stations x classes] */
    public Matrix avgOrbitSize;

    // ==================== Response Time Samples for CDF Computation ====================

    /**
     * Response time samples for empirical CDF computation [stations][classes].
     * Each element is a list of individual response time observations.
     */
    public java.util.List<Double>[][] respTimeSamples;

    /**
     * Passage time samples for empirical CDF computation [stations][classes].
     * Passage times track time from arrival to departure including routing delays.
     */
    public java.util.List<Double>[][] passTimeSamples;

    // ==================== Transient Probability Results ====================

    /**
     * Time points for transient probability analysis.
     */
    public Matrix tranProbT;

    /**
     * Transient state probabilities over time [numTimePoints x numStates].
     * Each row contains probability distribution at that time point.
     */
    public Matrix tranProbPit;

    /**
     * Transient state space matrix [numStates x stateVectorLength].
     * Defines the state corresponding to each column in tranProbPit.
     */
    public Matrix tranProbStateSpace;

    /**
     * Transient aggregated state probabilities per node [nodeIdx] -> [numTimePoints x numAggrStates].
     */
    public Map<Integer, Matrix> tranProbAggrPit;

    /**
     * Transient aggregated state space per node [nodeIdx] -> [numAggrStates x numClasses].
     */
    public Map<Integer, Matrix> tranProbAggrStateSpace;

    /**
     * Steady-state expected Markov reward per reward name, E[r] = (1/T) integral r(X(t)) dt.
     * Populated only when the model defines rewards via model.setReward.
     */
    public Map<String, Double> avgReward;

    /**
     * Single-run transient reward trajectory r(X(t)) per reward name.
     */
    public Map<String, double[]> tranReward;

    /**
     * Time points associated with {@link #tranReward}.
     */
    public Matrix rewardTime;

    /**
     * Names of the reward functions, in registration order.
     */
    public java.util.List<String> rewardNames;

    /**
     * Exact joint-state residence-time histogram (populated when
     * {@code LDESOptions.exportStateHistogram} is set). Each row is an aggregated
     * system state in the CTMC {@code stateSpaceAggr} layout (station-major,
     * class-minor, {@code nstations*nclasses} columns).
     */
    public Matrix stateHistogramSpace;

    /**
     * Residence time per state in {@link #stateHistogramSpace} (column vector).
     * Divide by its sum to obtain steady-state state probabilities.
     */
    public Matrix stateHistogramTime;

    /**
     * Integer joint-state trajectory (one aggregated state per row, same layout as
     * {@link #stateHistogramSpace}) sampled along the simulated path, downsampled to
     * a bounded number of points. Used by MATLAB/native Python for getTranReward.
     */
    public Matrix stateTrajectorySpace;

    /**
     * Time points associated with {@link #stateTrajectorySpace} (column vector).
     */
    public Matrix stateTrajectoryTime;

    /**
     * Initializes the respTimeSamples array with empty lists for each (station, class) pair.
     *
     * @param nStations number of stations
     * @param nClasses  number of classes
     */
    @SuppressWarnings("unchecked")
    public void initRespTimeSamples(int nStations, int nClasses) {
        respTimeSamples = new java.util.List[nStations][nClasses];
        for (int i = 0; i < nStations; i++) {
            for (int j = 0; j < nClasses; j++) {
                respTimeSamples[i][j] = new java.util.ArrayList<>();
            }
        }
    }

    /**
     * Constructs an empty LDESResult to allow field population from different sources.
     */
    public LDESResult() {
        this.varianceReductionMethod = "none";
    }

    /**
     * Constructs a LDESResult with the specified performance metrics and state information.
     *
     * @param QN queue lengths [stations x classes]
     * @param UN utilizations [stations x classes]
     * @param RN response times [stations x classes]
     * @param TN throughputs [stations x classes]
     * @param CN visit counts [stations x classes]
     * @param XN arrival rates [stations x classes]
     * @param tranSysState transient system states indexed by time
     * @param tranSync transient synchronization matrix
     * @param sn network structure information
     */
    public LDESResult(Matrix QN, Matrix UN, Matrix RN, Matrix TN, Matrix CN, Matrix XN,
                     Map<Integer, Matrix> tranSysState, Matrix tranSync, NetworkStruct sn) {
        this.QN = QN;
        this.UN = UN;
        this.RN = RN;
        this.TN = TN;
        this.CN = CN;
        this.XN = XN;
        this.tranSysState = tranSysState;
        this.tranSync = tranSync;
        this.sn = sn;
    }
}
