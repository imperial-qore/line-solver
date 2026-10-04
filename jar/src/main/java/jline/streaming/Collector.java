package jline.streaming;

import jline.lang.NetworkStruct;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Collects and aggregates metrics from SSA simulation for streaming to line-est.
 * Supports two modes:
 * - SAMPLED: Push metrics every N simulation events
 * - TIME_WINDOW: Accumulate time-weighted metrics and push averages at window end
 *
 * This collector is called from Solver_ssa.kt's save_log() function during simulation.
 */
public class Collector {

    private static final Logger logger = Logger.getLogger(Collector.class.getName());

    private final StreamingOptions options;
    private final OtlpMetricsClient grpcClient;
    private final HttpMetricsClient httpClient;
    private final NetworkStruct sn;

    // For SAMPLED mode
    private int eventsSinceLastPush = 0;

    // For TIME_WINDOW mode
    private double windowStartTime = 0.0;
    private final Map<String, Double> queueLengthAccum;    // station_class -> sum(qLen * dt)
    private final Map<String, Double> timeAccum;           // station_class -> sum(dt)
    private final Map<String, Double> throughputAccum;     // station_class -> sum(tput * dt)
    private final Map<String, Double> arrivalRateAccum;    // station_class -> sum(arvRate * dt)

    // Tracking for cumulative metrics
    private double lastPushTime = 0.0;
    private boolean initialized = false;

    /**
     * Create a new streaming metrics collector.
     * @param options Streaming configuration options
     * @param sn Network structure for station/class information
     */
    public Collector(StreamingOptions options, NetworkStruct sn) {
        this.options = options;
        this.sn = sn;

        // Initialize the appropriate client based on transport type
        if (options.transport == StreamingOptions.TransportType.HTTP) {
            this.httpClient = new HttpMetricsClient(options.endpoint, options.serviceName, options.failOnConnectionError);
            this.grpcClient = null;
        } else {
            this.grpcClient = new OtlpMetricsClient(options.endpoint, options.serviceName, options.failOnConnectionError);
            this.httpClient = null;
        }

        this.queueLengthAccum = new HashMap<String, Double>();
        this.timeAccum = new HashMap<String, Double>();
        this.throughputAccum = new HashMap<String, Double>();
        this.arrivalRateAccum = new HashMap<String, Double>();

        logger.log(Level.INFO, "Collector initialized: mode={0}, transport={1}, endpoint={2}",
                new Object[] { options.mode, options.transport, options.endpoint });
    }

    /**
     * Send metrics using the appropriate transport client.
     * @param metrics List of metrics to send
     * @return true if successful
     */
    private boolean sendMetrics(List<SSAMetricPoint> metrics) {
        if (httpClient != null) {
            return httpClient.sendMetrics(metrics);
        } else if (grpcClient != null) {
            return grpcClient.sendMetrics(metrics);
        }
        return false;
    }

    /**
     * Shutdown the appropriate transport client.
     */
    private void shutdownClient() {
        if (httpClient != null) {
            httpClient.shutdown();
        }
        if (grpcClient != null) {
            grpcClient.shutdown();
        }
    }

    /**
     * Record state observation from simulation loop.
     * Called from save_log() in Solver_ssa.kt.
     *
     * @param simulationTime Current simulation time
     * @param dt Time delta since last event
     * @param nirState Marginal state matrix (queue lengths per station/class as flat vector)
     * @param depRates Departure rates matrix [nclasses x nstateful]
     * @param arvRates Arrival rates matrix [nclasses x nstateful]
     */
    public void recordState(double simulationTime, double dt,
                            Matrix nirState, Matrix depRates, Matrix arvRates) {

        if (!initialized) {
            windowStartTime = simulationTime;
            initialized = true;
        }

        if (options.mode == StreamingOptions.StreamMode.SAMPLED) {
            eventsSinceLastPush++;
            if (eventsSinceLastPush >= options.sampleFrequency) {
                pushCurrentMetrics(simulationTime, nirState, depRates, arvRates);
                eventsSinceLastPush = 0;
            }
        } else {
            // TIME_WINDOW mode
            accumulateForWindow(simulationTime, dt, nirState, depRates, arvRates);
            if (simulationTime - windowStartTime >= options.timeWindowSeconds) {
                pushWindowAverages(simulationTime);
                windowStartTime = simulationTime;
                clearAccumulators();
            }
        }
    }

    /**
     * Push current instantaneous metrics (for SAMPLED mode).
     */
    private void pushCurrentMetrics(double time, Matrix nirState,
                                    Matrix depRates, Matrix arvRates) {
        List<SSAMetricPoint> metrics = new ArrayList<SSAMetricPoint>();
        long timestampNanos = System.currentTimeMillis() * 1_000_000L;

        int nstations = sn.nstations;
        int nclasses = sn.nclasses;

        for (int ist = 0; ist < nstations; ist++) {
            for (int k = 0; k < nclasses; k++) {
                String stationName = getStationName(ist);
                String className = getClassName(k);

                // Queue Length
                if (options.includeQueueLength && nirState != null) {
                    int idx = ist * nclasses + k;
                    if (idx < nirState.length()) {
                        double qLen = nirState.get(idx);
                        metrics.add(SSAMetricPoint.queueLength(ist, k, stationName, className, qLen, timestampNanos));
                    }
                }

                // Get stateful index for this station
                int isf = -1;
                if (sn.stationToStateful != null && ist < sn.stationToStateful.getNumCols()) {
                    isf = (int) sn.stationToStateful.get(ist);
                }

                // Throughput (departure rates)
                if (options.includeThroughput && depRates != null && isf >= 0) {
                    if (k < depRates.getNumRows() && isf < depRates.getNumCols()) {
                        double tput = depRates.get(k, isf);
                        metrics.add(SSAMetricPoint.throughput(ist, k, stationName, className, tput, timestampNanos));
                    }
                }

                // Arrival Rate
                if (options.includeArrivalRate && arvRates != null && isf >= 0) {
                    if (k < arvRates.getNumRows() && isf < arvRates.getNumCols()) {
                        double arvRate = arvRates.get(k, isf);
                        metrics.add(SSAMetricPoint.arrivalRate(ist, k, stationName, className, arvRate, timestampNanos));
                    }
                }

                // Response Time (Little's Law: R = Q / T)
                if (options.includeResponseTime && nirState != null && depRates != null && isf >= 0) {
                    int idx = ist * nclasses + k;
                    if (idx < nirState.length() && k < depRates.getNumRows() && isf < depRates.getNumCols()) {
                        double qLen = nirState.get(idx);
                        double tput = depRates.get(k, isf);
                        if (tput > 1e-10) {
                            double respTime = qLen / tput;
                            metrics.add(SSAMetricPoint.responseTime(ist, k, stationName, className, respTime, timestampNanos));
                        }
                    }
                }

                // Utilization (U = lambda * S / c)
                if (options.includeUtilization && arvRates != null && isf >= 0) {
                    if (k < arvRates.getNumRows() && isf < arvRates.getNumCols()) {
                        double arvRate = arvRates.get(k, isf);
                        // Get service rate and number of servers from sn
                        double util = computeUtilization(ist, k, arvRate);
                        if (!Double.isNaN(util)) {
                            metrics.add(SSAMetricPoint.utilization(ist, k, stationName, className, util, timestampNanos));
                        }
                    }
                }
            }
        }

        if (!metrics.isEmpty()) {
            boolean success = sendMetrics(metrics);
            if (!success) {
                logger.log(Level.WARNING, "Failed to push {0} metrics at time {1}",
                        new Object[] { metrics.size(), time });
            }
        }

        lastPushTime = time;
    }

    /**
     * Accumulate metrics for time-window averaging.
     */
    private void accumulateForWindow(double time, double dt, Matrix nirState,
                                     Matrix depRates, Matrix arvRates) {
        int nstations = sn.nstations;
        int nclasses = sn.nclasses;

        for (int ist = 0; ist < nstations; ist++) {
            for (int k = 0; k < nclasses; k++) {
                String key = ist + "_" + k;

                // Initialize accumulators if needed
                if (!queueLengthAccum.containsKey(key)) {
                    queueLengthAccum.put(key, 0.0);
                    timeAccum.put(key, 0.0);
                    throughputAccum.put(key, 0.0);
                    arrivalRateAccum.put(key, 0.0);
                }

                // Accumulate time
                double prevTime = timeAccum.get(key);
                timeAccum.put(key, prevTime + dt);

                // Queue length (time-weighted)
                if (nirState != null) {
                    int idx = ist * nclasses + k;
                    if (idx < nirState.length()) {
                        double qLen = nirState.get(idx);
                        double prevQLen = queueLengthAccum.get(key);
                        queueLengthAccum.put(key, prevQLen + qLen * dt);
                    }
                }

                // Get stateful index
                int isf = -1;
                if (sn.stationToStateful != null && ist < sn.stationToStateful.getNumCols()) {
                    isf = (int) sn.stationToStateful.get(ist);
                }

                // Throughput (time-weighted)
                if (depRates != null && isf >= 0) {
                    if (k < depRates.getNumRows() && isf < depRates.getNumCols()) {
                        double tput = depRates.get(k, isf);
                        double prevTput = throughputAccum.get(key);
                        throughputAccum.put(key, prevTput + tput * dt);
                    }
                }

                // Arrival rate (time-weighted)
                if (arvRates != null && isf >= 0) {
                    if (k < arvRates.getNumRows() && isf < arvRates.getNumCols()) {
                        double arvRate = arvRates.get(k, isf);
                        double prevArvRate = arrivalRateAccum.get(key);
                        arrivalRateAccum.put(key, prevArvRate + arvRate * dt);
                    }
                }
            }
        }
    }

    /**
     * Push time-weighted average metrics for the window.
     */
    private void pushWindowAverages(double time) {
        List<SSAMetricPoint> metrics = new ArrayList<SSAMetricPoint>();
        long timestampNanos = System.currentTimeMillis() * 1_000_000L;

        int nstations = sn.nstations;
        int nclasses = sn.nclasses;

        for (int ist = 0; ist < nstations; ist++) {
            for (int k = 0; k < nclasses; k++) {
                String key = ist + "_" + k;
                String stationName = getStationName(ist);
                String className = getClassName(k);

                Double totalTime = timeAccum.get(key);
                if (totalTime == null || totalTime <= 0) {
                    continue;
                }

                // Average queue length
                if (options.includeQueueLength) {
                    Double sumQLen = queueLengthAccum.get(key);
                    if (sumQLen != null) {
                        double avgQLen = sumQLen / totalTime;
                        metrics.add(SSAMetricPoint.queueLength(ist, k, stationName, className, avgQLen, timestampNanos));
                    }
                }

                // Average throughput
                if (options.includeThroughput) {
                    Double sumTput = throughputAccum.get(key);
                    if (sumTput != null) {
                        double avgTput = sumTput / totalTime;
                        metrics.add(SSAMetricPoint.throughput(ist, k, stationName, className, avgTput, timestampNanos));
                    }
                }

                // Average arrival rate
                if (options.includeArrivalRate) {
                    Double sumArvRate = arrivalRateAccum.get(key);
                    if (sumArvRate != null) {
                        double avgArvRate = sumArvRate / totalTime;
                        metrics.add(SSAMetricPoint.arrivalRate(ist, k, stationName, className, avgArvRate, timestampNanos));
                    }
                }

                // Response time from Little's Law
                if (options.includeResponseTime) {
                    Double sumQLen = queueLengthAccum.get(key);
                    Double sumTput = throughputAccum.get(key);
                    if (sumQLen != null && sumTput != null && sumTput > 1e-10) {
                        double avgQLen = sumQLen / totalTime;
                        double avgTput = sumTput / totalTime;
                        if (avgTput > 1e-10) {
                            double respTime = avgQLen / avgTput;
                            metrics.add(SSAMetricPoint.responseTime(ist, k, stationName, className, respTime, timestampNanos));
                        }
                    }
                }

                // Utilization
                if (options.includeUtilization) {
                    Double sumArvRate = arrivalRateAccum.get(key);
                    if (sumArvRate != null) {
                        double avgArvRate = sumArvRate / totalTime;
                        double util = computeUtilization(ist, k, avgArvRate);
                        if (!Double.isNaN(util)) {
                            metrics.add(SSAMetricPoint.utilization(ist, k, stationName, className, util, timestampNanos));
                        }
                    }
                }
            }
        }

        if (!metrics.isEmpty()) {
            boolean success = sendMetrics(metrics);
            if (!success) {
                logger.log(Level.WARNING, "Failed to push {0} window-averaged metrics at time {1}",
                        new Object[] { metrics.size(), time });
            }
        }

        lastPushTime = time;
    }

    /**
     * Clear window accumulators.
     */
    private void clearAccumulators() {
        queueLengthAccum.clear();
        timeAccum.clear();
        throughputAccum.clear();
        arrivalRateAccum.clear();
    }

    /**
     * Compute utilization for a station/class.
     * U = lambda * S / c (arrival rate * service time / servers)
     */
    private double computeUtilization(int station, int jobClass, double arrivalRate) {
        try {
            // Get number of servers
            int nservers = 1;
            if (sn.nservers != null && station < sn.nservers.length()) {
                nservers = Math.max(1, (int) sn.nservers.get(station));
            }

            // Get service rate (mu = 1/S)
            double serviceRate = 1.0;
            if (sn.rates != null && station < sn.rates.getNumRows() && jobClass < sn.rates.getNumCols()) {
                double mu = sn.rates.get(station, jobClass);
                if (mu > 0 && !Double.isInfinite(mu)) {
                    serviceRate = mu;
                }
            }

            // U = lambda / (c * mu)
            double util = arrivalRate / (nservers * serviceRate);
            return Math.min(1.0, Math.max(0.0, util));  // Clamp to [0, 1]
        } catch (Exception e) {
            return Double.NaN;
        }
    }

    /**
     * Get station name from network structure.
     */
    private String getStationName(int stationIndex) {
        try {
            if (sn.nodenames != null && stationIndex < sn.nodenames.size()) {
                return sn.nodenames.get(stationIndex);
            }
        } catch (Exception e) {
            // Ignore
        }
        return null;
    }

    /**
     * Get class name from network structure.
     */
    private String getClassName(int classIndex) {
        try {
            if (sn.classnames != null && classIndex < sn.classnames.size()) {
                return sn.classnames.get(classIndex);
            }
        } catch (Exception e) {
            // Ignore
        }
        return null;
    }

    /**
     * Flush any pending metrics (called at end of simulation).
     */
    public void flush(double finalTime) {
        if (options.mode == StreamingOptions.StreamMode.TIME_WINDOW) {
            // Push any remaining accumulated metrics
            if (!timeAccum.isEmpty()) {
                pushWindowAverages(finalTime);
            }
        }
    }

    /**
     * Shutdown the collector and release resources.
     */
    public void shutdown() {
        logger.log(Level.INFO, "Shutting down Collector");
        shutdownClient();
    }

    /**
     * Get the streaming options.
     */
    public StreamingOptions getOptions() {
        return options;
    }

    /**
     * Get the last push time.
     */
    public double getLastPushTime() {
        return lastPushTime;
    }
}
