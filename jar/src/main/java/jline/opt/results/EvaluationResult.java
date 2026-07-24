package jline.opt.results;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Result from evaluating a LINE network model via SolverAUTO, holding the
 * per-(station, class) performance metrics the optimizer reads. Mirrors the
 * native-Python {@code line_solver.opt.results.EvaluationResult}: per-station,
 * per-class metrics are stored as nested station -&gt; class -&gt; value maps
 * (the Java analogue of Python's tuple-keyed dicts), utilizations per station,
 * and system (end-to-end) metrics per chain/class.
 */
public class EvaluationResult {

    public boolean feasible = true;
    // metric kind -> station -> class -> value
    public final Map<String, Map<String, Double>> responseTimes = new LinkedHashMap<String, Map<String, Double>>();
    public final Map<String, Map<String, Double>> throughputs = new LinkedHashMap<String, Map<String, Double>>();
    public final Map<String, Map<String, Double>> queueLengths = new LinkedHashMap<String, Map<String, Double>>();
    public final Map<String, Double> utilizations = new LinkedHashMap<String, Double>();
    public final Map<String, Double> systemResponseTimes = new LinkedHashMap<String, Double>();
    public final Map<String, Double> systemThroughputs = new LinkedHashMap<String, Double>();
    public double solveTime = 0.0;
    public String solverUsed = "";
    /** Analytic sensitivities when product-form; null otherwise. */
    public SensitivityData sensitivities = null;

    private static void put(Map<String, Map<String, Double>> m, String station,
                            String jobclass, double value) {
        Map<String, Double> row = m.get(station);
        if (row == null) {
            row = new LinkedHashMap<String, Double>();
            m.put(station, row);
        }
        row.put(jobclass, value);
    }

    public void setResponseTime(String station, String jobclass, double v) {
        put(responseTimes, station, jobclass, v);
    }

    public void setThroughput(String station, String jobclass, double v) {
        put(throughputs, station, jobclass, v);
    }

    public void setQueueLength(String station, String jobclass, double v) {
        put(queueLengths, station, jobclass, v);
    }

    public double getResponseTime(String station, String jobclass) {
        Map<String, Double> row = responseTimes.get(station);
        if (row == null || !row.containsKey(jobclass)) {
            return Double.POSITIVE_INFINITY;
        }
        return row.get(jobclass);
    }

    /** Aggregate (mean) response time across classes at a station. */
    public double getResponseTime(String station) {
        Map<String, Double> row = responseTimes.get(station);
        if (row == null || row.isEmpty()) {
            return Double.POSITIVE_INFINITY;
        }
        double total = 0.0;
        for (double v : row.values()) {
            total += v;
        }
        return total / row.size();
    }

    public double getThroughput(String station, String jobclass) {
        Map<String, Double> row = throughputs.get(station);
        if (row == null || !row.containsKey(jobclass)) {
            return 0.0;
        }
        return row.get(jobclass);
    }

    public double getThroughput(String station) {
        Map<String, Double> row = throughputs.get(station);
        if (row == null) {
            return 0.0;
        }
        double total = 0.0;
        for (double v : row.values()) {
            total += v;
        }
        return total;
    }

    public double getUtilization(String station) {
        Double v = utilizations.get(station);
        return v == null ? 0.0 : v;
    }

    public double getQueueLength(String station, String jobclass) {
        Map<String, Double> row = queueLengths.get(station);
        if (row == null || !row.containsKey(jobclass)) {
            return 0.0;
        }
        return row.get(jobclass);
    }

    public double getQueueLength(String station) {
        Map<String, Double> row = queueLengths.get(station);
        if (row == null) {
            return 0.0;
        }
        double total = 0.0;
        for (double v : row.values()) {
            total += v;
        }
        return total;
    }

    public double getSystemResponseTime(String jobclass) {
        Double v = systemResponseTimes.get(jobclass);
        return v == null ? Double.POSITIVE_INFINITY : v;
    }

    /**
     * Throughput-weighted mean end-to-end response time across chains, matching
     * total-jobs / total-throughput by Little's law.
     */
    public double getSystemResponseTime() {
        if (systemResponseTimes.isEmpty()) {
            return Double.POSITIVE_INFINITY;
        }
        double totalTput = 0.0;
        double weighted = 0.0;
        for (Map.Entry<String, Double> e : systemResponseTimes.entrySet()) {
            double tput = systemThroughputs.containsKey(e.getKey())
                    ? systemThroughputs.get(e.getKey()) : 0.0;
            weighted += e.getValue() * tput;
            totalTput += tput;
        }
        if (totalTput > 0) {
            return weighted / totalTput;
        }
        double sum = 0.0;
        for (double v : systemResponseTimes.values()) {
            sum += v;
        }
        return sum / systemResponseTimes.size();
    }

    public double getSystemThroughput(String jobclass) {
        Double v = systemThroughputs.get(jobclass);
        return v == null ? 0.0 : v;
    }

    public double getSystemThroughput() {
        double sum = 0.0;
        for (double v : systemThroughputs.values()) {
            sum += v;
        }
        return sum;
    }
}
