package jline.streaming;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * Internal representation of a metric point from SSA simulation.
 * Used to collect metrics before sending to OTLP receiver.
 *
 * Metric names are aligned with line-est MetricType enum:
 * - queue_length -> QLen
 * - utilization -> Util
 * - throughput -> Tput
 * - response_time -> RespT
 * - arrival_rate -> ArvR
 */
public class SSAMetricPoint implements Serializable {

    /** Metric name (e.g., "queue_length", "throughput") */
    public final String metricName;

    /** Metric value */
    public final double value;

    /** Timestamp in nanoseconds since epoch */
    public final long timestampNanos;

    /** Labels for this metric (station, class, etc.) */
    public final Map<String, String> labels;

    /**
     * Create a new metric point.
     * @param metricName Name of the metric
     * @param value Metric value
     * @param timestampNanos Timestamp in nanoseconds
     * @param labels Additional labels (can be null)
     */
    public SSAMetricPoint(String metricName, double value, long timestampNanos, Map<String, String> labels) {
        this.metricName = metricName;
        this.value = value;
        this.timestampNanos = timestampNanos;
        this.labels = labels != null ? labels : new HashMap<String, String>();
    }

    /**
     * Create a queue length metric point.
     * @param station Station index
     * @param jobClass Job class index
     * @param value Queue length value
     * @param timestampNanos Timestamp in nanoseconds
     * @return New SSAMetricPoint
     */
    public static SSAMetricPoint queueLength(int station, int jobClass, double value, long timestampNanos) {
        Map<String, String> labels = new HashMap<String, String>();
        labels.put("station", String.valueOf(station));
        labels.put("class", String.valueOf(jobClass));
        return new SSAMetricPoint("queue_length", value, timestampNanos, labels);
    }

    /**
     * Create a queue length metric point with station/class names.
     * @param station Station index
     * @param jobClass Job class index
     * @param stationName Station name
     * @param className Class name
     * @param value Queue length value
     * @param timestampNanos Timestamp in nanoseconds
     * @return New SSAMetricPoint
     */
    public static SSAMetricPoint queueLength(int station, int jobClass, String stationName,
                                              String className, double value, long timestampNanos) {
        Map<String, String> labels = new HashMap<String, String>();
        labels.put("station", String.valueOf(station));
        labels.put("class", String.valueOf(jobClass));
        if (stationName != null) {
            labels.put("station_name", stationName);
        }
        if (className != null) {
            labels.put("class_name", className);
        }
        return new SSAMetricPoint("queue_length", value, timestampNanos, labels);
    }

    /**
     * Create a utilization metric point.
     * @param station Station index
     * @param jobClass Job class index
     * @param value Utilization value (0-1)
     * @param timestampNanos Timestamp in nanoseconds
     * @return New SSAMetricPoint
     */
    public static SSAMetricPoint utilization(int station, int jobClass, double value, long timestampNanos) {
        Map<String, String> labels = new HashMap<String, String>();
        labels.put("station", String.valueOf(station));
        labels.put("class", String.valueOf(jobClass));
        return new SSAMetricPoint("utilization", value, timestampNanos, labels);
    }

    /**
     * Create a utilization metric point with station/class names.
     * @param station Station index
     * @param jobClass Job class index
     * @param stationName Station name
     * @param className Class name
     * @param value Utilization value (0-1)
     * @param timestampNanos Timestamp in nanoseconds
     * @return New SSAMetricPoint
     */
    public static SSAMetricPoint utilization(int station, int jobClass, String stationName,
                                              String className, double value, long timestampNanos) {
        Map<String, String> labels = new HashMap<String, String>();
        labels.put("station", String.valueOf(station));
        labels.put("class", String.valueOf(jobClass));
        if (stationName != null) {
            labels.put("station_name", stationName);
        }
        if (className != null) {
            labels.put("class_name", className);
        }
        return new SSAMetricPoint("utilization", value, timestampNanos, labels);
    }

    /**
     * Create a throughput metric point.
     * @param station Station index
     * @param jobClass Job class index
     * @param value Throughput value (departures/time)
     * @param timestampNanos Timestamp in nanoseconds
     * @return New SSAMetricPoint
     */
    public static SSAMetricPoint throughput(int station, int jobClass, double value, long timestampNanos) {
        Map<String, String> labels = new HashMap<String, String>();
        labels.put("station", String.valueOf(station));
        labels.put("class", String.valueOf(jobClass));
        return new SSAMetricPoint("throughput", value, timestampNanos, labels);
    }

    /**
     * Create a throughput metric point with station/class names.
     * @param station Station index
     * @param jobClass Job class index
     * @param stationName Station name
     * @param className Class name
     * @param value Throughput value (departures/time)
     * @param timestampNanos Timestamp in nanoseconds
     * @return New SSAMetricPoint
     */
    public static SSAMetricPoint throughput(int station, int jobClass, String stationName,
                                             String className, double value, long timestampNanos) {
        Map<String, String> labels = new HashMap<String, String>();
        labels.put("station", String.valueOf(station));
        labels.put("class", String.valueOf(jobClass));
        if (stationName != null) {
            labels.put("station_name", stationName);
        }
        if (className != null) {
            labels.put("class_name", className);
        }
        return new SSAMetricPoint("throughput", value, timestampNanos, labels);
    }

    /**
     * Create a response time metric point.
     * @param station Station index
     * @param jobClass Job class index
     * @param value Response time value
     * @param timestampNanos Timestamp in nanoseconds
     * @return New SSAMetricPoint
     */
    public static SSAMetricPoint responseTime(int station, int jobClass, double value, long timestampNanos) {
        Map<String, String> labels = new HashMap<String, String>();
        labels.put("station", String.valueOf(station));
        labels.put("class", String.valueOf(jobClass));
        return new SSAMetricPoint("response_time", value, timestampNanos, labels);
    }

    /**
     * Create a response time metric point with station/class names.
     * @param station Station index
     * @param jobClass Job class index
     * @param stationName Station name
     * @param className Class name
     * @param value Response time value
     * @param timestampNanos Timestamp in nanoseconds
     * @return New SSAMetricPoint
     */
    public static SSAMetricPoint responseTime(int station, int jobClass, String stationName,
                                               String className, double value, long timestampNanos) {
        Map<String, String> labels = new HashMap<String, String>();
        labels.put("station", String.valueOf(station));
        labels.put("class", String.valueOf(jobClass));
        if (stationName != null) {
            labels.put("station_name", stationName);
        }
        if (className != null) {
            labels.put("class_name", className);
        }
        return new SSAMetricPoint("response_time", value, timestampNanos, labels);
    }

    /**
     * Create an arrival rate metric point.
     * @param station Station index
     * @param jobClass Job class index
     * @param value Arrival rate value (arrivals/time)
     * @param timestampNanos Timestamp in nanoseconds
     * @return New SSAMetricPoint
     */
    public static SSAMetricPoint arrivalRate(int station, int jobClass, double value, long timestampNanos) {
        Map<String, String> labels = new HashMap<String, String>();
        labels.put("station", String.valueOf(station));
        labels.put("class", String.valueOf(jobClass));
        return new SSAMetricPoint("arrival_rate", value, timestampNanos, labels);
    }

    /**
     * Create an arrival rate metric point with station/class names.
     * @param station Station index
     * @param jobClass Job class index
     * @param stationName Station name
     * @param className Class name
     * @param value Arrival rate value (arrivals/time)
     * @param timestampNanos Timestamp in nanoseconds
     * @return New SSAMetricPoint
     */
    public static SSAMetricPoint arrivalRate(int station, int jobClass, String stationName,
                                              String className, double value, long timestampNanos) {
        Map<String, String> labels = new HashMap<String, String>();
        labels.put("station", String.valueOf(station));
        labels.put("class", String.valueOf(jobClass));
        if (stationName != null) {
            labels.put("station_name", stationName);
        }
        if (className != null) {
            labels.put("class_name", className);
        }
        return new SSAMetricPoint("arrival_rate", value, timestampNanos, labels);
    }

    @Override
    public String toString() {
        return "SSAMetricPoint{" +
                "metricName='" + metricName + '\'' +
                ", value=" + value +
                ", timestampNanos=" + timestampNanos +
                ", labels=" + labels +
                '}';
    }
}
