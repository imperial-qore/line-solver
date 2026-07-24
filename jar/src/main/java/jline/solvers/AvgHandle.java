package jline.solvers;

import jline.lang.JobClass;
import jline.lang.Metric;
import jline.lang.nodes.Station;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Handle for managing performance metrics organized by station and job class.
 *
 * <p>AvgHandle provides a structured way to store and access performance metrics
 * (such as throughput, utilization, queue length, response time) that are computed
 * by LINE solvers. The metrics are organized in a two-level hierarchy: first by
 * station (queue, delay, etc.) and then by job class.
 *
 * <p>This class acts as a container that maps each (station, job class) pair to
 * its corresponding performance metric. It provides methods for storing, retrieving,
 * and managing these metrics efficiently.
 *
 * <p>Example usage:
 * <pre>
 * AvgHandle handle = new AvgHandle();
 * handle.put(queue1, classA, new Metric(MetricType.Throughput, 2.5));
 * Metric throughput = handle.get(queue1, classA);
 * </pre>
 *
 * @see Metric
 * @see Station
 * @see JobClass
 * @see NetworkAvgTable
 */
public class AvgHandle {

    /**
     * Two-level map: Station -> JobClass -> Metric
     */
    private final Map<Station, Map<JobClass, Metric>> data;

    /**
     * Constructs an empty AvgHandle.
     */
    public AvgHandle() {
        this.data = new HashMap<>();
    }

    /**
     * Helper method to find a Station in the data map by its station index.
     *
     * @param stationIdx the station index to find
     * @return the Station object with the matching index, or null if not found
     */
    private Station findStationByIndex(int stationIdx) {
        for (Station station : data.keySet()) {
            if (station != null && station.getStationIdx() == stationIdx) {
                return station;
            }
        }
        return null;
    }

    /**
     * Helper method to find a JobClass in a station's metrics by its job class index.
     *
     * @param stationMetrics the map of job class metrics for a station
     * @param jobClassIdx the job class index to find
     * @return the JobClass object with the matching index, or null if not found
     */
    private JobClass findJobClassByIndex(Map<JobClass, Metric> stationMetrics, int jobClassIdx) {
        for (JobClass jobClass : stationMetrics.keySet()) {
            if (jobClass.getIndex() == jobClassIdx) {
                return jobClass;
            }
        }
        return null;
    }

    /**
     * Retrieves the metric for a specific station and job class.
     *
     * @param station  the station to query
     * @param jobClass the job class to query
     * @return the metric for the given station and job class, or null if not found
     */
    public Metric get(Station station, JobClass jobClass) {
        // Handle null station (system-wide metrics)
        if (station == null) {
            Map<JobClass, Metric> stationMetrics = data.get(null);
            if (stationMetrics != null) {
                JobClass foundJobClass = findJobClassByIndex(stationMetrics, jobClass.getIndex());
                if (foundJobClass != null) {
                    return stationMetrics.get(foundJobClass);
                }
            }
            return null;
        }

        Station foundStation = findStationByIndex(station.getStationIdx());
        if (foundStation != null) {
            Map<JobClass, Metric> stationMetrics = data.get(foundStation);
            if (stationMetrics != null) {
                JobClass foundJobClass = findJobClassByIndex(stationMetrics, jobClass.getIndex());
                if (foundJobClass != null) {
                    return stationMetrics.get(foundJobClass);
                }
            }
        }
        return null;
    }

    /**
     * Retrieves all metrics for a specific station across all job classes.
     *
     * @param station the station to query
     * @return a map from job class to metric for the given station, empty if station not found
     */
    public Map<JobClass, Metric> get(Station station) {
        // Handle null station (system-wide metrics)
        if (station == null) {
            return data.getOrDefault(null, new HashMap<>());
        }

        Station foundStation = findStationByIndex(station.getStationIdx());
        if (foundStation != null) {
            return data.getOrDefault(foundStation, new HashMap<>());
        }
        return new HashMap<>();
    }

    /**
     * Retrieves the complete data structure containing all metrics.
     *
     * @return the full nested map: Station -> JobClass -> Metric
     */
    public Map<Station, Map<JobClass, Metric>> get() {
        return data;
    }

    /**
     * Checks if a metric exists for the specified station and job class.
     *
     * @param station  the station to check
     * @param jobClass the job class to check
     * @return true if a metric exists for the given station and job class, false otherwise
     */
    public boolean hasMetric(Station station, JobClass jobClass) {
        // Handle null station (system-wide metrics)
        if (station == null) {
            Map<JobClass, Metric> stationMetrics = data.get(null);
            if (stationMetrics != null) {
                for (Map.Entry<JobClass, Metric> classEntry : stationMetrics.entrySet()) {
                    if (classEntry.getKey().getIndex() == jobClass.getIndex()) {
                        return true;
                    }
                }
            }
            return false;
        }

        // Use index-based lookup for both station and jobClass
        for (Map.Entry<Station, Map<JobClass, Metric>> stationEntry : data.entrySet()) {
            Station entryStation = stationEntry.getKey();
            if (entryStation != null && entryStation.getStationIdx() == station.getStationIdx()) {
                Map<JobClass, Metric> stationMetrics = stationEntry.getValue();
                for (Map.Entry<JobClass, Metric> classEntry : stationMetrics.entrySet()) {
                    if (classEntry.getKey().getIndex() == jobClass.getIndex()) {
                        return true;
                    }
                }
                return false;
            }
        }
        return false;
    }

    /**
     * Checks if this handle contains no metrics.
     *
     * @return true if no metrics are stored, false otherwise
     */
    public boolean isEmpty() {
        return data.isEmpty();
    }

    /**
     * Returns the set of all stations that have metrics stored.
     *
     * @return a set containing all stations with stored metrics
     */
    public Set<Station> keySet() {
        return data.keySet();
    }

    /**
     * Stores a metric for a specific station and job class.
     *
     * @param station  the station associated with the metric
     * @param jobClass the job class associated with the metric
     * @param metric   the metric value to store
     */
    public void put(Station station, JobClass jobClass, Metric metric) {
        // Handle null station (system-wide metrics) - store directly with null key
        Station targetStation;
        if (station == null) {
            targetStation = null;
        } else {
            // Find existing station by index, or use the provided station if not found
            targetStation = findStationByIndex(station.getStationIdx());
            if (targetStation == null) {
                targetStation = station;
            }
        }

        // Get or create the station's metrics map
        Map<JobClass, Metric> stationMetrics = data.computeIfAbsent(targetStation, k -> new HashMap<>());

        // Find existing job class by index, or use the provided job class if not found
        JobClass targetJobClass = findJobClassByIndex(stationMetrics, jobClass.getIndex());
        if (targetJobClass == null) {
            targetJobClass = jobClass;
        }
        
        stationMetrics.put(targetJobClass, metric);
    }

    /**
     * Removes the metric for a specific station and job class.
     *
     * <p>If this removal leaves the station with no job classes, the station
     * entry is also removed from the data structure.
     *
     * @param station  the station to remove from
     * @param jobClass the job class to remove
     */
    public void remove(Station station, JobClass jobClass) {
        // Handle null station (system-wide metrics)
        Station foundStation;
        if (station == null) {
            foundStation = null;
        } else {
            foundStation = findStationByIndex(station.getStationIdx());
        }

        if (foundStation != null || station == null) {
            Map<JobClass, Metric> jobClassMetricMap = data.get(foundStation);
            if (jobClassMetricMap != null) {
                JobClass foundJobClass = findJobClassByIndex(jobClassMetricMap, jobClass.getIndex());
                if (foundJobClass != null) {
                    jobClassMetricMap.remove(foundJobClass);
                    if (jobClassMetricMap.isEmpty()) {
                        data.remove(foundStation);
                    }
                }
            }
        }
    }

}
