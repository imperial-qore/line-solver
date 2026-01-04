package jline.lang.nodeparam;

import jline.lang.NodeParam;

import java.util.ArrayList;
import java.util.List;

/**
 * Parameter container for logger nodes in queueing networks.
 * 
 * <p>This class configures logger nodes that record job passage events, service times,
 * and other performance metrics to files. Logger nodes are essential for trace generation,
 * performance monitoring, and validation of analytical results against simulation.</p>
 * 
 * <p>Logging capabilities include:
 * <ul>
 *   <li>Job event logging with timestamps and job identifiers</li>
 *   <li>Service time measurements for same-class and cross-class interactions</li>
 *   <li>Configurable output formats and file destinations</li>
 *   <li>Real-time performance monitoring</li>
 *   <li>Trace generation for replay and analysis</li>
 * </ul>
 * </p>
 * 
 * @see jline.lang.nodes.Logger
 * @since 1.0
 */
public class LoggerNodeParam extends NodeParam {
    
    /** List of output file names for different logging streams */
    public List<String> fileName;
    
    /** Base file path for logger output files */
    public String filePath;
    
    /** Whether to log simulation start time as reference */
    public boolean startTime;
    
    /** Whether to include timestamps in log entries */
    public boolean timestamp;
    
    /** Whether to include logger node name in log entries */
    public boolean loggerName;
    
    /** Whether to include unique job identifiers in log entries */
    public boolean jobID;
    
    /** Whether to include job class information in log entries */
    public boolean jobClass;
    
    /** Whether to measure and log service times for same-class jobs */
    public boolean timeSameClass;
    
    /** Whether to measure and log service times for any-class interactions */
    public boolean timeAnyClass;
    
    /**
     * Constructs a new LoggerNodeParam with default logging configuration.
     * Initializes fileName list with one empty entry.
     */
    public LoggerNodeParam() {
        this.fileName = new ArrayList<>();
        this.fileName.add(""); // Initialize with one empty string
    }

    /**
     * Checks if this logger parameter container is empty (no logging configuration set).
     * 
     * @return true if no logging parameters are configured, false otherwise
     */
    @Override
    public boolean isEmpty() {
        return (fileName == null || fileName.isEmpty()) &&
                filePath == null &&
                !startTime &&
                !timestamp &&
                !loggerName &&
                !jobID &&
                !jobClass &&
                !timeSameClass &&
                !timeAnyClass;
    }
}
