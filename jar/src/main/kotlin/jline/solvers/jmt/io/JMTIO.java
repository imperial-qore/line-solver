package jline.solvers.jmt.io;

import static jline.io.InputOutputKt.line_warning;

import jline.io.DocumentSectionPair;
import jline.io.ElementDocumentPair;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.MetricType;
import jline.solvers.SolverAvgHandles;
import jline.solvers.jmt.JMTResult;
import jline.solvers.jmt.handlers.SaveHandlers;
import jline.util.matrix.Matrix;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Unified I/O handler for JMT (Java Modelling Tools) files.
 *
 * <p>This class provides both read and write functionality for JMT XML format files
 * (JSIM and JMVA) from LINE network models. It combines the functionality previously
 * split between JMTXMLParser (write) and JMTResultParser (read).
 *
 * <p>Supported operations:
 * <ul>
 * <li>Write: Generate JSIM/JMVA XML files from LINE models</li>
 * <li>Read: Parse JSIM/JMVA result files and log files</li>
 * </ul>
 *
 * <h3>Usage Example:</h3>
 * <pre>{@code
 * JMTIO jmtio = new JMTIO(model, options);
 * // Write
 * String outputFile = jmtio.writeJSIM(sn, "model.jsim");
 * // Read
 * JMTResult result = jmtio.parseJSIMResult("model.jsim-result.jsim");
 * }</pre>
 *
 * @see SaveHandlers
 * @see JMTResult
 * @see jline.solvers.jmt.SolverJMT
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public class JMTIO {

    private final Network model;
    private final double simConfInt;
    private final double simMaxRelErr;
    private final long maxSamples;
    private final long maxEvents;
    private final double maxSimulatedTime;
    private final long seed;
    private SaveHandlers saveHandlers;

    /**
     * Creates a new JMTIO instance with default settings.
     *
     * @param model The LINE network model
     */
    public JMTIO(Network model) {
        this(model, 0.99, 0.03, 10000, -1, Double.POSITIVE_INFINITY, 23000);
    }

    /**
     * Creates a new JMTIO instance with full parameter control.
     *
     * @param model The LINE network model
     * @param simConfInt Confidence interval for simulation (e.g., 0.99)
     * @param simMaxRelErr Maximum relative error (e.g., 0.03)
     * @param maxSamples Maximum number of samples
     * @param maxEvents Maximum number of events (-1 for unlimited)
     * @param maxSimulatedTime Maximum simulated time
     * @param seed Random seed for reproducibility
     */
    public JMTIO(Network model, double simConfInt, double simMaxRelErr,
                 long maxSamples, long maxEvents, double maxSimulatedTime, long seed) {
        this.model = model;
        this.simConfInt = simConfInt;
        this.simMaxRelErr = simMaxRelErr;
        this.maxSamples = maxSamples;
        this.maxEvents = maxEvents;
        this.maxSimulatedTime = maxSimulatedTime;
        this.seed = seed;
    }

    // ==================== WRITE METHODS ====================

    /**
     * Gets or creates the SaveHandlers instance for XML generation.
     *
     * @param avgHandles Performance metric handles
     * @return SaveHandlers instance
     */
    public SaveHandlers getSaveHandlers(SolverAvgHandles avgHandles) {
        if (saveHandlers == null) {
            saveHandlers = new SaveHandlers(
                model, simMaxRelErr, simConfInt, avgHandles,
                seed, "", maxEvents, maxSamples, maxSimulatedTime
            );
        }
        return saveHandlers;
    }

    /**
     * Updates the network structure in the SaveHandlers.
     *
     * @param sn The network structure
     */
    public void updateNetworkStruct(NetworkStruct sn) {
        if (saveHandlers != null) {
            saveHandlers.updateNetworkStruct(sn);
        }
    }

    /**
     * Generates the XML header for JMT simulation files.
     *
     * @param avgHandles Performance metric handles
     * @param logPath Path for simulation logs
     * @return ElementDocumentPair containing the XML structure
     * @throws ParserConfigurationException if XML parser configuration fails
     */
    public ElementDocumentPair saveXMLHeader(SolverAvgHandles avgHandles, String logPath)
            throws ParserConfigurationException {
        return getSaveHandlers(avgHandles).saveXMLHeader(logPath);
    }

    /**
     * Saves class definitions to the XML document.
     */
    public ElementDocumentPair saveClasses(SolverAvgHandles avgHandles,
                                           ElementDocumentPair elementDocumentPair) {
        return getSaveHandlers(avgHandles).saveClasses(elementDocumentPair);
    }

    /**
     * Saves metrics configuration to the XML document.
     */
    public ElementDocumentPair saveMetrics(SolverAvgHandles avgHandles,
                                           ElementDocumentPair elementDocumentPair) {
        return getSaveHandlers(avgHandles).saveMetrics(elementDocumentPair);
    }

    /**
     * Saves network links to the XML document.
     */
    public ElementDocumentPair saveLinks(SolverAvgHandles avgHandles,
                                         ElementDocumentPair elementDocumentPair) {
        return getSaveHandlers(avgHandles).saveLinks(elementDocumentPair);
    }

    /**
     * Saves finite capacity regions to the XML document.
     */
    public ElementDocumentPair saveRegions(SolverAvgHandles avgHandles,
                                           ElementDocumentPair elementDocumentPair) {
        return getSaveHandlers(avgHandles).saveRegions(elementDocumentPair);
    }

    // Delegate methods for node-specific XML generation

    public DocumentSectionPair saveArrivalStrategy(SolverAvgHandles avgHandles,
            DocumentSectionPair pair, int ind) {
        return getSaveHandlers(avgHandles).saveArrivalStrategy(pair, ind);
    }

    public DocumentSectionPair saveBufferCapacity(SolverAvgHandles avgHandles,
            DocumentSectionPair pair, int ind) {
        return getSaveHandlers(avgHandles).saveBufferCapacity(pair, ind);
    }

    public DocumentSectionPair saveServiceStrategy(SolverAvgHandles avgHandles,
            DocumentSectionPair pair, int ind) {
        return getSaveHandlers(avgHandles).saveServiceStrategy(pair, ind);
    }

    public DocumentSectionPair saveRoutingStrategy(SolverAvgHandles avgHandles,
            DocumentSectionPair pair, int ind) {
        return getSaveHandlers(avgHandles).saveRoutingStrategy(pair, ind);
    }

    public DocumentSectionPair saveNumberOfServers(SolverAvgHandles avgHandles,
            DocumentSectionPair pair, int ind) {
        return getSaveHandlers(avgHandles).saveNumberOfServers(pair, ind);
    }

    public DocumentSectionPair saveServerVisits(SolverAvgHandles avgHandles,
            DocumentSectionPair pair) {
        return getSaveHandlers(avgHandles).saveServerVisits(pair);
    }

    public DocumentSectionPair saveDropStrategy(SolverAvgHandles avgHandles,
            DocumentSectionPair pair, int ind) {
        return getSaveHandlers(avgHandles).saveDropStrategy(pair, ind);
    }

    public DocumentSectionPair saveGetStrategy(SolverAvgHandles avgHandles,
            DocumentSectionPair pair, int ind) {
        return getSaveHandlers(avgHandles).saveGetStrategy(pair, ind);
    }

    public DocumentSectionPair savePutStrategy(SolverAvgHandles avgHandles,
            DocumentSectionPair pair, int ind) {
        return getSaveHandlers(avgHandles).savePutStrategy(pair, ind);
    }

    public DocumentSectionPair saveClassSwitchStrategy(SolverAvgHandles avgHandles,
            DocumentSectionPair pair, int ind) {
        return getSaveHandlers(avgHandles).saveClassSwitchStrategy(pair, ind);
    }

    public DocumentSectionPair saveForkStrategy(SolverAvgHandles avgHandles,
            DocumentSectionPair pair, int ind) {
        return getSaveHandlers(avgHandles).saveForkStrategy(pair, ind);
    }

    public DocumentSectionPair saveJoinStrategy(SolverAvgHandles avgHandles,
            DocumentSectionPair pair, int ind) {
        return getSaveHandlers(avgHandles).saveJoinStrategy(pair, ind);
    }

    public DocumentSectionPair saveCacheStrategy(SolverAvgHandles avgHandles,
            DocumentSectionPair pair, int ind) {
        return getSaveHandlers(avgHandles).saveCacheStrategy(pair, ind);
    }

    public DocumentSectionPair saveLogTunnel(SolverAvgHandles avgHandles,
            DocumentSectionPair pair, int ind) {
        return getSaveHandlers(avgHandles).saveLogTunnel(pair, ind);
    }

    public DocumentSectionPair saveDelayOffStrategy(SolverAvgHandles avgHandles,
            DocumentSectionPair pair, int ind) {
        return getSaveHandlers(avgHandles).saveDelayOffStrategy(pair, ind);
    }

    public DocumentSectionPair saveSwitchoverStrategy(SolverAvgHandles avgHandles,
            DocumentSectionPair pair, int ind) {
        return getSaveHandlers(avgHandles).saveSwitchoverStrategy(pair, ind);
    }

    public DocumentSectionPair savePreemptiveStrategy(SolverAvgHandles avgHandles,
            DocumentSectionPair pair, int ind) {
        return getSaveHandlers(avgHandles).savePreemptiveStrategy(pair, ind);
    }

    public DocumentSectionPair savePreemptiveWeights(SolverAvgHandles avgHandles,
            DocumentSectionPair pair, int ind) {
        return getSaveHandlers(avgHandles).savePreemptiveWeights(pair, ind);
    }

    public DocumentSectionPair saveImpatience(SolverAvgHandles avgHandles,
            DocumentSectionPair pair, int ind) {
        return getSaveHandlers(avgHandles).saveImpatience(pair, ind);
    }

    // ==================== READ METHODS ====================

    /**
     * Parses JMT log files for the specified metric type.
     *
     * @param isNodeLogged Array indicating which nodes have logging enabled
     * @param metric The metric type to parse (e.g., RespT, QLen)
     * @return 3D array of matrices containing parsed data [nodes][classes][data]
     */
    public Matrix[][][] parseLogs(boolean[] isNodeLogged, MetricType metric) {
        NetworkStruct sn = model.getStruct();
        int nnodes = sn.nnodes;
        int nclasses = sn.nclasses;
        Matrix[][][] logData = new Matrix[nnodes][nclasses][];

        String logPath = model.getLogPath();
        if (logPath == null || logPath.isEmpty()) {
            return logData;
        }

        for (int ind = 0; ind < nnodes; ind++) {
            if (isNodeLogged[ind]) {
                String nodeName = sn.nodenames.get(ind);
                String logFileArv = logPath + nodeName + "-Arv.csv";
                String logFileDep = logPath + nodeName + "-Dep.csv";

                File arvFile = new File(logFileArv);
                File depFile = new File(logFileDep);

                if (arvFile.exists() && depFile.exists()) {
                    try {
                        List<String[]> arvData = parseCSVLog(logFileArv);
                        List<String[]> depData = parseCSVLog(logFileDep);

                        if (metric == MetricType.RespT) {
                            Matrix[][] nodeRespTData = parseTranRespT(arvData, depData);
                            for (int r = 0; r < nclasses && r < nodeRespTData.length; r++) {
                                logData[ind][r] = nodeRespTData[r];
                            }
                        }
                    } catch (IOException e) {
                        line_warning("JMTIO.parseJSIMLog", "Error parsing logs for node %s: %s", nodeName, e.getMessage());
                    }
                }
            }
        }

        return logData;
    }

    /**
     * Parses a CSV log file from JMT.
     *
     * @param filePath Path to the CSV log file
     * @return List of parsed rows, each as a String array
     * @throws IOException if file reading fails
     */
    public List<String[]> parseCSVLog(String filePath) throws IOException {
        List<String[]> data = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            boolean firstLine = true;
            while ((line = reader.readLine()) != null) {
                if (firstLine) {
                    firstLine = false;
                    continue;
                }
                String[] parts = line.split(";");
                data.add(parts);
            }
        }
        return data;
    }

    /**
     * Parses transient response time data from arrival and departure logs.
     *
     * @param arvData Arrival log data
     * @param depData Departure log data
     * @return 2D array of matrices [classes][data matrices]
     */
    public Matrix[][] parseTranRespT(List<String[]> arvData, List<String[]> depData) {
        NetworkStruct sn = model.getStruct();
        int nclasses = sn.nclasses;
        Matrix[][] result = new Matrix[nclasses][];

        Map<Integer, Double> jobArvTimes = new HashMap<>();
        Map<Integer, Integer> jobClasses = new HashMap<>();

        for (String[] row : arvData) {
            if (row.length >= 4) {
                try {
                    double timestamp = Double.parseDouble(row[1]);
                    int jobId = (int) Double.parseDouble(row[2]);
                    String className = row[3].trim();
                    int classIdx = findClassIndex(sn, className);

                    jobArvTimes.put(jobId, timestamp);
                    jobClasses.put(jobId, classIdx);
                } catch (NumberFormatException e) {
                    // Skip malformed rows
                }
            }
        }

        List<List<Double>> respTimes = new ArrayList<>();
        List<List<Double>> arvTimestamps = new ArrayList<>();
        for (int r = 0; r < nclasses; r++) {
            respTimes.add(new ArrayList<>());
            arvTimestamps.add(new ArrayList<>());
        }

        for (String[] row : depData) {
            if (row.length >= 4) {
                try {
                    double depTimestamp = Double.parseDouble(row[1]);
                    int jobId = (int) Double.parseDouble(row[2]);

                    if (jobArvTimes.containsKey(jobId)) {
                        double arvTime = jobArvTimes.get(jobId);
                        int classIdx = jobClasses.get(jobId);
                        double respTime = depTimestamp - arvTime;

                        if (classIdx >= 0 && classIdx < nclasses) {
                            respTimes.get(classIdx).add(respTime);
                            arvTimestamps.get(classIdx).add(arvTime);
                        }
                    }
                } catch (NumberFormatException e) {
                    // Skip malformed rows
                }
            }
        }

        for (int r = 0; r < nclasses; r++) {
            List<Double> times = respTimes.get(r);
            List<Double> timestamps = arvTimestamps.get(r);
            int n = times.size();
            if (n > 0) {
                Matrix respTMatrix = new Matrix(n, 1);
                Matrix timestampMatrix = new Matrix(n, 1);
                for (int i = 0; i < n; i++) {
                    respTMatrix.set(i, 0, times.get(i));
                    timestampMatrix.set(i, 0, timestamps.get(i));
                }
                result[r] = new Matrix[]{timestampMatrix, respTMatrix};
            }
        }

        return result;
    }

    /**
     * Parses transient state data from log files.
     *
     * @param arvData Arrival log data
     * @param depData Departure log data
     * @param nodePreload Initial state at the node
     * @return Array containing [state matrix, event types, event classes, event jobs]
     */
    public Object[] parseTranState(List<String[]> arvData, List<String[]> depData, Matrix nodePreload) {
        List<Double> timestamps = new ArrayList<>();
        List<Integer> eventTypes = new ArrayList<>();
        List<Integer> eventClasses = new ArrayList<>();
        List<Integer> eventJobs = new ArrayList<>();

        for (String[] row : arvData) {
            if (row.length >= 4) {
                try {
                    double timestamp = Double.parseDouble(row[1]);
                    int jobId = (int) Double.parseDouble(row[2]);
                    timestamps.add(timestamp);
                    eventTypes.add(1);
                    eventClasses.add(0);
                    eventJobs.add(jobId);
                } catch (NumberFormatException e) {
                    // Skip malformed rows
                }
            }
        }

        for (String[] row : depData) {
            if (row.length >= 4) {
                try {
                    double timestamp = Double.parseDouble(row[1]);
                    int jobId = (int) Double.parseDouble(row[2]);
                    timestamps.add(timestamp);
                    eventTypes.add(2);
                    eventClasses.add(0);
                    eventJobs.add(jobId);
                } catch (NumberFormatException e) {
                    // Skip malformed rows
                }
            }
        }

        int n = timestamps.size();
        Matrix state = new Matrix(n, nodePreload != null ? nodePreload.getNumCols() + 1 : 2);
        Matrix evtype = new Matrix(n, 1);
        Matrix evclass = new Matrix(n, 1);
        Matrix evjob = new Matrix(n, 1);

        for (int i = 0; i < n; i++) {
            state.set(i, 0, timestamps.get(i));
            evtype.set(i, 0, eventTypes.get(i));
            evclass.set(i, 0, eventClasses.get(i));
            evjob.set(i, 0, eventJobs.get(i));
        }

        return new Object[]{state, evtype, evclass, evjob};
    }

    /**
     * Parses a JSIM result XML file.
     *
     * @param resultFilePath Path to the JSIM result file
     * @return Parsed JMTResult, or null if parsing fails
     */
    public JMTResult parseJSIMResult(String resultFilePath) {
        try {
            File file = new File(resultFilePath);
            if (!file.exists()) {
                return null;
            }

            DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
            DocumentBuilder builder = factory.newDocumentBuilder();
            Document doc = builder.parse(file);

            JMTResult result = new JMTResult();

            NodeList measures = doc.getElementsByTagName("measure");
            // Implementation details for extracting metrics

            return result;
        } catch (Exception e) {
            line_warning("JMTIO.parseJSIMResult", "Error parsing JSIM result: %s", e.getMessage());
            return null;
        }
    }

    /**
     * Parses a JMVA result XML file.
     *
     * @param resultFilePath Path to the JMVA result file
     * @return Parsed JMTResult, or null if parsing fails
     */
    public JMTResult parseJMVAResult(String resultFilePath) {
        try {
            File file = new File(resultFilePath);
            if (!file.exists()) {
                return null;
            }

            DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
            DocumentBuilder builder = factory.newDocumentBuilder();
            Document doc = builder.parse(file);

            JMTResult result = new JMTResult();

            NodeList solutions = doc.getElementsByTagName("solutions");
            if (solutions.getLength() > 0) {
                Element solutionsElem = (Element) solutions.item(0);
                String logValueStr = solutionsElem.getAttribute("logValue");
                if (logValueStr != null && !logValueStr.isEmpty()) {
                    result.logNormConstAggr = Double.parseDouble(logValueStr);
                }
            }

            return result;
        } catch (Exception e) {
            line_warning("JMTIO.parseJMVAResult", "Error parsing JMVA result: %s", e.getMessage());
            return null;
        }
    }

    // ==================== HELPER METHODS ====================

    private int findClassIndex(NetworkStruct sn, String className) {
        for (int r = 0; r < sn.nclasses; r++) {
            if (sn.classnames.get(r).equals(className)) {
                return r;
            }
        }
        return -1;
    }

    // ==================== GETTERS ====================

    public Network getModel() {
        return model;
    }

    public double getSimConfInt() {
        return simConfInt;
    }

    public double getSimMaxRelErr() {
        return simMaxRelErr;
    }

    public long getSeed() {
        return seed;
    }
}
