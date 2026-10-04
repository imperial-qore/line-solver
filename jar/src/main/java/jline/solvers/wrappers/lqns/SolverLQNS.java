/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.wrappers.lqns;

import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.io.QN2LQN;
import jline.io.Ret.ProbabilityResult;
import jline.io.Ret.SampleResult;
import jline.solvers.wrappers.lqns.qnsolver.Solver_qns_analyzer;
import jline.GlobalConstants;
import jline.lang.constant.CallType;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.VerboseLevel;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.LayeredNetworkElement;
import jline.lang.layered.LayeredNetworkStruct;
import jline.solvers.*;
import jline.util.matrix.Matrix;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.*;
import java.util.concurrent.TimeUnit;

import static jline.io.InputOutput.*;
import static jline.io.SysUtils.lineTempName;
import static jline.api.sn.SnHasProductForm.snHasProductForm;
import static jline.api.sn.SnHasOpenClasses.snHasOpenClasses;

/**
 * Solver interface to the LQNS external tools: lqns and lqsim for a Layered Queueing
 * Network, and qnsolver for a flat {@link Network}.
 *
 * <p>On a LayeredNetwork the methods are default, lqns, srvn, exactmva, srvn.exactmva,
 * sim, lqsim and lqns.default. On a Network they are default and qns, which leave the
 * multiserver approximation to the default, and qns.conway, qns.rolia, qns.zhou,
 * qns.suri, qns.reiser and qns.schmidt, which name it. A product-form or open Network
 * is written as a JMVA document for qnsolver; a closed non-product-form one is
 * converted by {@link jline.io.QN2LQN} and solved with lqns. This is the arrangement
 * SolverJMT has for JMVA. The layered table is {@link #getLNAvgTable()}.</p>
 * 
 * <p>SolverLQNS provides integration with the external LQNS (Layered Queueing Network Solver)
 * command-line tool developed at Carleton University. This solver enables high-performance
 * analysis of complex layered queueing networks through native C++ implementations.</p>
 * 
 * <p>Key LQNS integration capabilities:
 * <ul>
 *   <li>External tool integration via command-line interface</li>
 *   <li>LQN model serialization to LQNS XML format</li>
 *   <li>High-performance native solver algorithms</li>
 *   <li>Result parsing and integration back to LINE</li>
 *   <li>Support for complex software system models</li>
 * </ul>
 * </p>
 * 
 * <p><strong>Requirements:</strong> This solver requires the 'lqns' and 'lqsim' 
 * command-line tools to be installed and available in the system PATH. Tools can be 
 * obtained from: http://www.sce.carleton.ca/rads/lqns/</p>
 * 
 * @see jline.lang.layered.LayeredNetwork
 * @see LQNSOptions
 * @since 1.0
 */
public class SolverLQNS extends NetworkSolver {

    /** The layered model, or null when this solver holds a flat Network. */
    private LayeredNetwork lnModel;

    public SolverLQNS(LayeredNetwork lqnmodel) {
        this(lqnmodel, new SolverOptions(SolverType.LQNS));
    }

    public SolverLQNS(LayeredNetwork lqnmodel, SolverOptions options) {
        super(null, "SolverLQNS", options);
        this.lnModel = lqnmodel;
        // Solver.model is what the layered callers read; NetworkSolver.model stays null
        ((Solver) this).model = lqnmodel;
        //setOptions(Solver.parseOptions(varargin, defaultOptions()));
        if (!isAvailable() && !options.config.remote) {
            throw new RuntimeException(
                    "SolverLQNS requires the 'lqns' and 'lqsim' commands to be available in your system PATH.\n" +
                            "Obtain them from their authors at: http://www.sce.carleton.ca/rads/lqns/\n" +
                            "LINE ships no LQNS binary and does not redistribute one.\n\n" +
                            "Alternatively, point LINE at a host that already runs LQNS:\n" +
                            "     options.config.remote = true;\n" +
                            "     options.config.remote_url = \"http://localhost:8080\";\n\n");
        }
    }

    /**
     * A solver for a flat Network, whose method names are the qns ones. Unlike the
     * layered constructors it does not require lqns, which only a closed
     * non-product-form model reaches; the run checks for the binary it needs.
     */
    public SolverLQNS(Network model) {
        super(model, "SolverLQNS", networkDefaultOptions());
        this.result = new SolverResult();
    }

    public SolverLQNS(Network model, String method) {
        super(model, "SolverLQNS", networkDefaultOptions());
        this.options.method = method;
        this.result = new SolverResult();
    }

    public SolverLQNS(Network model, SolverOptions options) {
        super(model, "SolverLQNS", options);
        this.result = new SolverResult();
    }

    public SolverLQNS(Network model, Object... varargin) {
        super(model, "SolverLQNS", parseNetworkOptions(varargin));
        this.result = new SolverResult();
    }

    /** The layered model this solver holds, or null for a flat Network. */
    public LayeredNetwork getLayeredModel() {
        return this.lnModel;
    }

    /**
     * Options for a flat Network from name-value pairs: method, multiserver, timespan.
     */
    public static SolverOptions parseNetworkOptions(Object... varargin) {
        SolverOptions options = networkDefaultOptions();

        for (int i = 0; i < varargin.length; i += 2) {
            if (i + 1 < varargin.length && varargin[i] instanceof String) {
                String key = (String) varargin[i];
                Object value = varargin[i + 1];

                switch (key.toLowerCase()) {
                    case "method":
                        if (value instanceof String) {
                            options.method = (String) value;
                        }
                        break;
                    case "multiserver":
                        if (value instanceof String) {
                            options.config.multiserver = (String) value;
                        }
                        break;
                    case "timespan":
                        if (value instanceof double[]) {
                            options.timespan = (double[]) value;
                        }
                        break;
                }
            }
        }

        return options;
    }

    /**
     * The default options on a flat Network (the former SolverQNS defaults).
     */
    public static SolverOptions networkDefaultOptions() {
        SolverOptions options = new SolverOptions(SolverType.LQNS);
        options.method = "default";
        options.config.multiserver = "default";
        options.timespan = new double[]{GlobalConstants.Inf, GlobalConstants.Inf};
        return options;
    }

    /**
     * The method names on a flat Network: "default" and "qns" leave the
     * multiserver approximation to the default, "qns.NAME" selects NAME.
     */
    public static List<String> qnsMethods() {
        return Arrays.asList("default", "qns", "qns.conway", "qns.rolia", "qns.zhou",
                "qns.suri", "qns.reiser", "qns.schmidt");
    }

    /** True for "qns" and every "qns.NAME". */
    public static boolean isQnsMethod(String method) {
        if (method == null) {
            return false;
        }
        String m = method.toLowerCase();
        return m.equals("qns") || m.startsWith("qns.");
    }

    /**
     * The multiserver approximation a Network method selects: the part after
     * "qns.", or "default" for "default" and "qns".
     */
    public static String qnsMultiserver(String method) {
        if (method == null) {
            return "default";
        }
        String m = method.toLowerCase();
        return m.startsWith("qns.") ? m.substring(4) : "default";
    }

    /**
     * The layered feature set: what a layer of a LayeredNetwork may use.
     *
     * @return the feature set supported on the layered path
     */
    public static FeatureSet getLayeredFeatureSet() {
        FeatureSet s = new FeatureSet();
        String[] features = {
                "Sink",
                "Source",
                "Queue",
                "Coxian",
                "Erlang",
                "Exp",
                "HyperExp",
                "Buffer",
                "Server",
                "JobSink",
                "RandomSource",
                "ServiceTunnel",
                "SchedStrategy_PS",
                "SchedStrategy_FCFS",
                "ClosedClass"
        };
        s.setTrue(features);
        return s;
    }

    /**
     * The feature set on a flat Network (the qns methods). The layered path is
     * {@link #getLayeredFeatureSet()}.
     */
    public static FeatureSet getFeatureSet() {
        FeatureSet featSupported = new FeatureSet();

        // Node types
        featSupported.setTrue("Source");
        featSupported.setTrue("Sink");
        featSupported.setTrue("Queue");
        featSupported.setTrue("Delay");
        featSupported.setTrue("DelayStation");
        featSupported.setTrue("Router");
        featSupported.setTrue("ClassSwitch");
        featSupported.setTrue("Fork");
        featSupported.setTrue("Join");
        featSupported.setTrue("Forker");
        featSupported.setTrue("Joiner");
        featSupported.setTrue("Logger");

        // Service distributions
        featSupported.setTrue("Exp");
        featSupported.setTrue("HyperExp");
        featSupported.setTrue("Coxian");
        featSupported.setTrue("Cox2");
        featSupported.setTrue("APH");
        featSupported.setTrue("Erlang");
        featSupported.setTrue("Det");
        featSupported.setTrue("Gamma");
        featSupported.setTrue("Lognormal");
        featSupported.setTrue("MAP");
        featSupported.setTrue("MMPP2");
        featSupported.setTrue("Normal");
        featSupported.setTrue("PH");
        featSupported.setTrue("Pareto");
        featSupported.setTrue("Weibull");
        featSupported.setTrue("Uniform");
        featSupported.setTrue("Trace");
        featSupported.setTrue("Replayer");

        // Server types
        featSupported.setTrue("InfiniteServer");
        featSupported.setTrue("SharedServer");
        featSupported.setTrue("Server");
        featSupported.setTrue("Buffer");
        featSupported.setTrue("Dispatcher");
        featSupported.setTrue("JobSink");
        featSupported.setTrue("RandomSource");
        featSupported.setTrue("ServiceTunnel");
        featSupported.setTrue("LogTunnel");
        featSupported.setTrue("StatelessClassSwitcher");

        // Petri Net elements
        featSupported.setTrue("Linkage");
        featSupported.setTrue("Enabling");
        featSupported.setTrue("Timing");
        featSupported.setTrue("Firing");
        featSupported.setTrue("Storage");
        featSupported.setTrue("Place");
        featSupported.setTrue("Transition");

        // Scheduling strategies
        featSupported.setTrue("SchedStrategy_INF");
        featSupported.setTrue("SchedStrategy_PS");
        featSupported.setTrue("SchedStrategy_DPS");
        featSupported.setTrue("SchedStrategy_FCFS");
        featSupported.setTrue("SchedStrategy_GPS");
        featSupported.setTrue("SchedStrategy_SIRO");
        featSupported.setTrue("SchedStrategy_HOL");
        featSupported.setTrue("SchedStrategy_LCFS");
        featSupported.setTrue("SchedStrategy_LCFSPR");
        featSupported.setTrue("SchedStrategy_SEPT");
        featSupported.setTrue("SchedStrategy_LEPT");
        featSupported.setTrue("SchedStrategy_SJF");
        featSupported.setTrue("SchedStrategy_LJF");
        featSupported.setTrue("SchedStrategy_EXT");

        // Routing strategies
        featSupported.setTrue("RoutingStrategy_PROB");
        featSupported.setTrue("RoutingStrategy_RAND");
        featSupported.setTrue("RoutingStrategy_RROBIN");
        featSupported.setTrue("RoutingStrategy_WRROBIN");
        featSupported.setTrue("RoutingStrategy_SQ");

        // Job classes
        featSupported.setTrue("OpenClass");
        featSupported.setTrue("ClosedClass");

        // c-server stations: the JMVA document carries the count as an
        // <ldstation> and the LQN as a host multiplicity; "suri" and "schmidt"
        // refuse one on the qnsolver path, which stays structural (it is a
        // condition on options.config.multiserver, not on the model).
        // FiniteCapacity is NOT declared: neither document has a buffer, which
        // is what the binding-capacity gate in supportsModelMethod refuses.
        featSupported.setTrue("MultiServer");

        return featSupported;
    }

    /**
     * True if LQNS can run: a native lqns binary is on the PATH. LINE never runs
     * LQNS from a container image, because its licence is an evaluation
     * agreement that forbids redistribution, so the binary must be one the user
     * installed. To exercise a containerised build in the test suite, put a shim
     * on the PATH with {@code run-tests.sh --lqns-docker}.
     */
    public static boolean isAvailable() {
        return hasNativeLqns();
    }

    /** True if a native lqns binary is available on the PATH. */
    public static boolean hasNativeLqns() {
        Process process = null;
        try {
            // Use --help instead of --version because lqsim --version hangs waiting for input
            ProcessBuilder pb = new ProcessBuilder("lqns", "--help");
            pb.redirectErrorStream(true);
            process = pb.start();

            // Read output in a separate thread to prevent blocking
            final StringBuilder output = new StringBuilder();
            final Process finalProcess = process;
            Thread outputReader = new Thread(new Runnable() {
                @Override
                public void run() {
                    try {
                        BufferedReader reader = new BufferedReader(
                            new InputStreamReader(finalProcess.getInputStream()));
                        String line;
                        while ((line = reader.readLine()) != null) {
                            output.append(line).append("\n");
                        }
                        reader.close();
                    } catch (IOException e) {
                        // Ignore - process may have been killed
                    }
                }
            });
            outputReader.start();

            // Use timeout to prevent hanging on Windows/WSL
            boolean finished = process.waitFor(5, TimeUnit.SECONDS);
            if (!finished) {
                process.destroyForcibly();
                outputReader.interrupt();
                return false;
            }

            // Wait for output reader to finish
            outputReader.join(1000);

            // Check if output contains help info indicating the tool is available
            String outputStr = output.toString().toLowerCase();
            return outputStr.contains("lqns") || outputStr.contains("usage");
        } catch (IOException e) {
            return false;
        } catch (InterruptedException e) {
            return false;
        } finally {
            if (process != null) {
                process.destroyForcibly();
            }
        }
    }

    /** True if a native qnsolver binary is on the PATH (the qns methods on a Network). */
    public static boolean hasQnsolver() {
        return Solver_qns_analyzer.hasNativeQNSolver();
    }

    @Override
    public SolverResult getAvg() {
        if (this.lnModel != null) {
            return this.getEnsembleAvg();
        }
        return super.getAvg();
    }

    /**
     * The station table of a flat Network. A LayeredNetwork is not indexed by
     * station, so its table is {@link #getLNAvgTable()}.
     */
    @Override
    public NetworkAvgTable getAvgTable() {
        if (this.lnModel != null) {
            throw new RuntimeException("getAvgTable is indexed by station and this SolverLQNS holds a "
                    + "LayeredNetwork; use getLNAvgTable() for the layered table.");
        }
        return super.getAvgTable();
    }

    /** The per-element table of a LayeredNetwork. */
    public final LayeredNetworkAvgTable getLNAvgTable() {
        if (this.lnModel == null) {
            throw new RuntimeException("getLNAvgTable needs a LayeredNetwork and this SolverLQNS holds a "
                    + "flat Network; use getAvgTable().");
        }
        return jline.io.LineResultRecorder.around(this, "layered", () -> getLNAvgTableImpl());
    }

    /**
     * Body of {@link #getLNAvgTable()}, split out so {@link jline.io.LineResultRecorder}
     * sees what the getter RETURNED. The JAVA cross-codebase parity row is
     * measured from that rather than from what an example printed.
     */
    protected final LayeredNetworkAvgTable getLNAvgTableImpl() {
        LayeredSolverResult result = (LayeredSolverResult) this.getAvg();
        LayeredNetworkStruct lqn = this.getStruct();

        List<String> nodeNames = new ArrayList<>(lqn.names.values());
        List<String> nodeTypes = new ArrayList<>();
        for (int o = 0; o < nodeNames.size(); o++) {
            switch ((int) lqn.type.get(o)) {
                case LayeredNetworkElement.PROCESSOR:
                    nodeTypes.add("Processor");
                    break;
                case LayeredNetworkElement.TASK:
                    if (lqn.sched.get(o) == SchedStrategy.REF) {
                        nodeTypes.add("RefTask");
                    } else {
                        nodeTypes.add("Task");
                    }
                    break;
                case LayeredNetworkElement.ENTRY:
                    nodeTypes.add("Entry");
                    break;
                case LayeredNetworkElement.ACTIVITY:
                    nodeTypes.add("Activity");
                    break;
                case LayeredNetworkElement.CALL:
                    nodeTypes.add("Call");
                    break;
            }
        }
        LayeredNetworkAvgTable AvgTable = new LayeredNetworkAvgTable(result.QN.toList1D(), result.UN.toList1D(), result.RN.toList1D(), result.WN.toList1D(), result.AN.toList1D(), result.TN.toList1D());
        AvgTable.setNodeNames(nodeNames);
        AvgTable.setNodeTypes(nodeTypes);
        AvgTable.setOptions(this.options);
        return AvgTable;
    }

    protected SolverResult getEnsembleAvg() {
        try {
            this.runAnalyzer();
        } catch (IllegalAccessException | ParserConfigurationException e) {
            throw new RuntimeException(e);
        }

        LayeredSolverResult lqnsResult = (LayeredSolverResult) this.result.deepCopy();
        lqnsResult.QN = this.result.UN.copy();
        lqnsResult.UN = ((LayeredSolverResult) this.result).PN.copy();
        lqnsResult.RN = ((LayeredSolverResult) this.result).SN.copy();

        // UN is lqns' proc-utilization, verbatim for hosts, tasks and activities
        // and aggregated over the activity graph for entries, which lqns itself
        // reports as 0 in the activity-graph form. Both lqns and SolverLN report
        // the processor utilization summed over the host's servers, so no
        // rescaling by the host multiplicity applies.
        return lqnsResult;
    }

    public LayeredNetworkStruct getStruct() {
        if (this.lnModel == null) {
            throw new RuntimeException("getStruct returns the LayeredNetworkStruct and this SolverLQNS holds a "
                    + "flat Network; read its NetworkStruct from the model.");
        }
        return this.lnModel.getStruct();
    }

    public List<String> listValidMethods() {
        return listValidMethods(null);
    }

    /**
     * The method names for the model this solver holds: the layered ones on a
     * LayeredNetwork, {@link #qnsMethods()} on a flat Network. A non-null
     * argument asks about that Network instead.
     */
    public List<String> listValidMethods(Network model) {
        if (model != null || this.lnModel == null) {
            return qnsMethods();
        }
        return Arrays.asList("default", "lqns", "srvn", "exactmva", "srvn.exactmva", "sim", "lqsim", "lqns.default");
    }

    /**
     * 1-based position of the layered element called {@code name} whose kind is
     * {@code elemType}, or -1 when the struct carries no such element.
     *
     * The kind is part of the key, and has to be: a LINE-generated layered model
     * routinely gives a processor, its task and that task's entry the same name,
     * and lqn.names holds all three, so a name-only scan resolves to whichever
     * one it meets last and writes one element's .lqxo row into another's slot.
     * The result file states which kind each row describes -- it is the tag being
     * read -- so the ambiguity does not have to exist. Mirrors MATLAB
     * findlqnelem.m and the native-python _find_lqn_elem. See BUGS.md BUG-94.
     */
    private static int findLqnElem(LayeredNetworkStruct lqn, String name, int elemType) {
        int found = -1;
        for (Map.Entry<Integer, String> nameEntry : lqn.names.entrySet()) {
            int idx = nameEntry.getKey().intValue();
            if (!Objects.equals(name, nameEntry.getValue())) {
                continue;
            }
            if (lqn.type == null || idx >= lqn.type.getNumElements() || (int) lqn.type.get(idx) != elemType) {
                continue;
            }
            if (found < 0 || idx < found) {
                found = idx; // first declaration wins, as MATLAB find(...)(1) does
            }
        }
        return found;
    }

    public SolverResult parseXMLResults(String filename) throws IOException {
        SolverResult result = null;
        LayeredNetworkStruct lqn = this.getStruct();
        int numOfNodes = lqn.nidx;
        int numOfCalls = lqn.ncalls;
        Matrix AvgNodesUtilization = new Matrix(numOfNodes, 1);
        AvgNodesUtilization.fill(Double.NaN);
        Matrix AvgNodesPhase1Utilization = new Matrix(numOfNodes, 1);
        AvgNodesPhase1Utilization.fill(Double.NaN);
        Matrix AvgNodesPhase2Utilization = new Matrix(numOfNodes, 1);
        AvgNodesPhase2Utilization.fill(Double.NaN);
        Matrix AvgNodesPhase1ServiceTime = new Matrix(numOfNodes, 1);
        AvgNodesPhase1ServiceTime.fill(Double.NaN);
        Matrix AvgNodesPhase2ServiceTime = new Matrix(numOfNodes, 1);
        AvgNodesPhase2ServiceTime.fill(Double.NaN);
        Matrix AvgNodesThroughput = new Matrix(numOfNodes, 1);
        AvgNodesThroughput.fill(Double.NaN);
        Matrix AvgNodesProcWaiting = new Matrix(numOfNodes, 1);
        AvgNodesProcWaiting.fill(Double.NaN);
        Matrix AvgNodesProcUtilization = new Matrix(numOfNodes, 1);
        AvgNodesProcUtilization.fill(Double.NaN);
        Matrix AvgEdgesWaiting = new Matrix(numOfCalls, 1);
        AvgEdgesWaiting.fill(Double.NaN);

        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder dBuilder = null;
        try {
            dBuilder = dbFactory.newDocumentBuilder();
        } catch (ParserConfigurationException e) {
            throw new RuntimeException(e);
        }

        int dotIndex = filename.lastIndexOf('.');
        String resultFilename = filename.substring(0, dotIndex) + ".lqxo";
        File fXmlFile = new File(resultFilename);
        if (this.options.verbose == VerboseLevel.DEBUG) {
            System.out.println("Parsing LQNS result file: " + resultFilename);
            if (this.options.keep) {
                System.out.println("LQNS result file available at: " + resultFilename);
            }
        }
        Document doc = null;
        try {
            while (!fXmlFile.exists()) {
                TimeUnit.MILLISECONDS.sleep(1);
            }
            doc = dBuilder.parse(fXmlFile.toURI().toString());
        } catch (SAXException | InterruptedException e) {
            throw new RuntimeException(e);
        }
        doc.getDocumentElement().normalize();

        // Parse number of lqns iterations
        NodeList solverParams = doc.getElementsByTagName("solver-params");
        int iterations = 0;
        for (int i = 0; i < solverParams.getLength(); i++) {
            Node solverParamNode = solverParams.item(i);
            if (solverParamNode.getNodeType() == Node.ELEMENT_NODE) {
                Element solverParamElement = (Element) solverParamNode;
                NodeList resultGeneralList = solverParamElement.getElementsByTagName("result-general");
                if (resultGeneralList.getLength() > 0) {
                    Element resultGeneralElement = (Element) resultGeneralList.item(0);
                    String iterationsStr = resultGeneralElement.getAttribute("iterations");
                    try {
                        iterations = Integer.parseInt(iterationsStr);
                    } catch (NumberFormatException e) {
                        line_warning(mfilename(new Object() {
                        }), "Error parsing iterations value: " + iterationsStr);
                    }
                }
            }
        }

        NodeList procList = doc.getElementsByTagName("processor");

        for (int i = 0; i < procList.getLength(); i++) {
            Node procNode = procList.item(i);

            if (procNode.getNodeType() == Node.ELEMENT_NODE) {
                Element procElement = (Element) procNode;
                String procName = procElement.getAttribute("name");

                // Find the position of procName
                int procPos = findLqnElem(lqn, procName, LayeredNetworkElement.HOST);
                NodeList procResultList = procElement.getElementsByTagName("result-processor");
                if (procResultList.getLength() > 0 && procPos >= 0) {
                    Element procResultElement = (Element) procResultList.item(0);
                    String utilizationStr = procResultElement.getAttribute("utilization");
                    double uRes = 0.0;
                    if (!utilizationStr.isEmpty()) {
                        uRes = Double.parseDouble(utilizationStr);
                    }
                    AvgNodesProcUtilization.set(procPos, uRes);
                }

                // Assuming procElement is already defined as an Element and lqn.names is a List or array.
                NodeList taskList = procElement.getElementsByTagName("task");
                for (int j = 0; j < taskList.getLength(); j++) {
                    // Element - Task
                    Element taskElement = (Element) taskList.item(j);
                    String taskName = taskElement.getAttribute("name");

                    // Find the position of taskName
                    int taskPos = findLqnElem(lqn, taskName, LayeredNetworkElement.TASK);
                    NodeList taskResult = taskElement.getElementsByTagName("result-task");
                    Element resultElement = (Element) taskResult.item(0);

                    double uRes = Double.parseDouble(resultElement.getAttribute("utilization"));
                    String p1uResStr = resultElement.getAttribute("phase1-utilization");
                    double p1uRes = Double.NaN;
                    if (!p1uResStr.isEmpty()) {
                        p1uRes = Double.parseDouble(p1uResStr);
                    }
                    String p2uResStr = resultElement.getAttribute("phase2-utilization");
                    double p2uRes = Double.NaN;
                    if (!p2uResStr.isEmpty()) {
                        p2uRes = Double.parseDouble(p2uResStr);
                    }
                    double tRes = Double.parseDouble(resultElement.getAttribute("throughput"));
                    double puRes = Double.parseDouble(resultElement.getAttribute("proc-utilization"));
                    if (taskPos >= 0) {
                        AvgNodesUtilization.set(taskPos, uRes);
                        AvgNodesPhase1Utilization.set(taskPos, p1uRes);
                        AvgNodesThroughput.set(taskPos, tRes);
                        AvgNodesProcUtilization.set(taskPos, puRes);
                        AvgNodesPhase2Utilization.set(taskPos, p2uRes);
                    }

                    NodeList entryList = doc.getElementsByTagName("entry");
                    for (int k = 0; k < entryList.getLength(); k++) {
                        Element entryElement = (Element) entryList.item(k);
                        String entryName = entryElement.getAttribute("name");
                        // Find the position of entryName
                        int entryPos = findLqnElem(lqn, entryName, LayeredNetworkElement.ENTRY);
                        NodeList entryResult = entryElement.getElementsByTagName("result-entry");
                        Element firstEntryResult = (Element) entryResult.item(0);

                        uRes = Double.parseDouble(firstEntryResult.getAttribute("utilization"));
                        p1uResStr = firstEntryResult.getAttribute("phase1-utilization");
                        p1uRes = Double.NaN;
                        if (!p1uResStr.isEmpty()) {
                            p1uRes = Double.parseDouble(p1uResStr);
                        }
                        p2uResStr = firstEntryResult.getAttribute("phase2-utilization");
                        p2uRes = Double.NaN;
                        if (!p2uResStr.isEmpty()) {
                            p2uRes = Double.parseDouble(p2uResStr);
                        }
                        String p1stResStr = firstEntryResult.getAttribute("phase1-service-time");
                        double p1stRes = Double.NaN;
                        if (!p1stResStr.isEmpty()) {
                            p1stRes = Double.parseDouble(p1stResStr);
                        }
                        String p2stResStr = firstEntryResult.getAttribute("phase2-service-time");
                        double p2stRes = Double.NaN;
                        if (!p2stResStr.isEmpty()) {
                            p2stRes = Double.parseDouble(p2stResStr);
                        }
                        tRes = Double.parseDouble(firstEntryResult.getAttribute("throughput"));
                        puRes = Double.parseDouble(firstEntryResult.getAttribute("proc-utilization"));

                        /* If phase-specific attributes are missing from result-entry,
                           extract them from entry-phase-activities result-activity elements */
                        NodeList entryPhaseActsList = entryElement.getElementsByTagName("entry-phase-activities");
                        if (entryPhaseActsList.getLength() > 0) {
                            Element entryPhaseActsElement = (Element) entryPhaseActsList.item(0);
                            NodeList phaseActList = entryPhaseActsElement.getElementsByTagName("activity");
                            for (int pa = 0; pa < phaseActList.getLength(); pa++) {
                                Element phaseActElement = (Element) phaseActList.item(pa);
                                String phaseStr = phaseActElement.getAttribute("phase");
                                if (!phaseStr.isEmpty()) {
                                    int phase = Integer.parseInt(phaseStr);
                                    NodeList phaseActResult = phaseActElement.getElementsByTagName("result-activity");
                                    if (phaseActResult.getLength() > 0) {
                                        Element phaseResult = (Element) phaseActResult.item(0);
                                        if (phase == 1) {
                                            if (Double.isNaN(p1stRes)) {
                                                String stStr = phaseResult.getAttribute("service-time");
                                                if (!stStr.isEmpty()) p1stRes = Double.parseDouble(stStr);
                                            }
                                            if (Double.isNaN(p1uRes)) {
                                                String uStr = phaseResult.getAttribute("utilization");
                                                if (!uStr.isEmpty()) p1uRes = Double.parseDouble(uStr);
                                            }
                                        } else if (phase == 2) {
                                            if (Double.isNaN(p2stRes)) {
                                                String stStr = phaseResult.getAttribute("service-time");
                                                if (!stStr.isEmpty()) p2stRes = Double.parseDouble(stStr);
                                            }
                                            if (Double.isNaN(p2uRes)) {
                                                String uStr = phaseResult.getAttribute("utilization");
                                                if (!uStr.isEmpty()) p2uRes = Double.parseDouble(uStr);
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if (entryPos >= 0) {
                            AvgNodesUtilization.set(entryPos, uRes);
                            AvgNodesPhase1Utilization.set(entryPos, p1uRes);
                            AvgNodesPhase2Utilization.set(entryPos, p2uRes);
                            AvgNodesPhase1ServiceTime.set(entryPos, p1stRes);
                            AvgNodesPhase2ServiceTime.set(entryPos, p2stRes);
                            AvgNodesThroughput.set(entryPos, tRes);
                            AvgNodesProcUtilization.set(entryPos, puRes);
                        }
                    }
                }

                NodeList taskActsList = doc.getElementsByTagName("task-activities");
                if (taskActsList.getLength() > 0) {
                    for (int t = 0; t < taskActsList.getLength(); t++) {
                        Element taskActsElement = (Element) taskActsList.item(t);
                        NodeList actList = taskActsElement.getElementsByTagName("activity");
                        for (int l = 0; l < actList.getLength(); l++) {
                            Element actElement = (Element) actList.item(l);
                            if (actElement.getParentNode().getNodeName().equals("task-activities")) {
                                String actName = actElement.getAttribute("name");

                                // Find the position of actName
                                int actPos = findLqnElem(lqn, actName, LayeredNetworkElement.ACTIVITY);
                                if (actPos < 0) {
                                    continue;
                                }
                                NodeList actResult = actElement.getElementsByTagName("result-activity");
                                double uRes = Double.parseDouble(actResult.item(0).getAttributes().getNamedItem("utilization").getNodeValue());
                                double stRes = Double.parseDouble(actResult.item(0).getAttributes().getNamedItem("service-time").getNodeValue());
                                double tRes = Double.parseDouble(actResult.item(0).getAttributes().getNamedItem("throughput").getNodeValue());
                                String pwResStr = actResult.item(0).getAttributes().getNamedItem("proc-waiting").getNodeValue();
                                double pwRes = Double.NaN;
                                if (!pwResStr.isEmpty()) {
                                    pwRes = Double.parseDouble(pwResStr);
                                }
                                double mypuRes = Double.parseDouble(actResult.item(0).getAttributes().getNamedItem("proc-utilization").getNodeValue());

                                AvgNodesUtilization.set(actPos, uRes);
                                AvgNodesPhase1ServiceTime.set(actPos, stRes);
                                AvgNodesThroughput.set(actPos, tRes);
                                AvgNodesProcWaiting.set(actPos, pwRes);
                                AvgNodesProcUtilization.set(actPos, mypuRes);

                                String actID = lqn.names.get(actPos);
                                // Synchronous calls
                                NodeList synchCalls = actElement.getElementsByTagName("synch-call");
                                for (int m = 0; m < synchCalls.getLength(); m++) {
                                    Element callElement = (Element) synchCalls.item(m);
                                    String destName = callElement.getAttribute("dest");
                                    int destPos = findLqnElem(lqn, destName, LayeredNetworkElement.ENTRY);
                                    if (destPos < 0) {
                                        continue;
                                    }
                                    String destID = lqn.names.get(destPos);
                                    int callPos = -1;

                                    for (Map.Entry<Integer, String> entry : lqn.callnames.entrySet()) {
                                        if (Objects.equals(actID + "=>" + destID, entry.getValue())) {
                                            callPos = entry.getKey().intValue();
                                        }
                                    }
                                    if (callPos < 0) {
                                        continue;
                                    }
                                    NodeList callResult = callElement.getElementsByTagName("result-call");
                                    double wRes = Double.parseDouble(((Element) callResult.item(0)).getAttribute("waiting"));
                                    AvgEdgesWaiting.set(callPos, wRes);
                                }

                                // Asynchronous calls
                                NodeList asynchCalls = actElement.getElementsByTagName("asynch-call");
                                for (int m = 0; m < asynchCalls.getLength(); m++) {
                                    Element callElement = (Element) asynchCalls.item(m);
                                    String destName = callElement.getAttribute("dest");
                                    int destPos = findLqnElem(lqn, destName, LayeredNetworkElement.ENTRY);
                                    if (destPos < 0) {
                                        continue;
                                    }
                                    String destID = lqn.names.get(destPos);
                                    int callPos = -1;
                                    for (Map.Entry<Integer, String> entry : lqn.callnames.entrySet()) {
                                        if (Objects.equals(actID + "->" + destID, entry.getValue())) {
                                            callPos = entry.getKey().intValue();
                                        }
                                    }
                                    if (callPos < 0) {
                                        continue;
                                    }
                                    NodeList callResult = callElement.getElementsByTagName("result-call");
                                    double wRes = Double.parseDouble(((Element) callResult.item(0)).getAttribute("waiting"));
                                    AvgEdgesWaiting.set(callPos, wRes);
                                }
                            }
                        }
                    }
                }
                /* Also parse activity results from entry-phase-activities (PH1PH2 format) */
                NodeList entryPhaseActsListAll = doc.getElementsByTagName("entry-phase-activities");
                for (int ep = 0; ep < entryPhaseActsListAll.getLength(); ep++) {
                    Element epElement = (Element) entryPhaseActsListAll.item(ep);

                    // Find parent entry element to get throughput/proc-utilization
                    // (LQNS omits these attributes from entry-phase result-activity elements)
                    double parentEntryTput = Double.NaN;
                    double parentEntryPU = Double.NaN;
                    org.w3c.dom.Node parentNode = epElement.getParentNode();
                    if (parentNode instanceof Element) {
                        Element parentEntry = (Element) parentNode;
                        NodeList parentResults = parentEntry.getElementsByTagName("result-entry");
                        if (parentResults.getLength() > 0) {
                            Element parentResult = (Element) parentResults.item(0);
                            String ptStr = parentResult.getAttribute("throughput");
                            if (ptStr != null && !ptStr.isEmpty()) parentEntryTput = Double.parseDouble(ptStr);
                            String ppuStr = parentResult.getAttribute("proc-utilization");
                            if (ppuStr != null && !ppuStr.isEmpty()) parentEntryPU = Double.parseDouble(ppuStr);
                        }
                    }

                    NodeList epActList = epElement.getElementsByTagName("activity");
                    for (int l = 0; l < epActList.getLength(); l++) {
                        Element actElement = (Element) epActList.item(l);
                        String actName = actElement.getAttribute("name");

                        int actPos = findLqnElem(lqn, actName, LayeredNetworkElement.ACTIVITY);
                        if (actPos < 0) continue;

                        NodeList actResult = actElement.getElementsByTagName("result-activity");
                        if (actResult.getLength() == 0) continue;
                        Element resultEl = (Element) actResult.item(0);

                        String uResStr = resultEl.getAttribute("utilization");
                        if (!uResStr.isEmpty()) AvgNodesUtilization.set(actPos, Double.parseDouble(uResStr));
                        String stResStr = resultEl.getAttribute("service-time");
                        if (!stResStr.isEmpty()) AvgNodesPhase1ServiceTime.set(actPos, Double.parseDouble(stResStr));
                        String tResStr = resultEl.getAttribute("throughput");
                        if (!tResStr.isEmpty()) {
                            AvgNodesThroughput.set(actPos, Double.parseDouble(tResStr));
                        } else if (!Double.isNaN(parentEntryTput)) {
                            AvgNodesThroughput.set(actPos, parentEntryTput);
                        }
                        String pwResStr = resultEl.getAttribute("proc-waiting");
                        if (pwResStr != null && !pwResStr.isEmpty()) AvgNodesProcWaiting.set(actPos, Double.parseDouble(pwResStr));
                        String puResStr = resultEl.getAttribute("proc-utilization");
                        String hdStr = actElement.getAttribute("host-demand-mean");
                        if (!puResStr.isEmpty()) {
                            AvgNodesProcUtilization.set(actPos, Double.parseDouble(puResStr));
                        } else if (!Double.isNaN(parentEntryTput) && hdStr != null && !hdStr.isEmpty()) {
                            // see _kb/12-interfaces-and-docs.md (Wrappers: JAR subprocess-bridge notes: LQNS entry-phase utilization fallback)
                            AvgNodesProcUtilization.set(actPos, parentEntryTput * Double.parseDouble(hdStr));
                        } else if (!Double.isNaN(parentEntryPU) && epActList.getLength() == 1) {
                            // No throughput/demand available: the entry total is safe
                            // only when the entry has a single phase activity.
                            AvgNodesProcUtilization.set(actPos, parentEntryPU);
                        }

                        /* Parse synch-call waiting times */
                        String actID = lqn.names.get(actPos);
                        NodeList synchCalls = actElement.getElementsByTagName("synch-call");
                        for (int m = 0; m < synchCalls.getLength(); m++) {
                            Element callElement = (Element) synchCalls.item(m);
                            String destName = callElement.getAttribute("dest");
                            int destPos = findLqnElem(lqn, destName, LayeredNetworkElement.ENTRY);
                            if (destPos < 0) {
                                continue;
                            }
                            String destID = lqn.names.get(destPos);
                            int callPos = 0;
                            for (Map.Entry<Integer, String> entry : lqn.callnames.entrySet()) {
                                if (Objects.equals(actID + "->" + destID, entry.getValue())) {
                                    callPos = entry.getKey().intValue();
                                }
                            }
                            if (callPos > 0) {
                                NodeList callResult = callElement.getElementsByTagName("result-call");
                                if (callResult.getLength() > 0) {
                                    double wRes = Double.parseDouble(((Element) callResult.item(0)).getAttribute("waiting"));
                                    AvgEdgesWaiting.set(callPos - 1, wRes);
                                }
                            }
                        }
                    }
                }
            }
        }
        this.result = new LayeredSolverResult();

        // Processor utilization of an entry, aggregated from its activity graph.
        // lqns credits host work to whichever level carries the host demand: in
        // the activity-graph form an entry declares none, so lqns reports
        // result-entry proc-utilization as a literal 0 and the work sits on the
        // result-activity rows. The entry value is then the sum over the
        // activities reachable from the entry within its own task, which is what
        // lqn.actsof holds. In PH1PH2 form the same sum runs over the phase
        // activities and reproduces the value lqns reports there, so no form test
        // is needed. An entry with no activities, or any activity lqns left
        // unreported, keeps the raw attribute rather than a partial sum.
        for (int eoff = 0; eoff < lqn.nentries; eoff++) {
            int eidx = lqn.eshift + eoff;
            List<Integer> acts = lqn.actsof.get(eidx);
            if (acts == null || acts.isEmpty()) {
                continue;
            }
            double sumPU = 0;
            boolean complete = true;
            for (int a = 0; a < acts.size(); a++) {
                double pu = AvgNodesProcUtilization.get(acts.get(a));
                if (Double.isNaN(pu)) {
                    complete = false;
                    break;
                }
                sumPU += pu;
            }
            if (complete) {
                AvgNodesProcUtilization.set(eidx, sumPU);
            }
        }

        // Phase-1 service time of an entry lqns never invoked.
        // lqns omits phase1-service-time from result-entry exactly when the entry's
        // throughput is zero: nothing was served, so there is no per-invocation mean
        // to report. LINE then carried a NaN where the table says an entry HAS a
        // response time and every other solver reports one, breaking the NaN mask --
        // see _kb/06-solver-catalog.md. The value is taken from the activity rows, and
        // ONLY where they are unanimous: if every activity reachable from the entry
        // reports a zero service time then every aggregation law agrees on zero -- the
        // serial sum, the branch-weighted mean of an OrFork, the order statistic of an
        // AndFork -- so the derivation does not depend on which one applies.
        // It is deliberately NOT generalised the way ProcUtilization is above.
        // Utilizations add over an activity graph; response times do not. Measured over
        // the example corpus, sum(actsof) reproduces phase1-service-time on serial
        // chains only and misses it wherever the graph branches (lqn_workflows `Entry`:
        // 12.5667 reported against 8.5667 summed, lqn_fork_open_arrival `SE`: 0.841667
        // against 1.0), so a summed fallback would answer with a number lqns
        // contradicts. An entry whose activities are unreported, absent, or not all
        // zero keeps NaN.
        for (int eoff = 0; eoff < lqn.nentries; eoff++) {
            int eidx = lqn.eshift + eoff;
            if (!Double.isNaN(AvgNodesPhase1ServiceTime.get(eidx))) {
                continue;
            }
            List<Integer> acts = lqn.actsof.get(eidx);
            if (acts == null || acts.isEmpty()) {
                continue;
            }
            boolean allZero = true;
            for (int a = 0; a < acts.size(); a++) {
                double st = AvgNodesPhase1ServiceTime.get(acts.get(a));
                if (Double.isNaN(st) || st != 0.0) {
                    allZero = false;
                    break;
                }
            }
            if (allZero) {
                AvgNodesPhase1ServiceTime.set(eidx, 0.0);
            }
        }

        ((LayeredSolverResult) this.result).PN = new Matrix(AvgNodesProcUtilization);
        ((LayeredSolverResult) this.result).SN = new Matrix(AvgNodesPhase1ServiceTime);
        this.result.TN = new Matrix(AvgNodesThroughput);
        this.result.UN = new Matrix(AvgNodesUtilization);
        this.result.RN = new Matrix(AvgNodesProcWaiting);
        this.result.RN.fill(Double.NaN);
        this.result.QN = new Matrix(AvgNodesProcWaiting);
        this.result.QN.fill(Double.NaN);
        this.result.AN = new Matrix(lqn.nidx, 1);
        this.result.AN.fill(Double.NaN);
        this.result.WN = new Matrix(lqn.nidx, 1);
        this.result.WN.fill(Double.NaN);

        // Store raw metrics for getRawAvgTables (aligned with MATLAB RawAvg structure)
        LayeredSolverResult layeredResult = (LayeredSolverResult) this.result;
        layeredResult.rawUtilization = new Matrix(AvgNodesUtilization);
        layeredResult.rawPhase1Utilization = new Matrix(AvgNodesPhase1Utilization);
        layeredResult.rawPhase2Utilization = new Matrix(AvgNodesPhase2Utilization);
        layeredResult.rawPhase1ServiceTime = new Matrix(AvgNodesPhase1ServiceTime);
        layeredResult.rawPhase2ServiceTime = new Matrix(AvgNodesPhase2ServiceTime);
        layeredResult.rawThroughput = new Matrix(AvgNodesThroughput);
        layeredResult.rawProcWaiting = new Matrix(AvgNodesProcWaiting);
        layeredResult.rawProcUtilization = new Matrix(AvgNodesProcUtilization);
        layeredResult.rawEdgesWaiting = new Matrix(AvgEdgesWaiting);

        return result;
    }

    public void runAnalyzer() throws IllegalAccessException, ParserConfigurationException {
        runAnalyzer(this.options);
    }

    public void runAnalyzer(SolverOptions options)
            throws IllegalAccessException, ParserConfigurationException {

        if (this.lnModel == null) {
            if (options != null) {
                this.options = options;
            }
            runAnalyzerNetwork();
            return;
        }
        long t0 = System.nanoTime();
        if (options == null) options = this.options;
        // Propagate solver verbose level to global
        if (options != null) {
            GlobalConstants.Verbose = options.verbose;
        }
        jline.io.InputOutput.line_ack(options.verbose, "LQNS");
        // Solver console: on a LayeredNetwork getAvg goes to getEnsembleAvg, so
        // the run never reaches NetworkSolver.getAvg where every other solver's
        // run is opened. It opens its own here and closes it in the finally below.
        jline.io.LineConsole.beginRun(this, options);
        try {
            runAnalyzerBody(options, t0);
        } finally {
            jline.io.LineConsole.closeRun(this);
        }
    }

    private void runAnalyzerBody(SolverOptions options, long t0)
            throws IllegalAccessException, ParserConfigurationException {
        line_debug(options.verbose, String.format("LQNS solver starting: method=%s, multiserver=%s",
            options.method, options.config.multiserver != null ? options.config.multiserver : "default"));

        /* --- write .lqnx ------------------------------------------------ */
        String dirPath = null;
        try {
            dirPath = lineTempName("lqns", false);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        String fileName = dirPath + File.separator + "model.lqnx";
        jline.io.LineConsole.step("writing the LQN model to %s", fileName);
        this.lnModel.writeXML(fileName, false);
        jline.io.LineConsole.step("running the lqns binary as a subprocess");

        resetRandomGeneratorSeed(options.seed);

        /* --- CLI switches ------------------------------------------------ */
        String verboseFlag =
                (options.verbose == VerboseLevel.SILENT) ? "-a -w" : "";

        /* multiserver policy (Praqma) */
        String praqmaFlag = "";
        if (!"lqsim".equalsIgnoreCase(options.method)) {
            String pol = options.config.multiserver;          // default 'rolia'
            if (pol == null || pol.isEmpty() || "default".equals(pol)) pol = "rolia";
            praqmaFlag = "-Pmultiserver=" + pol;
        }

        /* build command */
        String cmd;
        int iterMax = options.iter_max;       /* iteration limit if needed */
        String common = verboseFlag + " " + praqmaFlag + " -Pstop-on-message-loss=false -x " + fileName;

        switch (options.method.toLowerCase()) {
            case "srvn":
                cmd = "lqns " + common + " -Playering=srvn";
                line_debug(options.verbose, "LQNS: using SRVN layering method");
                break;
            case "exactmva":
                cmd = "lqns " + common + " -Pmva=exact";
                line_debug(options.verbose, "LQNS: using exact MVA method");
                break;
            case "srvn.exactmva":
                cmd = "lqns " + common + " -Playering=srvn -Pmva=exact";
                line_debug(options.verbose, "LQNS: using SRVN layering with exact MVA");
                break;
            case "sim":
            case "lqsim":
                cmd = "lqsim " + common + " -A " + options.samples + " ";
                line_debug(options.verbose, String.format("LQNS: using LQSIM simulation, samples=%d", options.samples));
                break;
            case "lqns.default":
                cmd = "lqns " + verboseFlag + " " + praqmaFlag + " -x " + fileName;
                line_debug(options.verbose, "LQNS: using lqns.default method");
                break;
            case "default":
            case "lqns":
            default:
                cmd = "lqns " + common;
                line_debug(options.verbose, "LQNS: using default analytical method");
                break;
        }

        if (options.verbose == VerboseLevel.DEBUG) {
            System.out.println("\nLQNS command: " + cmd);
        }

        /* --- run --------------------------------------------------------- */
        line_debug(options.verbose, String.format("LQNS: executing command: %s", cmd));
        // Check for remote execution
        if (options.config.remote) {
            line_debug(options.verbose, String.format("LQNS: using remote execution at %s", options.config.remote_url));
            if (options.verbose == VerboseLevel.DEBUG) {
                System.out.println("\nUsing remote LQNS at: " + options.config.remote_url);
            }
            runRemoteLQNS(fileName, options);
        } else {
            // Local execution
            try {
                ProcessBuilder pb = new ProcessBuilder(cmd.split("\\s+"));
                pb.redirectErrorStream(true); // Merge stderr into stdout
                final Process p = pb.start();

                // see _kb/12-interfaces-and-docs.md (Wrappers: JAR subprocess-bridge notes: LQNS/LQSIM watchdog)
                final boolean[] killedByTimeout = {false};
                Thread watchdog = null;
                if (options != null && Double.isFinite(options.timeout) && options.timeout > 0) {
                    final long budgetMillis = (long) Math.ceil(options.timeout * 1000.0);
                    watchdog = new Thread(new Runnable() {
                        public void run() {
                            try {
                                if (!p.waitFor(budgetMillis, TimeUnit.MILLISECONDS)) {
                                    killedByTimeout[0] = true;
                                    p.destroyForcibly();
                                }
                            } catch (InterruptedException ie) {
                                // stop watching
                            }
                        }
                    });
                    watchdog.setDaemon(true);
                    watchdog.start();
                }

                // Consume output to prevent buffer blocking (especially on Windows)
                BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
                StringBuilder output = new StringBuilder();
                String line;
                while ((line = reader.readLine()) != null) {
                    output.append(line).append("\n");
                }

                int status = p.waitFor();
                if (watchdog != null) {
                    watchdog.interrupt();
                }
                if (killedByTimeout[0]) {
                    if (this.result != null) {
                        this.result.timedOut = true;
                    }
                    line_warning(mfilename(new Object() {
                            }),
                            "LQNS/LQSIM exceeded the wall-clock time budget (options.timeout=" + options.timeout + "s) and was terminated; returning an empty result.");
                    return;
                }
                if (status != 0) {
                    // Log captured output as error
                    String[] lines = output.toString().split("\n");
                    for (String errLine : lines) {
                        if (!errLine.isEmpty()) {
                            line_error(mfilename(new Object(){}), errLine);
                        }
                    }
                    line_error(mfilename(new Object() {
                            }),
                            "LQNS/LQSIM did not terminate correctly.");
                }
            } catch (IOException | InterruptedException e) {
                throw new RuntimeException(e);
            }
        }

        /* --- parse results ---------------------------------------------- */
        jline.io.LineConsole.step("parsing the lqns XML results");
        try {
            parseXMLResults(fileName);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        if (!options.keep) {
            File tempDir = new File(dirPath);
            if (tempDir.exists() && !tempDir.delete()) {
                line_warning(mfilename(new Object() {}), "Could not delete temporary directory: " + dirPath);
            }
        }

        long t1 = System.nanoTime();
        result.runtime = (t1 - t0) / 1000000000.0;
    }



    /**
     * Execute LQNS via remote REST API.
     *
     * @param fileName Path to the LQNX model file
     * @param options Solver options containing remote configuration
     */
    private void runRemoteLQNS(String fileName, SolverOptions options) {
        try {
            // Read LQNX model file
            String modelContent = new String(Files.readAllBytes(new File(fileName).toPath()), StandardCharsets.UTF_8);

            // Determine endpoint based on method
            String baseUrl = options.config.remote_url;
            if (!baseUrl.endsWith("/")) {
                baseUrl += "/";
            }
            String endpoint;
            if ("lqsim".equalsIgnoreCase(options.method) || "sim".equalsIgnoreCase(options.method)) {
                endpoint = baseUrl + "api/v1/solve/lqsim";
            } else {
                endpoint = baseUrl + "api/v1/solve/lqns";
            }

            // Build JSON request
            StringBuilder jsonBuilder = new StringBuilder();
            jsonBuilder.append("{\"model\":{\"content\":");
            jsonBuilder.append(escapeJson(modelContent));
            jsonBuilder.append(",\"base64\":false},\"options\":{\"include_raw_output\":true");

            // Add pragmas based on method
            if (!"lqsim".equalsIgnoreCase(options.method) && !"sim".equalsIgnoreCase(options.method)) {
                String multiserver = options.config.multiserver;
                if (multiserver == null || multiserver.isEmpty() || "default".equals(multiserver)) {
                    multiserver = "rolia";
                }
                jsonBuilder.append(",\"pragmas\":{\"multiserver\":\"").append(multiserver).append("\"");
                jsonBuilder.append(",\"stop_on_message_loss\":false");

                // Method-specific pragmas
                if ("srvn".equalsIgnoreCase(options.method)) {
                    jsonBuilder.append(",\"layering\":\"srvn\"");
                } else if ("exactmva".equalsIgnoreCase(options.method)) {
                    jsonBuilder.append(",\"mva\":\"exact-mva\"");
                } else if ("srvn.exactmva".equalsIgnoreCase(options.method)) {
                    jsonBuilder.append(",\"layering\":\"srvn\",\"mva\":\"exact-mva\"");
                }
                jsonBuilder.append("}");
            } else {
                // LQSIM options
                jsonBuilder.append(",\"blocks\":30");
                if (options.samples > 0) {
                    jsonBuilder.append(",\"run_time\":").append(options.samples);
                }
            }
            jsonBuilder.append("}}");

            String jsonRequest = jsonBuilder.toString();

            // Make HTTP request
            URL url = new URL(endpoint);
            HttpURLConnection conn = (HttpURLConnection) url.openConnection();
            conn.setRequestMethod("POST");
            conn.setRequestProperty("Content-Type", "application/json");
            conn.setRequestProperty("Accept", "application/json");
            conn.setConnectTimeout(30000);  // 30 seconds
            conn.setReadTimeout(300000);    // 5 minutes
            conn.setDoOutput(true);

            try (OutputStream os = conn.getOutputStream()) {
                byte[] input = jsonRequest.getBytes(StandardCharsets.UTF_8);
                os.write(input, 0, input.length);
            }

            int responseCode = conn.getResponseCode();
            if (responseCode != 200) {
                BufferedReader br = new BufferedReader(new InputStreamReader(conn.getErrorStream(), StandardCharsets.UTF_8));
                StringBuilder errorResponse = new StringBuilder();
                String line;
                while ((line = br.readLine()) != null) {
                    errorResponse.append(line);
                }
                br.close();
                throw new RuntimeException("Remote LQNS returned error: " + responseCode + " - " + errorResponse);
            }

            // Read response
            BufferedReader br = new BufferedReader(new InputStreamReader(conn.getInputStream(), StandardCharsets.UTF_8));
            StringBuilder response = new StringBuilder();
            String line;
            while ((line = br.readLine()) != null) {
                response.append(line);
            }
            br.close();

            // Parse JSON response to extract LQXO content
            String jsonResponse = response.toString();
            String lqxoContent = extractLqxoFromJson(jsonResponse);

            if (lqxoContent != null && !lqxoContent.isEmpty()) {
                // Write LQXO to file for parsing by existing code
                String lqxoFile = fileName.replace(".lqnx", ".lqxo");
                Files.write(new File(lqxoFile).toPath(), lqxoContent.getBytes(StandardCharsets.UTF_8));
            } else {
                // Check for error in response
                if (jsonResponse.contains("\"status\":\"error\"") || jsonResponse.contains("\"status\":\"failed\"")) {
                    String error = extractErrorFromJson(jsonResponse);
                    throw new RuntimeException("Remote solver returned error: " + error);
                }
                throw new RuntimeException("Remote solver did not return LQXO output");
            }

        } catch (IOException e) {
            throw new RuntimeException("Remote LQNS execution failed: " + e.getMessage(), e);
        }
    }

    /**
     * Escape string for JSON encoding.
     */
    private String escapeJson(String s) {
        if (s == null) return "null";
        StringBuilder sb = new StringBuilder("\"");
        for (char c : s.toCharArray()) {
            switch (c) {
                case '"': sb.append("\\\""); break;
                case '\\': sb.append("\\\\"); break;
                case '\b': sb.append("\\b"); break;
                case '\f': sb.append("\\f"); break;
                case '\n': sb.append("\\n"); break;
                case '\r': sb.append("\\r"); break;
                case '\t': sb.append("\\t"); break;
                default:
                    if (c < ' ') {
                        sb.append(String.format("\\u%04x", (int) c));
                    } else {
                        sb.append(c);
                    }
            }
        }
        sb.append("\"");
        return sb.toString();
    }

    /**
     * Extract LQXO content from JSON response.
     */
    private String extractLqxoFromJson(String json) {
        // Simple extraction - look for "lqxo":"..." pattern
        String marker = "\"lqxo\":\"";
        int start = json.indexOf(marker);
        if (start < 0) {
            marker = "\"lqxo\": \"";
            start = json.indexOf(marker);
        }
        if (start < 0) return null;

        start += marker.length();
        StringBuilder sb = new StringBuilder();
        boolean escape = false;
        for (int i = start; i < json.length(); i++) {
            char c = json.charAt(i);
            if (escape) {
                switch (c) {
                    case 'n': sb.append('\n'); break;
                    case 'r': sb.append('\r'); break;
                    case 't': sb.append('\t'); break;
                    case '\\': sb.append('\\'); break;
                    case '"': sb.append('"'); break;
                    default: sb.append(c);
                }
                escape = false;
            } else if (c == '\\') {
                escape = true;
            } else if (c == '"') {
                break;
            } else {
                sb.append(c);
            }
        }
        return sb.toString();
    }

    /**
     * Extract error message from JSON response.
     */
    private String extractErrorFromJson(String json) {
        String marker = "\"error\":\"";
        int start = json.indexOf(marker);
        if (start < 0) {
            marker = "\"error\": \"";
            start = json.indexOf(marker);
        }
        if (start < 0) return "Unknown error";

        start += marker.length();
        int end = json.indexOf("\"", start);
        if (end < 0) return "Unknown error";

        return json.substring(start, end);
    }

    /**
     * The lqsim simulator is stochastic; the analytical lqns/srvn methods are
     * deterministic.
     *
     * @param method the method name to classify
     * @return true if the method returns stochastic estimates
     */
    @Override
    public boolean isStochasticMethod(String method) {
        return method != null && (method.equalsIgnoreCase("sim") || method.equalsIgnoreCase("lqsim"));
    }

    /**
     * Coarse gate. On a flat Network it tests the model against
     * {@link #getFeatureSet()}; on a LayeredNetwork each layer against
     * {@link #getLayeredFeatureSet()}.
     */
    @Override
    public boolean supports(Network model) {
        if (this.lnModel == null) {
            Network net = (model != null) ? model : this.model;
            if (net == null) {
                return false;
            }
            return FeatureSet.supports(getFeatureSet(), net.getUsedLangFeatures());
        }
        LayeredNetwork lqnModel = this.lnModel;

        // Get the feature set supported by LQNS
        FeatureSet featSupported = getLayeredFeatureSet();

        // Check if the model is an ensemble (has multiple layers)
        int numLayers = lqnModel.getNumberOfLayers();
        if (numLayers > 0) {
            // For layered networks with multiple layers, check each layer's features
            List<Network> ensemble = lqnModel.getEnsemble();
            if (ensemble != null) {
                for (Network layer : ensemble) {
                    FeatureSet featUsed = layer.getUsedLangFeatures();
                    if (!FeatureSet.supports(featSupported, featUsed)) {
                        return false;
                    }
                }
            }
            return true;
        }

        // For models without layers, since LayeredNetwork doesn't implement getUsedLangFeatures,
        // we assume it's supported if it's a valid LayeredNetwork instance
        return true;
    }

    /**
     * The per-method envelope: the Network feature set for the qns methods, and
     * none on a LayeredNetwork, whose gate is the coarse layered test.
     */
    @Override
    public FeatureSet getMethodFeatureSet(String method) {
        if (this.lnModel != null) {
            return null;
        }
        return getFeatureSet();
    }

    /**
     * Structural finite-capacity gate.
     *
     * <p>On a LayeredNetwork the gate refuses the qns names, which need a flat
     * Network, and otherwise applies the coarse layered test.
     *
     * <p>On a Network: NOTHING on the qnsolver path reads sn.cap or sn.classcap -- the model is
     * written out for {@code qnsolver}, whose MVA-family algorithms have no
     * representation of a finite buffer -- so a capped station was solved as an
     * unbounded one and the table reported the unconstrained answer under this
     * solver's name. There is no registry feature name for plain capacity, hence
     * the structural test; SolverMVA, SolverNC, SolverAG and SolverFluid gate
     * the same way through the same helper.
     *
     * <p>Without it SolverAUTO.listValidMethods offered all eight qns method names on
     * the BAS-blocking model of cqn_bas_blocking.
     *
     * @param method the concrete method name
     * @return empty string if supported, else the offending reason
     */
    @Override
    public String supportsModelMethod(String method) {
        if (this.lnModel != null) {
            if (isQnsMethod(method)) {
                return "SolverLQNS: the method '" + method + "' reaches qnsolver, which takes a flat "
                        + "Network; on a LayeredNetwork use default, lqns, srvn, exactmva, srvn.exactmva, "
                        + "sim, lqsim or lqns.default.";
            }
            return supports((Network) null) ? "" : "Some features are not supported by the chosen solver.";
        }
        // A BINDING FINITE BUFFER FIRST, ahead of the base gate. Nothing under
        // the qnsolver path reads sn.cap or sn.classcap -- neither the JMVA document
        // qnsolver reads nor the LQN QN2LQN writes has a buffer -- and this
        // refusal names the station, the cap and the way out, where the feature
        // envelope can only say "(feature: FiniteCapacity)". Once FiniteCapacity
        // became a registry name on 2026-09-05 the base gate started answering
        // first and the useful sentence became unreachable.
        if (this.model != null) {
            String capReason = NetworkSolver.bindingCapacityReason(this.model,
                    this.model.getStruct(false), "SolverLQNS");
            if (capReason != null && !capReason.isEmpty()) {
                return capReason;
            }
        }
        String reason = super.supportsModelMethod(method);
        if (reason != null && !reason.isEmpty()) {
            return reason;
        }
        if (this.model == null) {
            return "";
        }
        NetworkStruct snLocal = this.model.getStruct(false);
        String immfeed = qnsImmfeedRefusal(snLocal);
        if (!immfeed.isEmpty()) {
            return immfeed;
        }
        if (this.model.hasProductFormSolution() || this.model.hasOpenClasses()) {
            String ms = qnsMultiserverRefusal(snLocal, qnsMultiserver(method));
            if (!ms.isEmpty()) {
                return ms;
            }
            // "jmva" as the engine: SolverLQNS carries its own immediate-feedback
            // sentence below and must not be handed SolverJMT's.
            String jmva = jline.solvers.wrappers.jmt.SolverJMT.jmtMethodRefusal(snLocal, method, null, "jmva");
            if (jmva != null && !jmva.isEmpty()) {
                return jmva;
            }
        }
        return "";
    }

    /**
     * Why the qns methods of SolverLQNS cannot serve a model with immediate feedback, or "" when the
     * model has none.
     *
     * <p>Immediate feedback (sn.immfeed) keeps a self-looping job on its server
     * instead of re-queueing it, and neither qns path can state that:
     * the JMVA document qnsolver reads carries a mean demand and a visit count
     * per chain, and the LQN QN2LQN writes turns the routing into OR-fork
     * precedences of pseudo-activities on the reference task, where a repeated
     * visit is a new call. Either would answer for re-queueing under this
     * solver's name.
     *
     * <p>ONE PREDICATE, TWO CALLERS: {@link #supportsModelMethod} (the gate,
     * hence model.help and SolverAUTO) and {@link #runAnalyzer} (the run, for a
     * caller with enableChecks off). SolverJMT keeps its own wording in
     * jmtMethodRefusal. Mirrors matlab/src/solvers/wrappers/LQNS/qns_immfeed_refusal.m.
     *
     * @param sn the network struct
     * @return the refusal, or "" when the model carries no immediate feedback
     */
    public static String qnsImmfeedRefusal(NetworkStruct sn) {
        // SnHasImmfeed and not any(sn.immfeed): the CLASS spelling marks every
        // station, so the raw matrix REFUSES a model no self-loop can exercise.
        if (sn == null || !jline.api.sn.SnHasImmfeed.snHasImmfeed(sn)) {
            return "";
        }
        return "SolverLQNS does not support immediate feedback (sn.immfeed): neither the JMVA "
                + "document qnsolver reads nor the LQN QN2LQN writes can keep a self-looping job "
                + "on its server. Use SolverCTMC or SolverSSA, whose state space carries the "
                + "self-loop.";
    }

    /**
     * Whether qnsolver's own -m switch offers this multiserver approximation.
     *
     * <p>THE RULE IS INSIDE THE MULTISERVER BRANCH, and that is not a detail.
     * Without a multiserver station the reference emits no -m at all and answers
     * under the caller's method name, so refusing "suri" there would refuse a
     * model this solver does solve.
     *
     * <p>"qnsolver -m" accepts conway, reiser, rolia and zhou. "suri" and
     * "schmidt" are LQNS approximations, reachable only on the
     * non-product-form closed lqns branch, and qnsolver has no flag for
     * either. Mirrors matlab/src/solvers/wrappers/LQNS/qns_multiserver_refusal.m
     * and the C++ is_qnsolver_multiserver.
     *
     * @param sn     the network struct
     * @param method the multiserver approximation, i.e. {@link #qnsMultiserver} of the method
     * @return the refusal, or "" when the pair is served
     */
    public static String qnsMultiserverRefusal(NetworkStruct sn, String method) {
        if (sn == null || method == null || method.isEmpty()) {
            return "";
        }
        boolean multiserver = false;
        if (sn.nservers != null) {
            for (int i = 0; i < sn.nservers.length() && !multiserver; i++) {
                double c = sn.nservers.get(i);
                if (c > 1 && !Double.isInfinite(c)) {
                    multiserver = true;
                }
            }
        }
        if (!multiserver) {
            // No multiserver station, so no -m flag is emitted and every method
            // name is served by the plain invocation.
            return "";
        }
        String ms = method.toLowerCase();
        if (ms.equals("default") || ms.equals("conway") || ms.equals("reiser")
                || ms.equals("rolia") || ms.equals("zhou")) {
            return "";
        }
        return "SolverLQNS: the multiserver approximation '" + ms + "' is one LQNS offers and "
                + "qnsolver does not: 'qnsolver -m' accepts conway, reiser, rolia and zhou only; "
                + "suri and schmidt are available only on the non-product-form closed SolverLQNS "
                + "branch.";
    }

    /**
     * Get raw average tables with detailed metrics for nodes and calls.
     * This method provides more detailed metrics than getAvgTable(), including
     * phase-specific utilization and service times, processor metrics, and call waiting times.
     * Returns two tables aligned with MATLAB's getRawAvgTables:
     * - NodeAvgTable: Node, NodeType, Utilization, Phase1Utilization, Phase2Utilization,
     *                 Phase1ServiceTime, Phase2ServiceTime, Throughput, ProcWaiting, ProcUtilization
     * - CallAvgTable: SourceNode, TargetNode, Type, Waiting
     *
     * @return Array containing [NodeAvgTable, CallAvgTable] with detailed metrics
     */
    public LayeredNetworkAvgTable[] getRawAvgTables() {
        // Ensure we have results
        if (this.result == null) {
            try {
                this.runAnalyzer();
            } catch (IllegalAccessException | ParserConfigurationException e) {
                throw new RuntimeException("Failed to run analyzer for raw average tables", e);
            }
        }

        LayeredNetworkStruct lqn = this.getStruct();
        LayeredSolverResult layeredResult = (LayeredSolverResult) this.result;

        // Node metrics table
        List<String> nodeNames = new ArrayList<>(lqn.names.values());
        List<String> nodeTypes = new ArrayList<>();

        // Build node types list (aligned with MATLAB)
        for (int o = 0; o < nodeNames.size(); o++) {
            switch ((int) lqn.type.get(o)) {
                case LayeredNetworkElement.PROCESSOR:
                    nodeTypes.add("Processor");
                    break;
                case LayeredNetworkElement.TASK:
                    nodeTypes.add("Task");
                    break;
                case LayeredNetworkElement.ENTRY:
                    nodeTypes.add("Entry");
                    break;
                case LayeredNetworkElement.ACTIVITY:
                    nodeTypes.add("Activity");
                    break;
                case LayeredNetworkElement.CALL:
                    nodeTypes.add("Call");
                    break;
                default:
                    nodeTypes.add("Unknown");
                    break;
            }
        }

        // Extract raw metrics from stored results
        List<Double> utilization = layeredResult.rawUtilization.toList1D();
        List<Double> phase1Utilization = layeredResult.rawPhase1Utilization.toList1D();
        List<Double> phase2Utilization = layeredResult.rawPhase2Utilization.toList1D();
        List<Double> phase1ServiceTime = layeredResult.rawPhase1ServiceTime.toList1D();
        List<Double> phase2ServiceTime = layeredResult.rawPhase2ServiceTime.toList1D();
        List<Double> throughput = layeredResult.rawThroughput.toList1D();
        List<Double> procWaiting = layeredResult.rawProcWaiting.toList1D();
        List<Double> procUtilization = layeredResult.rawProcUtilization.toList1D();

        // Create enhanced node table with detailed metrics
        DetailedLayeredNetworkAvgTable nodeTable = new DetailedLayeredNetworkAvgTable(
                utilization, phase1Utilization, phase2Utilization,
                phase1ServiceTime, phase2ServiceTime, throughput,
                procWaiting, procUtilization
        );
        nodeTable.setNodeNames(nodeNames);
        nodeTable.setNodeTypes(nodeTypes);
        nodeTable.setOptions(this.options);

        // Call metrics table
        DetailedLayeredNetworkAvgTable callTable;
        if (lqn.ncalls == 0) {
            // Empty call table
            callTable = new DetailedLayeredNetworkAvgTable(
                    new ArrayList<>(), new ArrayList<>(), new ArrayList<>(),
                    new ArrayList<>(), new ArrayList<>(), new ArrayList<>(),
                    new ArrayList<>(), new ArrayList<>()
            );
        } else {
            // Extract call information from lqn structure
            List<String> sourceNodes = new ArrayList<>();
            List<String> targetNodes = new ArrayList<>();
            List<String> callTypeStrings = new ArrayList<>();
            List<Double> waitingTimes = layeredResult.rawEdgesWaiting.toList1D();

            // Build call table data (callpair rows are 0-based call indices)
            for (int i = 0; i < lqn.ncalls; i++) {
                int sourceIdx = (int) lqn.callpair.get(i, 0);
                int targetIdx = (int) lqn.callpair.get(i, 1);
                sourceNodes.add(lqn.names.get(sourceIdx));
                targetNodes.add(lqn.names.get(targetIdx));

                // Get call type from lqn.calltype
                CallType callType = lqn.calltype.get(i);
                if (callType == CallType.SYNC) {
                    callTypeStrings.add("Synchronous");
                } else if (callType == CallType.ASYNC) {
                    callTypeStrings.add("Asynchronous");
                } else if (callType == CallType.FWD) {
                    callTypeStrings.add("Forwarding");
                } else {
                    callTypeStrings.add("Unknown");
                }
            }

            callTable = new DetailedLayeredNetworkAvgTable(
                    waitingTimes, new ArrayList<>(), new ArrayList<>(),
                    new ArrayList<>(), new ArrayList<>(), new ArrayList<>(),
                    new ArrayList<>(), new ArrayList<>()
            );
            callTable.setCallData(sourceNodes, targetNodes, callTypeStrings);
        }
        callTable.setOptions(this.options);

        return new LayeredNetworkAvgTable[]{nodeTable, callTable};
    }

    /**
     * Enhanced LayeredNetworkAvgTable with detailed metrics support
     */
    public static class DetailedLayeredNetworkAvgTable extends LayeredNetworkAvgTable {
        private List<Double> phase1Utilization;
        private List<Double> phase2Utilization; 
        private List<Double> phase1ServiceTime;
        private List<Double> phase2ServiceTime;
        private List<Double> procWaiting;
        private List<Double> procUtilization;
        
        // For call tables
        private List<String> sourceNodes;
        private List<String> targetNodes;
        private List<String> callTypes;

        public DetailedLayeredNetworkAvgTable(
                List<Double> utilization,
                List<Double> phase1Utilization, 
                List<Double> phase2Utilization,
                List<Double> phase1ServiceTime,
                List<Double> phase2ServiceTime,
                List<Double> throughput,
                List<Double> procWaiting,
                List<Double> procUtilization) {
            super(new ArrayList<>(), utilization, new ArrayList<>(), 
                  procWaiting, new ArrayList<>(), throughput);
            this.phase1Utilization = phase1Utilization;
            this.phase2Utilization = phase2Utilization;
            this.phase1ServiceTime = phase1ServiceTime;
            this.phase2ServiceTime = phase2ServiceTime;
            this.procWaiting = procWaiting;
            this.procUtilization = procUtilization;
        }

        public void setCallData(List<String> sourceNodes, List<String> targetNodes, List<String> callTypes) {
            this.sourceNodes = sourceNodes;
            this.targetNodes = targetNodes;
            this.callTypes = callTypes;
        }

        public List<Double> getPhase1Utilization() { return phase1Utilization; }
        public List<Double> getPhase2Utilization() { return phase2Utilization; }
        public List<Double> getPhase1ServiceTime() { return phase1ServiceTime; }
        public List<Double> getPhase2ServiceTime() { return phase2ServiceTime; }
        public List<Double> getProcWaiting() { return procWaiting; }
        public List<Double> getProcUtilization() { return procUtilization; }
        public List<String> getSourceNodes() { return sourceNodes; }
        public List<String> getTargetNodes() { return targetNodes; }
        public List<String> getCallTypes() { return callTypes; }

        @Override
        public void print(SolverOptions options, boolean printZeros) {
            if (sourceNodes != null) {
                // Print call table
                if (options.verbose != VerboseLevel.SILENT) {
                    String[] headers = {"SourceNode", "TargetNode", "Type", "Waiting"};
                    List<String[]> rows = new ArrayList<>();
                    for (int i = 0; i < sourceNodes.size(); i++) {
                        rows.add(new String[]{
                                sourceNodes.get(i), targetNodes.get(i),
                                callTypes.get(i), formatValue(getQLen().get(i))});
                    }
                    printFormattedTable(headers, rows);
                }
            } else {
                // Print enhanced node table with detailed metrics
                if (options.verbose != VerboseLevel.SILENT) {
                    String[] headers = {"Node", "NodeType", "Util", "P1Util", "P2Util",
                            "P1SvcT", "P2SvcT", "Tput", "ProcWait", "ProcUtil"};
                    List<String[]> rows = new ArrayList<>();
                    for (int i = 0; i < getNodeTypes().size(); i++) {
                        if (printZeros || hasNonZeroValues(i)) {
                            rows.add(new String[]{
                                    getNodeNames().get(i), getNodeTypes().get(i),
                                    formatValue(getUtil().get(i)),
                                    formatValue(phase1Utilization.get(i)),
                                    formatValue(phase2Utilization.get(i)),
                                    formatValue(phase1ServiceTime.get(i)),
                                    formatValue(phase2ServiceTime.get(i)),
                                    formatValue(getTput().get(i)),
                                    formatValue(procWaiting.get(i)),
                                    formatValue(procUtilization.get(i))});
                        }
                    }
                    printFormattedTable(headers, rows);
                }
            }
        }

        private boolean hasNonZeroValues(int i) {
            return getUtil().get(i) > GlobalConstants.Zero ||
                   getTput().get(i) > GlobalConstants.Zero ||
                   procUtilization.get(i) > GlobalConstants.Zero;
        }

        private String formatValue(double value) {
            if (Double.isNaN(value) || value == 0.0) {
                return "0";
            }
            return fmtValue(value, 5);
        }
    }

    /**
     * Returns the default solver options for the LQNS solver.
     *
     * @return Default solver options with SolverType.LQNS
     */
    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.LQNS);
    }

    /**
     * The qns methods: solve a flat Network with qnsolver, or with lqns through
     * QN2LQN when the model is closed and not product-form.
     */
    private void runAnalyzerNetwork() {
        long startTime = System.nanoTime();

        if (this.model == null) {
            throw new RuntimeException("Model is not provided");
        }

        if (this.options == null) {
            this.options = networkDefaultOptions();
        }
        // Propagate solver verbose level to global
        GlobalConstants.Verbose = options.verbose;

        if (this.sn == null) {
            this.sn = this.model.getStruct(false);
        }
        // The gate's own sentence for a caller running with enableChecks off:
        // neither path can keep a self-looping job on its server.
        String immfeedReason = qnsImmfeedRefusal(this.sn);
        if (!immfeedReason.isEmpty()) {
            throw new RuntimeException(immfeedReason);
        }
        jline.io.InputOutput.line_ack(options.verbose, "QNS");
        line_debug(options.verbose, String.format("SolverLQNS qns starting: method=%s, multiserver=%s, nstations=%d, nclasses=%d",
            options.method, options.config.multiserver, sn.nstations, sn.nclasses));

        // "qns.NAME" selects the multiserver approximation NAME; "default" and
        // "qns" leave it to the default. The engine below speaks the bare name,
        // and the caller's method name is restored once it has run.
        String requested = this.options.method;
        String method = qnsMultiserver(requested);
        this.options.method = method;
        try {
            runQnsEngine(method, startTime);
        } finally {
            this.options.method = requested;
        }
        // Report the engine's label under the qns names: qns.NAME, default/qns.NAME
        String label = this.result.method;
        if (label != null) {
            String prefix = "";
            if (label.startsWith("default/")) {
                prefix = "default/";
                label = label.substring(prefix.length());
            }
            if (!label.equals("default")) {
                label = "qns." + label;
            }
            this.result.method = prefix + label;
        }
    }

    private void runQnsEngine(String method, long startTime) {
        switch (method) {
            case "conway":
                this.options.config.multiserver = "conway";
                break;
            case "rolia":
                this.options.config.multiserver = "rolia";
                break;
            case "zhou":
                this.options.config.multiserver = "zhou";
                break;
            case "suri":
                this.options.config.multiserver = "suri";
                break;
            case "reiser":
                this.options.config.multiserver = "reiser";
                break;
            case "schmidt":
                this.options.config.multiserver = "schmidt";
                break;
            case "default":
                this.options.config.multiserver = "rolia";
                break;
        }

        boolean isProductForm = snHasProductForm(sn);
        boolean isOpen = snHasOpenClasses(sn);

        if (isProductForm || isOpen) {
            // Product-form or open: use qnsolver directly (QN2LQN does not support Source/Sink)
            if (!Solver_qns_analyzer.isQNSolverAvailable()) {
                throw new RuntimeException("SolverLQNS qns methods require the external 'qnsolver' tool for product-form and open networks. "
                        + "Obtain it from its authors at http://www.sce.carleton.ca/rads/lqns/; LINE ships no copy and "
                        + "runs none from a container image.");
            }
            line_debug(options.verbose, "SolverLQNS qns: product-form or open model, using qnsolver directly");
            Solver_qns_analyzer analyzer = new Solver_qns_analyzer(this);
            SolverResult result = analyzer.runAnalyzer();

            this.setAvgResults(result.QN, result.UN, result.RN, result.TN,
                    result.AN, result.WN, result.CN, result.XN,
                    result.runtime, result.method, result.iter);
        } else {
            // Non-product-form closed: convert to LQN and solve via SolverLQNS
            line_debug(options.verbose, "SolverLQNS qns: non-product-form closed model, converting to LQN and solving with lqns");
            LayeredNetwork lqnmodel = QN2LQN.convert(this.model);
            SolverOptions lqnsoptions = SolverLQNS.defaultOptions();
            lqnsoptions.verbose = this.options.verbose;

            // Map multiserver method to LQNS options
            String actualMethod = method;
            switch (method) {
                case "conway":
                    lqnsoptions.config.multiserver = "conway";
                    break;
                case "rolia":
                    lqnsoptions.config.multiserver = "rolia";
                    break;
                case "zhou":
                    lqnsoptions.config.multiserver = "zhou";
                    break;
                case "suri":
                    lqnsoptions.config.multiserver = "suri";
                    break;
                case "reiser":
                    lqnsoptions.config.multiserver = "reiser";
                    break;
                case "schmidt":
                    lqnsoptions.config.multiserver = "schmidt";
                    break;
                case "default":
                    lqnsoptions.config.multiserver = "rolia";
                    actualMethod = "rolia";
                    break;
            }

            LayeredNetworkAvgTable avgTable = new SolverLQNS(lqnmodel, lqnsoptions).getLNAvgTable();
            LayeredNetworkStruct lqn = lqnmodel.getStruct();

            int M = this.sn.nstations;
            int K = this.sn.nclasses;
            Matrix QN = new Matrix(M, K);
            Matrix UN = new Matrix(M, K);
            Matrix RN = new Matrix(M, K);
            Matrix TN = new Matrix(M, K);
            Matrix WN = new Matrix(M, K);

            List<Double> qlenList = avgTable.getQLen();
            List<Double> utilList = avgTable.getUtil();
            List<Double> respTList = avgTable.getRespT();
            List<Double> residTList = avgTable.getResidT();
            List<Double> tputList = avgTable.getTput();

            // see _kb/12-interfaces-and-docs.md (Wrappers: JAR subprocess-bridge notes: LQNS entry-phase utilization fallback)
            int entryOffset = lqn.eshift + sn.nchains;

            for (int r = 0; r < K; r++) {
                for (int i = 0; i < M; i++) {
                    // MATLAB: t = lqn.ashift + r + (i-1)*nclasses (1-based)
                    // Java: ashift is 0-based count of hosts+tasks+entries, list is 0-based
                    int t = lqn.ashift + r + i * K;
                    int e = entryOffset + r + i * K;
                    if (t < qlenList.size()) {
                        QN.set(i, r, qlenList.get(t));
                        // Use entry-level Util/Tput when activity-level values are NaN
                        double util = utilList.get(t);
                        double tput = tputList.get(t);
                        if (Double.isNaN(util) && e < utilList.size()) {
                            util = utilList.get(e);
                        }
                        if (Double.isNaN(tput) && e < tputList.size()) {
                            tput = tputList.get(e);
                        }
                        // see _kb/12-interfaces-and-docs.md (Wrappers: JAR subprocess-bridge notes: LQNS entry-phase utilization fallback)
                        // lqns sums the utilization over the host's servers, a Network station reports it per server
                        double nservers = this.sn.nservers.get(i);
                        if (!Double.isInfinite(nservers) && nservers > 0) {
                            util = util / nservers;
                        }
                        UN.set(i, r, util);
                        RN.set(i, r, respTList.get(t));
                        WN.set(i, r, residTList.get(t));
                        TN.set(i, r, tput);
                    }
                }
            }

            // Compute arrival rates from throughputs
            Matrix AN = new Matrix(M, K);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) {
                    AN.set(i, r, TN.get(i, r));
                }
            }

            double runtime = (System.nanoTime() - startTime) / 1000000000.0;

            // Handle default method naming
            if ("default".equals(method)) {
                actualMethod = "default/" + actualMethod;
            }

            int C = this.sn.nchains;
            Matrix CN = new Matrix(1, C);
            Matrix XN = new Matrix(1, C);
            this.setAvgResults(QN, UN, RN, TN, AN, WN, CN, XN,
                    runtime, actualMethod, 0);
        }
    }

    // Probability methods - SolverLQNS does not support detailed state probability analysis

    @Override
    public ProbabilityResult getProbNormConstAggr() {
        throw new RuntimeException("getProbNormConstAggr not supported by SolverLQNS");
    }

    @Override
    public ProbabilityResult getProb(int node, Matrix state) {
        throw new RuntimeException("getProb not supported by SolverLQNS");
    }

    @Override
    public ProbabilityResult getProb(int node) {
        throw new RuntimeException("getProb not supported by SolverLQNS");
    }

    @Override
    public ProbabilityResult getProbSys() {
        throw new RuntimeException("getProbSys not supported by SolverLQNS");
    }

    @Override
    public ProbabilityResult getProbAggr(int node, Matrix state_a) {
        throw new RuntimeException("getProbAggr not supported by SolverLQNS");
    }

    @Override
    public ProbabilityResult getProbAggr(int node) {
        throw new RuntimeException("getProbAggr not supported by SolverLQNS");
    }

    @Override
    public ProbabilityResult getProbSysAggr() {
        throw new RuntimeException("getProbSysAggr not supported by SolverLQNS");
    }

    // Sampling methods - SolverLQNS does not support sampling

    @Override
    public SampleResult sample(int node, int numEvents) {
        throw new RuntimeException("sample not supported by SolverLQNS");
    }

    @Override
    public SampleResult sampleAggr(int node, int numEvents) {
        throw new RuntimeException("sampleAggr not supported by SolverLQNS");
    }

    @Override
    public SampleResult sampleSys(int numEvents) {
        throw new RuntimeException("sampleSys not supported by SolverLQNS");
    }

    @Override
    public SampleResult sampleSysAggr(int numEvents) {
        throw new RuntimeException("sampleSysAggr not supported by SolverLQNS");
    }

    // Distribution methods are NOT overridden: the reference SolverLQNS declares
    // none, so getCdfRespT falls through to the NetworkSolver exponential
    // fallback and the transient/passage getters to the base refusals, exactly
    // as MATLAB's inheritance resolves them.

}
