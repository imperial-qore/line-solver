/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.jmt;

import static jline.GlobalConstants.Inf;
import static jline.GlobalConstants.NegInf;
import static jline.io.InputOutputKt.line_warning;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.io.*;
import jline.lang.FeatureSet;
import jline.lang.JobClass;
import jline.lang.Metric;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.RoutingMatrix;
import jline.lang.constant.*;
import jline.lang.nodes.Cache;
import jline.lang.nodes.ClassSwitch;
import jline.lang.nodes.Logger;
import jline.lang.nodes.Node;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Station;
import jline.lang.nodes.StatefulNode;
import jline.lang.sections.Section;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.solvers.AvgHandle;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.io.Ret.DistributionResult;
import jline.io.Ret.ProbabilityResult;
import jline.io.Ret.SampleResult;
import jline.solvers.jmt.handlers.SaveHandlers;
import jline.util.RandomManager;
import jline.util.Utils;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.NodeList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

import static java.lang.Double.NaN;
import static jline.api.sn.SnGetArvRFromTputKt.snGetArvRFromTput;
import static jline.api.sn.SnGetDemandsChainKt.snGetDemandsChain;
import static jline.io.InputOutputKt.*;
import static jline.io.SysUtilsKt.jmtGetPath;
import static jline.util.Utils.isInf;

/**
 * Solver interface to the Java Modelling Tools (JMT) simulation engine.
 * 
 * <p>SolverJMT provides integration with the JMT discrete-event simulation toolkit
 * for analyzing queueing networks through simulation. JMT offers powerful simulation
 * capabilities for complex network topologies and general service distributions that
 * may not be analytically tractable.</p>
 * 
 * <p>Key JMT solver capabilities:
 * <ul>
 *   <li>Discrete-event simulation via JMT engine</li>
 *   <li>Complex network topology support (fork-join, finite capacity, etc.)</li>
 *   <li>General service and interarrival time distributions</li>
 *   <li>Statistical analysis with confidence intervals</li>
 *   <li>Transient and steady-state performance metrics</li>
 *   <li>Model export to JMT JSIMG format</li>
 * </ul>
 * </p>
 * 
 * <p><strong>Requirements:</strong> This solver requires JMT.jar to be available
 * in the classpath. The solver can operate with or without the external JMT GUI
 * application installed.</p>
 * 
 * @see JMTResult
 * @see JMTOptions
 * @since 1.0
 */
public class SolverJMT extends NetworkSolver {
    public static final String FILE_FORMAT = "jsimg";
    public static final String JSIMG_PATH = "";
    public static final String XSI_NO_NAMESPACE_SCHEMA_LOCATION = "Archive.xsd";
    private SaveHandlers saveHandlers;
    private String jmtPath;
    private String filePath;
    private String fileName;
    private String lastCommandOutput = "";
    private double maxSimulatedTime;
    private long maxSamples;
    private long maxEvents;
    private long simulationTimeoutSeconds;  // Timeout in seconds for JMT subprocess (0 = no timeout)
    private long seed;
    private double simConfInt;
    private double simMaxRelErr;
    private JMTResult jmtResult;

    public SolverJMT(Network model) {
        this(model, SolverJMT.defaultOptions());
    }

    public SolverJMT(Network model, SolverOptions options) {
        super(model, "SolverJMT", options);
        this.simConfInt = 0.99;
        this.simMaxRelErr = 0.03;
        this.maxEvents = -1;
        this.simulationTimeoutSeconds = 0;  // No timeout by default
        this.jmtPath = jmtGetPath();
        this.result = new JMTResult();
        // Initialize seed and maxSamples from options
        this.seed = options.seed;
        this.maxSamples = options.samples;
    }

    public SolverJMT(Network model, SolverOptions options, String jmtPath) {
        super(model, "SolverJMT", options);
        this.simConfInt = 0.99;
        this.simMaxRelErr = 0.03;
        this.maxEvents = -1;
        this.simulationTimeoutSeconds = 0;  // No timeout by default
        this.jmtPath = jmtGetPath(jmtPath);
        this.result = new JMTResult();
        // Initialize seed and maxSamples from options
        this.seed = options.seed;
        this.maxSamples = options.samples;
    }

    public SolverJMT(Network model, String jmtPath) {
        this(model, SolverJMT.defaultOptions(), jmtPath);
    }

    public SolverJMT(Network model, Object... varargin) {
        this(model, SolverJMT.defaultOptions());
        this.options = SolverJMT.parseOptions(this.options, varargin);
        this.result = new JMTResult();
    }

    private static boolean anyElementIsInfinity(double[] array) {
        for (double val : array) {
            if (isInf(val)) {
                return true;
            }
        }
        return false;
    }

    private static int countNodesWithType(List<NodeType> nodetypes, NodeType type) {
        int count = 0;
        for (NodeType nodeType : nodetypes) {
            if (nodeType == type) count++;
        }
        return count;
    }

    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.JMT);
    }

    public static FeatureSet getFeatureSet() {
        FeatureSet featSupported = new FeatureSet();
        featSupported.setTrue(new String[]{
                "Sink",
                "Source",
                "Router",
                "ClassSwitch",
                "Delay",
                "DelayStation",
                "Queue",
                "Fork",
                "Join",
                "Forker",
                "Joiner",
                "Logger",
                "Coxian",
                "Cox2",
                "APH",
                "Erlang",
                "Exp",
                "HyperExp",
                "Det",
                "Gamma",
                "Lognormal",
                "MAP",
                "MMPP2",
                "Normal",
                "PH",
                "Pareto",
                "Weibull",
                "Replayer",
                "Uniform",
                "StatelessClassSwitcher",
                "InfiniteServer",
                "SharedServer",
                "Buffer",
                "Dispatcher",
                "Server",
                "JobSink",
                "RandomSource",
                "ServiceTunnel",
                "LogTunnel",
                "Buffer",
                "Linkage",
                "Enabling",
                "Timing",
                "Firing",
                "Storage",
                "Place",
                "Transition",
                "SchedStrategy_INF",
                "SchedStrategy_PS",
                "SchedStrategy_DPS",
                "SchedStrategy_FCFS",
                "SchedStrategy_GPS",
                "SchedStrategy_SIRO",
                "SchedStrategy_HOL",
                "SchedStrategy_LCFS",
                "SchedStrategy_LCFSPR",
                "SchedStrategy_SEPT",
                "SchedStrategy_LEPT",
                "SchedStrategy_SJF",
                "SchedStrategy_LJF",
                "SchedStrategy_SRPT",
                "SchedStrategy_SRPTPRIO",
                "SchedStrategy_POLLING",
                "RoutingStrategy_PROB",
                "RoutingStrategy_RAND",
                "RoutingStrategy_RROBIN",
                "RoutingStrategy_WRROBIN",
                "RoutingStrategy_KCHOICES",
                "SchedStrategy_EXT",
                "ClosedClass", "SelfLoopingClass",
                "OpenClass"
        });
        return featSupported;
    }

    private static void setAlgTypeName(Element algTypeElement, Matrix nservers, String method, String name) {
        double maxFiniteValue = NegInf;  // initial value set to negative infinity

        for (int i = 0; i < nservers.getNumRows(); i++) {
            if (!Utils.isInf(nservers.get(i, 0))) {
                double currentValue = nservers.get(i, 0);
                if (currentValue > maxFiniteValue) {
                    maxFiniteValue = currentValue;
                }
            }
        }

        if (maxFiniteValue > 1) {
            throw new RuntimeException(method + " does not support multi-server stations.");
        } else {
            algTypeElement.setAttribute("name", name);
        }
    }

    public static void viewModel(String jmtPath, String filename, ViewMode viewMode) {
        viewModel(jmtPath, filename, viewMode, VerboseLevel.STD);
    }

    public static void viewModel(String jmtPath, String filename, ViewMode viewMode, VerboseLevel verboseLevel) {
        Path path = Paths.get(filename).getParent();
        if (path == null) {
            filename = Paths.get(java.lang.System.getProperty("user.dir"), filename).toString();
        }

        boolean suppressOutput = verboseLevel != VerboseLevel.DEBUG;
        String redirectOutput = "";
        if (suppressOutput) {
            redirectOutput = " > /dev/null 2>&1";
            if (java.lang.System.getProperty("os.name").startsWith("Windows")) {
                redirectOutput = " > nul 2>&1";
            }
        }

        String cmd = String.format(
                "java -cp %s jmt.commandline.Jmt %s %s %s",
                jmtPath, viewMode.toString().toLowerCase(), filename, redirectOutput
        );

        if (verboseLevel == VerboseLevel.DEBUG) {
            java.lang.System.out.println("JMT view model command: " + cmd);
        }
        
        if (verboseLevel == VerboseLevel.DEBUG) {
            String output = SysUtilsKt.system(cmd);
            if (!output.isEmpty()) {
                java.lang.System.out.println("JMT view model command output: " + output);
            }
        } else {
            // Suppress system command output unless in DEBUG mode
            java.io.ByteArrayOutputStream devNull = new java.io.ByteArrayOutputStream();
            java.io.PrintStream nullStream = new java.io.PrintStream(devNull);
            java.io.PrintStream originalOut = System.out;
            java.io.PrintStream originalErr = System.err;
            
            try {
                System.setOut(nullStream);
                System.setErr(nullStream);
                SysUtilsKt.system(cmd);
            } finally {
                System.setOut(originalOut);
                System.setErr(originalErr);
                nullStream.close();
            }
        }
    }

    public static void viewModel(String filename, ViewMode viewMode) {
        viewModel(jmtGetPath(), filename, viewMode);
    }

    public static String writeJMVA(NetworkStruct sn, String outputFileName, SolverOptions options) {
        try {
            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document mvaDoc = dBuilder.newDocument();

            Element mvaElem = mvaDoc.createElement("model");
            mvaDoc.appendChild(mvaElem);
            mvaElem.setAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
            mvaElem.setAttribute("xsi:noNamespaceSchemaLocation", "JMTmodel.xsd");

            Element algTypeElement = mvaDoc.createElement("algType");

            switch (options.method) {
                case "jmva.recal":
                    setAlgTypeName(algTypeElement, sn.nservers, options.method, "RECAL");
                    break;
                case "jmva.comom":
                    setAlgTypeName(algTypeElement, sn.nservers, options.method, "CoMoM");
                    break;
                case "jmva.chow":
                    setAlgTypeName(algTypeElement, sn.nservers, options.method, "Chow");
                    break;
                case "jmva.bs":
                case "jmva.amva":
                    setAlgTypeName(algTypeElement, sn.nservers, options.method, "Bard-Schweitzer");
                    break;
                case "jmva.aql":
                    setAlgTypeName(algTypeElement, sn.nservers, options.method, "AQL");
                    break;
                case "jmva.lin":
                    setAlgTypeName(algTypeElement, sn.nservers, options.method, "Linearizer");
                    break;
                case "jmva.dmlin":
                    setAlgTypeName(algTypeElement, sn.nservers, options.method, "De Souza-Muntz Linearizer");
                    break;
                //case "jmva.ls":
                //    algTypeElement.setAttribute("name", "Logistic Sampling");
                //    break;
                default:
                    algTypeElement.setAttribute("name", "MVA");
                    break;
            }

            algTypeElement.setAttribute("tolerance", "1.0E-7");
            algTypeElement.setAttribute("maxSamples", String.valueOf(options.samples));

            int M = sn.nstations;    //number of stations
            Matrix NK = sn.njobs.transpose();  // initial population per class
            int C = sn.nchains;
            Matrix SCV = sn.scv;
            Matrix ST = new Matrix(sn.rates.getNumRows(), sn.rates.getNumCols());

            for (int i = 0; i < sn.rates.getNumRows(); i++) {
                for (int j = 0; j < sn.rates.getNumCols(); j++) {

                    double currentRate = sn.rates.get(i, j);
                    if (Double.isNaN(currentRate)) {
                        ST.set(i, j, 0);
                    } else {
                        ST.set(i, j, 1.0 / currentRate);
                    }

                    if (Double.isNaN(SCV.get(i, j))) {
                        SCV.set(i, j, 1);
                    }
                }
            }

            Ret.snGetDemands snGetDemandsChainReturn = snGetDemandsChain(sn);

            Element parametersElem = mvaDoc.createElement("parameters");
            Element classesElem = mvaDoc.createElement("classes");
            classesElem.setAttribute("number", String.valueOf(sn.nchains));
            Element stationsElem = mvaDoc.createElement("stations");
            int numberOfStations = sn.nstations - countNodesWithType(sn.nodetype, NodeType.Source);
            stationsElem.setAttribute("number", String.valueOf(numberOfStations));
            Element refStationsElem = mvaDoc.createElement("ReferenceStation");
            refStationsElem.setAttribute("number", String.valueOf(sn.nchains));
            Element algParamsElem = mvaDoc.createElement("algParams");

            boolean[] sourceid = new boolean[sn.nodetype.size()];
            for (int i = 0; i < sn.nodetype.size(); i++) {
                sourceid[i] = sn.nodetype.get(i) == NodeType.Source;
            }

            for (int c = 0; c < sn.nchains; c++) {
                Element classElem;
                double sumOfNJobs = 0.0;
                // Check if chains matrix has enough rows
                if (c < sn.chains.getNumRows()) {
                    // Iterate over all classes (columns) to find which classes belong to this chain
                    for (int k = 0; k < sn.chains.getNumCols(); k++) {
                        // sn.chains.get(c, k) > 0 means class k belongs to chain c
                        if (sn.chains.get(c, k) > 0 && k < sn.njobs.length()) {
                            sumOfNJobs += sn.njobs.get(k);
                        }
                    }
                }
                if (!Utils.isInf(sumOfNJobs) && !Double.isNaN(sumOfNJobs)) {
                    classElem = mvaDoc.createElement("closedclass");
                    classElem.setAttribute("population", String.valueOf(snGetDemandsChainReturn.Nchain.get(c)));
                    classElem.setAttribute("name", String.format("Chain%02d", c + 1));
                } else {
                    double rateSum = 0.0;
                    for (int i = 0; i < sourceid.length; i++) {
                        if (sourceid[i]) {
                            // Check if chains matrix has enough rows
                            if (c < sn.chains.getNumRows()) {
                                // Iterate over all classes (columns) to find which classes belong to this chain
                                for (int k = 0; k < sn.chains.getNumCols(); k++) {
                                    // sn.chains.get(c, k) > 0 means class k belongs to chain c
                                    if (sn.chains.get(c, k) > 0 && k < sn.rates.getNumCols()) {
                                        rateSum += sn.rates.get(i, k);
                                    }
                                }
                            }
                        }
                    }
                    classElem = mvaDoc.createElement("openclass");
                    classElem.setAttribute("rate", String.valueOf(rateSum));
                    classElem.setAttribute("name", String.format("Chain%02d", c + 1));
                }
                classesElem.appendChild(classElem);
            }

            boolean[] isLoadDep = new boolean[sn.nstations];
            for (int i = 0; i < sn.nstations; i++) {
                Element statElem = null;
                NodeType currentNodeType = sn.nodetype.get((int) sn.stationToNode.get(i));
                switch (currentNodeType) {
                    case Delay:
                        statElem = mvaDoc.createElement("delaystation");
                        statElem.setAttribute("name", sn.nodenames.get((int) sn.stationToNode.get(i)));
                        break;
                    case Queue:
                        if (sn.nservers.get(i) == 1) {
                            isLoadDep[i] = false;
                            statElem = mvaDoc.createElement("listation");
                        } else {
                            isLoadDep[i] = true;
                            statElem = mvaDoc.createElement("ldstation");
                        }
                        statElem.setAttribute("name", sn.nodenames.get((int) sn.stationToNode.get(i)));
                        statElem.setAttribute("servers", "1");
                        break;
                    default:
                        continue;
                }

                Element srvTimesElem = mvaDoc.createElement("servicetimes");
                for (int c = 0; c < sn.nchains; c++) {
                    if (isLoadDep[i]) {
                        Element statSrvTimeElem = mvaDoc.createElement("servicetimes");
                        statSrvTimeElem.setAttribute("customerclass", String.format("Chain%02d", c + 1));
                        String ldSrvString = String.valueOf(snGetDemandsChainReturn.STchain.get(i, c));
                        // For open models (Inf population), use nservers as cutoff
                        // since service time is constant at S/c for n >= c
                        int ldLimit;
                        if (anyElementIsInfinity(NK.toArray1D())) {
                            ldLimit = (int) sn.nservers.get(i);
                        } else {
                            ldLimit = (int) Arrays.stream(NK.toArray1D()).sum();
                        }

                        for (int n = 2; n <= ldLimit; n++) {
                            ldSrvString = String.format("%s;%s", ldSrvString, snGetDemandsChainReturn.STchain.get(i, c) / FastMath.min(n, sn.nservers.get(i)));
                        }

                        statSrvTimeElem.appendChild(mvaDoc.createTextNode(ldSrvString));
                        srvTimesElem.appendChild(statSrvTimeElem);
                    } else {
                        Element statSrvTimeElem = mvaDoc.createElement("servicetime");
                        statSrvTimeElem.setAttribute("customerclass", String.format("Chain%02d", c + 1));
                        statSrvTimeElem.appendChild(mvaDoc.createTextNode(String.valueOf(snGetDemandsChainReturn.STchain.get(i, c))));
                        srvTimesElem.appendChild(statSrvTimeElem);
                    }
                }
                statElem.appendChild(srvTimesElem);
                Element visitsElem = mvaDoc.createElement("visits");
                for (int c = 0; c < sn.nchains; c++) {
                    Element statVisitElem = mvaDoc.createElement("visit");
                    statVisitElem.setAttribute("customerclass", String.format("Chain%02d", c + 1));

                    double val;
                    if (snGetDemandsChainReturn.STchain.get(i, c) > 0) {
                        val = snGetDemandsChainReturn.Dchain.get(i, c) / snGetDemandsChainReturn.STchain.get(i, c);
                    } else {
                        val = 0;
                    }

                    statVisitElem.appendChild(mvaDoc.createTextNode(String.valueOf(val)));
                    visitsElem.appendChild(statVisitElem);
                }
                statElem.appendChild(visitsElem);
                stationsElem.appendChild(statElem);
            }

            int[] refstatchain = new int[C];
            for (int c = 0; c < sn.nchains; c++) {
                Matrix inchain = sn.inchain.get(c);
                refstatchain[c] = (int) sn.refstat.get((int) inchain.get(0));
            }
            for (int c = 0; c < sn.nchains; c++) {
                Element classRefElem = mvaDoc.createElement("Class");
                classRefElem.setAttribute("name", String.format("Chain%02d", c + 1));
                // For open chains, the refstat is Source which is excluded from JMVA output.
                // Use the first non-Source station as the reference station instead.
                int refIdx = refstatchain[c];
                NodeType refNodeType = sn.nodetype.get((int) sn.stationToNode.get(refIdx));
                if (refNodeType == NodeType.Source) {
                    for (int i = 0; i < sn.nstations; i++) {
                        NodeType nt = sn.nodetype.get((int) sn.stationToNode.get(i));
                        if (nt != NodeType.Source && nt != NodeType.Sink) {
                            refIdx = i;
                            break;
                        }
                    }
                }
                classRefElem.setAttribute("refStation", sn.nodenames.get((int) sn.stationToNode.get(refIdx)));
                refStationsElem.appendChild(classRefElem);
            }

            Element compareAlgsElem = mvaDoc.createElement("compareAlgs");
            compareAlgsElem.setAttribute("value", "false");
            algParamsElem.appendChild(algTypeElement);
            algParamsElem.appendChild(compareAlgsElem);

            parametersElem.appendChild(classesElem);
            parametersElem.appendChild(stationsElem);
            parametersElem.appendChild(refStationsElem);
            mvaElem.appendChild(parametersElem);
            mvaElem.appendChild(algParamsElem);
            writeXML(outputFileName, mvaDoc);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return outputFileName;
    }


    private List<Integer> findInd(String term, List<String> list) {
        List<Integer> res = new ArrayList<>();
        for (int i = 0; i < list.size(); i++) {
            String nodeName = list.get(i);
            if (nodeName.equalsIgnoreCase(term)) {
                res.add(i);
            }
        }
        return res;
    }

    public DistributionResult getCdfRespT() {
        AvgHandle R = getAvgRespTHandles();
        return getCdfRespT(R);
    }

    public DistributionResult getCdfRespT(AvgHandle RH) {
        if (GlobalConstants.DummyMode) {
            DistributionResult result = new DistributionResult();
            NetworkStruct sn = this.getStruct();
            result.numStations = sn.nstations;
            result.numClasses = sn.nclasses;
            result.distributionType = "ResponseTime";
            result.isTransient = false;
            result.cdfData = new ArrayList<>();
            for (int i = 0; i < sn.nstations; i++) {
                List<Matrix> stationData = new ArrayList<>();
                for (int r = 0; r < sn.nclasses; r++) {
                    stationData.add(null);
                }
                result.cdfData.add(stationData);
            }
            return result;
        }
        
        NetworkStruct sn = this.getStruct();
        
        // Get steady-state queue lengths to initialize the CDF model
        Matrix QN = this.getAvgQLen();
        Matrix n = QN.copy();
        
        // Adjust job numbers based on network constraints (following dev version logic)
        for (int r = 0; r < sn.nclasses; r++) {
            if (Double.isInfinite(sn.njobs.get(r))) {
                // Open class - use floor of queue lengths
                for (int i = 0; i < sn.nstations; i++) {
                    n.set(i, r, Math.floor(QN.get(i, r)));
                }
            } else {
                // Closed class - ensure total population equals njobs
                for (int i = 0; i < sn.nstations; i++) {
                    n.set(i, r, Math.floor(QN.get(i, r)));
                }
                double totalJobs = n.sumCols(r);
                if (totalJobs < sn.njobs.get(r)) {
                    // Find bottleneck station (maxpos equivalent) and add remaining jobs
                    int maxIdx = 0;
                    double maxVal = n.get(0, r);
                    for (int i = 1; i < sn.nstations; i++) {
                        if (n.get(i, r) > maxVal) {
                            maxVal = n.get(i, r);
                            maxIdx = i;
                        }
                    }
                    n.set(maxIdx, r, n.get(maxIdx, r) + sn.njobs.get(r) - totalJobs);
                }
            }
        }
        
        try {
            // Copy model to avoid modifying original (following dev version pattern)
            Network cdfmodel = this.model.copy();
            cdfmodel.resetNetwork();
            cdfmodel.reset();
            
            // Set up logging configuration based on R handles (following dev version logic)
            boolean[][] isNodeClassLogged = new boolean[cdfmodel.getNumberOfNodes()][cdfmodel.getNumberOfClasses()];
            
            for (int i = 0; i < cdfmodel.getNumberOfStations(); i++) {
                for (int r = 0; r < cdfmodel.getNumberOfClasses(); r++) {
                    Station station = cdfmodel.getStations().get(i);
                    JobClass jobClass = cdfmodel.getJobClasses().get(r);
                    if (RH != null && RH.hasMetric(station, jobClass)) {
                        Metric metric = RH.get(station, jobClass);
                        if (metric != null && !metric.isDisabled) {
                            String stationName = station.getName();
                            int ni = this.model.getNodeIndex(stationName);
                            isNodeClassLogged[ni][r] = true;
                        }
                    }
                }
            }
            
            // Convert boolean[][] to boolean[] for isNodeLogged
            boolean[] isNodeLogged = new boolean[cdfmodel.getNumberOfNodes()];
            for (int i = 0; i < cdfmodel.getNumberOfNodes(); i++) {
                boolean nodeLogged = false;
                for (int r = 0; r < cdfmodel.getNumberOfClasses(); r++) {
                    if (isNodeClassLogged[i][r]) {
                        nodeLogged = true;
                        break;
                    }
                }
                isNodeLogged[i] = nodeLogged;
            }
            
            // Set up logging path
            String logPath = SysUtilsKt.lineTempName("jmt_cdf_logs");
            
            // Create routing matrix from the original routing information
            RoutingMatrix Plinked = new RoutingMatrix(cdfmodel, cdfmodel.getClasses(), cdfmodel.getNodes());

            // Copy the original routing values from sn.rtorig
            // Need to find corresponding classes in cdfmodel since it's a copy
            for (JobClass origFromClass : sn.rtorig.keySet()) {
                for (JobClass origToClass : sn.rtorig.get(origFromClass).keySet()) {
                    // Find corresponding classes in cdfmodel by index
                    JobClass fromClass = cdfmodel.getClasses().get(origFromClass.getIndex() - 1);
                    JobClass toClass = cdfmodel.getClasses().get(origToClass.getIndex() - 1);
                    Matrix routingMatrix = sn.rtorig.get(origFromClass).get(origToClass);
                    Plinked.set(fromClass, toClass, routingMatrix);
                }
            }

            // Link and log (following dev version pattern)
            cdfmodel.linkAndLog(Plinked, isNodeLogged, logPath);
            
            // For CDF analysis, ensure closed class jobs start at reference stations
            // This is critical for correct transient analysis
            // NOTE: After linkAndLog, we need to use the updated network structure
            NetworkStruct cdfSn = cdfmodel.getStruct(false);
            Matrix n_ref = new Matrix(cdfSn.nstations, cdfSn.nclasses);
            n_ref.zero();
            
            // Map original stations to new stations after linkAndLog
            for (int r = 0; r < cdfSn.nclasses; r++) {
                if (!Double.isInfinite(cdfSn.njobs.get(r))) { // closed class
                    // Find the reference station in the new network
                    String refStatName = this.model.getStations().get((int)sn.refstat.get(r, 0)).getName();
                    int newRefStat = -1;
                    for (int i = 0; i < cdfmodel.getNumberOfStations(); i++) {
                        if (cdfmodel.getStations().get(i).getName().equals(refStatName)) {
                            newRefStat = i;
                            break;
                        }
                    }
                    if (newRefStat >= 0) {
                        n_ref.set(newRefStat, r, cdfSn.njobs.get(r));
                    }
                } else { // open class - use the steady-state distribution
                    // Map original stations to new stations
                    for (int i = 0; i < sn.nstations; i++) {
                        String statName = this.model.getStations().get(i).getName();
                        for (int j = 0; j < cdfmodel.getNumberOfStations(); j++) {
                            if (cdfmodel.getStations().get(j).getName().equals(statName)) {
                                n_ref.set(j, r, n.get(i, r));
                                break;
                            }
                        }
                    }
                }
            }
            cdfmodel.initFromMarginal(n_ref);

            // Run simulation to generate log data
            SolverJMT cdfSolver = new SolverJMT(cdfmodel, this.options);
            cdfSolver.getAvg();

            // Parse logs following dev version approach
            // After linkAndLog, cdfmodel may have more nodes than the original
            // We need to map which nodes in the transformed model should be logged
            boolean[] cdfIsNodeLogged = new boolean[cdfmodel.getNumberOfNodes()];
            
            // isNodeLogged was created for the original cdfmodel before linkAndLog
            // We need to find which nodes in the transformed model correspond to logged nodes
            for (int i = 0; i < Math.min(isNodeLogged.length, this.model.getNumberOfNodes()); i++) {
                if (isNodeLogged[i]) {
                    String nodeName = this.model.getNodes().get(i).getName();
                    // Find this node in the transformed cdfmodel
                    for (int j = 0; j < cdfmodel.getNumberOfNodes(); j++) {
                        if (cdfmodel.getNodes().get(j).getName().equals(nodeName)) {
                            cdfIsNodeLogged[j] = true;
                            break;
                        }
                    }
                }
            }
            Matrix[][][] logData = parseLogs(cdfmodel, cdfIsNodeLogged, MetricType.RespT);
            
            // Create distribution result
            DistributionResult result = new DistributionResult();
            result.numStations = sn.nstations;
            result.numClasses = sn.nclasses;
            result.distributionType = "ResponseTime";
            result.isTransient = false;
            
            // Convert from nodes in logData to stations
            // The logData is indexed by cdfmodel nodes, not original model nodes
            List<List<Matrix>> cdfData = new ArrayList<>();
            List<Station> origStations = this.model.getStations();
            for (int i = 0; i < origStations.size(); i++) {
                List<Matrix> stationData = new ArrayList<>();
                // Find the corresponding node in cdfmodel
                String stationName = origStations.get(i).getName();
                int cdfNodeIndex = -1;
                for (int j = 0; j < cdfmodel.getNumberOfNodes(); j++) {
                    if (cdfmodel.getNodes().get(j).getName().equals(stationName)) {
                        cdfNodeIndex = j;
                        break;
                    }
                }
                
                for (int r = 0; r < this.model.getNumberOfClasses(); r++) {
                    if (cdfNodeIndex >= 0 && cdfNodeIndex < logData.length && 
                        cdfIsNodeLogged[cdfNodeIndex] && logData[cdfNodeIndex] != null && 
                        r < logData[cdfNodeIndex].length && 
                        logData[cdfNodeIndex][r] != null && logData[cdfNodeIndex][r].length > 0) {
                        Matrix respTimes = logData[cdfNodeIndex][r][0]; // Response times
                        if (respTimes != null && respTimes.getNumRows() > 0) {
                            // Create empirical CDF from response times (ecdf equivalent)
                            Matrix[] cdfArray = createEmpiricalCDF(respTimes);
                            if (cdfArray != null && cdfArray.length == 2) {
                                // Store as [F, X] like dev version
                                Matrix combinedCDF = new Matrix(cdfArray[0].getNumRows(), 2);
                                for (int k = 0; k < cdfArray[0].getNumRows(); k++) {
                                    combinedCDF.set(k, 0, cdfArray[0].get(k, 0)); // F values
                                    combinedCDF.set(k, 1, cdfArray[1].get(k, 0)); // X values
                                }
                                stationData.add(combinedCDF);
                            } else {
                                stationData.add(null);
                            }
                        } else {
                            stationData.add(null);
                        }
                    } else {
                        stationData.add(null);
                    }
                }
                cdfData.add(stationData);
            }
            
            result.cdfData = cdfData;
            return result;
            
        } catch (Exception e) {
            line_warning("SolverJMT.getCdfRespT", "Error: %s", e.getMessage());
            return null;
        }
    }

    private Matrix[] createEmpiricalCDF(Matrix data) {
        // Empirical CDF creation equivalent to MATLAB's ecdf function
        double[] array = data.toArray1D();
        if (array.length == 0) {
            return new Matrix[]{new Matrix(0, 1), new Matrix(0, 1)};
        }
        
        // Sort the data
        java.util.Arrays.sort(array);
        
        // Find unique values and their counts
        java.util.List<Double> uniqueValues = new java.util.ArrayList<>();
        java.util.List<Integer> counts = new java.util.ArrayList<>();
        
        double currentValue = array[0];
        int currentCount = 1;
        
        for (int i = 1; i < array.length; i++) {
            if (array[i] == currentValue) {
                currentCount++;
            } else {
                uniqueValues.add(currentValue);
                counts.add(currentCount);
                currentValue = array[i];
                currentCount = 1;
            }
        }
        // Add the last group
        uniqueValues.add(currentValue);
        counts.add(currentCount);
        
        int n = array.length;
        int numUnique = uniqueValues.size();
        Matrix F = new Matrix(numUnique + 1, 1); // CDF values (include 0 at start)
        Matrix X = new Matrix(numUnique + 1, 1); // Data points
        
        // Start with F(0) = 0 at the first data point
        F.set(0, 0, 0.0);
        X.set(0, 0, uniqueValues.get(0));
        
        int cumulativeCount = 0;
        for (int i = 0; i < numUnique; i++) {
            cumulativeCount += counts.get(i);
            F.set(i + 1, 0, (double) cumulativeCount / n);
            X.set(i + 1, 0, uniqueValues.get(i));
        }
        
        return new Matrix[]{F, X};
    }

    public String getFileName() {
        return fileName;
    }

    public void setFileName(String fileName) {
        this.fileName = fileName;
    }

    private SaveHandlers getSaveHandlers() {
        if (saveHandlers == null) {
            saveHandlers = new SaveHandlers(
                getModel(),
                getSimMaxRelErr(),
                getSimConfInt(),
                getAvgHandles(),
                getSeed(),
                getFileName(),
                getMaxEvents(),
                getMaxSamples(),
                getMaxSimulatedTime()
            );
        }
        return saveHandlers;
    }

    public String getFilePath() {
        return filePath;
    }

    public void setFilePath(String filePath) {
        this.filePath = filePath;
    }

    public String getJMVATempPath() {
        if (this.filePath == null || this.fileName == null) {
            try {
                this.filePath = SysUtilsKt.lineTempName("jmva");
            } catch (IOException ioe) {
                ioe.printStackTrace();
                throw new RuntimeException("Unable to get JMVA temp path");
            }
            this.fileName = "model";
        }
        String fname = this.fileName + ".jmva";
        return filePath + File.separator + fname;
    }

    public String getJSIMTempPath() {
        if (this.filePath == null || this.fileName == null) {
            try {
                this.filePath = SysUtilsKt.lineTempName("jsim");
            } catch (IOException ioe) {
                ioe.printStackTrace();
                throw new RuntimeException("Unable to get JSIM temp path");
            }
            this.fileName = "jmodel";
        }
        String fname = this.fileName + ".jsim";
        return filePath + File.separator + fname;
    }

    public String getJmtJarPath() {
        return jmtPath;
    }

    public void setJmtJarPath(String path) {
        this.jmtPath = path;
    }

    public long getMaxEvents() {
        return maxEvents;
    }

    public void setMaxEvents(long maxEvents) {
        this.maxEvents = maxEvents;
    }

    public long getMaxSamples() {
        return maxSamples;
    }

    public void setMaxSamples(long maxSamples) {
        this.maxSamples = maxSamples;
    }

    public double getMaxSimulatedTime() {
        return maxSimulatedTime;
    }

    public void setMaxSimulatedTime(double maxSimulatedTime) {
        this.maxSimulatedTime = maxSimulatedTime;
    }

    public long getSimulationTimeoutSeconds() {
        return simulationTimeoutSeconds;
    }

    public void setSimulationTimeoutSeconds(long timeoutSeconds) {
        this.simulationTimeoutSeconds = timeoutSeconds;
    }

    public double getProbAggr(Node node, Matrix state_a) {
        double Pr = NaN;
        if (GlobalConstants.DummyMode) {
            return Pr;
        }

        try {
            NetworkStruct sn = getStruct();
            SampleResult stationStateAggr = this.sampleAggr(node);

            // Validate sample result
            if (stationStateAggr == null) {
                line_warning(mfilename(new Object[]{}), "JMT getProbAggr: no sample result available, returning NaN.");
                return Pr;
            }

            // Get the state matrix from the sample result
            Matrix stateMatrix = stationStateAggr.getStateMatrix();

            // Validate state matrix
            if (stateMatrix == null || stateMatrix.getNumRows() == 0) {
                line_warning(mfilename(new Object[]{}), "JMT getProbAggr: empty state matrix, returning 0.");
                return 0.0;
            }

            // Validate state_a
            if (state_a == null || state_a.getNumRows() == 0) {
                line_warning(mfilename(new Object[]{}), "JMT getProbAggr: target state not set, returning 0.");
                return 0.0;
            }

            // Check dimension compatibility
            if (state_a.getNumCols() != stateMatrix.getNumCols()) {
                line_warning(mfilename(new Object[]{}), "JMT getProbAggr: state dimensions mismatch (target: " +
                    state_a.getNumCols() + " cols, samples: " + stateMatrix.getNumCols() + " cols), returning 0.");
                return 0.0;
            }

            // Find rows that match the requested state
            List<Integer> rows = Matrix.findRows(stateMatrix, state_a);

            // Get time points and calculate time differences
            Matrix t = stationStateAggr.t;
            if (t == null || t.getNumRows() == 0) {
                line_warning(mfilename(new Object[]{}), "JMT getProbAggr: no time data available, returning 0.");
                return 0.0;
            }

            Matrix dt = new Matrix(t.getNumRows(), 1);

            // Calculate time differences: dt = diff(t) with last element as 0
            for (int i = 0; i < t.getNumRows() - 1; i++) {
                dt.set(i, 0, t.get(i + 1, 0) - t.get(i, 0));
            }
            dt.set(t.getNumRows() - 1, 0, 0.0); // Last element is 0

            // Calculate probability as sum of time spent in matching states divided by total time
            double numerator = 0.0;
            double denominator = 0.0;

            for (int i = 0; i < dt.getNumRows(); i++) {
                double timeSpent = dt.get(i, 0);
                denominator += timeSpent;
                if (rows.contains(i)) {
                    numerator += timeSpent;
                }
            }

            if (denominator > 0) {
                Pr = numerator / denominator;
            } else {
                Pr = 0.0;
            }

        } catch (IOException e) {
            line_error(mfilename(new Object[]{}), "IOException in getProbAggr(): " + e.getMessage());
        } catch (Exception e) {
            line_error(mfilename(new Object[]{}), "Exception in getProbAggr(): " + e.getMessage());
        }
        return Pr;
    }

    public double getProbAggr(Node node) {
        double Pr = NaN;
        if (GlobalConstants.DummyMode) {
            return Pr;
        }
        try {
            Matrix state_a = sn.state.get(this.model.getStations().get((int) sn.stationToStateful.get((int) sn.nodeToStation.get(node.getNodeIndex()))));
            return getProbAggr(node, state_a);
        } catch (Exception e) {
            line_error(mfilename(new Object[]{}), "Exception in getProbAggr(): " + e.getMessage());
            return Pr;
        }
    }

    /*
     * This method returns the normalizing constant of the joint distribution of the system state and the aggregated state.
     * The normalizing constant is the sum of the probabilities of all possible states of the system.
     * @return the normalizing constant of the joint distribution of the system state and the aggregated state
     * @throws ParserConfigurationException
     */
    public ProbabilityResult getProbNormConstAggr() {
        Double lNormConst;
        if (GlobalConstants.DummyMode) {
            lNormConst = NaN;
            return new ProbabilityResult(lNormConst, true);
        }
        //case "jmva.ls":
        try {
            switch (options.method) {
                case "jmva":
                case "jmva.recal":
                case "jmva.comom":
                    this.runAnalyzer();
                    lNormConst = this.jmtResult.logNormConstAggr;
                    break;
                default:
                    lNormConst = NaN;
                    line_error(mfilename(new Object[]{}), "Selected solver method does not compute normalizing constants. Choose either jmva.recal, jmva.comom, or jmva.ls.");
            }
        } catch (Exception e) {
            lNormConst = NaN;
        }
        return new ProbabilityResult(lNormConst, true);
    }

    public SolverResult getResults() {
        SolverOptions options = this.options;
        SolverResult solverResult = new JMTResult();
        switch (options.method) {
            case "jsim":
            case "default":
                this.jmtResult = getResultsJSIM();
                break;
            default:
                this.jmtResult = getResultsJMVA();
                break;
        }


        NetworkStruct sn = getStruct();
        int numOfNodes = sn.nnodes;
        int numOfCache = 0;
        List<Cache> caches = new ArrayList<>();
        List<Integer> cacheNodeIndices = new ArrayList<>();
        for (int r = 0; r < numOfNodes; r++) {
            Node node = model.getNodes().get(r);
            if (node instanceof Cache) {
                numOfCache++;
                caches.add((Cache) node);
                cacheNodeIndices.add(r);
            }
        }

        solverResult.QN = new Matrix(sn.nstations, sn.nclasses);
        solverResult.UN = new Matrix(sn.nstations, sn.nclasses);
        solverResult.RN = new Matrix(sn.nstations, sn.nclasses);
        solverResult.TN = new Matrix(sn.nstations, sn.nclasses);
        Matrix cacheTN = new Matrix(numOfCache, sn.nclasses);
        // Node-indexed matrices for cache metrics (will be stored in JMTResult)
        Matrix nodeCacheTN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix nodeCacheAN = new Matrix(sn.nnodes, sn.nclasses);
        solverResult.AN = new Matrix(sn.nstations, sn.nclasses);
        solverResult.WN = new Matrix(sn.nstations, sn.nclasses);
        solverResult.TardN = new Matrix(sn.nstations, sn.nclasses);
        solverResult.SysTardN = new Matrix(1, sn.nclasses);
        solverResult.CN = new Matrix(1, sn.nchains);  // Mean system response times per chain
        solverResult.XN = new Matrix(1, sn.nchains);  // Mean system throughputs per chain

        // Initialize FCR metric matrices
        if (sn.nregions > 0) {
            solverResult.QNfcr = new Matrix(sn.nregions, sn.nclasses);
            solverResult.UNfcr = new Matrix(sn.nregions, sn.nclasses);
            solverResult.RNfcr = new Matrix(sn.nregions, sn.nclasses);
            solverResult.TNfcr = new Matrix(sn.nregions, sn.nclasses);
            solverResult.ANfcr = new Matrix(sn.nregions, sn.nclasses);
            solverResult.WNfcr = new Matrix(sn.nregions, sn.nclasses);
            // Initialize with zeros (will be populated from JMT results)
            solverResult.QNfcr.fill(0.0);
            solverResult.UNfcr.fill(Double.NaN);  // JMT doesn't provide Util for regions
            solverResult.RNfcr.fill(0.0);
            solverResult.TNfcr.fill(0.0);
            solverResult.ANfcr.fill(Double.NaN);  // JMT doesn't provide ArvR for regions
            solverResult.WNfcr.fill(0.0);
        }

        for (int m = 0; m < this.jmtResult.metrics.size(); m++) {
            Metric metric = this.jmtResult.metrics.get(m);
            String stationName = metric.getStationName();

            // Check if this is a region (FCR) metric
            if ("region".equals(metric.getNodeType())) {
                // Parse FCR index from station name "FCRegion1", "FCRegion2", etc.
                if (stationName != null && stationName.startsWith("FCRegion")) {
                    try {
                        int fcrIndex = Integer.parseInt(stationName.substring(8)) - 1;
                        if (fcrIndex >= 0 && fcrIndex < sn.nregions) {
                            // FCR metrics are aggregate (not per-class), so apply to all classes
                            switch (metric.getMetricType()) {
                                case QLen:
                                    for (int r = 0; r < sn.nclasses; r++) {
                                        solverResult.QNfcr.set(fcrIndex, r, metric.getMeanValue() / sn.nclasses);
                                    }
                                    break;
                                case RespT:
                                    for (int r = 0; r < sn.nclasses; r++) {
                                        solverResult.RNfcr.set(fcrIndex, r, metric.getMeanValue());
                                    }
                                    break;
                                case ResidT:
                                    for (int r = 0; r < sn.nclasses; r++) {
                                        solverResult.WNfcr.set(fcrIndex, r, metric.getMeanValue());
                                    }
                                    break;
                                case Tput:
                                    for (int r = 0; r < sn.nclasses; r++) {
                                        solverResult.TNfcr.set(fcrIndex, r, metric.getMeanValue() / sn.nclasses);
                                    }
                                    break;
                                default:
                                    break;
                            }
                        }
                    } catch (NumberFormatException e) {
                        // Ignore malformed FCR names
                    }
                }
                continue;  // Skip further processing for region metrics
            }
            List<Integer> indList = findInd(stationName, sn.nodenames);
            List<Integer> istStations = new ArrayList<>();
            for (int ind : indList) {
                istStations.add((int) sn.nodeToStation.get(ind));
            }
            List<Integer> rList = new ArrayList<>();
            String metricClass = metric.getClassName();
            
            // Check if this is a JMVA result (chain-based)
            if (metricClass != null && metricClass.startsWith("Chain")) {
                // Extract chain number from "Chain01", "Chain02", etc.
                try {
                    int chainIdx = Integer.parseInt(metricClass.substring(5)) - 1;
                    if (chainIdx >= 0 && chainIdx < sn.nchains) {
                        // Get all classes in this chain
                        Matrix inchain = sn.inchain.get(chainIdx);
                        for (int i = 0; i < inchain.length(); i++) {
                            rList.add((int) inchain.get(i));
                        }
                    }
                } catch (NumberFormatException | StringIndexOutOfBoundsException e) {
                    // Fall back to exact match
                    for (int i = 0; i < sn.classnames.size(); i++) {
                        String className = sn.classnames.get(i);
                        if (metricClass.equalsIgnoreCase(className)) {
                            rList.add(i);
                        }
                    }
                }
            } else {
                // Regular class-based matching
                for (int i = 0; i < sn.classnames.size(); i++) {
                    String className = sn.classnames.get(i);
                    if (metricClass != null && metricClass.equalsIgnoreCase(className)) {
                        rList.add(i);
                    }
                }
            }
            boolean open = true;
            for (int r : rList) {
                if (!Utils.isInf(sn.njobs.get(r)) && !Double.isNaN(sn.njobs.get(r))) {
                    open = false;
                }
            }
            int sumJobs = 0;
            for (int r : rList) {
                sumJobs += sn.njobs.get(r);
            }
            
            // Debug: check if rList is empty for JMVA
            if (rList.isEmpty()) {
                line_warning("SolverJMT", "rList is empty for metric class: %s", metricClass);
            }

            // For closed classes, filter metrics with insufficient analyzed samples (matches MATLAB getResults.m)
            // A class is considered recurrent only if analyzedSamples > total jobs in its chain
            if (!open && metric.getAnalyzedSamples() <= sumJobs) {
                continue;  // Leave as 0 (starved class)
            }

            switch (metric.getMetricType()) {
                case QLen:
                    solverResult.QN = setValues(solverResult.QN, istStations, rList, metric.getMeanValue());
                    break;
                case Util:
                    solverResult.UN = setValues(solverResult.UN, istStations, rList, metric.getMeanValue());
                    break;
                case RespT:
                    solverResult.RN = setValues(solverResult.RN, istStations, rList, metric.getMeanValue());
                    break;
                case ResidT:
                    solverResult.WN = setValues(solverResult.WN, istStations, rList, metric.getMeanValue());
                    break;
                case ArvR:
                    if (istStations.get(0) >= 0) {
                        // Regular station
                        solverResult.AN = setValues(solverResult.AN, istStations, rList, metric.getMeanValue());
                    } else {
                        // This is a cache node - store in nodeCacheAN using the actual node index
                        // Validate bounds before setting
                        if (!indList.isEmpty() && indList.get(0) < sn.nnodes) {
                            nodeCacheAN = setValues(nodeCacheAN, indList, rList, metric.getMeanValue());
                        }
                    }
                    break;
                case Tput:
                    if (istStations.get(0) >= 0) {
                        // Regular station
                        solverResult.TN = setValues(solverResult.TN, istStations, rList, metric.getMeanValue());
                    } else {
                        // This is a cache node - store in both cacheTN and nodeCacheTN
                        // Validate bounds before setting
                        if (!indList.isEmpty() && indList.get(0) < sn.nnodes) {
                            int cacheIdx = 0;
                            for (int idx = 0; idx < cacheNodeIndices.size(); idx++) {
                                if (cacheNodeIndices.get(idx).equals(indList.get(0))) {
                                    cacheIdx = idx;
                                    break;
                                }
                            }
                            List<Integer> cacheIdxList = new ArrayList<>();
                            cacheIdxList.add(cacheIdx);
                            cacheTN = setValues(cacheTN, cacheIdxList, rList, metric.getMeanValue());
                            // Also store in node-indexed matrix using the actual node index
                            nodeCacheTN = setValues(nodeCacheTN, indList, rList, metric.getMeanValue());
                        }
                    }
                    break;
                case Tard:
                    solverResult.TardN = setValues(solverResult.TardN, istStations, rList, metric.getMeanValue());
                    break;
                case SysTard:
                    // System tardiness is a 1 x classes metric (no station index)
                    for (int r : rList) {
                        solverResult.SysTardN.set(0, r, metric.getMeanValue());
                    }
                    break;
            }
        }

        Matrix hitProb = Matrix.zeros(numOfCache, sn.nclasses);
        for (int i = 0; i < numOfCache; i++) {
            for (int j = 0; j < sn.nclasses / 3; j++) {
                double total = solverResult.TN.get(i, 3 * j);
                if (j == 0) {
                    hitProb.set(i, 0, cacheTN.get(i, 1) / total);
                } else {
                    double hit = cacheTN.get(i, 3 * j + 1);
                    hitProb.set(i, j + 2, cacheTN.get(i, 3 * j + 1) / total);
                }
            }
            caches.get(i).setResultHitProb(hitProb);
        }


        this.result = solverResult;
        solverResult.method = this.options.method != null ? this.options.method : "default";

        // Store cache metrics in JMTResult for use in getAvgNode()
        if (solverResult instanceof JMTResult) {
            JMTResult jmtResult = (JMTResult) solverResult;
            jmtResult.cacheTN = nodeCacheTN;
            jmtResult.cacheAN = nodeCacheAN;
            jmtResult.cacheNodeIndices = cacheNodeIndices;
        }

        return solverResult;
    }

    public JMTResult getResultsJMVA() {
        JMTResult result = new JMTResult();
        String fileName = this.getFileName() + ".jmva-result.jmva";
        Path filePath = Paths.get(getFilePath(), fileName);
        File file = new File(filePath.toUri());
        if (file.exists() && !file.isDirectory()) {
            DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
            try {
                DocumentBuilder builder = factory.newDocumentBuilder();
                Document doc = builder.parse(file);
                
                // Get normalizing constant if available
                NodeList normConstList = doc.getElementsByTagName("normconst");
                if (normConstList.getLength() > 0) {
                    org.w3c.dom.Node normConstNode = normConstList.item(0);
                    if (normConstNode.getAttributes().getNamedItem("logValue") != null) {
                        String logValueStr = normConstNode.getAttributes().getNamedItem("logValue").getNodeValue();
                        result.logNormConstAggr = Double.parseDouble(logValueStr);
                    }
                }
                
                // Parse station results
                NodeList stationResults = doc.getElementsByTagName("stationresults");
                for (int i = 0; i < stationResults.getLength(); i++) {
                    org.w3c.dom.Node stationNode = stationResults.item(i);
                    String stationName = stationNode.getAttributes().getNamedItem("station").getNodeValue();
                    
                    // Parse class results within each station
                    NodeList classResults = ((Element)stationNode).getElementsByTagName("classresults");
                    for (int j = 0; j < classResults.getLength(); j++) {
                        org.w3c.dom.Node classNode = classResults.item(j);
                        String customerClass = classNode.getAttributes().getNamedItem("customerclass").getNodeValue();
                        
                        // Parse measures within each class
                        NodeList measures = ((Element)classNode).getElementsByTagName("measure");
                        for (int k = 0; k < measures.getLength(); k++) {
                            NamedNodeMap attributes = measures.item(k).getAttributes();
                            Metric metric = new Metric(attributes);
                            // Set the station and class from parent nodes
                            metric.setStationName(stationName);
                            metric.setClassName(customerClass);
                            result.metrics.add(metric);
                        }
                    }
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        } else {
            String errorMsg = "JMT did not output a result file, the analysis has likely failed.";
            if (this.lastCommandOutput != null && !this.lastCommandOutput.isEmpty()) {
                errorMsg += " JMT output: " + this.lastCommandOutput;
            }
            line_error(mfilename(new Object[]{}), errorMsg);
        }
        return result;
    }

    public JMTResult getResultsJSIM() {
        JMTResult result = new JMTResult();
        String fileName = this.getFileName() + ".jsim-result.jsim";
        Path filePath = Paths.get(getFilePath(), fileName);
        File file = new File(filePath.toUri());
        if (file.exists() && !file.isDirectory()) {
            DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
            try {
                DocumentBuilder builder = factory.newDocumentBuilder();
                Document doc = builder.parse(file);
                NodeList measure = doc.getElementsByTagName("measure");
                for (int i = 0; i < measure.getLength(); i++) {
                    NamedNodeMap attributes = measure.item(i).getAttributes();
                    Metric metric = new Metric(attributes);
                    result.metrics.add(metric);
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        } else {
            String errorMsg = "JMT did not output a result file, the simulation has likely failed.";
            if (this.lastCommandOutput != null && !this.lastCommandOutput.isEmpty()) {
                errorMsg += " JMT output: " + this.lastCommandOutput;
            }
            line_error(mfilename(new Object[]{}), errorMsg);
        }
        return result;
    }

    /**
     * Computes average performance metrics at steady-state for all nodes.
     * This method overrides NetworkSolver.getAvgNode() to use JMT simulation values
     * for cache node throughputs and arrival rates instead of computing them from routing probabilities.
     *
     * @return solver result containing node-level average metrics
     */
    @Override
    public SolverResult getAvgNode() {
        // Call parent implementation first
        SolverResult noderesult = super.getAvgNode();

        // If we have cached metrics from JMT simulation, use those for cache nodes
        if (this.result instanceof JMTResult) {
            JMTResult jmtResult = (JMTResult) this.result;
            if (jmtResult.cacheNodeIndices != null && !jmtResult.cacheNodeIndices.isEmpty()) {
                // Override the computed cache throughputs with the actual JMT values
                if (jmtResult.cacheTN != null) {
                    for (int cacheIdx : jmtResult.cacheNodeIndices) {
                        for (int r = 0; r < sn.nclasses; r++) {
                            double cachedValue = jmtResult.cacheTN.get(cacheIdx, r);
                            if (cachedValue > 0) {
                                noderesult.TN.set(cacheIdx, r, cachedValue);
                            }
                        }
                    }
                }
                // Override the computed cache arrival rates with the actual JMT values
                if (jmtResult.cacheAN != null) {
                    for (int cacheIdx : jmtResult.cacheNodeIndices) {
                        for (int r = 0; r < sn.nclasses; r++) {
                            double cachedValue = jmtResult.cacheAN.get(cacheIdx, r);
                            if (cachedValue > 0) {
                                noderesult.AN.set(cacheIdx, r, cachedValue);
                            }
                        }
                    }
                }

                // For cache models, fix arrival rates at Sink and ClassSwitch nodes
                // These arrival rates should equal cache output throughputs, not routing-based computation
                for (int ind = 0; ind < sn.nnodes; ind++) {
                    Node node = model.getNodes().get(ind);
                    if (node instanceof Sink || node instanceof ClassSwitch) {
                        // Get the cache node to read throughputs from
                        for (int cacheIdx : jmtResult.cacheNodeIndices) {
                            Cache cacheNode = (Cache) model.getNodes().get(cacheIdx);
                            // For each class, if it's a hit or miss class, set arrival rate = cache throughput
                            // For init class (non-hit/miss), set arrival rate to 0
                            for (int r = 0; r < sn.nclasses; r++) {
                                // Check if this class is a hit class or miss class from the cache
                                boolean isHitOrMissClass = false;
                                for (int col = 0; col < cacheNode.getHitClass().getNumCols(); col++) {
                                    if ((int) cacheNode.getHitClass().get(col) == r ||
                                        (int) cacheNode.getMissClass().get(col) == r) {
                                        isHitOrMissClass = true;
                                        break;
                                    }
                                }
                                if (isHitOrMissClass) {
                                    // Arrival rate at Sink/ClassSwitch = throughput at Cache for this class
                                    noderesult.AN.set(ind, r, noderesult.TN.get(cacheIdx, r));
                                } else if (node instanceof Sink) {
                                    // InitClass doesn't arrive at Sink (gets converted at Cache)
                                    noderesult.AN.set(ind, r, 0.0);
                                }
                                // For ClassSwitch nodes, InitClass keeps its computed arrival rate
                            }
                        }
                    }
                }
            }
        }

        return noderesult;
    }

    public long getSeed() {
        return seed;
    }

    public void setSeed(int seed) {
        this.seed = seed;
    }

    public double getSimConfInt() {
        return simConfInt;
    }

    public void setSimConfInt(double simConfInt) {
        this.simConfInt = simConfInt;
    }

    public double getSimMaxRelErr() {
        return simMaxRelErr;
    }

    public void setSimMaxRelErr(double simMaxRelErr) {
        this.simMaxRelErr = simMaxRelErr;
    }

    public NetworkStruct getStruct() {
        if (this.sn == null)
            this.sn = this.model.getStruct(true);
        return this.sn;
    }

    public DistributionResult getTranCdfPassT() {
        NetworkStruct sn = getStruct();
        Matrix RD = new Matrix(sn.nstations, sn.nclasses);
        if (GlobalConstants.DummyMode) {
            return new DistributionResult(sn.nstations, sn.nclasses, "passage_time");
        }
        AvgHandle R = getAvgRespTHandles();
        return getTranCdfPassT(R);
    }

    public DistributionResult getTranCdfPassT(AvgHandle R) {
        if (GlobalConstants.DummyMode) {
            return new DistributionResult();
        }

        NetworkStruct sn = this.model.getStruct(false);
        
        if (R == null) {
            R = getAvgRespTHandles();
        }
        
        DistributionResult result = new DistributionResult(sn.nstations, sn.nclasses, "passage_time");
        
        // Determine which nodes and classes need to be logged
        boolean[][] isNodeClassLogged = new boolean[this.model.getNumberOfNodes()][this.model.getNumberOfClasses()];
        for (int i = 0; i < this.model.getNumberOfStations(); i++) {
            for (int r = 0; r < this.model.getNumberOfClasses(); r++) {
                Station station = this.model.getStations().get(i);
                JobClass jobClass = this.model.getJobClasses().get(r);
                if (R.hasMetric(station, jobClass)) {
                    Metric metric = R.get(station, jobClass);
                    if (metric != null && !metric.isDisabled) {
                        int ni = this.model.getNodeIndex(station);
                        isNodeClassLogged[ni][r] = true;
                    }
                }
            }
        }
        
        // Determine which nodes need to be logged
        boolean[] isNodeLogged = new boolean[this.model.getNumberOfNodes()];
        for (int ni = 0; ni < this.model.getNumberOfNodes(); ni++) {
            for (int r = 0; r < this.model.getNumberOfClasses(); r++) {
                if (isNodeClassLogged[ni][r]) {
                    isNodeLogged[ni] = true;
                    break;
                }
            }
        }
        
        try {
            // Set up logging path
            String logPath = SysUtilsKt.lineTempName("jmt_logs");
            
            // Link and log the model
            RoutingMatrix P = this.model.initRoutingMatrix();
            this.model.linkAndLog(P, isNodeLogged, logPath);
            
            // Run simulation to generate log data
            this.getAvg();
            
            // Parse logs and compute CDFs (using same parsing as response time)
            Matrix[][][] RD = parseLogs(this.model, isNodeLogged, MetricType.RespT);
            
            // Convert from nodes to stations and compute empirical CDFs
            for (int i = 0; i < this.model.getNumberOfStations(); i++) {
                int ni = this.model.getNodeIndex(this.model.getStations().get(i));
                for (int r = 0; r < this.model.getNumberOfClasses(); r++) {
                    if (isNodeClassLogged[ni][r] && RD[ni][r] != null) {
                        // Compute empirical CDF using the response time data
                        double[] passageTimes = RD[ni][r][0].toArray1D();
                        if (passageTimes.length > 0) {
                            Matrix cdf = computeEmpiricalCDF(passageTimes);
                            result.setCdf(i, r, cdf);
                        }
                    }
                }
            }
            
        } catch (Exception e) {
            throw new RuntimeException("Error in getTranCdfPassT: " + e.getMessage(), e);
        }
        
        return result;
//        int numberOfClasses = cdfmodel.getNumberOfClasses();
//        boolean[][] isNodeClassLogged = new boolean[numberOfNodes][numberOfClasses];
//        for (int i = 0; i < numberOfNodes; i++) {
//            for (int r = 0; r < numberOfClasses; r++) {
//                isNodeClassLogged[i][r] = false;
//            }
//        }
//        for (int i = 0; i < numberOfNodes; i++) {
//            for (int r = 0; r < numberOfClasses; r++) {
//                if (!R.get(sn.stations.get((int) sn.nodeToStation.get(i))).get(sn.jobclasses.get(r)).isDisabled){
//                    int ni = this.model.getNodeIndex(cdfmodel.getStationFromIndex(i));
//                    isNodeClassLogged[ni][r] = true;
//                }
//            }
//        }
    }

    public DistributionResult getTranCdfRespT() {
        AvgHandle R = getAvgRespTHandles();
        return getCdfRespT(R);
    }

    public DistributionResult getTranCdfRespT(AvgHandle R) {
        if (GlobalConstants.DummyMode) {
            return new DistributionResult();
        }

        NetworkStruct sn = this.model.getStruct(false);
        
        if (R == null) {
            R = getAvgRespTHandles();
        }
        
        DistributionResult result = new DistributionResult(sn.nstations, sn.nclasses, "response_time");
        
        // Determine which nodes and classes need to be logged
        boolean[][] isNodeClassLogged = new boolean[this.model.getNumberOfNodes()][this.model.getNumberOfClasses()];
        for (int i = 0; i < this.model.getNumberOfStations(); i++) {
            for (int r = 0; r < this.model.getNumberOfClasses(); r++) {
                Station station = this.model.getStations().get(i);
                JobClass jobClass = this.model.getJobClasses().get(r);
                if (R.hasMetric(station, jobClass)) {
                    Metric metric = R.get(station, jobClass);
                    if (metric != null && !metric.isDisabled) {
                        int ni = this.model.getNodeIndex(station);
                        isNodeClassLogged[ni][r] = true;
                    }
                }
            }
        }
        
        // Determine which nodes need to be logged
        boolean[] isNodeLogged = new boolean[this.model.getNumberOfNodes()];
        for (int ni = 0; ni < this.model.getNumberOfNodes(); ni++) {
            for (int r = 0; r < this.model.getNumberOfClasses(); r++) {
                if (isNodeClassLogged[ni][r]) {
                    isNodeLogged[ni] = true;
                    break;
                }
            }
        }
        
        try {
            // Set up logging path
            String logPath = SysUtilsKt.lineTempName("jmt_logs");
            
            // Link and log the model (in a transient version of the model)
            // Use the original routing matrix from the network structure
            RoutingMatrix P = new RoutingMatrix(this.model, this.model.getClasses(), this.model.getNodes());
            this.model.linkAndLog(P, isNodeLogged, logPath);
            
            // Run simulation to generate log data
            this.getAvg();
            
            // Parse logs and compute CDFs
            Matrix[][][] RD = parseLogs(this.model, isNodeLogged, MetricType.RespT);
            
            // Convert from nodes to stations and compute empirical CDFs
            for (int i = 0; i < this.model.getNumberOfStations(); i++) {
                int ni = this.model.getNodeIndex(this.model.getStations().get(i));
                for (int r = 0; r < this.model.getNumberOfClasses(); r++) {
                    if (isNodeClassLogged[ni][r] && RD[ni][r] != null) {
                        // Compute empirical CDF using the response time data
                        double[] respTimes = RD[ni][r][0].toArray1D();
                        if (respTimes.length > 0) {
                            Matrix cdf = computeEmpiricalCDF(respTimes);
                            // Set CDF in result structure - format as [F, X] where F is CDF values, X is data values
                            result.setCdf(i, r, cdf);
                        }
                    }
                }
            }
            
        } catch (Exception e) {
            throw new RuntimeException("Error in getTranCdfRespT: " + e.getMessage(), e);
        }
        
        return result;
    }
    
    /**
     * Computes transient average station metrics over the specified time interval.
     * This method overrides NetworkSolver.getTranAvg() to provide JMT-specific transient analysis.
     */
    @Override
    public void getTranAvg() {
        if (Double.isInfinite(this.options.timespan[1])) {
            throw new RuntimeException("Transient analysis requires finite timespan. Please specify timespan option, e.g., SolverJMT(model, \"timespan\", new double[]{0, 10}).");
        }
        
        // Set up transient handles
        this.tranHandles = model.getTranHandles();
        NetworkStruct sn = this.model.getStruct(true);
        
        // Prepare logging configuration for transient metrics
        boolean[] isNodeLogged = new boolean[this.model.getNumberOfNodes()];
        boolean[][] isNodeClassLogged = new boolean[this.model.getNumberOfNodes()][this.model.getNumberOfClasses()];
        
        // Enable logging for all active handles
        for (int i = 0; i < this.model.getNumberOfStations(); i++) {
            for (int r = 0; r < this.model.getNumberOfClasses(); r++) {
                Station station = this.model.getStations().get(i);
                JobClass jobClass = this.model.getJobClasses().get(r);
                
                // Check if any metric is requested for this station-class pair
                boolean hasActiveMetric = false;
                if (this.tranHandles.Qt.hasMetric(station, jobClass) && !this.tranHandles.Qt.get(station, jobClass).isDisabled) {
                    hasActiveMetric = true;
                }
                if (this.tranHandles.Ut.hasMetric(station, jobClass) && !this.tranHandles.Ut.get(station, jobClass).isDisabled) {
                    hasActiveMetric = true;
                }
                if (this.tranHandles.Tt.hasMetric(station, jobClass) && !this.tranHandles.Tt.get(station, jobClass).isDisabled) {
                    hasActiveMetric = true;
                }
                
                if (hasActiveMetric) {
                    int ni = this.model.getNodeIndex(station);
                    isNodeLogged[ni] = true;
                    isNodeClassLogged[ni][r] = true;
                }
            }
        }
        
        try {
            // Set up logging path for transient analysis
            String logPath = SysUtilsKt.lineTempName("jmt_logs_tran");
            
            // Set initial method if not specified
            if (this.options.method.equals("default") || this.options.method.isEmpty()) {
                this.options.method = "jsim.standard";
            }
            
            // Link and log the model for transient data collection
            RoutingMatrix P = this.model.initRoutingMatrix();
            this.model.linkAndLog(P, isNodeLogged, logPath);
            
            // Run simulation to generate transient data
            this.runAnalyzer();
            
            // Parse transient data from simulation logs
            Matrix[][][] tranQData = parseTransientLogs(this.model, isNodeLogged, isNodeClassLogged, MetricType.QLen);
            Matrix[][][] tranUData = parseTransientLogs(this.model, isNodeLogged, isNodeClassLogged, MetricType.Util);
            Matrix[][][] tranTData = parseTransientLogs(this.model, isNodeLogged, isNodeClassLogged, MetricType.Tput);
            
            // Convert node-based data to station-based results
            Matrix[][] Qt = new Matrix[1][sn.nstations * sn.nclasses];
            Matrix[][] Ut = new Matrix[1][sn.nstations * sn.nclasses];
            Matrix[][] Tt = new Matrix[1][sn.nstations * sn.nclasses];
            
            for (int i = 0; i < sn.nstations; i++) {
                int ni = (int) sn.stationToNode.get(i);
                for (int r = 0; r < sn.nclasses; r++) {
                    int index = i * sn.nclasses + r;
                    
                    // Queue length transients
                    if (isNodeClassLogged[ni][r] && tranQData[ni][r] != null && tranQData[ni][r].length > 0) {
                        Qt[0][index] = tranQData[ni][r][0]; // Time series data
                    } else {
                        Qt[0][index] = new Matrix(0, 2); // Empty matrix [value, time] format
                    }
                    
                    // Utilization transients
                    if (isNodeClassLogged[ni][r] && tranUData[ni][r] != null && tranUData[ni][r].length > 0) {
                        Ut[0][index] = tranUData[ni][r][0];
                    } else {
                        Ut[0][index] = new Matrix(0, 2);
                    }
                    
                    // Throughput transients
                    if (isNodeClassLogged[ni][r] && tranTData[ni][r] != null && tranTData[ni][r].length > 0) {
                        Tt[0][index] = tranTData[ni][r][0];
                    } else {
                        Tt[0][index] = new Matrix(0, 2);
                    }
                }
            }
            
            // Store transient results using the parent method
            this.setTranAvgResults(Qt, Ut, new Matrix[0][0], Tt, new Matrix[0][0], new Matrix[0][0], this.result.runtime);
            
        } catch (Exception e) {
            throw new RuntimeException("Error in getTranAvg: " + e.getMessage(), e);
        }
    }
    
    /**
     * Returns transient queue length results from the last getTranAvg() call.
     * @return Matrix array with transient queue length data [stations x classes]
     */
    public Matrix[][] getTranQLen() {
        if (this.result == null || this.result.QNt == null || this.result.QNt.length == 0) {
            NetworkStruct sn = this.model.getStruct(false);
            return new Matrix[sn.nstations][sn.nclasses];
        }
        
        NetworkStruct sn = this.model.getStruct(false);
        Matrix[][] resultMatrix = new Matrix[sn.nstations][sn.nclasses];
        
        // The transient results are stored as Qt[time_index][station*class + class_index]
        // Convert back to [stations x classes] format
        for (int i = 0; i < sn.nstations; i++) {
            for (int r = 0; r < sn.nclasses; r++) {
                int index = i * sn.nclasses + r;
                if (index < this.result.QNt[0].length) {
                    resultMatrix[i][r] = this.result.QNt[0][index];
                }
            }
        }
        
        return resultMatrix;
    }
    
    /**
     * Returns transient utilization results from the last getTranAvg() call.
     * @return Matrix array with transient utilization data [stations x classes]
     */
    public Matrix[][] getTranUtil() {
        if (this.result == null || this.result.UNt == null || this.result.UNt.length == 0) {
            NetworkStruct sn = this.model.getStruct(false);
            return new Matrix[sn.nstations][sn.nclasses];
        }
        
        NetworkStruct sn = this.model.getStruct(false);
        Matrix[][] resultMatrix = new Matrix[sn.nstations][sn.nclasses];
        
        for (int i = 0; i < sn.nstations; i++) {
            for (int r = 0; r < sn.nclasses; r++) {
                int index = i * sn.nclasses + r;
                if (index < this.result.UNt[0].length) {
                    resultMatrix[i][r] = this.result.UNt[0][index];
                }
            }
        }
        
        return resultMatrix;
    }
    
    /**
     * Returns transient throughput results from the last getTranAvg() call.
     * @return Matrix array with transient throughput data [stations x classes]
     */
    public Matrix[][] getTranTput() {
        if (this.result == null || this.result.TNt == null || this.result.TNt.length == 0) {
            NetworkStruct sn = this.model.getStruct(false);
            return new Matrix[sn.nstations][sn.nclasses];
        }
        
        NetworkStruct sn = this.model.getStruct(false);
        Matrix[][] resultMatrix = new Matrix[sn.nstations][sn.nclasses];
        
        for (int i = 0; i < sn.nstations; i++) {
            for (int r = 0; r < sn.nclasses; r++) {
                int index = i * sn.nclasses + r;
                if (index < this.result.TNt[0].length) {
                    resultMatrix[i][r] = this.result.TNt[0][index];
                }
            }
        }
        
        return resultMatrix;
    }

    /**
     * Parses JMT log files to extract performance metrics data.
     * 
     * @param model The network model
     * @param isNodeLogged Array indicating which nodes are logged
     * @param metric The metric type to parse
     * @return 3D array of matrices containing parsed data [nodes][classes][data]
     */
    private Matrix[][][] parseLogs(Network model, boolean[] isNodeLogged, MetricType metric) {
        NetworkStruct sn = model.getStruct(false);
        int nclasses = sn.nclasses;
        Matrix[][][] logData = new Matrix[sn.nnodes][nclasses][];

        String logPath = model.getLogPath();

        for (int ind = 0; ind < sn.nnodes; ind++) {
            boolean isStateful = sn.isstateful.get(ind) == 1.0;
            boolean isLogged = ind < isNodeLogged.length && isNodeLogged[ind];

            if (isStateful && isLogged) {
                String nodeName = model.getNodeNames().get(ind);
                String logFileArv = logPath + "/" + nodeName + "-Arv.csv";
                String logFileDep = logPath + "/" + nodeName + "-Dep.csv";

                boolean arvExists = new java.io.File(logFileArv).exists();
                boolean depExists = new java.io.File(logFileDep).exists();

                try {
                    if (arvExists && depExists) {
                        
                        // Parse arrival data
                        java.util.List<String[]> arvData = parseCSVLog(logFileArv);
                        java.util.List<String[]> depData = parseCSVLog(logFileDep);
                        
                        if (metric == MetricType.RespT) {
                            // Parse response time data
                            Matrix[][] nodeRespTData = parseTranRespT(arvData, depData, model);
                            for (int r = 0; r < Math.min(nclasses, nodeRespTData.length); r++) {
                                if (nodeRespTData[r] != null) {
                                    logData[ind][r] = nodeRespTData[r];
                                }
                            }
                        } else if (metric == MetricType.QLen) {
                            // Parse queue length data from arrival/departure logs
                            Matrix[][] nodeQLenData = parseTranQLen(arvData, depData, model);
                            for (int r = 0; r < Math.min(nclasses, nodeQLenData.length); r++) {
                                if (nodeQLenData[r] != null) {
                                    logData[ind][r] = nodeQLenData[r];
                                }
                            }
                        }
                    }
                } catch (Exception e) {
                    // Log the error but continue with other nodes
                    line_warning("SolverJMT.getLog", "Error parsing logs for node %d: %s", ind, e.getMessage());
                }
            }
        }

        return logData;
    }

    /**
     * Parses CSV log files with semicolon delimiter.
     * 
     * @param filePath Path to the CSV file
     * @return List of parsed rows
     */
    private java.util.List<String[]> parseCSVLog(String filePath) throws java.io.IOException {
        java.util.List<String[]> data = new java.util.ArrayList<>();
        try (java.io.BufferedReader br = new java.io.BufferedReader(new java.io.FileReader(filePath))) {
            String line;
            boolean firstLine = true;
            while ((line = br.readLine()) != null) {
                if (firstLine) {
                    firstLine = false;
                    continue; // Skip header
                }
                String[] values = line.split(";");
                data.add(values);
            }
        }
        return data;
    }
    
    /**
     * Parses transient response time data from arrival and departure logs.
     * 
     * @param arvData Arrival log data
     * @param depData Departure log data
     * @param model The network model
     * @return Array of matrices containing response time data per class
     */
    private Matrix[][] parseTranRespT(java.util.List<String[]> arvData, java.util.List<String[]> depData, Network model) {
        // Following MATLAB approach: combine all events and sort by jobID then timestamp
        java.util.List<JobEvent> allEvents = new java.util.ArrayList<>();
        
        // Parse arrival data (type = +1)
        for (String[] row : arvData) {
            if (row.length >= 4) {
                double timestamp = Double.parseDouble(row[1]);
                int jobId = (int) Double.parseDouble(row[2]);
                String className = row[3];
                int classId = getClassIndex(model, className);
                allEvents.add(new JobEvent(timestamp, jobId, classId, 1)); // +1 for arrival
            }
        }
        
        // Parse departure data (type = -1)
        for (String[] row : depData) {
            if (row.length >= 4) {
                double timestamp = Double.parseDouble(row[1]);
                int jobId = (int) Double.parseDouble(row[2]);
                String className = row[3];
                int classId = getClassIndex(model, className);
                allEvents.add(new JobEvent(timestamp, jobId, classId, -1)); // -1 for departure
            }
        }
        
        // Sort by jobId, then by timestamp
        allEvents.sort((a, b) -> {
            if (a.jobId != b.jobId) return Integer.compare(a.jobId, b.jobId);
            return Double.compare(a.timestamp, b.timestamp);
        });
        
        // Process events per job ID to compute response times
        int nclasses = model.getNumberOfClasses();
        java.util.List<java.util.List<Double>> classRespTimes = new java.util.ArrayList<>();
        for (int i = 0; i < nclasses; i++) {
            classRespTimes.add(new java.util.ArrayList<>());
        }
        
        int currentJobId = -1;
        java.util.List<JobEvent> currentJobEvents = new java.util.ArrayList<>();
        
        for (JobEvent event : allEvents) {
            if (event.jobId != currentJobId) {
                // Process previous job's events
                if (!currentJobEvents.isEmpty()) {
                    processJobEvents(currentJobEvents, classRespTimes);
                }
                currentJobId = event.jobId;
                currentJobEvents = new java.util.ArrayList<>();
            }
            currentJobEvents.add(event);
        }
        // Process last job's events
        if (!currentJobEvents.isEmpty()) {
            processJobEvents(currentJobEvents, classRespTimes);
        }
        
        // Convert to result format
        Matrix[][] classRespTData = new Matrix[nclasses][];
        for (int r = 0; r < nclasses; r++) {
            if (!classRespTimes.get(r).isEmpty()) {
                Matrix respTMatrix = new Matrix(classRespTimes.get(r).size(), 1);
                for (int i = 0; i < classRespTimes.get(r).size(); i++) {
                    respTMatrix.set(i, 0, classRespTimes.get(r).get(i));
                }
                classRespTData[r] = new Matrix[]{respTMatrix};
            }
        }
        
        return classRespTData;
    }
    
    // Helper method to process events for a single job
    private void processJobEvents(java.util.List<JobEvent> events, java.util.List<java.util.List<Double>> classRespTimes) {
        // Find first arrival and last departure
        int firstArrival = -1;
        int lastDeparture = -1;
        
        for (int i = 0; i < events.size(); i++) {
            if (events.get(i).type > 0 && firstArrival < 0) {
                firstArrival = i;
            }
            if (events.get(i).type < 0) {
                lastDeparture = i;
            }
        }
        
        if (firstArrival >= 0 && lastDeparture >= firstArrival) {
            // Process events between first arrival and last departure
            java.util.List<JobEvent> relevantEvents = events.subList(firstArrival, lastDeparture + 1);
            
            // Pair arrivals with departures
            for (int i = 0; i < relevantEvents.size() - 1; i += 2) {
                if (i + 1 < relevantEvents.size() && 
                    relevantEvents.get(i).type > 0 && 
                    relevantEvents.get(i + 1).type < 0) {
                    // We have an arrival followed by a departure
                    double respTime = relevantEvents.get(i + 1).timestamp - relevantEvents.get(i).timestamp;
                    int classId = relevantEvents.get(i).classId;
                    if (respTime >= 0 && classId >= 0 && classId < classRespTimes.size()) {
                        classRespTimes.get(classId).add(respTime);
                    }
                }
            }
        }
    }
    
    // Helper class for job events
    private static class JobEvent {
        double timestamp;
        int jobId;
        int classId;
        int type; // +1 for arrival, -1 for departure
        
        JobEvent(double timestamp, int jobId, int classId, int type) {
            this.timestamp = timestamp;
            this.jobId = jobId;
            this.classId = classId;
            this.type = type;
        }
    }
    
    /**
     * Helper method to get class index from class name.
     * 
     * @param model The network model
     * @param className The class name
     * @return The class index
     */
    private int getClassIndex(Network model, String className) {
        java.util.List<String> classNames = model.getClassNames();
        for (int i = 0; i < classNames.size(); i++) {
            if (classNames.get(i).equals(className)) {
                return i;
            }
        }
        return 0; // Default to first class if not found
    }

    /**
     * Parses transient queue length data from arrival and departure logs.
     * Computes queue length over time by tracking arrivals (+1) and departures (-1).
     *
     * @param arvData Arrival log data
     * @param depData Departure log data
     * @param model The network model
     * @return Array of matrices containing queue length data per class [time, qlen]
     */
    private Matrix[][] parseTranQLen(java.util.List<String[]> arvData, java.util.List<String[]> depData, Network model) {
        int nclasses = model.getNumberOfClasses();

        // Create lists to hold (timestamp, qlen_change) events per class
        java.util.List<java.util.List<double[]>> classEvents = new java.util.ArrayList<>();
        for (int r = 0; r < nclasses; r++) {
            classEvents.add(new java.util.ArrayList<>());
        }

        // Parse arrival data (queue length change = +1)
        for (String[] row : arvData) {
            if (row.length >= 4) {
                double timestamp = Double.parseDouble(row[1]);
                String className = row[3];
                int classId = getClassIndex(model, className);
                if (classId >= 0 && classId < nclasses) {
                    classEvents.get(classId).add(new double[]{timestamp, 1.0});
                }
            }
        }

        // Parse departure data (queue length change = -1)
        for (String[] row : depData) {
            if (row.length >= 4) {
                double timestamp = Double.parseDouble(row[1]);
                String className = row[3];
                int classId = getClassIndex(model, className);
                if (classId >= 0 && classId < nclasses) {
                    classEvents.get(classId).add(new double[]{timestamp, -1.0});
                }
            }
        }

        // Build result matrices for each class
        Matrix[][] classQLenData = new Matrix[nclasses][];

        for (int r = 0; r < nclasses; r++) {
            java.util.List<double[]> events = classEvents.get(r);

            if (events.isEmpty()) {
                continue;
            }

            // Sort events by timestamp
            events.sort((a, b) -> Double.compare(a[0], b[0]));

            // Compute cumulative queue length over time
            java.util.List<double[]> timeQLenPairs = new java.util.ArrayList<>();
            double currentQLen = 0.0;

            for (double[] event : events) {
                double timestamp = event[0];
                double change = event[1];
                currentQLen += change;
                if (currentQLen < 0) currentQLen = 0; // Safety check
                timeQLenPairs.add(new double[]{timestamp, currentQLen});
            }

            if (!timeQLenPairs.isEmpty()) {
                // Create time matrix and qlen matrix
                Matrix timeMatrix = new Matrix(timeQLenPairs.size(), 1);
                Matrix qlenMatrix = new Matrix(timeQLenPairs.size(), 1);

                for (int i = 0; i < timeQLenPairs.size(); i++) {
                    timeMatrix.set(i, 0, timeQLenPairs.get(i)[0]);
                    qlenMatrix.set(i, 0, timeQLenPairs.get(i)[1]);
                }

                classQLenData[r] = new Matrix[]{timeMatrix, qlenMatrix};
            }
        }

        return classQLenData;
    }

    /**
     * Parses transient metrics from JMT simulation logs for the specified metric type.
     * This method extracts time-series data for transient analysis.
     * 
     * @param model The network model
     * @param isNodeLogged Array indicating which nodes are logged
     * @param isNodeClassLogged 2D array indicating which node-class pairs are logged  
     * @param metricType The type of metric to parse (QLen, Util, Tput)
     * @return 3D array of matrices [node][class][time_series_data]
     */
    private Matrix[][][] parseTransientLogs(Network model, boolean[] isNodeLogged, boolean[][] isNodeClassLogged, MetricType metricType) {
        NetworkStruct sn = model.getStruct(false);
        int nnodes = model.getNumberOfNodes();  // Use model.getNumberOfNodes() to match array size
        int nclasses = sn.nclasses;
        Matrix[][][] tranData = new Matrix[nnodes][nclasses][];
        
        // Iterate over all nodes that have logging enabled
        for (int ni = 0; ni < nnodes && ni < isNodeLogged.length; ni++) {
            if (!isNodeLogged[ni]) continue;
            
            String nodeName = model.getNodeNames().get(ni);
            String arvLogFile = model.getLogPath() + "/" + nodeName + "-Arv.csv";
            String depLogFile = model.getLogPath() + "/" + nodeName + "-Dep.csv";
            
            try {
                java.io.File arvFile = new java.io.File(arvLogFile);
                java.io.File depFile = new java.io.File(depLogFile);
                
                if (arvFile.exists() && depFile.exists()) {
                    // Parse both arrival and departure logs
                    java.util.List<String[]> arvData = parseCSVLog(arvLogFile);
                    java.util.List<String[]> depData = parseCSVLog(depLogFile);
                    
                    // Process transient data based on metric type
                    switch (metricType) {
                        case QLen:
                            tranData[ni] = parseTransientQueueLength(arvData, depData, model, isNodeClassLogged[ni]);
                            break;
                        case Util:
                            tranData[ni] = parseTransientUtilization(arvData, depData, model, isNodeClassLogged[ni]);
                            break;
                        case Tput:
                            tranData[ni] = parseTransientThroughput(depData, model, isNodeClassLogged[ni]);
                            break;
                        default:
                            // Unsupported metric type, leave empty
                            break;
                    }
                }
            } catch (Exception e) {
                // Log error but continue with other nodes
                line_warning("SolverJMT", "Could not parse transient logs for node %d, metric %s: %s", ni, metricType, e.getMessage());
            }
        }
        
        return tranData;
    }
    
    /**
     * Parses transient queue length data from arrival and departure logs.
     */
    private Matrix[][] parseTransientQueueLength(java.util.List<String[]> arvData, java.util.List<String[]> depData, Network model, boolean[] isClassLogged) {
        int nclasses = model.getNumberOfClasses();
        Matrix[][] qlenData = new Matrix[nclasses][];
        
        try {
            // Create time-ordered event list for all classes
            java.util.List<TransientEvent> allEvents = new java.util.ArrayList<>();
            
            // Add arrivals (+1 for queue length)
            for (String[] row : arvData) {
                if (row.length >= 4) {
                    double timestamp = Double.parseDouble(row[1]);
                    String className = row[3];
                    int classId = getClassIndex(model, className);
                    if (classId < isClassLogged.length && isClassLogged[classId]) {
                        allEvents.add(new TransientEvent(timestamp, classId, 1));
                    }
                }
            }
            
            // Add departures (-1 for queue length)
            for (String[] row : depData) {
                if (row.length >= 4) {
                    double timestamp = Double.parseDouble(row[1]);
                    String className = row[3];
                    int classId = getClassIndex(model, className);
                    if (classId < isClassLogged.length && isClassLogged[classId]) {
                        allEvents.add(new TransientEvent(timestamp, classId, -1));
                    }
                }
            }
            
            // Sort events by timestamp
            allEvents.sort(java.util.Comparator.comparingDouble(e -> e.timestamp));
            
            // Apply timespan filter
            double startTime = this.options.timespan[0];
            double endTime = this.options.timespan[1];
            
            // Track queue lengths per class
            int[] queueLengths = new int[nclasses];
            java.util.List<java.util.List<Double>> timeSeries = new java.util.ArrayList<>();
            java.util.List<java.util.List<Double>> valueSeries = new java.util.ArrayList<>();
            
            for (int r = 0; r < nclasses; r++) {
                timeSeries.add(new java.util.ArrayList<>());
                valueSeries.add(new java.util.ArrayList<>());
                if (isClassLogged[r]) {
                    // Initialize with zero at start time
                    timeSeries.get(r).add(startTime);
                    valueSeries.get(r).add(0.0);
                }
            }
            
            // Process events to build queue length time series
            for (TransientEvent event : allEvents) {
                if (event.timestamp >= startTime && event.timestamp <= endTime) {
                    int r = event.classId;
                    if (r >= 0 && r < nclasses && isClassLogged[r]) {
                        queueLengths[r] = Math.max(0, queueLengths[r] + event.delta);
                        timeSeries.get(r).add(event.timestamp);
                        valueSeries.get(r).add((double) queueLengths[r]);
                    }
                }
            }
            
            // Convert to matrix format [value, time] as expected by LINE
            for (int r = 0; r < nclasses; r++) {
                if (isClassLogged[r] && !timeSeries.get(r).isEmpty()) {
                    int numPoints = timeSeries.get(r).size();
                    Matrix tranMatrix = new Matrix(numPoints, 2);
                    for (int i = 0; i < numPoints; i++) {
                        tranMatrix.set(i, 0, valueSeries.get(r).get(i)); // value
                        tranMatrix.set(i, 1, timeSeries.get(r).get(i));   // time
                    }
                    qlenData[r] = new Matrix[]{tranMatrix};
                }
            }

        } catch (Exception e) {
            line_warning("SolverJMT.parseTransientQueueLength", "Error parsing queue length transients: %s", e.getMessage());
        }

        return qlenData;
    }

    /**
     * Parses transient utilization data (simplified - assumes server utilization proportional to queue occupancy).
     */
    private Matrix[][] parseTransientUtilization(java.util.List<String[]> arvData, java.util.List<String[]> depData, Network model, boolean[] isClassLogged) {
        // For simplicity, approximate utilization as min(queueLength, 1) for single server stations
        // This is a simplification - proper utilization would require service time tracking
        Matrix[][] qlenData = parseTransientQueueLength(arvData, depData, model, isClassLogged);
        
        if (qlenData != null) {
            for (int r = 0; r < qlenData.length; r++) {
                if (qlenData[r] != null && qlenData[r].length > 0) {
                    Matrix qMatrix = qlenData[r][0];
                    Matrix uMatrix = new Matrix(qMatrix.getNumRows(), qMatrix.getNumCols());
                    for (int i = 0; i < qMatrix.getNumRows(); i++) {
                        double qlen = qMatrix.get(i, 0);
                        double util = Math.min(qlen, 1.0); // Simplified utilization
                        uMatrix.set(i, 0, util); // value
                        uMatrix.set(i, 1, qMatrix.get(i, 1)); // time
                    }
                    qlenData[r][0] = uMatrix;
                }
            }
        }
        
        return qlenData;
    }
    
    /**
     * Parses transient throughput data from departure events.
     */
    private Matrix[][] parseTransientThroughput(java.util.List<String[]> depData, Network model, boolean[] isClassLogged) {
        int nclasses = model.getNumberOfClasses();
        Matrix[][] tputData = new Matrix[nclasses][];
        
        try {
            double startTime = this.options.timespan[0];
            double endTime = this.options.timespan[1];
            double windowSize = Math.max(1.0, (endTime - startTime) / 100.0); // Adaptive window size
            
            // Count departures per class in time windows
            for (int r = 0; r < nclasses; r++) {
                if (!isClassLogged[r]) continue;
                
                java.util.List<Double> timePoints = new java.util.ArrayList<>();
                java.util.List<Double> tputValues = new java.util.ArrayList<>();
                
                // Create time windows
                for (double t = startTime; t < endTime; t += windowSize) {
                    double windowEnd = Math.min(t + windowSize, endTime);
                    int count = 0;
                    
                    // Count departures in this window for class r
                    for (String[] row : depData) {
                        if (row.length >= 4) {
                            double timestamp = Double.parseDouble(row[1]);
                            String className = row[3];
                            int classId = getClassIndex(model, className);
                            
                            if (classId == r && timestamp >= t && timestamp < windowEnd) {
                                count++;
                            }
                        }
                    }
                    
                    double throughput = count / windowSize; // Jobs per time unit
                    timePoints.add(t + windowSize / 2); // Middle of window
                    tputValues.add(throughput);
                }
                
                // Convert to matrix format
                if (!timePoints.isEmpty()) {
                    Matrix tranMatrix = new Matrix(timePoints.size(), 2);
                    for (int i = 0; i < timePoints.size(); i++) {
                        tranMatrix.set(i, 0, tputValues.get(i)); // value  
                        tranMatrix.set(i, 1, timePoints.get(i));  // time
                    }
                    tputData[r] = new Matrix[]{tranMatrix};
                }
            }

        } catch (Exception e) {
            line_warning("SolverJMT.parseTransientThroughput", "Error parsing throughput transients: %s", e.getMessage());
        }

        return tputData;
    }

    /**
     * Helper class for transient event processing.
     */
    private static class TransientEvent {
        double timestamp;
        int classId;
        int delta; // +1 for arrival, -1 for departure
        
        TransientEvent(double timestamp, int classId, int delta) {
            this.timestamp = timestamp;
            this.classId = classId;
            this.delta = delta;
        }
    }
    
    /**
     * Computes empirical cumulative distribution function from data.
     * 
     * @param data Input data array
     * @return Matrix with CDF values [F, X] where F is cumulative probability, X is sorted data
     */
    private Matrix computeEmpiricalCDF(double[] data) {
        if (data.length == 0) {
            return new Matrix(2, 0);
        }
        
        // Sort the data
        java.util.Arrays.sort(data);
        
        // Remove duplicates and compute CDF
        java.util.List<Double> uniqueValues = new java.util.ArrayList<>();
        java.util.List<Double> cdfValues = new java.util.ArrayList<>();
        
        uniqueValues.add(data[0]);
        double currentValue = data[0];
        int currentCount = 1;
        
        for (int i = 1; i < data.length; i++) {
            if (data[i] != currentValue) {
                cdfValues.add((double) currentCount / data.length);
                uniqueValues.add(data[i]);
                currentValue = data[i];
                currentCount = 1;
            } else {
                currentCount++;
            }
        }
        // Add the last point
        cdfValues.add(1.0);
        
        // Convert to matrix format [F; X]
        Matrix result = new Matrix(2, uniqueValues.size());
        for (int i = 0; i < uniqueValues.size(); i++) {
            result.set(0, i, cdfValues.get(i));
            result.set(1, i, uniqueValues.get(i));
        }
        
        return result;
    }


    protected boolean hasAvgResults() {
        return hasResults();
    }

    public void jsimgView() {
        jsimgView(jmtGetPath(), this.options);
    }

    public void jsimgView(SolverOptions options) {
        jsimgView(jmtGetPath(), options);
    }

    public void jsimgView(String jmtPath, SolverOptions options) {
        if (this.enableChecks && !supports(this.model)) {
            line_error(mfilename(new Object() {
            }), "This model contains features not supported by the solver.");
        }

        if (options == null) {
            options = defaultOptions();
        }

        if (options.samples == 0) {
            options.samples = 10000;
        } else if (options.samples < 5000) {
            // line_warning
            line_error(mfilename(new Object[]{}), "JMT requires at least 5000 samples for each metric. Setting the samples to 5000.");
            options.samples = 5000;
        }

        // set seed and maxSamples
        this.seed = options.seed;
        RandomManager.setMasterSeed(options.seed);
        this.maxSamples = options.samples;

        NetworkStruct sn = getStruct();
        try {
            this.writeJSIM(sn);
        } catch (ParserConfigurationException e) {
            line_error(mfilename(new Object() {
            }), "XML parsing error.");
        }

        String fileName = this.getFilePath() + File.separator + this.getFileName() + ".jsim";
        if (options.verbose!=VerboseLevel.SILENT){
        java.lang.System.out.println("JMT Model: " + fileName);
        java.lang.System.out.flush();
    }

        viewModel(jmtPath, fileName, ViewMode.JSIMG, options.verbose);
    }

    public void jsimgView(String jmtPath) {
        jsimgView(jmtPath, SolverJMT.defaultOptions());
    }

    public void jsimwView(String jmtPath) {
        jsimwView(jmtPath, SolverJMT.defaultOptions());
    }

    public void jsimwView(String jmtPath, SolverOptions options) {
        if (this.enableChecks && !supports(this.model)) {
            line_error(mfilename(new Object() {
            }), "This model contains features not supported by the solver.");
        }

        if (options.samples < 5000) {
            line_error(mfilename(new Object[]{}), "JMT requires at least 5000 samples for each metric. Setting the samples to 5000.");
            options.samples = 5000;
        }

        this.seed = options.seed;
        RandomManager.setMasterSeed(options.seed);
        this.maxSamples = options.samples;

        NetworkStruct sn = getStruct();
        try {
            this.writeJSIM(sn);
        } catch (ParserConfigurationException e) {
            line_error(mfilename(new Object() {
            }), "XML parsing error.");
        }

        String fileName = this.getFilePath() + File.separator + this.getFileName() + ".jsim";
        //java.lang.System.out.println("JMT Model: " + fileName);

        viewModel(jmtPath, fileName, ViewMode.JSIMW, options.verbose);
    }

    public void jsimwView() throws ParserConfigurationException {
        jsimwView(jmtGetPath(), this.options);
    }

    public List<String> listValidMethods() {
        return listValidMethods(null);
    }

    public List<String> listValidMethods(Network model) {

        return Arrays.asList(
                "default",
                "jsim",
                "jmva",
                "jmva.amva",
                "jmva.mva",
                "jmva.recal",
                "jmva.comom",
                "jmva.chow",
                "jmva.bs",
                "jmva.aql",
                "jmva.lin",
                "jmva.dmlin"); // "jmva.ls"
    }

    public double probSysStateAggr() {
        if (GlobalConstants.DummyMode) {
            return NaN;
        }
        
        try {
            NetworkStruct sn = this.getStruct();
            
            // Get system state samples
            SampleResult tranSysStateAggr = this.sampleSysAggr();
            if (tranSysStateAggr == null || tranSysStateAggr.state == null) {
                return 0.0;
            }
            
            // Extract time and state data
            Matrix timeData = tranSysStateAggr.t;
            Object stateDataObj = tranSysStateAggr.state;
            
            if (timeData == null || stateDataObj == null) {
                return 0.0;
            }
            
            // Cast state data to Matrix (assuming single node sampling)
            Matrix stateData;
            if (stateDataObj instanceof Matrix) {
                stateData = (Matrix) stateDataObj;
            } else {
                return 0.0; // Cannot handle list format for this operation
            }
            
            // Convert time differences (like MATLAB: TSS(:,1)=[diff(TSS(:,1));0])
            Matrix timeDiffs = new Matrix(timeData.getNumRows(), 1);
            for (int i = 0; i < timeData.getNumRows() - 1; i++) {
                timeDiffs.set(i, 0, timeData.get(i + 1, 0) - timeData.get(i, 0));
            }
            timeDiffs.set(timeDiffs.getNumRows() - 1, 0, 0.0); // Last element is 0
            
            // Get current network state in marginal form
            Matrix currentStateMarginal = getCurrentStateMarginal(sn);
            
            // Find rows in stateData that match the current state
            double totalTimeInState = 0.0;
            double totalTime = timeDiffs.elementSum();
            
            for (int i = 0; i < stateData.getNumRows(); i++) {
                boolean stateMatches = true;
                for (int j = 0; j < currentStateMarginal.length() && j < stateData.getNumCols(); j++) {
                    if (Math.abs(stateData.get(i, j) - currentStateMarginal.get(j)) > 1e-10) {
                        stateMatches = false;
                        break;
                    }
                }
                
                if (stateMatches) {
                    totalTimeInState += timeDiffs.get(i, 0);
                }
            }
            
            if (totalTime > 0) {
                return totalTimeInState / totalTime;
            } else {
                return 0.0;
            }
            
        } catch (Exception e) {
            line_warning("SolverJMT", "The state was not seen during the simulation.");
            return 0.0;
        }
    }
    
    private Matrix getCurrentStateMarginal(NetworkStruct sn) {
        // Convert current network state to marginal representation
        // This is a simplified implementation - may need refinement
        Matrix marginal = new Matrix(sn.nstateful * sn.nclasses, 1);
        int idx = 0;
        
        for (int isf = 0; isf < sn.nstateful; isf++) {
            int ind = (int) sn.statefulToNode.get(isf);
            // Get marginal state for this stateful node
            if (sn.state != null && sn.state.size() > isf && sn.state.get(isf) != null) {
                Matrix nodeState = sn.state.get(isf);
                State.StateMarginalStatistics stats = ToMarginal.toMarginal(sn, ind, nodeState, null, null, null, null, null);
                if (stats != null && stats.nir != null) {
                    for (int r = 0; r < Math.min(sn.nclasses, stats.nir.getNumCols()); r++) {
                        if (idx < marginal.length()) {
                            marginal.set(idx++, 0, stats.nir.get(0, r));
                        }
                    }
                }
            } else {
                // Fill with zeros if state not available
                for (int r = 0; r < sn.nclasses && idx < marginal.length(); r++) {
                    marginal.set(idx++, 0, 0.0);
                }
            }
        }
        
        return marginal;
    }

    @Override
    public void runAnalyzer() throws ParserConfigurationException {
        long startTime = java.lang.System.nanoTime();
//        if (this.enableChecks && !SolverJMT.supports(this.model)) {
//            throw new RuntimeException("This model contains features not supported by the solver.");
//        }

        if (this.model == null)
            throw new RuntimeException("Model is not provided");
        if (this.options == null)
            this.options = new SolverOptions(SolverType.JMT);
        if (options.verbose == null) {
            options.verbose = VerboseLevel.values()[0];
        }
        if (options.samples == 0) {
            options.samples = 10000;
        } else if (options.samples < 5000) {
            //if (!options.method.equalsIgnoreCase("jmva.ls")) {
            line_warning(mfilename(new Object() {
            }), String.format("JMT requires at least 5000 samples for each metric, the current value is %d. Starting the simulation with 5000 samples.%n", options.samples));
            //}
            options.samples = 5000;
        }
        if (options.seed == 0) {
            options.seed = RandomManager.generateRandomSeed();
        }
        this.seed = options.seed;
        RandomManager.setMasterSeed(options.seed);
        if (options.timespan == null) {
            options.timespan = new double[]{0.0, Inf};
        } else {
            this.maxSimulatedTime = options.timespan[1];
        }

        if (!this.model.hasInitState()) {
            this.model.initDefault();
        }
        this.maxSamples = options.samples;

        // Initialize cache hit/miss probabilities with default values if not already set.
        // This prevents NaN values in the routing matrix which JMT cannot handle.
        // Without this, JMT would require running MVA or another solver first to compute
        // the actual cache hit/miss probabilities.
        // Default hit probability is estimated as cache_size / number_of_items.
        for (Node node : this.model.getNodes()) {
            if (node instanceof Cache) {
                Cache cacheNode = (Cache) node;
                if (cacheNode.getCacheServer().actualHitProb.isEmpty()) {
                    Matrix hitClass = cacheNode.getHitClass();
                    Matrix hitProb = new Matrix(1, hitClass.length());
                    Matrix missProb = new Matrix(1, hitClass.length());
                    // Compute cache size as sum of all level capacities
                    double cacheSize = cacheNode.getItemLevelCap().elementSum();
                    int numItems = cacheNode.getNumberOfItems();
                    // Initial hit probability estimate: cache_size / number_of_items
                    double defaultHitProb = Math.min(cacheSize / numItems, 1.0);
                    for (int k = 0; k < hitClass.length(); k++) {
                        if (hitClass.get(k) >= 0) {
                            hitProb.set(k, defaultHitProb);
                            missProb.set(k, 1.0 - defaultHitProb);
                        }
                    }
                    cacheNode.setResultHitProb(hitProb);
                    cacheNode.setResultMissProb(missProb);
                    // Force struct refresh since cache probabilities affect routing
                    this.model.refreshStruct(true);
                    this.sn = null; // Clear cached struct to force refresh
                }
            }
        }

        NetworkStruct sn = this.getStruct();
        String fname = "";
        String cmd = "";
        String cmdOutput = "";
        long runTime = 0;
        SolverResult solverResult = null;


        switch (options.method) {
            case "jsim":
            case "default":
                fname = this.writeJSIM(sn);
                cmd = String.format(
                        "java -cp %s jmt.commandline.Jmt sim %s -seed %s --illegal-access=permit",
                        this.jmtPath, fname, options.seed
                );

                if (options.verbose != VerboseLevel.SILENT) {
                    java.lang.System.out.println("JMT Model: " + fname);
                }
                if (options.verbose == VerboseLevel.DEBUG) {
                    java.lang.System.out.println("JMT Command: " + cmd);
                }
                if (options.verbose == VerboseLevel.DEBUG) {
                    cmdOutput = SysUtilsKt.system(cmd, simulationTimeoutSeconds);
                } else {
                    // Suppress system command output unless in DEBUG mode
                    java.io.ByteArrayOutputStream devNull = new java.io.ByteArrayOutputStream();
                    java.io.PrintStream nullStream = new java.io.PrintStream(devNull);
                    java.io.PrintStream originalOut = System.out;
                    java.io.PrintStream originalErr = System.err;

                    try {
                        System.setOut(nullStream);
                        System.setErr(nullStream);
                        cmdOutput = SysUtilsKt.system(cmd, simulationTimeoutSeconds);
                    } finally {
                        System.setOut(originalOut);
                        System.setErr(originalErr);
                        nullStream.close();
                    }
                }

                this.lastCommandOutput = cmdOutput; // Store for error reporting
                // Check for timeout and warn user
                if (cmdOutput.startsWith("TIMEOUT:")) {
                    line_warning(mfilename(new Object[]{}), "JMT simulation timed out. Consider reducing samples or using a different solver.");
                }
                if (options.verbose != VerboseLevel.SILENT && !cmdOutput.isEmpty()) {
                    java.lang.System.out.println("JMT Command output: " + cmdOutput);
                }
                runTime = java.lang.System.nanoTime() - startTime;

                solverResult = getResults();
                if (!options.keep) {
//                    try {
//                        //JMT.removeDirectory(Paths.get(this.getFilePath()));
//                    } catch (IOException ioe) {
//                        ioe.printStackTrace();
//                    }
                }
                solverResult.runtime = runTime / 1000000000.0;
                this.result = solverResult;
                
                if (this.options.verbose != VerboseLevel.SILENT) {
                    System.out.printf(
                            "%s analysis [method: %s, lang: %s, env: %s] completed in %fs.\n",
                            this.name.replaceFirst("^Solver", ""),
                            this.result.method,
                            "java",
                            System.getProperty("java.version"),
                            this.result.runtime
                    );
                    System.out.flush();
                }
                break;
            case "closing":
                // Closing simulation for transient analysis
                NetworkStruct snClosing = this.getStruct();
                double initSeed = options.seed;
                double[] initTimeSpan = options.timespan.clone();
                
                // Set timespan for closing simulation
                options.timespan[0] = options.timespan[1];
                
                if (Double.isFinite(options.timespan[1])) {
                    // Perform multiple simulation runs with different seeds
                    for (int it = 0; it < options.iter_max; it++) {
                        options.seed = (int)(initSeed + it);
                        
                        // Run transient analysis (simplified implementation)
                        fname = this.writeJSIM(snClosing);
                        cmd = String.format(
                                "java -cp %s jmt.commandline.Jmt sim %s -seed %s --illegal-access=permit",
                                this.jmtPath, fname, options.seed
                        );
                        
                        if (options.verbose != VerboseLevel.SILENT) {
                            java.lang.System.out.println("JMT Closing Simulation " + (it + 1) + "/" + options.iter_max);
                        }

                        if (options.verbose == VerboseLevel.DEBUG) {
                            SysUtilsKt.system(cmd, simulationTimeoutSeconds);
                        } else {
                            // Suppress system command output unless in DEBUG mode
                            java.io.ByteArrayOutputStream devNull = new java.io.ByteArrayOutputStream();
                            java.io.PrintStream nullStream = new java.io.PrintStream(devNull);
                            java.io.PrintStream originalOut = System.out;
                            java.io.PrintStream originalErr = System.err;

                            try {
                                System.setOut(nullStream);
                                System.setErr(nullStream);
                                SysUtilsKt.system(cmd, simulationTimeoutSeconds);
                            } finally {
                                System.setOut(originalOut);
                                System.setErr(originalErr);
                                nullStream.close();
                            }
                        }
                    }
                    
                    // Restore original options
                    options.seed = (int)initSeed;
                    options.timespan = initTimeSpan;
                    
                    // Parse results from the final simulation
                    SolverResult result = this.getResults();
                    this.setAvgResults(result.QN, result.UN, result.RN, result.TN, result.AN, result.WN, result.CN, result.XN, 
                                     result.runtime, result.method, result.iter);
                } else {
                    throw new RuntimeException("Closing method requires finite timespan[1]");
                }

//                SolverOptions options = this.options;
//                int initSeed = this.options.seed;
//                if (options.timespan.length > 1){
//                    double[] initTimeSpan = this.options.timespan;
//                    this.options.timespan[0] = this.options.timespan[1];
//                    if (Double.isFinite(this.options.timespan[1])){
//                        double [] tu;
//                        for (int it = 0; it < options.iter_max; it++) {
//                            this.options.seed = initSeed + it - 1;
//                            TranSysStateAggr{it} = sampleSysAggr(self);
//                            if isempty(tu)
//                            tu = TranSysStateAggr{it}.t;
//                        else
//                            % we need to limit the time series at the minimum
//                            % as otherwise the predictor of the state cannot
//                            % take into account constraints that exist on the
//                            % state space
//                            tumax = min(max(tu),max(TranSysStateAggr{it}.t));
//                            tu = union(tu, TranSysStateAggr{it}.t);
//                            tu = tu(tu<=tumax);
//                            end
//                        }
//                        Matrix QNt = new Matrix(sn.nstations, sn.nclasses, tu.length, 2);
//                        Matrix UNt = new Matrix(sn.nstations, sn.nclasses, tu.length, 2);
//                        Matrix TNt = new Matrix(sn.nstations, sn.nclasses, tu.length, 2);
//
//                    }
//                }

            case "jmva":
            case "jmva.amva":
            case "jmva.mva":
            case "jmva.recal":
            case "jmva.comom":
            case "jmva.chow":
            case "jmva.bs":
            case "jmva.aql":
            case "jmva.lin":
            case "jmva.dmlin":
            case "jmva.ls":
            case "jmt.jmva":
            case "jmt.jmva.mva":
            case "jmt.jmva.amva":
            case "jmt.jmva.recal":
            case "jmt.jmva.comom":
            case "jmt.jmva.chow":
            case "jmt.jmva.bs":
            case "jmt.jmva.aql":
            case "jmt.jmva.lin":
            case "jmt.jmva.dmlin":
            case "jmt.jmva.ls":
                fname = writeJMVA(sn, getJMVATempPath(), this.options);
                cmd = String.format(
                        "java -cp %s jmt.commandline.Jmt mva %s -seed %s --illegal-access=permit",
                        this.jmtPath, fname, this.options.seed
                );

                if (this.options.verbose != VerboseLevel.SILENT) {
                    java.lang.System.out.println("JMT Model: " + fname);
                }
                if (this.options.verbose == VerboseLevel.DEBUG) {
                    java.lang.System.out.println("JMT Command: " + cmd);
                }
                if (options.verbose == VerboseLevel.DEBUG) {
                    cmdOutput = SysUtilsKt.system(cmd, simulationTimeoutSeconds);
                } else {
                    // Suppress system command output unless in DEBUG mode
                    java.io.ByteArrayOutputStream devNull = new java.io.ByteArrayOutputStream();
                    java.io.PrintStream nullStream = new java.io.PrintStream(devNull);
                    java.io.PrintStream originalOut = System.out;
                    java.io.PrintStream originalErr = System.err;

                    try {
                        System.setOut(nullStream);
                        System.setErr(nullStream);
                        cmdOutput = SysUtilsKt.system(cmd, simulationTimeoutSeconds);
                    } finally {
                        System.setOut(originalOut);
                        System.setErr(originalErr);
                        nullStream.close();
                    }
                }

                this.lastCommandOutput = cmdOutput; // Store for error reporting
                // Check for timeout and warn user
                if (cmdOutput.startsWith("TIMEOUT:")) {
                    line_warning(mfilename(new Object[]{}), "JMT MVA timed out. Consider using a different solver.");
                }
                if (options.verbose != VerboseLevel.SILENT && !cmdOutput.isEmpty()) {
                    java.lang.System.out.println("JMT Command output: " + cmdOutput);
                }
                runTime = java.lang.System.nanoTime() - startTime;

                solverResult = getResults();
                if (solverResult != null) {
                    solverResult.runtime = runTime / 1000000000.0;
                }
                
                if (this.options.verbose != VerboseLevel.SILENT && this.result != null) {
                    System.out.printf(
                            "%s analysis [method: %s, lang: %s, env: %s] completed in %fs.\n",
                            this.name.replaceFirst("^Solver", ""),
                            this.result.method,
                            "java",
                            System.getProperty("java.version"),
                            this.result.runtime
                    );
                    System.out.flush();
                }
                if (!this.options.keep) {
//                    try {
//                        //JMT.removeDirectory(Paths.get(this.getFilePath()));
//                    } catch (IOException ioe) {
//                        ioe.printStackTrace();
//                    }
                }
                sn = this.model.getStruct();
                AvgHandle TH = getAvgTputHandles();
                Matrix AN = snGetArvRFromTput(sn, solverResult.TN, TH);
                this.setAvgResults(solverResult.QN, solverResult.UN, solverResult.RN, solverResult.TN, AN, solverResult.WN, solverResult.CN, solverResult.XN, solverResult.runtime, options.method, 1);
                break;
            default:
                line_warning(mfilename(new Object() {
                }), "Warning: This solver does not support the specified method. Setting to default.");
                this.options.method = "default";
                runAnalyzer();
        }
    }

    public SampleResult sampleAggr(Node node, int numEvents, boolean markActivePassive) throws IOException {
        if (GlobalConstants.DummyMode) {
            return null;
        }
        
        NetworkStruct sn = this.getStruct();
        
        if (node == null) {
            throw new RuntimeException("sampleAggr requires to specify a node.");
        }
        
        if (numEvents > 0) {
            line_warning("SolverJMT", "JMT does not allow to fix the number of events for individual nodes. The number of returned events may be inaccurate.");
        }
        
        
        try {
            // Create a temp model (following the pattern used in MATLAB)
            Network modelCopy = this.model.copy();
            modelCopy.resetNetwork();
            
            // Determine the nodes to log - match MATLAB implementation
            boolean[][] isNodeClassLogged = new boolean[modelCopy.getNumberOfNodes()][modelCopy.getNumberOfClasses()];
            // Use the node index from the copied model, not the original
            int copyInd = modelCopy.getNodeIndex(node.getName());
            if (copyInd >= 0) {
                for (int r = 0; r < modelCopy.getNumberOfClasses(); r++) {
                    isNodeClassLogged[copyInd][r] = true;
                }
            }
            
            // Apply logging to the copied model (use original routing matrix like MATLAB sn.rtorig)
            RoutingMatrix Plinked = new RoutingMatrix(modelCopy, modelCopy.getClasses(), modelCopy.getNodes());
            // Copy the original routing values from sn.rtorig
            for (JobClass fromClass : sn.rtorig.keySet()) {
                for (JobClass toClass : sn.rtorig.get(fromClass).keySet()) {
                    Plinked.set(fromClass.getIndex() - 1, toClass.getIndex() - 1, sn.rtorig.get(fromClass).get(toClass));
                }
            }
            boolean[] isNodeLogged = new boolean[modelCopy.getNumberOfNodes()];
            for (int i = 0; i < modelCopy.getNumberOfNodes(); i++) {
                boolean nodeLogged = false;
                for (int j = 0; j < modelCopy.getNumberOfClasses(); j++) {
                    if (isNodeClassLogged[i][j]) {
                        nodeLogged = true;
                        break;
                    }
                }
                isNodeLogged[i] = nodeLogged;
            }
            
            String logPath = SysUtilsKt.lineTempName("jmt_sample_aggr_logs");
            modelCopy.linkAndLog(Plinked, isNodeLogged, logPath);
            
            // Force the struct to be rebuilt after linkAndLog resets it
            modelCopy.refreshStruct(false);
            
            // Simulate the model copy and retrieve log data
            SolverJMT solverjmt = new SolverJMT(modelCopy, this.getOptions());
            if (numEvents > 0) {
                // Use a more conservative multiplier to avoid excessive simulation time
                long maxEvents = Math.min(numEvents * 2, numEvents + 10000);
                solverjmt.setMaxEvents(maxEvents);
            } else {
                solverjmt.setMaxEvents(-1);
                numEvents = this.options.samples;
            }
            // Set a timeout for the sampling simulation (60 seconds by default)
            solverjmt.setSimulationTimeoutSeconds(60);
            solverjmt.getAvg(); // log data
            
            Matrix[][][] logData = parseLogs(modelCopy, isNodeLogged, MetricType.QLen);
            
            // Convert from nodes in logData to stations - match MATLAB logic
            // Get the struct that was already refreshed after linkAndLog
            NetworkStruct copySn = modelCopy.getStruct(true);
            // Get node index from both models
            int copyNodeInd = modelCopy.getNodeIndex(node.getName());
            int nodeInd = this.model.getNodeIndex(node.getName());
            
            // Get the stateful index - this maps the node to its stateful node index
            // Handle potential mismatch between original and copied model
            int isf = -1;
            
            
            // Try to get the stateful index safely
            // The nodeToStateful matrix should have one row per node, but sometimes it's not properly initialized
            // Check both models to find the correct mapping
            if (copyNodeInd >= 0 && copyNodeInd < copySn.nnodes) {
                // Check if this is a stateful node
                if (copySn.isstateful != null && copyNodeInd < copySn.isstateful.getNumRows() && 
                    copySn.isstateful.get(copyNodeInd, 0) == 1.0) {
                    
                    // Get the stateful index
                    // nodeToStateful is a row vector (1 x nnodes), not a column vector!
                    if (copySn.nodeToStateful != null && copyNodeInd < copySn.nodeToStateful.getNumCols()) {
                        isf = (int) copySn.nodeToStateful.get(0, copyNodeInd);  // Row 0, column copyNodeInd
                    } else {
                        // Fallback: use the original model's mapping
                        if (sn.nodeToStateful != null && nodeInd < sn.nodeToStateful.getNumCols()) {
                            isf = (int) sn.nodeToStateful.get(0, nodeInd);  // Row 0, column nodeInd
                        }
                    }
                }
            }
            
            // If we couldn't get a valid stateful index, the node cannot be sampled
            if (isf < 0) {
                throw new IOException("Node '" + node.getName() + "' cannot be sampled (not a stateful node or missing state mapping)");
            }
            
            Matrix t = null;
            Matrix[] nir = new Matrix[sn.nclasses];
            List<List<Object>> eventLists = new ArrayList<List<Object>>();
            for (int r = 0; r < sn.nclasses; r++) {
                eventLists.add(new ArrayList<Object>());
            }
            
            for (int r = 0; r < sn.nclasses; r++) {
                if (copyNodeInd >= logData.length || logData[copyNodeInd][r] == null || logData[copyNodeInd][r].length == 0) {
                    nir[r] = new Matrix(0, 0);
                } else {
                    if (copyNodeInd < isNodeClassLogged.length && r < isNodeClassLogged[copyNodeInd].length && isNodeClassLogged[copyNodeInd][r]) {
                        Matrix logMatrix = logData[copyNodeInd][r][0]; // Get first matrix from array
                        if (logMatrix != null && !logMatrix.isEmpty() && logMatrix.getNumCols() >= 2) {
                            // Get unique timestamps (matching MATLAB's unique behavior)
                            // Assuming log matrix has columns: [time, qlen, ...]
                            int timeColIdx = 0;
                            int qlenColIdx = 1;

                            List<Double> timeValues = new ArrayList<Double>();
                            List<Double> qlenValues = new ArrayList<Double>();

                            // Extract time and qlen data
                            for (int i = 0; i < logMatrix.getNumRows(); i++) {
                                timeValues.add(logMatrix.get(i, timeColIdx));
                                qlenValues.add(logMatrix.get(i, qlenColIdx));
                            }
                            
                            // Simple unique implementation - get first occurrence of each unique timestamp
                            List<Double> uniqueTimes = new ArrayList<Double>();
                            List<Double> uniqueQLens = new ArrayList<Double>();
                            Set<Double> seen = new HashSet<Double>();
                            
                            for (int i = 0; i < timeValues.size(); i++) {
                                double time = timeValues.get(i);
                                if (!seen.contains(time)) {
                                    seen.add(time);
                                    uniqueTimes.add(time);
                                    uniqueQLens.add(qlenValues.get(i));
                                }
                            }
                            
                            if (!uniqueTimes.isEmpty()) {
                                // Convert to matrix and shift time series like MATLAB
                                Matrix timeMatrix = new Matrix(uniqueTimes.size(), 1);
                                Matrix qlenMatrix = new Matrix(uniqueQLens.size(), 1);
                                
                                for (int i = 0; i < uniqueTimes.size(); i++) {
                                    timeMatrix.set(i, 0, uniqueTimes.get(i));
                                    qlenMatrix.set(i, 0, uniqueQLens.get(i));
                                }
                                
                                // Shift time series: t = [t(2:end); t(end)]
                                if (timeMatrix.getNumRows() > 1) {
                                    Matrix shiftedTime = new Matrix(timeMatrix.getNumRows(), 1);
                                    for (int i = 0; i < timeMatrix.getNumRows() - 1; i++) {
                                        shiftedTime.set(i, 0, timeMatrix.get(i + 1, 0));
                                    }
                                    shiftedTime.set(timeMatrix.getNumRows() - 1, 0, timeMatrix.get(timeMatrix.getNumRows() - 1, 0));
                                    t = shiftedTime;
                                    nir[r] = qlenMatrix;
                                }
                            }
                        }
                    }
                }
            }
            
            // Handle timespan filtering if finite
            if (t != null && Double.isFinite(this.options.timespan[1])) {
                int stopAt = -1;
                for (int i = 0; i < t.getNumRows(); i++) {
                    if (t.get(i, 0) > this.options.timespan[1]) {
                        stopAt = i;
                        break;
                    }
                }
                if (stopAt >= 1) {
                    Matrix newT = new Matrix(stopAt, 1);
                    for (int i = 0; i < stopAt; i++) {
                        newT.set(i, 0, t.get(i, 0));
                    }
                    t = newT;
                    
                    for (int r = 0; r < nir.length; r++) {
                        if (nir[r] != null && nir[r].getNumRows() > stopAt) {
                            Matrix newNir = new Matrix(stopAt, 1);
                            for (int i = 0; i < stopAt; i++) {
                                newNir.set(i, 0, nir[r].get(i, 0));
                            }
                            nir[r] = newNir;
                        }
                    }
                }
            }
            
            // Warning if insufficient events
            if (t != null && t.getNumRows() < 1 + numEvents) {
                line_warning("SolverJMT", "LINE could not estimate correctly the JMT simulation length to return the desired number of events at the specified node. Try to re-run increasing the number of events.");
            }
            
            // Create station state aggregate result
            SampleResult stationStateAggr = new SampleResult();
            stationStateAggr.nodeIndex = nodeInd;
            stationStateAggr.isAggregate = true;
            
            if (t != null) {
                // Limit to requested number of events + 1
                int maxLength = Math.min(t.getNumRows(), 1 + numEvents);
                Matrix finalT = new Matrix(maxLength, 1);
                for (int i = 0; i < maxLength; i++) {
                    finalT.set(i, 0, t.get(i, 0));
                }
                
                // Adjust time: [0; t(1:end-1)]
                Matrix adjustedT = new Matrix(maxLength, 1);
                adjustedT.set(0, 0, 0.0);
                for (int i = 1; i < maxLength; i++) {
                    adjustedT.set(i, 0, finalT.get(i - 1, 0));
                }
                stationStateAggr.t = adjustedT;
                
                // Combine class data into state matrix
                Matrix stateMatrix = new Matrix(maxLength, sn.nclasses);
                for (int r = 0; r < sn.nclasses; r++) {
                    if (nir[r] != null) {
                        int copyLength = Math.min(maxLength, nir[r].getNumRows());
                        for (int i = 0; i < copyLength; i++) {
                            stateMatrix.set(i, r, nir[r].get(i, 0));
                        }
                    }
                }
                stationStateAggr.state = stateMatrix;
            } else {
                // No data case
                stationStateAggr.t = new Matrix(1, 1);
                stationStateAggr.t.set(0, 0, 0.0);
                stationStateAggr.state = new Matrix(1, sn.nclasses);
            }
            
            // TODO: Process events properly - for now return empty
            stationStateAggr.event = new Matrix(0, 0);
            
            return stationStateAggr;
            
        } catch (Exception e) {
            line_warning("SolverJMT.sampleAggr", "Error: %s", e.getMessage());
            throw new IOException("Failed to sample node state: " + e.getMessage(), e);
        }
    }

    public SampleResult sampleAggr(Node node, int numEvents) throws IOException {
        return sampleAggr(node, numEvents, false);
    }

    public SampleResult sampleAggr(Node node) throws IOException {
        line_warning("SolverJMT", "JMT does not allow to fix the number of events for individual nodes. The number of returned events may be inaccurate.");
        return sampleAggr(node, this.options.samples, false);
    }

    public SampleResult sampleSysAggr(long numEvents, boolean markActivePassive) {
        if (GlobalConstants.DummyMode) {
            return null;
        }

        try {
            NetworkStruct sn = this.getStruct();
            numEvents = numEvents - 1; // Include initialization as an event

            // Use the existing model for sampling
            Network modelCopy = this.model.copy();
            modelCopy.resetNetwork();

            // Set up logging for all non-source stations
            boolean[][] isNodeClassLogged = new boolean[modelCopy.getNumberOfNodes()][modelCopy.getNumberOfClasses()];

            for (int i = 0; i < modelCopy.getNumberOfStations(); i++) {
                int nodeIndex = this.model.getNodeIndex(modelCopy.getStationNames().get(i));
                if (sn.nodetype.get(nodeIndex) != NodeType.Source) {
                    for (int r = 0; r < modelCopy.getNumberOfClasses(); r++) {
                        isNodeClassLogged[nodeIndex][r] = true;
                    }
                }
            }

            // Set up routing and logging (use original routing matrix like MATLAB sn.rtorig)
            String logPath = SysUtilsKt.lineTempName("jmt_sys_sample_logs");
            RoutingMatrix P = new RoutingMatrix(modelCopy, modelCopy.getClasses(), modelCopy.getNodes());
            // Copy the original routing values from sn.rtorig
            for (JobClass fromClass : sn.rtorig.keySet()) {
                for (JobClass toClass : sn.rtorig.get(fromClass).keySet()) {
                    P.set(fromClass.getIndex() - 1, toClass.getIndex() - 1, sn.rtorig.get(fromClass).get(toClass));
                }
            }

            // Convert boolean[][] to boolean[] for isNodeLogged
            boolean[] isNodeLogged = new boolean[modelCopy.getNumberOfNodes()];
            for (int i = 0; i < modelCopy.getNumberOfNodes(); i++) {
                boolean nodeLogged = false;
                for (int r = 0; r < modelCopy.getNumberOfClasses(); r++) {
                    if (isNodeClassLogged[i][r]) {
                        nodeLogged = true;
                        break;
                    }
                }
                isNodeLogged[i] = nodeLogged;
            }

            modelCopy.linkAndLog(P, isNodeLogged, logPath);

            // Create solver options for the copy
            SolverOptions copyOptions = this.options.copy();
            copyOptions.samples = (int)numEvents;

            // Create solver for the copy and run simulation
            SolverJMT sampleSolver = new SolverJMT(modelCopy, copyOptions);
            // Set max events proportional to samples but capped for performance
            // Use a more conservative multiplier to avoid excessive simulation time
            long maxEvents = Math.min(numEvents * 2, numEvents + 10000);
            sampleSolver.setMaxEvents(maxEvents);
            // Set a timeout for the sampling simulation (60 seconds by default)
            // This prevents getProbAggr/getProbSysAggr from hanging indefinitely
            sampleSolver.setSimulationTimeoutSeconds(60);
            sampleSolver.runAnalyzer(); // Generate log data

            // Process log data with efficient approach - parse all events together
            NetworkStruct sampleSn = modelCopy.getStruct();
            // Get initial state aggregation like MATLAB's sn_get_state_aggr
            java.util.Map<StatefulNode, Matrix> initialStateAggr = jline.api.sn.SnGetStateAggrKt.snGetStateAggr(sampleSn);
            SampleResult result = parseSystemStateFromLogs(modelCopy, isNodeLogged, sampleSn, initialStateAggr);
            return result;

        } catch (Exception e) {
            line_warning("SolverJMT.sampleSysAggr", "Error: %s", e.getMessage());
            return null;
        }
    }

    /**
     * Efficiently parses system state from logs by processing all events together.
     * This approach avoids the expensive time series merging and interpolation.
     *
     * @param model The network model
     * @param isNodeLogged Array indicating which nodes are logged
     * @param sn Network structure
     * @param initialStateAggr Map of stateful node to initial state (nir matrix)
     */
    private SampleResult parseSystemStateFromLogs(Network model, boolean[] isNodeLogged, NetworkStruct sn,
                                                   java.util.Map<StatefulNode, Matrix> initialStateAggr) {
        String logPath = model.getLogPath();
        int nstations = sn.nstations;
        int nclasses = sn.nclasses;

        // Limit the number of events to read to avoid excessive memory and time consumption
        // With 3 stations and 2 files per station (arv/dep), this allows up to 60000 total events
        final int MAX_EVENTS_PER_FILE = 10000;

        // Collect all events from all logged stations
        java.util.List<SystemEvent> allEvents = new java.util.ArrayList<>();

        for (int ist = 0; ist < nstations; ist++) {
            int nodeIndex = (int) sn.stationToNode.get(ist);
            if (nodeIndex >= isNodeLogged.length || !isNodeLogged[nodeIndex]) continue;
            if (sn.nodetype.get(nodeIndex) == NodeType.Source) continue;

            String nodeName = model.getNodeNames().get(nodeIndex);
            String logFileArv = logPath + "/" + nodeName + "-Arv.csv";
            String logFileDep = logPath + "/" + nodeName + "-Dep.csv";

            try {
                java.io.File arvFile = new java.io.File(logFileArv);
                java.io.File depFile = new java.io.File(logFileDep);

                if (arvFile.exists() && depFile.exists()) {
                    int eventsRead = 0;
                    // Parse arrivals
                    try (java.io.BufferedReader br = new java.io.BufferedReader(new java.io.FileReader(arvFile))) {
                        String line;
                        br.readLine(); // Skip header
                        while ((line = br.readLine()) != null && eventsRead < MAX_EVENTS_PER_FILE) {
                            String[] parts = line.split(";");
                            if (parts.length >= 4) {
                                double timestamp = Double.parseDouble(parts[1]);
                                String className = parts[3];
                                int classId = getClassIndex(model, className);
                                allEvents.add(new SystemEvent(timestamp, ist, classId, 1)); // +1 for arrival
                                eventsRead++;
                            }
                        }
                    }

                    eventsRead = 0;
                    // Parse departures
                    try (java.io.BufferedReader br = new java.io.BufferedReader(new java.io.FileReader(depFile))) {
                        String line;
                        br.readLine(); // Skip header
                        while ((line = br.readLine()) != null && eventsRead < MAX_EVENTS_PER_FILE) {
                            String[] parts = line.split(";");
                            if (parts.length >= 4) {
                                double timestamp = Double.parseDouble(parts[1]);
                                String className = parts[3];
                                int classId = getClassIndex(model, className);
                                allEvents.add(new SystemEvent(timestamp, ist, classId, -1)); // -1 for departure
                                eventsRead++;
                            }
                        }
                    }
                }
            } catch (Exception e) {
                line_warning("SolverJMT.parseSystemStateFromLogs", "Error parsing logs for station %d: %s", ist, e.getMessage());
            }
        }

        if (allEvents.isEmpty()) {
            return null;
        }

        // Sort all events by timestamp
        allEvents.sort((a, b) -> Double.compare(a.timestamp, b.timestamp));

        // Track current state for all stations/classes - initialize from preload
        int[][] currentState = new int[nstations][nclasses];
        // Initialize from initial state aggregation (like MATLAB's nodePreload)
        if (initialStateAggr != null) {
            for (int ist = 0; ist < nstations; ist++) {
                int isf = (int) sn.stationToStateful.get(ist);
                StatefulNode statefulNode = sn.stateful.get(isf);
                if (statefulNode != null && initialStateAggr.containsKey(statefulNode)) {
                    Matrix nirMatrix = initialStateAggr.get(statefulNode);
                    // nirMatrix is a row vector with nclasses columns
                    for (int r = 0; r < nclasses && r < nirMatrix.length(); r++) {
                        currentState[ist][r] = (int) nirMatrix.get(r);
                    }
                }
            }
        }

        // Collect sampled states - group events by timestamp to avoid
        // intermediate invalid states when multiple events happen at same time
        java.util.List<double[]> sampledTimes = new java.util.ArrayList<>();
        java.util.List<int[][]> sampledStates = new java.util.ArrayList<>();

        // Group events by timestamp and process each group together
        int eventIdx = 0;
        while (eventIdx < allEvents.size()) {
            double currentTimestamp = allEvents.get(eventIdx).timestamp;

            // Process all events at this timestamp
            while (eventIdx < allEvents.size() && allEvents.get(eventIdx).timestamp == currentTimestamp) {
                SystemEvent event = allEvents.get(eventIdx);
                if (event.classId >= 0 && event.classId < nclasses &&
                    event.stationIndex >= 0 && event.stationIndex < nstations) {
                    currentState[event.stationIndex][event.classId] += event.change;
                    if (currentState[event.stationIndex][event.classId] < 0) {
                        currentState[event.stationIndex][event.classId] = 0;
                    }
                }
                eventIdx++;
            }

            // Sample state after processing all events at this timestamp
            sampledTimes.add(new double[]{currentTimestamp});
            int[][] stateCopy = new int[nstations][nclasses];
            for (int i = 0; i < nstations; i++) {
                stateCopy[i] = currentState[i].clone();
            }
            sampledStates.add(stateCopy);
        }

        // Convert to result matrices
        int numSamples = sampledTimes.size();
        if (numSamples == 0) {
            return null;
        }

        Matrix timeMatrix = new Matrix(numSamples, 1);
        Matrix stateMatrix = new Matrix(numSamples, nstations * nclasses);

        for (int i = 0; i < numSamples; i++) {
            timeMatrix.set(i, 0, sampledTimes.get(i)[0]);
            int[][] state = sampledStates.get(i);
            for (int ist = 0; ist < nstations; ist++) {
                for (int r = 0; r < nclasses; r++) {
                    stateMatrix.set(i, ist * nclasses + r, state[ist][r]);
                }
            }
        }

        SampleResult result = new SampleResult();
        result.t = timeMatrix;
        result.state = stateMatrix;
        result.isAggregate = true;
        result.event = new Matrix(0, 0);

        return result;
    }

    // Helper class for system events
    private static class SystemEvent {
        double timestamp;
        int stationIndex;
        int classId;
        int change; // +1 for arrival, -1 for departure

        SystemEvent(double timestamp, int stationIndex, int classId, int change) {
            this.timestamp = timestamp;
            this.stationIndex = stationIndex;
            this.classId = classId;
            this.change = change;
        }
    }

    private Matrix mergeTimeSeries(Matrix time1, Matrix time2) {
        // Simple merge - combine unique time points and sort
        java.util.Set<Double> timeSet = new java.util.TreeSet<>();
        
        for (int i = 0; i < time1.getNumRows(); i++) {
            timeSet.add(time1.get(i, 0));
        }
        for (int i = 0; i < time2.getNumRows(); i++) {
            timeSet.add(time2.get(i, 0));
        }
        
        Matrix merged = new Matrix(timeSet.size(), 1);
        int idx = 0;
        for (Double time : timeSet) {
            merged.set(idx++, 0, time);
        }
        
        return merged;
    }
    
    private double interpolateState(Matrix timeData, Matrix stateData, int classIdx, double targetTime) {
        // Simple previous-value interpolation (step function)
        if (timeData == null || stateData == null || timeData.getNumRows() == 0) {
            return 0.0;
        }
        
        // Find the last time point <= targetTime
        int lastIdx = -1;
        for (int i = 0; i < timeData.getNumRows(); i++) {
            if (timeData.get(i, 0) <= targetTime) {
                lastIdx = i;
            } else {
                break;
            }
        }
        
        if (lastIdx >= 0 && lastIdx < stateData.getNumRows() && classIdx < stateData.getNumCols()) {
            return stateData.get(lastIdx, classIdx);
        } else {
            return 0.0;
        }
    }

    public SampleResult sampleSysAggr(long numEvents) {
        boolean markActivePassive = false;
        return sampleSysAggr(numEvents, markActivePassive);
    }

    public SampleResult sampleSysAggr() {
        long numEvents = this.options.samples;
        return sampleSysAggr(numEvents);
    }

    /**
     * Gets transient probability for a specific node's aggregated state.
     * Currently not fully implemented - returns empty result.
     *
     * @param node The node of interest
     * @return TransientProbabilityResult containing time points and probabilities
     */
    public JMTResult.TransientProbabilityResult getTranProbAggr(Node node) {
        if (GlobalConstants.DummyMode) {
            return new JMTResult.TransientProbabilityResult();
        }
        
        if (node == null) {
            throw new RuntimeException("getTranProbAggr requires to specify a node.");
        }
        
        if (!Double.isFinite(this.options.timespan[1])) {
            throw new RuntimeException("getTranProbAggr in SolverJMT requires to specify a finite timespan T, e.g., SolverJMT(model, \"timespan\", new double[]{0, T}).");
        }
        
        // Currently not implemented in MATLAB either - throws "Method not implemented yet"
        throw new RuntimeException("Method not implemented yet.");
    }

    /**
     * Gets probability of the current system state in aggregated form.
     * Uses simulation sampling to estimate the probability.
     *
     * @return Probability of the current system state
     */
    public ProbabilityResult getProbSysAggr() {
        if (GlobalConstants.DummyMode) {
            return new ProbabilityResult(Double.NaN);
        }

        try {
            NetworkStruct sn = getStruct();

            // Get system state samples
            SampleResult tranSysStateAggr = this.sampleSysAggr();
            if (tranSysStateAggr == null || tranSysStateAggr.state == null) {
                line_warning("SolverJMT.getProbSysAggr", "Unable to extract state samples from JMT simulation. This feature requires simulation logging which may not be supported for all model types.");
                return new ProbabilityResult(0.0);
            }

            // Build time-state matrix similar to MATLAB's TSS
            int numSamples = tranSysStateAggr.t.length();
            int numStatefulNodes = (int) sn.nstateful;
            int numClasses = (int) sn.nclasses;

            // Calculate current system state in aggregated form (nir format)
            Matrix currentNir = new Matrix(numStatefulNodes, numClasses);
            for (int isf = 0; isf < numStatefulNodes; isf++) {
                int ind = (int) sn.statefulToNode.get(isf);
                StatefulNode statefulNode = sn.stateful.get(isf);
                Matrix nodeState = sn.state.get(statefulNode);
                if (nodeState != null) {
                    State.StateMarginalStatistics stats = ToMarginal.toMarginal(sn, ind, nodeState, null, null, null, null, null);
                    if (stats != null && stats.nir != null) {
                        for (int r = 0; r < numClasses; r++) {
                            currentNir.set(isf, r, stats.nir.get(0, r));
                        }
                    }
                }
            }

            // Calculate time differences (duration each state was observed)
            Matrix timeDiffs = new Matrix(numSamples, 1);
            for (int i = 0; i < numSamples - 1; i++) {
                timeDiffs.set(i, 0, tranSysStateAggr.t.get(i + 1, 0) - tranSysStateAggr.t.get(i, 0));
            }
            timeDiffs.set(numSamples - 1, 0, 0.0); // Last sample has 0 duration

            double totalTime = timeDiffs.elementSum();
            double matchingTime = 0.0;

            // state is a Matrix with rows=samples, columns=station*class combinations
            Matrix stateMatrix = (Matrix) tranSysStateAggr.state;

            // Find samples that match current system state
            for (int sample = 0; sample < numSamples; sample++) {
                boolean matches = true;

                // Check if this sample matches current state
                for (int isf = 0; isf < numStatefulNodes && matches; isf++) {
                    for (int r = 0; r < numClasses && matches; r++) {
                        double sampleValue = stateMatrix.get(sample, isf * numClasses + r);
                        double currentValue = currentNir.get(isf, r);
                        if (Math.abs(sampleValue - currentValue) > 1e-10) {
                            matches = false;
                        }
                    }
                }

                if (matches) {
                    matchingTime += timeDiffs.get(sample, 0);
                }
            }

            if (totalTime > 0) {
                double prob = matchingTime / totalTime;
                return new ProbabilityResult(prob);
            } else {
                line_warning("SolverJMT", "The state was not seen during the simulation.");
                return new ProbabilityResult(0.0);
            }

        } catch (Exception e) {
            line_warning("SolverJMT.getProbSysAggr", "Error: %s", e.getMessage());
            return new ProbabilityResult(0.0);
        }
    }

    // XML serialization methods have been moved to SaveHandlers class
    // Use getSaveHandlers().saveXxx() for all XML generation operations

    private Matrix setValues(Matrix matrix, List<Integer> istStations, List<Integer> rList, double value) {
        for (int i : istStations) {
            for (int r : rList) {
                matrix.set(i, r, value);
            }
        }
        return matrix;
    }

    @Override
    public boolean supports(Network model) {
        FeatureSet featUsed = model.getUsedLangFeatures();
        FeatureSet featSupported = SolverJMT.getFeatureSet();
        return FeatureSet.supports(featSupported, featUsed);
    }

    /*
     * This method writes the model to a JSIM file and returns the path to the file.
     * The file is written to the directory specified by the logPath attribute of the model.
     * The file name is the name of the model.
     * The file is written in the JSIM format.
     * @return the path to the JSIM file
     */
    public String writeJSIM(NetworkStruct sn, String outputFileName) throws ParserConfigurationException {
        // Update SaveHandlers with the current network structure
        getSaveHandlers().updateNetworkStruct(sn);
        
        ElementDocumentPair xml = getSaveHandlers().saveXMLHeader(this.model.getLogPath());
        xml = getSaveHandlers().saveClasses(xml);
        int numOfClasses = sn.nclasses;
        int numOfNodes = sn.nnodes;

        for (int i = 0; i < numOfNodes; i++) {
            Node currentNode = this.model.getNodes().get(i);
            Element node = xml.simDoc.createElement("node");
            node.setAttribute("name", currentNode.getName());
            List<Section> nodeSections = Arrays.asList(currentNode.getInput(), currentNode.getServer(), currentNode.getOutput());
            for (int j = 0; j < nodeSections.size(); j++) {
                Element xml_section = xml.simDoc.createElement("section");
                Section currentSection = nodeSections.get(j);
                // For Logger nodes, we need to include Generic sections too
                boolean isLogger = currentNode instanceof Logger;
                if (currentSection != null && (isLogger || !currentSection.getClassName().startsWith("Generic "))) {
                    xml_section.setAttribute("className", currentSection.getClassName());

                    // Override className for preemptive strategies - JMT requires PreemptiveServer
                    if (currentSection.getClassName().equals("Server") && currentNode instanceof jline.lang.nodes.Queue) {
                        jline.lang.nodes.Queue queue = (jline.lang.nodes.Queue) currentNode;
                        SchedStrategy sched = queue.getSchedStrategy();
                        if (sched == SchedStrategy.SRPT || sched == SchedStrategy.SRPTPRIO ||
                            sched == SchedStrategy.LCFSPR || sched == SchedStrategy.LCFSPRPRIO ||
                            sched == SchedStrategy.LCFSPI || sched == SchedStrategy.LCFSPIPRIO ||
                            sched == SchedStrategy.FCFSPR || sched == SchedStrategy.FCFSPRPRIO ||
                            sched == SchedStrategy.FCFSPI || sched == SchedStrategy.FCFSPIPRIO) {
                            xml_section.setAttribute("className", "PreemptiveServer");
                        }
                    }

                    DocumentSectionPair simXML = new DocumentSectionPair(xml.simDoc, xml_section);
                    switch (currentSection.getClassName()) {
                        case "Buffer":
                            simXML.section.setAttribute("className", "Queue"); // overwrite with JMT class name
                            simXML = getSaveHandlers().saveBufferCapacity(simXML, i);
                            simXML = getSaveHandlers().saveDropStrategy(simXML, i); // unfinished
                            simXML = getSaveHandlers().saveGetStrategy(simXML, i);
                            simXML = getSaveHandlers().savePutStrategy(simXML, i);
                            simXML = getSaveHandlers().saveImpatience(simXML, i);
                            break;
                        case "Server":
                        case "jline.Server":
                        case "PreemptiveServer":
                            simXML = getSaveHandlers().saveNumberOfServers(simXML, i);
                            simXML = getSaveHandlers().saveServerVisits(simXML);
                            simXML = getSaveHandlers().saveServiceStrategy(simXML, i);
                            simXML = getSaveHandlers().saveDelayOffStrategy(simXML, i);
                            break;
                        case "PollingServer":
                            simXML = getSaveHandlers().setPollingServerClassName(simXML, i);
                            simXML = getSaveHandlers().saveNumberOfServers(simXML, i);
                            simXML = getSaveHandlers().saveServerVisits(simXML);
                            simXML = getSaveHandlers().saveServiceStrategy(simXML, i);
                            simXML = getSaveHandlers().saveSwitchoverStrategy(simXML, i);
                            break;
                        case "SharedServer":
                            simXML.section.setAttribute("className", "PSServer"); // overwrite with JMT class name
                            simXML = getSaveHandlers().saveNumberOfServers(simXML, i);
                            simXML = getSaveHandlers().saveServerVisits(simXML);
                            simXML = getSaveHandlers().saveServiceStrategy(simXML, i);
                            simXML = getSaveHandlers().saveDelayOffStrategy(simXML, i);
                            simXML = getSaveHandlers().savePreemptiveStrategy(simXML, i);
                            simXML = getSaveHandlers().savePreemptiveWeights(simXML, i);
                            break;
                        case "InfiniteServer":
                        case "jline.InfiniteServer":
                            simXML.section.setAttribute("className", "Delay"); // overwrite with JMT class name
                            simXML = getSaveHandlers().saveServiceStrategy(simXML, i);
                            break;
                        case "RandomSource":
                        case "jline.RandomSource":
                            simXML = getSaveHandlers().saveArrivalStrategy(simXML, i);
                            break;
                        case "Dispatcher":
                        case "jline.Dispatcher":
                        case "ClassSwitchDispatcher":
                            simXML.section.setAttribute("className", "Router"); // overwrite with JMT class name
                            simXML = getSaveHandlers().saveRoutingStrategy(simXML, i);
                            break;
                        case "StatelessClassSwitcher":
                        case "jline.StatelessClassSwitcher":
                            simXML.section.setAttribute("className", "ClassSwitch"); // overwrite with JMT class name
                            simXML = getSaveHandlers().saveClassSwitchStrategy(simXML, i);
                            break;
                        case "Cache":
                            simXML.section.setAttribute("className", "Cache"); // overwrite with JMT class name
                            //System.out.println("simXML: " + simXML + "; i:" + i);
                            simXML = getSaveHandlers().saveCacheStrategy(simXML, i);
                            break;
                        case "LogTunnel":
                            simXML = getSaveHandlers().saveLogTunnel(simXML, i);
                            break;
                        case "Joiner":
                            simXML.section.setAttribute("className", "Join"); // overwrite with JMT class name
                            simXML = getSaveHandlers().saveJoinStrategy(simXML, i);
                            break;
                        case "Forker":
                            simXML.section.setAttribute("className", "Fork"); // overwrite with JMT class name
                            simXML = getSaveHandlers().saveForkStrategy(simXML, i);
                            break;
                        case "Storage":
                            simXML.section.setAttribute("className", "Storage"); // overwrite with JMT class name
                            simXML = getSaveHandlers().saveTotalCapacity(simXML, i);
                            simXML = getSaveHandlers().savePlaceCapacities(simXML, i);
                            simXML = getSaveHandlers().saveDropRule(simXML, i);
                            simXML = getSaveHandlers().saveGetStrategy(simXML);
                            simXML = getSaveHandlers().savePutStrategies(simXML, i);
                            break;
                        case "Enabling":
                            simXML.section.setAttribute("className", "Enabling"); // overwrite with JMT class name
                            simXML = getSaveHandlers().saveEnablingConditions(simXML, i);
                            simXML = getSaveHandlers().saveInhibitingConditions(simXML, i);
                            break;
                        case "Firing":
                            simXML.section.setAttribute("className", "Firing"); // overwrite with JMT class name
                            simXML = getSaveHandlers().saveFiringOutcomes(simXML, i);
                            break;
                        case "Timing":
                            simXML.section.setAttribute("className", "Timing"); // overwrite with JMT class name
                            simXML = getSaveHandlers().saveModeNames(simXML, i);
                            simXML = getSaveHandlers().saveNumbersOfServers(simXML, i);
                            simXML = getSaveHandlers().saveTimingStrategies(simXML, i);
                            simXML = getSaveHandlers().saveFiringPriorities(simXML, i);
                            simXML = getSaveHandlers().saveFiringWeights(simXML, i);
                            break;
                        case "Generic Input":
                            // For Logger nodes, create a Queue section for the input
                            if (currentNode instanceof Logger) {
                                simXML.section.setAttribute("className", "Queue");
                                simXML = getSaveHandlers().saveBufferCapacity(simXML, i);
                                simXML = getSaveHandlers().saveDropStrategy(simXML, i);
                                simXML = getSaveHandlers().saveGetStrategy(simXML, i);
                                simXML = getSaveHandlers().savePutStrategy(simXML, i);
                                simXML = getSaveHandlers().saveImpatience(simXML, i);
                            }
                            break;
                        case "Generic Output":
                            // For Logger nodes, output section is handled by Router section below
                            if (currentNode instanceof Logger) {
                                continue; // Skip adding this section, Router will be added separately
                            }
                            break;
                    }
                    node.appendChild(simXML.section);
                }
            }
            
            // Special handling for Logger nodes - they need Router sections even with Generic Output
            if (currentNode instanceof Logger) {
                Element routerSection = xml.simDoc.createElement("section");
                routerSection.setAttribute("className", "Router");
                DocumentSectionPair routerXML = new DocumentSectionPair(xml.simDoc, routerSection);
                routerXML = getSaveHandlers().saveRoutingStrategy(routerXML, i);
                node.appendChild(routerXML.section);
            }
            
            xml.simElem.appendChild(node);
        }
        xml = getSaveHandlers().saveMetrics(xml);
        xml = getSaveHandlers().saveLinks(xml);
        xml = getSaveHandlers().saveRegions(xml);

        boolean hasReferenceNodes = false;
        Element preloadNode = xml.simDoc.createElement("preload");
        Map<StatefulNode, Matrix> s0 = sn.state;
        int numOfStations = sn.nstations;

        for (int i = 0; i < numOfStations; i++) {
            boolean isReferenceNode = false;
            int nodeIndex = (int) sn.stationToNode.get(i);
            int isf = (int) sn.stationToStateful.get(i);
            Element stationPopulationsNode = null;

            if (sn.nodetype.get(nodeIndex) != NodeType.Source && sn.nodetype.get(nodeIndex) != NodeType.Join) {
                State.StateMarginalStatistics sms = ToMarginal.toMarginal(sn, nodeIndex, s0.get(this.model.getStatefulNodes().get(isf)), null, null, null, null, null);
                stationPopulationsNode = xml.simDoc.createElement("stationPopulations");
                stationPopulationsNode.setAttribute("stationName", sn.nodenames.get(nodeIndex));

                // TODO: here we assume that the current state is the one on top of the state data structure
                // however it could be that stateprior places the mass on another (or multiple other) states
                for (int r = 0; r < numOfClasses; r++) {
                    Element classPopulationNode = xml.simDoc.createElement("classPopulation");

                    if (Double.isInfinite(sn.njobs.get(r)) || sn.njobs.get(r) == Integer.MAX_VALUE) {
                        isReferenceNode = true;
                        classPopulationNode.setAttribute("population", String.valueOf(Math.round(sms.nir.get(0, r))));
                        classPopulationNode.setAttribute("refClass", sn.classnames.get(r));
                        stationPopulationsNode.appendChild(classPopulationNode);
                    } else {
                        isReferenceNode = true;
                        classPopulationNode.setAttribute("population", String.valueOf(Math.round(sms.nir.get(0, r))));
                        classPopulationNode.setAttribute("refClass", sn.classnames.get(r));
                        stationPopulationsNode.appendChild(classPopulationNode);
                    }
                }
            }

            if (isReferenceNode) {
                preloadNode.appendChild(stationPopulationsNode);
            }
            hasReferenceNodes = hasReferenceNodes || isReferenceNode;
        }

        if (hasReferenceNodes) {
            xml.simElem.appendChild(preloadNode);
        }
        try {
            InputOutputKt.writeXML(outputFileName, xml.simDoc);
        } catch (Exception e) {
            e.printStackTrace();
            try {
                InputOutputKt.writeXML(outputFileName, xml.simDoc);
            } catch (Exception retryException) {
                retryException.printStackTrace();
            }
        }
        return outputFileName;
    }

    public String writeJSIM(NetworkStruct sn) throws ParserConfigurationException {
        String outputFileName = getJSIMTempPath();
        return writeJSIM(sn, outputFileName);
    }

    /**
     * Writes queueing network model to JMT JSIMG format.
     * Delegates to the standalone QN2JSIMG class in the io package.
     * @param sn the network structure
     * @param outputFileName the output file name
     * @return the path to the JSIM file
     * @throws ParserConfigurationException if XML parsing fails
     */
    public String QN2JSIMG(NetworkStruct sn, String outputFileName) throws ParserConfigurationException {
        return jline.io.QN2JSIMG.writeJSIM(this.model, sn, outputFileName, getSaveHandlers());
    }

    /**
     * Writes queueing network model to JMT JSIMG format.
     * Delegates to the standalone QN2JSIMG class in the io package.
     * @param sn the network structure
     * @return the path to the JSIM file
     * @throws ParserConfigurationException if XML parsing fails
     */
    public String QN2JSIMG(NetworkStruct sn) throws ParserConfigurationException {
        return jline.io.QN2JSIMG.writeJSIM(this.model, sn, null, getSaveHandlers());
    }


    public enum ViewMode {
        JSIMW,
        JSIMG
    }
    
    public static class EventInfo {
        public int node;
        public int jobclass;
        public double t;
    }
}
