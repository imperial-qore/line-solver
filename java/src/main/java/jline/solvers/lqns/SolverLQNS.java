package jline.solvers.lqns;

import jline.lang.constant.GlobalConstants;
import jline.lang.constant.SolverType;
import jline.lang.constant.VerboseLevel;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.LayeredNetworkElement;
import jline.lang.layered.LayeredNetworkStruct;
import jline.solvers.*;

import java.io.*;
import java.util.*;
import java.util.concurrent.TimeUnit;
import javax.xml.parsers.*;

import jline.util.Matrix;
import org.w3c.dom.*;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

import static jline.io.SysUtils.lineTempName;

// LayeredNetworkSolver not available in Java
public class SolverLQNS extends Solver {
    public SolverLQNS(LayeredNetwork lqnmodel) {
        this(lqnmodel, new SolverOptions(SolverType.LQNS));
    }

    public SolverLQNS(LayeredNetwork lqnmodel, SolverOptions options) {
        super(lqnmodel, "SolverLQNS", options);
        //setOptions(Solver.parseOptions(varargin, defaultOptions()));
        if (!isAvailable()) {
            throw new RuntimeException("SolverLQNS requires the lqns and lqsim commands to be available on the system path. Please visit: http://www.sce.carleton.ca/rads/lqns/");
        }
    }

    protected void runAnalyzer() throws IllegalAccessException, ParserConfigurationException {
        this.runAnalyzer(this.options);
    }

    protected void runAnalyzer(SolverOptions options) throws IllegalAccessException, ParserConfigurationException {
        // TODO: add runtime
        long startTimeMillis = System.currentTimeMillis();
        if (options == null) {
            options = this.options;  // Implement getOptions to retrieve default or previously set options
        }

        String dirpath = null;
        try {
            dirpath = lineTempName("lqns");
        } catch (IOException e) {
            return;
        }
        String filename = dirpath + File.separator + "model.lqnx";
        ((LayeredNetwork) this.model).writeXML(filename, false);

        this.resetRandomGeneratorSeed(options.seed);

        String verbose = options.verbose == VerboseLevel.SILENT ? "" : "-a -w";
        String multiserver_praqma = ""; //getMultiserverPraqma(options);  // TODO: Implement this method based on the switch case in MATLAB

        String cmd = "lqns " + verbose + " " + multiserver_praqma + " " + filename;

        if (options.verbose == VerboseLevel.DEBUG) {
            System.out.println("\nLQNS command: " + cmd);
        }
        try {
            Process process = Runtime.getRuntime().exec(cmd);
            int status = process.waitFor();
        } catch (IOException e) {
            return;
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }
        try {
            this.parseXMLResults(filename);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        if (!options.keep) {
            // Remove the directory
            new File(dirpath).delete();
        }
        long endTimeMillis = System.currentTimeMillis();
        result.runtime = (endTimeMillis - startTimeMillis) / 1000.0;

    }

    public LayeredNetworkStruct getStruct() {
        return ((LayeredNetwork) this.model).getStruct();
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
                        System.err.println("Error parsing iterations value: " + iterationsStr);
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
                int procPos = 0;
                for (Map.Entry<Integer, String> entry : lqn.names.entrySet()) {
                    if (Objects.equals(procName, entry.getValue())) {
                        procPos = entry.getKey().intValue();
                    }
                }
                NodeList procResultList = procElement.getElementsByTagName("result-processor");
                if (procResultList.getLength() > 0) {
                    Element procResultElement = (Element) procResultList.item(0);
                    String utilizationStr = procResultElement.getAttribute("utilization");
                    double uRes = 0.0;
                    if (!utilizationStr.isEmpty()) {
                        uRes = Double.parseDouble(utilizationStr);
                    }
                    AvgNodesProcUtilization.set(procPos - 1, uRes);
                }

                // Assuming procElement is already defined as an Element and lqn.names is a List or array.
                NodeList taskList = procElement.getElementsByTagName("task");
                for (int j = 0; j < taskList.getLength(); j++) {
                    // Element - Task
                    Element taskElement = (Element) taskList.item(j);
                    String taskName = taskElement.getAttribute("name");

                    // Find the position of taskName
                    int taskPos = 0;
                    for (Map.Entry<Integer, String> entry : lqn.names.entrySet()) {
                        if (Objects.equals(taskName, entry.getValue())) {
                            taskPos = entry.getKey().intValue();
                        }
                    }
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
                    AvgNodesUtilization.set(taskPos - 1, uRes);
                    AvgNodesPhase1Utilization.set(taskPos - 1, p1uRes);
                    AvgNodesThroughput.set(taskPos - 1, tRes);
                    AvgNodesProcUtilization.set(taskPos - 1, puRes);
                    AvgNodesPhase2Utilization.set(taskPos - 1, p2uRes);

                    NodeList entryList = doc.getElementsByTagName("entry");
                    for (int k = 0; k < entryList.getLength(); k++) {
                        Element entryElement = (Element) entryList.item(k);
                        String entryName = entryElement.getAttribute("name");
                        // Find the position of entryName
                        int entryPos = 0;
                        for (Map.Entry<Integer, String> entry : lqn.names.entrySet()) {
                            if (Objects.equals(entryName, entry.getValue())) {
                                entryPos = entry.getKey().intValue();
                            }
                        }
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
                            p2stRes = Double.parseDouble(p2uResStr);
                        }
                        tRes = Double.parseDouble(firstEntryResult.getAttribute("throughput"));
                        puRes = Double.parseDouble(firstEntryResult.getAttribute("proc-utilization"));

                        AvgNodesUtilization.set(entryPos - 1, uRes);
                        AvgNodesPhase1Utilization.set(entryPos - 1, p1uRes);
                        AvgNodesPhase2Utilization.set(entryPos - 1, p2uRes);
                        AvgNodesPhase1ServiceTime.set(entryPos - 1, p1stRes);
                        AvgNodesPhase2ServiceTime.set(entryPos - 1, p2stRes);
                        AvgNodesThroughput.set(entryPos - 1, tRes);
                        AvgNodesProcUtilization.set(entryPos - 1, puRes);
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
                                // TODO: there may be a performance bug here as actName seems to visit more than once
                                //  the same activity

                                // Find the position of actName
                                int actPos = 0;
                                for (Map.Entry<Integer, String> entry : lqn.names.entrySet()) {
                                    if (Objects.equals(actName, entry.getValue())) {
                                        actPos = entry.getKey().intValue();
                                    }
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

                                AvgNodesUtilization.set(actPos - 1, uRes);
                                AvgNodesPhase1ServiceTime.set(actPos - 1, stRes);
                                AvgNodesThroughput.set(actPos - 1, tRes);
                                AvgNodesProcWaiting.set(actPos - 1, pwRes);
                                AvgNodesProcUtilization.set(actPos - 1, mypuRes);

                                String actID = lqn.names.get(actPos);
                                // Synchronous calls
                                NodeList synchCalls = actElement.getElementsByTagName("synch-call");
                                for (int m = 0; m < synchCalls.getLength(); m++) {
                                    Element callElement = (Element) synchCalls.item(m);
                                    String destName = callElement.getAttribute("dest");
                                    int destPos = 0;
                                    for (Map.Entry<Integer, String> entry : lqn.names.entrySet()) {
                                        if (Objects.equals(destName, entry.getValue())) {
                                            destPos = entry.getKey().intValue();
                                        }
                                    }
                                    String destID = lqn.names.get(destPos);
                                    int callPos = 0;

                                    for (Map.Entry<Integer, String> entry : lqn.callnames.entrySet()) {
                                        if (Objects.equals(Arrays.asList(actID + "=>" + destID).get(0), entry.getValue())) {
                                            callPos = entry.getKey().intValue();
                                        }
                                    }
                                    NodeList callResult = callElement.getElementsByTagName("result-call");
                                    double wRes = Double.parseDouble(((Element) callResult.item(0)).getAttribute("waiting"));
                                    AvgEdgesWaiting.set(callPos - 1, wRes);
                                }

                                // Asynchronous calls
                                NodeList asynchCalls = actElement.getElementsByTagName("asynch-call");
                                for (int m = 0; m < asynchCalls.getLength(); m++) {
                                    Element callElement = (Element) asynchCalls.item(m);
                                    String destName = callElement.getAttribute("dest");
                                    int destPos = 0;
                                    for (Map.Entry<Integer, String> entry : lqn.names.entrySet()) {
                                        if (Objects.equals(destName, entry.getValue())) {
                                            destPos = entry.getKey().intValue();
                                        }
                                    }
                                    String destID = lqn.names.get(destPos);
                                    int callPos = 0;
                                    for (Map.Entry<Integer, String> entry : lqn.callnames.entrySet()) {
                                        if (Objects.equals(Arrays.asList(actID + "->" + destID).get(0), entry.getValue())) {
                                            callPos = entry.getKey().intValue();
                                        }
                                    }
                                    NodeList callResult = callElement.getElementsByTagName("result-call");
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

        ((LayeredSolverResult) this.result).PN = new Matrix(AvgNodesProcUtilization);
        ((LayeredSolverResult) this.result).SN = new Matrix(AvgNodesPhase1ServiceTime);
        ((LayeredSolverResult) this.result).TN = new Matrix(AvgNodesThroughput);
        ((LayeredSolverResult) this.result).UN = new Matrix(AvgNodesUtilization);
        ((LayeredSolverResult) this.result).RN = new Matrix(AvgNodesProcWaiting);
        ((LayeredSolverResult) this.result).RN.fill(Double.NaN);
        ((LayeredSolverResult) this.result).QN = new Matrix(AvgNodesProcWaiting);
        ((LayeredSolverResult) this.result).QN.fill(Double.NaN);
        ((LayeredSolverResult) this.result).AN = new Matrix(lqn.nidx, 1);
        ((LayeredSolverResult) this.result).WN = new Matrix(lqn.nidx, 1);
        return result;
    }

    // Static methods
    public static List<String> listValidMethods() {
        // Implementation of listValidMethods
        return Arrays.asList("default", "lqns", "srvn", "exactmva", "srvn.exactmva", "sim", "lqsim", "lqnsdefault");
    }

    public static boolean isAvailable() {
        boolean isAvailable = false;
        try {
            Process process = Runtime.getRuntime().exec("lqsim -V -H"); //lqns seems to return exit code 4
            process.waitFor();
            if (process.exitValue() == 0) {
                isAvailable = true;
            }
        } catch (IOException e) {
            return isAvailable;
        } catch (InterruptedException e) {
            return isAvailable;
        }
        return isAvailable;
    }

    public SolverResult getAvg() {
        return this.getEnsembleAvg();
    }

    protected SolverResult getEnsembleAvg() {
        // TODO: add useLQNSnaming

        boolean useLQNSnaming = false;
        try {
            this.runAnalyzer();
        } catch (IllegalAccessException | ParserConfigurationException e) {
            throw new RuntimeException(e);
        }

        LayeredSolverResult lqnsResult = (LayeredSolverResult) this.result.deepCopy();
        if (!useLQNSnaming) {
            lqnsResult.QN = ((LayeredSolverResult) this.result).UN.clone();
            lqnsResult.UN = ((LayeredSolverResult) this.result).PN.clone();
            lqnsResult.RN = ((LayeredSolverResult) this.result).SN.clone();
        }

        return lqnsResult;
    }

    public final LayeredNetworkAvgTable getAvgTable() {
        LayeredSolverResult result = (LayeredSolverResult) this.getAvg();
        LayeredNetworkStruct lqn = this.getStruct();

        List<String> nodeNames = new ArrayList<>(lqn.names.values());
        List<String> nodeTypes = new ArrayList<>();
        for (int o = 0; o < nodeNames.size(); o++) {
            switch ((int) lqn.type.get(1 + o)) {
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
            }
        }
        LayeredNetworkAvgTable AvgTable = new LayeredNetworkAvgTable(result.QN.toList1D(), result.UN.toList1D(), result.RN.toList1D(), result.WN.toList1D(), result.TN.toList1D());
        AvgTable.setNodeNames(nodeNames);
        AvgTable.setNodeTypes(nodeTypes);
        AvgTable.setOptions(this.options);
        return AvgTable;
    }

    public static void main(String[] args) throws Exception {
        GlobalConstants.getInstance().setVerbose(VerboseLevel.DEBUG);
        LayeredNetwork lqnmodel = jline.examples.LayeredExamples.test0();
        SolverLQNS s = new SolverLQNS(lqnmodel);
        LayeredNetworkAvgTable AvgTable = s.getAvgTable();
        AvgTable.printTable();
    }
}
