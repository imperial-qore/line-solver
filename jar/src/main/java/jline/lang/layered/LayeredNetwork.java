/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;

import jline.GlobalConstants;
import static jline.GlobalConstants.Inf;

import jline.cli.LineDockerClient;
import jline.io.SysUtils;
import jline.io.ModelVisualizer;
import jline.lang.Copyable;
import jline.lang.Ensemble;
import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.CallType;
import jline.lang.constant.ActivityPrecedenceType;
import jline.lang.constant.ProcessType;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.*;
import jline.util.matrix.MatrixCell;
import jline.solvers.ln.SolverLN;
import jline.util.graph.DirectedGraph;
import jline.util.Utils;
import jline.util.matrix.Matrix;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.ArrayList;
import java.util.List;
import java.util.Collections;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

import static jline.api.lsn.LsnMaxMultiplicity.lsnMaxMultiplicity;
import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.line_warning;
import static jline.io.InputOutput.mfilename;
import static jline.io.SysUtils.jlqnGetPath;
import static jline.io.SysUtils.lineTempName;
import static jline.lang.constant.ActivityPrecedenceType.*;
import jline.VerboseLevel;

/**
 * LayeredNetwork represents a layered queueing network (LQN) model for performance analysis
 * of distributed and multi-tiered software systems.
 * 
 * <p>Layered queueing networks extend traditional queueing networks by modeling software
 * systems where processes can both serve requests and make requests to other processes,
 * creating a layered architecture. This is particularly useful for analyzing:
 * <ul>
 * <li>Multi-tier web applications (web server, application server, database)</li>
 * <li>Service-oriented architectures and microservices</li>
 * <li>Client-server systems with nested service calls</li>
 * <li>Cloud computing and distributed systems</li>
 * </ul>
 * 
 * <p>The network consists of:
 * <ul>
 * <li><b>Hosts:</b> Physical or virtual processors that execute tasks</li>
 * <li><b>Tasks:</b> Software processes that can serve requests and make calls</li>
 * <li><b>Entries:</b> Service interfaces exposed by tasks</li>
 * <li><b>Activities:</b> Individual processing steps within tasks</li>
 * <li><b>Precedences:</b> Execution order and control flow relationships</li>
 * </ul>
 * 
 * <p>LayeredNetwork supports both open models (with external arrivals) and closed models
 * (with fixed population), and can model complex interactions including synchronous calls,
 * asynchronous messaging, fork-join parallelism, and probabilistic routing.
 * 
 * <p>The model can be solved using various algorithms including Mean Value Analysis (MVA),
 * simulation, and matrix-analytic methods to obtain performance metrics such as response
 * times, throughputs, and resource utilizations.
 * 
 * @see Task
 * @see Entry
 * @see Activity
 * @see ActivityPrecedence
 * @see Host
 * @see Processor
 */
public class LayeredNetwork extends Ensemble implements Copyable {

    /**
     * Which solvers and solver methods can analyze THIS layered model.
     *
     * <pre>
     * model.findSolver()             every (solver, method) pair that runs
     * model.findSolver("tran", false)  ... that returns transients
     * model.findSolver("", true)     also the pairs that are refused, and why
     * </pre>
     *
     * <p>One row per pair; see {@link jline.solvers.auto.SolverCandidate} for
     * the columns and {@code SolverCandidate.toTable} to print them.
     *
     * <p>The families in play are the layered ones, "ln" and "lqns": the flat
     * Network families describe what they accept INSIDE a layer, so answering
     * with them would answer a question that was not asked.
     *
     * <p>{@link #findMethod()} and {@link #help()} are aliases.
     *
     * @return one row per runnable (family, method) pair
     */
    public java.util.List<jline.solvers.auto.SolverCandidate> findSolver() {
        return findSolver("", false);
    }

    /**
     * Which solvers and solver methods can analyze this layered model, narrowed
     * to one measure and optionally including the refused pairs.
     *
     * @param metric  a measure group ("tran") or the accessor that returns it
     *                ("getTranAvg"); "" or "any" keeps every pair
     * @param showAll keep the refused pairs too, with the reason each was refused
     * @return the matching rows
     */
    public java.util.List<jline.solvers.auto.SolverCandidate> findSolver(String metric,
                                                                        boolean showAll) {
        return jline.solvers.auto.SolverAUTO.findSolverLayered(this, metric, showAll);
    }

    /**
     * Alias of {@link #findSolver()}.
     *
     * @return one row per runnable (family, method) pair
     */
    public java.util.List<jline.solvers.auto.SolverCandidate> findMethod() {
        return findSolver("", false);
    }

    /**
     * Alias of {@link #findSolver(String, boolean)}.
     *
     * @param metric  a measure group or the accessor that returns it
     * @param showAll keep the refused pairs too
     * @return the matching rows
     */
    public java.util.List<jline.solvers.auto.SolverCandidate> findMethod(String metric,
                                                                        boolean showAll) {
        return findSolver(metric, showAll);
    }

    /**
     * Alias of {@link #findSolver()}: what can this model be solved with?
     *
     * @return one row per runnable (family, method) pair
     */
    public java.util.List<jline.solvers.auto.SolverCandidate> help() {
        return findSolver("", false);
    }

    /**
     * Alias of {@link #findSolver(String, boolean)}.
     *
     * @param metric  a measure group or the accessor that returns it
     * @param showAll keep the refused pairs too
     * @return the matching rows
     */
    public java.util.List<jline.solvers.auto.SolverCandidate> help(String metric,
                                                                  boolean showAll) {
        return findSolver(metric, showAll);
    }

    private final Param param;
    protected Map<Integer, Host> hosts;
    protected Map<Integer, Task> tasks;
    protected Map<Integer, Task> reftasks;
    protected Map<Integer, Activity> activities;
    protected Map<Integer, Entry> entries;
    protected Map<Integer, LayeredNetworkElement> nodes;
    //private final Aux aux;
    private Matrix lqnGraph;
    private Matrix taskGraph;
    private LayeredNetworkStruct lsn;
    private FeatureSet usedFeatures;
    /** Per-layer warm-start marginals from the last initFromMarginal call. */
    private List<Matrix> initMarginalBlocks;

    /**
     * Creates a new layered queueing network with the specified name.
     * 
     * @param name the name of the layered network
     */
    public LayeredNetwork(String name) {
        super(name);
//        this.aux = new Aux();
//        this.lqnGraph = new JLineMatrix(0,0,0);
//        this.taskGraph = new JLineMatrix(0,0,0);
        this.ensemble = new ArrayList<>();
        this.hosts = new HashMap<>();
        this.activities = new HashMap<>();
        this.tasks = new HashMap<>();
        this.reftasks = new HashMap<>();
        this.entries = new HashMap<>();
        this.nodes = new HashMap<>();
        this.lsn = null;
        this.param = new Param();
        this.param.Nodes.RespT = 0;
        this.param.Nodes.Tput = 0;
        this.param.Nodes.Util = 0;
        this.param.Edges.RespT = 0;
        this.param.Edges.Tput = 0;
    }

    public Map<Integer, Host> getHosts() {
        return hosts;
    }

    public Map<Integer, Task> getTasks() {
        return tasks;
    }

    public Map<Integer, Entry> getEntries() {
        return entries;
    }

    public Map<Integer, Activity> getActivities() {
        return activities;
    }

    /**
     * Loads a layered queueing network from an XML file.
     * 
     * @param filename the path to the XML file to load
     * @param verbose if true, enables verbose output during loading
     * @return the loaded LayeredNetwork instance
     */
    public static LayeredNetwork load(String filename, boolean verbose) {
        return parseXML(filename, verbose);
    }

    /**
     * Loads a layered queueing network from an XML file with default verbose setting.
     * 
     * @param filename the path to the XML file to load
     * @return the loaded LayeredNetwork instance
     */
    public static LayeredNetwork load(String filename) {
        return load(filename, false);
    }

    /**
     * Parses a layered queueing network from an XML file with default verbose setting.
     * 
     * @param filename the path to the XML file to parse
     * @return the parsed LayeredNetwork instance
     */
    public static LayeredNetwork parseXML(String filename) {
        return parseXML(filename, false);
    }

    /**
     * Parses a layered queueing network from an XML file.
     * 
     * @param filename the path to the XML file to parse
     * @param verbose if true, enables verbose output during parsing
     * @return the parsed LayeredNetwork instance
     */
    public static LayeredNetwork parseXML(String filename, boolean verbose) {

        LayeredNetwork myLN = new LayeredNetwork(filename.replace("_", "\\_"));

        // File validation like MATLAB version
        File file = new File(filename);
        if (!file.exists()) {
            line_error(mfilename(new Object() {}), "File cannot be found. Verify the current directory and the specified filename.");
        }

        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder dBuilder = null;
        try {
            dBuilder = dbFactory.newDocumentBuilder();
        } catch (ParserConfigurationException e) {
            line_error(mfilename(new Object() {}), "XML parser configuration error: " + e.getMessage());
        }

        Document doc = null;
        try {
            // Handle file path resolution like MATLAB version
            if (new File(filename).getParent() == null) {
                // Use which() equivalent or current directory
                doc = dBuilder.parse(new File(System.getProperty("user.dir"), filename));
            } else {
                doc = dBuilder.parse(filename);
            }
        } catch (SAXException e) {
            line_error(mfilename(new Object() {}), "XML parsing error: " + e.getMessage());
        } catch (IOException e) {
            line_error(mfilename(new Object() {}), "File I/O error: " + e.getMessage());
        }

        doc.getDocumentElement().normalize();
        validateInputModel(doc);

        if (verbose) {
            System.out.println("Parsing LQN file" + filename);
            System.out.println("Root element:" + doc.getDocumentElement().getNodeName());
        }

        Map<Integer, String> hosts = new HashMap<>(); //list of hosts - Proc
        Map<Integer, Map<String, List<Integer>>> tasks = new HashMap<>(); //list of tasks - Task, ProcID
        Map<Integer, Map<String, List<Integer>>> entries = new HashMap<>(); //list of entries - Entry, TaskID, ProcID
        Map<Integer, Map<String, List<Integer>>> activities = new HashMap<>(); //list of activities - Act, TaskID, ProcID
        int procID = 1;
        int taskID = 1;
        int entryID = 1;
        int actID = 1;
        Map<Integer, Processor> procObj = new HashMap<>();
        Map<Integer, Task> taskObj = new HashMap<>();
        Map<Integer, Entry> entryObj = new HashMap<>();
        Map<Integer, Activity> actObj = new HashMap<>();

        NodeList procList = doc.getElementsByTagName("processor");

        for (int i = 0; i < procList.getLength(); i++) {
            Element procElement = (Element) procList.item(i);
            String name = procElement.getAttribute("name");
            String scheduling = procElement.getAttribute("scheduling");
            // LQN schema default for a processor with no scheduling attribute is FCFS.
            if (scheduling == null || scheduling.isEmpty()) {
                scheduling = "fcfs";
            }

            String multiplicityString = procElement.getAttribute("multiplicity");
            double multiplicity = multiplicityString.isEmpty() ? 1.0 : Double.parseDouble(multiplicityString);

            String replicationString = procElement.getAttribute("replication");
            double replication = replicationString.isEmpty() ? 1.0 : Double.parseDouble(replicationString);

            if (scheduling.equals("inf")) {
                // Override finite multiplicity to infinity for INF scheduling processors
                multiplicity = Inf;
            } else if (Double.isNaN(multiplicity)) {
                multiplicity = 1;
            }

            String quantumString = procElement.getAttribute("quantum");
            double quantum = quantumString.isEmpty() ? 0.001 : Double.parseDouble(quantumString);

            String speedFactorString = procElement.getAttribute("speed-factor");
            double speedFactor = speedFactorString.isEmpty() ? 1.0 : Double.parseDouble(speedFactorString);

            Processor newProc = new Processor(myLN, name, (int) multiplicity, SchedStrategy.fromText(scheduling), quantum, speedFactor);
            newProc.setReplication((int) replication);
            procObj.put(procObj.size(), newProc);

            NodeList taskList = procElement.getElementsByTagName("task");

            for (int j = 0; j < taskList.getLength(); j++) {
                Element taskElement = (Element) taskList.item(j);
                String tName = taskElement.getAttribute("name");
                String tScheduling = taskElement.getAttribute("scheduling");
                // LQN schema default for a task with no scheduling attribute is FCFS.
                if (tScheduling == null || tScheduling.isEmpty()) {
                    tScheduling = "fcfs";
                }

                String treplicationString = taskElement.getAttribute("replication");
                double tReplication = treplicationString.isEmpty() ? 1.0 : Double.parseDouble(treplicationString);

                String tmultiplicityString = taskElement.getAttribute("multiplicity");
                double tMultiplicity = tmultiplicityString.isEmpty() ? 1.0 : Double.parseDouble(tmultiplicityString);

                // Override finite multiplicity to infinity for INF scheduling tasks (matches MATLAB behavior)
                if ("inf".equalsIgnoreCase(tScheduling)) {
                    if (Double.isFinite(tMultiplicity)) {
                        line_warning(mfilename(new Object(){}), "A finite multiplicity is specified for a task with inf scheduling. Remove it or set it to inf.");
                    }
                    tMultiplicity = Inf;
                } else if (Double.isNaN(tMultiplicity)) {
                    tMultiplicity = 1;
                }

                String tthinkTimeMeanString = taskElement.getAttribute("think-time");
                double tThinkTimeMean = tthinkTimeMeanString.isEmpty() ? 0.0 : Double.parseDouble(tthinkTimeMeanString);

                Distribution thinkTime;
                if (tThinkTimeMean <= 0.0) {
                    thinkTime = Immediate.getInstance();
                } else {
                    thinkTime = Exp.fitMean(tThinkTimeMean);
                }
                /* LINE .lqnx dialect: the presence of <cache> is what makes this a CacheTask.
                 * itemLevelCap is read back as the ARRAY it is, one entry per <level>. */
                Element cacheElement = firstDirectChild(taskElement, "cache");
                Task newTask;
                if (cacheElement != null) {
                    String itemsStr = cacheElement.getAttribute("items");
                    int cacheItems = itemsStr.isEmpty() ? 1 : Integer.parseInt(itemsStr);
                    String replStr = cacheElement.getAttribute("replacement");
                    ReplacementStrategy repl = replStr.isEmpty()
                            ? ReplacementStrategy.FIFO : ReplacementStrategy.fromText(replStr);
                    List<Integer> caps = new ArrayList<Integer>();
                    NodeList levelList = cacheElement.getElementsByTagName("level");
                    for (int lv = 0; lv < levelList.getLength(); lv++) {
                        String capStr = ((Element) levelList.item(lv)).getAttribute("capacity");
                        caps.add(capStr.isEmpty() ? 1 : Integer.parseInt(capStr));
                    }
                    if (caps.isEmpty()) {
                        caps.add(1);
                    }
                    int[] levelCaps = new int[caps.size()];
                    for (int lv = 0; lv < caps.size(); lv++) {
                        levelCaps[lv] = caps.get(lv);
                    }
                    CacheTask newCacheTask = new CacheTask(myLN, tName, cacheItems, levelCaps, repl,
                            (int) tMultiplicity, SchedStrategy.fromText(tScheduling));
                    if (Boolean.parseBoolean(cacheElement.getAttribute("retrieval"))) {
                        newCacheTask.setRetrieval(true);
                    }
                    newCacheTask.setThinkTime(thinkTime);
                    newTask = newCacheTask;
                } else {
                    newTask = new Task(myLN, tName, (int) tMultiplicity, SchedStrategy.fromText(tScheduling), thinkTime);
                }
                newTask.setReplication((int) replication);

                /* LINE .lqnx dialect: setup and delay-off times. Reconstructed through the
                 * model's own setters so a mean plus an SCV rebuilds the family the setter
                 * would have built (Exp at SCV 1, an APH fit otherwise). */
                Element setupElement = firstDirectChild(taskElement, "setup");
                if (setupElement != null) {
                    newTask.setSetupTime(distFromMeanAndSCV(setupElement));
                }
                Element delayOffElement = firstDirectChild(taskElement, "delay-off");
                if (delayOffElement != null) {
                    newTask.setDelayOffTime(distFromMeanAndSCV(delayOffElement));
                }

                // Parse priority attribute if present
                String tPriorityString = taskElement.getAttribute("priority");
                if (!tPriorityString.isEmpty()) {
                    newTask.setPriority(Integer.parseInt(tPriorityString));
                }

                // Parse fan-in element if present (used for replication load distribution)
                NodeList fanInList = taskElement.getElementsByTagName("fan-in");
                if (fanInList.getLength() > 0) {
                    Element fanInElement = (Element) fanInList.item(0);
                    String source = fanInElement.getAttribute("source");
                    String value = fanInElement.getAttribute("value");
                    if (!source.isEmpty() && !value.isEmpty()) {
                        newTask.setFanIn(source, Integer.parseInt(value));
                    }
                }

                // Parse fan-out elements if present (used for replication load distribution)
                NodeList fanOutList = taskElement.getElementsByTagName("fan-out");
                for (int fo = 0; fo < fanOutList.getLength(); fo++) {
                    Element fanOutElement = (Element) fanOutList.item(fo);
                    String dest = fanOutElement.getAttribute("dest");
                    String value = fanOutElement.getAttribute("value");
                    if (!dest.isEmpty() && !value.isEmpty()) {
                        newTask.setFanOut(dest, Integer.parseInt(value));
                    }
                }

                newTask.on(newProc);  // Assign task to its processor
                taskObj.put(taskObj.size(), newTask);

                NodeList entryList = taskElement.getElementsByTagName("entry");
                for (int k = 0; k < entryList.getLength(); k++) {
                    Element entryElement = (Element) entryList.item(k);
                    String eName = entryElement.getAttribute("name");
                    /* LINE .lqnx dialect: the presence of <item-entry> is what makes this an
                     * ItemEntry. Its popularity is rebuilt from the flat parameter list, split
                     * on the declared cardinality. */
                    Element itemElement = firstDirectChild(entryElement, "item-entry");
                    Entry newEntry;
                    if (itemElement != null) {
                        String cardStr = itemElement.getAttribute("cardinality");
                        int cardinality = cardStr.isEmpty() ? 1 : Integer.parseInt(cardStr);
                        Distribution popularity = readAccessPopularity(itemElement, cardinality);
                        if (popularity == null) {
                            // An ItemEntry with no popularity is still an ItemEntry; give it the
                            // uniform law over its items rather than degrading it to a plain Entry.
                            Matrix uniform = new Matrix(1, cardinality);
                            for (int q = 0; q < cardinality; q++) {
                                uniform.set(0, q, 1.0 / cardinality);
                            }
                            popularity = new DiscreteSampler(uniform);
                        }
                        newEntry = new ItemEntry(myLN, eName, cardinality, popularity);
                    } else {
                        newEntry = new Entry(myLN, eName);
                    }

                    // Parse entry type attribute if present
                    String eType = entryElement.getAttribute("type");
                    if (!eType.isEmpty()) {
                        newEntry.setType(eType);
                    }

                    String openArrivalRateMeanString = entryElement.getAttribute("open-arrival-rate");
                    if (!openArrivalRateMeanString.isEmpty()) {
                        double openArrivalRate = Double.parseDouble(openArrivalRateMeanString);
                        newEntry.setArrival(Exp.fitMean(1.0 / openArrivalRate));
                    }

                    entryObj.put(entryObj.size(), newEntry);

                    // Assign entry to its parent task
                    newTask.addEntry(newEntry);
                    newEntry.parent = newTask;

                    // Parse forwarding calls
                    NodeList forwardingList = entryElement.getElementsByTagName("forwarding");
                    for (int fw = 0; fw < forwardingList.getLength(); fw++) {
                        Element fwdElement = (Element) forwardingList.item(fw);
                        String destName = fwdElement.getAttribute("dest");
                        double prob = fwdElement.hasAttribute("prob") ?
                            Double.parseDouble(fwdElement.getAttribute("prob")) : 1.0;
                        newEntry.forward(destName, prob);
                    }

                    NodeList entryPhaseActsList = entryElement.getElementsByTagName("entry-phase-activities");
                    if (entryPhaseActsList.getLength() > 0) {
                        Element entryPhaseActsElement = (Element) entryPhaseActsList.item(0);
                        NodeList actList = entryPhaseActsElement.getElementsByTagName("activity");
                        Map<Integer, String> nameList = new HashMap<>();
                        // the lowest phase present is bound to the entry; an entry
                        // declaring only phase 2 would otherwise have no bound
                        // activity and its host demand would never be routed
                        int minPhase = Integer.MAX_VALUE;
                        for (int l = 0; l < actList.getLength(); l++) {
                            int ph = (int) Double.parseDouble(((Element) actList.item(l)).getAttribute("phase"));
                            if (ph < minPhase) {
                                minPhase = ph;
                            }
                        }
                        for (int l = 0; l < actList.getLength(); l++) {
                            Element actElement = (Element) actList.item(l);
                            double phase = Double.parseDouble(actElement.getAttribute("phase"));
                            nameList.put((int) phase, actElement.getAttribute("name"));

                            String hostDemandMeanString = actElement.getAttribute("host-demand-mean");
                            double hostDemandMean = hostDemandMeanString.isEmpty() ? 0.0 : Double.parseDouble(hostDemandMeanString);

                            String hostDemandSCVString = actElement.getAttribute("host-demand-cvsq");
                            double hostDemandSCV = hostDemandSCVString.isEmpty() ? 1.0 : Double.parseDouble(hostDemandSCVString);

                            Distribution hostDemand = Immediate.getInstance();
                            if (hostDemandMean > 0) {
                                if (hostDemandSCV <= 0) {
                                    hostDemand = new Det(hostDemandMean);
                                } else if (hostDemandSCV == 1) {
                                    hostDemand = Exp.fitMean(hostDemandMean);
                                } else {
                                    hostDemand = APH.fitMeanAndSCV(hostDemandMean, hostDemandSCV);
                                }
                            }
                            String boundToEntry;
                            if ((int) phase == minPhase) {
                                boundToEntry = newEntry.getName();
                            } else {
                                boundToEntry = "";
                            }

                            String callOrder = actElement.getAttribute("call-order");
                            Activity newAct = new Activity(myLN, nameList.get((int) phase), hostDemand, boundToEntry, callOrder);
                            newAct.setPhase((int) phase);  // Store phase number from XML

                            // Parse activity think-time
                            String actThinkTimeMeanString = actElement.getAttribute("think-time");
                            if (!actThinkTimeMeanString.isEmpty()) {
                                double actThinkTimeMean = Double.parseDouble(actThinkTimeMeanString);
                                if (actThinkTimeMean > 0.0) {
                                    newAct.setThinkTime(actThinkTimeMean);
                                }
                            }

                            actObj.put(actObj.size(), newAct);

                            NodeList synchCalls = actElement.getElementsByTagName("synch-call");
                            for (int m = 0; m < synchCalls.getLength(); m++) {
                                Element callElement = (Element) synchCalls.item(m);
                                String dest = callElement.getAttribute("dest");
                                double mean = Double.parseDouble(callElement.getAttribute("calls-mean"));
                                newAct.synchCall(dest, mean);
                            }

                            NodeList asynchCalls = actElement.getElementsByTagName("asynch-call");
                            for (int m = 0; m < asynchCalls.getLength(); m++) {
                                Element callElement = (Element) asynchCalls.item(m);
                                String dest = callElement.getAttribute("dest");
                                double mean = Double.parseDouble(callElement.getAttribute("calls-mean"));
                                newAct.asynchCall(dest, mean);
                            }

                            parseCallGroups(actElement, newAct);

                            Map<String, List<Integer>> tempMap = new HashMap<>();
                            List<Integer> tempList = new ArrayList<>();
                            tempList.add(taskID);
                            tempList.add(procID);
                            tempMap.put(newAct.getName(), tempList);
                            activities.put(activities.size(), tempMap);
                            newTask.addActivity(newAct);
                            newAct.setParent(newTask);
                            actID++;
                        }


                        List<Integer> phases = new ArrayList<>(nameList.keySet());
                        Collections.sort(phases);
                        for (int l = 0; l < phases.size() - 1; l++) {
                            ActivityPrecedence newPrec = new ActivityPrecedence(Collections.singletonList(nameList.get(phases.get(l))), Collections.singletonList(nameList.get(phases.get(l + 1))));
                            newTask.addPrecedence(newPrec);
                        }

//                        if (!nameList.isEmpty()) {
//                            newEntry.replyActivity.put(1, nameList.get(1));
//                        }

                        Map<String, List<Integer>> tempMap = new HashMap<>();
                        List<Integer> tempList = new ArrayList<>();
                        tempList.add(taskID);
                        tempList.add(procID);
                        tempMap.put(newTask.getName(), tempList);
                        entries.put(entries.size(), tempMap);
                        // Note: newTask.addEntry(newEntry) already called at line 300
                        // Note: newEntry.parent = newTask already set at line 301
                        entryID++;
                    }
                }  // End of entry loop

                    NodeList taskActsList = taskElement.getElementsByTagName("task-activities");
                    if (taskActsList.getLength() > 0) {
                        Element taskActsElement = (Element) taskActsList.item(0);
                        NodeList actList = taskActsElement.getElementsByTagName("activity");
                        for (int l = 0; l < actList.getLength(); l++) {
                            Element actElement = (Element) actList.item(l);
                            if (actElement.getParentNode().getNodeName().equals("task-activities")) {
                                String actName = actElement.getAttribute("name");

                                String hostDemandMeanString = actElement.getAttribute("host-demand-mean");
                                double hostDemandMean = hostDemandMeanString.isEmpty() ? 0.0 : Double.parseDouble(hostDemandMeanString);

                                String hostDemandSCVString = actElement.getAttribute("host-demand-cvsq");
                                double hostDemandSCV = hostDemandSCVString.isEmpty() ? 1.0 : Double.parseDouble(hostDemandSCVString);

                                Distribution hostDemand = Immediate.getInstance();
                                if (Double.isNaN(hostDemandSCV)) {
                                    hostDemandSCV = 1.0;
                                }
                                if (hostDemandMean <= 0.0) {
                                    hostDemand = Immediate.getInstance();
                                } else {
                                    if (hostDemandSCV <= 0.0) {
                                        hostDemand = new Det(hostDemandMean);
                                    } else if (hostDemandSCV < 1.0) {
                                        hostDemand = APH.fitMeanAndSCV(hostDemandMean, hostDemandSCV);
                                    } else if (hostDemandSCV == 1.0) {
                                        hostDemand = Exp.fitMean(hostDemandMean);
                                    } else {
                                        hostDemand = HyperExp.fitMeanAndSCV(hostDemandMean, hostDemandSCV);
                                    }
                                }
                                String boundToEntry = actElement.getAttribute("bound-to-entry");
                                String callOrder = actElement.getAttribute("call-order");
                                Activity newAct = new Activity(myLN, actName, hostDemand, boundToEntry, callOrder);

                                // Parse activity think-time
                                String actThinkTimeMeanString = actElement.getAttribute("think-time");
                                if (!actThinkTimeMeanString.isEmpty()) {
                                    double actThinkTimeMean = Double.parseDouble(actThinkTimeMeanString);
                                    if (actThinkTimeMean > 0.0) {
                                        newAct.setThinkTime(actThinkTimeMean);
                                    }
                                }

                                actObj.put(actObj.size(), newAct);

                                NodeList synchCalls = actElement.getElementsByTagName("synch-call");
                                for (int m = 0; m < synchCalls.getLength(); m++) {
                                    Element callElement = (Element) synchCalls.item(m);
                                    String dest = callElement.getAttribute("dest");
                                    double mean = Double.parseDouble(callElement.getAttribute("calls-mean"));
                                    newAct.synchCall(dest, mean);
                                }

                                NodeList asynchCalls = actElement.getElementsByTagName("asynch-call");
                                for (int m = 0; m < asynchCalls.getLength(); m++) {
                                    Element callElement = (Element) asynchCalls.item(m);
                                    String dest = callElement.getAttribute("dest");
                                    double mean = Double.parseDouble(callElement.getAttribute("calls-mean"));
                                    newAct.asynchCall(dest, mean);
                                }

                                parseCallGroups(actElement, newAct);

                                Map<String, List<Integer>> tempMap = new HashMap<>();
                                List<Integer> tempList = new ArrayList<>();
                                tempList.add(taskID);
                                tempList.add(procID);
                                tempMap.put(newAct.getName(), tempList);
                                activities.put(activities.size(), tempMap);
                                newTask.addActivity(newAct);
                                newAct.setParent(newTask);
                                actID++;
                            }
                        }

                        NodeList precList = taskActsElement.getElementsByTagName("precedence");
                        for (int l = 0; l < precList.getLength(); l++) {
                            Element precElement = (Element) precList.item(l);

                            // Pre-activity parsing
                            String[] preTypes = {ActivityPrecedenceType.PRE_SEQ, ActivityPrecedenceType.PRE_AND, ActivityPrecedenceType.PRE_OR};
                            NodeList preList = null;
                            String preType = "";
                            for (String type : preTypes) {
                                preType = type;
                                preList = precElement.getElementsByTagName(preType);
                                if (preList.getLength() > 0) {
                                    break;
                                }
                            }

                            Element preElement = (Element) preList.item(0);
                            Matrix preParams = null;
                            NodeList preActList = preElement.getElementsByTagName("activity");
                            List<String> preActs = new ArrayList<>();
                            
                            if (preType.equals(ActivityPrecedenceType.PRE_OR)) {
                                // see _kb/04-networkstruct.md (LayeredNetwork.getStruct() graph/validation rules) for rationale
                                for (int m = 0; m < preActList.getLength(); m++) {
                                    Element preActElement = (Element) preActList.item(m);
                                    preActs.add(preActElement.getAttribute("name"));
                                }
                            } else if (preType.equals(ActivityPrecedenceType.PRE_AND)) {
                                for (int m = 0; m < preActList.getLength(); m++) {
                                    Element preActElement = (Element) preActList.item(m);
                                    preActs.add(preActElement.getAttribute("name"));
                                }
                                String quorumStr = preElement.getAttribute("quorum");
                                if (!quorumStr.isEmpty()) {
                                    preParams = Matrix.singleton(Double.parseDouble(quorumStr));
                                }
                            } else {
                                Element preActElement = (Element) preActList.item(0);
                                preActs.add(preActElement.getAttribute("name"));
                            }

                            // Post-activity parsing. The post side is minOccurs="0" in
                            // lqn-core.xsd: a precedence carrying only a pre element declares a
                            // TERMINAL activity and no successor, so it contributes no edge.
                            String[] postTypes = {ActivityPrecedenceType.POST_SEQ, ActivityPrecedenceType.POST_AND, ActivityPrecedenceType.POST_OR, ActivityPrecedenceType.POST_LOOP, ActivityPrecedenceType.POST_CACHE};
                            NodeList postList = null;
                            String postType = null;
                            boolean hasPost = false;
                            for (String type : postTypes) {
                                postType = type;
                                postList = precElement.getElementsByTagName(postType);
                                if (postList.getLength() > 0) {
                                    hasPost = true;
                                    break;
                                }
                            }
                            if (!hasPost) {
                                continue;
                            }

                            Element postElement = (Element) postList.item(0);
                            NodeList postActList = postElement.getElementsByTagName("activity");
                            List<String> postActs = new ArrayList<>();
                            Matrix postParams = null;
                            
                            if (postType.equals(ActivityPrecedenceType.POST_OR)) {
                                postParams = new Matrix(1, postActList.getLength());
                                for (int m = 0; m < postActList.getLength(); m++) {
                                    Element postActElement = (Element) postActList.item(m);
                                    postActs.add(postActElement.getAttribute("name"));
                                    postParams.set(0, m, Double.parseDouble(postActElement.getAttribute("prob")));
                                }
                            } else if (postType.equals(ActivityPrecedenceType.POST_LOOP)) {
                                postParams = new Matrix(1, postActList.getLength());
                                for (int m = 0; m < postActList.getLength(); m++) {
                                    Element postActElement = (Element) postActList.item(m);
                                    postActs.add(postActElement.getAttribute("name"));
                                    postParams.set(0, m, Double.parseDouble(postActElement.getAttribute("count")));
                                }
                                postActs.add(postElement.getAttribute("end"));
                            } else if (postType.equals(ActivityPrecedenceType.POST_CACHE)) {
                                // cache-result is explicit; a file written without it (or by
                                // another codebase) still loads, on document order.
                                String hitAct = null, missAct = null;
                                List<String> unlabelled = new ArrayList<>();
                                for (int m = 0; m < postActList.getLength(); m++) {
                                    Element postActElement = (Element) postActList.item(m);
                                    String actName = postActElement.getAttribute("name");
                                    String result = postActElement.getAttribute("cache-result");
                                    if ("hit".equalsIgnoreCase(result)) {
                                        hitAct = actName;
                                    } else if ("miss".equalsIgnoreCase(result)) {
                                        missAct = actName;
                                    } else {
                                        unlabelled.add(actName);
                                    }
                                }
                                for (String actName : unlabelled) {
                                    if (hitAct == null) {
                                        hitAct = actName;
                                    } else if (missAct == null) {
                                        missAct = actName;
                                    } else {
                                        postActs.add(actName);
                                    }
                                }
                                if (hitAct != null) {
                                    postActs.add(0, hitAct);
                                }
                                if (missAct != null) {
                                    postActs.add(Math.min(1, postActs.size()), missAct);
                                }
                            } else {
                                for (int m = 0; m < postActList.getLength(); m++) {
                                    Element postActElement = (Element) postActList.item(m);
                                    postActs.add(postActElement.getAttribute("name"));
                                }
                            }
                            
                            ActivityPrecedence newPrec = new ActivityPrecedence(preActs, postActs, preType, postType, preParams, postParams);
                            newTask.addPrecedence(newPrec);
                        }
                        NodeList replyList = taskActsElement.getElementsByTagName("reply-entry");
                        for (int l = 0; l < replyList.getLength(); l++) {
                            Element replyElement = (Element) replyList.item(l);
                            String replyName = replyElement.getAttribute("name");
                            
                            // Find the entry object by name
                            Entry replyEntry = null;
                            for (Entry entry : entryObj.values()) {
                                if (entry.getName().equals(replyName)) {
                                    replyEntry = entry;
                                    break;
                                }
                            }
                            
                            if (replyEntry != null) {
                                NodeList replyActList = replyElement.getElementsByTagName("reply-activity");
                                for (int m = 0; m < replyActList.getLength(); m++) {
                                    Element replyActElement = (Element) replyActList.item(m);
                                    String replyActName = replyActElement.getAttribute("name");
                                    replyEntry.replyActivity.put(replyEntry.replyActivity.size() + 1, replyActName);
                                }
                            }
                        }
                    }
            }
        }
        return myLN;
    }

    /**
     * Rejects a structurally inconsistent LQN document.
     *
     * Run on the parsed document before any object is built, so that a defective
     * input is named at its source instead of surfacing as a downstream failure.
     * The same checks, in the same order and with the same messages, are applied
     * by the MATLAB, Python and C++ readers.
     *
     * @param doc the parsed LQN document
     */
    private static void validateInputModel(Document doc) {
        final double tol = 1e-6;
        List<String> procNames = new ArrayList<String>();
        List<String> taskNames = new ArrayList<String>();
        List<String> entryNames = new ArrayList<String>();
        List<String> entryOwner = new ArrayList<String>(); // task owning entryNames.get(k)
        List<Boolean> isRefEntry = new ArrayList<Boolean>();
        List<String> callDests = new ArrayList<String>();
        List<String> replyEntries = new ArrayList<String>();
        boolean hasRefTask = false;
        boolean hasOpenArrival = false;

        NodeList procList = doc.getElementsByTagName("processor");
        for (int i = 0; i < procList.getLength(); i++) {
            Element procElement = (Element) procList.item(i);
            String procName = procElement.getAttribute("name");
            if (procNames.contains(procName)) {
                line_error(mfilename(new Object() {}), String.format("Duplicate processor name \"%s\".", procName));
            }
            procNames.add(procName);

            NodeList taskList = procElement.getElementsByTagName("task");
            for (int j = 0; j < taskList.getLength(); j++) {
                Element taskElement = (Element) taskList.item(j);
                String taskName = taskElement.getAttribute("name");
                if (taskNames.contains(taskName)) {
                    line_error(mfilename(new Object() {}), String.format("Duplicate task name \"%s\".", taskName));
                }
                taskNames.add(taskName);
                boolean isRef = "ref".equalsIgnoreCase(taskElement.getAttribute("scheduling"));
                hasRefTask = hasRefTask || isRef;

                NodeList entryList = taskElement.getElementsByTagName("entry");
                if (entryList.getLength() == 0) {
                    line_error(mfilename(new Object() {}), String.format("Task \"%s\" has no entries.", taskName));
                }
                for (int k = 0; k < entryList.getLength(); k++) {
                    Element entryElement = (Element) entryList.item(k);
                    String entryName = entryElement.getAttribute("name");
                    if (entryNames.contains(entryName)) {
                        line_error(mfilename(new Object() {}), String.format("Duplicate entry name \"%s\".", entryName));
                    }
                    entryNames.add(entryName);
                    entryOwner.add(taskName);
                    isRefEntry.add(isRef);

                    String openArrivalRateString = entryElement.getAttribute("open-arrival-rate");
                    if (!openArrivalRateString.isEmpty()) {
                        double openArrivalRate = parseProb(openArrivalRateString, Double.NaN);
                        if (openArrivalRate > 0) {
                            hasOpenArrival = true;
                            if (isRef) {
                                line_error(mfilename(new Object() {}), String.format("Entry \"%s\" belongs to reference task \"%s\" and cannot have open arrivals.", entryName, taskName));
                            }
                        }
                    }

                    NodeList fwdList = entryElement.getElementsByTagName("forwarding");
                    if (isRef && fwdList.getLength() > 0) {
                        line_error(mfilename(new Object() {}), String.format("Entry \"%s\" belongs to reference task \"%s\" and cannot forward requests.", entryName, taskName));
                    }
                    double fwdTotal = 0.0;
                    for (int fw = 0; fw < fwdList.getLength(); fw++) {
                        Element fwdElement = (Element) fwdList.item(fw);
                        double prob = parseProb(fwdElement.getAttribute("prob"), 1.0);
                        if (Double.isNaN(prob) || prob < 0.0 || prob > 1.0) {
                            line_error(mfilename(new Object() {}), String.format("Forwarding from entry \"%s\" to entry \"%s\" has an invalid probability of %s.", entryName, fwdElement.getAttribute("dest"), fmtNum(prob)));
                        }
                        fwdTotal += prob;
                    }
                    if (fwdTotal > 1.0 + tol) {
                        line_error(mfilename(new Object() {}), String.format("Entry \"%s\" has a total forwarding probability of %s.", entryName, fmtNum(fwdTotal)));
                    }
                }

                // activity names are unique within their task; a name under a pre or post list is a reference, not a declaration
                List<String> actNames = new ArrayList<String>();
                NodeList actList = taskElement.getElementsByTagName("activity");
                for (int l = 0; l < actList.getLength(); l++) {
                    Element actElement = (Element) actList.item(l);
                    String parentTag = actElement.getParentNode().getNodeName();
                    if (!"task-activities".equals(parentTag) && !"entry-phase-activities".equals(parentTag)) {
                        continue;
                    }
                    String actName = actElement.getAttribute("name");
                    if (actNames.contains(actName)) {
                        line_error(mfilename(new Object() {}), String.format("Duplicate activity name \"%s\" in task \"%s\".", actName, taskName));
                    }
                    actNames.add(actName);
                }

                NodeList synchCalls = taskElement.getElementsByTagName("synch-call");
                for (int m = 0; m < synchCalls.getLength(); m++) {
                    callDests.add(((Element) synchCalls.item(m)).getAttribute("dest"));
                }
                NodeList asynchCalls = taskElement.getElementsByTagName("asynch-call");
                for (int m = 0; m < asynchCalls.getLength(); m++) {
                    callDests.add(((Element) asynchCalls.item(m)).getAttribute("dest"));
                }
                NodeList taskFwdList = taskElement.getElementsByTagName("forwarding");
                for (int fw = 0; fw < taskFwdList.getLength(); fw++) {
                    callDests.add(((Element) taskFwdList.item(fw)).getAttribute("dest"));
                }

                NodeList orList = taskElement.getElementsByTagName("post-OR");
                for (int l = 0; l < orList.getLength(); l++) {
                    NodeList branchList = ((Element) orList.item(l)).getElementsByTagName("activity");
                    double branchTotal = 0.0;
                    for (int m = 0; m < branchList.getLength(); m++) {
                        Element branchElement = (Element) branchList.item(m);
                        double prob = parseProb(branchElement.getAttribute("prob"), 1.0);
                        if (Double.isNaN(prob) || prob < 0.0 || prob > 1.0) {
                            line_error(mfilename(new Object() {}), String.format("Activity \"%s\" in task \"%s\" has an invalid branch probability of %s.", branchElement.getAttribute("name"), taskName, fmtNum(prob)));
                        }
                        branchTotal += prob;
                    }
                    if (Math.abs(branchTotal - 1.0) > tol) {
                        line_error(mfilename(new Object() {}), String.format("Branch probabilities of an OR-fork in task \"%s\" sum to %s instead of 1.", taskName, fmtNum(branchTotal)));
                    }
                }

                NodeList replyList = taskElement.getElementsByTagName("reply-entry");
                for (int l = 0; l < replyList.getLength(); l++) {
                    replyEntries.add(((Element) replyList.item(l)).getAttribute("name"));
                }
            }
        }

        for (int c = 0; c < callDests.size(); c++) {
            int idx = entryNames.indexOf(callDests.get(c));
            if (idx >= 0 && isRefEntry.get(idx)) {
                line_error(mfilename(new Object() {}), String.format("Entry \"%s\" belongs to reference task \"%s\" and cannot receive requests.", entryNames.get(idx), entryOwner.get(idx)));
            }
        }

        for (int r = 0; r < replyEntries.size(); r++) {
            int idx = entryNames.indexOf(replyEntries.get(r));
            if (idx >= 0 && isRefEntry.get(idx)) {
                line_error(mfilename(new Object() {}), String.format("Entry \"%s\" belongs to reference task \"%s\" and cannot be replied to.", entryNames.get(idx), entryOwner.get(idx)));
            }
        }

        if (!hasRefTask && !hasOpenArrival) {
            line_error(mfilename(new Object() {}), "The model has no reference task and no open arrivals.");
        }
    }

    /**
     * Reads a numeric attribute of the input document, as the other readers do.
     *
     * @param s     the attribute value, empty when the attribute is absent
     * @param dflt  the value an absent attribute stands for
     * @return the parsed value, or NaN when the text is not a number
     */
    private static double parseProb(String s, double dflt) {
        if (s == null || s.isEmpty()) {
            return dflt;
        }
        try {
            return Double.parseDouble(s);
        } catch (NumberFormatException e) {
            return Double.NaN;
        }
    }

    /**
     * Formats a number the way the MATLAB, Python and C++ readers do, so that the
     * validation messages agree across the codebases.
     *
     * @param v the value to format
     * @return the shortest general-format rendering of v
     */
    private static String fmtNum(double v) {
        String s = String.format(java.util.Locale.US, "%g", v);
        if (s.indexOf('.') >= 0 && s.indexOf('e') < 0 && s.indexOf('E') < 0) {
            s = s.replaceAll("0+$", "");
            if (s.endsWith(".")) {
                s = s.substring(0, s.length() - 1);
            }
        }
        return s;
    }

    /**
     * Reads a layered queueing network from an XML file with default verbose setting.
     *
     * @param filename the path to the XML file to read
     * @return the read LayeredNetwork instance
     */
    public static LayeredNetwork readXML(String filename) {
        return parseXML(filename, false);
    }

    /**
     * Reads a layered queueing network from an XML file.
     * 
     * @param filename the path to the XML file to read
     * @param verbose if true, enables verbose output during reading
     * @return the read LayeredNetwork instance
     */
    public static LayeredNetwork readXML(String filename, boolean verbose) {
        return parseXML(filename, verbose);
    }

    /**
     * Views a layered queueing network model using the default JLQN path.
     * 
     * @param filename the path to the model file to view
     */
    public static void viewModel(String filename) {
        viewModel(jlqnGetPath(), filename);
    }

    /**
     * Views a layered queueing network model using the specified JLQN path.
     * 
     * @param jlqnPath the path to the JLQN executable
     * @param filename the path to the model file to view
     */
    public static void viewModel(String jlqnPath, String filename) {
        Path path = Paths.get(filename).getParent();
        if (path == null) {
            filename = Paths.get(java.lang.System.getProperty("user.dir"), filename).toString();
        }

        // No shell is involved (SysUtils.system spawns the argv directly), so a
        // redirection appended here would reach JLQN as extra arguments. The
        // launcher is resolved rather than assumed to be on the PATH.
        String cmd = String.format(
                "\"%s\" -cp \"%s\" jlqn.commandline.Jlqn \"%s\"",
                SysUtils.javaLauncher(), jlqnPath, filename
        );

        java.lang.System.out.println("JLQN view model command: " + cmd);
        String output = SysUtils.system(cmd);
        java.lang.System.out.println("JLQN view model command output: " + output);
    }

    private static void writeJLQNActivity(Document document, Element activities, String name, String task, String hostDemandMean) {
        Element activity = document.createElement("activity");
        activity.setAttribute("name", name);
        activity.setAttribute("task", task);
        activity.setAttribute("host_demand_mean", hostDemandMean);
        activities.appendChild(activity);
    }

    private static void writeJLQNCall(Document document, Element calls, String name, String activity, String entry, CallType type, double callMean) {
        Element call = document.createElement("call");
        call.setAttribute("activity", activity);
        call.setAttribute("entry", entry);
        call.setAttribute("mean_repeat", "" + callMean);
        call.setAttribute("name", name);
        switch (type) {
            case ASYNC:
                call.setAttribute("type", "Asynchronous");
                break;
            case SYNC:
                call.setAttribute("type", "Synchronous");
                break;
            case FWD:
                call.setAttribute("type", "Forwarding");
                break;
            default:
                // no-op
        }
        calls.appendChild(call);
    }

    private static void writeJLQNEntry(Document document, Element entries, String name, String bndto_activity, String reply_activity, String task,
                                       Map<Integer, String> forwardingDests, Matrix forwardingProbs) {
        Element entry = document.createElement("entry");
        entry.setAttribute("arrival_rate", "0.0");
        entry.setAttribute("bound_to_activity", bndto_activity);
        entry.setAttribute("name", name);
        entry.setAttribute("priority", "0");
        entry.setAttribute("reply_to_activity", reply_activity);
        entry.setAttribute("task", task);

        // Emit per-entry <forwarding> children (matches MATLAB writeXML.m and the
        // existing reader at parseXML/parse path).
        if (forwardingDests != null) {
            for (int fw = 0; fw < forwardingDests.size(); fw++) {
                String destName = forwardingDests.get(fw);
                if (destName == null) {
                    continue;
                }
                Element fwdElement = document.createElement("forwarding");
                fwdElement.setAttribute("dest", destName);
                double prob = (forwardingProbs != null && fw < forwardingProbs.getNumCols())
                    ? forwardingProbs.get(0, fw) : 1.0;
                fwdElement.setAttribute("prob", String.valueOf(prob));
                entry.appendChild(fwdElement);
            }
        }

        entries.appendChild(entry);
    }

    private static void writeJLQNPrecedence(Document document, Element precedence, String activity, String params, String type) {
        Element precedenceActivity = document.createElement("precedence-activity");
        precedenceActivity.setAttribute("activity", activity);
        if (type.compareTo("end") == 0) {
            precedenceActivity.setAttribute("params", "1.0");
        } else {
            precedenceActivity.setAttribute("params", params);
        }
        precedenceActivity.setAttribute("type", type);
        precedence.appendChild(precedenceActivity);
    }

    private static void writeJLQNProcessor(Document document, Element processors, String name, String scheduling, int multiplicity) {
        Element processor = document.createElement("processor");
        if (multiplicity == Integer.MAX_VALUE) {
            processor.setAttribute("multiplicity", "-1");
        } else {
            processor.setAttribute("multiplicity", "" + multiplicity);
        }
        processor.setAttribute("name", name);
        processor.setAttribute("quantum", "0.0");
        processor.setAttribute("replicas", "1");
        processor.setAttribute("scheduling", scheduling);
        processor.setAttribute("speed_factor", "1.0");
        processors.appendChild(processor);
    }

    private static void writeJLQNTask(Document document, Element tasks, Task t) {
        Element task = document.createElement("task");
        task.setAttribute("multiplicity", "" + t.multiplicity);
        task.setAttribute("name", t.getName());
        // The reader parses this attribute back into Task.priority, so writing a
        // constant 0 would silently drop the priority of a prioritised task.
        task.setAttribute("priority", "" + t.getPriority());
        task.setAttribute("processor", t.parent.getName());
        task.setAttribute("replicas", "1");
        task.setAttribute("scheduling", t.scheduling.toString());
        task.setAttribute("think_time_mean", "" + t.thinkTime.getMean());
        tasks.appendChild(task);
    }

    /**
     * Generates the graph representation of the layered network.
     * This method is currently not implemented.
     */
    public void generateGraph() {
    }

    @Override
    public List<Network> getEnsemble() {
        if (this.ensemble.isEmpty()) {
            SolverLN solver = new SolverLN(this);
            this.ensemble = solver.getEnsemble();
        }
        return this.ensemble;
    }

    /**
     * Gets the list of network layers in this layered network.
     * 
     * @return the list of network layers
     */
    public List<Network> getLayers() {
        return getEnsemble();
    }

    // ------------------------------------------------------------------
    // Aggregate flat interface used by SolverENV (mirrors the MATLAB
    // @LayeredNetwork additions and the native-Python LayeredNetwork).
    //
    // An LQN is exposed as the block-diagonal union of its layer networks:
    // aggregate station/class/node counts are the sums over layers, and
    // matrices are laid out block-diagonally in the exact ordering produced by
    // SolverLN.getTranAvg. All these methods share that layout so SolverENV can
    // treat an LN stage through the same interface as a flat Network.
    // ------------------------------------------------------------------

    /** Aggregate (sum over layers) number of stations. */
    public int getNumberOfStations() {
        int m = 0;
        for (Network layer : getEnsemble()) m += layer.getNumberOfStations();
        return m;
    }

    /** Aggregate (sum over layers) number of classes. */
    public int getNumberOfClasses() {
        int k = 0;
        for (Network layer : getEnsemble()) k += layer.getNumberOfClasses();
        return k;
    }

    /** Aggregate (sum over layers) number of nodes. */
    public int getNumberOfNodes() {
        int n = 0;
        for (Network layer : getEnsemble()) n += layer.getNumberOfNodes();
        return n;
    }

    /** Aggregate (sum over layers) number of stateful nodes. */
    public int getNumberOfStatefulNodes() {
        int n = 0;
        for (Network layer : getEnsemble()) n += layer.getNumberOfStatefulNodes();
        return n;
    }

    /** Aggregate list of stations over the layer networks (block order). */
    public List<jline.lang.nodes.Station> getStations() {
        List<jline.lang.nodes.Station> out = new ArrayList<jline.lang.nodes.Station>();
        for (Network layer : getEnsemble()) out.addAll(layer.getStations());
        return out;
    }

    /**
     * Split the aggregate (M x K) marginal queue-length matrix into per-layer
     * blocks and warm-start each layer network. State-continuity mechanism used
     * by SolverENV across environment switches.
     */
    public void initFromMarginal(Matrix n) {
        List<Network> layers = getEnsemble();
        List<Matrix> blocks = new ArrayList<Matrix>();
        int roff = 0, coff = 0;
        for (Network layer : layers) {
            int me = layer.getNumberOfStations();
            int ke = layer.getNumberOfClasses();
            Matrix block = new Matrix(me, ke);
            for (int i = 0; i < me; i++)
                for (int j = 0; j < ke; j++)
                    block.set(i, j, n.get(roff + i, coff + j));
            layer.initFromMarginal(block);
            blocks.add(block);
            roff += me;
            coff += ke;
        }
        this.initMarginalBlocks = blocks;
    }

    /**
     * Per-layer blocks of the warm start supplied by the last
     * {@link #initFromMarginal} call, or null if none was supplied. The layered
     * fixed point hard-resets every layer state when it detects convergence,
     * which discards the warm start, so SolverLN replays these blocks just
     * before running the layer transients; that is what lets a stage of
     * SolverENV resume from the marginal handed over at the environment switch.
     *
     * @return the per-layer warm-start marginals, or null
     */
    public List<Matrix> getInitMarginalBlocks() {
        return this.initMarginalBlocks;
    }

    /**
     * Build the aggregate, block-diagonal {@link NetworkStruct} for this LQN by
     * concatenating the per-layer structs. This lets SolverENV consume an LN
     * stage through its usual NetworkStruct-based code path (sn[e]).
     */
    public NetworkStruct getEnvStruct() {
        List<Network> layers = getEnsemble();
        List<NetworkStruct> ls = new ArrayList<NetworkStruct>();
        int M = 0, K = 0, nch = 0, ncj = 0, nn = 0, nst = 0;
        for (Network layer : layers) {
            NetworkStruct s = layer.getStruct(true);
            ls.add(s);
            M += s.nstations;
            K += s.nclasses;
            nch += s.nchains;
            ncj += s.nclosedjobs;
            nn += s.nnodes;
            nst += s.nstateful;
        }
        NetworkStruct agg = new NetworkStruct();
        agg.nstations = M;
        agg.nclasses = K;
        agg.nchains = nch;
        agg.nclosedjobs = ncj;
        agg.nnodes = nn;
        agg.nstateful = nst;
        agg.stations = new ArrayList<jline.lang.nodes.Station>();
        agg.jobclasses = new ArrayList<jline.lang.JobClass>();
        agg.sched = new HashMap<jline.lang.nodes.Station, SchedStrategy>();
        agg.nservers = new Matrix(M, 1);
        agg.rates = new Matrix(M, K);
        agg.njobs = new Matrix(1, K);
        agg.chains = new Matrix(Math.max(nch, 1), Math.max(K, 1));
        agg.refstat = new Matrix(K, 1);  // refstat is (nclasses x 1), unlike njobs
        int roff = 0, coff = 0, choff = 0;
        for (NetworkStruct s : ls) {
            int me = s.nstations, ke = s.nclasses, che = s.nchains;
            for (int i = 0; i < me; i++) agg.stations.add(s.stations.get(i));
            for (int j = 0; j < ke; j++) agg.jobclasses.add(s.jobclasses.get(j));
            if (s.sched != null) agg.sched.putAll(s.sched);
            for (int i = 0; i < me; i++)
                agg.nservers.set(roff + i, 0, s.nservers.get(i, 0));
            for (int i = 0; i < me; i++)
                for (int j = 0; j < ke; j++)
                    agg.rates.set(roff + i, coff + j, s.rates.get(i, j));
            for (int j = 0; j < ke; j++)
                agg.njobs.set(0, coff + j, s.njobs.get(0, j));
            int cr = (s.chains != null) ? Math.min(che, s.chains.getNumRows()) : 0;
            int cc = (s.chains != null) ? Math.min(ke, s.chains.getNumCols()) : 0;
            for (int c = 0; c < cr; c++)
                for (int j = 0; j < cc; j++)
                    agg.chains.set(choff + c, coff + j, s.chains.get(c, j));
            for (int j = 0; j < ke; j++)
                agg.refstat.set(coff + j, 0, s.refstat.get(j, 0) + roff);
            roff += me;
            coff += ke;
            choff += che;
        }
        return agg;
    }

    /**
     * Retrieves a network element by its name.
     * 
     * @param nodeName the name of the node to find
     * @return the LayeredNetworkElement with the specified name, or null if not found
     */
    public LayeredNetworkElement getNodeByName(String nodeName) {
        List<String> nodenames = this.getNodeNames();
        for (int idx = 0; idx < nodenames.size(); idx++) {
            if (nodenames.get(idx).equals(nodeName)) {
                return this.nodes.get(idx);
            }
        }
        return null;
    }

    /**
     * Gets the index of a network element in the node collection.
     * 
     * @param node the network element to find the index for
     * @return the index of the node, or -1 if not found
     */
    public Integer getNodeIndex(LayeredNetworkElement node) {
        List<String> nodenames = this.getNodeNames();
        String nodeName = node.getName();
        for (int idx = 0; idx < nodenames.size(); idx++) {
            if (nodenames.get(idx).equals(nodeName)) {
                return idx;
            }
        }
        return -1;
    }

    /**
     * Gets the names of all network elements (hosts, tasks, entries, activities).
     * 
     * @return a list containing all node names in the network
     */
    public List<String> getNodeNames() {
        List<String> nodenames = new ArrayList<>();
        for (int h = 0; h < this.hosts.size(); h++) {
            nodenames.add(this.hosts.get(h).getName());
        }
        for (int t = 0; t < this.tasks.size(); t++) {
            nodenames.add(this.tasks.get(t).getName());
        }
        for (int e = 0; e < this.entries.size(); e++) {
            nodenames.add(this.entries.get(e).getName());
        }
        for (int a = 0; a < this.activities.size(); a++) {
            nodenames.add(this.activities.get(a).getName());
        }
        return nodenames;
    }

    /**
     * Gets the number of layers in the layered network.
     *
     * @return the number of layers
     */
    public int getNumberOfLayers() {
        return getNumberOfModels();
    }

    /**
     * Gets the number of models in the layered network ensemble.
     * 
     * @return the number of models
     */
    public int getNumberOfModels() {
        if (this.ensemble.isEmpty()) {
            getEnsemble();
        }
        return this.ensemble.size();
    }

    /**
     * Gets the structural representation of the layered network.
     * 
     * @return the LayeredNetworkStruct containing the network structure
     */
    public LayeredNetworkStruct getStruct() {
        return this.getStruct(false);
    }

    /**
     * Gets the structural representation of the layered network.
     * 
     * @param regenerate if true, forces regeneration of the structure
     * @return the LayeredNetworkStruct containing the network structure
     */
    public LayeredNetworkStruct getStruct(boolean regenerate) {
        if (!regenerate & this.lsn != null) {
            return this.lsn;
        }
        LayeredNetworkStruct lsn = new LayeredNetworkStruct();

        lsn.nidx = 0;
        lsn.hshift = 0;
        lsn.nhosts = this.hosts.size();
        lsn.ntasks = this.tasks.size();
        lsn.nentries = this.entries.size();
        lsn.nacts = this.activities.size();
        lsn.tshift = lsn.nhosts;
        lsn.eshift = lsn.nhosts + lsn.ntasks;
        lsn.ashift = lsn.nhosts + lsn.ntasks + lsn.nentries;
        lsn.cshift = lsn.nhosts + lsn.ntasks + lsn.nentries + lsn.nacts;

        // analyze static properties
        lsn.nidx = lsn.nhosts + lsn.ntasks + lsn.nentries + lsn.nacts;
        // Element indices are 0-based: local h,t,e,a run 0..n-1 and the absolute
        // index is shift+local, so 0..nidx-1. See _kb/04-networkstruct.md
        int idx = 0;

        lsn.tasksof = new HashMap<>(lsn.nhosts);
        lsn.entriesof = new HashMap<>(lsn.nhosts + lsn.ntasks);
        lsn.actsof = new HashMap<>(lsn.nhosts + lsn.ntasks + lsn.nentries);

        lsn.callsof = new HashMap<>(lsn.nacts);

        lsn.hostdem = new HashMap<>();
        lsn.hostdem_type = new HashMap<>();
        lsn.hostdem_params = new HashMap<>();
        lsn.hostdem_mean = new HashMap<>();
        lsn.hostdem_scv = new HashMap<>();
        lsn.hostdem_proc = new HashMap<>();

        lsn.think = new HashMap<>();
        lsn.think_type = new HashMap<>();
        lsn.think_params = new HashMap<>();
        lsn.think_mean = new HashMap<>();
        lsn.think_scv = new HashMap<>();
        lsn.think_proc = new HashMap<>();

        lsn.actthink = new HashMap<>();
        lsn.actthink_type = new HashMap<>();
        lsn.actthink_params = new HashMap<>();
        lsn.actthink_mean = new HashMap<>();
        lsn.actthink_scv = new HashMap<>();
        lsn.actthink_proc = new HashMap<>();

        lsn.sched = new HashMap<>();

        lsn.names = new HashMap<>();
        lsn.hashnames = new HashMap<>();

        lsn.mult = new Matrix(1, lsn.nhosts + lsn.ntasks, lsn.nhosts + lsn.ntasks);
        lsn.maxmult = new Matrix(1, lsn.nhosts + lsn.ntasks, lsn.nhosts + lsn.ntasks);

        lsn.repl = new Matrix(1, lsn.nhosts + lsn.ntasks, lsn.nhosts + lsn.ntasks);
        // Task scheduling priority, read by the priority disciplines; 0 elsewhere.
        lsn.prio = new Matrix(1, lsn.nhosts + lsn.ntasks, lsn.nhosts + lsn.ntasks);
        lsn.type = new Matrix(1, lsn.nidx, lsn.nidx);
        lsn.graph = new Matrix(lsn.nidx, lsn.nidx, lsn.nidx * lsn.nidx);
        lsn.dag = new Matrix(lsn.nidx, lsn.nidx, lsn.nidx * lsn.nidx);
        lsn.replygraph = new Matrix(lsn.nacts, lsn.nentries, lsn.nentries * lsn.nacts);
        lsn.actphase = new Matrix(1, lsn.nacts, lsn.nacts);  // Phase for each activity (default=1)
        for (int a = 0; a < lsn.nacts; a++) {
            lsn.actphase.set(0, a, 1.0);  // Default phase is 1
        }

        lsn.nitems = new Matrix(1, lsn.nidx, lsn.ntasks + lsn.nacts);

        lsn.itemcap = new HashMap<>();
        lsn.itemproc = new HashMap<>();
        lsn.itemproc_type = new HashMap<>();
        lsn.itemproc_params = new HashMap<>();
        lsn.itemproc_mean = new HashMap<>();
        lsn.itemproc_scv = new HashMap<>();
        lsn.itemproc_proc = new HashMap<>();


        lsn.iscache = new Matrix(1, lsn.nidx, lsn.nhosts + lsn.ntasks);
        lsn.hasretrieval = new Matrix(1, lsn.nidx, lsn.nhosts + lsn.ntasks);
        lsn.replacestrat = new Matrix(1, lsn.nidx, lsn.nhosts + lsn.ntasks);

        lsn.setuptime = new HashMap<>();
        lsn.setuptime_type = new HashMap<>();
        lsn.setuptime_params = new HashMap<>();
        lsn.setuptime_mean = new HashMap<>();
        lsn.setuptime_scv = new HashMap<>();
        lsn.setuptime_proc = new HashMap<>();

        lsn.delayofftime = new HashMap<>();
        lsn.delayofftime_type = new HashMap<>();
        lsn.delayofftime_params = new HashMap<>();
        lsn.delayofftime_mean = new HashMap<>();
        lsn.delayofftime_scv = new HashMap<>();
        lsn.delayofftime_proc = new HashMap<>();

        lsn.hassetup = new Matrix(1, lsn.nidx, lsn.nidx);

        // Admission constraints on the layer station of a host or task -- see _kb/04-networkstruct.md
        lsn.lincon = new HashMap<>();

        // Service-rate dependences on the layer station of a host or task -- see _kb/04-networkstruct.md
        lsn.lldscaling = new HashMap<>();
        lsn.cdscaling = new HashMap<>();
        lsn.cdscalingpeak = new HashMap<>();
        lsn.jdscaling = new HashMap<>();
        lsn.jdscalingpeak = new HashMap<>();
        lsn.pools = new HashMap<>();

        lsn.arrival = new HashMap<>();
        lsn.arrival_type = new HashMap<>();
        lsn.arrival_params = new HashMap<>();
        lsn.arrival_mean = new HashMap<>();
        lsn.arrival_scv = new HashMap<>();
        lsn.arrival_proc = new HashMap<>();

        // A host has no parent. -1 is the unset sentinel: 0 is now the first host,
        // so it can no longer double as "absent" the way it did when idx was 1-based.
        lsn.parent = new Matrix(1, lsn.nidx, lsn.nidx);
        for (int i = 0; i < lsn.nidx; i++) {
            lsn.parent.set(0, i, -1);
        }

        for (int i = 0; i < lsn.nhosts; i++) {
            lsn.sched.put(idx, this.hosts.get(i).scheduling);
            // If multiplicity is Integer.MAX_VALUE, store as Inf (same convention as tasks below)
            int hostMult = this.hosts.get(i).multiplicity;
            lsn.mult.set(0, idx, hostMult == Integer.MAX_VALUE ? Inf : hostMult);
            lsn.repl.set(0, idx, this.hosts.get(i).replication);
            lsn.names.put(idx, this.hosts.get(i).getName());
            lsn.hashnames.put(idx, "P:" + lsn.names.get(idx));
            lsn.type.set(0, idx, LayeredNetworkElement.HOST);
            if (this.hosts.get(i).hasLinearConstraints()) {
                lsn.lincon.put(idx, this.hosts.get(i).getLinearConstraints());
            }
            idx = idx + 1;
        }

        for (int i = 0; i < lsn.ntasks; i++) {

            lsn.sched.put(idx, this.tasks.get(i).scheduling);
            lsn.hostdem.put(idx, Immediate.getInstance());
            DistParams hostdemParams = extractDistParams(Immediate.getInstance());
            lsn.hostdem_type.put(idx, hostdemParams.type);
            lsn.hostdem_params.put(idx, hostdemParams.params);
            lsn.hostdem_mean.put(idx, hostdemParams.mean);
            lsn.hostdem_scv.put(idx, hostdemParams.scv);
            lsn.hostdem_proc.put(idx, hostdemParams.proc);
            lsn.think.put(idx, this.tasks.get(i).thinkTime);
            DistParams thinkParams = extractDistParams(this.tasks.get(i).thinkTime);
            lsn.think_type.put(idx, thinkParams.type);
            lsn.think_params.put(idx, thinkParams.params);
            lsn.think_mean.put(idx, thinkParams.mean);
            lsn.think_scv.put(idx, thinkParams.scv);
            lsn.think_proc.put(idx, thinkParams.proc);
            // If multiplicity is Integer.MAX_VALUE, store as Inf (handles INF scheduling tasks)
            int taskMult = this.tasks.get(i).multiplicity;
            lsn.mult.set(0, idx, taskMult == Integer.MAX_VALUE ? Inf : taskMult);
            lsn.repl.set(0, idx, this.tasks.get(i).replication);
            lsn.prio.set(0, idx, this.tasks.get(i).getPriority());
            lsn.names.put(idx, this.tasks.get(i).getName());

            if (lsn.sched.get(idx) == SchedStrategy.REF) {
                lsn.hashnames.put(idx, "R:" + lsn.names.get(idx));
            } else {
                lsn.hashnames.put(idx, "T:" + lsn.names.get(idx));
            }

            if (this.tasks.get(i) instanceof CacheTask) {
                lsn.iscache.set(0, idx, 1);
                lsn.nitems.set(0, idx, ((CacheTask) this.tasks.get(i)).items);
                lsn.itemcap.put(idx, ((CacheTask) this.tasks.get(i)).getItemLevelCap());  // Now stores int[] for multi-level support
                lsn.replacestrat.set(0, idx, ((CacheTask) this.tasks.get(i)).replacestrategy.ordinal());
                if (((CacheTask) this.tasks.get(i)).hasRetrieval()) {
                    lsn.hasretrieval.set(0, idx, 1);
                }
                lsn.hashnames.put(idx, "C:" + lsn.names.get(idx));
            } else if (this.tasks.get(i).hasSetupDelayoff()) {
                // Task has setup/delayoff configured (not just SetupTask)
                Distribution setupDist = this.tasks.get(i).getSetupTime();
                lsn.setuptime.put(idx, setupDist);
                DistParams setupParams = extractDistParams(setupDist);
                lsn.setuptime_type.put(idx, setupParams.type);
                lsn.setuptime_params.put(idx, setupParams.params);
                lsn.setuptime_mean.put(idx, setupParams.mean);
                lsn.setuptime_scv.put(idx, setupParams.scv);
                lsn.setuptime_proc.put(idx, setupParams.proc);
                Distribution delayoffDist = this.tasks.get(i).getDelayOffTime();
                lsn.delayofftime.put(idx, delayoffDist);
                DistParams delayoffParams = extractDistParams(delayoffDist);
                lsn.delayofftime_type.put(idx, delayoffParams.type);
                lsn.delayofftime_params.put(idx, delayoffParams.params);
                lsn.delayofftime_mean.put(idx, delayoffParams.mean);
                lsn.delayofftime_scv.put(idx, delayoffParams.scv);
                lsn.delayofftime_proc.put(idx, delayoffParams.proc);
                lsn.hashnames.put(idx, "T:" + lsn.names.get(idx));
                lsn.hassetup.set(0, idx, 1);  // Store at idx to match setuptime/delayofftime indexing
            }

            int pidx = -1;
            Task currentTask = this.tasks.get(i);
            if (currentTask.parent == null) {
                line_error(mfilename(new Object() {}), "Task " + currentTask.getName() + " has no parent processor assigned during XML parsing");
            }
            for (int id = 0; id < this.hosts.size(); id++) {
                if (this.hosts.get(id).getName().equals(currentTask.parent.getName())) {
                    pidx = lsn.hshift + id;
                    break;
                }
            }

            lsn.parent.set(0, idx, pidx);
            lsn.graph.set(idx, pidx, 1);

            lsn.type.set(0, idx, LayeredNetworkElement.TASK);
            if (this.tasks.get(i).hasLinearConstraints()) {
                lsn.lincon.put(idx, this.tasks.get(i).getLinearConstraints());
            }
            idx++;
        }

        // Adjust task replication to account for host processor replication.
        // In LQN, task repl >= host repl. If task repl is 1 (default), inherit host repl.
        for (int t = 0; t < lsn.ntasks; t++) {
            int tidx = lsn.tshift + t;
            int pidx2 = (int) lsn.parent.get(0, tidx);
            double hostRepl = lsn.repl.get(0, pidx2);
            double taskRepl = lsn.repl.get(0, tidx);
            lsn.repl.set(0, tidx, Math.max(taskRepl, hostRepl));
        }

        // Build fan-out matrix from Task objects' fanOutMap
        // fanout(source_task_idx, dest_task_idx) = fan-out value (0 means not set)
        lsn.fanout = new Matrix(lsn.nidx, lsn.nidx);
        // Build task name -> index mapping
        Map<String, Integer> taskNameToIdx = new HashMap<>();
        for (int t = 0; t < lsn.ntasks; t++) {
            int tidx = lsn.tshift + t;
            taskNameToIdx.put(this.tasks.get(t).getName(), tidx);
        }
        // Populate fanout matrix
        for (int t = 0; t < lsn.ntasks; t++) {
            int tidx = lsn.tshift + t;
            Map<String, Integer> fanOutMap = this.tasks.get(t).getFanOutMap();
            for (Map.Entry<String, Integer> entry : fanOutMap.entrySet()) {
                Integer destIdx = taskNameToIdx.get(entry.getKey());
                if (destIdx != null) {
                    lsn.fanout.set(tidx, destIdx, entry.getValue());
                }
            }
        }

        for (int p = lsn.hshift; p < lsn.hshift + lsn.nhosts; p++) {

            if (!lsn.tasksof.containsKey(p)) {
                lsn.tasksof.put(p, new ArrayList<>());
            }
            for (int id = 0; id < lsn.parent.length(); id++) {

                if (lsn.parent.get(0, id) == p) {
                    lsn.tasksof.get(p).add(id);
                }
            }
        }

        for (int e = 0; e < lsn.nentries; e++) {
            idx = lsn.eshift + e;

            lsn.names.put(idx, this.entries.get(e).getName());

            // Extract open arrival distribution if present
            Entry entry = this.entries.get(e);
            if (entry.getArrival() != null) {
                Distribution arrivalDist = entry.getArrival();
                lsn.arrival.put(idx, arrivalDist);
                DistParams arrivalParams = extractDistParams(arrivalDist);
                lsn.arrival_type.put(idx, arrivalParams.type);
                lsn.arrival_params.put(idx, arrivalParams.params);
                lsn.arrival_mean.put(idx, arrivalParams.mean);
                lsn.arrival_scv.put(idx, arrivalParams.scv);
                lsn.arrival_proc.put(idx, arrivalParams.proc);
            }

            lsn.hashnames.put(idx, "E:" + lsn.names.get(idx));

            if (this.entries.get(e) instanceof ItemEntry) {
                lsn.hashnames.put(idx, "I:" + lsn.names.get(idx));
                lsn.nitems.set(0, idx, ((ItemEntry) this.entries.get(e)).getCardinality());
                Distribution popularity = ((ItemEntry) this.entries.get(e)).getPopularity();
                if (popularity != null) {
                    lsn.itemproc.put(idx, (DiscreteDistribution) popularity);
                    DistParams itemParams = extractDistParams(popularity);
                    lsn.itemproc_type.put(idx, itemParams.type);
                    lsn.itemproc_params.put(idx, itemParams.params);
                    lsn.itemproc_mean.put(idx, itemParams.mean);
                    lsn.itemproc_scv.put(idx, itemParams.scv);
                    lsn.itemproc_proc.put(idx, itemParams.proc);
                }
            } else {
                // Note: replygraph population from explicit repliesTo() is deferred
                // until after activity hashnames are populated (see below)
            }

            lsn.hostdem.put(idx, Immediate.getInstance());
            DistParams entryHostdemParams = extractDistParams(Immediate.getInstance());
            lsn.hostdem_type.put(idx, entryHostdemParams.type);
            lsn.hostdem_params.put(idx, entryHostdemParams.params);
            lsn.hostdem_mean.put(idx, entryHostdemParams.mean);
            lsn.hostdem_scv.put(idx, entryHostdemParams.scv);
            lsn.hostdem_proc.put(idx, entryHostdemParams.proc);
            int tidx = 0;
            Entry currentEntry = this.entries.get(e);
            if (currentEntry.parent == null) {
                line_error(mfilename(new Object() {}), "Entry " + currentEntry.getName() + " has no parent task assigned during XML parsing");
            }
            for (int id = 0; id < this.tasks.size(); id++) {
                if (currentEntry.parent.getName().equals(this.tasks.get(id).getName())) {
                    tidx = lsn.tshift + id;
                    break;
                }
            }
            lsn.parent.set(0, idx, tidx);
            lsn.graph.set(tidx, idx, 1);
            if (!lsn.entriesof.containsKey(tidx)) {
                lsn.entriesof.put(tidx, new ArrayList<>());

            }
            lsn.entriesof.get(tidx).add(idx);

            lsn.type.set(0, idx, LayeredNetworkElement.ENTRY);
            idx++;
        }

        // Admission constraint columns are only resolvable once tasksof/entriesof exist
        for (int cidx = 0; cidx < lsn.nhosts + lsn.ntasks; cidx++) {
            LayeredNetworkElement elem;
            List<Integer> colIdx;
            String colwhat;
            if (cidx < lsn.tshift) {
                elem = this.hosts.get(cidx - lsn.hshift);
                colIdx = lsn.tasksof.get(cidx);
                colwhat = "tasks on this host";
            } else {
                elem = this.tasks.get(cidx - lsn.tshift);
                colIdx = lsn.entriesof.get(cidx);
                colwhat = "entries of this task";
            }
            if (colIdx == null) {
                colIdx = new ArrayList<>();
            }
            int ncols = colIdx.size();
            // Rate dependences share the operand order of the constraint columns
            if (elem.lldScaling != null) {
                lsn.lldscaling.put(cidx, elem.lldScaling);
            }
            if (elem.lcdScaling != null) {
                lsn.cdscaling.put(cidx, elem.lcdScaling);
                lsn.cdscalingpeak.put(cidx, expandPeak(elem.lcdScalingPeak, ncols, lsn.names.get(cidx), colwhat, "Class"));
            }
            if (elem.ljdScaling != null) {
                lsn.jdscaling.put(cidx, elem.ljdScaling);
                lsn.jdscalingpeak.put(cidx, expandPeak(elem.ljdScalingPeak, ncols, lsn.names.get(cidx), colwhat, "Joint"));
            }
            // Compatibility pools name the operands they may serve, so the names
            // become columns only here, on the same operand order as the
            // constraints below.
            if (!elem.serverPools.isEmpty()) {
                List<String> poolCols = new ArrayList<>();
                for (int j = 0; j < ncols; j++) {
                    poolCols.add(lsn.names.get(colIdx.get(j)));
                }
                int npools = elem.serverPools.size();
                Matrix compat = new Matrix(npools, ncols);
                Matrix counts = new Matrix(1, npools);
                Matrix rates = new Matrix(1, npools);
                List<String> poolNames = new ArrayList<>();
                for (int t = 0; t < npools; t++) {
                    LayeredNetworkElement.ServerPool pool = elem.serverPools.get(t);
                    poolNames.add(pool.name);
                    counts.set(0, t, pool.count);
                    rates.set(0, t, pool.rate);
                    for (int k = 0; k < pool.compatible.size(); k++) {
                        int posCol = poolCols.indexOf(pool.compatible.get(k));
                        if (posCol < 0) {
                            throw new IllegalArgumentException("Server pool '" + pool.name + "' on "
                                    + lsn.names.get(cidx) + " names " + pool.compatible.get(k)
                                    + ", which is not one of the " + colwhat + ".");
                        }
                        compat.set(t, posCol, 1);
                    }
                }
                // An operand no pool can serve would be served at rate zero and
                // never complete, so it is a declaration error, not an empty column.
                for (int j = 0; j < ncols; j++) {
                    boolean served = false;
                    for (int t = 0; t < npools; t++) {
                        if (compat.get(t, j) != 0) {
                            served = true;
                        }
                    }
                    if (!served) {
                        throw new IllegalArgumentException(poolCols.get(j) + " on "
                                + lsn.names.get(cidx) + " is compatible with no server pool, so it "
                                + "can never be served.");
                    }
                }
                // The pools describe HOW the declared servers are shared, not how
                // many there are, so the two statements have to agree. Letting them
                // diverge would leave the layer station sized by the multiplicity
                // and scaled by a peak taken over a different number of servers,
                // reporting a utilization against a denominator never declared.
                double totalServers = 0;
                for (int t = 0; t < npools; t++) {
                    totalServers += counts.get(0, t);
                }
                double multc = lsn.mult.get(0, cidx);
                if (!Double.isInfinite(multc) && totalServers != multc) {
                    throw new IllegalArgumentException("Server pools on " + lsn.names.get(cidx)
                            + " hold " + totalServers + " servers but its multiplicity is " + multc
                            + "; the pools partition the declared servers, so the two must agree.");
                }
                lsn.pools.put(cidx, new LayeredNetworkStruct.ServerPools(poolNames, counts, rates, compat));
            }
            Matrix[] pos = lsn.lincon.get(cidx);
            if (pos != null && pos[0] != null && pos[0].getNumCols() != ncols) {
                throw new IllegalArgumentException("Admission constraint on " + lsn.names.get(cidx) + " has "
                        + pos[0].getNumCols() + " columns but there are " + ncols + " " + colwhat + ".");
            }
            if (elem.linConRows.isEmpty()) {
                continue;
            }
            // resolve rows declared by operand name against this server's columns
            List<String> colNames = new ArrayList<>();
            for (int j = 0; j < ncols; j++) {
                colNames.add(lsn.names.get(colIdx.get(j)));
            }
            int nnamed = elem.linConRows.size();
            Matrix Anamed = new Matrix(nnamed, ncols);
            Matrix bnamed = new Matrix(nnamed, 1);
            for (int r = 0; r < nnamed; r++) {
                LayeredNetworkElement.LinConRow namedRow = elem.linConRows.get(r);
                for (int k = 0; k < namedRow.names.size(); k++) {
                    int position = colNames.indexOf(namedRow.names.get(k));
                    if (position < 0) {
                        throw new IllegalArgumentException("Admission constraint on " + lsn.names.get(cidx) + " names "
                                + namedRow.names.get(k) + ", which is not one of the " + colwhat + ".");
                    }
                    Anamed.set(r, position, namedRow.coeffs[k]);
                }
                bnamed.set(r, 0, namedRow.cap);
            }
            if (pos == null || pos[0] == null || pos[1] == null) {
                lsn.lincon.put(cidx, new Matrix[]{Anamed, bnamed});
            } else {
                lsn.lincon.put(cidx, new Matrix[]{Matrix.concatRows(pos[0], Anamed, null),
                        Matrix.concatRows(pos[1], bnamed, null)});
            }
        }

        for (int a = 0; a < lsn.nacts; a++) {

            lsn.names.put(idx, this.activities.get(a).getName());
            lsn.hashnames.put(idx, "A:" + lsn.names.get(idx));
            Distribution actHostDemand = this.activities.get(a).hostDemand;
            lsn.hostdem.put(idx, actHostDemand);
            DistParams actHostdemParams = extractDistParams(actHostDemand);
            lsn.hostdem_type.put(idx, actHostdemParams.type);
            lsn.hostdem_params.put(idx, actHostdemParams.params);
            lsn.hostdem_mean.put(idx, actHostdemParams.mean);
            lsn.hostdem_scv.put(idx, actHostdemParams.scv);
            lsn.hostdem_proc.put(idx, actHostdemParams.proc);
            Distribution actThinkTime = this.activities.get(a).getThinkTime();
            lsn.actthink.put(idx, actThinkTime);
            DistParams actThinkParams = extractDistParams(actThinkTime);
            lsn.actthink_type.put(idx, actThinkParams.type);
            lsn.actthink_params.put(idx, actThinkParams.params);
            lsn.actthink_mean.put(idx, actThinkParams.mean);
            lsn.actthink_scv.put(idx, actThinkParams.scv);
            lsn.actthink_proc.put(idx, actThinkParams.proc);
            int tidx = 0;
            for (int id = 0; id < this.tasks.size(); id++) {
                if (this.activities.get(a).parent != null && this.activities.get(a).parent.getName().equals(this.tasks.get(id).getName())) {
                    tidx = lsn.tshift + id;
                    break;
                }
            }
            lsn.parent.set(0, idx, tidx);
            if (!lsn.actsof.containsKey(tidx)) {
                lsn.actsof.put(tidx, new ArrayList<>());
            }
            lsn.actsof.get(tidx).add(idx);
            lsn.type.set(0, idx, LayeredNetworkElement.ACTIVITY);
            // Store activity phase (1 or 2)
            lsn.actphase.set(0, a, this.activities.get(a).getPhase());
            idx++;
        }

        // see _kb/04-networkstruct.md (LayeredNetwork.getStruct() graph/validation rules) for rationale
        for (int e = 0; e < lsn.nentries; e++) {
            int eidx = lsn.eshift + e;
            for (String ra : this.entries.get(e).replyActivity.values()) {
                int ractidx = Utils.findString(lsn.hashnames, "A:" + ra);
                if (ractidx >= 0) {
                    lsn.replygraph.set(ractidx - lsn.ashift, eidx - lsn.eshift, 1);
                }
            }
        }

        lsn.graph.set(lsn.nidx - 1, lsn.nidx - 1, 0);

        Map<Integer, Task> tasks = this.tasks;
        // pre-incremented before each use, so the first call lands on index 0
        int cidx = -1;

        lsn.calltype = new HashMap<>();
        lsn.iscaller = new Matrix(lsn.nidx, lsn.nidx, (lsn.ntasks + lsn.nacts) * (lsn.ntasks + lsn.nentries));
        lsn.issynccaller = new Matrix(lsn.nidx, lsn.nidx, (lsn.ntasks + lsn.nacts) * (lsn.ntasks + lsn.nentries));
        lsn.isasynccaller = new Matrix(lsn.nidx, lsn.nidx, (lsn.ntasks + lsn.nacts) * (lsn.ntasks + lsn.nentries));
        lsn.callpair = new Matrix(lsn.nidx, 2, lsn.nidx * 2);
        lsn.callproc = new HashMap<>();
        lsn.callproc_type = new HashMap<>();
        lsn.callproc_params = new HashMap<>();
        lsn.callproc_mean = new HashMap<>();
        lsn.callproc_scv = new HashMap<>();
        lsn.callproc_proc = new HashMap<>();
        lsn.callnames = new HashMap<>();
        lsn.callhashnames = new HashMap<>();
        lsn.taskgraph = new Matrix(lsn.ntasks + lsn.tshift, lsn.ntasks + lsn.tshift, (lsn.ntasks + lsn.tshift) * (lsn.ntasks + lsn.tshift));
        lsn.actpretype = new Matrix(1, lsn.nidx, lsn.nacts);
        lsn.actposttype = new Matrix(1, lsn.nidx, lsn.nacts);
        lsn.actquorum = new Matrix(1, lsn.nidx, lsn.nacts);

        Matrix loop_back_edges = new Matrix(lsn.nidx, lsn.nidx, lsn.nidx * lsn.nidx);
        List<int[]> loopInfoList = new ArrayList<>();

        // Track boundToEntry mappings to validate uniqueness
        Map<String, String> entryToActivityMap = new HashMap<>();

        // Initialize callsof for all activities to prevent null access in SolverLN
        for (int i = 0; i < this.activities.size(); i++) {
            Activity activity = this.activities.get(i);
            int aidx = Utils.findString(lsn.hashnames, "A:" + activity.getName());
            if (aidx >= 0 && !lsn.callsof.containsKey(aidx)) {
                lsn.callsof.put(aidx, new ArrayList<>());
            }
        }

        for (int t = 0; t < lsn.ntasks; t++) {
            int tidx = lsn.tshift + t;

            for (int a = 0; a < tasks.get(t).activities.size(); a++) {
                int aidx = Utils.findString(lsn.hashnames, "A:" + tasks.get(t).activities.get(a).getName());
                if (aidx >= 0) {
                    lsn.callsof.put(aidx, new ArrayList<>());
                }

                String boundToEntry = tasks.get(t).activities.get(a).boundToEntry;
                int eidx = Utils.findString(lsn.hashnames, "E:" + boundToEntry);
                if (eidx < 0) {
                    eidx = Utils.findString(lsn.hashnames, "I:" + boundToEntry);
                }
                if (eidx >= 0) {
                    lsn.graph.set(eidx, aidx, 1);
                    
                    // Check if this entry is already bound to another activity
                    String activityName = tasks.get(t).activities.get(a).getName();
                    if (entryToActivityMap.containsKey(boundToEntry)) {
                        line_error(mfilename(new Object() {
                        }), "Multiple activities (" + entryToActivityMap.get(boundToEntry) + ", " + activityName + ") are bound to the same entry: " + boundToEntry);
                    } else {
                        entryToActivityMap.put(boundToEntry, activityName);
                    }
                }

                for (int s = 0; s < tasks.get(t).activities.get(a).syncCallDests.size(); s++) {

                    int target_eidx = Utils.findString(lsn.hashnames, "E:" + tasks.get(t).activities.get(a).syncCallDests.get(s));

                    if (target_eidx < 0) {
                        target_eidx = Utils.findString(lsn.hashnames, "I:" + tasks.get(t).activities.get(a).syncCallDests.get(s));
                    }
                    int target_tidx = (int) lsn.parent.get(target_eidx);
                    cidx++;

                    lsn.calltype.put(cidx, CallType.SYNC);
                    lsn.callpair.set(cidx, 0, aidx);
                    lsn.callpair.set(cidx, 1, target_eidx);
                    
                    // Check for self-call: if activity is bound to an entry and calls the same entry
                    boolean isSelfCall = false;
                    if (tasks.get(t).activities.get(a).boundToEntry != null) {
                        String boundEntryName = tasks.get(t).activities.get(a).boundToEntry;
                        String targetEntryName = lsn.hashnames.get(target_eidx);
                        // boundToEntry stores just the name, but hashnames uses "E:" prefix
                        if (targetEntryName != null && targetEntryName.startsWith("E:") && 
                            boundEntryName.equals(targetEntryName.substring(2))) {
                            isSelfCall = true;
                        }
                    }
                    
                    // Check for self-call first (more specific error)
                    if (isSelfCall) {
                        line_error(mfilename(new Object() {
                        }), "An entry calls itself, which creates an invalid self-referencing cycle.");
                    } else if (tidx == target_tidx) {
                        line_error(mfilename(new Object() {
                        }), "An entry on a task cannot call another entry on the same task.");
                    }
                    lsn.callnames.put(cidx, lsn.names.get(aidx) + "=>" + lsn.names.get(target_eidx));
                    lsn.callhashnames.put(cidx, lsn.hashnames.get(aidx) + "=>" + lsn.hashnames.get(target_eidx));
                    Distribution syncCallDist = callCountDist(tasks.get(t).activities.get(a).syncCallMeans.get(s));
                    lsn.callproc.put(cidx, syncCallDist);
                    DistParams syncCallParams = extractDistParams(syncCallDist);
                    lsn.callproc_type.put(cidx, syncCallParams.type);
                    lsn.callproc_params.put(cidx, syncCallParams.params);
                    lsn.callproc_mean.put(cidx, syncCallParams.mean);
                    lsn.callproc_scv.put(cidx, syncCallParams.scv);
                    lsn.callproc_proc.put(cidx, syncCallParams.proc);

                    lsn.callsof.get(aidx).add(cidx);
                    lsn.iscaller.set(aidx, target_tidx, 1);
                    lsn.iscaller.set(aidx, target_eidx, 1);
                    lsn.iscaller.set(tidx, target_tidx, 1);//1 -> true
                    lsn.iscaller.set(tidx, target_eidx, 1);
                    lsn.issynccaller.set(tidx, target_tidx, 1);
                    lsn.issynccaller.set(tidx, target_eidx, 1);
                    lsn.issynccaller.set(aidx, target_eidx, 1);
                    lsn.issynccaller.set(aidx, target_tidx, 1);
                    lsn.taskgraph.set(tidx, target_tidx, 1);
                    lsn.graph.set(aidx, target_eidx, 1);
                }

                for (Activity.CallGroup grp : tasks.get(t).activities.get(a).getSyncCallGroups()) {
                    List<Integer> gtargets = new ArrayList<>();
                    for (String dest : grp.dests) {
                        int geidx = Utils.findString(lsn.hashnames, "E:" + dest);
                        if (geidx < 0) {
                            geidx = Utils.findString(lsn.hashnames, "I:" + dest);
                        }
                        if (geidx >= 0) {
                            gtargets.add(geidx);
                        }
                    }
                    if (gtargets.size() >= 2) {
                        lsn.callgroups.add(new LayeredNetworkStruct.CallGroupStruct(
                                aidx, grp.strategy, gtargets));
                    }
                }

                for (int s = 0; s < tasks.get(t).activities.get(a).asyncCallDests.size(); s++) {
                    String target_entry_name = tasks.get(t).activities.get(a).asyncCallDests.get(s);
                    int target_eidx = Utils.findString(lsn.hashnames, "E:" + target_entry_name);
                    if (target_eidx < 0) {
                        target_eidx = Utils.findString(lsn.hashnames, "I:" + target_entry_name);
                    }
                    // Validate that the target entry exists
                    if (target_eidx < 0) {
                        line_error(mfilename(new Object() {}), "Activity \"" + tasks.get(t).activities.get(a).getName() + "\" has an async call to non-existent entry \"" + target_entry_name + "\".");
                    }
                    int target_tidx = (int) lsn.parent.get(target_eidx);
                    // Check for self-referential async calls (task calling itself)
                    if (tidx == target_tidx) {
                        line_error(mfilename(new Object() {}), "Activity \"" + tasks.get(t).activities.get(a).getName() + "\" in task \"" + tasks.get(t).getName() + "\" has an async call to an entry on the same task. Async self-calls are not supported.");
                    }
                    cidx++;

                    lsn.calltype.put(cidx, CallType.ASYNC);
                    lsn.callpair.set(cidx, 0, aidx);
                    lsn.callpair.set(cidx, 1, target_eidx);
                    lsn.callnames.put(cidx, lsn.names.get(aidx) + "->" + lsn.names.get(target_eidx));
                    lsn.callhashnames.put(cidx, lsn.hashnames.get(aidx) + "->" + lsn.hashnames.get(target_eidx));
                    Distribution asyncCallDist = callCountDist(tasks.get(t).activities.get(a).asyncCallMeans.get(s));
                    lsn.callproc.put(cidx, asyncCallDist);
                    DistParams asyncCallParams = extractDistParams(asyncCallDist);
                    lsn.callproc_type.put(cidx, asyncCallParams.type);
                    lsn.callproc_params.put(cidx, asyncCallParams.params);
                    lsn.callproc_mean.put(cidx, asyncCallParams.mean);
                    lsn.callproc_scv.put(cidx, asyncCallParams.scv);
                    lsn.callproc_proc.put(cidx, asyncCallParams.proc);
                    lsn.callsof.get(aidx).add(cidx);
                    lsn.iscaller.set(aidx, target_tidx, 1);
                    lsn.iscaller.set(aidx, target_eidx, 1);
                    lsn.iscaller.set(tidx, target_tidx, 1);//1 -> true
                    lsn.iscaller.set(tidx, target_eidx, 1);
                    lsn.isasynccaller.set(tidx, target_tidx, 1);
                    lsn.isasynccaller.set(tidx, target_eidx, 1);
                    lsn.isasynccaller.set(aidx, target_eidx, 1);
                    lsn.isasynccaller.set(aidx, target_tidx, 1);
                    lsn.taskgraph.set(tidx, target_tidx, 1);
                    lsn.graph.set(aidx, target_eidx, 1);
                }
            }
        }

        /* ========== Process forwarding calls from entries ========== */
        for (int e = 0; e < this.entries.size(); e++) {
            Entry entry = this.entries.get(e);
            int eidx = Utils.findString(lsn.hashnames, "E:" + entry.getName());
            if (eidx < 0) {
                eidx = Utils.findString(lsn.hashnames, "I:" + entry.getName());
            }
            if (eidx < 0) {
                continue;
            }
            int source_tidx = (int) lsn.parent.get(eidx);

            for (int fw = 0; fw < entry.getForwardingDests().size(); fw++) {
                String target_entry_name = entry.getForwardingDests().get(fw);
                int target_eidx = Utils.findString(lsn.hashnames, "E:" + target_entry_name);
                if (target_eidx < 0) {
                    target_eidx = Utils.findString(lsn.hashnames, "I:" + target_entry_name);
                }
                if (target_eidx < 0) {
                    line_error(mfilename(new Object() {}), "Entry \"" + entry.getName() + "\" forwards to non-existent entry \"" + target_entry_name + "\".");
                }
                int target_tidx = (int) lsn.parent.get(target_eidx);

                // Validate: cannot forward to same task
                if (source_tidx == target_tidx) {
                    line_error(mfilename(new Object() {}), "Entry \"" + entry.getName() + "\" cannot forward to entry \"" + target_entry_name + "\" on the same task.");
                }

                cidx++;
                lsn.calltype.put(cidx, CallType.FWD);
                lsn.callpair.set(cidx, 0, eidx);
                lsn.callpair.set(cidx, 1, target_eidx);
                lsn.callnames.put(cidx, lsn.names.get(eidx) + "~>" + lsn.names.get(target_eidx));
                lsn.callhashnames.put(cidx, lsn.hashnames.get(eidx) + "~>" + lsn.hashnames.get(target_eidx));

                // Forwarding probability (not using callproc as forwarding is deterministic choice)
                double fwdProb = entry.getForwardingProbs().get(fw);
                Distribution fwdCallDist = callCountDist(fwdProb);
                lsn.callproc.put(cidx, fwdCallDist); // Store as mean calls
                DistParams fwdCallParams = extractDistParams(fwdCallDist);
                lsn.callproc_type.put(cidx, fwdCallParams.type);
                lsn.callproc_params.put(cidx, fwdCallParams.params);
                lsn.callproc_mean.put(cidx, fwdCallParams.mean);
                lsn.callproc_scv.put(cidx, fwdCallParams.scv);
                lsn.callproc_proc.put(cidx, fwdCallParams.proc);

                // Update task graph to reflect forwarding relationship (matches MATLAB getStruct.m)
                lsn.taskgraph.set(source_tidx, target_tidx, 1);
                lsn.graph.set(eidx, target_eidx, 1);

                // see _kb/04-networkstruct.md (LayeredNetwork.getStruct() graph/validation rules) for rationale
            }
        }

        lsn.ncalls = cidx;  // Update total number of calls

        for (int t = 0; t < lsn.ntasks; t++) {
            int tidx = lsn.tshift + t;

            for (int a = 0; a < tasks.get(t).activities.size(); a++) {

                /* precedence handling ------------------------------------ */
                for (int ap = 0; ap < tasks.get(t).precedences.size(); ap++) {
                    String pretype = tasks.get(t).precedences.get(ap).preType;
                    String posttype = tasks.get(t).precedences.get(ap).postType;
                    List<String> preacts = tasks.get(t).precedences.get(ap).preActs;
                    List<String> postacts = tasks.get(t).precedences.get(ap).postActs;

                    // Validate PRE_AND activities exist before processing
                    if (pretype.equals(PRE_AND)) {
                        if (preacts.isEmpty()) {
                            line_error(mfilename(new Object() {}), "PRE_AND precedence in task \"" + tasks.get(t).getName() + "\" has no pre activities.");
                        }
                        for (int prea = 0; prea < preacts.size(); prea++) {
                            int preaidx = Utils.findString(lsn.hashnames, "A:" + preacts.get(prea));
                            if (preaidx < 0) {
                                line_error(mfilename(new Object() {}), "PRE_AND precedence references non-existent activity \"" + preacts.get(prea) + "\" in task \"" + tasks.get(t).getName() + "\".");
                            }
                            if (preaidx >= 0 && lsn.parent.get(preaidx) != tidx) {
                                line_error(mfilename(new Object() {}), "PRE_AND precedence in task \"" + tasks.get(t).getName() + "\" references activity \"" + preacts.get(prea) + "\" from a different task.");
                            }
                        }

                        // see _kb/04-networkstruct.md (LayeredNetwork.getStruct() graph/validation rules) for rationale
                        Set<String> joinEntries = new HashSet<String>();
                        for (String pa : preacts) {
                            for (Activity act : tasks.get(t).activities) {
                                if (act.getName().equals(pa)) {
                                    if (act.boundToEntry != null && !act.boundToEntry.isEmpty()) {
                                        joinEntries.add(act.boundToEntry);
                                    }
                                    break;
                                }
                            }
                        }
                        if (joinEntries.size() > 1) {
                            line_error(mfilename(new Object() {}), "AND-join in task \""
                                + tasks.get(t).getName() + "\" synchronizes activities bound to "
                                + "different entries " + joinEntries + ": cross-entry synchronization "
                                + "(LQN synchronization server) is not supported.");
                        }
                    }

                    // Validate POST_AND activities exist before processing
                    if (posttype.equals(POST_AND)) {
                        if (postacts.isEmpty()) {
                            line_error(mfilename(new Object() {}), "POST_AND precedence in task \"" + tasks.get(t).getName() + "\" has no post activities.");
                        }
                        for (int posta = 0; posta < postacts.size(); posta++) {
                            int postaidx = Utils.findString(lsn.hashnames, "A:" + postacts.get(posta));
                            if (postaidx < 0) {
                                line_error(mfilename(new Object() {}), "POST_AND precedence references non-existent activity \"" + postacts.get(posta) + "\" in task \"" + tasks.get(t).getName() + "\".");
                            }
                            if (postaidx >= 0 && lsn.parent.get(postaidx) != tidx) {
                                line_error(mfilename(new Object() {}), "POST_AND precedence in task \"" + tasks.get(t).getName() + "\" references activity \"" + postacts.get(posta) + "\" from a different task.");
                            }
                        }
                    }

                    // see _kb/04-networkstruct.md (LayeredNetwork.getStruct() graph/validation rules) for rationale
                    int quorumCount = pretype.equals(PRE_AND)
                            ? ActivityPrecedence.getQuorumCount(tasks.get(t).precedences.get(ap).preParams, preacts.size())
                            : 0;

                    for (int prea = 0; prea < preacts.size(); prea++) {
                        int preaidx = Utils.findString(lsn.hashnames, "A:" + tasks.get(t).precedences.get(ap).preActs.get(prea));
                        double preParam = 1.0;

                        switch (posttype) {
                            case POST_OR:
                                for (int posta = 0; posta < postacts.size(); posta++) {
                                    int postaidx = Utils.findString(lsn.hashnames, "A:" + tasks.get(t).precedences.get(ap).postActs.get(posta));
                                    Matrix probs = tasks.get(t).precedences.get(ap).postParams;
                                    double postParam = probs.get(posta);
                                    lsn.graph.set(preaidx, postaidx, preParam * postParam);
                                    lsn.actpretype.set(0, preaidx, ActivityPrecedence.getPrecedenceId(tasks.get(t).precedences.get(ap).preType));
                                    lsn.actposttype.set(0, postaidx, ActivityPrecedence.getPrecedenceId(tasks.get(t).precedences.get(ap).postType));
                                }
                                break;
                            case POST_AND:
                                for (int posta = 0; posta < postacts.size(); posta++) {
                                    int postaidx = Utils.findString(lsn.hashnames, "A:" + tasks.get(t).precedences.get(ap).postActs.get(posta));
                                    lsn.graph.set(preaidx, postaidx, 1.0);
                                    lsn.actpretype.set(0, preaidx, ActivityPrecedence.getPrecedenceId(tasks.get(t).precedences.get(ap).preType));
                                    lsn.actposttype.set(0, postaidx, ActivityPrecedence.getPrecedenceId(tasks.get(t).precedences.get(ap).postType));
                                }
                                break;
                            case POST_LOOP:
                                Matrix counts = tasks.get(t).precedences.get(ap).postParams;
                                int loopEntryAidx = Utils.findString(lsn.hashnames, "A:" + tasks.get(t).precedences.get(ap).preActs.get(0));
                                int loopStartAidx = Utils.findString(lsn.hashnames, "A:" + tasks.get(t).precedences.get(ap).postActs.get(0));
                                int loopEndAidx = Utils.findString(lsn.hashnames, "A:" + tasks.get(t).precedences.get(ap).postActs.get(postacts.size() - 1));

                                if (counts.value() < 1) {
                                    // When expected iterations < 1, we may skip loop entirely
                                    // E[iterations] = counts means P(enter loop) = counts
                                    lsn.graph.set(loopEntryAidx, loopStartAidx, counts.value());
                                    lsn.graph.set(loopEntryAidx, loopEndAidx, 1.0 - counts.value());
                                    // Process activities inside the loop as serial chain
                                    int curAidx = loopStartAidx;
                                    for (int posta = 1; posta < postacts.size() - 1; posta++) {
                                        int postaidx = Utils.findString(lsn.hashnames, "A:" + tasks.get(t).precedences.get(ap).postActs.get(posta));
                                        lsn.graph.set(curAidx, postaidx, 1.0);
                                        lsn.actposttype.set(0, postaidx, ActivityPrecedence.getPrecedenceId(tasks.get(t).precedences.get(ap).postType));
                                        curAidx = postaidx;
                                    }
                                    // After loop body, always exit to end (no looping back)
                                    lsn.graph.set(curAidx, loopEndAidx, 1.0);
                                    lsn.actposttype.set(0, loopStartAidx, ActivityPrecedence.getPrecedenceId(tasks.get(t).precedences.get(ap).postType));
                                } else {
                                    // When expected iterations >= 1, always enter loop
                                    // E[iterations] = 1/(1-p) = counts => p = 1 - 1/counts
                                    int curAidx = loopEntryAidx;
                                    for (int posta = 0; posta < postacts.size() - 1; posta++) {
                                        int postaidx = Utils.findString(lsn.hashnames, "A:" + tasks.get(t).precedences.get(ap).postActs.get(posta));
                                        lsn.graph.set(curAidx, postaidx, 1.0);
                                        lsn.actposttype.set(0, postaidx, ActivityPrecedence.getPrecedenceId(tasks.get(t).precedences.get(ap).postType));
                                        curAidx = postaidx;
                                    }
                                    loop_back_edges.set(curAidx, loopStartAidx, 1);
                                    lsn.graph.set(curAidx, loopStartAidx, 1.0 - 1.0 / counts.value());
                                    lsn.graph.set(curAidx, loopEndAidx, 1.0 / counts.value());
                                    loopInfoList.add(new int[]{loopStartAidx, loopEndAidx});
                                }
                                lsn.actposttype.set(0, loopEndAidx, ActivityPrecedence.getPrecedenceId(tasks.get(t).precedences.get(ap).postType));
                                break;
                            default:
                                for (int posta = 0; posta < postacts.size(); posta++) {
                                    int postaidx = Utils.findString(lsn.hashnames, "A:" + tasks.get(t).precedences.get(ap).postActs.get(posta));
                                    double postParam = 1;
                                    lsn.graph.set(preaidx, postaidx, preParam * postParam);
                                    lsn.actpretype.set(0, preaidx, ActivityPrecedence.getPrecedenceId(tasks.get(t).precedences.get(ap).preType));
                                    lsn.actposttype.set(0, postaidx, ActivityPrecedence.getPrecedenceId(tasks.get(t).precedences.get(ap).postType));
                                    if (quorumCount > 0) {
                                        lsn.actquorum.set(0, postaidx, quorumCount);
                                    }
                                }
                        }
                    }
                }
            }

        }

        /* Compute entry-to-activity reachability within the same task */
        for (int eoff = 0; eoff < lsn.nentries; eoff++) {
            int eidx = lsn.eshift + eoff; // global entry index
            int tidx = (int) lsn.parent.get(eidx);
            boolean[] visited = new boolean[lsn.nidx];
            Deque<Integer> stack = new ArrayDeque<>();
            stack.push(eidx);
            visited[eidx] = true;
            while (!stack.isEmpty()) {
                int v = stack.pop();
                for (int nbr = 0; nbr < lsn.nidx; nbr++) {
                    if (lsn.graph.get(v, nbr) != 0 && !visited[nbr]) {
                        visited[nbr] = true;
                        stack.push(nbr);
                    }
                }
            }
            List<Integer> acts = new ArrayList<>();
            for (int n = 0; n < lsn.nidx; n++) {
                if (visited[n] && lsn.type.get(n) == LayeredNetworkElement.ACTIVITY && lsn.parent.get(n) == tidx) {
                    acts.add(n);
                }
            }
            lsn.actsof.put(eidx, acts);
        }

        for (int t = 0; t < lsn.ntasks; t++) {
            int tidx = lsn.tshift + t;
            for (int aidx : lsn.actsof.getOrDefault(tidx, new ArrayList<>())) {
                List<Integer> postaidxs = new ArrayList<>();
                for (int col = 0; col < lsn.graph.getNumCols(); col++) {
                    if (lsn.graph.get(aidx, col) != 0) {
                        postaidxs.add(col);
                    }
                }

                boolean isreply = true;
                for (int postaidx : postaidxs) {
                    for (int acts : lsn.actsof.get(tidx)) {
                        if (acts == postaidx) {
                            isreply = false;
                            break;
                        }
                    }
                }
                if (isreply) {
                    int parentidx = aidx;
                    while (lsn.type.get(parentidx) != LayeredNetworkElement.ENTRY) {
                        List<Integer> ancestors = new ArrayList<>();
                        for (int row = 0; row < lsn.graph.getNumRows(); row++) {
                            if (lsn.graph.get(row, parentidx) != 0) {
                                ancestors.add(row);
                            }
                        }
                        if (ancestors.isEmpty()) {
                            // No ancestors found, break out of the loop
                            break;
                        }
                        parentidx = ancestors.get(0);
                    }
                    if (lsn.type.get(parentidx) == LayeredNetworkElement.ENTRY) {
                        lsn.replygraph.set(aidx - lsn.ashift, parentidx - lsn.eshift, 1);
                    }
                }
            }
        }

        lsn.ncalls = lsn.calltype.size();
        List<Integer> toRemove = new ArrayList<>();
        for (
                int i = lsn.ncalls; i < lsn.callpair.getNumRows(); i++) {
            toRemove.add(i);
        }
        lsn.callpair.removeRows(toRemove);

        List<Integer> tidxs = new ArrayList<>();
        for (
                int i = 0; i < lsn.sched.size(); i++) {
            if (lsn.sched.get(i) == SchedStrategy.REF) {
                tidxs.add(i);
            }
        }

        for (int tidx : tidxs) {
            if (lsn.type.get(tidx) == LayeredNetworkElement.TASK) {
                List<Integer> callers = new ArrayList<>();
                for (int row = 0; row < lsn.taskgraph.getNumRows(); row++) {
                    if (lsn.taskgraph.get(row, tidx) != 0) {
                        callers.add(row);
                    }
                }
                List<Integer> callers_inf = new ArrayList<>();
                for (int caller : callers) {
                    if (lsn.mult.get(caller) < 0) {
                        callers_inf.add(1);
                    } else {
                        callers_inf.add(0);
                    }
                }
            }
        }

        lsn.isref = new Matrix(1, lsn.nhosts + lsn.ntasks, lsn.ntasks);
        for (int col = 0; col < lsn.sched.size(); col++) {
            if (lsn.sched.get(col) == SchedStrategy.REF) {
                lsn.isref.set(0, col, 1);
            }
        }

        // Create schedid matrix (scheduling strategy ordinal values) for Python compatibility
        lsn.schedid = new Matrix(1, lsn.nhosts + lsn.ntasks, lsn.nhosts + lsn.ntasks);
        for (int col = 0; col < lsn.sched.size(); col++) {
            SchedStrategy strategy = lsn.sched.get(col);
            if (strategy != null) {
                lsn.schedid.set(0, col, strategy.ordinal());
            }
        }

        // Create replacement alias for replacestrat (for Python compatibility)
        lsn.replacement = lsn.replacestrat;

        // Only process cache items if there are any
        if (lsn.nitems.getNumCols() > 0) {
            for (int i = 0; i < lsn.nitems.getNumCols(); i++) {
                if (lsn.nitems.get(0, i) > 0) {
                    lsn.iscache.set(0, i, 1);
                }
            }
        }

        // dag vs graph differences: see _kb/04-networkstruct.md (LayeredNetworkStruct field tables)
        Matrix dag = lsn.graph.copy();
        // Reverse edges from TASK to ENTRY for non-reference tasks
        // This enables proper flow propagation in lsn_max_multiplicity
        for (int i = 0; i < lsn.nidx; i++) {
            if (lsn.type.get(i) == LayeredNetworkElement.TASK &&
                    lsn.isref.get(0, i) == 0) {
                for (int j = 0; j < lsn.nidx; j++) {
                    if (lsn.type.get(j) == LayeredNetworkElement.ENTRY && dag.get(i, j) != 0) {
                        dag.set(i, j, 0);
                        dag.set(j, i, 1);
                    }
                }
            }
        }
        for (int r = 0; r < lsn.nidx; r++) {
            for (int c = 0; c < lsn.nidx; c++) {
                if (loop_back_edges.get(r, c) != 0) {
                    dag.set(r, c, 0);
                }
            }
        }

        lsn.dag = dag;

        int newRow = lsn.taskgraph.getNumCols() - lsn.nhosts;
        int newCol = lsn.taskgraph.getNumRows() - lsn.nhosts;
        Matrix taskgraphSection = new Matrix(newRow, newCol);
        for (
                int r = 0;
                r < newRow; r++) {
            for (int c = 0; c < newCol; c++) {
                taskgraphSection.set(r, c, lsn.taskgraph.get(lsn.tshift + r, lsn.tshift + c));
            }
        }

        Matrix labels = new Matrix(1, newCol);
        labels.zero();
        int ccc = 0;
        List<Double> roots = new ArrayList<>();
        List<Double> vectorList = new ArrayList<>();
        List<Double> vectorListNew;
        // graph_connected_components
        while (!labels.findZero().

                isEmpty()) {
            Matrix ccZero = labels.findZero();
            double fue = ccZero.value(); // first unexplored vertex
            roots.add(fue);
            vectorList.add(fue);
            ccc++;
            labels.set(0, (int) fue, ccc);
            while (!vectorList.isEmpty()) {
                vectorListNew = new ArrayList<>();
                for (int lc = 0; lc < vectorList.size(); lc++) {
                    int point = vectorList.get(lc).intValue();
                    Matrix connectedPoint = Matrix.extractRows(taskgraphSection, point, point + 1, null);
                    connectedPoint = connectedPoint.find();
                    List<Double> labelConnectedPoints = new ArrayList<>();
                    for (Double cp : connectedPoint.toList1D()) {
                        labelConnectedPoints.add(labels.get(0, cp.intValue()));
                    }
                    List<Double> intermediate = new ArrayList<>();
                    for (int lcpIdx = 0; lcpIdx < labelConnectedPoints.size(); lcpIdx++) {
                        if (labelConnectedPoints.get(lcpIdx) == 0) {
                            intermediate.add((double) lcpIdx);
                        }
                    }
                    List<Double> cp1 = new ArrayList<>();
                    for (Double icp : intermediate) {
                        cp1.add(connectedPoint.get(icp.intValue()));
                    }
                    for (Double cp1value : cp1) {
                        labels.set(0, cp1value.intValue(), ccc);
                    }
                    vectorListNew.addAll(cp1);
                }
                vectorList = new ArrayList<>(vectorListNew);
            }

        }
        //  [conncomps, roots]=graph_connected_components(lsn.taskgraph(lsn.nhosts+1:end, lsn.nhosts+1:end));
        // Map each component label to the absolute index of that component's root task.
        // The result gets its own matrix rather than overwriting `labels` in place: a
        // rewritten entry holds an absolute index, which can coincide with a later
        // component label and be relabelled a second time.
        lsn.conntasks = new Matrix(1, newCol);
        for (
                int r = 1; r < roots.size() + 1; r++) {
            for (int ctidx = 0; ctidx < labels.length(); ctidx++) {
                if (labels.get(0, ctidx) == r) {
                    lsn.conntasks.set(0, ctidx, lsn.tshift + roots.get(r - 1));
                }
            }
        }


        dag = lsn.dag.copy();
        for (int i = 0; i < dag.getNumRows(); i++) {
            for (int j = 0; j < dag.getNumCols(); j++) {
                if (loop_back_edges.get(i, j) != 0) {
                    dag.set(i, j, 0);              // cut the back edge
                }
            }
        }
        lsn.dag = dag;

        // see _kb/04-networkstruct.md (Python native layered.py: LayeredNetwork.getStruct() port notes) for rationale
        for (int tidx = 0; tidx < lsn.nhosts + lsn.ntasks; tidx++) {
            if (lsn.sched.get(tidx) == SchedStrategy.INF &&
                lsn.type.get(tidx) == LayeredNetworkElement.TASK) {
                // Sum caller multiplicities (Inf + finite = Inf, matching MATLAB behavior)
                double totalCallerMult = 0;
                for (int row = 0; row < lsn.taskgraph.getNumRows(); row++) {
                    if (lsn.taskgraph.get(row, tidx) != 0) {
                        totalCallerMult += lsn.mult.get(0, row);
                    }
                }
                if (totalCallerMult > 0) {
                    lsn.mult.set(0, tidx, totalCallerMult);
                }
                // Otherwise keep the original value (Inf for INF scheduling with no callers)
            }
        }

        // Compute bounds on multiplicies for host processors and non-ref tasks
        if (DirectedGraph.isDAG(dag)) {               // topological-sort based test
            lsn.maxmult = lsnMaxMultiplicity(lsn).transpose();
        } else {
            line_error(mfilename(new Object() {
            }), "A cycle exists in an activity graph.");
        }

        // every entry must have a boundTo activity: getStruct.m's guard, absent here
        // until 2026-08-15. An entry with an empty <entry-phase-activities> reaches no
        // activity, so it has no service and no reply; it used to build a struct and
        // report a row of NaN instead of being named. Checked before the reply guard
        // below, the order getStruct.m refuses them in.
        for (int e = 0; e < lsn.nentries; e++) {
            int eidx = lsn.eshift + e;
            boolean bound = false;
            for (int succ = lsn.ashift; succ < lsn.nidx; succ++) {
                if (lsn.graph.get(eidx, succ) != 0) {
                    bound = true;
                    break;
                }
            }
            if (!bound) {
                line_error(mfilename(new Object() {
                }), "An entry does not have any boundTo activity.");
            }
        }

        // non-terminal reply activity validity: see _kb/06-solver-catalog.md ("Activity-graph validity and .lqnx writer rules")
        for (int a = 0; a < lsn.replygraph.getNumRows(); a++) {
            for (int b = 0; b < lsn.replygraph.getNumCols(); b++) {
                if (lsn.replygraph.get(a, b) > 0) {          // activity 'a' replies
                    int aidx = lsn.ashift + a;                // global activity index
                    for (int succ = 0; succ < lsn.graph.getNumCols(); succ++) {
                        if (lsn.graph.get(aidx, succ) != 0 && succ >= lsn.ashift) {
                            int succActIdx = succ - lsn.ashift;   // 0-based activity number
                            if (succActIdx >= 0 && succActIdx < lsn.nacts) {
                                int succPhase = (int) lsn.actphase.get(0, succActIdx);
                                if (succPhase != 2) {             // phase-1 successor => invalid
                                    line_error(mfilename(new Object() {
                                    }), "Unsupported replyTo in non-terminal activity.");
                                }
                            }
                        }
                    }
                }
            }
        }

        // Call check: an entry may be called either synchronously or asynchronously, but not both.
        if (lsn.callpair != null && lsn.callpair.getNumRows() > 0) {

            /* Remember the first call-kind we encounter for each entry */
            Map<Integer, CallType> entryCallKind = new HashMap<>();

            int nCalls = lsn.callpair.getNumRows();
            for (int iter_cidx = 0; iter_cidx < nCalls; iter_cidx++) {

                int targetEidx = (int) lsn.callpair.get(iter_cidx, 1);   // callee entry column
                CallType kind = lsn.calltype.get(iter_cidx);            // SYNC, ASYNC, or FWD

                // Skip FWD calls - they're transformed to SYNC pseudo calls
                if (kind == CallType.FWD) {
                    continue;
                }

                CallType prev = entryCallKind.putIfAbsent(targetEidx, kind);
                if (prev != null && prev != kind) {
                    line_error(mfilename(new Object() {
                    }), "An entry is called both synchronously and asynchronously.");
                }
            }
        }

        lsn.nidx = lsn.nhosts + lsn.ntasks + lsn.nentries + lsn.nacts;
        this.lsn = lsn.copy();
        return lsn;
    }

    /**
     * Initializes the layered network by generating the graph and setting default parameters.
     */
    public void init() {
        this.generateGraph();
        this.initDefault();
        this.param.Nodes.RespT = 0;
        this.param.Nodes.Tput = 0;
        this.param.Nodes.Util = 0;
        this.param.Nodes.QLen = 0;
        this.param.Edges.RespT = 0;
        this.param.Edges.Tput = 0;
        this.param.Edges.QLen = 0;
    }

    /**
     * Initializes the layered network with default settings.
     * Currently this method contains placeholder implementation.
     */
    public void initDefault() {
        // TODO: is it necessary to have a version where, per LINE, nodes can be passed in as a parameter?
        // Check in particular how in dev/ that particular version is then used
    }

    /**
     * Resets the layered network to its initial state.
     * 
     * @param isHard if true, performs a hard reset clearing all network elements;
     *               if false, performs a soft reset clearing only the ensemble
     */
    public void reset(boolean isHard) {
        this.ensemble = new ArrayList<>();
        // The warm start is indexed by layer, so it cannot outlive the layers.
        this.initMarginalBlocks = null;
        if (isHard) {
            this.hosts = new HashMap<>();
            this.activities = new HashMap<>();
            this.tasks = new HashMap<>();
            this.reftasks = new HashMap<>();
            this.entries = new HashMap<>();
            this.nodes = new HashMap<>();
        }
    }

    /**
     * Sends the model to a local server using the default IP address (127.0.0.1).
     * 
     * @param outputPath the path where the model output should be stored
     * @param portNumber the port number for the server connection
     */
    public void sendModel(String outputPath, String portNumber) {
        this.sendModel(outputPath, "127.0.0.1", portNumber);
    }

    /**
     * Sends the model to a server at the specified IP address and port.
     * 
     * @param outputPath the path where the model output should be stored
     * @param ipNumber the IP address of the server
     * @param portNumber the port number for the server connection
     */
    public void sendModel(String outputPath, String ipNumber, String portNumber) {
        String filePath = null;
        try {
            filePath = lineTempName("layered");
        } catch (IOException e) {
            return;
        }
        try {
            this.writeXML(filePath + "/model.xml", true);
            LineDockerClient.sendModel(filePath + "/model.xml", outputPath, ipNumber, portNumber);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Prints a summary of the layered network structure to standard output.
     */
    public void summary() {
        LayeredNetworkStruct this_lqn = getStruct();
        this_lqn.print();
    }

    /**
     * Views the layered network model using ModelVisualizer.
     */
    public void view() {
        plot();
    }

    /**
     * Displays an interactive graph visualization of this layered network using JUNG.
     *
     * <p>The visualization shows:
     * <ul>
     * <li>Hosts as red pyramids</li>
     * <li>Tasks as red parallelograms</li>
     * <li>Entries as red rectangles</li>
     * <li>Activities as red circles</li>
     * </ul>
     *
     * <p>The window includes menus for:
     * <ul>
     * <li>Switching between different layout algorithms (Hierarchical, Circle, Force-Directed, etc.)</li>
     * <li>Switching mouse modes (Pan/Zoom vs Pick/Move)</li>
     * </ul>
     */
    public void plot() {
        plot(getName(), 800, 600);
    }

    /**
     * Displays an interactive graph visualization with custom title.
     *
     * @param title the window title
     */
    public void plot(String title) {
        plot(title, 800, 600);
    }

    /**
     * Displays an interactive graph visualization with custom title and dimensions.
     *
     * @param title  the window title
     * @param width  the window width in pixels
     * @param height the window height in pixels
     */
    public void plot(String title, int width, int height) {
        ModelVisualizer visualizer = new ModelVisualizer(this);
        visualizer.buildGraph();
        visualizer.show(title, width, height);
    }

    /**
     * Writes the layered network to a JLQN file with default naming.
     * 
     * @param filename the path to write the JLQN file
     * @throws Exception if there's an error writing the file
     */
    public void writeJLQN(String filename) throws Exception {
        writeJLQN(filename, false);
    }

    /**
     * Writes the layered network to a JLQN file.
     * 
     * @param filename the path to write the JLQN file
     * @param abstractNames if true, uses abstract names in the output
     * @throws Exception if there's an error writing the file
     */
    public void writeJLQN(String filename, boolean abstractNames) throws Exception {
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        DocumentBuilder builder = factory.newDocumentBuilder();

        // Create a new document
        Document document = builder.newDocument();

        // Create XML declaration
        DOMImplementation domImpl = document.getImplementation();
        document.appendChild(domImpl.createDocumentType("jlqn", null, "JLQNmodel.xsd"));

        // Root element
        Element root = document.createElement("jlqn");
        root.setAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
        root.setAttribute("xsi:noNamespaceSchemaLocation", "JLQNmodel.xsd");
        document.appendChild(root);

        // Parameters element
        Element parameters = document.createElement("parameters");
        root.appendChild(parameters);

        // Processors
        Element processors = document.createElement("processors");
        processors.setAttribute("number", "" + this.hosts.size());
        parameters.appendChild(processors);

        for (int h = 0; h < this.hosts.size(); h++) {
            Host H = this.hosts.get(h);
            writeJLQNProcessor(document, processors, H.getName(), H.scheduling.toString(), H.multiplicity);
        }

        // Tasks
        Element tasks = document.createElement("tasks");
        tasks.setAttribute("number", "" + this.tasks.size());
        parameters.appendChild(tasks);

        for (int t = 0; t < this.tasks.size(); t++) {
            Task T = this.tasks.get(t);
            writeJLQNTask(document, tasks, T);
        }

        // Entries
        Element entries = document.createElement("entries");
        entries.setAttribute("number", "" + this.entries.size());
        parameters.appendChild(entries);

        for (int e = 0; e < this.entries.size(); e++) {
            Entry E = this.entries.get(e);
            String E_boundTo = "";
            String E_replyTo = "";
            for (int a = 0; a < this.activities.size(); a++) {
                Activity A = this.activities.get(a);
                if (A.boundToEntry.compareTo(E.getName()) == 0) {
                    E_boundTo = A.getName();
                    // Get the first reply activity value (if any exist)
                    if (!E.replyActivity.isEmpty()) {
                        E_replyTo = E.replyActivity.values().iterator().next();
                    }
                }
            }
            writeJLQNEntry(document, entries, E.getName(), E_boundTo, E_replyTo, E.parent.getName(),
                E.getForwardingDests(), E.getForwardingProbs());
        }

        // Activities
        Element activities = document.createElement("activities");
        activities.setAttribute("number", "" + this.activities.size());
        parameters.appendChild(activities);

        for (int t = 0; t < this.activities.size(); t++) {
            Activity A = this.activities.get(t);
            writeJLQNActivity(document, activities, A.getName(), A.parent.getName(), "" + A.hostDemand.getMean());
        }

        // Calls
        Element calls = document.createElement("calls");
        int ncalls = 0;
        for (int t = 0; t < this.activities.size(); t++) {
            Activity A = this.activities.get(t);
            ncalls += A.syncCallDests.size() + A.asyncCallDests.size();
        }
        calls.setAttribute("number", "" + ncalls);
        parameters.appendChild(calls);

        for (int t = 0; t < this.activities.size(); t++) {
            Activity A = this.activities.get(t);
            for (int c = 0; c < A.syncCallDests.size(); c++) {
                String call_dest = A.syncCallDests.get(c);
                writeJLQNCall(document, calls, A.getName() + "=>" + call_dest, A.getName(), call_dest, CallType.SYNC, A.syncCallMeans.get(c));
            }
            for (int c = 0; c < A.asyncCallDests.size(); c++) {
                String call_dest = A.asyncCallDests.get(c);
                writeJLQNCall(document, calls, A.getName() + "->" + call_dest, A.getName(), call_dest, CallType.ASYNC, A.asyncCallMeans.get(c));
            }
            // Forwarding calls are per-entry (emitted as <forwarding> children of
            // <entry> in writeJLQNEntry), not per-activity, so nothing to do here.
        }

        // Precedences
        Element precedences = document.createElement("precedences");
        int nprec = 0;
        for (int t = 0; t < this.tasks.size(); t++) {
            nprec += this.tasks.get(t).precedences.size();
        }
        precedences.setAttribute("number", "" + nprec);
        parameters.appendChild(precedences);

        for (int t = 0; t < this.tasks.size(); t++) {
            for (int p = 0; p < this.tasks.get(t).precedences.size(); p++) {
                ActivityPrecedence P = this.tasks.get(t).precedences.get(p);
                Element precedence = document.createElement("precedence");
                if (P.preType.equals(PRE_SEQ) && P.postType.equals(POST_LOOP)) {
                    precedence.setAttribute("type", "loop");
                } else if (P.preType.equals(PRE_OR)) {
                    precedence.setAttribute("type", "or-join");
                } else if (P.preType.equals(PRE_AND)) {
                    precedence.setAttribute("type", "and-join");
                } else if (P.postType.equals(POST_OR)) {
                    precedence.setAttribute("type", "or-fork");
                } else if (P.postType.equals(POST_AND)) {
                    precedence.setAttribute("type", "and-fork");
                } else {
                    precedence.setAttribute("type", "seq");
                }
                Matrix quorum = null;
                if (P.preParams != null) {
                    if (P.preType.equals(PRE_AND) && P.postType.equals(POST_SEQ)) {
                        quorum = P.preParams;
                        for (int a = 0; a < P.preActs.size(); a++)
                            writeJLQNPrecedence(document, precedence, P.preActs.get(a), "1.0", "pre");
                    } else {
                        for (int a = 0; a < P.preActs.size(); a++)
                            writeJLQNPrecedence(document, precedence, P.preActs.get(a), "" + P.preParams.get(0), "pre");
                    }
                } else {
                    for (int a = 0; a < P.preActs.size(); a++)
                        writeJLQNPrecedence(document, precedence, P.preActs.get(a), "1.0", "pre");
                }
                boolean isLoopFirst = true;
                String loopActType;
                if (P.postParams != null) {
                    for (int a = 0; a < P.postActs.size(); a++) {
                        if (P.preType.equals(PRE_SEQ) && P.postType.equals(POST_LOOP)) {
                            if (isLoopFirst) {
                                loopActType = "post";
                                isLoopFirst = false;
                            } else {
                                loopActType = "end";
                            }
                            writeJLQNPrecedence(document, precedence, P.postActs.get(a), "" + P.postParams.get(0), loopActType);
                        } else {
                            writeJLQNPrecedence(document, precedence, P.postActs.get(a), "" + P.postParams.get(a), "post");
                        }
                    }
                } else {
                    for (int a = 0; a < P.postActs.size(); a++)
                        if (P.preType.equals(PRE_SEQ) && P.postType.equals(POST_LOOP)) {
                            if (isLoopFirst) {
                                loopActType = "post";
                                isLoopFirst = false;
                            } else {
                                loopActType = "end";
                            }
                            writeJLQNPrecedence(document, precedence, P.postActs.get(a), "1.0", loopActType);
                        } else if (P.preType.equals(PRE_AND) && P.postType.equals(POST_SEQ)) {
                            // for now JLQN quorum supports only the total number of jobs to wait
                            writeJLQNPrecedence(document, precedence, P.postActs.get(a), "" + quorum.get(0), "pre"); // Use first element to match MATLAB behavior
                        } else {
                            writeJLQNPrecedence(document, precedence, P.postActs.get(a), "1.0", "post");
                        }
                }
                precedences.appendChild(precedence);
            }
        }

        TransformerFactory transformerFactory = TransformerFactory.newInstance();
        Transformer transformer = transformerFactory.newTransformer();
        transformer.setOutputProperty(OutputKeys.INDENT, "yes");
        transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "4"); // Set indent amount to 4 spaces

        DOMSource domSource = new DOMSource(document);
        StringWriter writer = new StringWriter();
        StreamResult result = new StreamResult(writer);
        transformer.transform(domSource, result);

        File file = new File(filename);
        StreamResult fileResult = new StreamResult(file);
        transformer.transform(domSource, fileResult);
        if (GlobalConstants.Verbose != VerboseLevel.SILENT) {
            System.out.println("JLQN file saved as " + filename);
        }
    }

    /**
     * The first DIRECT child of {@code parent} with the given tag, or null.
     *
     * <p>{@code getElementsByTagName} searches every descendant, which is wrong for
     * the dialect's task- and entry-level elements: a {@code <cache>} belongs to the
     * task that declares it, not to any element nested below it.</p>
     *
     * @param parent the element to search
     * @param tag    the child tag name
     * @return the first matching direct child, or null when there is none
     */
    private static Element firstDirectChild(Element parent, String tag) {
        NodeList children = parent.getChildNodes();
        for (int i = 0; i < children.getLength(); i++) {
            org.w3c.dom.Node child = children.item(i);
            if (child instanceof Element && tag.equals(child.getNodeName())) {
                return (Element) child;
            }
        }
        return null;
    }

    /**
     * Rebuilds a distribution from a {@code mean}/{@code scv} pair, choosing the same
     * family {@code Task.setSetupTime(double)} would: exponential at unit SCV,
     * deterministic at zero, an APH fit otherwise.
     *
     * @param element an element carrying {@code mean} and optionally {@code scv}
     * @return the reconstructed distribution
     */
    private static Distribution distFromMeanAndSCV(Element element) {
        String meanStr = element.getAttribute("mean");
        double mean = meanStr.isEmpty() ? 0.0 : Double.parseDouble(meanStr);
        String scvStr = element.getAttribute("scv");
        double scv = scvStr.isEmpty() ? 1.0 : Double.parseDouble(scvStr);
        if (mean <= GlobalConstants.FineTol) {
            return Immediate.getInstance();
        }
        if (scv <= 0) {
            return new Det(mean);
        }
        if (Math.abs(scv - 1.0) <= GlobalConstants.FineTol) {
            return Exp.fitMean(mean);
        }
        return APH.fitMeanAndSCV(mean, scv);
    }

    /**
     * Writes an {@code ItemEntry}'s access popularity as the LINE .lqnx dialect
     * encodes it: the distribution CLASS NAME plus one {@code <parameter value>}
     * child per constructor argument, in constructor order, mirroring dist2json.
     *
     * <p>A vector-valued argument contributes one child per element. The reader
     * splits them using the {@code cardinality} on the parent {@code <item-entry>},
     * so a {@code DiscreteSampler} written as n parameters is p over the default
     * support 1..n, and one written as 2n parameters is p followed by an explicit
     * support x. That keeps the flat parameter list unambiguous without adding an
     * element the spec does not have.</p>
     *
     * @param doc        the document being built
     * @param itemElement the {@code <item-entry>} element to append to
     * @param popularity the popularity distribution, or null when the entry has none
     * @param entryName  entry name, for the diagnostic on an unsupported class
     */
    private static void writeAccessPopularity(Document doc, Element itemElement,
                                              Distribution popularity, String entryName) {
        if (popularity == null) {
            return;
        }
        Element popElement = doc.createElement("access-popularity");
        itemElement.appendChild(popElement);
        popElement.setAttribute("name", popularity.getName());
        List<Double> values = new ArrayList<Double>();
        if (popularity instanceof DiscreteSampler) {
            Matrix pMat = (Matrix) popularity.getParam(1).getValue();
            Matrix xMat = (Matrix) popularity.getParam(2).getValue();
            boolean defaultSupport = xMat != null && xMat.length() == pMat.length();
            for (int k = 0; defaultSupport && k < xMat.length(); k++) {
                if (Math.abs(xMat.get(k) - (k + 1)) > GlobalConstants.FineTol) {
                    defaultSupport = false;
                }
            }
            for (int k = 0; k < pMat.length(); k++) {
                values.add(pMat.get(k));
            }
            // The default support is reconstructible, so it is not written; an explicit
            // one is, and the doubled length is what tells the reader which it has.
            if (!defaultSupport) {
                for (int k = 0; k < xMat.length(); k++) {
                    values.add(xMat.get(k));
                }
            }
        } else if (popularity instanceof Zipf) {
            // Zipf stores params 1=p, 2=x, 3=s, 4=n, so its CONSTRUCTOR order is (s, n)
            // = params 3 and 4. Writing p and x instead would be lossy: the reader
            // cannot recover s from them. Always exactly two parameters, whatever the
            // cardinality, because a Zipf DERIVES p and x from (s, n).
            values.add(((Number) popularity.getParam(3).getValue()).doubleValue());
            values.add(((Number) popularity.getParam(4).getValue()).doubleValue());
        } else {
            line_error(mfilename(new Object() {}),
                    "Entry " + entryName + " has an access popularity of class "
                    + popularity.getName() + ", which the .lqnx dialect does not encode yet "
                    + "(it carries DiscreteSampler and Zipf). Writing it without its parameters "
                    + "would produce a file describing a different model, so it is refused here "
                    + "rather than silently degraded.");
            return;
        }
        for (int k = 0; k < values.size(); k++) {
            Element paramElement = doc.createElement("parameter");
            popElement.appendChild(paramElement);
            paramElement.setAttribute("value", Double.toString(values.get(k)));
        }
    }

    /**
     * Rebuilds an access popularity written by {@link #writeAccessPopularity}.
     *
     * @param itemElement the {@code <item-entry>} element
     * @param cardinality the item cardinality, used to split a flat parameter list
     * @return the distribution, or null when the element carries none
     */
    private static Distribution readAccessPopularity(Element itemElement, int cardinality) {
        NodeList popList = itemElement.getElementsByTagName("access-popularity");
        if (popList.getLength() == 0) {
            return null;
        }
        Element popElement = (Element) popList.item(0);
        String className = popElement.getAttribute("name");
        NodeList paramList = popElement.getElementsByTagName("parameter");
        List<Double> values = new ArrayList<Double>();
        for (int k = 0; k < paramList.getLength(); k++) {
            values.add(Double.parseDouble(((Element) paramList.item(k)).getAttribute("value")));
        }
        if ("Zipf".equals(className)) {
            // Always exactly (s, n); the cardinality split used for DiscreteSampler does
            // NOT apply, because a Zipf derives p and x from those two. Any other count
            // is refused by name rather than guessed at.
            if (values.size() != 2) {
                line_error(mfilename(new Object() {}),
                        "An access popularity of class Zipf carries exactly two parameters, "
                        + "the shape s then the item count n, but this one carries "
                        + values.size() + ".");
                return null;
            }
            return new Zipf(values.get(0), (int) Math.round(values.get(1)));
        }
        if ("DiscreteSampler".equals(className)) {
            if (values.isEmpty()) {
                return null;
            }
            int n = (values.size() == 2 * cardinality) ? cardinality : values.size();
            Matrix pMat = new Matrix(1, n);
            for (int k = 0; k < n; k++) {
                pMat.set(0, k, values.get(k));
            }
            if (values.size() == 2 * cardinality) {
                Matrix xMat = new Matrix(1, n);
                for (int k = 0; k < n; k++) {
                    xMat.set(0, k, values.get(n + k));
                }
                return new DiscreteSampler(pMat, xMat);
            }
            return new DiscreteSampler(pMat);
        }
        line_error(mfilename(new Object() {}),
                "Access popularity of class " + className + " is not one the .lqnx dialect "
                + "decodes (it carries DiscreteSampler and Zipf).");
        return null;
    }

    /**
     * Writes the layered network to an XML file with default naming.
     *
     * @param filename the path to write the XML file
     */
    public void writeXML(String filename) {
        writeXML(filename, false);
    }

    /**
     * Writes the layered network to an XML file.
     * 
     * @param filename the path to write the XML file
     * @param abstractNames if true, uses abstract names in the output
     */
    public void writeXML(String filename, boolean abstractNames) {
        // First pass: build a name->hash map (optionally with abstract names)
        Map<String, String> nodeHashMap = new HashMap<>();

        int tctr = 0, ectr = 0, actr = 0;

        /* hosts is a Map<Integer,Host> so keySet() is fine here */
        List<Integer> hostIds = new ArrayList<>(this.hosts.keySet());
        Collections.sort(hostIds);                      // deterministic order

        for (Integer pId : hostIds) {
            Processor curProc = (Processor) this.hosts.get(pId);
            nodeHashMap.put(curProc.getName(),
                    abstractNames ? "P" + pId : curProc.getName());

            /* -------- tasks are stored in a List inside each processor -------- */
            for (int t = 0; t < curProc.tasks.size(); t++) {
                Task curTask = curProc.tasks.get(t);
                tctr++;
                nodeHashMap.put(curTask.getName(),
                        abstractNames ? "T" + tctr : curTask.getName());

                /* entries list */
                for (int e = 0; e < curTask.entries.size(); e++) {
                    Entry curEntry = curTask.entries.get(e);
                    ectr++;
                    nodeHashMap.put(curEntry.getName(),
                            abstractNames ? "E" + ectr : curEntry.getName());
                }

                /* activities list */
                for (int a = 0; a < curTask.activities.size(); a++) {
                    Activity curAct = curTask.activities.get(a);
                    actr++;
                    nodeHashMap.put(curAct.getName(),
                            abstractNames ? "A" + actr : curAct.getName());
                }
            }
        }

        /* ------------------------------------------------------------------
         * XML document root
         * ------------------------------------------------------------------ */
        DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder docBuilder;
        try {
            docBuilder = docFactory.newDocumentBuilder();
        } catch (ParserConfigurationException e) {
            throw new RuntimeException(e);
        }
        Document doc = docBuilder.newDocument();
        Element rootElement = doc.createElement("lqn-model");
        doc.appendChild(rootElement);
        rootElement.setAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
        rootElement.setAttribute("xsi:noNamespaceSchemaLocation", "lqn.xsd");
        rootElement.setAttribute("name", this.getName());

        // rendezvous-entry gating for reply-entry emission: see _kb/06-solver-catalog.md ("Activity-graph validity and .lqnx writer rules")
        Set<String> rendezvousEntries = new HashSet<>();
        for (Activity act : this.activities.values()) {
            rendezvousEntries.addAll(act.syncCallDests.values());
        }
        for (Entry ent : this.entries.values()) {
            rendezvousEntries.addAll(ent.getForwardingDests().values());
        }

        // Second pass: emit processors, tasks, entries, activities
        for (Integer pId : hostIds) {
            Processor curProc = (Processor) this.hosts.get(pId);
            Element procElement = doc.createElement("processor");
            rootElement.appendChild(procElement);
            procElement.setAttribute("name", nodeHashMap.get(curProc.getName()));
            procElement.setAttribute("scheduling", lqnSchedText(curProc.scheduling));
            if (curProc.replication > 1) {
                procElement.setAttribute("replication", Integer.toString(curProc.replication));
            }
            if (!curProc.scheduling.equals(SchedStrategy.INF)) {
                procElement.setAttribute("multiplicity", Integer.toString(curProc.multiplicity));
            }
            if (curProc.scheduling.equals(SchedStrategy.PS) ||
                    curProc.scheduling.equals(SchedStrategy.PSPRIO) ||
                    curProc.scheduling.equals(SchedStrategy.GPS)) {
                // PS, PSPRIO and CFS (GPS) processors require a scheduling quantum.
                procElement.setAttribute("quantum", Double.toString(curProc.quantum));
            }
            procElement.setAttribute("speed-factor", Double.toString(curProc.speedFactor));

            /* ----------------------------- TASKS -------------------------- */
            for (int t = 0; t < curProc.tasks.size(); t++) {
                Task curTask = curProc.tasks.get(t);
                Element taskElement = doc.createElement("task");
                procElement.appendChild(taskElement);
                taskElement.setAttribute("name", nodeHashMap.get(curTask.getName()));
                taskElement.setAttribute("scheduling", lqnSchedText(curTask.scheduling));
                if (curTask.replication > 1) {
                    taskElement.setAttribute("replication", Integer.toString(curTask.replication));
                }
                if (!curTask.scheduling.equals(SchedStrategy.INF)) {
                    taskElement.setAttribute("multiplicity", Integer.toString(curTask.multiplicity));
                }
                // see _kb/11-conventions-and-gotchas.md (.lqnx cannot carry a non-reference task's think time) for rationale
                if (curTask.scheduling.equals(SchedStrategy.REF)) {
                    taskElement.setAttribute("think-time", Double.toString(curTask.thinkTimeMean));
                } else if (curTask.thinkTimeMean > 0) {
                    line_warning(mfilename(new Object(){}),
                            "Task %s is not a reference task, so its think time (%g) is not written: "
                            + "the LQN XML schema accepts think-time on reference tasks only.",
                            curTask.getName(), curTask.thinkTimeMean);
                }

                /* fan-out/fan-in are parsed by this reader, by MATLAB parseXML and by the cpp
                 * lqn_reader, and were written by no codebase, so a replicated model lost its
                 * call multiplicities on every round trip. lqn-core.xsd (TaskType) places them
                 * before the entries. */
                for (Map.Entry<String, Integer> fo : curTask.getFanOutMap().entrySet()) {
                    Element fanOutElement = doc.createElement("fan-out");
                    taskElement.appendChild(fanOutElement);
                    fanOutElement.setAttribute("dest", nodeHashMap.getOrDefault(fo.getKey(), fo.getKey()));
                    fanOutElement.setAttribute("value", Integer.toString(fo.getValue()));
                }
                if (!curTask.getFanInSource().isEmpty() && curTask.getFanInValue() > 0) {
                    Element fanInElement = doc.createElement("fan-in");
                    taskElement.appendChild(fanInElement);
                    fanInElement.setAttribute("source",
                            nodeHashMap.getOrDefault(curTask.getFanInSource(), curTask.getFanInSource()));
                    fanInElement.setAttribute("value", Integer.toString(curTask.getFanInValue()));
                }

                /* LINE .lqnx dialect: a CacheTask and a task's setup/delay-off times. The
                 * base schema carries neither, so writing a plain <task> emitted a VALID
                 * LQN of a DIFFERENT model -- lcq_threehosts then read a cache hit ratio of
                 * exactly 0.5 (the unweighted POST_CACHE branch) and lqn_setup a processor
                 * utilization of 0.75 against a true 0.43577. See LQNX_CACHE_SPEC. */
                if (curTask instanceof CacheTask) {
                    CacheTask cacheTask = (CacheTask) curTask;
                    Element cacheElement = doc.createElement("cache");
                    taskElement.appendChild(cacheElement);
                    cacheElement.setAttribute("items", Integer.toString(cacheTask.getItems()));
                    cacheElement.setAttribute("replacement", cacheTask.getReplacestrategy().name());
                    cacheElement.setAttribute("retrieval", Boolean.toString(cacheTask.hasRetrieval()));
                    // itemLevelCap is an ARRAY: one <level> per cache list, in order. Never
                    // collapsed to a scalar, because a multi-list cache is the normal case.
                    int[] levelCaps = cacheTask.getItemLevelCap();
                    for (int lv = 0; lv < levelCaps.length; lv++) {
                        Element levelElement = doc.createElement("level");
                        cacheElement.appendChild(levelElement);
                        levelElement.setAttribute("capacity", Integer.toString(levelCaps[lv]));
                    }
                }
                if (curTask.getSetupTime() != null && curTask.getSetupTimeMean() > GlobalConstants.FineTol) {
                    Element setupElement = doc.createElement("setup");
                    taskElement.appendChild(setupElement);
                    setupElement.setAttribute("mean", Double.toString(curTask.getSetupTimeMean()));
                    setupElement.setAttribute("scv", Double.toString(curTask.getSetupTimeSCV()));
                }
                if (curTask.getDelayOffTime() != null && curTask.getDelayOffTimeMean() > GlobalConstants.FineTol) {
                    Element delayOffElement = doc.createElement("delay-off");
                    taskElement.appendChild(delayOffElement);
                    delayOffElement.setAttribute("mean", Double.toString(curTask.getDelayOffTimeMean()));
                    delayOffElement.setAttribute("scv", Double.toString(curTask.getDelayOffTimeSCV()));
                }

                /* Track activities that are exported as entry-phase-activities (not to be duplicated in task-activities) */
                Set<String> phaseActivityNames = new HashSet<>();

                /* Build set of activities that participate in precedences - these MUST be in task-activities */
                Set<String> activitiesInPrecedences = new HashSet<>();
                for (ActivityPrecedence prec : curTask.precedences) {
                    activitiesInPrecedences.addAll(prec.preActs);
                    activitiesInPrecedences.addAll(prec.postActs);
                }

                /* Track entries whose type was overridden to NONE due to precedence activities */
                Set<String> overriddenToNoneEntries = new HashSet<>();

                /* ----------- ENTRIES -------------- */
                for (int e = 0; e < curTask.entries.size(); e++) {
                    Entry curEntry = curTask.entries.get(e);
                    Element entryElement = doc.createElement("entry");
                    taskElement.appendChild(entryElement);
                    entryElement.setAttribute("name", nodeHashMap.get(curEntry.getName()));

                    // open-arrival-rate emission rationale: see _kb/06-solver-catalog.md ("Activity-graph validity and .lqnx writer rules")
                    if (curEntry.getArrival() != null) {
                        double arrMean = curEntry.getArrival().getMean();
                        if (arrMean > 0 && !Double.isInfinite(arrMean)) {
                            entryElement.setAttribute("open-arrival-rate", Double.toString(1.0 / arrMean));
                        }
                    }

                    /* LINE .lqnx dialect: the presence of <item-entry> is what makes this an
                     * ItemEntry on read. See LQNX_CACHE_SPEC. */
                    if (curEntry instanceof ItemEntry) {
                        ItemEntry itemEntry = (ItemEntry) curEntry;
                        Element itemElement = doc.createElement("item-entry");
                        entryElement.appendChild(itemElement);
                        itemElement.setAttribute("cardinality", Integer.toString(itemEntry.getCardinality()));
                        writeAccessPopularity(doc, itemElement, itemEntry.getPopularity(), curEntry.getName());
                    }

                    /* Get entry type from the entry itself */
                    String entryType = curEntry.getType();
                    if (entryType == null || entryType.isEmpty()) {
                        entryType = "NONE";
                    }

                    // PH1PH2 collapse-to-NONE rule: see _kb/06-solver-catalog.md ("Activity-graph validity and .lqnx writer rules")
                    if (entryType.equals("PH1PH2")) {
                        boolean hasMultiPhase = false;
                        for (Activity act : curTask.activities) {
                            if (act.getPhase() > 1) {
                                hasMultiPhase = true;
                                break;
                            }
                        }
                        if (!hasMultiPhase) {
                            /* Check if any activity is bound to this entry */
                            boolean hasBoundActivities = false;
                            for (Activity act : curTask.activities) {
                                if (act.boundToEntry.equals(curEntry.getName())) {
                                    hasBoundActivities = true;
                                    break;
                                }
                            }
                            if (hasBoundActivities && !activitiesInPrecedences.isEmpty()) {
                                entryType = "NONE";
                                overriddenToNoneEntries.add(curEntry.getName());
                            }
                        }
                    }
                    entryElement.setAttribute("type", entryType);

                    // Collect phase activities for PH1PH2 entries: phase-1 is bound to the entry, later phases found via the phase-sequence precedence chain
                    List<Activity> phaseActivities = new ArrayList<>();
                    if (entryType.equals("PH1PH2")) {
                        /* Find the phase-1 activity bound to this entry */
                        Activity phase1Act = null;
                        for (Activity act : curTask.activities) {
                            if (act.boundToEntry.equals(curEntry.getName())) {
                                phase1Act = act;
                                phaseActivities.add(act);
                                phaseActivityNames.add(act.getName());
                                break;
                            }
                        }
                        /* Walk precedence chain to find subsequent phases */
                        if (phase1Act != null) {
                            String currentName = phase1Act.getName();
                            boolean found = true;
                            while (found) {
                                found = false;
                                for (ActivityPrecedence prec : curTask.precedences) {
                                    if (prec.preActs.size() == 1 && prec.postActs.size() == 1 &&
                                        prec.preActs.get(0).equals(currentName)) {
                                        String nextName = prec.postActs.get(0);
                                        /* Find the activity object */
                                        for (Activity act : curTask.activities) {
                                            if (act.getName().equals(nextName) && act.getPhase() > 0) {
                                                phaseActivities.add(act);
                                                phaseActivityNames.add(act.getName());
                                                currentName = nextName;
                                                found = true;
                                                break;
                                            }
                                        }
                                        break;
                                    }
                                }
                            }
                        }
                    }

                    /* Export entry-phase-activities for PH1PH2 entries */
                    if (!phaseActivities.isEmpty()) {
                        Element entryPhaseActsElement = doc.createElement("entry-phase-activities");
                        entryElement.appendChild(entryPhaseActsElement);

                        for (Activity phaseAct : phaseActivities) {
                            Element phaseActElement = doc.createElement("activity");
                            entryPhaseActsElement.appendChild(phaseActElement);
                            phaseActElement.setAttribute("name", nodeHashMap.get(phaseAct.getName()));
                            phaseActElement.setAttribute("phase", Integer.toString(phaseAct.getPhase()));
                            phaseActElement.setAttribute("host-demand-mean", Double.toString(phaseAct.hostDemandMean));
                            if (phaseAct.hostDemandSCV > 0) {
                                phaseActElement.setAttribute("host-demand-cvsq", Double.toString(phaseAct.hostDemandSCV));
                            }
                            if (phaseAct.thinkTimeMean > GlobalConstants.Zero) {
                                phaseActElement.setAttribute("think-time", Double.toString(phaseAct.thinkTimeMean));
                            }

                            /* synchronous calls for phase activity */
                            for (int sc = 0; sc < phaseAct.syncCallDests.size(); sc++) {
                                Element syncCallElement = doc.createElement("synch-call");
                                phaseActElement.appendChild(syncCallElement);
                                syncCallElement.setAttribute("dest", nodeHashMap.get(phaseAct.syncCallDests.get(sc)));
                                syncCallElement.setAttribute("calls-mean", Double.toString(phaseAct.syncCallMeans.get(sc)));
                            }

                            /* asynchronous calls for phase activity */
                            for (int ac = 0; ac < phaseAct.asyncCallDests.size(); ac++) {
                                Element asyncCallElement = doc.createElement("asynch-call");
                                phaseActElement.appendChild(asyncCallElement);
                                asyncCallElement.setAttribute("dest", nodeHashMap.get(phaseAct.asyncCallDests.get(ac)));
                                asyncCallElement.setAttribute("calls-mean", Double.toString(phaseAct.asyncCallMeans.get(ac)));
                            }

                            writeCallGroups(doc, phaseActElement, phaseAct, nodeHashMap);
                        }
                    }

                    /* forwarding calls */
                    for (int fw = 0; fw < curEntry.getForwardingDests().size(); fw++) {
                        Element fwdElement = doc.createElement("forwarding");
                        entryElement.appendChild(fwdElement);
                        fwdElement.setAttribute("dest", nodeHashMap.get(curEntry.getForwardingDests().get(fw)));
                        fwdElement.setAttribute("prob", Double.toString(curEntry.getForwardingProbs().get(fw)));
                    }
                }

                /* ----------- ACTIVITIES ----------- */
                /* Check if any activities will go into task-activities (not already in entry-phase-activities) */
                boolean hasTaskActivities = false;
                for (int a = 0; a < curTask.activities.size(); a++) {
                    if (!phaseActivityNames.contains(curTask.activities.get(a).getName())) {
                        hasTaskActivities = true;
                        break;
                    }
                }
                /* Skip empty task-activities element for pure PH1PH2 tasks */
                if (!hasTaskActivities && curTask.precedences.isEmpty()) {
                    continue;  // No task-level activities, precedences, or reply-entries needed
                }

                Element taskActsElement = doc.createElement("task-activities");
                taskElement.appendChild(taskActsElement);

                for (int a = 0; a < curTask.activities.size(); a++) {
                    Activity curAct = curTask.activities.get(a);

                    /* Skip activities that were already exported in entry-phase-activities */
                    if (phaseActivityNames.contains(curAct.getName())) {
                        continue;
                    }

                    Element actElement = doc.createElement("activity");
                    taskActsElement.appendChild(actElement);
                    actElement.setAttribute("host-demand-mean", Double.toString(curAct.hostDemandMean));
                    actElement.setAttribute("host-demand-cvsq", Double.toString(curAct.hostDemandSCV));
                    if (!curAct.boundToEntry.isEmpty()) {
                        actElement.setAttribute("bound-to-entry", nodeHashMap.get(curAct.boundToEntry));
                    }
                    actElement.setAttribute("call-order", curAct.callOrder);
                    actElement.setAttribute("name", nodeHashMap.get(curAct.getName()));
                    if (curAct.thinkTimeMean > GlobalConstants.Zero) {
                        actElement.setAttribute("think-time", Double.toString(curAct.thinkTimeMean));
                    }

                    /* synchronous calls */
                    for (int sc = 0; sc < curAct.syncCallDests.size(); sc++) {
                        Element syncCallElement = doc.createElement("synch-call");
                        actElement.appendChild(syncCallElement);
                        syncCallElement.setAttribute("dest", nodeHashMap.get(curAct.syncCallDests.get(sc)));
                        syncCallElement.setAttribute("calls-mean", Double.toString(curAct.syncCallMeans.get(sc)));
                    }

                    /* asynchronous calls */
                    for (int ac = 0; ac < curAct.asyncCallDests.size(); ac++) {
                        Element asyncCallElement = doc.createElement("asynch-call");
                        actElement.appendChild(asyncCallElement);
                        asyncCallElement.setAttribute("dest", nodeHashMap.get(curAct.asyncCallDests.get(ac)));
                        asyncCallElement.setAttribute("calls-mean", Double.toString(curAct.asyncCallMeans.get(ac)));
                    }

                    writeCallGroups(doc, actElement, curAct, nodeHashMap);
                }

                /* ----------- PRECEDENCES ---------- */
                for (int ap = 0; ap < curTask.precedences.size(); ap++) {
                    ActivityPrecedence curActPrec = curTask.precedences.get(ap);

                    /* Skip precedences that only involve phase activities (PH1PH2 entries handle sequencing internally) */
                    boolean allPhase = true;
                    for (String act : curActPrec.preActs) {
                        if (!phaseActivityNames.contains(act)) { allPhase = false; break; }
                    }
                    if (allPhase) {
                        for (String act : curActPrec.postActs) {
                            if (!phaseActivityNames.contains(act)) { allPhase = false; break; }
                        }
                    }
                    if (allPhase) continue;

                    Element actPrecElement = doc.createElement("precedence");
                    taskActsElement.appendChild(actPrecElement);

                    Element preElement = doc.createElement(curActPrec.preType);
                    actPrecElement.appendChild(preElement);
                    if (curActPrec.preType.equals(PRE_AND)) {
                        // Emit the quorum attribute only for a genuine quorum, i.e. k < n. A join
                        // that waits for all predecessors is the LQN default and carries no attribute.
                        int quorum = ActivityPrecedence.getQuorumCount(curActPrec.preParams, curActPrec.preActs.size());
                        if (quorum > 0 && quorum != curActPrec.preActs.size()) {
                            preElement.setAttribute("quorum", Integer.toString(quorum));
                        }
                    }
                    for (String pra : curActPrec.preActs) {
                        Element preActElement = doc.createElement("activity");
                        preElement.appendChild(preActElement);
                        preActElement.setAttribute("name", nodeHashMap.get(pra));
                    }

                    Element postElement = doc.createElement(curActPrec.postType);
                    actPrecElement.appendChild(postElement);
                    if (curActPrec.postType.equals(POST_OR)) {
                        for (int i = 0; i < curActPrec.postActs.size(); i++) {
                            Element postActElement = doc.createElement("activity");
                            postElement.appendChild(postActElement);
                            postActElement.setAttribute("name", nodeHashMap.get(curActPrec.postActs.get(i)));
                            postActElement.setAttribute("prob", Double.toString(curActPrec.postParams.get(i)));
                        }
                    } else if (curActPrec.postType.equals(POST_LOOP)) {
                        for (int i = 0; i < curActPrec.postActs.size() - 1; i++) {
                            Element postActElement = doc.createElement("activity");
                            postElement.appendChild(postActElement);
                            postActElement.setAttribute("name", nodeHashMap.get(curActPrec.postActs.get(i)));
                            postActElement.setAttribute("count", Double.toString(curActPrec.postParams.get(i)));
                        }
                        postElement.setAttribute("end",
                                curActPrec.postActs.get(curActPrec.postActs.size() - 1));
                    } else if (curActPrec.postType.equals(ActivityPrecedenceType.POST_CACHE)) {
                        // cache-result is explicit rather than positional, so a reader never has
                        // to infer the hit branch from document order. Readers still fall back to
                        // that order when the attribute is absent, so older files keep loading.
                        for (int i = 0; i < curActPrec.postActs.size(); i++) {
                            Element postActElement = doc.createElement("activity");
                            postElement.appendChild(postActElement);
                            postActElement.setAttribute("name", nodeHashMap.get(curActPrec.postActs.get(i)));
                            postActElement.setAttribute("cache-result", i == 0 ? "hit" : "miss");
                        }
                    } else {
                        for (String poa : curActPrec.postActs) {
                            Element postActElement = doc.createElement("activity");
                            postElement.appendChild(postActElement);
                            postActElement.setAttribute("name", nodeHashMap.get(poa));
                        }
                    }
                }

                /* ----------- REPLY-ACTIVITIES ----- */
                /* Build a set of activity names in this task for validation */
                Set<String> taskActivityNames = new HashSet<>();
                for (Activity act : curTask.activities) {
                    taskActivityNames.add(act.getName());
                }

                for (Entry curEntry : curTask.entries) {
                    /* Skip reference tasks - they don't reply */
                    if (curTask.scheduling.equals(SchedStrategy.REF)) {
                        continue;
                    }

                    /* Skip entries that do not accept rendezvous (only async-called or
                       uncalled) - emitting a reply for them is an invalid LQN model. */
                    if (!rendezvousEntries.contains(curEntry.getName())) {
                        continue;
                    }

                    /* Skip PH1PH2 entries - they have implicit replies and don't need reply-entry elements
                       But process entries that were overridden to NONE due to precedence activities */
                    String entryType = curEntry.getType();
                    if (entryType != null && entryType.equals("PH1PH2") &&
                        !overriddenToNoneEntries.contains(curEntry.getName())) {
                        continue;
                    }

                    /* First, clear any cross-task reply activities that may have been incorrectly set */
                    Map<Integer, String> validReplyMap = new HashMap<>();
                    for (Map.Entry<Integer, String> raEntry : curEntry.replyActivity.entrySet()) {
                        if (taskActivityNames.contains(raEntry.getValue()) && !phaseActivityNames.contains(raEntry.getValue())) {
                            validReplyMap.put(raEntry.getKey(), raEntry.getValue());
                        }
                    }
                    curEntry.replyActivity.clear();
                    curEntry.replyActivity.putAll(validReplyMap);

                    /* If the entry has no reply-activities yet, infer them from lsn.replygraph */
                    if (curEntry.replyActivity.isEmpty()) {
                        LayeredNetworkStruct lsn = this.getStruct();

                        /* locate the column of this entry inside lsn.replygraph   */
                        int eidx = -1;
                        for (Map.Entry<Integer, Entry> kv : this.entries.entrySet()) {
                            if (kv.getValue() == curEntry) {
                                eidx = kv.getKey();                  // entry-local index, as replygraph columns are
                                break;
                            }
                        }

                        /* add every activity that replies to this entry, but only if it belongs to this task and is not a phase activity */
                        if (eidx >= 0 && eidx < lsn.replygraph.getNumCols()) {
                            for (int ra = 0; ra < lsn.replygraph.getNumRows(); ra++) {
                                if (lsn.replygraph.get(ra, eidx) != 0) {
                                    String actName = lsn.names.get(lsn.ashift + ra); // replygraph rows are activity-local, matching names at ashift+ra
                                    // Only add if the activity belongs to this task and is not already in entry-phase-activities
                                    if (taskActivityNames.contains(actName) && !phaseActivityNames.contains(actName)) {
                                        Integer actIndex = findActivityIndexByName(actName);
                                        if (actIndex != null) {
                                            curEntry.replyActivity.put(actIndex, actName);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    /* Emit <reply-entry> only if at least one valid reply-activity exists */
                    if (!curEntry.replyActivity.isEmpty()) {
                        Element entryReplyElement = doc.createElement("reply-entry");
                        taskActsElement.appendChild(entryReplyElement);
                        entryReplyElement.setAttribute("name",
                                nodeHashMap.get(curEntry.getName()));

                        for (String ra : curEntry.replyActivity.values()) {
                            Element entryReplyActElement = doc.createElement("reply-activity");
                            entryReplyElement.appendChild(entryReplyActElement);
                            entryReplyActElement.setAttribute("name", nodeHashMap.get(ra));
                        }
                    }
                }
            }
        }

        // write XML to file
        try {
            TransformerFactory transformerFactory = TransformerFactory.newInstance();
            Transformer transformer = transformerFactory.newTransformer();
            transformer.setOutputProperty(OutputKeys.INDENT, "yes");
            transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");
            transformer.setOutputProperty(OutputKeys.VERSION, "1.0");
            transformer.setOutputProperty(OutputKeys.ENCODING, "UTF-8");
            transformer.setOutputProperty(OutputKeys.STANDALONE, "no");
            DOMSource source = new DOMSource(doc);
            StreamResult result = new StreamResult(new File(filename));
            transformer.transform(source, result);
        } catch (TransformerException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * LINE dialect &lt;call-group&gt;: which of the synch-calls written above one
     * dispatcher issues, and under which strategy. The member calls stay ordinary
     * synch-calls, so a reader that ignores this element still sees the same
     * aggregate call means -- which is what lqns and lqsim, having no dispatcher,
     * should see.
     */
    private static void writeCallGroups(Document doc, Element actElement, Activity act,
                                        Map<String, String> nodeHashMap) {
        for (Activity.CallGroup grp : act.getSyncCallGroups()) {
            Element grpElement = doc.createElement("call-group");
            actElement.appendChild(grpElement);
            grpElement.setAttribute("strategy", callGroupStrategyName(grp.strategy));
            for (String dest : grp.dests) {
                Element destElement = doc.createElement("dest");
                grpElement.appendChild(destElement);
                destElement.setAttribute("name", nodeHashMap.get(dest));
            }
        }
    }

    /**
     * RoutingStrategy -&gt; the wire enum name, spelled as the JSON interchange spells
     * it. Only the two strategies a call group can be built with are named: WRROBIN
     * would need per-target weights the group API does not take, and the remaining
     * strategies are not dispatch policies at all, so an unnamed one is an error
     * rather than a silent PROB.
     */
    private static String callGroupStrategyName(RoutingStrategy strategy) {
        if (strategy == RoutingStrategy.RROBIN) {
            return "RROBIN";
        }
        if (strategy == RoutingStrategy.JSQ) {
            return "JSQ";
        }
        throw new IllegalArgumentException("Call groups carry RROBIN or JSQ; routing strategy "
                + strategy + " cannot be written to .lqnx");
    }

    /**
     * Wire enum name -&gt; RoutingStrategy, the inverse of callGroupStrategyName.
     */
    private static RoutingStrategy callGroupStrategyOf(String name, String actName) {
        String key = name == null ? "" : name.trim().toUpperCase();
        if (key.equals("RROBIN")) {
            return RoutingStrategy.RROBIN;
        }
        if (key.equals("JSQ")) {
            return RoutingStrategy.JSQ;
        }
        throw new IllegalArgumentException("Activity \"" + actName + "\" declares a call group with"
                + " an unrecognized strategy \"" + name + "\"; the dialect spells them RROBIN and JSQ");
    }

    /**
     * Reads the LINE dialect &lt;call-group&gt; children of an activity element into ACT.
     *
     * The member calls are ordinary synch-call elements and have already been read,
     * so only the grouping is recorded; issuing them again would double the call rate.
     */
    private static void parseCallGroups(Element actElement, Activity act) {
        NodeList children = actElement.getChildNodes();
        for (int i = 0; i < children.getLength(); i++) {
            if (!(children.item(i) instanceof Element)) {
                continue;
            }
            Element grpElement = (Element) children.item(i);
            if (!"call-group".equals(grpElement.getNodeName())) {
                continue;
            }
            RoutingStrategy strategy =
                    callGroupStrategyOf(grpElement.getAttribute("strategy"), act.getName());
            List<String> dests = new ArrayList<>();
            NodeList destNodes = grpElement.getChildNodes();
            for (int j = 0; j < destNodes.getLength(); j++) {
                if (!(destNodes.item(j) instanceof Element)) {
                    continue;
                }
                Element destElement = (Element) destNodes.item(j);
                if ("dest".equals(destElement.getNodeName())) {
                    dests.add(destElement.getAttribute("name"));
                }
            }
            act.recordCallGroup(strategy, dests);
        }
    }

    /**
     * LQN-valid scheduling name for a processor/task in writeXML. LINE maps the LQN
     * "cfs" (completely fair scheduling) discipline onto GPS; write it back as "cfs"
     * so lqns/lqsim accept the round-tripped model (toText(GPS)="gps" is not valid LQN).
     * The same holds for "pri", the lqns spelling of preemptive priority resume,
     * which LINE holds as FCFSPRPRIO (toText gives "fcfsprprio", not valid LQN).
     */
    private static String lqnSchedText(SchedStrategy scheduling) {
        if (scheduling == SchedStrategy.GPS) {
            return "cfs";
        }
        if (scheduling == SchedStrategy.FCFSPRPRIO) {
            return "pri";
        }
        return SchedStrategy.toText(scheduling);
    }

    /**
     * Initialize the used features for each model in the ensemble
     */
    public void initUsedFeatures() {
        if (this.ensemble.isEmpty()) {
            getEnsemble();
        }
        // Initialize usedFeatures for each model in the ensemble
        // This corresponds to the MATLAB implementation
        this.usedFeatures = new FeatureSet();
    }

    /**
     * Set a language feature as used
     * @param feature the name of the feature to mark as used
     */
    public void setUsedFeatures(String feature) {
        if (this.usedFeatures == null) {
            this.usedFeatures = new FeatureSet();
        }
        this.usedFeatures.setTrue(feature);
    }

    /**
     * Get the used language features by analyzing the layered network structure
     * @return FeatureSet containing all features used in this layered network
     */
    /**
     * Marks a host demand or think time under its registry name.
     *
     * Reads getFeatureName, as Network.getUsedLangFeatures does, and applies the
     * same MarkedMAP / MarkedMMPP to MMAP normalization. The previous
     * getClass().getSimpleName() was the JAVA CLASS name, which coincides with
     * the registry name for most distributions and not for all of them -- and
     * setTrue line_errors on a name it does not know, so an LQN whose host
     * demand was one of the others crashed the feature scan instead of being
     * gated by it.
     *
     * @param distribution the host demand or think time to mark
     */
    private void markDistribution(Distribution distribution) {
        String feature = distribution.getFeatureName();
        if ("MarkedMAP".equals(feature) || "MarkedMMPP".equals(feature)) {
            feature = "MMAP";
        }
        this.usedFeatures.setTrue(feature);
    }

    /**
     * Marks one side of an activity precedence, skipping the types that gate
     * nothing (POST_LOOP, and a null side on a precedence built without one).
     *
     * @param precedenceType the preType or postType string
     */
    private void markPrecedence(String precedenceType) {
        if (precedenceType == null) {
            return;
        }
        String feature = ActivityPrecedenceType.toFeature(precedenceType);
        if (feature.length() > 0) {
            this.usedFeatures.setTrue(feature);
        }
    }

    public FeatureSet getUsedLangFeatures() {
        if (this.usedFeatures == null) {
            this.usedFeatures = new FeatureSet();
        }

        // Analyze hosts (processors)
        for (Host host : this.hosts.values()) {
            // Host and Processor are the two registry spellings of the same
            // construct, and getLNFeatureSet declares both.
            this.usedFeatures.setTrue("Host");
            this.usedFeatures.setTrue("Processor");
            // Mark processor scheduling features
            switch (host.scheduling) {
                case INF:
                    this.usedFeatures.setTrue("SchedStrategy_INF");
                    break;
                case FCFS:
                    this.usedFeatures.setTrue("SchedStrategy_FCFS");
                    break;
                case PS:
                    this.usedFeatures.setTrue("SchedStrategy_PS");
                    break;
                case DPS:
                    this.usedFeatures.setTrue("SchedStrategy_DPS");
                    break;
                case GPS:
                    this.usedFeatures.setTrue("SchedStrategy_GPS");
                    break;
                case PSPRIO:
                    this.usedFeatures.setTrue("SchedStrategy_PSPRIO");
                    break;
                case DPSPRIO:
                    this.usedFeatures.setTrue("SchedStrategy_DPSPRIO");
                    break;
                case GPSPRIO:
                    this.usedFeatures.setTrue("SchedStrategy_GPSPRIO");
                    break;
                case POLLING:
                    this.usedFeatures.setTrue("SchedStrategy_POLLING");
                    break;
                case LCFS:
                    this.usedFeatures.setTrue("SchedStrategy_LCFS");
                    break;
                case SIRO:
                    this.usedFeatures.setTrue("SchedStrategy_SIRO");
                    break;
                case HOL:
                    this.usedFeatures.setTrue("SchedStrategy_HOL");
                    break;
                default:
                    break;
            }
        }

        // Analyze tasks
        for (Task task : this.tasks.values()) {
            this.usedFeatures.setTrue("Task");
            // Mark task scheduling features
            switch (task.scheduling) {
                case REF:
                    // A reference task IS a scheduling discipline, and the only
                    // one every LQN carries, so leaving it unmarked made the
                    // SchedStrategy_REF entry unreachable. SolverLDES, the sole
                    // consumer of this set (getLNFeatureSet), already declares it.
                    this.usedFeatures.setTrue("SchedStrategy_REF");
                    break;
                case INF:
                    this.usedFeatures.setTrue("SchedStrategy_INF");
                    break;
                case FCFS:
                    this.usedFeatures.setTrue("SchedStrategy_FCFS");
                    break;
                case LCFS:
                    this.usedFeatures.setTrue("SchedStrategy_LCFS");
                    break;
                case SIRO:
                    this.usedFeatures.setTrue("SchedStrategy_SIRO");
                    break;
                case HOL:
                    this.usedFeatures.setTrue("SchedStrategy_HOL");
                    break;
                case PS:
                    // Recorded for completeness. The featset carries one
                    // SchedStrategy_PS name for hosts and tasks alike, so a PS task
                    // passes the gate and is rejected imperatively by the engine
                    // that cannot serve it (a task holds threads, it does not
                    // divide them).
                    this.usedFeatures.setTrue("SchedStrategy_PS");
                    break;
                default:
                    break;
            }

            // Analyze think time distribution
            if (task.thinkTime != null && !(task.thinkTime instanceof Immediate)) {
                markDistribution(task.thinkTime);
            }

            // Check if it's a cache task
            if (task instanceof CacheTask) {
                this.usedFeatures.setTrue("Cache");
                this.usedFeatures.setTrue("CacheTask");
                this.usedFeatures.setTrue(
                        ReplacementStrategy.toFeature(((CacheTask) task).replacestrategy));
            }
        }

        // Analyze activities
        for (Activity activity : this.activities.values()) {
            this.usedFeatures.setTrue("Activity");
            // Analyze host demand distribution
            if (activity.hostDemand != null && !(activity.hostDemand instanceof Immediate)) {
                markDistribution(activity.hostDemand);
            }

            // Mark call features
            if (!activity.syncCallDests.isEmpty()) {
                this.usedFeatures.setTrue("SyncCall");
            }
            if (!activity.asyncCallDests.isEmpty()) {
                this.usedFeatures.setTrue("AsyncCall");
            }
        }

        // Analyze entries
        for (Entry entry : this.entries.values()) {
            this.usedFeatures.setTrue("Entry");
            if (entry instanceof ItemEntry) {
                // Item-based entries for cache modeling
                this.usedFeatures.setTrue("Cache");
                this.usedFeatures.setTrue("ItemEntry");
            }
        }

        // Analyze activity precedences. The six ActivityPrecedence_ entries are
        // declared by getLNFeatureSet and were emitted by nothing, so an OR-fork
        // reached a solver that does not implement one exactly as a sequence did.
        for (Task task : this.tasks.values()) {
            for (ActivityPrecedence precedence : task.getPrecedences()) {
                markPrecedence(precedence.getPreType());
                markPrecedence(precedence.getPostType());
            }
        }

        return this.usedFeatures;
    }

    /**
     * Validates the LayeredNetwork configuration.
     *
     * Ensures that if entries are defined, activities are also defined to serve those entries.
     * An entry without a bound activity is not a functional LQN model.
     *
     * @throws RuntimeException if entries exist but no activities are defined
     */
    public void sanitize() {
        int numEntries = this.entries.size();
        int numActivities = this.activities.size();

        if (numEntries > 0 && numActivities == 0) {
            line_error(mfilename(new Object() {}),
                String.format("LayeredNetwork '%s' has %d entry(ies) but no activities. " +
                    "Entries must be bound to activities to form a valid LQN model. " +
                    "Use activity.boundTo(entry) to establish the binding.",
                    this.getName(), numEntries));
        }
    }

    private class Param {

        public Nodes Nodes;
        public Edges Edges;

        public Param() {
            this.Nodes = new Nodes();
            this.Edges = new Edges();
        }
    }

    private class Nodes {
        double QLen;
        double RespT;
        double Tput;
        double Util;
    }

    private class Edges {
        double QLen;
        double RespT;
        double Tput;
    }

    /**
     * Helper method to find an activity's index by name
     * @param activityName The name of the activity to find
     * @return The index of the activity in the model, or null if not found
     */
    private Integer findActivityIndexByName(String activityName) {
        for (Map.Entry<Integer, Activity> entry : this.activities.entrySet()) {
            if (entry.getValue().getName().equals(activityName)) {
                return entry.getKey();
            }
        }
        return null;
    }

    /**
     * Per-operand peak rate scaling of a class- or joint-dependence declaration,
     * broadcasting a 1x1 declaration onto the server's operands.
     */
    private static Matrix expandPeak(Matrix peak, int ncols, String elemname, String colwhat, String what) {
        if (peak.length() == 1) {
            Matrix out = new Matrix(1, ncols);
            for (int j = 0; j < ncols; j++) {
                out.set(0, j, peak.get(0));
            }
            return out;
        }
        if (peak.length() != ncols) {
            throw new IllegalArgumentException(what + "-dependence peak rate on " + elemname + " has "
                    + peak.length() + " entries but there are " + ncols + " " + colwhat + ".");
        }
        Matrix out = new Matrix(1, ncols);
        for (int j = 0; j < ncols; j++) {
            out.set(0, j, peak.get(j));
        }
        return out;
    }

    /**
     * Helper class to hold extracted distribution parameters
     */
    private static class DistParams {
        ProcessType type;
        Matrix params;
        Double mean;
        Double scv;
        MatrixCell proc;

        DistParams(ProcessType type, Matrix params, Double mean, Double scv, MatrixCell proc) {
            this.type = type;
            this.params = params;
            this.mean = mean;
            this.scv = scv;
            this.proc = proc;
        }
    }

    /**
     * Distribution of the number of calls issued per invocation.
     *
     * A mean below 1 is a call that either happens or does not, hence Bernoulli.
     * Geometric(1/m) is undefined there: its parameter would exceed 1 and its SCV
     * (1-p) would come out negative.
     *
     * @param meanCalls mean number of calls
     * @return the call-count distribution
     */
    public static Distribution callCountDist(double meanCalls) {
        if (Double.isNaN(meanCalls) || meanCalls <= GlobalConstants.FineTol) {
            return new Immediate();
        } else if (meanCalls < 1.0) {
            return new Bernoulli(meanCalls);
        }
        return new Geometric(1.0 / meanCalls);
    }

    /**
     * Extract primitive parameters from a Distribution object
     *
     * @param dist The distribution to extract parameters from
     * @return DistParams containing type, params, mean, scv, and process representation
     */
    private static DistParams extractDistParams(Distribution dist) {
        if (dist == null) {
            return new DistParams(ProcessType.DISABLED, null, Double.NaN, Double.NaN, null);
        }

        ProcessType dtype = ProcessType.fromDistribution(dist);
        Double mean = Double.NaN;
        Double scv = Double.NaN;
        Matrix params = null;
        MatrixCell proc = null;

        // Get mean and SCV
        try {
            mean = dist.getMean();
        } catch (Exception e) {
            mean = Double.NaN;
        }

        try {
            scv = dist.getSCV();
        } catch (Exception e) {
            scv = Double.NaN;
        }

        // Get process representation if available (for Markovian distributions)
        try {
            if (dist instanceof ContinuousDistribution) {
                proc = ((ContinuousDistribution) dist).getProcess();
            }
        } catch (Exception e) {
            proc = null;
        }

        // Extract parameters based on distribution type
        try {
            switch (dtype) {
                case DISABLED:
                case IMMEDIATE:
                    params = null;
                    if (dtype == ProcessType.IMMEDIATE) {
                        mean = 0.0;
                        scv = 0.0;
                    }
                    break;
                case EXP:
                    params = Matrix.singleton((Double) dist.getParam(1).getValue());
                    break;
                case ERLANG:
                    params = new Matrix(2, 1, 2);
                    params.set(0, 0, (Double) dist.getParam(1).getValue());
                    params.set(1, 0, (Double) dist.getParam(2).getValue());
                    break;
                case HYPEREXP:
                    params = new Matrix(3, 1, 3);
                    params.set(0, 0, (Double) dist.getParam(1).getValue());
                    params.set(1, 0, (Double) dist.getParam(2).getValue());
                    params.set(2, 0, (Double) dist.getParam(3).getValue());
                    break;
                case GEOMETRIC:
                    params = Matrix.singleton((Double) dist.getParam(1).getValue());
                    break;
                case DET:
                    params = Matrix.singleton((Double) dist.getParam(1).getValue());
                    break;
                default:
                    // For complex types, store minimal params
                    params = null;
                    break;
            }
        } catch (Exception e) {
            params = null;
        }

        return new DistParams(dtype, params, mean, scv, proc);
    }

    /**
     * Reconstruct a Distribution object from primitive parameters.
     * This is the inverse of extractDistParams and enables migration away from
     * storing Distribution objects in LayeredNetworkStruct.
     *
     * @param type The ProcessType of the distribution
     * @param params The distribution parameters (may be null for some types)
     * @param mean The mean of the distribution
     * @param scv The squared coefficient of variation
     * @param proc The process representation (for PH, MAP types)
     * @return A Distribution object reconstructed from the primitives
     */
    public static Distribution reconstructDistribution(ProcessType type, Matrix params,
            Double mean, Double scv, MatrixCell proc) {
        if (type == null || type == ProcessType.DISABLED) {
            return null;
        }

        switch (type) {
            case IMMEDIATE:
                return Immediate.getInstance();
            case EXP:
                if (params != null && params.length() > 0) {
                    return new Exp(params.get(0, 0));
                } else if (mean != null && !Double.isNaN(mean) && mean > 0) {
                    return new Exp(1.0 / mean);
                }
                return Immediate.getInstance();
            case ERLANG:
                if (params != null && params.length() >= 2) {
                    return new Erlang(params.get(0, 0), (int) params.get(1, 0));
                } else if (mean != null && scv != null && !Double.isNaN(mean) && !Double.isNaN(scv) && scv > 0) {
                    int k = (int) Math.max(1, Math.round(1.0 / scv));
                    return new Erlang(k / mean, k);
                }
                return new Exp(mean != null && mean > 0 ? 1.0 / mean : 1.0);
            case HYPEREXP:
                if (params != null && params.length() >= 3) {
                    return new HyperExp(params.get(0, 0), params.get(1, 0), params.get(2, 0));
                }
                return new Exp(mean != null && mean > 0 ? 1.0 / mean : 1.0);
            case DET:
                if (params != null && params.length() > 0) {
                    return new Det(params.get(0, 0));
                } else if (mean != null && !Double.isNaN(mean)) {
                    return new Det(mean);
                }
                return new Det(1.0);
            case GEOMETRIC:
                if (params != null && params.length() > 0) {
                    return new Geometric(params.get(0, 0));
                }
                return new Geometric(0.5);
            case PH:
            case APH:
                if (proc != null && proc.size() >= 2) {
                    Matrix alpha = proc.get(0);
                    Matrix T = proc.get(1);
                    return new PH(alpha, T);
                }
                return new Exp(mean != null && mean > 0 ? 1.0 / mean : 1.0);
            case MAP:
                if (proc != null && proc.size() >= 2) {
                    Matrix D0 = proc.get(0);
                    Matrix D1 = proc.get(1);
                    return new MAP(D0, D1);
                }
                return new Exp(mean != null && mean > 0 ? 1.0 / mean : 1.0);
            case COXIAN:
                if (proc != null && proc.size() >= 2) {
                    Matrix alpha = proc.get(0);
                    Matrix T = proc.get(1);
                    return new Coxian(alpha, T);
                }
                return new Exp(mean != null && mean > 0 ? 1.0 / mean : 1.0);
            case COX2:
                // Cox2 takes (mu1, mu2, phi1) - extract from params if available
                if (params != null && params.length() >= 3) {
                    return new Cox2(params.get(0, 0), params.get(1, 0), params.get(2, 0));
                } else if (proc != null && proc.size() >= 2) {
                    // Fall back to Coxian representation
                    Matrix alpha = proc.get(0);
                    Matrix T = proc.get(1);
                    return new Coxian(alpha, T);
                }
                return new Exp(mean != null && mean > 0 ? 1.0 / mean : 1.0);
            default:
                // For unsupported types, fall back to Exp with given mean
                if (mean != null && !Double.isNaN(mean) && mean > 0) {
                    return new Exp(1.0 / mean);
                }
                return Immediate.getInstance();
        }
    }

    /**
     * Reconstruct a DiscreteDistribution object from primitive parameters.
     *
     * @param type The ProcessType of the distribution
     * @param params The distribution parameters (may be null for some types)
     * @param mean The mean of the distribution
     * @return A DiscreteDistribution object reconstructed from the primitives
     */
    public static DiscreteDistribution reconstructDiscreteDistribution(ProcessType type,
            Matrix params, Double mean) {
        if (type == null || type == ProcessType.DISABLED) {
            return null;
        }

        switch (type) {
            case GEOMETRIC:
                if (params != null && params.length() > 0) {
                    return new Geometric(params.get(0, 0));
                }
                return new Geometric(0.5);
            case BINOMIAL:
                if (params != null && params.length() >= 2) {
                    return new Binomial((int) params.get(0, 0), params.get(1, 0));
                }
                return new Binomial(1, 0.5);
            case POISSON:
                if (params != null && params.length() > 0) {
                    return new Poisson(params.get(0, 0));
                } else if (mean != null && !Double.isNaN(mean)) {
                    return new Poisson(mean);
                }
                return new Poisson(1.0);
            default:
                return null;
        }
    }
}
