/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;

import jline.GlobalConstants;
import static jline.GlobalConstants.Inf;

import jline.cli.LineDockerClient;
import jline.io.SysUtilsKt;
import jline.io.ModelVisualizer;
import jline.lang.Copyable;
import jline.lang.Ensemble;
import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.constant.CallType;
import jline.lang.constant.ActivityPrecedenceType;
import jline.lang.constant.ProcessType;
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

import static jline.api.lsn.LsnMaxMultiplicityKt.lsnMaxMultiplicity;
import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;
import static jline.io.SysUtilsKt.jlqnGetPath;
import static jline.io.SysUtilsKt.lineTempName;
import static jline.lang.constant.ActivityPrecedenceType.*;

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

                String treplicationString = taskElement.getAttribute("replication");
                double tReplication = treplicationString.isEmpty() ? 1.0 : Double.parseDouble(treplicationString);

                String tmultiplicityString = taskElement.getAttribute("multiplicity");
                double tMultiplicity = tmultiplicityString.isEmpty() ? 1.0 : Double.parseDouble(tmultiplicityString);

                String tthinkTimeMeanString = taskElement.getAttribute("think_time");
                double tThinkTimeMean = tthinkTimeMeanString.isEmpty() ? 0.0 : Double.parseDouble(tthinkTimeMeanString);

                Distribution thinkTime;
                if (tThinkTimeMean <= 0.0) {
                    thinkTime = Immediate.getInstance();
                } else {
                    thinkTime = Exp.fitMean(tThinkTimeMean);
                }
                Task newTask = new Task(myLN, tName, (int) tMultiplicity, SchedStrategy.fromText(tScheduling), thinkTime);
                newTask.setReplication((int) replication);

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

                newTask.on(newProc);  // Assign task to its processor
                taskObj.put(taskObj.size(), newTask);

                NodeList entryList = taskElement.getElementsByTagName("entry");
                for (int k = 0; k < entryList.getLength(); k++) {
                    Element entryElement = (Element) entryList.item(k);
                    String eName = entryElement.getAttribute("name");
                    Entry newEntry = new Entry(myLN, eName);

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
                            if (phase == 1) {
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


                        for (int l = 0; l < nameList.size() - 1; l++) {

                            ActivityPrecedence newPrec = new ActivityPrecedence(Collections.singletonList(nameList.get(l)), Collections.singletonList(nameList.get(l + 1)));
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
                                preParams = new Matrix(1, preActList.getLength());
                                for (int m = 0; m < preActList.getLength(); m++) {
                                    Element preActElement = (Element) preActList.item(m);
                                    preActs.add(preActElement.getAttribute("name"));
                                    preParams.set(0, m, Double.parseDouble(preActElement.getAttribute("prob")));
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

                            // Post-activity parsing
                            String[] postTypes = {ActivityPrecedenceType.POST_SEQ, ActivityPrecedenceType.POST_AND, ActivityPrecedenceType.POST_OR, ActivityPrecedenceType.POST_LOOP, ActivityPrecedenceType.POST_CACHE};
                            NodeList postList = null;
                            String postType = null;
                            for (String type : postTypes) {
                                postType = type;
                                postList = precElement.getElementsByTagName(postType);
                                if (postList.getLength() > 0) break;
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
        }
        return myLN;
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

        String redirectOutput = " > /dev/null";
        if (java.lang.System.getProperty("os.name").startsWith("Windows")) {
            redirectOutput = " > nul 2>&1";
        }

        String cmd = String.format(
                "java -cp %s jlqn.commandline.Jlqn %s %s",
                jlqnPath, filename, redirectOutput
        );

        java.lang.System.out.println("JLQN view model command: " + cmd);
        String output = SysUtilsKt.system(cmd);
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

    private static void writeJLQNEntry(Document document, Element entries, String name, String bndto_activity, String reply_activity, String task) {
        Element entry = document.createElement("entry");
        entry.setAttribute("arrival_rate", "0.0");
        entry.setAttribute("bound_to_activity", bndto_activity);
        entry.setAttribute("name", name);
        entry.setAttribute("priority", "0");
        entry.setAttribute("reply_to_activity", reply_activity);
        entry.setAttribute("task", task);
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
        task.setAttribute("priority", "0");
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
        int idx = 1;

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

        lsn.mult = new Matrix(1, lsn.nhosts + lsn.ntasks + 1, lsn.nhosts + lsn.ntasks);
        lsn.maxmult = new Matrix(1, lsn.nhosts + lsn.ntasks + 1, lsn.nhosts + lsn.ntasks);

        lsn.repl = new Matrix(1, lsn.nhosts + lsn.ntasks + 1, lsn.nhosts + lsn.ntasks);
        lsn.type = new Matrix(1, lsn.nidx + 1, lsn.nidx);
        lsn.graph = new Matrix(lsn.nidx + 1, lsn.nidx + 1, lsn.nidx * lsn.nidx);
        lsn.dag = new Matrix(lsn.nidx + 1, lsn.nidx + 1, lsn.nidx * lsn.nidx);
        lsn.replygraph = new Matrix(lsn.nacts + 1, lsn.nentries + 1, lsn.nentries * lsn.nacts);
        lsn.actphase = new Matrix(1, lsn.nacts + 1, lsn.nacts);  // Phase for each activity (default=1)
        for (int a = 1; a <= lsn.nacts; a++) {
            lsn.actphase.set(0, a, 1.0);  // Default phase is 1
        }

        lsn.nitems = new Matrix(1, lsn.nidx + 1, lsn.ntasks + lsn.nacts);

        lsn.itemcap = new HashMap<>();
        lsn.itemproc = new HashMap<>();
        lsn.itemproc_type = new HashMap<>();
        lsn.itemproc_params = new HashMap<>();
        lsn.itemproc_mean = new HashMap<>();
        lsn.itemproc_scv = new HashMap<>();
        lsn.itemproc_proc = new HashMap<>();


        lsn.iscache = new Matrix(1, lsn.nidx + 1, lsn.nhosts + lsn.ntasks);
        lsn.replacestrat = new Matrix(1, lsn.nidx + 1, lsn.nhosts + lsn.ntasks);

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

        lsn.isfunction = new Matrix(1, lsn.nidx + 1, lsn.nidx);

        lsn.arrival = new HashMap<>();
        lsn.arrival_type = new HashMap<>();
        lsn.arrival_params = new HashMap<>();
        lsn.arrival_mean = new HashMap<>();
        lsn.arrival_scv = new HashMap<>();
        lsn.arrival_proc = new HashMap<>();

        lsn.parent = new Matrix(1, lsn.nidx + 1, lsn.nidx);

        for (int i = 0; i < lsn.nhosts; i++) {
            lsn.sched.put(idx, this.hosts.get(i).scheduling);
            lsn.mult.set(0, idx, this.hosts.get(i).multiplicity);
            lsn.repl.set(0, idx, this.hosts.get(i).replication);
            lsn.names.put(idx, this.hosts.get(i).getName());
            lsn.hashnames.put(idx, "P:" + lsn.names.get(idx));
            lsn.type.set(0, idx, LayeredNetworkElement.HOST);
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
            lsn.mult.set(0, idx, this.tasks.get(i).multiplicity);
            lsn.repl.set(0, idx, this.tasks.get(i).replication);
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
                lsn.hashnames.put(idx, "C:" + lsn.names.get(idx));
            } else if (this.tasks.get(i).hasSetupDelayoff()) {
                // Task has setup/delayoff configured (not just FunctionTask)
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
                lsn.hashnames.put(idx, "F:" + lsn.names.get(idx));
                lsn.isfunction.set(0, idx - 1, 1);  // Matrix uses 0-based indexing
            }

            int pidx = 0;
            Task currentTask = this.tasks.get(i);
            if (currentTask.parent == null) {
                line_error(mfilename(new Object() {}), "Task " + currentTask.getName() + " has no parent processor assigned during XML parsing");
            }
            for (int id = 0; id < this.hosts.size(); id++) {
                if (this.hosts.get(id).getName().equals(currentTask.parent.getName())) {
                    pidx = id + 1;
                    break;
                }
            }

            lsn.parent.set(0, idx, pidx);
            lsn.graph.set(idx, pidx, 1);

            lsn.type.set(0, idx, LayeredNetworkElement.TASK);
            idx++;
        }


        for (int p = 1; p <= lsn.nhosts; p++) {

            if (!lsn.tasksof.containsKey(p)) {
                lsn.tasksof.put(p, new ArrayList<>());
            }
            for (int id = 1; id < lsn.parent.length(); id++) {

                if (lsn.parent.get(0, id) == p) {
                    lsn.tasksof.get(p).add(id);
                }
            }
        }

        for (int e = 0; e < lsn.nentries; e++) {
            idx = lsn.eshift + e + 1;

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
                for (String ra : this.entries.get(e).replyActivity.values()) {
                    int ractidx = Utils.findString(lsn.hashnames, "A:" + ra);
                    if (ractidx > 0) {
                        lsn.replygraph.set(ractidx - lsn.ashift, idx - lsn.eshift, 1);
                    }
                }
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
                    tidx = lsn.nhosts + id + 1;
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
                    tidx = lsn.nhosts + id + 1;
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
            lsn.actphase.set(0, a + 1, this.activities.get(a).getPhase());
            idx++;
        }


        lsn.graph.set(lsn.nidx, lsn.nidx, 0);

        Map<Integer, Task> tasks = this.tasks;
        int cidx = 0;

        lsn.calltype = new HashMap<>();
        lsn.iscaller = new Matrix(lsn.nidx + 1, lsn.nidx + 1, (lsn.ntasks + lsn.nacts) * (lsn.ntasks + lsn.nentries));
        lsn.issynccaller = new Matrix(lsn.nidx + 1, lsn.nidx + 1, (lsn.ntasks + lsn.nacts) * (lsn.ntasks + lsn.nentries));
        lsn.isasynccaller = new Matrix(lsn.nidx + 1, lsn.nidx + 1, (lsn.ntasks + lsn.nacts) * (lsn.ntasks + lsn.nentries));
        lsn.callpair = new Matrix(lsn.nidx + 1, 3, lsn.nidx * 3);
        lsn.callproc = new HashMap<>();
        lsn.callproc_type = new HashMap<>();
        lsn.callproc_params = new HashMap<>();
        lsn.callproc_mean = new HashMap<>();
        lsn.callproc_scv = new HashMap<>();
        lsn.callproc_proc = new HashMap<>();
        lsn.callnames = new HashMap<>();
        lsn.callhashnames = new HashMap<>();
        lsn.taskgraph = new Matrix(lsn.ntasks + lsn.tshift + 1, lsn.ntasks + lsn.tshift + 1, (lsn.ntasks + lsn.tshift) * (lsn.ntasks + lsn.tshift));
        lsn.actpretype = new Matrix(1, lsn.nidx + 1, lsn.nacts);
        lsn.actposttype = new Matrix(1, lsn.nidx + 1, lsn.nacts);

        Matrix loop_back_edges = new Matrix(lsn.nidx + 1, lsn.nidx + 1, lsn.nidx * lsn.nidx);

        // Track boundToEntry mappings to validate uniqueness
        Map<String, String> entryToActivityMap = new HashMap<>();

        // Initialize callsof for all activities to prevent null access in SolverLN
        for (int i = 0; i < this.activities.size(); i++) {
            Activity activity = this.activities.get(i);
            int aidx = Utils.findString(lsn.hashnames, "A:" + activity.getName());
            if (aidx > 0 && !lsn.callsof.containsKey(aidx)) {
                lsn.callsof.put(aidx, new ArrayList<>());
            }
        }

        for (int t = 0; t < lsn.ntasks; t++) {
            int tidx = lsn.tshift + t + 1;

            for (int a = 0; a < tasks.get(t).activities.size(); a++) {
                int aidx = Utils.findString(lsn.hashnames, "A:" + tasks.get(t).activities.get(a).getName());
                if (aidx > 0) {
                    lsn.callsof.put(aidx, new ArrayList<>());
                }

                String boundToEntry = tasks.get(t).activities.get(a).boundToEntry;
                int eidx = Utils.findString(lsn.hashnames, "E:" + boundToEntry);
                if (eidx <= 0) {
                    eidx = Utils.findString(lsn.hashnames, "I:" + boundToEntry);
                }
                if (eidx > 0) {
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

                    if (target_eidx <= 0) {
                        target_eidx = Utils.findString(lsn.hashnames, "I:" + tasks.get(t).activities.get(a).syncCallDests.get(s));
                    }
                    int target_tidx = (int) lsn.parent.get(target_eidx);
                    cidx++;

                    lsn.calltype.put(cidx, CallType.SYNC);
                    lsn.callpair.set(cidx, 1, aidx);
                    lsn.callpair.set(cidx, 2, target_eidx);
                    
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
                    Geometric syncCallDist = new Geometric(1.0 / tasks.get(t).activities.get(a).syncCallMeans.get(s));
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

                for (int s = 0; s < tasks.get(t).activities.get(a).asyncCallDests.size(); s++) {
                    String target_entry_name = tasks.get(t).activities.get(a).asyncCallDests.get(s);
                    int target_eidx = Utils.findString(lsn.hashnames, "E:" + target_entry_name);
                    if (target_eidx <= 0) {
                        target_eidx = Utils.findString(lsn.hashnames, "I:" + target_entry_name);
                    }
                    // Validate that the target entry exists
                    if (target_eidx <= 0) {
                        line_error(mfilename(new Object() {}), "Activity \"" + tasks.get(t).activities.get(a).getName() + "\" has an async call to non-existent entry \"" + target_entry_name + "\".");
                    }
                    int target_tidx = (int) lsn.parent.get(target_eidx);
                    // Check for self-referential async calls (task calling itself)
                    if (tidx == target_tidx) {
                        line_error(mfilename(new Object() {}), "Activity \"" + tasks.get(t).activities.get(a).getName() + "\" in task \"" + tasks.get(t).getName() + "\" has an async call to an entry on the same task. Async self-calls are not supported.");
                    }
                    cidx++;

                    lsn.calltype.put(cidx, CallType.ASYNC);
                    lsn.callpair.set(cidx, 1, aidx);
                    lsn.callpair.set(cidx, 2, target_eidx);
                    lsn.callnames.put(cidx, lsn.names.get(aidx) + "->" + lsn.names.get(target_eidx));
                    lsn.callhashnames.put(cidx, lsn.hashnames.get(aidx) + "->" + lsn.hashnames.get(target_eidx));
                    Geometric asyncCallDist = new Geometric(1.0 / tasks.get(t).activities.get(a).asyncCallMeans.get(s));
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
            if (eidx <= 0) {
                eidx = Utils.findString(lsn.hashnames, "I:" + entry.getName());
            }
            if (eidx <= 0) {
                continue;
            }
            int source_tidx = (int) lsn.parent.get(eidx);

            for (int fw = 0; fw < entry.getForwardingDests().size(); fw++) {
                String target_entry_name = entry.getForwardingDests().get(fw);
                int target_eidx = Utils.findString(lsn.hashnames, "E:" + target_entry_name);
                if (target_eidx <= 0) {
                    target_eidx = Utils.findString(lsn.hashnames, "I:" + target_entry_name);
                }
                if (target_eidx <= 0) {
                    line_error(mfilename(new Object() {}), "Entry \"" + entry.getName() + "\" forwards to non-existent entry \"" + target_entry_name + "\".");
                }
                int target_tidx = (int) lsn.parent.get(target_eidx);

                // Validate: cannot forward to same task
                if (source_tidx == target_tidx) {
                    line_error(mfilename(new Object() {}), "Entry \"" + entry.getName() + "\" cannot forward to entry \"" + target_entry_name + "\" on the same task.");
                }

                cidx++;
                lsn.calltype.put(cidx, CallType.FWD);
                lsn.callpair.set(cidx, 1, eidx);
                lsn.callpair.set(cidx, 2, target_eidx);
                lsn.callnames.put(cidx, lsn.names.get(eidx) + "~>" + lsn.names.get(target_eidx));
                lsn.callhashnames.put(cidx, lsn.hashnames.get(eidx) + "~>" + lsn.hashnames.get(target_eidx));

                // Forwarding probability (not using callproc as forwarding is deterministic choice)
                double fwdProb = entry.getForwardingProbs().get(fw);
                Geometric fwdCallDist = new Geometric(1.0 / fwdProb);
                lsn.callproc.put(cidx, fwdCallDist); // Store as mean calls
                DistParams fwdCallParams = extractDistParams(fwdCallDist);
                lsn.callproc_type.put(cidx, fwdCallParams.type);
                lsn.callproc_params.put(cidx, fwdCallParams.params);
                lsn.callproc_mean.put(cidx, fwdCallParams.mean);
                lsn.callproc_scv.put(cidx, fwdCallParams.scv);
                lsn.callproc_proc.put(cidx, fwdCallParams.proc);

                // Note: Forwarding is a reply redirection mechanism that does NOT create
                // task dependencies or activity precedence. The forwarding call relationship
                // is captured in calltype, callpair, and callnames, which is sufficient for
                // solver processing. We do NOT update taskgraph or graph for forwarding.
                //
                // lsn.taskgraph.set(source_tidx, target_tidx, 1);  // DISABLED - causes cycles
            }
        }

        /* ========== Transform forwarding to pseudo rendezvous calls ========== */
        // LQNS transforms forwarding by adding pseudo rendezvous calls from original
        // callers to forwarding targets. This ensures forwarding targets receive throughput.
        //
        // For each forwarding call from entry A to entry B with probability P:
        //   - Find all calls TO entry A
        //   - For each caller C with rate R, add a pseudo call from C to B with rate R*P
        //
        // This implements the same transformation as LQNS's addForwardingRendezvous().

        int initial_cidx = cidx;  // Save initial count before transformation

        for (int fwd_call_idx = 1; fwd_call_idx <= initial_cidx; fwd_call_idx++) {
            if (lsn.calltype.get(fwd_call_idx) != CallType.FWD) {
                continue;
            }

            int fwd_source_eidx = (int) lsn.callpair.get(fwd_call_idx, 1);  // Entry that does forwarding
            int fwd_target_eidx = (int) lsn.callpair.get(fwd_call_idx, 2);  // Entry receiving forward
            double fwd_prob = lsn.callproc.get(fwd_call_idx).getMean();  // Forwarding probability (stored as mean)

            // Find all SYNC calls TO the forwarding source entry
            for (int caller_call_idx = 1; caller_call_idx <= initial_cidx; caller_call_idx++) {
                if (lsn.calltype.get(caller_call_idx) != CallType.SYNC) {
                    continue;
                }

                int call_target = (int) lsn.callpair.get(caller_call_idx, 2);
                if (call_target == fwd_source_eidx) {
                    // Found a call TO the forwarding source
                    int caller_aidx = (int) lsn.callpair.get(caller_call_idx, 1);  // Activity making the call
                    double call_rate = lsn.callproc.get(caller_call_idx).getMean();

                    // Create pseudo rendezvous from caller to forwarding target
                    cidx++;
                    lsn.calltype.put(cidx, CallType.SYNC);
                    lsn.callpair.set(cidx, 1, caller_aidx);
                    lsn.callpair.set(cidx, 2, fwd_target_eidx);
                    lsn.callnames.put(cidx, lsn.names.get(caller_aidx) + "->" + lsn.names.get(fwd_target_eidx) + "[fwd]");
                    lsn.callhashnames.put(cidx, lsn.hashnames.get(caller_aidx) + "->" + lsn.hashnames.get(fwd_target_eidx) + "[fwd]");

                    // Pseudo call rate = original call rate * forwarding probability
                    double pseudo_rate = call_rate * fwd_prob;
                    lsn.callproc.put(cidx, new Geometric(pseudo_rate));

                    // CRITICAL: Update all matrices so the pseudo call is recognized during submodel construction
                    // (same as regular SYNC calls at lines 1123-1133)
                    int caller_tidx = (int) lsn.parent.get(0, caller_aidx);
                    int target_tidx = (int) lsn.parent.get(0, fwd_target_eidx);

                    lsn.callsof.get(caller_aidx).add(cidx);
                    lsn.iscaller.set(caller_aidx, target_tidx, 1);
                    lsn.iscaller.set(caller_aidx, fwd_target_eidx, 1);
                    lsn.iscaller.set(caller_tidx, target_tidx, 1);
                    lsn.iscaller.set(caller_tidx, fwd_target_eidx, 1);
                    lsn.issynccaller.set(caller_tidx, target_tidx, 1);
                    lsn.issynccaller.set(caller_tidx, fwd_target_eidx, 1);
                    lsn.issynccaller.set(caller_aidx, fwd_target_eidx, 1);
                    lsn.issynccaller.set(caller_aidx, target_tidx, 1);
                    lsn.taskgraph.set(caller_tidx, target_tidx, 1);
                    lsn.graph.set(caller_aidx, fwd_target_eidx, 1);
                }
            }
        }

        lsn.ncalls = cidx;  // Update total number of calls

        for (int t = 0; t < lsn.ntasks; t++) {
            int tidx = lsn.tshift + t + 1;

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
                            if (preaidx <= 0) {
                                line_error(mfilename(new Object() {}), "PRE_AND precedence references non-existent activity \"" + preacts.get(prea) + "\" in task \"" + tasks.get(t).getName() + "\".");
                            }
                            if (preaidx > 0 && lsn.parent.get(preaidx) != tidx) {
                                line_error(mfilename(new Object() {}), "PRE_AND precedence in task \"" + tasks.get(t).getName() + "\" references activity \"" + preacts.get(prea) + "\" from a different task.");
                            }
                        }
                    }

                    // Validate POST_AND activities exist before processing
                    if (posttype.equals(POST_AND)) {
                        if (postacts.isEmpty()) {
                            line_error(mfilename(new Object() {}), "POST_AND precedence in task \"" + tasks.get(t).getName() + "\" has no post activities.");
                        }
                        for (int posta = 0; posta < postacts.size(); posta++) {
                            int postaidx = Utils.findString(lsn.hashnames, "A:" + postacts.get(posta));
                            if (postaidx <= 0) {
                                line_error(mfilename(new Object() {}), "POST_AND precedence references non-existent activity \"" + postacts.get(posta) + "\" in task \"" + tasks.get(t).getName() + "\".");
                            }
                            if (postaidx > 0 && lsn.parent.get(postaidx) != tidx) {
                                line_error(mfilename(new Object() {}), "POST_AND precedence in task \"" + tasks.get(t).getName() + "\" references activity \"" + postacts.get(posta) + "\" from a different task.");
                            }
                        }
                    }

                    for (int prea = 0; prea < preacts.size(); prea++) {
                        int preaidx = Utils.findString(lsn.hashnames, "A:" + tasks.get(t).precedences.get(ap).preActs.get(prea));
                        double preParam = 1.0;
                        if (pretype.equals(PRE_AND)) {
                            Matrix quorum = tasks.get(t).precedences.get(ap).preParams;
                            if (!quorum.isEmpty()) {
                                preParam = quorum.get(0) / preacts.size(); // Use first element to match MATLAB behavior
                            }
                        }

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
                                lsn.actposttype.set(0, loopEndAidx, ActivityPrecedence.getPrecedenceId(tasks.get(t).precedences.get(ap).postType));
                                break;
                            default:
                                for (int posta = 0; posta < postacts.size(); posta++) {
                                    int postaidx = Utils.findString(lsn.hashnames, "A:" + tasks.get(t).precedences.get(ap).postActs.get(posta));
                                    double postParam = 1;
                                    lsn.graph.set(preaidx, postaidx, preParam * postParam);
                                    lsn.actpretype.set(0, preaidx, ActivityPrecedence.getPrecedenceId(tasks.get(t).precedences.get(ap).preType));
                                    lsn.actposttype.set(0, postaidx, ActivityPrecedence.getPrecedenceId(tasks.get(t).precedences.get(ap).postType));
                                }
                        }
                    }
                }
            }

        }

        /* Compute entry-to-activity reachability within the same task */
        for (int eoff = 0; eoff < lsn.nentries; eoff++) {
            int eidx = lsn.eshift + eoff + 1; // global entry index
            int tidx = (int) lsn.parent.get(eidx);
            boolean[] visited = new boolean[lsn.nidx + 1];
            Deque<Integer> stack = new ArrayDeque<>();
            stack.push(eidx);
            visited[eidx] = true;
            while (!stack.isEmpty()) {
                int v = stack.pop();
                for (int nbr = 1; nbr <= lsn.nidx; nbr++) {
                    if (lsn.graph.get(v, nbr) != 0 && !visited[nbr]) {
                        visited[nbr] = true;
                        stack.push(nbr);
                    }
                }
            }
            List<Integer> acts = new ArrayList<>();
            for (int n = 1; n <= lsn.nidx; n++) {
                if (visited[n] && lsn.type.get(n) == LayeredNetworkElement.ACTIVITY && lsn.parent.get(n) == tidx) {
                    acts.add(n);
                }
            }
            lsn.actsof.put(eidx, acts);
        }

        for (int t = 0; t < lsn.ntasks; t++) {
            int tidx = lsn.tshift + t + 1;
            for (int aidx : lsn.actsof.getOrDefault(tidx, new ArrayList<>())) {
                List<Integer> postaidxs = new ArrayList<>();
                // Start from col=1 to skip the 0-padding column (Java uses 1-based indexing in matrices)
                for (int col = 1; col < lsn.graph.getNumCols(); col++) {
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
                        // Start from row=1 to skip the 0-padding row (Java uses 1-based indexing in matrices)
                        for (int row = 1; row < lsn.graph.getNumRows(); row++) {
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
                int i = lsn.ncalls + 1; i < lsn.callpair.getNumRows(); i++) {
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
                // Start from row=1 to skip the 0-padding row (Java uses 1-based indexing in matrices)
                for (int row = 1; row < lsn.taskgraph.getNumRows(); row++) {
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

        lsn.isref = new Matrix(1, lsn.nhosts + lsn.ntasks + 1, lsn.ntasks);
        for (int col = 1; col <= lsn.sched.size(); col++) {
            if (lsn.sched.get(col) == SchedStrategy.REF) {
                lsn.isref.set(0, col, 1);
            }
        }

        // Create schedid matrix (scheduling strategy ordinal values) for Python compatibility
        lsn.schedid = new Matrix(1, lsn.nhosts + lsn.ntasks + 1, lsn.nhosts + lsn.ntasks);
        for (int col = 1; col <= lsn.sched.size(); col++) {
            SchedStrategy strategy = lsn.sched.get(col);
            if (strategy != null) {
                lsn.schedid.set(0, col, strategy.ordinal());
            }
        }

        // Create replacement alias for replacestrat (for Python compatibility)
        lsn.replacement = lsn.replacestrat;

        // Only process cache items if there are any
        if (lsn.nitems.getNumCols() > 1) {
            for (int i = 1; i < lsn.nitems.getNumCols(); i++) {
                if (lsn.nitems.get(0, i) > 0) {
                    lsn.iscache.set(0, i, 1);
                }
            }
        }

        /*
         * The dag differs from the graph:
         * - dag swaps the direction of entry-task edges
         * - dag removes loop edges
         */
        Matrix dag = lsn.graph.copy();
        for (int i = 1; i <= lsn.nidx; i++) {
            if (lsn.type.get(i) == LayeredNetworkElement.TASK &&
                    lsn.isref.get(0, i) == 0) {
                for (int j = 1; j <= lsn.nidx; j++) {
                    if (lsn.type.get(j) == LayeredNetworkElement.ENTRY && dag.get(i, j) != 0) {
                        dag.set(i, j, 0);
                        dag.set(j, i, 1);
                    }
                }
            }
        }
        for (int r = 1; r <= lsn.nidx; r++) {
            for (int c = 1; c <= lsn.nidx; c++) {
                if (loop_back_edges.get(r, c) != 0) {
                    dag.set(r, c, 0);
                }
            }
        }

        lsn.dag = dag;

        int newRow = lsn.taskgraph.getNumCols() - lsn.nhosts - 1;
        int newCol = lsn.taskgraph.getNumRows() - lsn.nhosts - 1;
        Matrix taskgraphSection = new Matrix(newRow, newCol);
        for (
                int r = 0;
                r < newRow; r++) {
            for (int c = 0; c < newCol; c++) {
                taskgraphSection.set(r, c, lsn.taskgraph.get(lsn.nhosts + r + 1, lsn.nhosts + c + 1));
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
            roots.add(fue + 1);
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
        lsn.conntasks = labels;
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

        // Compute bounds on multiplicies for host processors and non-ref tasks
        if (DirectedGraph.isDAG(dag)) {               // topological-sort based test
            lsn.maxmult = lsnMaxMultiplicity(lsn).transpose();
        } else {
            line_error(mfilename(new Object() {
            }), "A cycle exists in an activity graph.");
        }

        /* Check for non-terminal reply activities */
        /* An activity that replies to an entry should not have Phase 1 successor activities */
        /* Phase 2 successors are allowed (post-reply processing) */
        for (int a = 0; a < lsn.replygraph.getNumRows(); a++) {
            for (int b = 0; b < lsn.replygraph.getNumCols(); b++) {
                if (lsn.replygraph.get(a, b) > 0) {          // activity 'a' replies
                    int aidx = lsn.ashift + a;                // global activity index
                    for (int succ = 0; succ < lsn.graph.getNumCols(); succ++) {
                        if (lsn.graph.get(aidx, succ) != 0) { // successor activity
                            if (succ > lsn.eshift + lsn.nentries) {
                                int succActIdx = succ - lsn.ashift;  // convert to activity array index (0-based)
                                int succPhase = (int) lsn.actphase.get(0, succActIdx + 1);  // actphase uses 1-based column index
                                if (succPhase == 1) {
                                    line_error(mfilename(new Object() {
                                    }), "Unsupported replyTo in non-terminal activity.");
                                }
                            }
                        }
                    }
                }
            }
        }

        /* ------------------------------------------------------------------ */
        /* 4.  Call check: an entry may be called either synchronously or      */
        /*     asynchronously, but not both.                                   */
        /* ------------------------------------------------------------------ */
        if (lsn.callpair != null && lsn.callpair.getNumRows() > 0) {

            /* Remember the first call-kind we encounter for each entry */
            Map<Integer, CallType> entryCallKind = new HashMap<>();

            int nCalls = lsn.callpair.getNumRows();
            for (int iter_cidx = 0; iter_cidx < nCalls; iter_cidx++) {

                int targetEidx = (int) lsn.callpair.get(iter_cidx, 2);   // 3-rd column (0-based index)
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
            writeJLQNEntry(document, entries, E.getName(), E_boundTo, E_replyTo, E.parent.getName());
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
            ncalls += A.syncCallDests.size() + A.asyncCallDests.size(); // TODO: add forwarding calls
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
            // TODO: forwarding calls
//            for (int c=0; c<A.fwdCallDests.size(); c++) {
//                String call_dest = A.fwdCallDests.get(c);
//                writeJLQNCall(document, calls, A.getName()+"->"+call_dest, A.getName(), call_dest, CallType.FWD);
//            }
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
        System.out.println("JLQN file saved as " + filename);
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
        /* ------------------------------------------------------------------
         * First pass – build a name→hash map (optionally with abstract names)
         * ------------------------------------------------------------------ */
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

        /* ------------------------------------------------------------------
         * Second pass – emit processors, tasks, entries, activities
         * ------------------------------------------------------------------ */
        for (Integer pId : hostIds) {
            Processor curProc = (Processor) this.hosts.get(pId);
            Element procElement = doc.createElement("processor");
            rootElement.appendChild(procElement);
            procElement.setAttribute("name", nodeHashMap.get(curProc.getName()));
            procElement.setAttribute("scheduling", SchedStrategy.toText(curProc.scheduling));
            if (curProc.replication > 1) {
                procElement.setAttribute("replication", Integer.toString(curProc.replication));
            }
            if (!curProc.scheduling.equals(SchedStrategy.INF)) {
                procElement.setAttribute("multiplicity", Integer.toString(curProc.multiplicity));
            }
            if (curProc.scheduling.equals(SchedStrategy.PS) ||
                    curProc.scheduling.equals(SchedStrategy.PSPRIO)) {
                procElement.setAttribute("quantum", Double.toString(curProc.quantum));
            }
            procElement.setAttribute("speed-factor", Double.toString(curProc.speedFactor));

            /* ----------------------------- TASKS -------------------------- */
            for (int t = 0; t < curProc.tasks.size(); t++) {
                Task curTask = curProc.tasks.get(t);
                Element taskElement = doc.createElement("task");
                procElement.appendChild(taskElement);
                taskElement.setAttribute("name", nodeHashMap.get(curTask.getName()));
                taskElement.setAttribute("scheduling", SchedStrategy.toText(curTask.scheduling));
                if (curTask.replication > 1) {
                    taskElement.setAttribute("replication", Integer.toString(curTask.replication));
                }
                if (!curTask.scheduling.equals(SchedStrategy.INF)) {
                    taskElement.setAttribute("multiplicity", Integer.toString(curTask.multiplicity));
                }
                if (curTask.scheduling.equals(SchedStrategy.REF)) {
                    taskElement.setAttribute("think-time", Double.toString(curTask.thinkTimeMean));
                }

                /* ----------- ENTRIES -------------- */
                for (int e = 0; e < curTask.entries.size(); e++) {
                    Entry curEntry = curTask.entries.get(e);
                    Element entryElement = doc.createElement("entry");
                    taskElement.appendChild(entryElement);
                    entryElement.setAttribute("name", nodeHashMap.get(curEntry.getName()));
                    entryElement.setAttribute("type", "NONE");

                    /* forwarding calls */
                    for (int fw = 0; fw < curEntry.getForwardingDests().size(); fw++) {
                        Element fwdElement = doc.createElement("forwarding");
                        entryElement.appendChild(fwdElement);
                        fwdElement.setAttribute("dest", nodeHashMap.get(curEntry.getForwardingDests().get(fw)));
                        fwdElement.setAttribute("prob", Double.toString(curEntry.getForwardingProbs().get(fw)));
                    }
                }

                /* ----------- ACTIVITIES ----------- */
                Element taskActsElement = doc.createElement("task-activities");
                taskElement.appendChild(taskActsElement);

                for (int a = 0; a < curTask.activities.size(); a++) {
                    Activity curAct = curTask.activities.get(a);
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
                }

                /* ----------- PRECEDENCES ---------- */
                for (int ap = 0; ap < curTask.precedences.size(); ap++) {
                    ActivityPrecedence curActPrec = curTask.precedences.get(ap);
                    Element actPrecElement = doc.createElement("precedence");
                    taskActsElement.appendChild(actPrecElement);

                    Element preElement = doc.createElement(curActPrec.preType);
                    actPrecElement.appendChild(preElement);
                    if (curActPrec.preType.equals(PRE_AND) &&
                            curActPrec.preParams != null &&
                            !curActPrec.preParams.isEmpty()) {
                        int quorum = (int) curActPrec.preParams.get(0); // Use first element to match MATLAB behavior
                        if (quorum != curActPrec.preActs.size()) {
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
                    } else {
                        for (String poa : curActPrec.postActs) {
                            Element postActElement = doc.createElement("activity");
                            postElement.appendChild(postActElement);
                            postActElement.setAttribute("name", nodeHashMap.get(poa));
                        }
                    }
                }

                /* ----------- REPLY-ACTIVITIES ----- */
                for (Entry curEntry : curTask.entries) {

                    /* If the entry has no reply-activities yet, infer them from lsn.replygraph */
                    if (curEntry.replyActivity.isEmpty()) {
                        LayeredNetworkStruct lsn = this.getStruct();

                        /* locate the zero-based column of this entry inside lsn.replygraph   */
                        int eidx = -1;
                        for (Map.Entry<Integer, Entry> kv : this.entries.entrySet()) {
                            if (kv.getValue() == curEntry) {
                                eidx = kv.getKey() - 1;              // convert 1-based → 0-based
                                break;
                            }
                        }

                        /* add every activity that replies to this entry */
                        if (eidx >= 0) {
                            for (int ra = 0; ra < lsn.replygraph.getNumRows(); ra++) {
                                if (lsn.replygraph.get(ra, eidx) != 0) {
                                    String actName = lsn.names.get(lsn.ashift + ra + 1); // back to 1-based
                                    // Find the activity index by name
                                    Integer actIndex = findActivityIndexByName(actName);
                                    if (actIndex != null) {
                                        curEntry.replyActivity.put(actIndex, actName);
                                    }
                                }
                            }
                        }
                    }

                    /* Emit <reply-entry> only if at least one reply-activity exists */
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

        /* ------------------------------------------------------------------
         * write XML to file
         * ------------------------------------------------------------------ */
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
    public FeatureSet getUsedLangFeatures() {
        if (this.usedFeatures == null) {
            this.usedFeatures = new FeatureSet();
        }

        // Analyze hosts (processors)
        for (Host host : this.hosts.values()) {
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
                default:
                    break;
            }
        }

        // Analyze tasks
        for (Task task : this.tasks.values()) {
            // Mark task scheduling features
            switch (task.scheduling) {
                case REF:
                    // Reference tasks are special in layered networks
                    break;
                case INF:
                    this.usedFeatures.setTrue("SchedStrategy_INF");
                    break;
                case FCFS:
                    this.usedFeatures.setTrue("SchedStrategy_FCFS");
                    break;
                default:
                    break;
            }

            // Analyze think time distribution
            if (task.thinkTime != null && !(task.thinkTime instanceof Immediate)) {
                String distName = task.thinkTime.getClass().getSimpleName();
                this.usedFeatures.setTrue(distName);
            }

            // Check if it's a cache task
            if (task instanceof CacheTask) {
                this.usedFeatures.setTrue("Cache");
                CacheTask cacheTask = (CacheTask) task;
                // Mark replacement strategy features if available
                // Note: replacement strategy handling would go here
            }
        }

        // Analyze activities
        for (Activity activity : this.activities.values()) {
            // Analyze host demand distribution
            if (activity.hostDemand != null && !(activity.hostDemand instanceof Immediate)) {
                String distName = activity.hostDemand.getClass().getSimpleName();
                this.usedFeatures.setTrue(distName);
            }

            // Mark call features
            if (!activity.syncCallDests.isEmpty()) {
                // Has synchronous calls
            }
            if (!activity.asyncCallDests.isEmpty()) {
                // Has asynchronous calls
            }
        }

        // Analyze entries
        for (Entry entry : this.entries.values()) {
            if (entry instanceof ItemEntry) {
                // Item-based entries for cache modeling
                this.usedFeatures.setTrue("Cache");
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
