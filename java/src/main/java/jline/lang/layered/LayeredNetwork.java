package jline.lang.layered;

import jline.lang.Ensemble;
import jline.lang.Network;
import jline.lang.constant.CallType;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.DiscreteDistribution;
import jline.lang.distributions.Distribution;
import jline.lang.distributions.Geometric;
import jline.lang.distributions.Immediate;
import jline.solvers.ln.SolverLN;
import jline.util.Matrix;
import jline.util.Utils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.*;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import java.io.*;
import java.util.*;

import static jline.io.SysUtils.lineTempName;
import static jline.lang.constant.ActivityPrecedenceType.*;

/**
 * A layered queueing network model
 */
public class LayeredNetwork extends Ensemble {

    //private final Aux aux;
    private Matrix lqnGraph;
    private Matrix taskGraph;
    private Param param;

    protected Map<Integer, Host> hosts;
    protected Map<Integer, Task> tasks;
    protected Map<Integer, Task> reftasks;
    protected Map<Integer, Activity> activities;
    protected Map<Integer, Entry> entries;

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
        this.param = new Param();
        this.param.Nodes.RespT = 0;
        this.param.Nodes.Tput = 0;
        this.param.Nodes.Util = 0;
        this.param.Edges.RespT = 0;
        this.param.Edges.Tput = 0;
    }

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

    public void generateGraph() {
    }

    public void initDefault() {
    }

    public void sendModel(String outputPath, String portNumber) {
        this.sendModel(outputPath, "127.0.0.1", portNumber);
    }

    public void sendModel(String outputPath, String ipNumber, String portNumber) {
        String filePath = null;
        try {
            filePath = lineTempName("layered");
        } catch (IOException e) {
            return;
        }
        try {
            this.writeXML(filePath + "/model.xml", true);
            File modelXML = new File(filePath + "/model.xml");
            File out = new File(outputPath);
            ProcessBuilder processBuilder = new ProcessBuilder("java", "-jar", filePath + "/lineclient.jar", ipNumber, portNumber, "-i", "lqnx");
            processBuilder.redirectErrorStream(true);
            processBuilder.redirectInput(modelXML);
            processBuilder.redirectOutput(out);
            Process process = processBuilder.start();
            process.waitFor();
            InputStream inputStream = process.getInputStream();
            BufferedReader input = new BufferedReader(new InputStreamReader(inputStream));
            String ss = null;
            while ((ss = input.readLine()) != null) {
                System.out.println(ss);
            }
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }


    public void writeXML(String filename, boolean abstractNames) {
        Map<String, String> nodeHashMap = new HashMap<>();

        int tctr = 0, ectr = 0, actr = 0;

        if (abstractNames) {
            for (int p = 0; p < this.hosts.size(); p++) {
                Processor curProc = (Processor) this.hosts.get(p);
                nodeHashMap.put(curProc.getName(), "P" + p);
                for (int t = 0; t < curProc.tasks.size(); t++) {
                    Task curTask = curProc.tasks.get(t);
                    tctr++;
                    nodeHashMap.put(curTask.getName(), "T" + tctr);
                    for (int e = 0; e < curTask.entries.size(); e++) {
                        Entry curEntry = curTask.entries.get(e);
                        ectr++;
                        nodeHashMap.put(curEntry.getName(), "E" + ectr);
                    }
                    for (int a = 0; a < curTask.activities.size(); a++) {
                        Activity curAct = curTask.activities.get(a);
                        actr++;
                        nodeHashMap.put(curAct.getName(), "A" + actr);
                    }
                }
            }
        } else {
            for (int p = 0; p < this.hosts.size(); p++) {
                Processor curProc = (Processor) this.hosts.get(p);
                nodeHashMap.put(curProc.getName(), curProc.getName());
                for (int t = 0; t < curProc.tasks.size(); t++) {
                    Task curTask = curProc.tasks.get(t);
                    nodeHashMap.put(curTask.getName(), curTask.getName());
                    for (int e = 0; e < curTask.entries.size(); e++) {
                        Entry curEntry = curTask.entries.get(e);
                        nodeHashMap.put(curEntry.getName(), curEntry.getName());
                    }
                    for (int a = 0; a < curTask.activities.size(); a++) {
                        Activity curAct = curTask.activities.get(a);
                        nodeHashMap.put(curAct.getName(), curAct.getName());
                    }
                }
            }
        }

        String precision = "%10.15e";
        DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder docBuilder = null;
        try {
            docBuilder = docFactory.newDocumentBuilder();
        } catch (ParserConfigurationException e) {
            e.printStackTrace();
        }
        Document doc = docBuilder.newDocument();
        Element rootElement = doc.createElement("lqn-model");
        doc.appendChild(rootElement);
        rootElement.setAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
        rootElement.setAttribute("xsi:noNamespaceSchemaLocation", "lqn.xsd");
        rootElement.setAttribute("name", this.getName());

        for (int p = 0; p < this.hosts.size(); p++) {
            Processor curProc = (Processor) this.hosts.get(p);
            Element procElement = doc.createElement("processor");
            rootElement.appendChild(procElement);
            procElement.setAttribute("name", nodeHashMap.get(curProc.getName()));
            procElement.setAttribute("scheduling", SchedStrategy.toText(curProc.scheduling));
            if (curProc.replication > 1) {
                procElement.setAttribute("replication", Integer.toString(curProc.replication));
            }
            if (!curProc.scheduling.equals(SchedStrategy.INF)) {
                String mult = Integer.toString(curProc.multiplicity);
                procElement.setAttribute("multiplicity", mult);
            }
            if (curProc.scheduling.equals(SchedStrategy.PS)) {
                procElement.setAttribute("quantum", Double.toString(curProc.quantum));
            }
            procElement.setAttribute("speed-factor", Double.toString(curProc.speedFactor));
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
                for (int e = 0; e < curTask.entries.size(); e++) {
                    Entry curEntry = curTask.entries.get(e);
                    Element entryElement = doc.createElement("entry");
                    taskElement.appendChild(entryElement);
                    entryElement.setAttribute("name", nodeHashMap.get(curEntry.getName()));
                    entryElement.setAttribute("type", "NONE");
                }
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

                    for (int sc = 0; sc < curAct.syncCallDests.size(); sc++) {

                        Element syncCallElement = doc.createElement("synch-call");
                        actElement.appendChild(syncCallElement);
                        syncCallElement.setAttribute("dest", nodeHashMap.get(curAct.syncCallDests.get(sc)));
                        syncCallElement.setAttribute("calls-mean", Double.toString(curAct.syncCallMeans.get(sc)));
                    }
                    for (int ac = 0; ac < curAct.asyncCallDests.size(); ac++) {
                        Element asyncCallElement = doc.createElement("asynch-call");
                        actElement.appendChild(asyncCallElement);
                        asyncCallElement.setAttribute("dest", nodeHashMap.get(curAct.asyncCallDests.get(ac)));
                        asyncCallElement.setAttribute("calls-mean", Double.toString(curAct.asyncCallMeans.get(ac)));
                    }
                }
                for (int ap = 0; ap < curTask.precedences.size(); ap++) {
                    ActivityPrecedence curActPrec = curTask.precedences.get(ap);
                    Element actPrecElement = doc.createElement("precedence");
                    taskActsElement.appendChild(actPrecElement);

                    Element preElement = doc.createElement(curActPrec.preType);
                    actPrecElement.appendChild(preElement);
                    if (curActPrec.preType.equals(PRE_AND) && !curActPrec.preParams.isEmpty()) {
                        preElement.setAttribute("quorum", Double.toString(curActPrec.preParams.get(1)));
                    }
                    for (int pra = 0; pra < curActPrec.preActs.size(); pra++) {
                        Element preActElement = doc.createElement("activity");
                        preElement.appendChild(preActElement);
                        preActElement.setAttribute("name", nodeHashMap.get(curActPrec.preActs.get(pra)));
                    }

                    Element postElement = doc.createElement(curActPrec.postType);
                    actPrecElement.appendChild(postElement);
                    if (curActPrec.postType.equals(POST_OR)) {
                        for (int poa = 0; poa < curActPrec.postActs.size(); poa++) {
                            Element postActElement = doc.createElement("activity");
                            postElement.appendChild(postActElement);
                            postActElement.setAttribute("name", nodeHashMap.get(curActPrec.postActs.get(poa)));
                            postActElement.setAttribute("prob", Double.toString(curActPrec.postParams.get(poa)));
                        }
                    } else if (curActPrec.postType.equals(POST_LOOP)) {
                        for (int poa = 0; poa < curActPrec.postActs.size() - 1; poa++) {
                            Element postActElement = doc.createElement("activity");
                            postElement.appendChild(postActElement);
                            postActElement.setAttribute("name", nodeHashMap.get(curActPrec.postActs.get(poa)));
                            postActElement.setAttribute("count", Double.toString(curActPrec.postParams.get(poa)));
                        }
                        postElement.setAttribute("end", curActPrec.postActs.get(curActPrec.postActs.size() - 1));//{end}?
                    } else {
                        for (int poa = 0; poa < curActPrec.postActs.size(); poa++) {
                            Element postActElement = doc.createElement("activity");
                            postElement.appendChild(postActElement);
                            postActElement.setAttribute("name", nodeHashMap.get(curActPrec.postActs.get(poa)));
                        }
                    }
                }
                if (curTask.scheduling != SchedStrategy.REF) {
                    for (int e = 0; e < curTask.entries.size(); e++) {
                        Entry curEntry = curTask.entries.get(e);
                        if (!curEntry.replyActivity.isEmpty()) {
                            Element entryReplyElement = doc.createElement("reply-entry");
                            taskActsElement.appendChild(entryReplyElement);
                            entryReplyElement.setAttribute("name", nodeHashMap.get(curEntry.getName()));
                            for (int r = 0; r < curEntry.replyActivity.size(); r++) {
                                Element entryReplyActElement = doc.createElement("reply-activity");
                                entryReplyElement.appendChild(entryReplyActElement);
                                entryReplyActElement.setAttribute("name", nodeHashMap.get(curEntry.replyActivity.get(r)));
                            }
                        }
                    }
                }
            }
        }
        TransformerFactory transformerFactory = TransformerFactory.newInstance();
        try {
            Transformer transformer = transformerFactory.newTransformer();
            transformer.setOutputProperty(OutputKeys.INDENT, "yes");
            DOMSource source = new DOMSource(doc);
            StreamResult result = new StreamResult(new File(filename));//todo
            transformer.transform(source, result);
        } catch (TransformerConfigurationException e) {
            e.printStackTrace();
        } catch (TransformerException e) {
            e.printStackTrace();
        }
    }


    public LayeredNetworkStruct getStruct() {

        LayeredNetworkStruct lqn = new LayeredNetworkStruct();

        lqn.nidx = 0;
        lqn.hshift = 0;
        lqn.cshift = 0;
        lqn.nhosts = this.hosts.size();
        lqn.ntasks = this.tasks.size();
        lqn.nentries = this.entries.size();
        lqn.nacts = this.activities.size();
        lqn.tshift = lqn.nhosts;
        lqn.eshift = lqn.nhosts + lqn.ntasks;
        lqn.ashift = lqn.nhosts + lqn.ntasks + lqn.nentries;

        // analyze static properties
        lqn.nidx = lqn.nhosts + lqn.ntasks + lqn.nentries + lqn.nacts;
        int idx = 1;

        lqn.tasksof = new HashMap<>(lqn.nhosts);
        lqn.entriesof = new HashMap<>(lqn.nhosts + lqn.ntasks);
        lqn.actsof = new HashMap<>(lqn.nhosts + lqn.ntasks);
        lqn.callsof = new HashMap<>(lqn.nacts);

        lqn.hostdem = new HashMap<>();
        lqn.think = new HashMap<>();
        lqn.sched = new HashMap<>();
        lqn.schedid = new Matrix(1, lqn.nhosts + lqn.ntasks + 1, lqn.nhosts + lqn.ntasks);

        lqn.names = new HashMap<>();
        lqn.hashnames = new HashMap<>();

        lqn.mult = new Matrix(1, lqn.nhosts + lqn.ntasks + 1, lqn.nhosts + lqn.ntasks);
        lqn.repl = new Matrix(1, lqn.nhosts + lqn.ntasks + 1, lqn.nhosts + lqn.ntasks);
        lqn.type = new Matrix(1, lqn.nidx + 1, lqn.nidx);
        lqn.graph = new Matrix(lqn.nidx + 1, lqn.nidx + 1, lqn.nidx * lqn.nidx);
        lqn.replygraph = new Matrix(lqn.nidx + 1, lqn.nidx + 1, lqn.nentries * lqn.nacts);


        lqn.nitems = new Matrix(1, lqn.nhosts + lqn.ntasks + lqn.nentries + 1, lqn.ntasks + lqn.nacts);

        lqn.itemcap = new HashMap<>();
        lqn.itemproc = new HashMap<>();


        lqn.iscache = new Matrix(1, lqn.nhosts + lqn.ntasks + 1, lqn.nhosts + lqn.ntasks);
        lqn.replacement = new Matrix(1, lqn.nhosts + lqn.ntasks + 1, lqn.nhosts + lqn.ntasks);


        lqn.parent = new Matrix(1, lqn.nidx + 1, lqn.nidx);

        for (int i = 0; i < lqn.nhosts; i++) {

            lqn.sched.put(idx, this.hosts.get(i).scheduling);
            lqn.schedid.set(0, idx, SchedStrategy.toID(lqn.sched.get(idx)));
            lqn.mult.set(0, idx, this.hosts.get(i).multiplicity);
            lqn.repl.set(0, idx, this.hosts.get(i).replication);
            lqn.names.put(idx, this.hosts.get(i).getName());
            lqn.hashnames.put(idx, "P:" + lqn.names.get(idx));
            lqn.type.set(0, idx, LayeredNetworkElement.HOST);
            idx = idx + 1;
        }

        for (int i = 0; i < lqn.ntasks; i++) {

            lqn.sched.put(idx, this.tasks.get(i).scheduling);
            lqn.schedid.set(0, idx, SchedStrategy.toID(lqn.sched.get(idx)));
            lqn.hostdem.put(idx, new Immediate());
            lqn.think.put(idx, this.tasks.get(i).thinkTime);
            lqn.mult.set(0, idx, this.tasks.get(i).multiplicity);
            lqn.repl.set(0, idx, this.tasks.get(i).replication);
            lqn.names.put(idx, this.tasks.get(i).getName());

            if (lqn.schedid.get(0, idx) == SchedStrategy.toID(SchedStrategy.REF)) {
                lqn.hashnames.put(idx, "R:" + lqn.names.get(idx));
            } else {
                lqn.hashnames.put(idx, "T:" + lqn.names.get(idx));
            }

            if (this.tasks.get(i) instanceof CacheTask) {
                lqn.nitems.set(0, idx, ((CacheTask) this.tasks.get(i)).items);
                lqn.itemcap.put(idx, ((CacheTask) this.tasks.get(i)).itemLevelCap);
                lqn.replacement.set(0, idx, ((CacheTask) this.tasks.get(i)).replacementPolicy);
                lqn.hashnames.put(idx, "C:" + lqn.names.get(idx));
            }

            int pidx = 0;
            for (int id = 0; id < this.hosts.size(); id++) {
                if (this.hosts.get(id).getName().equals(this.tasks.get(i).parent.getName())) {
                    pidx = id + 1;
                    break;
                }
            }

            lqn.parent.set(0, idx, pidx);
            lqn.graph.set(idx, pidx, 1);//col,row?

            lqn.type.set(0, idx, LayeredNetworkElement.TASK);
            idx++;
        }


        for (int p = 1; p <= lqn.nhosts; p++) {

            if (!lqn.tasksof.containsKey(p)) {
                lqn.tasksof.put(p, new ArrayList<>());
            }
            for (int id = 1; id < lqn.parent.length(); id++) {

                if (lqn.parent.get(0, id) == p) {
                    lqn.tasksof.get(p).add(id);
                }
            }
        }

        for (int e = 0; e < lqn.nentries; e++) {

            lqn.names.put(idx, this.entries.get(e).getName());

            lqn.hashnames.put(idx, "E:" + lqn.names.get(idx));

            if (this.entries.get(e) instanceof ItemEntry) {
                lqn.hashnames.put(idx, "I:" + lqn.names.get(idx));
                lqn.nitems.set(0, idx, ((ItemEntry) this.entries.get(e)).getCardinality());
                lqn.itemproc.put(idx, (DiscreteDistribution) ((ItemEntry) this.entries.get(e)).getPopularity());
            }
            lqn.hostdem.put(idx, new Immediate());
            int tidx = 0;
            for (int id = 0; id < this.tasks.size(); id++) {
                if (this.entries.get(e).parent.getName().equals(this.tasks.get(id).getName())) {
                    tidx = lqn.nhosts + id + 1;
                    break;
                }
            }
            lqn.parent.set(0, idx, tidx);
            lqn.graph.set(tidx, idx, 1);
            if (!lqn.entriesof.containsKey(tidx)) {
                lqn.entriesof.put(tidx, new ArrayList<>());

            }
            lqn.entriesof.get(tidx).add(idx);

            lqn.type.set(0, idx, LayeredNetworkElement.ENTRY);
            idx++;
        }

        for (int a = 0; a < lqn.nacts; a++) {

            lqn.names.put(idx, this.activities.get(a).getName());
            lqn.hashnames.put(idx, "A:" + lqn.names.get(idx));
            lqn.hostdem.put(idx, this.activities.get(a).hostDemand);
            int tidx = 0;
            for (int id = 0; id < this.tasks.size(); id++) {
                if (this.activities.get(a).parent.getName().equals(this.tasks.get(id).getName())) {
                    tidx = lqn.nhosts + id + 1;
                    break;
                }
            }
            lqn.parent.set(0, idx, tidx);
            if (!lqn.actsof.containsKey(tidx)) {
                lqn.actsof.put(tidx, new ArrayList<>());
            }
            lqn.actsof.get(tidx).add(idx);
            lqn.type.set(0, idx, LayeredNetworkElement.ACTIVITY);
            idx++;
        }


        lqn.graph.set(lqn.nidx, lqn.nidx, 0);

        Map<Integer, Task> tasks = this.tasks;
        int cidx = 0;

        lqn.calltype = new HashMap<>();
        lqn.iscaller = new Matrix(lqn.nidx + 1, lqn.nidx + 1, (lqn.ntasks + lqn.nacts) * (lqn.ntasks + lqn.nentries));
        lqn.issynccaller = new Matrix(lqn.nidx + 1, lqn.nidx + 1, (lqn.ntasks + lqn.nacts) * (lqn.ntasks + lqn.nentries));
        lqn.isasynccaller = new Matrix(lqn.nidx + 1, lqn.nidx + 1, (lqn.ntasks + lqn.nacts) * (lqn.ntasks + lqn.nentries));
        lqn.callpair = new Matrix(lqn.nacts * lqn.nentries + 1, 3, lqn.nacts * lqn.nentries);
        lqn.callproc = new HashMap<>();
        lqn.callnames = new HashMap<>();
        lqn.callhashnames = new HashMap<>();
        lqn.taskgraph = new Matrix(lqn.ntasks + lqn.tshift + 1, lqn.ntasks + lqn.tshift + 1, (lqn.ntasks + lqn.tshift) * (lqn.ntasks + lqn.tshift));
        lqn.actpretype = new Matrix(1, lqn.nidx + 1, lqn.nacts);
        lqn.actposttype = new Matrix(1, lqn.nidx + 1, lqn.nacts);

        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t + 1;

            for (int a = 0; a < tasks.get(t).activities.size(); a++) {
                int aidx = Utils.findString(lqn.hashnames, "A:" + tasks.get(t).activities.get(a).getName());
                lqn.callsof.put(aidx, new ArrayList<>());

                String boundToEntry = tasks.get(t).activities.get(a).boundToEntry;
                int eidx = Utils.findString(lqn.hashnames, "E:" + boundToEntry);
                if (eidx < 0) {
                    eidx = Utils.findString(lqn.hashnames, "I:" + boundToEntry);
                }
                if (eidx > 0) {
                    lqn.graph.set(eidx, aidx, 1);
                }

                for (int s = 0; s < tasks.get(t).activities.get(a).syncCallDests.size(); s++) {

                    int target_eidx = Utils.findString(lqn.hashnames, "E:" + tasks.get(t).activities.get(a).syncCallDests.get(s));

                    if (target_eidx <= 0) {
                        target_eidx = Utils.findString(lqn.hashnames, "I:" + tasks.get(t).activities.get(a).syncCallDests.get(s));
                    }
                    int target_tidx = (int) lqn.parent.get(target_eidx);
                    cidx++;

                    lqn.calltype.put(cidx, CallType.SYNC);
                    lqn.callpair.set(cidx, 1, aidx);
                    lqn.callpair.set(cidx, 2, target_eidx);
                    lqn.callnames.put(cidx, lqn.names.get(aidx) + "=>" + lqn.names.get(target_eidx));
                    lqn.callhashnames.put(cidx, lqn.hashnames.get(aidx) + "=>" + lqn.hashnames.get(target_eidx));
                    lqn.callproc.put(cidx, new Geometric(1 / tasks.get(t).activities.get(a).syncCallMeans.get(s)));

                    lqn.callsof.get(aidx).add(cidx);
                    lqn.iscaller.set(aidx, target_tidx, 1);
                    lqn.iscaller.set(aidx, target_eidx, 1);
                    lqn.iscaller.set(tidx, target_tidx, 1);//1 -> true
                    lqn.iscaller.set(tidx, target_eidx, 1);
                    lqn.issynccaller.set(tidx, target_tidx, 1);
                    lqn.issynccaller.set(tidx, target_eidx, 1);
                    lqn.issynccaller.set(aidx, target_eidx, 1);
                    lqn.issynccaller.set(aidx, target_tidx, 1);
                    lqn.taskgraph.set(tidx, target_tidx, 1);
                    lqn.graph.set(aidx, target_eidx, 1);
                }

                for (int s = 0; s < tasks.get(t).activities.get(a).asyncCallDests.size(); s++) {
                    int target_eidx = Utils.findString(lqn.hashnames, "E:" + tasks.get(t).activities.get(a).asyncCallDests.get(s));
                    int target_tidx = (int) lqn.parent.get(target_eidx);
                    cidx++;

                    lqn.calltype.put(cidx, CallType.ASYNC);
                    lqn.callpair.set(cidx, 1, aidx);
                    lqn.callpair.set(cidx, 2, target_eidx);
                    lqn.callnames.put(cidx, lqn.names.get(aidx) + "->" + lqn.names.get(target_eidx));
                    lqn.callhashnames.put(cidx, lqn.hashnames.get(aidx) + "->" + lqn.hashnames.get(target_eidx));
                    lqn.callproc.put(cidx, new Geometric(1 / tasks.get(t).activities.get(a).syncCallMeans.get(s)));
                    lqn.callsof.get(aidx).add(cidx);
                    lqn.iscaller.set(tidx, target_tidx, 1);//1 -> true
                    lqn.iscaller.set(tidx, target_eidx, 1);
                    lqn.iscaller.set(aidx, target_eidx, 1);
                    lqn.iscaller.set(aidx, target_tidx, 1);
                    lqn.isasynccaller.set(tidx, target_tidx, 1);
                    lqn.isasynccaller.set(tidx, target_eidx, 1);
                    lqn.isasynccaller.set(aidx, target_eidx, 1);
                    lqn.isasynccaller.set(aidx, target_tidx, 1);
                    lqn.taskgraph.set(tidx, target_tidx, 1);
                    lqn.graph.set(aidx, target_eidx, 1);
                }

                for (int ap = 0; ap < tasks.get(t).precedences.size(); ap++) {
                    String pretype = tasks.get(t).precedences.get(ap).preType;
                    String posttype = tasks.get(t).precedences.get(ap).postType;
                    List<String> preacts = tasks.get(t).precedences.get(ap).preActs;
                    List<String> postacts = tasks.get(t).precedences.get(ap).postActs;
                    for (int prea = 0; prea < preacts.size(); prea++) {
                        int preaidx = Utils.findString(lqn.hashnames, "A:" + tasks.get(t).precedences.get(ap).preActs.get(prea));
                        double preParam = 1;
                        if (pretype.equals(PRE_AND)) {
                            Matrix quorum = tasks.get(t).precedences.get(ap).preParams;
                            if (!quorum.isEmpty()) {
                                preParam = quorum.get(0) / preacts.size();
                            }
                        }

                        switch (posttype) {
                            case POST_OR:
                                for (int posta = 0; posta < postacts.size(); posta++) {
                                    int postaidx = Utils.findString(lqn.hashnames, "A:" + tasks.get(t).precedences.get(ap).postActs.get(posta));
                                    Matrix probs = tasks.get(t).precedences.get(ap).postParams;
                                    double postParam = probs.get(posta);
                                    lqn.graph.set(preaidx, postaidx, preParam * postParam);
                                    lqn.actpretype.set(0, preaidx, ActivityPrecedence.getPrecedenceId(tasks.get(t).precedences.get(ap).preType));//TODO:String to ASCII array
                                    lqn.actposttype.set(0, postaidx, ActivityPrecedence.getPrecedenceId(tasks.get(t).precedences.get(ap).postType));
                                }
                            case POST_LOOP:
                                Matrix counts = tasks.get(t).precedences.get(ap).postParams;
                                int enda = postacts.size();
                                int endaidx = Utils.findString(lqn.hashnames, "A:" + tasks.get(t).precedences.get(ap).postActs.get(enda));
                                for (int posta = 0; posta < postacts.size() - 1; posta++) {
                                    int postaidx = Utils.findString(lqn.hashnames, "A:" + tasks.get(t).precedences.get(ap).postActs.get(posta));
                                    double postParam = 1.0 / (postacts.size() - 1);
                                    lqn.graph.set(preaidx, postaidx, preParam * postParam);
                                    lqn.graph.set(postaidx, postaidx, 1 - 1.0 / counts.length());
                                    lqn.graph.set(postaidx, endaidx, 1.0 / counts.length());
                                    lqn.actposttype.set(0, endaidx, ActivityPrecedence.getPrecedenceId(tasks.get(t).precedences.get(ap).postType));
                                }
                            default:
                                for (int posta = 0; posta < postacts.size(); posta++) {
                                    int postaidx = Utils.findString(lqn.hashnames, "A:" + tasks.get(t).precedences.get(ap).postActs.get(posta));
                                    double postParam = 1;
                                    lqn.graph.set(preaidx, postaidx, preParam * postParam);
                                    lqn.actpretype.set(0, preaidx, ActivityPrecedence.getPrecedenceId(tasks.get(t).precedences.get(ap).preType));
                                    lqn.actposttype.set(0, postaidx, ActivityPrecedence.getPrecedenceId(tasks.get(t).precedences.get(ap).postType));
                                }
                        }
                    }
                }
            }

        }

        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t + 1;
            for (int aidx : lqn.actsof.get(tidx)) {
                List<Integer> postaidxs = new ArrayList<>();
                for (int col = 0; col < lqn.graph.numCols; col++) {
                    if (lqn.graph.get(aidx, col) != 0) {
                        postaidxs.add(col);
                    }
                }

                boolean isreply = true;
                for (int postaidx : postaidxs) {
                    for (int acts : lqn.actsof.get(tidx)) {
                        if (acts == postaidx) {
                            isreply = false;
                            break;
                        }
                    }
                }
                if (isreply) {

                    int parentidx = aidx;
                    while (lqn.type.get(parentidx) != LayeredNetworkElement.ENTRY) {
                        List<Integer> ancestors = new ArrayList<>();
                        for (int row = 0; row < lqn.graph.numRows; row++) {
                            if (lqn.graph.get(row, parentidx) != 0) {
                                ancestors.add(row);
                            }
                        }
                        parentidx = ancestors.get(0);
                    }
                    if (lqn.type.get(parentidx) == LayeredNetworkElement.ENTRY) {
                        lqn.replygraph.set(aidx, parentidx, 1);
                    }
                }
            }
        }
        lqn.ncalls = lqn.calltype.size();

        List<Integer> tidxs = new ArrayList<>();
        for (int i = 0; i < lqn.schedid.length(); i++) {
            if (((int) lqn.schedid.get(i)) == 0) {
                tidxs.add(i);
            }
        }

        for (int tidx : tidxs) {
            if (lqn.type.get(tidx) == LayeredNetworkElement.TASK) {
                List<Integer> callers = new ArrayList<>();
                for (int row = 0; row < lqn.taskgraph.numRows; row++) {
                    if (lqn.taskgraph.get(row, tidx) != 0) {
                        callers.add(row);
                    }
                }
                List<Integer> callers_inf = new ArrayList<>();
                for (int caller : callers) {
                    if (lqn.mult.get(caller) < 0) {
                        callers_inf.add(1);
                    } else {
                        callers_inf.add(0);
                    }
                }
            }
        }
        lqn.isref = new Matrix(1, lqn.nhosts + lqn.ntasks + 1, lqn.ntasks);
        for (int col = 1; col < lqn.schedid.numCols; col++) {
            if (lqn.schedid.get(0, col) == SchedStrategy.toID(SchedStrategy.REF)) {
                lqn.isref.set(0, col, 1);

            }
        }

        lqn.iscache = new Matrix(1, lqn.nhosts + lqn.ntasks + 1, lqn.ntasks);
        for (int i = 1; i < lqn.nitems.length(); i++) {
            if (lqn.nitems.get(i) > 0) {
                lqn.iscache.set(0, i, 1);
            }
        }

        // TODO:
        // lqn.refset = zeros(lqn.nidx,1);
        //  [conncomps, roots]=graph_connected_components(lqn.taskgraph(lqn.nhosts+1:end, lqn.nhosts+1:end));
        // lqn.conntasks = conncomps;
        // for r=1:length(roots)
        // lqn.conntasks(find(lqn.conntasks == r)) = lqn.tshift+roots(r);
        // end

        return lqn;
    }


    public int getNumberOfLayers() {
        return getNumberOfModels();
    }

    public int getNumberOfModels() {
        if (this.ensemble.isEmpty()) {
            getEnsemble();
        }
        return this.ensemble.size();
    }

    @Override
    public List<Network> getEnsemble() {
        if (this.ensemble.isEmpty()) {
            SolverLN solver = new SolverLN(this);
            this.ensemble = solver.getEnsemble();
        }
        return this.ensemble;
    }

    public List<Network> getLayers() {
        return getEnsemble();
    }

    public void summary() {
        LayeredNetworkStruct this_lqn = getStruct();
        this_lqn.print();
    }


    public static LayeredNetwork parseXML(String filename, int verbose) {

        LayeredNetwork myLN = new LayeredNetwork(filename.replace("_", "\\_"));

        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder dBuilder = null;
        try {
            dBuilder = dbFactory.newDocumentBuilder();
        } catch (ParserConfigurationException e) {
            e.printStackTrace();
        }

        Document doc = null;
        try {
            doc = dBuilder.parse(filename);
        } catch (SAXException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        doc.getDocumentElement().normalize();

        if (verbose > 0) {
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
            String scheduling = procElement.getAttribute("scheduling");//TODO:to ScheduleStrat
            double multiplicity = Double.parseDouble(procElement.getAttribute("multiplicity"));
            double replication = Double.parseDouble(procElement.getAttribute("replication"));

//           if(String.valueOf(scheduling).equals("inf")){//TODO inf/finite/nan in java
//               //if isfinite(multiplicity) //TODO: inf in java?
//               // line_warning(mfilename,'A finite multiplicity is specified for a host processor with inf scheduling. Remove it or set it to inf.');
//               //multiplicity = Inf
//           }else if(isNAN(multiplicity)){
//               multiplicity = 1;
//           }

            double quantum = Double.parseDouble(procElement.getAttribute("quantum"));
            double speedFactor = Double.parseDouble(procElement.getAttribute("speedFactor"));

            Processor newProc = new Processor(myLN, name, (int) multiplicity, SchedStrategy.valueOf(scheduling), quantum, speedFactor);
            newProc.setReplication((int) replication);
            procObj.put(procObj.size(), newProc);

            NodeList taskList = procElement.getElementsByTagName("task");

            for (int j = 0; j < taskList.getLength(); j++) {
                Element taskElement = (Element) taskList.item(j);
                String tName = taskElement.getAttribute("name");
                String tScheduling = taskElement.getAttribute("scheduling");
                double tReplication = Double.parseDouble(taskElement.getAttribute("replication"));
                double tMultiplicity = Double.parseDouble(taskElement.getAttribute("multiplicity"));
                double tThinkTimeMean = Double.parseDouble(taskElement.getAttribute("think_time"));

                Task newTask = new Task(myLN, tName, (int) tMultiplicity, SchedStrategy.valueOf(tScheduling), new Immediate());
                //todo: if tThinkTimeMean>0 thinkTime = Exp.fitMean(thinkTimeMean);
                newTask.setReplication((int) replication);
                taskObj.put(taskObj.size(), newTask);

                NodeList entryList = taskElement.getElementsByTagName("entry");
                for (int k = 0; k < entryList.getLength(); k++) {
                    Element entryElement = (Element) entryList.item(k);
                    String eName = entryElement.getAttribute("name");
                    Entry newEntry = new Entry(myLN, eName);
                    double openArrivalRate = Double.parseDouble(entryElement.getAttribute("open-arrival-rate"));
                    entryObj.put(entryObj.size(), newEntry);

                    NodeList entryPhaseActsList = entryElement.getElementsByTagName("entry-phase-activities");
                    if (entryPhaseActsList.getLength() > 0) {
                        Element entryPhaseActsElement = (Element) entryPhaseActsList.item(0);
                        NodeList actList = entryPhaseActsElement.getElementsByTagName("activity");
                        Map<Integer, String> nameList = new HashMap<>();
                        for (int l = 0; l < actList.getLength(); l++) {
                            Element actElement = (Element) actList.item(l);
                            double phase = Double.parseDouble(actElement.getAttribute("phase"));
                            nameList.put((int) phase, actElement.getAttribute("name"));
                            double hostDemandMean = Double.parseDouble(actElement.getAttribute("host-demand-mean"));
                            double hostDemandSCV = Double.parseDouble(actElement.getAttribute("host-demand-cvsq"));
                            Distribution hostDemand = new Immediate();
                            if (hostDemandMean <= 0) {
                                hostDemand = new Immediate();//todo
                            } else {
//                              if(hostDemandSCV<=0){
//                                  hostDemand = Det(hostDemandMean);
//                              }else if(hostDemandSCV<1){
//                                  hostDemand = APH.fitMeanAndSCV(hostDemandMean, hostDemandSCV);
//                              }else if(((int) hostDemand) ==1){
//                                  hostDemand = Exp.fitMean(hostDemandMean);
//                              }else{
//                                  hostDemand = HyperExp.fitMeanAndSCV(hostDemandMean,hostDemandSCV);
//                              }
                            }
                            String boundToEntry;
                            if (phase == 1) {
                                boundToEntry = newEntry.getName();
                            } else {
                                boundToEntry = "";
                            }

                            String callOrder = actElement.getAttribute("call-order");
                            Activity newAct = new Activity(myLN, nameList.get(phase), hostDemand, boundToEntry, callOrder);
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

                        if (!nameList.isEmpty()) {
                            newEntry.replyActivity.put(1, nameList.get(1));
                        }

                        Map<String, List<Integer>> tempMap = new HashMap<>();
                        List<Integer> tempList = new ArrayList<>();
                        tempList.add(taskID);
                        tempList.add(procID);
                        tempMap.put(newTask.getName(), tempList);
                        entries.put(entries.size(), tempMap);
                        newTask.addEntry(newEntry);
                        newEntry.parent = newTask;
                        entryID++;
                    }

                    NodeList taskActsList = taskElement.getElementsByTagName("task-activities");
                    if (taskActsList.getLength() > 0) {
                        Element taskActsElement = (Element) taskActsList.item(0);
                        NodeList actList = taskActsElement.getElementsByTagName("activity");
                        for (int l = 0; l < actList.getLength() - 1; l++) {
                            Element actElement = (Element) actList.item(l);
                            if (actElement.getParentNode().getNodeName().equals("task-activities")) {
                                String actName = actElement.getAttribute("name");
                                double hostDemandMean = Double.parseDouble(actElement.getAttribute("host-demand-mean"));
                                double hostDemandSCV = Double.parseDouble(actElement.getAttribute("host-demand-cvsq"));
                                Distribution hostDemand = new Immediate();
                                //if(isnan(hostDemandSCV)) hostDemandSCV = 1;
                                if (hostDemandMean <= 0) {
                                    hostDemand = new Immediate();//todo
                                } else {
//                                  if(hostDemandSCV<=0){
//                                      hostDemand = Det(hostDemandMean);
//                                  }else if(hostDemandSCV<1){
//                                      hostDemand = APH.fitMeanAndSCV(hostDemandMean, hostDemandSCV);
//                                  }else if(((int) hostDemand) ==1){
//                                      hostDemand = Exp.fitMean(hostDemandMean);
//                                  }else{
//                                      hostDemand = HyperExp.fitMeanAndSCV(hostDemandMean,hostDemandSCV);
//                                  }
                                }
                                String boundToEntry = actElement.getAttribute("bound-to-entry");
                                String callOrder = actElement.getAttribute("call-order");
                                Activity newAct = new Activity(myLN, actName, hostDemand, boundToEntry, callOrder);
                                actObj.put(actObj.size(), newAct);

                                NodeList synchCalls = actElement.getElementsByTagName("synch-call");
                                for (int m = 0; m < synchCalls.getLength() - 1; m++) {
                                    Element callElement = (Element) synchCalls.item(m);
                                    String dest = callElement.getAttribute("dest");
                                    double mean = Double.parseDouble(callElement.getAttribute("calls-mean"));
                                    newAct.synchCall(dest, mean);
                                }

                                NodeList asynchCalls = actElement.getElementsByTagName("asynch-call");
                                for (int m = 0; m < asynchCalls.getLength() - 1; m++) {
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
                        for (int l = 0; l < precList.getLength() - 1; l++) {
                            Element precElement = (Element) precList.item(l);

                            Map<Integer, String> preTypes = new HashMap<>();
                            preTypes.put(0, PRE_SEQ);
                            preTypes.put(1, PRE_AND);
                            preTypes.put(2, PRE_OR);
                            NodeList preList = null;
                            String preType = "";
                            for (int m = 0; m < preTypes.size(); m++) {
                                preType = preTypes.get(m);
                                preList = precElement.getElementsByTagName(preType);
                                if (preList.getLength() > 0) {
                                    break;
                                }
                            }

                            Element preElement = (Element) preList.item(0);
                            Matrix preParams = new Matrix(0, 0, 0);
                            NodeList preActList = preElement.getElementsByTagName("activity");
                            Map<Integer, String> preActs = new HashMap<>();
                            if (preType.equals(PRE_OR)) {
                                for (int m = 0; m < preActList.getLength() - 1; m++) {
                                    Element preActElement = (Element) preActList.item(m);
                                    preActs.put(m + 1, preActElement.getAttribute("name"));
                                    preParams.set(0, m + 1, Double.parseDouble(preActElement.getAttribute("prob")));
                                }
                            } else if (preType.equals(PRE_AND)) {
                                for (int m = 0; m < preActList.getLength() - 1; m++) {
                                    Element preActElement = (Element) preActList.item(m);
                                    preActs.put(m + 1, preActElement.getAttribute("name"));
                                }
                            } else {
                                Element preActElement = (Element) preActList.item(0);
                                preActs.put(1, preActElement.getAttribute("name"));
                            }
                            if (preParams.isEmpty()) preParams = null;

                            Map<Integer, String> postTypes = new HashMap<>();
                            postTypes.put(0, POST_SEQ);
                            postTypes.put(1, POST_AND);
                            postTypes.put(2, POST_OR);
                            postTypes.put(3, POST_LOOP);
                            NodeList postList = null;
                            String postType = null;
                            for (int m = 0; m < postTypes.size(); m++) {
                                postType = postTypes.get(m);
                                postList = precElement.getElementsByTagName(postType);
                                if (postList.getLength() > 0) break;
                            }
                            Element postElement = (Element) postList.item(0);
                            NodeList postActList = postElement.getElementsByTagName("activity");
                            Map<Integer, String> postActs = new HashMap<>();
                            Matrix postParams = new Matrix(0, 0, 0);
                            if (postType.equals(POST_OR)) {
                                for (int m = 0; m < postActList.getLength() - 1; m++) {
                                    Element postActElement = (Element) postActList.item(m);
                                    postActs.put(m + 1, postActElement.getAttribute("name"));
                                    preParams.set(0, m + 1, Double.parseDouble(postActElement.getAttribute("prob")));
                                }
                            } else if (postType.equals(POST_LOOP)) {
                                for (int m = 0; m < postActList.getLength() - 1; m++) {
                                    Element postActElement = (Element) postActList.item(m);
                                    postActs.put(m + 1, postActElement.getAttribute("name"));
                                    postParams.set(0, m + 1, Double.parseDouble(postActElement.getAttribute("count")));
                                }
                                postActs.put(postActs.size(), postElement.getAttribute("end"));
                            } else {
                                for (int m = 0; m < postActList.getLength() - 1; m++) {
                                    Element postActElement = (Element) postActList.item(m);
                                    postActs.put(m + 1, postActElement.getAttribute("name"));
                                }
                            }
                            //ActivityPrecedence newPrec = new ActivityPrecedence(preActs, postActs, preType, postType, preParams, postParams);
                            //newTask.addPrecedence(newPrec);
                        }
                        NodeList replyList = taskActsElement.getElementsByTagName("reply-entry");
                        for (int l = 0; l < replyList.getLength() - 1; l++) {
                            Element replyElement = (Element) replyList.item(l);
                            String replyName = replyElement.getAttribute("name");
                        }
                    }
                }
            }
        }
        return myLN;
    }

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
}
